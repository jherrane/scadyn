module shapebeam
use T_matrix
use bessel
use gaussquad

implicit none
contains

!******************************************************************************
! The routine computes the SVWF expansion coefficients 
! for a time harmonic x-polarized planewave
! propagating +z-direction  with the wave number k
subroutine pw(Nmax, k, a_nm, b_nm)
complex(dp), dimension((Nmax+1)**2-1) :: a_nm, b_nm
real(dp) :: k, C
integer :: Nmax, ind, n, m

ind = 0

a_nm = dcmplx(0.0,0.0)
b_nm = dcmplx(0.0,0.0)

do n = 1,Nmax
   do m = -n, n 
      ind = ind + 1
      
      C = sqrt(pi*(2*n+1d0))

      if(abs(m) == 1) then 
         a_nm(ind) = sign(1,m)*dcmplx(0.0,1.0)**(n+1.0) * C
         b_nm(ind) = dcmplx(0.0,1.0)**(n+1.0) * C
      end if
   
   end do
end do

end subroutine pw

!******************************************************************************
! The routine computes the SVWF expansion coefficients 
! for a time harmonic x-and y-polarized planewave
! propagating +z-direction  with the wave number k
subroutine pw2(Nmax, k, a_nm, b_nm, a_nm2, b_nm2)
complex(dp), dimension((Nmax+1)**2-1) :: a_nm, b_nm, a_nm2, b_nm2
integer :: Nmax, n, las, nm_in
real(dp) :: k

complex(dp), dimension(:), allocatable :: rotD
integer, dimension(:,:), allocatable :: indD

call pw(Nmax, k, a_nm, b_nm)

las = 0
do n = 1,Nmax
   las = las + (2*n+1)**2
end do

allocate(rotD(las))
allocate(indD(las,2))

call sph_rotation_sparse(0.0d0, -pi/2.0d0, Nmax, rotD, indD)
a_nm2 = sparse_matmul(rotD,indD,a_nm,nm_in)
b_nm2 = sparse_matmul(rotD,indD,b_nm,nm_in)

end subroutine pw2

!******************************************************************************
! The routine computes the SVWF expansion coefficients 
! for a time harmonic spherically polarized planewave (x+iy)
! propagating +z-direction  with the wave number k (Moreira et al 2016)
subroutine pw_sph(Nmax, k, a_nm, b_nm)
complex(dp), dimension((Nmax+1)**2-1) :: a_nm, b_nm
real(dp) :: k, C
integer :: Nmax, ind, n, m

ind = 0

a_nm = dcmplx(0.0,0.0)
b_nm = dcmplx(0.0,0.0)

do n = 1,Nmax
   do m = -n, n 
      ind = ind + 1
      
      C = sqrt(pi*(2*n+1d0))

      if(m == 1) then 
         a_nm(ind) = dcmplx(0.0,1.0)**n * C
         b_nm(ind) = -dcmplx(0.0,1.0)**(n+1.0) * C
      end if
   
   end do
end do

end subroutine pw_sph

!******************************************************************************
! The Gaussian amplitude profile in localized approximation (kw0>=5) a'la
! MSTM 3.0 (Mackowski et al 2013)
subroutine gaussian_beams(matrices,mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: width
integer :: i 

width = 5d0/(maxval(mesh%ki))/sqrt(2d0)

do i = 1,matrices%bars
    call gaussian_beam_shape(matrices,mesh,i,matrices%Nmaxs(i),width)
end do

end subroutine gaussian_beams

!******************************************************************************
! The Gaussian amplitude profile in localized approximation (kw0>=5) a'la
! MSTM 3.0 (Mackowski et al 2013)
subroutine gaussian_beam_shape(matrices, mesh, i, Nmax, width)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: C, gn, kw0, width
integer :: n, m, ind, i, Nmax

kw0 = mesh%ki(i)*width

if(kw0<5d0) then
    write(*,'(A)'), "    Problem: width of Gaussian beam at focal point is smaller than wavelength!"
    write(*,'(2(A, ES9.3))'), "     Wavelength is ", 2d0*pi/mesh%ki(i), ", width is ", width
end if

ind = 0
do n = 1,Nmax
    gn = dexp(-((dble(n)+.5d0)/kw0)**2)
    do m = -n,n
        ind = ind + 1
        matrices%as(ind,i) = gn*matrices%as(ind,i)
        matrices%bs(ind,i) = gn*matrices%bs(ind,i)
    end do 
end do

end subroutine gaussian_beam_shape

!******************************************************************************

subroutine laguerre_gaussian_beams(matrices, mesh, p, l)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: width
integer :: i, p, l

width = 0.3d0/(maxval(mesh%ki))

do i = 1,matrices%bars
    call laguerre_gauss_num(matrices, mesh, i, p, l, width)
end do

end subroutine laguerre_gaussian_beams

!******************************************************************************
! Axisymmetric LG-beam. In here, instead of some other routines.
subroutine laguerre_gauss_num(matrices, mesh, whichWL, n, m, w0)
type (data) :: matrices
type (mesh_struct) :: mesh
integer :: whichWL, nmax, i, j, m, n, ind
real(dp) :: w0, r(3), theta, phi, x, f, n_j, norm
complex(dp), dimension(:), allocatable :: a_jm, b_jm
complex(dp) :: Enm, dEnm, Yjm
real(dp), dimension(:,:), allocatable :: PP
real(dp), dimension(:), allocatable :: w
real(dp), dimension(:), allocatable :: Pjm
CHARACTER(LEN=80) :: fname = 'points.h5'

call integ_points(PP,w,20,20)
! call sample_points(PP,w,20,20)
! call write2file(dcmplx(PP),fname)

f = 1d0/(mesh%ki(whichWL)*w0)
nmax = matrices%Nmaxs(whichWL)

allocate(a_jm((nmax+1)**2-1), b_jm((nmax+1)**2-1))
a_jm = dcmplx(0d0)
b_jm = dcmplx(0d0)

ind = 0
do j = 1,nmax
   if(allocated(Pjm)) then 
      deallocate(Pjm)
   end if
   allocate(Pjm(j+1))

   ind = j*(j+1)+m

   ! Theta-phi-integration is done using Gaussian quadratures in the i-loop
   do i = 1, size(PP,2)
      theta = PP(2,i)
      phi = PP(3,i)

      ! Calculate the normalized spherical harmonics
      call legendre2(j,cos(theta),Pjm) ! The same as legendre of MATLAB
      if(m+1>j)then
         Yjm=dcmplx(0d0)
      else
         Yjm = Pjm(m+1)*exp(i1*m*phi)* &
         sqrt((2d0*j+1)/(4d0*pi)*factorial(j-m)/factorial(j+m))
      end if
      n_j = 1/(sqrt(dble(j*(j+1))))

      x = sin(theta)/f/sqrt(2d0)
      Enm = E_nm(x,f,n,m)
      dEnm = dE_nm(theta,f,n,m)

      a_jm(ind) = a_jm(ind) + 2d0*n_j*i1**j*dconjg(Yjm)*exp(i1*m*phi)* &
         ((2d0*sin(phi)*sin(theta)**2 -(i1*m)*cos(phi))*Enm &
         -sin(phi)*sin(theta)*cos(theta)*dEnm)*w(i)

      b_jm(ind) = b_jm(ind) - 2d0*n_j*i1**j *dconjg(Yjm)*exp(i1*m*phi)* &
         (sin(theta)*cos(phi)*dEnm - i1*m*cos(theta)*sin(phi)*Enm)*w(i)
   end do
end do

! Normalization so that sum_jm(a_jm + b_jm) = 1
norm = 0d0
do i = 1,size(a_jm)
 norm = norm + abs(a_jm(i)) + abs(b_jm(i))
end do

matrices%as(1:(nmax+1)**2-1,whichWL) = a_jm/norm
matrices%bs(1:(nmax+1)**2-1,whichWL) = b_jm/norm

end subroutine laguerre_gauss_num

!******************************************************************************

subroutine gaussi(P, w, order)
real(dp), dimension(:), allocatable :: P
real(dp), dimension(:), allocatable :: w
integer, intent(in) :: order

allocate(P(order), w(order))

call cpquad(order,dble(0),"Legendre",w,P)

P = P/2.0d0+0.5d0
w = w/2.0d0

end subroutine gaussi

!******************************************************************************

subroutine integ_points(P,w,M,N)
real(dp), dimension(:,:), allocatable :: P
real(dp), dimension(:), allocatable :: w
integer :: M, N
real(dp), dimension(:), allocatable :: P_theta, w_theta 
integer :: i1, i2, j1
real(dp) :: phi
call gaussi(P_theta, w_theta, M)

allocate(P(3,M*N),w(M*N))
P_theta = pi * P_theta
w_theta = pi * w_theta

j1 = 1
do i1 = 1,N
   phi = dble(i1-1) * 2*pi / dble(N)
   do i2 = 1,M
      P(1,j1) = cos(phi) * sin(P_theta(i2))
      P(2,j1) = sin(phi) * sin(P_theta(i2))
      P(3,j1) = cos(P_theta(i2))
      w(j1) = W_theta(i2) / dble(N) * sin(P_theta(i2)) * 2 * pi
      
      j1 = j1 + 1
   end do
end do


end subroutine integ_points

!******************************************************************************

function E_nm(x, f, n, m) result(Enm)
integer :: n, m
real(dp) :: x, f
real(dp) :: L
complex(dp) :: Enm

L = L_pl(n, m, x**2)
if(n<0) then
   Enm = 0
else
   Enm = (x**m/(i1**(2*n+m+1)*2*f**2))*L*exp(-0.5d0*x**2)
end if

end function E_nm

!******************************************************************************

function dE_nm(theta, f, n, m) result(dEnm)
integer :: n, m
real(dp) :: theta, f, x
real(dp), dimension(1) :: LL
complex(dp) :: dEnm

x = sin(theta)/sqrt(2d0)/f

dEnm = ((m+x**2)/tan(theta))*E_nm(x,f,n,m) &
       - i1*(2*x**2/tan(theta))*E_nm(x,f,n-1,m+1)
if(theta < 1d-9) dEnm = dcmplx(0d0)

end function dE_nm

!******************************************************************************
! The generalized Laguerre polynomial as a summation, same as in MATLAB
function L_pl(p,l,x) result(Lpl)
real(dp) :: x, Lpl
integer :: n, m, p, l, j

Lpl = 1d0
Lpl = Lpl*choose(p+l,p)

do m = 1,p
  Lpl = Lpl + (-1)**m/factorial(m)*choose(p+l,p-m)*x**m
enddo

end function L_pl

!******************************************************************************
! Axisymmetric LG-beam
subroutine laguerre_gauss_farfield(matrices, mesh, i, p, l, w0)
type (data) :: matrices
type (mesh_struct) :: mesh
integer :: i, nmax, total_modes, iii, jjj, ntheta, nphi, p, l, tp, info, lwork, ind
integer, dimension(:), allocatable :: nn, mm, nn_old
complex(dp) :: x, y, BCP(9)
real(dp) :: k, w0, offset(3), zero_rejection_level, truncation_angle
real(dp), dimension(:), allocatable :: theta, phi, rw, LL, work
complex(dp), dimension(:), allocatable :: beam_envelope, Ex, Ey, &
Etheta, Ephi, e_field, expansion_coefficients, a, b, a_nm, b_nm
complex(dp), dimension(:,:), allocatable :: coefficient_matrix

k = mesh%ki(i)
nmax = matrices%Nmaxs(i)
allocate(a_nm((nmax+1)**2-1), b_nm((nmax+1)**2-1))
truncation_angle = 90
x = dcmplx(1d0,0d0)
y = dcmplx(0d0,0d0)

! total_modes = nmax**2 + 2*nmax
total_modes = 2*nmax
allocate(nn(total_modes), mm(total_modes), nn_old(total_modes))
do iii = 1,total_modes
    nn(iii) = ceiling(dble(iii)/2d0)
    mm(iii) = (-1d0)**(iii) + l
end do
nn_old = nn
nn = PACK(nn, nn_old >= abs(mm))
mm = PACK(mm, nn_old >= abs(mm))
ntheta = 2*(nmax+1)
nphi = 3
tp = ntheta*nphi
allocate(theta(tp), phi(tp))
call angular_grid(ntheta,nphi,theta,phi)
allocate(e_field(2*tp), rw(tp), LL(tp), beam_envelope(tp), Ex(tp), Ey(tp), &
    Etheta(tp), Ephi(tp))

do iii = 1,tp
 rw(i) = k**2*w0**2*(dtan(theta(iii)))**2/2
end do 

LL = laguerre(p, abs(l), rw)

do iii = 1,tp
 beam_envelope(i) = rw(iii)**(abs(l/2))*LL(iii)*zexp(-rw(iii)/2 + &
    dcmplx(0d0,1d0)*l*phi(iii))
 if(theta(iii)<pi*(180-truncation_angle)/180) beam_envelope(iii) = 0
end do 

! if(vlen(offset)>1d-9)
 ! Phase shift is exp(-i*k * offset.rhat)
 ! rhat = rtpv2xyzv( ones(size(theta)), zeros(size(theta)), zeros(size(theta)), ones(size(theta)), theta, phi );
 ! [offset,rhat] = matchsize(offset,rhat);
 ! phase_shift = exp( -i * k * dot(offset,rhat,2) );
 ! beam_envelope = beam_envelope .* phase_shift;
! endif

Ex = x*beam_envelope
Ey = y*beam_envelope
do iii = 1,tp
 Etheta(iii) = -Ex(i)*dcos(phi(iii)) - Ey(iii)*dsin(phi(iii))
 Ephi = -Ex(iii)*dsin(phi(iii)) + Ey(iii)*dcos(phi(iii))
end do

e_field = [Etheta,Ephi]
allocate(coefficient_matrix(2*tp,2*size(nn,1)), expansion_coefficients(2*size(nn,1)))

do iii = 1,size(nn,1)
 do jjj = 1,tp
  BCP = vsh(nn(iii), mm(iii), theta(jjj), phi(jjj) )
  coefficient_matrix(jjj,iii) = BCP(5) * dcmplx(0d0,1d0)**(nn(iii)+1)/dsqrt(dble(nn(iii))*(nn(iii)+1))
  coefficient_matrix(tp + jjj,iii) = BCP(6) * dcmplx(0d0,1d0)**(nn(iii)+1)/dsqrt(dble(nn(iii))*(nn(iii)+1))
  coefficient_matrix(jjj,iii+size(nn,1)) = BCP(2) * dcmplx(0d0,1d0)**(nn(iii))/dsqrt(dble(nn(iii))*(nn(iii)+1))
  coefficient_matrix(tp + jjj,iii+size(nn,1)) = BCP(3) * dcmplx(0d0,1d0)**(nn(iii))/dsqrt(dble(nn(iii))*(nn(iii)+1))
 end do
end do

lwork = 64*min(2*tp, 2*size(nn,1))
allocate(work(lwork))
call zgels('N', 2*tp, 2*size(nn,1), 1, coefficient_matrix, 2*tp, e_field, &
    2*tp, work, lwork, info)

expansion_coefficients = e_field(1:2*size(nn,1))
allocate(a(size(nn,1)), b(size(nn,1)))
a = expansion_coefficients(1:size(nn,1))
b = expansion_coefficients(size(nn,1)+1:2*size(nn,1))

do iii = 1,size(nn,1)
  ind = nn(iii)*(nn(iii)+1) + mm(iii)
  a_nm(ind) = a(iii)
  b_nm(ind) = b(iii)
end do

matrices%as(1:(nmax+1)**2-1,i) = a_nm
matrices%bs(1:(nmax+1)**2-1,i) = b_nm

end subroutine laguerre_gauss_farfield

!******************************************************************************

function LG_mode(p, l, r, phi) result(res)
integer :: p, l, fp, fpl
real(dp), dimension(:), allocatable :: r2
real(dp) :: r, phi
complex(dp) :: LG(1), res

fp = factorial(p)
fpl = factorial(p+abs(l))
allocate(r2(1))
r2(1) = 2*r**2
LG = sqrt(2d0*fp/(pi*fpl)) * &
(sqrt(2d0)*r)**abs(l)*laguerre(p,abs(l),r2)*exp(dcmplx(0,1)*l*phi) * &
 exp(-r**2) * exp(dcmplx(0,1) * (2*p + abs(l) + 1) * pi/2)

res = LG(1)

end function LG_mode

!******************************************************************************

subroutine fields_out(matrices,mesh,which,n)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: n,nn, i, j, ind, which
real(dp) :: lim
real(dp), allocatable :: z(:), y(:), grid(:,:)
complex(dp) :: F(3), G(3) 
complex(dp), dimension(:,:), allocatable :: E  

nn = n*n

grid = field_grid(matrices,mesh)

allocate(E(3,nn))
matrices%field_points = grid
do i = 1,nn
    call calc_fields(matrices%as(:,which),matrices%bs(:,which),&
    dcmplx(mesh%ki(which)), grid(:,i),F,G,0)
    E(:,i) = F
    ! E(:,i) = crossCC(F,G)/mu
end do
matrices%E_field = E

end subroutine fields_out
   
!******************************************************************************

subroutine scat_fields_out(matrices, mesh, which, n)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: n, nn, i, j, ind, which
real(dp) :: lim
real(dp), allocatable :: grid(:,:)
complex(dp), dimension(:) , allocatable :: p, q, p90, q90
complex(dp) :: F(3), G(3) 
complex(dp), dimension(:,:), allocatable :: E  

nn = n*n

grid = field_grid(matrices,mesh)

! Ensure that directions are ok. They might already be...
matrices%Rkt = eye(3)
allocate(E(3,nn))
call rot_setup(matrices)
call scattered_fields(matrices,1d0,p,q,p90,q90,which)

matrices%field_points = grid
do i = 1,nn
    call calc_fields(p,q,dcmplx(mesh%ki(which)), grid(:,i),F,G,1)
    E(:,i) = F
    ! E(:,i) = crossCC(F,G)/mu
end do

matrices%E_field = E


end subroutine scat_fields_out

!******************************************************************************

function field_grid(matrices,mesh) result(grid)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: n, nn, i, j, ind, which
real(dp) :: lim
real(dp), allocatable :: z(:), y(:), grid(:,:)
complex(dp) :: F(3), G(3) 
complex(dp), dimension(:,:), allocatable :: E  

lim = 3.5d0
n = 100
nn = n*n
allocate(z(n),y(n),grid(3,nn))

call linspace(-lim,lim,n,z)
call linspace(-lim,lim,n,y)
grid = 0d0
do i = 1,n
    do j = 1,n
        ind = n*(j-1)+i 
        grid(1,ind) = z(i)*mesh%a 
        grid(2,ind) = y(j)*mesh%a
    end do
end do

! do i = 1,size(grid,2)
!   write(*, '(3ES11.3)') grid(:,i)
! end do

end function field_grid

!******************************************************************************

subroutine write_fields(matrices,mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: width
integer :: i, n
character(LEN=80) :: gridname, fieldname, scatfield
i = 1
n = 100
gridname = 'grid_xyz.h5'
fieldname = 'E_field.h5'
scatfield = 'E_scat.h5'

call fields_out(matrices,mesh,i,n)
call write2file(dcmplx(matrices%field_points),gridname)
call write2file(matrices%E_field,fieldname)

call scat_fields_out(matrices,mesh,i,n)
call write2file(matrices%E_field,scatfield)

end subroutine write_fields

end module shapebeam