module shapebeam
use T_matrix

implicit none
contains

!*******************************************************************************
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

!*******************************************************************************
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

!*******************************************************************************
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

!*******************************************************************************
! The Gaussian amplitude profile in localized approximation (kw0>=5) a'la
! MSTM 3.0 (Mackowski et al 2013)
subroutine gaussian_beams(matrices,mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: width
integer :: i 

width = 5d0/(minval(mesh%ki))

do i = 1,matrices%bars
    call gaussian_beam_shape(matrices,mesh,i,matrices%Nmaxs(i),width)
end do

end subroutine gaussian_beams

!*******************************************************************************
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

!*******************************************************************************
! Calculate x-polarized bessel beams propagating along the z-axis.
subroutine bessel_beams(matrices,mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: i 

do i = 1,matrices%bars
    call bessel_beam_z(matrices,mesh,i,matrices%Nmaxs(i), 3, 3, 3, 3)
end do

end subroutine bessel_beams

!*******************************************************************************
! Routine for the VSWF expansion of a Bessel beam at origin
subroutine bessel_beam_z(matrices, mesh, i, Nmax, lmx, M, S, TM)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: M        ! 
integer :: S        ! 
integer :: lmax, lmx ! 
integer :: Nmax     ! 
integer :: TM       ! 
integer :: i        ! 

real(dp) :: k       ! The wavenumber
real(dp) :: kz      ! z-directed part of the wavenumber
real(dp) :: gamma   ! Transverse wavenumber
complex(dp), dimension(:), allocatable :: GTE  ! The final expansion, first half (a)
complex(dp), dimension(:), allocatable :: GTM  ! The final expansion, second half (b)
real(dp) :: rho, gr, cth, sth, cph, sph, zz
real(dp), dimension(:), allocatable :: Qlm, dQlm, Jn, dn 
integer, dimension(:), allocatable :: ll, mm
integer :: ls, le, mo, l
complex(dp) :: A0, psi

lmax = lmx*(lmx+2)+1
allocate(Qlm(lmax), dQlm(lmax), Jn(lmax), Dn(lmax))
allocate(ll(lmax), mm(lmax))
allocate(GTE(lmax),GTM(lmax))

ls = abs(M-lmax-1)
le = abs(M+lmax+1)
if(ls>le) then
    le = ls
    ls = 0
endif

call vswf_qlm(cth,lmax,Qlm,dQlm,ll,mm)
call bess_cyl(le,gr,Jn,Dn,Nmax)

do l = 1,lmax
    do m = -l,l
        mo = M-S*m
        A0 = 4d0*pi*I**(l-S*m)/sqrt(dble(l*(l+1)))
        if(mo<0) then
            psi = Jn(abs(mo))*(cph+S*I*sph)**(mo)*1*(-1)**abs(mo) !*exp(I*zz=0)
        else
            psi = Jn(abs(mo))*(cph+S*I*sph)**(mo)*1 !*exp(I*zz=0)
        end if

        if(TM == 1)then
            GTE(jlm(l,m)) =   A0*psi*m*Qlm(jlm(l,m))
            GTM(jlm(l,m)) = I*A0*psi*sth**2*dQlm(jlm(l,m))
        else
            GTE(jlm(l,m)) = I*A0*psi*sth**2*dQlm(jlm(l,m))
            GTM(jlm(l,m)) =  -A0*psi*m*Qlm(jlm(l,m))
        end if
    end do
end do

! Write result to matrices%a etc

end subroutine bessel_beam_z

!*******************************************************************************

subroutine vswf_qlm(x,lmax,Qlm,dQlm,lo,mo)
real(dp), intent(inout) :: x
real(dp), dimension(:), allocatable, intent(inout) :: Qlm, dQlm
integer, dimension(:), allocatable, intent(inout) :: lo, mo
integer :: l, m, lmax
real(dp) :: kl, cth, sth, cs2

cth = x
sth = sqrt(1d0-cth**2)
cs2 = sth**(-2)

! Qlm - First 4 terms
Qlm(jlm(0,0))=1/sqrt(4d0*pi)
Qlm(jlm(1, 1))=-gammaQ(1)*sth*Qlm(jlm(0,0))! Q11
Qlm(jlm(1, 0))=sqrt(3d0)*cth*Qlm(jlm(0,0))  ! Q10
Qlm(jlm(1,-1))=-Qlm(jlm(1,1))               ! Q11*(-1)

! dQlm - First 4 terms
dQlm(jlm(0, 0))=0.0
dQlm(jlm(1, 1))= -cth*Qlm(jlm(1, 1))*cs2
dQlm(jlm(1, 0))=-(cth*Qlm(jlm(1, 0))-sqrt(3d0)*Qlm(jlm(0,0)))*cs2
dQlm(jlm(1,-1))= -cth*Qlm(jlm(1,-1))*cs2

! lo - First 4 terms
lo(jlm(0, 0)) = 0
lo(jlm(1,-1)) = 1
lo(jlm(1, 0)) = 1
lo(jlm(1, 1)) = 1

! mo - First 4 terms
mo(jlm(0, 0)) =  0
mo(jlm(1,-1)) = -1
mo(jlm(1, 0)) =  0
mo(jlm(1, 1)) =  1

do l = 2,lmax
    ! Qlm positives extremes 
    Qlm(jlm(l, l  ))=-gammaQ(l)*sth*Qlm(jlm(l-1,l-1))
    Qlm(jlm(l, l-1))= deltaQ(l)*cth*Qlm(jlm(l-1,l-1))
    ! dQlm positives extremes 
    Qlm(jlm(l,-l  ))=(-1)**l*Qlm(jlm(l,l  ))
    Qlm(jlm(l,-l+1))=(-1)**(l-1)*Qlm(jlm(l,l-1))
                  
    kl=1d0/(sqrt(dble(l)*(l+1)))

    do m = l,0,-1
        lo(jlm(l, m)) =  l
        lo(jlm(l,-m)) =  l
        mo(jlm(l, m)) =  m
        mo(jlm(l,-m)) = -m

        if(m<l-1) then
            Qlm(jlm(l, m))=alfaQ(l,m)*cth*Qlm(jlm(l-1,m)) &
               -betaQ(l,m)*Qlm(jlm(l-2,m)) 
            Qlm(jlm(l,-m))=(-1)**m*Qlm(jlm(l, m))
        end if

        dQlm(jlm(l, m)) = -l*cth*cs2*Qlm(jlm(l, m))+ &
         sqrt((2d0*l+1.)/(2*l-1.))*sqrt(dble(l+m)*(l-m))*cs2*Qlm(jlm(l-1, m))

        dQlm(jlm(l,-m)) = -l*cth*cs2*Qlm(jlm(l,-m))+ &
         sqrt((2d0*l+1.)/(2*l-1.))*sqrt(dble(l+m)*(l-m))*cs2*Qlm(jlm(l-1,-m))
    end do
end do 

end subroutine vswf_qlm

!*******************************************************************************
! Array of cylindrical Bessel functions
subroutine bess_cyl(nmx, x, Jn, Dn, Nmax)
real(dp), intent(inout) :: x
real(dp), dimension(:), allocatable, intent(in) :: Jn, Dn 
integer, intent(inout) :: nmx
integer, intent(in) :: Nmax
integer :: l, m
real(dp) :: kl, cth, sth, cs2

! // Array of cylindrical Bessel func
! void bess_cyl(/* FUNCTION */
!       int *nmax, 
!       double *x, 
!       double *Jn, 
!       double *Dn, 
!       int *NMAX
!    ){
! //--------------------------------------
!    Jn[*nmax]=0.0;  // n-th Bessel func
!    Dn[*nmax]=1.0;  // Derivative of the n-th Bessel func
!    int dig=15;     // Number of digits of precision
!    int n,j;
!    bess_csv(nmax,&dig,x,&Jn[*nmax],&Dn[*nmax]);
!         //------------------------------------------------------------------------------
!         //Starting values for downward recurrence, Cylindrical Bessel func
!         //For calculation of J_n(x), one can calculate by downward recurrence from
!         //   J_{n-dig}(x), where dig is the number of precision in J_n(x).
!         void bess_csv(/* FUNCTION */
!               int *nmax,
!               int *dig,
!               double *x,
!               double *JN,
!               double *DN
!            ){
!         //--------------------------------------
!            double Jn[*dig], Dn[*dig];
!            Jn[*dig]=0.0;  // n-th Bessel func
!            Dn[*dig]=1.0;  // Derivative of the n-th Bessel func
!            int n,j;
!            for(n=*dig;n>0;n--){
!               // Downward Recursion
!               Jn[n-1]=lcfe_afs(*nmax+n  ,*x)*Jn[n  ]+Dn[n];
!               Dn[n-1]=lcfe_afs(*nmax+n-1,*x)*Jn[n-1]-Jn[n];
!               //Renormalization may be required
!               if(fabs(Jn[n-1])>1e+2){
!                  for(j=*dig;j>=n-1;j--){
!                     Dn[j]=Dn[j]/Jn[n-1];
!                     Jn[j]=Jn[j]/Jn[n-1];
!                  }
!               }
!            }
!            *JN=Jn[0];
!            *DN=Dn[0];
!         }
!    //lcfe_cbl(nmax,x,NMAX,&Dn[*nmax]);
!    for(n=*nmax;n>0;n--){
!       // Downward Recursion
!       Jn[n-1]=lcfe_afs(n  ,*x)*Jn[n  ]+Dn[n];
!       Dn[n-1]=lcfe_afs(n-1,*x)*Jn[n-1]-Jn[n];
!       //Renormalization may be required
!       if(fabs(Jn[n-1])>1e2){
!          for(j=*nmax;j>=n-1;j--){
!             Dn[j]=Dn[j]/Jn[n-1];
!             Jn[j]=Jn[j]/Jn[n-1];
!          }
!       }
!    }
!    // Normalization factor for the func
!    double pf=1.0;
!    // Using j0(x) and j1(x) from math.h,
!    // otherwise it should be entered from R
!    // or calculated internally.
!    double J0=j0(*x);
!    double J1=j1(*x);
!    if(fabs(J0)>1e-10){
!       //pf=*J0/Jn[0];
!       pf=j0(*x)/Jn[0];
!    }else{
!       //pf=*J1/Jn[1];
!       pf=j1(*x)/Jn[1];
!    }
!    // Normalization factor for the derivative
!    double pd=1.0;
!    if(fabs(J1)>1e-10){
!       pd=-(J1)/Dn[0];
!    }else{
!       pd=0.5*(Jn[0]-Jn[2])/Dn[1];
!    }
!    // Normalizations
!    for(n=0;n<=*nmax;n++){
!       Jn[n]=Jn[n]*pf;
!       Dn[n]=Dn[n]*pd;
!    }
! }

end subroutine bess_cyl


!*******************************************************************************

function alfaQ(l,m) result(a)
integer :: l, m
real(dp) :: a
a = sqrt(((2d0*l-1)*(2d0*l+1))/((l-m)*(l+m)))
end function alfaQ

!*******************************************************************************

function betaQ(l,m) result(b)
integer :: l, m
real(dp) :: b
b = sqrt((2d0*l+1)/(2*l-3))*sqrt((dble(l+m-1)*(l-m-1))/((l-m)*(l+m)))
end function betaQ

!*******************************************************************************

function gammaQ(l) result(g)
integer :: l
real(dp) :: g
g = sqrt((2d0*l+1)/(2d0*l))
end function gammaQ

!*******************************************************************************

function deltaQ(l) result(d)
integer :: l
real(dp) :: d
d = sqrt(2d0*l+1)
end function deltaQ

!*******************************************************************************

function jlm(l,m) result(j)
integer :: l, m, j
if(abs(m)>l) then
    j = 0
else
    j = l*(l+1)+m
endif
end function jlm

!*******************************************************************************

subroutine laguerre_gaussian_beams(matrices, mesh, p, l)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: width
integer :: i , p, l

width = 5d0/(minval(mesh%ki))

do i = 1,matrices%bars
    call laguerre_gauss_farfield(matrices, mesh, i, p, l, width)
end do

end subroutine laguerre_gaussian_beams

!*******************************************************************************
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

!*******************************************************************************

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

!*******************************************************************************

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
end do
matrices%E_field = E

end subroutine fields_out

!*******************************************************************************

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
end do

matrices%E_field = E


end subroutine scat_fields_out

!*******************************************************************************

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
        grid(3,ind) = z(i)*mesh%a 
        grid(2,ind) = y(j)*mesh%a
    end do
end do

end function field_grid

!*******************************************************************************

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