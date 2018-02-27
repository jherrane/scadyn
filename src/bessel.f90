module bessel
use T_matrix

implicit none
contains

!*******************************************************************************
! Routine for the VSWF expansion of a Bessel beam at origin
subroutine bessel_beam_z(matrices, mesh, i, lmx)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: MMM       ! Order of the Bessel Beam.
integer :: NNN
integer :: S         ! Chirality
integer :: lmax, lmx ! 
integer :: Nmax      ! 
integer :: TM        ! 1 if transverse magnetic mode is considered, 0 otherwise (?)
integer :: i         ! 

real(dp) :: k       ! The wavenumber
real(dp) :: kx, ky, kz      ! x,y,z-directed parts of the wavenumber
real(dp) :: gamma   ! Transverse wavenumber
complex(dp) :: GTE(0:lmx*(lmx+2)), GTE2(0:lmx*(lmx+2))  ! The final expansion, first half (a)
complex(dp) :: GTM(0:lmx*(lmx+2)), GTM2(0:lmx*(lmx+2))  ! The final expansion, second half (b)
complex(8), dimension(:), allocatable :: rotD
integer, dimension(:,:), allocatable :: indD
real(dp) :: rho, gr, cth, sth, cph, sph, zz, lambda, a, b, c, x, y, J0, J1, pf, pd
real(dp) :: Jn(0:lmx*(lmx+2)), Dn(0:lmx*(lmx+2))
real(dp) :: Qlm(0:lmx*(lmx+2)+1), dQlm(0:lmx*(lmx+2)+1)
real(dp) :: Pmn(0:lmx, 0:lmx), DPmn(0:lmx, 0:lmx)
integer, dimension(:), allocatable :: ll, mm
integer :: ls, le, mo, l, m, iii, las, nm
complex(dp) :: Wlsm, psi, i1, gte1, gtm1, T(3,3)

i1 = dcmplx(0d0,1d0)
MMM = 3
NNN = 5
S = -1
TM = 0
x = 1d-7
y = 1d-7

lmax = lmx*(lmx+2)

allocate(ll(lmax), mm(lmax))

! Transform matrix from circular to cartesian components
! T = reshape(0.5d0*dcmplx([sqrt(2d0),0d0,sqrt(2d0),-i1*sqrt(2d0),&
   ! 0d0,i1*sqrt(2d0),0d0,1d0,0d0]),[3,3])
T = eye(3)
T = transpose(T)

! Waveguide
lambda = 2d0*pi/mesh%ki(i)
a = 8d0*lambda
b = 8d0*lambda
c = 6d0*lambda
k = mesh%ki(i)
kx = MMM*pi/a
ky = NNN*pi/b 
gamma = sqrt(kx**2+ky**2)
kz = sqrt(k**2-gamma**2)

rho = sqrt(x**2+y**2)

cth = kz/k
gr = gamma*rho
cph = x/rho
sph = y/rho
sth = gamma/k
zz = kz*0d0

ls = abs(MMM-lmax-1)
le = abs(MMM+lmax+1)

if(ls>le) then
    le = ls
    ls = 0
endif

Jn = 0d0
Dn = 0d0

call vswf_qlm(cth,lmx,Qlm,dQlm)
call bess_cyl(lmax,gr,Jn,Dn)

Nmax = matrices%Nmaxs(i)
nm = (Nmax + 1)**2 - 1
las  =  ( Nmax  + 1 ) * ( 2 * Nmax + 1 ) * ( 2 * Nmax + 3 ) / 3 - 1
if(allocated(rotD)) deallocate(rotD,indD)
allocate(rotD(las))
allocate(indD(las,2))
call sph_rotation_sparse_gen2(eye(3), Nmax, rotD, indD)

do l = 1,lmx 
   do m = -l,l
      mo = MMM-S*m
      Wlsm = 4d0*pi*i1**(l-S*m)/sqrt(dble(l*(l+1)))
      if(mo<0) then
         psi = Jn(abs(mo))*(cph+S*i1*sph)**(mo)*(-1)**abs(mo)*exp(I*zz)
      else
         psi = Jn(abs(mo))*(cph+S*i1*sph)**(mo)*exp(I*zz)
      end if

      if(TM == 1)then
         gte1 =   Wlsm*psi*m*Qlm(jlm(l,m))
         gtm1 = i1*Wlsm*psi*sth**2*dQlm(jlm(l,m))
      else
         gte1 = i1*Wlsm*psi*sth**2*dQlm(jlm(l,m))
         gtm1 =  -Wlsm*psi*m*Qlm(jlm(l,m))
      end if

      GTE(jlm(l,m)) = gte1
      GTM(jlm(l,m)) = gtm1
   end do
end do

GTE2 = sparse_matmul(rotD, indD, GTE, nm)
GTM2 = sparse_matmul(rotD, indD, GTM, nm)

! Write result to matrices%a etc
matrices%as(1:lmax,i) = GTE2(1:lmax)
matrices%bs(1:lmax,i) = GTM2(1:lmax)

end subroutine bessel_beam_z

!*******************************************************************************
! The associated Legendre functions Qlm and their derivatives
subroutine vswf_qlm(x,lmax,Qlm,dQlm)
real(dp), intent(inout) :: x
real(dp) :: Qlm(0:lmax*(lmax+2)+1), dQlm(0:lmax*(lmax+2)+1)
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
dQlm(jlm(0,0))=0d0
dQlm(jlm(1,1))= -cth*Qlm(jlm(1, 1))*cs2
dQlm(jlm(1,0))=-(cth*Qlm(jlm(1, 0))-sqrt(3d0)*Qlm(jlm(0,0)))*cs2
dQlm(jlm(1,-1))= -cth*Qlm(jlm(1,-1))*cs2


do l = 2,lmax
   ! Qlm positives extremes 
   Qlm(jlm(l, l  ))=-gammaQ(l)*sth*Qlm(jlm(l-1,l-1))
   Qlm(jlm(l, l-1))= deltaQ(l)*cth*Qlm(jlm(l-1,l-1))
   ! dQlm positives extremes 
   Qlm(jlm(l,-l  ))=(-1)**l*Qlm(jlm(l,l  ))
   Qlm(jlm(l,-l+1))=(-1)**(l-1)*Qlm(jlm(l,l-1))
               
   kl=1d0/(sqrt(dble(l)*(l+1)))

   do m = l,0,-1
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
subroutine bess_cyl(nmax, x, Jn, Dn)
integer, intent(in) :: nmax
real(dp), intent(in) :: x
real(dp), intent(inout) :: Jn(0:nmax*(nmax+2)), Dn(0:nmax*(nmax+2))
integer :: l, m, digits, n, j
real(dp) :: kl, cth, sth, cs2, pf, pd, J0, J1, j_n, d_n

digits = 15

call bess_csv(nmax, digits, x, j_n, d_n)
Jn(nmax) = j_n
Dn(nmax) = d_n

! Downward recursion
do n = nmax, 1, -1
   Jn(n-1) = dble(n)/x*Jn(n)+Dn(n)
   Dn(n-1) = dble(n-1)/x*Jn(n-1)-Jn(n)
   if(abs(Jn(n-1))>1d2) then
      do j = nmax, n-1, -1
         Dn(j) = Dn(j)/Jn(n-1)
         Jn(j) = Jn(j)/Jn(n-1)
      end do
   end if
end do

! Normalization factor for the function
pf = 1d0
J0 = bessel_jn(0,x)
J1 = bessel_jn(1,x)
if(abs(J0)>1d-10) then
   pf = bessel_jn(0,x)/Jn(0)
else
   pf = bessel_jn(1,x)/Jn(1)
end if

! Normalization factor for the derivative
pd = 1d0
if(abs(J1)>1d-10)then
   pd = -J1/Dn(0)
else  
   pd = 0.5d0*(Jn(0)-Jn(2))/Dn(0)
end if

! Normalizations
do n = 0, nmax
   Jn(n) = Jn(n)*pf
   Dn(n) = Dn(n)*pd
end do 

end subroutine bess_cyl

!******************************************************************************
! Starting values for downward recurrence, Cylindrical Bessel func
! For calculation of J_n(x), one can calculate by downward recurrence from
! J_{n-dig}(x), where dig is the number of precision in J_n(x).
subroutine bess_csv(nmax, dig, x, Jnn, Dnn)
integer, intent(in) :: nmax, dig
real(dp), intent(in) :: x
real(dp), intent(out) :: Jnn, Dnn
integer :: n, m
real(dp) :: Jn(dig), Dn(dig)

Jn = 0d0
Dn = 0d0

Jn(dig) = 0d0
Dn(dig) = 1d0

do n = dig, 2, -1
   Jn(n-1) = (dble(nmax+n)/x)*Jn(n)+Dn(n)
   Dn(n-1) = (dble(nmax+n-1)/x)*Jn(n-1)-Jn(n)
   if(abs(Jn(n-1))>1d2)then
      do m = dig, n-1, -1
         Dn(m) = Dn(m)/Jn(n-1)
         Jn(m) = Jn(m)/Jn(n-1)
      end do
   end if
end do 

Jnn = Jn(1)
Dnn = Dn(1)

end subroutine bess_csv

!******************************************************************************

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


end module bessel