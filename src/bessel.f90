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
complex(dp) :: GTE(0:lmx*(lmx+2))! The final expansion, first half (a)
complex(dp) :: GTM(0:lmx*(lmx+2)) ! The final expansion, second half (b)
real(dp) :: rho, gr, cth, sth, cph, sph, zz, lambda, a, b, c, x, y, J0, J1, pf, pd
real(dp) :: Jn(0:lmx*(lmx+2)), Dn(0:lmx*(lmx+2))
real(dp) :: Qlm(0:lmx*(lmx+2)+1), dQlm(0:lmx*(lmx+2)+1)
real(dp) :: Pmn(0:lmx, 0:lmx), DPmn(0:lmx, 0:lmx)
integer, dimension(:), allocatable :: ll, mm
integer :: ls, le, mo, l, m, iii, ind
complex(dp) :: Wlsm, psi, i1, gte1, gtm1, gtem(2)

i1 = dcmplx(0d0,1d0)
MMM = 3
NNN = 5
S = -1
TM = 0
x = 1d-7
y = 1d-7

lmax = lmx*(lmx+2)

allocate(ll(lmax), mm(lmax))

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

gr = gamma*rho

cth = kz/k
sth = gamma/k
cph = x/rho
sph = y/rho

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

ind = 0
do l = 1,lmx 
   do m = -l,l
      ind = ind + 1
      mo = MMM-S*m
      
      if(mo<0) then
         psi = bessel_jn(mo,gr)*(cph+S*i1*sph)**(mo)*(-1)**abs(mo)*exp(I*zz)
      else
         psi = bessel_jn(mo,gr)*(cph+S*i1*sph)**(mo)*exp(I*zz)
      end if

      if(l>0)then
         Wlsm = 4d0*pi*i1**(l-S*m)/sqrt(dble(l*(l+1)))
      else
         Wlsm = 0d0
      end if

      if(TM == 1)then
         gte1 =   Wlsm*psi*m*Qlm(jlm(l,m))
         gtm1 = i1*Wlsm*psi*sth**2*dQlm(jlm(l,m))
      else
         gte1 = i1*Wlsm*psi*sth**2*dQlm(jlm(l,m))
         gtm1 =  -Wlsm*psi*m*Qlm(jlm(l,m))
      end if

      gtem = [gte1, gtm1]

      ! Fixing signs and stuff
      if(mod(ind,2)/=0)then
         gtem = gtem/(-i1)**(ind)
      else
         gtem = gtem/(-i1)**(2*ind-1)
      end if

      GTE(jlm(l,m)) = gtem(1)
      GTM(jlm(l,m)) = gtem(2)
   end do
end do
! Write result to matrices%a etc
matrices%as(1:lmax,i) = GTE(1:lmax)
matrices%bs(1:lmax,i) = GTM(1:lmax)

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

!******************************************************************************

subroutine bbz_taylor(matrices, mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: i
integer :: Nmax, nm
complex(dp), dimension(:), allocatable :: a_in, b_in

do i = 1,matrices%bars
   if(allocated(a_in)) deallocate(a_in,b_in)
   Nmax = matrices%Nmaxs(i)
   nm = (Nmax+1)**2-1
   allocate(a_in(nm))
   allocate(b_in(nm))
   call bess(Nmax, dble(mesh%ki(i)), a_in, b_in)

   matrices%as(1:nm,i) = a_in
   matrices%bs(1:nm,i) = b_in
end do


end subroutine bbz_taylor

!******************************************************************************
! The routine computes the SVWF expansion coefficients 
! for a time harmonic x-polarized Bessel beams
! propagating +z-direction  with the wave number k
subroutine bess(Nmax, k, a_nm, b_nm)
complex(dp), dimension((Nmax+1)**2-1) :: a_nm, b_nm
integer :: Nmax
real(dp) :: k

integer :: ind, n, m, mm
real(dp) :: C, E0, rho, x, y, z, vec(3), r, theta, phi, &
pi_tilde, tau_tilde, cth
complex(dp) :: Un, Ip, Im
real(dp) :: Pnm(0:Nmax*(Nmax+2)+1), dPnm(0:Nmax*(Nmax+2)+1)

E0 = 1d0
ind = 0
x = 1d-7
y = 0d0
z = 1d-7
vec = cart2sph([x,y,z])
r = vec(1)
theta = vec(2)
phi = vec(3)
rho = k*sqrt(x**2 + y**2)*sin(theta)
cth = cos(theta)

call vswf_qlm(cth,Nmax,Pnm,dPnm)

do n = 1,Nmax
   Un = 4d0*pi*i1**(n)/(n*(n+1))

   do m = -n, n 
      ind = ind + 1
      mm = abs(m)

      C = sqrt((2d0*n+1.0d0)*factorial(n-mm)/factorial(n+mm)/(4d0*pi))
      pi_tilde = C*m/sin(theta)*Pnm(jlm(n,m))
      tau_tilde = C*dPnm(jlm(n,m))

      Ip = pi*(exp(i1*(m-1)*phi)*bessel_jn(1-m,rho) + &
      exp(i1*(m+1)*phi)*bessel_jn(-1-m,rho)) 
      Im = pi*(exp(i1*(m-1)*phi)*bessel_jn(1-m,rho) - &
      exp(i1*(m+1)*phi)*bessel_jn(-1-m,rho)) 

      a_nm(ind) = E0*Un*exp(i1*k*z*cos(theta))*(tau_tilde*Ip + pi_tilde*Im)
      b_nm(ind) = E0*Un*exp(i1*k*z*cos(theta))*(pi_tilde*Ip+tau_tilde*Im)

   end do 
end do

end subroutine bess

end module bessel