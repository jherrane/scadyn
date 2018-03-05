module bessel
use T_matrix

implicit none
contains

!******************************************************************************
! The routine replaces the incident beams with Bessel beams.
subroutine bessel_beams(matrices, mesh)
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
   call bessel_z(Nmax, dble(mesh%ki(i)), a_in, b_in)

   matrices%as(1:nm,i) = a_in
   matrices%bs(1:nm,i) = b_in
end do


end subroutine bessel_beams

!******************************************************************************
! The routine computes the SVWF expansion coefficients 
! for a time harmonic x-polarized Bessel beams
! propagating +z-direction  with the wave number k
subroutine bessel_z(Nmax, k, a_nm, b_nm)
integer :: Nmax, ind, n, m, mm
complex(dp), dimension((Nmax+1)**2-1) :: a_nm, b_nm
real(dp) :: Pnm(0:Nmax*(Nmax+2)+1), dPnm(0:Nmax*(Nmax+2)+1)
real(dp) :: k, C, E0, rho, x, y, z, vec(3), r, theta, phi, &
pi_tilde, tau_tilde, cth
complex(dp) :: Un, Ip, Im

E0 = 1d0
ind = 0
x = 1d-7
y = 0d-7
z = 1d-7
vec = cart2sph([x,y,z])
r = vec(1)
theta = vec(2)
phi = vec(3)
rho = k*sqrt(x**2 + y**2)*sin(theta)
cth = cos(theta)

call legendre_Pnm(cth,Nmax,Pnm,dPnm)

do n = 1,Nmax
   Un = 4d0*pi*i1**(n)/(n*(n+1))

   do m = -n, n 
      ind = ind + 1
      mm = abs(m)

      C = sqrt((2d0*n+1.0d0)*factorial(n-mm)/factorial(n+mm)/(4d0*pi))
      pi_tilde = C*m/sin(theta)*Pnm(jnm(n,m))
      tau_tilde = C*dPnm(jnm(n,m))

      Ip = pi*(exp(i1*(m-1)*phi)*bessel_jn(1-m,rho) + &
      exp(i1*(m+1)*phi)*bessel_jn(-1-m,rho)) 
      Im = pi*(exp(i1*(m-1)*phi)*bessel_jn(1-m,rho) - &
      exp(i1*(m+1)*phi)*bessel_jn(-1-m,rho)) 

      a_nm(ind) = E0*Un*exp(i1*k*z*cos(theta))*(tau_tilde*Ip + pi_tilde*Im)
      b_nm(ind) = E0*Un*exp(i1*k*z*cos(theta))*(pi_tilde*Ip+tau_tilde*Im)
   end do 
end do

end subroutine bessel_z

!*******************************************************************************
! The associated Legendre functions Pnm and their derivatives
subroutine legendre_Pnm(cth,Nmax,Pnm,dPnm)
real(dp), intent(in) :: cth
real(dp) :: Pnm(0:Nmax*(Nmax+2)+1), dPnm(0:Nmax*(Nmax+2)+1)
integer :: n, m, Nmax
real(dp) :: kn, sth, cs2, alphaP(1:Nmax, 0:Nmax), betaP(1:Nmax, 0:Nmax),&
gammaP(1:Nmax), deltaP(1:Nmax) 

sth = sqrt(1d0-cth**2)
cs2 = sth**(-2)

do n = 1, Nmax
   do m = n,0,-1
      alphaP(n,m) = sqrt(((2d0*n-1)*(2d0*n+1))/((n-m)*(n+m)))
      betaP(n,m) = sqrt((2d0*n+1)/(2*n-3))*sqrt((dble(n+m-1)*(n-m-1))/((n-m)*(n+m)))
   end do
   gammaP(n) = sqrt((2d0*n+1)/(2d0*n))
   deltaP(n) = sqrt(2d0*n+1)
end do

! Pnm - First 4 terms
Pnm(jnm(0 ,0)) = 1/sqrt(4d0*pi)
Pnm(jnm(1, 1)) = -gammaP(1)*sth*Pnm(jnm(0,0))! P11
Pnm(jnm(1, 0)) = sqrt(3d0)*cth*Pnm(jnm(0,0))  ! P10
Pnm(jnm(1,-1)) = -Pnm(jnm(1,1))               ! P11*(-1)

! dPnm - First 4 terms
dPnm(jnm(0, 0))  = 0d0
dPnm(jnm(1, 1))  = -cth*Pnm(jnm(1, 1))*cs2
dPnm(jnm(1, 0))  = -(cth*Pnm(jnm(1, 0))-sqrt(3d0)*Pnm(jnm(0,0)))*cs2
dPnm(jnm(1, -1)) = -cth*Pnm(jnm(1,-1))*cs2

do n = 2,Nmax
   ! Pnm positives extremes 
   Pnm(jnm(n, n  ))  = -gammaP(n)*sth*Pnm(jnm(n-1,n-1))
   Pnm(jnm(n, n-1))  = deltaP(n)*cth*Pnm(jnm(n-1,n-1))
   ! dPnm positives extremes 
   Pnm(jnm(n, -n  )) = (-1)**n*Pnm(jnm(n,n  ))
   Pnm(jnm(n, -n+1)) = (-1)**(n-1)*Pnm(jnm(n,n-1))
               
   kn=1d0/(sqrt(dble(n)*(n+1)))

   do m = n,0,-1
      if(m<n-1) then
         Pnm(jnm(n, m)) = alphaP(n,m)*cth*Pnm(jnm(n-1,m)) &
            -betaP(n,m)*Pnm(jnm(n-2,m)) 
         Pnm(jnm(n,-m)) = (-1)**m*Pnm(jnm(n, m))
      end if

      dPnm(jnm(n, m)) = -n*cth*cs2*Pnm(jnm(n, m))+ &
      sqrt((2d0*n+1.)/(2*n-1.))*sqrt(dble(n+m)*(n-m))*cs2*Pnm(jnm(n-1, m))

      dPnm(jnm(n,-m)) = -n*cth*cs2*Pnm(jnm(n,-m))+ &
      sqrt((2d0*n+1.)/(2*n-1.))*sqrt(dble(n+m)*(n-m))*cs2*Pnm(jnm(n-1,-m))
   end do 
end do 

end subroutine legendre_Pnm

!*******************************************************************************

function jnm(n,m) result(j)
integer :: n, m, j
if(abs(m)>l) then
    j = 0
else
    j = n*(n+1)+m
endif
end function jnm

end module bessel
