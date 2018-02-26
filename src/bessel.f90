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
complex(dp), dimension(:), allocatable :: GTE  ! The final expansion, first half (a)
complex(dp), dimension(:), allocatable :: GTM  ! The final expansion, second half (b)
real(dp) :: rho, gr, cth, sth, cph, sph, zz, lambda, a, b, c, x, y, J0, J1, pf, pd
real(dp), dimension(:), allocatable :: Qlm, dQlm, Jn, dn, Yn, DYn
integer, dimension(:), allocatable :: ll, mm
integer :: ls, le, mo, l, m
complex(dp) :: Wlsm, psi, i1

i1 = dcmplx(0d0,1d0)
MMM = 3
NNN = 5
S = -1
TM = 1
x = 1d-7
y = 1d-7

lmax = lmx*(lmx+2)+1

allocate(Qlm(lmax), dQlm(lmax))
allocate(ll(lmax), mm(lmax))
allocate(GTE(lmax),GTM(lmax))

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

ls = abs(MMM-lmax-1)
le = abs(MMM+lmax+1)

if(ls>le) then
    le = ls
    ls = 0
endif
allocate(Jn(le), Dn(le), Yn(le), DYn(le))
Jn = 0d0
Dn = 0d0

call lpmn(lmx, MMM, lmx, cth, Qlm, dQlm )
call jyna(lmax, gr, lmax, Jn, Dn, Yn, DYn )

! Normalizations
pf = 1d0
J0 = bessel_jn(0,x)
J1 = bessel_jn(1,x)
if(abs(J0)>1d-10)then
   pf = J0/Jn(1)
else
   pf = J1/Jn(2)
end if

pd = 1d0
if(abs(J1)>1d-10)then
   pd = -J1/Dn(1)
else
   pd = 0.5*(Jn(1)-Jn(3))/Dn(2)
end if

do l = 1, lmax
   Jn(l) = Jn(l)*pf
   Dn(l) = Dn(l)*pd
end do

do l = 1,lmx
    do m = -l,l
        mo = MMM-S*m
        Wlsm = 4d0*pi*i1**(l-S*m)/sqrt(dble(l*(l+1)))
        if(mo<0) then
            psi = Jn(abs(mo)+1)*(cph+S*i1*sph)**(mo)*1*(-1)**abs(mo) !*exp(I*zz=0)
        else
            psi = Jn(abs(mo)+1)*(cph+S*i1*sph)**(mo)*1 !*exp(I*zz=0)
        end if
        if(TM == 1)then
            GTE(jlm(l,m)) =   Wlsm*psi*m*Qlm(jlm(l,m))
            GTM(jlm(l,m)) = i1*Wlsm*psi*sth**2*dQlm(jlm(l,m))
        else
            GTE(jlm(l,m)) = i1*Wlsm*psi*sth**2*dQlm(jlm(l,m))
            GTM(jlm(l,m)) =  -Wlsm*psi*m*Qlm(jlm(l,m))
        end if
    end do
end do

! Write result to matrices%a etc
matrices%as(1:lmax-1,i) = GTE(1:lmax-1)
matrices%bs(1:lmax-1,i) = GTM(1:lmax-1)

end subroutine bessel_beam_z

!*******************************************************************************
! The associated Legendre functions Qlm and their derivatives
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
Qlm(jlm(1,1))=1/sqrt(4d0*pi)
Qlm(jlm(2, 2))=-gammaQ(1)*sth*Qlm(jlm(1,1))! Q11
Qlm(jlm(2, 1))=sqrt(3d0)*cth*Qlm(jlm(1,1))  ! Q10
Qlm(jlm(2,-2))=-Qlm(jlm(1,1))               ! Q11*(-1)

! dQlm - First 4 terms
dQlm(jlm(0,0))=0d0
dQlm(jlm(1,1))= -cth*Qlm(jlm(2, 2))*cs2
dQlm(jlm(1,0))=-(cth*Qlm(jlm(2, 1))-sqrt(3d0)*Qlm(jlm(1,1)))*cs2
dQlm(jlm(1,-1))= -cth*Qlm(jlm(2,-2))*cs2

! lo - First 4 terms
lo(jlm(0,0)) = 1
lo(jlm(1,-1)) = 2
lo(jlm(1, 0)) = 2
lo(jlm(1, 1)) = 2

! mo - First 4 terms
mo(jlm(0, 0)) =  1
mo(jlm(1,-1)) = -2
mo(jlm(1, 0)) =  1
mo(jlm(1, 1)) =  2

do l = 2,lmax
   ! Qlm positives extremes 
   Qlm(jlm(l, l  ))=-gammaQ(l)*sth*Qlm(jlm(l-1,l-1))
   Qlm(jlm(l, l-1))= deltaQ(l)*cth*Qlm(jlm(l-1,l-1))
   ! dQlm positives extremes 
   Qlm(jlm(l,-l  ))=(-1)**l*Qlm(jlm(l,l  ))
   Qlm(jlm(l,-l+1))=(-1)**(l-1)*Qlm(jlm(l,l-1))
               
   kl=1d0/(sqrt(dble(l)*(l+1)))

   do m = l,0,-1
      ! print*, 'jlm(',l,', ',m,') = ', jlm(l, m)
      ! print*, 'jlm(',l,', -',m,') = ', jlm(l, -m)
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
! print*, 'lo = ', lo

end subroutine vswf_qlm

! !*******************************************************************************
! ! Array of cylindrical Bessel functions
! subroutine bess_cyl(nmax, x, Jn, Dn)
! real(dp), intent(inout) :: x
! real(dp), dimension(:), allocatable, intent(inout) :: Jn, Dn 
! integer, intent(inout) :: nmax
! integer :: l, m, digits, n, j
! real(dp) :: kl, cth, sth, cs2, pf, pd, J0, J1, j_n, d_n

! digits = 15

! call bess_csv(nmax, digits, x, j_n, d_n)
! Jn(nmax) = j_n
! Dn(nmax) = d_n

! ! Downward recursion
! do n = nmax, 2, -1
!    Jn(n-1) = dble(n)/x*Jn(n)+Dn(n)
!    Dn(n-1) = dble(n-1)/x*Jn(n-1)-Jn(n)
!    if(abs(Jn(n-1))>1d2) then
!       do j = nmax, n-1, -1
!          Dn(j) = Dn(j)/Jn(n-1)
!          Jn(j) = Jn(j)/Jn(n-1)
!       end do
!    end if
! end do

! ! Normalization factor for the function
! pf = 1d0
! J0 = bessel_jn(0,x)
! J1 = bessel_jn(1,x)

! if(abs(J0)>1d-10) then
!    pf = bessel_jn(0,x)/Jn(1)
! else
!    pf = bessel_jn(1,x)/Jn(2)
! end if

! ! Normalization factor for the derivative
! pd = 1d0
! if(abs(J1)>1d-10)then
!    pd = -J1/Dn(1)
! else  
!    pd = 0.5d0*(Jn(1)-Jn(3))/Dn(2)
! end if

! ! Normalizations
! do n = 1, nmax
!    Jn(n) = Jn(n)*pf
!    Dn(n) = Dn(n)*pd
! end do 

! end subroutine bess_cyl

! !******************************************************************************
! ! Starting values for downward recurrence, Cylindrical Bessel func
! ! For calculation of J_n(x), one can calculate by downward recurrence from
! ! J_{n-dig}(x), where dig is the number of precision in J_n(x).
! subroutine bess_csv(nmax, dig, x, Jnn, Dnn)
! integer, intent(in) :: nmax, dig
! real(dp), intent(in) :: x
! real(dp), intent(out) :: Jnn, Dnn
! integer :: n, m
! real(dp) :: Jn(dig), Dn(dig)

! Jn = 0d0
! Dn = 0d0

! Jn(dig) = 0d0
! Dn(dig) = 1d0

! do n = dig, 2, -1
!    Jn(n-1) = (dble(nmax+n)/x)*Jn(n)+Dn(n)
!    Dn(n-1) = (dble(nmax+n-1)/x)*Jn(n-1)-Jn(n)
!    if(abs(Jn(n-1))>1d2)then
!       do m = dig, n-1, -1
!          Dn(m) = Dn(m)/Jn(n-1)
!          Jn(m) = Jn(m)/Jn(n-1)
!       end do
!    end if
! end do 

! Jnn = Jn(1)
! Dnn = Dn(1)

! end subroutine bess_csv

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
    j = 1
else
    j = l*(l+1)+m+1
endif
end function jlm


end module bessel