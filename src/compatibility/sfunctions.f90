module sfunctions
   use common
   implicit none
contains
!******************************************************************************
!## SFUNCTIONS ##
!     ==========================================================
!     Purpose: Compute spherical Bessel functions jn(z) & yn(z)
!              for a complex argument
!     Input :  z --- Complex argument
!              n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
!     Output:  CSJ(n) --- jn(z)
!              CSY(n) --- yn(z)
!              NM --- Highest order computed
!     Routines called:
!              MSTA1 and MSTA2 for computing the starting
!              point for backward recurrence
!     ==========================================================
   subroutine cspherebessel(n, z, csj, csy)
      implicit none
      integer :: n, nm, k, m
      real(8) :: a0
      complex(8) :: z, csj(0:n), csy(0:n), csa, csb, cs, cf0, cf1, cf
      a0 = cdabs(z)
      nm = n
      if (a0 .lt. 1.0d-60) then
         csj = (0.d0, 0.d0)
         csy = (-1.d300, 0.d0)
         csy(0) = (1.d0, 0.d0)
         return
      endif
      csj = (0.d0, 0.d0)
      csj(0) = cdsin(z)/z
      csj(1) = (csj(0) - cdcos(z))/z
      if (n .ge. 2) then
         csa = csj(0)
         csb = csj(1)
         m = msta1(a0, 200)
         if (m .lt. n) then
            nm = m
         else
            m = msta2(a0, n, 15)
         endif
         cf0 = dcmplx(0d0, 0d0)
         cf1 = dcmplx(1.0d0 - 100)
         do k = m, 0, -1
            cf = (2.0d0*k + 3.0d0)*cf1/z - cf0
            if (k .le. nm) csj(k) = cf
            cf0 = cf1
            cf1 = cf
         enddo
         if (cdabs(csa) .gt. cdabs(csb)) cs = csa/cf
         if (cdabs(csa) .le. cdabs(csb)) cs = csb/cf0
         do k = 0, min(nm, n)
            csj(k) = cs*csj(k)
         enddo
      endif
      csy = (1.d200, 0.d0)
      csy(0) = -cdcos(z)/z
      csy(1) = (csy(0) - cdsin(z))/z
      do k = 2, min(nm, n)
         if (cdabs(csj(k - 1)) .gt. cdabs(csj(k - 2))) then
            csy(k) = (csj(k)*csy(k - 1) - 1.0d0/(z*z))/csj(k - 1)
         else
            csy(k) = (csj(k)*csy(k - 2) - (2.0d0*k - 1.0d0)/z**3)/csj(k - 2)
         endif
      enddo
   end subroutine cspherebessel
!******************************************************************************
!     ===================================================
!     Purpose: Determine the starting point for backward
!              recurrence such that the magnitude of
!              Jn(x) at that point is about 10^(-MP)
!     Input :  x     --- Argument of Jn(x)
!              MP    --- Value of magnitude
!     Output:  MSTA1 --- Starting point
!     ===================================================
   integer function msta1(x, mp)
      implicit none
      integer :: mp, n0, n1, it, nn
      real(8) :: x, a0, f1, f, f0
      a0 = dabs(x)
      n0 = int(1.1*a0) + 1
      f0 = envj(n0, a0) - dble(mp)
      n1 = n0 + 5
      f1 = envj(n1, a0) - dble(mp)
      do it = 1, 20
         nn = n1 - int((n1 - n0)/(1.0d0 - f0/f1))
         f = envj(nn, a0) - dble(mp)
         if (abs(nn - n1) .lt. 1) exit
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      enddo
      msta1 = nn
   end function msta1

!******************************************************************************
!     ===================================================
!     Purpose: Determine the starting point for backward
!              recurrence such that all Jn(x) has MP
!              significant digits
!     Input :  x  --- Argument of Jn(x)
!              n  --- Order of Jn(x)
!              MP --- Significant digit
!     Output:  MSTA2 --- Starting point
!     ===================================================
   integer function msta2(x, n, mp)
      implicit none
      integer :: n, mp, n0, n1, it, nn
      real(8) :: x, a0, hmp, ejn, obj, f0, f1, f
      a0 = dabs(x)
      hmp = 0.5d0*dble(mp)
      ejn = envj(n, a0)
      if (ejn .le. hmp) then
         obj = dble(mp)
         n0 = int(1.1*a0)
      else
         obj = hmp + ejn
         n0 = n
      endif
      f0 = envj(n0, a0) - obj
      n1 = n0 + 5
      f1 = envj(n1, a0) - obj
      do it = 1, 20
         nn = n1 - int((n1 - n0)/(1.0d0 - f0/f1))
         f = envj(nn, a0) - obj
         if (abs(nn - n1) .lt. 1) exit
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      enddo
      msta2 = nn + 10
   end function msta2

!******************************************************************************

   real(8) function envj(n, x)
      implicit none
      integer :: n
      real(8) :: x
      n = max(1, abs(n))
      envj = 0.5d0*dlog10(6.28d0*n) - dble(n)*dlog10(1.36d0*x/n)
   end function envj

!******************************************************************************

   subroutine sphj2(n, x, nm, sj, dj)
!*****************************************************************************80
!
!  SPHJ computes spherical Bessel functions jn(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) SJ(0:N), the values of jn(x).
!
!    Output, real ( kind = 8 ) DJ(0:N), the values of jn'(x).
!
      implicit none

      integer :: n

      real(dp) :: cs
      real(dp) ::  dj(0:n)
      real(dp) :: f
      real(dp) :: f0
      real(dp) :: f1
      integer :: k
      integer  :: m
      integer  :: nm
      real(dp) ::  sa
      real(dp) :: sb
      real(dp) :: sj(0:n)
      real(dp) :: x

      nm = n

      if (abs(x) <= 1.0D-100) then
         do k = 0, n
            sj(k) = 0.0D+00
            dj(k) = 0.0D+00
         end do
         sj(0) = 1.0D+00
         dj(1) = 0.3333333333333333D+00
         return
      end if

      sj(0) = sin(x)/x
      sj(1) = (sj(0) - cos(x))/x

      if (2 <= n) then

         sa = sj(0)
         sb = sj(1)
         m = msta1(x, 200)
         if (m < n) then
            nm = m
         else
            m = msta2(x, n, 15)
         end if

         f0 = 0.0D+00
         f1 = 1.0D+00-100
         do k = m, 0, -1
            f = (2.0D+00*k + 3.0D+00)*f1/x - f0
            if (k <= nm) then
               sj(k) = f
            end if
            f0 = f1
            f1 = f
         end do

         if (abs(sa) <= abs(sb)) then
            cs = sb/f0
         else
            cs = sa/f
         end if

         do k = 0, nm
            sj(k) = cs*sj(k)
         end do

      end if

      dj(0) = (cos(x) - sin(x)/x)/x
      do k = 1, nm
         dj(k) = sj(k - 1) - (k + 1.0D+00)*sj(k)/x
      end do

      return
   end subroutine sphj2

   subroutine sphy(n, x, nm, sy, dy)

!*****************************************************************************80
!
!! SPHY computes spherical Bessel functions yn(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) SY(0:N), DY(0:N), the values of yn(x) and yn'(x).
!
      implicit none

      integer(kind=4) n

      real(kind=8) dy(0:n)
      real(kind=8) f
      real(kind=8) f0
      real(kind=8) f1
      integer(kind=4) k
      integer(kind=4) nm
      real(kind=8) sy(0:n)
      real(kind=8) x

      nm = n

      if (x < 1.0D-60) then
         do k = 0, n
            sy(k) = -1.0D+300
            dy(k) = 1.0D+300
         end do
         return
      end if

      sy(0) = -cos(x)/x
      sy(1) = (sy(0) - sin(x))/x
      f0 = sy(0)
      f1 = sy(1)
      do k = 2, n
         f = (2.0D+00*k - 1.0D+00)*f1/x - f0
         sy(k) = f
         if (1.0D+300 <= abs(f)) then
            exit
         end if
         f0 = f1
         f1 = f
      end do

      nm = k - 1
      dy(0) = (sin(x) + cos(x)/x)/x
      do k = 1, nm
         dy(k) = sy(k - 1) - (k + 1.0D+00)*sy(k)/x
      end do

      return
   end subroutine sphy

   subroutine csphjy(n, z, nm, csj, cdj, csy, cdy)
!*****************************************************************************80
!
!! CSPHJY: spherical Bessel functions jn(z) and yn(z) for complex argument.
!
!  Discussion:
!
!    This procedure computes spherical Bessel functions jn(z) and yn(z)
!    and their derivatives for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of jn(z) and yn(z).
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, complex ( kind = 8 ) CSJ(0:N0, CDJ(0:N), CSY(0:N), CDY(0:N),
!    the values of jn(z), jn'(z), yn(z), yn'(z).
!
      implicit none

      integer(kind=4) n

      real(kind=8) a0
      complex(kind=8) csj(0:n)
      complex(kind=8) cdj(0:n)
      complex(kind=8) csy(0:n)
      complex(kind=8) cdy(0:n)
      complex(kind=8) cf
      complex(kind=8) cf0
      complex(kind=8) cf1
      complex(kind=8) cs
      complex(kind=8) csa
      complex(kind=8) csb
      integer(kind=4) k
      integer(kind=4) m
      integer(kind=4) nm
      complex(kind=8) z

      a0 = abs(z)
      nm = n

      if (a0 < 1.0D-60) then
         do k = 0, n
            csj(k) = dcmplx(0.0D+00)
            cdj(k) = dcmplx(0.0D+00)
            csy(k) = dcmplx(-1.0D+300)
            cdy(k) = dcmplx(1.0D+300)
         end do
         csj(0) = cmplx(1.0D+00, 0.0D+00, kind=8)
         cdj(1) = cmplx(0.333333333333333D+00, 0.0D+00, kind=8)
         return
      end if

      csj(0) = sin(z)/z
      csj(1) = (csj(0) - cos(z))/z

      if (2 <= n) then
         csa = csj(0)
         csb = csj(1)
         m = msta1(a0, 200)
         if (m < n) then
            nm = m
         else
            m = msta2(a0, n, 15)
         end if
         cf0 = dcmplx(0.0D+00)
         cf1 = dcmplx(1.0D+00-100)
         do k = m, 0, -1
            cf = (2.0D+00*k + 3.0D+00)*cf1/z - cf0
            if (k <= nm) then
               csj(k) = cf
            end if
            cf0 = cf1
            cf1 = cf
         end do

         if (abs(csa) <= abs(csb)) then
            cs = csb/cf0
         else
            cs = csa/cf
         end if

         do k = 0, nm
            csj(k) = cs*csj(k)
         end do

      end if

      cdj(0) = (cos(z) - sin(z)/z)/z
      do k = 1, nm
         cdj(k) = csj(k - 1) - (k + 1.0D+00)*csj(k)/z
      end do
      csy(0) = -cos(z)/z
      csy(1) = (csy(0) - sin(z))/z
      cdy(0) = (sin(z) + cos(z)/z)/z
      cdy(1) = (2.0D+00*cdy(0) - cos(z))/z

      do k = 2, nm
         if (abs(csj(k - 2)) < abs(csj(k - 1))) then
            csy(k) = (csj(k)*csy(k - 1) - 1.0D+00/(z*z))/csj(k - 1)
         else
            csy(k) = (csj(k)*csy(k - 2) &
                      - (2.0D+00*k - 1.0D+00)/z**3)/csj(k - 2)
         end if
      end do

      do k = 2, nm
         cdy(k) = csy(k - 1) - (k + 1.0D+00)*csy(k)/z
      end do

      return
   end subroutine csphjy

!******************************************************************************

   subroutine cspherebessel2(n, z, csj, csy)
      implicit none
      integer :: n, nm
      complex(8) :: z, csj(0:n), csy(0:n), cdj(0:n), cdy(0:n)

      call csphjy(n, z, nm, csj, cdj, csy, cdy)

   end subroutine cspherebessel2

!******************************************************************************
! Returns Spherical Hankel functions of the 1st kind
! of orders 0:N at x
! H = sqrt(pi/2/x)* besselh(n+1/2,x)
   function sphankel(N, x) result(H)
      integer :: N, nn
      complex(dp) :: x
      complex(dp) :: H(N + 1)

      H(:) = dcmplx(0.0d0, 0.0d0)

      H(1) = -dcmplx(0.0d0, 1.0d0)*exp(dcmplx(0.0d0, 1.0d0)*x)/x

      if (N > 0) then
         H(2) = (1.0d0/x - dcmplx(0.0d0, 1.0d0))*H(1)
      end if

      do nn = 2, N
         H(nn + 1) = (2.0d0*nn - 1.0d0)/x*H(nn) - H(nn - 1)
      end do

   end function sphankel

!******************************************************************************

   function spbessel(N, x) result(H)
      integer :: N
      complex(dp) :: x
      complex(dp) :: H(N + 1)

      complex(dp) :: j0, j1
      integer :: i1

      j0 = sin(x)/x
      j1 = sin(x)/x**2 - cos(x)/x

      H(1) = j0
      H(2) = j1

      do i1 = 2, N
         H(i1 + 1) = (2.0d0*i1 + 1.0d0)/x*H(i1) - H(i1 - 1)
      end do

   end function spbessel

!******************************************************************************

   subroutine sbesseljd(n, z, j, dj)
      implicit none
      integer :: n, nm, i
      complex(8) :: z, csj(0:n), csy(0:n), cdj(0:n), cdy(0:n), j(n), dj(n)

      call csphjy(n, z, nm, csj, cdj, csy, cdy)
      j = csj(1:n)
      do i = 1, n
         dj(i) = csj(i - 1) - i*csj(i)/z
      end do

   end subroutine sbesseljd

!******************************************************************************

   subroutine shankeljd(n, z, h, dh)
      implicit none
      integer :: n, i
      complex(8) :: z, hnkl(n + 2), h(n), dh(n)

      hnkl = sphankel(n + 1, z)
      h = hnkl(2:n + 1)

      do i = 1, n
         dh(i) = hnkl(i) - i*hnkl(i + 1)/z
      end do

   end subroutine shankeljd

!******************************************************************************
! Riccati-Bessselh up to order N
   function riccati_besselh(N, x) result(zeta)
      integer :: N
      complex(dp) :: x
      complex(dp) :: zeta(N + 1, 2)

      complex(dp) :: sph(N + 2)
      integer :: nn

      sph = sphankel(N + 1, x)

      do nn = 0, N
         zeta(nn + 1, 1) = -x*sph(nn + 1)
         zeta(nn + 1, 2) = -(nn + 1)*sph(nn + 1) + x*sph(nn + 2)
      end do

   end function riccati_besselh

!******************************************************************************

   subroutine riccati_bessel(N, x, zeta, d_zeta, psi, d_psi)
      integer :: N, i1
      complex(dp) :: x
      complex(dp), dimension(N + 1) :: zeta, d_zeta, psi, d_psi
      complex(dp), dimension(N + 2) :: sphj, sphy, sphh

      call cspherebessel(N + 1, x, sphj, sphy)

! spherical hankel of the 2nd kind
!sphh = sphj + dcmplx(0.0, 1.0) * sphy
      sphh = sphankel(N + 1, x)

      do i1 = 1, N + 1

         ! riccati-bessel and derivative
         psi(i1) = x*sphj(i1)
         d_psi(i1) = (i1)*sphj(i1) - x*sphj(i1 + 1)

         ! riccati-hankel and derivative
         zeta(i1) = -x*sphh(i1)
         d_zeta(i1) = -(i1)*sphh(i1) + x*sphh(i1 + 1)

      end do

   end subroutine riccati_bessel

!******************************************************************************

   subroutine legendre2(N, x, P_nm)
      integer :: N
      real(dp) :: x
      real(dp), dimension(N + 1) :: P_nm
      real(dp) :: PM(N + 2, N + 2), PD(N + 2, N + 2)
      integer :: i1, i2, LS
      real(dp) :: xq, xs

      PM(:, :) = 0.0d0
      PD(:, :) = 0.0d0

      PM(1, 1) = 1.0d0

      if (abs(x) == 1.0d0) then

         do i1 = 1, N
            PM(1, i1 + 1) = x**i1
            PD(1, i1 + 1) = 0.5d0*i1*(i1 + 1.0d0)*x**(i1 + 1.0d0)
         end do

         do i2 = 1, N
            do i1 = 1, N
               if (i1 == 1) then
                  PD(i1 + 1, i2 + 1) = 1.0d+300
               else if (i1 == 2) then
                  PD(i1 + 1, i2 + 1) = -0.25d0*(i2 + 2.0d0)*(i2 + 1.0d0)*i2*(i2 - 1.0d0)*x**(i2 + 1d0)
               end if
            end do
         end do

      end if

      LS = 1
      if (abs(x) > 1.0d0) then
         LS = -1
      end if

      xq = sqrt(LS*(1.0d0 - x*x))
      xs = LS*(1.0d0 - x*x)

      do i1 = 1, N
         PM(i1 + 1, i1 + 1) = -LS*(2.0d0*i1 - 1.0d0)*xq*PM(i1, i1)
      end do

      do i1 = 0, N
         PM(i1 + 1, i1 + 2) = (2.0d0*i1 + 1.0d0)*x*PM(i1 + 1, i1 + 1)
      end do

      do i1 = 0, N
         do i2 = i1 + 2, N
            PM(i1 + 1, i2 + 1) = ((2.0d0*i2 - 1.0d0)*x*PM(i1 + 1, i2) - &
                                  (i1 + i2 - 1.0d0)*PM(i1 + 1, i2 - 1))/dble(i2 - i1)

         end do
      end do

      PD(1, 1) = 0d0

      do i2 = 1, N
         PD(1, i2 + 1) = LS*i2*(PM(1, i2) - x*PM(1, i2 + 1))/xs
      end do

      do i1 = 1, N
         do i2 = 1, N
            PD(i1 + 1, i2 + 1) = LS*i1*x*PM(i1 + 1, i2 + 1)/xs + (i2 + i1)* &
                                 (i2 - i1 + 1.0d0)/xq*PM(i1, i2 + 1)
         end do
      end do

!print*, size(PM,1), size(PM,2)!, size(P_nm)
      P_nm = PM(1:N + 1, N + 1)

!do i1 = 1, N+1
!   print*, PM(i1,N+1)
!end do
   end subroutine legendre2

   function sphankel_a(N, x) result(H)
      integer :: N, nn
      complex(dp) :: x
      complex(dp) :: H(N + 1)

      H(:) = dcmplx(0.0, 0.0)

      do nn = 1, N
         H(nn) = dcmplx(0.0, -1.0)**(nn)*exp(dcmplx(0.0, 1.0)*x)/x
      end do

   end function sphankel_a

   function laguerre(p, l, x) result(Lpl)
      real(dp), dimension(:), allocatable :: x, Lpl
      integer :: n, m, p, l, j
      n = size(x, 1)
      allocate (Lpl(n))
      Lpl = 1d0
      Lpl = Lpl*dble(choose(p + l, p))

      do m = 1, p
         do j = 1, size(x, 1)
            Lpl(j) = Lpl(j) + (-1d0)**dble(m)/factorial(m)*dble(choose(p + l, p - m))*x(j)**dble(m)
         enddo
      enddo

   end function laguerre

   subroutine jyna(n, x, nm, bj, dj, by, dy)

!*****************************************************************************80
!
!! JYNA computes Bessel functions Jn(x) and Yn(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 April 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) BJ(0:N), DJ(0:N), BY(0:N), DY(0:N), the values
!    of Jn(x), Jn'(x), Yn(x), Yn'(x).
!
      implicit none

      integer(kind=4) n

      real(kind=8) bj(0:n)
      real(kind=8) bj0
      real(kind=8) bj1
      real(kind=8) bjk
      real(kind=8) by(0:n)
      real(kind=8) by0
      real(kind=8) by1
      real(kind=8) cs
      real(kind=8) dj(0:n)
      real(kind=8) dj0
      real(kind=8) dj1
      real(kind=8) dy(0:n)
      real(kind=8) dy0
      real(kind=8) dy1
      real(kind=8) f
      real(kind=8) f0
      real(kind=8) f1
      real(kind=8) f2
      integer(kind=4) k
      integer(kind=4) m
      integer(kind=4) nm
      real(kind=8) x

      nm = n

      if (x < 1.0D-100) then

         do k = 0, n
            bj(k) = 0.0D+00
            dj(k) = 0.0D+00
            by(k) = -1.0D+300
            dy(k) = 1.0D+300
         end do
         bj(0) = 1.0D+00
         dj(1) = 0.5D+00
         return

      end if

      call jy01b(x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1)
      bj(0) = bj0
      bj(1) = bj1
      by(0) = by0
      by(1) = by1
      dj(0) = dj0
      dj(1) = dj1
      dy(0) = dy0
      dy(1) = dy1

      if (n <= 1) then
         return
      end if

      if (n < int(0.9D+00*x)) then

         do k = 2, n
            bjk = 2.0D+00*(k - 1.0D+00)/x*bj1 - bj0
            bj(k) = bjk
            bj0 = bj1
            bj1 = bjk
         end do

      else

         m = msta1(x, 200)

         if (m < n) then
            nm = m
         else
            m = msta2(x, n, 15)
         end if

         f2 = 0.0D+00
         f1 = 1.0D-100
         do k = m, 0, -1
            f = 2.0D+00*(k + 1.0D+00)/x*f1 - f2
            if (k <= nm) then
               bj(k) = f
            end if
            f2 = f1
            f1 = f
         end do

         if (abs(bj1) < abs(bj0)) then
            cs = bj0/f
         else
            cs = bj1/f2
         end if

         do k = 0, nm
            bj(k) = cs*bj(k)
         end do

      end if

      do k = 2, nm
         dj(k) = bj(k - 1) - k/x*bj(k)
      end do

      f0 = by(0)
      f1 = by(1)
      do k = 2, nm
         f = 2.0D+00*(k - 1.0D+00)/x*f1 - f0
         by(k) = f
         f0 = f1
         f1 = f
      end do

      do k = 2, nm
         dy(k) = by(k - 1) - k*by(k)/x
      end do

      return
   end

   subroutine jy01b(x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1)

!*****************************************************************************80
!
!! JY01B computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
!    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
!
      implicit none

      real(kind=8) a0
      real(kind=8) bj0
      real(kind=8) bj1
      real(kind=8) by0
      real(kind=8) by1
      real(kind=8) dj0
      real(kind=8) dj1
      real(kind=8) dy0
      real(kind=8) dy1
      real(kind=8) p0
      real(kind=8) p1
      real(kind=8) pi
      real(kind=8) q0
      real(kind=8) q1
      real(kind=8) t
      real(kind=8) t2
      real(kind=8) ta0
      real(kind=8) ta1
      real(kind=8) x

      pi = 3.141592653589793D+00

      if (x == 0.0D+00) then

         bj0 = 1.0D+00
         bj1 = 0.0D+00
         dj0 = 0.0D+00
         dj1 = 0.5D+00
         by0 = -1.0D+300
         by1 = -1.0D+300
         dy0 = 1.0D+300
         dy1 = 1.0D+300
         return

      else if (x <= 4.0D+00) then

         t = x/4.0D+00
         t2 = t*t

         bj0 = (((((( &
                    -0.5014415D-03*t2 &
                    + 0.76771853D-02)*t2 &
                    - 0.0709253492D+00)*t2 &
                   + 0.4443584263D+00)*t2 &
                  - 1.7777560599D+00)*t2 &
                 + 3.9999973021D+00)*t2 &
                - 3.9999998721D+00)*t2 &
               + 1.0D+00

         bj1 = t*((((((( &
                       -0.1289769D-03*t2 &
                       + 0.22069155D-02)*t2 &
                       - 0.0236616773D+00)*t2 &
                      + 0.1777582922D+00)*t2 &
                     - 0.8888839649D+00)*t2 &
                    + 2.6666660544D+00)*t2 &
                   - 3.9999999710D+00)*t2 &
                  + 1.9999999998D+00)

         by0 = ((((((( &
                     -0.567433D-04*t2 &
                     + 0.859977D-03)*t2 &
                     - 0.94855882D-02)*t2 &
                    + 0.0772975809D+00)*t2 &
                   - 0.4261737419D+00)*t2 &
                  + 1.4216421221D+00)*t2 &
                 - 2.3498519931D+00)*t2 &
                + 1.0766115157D+00)*t2 &
               + 0.3674669052D+00

         by0 = 2.0D+00/pi*log(x/2.0D+00)*bj0 + by0

         by1 = (((((((( &
                      0.6535773D-03*t2 &
                      - 0.0108175626D+00)*t2 &
                      + 0.107657606D+00)*t2 &
                     - 0.7268945577D+00)*t2 &
                    + 3.1261399273D+00)*t2 &
                   - 7.3980241381D+00)*t2 &
                  + 6.8529236342D+00)*t2 &
                 + 0.3932562018D+00)*t2 &
                - 0.6366197726D+00)/x

         by1 = 2.0D+00/pi*log(x/2.0D+00)*bj1 + by1

      else

         t = 4.0D+00/x
         t2 = t*t
         a0 = sqrt(2.0D+00/(pi*x))

         p0 = (((( &
                 -0.9285D-05*t2 &
                 + 0.43506D-04)*t2 &
                 - 0.122226D-03)*t2 &
                + 0.434725D-03)*t2 &
               - 0.4394275D-02)*t2 &
              + 0.999999997D+00

         q0 = t*((((( &
                    0.8099D-05*t2 &
                    - 0.35614D-04)*t2 &
                    + 0.85844D-04)*t2 &
                   - 0.218024D-03)*t2 &
                  + 0.1144106D-02)*t2 &
                 - 0.031249995D+00)

         ta0 = x - 0.25D+00*pi
         bj0 = a0*(p0*cos(ta0) - q0*sin(ta0))
         by0 = a0*(p0*sin(ta0) + q0*cos(ta0))

         p1 = (((( &
                 0.10632D-04*t2 &
                 - 0.50363D-04)*t2 &
                 + 0.145575D-03)*t2 &
                - 0.559487D-03)*t2 &
               + 0.7323931D-02)*t2 &
              + 1.000000004D+00

         q1 = t*((((( &
                    -0.9173D-05*t2 &
                    + 0.40658D-04)*t2 &
                    - 0.99941D-04)*t2 &
                   + 0.266891D-03)*t2 &
                  - 0.1601836D-02)*t2 &
                 + 0.093749994D+00)

         ta1 = x - 0.75D+00*pi
         bj1 = a0*(p1*cos(ta1) - q1*sin(ta1))
         by1 = a0*(p1*sin(ta1) + q1*cos(ta1))

      end if

      dj0 = -bj1
      dj1 = bj0 - bj1/x
      dy0 = -by1
      dy1 = by0 - by1/x

      return
   end

   subroutine lpmn(mm, m, n, x, pm, pd)
!*****************************************************************************80
!
!! LPMN computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    19 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the leading dimension of PM and PD.
!    Input, integer ( kind = 4 ) M, the order of Pmn(x).
!    Input, integer ( kind = 4 ) N, the degree of Pmn(x).
!    Input, real ( kind = 8 ) X, the argument of Pmn(x).
!    Output, real ( kind = 8 ) PM(0:MM,0:N), PD(0:MM,0:N), the
!    values of Pmn(x) and Pmn'(x).

      implicit none
      integer(kind=4) mm
      integer(kind=4) n

      integer(kind=4) i
      integer(kind=4) j
      integer(kind=4) ls
      integer(kind=4) m
      real(kind=8) pd(0:mm, 0:n)
      real(kind=8) pm(0:mm, 0:n)
      real(kind=8) x
      real(kind=8) xq
      real(kind=8) xs

      do i = 0, n
         do j = 0, m
            pm(j, i) = 0.0D+00
            pd(j, i) = 0.0D+00
         end do
      end do

      pm(0, 0) = 1.0D+00

      if (abs(x) == 1.0D+00) then

         do i = 1, n
            pm(0, i) = x**i
            pd(0, i) = 0.5D+00*i*(i + 1.0D+00)*x**(i + 1)
         end do

         do j = 1, n
            do i = 1, m
               if (i == 1) then
                  pd(i, j) = 1.0D+300
               else if (i == 2) then
                  pd(i, j) = -0.25D+00*(j + 2)*(j + 1)*j &
                             *(j - 1)*x**(j + 1)
               end if
            end do
         end do

         return

      end if

      if (1.0D+00 < abs(x)) then
         ls = -1
      else
         ls = +1
      end if

      xq = sqrt(ls*(1.0D+00-x*x))
      xs = ls*(1.0D+00-x*x)
      do i = 1, m
         pm(i, i) = -ls*(2.0D+00*i - 1.0D+00)*xq*pm(i - 1, i - 1)
      end do

      do i = 0, m
         pm(i, i + 1) = (2.0D+00*i + 1.0D+00)*x*pm(i, i)
      end do

      do i = 0, m
         do j = i + 2, n
            pm(i, j) = ((2.0D+00*j - 1.0D+00)*x*pm(i, j - 1) - &
                        (i + j - 1.0D+00)*pm(i, j - 2))/(j - i)
         end do
      end do

      pd(0, 0) = 0.0D+00
      do j = 1, n
         pd(0, j) = ls*j*(pm(0, j - 1) - x*pm(0, j))/xs
      end do

      do i = 1, m
         do j = i, n
            pd(i, j) = ls*i*x*pm(i, j)/xs + (j + i) &
                       *(j - i + 1.0D+00)/xq*pm(i - 1, j)
         end do
      end do

      return
   end

!*****************************************************************************80
   
   subroutine elit ( hk, phi, fe, ee )
!! ELIT: complete and incomplete elliptic integrals F(k,phi) and E(k,phi).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    12 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) HK, the modulus, between 0 and 1.
!
!    Input, real ( kind = 8 ) PHI, the argument in degrees.
!
!    Output, real ( kind = 8 ) FE, EE, the values of F(k,phi) and E(k,phi).
!
      implicit none

      real ( kind = 8 ) a
      real ( kind = 8 ) a0
      real ( kind = 8 ) b
      real ( kind = 8 ) b0
      real ( kind = 8 ) c
      real ( kind = 8 ) ce
      real ( kind = 8 ) ck
      real ( kind = 8 ) d
      real ( kind = 8 ) d0
      real ( kind = 8 ) ee
      real ( kind = 8 ) fac
      real ( kind = 8 ) fe
      real ( kind = 8 ) g
      real ( kind = 8 ) hk
      integer ( kind = 4 ) n
      real ( kind = 8 ) phi
      real ( kind = 8 ) pi
      real ( kind = 8 ) r

      g = 0.0D+00
      pi = 3.14159265358979D+00
      a0 = 1.0D+00
      b0 = sqrt ( 1.0D+00 - hk * hk )
      d0 = ( pi / 180.0D+00 ) * phi
      r = hk * hk

      if ( hk == 1.0D+00 .and. phi == 90.0D+00 ) then

       fe = 1.0D+300
       ee = 1.0D+00

      else if ( hk == 1.0D+00 ) then

       fe = log ( ( 1.0D+00 + sin ( d0 ) ) / cos ( d0 ) )
       ee = sin ( d0 )

      else

       fac = 1.0D+00
       do n = 1, 40
         a = ( a0 + b0 ) /2.0D+00
         b = sqrt ( a0 * b0 )
         c = ( a0 - b0 ) / 2.0D+00
         fac = 2.0D+00 * fac
         r = r + fac * c * c
         if ( phi /= 90.0D+00 ) then
           d = d0 + atan ( ( b0 / a0 ) * tan ( d0 ) )
           g = g + c * sin( d )
           d0 = d + pi * int ( d / pi + 0.5D+00 )
         end if
         a0 = a
         b0 = b
         if ( c < 1.0D-07 ) then
           exit
         end if
       end do

       ck = pi / ( 2.0D+00 * a )
       ce = pi * ( 2.0D+00 - r ) / ( 4.0D+00 * a )
       if ( phi == 90.0D+00 ) then
         fe = ck
         ee = ce
       else
         fe = d / ( fac * a )
         ee = fe * ce / ck + g
       end if

      end if

      return
   end subroutine elit

!*****************************************************************************80

   subroutine jelp ( u, hk, esn, ecn, edn, eph )
!! JELP computes Jacobian elliptic functions SN(u), CN(u), DN(u).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U, the argument.
!
!    Input, real ( kind = 8 ) HK, the modulus, between 0 and 1.
!
!    Output, real ( kind = 8 ) ESN, ECN, EDN, EPH, the values of
!    sn(u), cn(u), dn(u), and phi (in degrees).
!
      implicit none

      real ( kind = 8 ) a
      real ( kind = 8 ) a0
      real ( kind = 8 ) b
      real ( kind = 8 ) b0
      real ( kind = 8 ) c
      real ( kind = 8 ) d
      real ( kind = 8 ) dn
      real ( kind = 8 ) ecn
      real ( kind = 8 ) edn
      real ( kind = 8 ) eph
      real ( kind = 8 ) esn
      real ( kind = 8 ) hk
      integer ( kind = 4 ) j
      integer ( kind = 4 ) n 
      real ( kind = 8 ) pi
      real ( kind = 8 ) r(40)
      real ( kind = 8 ) sa
      real ( kind = 8 ) t
      real ( kind = 8 ) u

      pi = 3.14159265358979D+00
      a0 = 1.0D+00
      b0 = sqrt ( 1.0D+00 - hk * hk )

      do n = 1, 40

       a = ( a0 + b0 ) / 2.0D+00
       b = sqrt ( a0 * b0 )
       c = ( a0 - b0 ) / 2.0D+00
       r(n) = c / a

       if ( c < 1.0D-07 ) then
         exit
       end if

       a0 = a
       b0 = b

      end do

      dn = 2.0D+00 ** n * a * u

      do j = n, 1, -1
       t = r(j) * sin ( dn )
       sa = atan ( t / sqrt ( abs ( 1.0D+00 - t * t )))
       d = 0.5D+00 * ( dn + sa )
       dn = d
      end do

      eph = d * 180.0D+00 / pi
      esn = sin ( d )
      ecn = cos ( d )
      edn = sqrt ( 1.0D+00 - hk * hk * esn * esn )

      return
   end

end module sfunctions
