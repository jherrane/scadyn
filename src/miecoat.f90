module miecoat
   use common
   use sfunctions
   implicit none

contains

!****************************************************************************80
!
! Bohren & Huffman's Mie coefficients for coated sphere
! All bessel functions computed by upward recurrence.
! Input:
!        X = 2*PI*RCORE*REFMED/WAVEL
!        Y = 2*PI*RMANT*REFMED/WAVEL
!        RFREL1 = REFCOR/REFMED
!        RFREL2 = REFMAN/REFMED
! where  REFCOR = complex refr.index of core)
!        REFMAN = complex refr.index of mantle)
!        REFMED = real refr.index of medium)
!        RCORE = radius of core
!        RMANT = radius of mantle
!        WAVEL = wavelength of light in ambient medium
   subroutine bhcoat(nstop, xx, yy, rrfrl1, rrfrl2, a_n, b_n)
      real(dp) :: xx, yy
      complex(dp) :: rrfrl1, rrfrl2
      integer :: iflag, n, nstop
      real(dp) :: chi0y, chi1y, chiy, del, psi0y, psi1y, psiy, x, y
      complex(dp) :: a_n(nstop), b_n(nstop)
      complex(dp) :: amess1, amess2, amess3, amess4, an, an1, ancap, &
                     bn, bn1, bncap, brack, &
                     chi0x2, chi0y2, chi1x2, chi1y2, chix2, chipx2, chipy2, chiy2, crack, &
                     d0x1, d0x2, d0y2, d1x1, d1x2, d1y2, dnbar, gnbar, ii, &
                     refrel, rfrel1, rfrel2, &
                     xi0y, xi1y, xiy, &
                     x1, x2, y2

!*** del is the inner sphere convergence criterion
      del = 1d-8
      ii = dcmplx(0d0, 1d0)
      x = xx
      y = yy
      rfrel1 = rrfrl1
      rfrel2 = rrfrl2

      x1 = rfrel1*x
      x2 = rfrel2*x
      y2 = rfrel2*y

      refrel = rfrel2/rfrel1

! Series truncation rule replaced
!ystop = y + 4.*y**0.3333 + 2.0
!nstop = ystop

      d0x1 = cos(x1)/sin(x1)
      d0x2 = cos(x2)/sin(x2)
      d0y2 = cos(y2)/sin(y2)

      psi0y = cos(y)
      psi1y = sin(y)
      chi0y = -sin(y)
      chi1y = cos(y)

      xi0y = psi0y-II*chi0y
      xi1y = psi1y-II*chi1y

      chi0y2 = -sin(y2)
      chi1y2 = cos(y2)
      chi0x2 = -sin(x2)
      chi1x2 = cos(x2)

! Inner sphere is done when flag = 1
      iflag = 0

      do n = 1, nstop
         psiy = (2d0*n - 1d0)*psi1y/y - psi0y
         chiy = (2d0*n - 1d0)*chi1y/y - chi0y
         xiy = psiy - II*chiy
         d1y2 = 1d0/(n/y2 - d0y2) - n/y2

         if (iflag /= 1) then
            d1x1 = 1d0/(n/x1 - d0x1) - n/x1
            d1x2 = 1d0/(n/x2 - d0x2) - n/x2
            chix2 = (2d0*n - 1d0)*chi1x2/x2 - chi0x2
            chiy2 = (2d0*n - 1d0)*chi1y2/y2 - chi0y2
            chipx2 = chi1x2 - n*chix2/x2
            chipy2 = chi1y2 - n*chiy2/y2
            ancap = refrel*d1x1 - d1x2
            ancap = ancap/(refrel*d1x1*chix2 - chipx2)
            ancap = ancap/(chix2*d1x2 - chipx2)
            brack = ancap*(chiy2*d1y2 - chipy2)
            bncap = refrel*d1x2 - d1x1
            bncap = bncap/(refrel*chipx2 - d1x1*chix2)
            bncap = bncap/(chix2*d1x2 - chipx2)
            crack = bncap*(chiy2*d1y2 - chipy2)
            amess1 = brack*chipy2
            amess2 = brack*chiy2
            amess3 = crack*chipy2
            amess4 = crack*chiy2

            if (ABS(amess1) <= del*ABS(d1y2) .OR. ABS(amess2) <= del &
                .OR. ABS(amess3) <= del*ABS(d1y2) .OR. ABS(amess4) <= del) then
               brack = (0d0, 0d0)
               crack = (0d0, 0d0)
               iflag = 1
            end if
         end if

         dnbar = d1y2 - brack*chipy2
         dnbar = dnbar/(1d0 - brack*chiy2)
         gnbar = d1y2 - crack*chipy2
         gnbar = gnbar/(1d0 - crack*chiy2)

         ! Output rule as in bhmie
         if (n > 1) then
            an1 = an
            bn1 = bn
            a_n(n - 1) = -an1
            b_n(n - 1) = -bn1
         end if

         an = (dnbar/rfrel2 + n/y)*psiy - psi1y
         an = an/((dnbar/rfrel2 + n/y)*xiy - xi1y)
         bn = (rfrel2*gnbar + n/y)*psiy - psi1y
         bn = bn/((rfrel2*gnbar + n/y)*xiy - xi1y)

         psi0y = psi1y
         psi1y = psiy
         chi0y = chi1y
         chi1y = chiy
         xi1y = psi1y-II*chi1y
         chi0x2 = chi1x2
         chi1x2 = chix2
         chi0y2 = chi1y2
         chi1y2 = chiy2
         d0x1 = d1x1
         d0x2 = d1x2
         d0y2 = d1y2
      end do

   end subroutine bhcoat

!****************************************************************************80
!* Mie scattering routine for coated spheres from BH
   subroutine coated_mie_coeff_nm(N, x, mr, a_nm, b_nm, c_nm, d_nm)
      integer :: N
      real(dp) :: x(2)
      complex(dp) :: mr(2)
      complex(dp), dimension((N + 1)**2 - 1) :: a_nm, b_nm, c_nm, d_nm

      complex(dp), dimension(N) :: a_n, b_n
      integer :: las, i1, m

      call BHCOAT(N + 1, x(1), x(2), mr(1), mr(2), a_n, b_n)

      las = 0; 
      do i1 = 1, N
         do m = -i1, i1
            las = las + 1; 
            a_nm(las) = a_n(i1); 
            b_nm(las) = b_n(i1); 
            c_nm(las) = dcmplx(0d0) !c_n(i1);
            d_nm(las) = dcmplx(0d0) !d_n(i1);
         end do
      end do

   end subroutine coated_mie_coeff_nm
end module miecoat
