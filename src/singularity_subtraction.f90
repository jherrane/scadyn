module singularity_subtraction
! Copyright (c) 2018 Johannes Markkanen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common
   implicit none

! Number of terms to be extracted
   integer, parameter :: Np = 2

contains

!***********************************************************************
! int_S N 1/R + R + R³ +...   dS,  S is a triangle, N is a shape function
   function analytical_integration(tri_coord, n_vec, r) result(intN)
      real(dp), intent(in) :: tri_coord(3, 3), n_vec(3), r(3)
      real(dp) :: intN(Np)

      real(dp) :: Kl(Np + 1), rmp1(3), slv(3), sl(3), ml(3), pl1mr(3), pl2mpl1(3), rmplm1(3), pl1mpl0(3), tl(3), Il(3)
      real(dp) :: w0, normslv, tl0, ll, sml, spl, Rpl, Rml, R0l, Bl, Km3
      integer :: ind(7), nn, m, n, l2, l1, l0, lm1, l
      ind = [3, 1, 2, 3, 1, 2, 3]

      Kl(:) = 0
      nn = 1

      rmp1 = r - tri_coord(:, 1)
      w0 = dot_product(rmp1, n_vec)

      do m = 0, Np
         n = 2*m - 1

         do l = 1, 3
            l2 = ind(l + 3)
            l1 = ind(l + 2)
            l0 = ind(l + 1)
            lm1 = ind(l)

            slv = tri_coord(:, l2) - tri_coord(:, l1)

            normslv = sqrt(dot_product(slv, slv))
            sl = slv/normslv
            ml = crossRR(sl, n_vec)

            pl1mr = tri_coord(:, l1) - r
            tl0 = dot_product(pl1mr, ml)

            pl2mpl1 = tri_coord(:, l2) - tri_coord(:, l1)

            ll = sqrt(dot_product(pl2mpl1, pl2mpl1))

            sml = dot_product(pl1mr, sl)
            spl = sml + ll

            rmplm1 = r - tri_coord(:, lm1)

            Rpl = sqrt(dot_product(rmplm1, rmplm1))
            Rml = sqrt(dot_product(pl1mr, pl1mr))

            R0l = sqrt(tl0*tl0 + w0*w0)

            pl1mpl0 = tri_coord(:, l1) - tri_coord(:, l0)

            !hl = dot_product(pl1mpl0,ml)

            !---------------------------------------------------------------!

            if (n == -1) then
               !mvek(:,l) = ml
               tl(l) = tl0
               !hhl(l) =hl

               if (abs(w0) > 1.0d-10) then ! Check this

                  Bl = atan((tl0*spl)/(R0l*R0l + abs(w0)*Rpl)) - atan((tl0*sml)/(R0l*R0l + abs(w0)*Rml)); 
                  Km3 = 1.0/abs(w0)*Bl

               else

                  Km3 = 0
                  Il(l) = 0
               end if

               if ((Rml + sml) > (Rpl - spl)) then
                  Il(l) = log(abs((Rpl + spl)/(Rml + sml)))
               else
                  Il(l) = log(abs((Rml - sml)/(Rpl - spl)))
               end if

               Kl(nn) = Kl(nn) + 1/(dble(n) + 2.0)*(dble(n)*w0*w0*Km3 + tl0*Il(l))

            else !----------------------------------------------------!

               Il(l) = 1/(dble(n) + 1.0)*(spl*Rpl**dble(n) - sml*Rml**dble(n) + dble(n)*R0l*R0l*Il(l))
            end if

         end do

         if (n .ne. -1) then
            Kl(nn) = 1/(dble(n) + 2.0)*(dble(n)*w0*w0*Kl(nn - 1) + (tl(1)*Il(1) + tl(2)*Il(2) + tl(3)*Il(3)))
         end if

         nn = nn + 1
      end do

      intN(:) = Kl(1:Np)

   end function analytical_integration

!****************************************************************************80
! Computes integral I(rf) = int_V G(rp, rf) drp
   function singularity_subtraction_int_V_G(rf, k, Pb_tet, wb_tet, B_coord) result(intN)
!type (element_data), intent(in) :: elem
      real(dp), intent(in) :: rf(3), k
      real(dp), dimension(:) :: wb_tet
      real(dp), dimension(:, :) :: Pb_tet, B_coord

      complex(dp) :: intN

      real(dp), parameter :: fact(6) = [1.0, -0.5, 1.0/24, -1.0/720, 1.0/40320, -1.0/3628800]
      real(dp), parameter :: ak(5) = [1.0/2, 1.0/4, 1.0/6, 1.0/8, 1.0/10]

      complex(dp) :: intR, II, G0
      real(dp) :: intR2(Np), intR4(Np)

      integer :: i1, l, f1, l2, j1
      real(dp) :: rp(3), RR, tet_n(3, 4), tri_coord(3, 3), tri_n(3), rr2(3), cut

      real(dp) :: f, a1, f2

      intR = dcmplx(0.0, 0.0)

!*****************************  Numerical part   ******************************
      cut = 1.0d-10*sum(wb_tet)**(1.0/3.0)

      do i1 = 1, size(wb_tet)
         rp = Pb_tet(:, i1)
         RR = norm(rp, rf)

         !if(RR < 1.0d-10) then ! Check this
         if (RR < cut) then
            II = dcmplx(0.0, k/(4.0*pi))

         else

            f = 0
            do l = 1, Np
               f = f + 1.0/(4*pi)*k**(2.0*(l - 1.0))*fact(l)*RR**(2.0*(l - 1.0) - 1.0)
            end do

            G0 = Gr(rf, rp, k)
            II = G0 - dcmplx(f, 0)

         end if

         intR = intR + II*wb_tet(i1)
      end do

!****************************  Analytical part  *********************************

      intR4(:) = 0

      tet_n = tetra_n_vectors(B_coord)

      do f1 = 1, 4

         tri_coord = tetra_face(B_coord, f1)
         tri_n = tet_n(:, f1)

         intR2 = analytical_integration(tri_coord, tri_n, rf)
         !print*,intR2(1,1)
         rr2 = tri_coord(:, 1) - rf
         a1 = dot_product(tri_n, rr2)

         do j1 = 1, Np
            intR4(j1) = intR4(j1) + a1*intR2(j1)*ak(j1)
         end do

      end do
!print*,intR4
      f2 = 0
      do l2 = 1, Np
         f2 = f2 + 1.0/(4*pi)*k**(2.0*(l2 - 1))*fact(l2)*intR4(l2)
      end do

      intN = intR + dcmplx(f2, 0)
!print*,f2

   end function singularity_subtraction_int_V_G

!****************************************************************************80
! Computes int(rf) = int_S GN dS anatytically
   function singularity_subtraction_int_S_G(rf, tri_coord, tri_n, k, Pb_tri, wb_tri) result(intN)
!type (element_data), intent(in) :: elem

      real(dp), intent(in) :: rf(3), tri_n(3), k
      real(dp), intent(in) :: tri_coord(3, 3)
      real(dp), dimension(:, :), intent(in) :: Pb_tri
      real(dp), dimension(:), intent(in) :: wb_tri

      real(dp), parameter :: fact(6) = [1.0, -0.5, 1.0/24, -1.0/720, 1.0/40320, -1.0/3628800]

      real(dp) :: rp(3), R, f, intR2(Np), f2, cut
      complex(dp) :: intN, intR, II, G0
      integer :: i1, l, l2

!_____________________Numerical part____________________________________!

      intR = 0
      cut = 1.0d-10*sum(wb_tri)**(1.0/2.0)
      do i1 = 1, size(wb_tri)

         rp = Pb_tri(:, i1)
         R = norm(rp, rf)
         !if(R< 0.0000000001) then ! Check this
         if (R < cut) then
            II = dcmplx(0, k/(4*pi))
         else
            f = 0
            do l = 1, Np
               f = f + 1/(4*pi)*k**(2.0*(dble(l) - 1.0))*fact(l)*R**(2.0*(dble(l) - 1) - 1)
            end do
            G0 = Gr(rf, rp, k)
            II = G0 - dcmplx(f, 0)
         end if

         intR = intR + II*wb_tri(i1)

      end do

!________________Analytical part_______________________________!

      intR2 = analytical_integration(tri_coord, tri_n, rf)

      f2 = 0.0

      do l2 = 1, Np
         f2 = f2 + 1/(4*pi)*k**(2.0*(dble(l2) - 1.0))*fact(l2)*(intR2(l2))
      end do

      intN = intR + dcmplx(f2, 0)

   end function singularity_subtraction_int_S_G

!****************************************************************************80
! Computes integral I(rf) = int_V G(rp, rf) drp
   function singularity_subtraction_int_V_gradG(rf, k, Pb_tet, wb_tet, B_coord) result(intN)
!type (element_data), intent(in) :: elem
      real(dp), intent(in) :: rf(3), k
      real(dp), dimension(:) :: wb_tet
      real(dp), dimension(:, :) :: Pb_tet, B_coord

      complex(dp) :: intN(3)

      real(dp), parameter :: fact(6) = [1.0, -0.5, 1.0/24, -1.0/720, 1.0/40320, -1.0/3628800]
      real(dp), parameter :: ak(5) = [1.0/2, 1.0/4, 1.0/6, 1.0/8, 1.0/10]

      complex(dp) :: intR(3), II(3), G0(3)
      real(dp) :: intR2(Np), intR4(3, Np)

      integer :: i1, l, f1, l2, j1
      real(dp) :: rp(3), RR, tet_n(3, 4), tri_coord(3, 3), tri_n(3), rr2(3), cut
      real(dp) :: f, a1, f2(3)

      intR(:) = dcmplx(0.0, 0.0)

!*****************************  Numerical part   ******************************
      cut = 1.0d-10*sum(wb_tet)**(1.0/3.0)

      do i1 = 1, size(wb_tet)
         rp = Pb_tet(:, i1)
         RR = norm(rp, rf)
         if (RR < cut) then
            II(:) = dcmplx(0.0, 0.0)! Check this
         else
            f = 0.0
            do l = 1, Np
               f = f + 1.0/(4*pi)*k**(2.0*(l - 1.0))*fact(l)*dble(2*(l - 1) - 1)*RR**(2.0*(l - 1.0) - 2.0)
            end do

            G0 = grad_G(rf, rp, k)
            II = G0 - dcmplx(f, 0.0)*(rf - rp)/RR
         end if
         intR = intR - II*wb_tet(i1)
      end do

!****************************  Analytical part  *********************************

      intR4(:, :) = 0.0

      tet_n = tetra_n_vectors(B_coord)

      do f1 = 1, 4
         tri_coord = tetra_face(B_coord, f1)
         tri_n = tet_n(:, f1)

         intR2 = analytical_integration(tri_coord, tri_n, rf)

         do j1 = 1, Np
            intR4(:, j1) = intR4(:, j1) + tri_n*intR2(j1)
         end do
      end do

      f2(:) = 0.0
      do l2 = 1, Np
         f2 = f2 + 1.0/(4*pi)*k**(2.0*(l2 - 1))*fact(l2)*intR4(:, l2)
      end do

      intN = intR + dcmplx(f2, 0.0)

   end function singularity_subtraction_int_V_gradG

end module singularity_subtraction
