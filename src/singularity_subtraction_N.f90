module singularity_subtraction_N
! Copyright (c) 2018 Johannes Markkanen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common
   implicit none

! Number of terms to be extracted
!integer, parameter :: Np = 1

contains

!****************************************************************************80
! int_S N 1/R + R + R³ +...   dS,  S is a triangle, N is a shape function
   function analytical_integration_N(tri_coord, n_vec, r, Np) result(intN)
      real(dp), intent(in) :: tri_coord(3, 3), n_vec(3), r(3)
      integer :: Np
      real(dp) :: intN(3, Np)
      real(dp) :: Kl(Np + 1), rmp1(3), slv(3), sl(3), ml(3), pl1mr(3), &
                  pl2mpl1(3), rmplm1(3), pl1mpl0(3), tl(3), Il(3), mvek(3, 3)
      real(dp) :: w0, normslv, tl0, ll, sml, spl, Rpl, Rml, R0l, Bl, &
                  Km3, hl, Inl(3, Np + 1), v1(3, Np + 1), d1, d2, d3, hhl(3)
      integer :: ind(7), nn, m, n, l2, l1, l0, lm1, l, i1, i2, nm
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

            hl = dot_product(pl1mpl0, ml)

            !---------------------------------------------------------------!

            if (n == -1) then
               mvek(:, l) = ml
               tl(l) = tl0
               hhl(l) = hl

               if (abs(w0) > 1.0d-10) then ! Check this

                  Bl = atan((tl0*spl)/(R0l*R0l + abs(w0)*Rpl)) - atan((tl0*sml)/(R0l*R0l + abs(w0)*Rml)); 
                  Km3 = 1.0/abs(w0)*Bl

               else

                  Km3 = 0.0
                  Il(l) = 0.0
               end if

               !if(abs(Rml + sml) > abs(Rpl-spl)) then
               !   Il(l) = log(((Rpl+spl)/(Rml+sml)))
               !else
               !   Il(l) = log(((Rml-sml)/(Rpl-spl)))
               !end if

               if (abs(Rml + sml) > 0.0) then
                  Il(l) = log(((Rpl + spl)/(Rml + sml)))
               else
                  Il(l) = 0.0
               end if

               Kl(nn) = Kl(nn) + 1/(dble(n) + 2.0)*(dble(n)*w0*w0*Km3 + tl0*Il(l))

            else !----------------------------------------------------!

               Il(l) = 1.0/(dble(n) + 1.0)*(spl*Rpl**dble(n) - sml*Rml**dble(n) + dble(n)*R0l*R0l*Il(l))

            end if

         end do

         if (n .ne. -1) then
            Kl(nn) = 1.0/(dble(n) + 2.0)*(dble(n)*w0*w0*Kl(nn - 1) + (tl(1)*Il(1) + tl(2)*Il(2) + tl(3)*Il(3)))
         end if

         Inl(:, nn) = Il
         nn = nn + 1
      end do

      v1(:, :) = 0.0
      do i1 = 1, Np + 1
         do i2 = 1, 3
            v1(:, i1) = v1(:, i1) + Inl(i2, i1)*mvek(:, i2)
         end do
      end do
!print*,v1
      do nm = 1, Np

         d1 = dot_product(mvek(:, 1), v1(:, nm + 1))
         d2 = dot_product(mvek(:, 2), v1(:, nm + 1))
         d3 = dot_product(mvek(:, 3), v1(:, nm + 1))

         intN(1, nm) = 1.0/hhl(1)*(tl(1)*Kl(nm) - 1.0/(2.0*dble(nm) - 1.0)*d1)
         intN(2, nm) = 1.0/hhl(2)*(tl(2)*Kl(nm) - 1.0/(2.0*dble(nm) - 1.0)*d2)
         intN(3, nm) = 1.0/hhl(3)*(tl(3)*Kl(nm) - 1.0/(2.0*dble(nm) - 1.0)*d3)
      end do

!print*,intN

   end function analytical_integration_N

!****************************************************************************80
! Computes integral I(rf) = int_V G(rp, rf) drp
   function singularity_subtraction_int_V_GN(rf, k, Pb_tet, wb_tet, B_coord, shape_tet) result(intN)
!type (element_data), intent(in) :: elem
      real(dp), intent(in) :: rf(3), k
      real(dp), dimension(:) :: wb_tet
      real(dp), dimension(:, :) :: Pb_tet, B_coord, shape_tet

      complex(dp) :: intN(4)
      integer, parameter :: Np = 1

      real(dp), parameter :: fact(6) = [1.0, -0.5, 1.0/24, -1.0/720, 1.0/40320, -1.0/3628800]
      real(dp), parameter :: ak(5) = [1.0/2, 1.0/4, 1.0/6, 1.0/8, 1.0/10]
      integer, parameter :: ind(12) = [2, 3, 4, 1, 4, 3, 1, 2, 4, 1, 3, 2]
      real(dp), parameter :: ak2(5) = [1.0/1, 1.0/3, 1.0/5, 1.0/7, 1.0/9]

      complex(dp) :: intR(4), II, G0
      real(dp) :: intR2(3, Np + 1), intR4(Np), intR6(4, Np), n_dN(4), intR5(4, Np)

      integer :: i1, l, f1, l2, j1, N_ind(3), j2
      real(dp) :: rp(3), RR, tet_n(3, 4), tri_coord(3, 3), tri_n(3), rr2(3), cut

      real(dp) :: f, a1, f2(4), rot(3, 6), dN(3, 4)

      intR(:) = dcmplx(0.0, 0.0)

!*****************************  Numerical part   ******************************
      cut = 1.0d-10*sum(wb_tet)**(1.0/3.0)

      do i1 = 1, size(wb_tet)
         rp = Pb_tet(:, i1)
         RR = norm(rp, rf)

         !if(RR < 1.0d-10) then ! Check this
         if (RR < cut) then
            II = dcmplx(0.0, k/(4.0*pi))

         else

            f = 0.0
            do l = 1, Np
               f = f + 1.0/(4*pi)*k**(2.0*(l - 1.0))*fact(l)*RR**(2.0*(l - 1.0) - 1.0)
            end do

            G0 = Gr(rf, rp, k)
            II = G0 - dcmplx(f, 0)

         end if

         intR = intR + II*shape_tet(:, i1)*wb_tet(i1)
      end do

!****************************  Analytical part  *********************************

      intR4(:) = 0.0
      intR5(:, :) = 0.0
      intR6(:, :) = 0.0
      tet_n = tetra_n_vectors(B_coord)

      do f1 = 1, 4

         N_ind(1) = ind(3*(f1 - 1) + 1)
         N_ind(2) = ind(3*(f1 - 1) + 2)
         N_ind(3) = ind(3*(f1 - 1) + 3)

         tri_coord = tetra_face(B_coord, f1)
         tri_n = tet_n(:, f1)

         call gradshape(rot, dN, B_coord)

         intR2 = analytical_integration_N(tri_coord, tri_n, rf, Np + 1)
         !print*,intR2(1,1)
         rr2 = tri_coord(:, 1) - rf
         a1 = dot_product(tri_n, rr2)

         do j1 = 1, Np
            intR4(j1) = intR4(j1) + a1*sum(intR2(:, j1))*ak(j1)

            intR6(N_ind, j1) = intR6(N_ind, j1) + intR2(:, j1)*a1*ak(j1)

         end do

         n_dN(1) = dot_product(tri_n, dN(:, 1))
         n_dN(2) = dot_product(tri_n, dN(:, 2))
         n_dN(3) = dot_product(tri_n, dN(:, 3))
         n_dN(4) = dot_product(tri_n, dN(:, 4))

         do j2 = 1, Np
            intR5(:, j2) = intR5(:, j2) + n_dN*sum(intR2(:, j2 + 1))*ak(j2)*ak2(j2)
         end do

      end do

      f2(:) = 0.0
      do l2 = 1, Np
         f2 = f2 + 1.0/(4*pi)*k**(2.0*(l2 - 1))*fact(l2)*(-intR5(:, l2) + intR6(:, l2))
      end do

      intN = intR + dcmplx(f2, 0)

   end function singularity_subtraction_int_V_GN

!********************************************************************************
! Computes int(rf) = int_S GN dS analytically
   function singularity_subtraction_int_S_GN(rf, tri_coord, tri_n, k, Pb_tri, wb_tri, tri_shape) result(intN)
!type (element_data), intent(in) :: elem

      real(dp), intent(in) :: rf(3), tri_n(3), k
      real(dp), intent(in) :: tri_coord(3, 3)
      real(dp), dimension(:, :), intent(in) :: Pb_tri, tri_shape
      real(dp), dimension(:), intent(in) :: wb_tri

      integer, parameter :: Np = 1

      real(dp), parameter :: fact(6) = [1.0, -0.5, 1.0/24, -1.0/720, 1.0/40320, -1.0/3628800]

      real(dp) :: rp(3), R, f, intR2(3, Np), f2(3), cut
      complex(dp) :: intN(3), intR(3), II, G0
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
            f = 0.0
            do l = 1, Np
               f = f + 1/(4*pi)*k**(2.0*(dble(l) - 1.0))*fact(l)*R**(2.0*(dble(l) - 1) - 1)
            end do
            G0 = Gr(rf, rp, k)
            II = G0 - dcmplx(f, 0)
         end if

         intR = intR + II*wb_tri(i1)*tri_shape(:, i1)

      end do

!________________Analytical part_______________________________!

      intR2 = analytical_integration_N(tri_coord, tri_n, rf, Np)

      f2(:) = 0.0

      do l2 = 1, Np
         f2 = f2 + 1/(4*pi)*k**(2.0*(dble(l2) - 1.0))*fact(l2)*(intR2(:, l2))
      end do

      intN = intR + dcmplx(f2, 0)

   end function singularity_subtraction_int_S_GN

end module singularity_subtraction_N
