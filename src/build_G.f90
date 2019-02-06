subroutine build_G(matrices2, mesh2)
! Copyright (c) 2018 Johannes Markkanen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common, matrices => matrices, mesh => mesh
   use possu
   implicit none

   type(mesh_struct) :: mesh2
   type(data) :: matrices2

   double complex, dimension(:, :, :), allocatable :: A
   double complex :: cmp
   integer :: m1, l1, k1, las
   double precision :: R

   allocate (A(2*mesh%Nx, 2*mesh%Ny, 2*Mesh%Nz))
   allocate (matrices%Fg(2*mesh%Nx, 2*mesh%Ny, 2*Mesh%Nz))
   A(:, :, :) = dcmplx(0.0, 0.0)
   las = 1

   do m1 = 1, mesh%Nz
      do l1 = 1, mesh%Ny
         do k1 = 1, mesh%Nx
            R = norm(mesh%nodes(:, 1), mesh%nodes(:, las))
            if (R > 0d0) then
               cmp = cdexp(dcmplx(0.0, mesh%k*R))/(4*PI*R)
            else
               cmp = dcmplx(0.0, 0.0)
            end if

            A(k1, l1, m1) = cmp

            if (k1 .ne. 1) then
               A(2*mesh%Nx - k1 + 2, l1, m1) = cmp
            end if
            if (k1 .ne. 1 .and. m1 .ne. 1) then
               A(2*mesh%Nx - k1 + 2, l1, 2*mesh%Nz - m1 + 2) = cmp
            end if
            if (k1 .ne. 1 .and. l1 .ne. 1) then
               A(2*mesh%Nx - k1 + 2, 2*mesh%Ny - l1 + 2, m1) = cmp
            end if
            if (k1 .ne. 1 .and. l1 .ne. 1 .and. m1 .ne. 1) then
               A(2*mesh%Nx - k1 + 2, 2*mesh%Ny - l1 + 2, 2*mesh%Nz - m1 + 2) = cmp
            end if
            if (m1 .ne. 1) then
               A(k1, l1, 2*mesh%Nz - m1 + 2) = cmp
            end if
            if (l1 .ne. 1) then
               A(k1, 2*mesh%Ny - l1 + 2, m1) = cmp
            end if
            if (l1 .ne. 1 .and. m1 .ne. 1) then
               A(k1, 2*mesh%Ny - l1 + 2, 2*mesh%Nz - m1 + 2) = cmp
            end if

            las = las + 1

         end do
      end do
   end do
   A(1, 1, 1) = dcmplx(0.0, 0.0)

   call FFT3d(A)
   matrices%Fg = A

end subroutine
