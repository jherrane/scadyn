module precorrection
! Copyright (c) 2018 Johannes Markkanen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common
   use integration_points
   use io
   use possu
   use geometry
   use sparse !lin
   use sparse_mat
   use singularity_subtraction_N
   use singularity_subtraction

   implicit none

contains

!****************************************************************************80

   function compute_pG(mesh, t_cube, b_cube) result(G)
      implicit none
      type(mesh_struct) :: mesh

      integer, intent(in) :: t_cube, b_cube

      double complex, dimension(:, :), allocatable :: G
      integer :: N, i1, i2
      double precision :: r(3), rp(3)

      N = size(mesh%elem_box, 1)

      allocate (G(N, N))
      G(:, :) = dcmplx(0.0, 0.0)

      do i1 = 1, N
         r = mesh%nodes(:, mesh%elem_box(i1, t_cube))
         do i2 = 1, N
            if (mesh%elem_box(i1, t_cube) .ne. mesh%elem_box(i2, b_cube)) then
               rp = mesh%nodes(:, mesh%elem_box(i2, b_cube))
               G(i1, i2) = Gr(r, rp, mesh%k)
            end if
         end do
      end do

   end function compute_pG

!****************************************************************************80

   subroutine precorrect_const(matrices, tet_T, tet_B, G, corr1, corr2)
      type(data), intent(in) :: matrices
      double complex, dimension(:, :), intent(in) :: G
      double complex, dimension(3, 3) :: corr1, corr2
      integer :: tet_T, tet_B
      double complex :: a, vec_x(size(G, 1)), vec_y(size(G, 1)), vec_z(size(G, 1))

      corr1(:, :) = dcmplx(0.0, 0.0)
      corr2(:, :) = dcmplx(0.0, 0.0)

      a = dot_product(conjg(matmul(G, matrices%S(:, tet_B))), matrices%S(:, tet_T))

      corr1(1, 1) = a
      corr1(2, 2) = a
      corr1(3, 3) = a

      vec_x = matmul(G, matrices%Sx(:, tet_B))
      vec_y = matmul(G, matrices%Sy(:, tet_B))
      vec_z = matmul(G, matrices%Sz(:, tet_B))

      corr2(1, 1) = dot_product(conjg(vec_x), matrices%Sx(:, tet_T))
      corr2(1, 2) = dot_product(conjg(vec_y), matrices%Sx(:, tet_T))
      corr2(1, 3) = dot_product(conjg(vec_z), matrices%Sx(:, tet_T))

      corr2(2, 1) = dot_product(conjg(vec_x), matrices%Sy(:, tet_T))
      corr2(2, 2) = dot_product(conjg(vec_y), matrices%Sy(:, tet_T))
      corr2(2, 3) = dot_product(conjg(vec_z), matrices%Sy(:, tet_T))

      corr2(3, 1) = dot_product(conjg(vec_x), matrices%Sz(:, tet_T))
      corr2(3, 2) = dot_product(conjg(vec_y), matrices%Sz(:, tet_T))
      corr2(3, 3) = dot_product(conjg(vec_z), matrices%Sz(:, tet_T))

   end subroutine precorrect_const

!****************************************************************************80

   subroutine precorrect_lin(matrices, tet_T, tet_B, G, corr1, corr2)
      type(data), intent(in) :: matrices
      double complex, dimension(:, :), intent(in) :: G
      double complex, dimension(12, 12) :: corr1, corr2
      integer :: tet_T, tet_B, i1, i2
      double complex :: a, vec_x(size(G, 1)), vec_y(size(G, 1)), vec_z(size(G, 1))

      corr1(:, :) = dcmplx(0.0, 0.0)
      corr2(:, :) = dcmplx(0.0, 0.0)

      do i1 = 1, 4
         do i2 = 1, 4
            a = dot_product(conjg(matmul(G, matrices%S(:, 4*(tet_B - 1) + i1))), matrices%S(:, 4*(tet_T - 1) + i2))
            corr1(i2, i1) = a
            corr1(4 + i2, 4 + i1) = a
            corr1(8 + i2, 8 + i1) = a
         end do
      end do

      do i1 = 1, 4

         vec_x = matmul(G, matrices%Sx(:, 4*(tet_B - 1) + i1))
         vec_y = matmul(G, matrices%Sy(:, 4*(tet_B - 1) + i1))
         vec_z = matmul(G, matrices%Sz(:, 4*(tet_B - 1) + i1))

         do i2 = 1, 4
            corr2(i2, i1) = dot_product(conjg(vec_x), matrices%Sx(:, 4*(tet_T - 1) + i2))
            corr2(i2, 4 + i1) = dot_product(conjg(vec_y), matrices%Sx(:, 4*(tet_T - 1) + i2))
            corr2(i2, 8 + i1) = dot_product(conjg(vec_z), matrices%Sx(:, 4*(tet_T - 1) + i2))

            corr2(4 + i2, i1) = dot_product(conjg(vec_x), matrices%Sy(:, 4*(tet_T - 1) + i2))
            corr2(4 + i2, 4 + i1) = dot_product(conjg(vec_y), matrices%Sy(:, 4*(tet_T - 1) + i2))
            corr2(4 + i2, 8 + i1) = dot_product(conjg(vec_z), matrices%Sy(:, 4*(tet_T - 1) + i2))

            corr2(8 + i2, i1) = dot_product(conjg(vec_x), matrices%Sz(:, 4*(tet_T - 1) + i2))
            corr2(8 + i2, 4 + i1) = dot_product(conjg(vec_y), matrices%Sz(:, 4*(tet_T - 1) + i2))
            corr2(8 + i2, 8 + i1) = dot_product(conjg(vec_z), matrices%Sz(:, 4*(tet_T - 1) + i2))

         end do
      end do

   end subroutine precorrect_lin

!****************************************************************************80

   subroutine compute_near_zone_interactions_const(matrices, mesh)
      type(mesh_struct), intent(in) :: mesh
      type(data) :: matrices

      integer(kind=8) :: mem
      integer :: t_cube, b_cube, i1, t1, b1, tet_T, tet_B
      integer :: T_nodes(4), B_nodes(4), bases(mesh%N_tet_cube), tests(mesh%N_tet_cube)
      integer, dimension(:), allocatable :: near_cubes
      double precision ::  vol_T, vol_B, B_coord(3, 4), T_coord(3, 4)
      double complex :: mat_block(3, 3), corr(3, 3), corr2(3, 3), ele_eps
      double complex, dimension(:, :), allocatable :: Gcorr
      double precision :: T_rot(3, 6), T_dN(3, 4), B_rot(3, 6), B_dN(3, 4)
      double precision, dimension(:, :), allocatable :: P0_tet, P0_tri
      double precision, dimension(:), allocatable :: w0_tet, w0_tri
      double complex :: alok(3, 3), alok1(3, 3), alok2(3, 3)

      integer :: blok, M_cube

      M_cube = (2*mesh%near_zone + 1)**3

      allocate (Gcorr(size(mesh%elem_box, 1), size(mesh%elem_box, 1)))

      allocate (near_cubes(M_cube))

      corr(:, :) = dcmplx(0.0, 0.0)
      corr2(:, :) = dcmplx(0.0, 0.0)

      call inttetra(P0_tet, w0_tet, 5)
      call inttri(P0_tri, w0_tri, 5)

      call allocate_Acorr(matrices%sp_mat, matrices%sp_ind, mesh)

      matrices%sp_mat(:, :, :) = cmplx(0.0, 0.0)
      matrices%sp_ind(:, :) = 0

      mem = sizeof(matrices%sp_mat)/1024/1024 + &
            sizeof(matrices%sp_ind)/1024/1024 + &
            sizeof(matrices%S)/1024/1024*4 + &
            sizeof(matrices%indS)/1024/1024 + &
            sizeof(matrices%Fg)/1024/1024*2 + &
            (5 + mesh%maxit)*8*2*3*size(mesh%elem, 2)/1024/1024 + &
            sizeof(mesh%node)/1024/1024 + &
            sizeof(mesh%elem)/1024/1024 + &
            sizeof(mesh%elem_box)/1024/1024 + &
            sizeof(mesh%nodes)/1024/1024 + &
            sizeof(mesh%tetras)/1024/1024

      print *, 'Estimated memory usage:', mem, 'Mb'

      blok = 1

      do t_cube = 1, mesh%N_cubes

         tests = mesh%tetras(:, t_cube)

         near_cubes = find_near_cubes(mesh, t_cube)

         do i1 = 1, M_cube
            if (tests(1) == 0) exit
            b_cube = near_cubes(i1)
            if (b_cube == 0) exit

            bases = mesh%tetras(:, b_cube)
            Gcorr = compute_pG(mesh, t_cube, b_cube)

            do t1 = 1, mesh%N_tet_cube
               tet_T = tests(t1)
               if (tet_T == 0) exit

               T_nodes = mesh%elem(:, tet_T)
               T_coord = mesh%node(:, T_nodes)
               vol_T = tetra_volume(T_coord)

               call gradshape(T_rot, T_dN, T_coord)

               ele_eps = mesh%eps(tet_T)

               do b1 = 1, mesh%N_tet_cube
                  tet_B = bases(b1)

                  if (tet_B == 0) exit
                  if (tet_B >= tet_T) then
                     B_nodes = mesh%elem(:, tet_B)
                     B_coord = mesh%node(:, B_nodes)
                     vol_B = tetra_volume(B_coord)

                     call gradshape(B_rot, B_dN, B_coord)

                     !_________-Correction term____________________________________

                     call precorrect_const(matrices, tet_T, tet_B, Gcorr, corr, corr2)
                     corr = corr*sqrt(vol_T*vol_B)
                     corr2 = corr2*sqrt(vol_T*vol_B)

                     !_____________________________________________________________

                     alok(:, :) = dcmplx(0.0, 0.0)
                     if (tet_T == tet_B) then
                        alok(1, 1) = 1.0/(ele_eps - 1.0)
                        alok(2, 2) = 1.0/(ele_eps - 1.0)
                        alok(3, 3) = 1.0/(ele_eps - 1.0)

                        !alok(1,1) = ele_eps / (ele_eps-1.0)
                        !alok(2,2) = ele_eps / (ele_eps-1.0)
                        !alok(3,3) = ele_eps / (ele_eps-1.0)
                     end if

                     call integrate_V_V_G(alok1, mesh, P0_tet, w0_tet, P0_tet, w0_tet, T_coord, B_coord)
                     call integrate_dV_dV_G(alok2, mesh, P0_tri, w0_tri, P0_tri, w0_tri, T_coord, B_coord)

                     !call integrate_nx_dV_dV_G(alok2,mesh,P0_tri,w0_tri,P0_tri,w0_tri ,T_coord,B_coord)

                     ! call integrate_dV_V_gradG(alok2,mesh,P0_tet,w0_tet,P0_tri,w0_tri ,T_coord,B_coord)

                     mat_block = alok + 1.0/sqrt(vol_B*vol_T)* &
                                 (-mesh%k**2*(alok1 - corr) + (alok2 - corr2))

                     !mat_block = alok + 1.0/sqrt(vol_B*vol_T) * &
!(-mesh%k**2 * (-corr) + (-alok2 - corr2) )

                     matrices%sp_mat(:, :, blok) = cmplx(mat_block)
                     matrices%sp_ind(1, blok) = tet_T
                     matrices%sp_ind(2, blok) = tet_B

                     blok = blok + 1

                  end if
               end do
            end do
         end do
         call print_bar(t_cube, mesh%N_cubes)
      end do

   end subroutine compute_near_zone_interactions_const

!****************************************************************************80

   subroutine compute_near_zone_interactions_lin(matrices, mesh)

      type(mesh_struct), intent(in) :: mesh
      type(data) :: matrices

      integer(kind=8) :: mem
      integer :: t_cube, b_cube, i1, t1, b1, tet_T, tet_B
      integer :: T_nodes(4), B_nodes(4), bases(mesh%N_tet_cube), tests(mesh%N_tet_cube)
      integer, dimension(:), allocatable :: near_cubes
      double precision ::  vol_T, vol_B, B_coord(3, 4), T_coord(3, 4)
      double complex :: mat_block(12, 12), corr(12, 12), corr2(12, 12), ele_eps
      double complex, dimension(:, :), allocatable :: Gcorr
      double precision :: T_rot(3, 6), T_dN(3, 4), B_rot(3, 6), B_dN(3, 4)
      double precision, dimension(:, :), allocatable :: P0_tet, P0_tri
      double precision, dimension(:), allocatable :: w0_tet, w0_tri
      double complex :: alok(12, 12), alok1(12, 12), alok2(12, 12), alok3(12, 12), alok4(12, 12), alok5(12, 12)
      double complex :: GNN(4, 4)
      integer :: ind(12, 2), blok, M_cube

      M_cube = (2*mesh%near_zone + 1)**3

      allocate (Gcorr(size(mesh%elem_box, 1), size(mesh%elem_box, 1)))

      allocate (near_cubes(M_cube))

      corr(:, :) = dcmplx(0.0, 0.0)
      corr2(:, :) = dcmplx(0.0, 0.0)

      call inttetra(P0_tet, w0_tet, 5)
      call inttri(P0_tri, w0_tri, 5)

      call allocate_Acorr(matrices%sp_mat, matrices%sp_ind, mesh)

      matrices%sp_mat(:, :, :) = cmplx(0.0, 0.0)
      matrices%sp_ind(:, :) = 0

      mem = sizeof(matrices%sp_mat)/1024/1024 + &
            sizeof(matrices%sp_ind)/1024/1024 + &
            sizeof(matrices%S)/1024/1024*4 + &
            sizeof(matrices%indS)/1024/1024 + &
            sizeof(matrices%Fg)/1024/1024*2 + &
            (5 + mesh%maxit)*8*2*12*size(mesh%elem, 2)/1024/1024 + &
            sizeof(mesh%node)/1024/1024 + &
            sizeof(mesh%elem)/1024/1024 + &
            sizeof(mesh%elem_box)/1024/1024 + &
            sizeof(mesh%nodes)/1024/1024 + &
            sizeof(mesh%tetras)/1024/1024

      print *, 'Estimated memory usage:', mem, 'Mb'

      blok = 1

      do t_cube = 1, mesh%N_cubes

         tests = mesh%tetras(:, t_cube)

         near_cubes = find_near_cubes(mesh, t_cube)

         do i1 = 1, M_cube
            if (tests(1) == 0) exit
            b_cube = near_cubes(i1)
            if (b_cube == 0) exit

            bases = mesh%tetras(:, b_cube)
            Gcorr = compute_pG(mesh, t_cube, b_cube)

            do t1 = 1, mesh%N_tet_cube
               tet_T = tests(t1)
               if (tet_T == 0) exit

               T_nodes = mesh%elem(:, tet_T)
               T_coord = mesh%node(:, T_nodes)
               vol_T = tetra_volume(T_coord)

               call gradshape(T_rot, T_dN, T_coord)

               ele_eps = mesh%eps(tet_T)

               do b1 = 1, mesh%N_tet_cube
                  tet_B = bases(b1)

                  if (tet_B == 0) exit
                  if (tet_B >= tet_T) then
                     B_nodes = mesh%elem(:, tet_B)
                     B_coord = mesh%node(:, B_nodes)
                     vol_B = tetra_volume(B_coord)

                     call gradshape(B_rot, B_dN, B_coord)

                     !_________-Correction term____________________________________

                     call precorrect_lin(matrices, tet_T, tet_B, Gcorr, corr, corr2)
                     corr = corr*sqrt(vol_T*vol_B)
                     corr2 = corr2*sqrt(vol_T*vol_B)

                     !_____________________________________________________________

                     alok(:, :) = dcmplx(0.0, 0.0)
                     if (tet_T == tet_B) then
                        call integrate_NN(alok, P0_tet, w0_tet, T_coord)
                     end if

                     call integrate_V_V_GNN(GNN, mesh, P0_tet, w0_tet, T_coord, B_coord)
                     call integrate_dV_dV_GNN(alok2, mesh, P0_tri, w0_tri, T_coord, B_coord)
                     call integrate_dV_V_GNN(alok4, mesh, P0_tet, w0_tet, P0_tri, w0_tri, T_coord, B_coord, B_dN)
                     call integrate_V_dV_GNN(alok5, mesh, P0_tet, w0_tet, P0_tri, w0_tri, T_coord, B_coord, T_dN)
                     call aloks(alok1, alok3, ind, GNN, tet_T, tet_B, T_dN, B_dN)

                     mat_block = 1.0/sqrt(vol_B*vol_T)*(alok/(ele_eps - 1.0) + &
                                                        (-mesh%k**2*(alok1 - corr) + (alok2 + alok3 - alok4 - alok5) - corr2))

                     matrices%sp_mat(:, :, blok) = cmplx(mat_block)
                     matrices%sp_ind(1, blok) = tet_T
                     matrices%sp_ind(2, blok) = tet_B

                     blok = blok + 1

                  end if
               end do
            end do
         end do
      end do

   end subroutine compute_near_zone_interactions_lin

!****************************************************************************80

   subroutine integrate_NN(alok, P0_tet, w0_tet, T_coord)
      double precision, dimension(:, :), intent(in) :: P0_tet
      double precision, dimension(:), intent(in) :: w0_tet
      double precision, intent(in) ::  T_coord(3, 4)
      double complex, dimension(12, 12) :: alok

      double complex, dimension(4, 4) :: intN
      double precision :: Pt(3, size(w0_tet)), wt(size(w0_tet))
      double precision :: shape_tet(4, size(w0_tet))
      integer :: i1

      shape_tet(1, :) = 1 - P0_tet(1, :) - P0_tet(2, :) - P0_tet(3, :)
      shape_tet(2, :) = P0_tet(1, :)
      shape_tet(3, :) = P0_tet(2, :)
      shape_tet(4, :) = P0_tet(3, :)

      call linmap_tet(Pt, wt, T_coord, P0_tet, w0_tet)
      alok(:, :) = dcmplx(0.0, 0.0)
      intN(:, :) = dcmplx(0.0, 0.0)
      do i1 = 1, size(wt)

         intN(:, 1) = intN(:, 1) + shape_tet(1, i1)*shape_tet(:, i1)*wt(i1)
         intN(:, 2) = intN(:, 2) + shape_tet(2, i1)*shape_tet(:, i1)*wt(i1)
         intN(:, 3) = intN(:, 3) + shape_tet(3, i1)*shape_tet(:, i1)*wt(i1)
         intN(:, 4) = intN(:, 4) + shape_tet(4, i1)*shape_tet(:, i1)*wt(i1)

      end do

      alok(1:4, 1:4) = intN
      alok(5:8, 5:8) = intN
      alok(9:12, 9:12) = intN

   end subroutine integrate_NN

!****************************************************************************80

   subroutine integrate_V_V_GNN(intN, mesh, P0_tet, w0_tet, T_coord, B_coord)
      type(mesh_struct), intent(in) :: mesh
      double precision, dimension(:, :), intent(in) :: P0_tet
      double precision, dimension(:), intent(in) :: w0_tet
      double precision, intent(in) ::  T_coord(3, 4), B_coord(3, 4)

      integer :: i1

      double precision :: Pb(3, size(w0_tet)), wb(size(w0_tet))
      double precision :: Pt(3, size(w0_tet)), wt(size(w0_tet))
      double precision :: shape_tet(4, size(w0_tet)), rf(3)

      double complex :: intN(4, 4), I(4)

      shape_tet(1, :) = 1 - P0_tet(1, :) - P0_tet(2, :) - P0_tet(3, :)
      shape_tet(2, :) = P0_tet(1, :)
      shape_tet(3, :) = P0_tet(2, :)
      shape_tet(4, :) = P0_tet(3, :)

      call linmap_tet(Pt, wt, T_coord, P0_tet, w0_tet)
      call linmap_tet(Pb, wb, B_coord, P0_tet, w0_tet)

      intN(:, :) = dcmplx(0, 0)

      do i1 = 1, size(wt)
         rf = Pt(:, i1)
         I = singularity_subtraction_int_V_GN(rf, mesh%k, Pb, wb, B_coord, shape_tet)

         intN(1, :) = intN(1, :) + I*shape_tet(1, i1)*wt(i1)
         intN(2, :) = intN(2, :) + I*shape_tet(2, i1)*wt(i1)
         intN(3, :) = intN(3, :) + I*shape_tet(3, i1)*wt(i1)
         intN(4, :) = intN(4, :) + I*shape_tet(4, i1)*wt(i1)

      end do

   end subroutine integrate_V_V_GNN

!****************************************************************************80

   subroutine integrate_dV_V_GNN(alok, mesh, P0_tet, w0_tet, P0_tri, w0_tri, T_coord, B_coord, B_dN)
      type(mesh_struct), intent(in) :: mesh
      double precision, dimension(:, :), intent(in) :: P0_tet, P0_tri
      double precision, dimension(:), intent(in) :: w0_tet, w0_tri
      double precision, intent(in) ::  T_coord(3, 4), B_coord(3, 4), B_dN(3, 4)
      double complex :: alok(12, 12)

      integer :: i1, i2, face_t, T_nodes(3)
      integer, parameter :: ind(12) = [2, 3, 4, 1, 4, 3, 1, 2, 4, 1, 3, 2]

      double precision :: Pb(3, size(w0_tet)), wb(size(w0_tet))
      double precision :: Pt(3, size(w0_tri)), wt(size(w0_tri))
      double precision :: shape_tet(4, size(w0_tet)), shape_tri(3, size(w0_tri)), rf(3)
      double precision :: T_tri_coord(3, 3), T_nvec(3), dd(12)

      double complex :: intN(4), I, II(12)

      shape_tet(1, :) = 1 - P0_tet(1, :) - P0_tet(2, :) - P0_tet(3, :)
      shape_tet(2, :) = P0_tet(1, :)
      shape_tet(3, :) = P0_tet(2, :)
      shape_tet(4, :) = P0_tet(3, :)

      shape_tri(1, :) = 1 - P0_tri(1, :) - P0_tri(2, :)
      shape_tri(2, :) = P0_tri(1, :)
      shape_tri(3, :) = P0_tri(2, :)

      dd = [B_dN(1, :), B_dN(2, :), B_dN(3, :)]

      call linmap_tet(Pb, wb, B_coord, P0_tet, w0_tet)

      II(:) = dcmplx(0.0, 0.0)

      do face_t = 1, 4
         T_nodes(1) = ind(3*(face_t - 1) + 1)
         T_nodes(2) = ind(3*(face_t - 1) + 2)
         T_nodes(3) = ind(3*(face_t - 1) + 3)

         T_tri_coord = T_coord(:, T_nodes)
         T_nvec = tri_n_vectors(T_tri_coord)
         call linmap_tri(Pt, wt, T_tri_coord, P0_tri, W0_tri)

         intN(:) = dcmplx(0.0, 0.0)
         do i1 = 1, size(wt)
            rf = Pt(:, i1)
            I = sum(singularity_subtraction_int_V_GN(rf, mesh%k, Pb, wb, B_coord, shape_tet))
            intN(T_nodes) = intN(T_nodes) + I*shape_tri(:, i1)*wt(i1)
         end do

         II(1:4) = II(1:4) + intN*T_nvec(1)
         II(5:8) = II(5:8) + intN*T_nvec(2)
         II(9:12) = II(9:12) + intN*T_nvec(3)

      end do

      do i2 = 1, 12
         alok(i2, :) = II(i2)*dd
      end do

   end subroutine integrate_dV_V_GNN

!****************************************************************************80

   subroutine integrate_dV_dV_GNN(alok, mesh, P0_tri, w0_tri, T_coord, B_coord)
      type(mesh_struct), intent(in) :: mesh
      double precision, dimension(:, :), intent(in) :: P0_tri
      double precision, dimension(:), intent(in) :: w0_tri
      double precision, intent(in) :: T_coord(3, 4), B_coord(3, 4)

      integer :: T_nodes(3), B_nodes(3), face_b, face_t, i1

      double precision :: Pt(3, size(w0_tri)), wt(size(w0_tri))
      double precision :: Pb(3, size(w0_tri)), wb(size(w0_tri)), shape_tri(3, size(w0_tri))
      integer, parameter :: ind(12) = [2, 3, 4, 1, 4, 3, 1, 2, 4, 1, 3, 2]
      double precision :: T_tri_coord(3, 3), B_tri_coord(3, 3)
      double precision :: T_nvec(3), B_nvec(3), rf(3)
      double complex :: Ix(4), Iy(4), Iz(4), IIx(4, 4), IIy(4, 4), IIz(4, 4), alok(12, 12)
      double complex :: intN(3)

      shape_tri(1, :) = 1.0 - P0_tri(1, :) - P0_tri(2, :)
      shape_tri(2, :) = P0_tri(1, :)
      shape_tri(3, :) = P0_tri(2, :)

      alok(:, :) = dcmplx(0.0, 0.0)

      do face_t = 1, 4
         T_nodes(1) = ind(3*(face_t - 1) + 1)
         T_nodes(2) = ind(3*(face_t - 1) + 2)
         T_nodes(3) = ind(3*(face_t - 1) + 3)

         !T_tri_coord = tetra_face(T_coord,face_t)
         T_tri_coord = T_coord(:, T_nodes)
         T_nvec = tri_n_vectors(T_tri_coord)
         call linmap_tri(Pt, wt, T_tri_coord, P0_tri, W0_tri)

         IIx(:, :) = dcmplx(0.0, 0.0)
         IIy(:, :) = dcmplx(0.0, 0.0)
         IIz(:, :) = dcmplx(0.0, 0.0)

         do i1 = 1, size(wt)
            rf = Pt(:, i1)

            !________________________________________________________
            Ix(:) = dcmplx(0.0, 0.0)
            Iy(:) = dcmplx(0.0, 0.0)
            Iz(:) = dcmplx(0.0, 0.0)

            do face_b = 1, 4

               B_nodes(1) = ind(3*(face_b - 1) + 1)
               B_nodes(2) = ind(3*(face_b - 1) + 2)
               B_nodes(3) = ind(3*(face_b - 1) + 3)

               !B_tri_coord = tetra_face(B_coord,face_b)
               B_tri_coord = B_coord(:, B_nodes)
               B_nvec = tri_n_vectors(B_tri_coord)
               call linmap_tri(Pb, wb, B_tri_coord, P0_tri, W0_tri)

               intN = singularity_subtraction_int_S_GN(rf, B_tri_coord, B_nvec, mesh%k, Pb, wb, shape_tri)

               Ix(B_nodes) = Ix(B_nodes) + intN*B_nvec(1)
               Iy(B_nodes) = Iy(B_nodes) + intN*B_nvec(2)
               Iz(B_nodes) = Iz(B_nodes) + intN*B_nvec(3)

            end do
            !_________________________________________________________

            IIx(T_nodes(1), :) = IIx(T_nodes(1), :) + Ix*shape_tri(1, i1)*wt(i1)
            IIx(T_nodes(2), :) = IIx(T_nodes(2), :) + Ix*shape_tri(2, i1)*wt(i1)
            IIx(T_nodes(3), :) = IIx(T_nodes(3), :) + Ix*shape_tri(3, i1)*wt(i1)

            IIy(T_nodes(1), :) = IIy(T_nodes(1), :) + Iy*shape_tri(1, i1)*wt(i1)
            IIy(T_nodes(2), :) = IIy(T_nodes(2), :) + Iy*shape_tri(2, i1)*wt(i1)
            IIy(T_nodes(3), :) = IIy(T_nodes(3), :) + Iy*shape_tri(3, i1)*wt(i1)

            IIz(T_nodes(1), :) = IIz(T_nodes(1), :) + Iz*shape_tri(1, i1)*wt(i1)
            IIz(T_nodes(2), :) = IIz(T_nodes(2), :) + Iz*shape_tri(2, i1)*wt(i1)
            IIz(T_nodes(3), :) = IIz(T_nodes(3), :) + Iz*shape_tri(3, i1)*wt(i1)

         end do

         alok(1:4, 1:4) = alok(1:4, 1:4) + IIx*T_nvec(1)
         alok(1:4, 5:8) = alok(1:4, 5:8) + IIy*T_nvec(1)
         alok(1:4, 9:12) = alok(1:4, 9:12) + IIz*T_nvec(1)

         alok(5:8, 1:4) = alok(5:8, 1:4) + IIx*T_nvec(2)
         alok(5:8, 5:8) = alok(5:8, 5:8) + IIy*T_nvec(2)
         alok(5:8, 9:12) = alok(5:8, 9:12) + IIz*T_nvec(2)

         alok(9:12, 1:4) = alok(9:12, 1:4) + IIx*T_nvec(3)
         alok(9:12, 5:8) = alok(9:12, 5:8) + IIy*T_nvec(3)
         alok(9:12, 9:12) = alok(9:12, 9:12) + IIz*T_nvec(3)

      end do

   end subroutine integrate_dV_dV_GNN

!**********************************************************************

   subroutine integrate_V_dV_GNN(alok, mesh, P0_tet, w0_tet, P0_tri, w0_tri, T_coord, B_coord, T_dN)
      type(mesh_struct), intent(in) :: mesh
      double precision, dimension(:, :), intent(in) :: P0_tri, P0_tet
      double precision, dimension(:), intent(in) :: w0_tri, w0_tet
      double precision, intent(in) :: T_coord(3, 4), B_coord(3, 4), T_dN(3, 4)

      integer :: B_nodes(3), face_b, i1, i2

      double precision :: Pt(3, size(w0_tet)), wt(size(w0_tet)), shape_tet(4, size(w0_tet))
      double precision :: Pb(3, size(w0_tri)), wb(size(w0_tri)), shape_tri(3, size(w0_tri))
      integer, parameter :: ind(12) = [2, 3, 4, 1, 4, 3, 1, 2, 4, 1, 3, 2]
      double precision :: B_tri_coord(3, 3)
      double precision :: B_nvec(3), rf(3), dd(12)
      double complex :: I(12), II(12), alok(12, 12)
      double complex :: intN(3)

      shape_tri(1, :) = 1.0 - P0_tri(1, :) - P0_tri(2, :)
      shape_tri(2, :) = P0_tri(1, :)
      shape_tri(3, :) = P0_tri(2, :)

      shape_tet(1, :) = 1 - P0_tet(1, :) - P0_tet(2, :) - P0_tet(3, :)
      shape_tet(2, :) = P0_tet(1, :)
      shape_tet(3, :) = P0_tet(2, :)
      shape_tet(4, :) = P0_tet(3, :)

      dd = [T_dN(1, :), T_dN(2, :), T_dN(3, :)]

      alok(:, :) = dcmplx(0.0, 0.0)
      II(:) = dcmplx(0.0, 0.0)

      call linmap_tet(Pt, wt, T_coord, P0_tet, w0_tet)

      do i1 = 1, size(wt)

         rf = Pt(:, i1)

         I(:) = dcmplx(0.0, 0.0)

         do face_b = 1, 4

            B_nodes(1) = ind(3*(face_b - 1) + 1)
            B_nodes(2) = ind(3*(face_b - 1) + 2)
            B_nodes(3) = ind(3*(face_b - 1) + 3)

            !B_tri_coord = tetra_face(B_coord,face_b)
            B_tri_coord = B_coord(:, B_nodes)
            B_nvec = tri_n_vectors(B_tri_coord)
            call linmap_tri(Pb, wb, B_tri_coord, P0_tri, W0_tri)

            intN = singularity_subtraction_int_S_GN(rf, B_tri_coord, B_nvec, mesh%k, Pb, wb, shape_tri)

            I(B_nodes) = I(B_nodes) + intN*B_nvec(1)
            I(B_nodes + 4) = I(B_nodes + 4) + intN*B_nvec(2)
            I(B_nodes + 8) = I(B_nodes + 8) + intN*B_nvec(3)

         end do

         II = II + I*wt(i1)

      end do

      do i2 = 1, 12
         alok(i2, :) = dd(i2)*II
      end do

   end subroutine integrate_V_dV_GNN

!**************************************************************************
! Local matrices
   subroutine aloks(alok, alok2, ind, GNN, test_T, basis_T, T_dN, B_dN)
      double complex :: alok(12, 12), alok2(12, 12)
      integer :: ind(12, 2), i2
      double complex :: GNN(4, 4)
      double precision :: T_dN(3, 4), B_dN(3, 4), ddT(12), ddB(12)
      integer :: test_T, basis_T, i1, nt, nb, ct, cb

      alok = dcmplx(0.0, 0.0)
      alok2 = dcmplx(0.0, 0.0)

      alok(1:4, 1:4) = GNN
      alok(5:8, 5:8) = GNN
      alok(9:12, 9:12) = GNN

      do i1 = 1, 12
         ind(i1, 1) = 12*(test_T - 1) + i1
         ind(i1, 2) = 12*(basis_T - 1) + i1
      end do

      ddT = [T_dN(1, :), T_dN(2, :), T_dN(3, :)]
      ddB = [B_dN(1, :), B_dN(2, :), B_dN(3, :)]

      do i2 = 1, 12
         alok2(i2, :) = ddT(i2)*ddB
      end do

      alok2 = alok2*sum(GNN)

   end subroutine

!****************************************************************************80

   subroutine integrate_V_V_G(intN2, mesh, P0_tet, w0_tet, P02_tet, w02_tet, T_coord, B_coord)
      type(mesh_struct), intent(in) :: mesh
      double precision, dimension(:, :), intent(in) :: P0_tet, P02_tet
      double precision, dimension(:), intent(in) :: w0_tet, w02_tet
      double precision, intent(in) ::  T_coord(3, 4), B_coord(3, 4)

      integer :: i1
      double precision :: Pt(3, size(w02_tet)), wt(size(w02_tet))
      double precision :: Pb(3, size(w0_tet)), wb(size(w0_tet))

      double precision :: rf(3), bt(3, 3)
      double complex :: intN, intN2(3, 3)

      bt(:, :) = 0
      bt(1, 1) = 1
      bt(2, 2) = 1
      bt(3, 3) = 1

      call linmap_tet(Pt, wt, T_coord, P02_tet, w02_tet)
      call linmap_tet(Pb, wb, B_coord, P0_tet, w0_tet)

      intN2(:, :) = dcmplx(0, 0)

      do i1 = 1, size(wt)

         rf = Pt(:, i1)
         intN = singularity_subtraction_int_V_G(rf, mesh%k, Pb, wb, B_coord)
         intN2 = intN2 + bt*intN*wt(i1)

      end do

   end subroutine integrate_V_V_G

!****************************************************************************80

   subroutine integrate_dV_dV_G(intN2, mesh, P0_tri, w0_tri, P02_tri, w02_tri, T_coord, B_coord)
      type(mesh_struct), intent(in) :: mesh
      double precision, dimension(:, :), intent(in) :: P0_tri, P02_tri
      double precision, dimension(:), intent(in) :: w0_tri, w02_tri
      double precision, intent(in) :: T_coord(3, 4), B_coord(3, 4)

      integer :: face_t, face_b, i1
      double precision :: Pt(3, size(w02_tri)), wt(size(w02_tri))
      double precision :: Pb(3, size(w0_tri)), wb(size(w0_tri))

      double precision :: T_tri_coord(3, 3), B_tri_coord(3, 3), bt(3, 3)
      double precision :: rf(3), T_nvec(3), B_nvec(3)
      double complex :: intN, intN2(3, 3)

      intN2(:, :) = dcmplx(0, 0)

      do face_t = 1, 4

         T_tri_coord = tetra_face(T_coord, face_t)
         T_nvec = tri_n_vectors(T_tri_coord)
         call linmap_tri(Pt, wt, T_tri_coord, P02_tri, W02_tri)

         do face_b = 1, 4
            B_tri_coord = tetra_face(B_coord, face_b)
            B_nvec = tri_n_vectors(B_tri_coord)
            call linmap_tri(Pb, wb, B_tri_coord, P0_tri, W0_tri)

            bt(:, 1) = T_nvec*B_nvec(1)
            bt(:, 2) = T_nvec*B_nvec(2)
            bt(:, 3) = T_nvec*B_nvec(3)

            do i1 = 1, size(wt)

               rf = Pt(:, i1)
               intN = singularity_subtraction_int_S_G(rf, B_tri_coord, B_nvec, mesh%k, Pb, wb)
               intN2 = intN2 + bt*intN*wt(i1)
            end do

         end do

      end do

   end subroutine integrate_dV_dV_G

end module precorrection
