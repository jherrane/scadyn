module geometry
! Copyright (c) 2018 Johannes Markkanen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common

   implicit none

contains

!****************************************************************************80
! Finds the index of the cube in which point is
   integer function point_in_cube(mesh, point)

      real(dp), dimension(3), intent(in) :: point
      type(mesh_struct) :: mesh

      integer :: Mx, My, Mz, Nx_box, Ny_box, Nz_box
      real(dp), dimension(3, (mesh%M_ex + 1)**3) :: nn
      real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax

      Mx = floor((point(1) - mesh%min_coord(1))/mesh%box_delta); 
      My = floor((point(2) - mesh%min_coord(2))/mesh%box_delta); 
      Mz = floor((point(3) - mesh%min_coord(3))/mesh%box_delta); 
      Nx_box = mesh%Nx_cube
      Ny_box = mesh%Ny_cube

      point_in_cube = Mx + Nx_box*My + &
                      Nx_box*Ny_box*Mz + 1; 
!  Check
      if (point_in_cube < 1) then
         print *, 'something is wrong (Point is not in any cube)'
      end if

      nn = mesh%nodes(:, mesh%etopol_box(:, point_in_cube))

      xmin = minval(nn(1, :)); 
      ymin = minval(nn(2, :)); 
      zmin = minval(nn(3, :)); 
      xmax = maxval(nn(1, :)); 
      ymax = maxval(nn(2, :)); 
      zmax = maxval(nn(3, :)); 
      if (point(1) < xmin .or. point(1) >= xmax) then
         print *, 'something is wrong x '
      end if

      if (point(2) < ymin .or. point(2) >= ymax) then
         print *, 'something is wrong y'
      end if

      if (point(3) < zmin .or. point(3) >= zmax) then
         print *, 'something is wrong z'
      end if

   end function point_in_cube

!****************************************************************************80
! Routine creates 3D grid for fft
! Output: nodes (3,N_nodes) Coordinates of grid points
   subroutine build_grid(mesh)
      type(mesh_struct) :: mesh

      real(dp), dimension(:, :), allocatable :: nodes
      real(dp), dimension(3) :: min_coord, max_coord
      real(dp) :: delta, min_x, min_y, min_z, max_x, max_y, max_z, PP(3, 4), vol
      integer :: las, m1, l1, k1, P(4), i1

      min_coord = 1.0*minval(mesh%coord, dim=2)
      max_coord = 1.0*maxval(mesh%coord, dim=2)

      mesh%min_coord = min_coord

      min_x = min_coord(1)
      min_y = min_coord(2)
      min_z = min_coord(3)

      max_x = max_coord(1)
      max_y = max_coord(2)
      max_z = max_coord(3)

      vol = 0.0
      do i1 = 1, mesh%N_tet
         P = mesh%etopol(:, i1)
         PP = mesh%coord(:, P)
         vol = vol + tetra_volume(PP)
      end do

      delta = mesh%grid_size*1.0/mesh%M_ex*(vol/dble(mesh%N_tet))**(1.0/3.0)

      mesh%Nx = int(ceiling(abs(max_coord(1) - min_coord(1))/delta)) + 1
      mesh%Ny = int(ceiling(abs(max_coord(2) - min_coord(2))/delta)) + 1
      mesh%Nz = int(ceiling(abs(max_coord(3) - min_coord(3))/delta)) + 1

      if (mesh%M_ex == 1 .or. mesh%M_ex == 2) then
         mesh%Nx = mesh%Nx + mod(mesh%Nx + 1, 2) + 1
         mesh%Ny = mesh%Ny + mod(mesh%Ny + 1, 2) + 1
         mesh%Nz = mesh%Nz + mod(mesh%Nz + 1, 2) + 1
      end if

      if (mesh%M_ex == 3) then
         mesh%Nx = mesh%Nx - mod(mesh%Nx - 1, 3) + 3
         mesh%Ny = mesh%Ny - mod(mesh%Ny - 1, 3) + 3
         mesh%Nz = mesh%Nz - mod(mesh%Nz - 1, 3) + 3
      end if

      allocate (nodes(3, mesh%Nx*mesh%Ny*mesh%Nz))

      las = 1; 
      do m1 = 1, mesh%Nz
         do l1 = 1, mesh%Ny
            do k1 = 1, mesh%Nx
               nodes(1, las) = min_x + delta*(k1 - 1); 
               nodes(2, las) = min_y + delta*(l1 - 1); 
               nodes(3, las) = min_z + delta*(m1 - 1); 
               las = las + 1; 
            end do
         end do
      end do

      allocate (mesh%nodes(3, mesh%Nx*mesh%Ny*mesh%Nz))

      mesh%delta = delta
      mesh%nodes = nodes

   end subroutine build_grid

!****************************************************************************80
! Routine creates sub-cubes
! Output: etopol_box ((M+1)^3,N_box) indeces of box's nodes
   subroutine build_box(mesh)

      type(mesh_struct) :: mesh

      integer, dimension(:, :), allocatable :: etopol_box
      integer, dimension(:, :), allocatable :: etopol_box2
      integer, dimension(:, :), allocatable :: etopol_box3
      integer :: las, x, y, z, Nx, Ny, Nz, M
      integer, dimension(8) :: co, loc
      integer, dimension(9) :: co2
      integer, dimension(27) :: co3, loc2
      integer :: M3(16), coM3(64), loc3(64)

      allocate (etopol_box(8, (mesh%Nx - 1)*(mesh%Ny - 1)*(mesh%Nz - 1)))
      allocate (etopol_box2(27, (mesh%Nx - 2)/2*(mesh%Ny - 2)/2*(mesh%Nz - 2)/2))
      allocate (etopol_box3(64, (mesh%Nx - 1)/3*(mesh%Ny - 1)/3*(mesh%Nz - 1)/3))

      Nx = mesh%Nx
      Ny = mesh%Ny
      Nz = mesh%Nz

      M = mesh%M_ex

      if (M == 1) then
         co = (/1, 2, Nx + 2, Nx + 1, Nx*Ny + 1, Nx*Ny + 2, Nx*Ny + Nx + 2, Nx*Ny + Nx + 1/)

         las = 0
         loc = co

         do z = 1, Nz - 1
            do y = 1, Ny - 1
               do x = 1, Nx - 1
                  etopol_box(:, las + 1) = loc
                  las = las + 1

                  if (x == (Nx - 1)) then
                     loc = loc + 1
                  end if

                  if (x == (Nx - 1) .and. y == (Ny - 1)) then
                     loc = loc + Nx
                  end if
                  loc = loc + 1
               end do
            end do
         end do

         allocate (mesh%etopol_box(8, (mesh%Nx - 1)*(mesh%Ny - 1)*(mesh%Nz - 1)))
         mesh%etopol_box = etopol_box
         mesh%N_cubes = (Nx - 1)*(Ny - 1)*(Nz - 1)
         mesh%Nx_cube = (Nx - 1)
         mesh%Ny_cube = (Ny - 1)
         mesh%Nz_cube = (Nz - 1)
      end if

      if (M == 2) then

         co2 = (/1, 2, 3, Nx + 1, Nx + 2, Nx + 3, 2*Nx + 1, 2*Nx + 2, 2*Nx + 3/)

         co3 = [co2, Nx*Ny + co2, 2*Nx*Ny + co2]

         las = 0; 
         loc2 = co3; 
         do z = 1, (Nz - 2)/M
            do y = 1, (Ny - 2)/M
               do x = 1, (Nx - 2)/M
                  etopol_box2(:, las + 1) = loc2
                  las = las + 1

                  if (x == (Nx - 2)/M) then
                     loc2 = loc2 + 2 + Nx
                  end if

                  if (x == (Nx - 2)/M .and. y == (Ny - 2)/M) then
                     loc2 = loc2 + Nx*Ny + Nx + Nx
                  end if
                  loc2 = loc2 + 2
               end do
            end do
         end do

         allocate (mesh%etopol_box(27, (mesh%Nx - 2)/2*(mesh%Ny - 2)/2*(mesh%Nz - 2)/2))
         mesh%etopol_box = etopol_box2
         mesh%Nx_cube = (Nx - 2)/2
         mesh%Ny_cube = (Ny - 2)/2
         mesh%Nz_cube = (Nz - 2)/2

         mesh%N_cubes = (Nx - 2)/2*(Ny - 2)/2*(Nz - 2)/2

      end if

      if (M == 3) then

         M3 = (/1, 2, 3, 4, Nx + 1, Nx + 2, Nx + 3, Nx + 4, 2*Nx + 1, 2*Nx + 2, 2*Nx + 3, 2*Nx + 4, &
                3*Nx + 1, 3*Nx + 2, 3*Nx + 3, 3*Nx + 4/)

         coM3 = [M3, Nx*Ny + M3, 2*Nx*Ny + M3, 3*Nx*Ny + M3]

         las = 0; 
         loc3 = coM3; 
         if (mod((Nx - 1), 3) .ne. 0) then
            print *, 'ERROR grid size'
            stop
         end if
         if (mod((Ny - 1), 3) .ne. 0) then
            print *, 'ERROR grid size'
            stop
         end if
         if (mod((Nz - 1), 3) .ne. 0) then
            print *, 'ERROR grid size'
            stop
         end if

         do z = 1, (Nz - 1)/M
            do y = 1, (Ny - 1)/M
               do x = 1, (Nx - 1)/M
                  etopol_box3(:, las + 1) = loc3
                  las = las + 1

                  if (x == (Nx - 1)/M) then
                     loc3 = loc3 + 1 + 2*Nx
                  end if

                  if (x == (Nx - 1)/M .and. y == (Ny - 1)/M) then
                     loc3 = loc3 + 2*Nx*Ny + Nx
                  end if
                  loc3 = loc3 + 3
               end do
            end do
         end do

         allocate (mesh%etopol_box(64, (mesh%Nx - 1)/3*(mesh%Ny - 1)/3*(mesh%Nz - 1)/3))
         mesh%etopol_box = etopol_box3
         mesh%N_cubes = (Nx - 1)/3*(Ny - 1)/3*(Nz - 1)/3
         mesh%Nx_cube = (Nx - 1)/3
         mesh%Ny_cube = (Ny - 1)/3
         mesh%Nz_cube = (Nz - 1)/3
      end if

      mesh%box_delta = mesh%delta*M
   end subroutine build_box

!****************************************************************************80
! Routine finds near-zone cubes
! Output: mesh.near_cubes(27,N_cubes)
   function find_near_cubes(mesh, cube) result(N)

      type(mesh_struct) :: mesh
      integer, intent(in) :: cube
      integer, dimension(:), allocatable ::  N

      integer :: Mx, My, Mz, max_cube, i1, i2, i3, las, x, y, z, xx, yy, zz, M, pM

      M = (2*mesh%near_zone + 1)
      pM = mesh%near_zone + 1
      allocate (N(M**3))
      N(:) = 0

      Mx = mesh%Nx_cube
      My = mesh%Ny_cube
      Mz = mesh%Nz_cube

      max_cube = Mx*My*Mz

      z = (cube - 1)/(Mx*My) + 1
!check
      if (z > Mz) then
         print *, 'Error in function find_near_cubes'
      end if

      y = (cube - (z - 1)*Mx*My - 1)/Mx + 1
      x = (cube - (z - 1)*Mx*My - (y - 1)*Mx - 1) + 1

      if (x + Mx*(y - 1) + Mx*My*(z - 1) .ne. cube) then
         print *, 'Error in find_near_cubes'
      end if

      las = 1
      do i1 = 1, M
         zz = z - pM + i1
         if (zz > 0 .and. zz < Mz + 1) then
            do i2 = 1, M
               yy = y - pM + i2
               if (yy > 0 .and. yy < My + 1) then
                  do i3 = 1, M
                     xx = x - pM + i3
                     if (xx > 0 .and. xx < Mx + 1) then
                        N(las) = xx + Mx*(yy - 1) + Mx*My*(zz - 1)
                        las = las + 1
                     end if
                  end do
               end if
            end do
         end if
      end do

   end function find_near_cubes

!****************************************************************************80

   subroutine tetras_in_cubes(mesh)

      type(mesh_struct) :: mesh

      integer :: tetra, cube, i1, i2, max
      integer, dimension(4) :: p
      real(dp), dimension(3) :: cp
      real(dp), dimension(3, 4) :: co
      integer, dimension(mesh%N_tet) :: tmp
      integer, dimension(mesh%N_cubes) :: tmp2
      integer, dimension(:, :), allocatable :: tetras

      do tetra = 1, mesh%N_tet

         p = mesh%etopol(:, tetra)
         co = mesh%coord(:, p)

         cp = (co(:, 1) + co(:, 2) + co(:, 3) + co(:, 4))/4 !Center of tetra

         cube = point_in_cube(mesh, cp)

         tmp(tetra) = cube

      end do

      tmp2(:) = 0
      do i1 = 1, mesh%N_tet
         tmp2(tmp(i1)) = tmp2(tmp(i1)) + 1
      end do

      max = maxval(tmp2)

      allocate (tetras(max, mesh%N_cubes))
      tetras(:, :) = 0

      tmp2(:) = 0

      do i2 = 1, mesh%N_tet
         tmp2(tmp(i2)) = tmp2(tmp(i2)) + 1
         tetras(tmp2(tmp(i2)), tmp(i2)) = i2

      end do

      allocate (mesh%tetras(max, mesh%N_cubes))
      mesh%tetras = tetras
      mesh%N_tet_cube = max

   end subroutine tetras_in_cubes

end module geometry
