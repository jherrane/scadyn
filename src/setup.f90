module setup
   use common
   use io
   use mie
   use translations
   use projection
   use octtree

   implicit none
contains

!****************************************************************************80

   subroutine init_geometry()
      type(data_struct), dimension(:), allocatable :: sphere
      type(level_struct), dimension(:), allocatable :: otree
      real(dp) :: ka, maxrad, maxrad_sph, vol
      integer :: i, Nspheres, sph, max_level, max_sph
      complex(dp) :: k

!* TETRAHEDRAL MESH *************************************************
      if (matrices%is_aggr .NE. 1 .AND. use_mie .NE. 1) then
         if (matrices%Tmat == 0) print *, 'Order of basis functions  =', mesh%order
         if (mesh%M_ex > 3 .or. mesh%M_ex < 1) then
            print *, 'ERROR: order should be 1, 2 or 3'
            stop
         end if
         if (mesh%order > 1 .or. mesh%M_ex < 0) then
            print *, 'ERROR: Expansion order should be 0 or 1 '
            stop
         end if

         print *, 'Reading mesh...'
         call read_mesh() ! io
         call int_points(mesh)

         if (matrices%Tmat == 0) print *, 'Initializing FFT... '
         if (allocated(mesh%nodes)) deallocate (mesh%nodes, mesh%etopol_box, mesh%tetras)
         call build_grid2(mesh) ! geometry
         call build_box(mesh) ! geometry
         call tetras_in_cubes(mesh) ! geometry

         if (matrices%Tmat == 0) then
            print *, '   Grid size            = ', (/mesh%Nx, mesh%Ny, mesh%Nz/)
            print *, '   Delta grid           = ', real(mesh%delta)
            print *, '   Number of cubes      = ', mesh%N_cubes
            print *, '   Delta cubes          = ', real(mesh%box_delta)
            print *, '   Elems. in cube (max) = ', mesh%N_tet_cube
            print *, 'Done'
         end if

         if (matrices%Tmat == 0) call construct_projectors() ! projection

! Check whether using maxval is good or not
         do i = 1, matrices%bars
            ka = mesh%ki(i)*(dble(maxval([mesh%Nx, mesh%Ny, mesh%Nz])) &
                             *mesh%delta)/2.0d0
            matrices%Nmaxs(i) = truncation_order(ka)
         end do

      else if (use_mie .NE. 1) then !* SPHERICAL AGGREGATE *******************
         call read_aggr()
         Nspheres = size(mesh%radius)
         allocate (sphere(Nspheres))

         do i = 1, matrices%bars
            k = dcmplx(mesh%ki(i), 0d0)
            maxrad = 0d0
            vol = 0d0
            do sph = 1, Nspheres
               sphere(sph)%cp = mesh%coord(:, sph)
               sphere(sph)%r = mesh%radius(sph)

               ka = dble(mesh%ki(i))*mesh%radius(sph)
               sphere(sph)%Nmax = truncation_order2(ka) ! Note

               sphere(sph)%Tmat_ind = 1

               maxrad_sph = sqrt(dot_product(mesh%coord(:, sph), mesh%coord(:, sph))) + mesh%radius(sph)
               if (maxrad < maxrad_sph) then
                  maxrad = maxrad_sph
               end if

               if (max_sph < sphere(sph)%Nmax) then
                  max_sph = sphere(sph)%Nmax
               end if

               vol = vol + 4.0/3.0*pi*(mesh%radius(sph))**3
            end do
            if (allocated(otree)) deallocate (otree)
            call create_octtree(sphere, otree, dble(mesh%ki(i)), max_level)
         end do
         mesh%maxrad = maxrad
         call int_points(mesh)

         do i = 1, matrices%bars
            ka = mesh%ki(i)*sqrt(3.0)/2.0*otree(1)%tree(1)%dl
            matrices%Nmaxs(i) = truncation_order2(ka)
         end do
      else ! MIE *************************************************************
         call int_points_mie(mesh)
         do i = 1, matrices%bars
            ka = mesh%ki(i)*mesh%a
            matrices%Nmaxs(i) = truncation_order(ka)
         end do

         allocate (mesh%param(1))
         mesh%param = dcmplx(matrices%refr**2 - matrices%refi**2, &
                             2d0*matrices%refr*matrices%refi)
         write (*, '(2(A,F5.3))') '    Refractive index     =   ', matrices%refr, ' + i', matrices%refi
         write (*, '(2(A,F5.3))') '    Dielectric constant  =   ', matrices%refr**2 - matrices%refi**2, &
            ' + i', 2d0*matrices%refr*matrices%refi
      end if
   end subroutine init_geometry

!****************************************************************************80

   subroutine construct_projectors()
      integer :: i, M_box
      real(dp) :: k_orig

      k_orig = mesh%k ! Keep the original just in case
      M_box = size(mesh%etopol_box, 1)

      if (.not. allocated(matrices%rhs)) then
         if (mesh%order == 0) then
            allocate (matrices%rhs(3*mesh%N_tet))
            allocate (matrices%x(3*mesh%N_tet))
            allocate (matrices%Ax(3*mesh%N_tet))

            allocate (matrices%listS(M_box, mesh%N_tet, matrices%bars))
            allocate (matrices%listSx(M_box, mesh%N_tet, matrices%bars))
            allocate (matrices%listSy(M_box, mesh%N_tet, matrices%bars))
            allocate (matrices%listSz(M_box, mesh%N_tet, matrices%bars))
            allocate (matrices%listindS(M_box, mesh%N_tet, matrices%bars))
         end if
         if (mesh%order == 1) then
            allocate (matrices%rhs(4*3*mesh%N_tet))
            allocate (matrices%x(4*3*mesh%N_tet))
            allocate (matrices%Ax(4*3*mesh%N_tet))

            allocate (matrices%listS(M_box, 4*mesh%N_tet, matrices%bars))
            allocate (matrices%listSx(M_box, 4*mesh%N_tet, matrices%bars))
            allocate (matrices%listSy(M_box, 4*mesh%N_tet, matrices%bars))
            allocate (matrices%listSz(M_box, 4*mesh%N_tet, matrices%bars))
            allocate (matrices%listindS(M_box, 4*mesh%N_tet, matrices%bars))
         end if
      end if
      print *, 'Constructing ', trim(mesh%projector), '-projectors... '

      do i = 1, size(mesh%ki)
         if (allocated(matrices%S)) then
            deallocate (matrices%S, matrices%Sx, matrices%Sy, matrices%Sz, matrices%indS)
         end if
         mesh%k = mesh%ki(i)
         call print_bar(i, size(mesh%ki))
         if (mesh%order == 1) call pfft_projection_lin(matrices, mesh)
         if (mesh%order == 0) call pfft_projection_const(matrices, mesh)
         matrices%listS(:, :, i) = matrices%S
         matrices%listSx(:, :, i) = matrices%Sx
         matrices%listSy(:, :, i) = matrices%Sy
         matrices%listSz(:, :, i) = matrices%Sz
         matrices%listindS(:, :, i) = matrices%indS
      end do

      mesh%k = k_orig

   end subroutine construct_projectors

!****************************************************************************80

   subroutine update_projections(i)
      integer :: i

      matrices%S = matrices%listS(:, :, i)
      matrices%Sx = matrices%listSx(:, :, i)
      matrices%Sy = matrices%listSy(:, :, i)
      matrices%Sz = matrices%listSz(:, :, i)
      matrices%indS = matrices%listindS(:, :, i)

   end subroutine update_projections

!****************************************************************************80

   subroutine allocate_inc_wave()
      complex(dp), dimension(:), allocatable :: a_in, b_in
      integer :: Nmax, i, las, nm

      Nmax = maxval(matrices%Nmaxs)

      allocate (matrices%as((Nmax + 1)**2 - 1, matrices%bars))
      allocate (matrices%bs((Nmax + 1)**2 - 1, matrices%bars))

      las = (Nmax + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1

      allocate (matrices%rotDs(las, matrices%bars))
      allocate (matrices%indDs(las, 2, matrices%bars))
      allocate (matrices%rotD90s(las, matrices%bars))
      allocate (matrices%indD90s(las, 2, matrices%bars))
      allocate (matrices%rotXs(las, matrices%bars))
      allocate (matrices%indXs(las, 2, matrices%bars))
      allocate (matrices%rotYs(las, matrices%bars))
      allocate (matrices%indYs(las, 2, matrices%bars))

      do i = 1, matrices%bars
         if (allocated(a_in)) deallocate (a_in, b_in)
         Nmax = matrices%Nmaxs(i)
         nm = (Nmax + 1)**2 - 1
         allocate (a_in(nm))
         allocate (b_in(nm))
         call planewave(Nmax, dcmplx(mesh%ki(i)), a_in, b_in)

         matrices%as(1:nm, i) = a_in
         matrices%bs(1:nm, i) = b_in
      end do

   end subroutine allocate_inc_wave

!****************************************************************************80

   subroutine update_inc_wave(i, a_nm, b_nm, a_nm90, b_nm90)
      complex(dp), dimension(:), allocatable :: a_in, b_in, a90, b90, &
                                                a, b, a_nm, b_nm, a_nm90, b_nm90
      complex(dp), dimension(:, :), allocatable :: Taa, Tab, Tba, Tbb
      complex(8), dimension(:), allocatable :: rotD, rotD90
      integer, dimension(:, :), allocatable :: indD, indD90
      integer :: Nmax, i, las, nm

      matrices%E0 = matrices%E_rel(i)*matrices%E*matrices%E0hat
      matrices%E90 = matrices%E_rel(i)*matrices%E*matrices%E90hat
      mesh%k = mesh%ki(i)

      Nmax = matrices%Nmaxs(i)
      las = (Nmax + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1
      nm = (Nmax + 1)**2 - 1
      rotD = matrices%rotDs(1:las, i)
      indD = matrices%indDs(1:las, :, i)
      rotD90 = matrices%rotD90s(1:las, i)
      indD90 = matrices%indD90s(1:las, :, i)

      Taa = matrices%Taai(1:nm, 1:nm, i)
      Tab = matrices%Tabi(1:nm, 1:nm, i)
      Tba = matrices%Tbai(1:nm, 1:nm, i)
      Tbb = matrices%Tbbi(1:nm, 1:nm, i)

      a_in = matrices%E_rel(i)*matrices%E*matrices%as(1:nm, i)
      b_in = matrices%E_rel(i)*matrices%E*matrices%bs(1:nm, i)

      a = sparse_matmul(rotD, indD, a_in, nm)
      b = sparse_matmul(rotD, indD, b_in, nm)
      a90 = sparse_matmul(rotD90, indD90, a_in, nm)
      b90 = sparse_matmul(rotD90, indD90, b_in, nm)

      a_nm = matmul(Taa, a) + matmul(Tab, b)
      b_nm = matmul(Tbb, b) + matmul(Tba, a)
      a_nm90 = matmul(Taa, a90) + matmul(Tab, b90)
      b_nm90 = matmul(Tbb, b90) + matmul(Tba, a90)

   end subroutine update_inc_wave

!****************************************************************************80

   subroutine rot_setup()
      real(dp), dimension(3, 3) ::  RT, RT90

      RT = transpose(matrices%R)
      RT90 = matmul(RT, matrices%R90_init)

      matrices%Rk = matrices%R
      matrices%khat = matmul(RT, [0d0, 0d0, 1d0])
      matrices%E0hat = dcmplx(matmul(RT, [1d0, 0d0, 0d0]))
      matrices%E90hat = dcmplx(matmul(RT, [0d0, 1d0, 0d0]))

      call gen_rotations(RT, RT90)

   end subroutine rot_setup

!****************************************************************************80

   subroutine fix_band()
      if (dabs(matrices%lambda1 - 2d0*pi/mesh%ki(1)) > 1d-7) then
         matrices%lambda1 = 2d0*pi/mesh%ki(1)
         matrices%lambda2 = 2d0*pi/mesh%ki(matrices%bars)
         if (matrices%waves == 'bnd') call calc_E_rel()
      end if

   end subroutine fix_band

!****************************************************************************80

   subroutine gen_rotations(R, R90)
      integer :: las, Nmax, i
      complex(8), dimension(:), allocatable :: rotD, rotD90, rotX, rotY
      integer, dimension(:, :), allocatable :: indD, indD90, indX, indY
      real(dp), dimension(3, 3) :: R, R90, Rx, Ry

      Rx = reshape([0d0, 0d0, 1d0, 0d0, 1d0, 0d0, -1d0, 0d0, 0d0], [3, 3])
      Ry = reshape([1d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0, -1d0, 0d0], [3, 3])

      do i = 1, matrices%bars
         Nmax = matrices%Nmaxs(i)
         las = (Nmax + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1

         if (allocated(rotD)) deallocate (rotD, indD, rotD90, indD90, rotX, indX, rotY, indY)
         allocate (rotD(las))
         allocate (indD(las, 2))
         allocate (rotD90(las))
         allocate (indD90(las, 2))
         allocate (rotX(las))
         allocate (indX(las, 2))
         allocate (rotY(las))
         allocate (indY(las, 2))

         call sph_rotation_sparse_gen(mat2euler(R), Nmax, rotD, indD)
         call sph_rotation_sparse_gen(mat2euler(R90), Nmax, rotD90, indD90)
         call sph_rotation_sparse_gen(mat2euler(Rx), Nmax, rotX, indX)
         call sph_rotation_sparse_gen(mat2euler(Ry), Nmax, rotY, indY)

         matrices%rotDs(1:las, i) = rotD
         matrices%indDs(1:las, :, i) = indD
         matrices%rotD90s(1:las, i) = rotD90
         matrices%indD90s(1:las, :, i) = indD90
         matrices%rotXs(1:las, i) = rotX
         matrices%indXs(1:las, :, i) = indX
         matrices%rotYs(1:las, i) = rotY
         matrices%indYs(1:las, :, i) = indY
      end do

   end subroutine gen_rotations

!****************************************************************************80

   subroutine int_points(mesh)
      type(mesh_struct) :: mesh
      real(dp), dimension(:, :), allocatable :: P
      real(dp), dimension(:), allocatable :: w
      real(dp) :: d, d2
      integer :: i1

      if (mesh%is_mesh == 1) then
         d = 0d0
         do i1 = 1, size(mesh%coord, 2)
            d2 = vlen(mesh%coord(:, i1))
            if (d2 > d) then
               d = d2
            end if
         end do
         d = 1.05d0*d
      else
         d = 1.05d0*(mesh%maxrad + mesh%radius(1))
      end if

      call sample_points(P, w, 20, 20)
      if (.not. allocated(mesh%P)) then
         allocate (mesh%P(size(P, 1), size(P, 2)), mesh%w(size(w)))
      end if

! These points should enclose the object
      mesh%P = P*d !* mesh%a
      mesh%w = w*d**2 !* mesh%a**2

   end subroutine int_points

!****************************************************************************80

   subroutine int_points_mie(mesh)
      type(mesh_struct) :: mesh
      real(dp), dimension(:, :), allocatable :: P
      real(dp), dimension(:), allocatable :: w
      real(dp) :: d

      d = 1.05d0*(mesh%a)

      call sample_points(P, w, 20, 20)
      if (.not. allocated(mesh%P)) then
         allocate (mesh%P(size(P, 1), size(P, 2)), mesh%w(size(w)))
      end if

! These points should enclose the object
      mesh%P = P*d !* mesh%a
      mesh%w = w*d**2 !* mesh%a**2

   end subroutine int_points_mie

!****************************************************************************80

   subroutine polarization()
      real(dp) :: k0(3), k1(3), R(3, 3)

      k0 = [0d0, 0d0, 1d0]
      k1 = matrices%khat
      R = rotate_a_to_b(k0, k1)

      matrices%E0hat = dcmplx(matmul(R, [1d0, 0d0, 0d0]))
      matrices%E90hat = dcmplx(matmul(R, [0d0, 1d0, 0d0]))

   end subroutine polarization

!****************************************************************************80
! B_lambda calculates the black body radiation intensity at
! wavelength lambda (m) and temperature T (K) using Planck's
! distribution.
   function B_lambda(lambda, T) result(I)
      real(dp), intent(in) :: T
      real(dp) :: I, lambda

      lambda = lambda
      I = (4d0*pi*hbar*cc**2d0/(lambda**5d0))/(exp(2d0*pi*hbar*cc/(lambda*k_b*T)) - 1d0)

   end function B_lambda

!****************************************************************************80
! find_k calculates reasonable values for n=bars wavenumbers
! in the approximate wavelength range lambda1-lambda2,
! with the priority that one wavelength in the range is at the
! maximum
   subroutine find_k()
      real(dp) :: dlmbda, lmax, fix
      real(dp), dimension(:), allocatable :: c, diff, absdif
      integer :: i, n, imin

! c(i) are the centers of histogram bars (wavelength)
      n = matrices%bars
      allocate (c(n))
      allocate (diff(n))
      allocate (absdif(n))

! Calculate lambda_max of distribution and width of bars
      lmax = b/matrices%temp
      dlmbda = (matrices%lambda2 - matrices%lambda1)/n

! Calculate c(i) and test vectors for adjusting the centers
      do i = 1, n
         c(i) = matrices%lambda1 + (2d0*i - 1d0)*dlmbda/2d0
         diff(i) = lmax - c(i)
         absdif(i) = abs(diff(i))
      end do

! Find the minimum distance from lambda_max
      imin = minloc(absdif, 1)

! fix = amount to slide c(i)
      if (imin > 0) then
         fix = diff(imin)
      else
         fix = 0d0
      end if

! Test whether or not the slide takes the minimum to negative wavelength
      if (diff(imin) < 0) then
         if (diff(imin) + c(1) < 0) then
            fix = -matrices%lambda1
         end if
      end if

! Calculate the final wavelengths and the corresponding wave numbers
      do i = 1, n
         c(i) = c(i) + fix
         mesh%ki(i) = 2d0*pi/c(i)
      end do

   end subroutine find_k

!****************************************************************************80
! Calculate the relative amplitudes in the field
   subroutine calc_E_rel()
      real(dp) :: lambda, N, sm, sm2, dlmbda, M
      integer :: i

      if (matrices%bars == 1) then
         matrices%E_rel = 1d0
      else
         dlmbda = (matrices%lambda2 - matrices%lambda1)/matrices%bars
         sm = 0.0d0
         sm2 = 0.0d0
         M = 0.5d0*epsilon*cc*(matrices%E)**2
         do i = 1, matrices%bars
            lambda = 2d0*pi/mesh%ki(i)
            sm2 = sm2 + B_lambda(lambda, matrices%temp)
         end do
         N = M/(sm2*dlmbda)
         do i = 1, matrices%bars
            lambda = 2d0*pi/mesh%ki(i)
            matrices%E_rel(i) = sqrt(2d0*N*B_lambda(lambda, matrices%temp)*dlmbda &
                                     /(epsilon*cc))
            sm = sm + matrices%E_rel(i)
         end do
         matrices%E_rel = matrices%E_rel/sm
      end if
       
   end subroutine calc_E_rel

!****************************************************************************80
! Setup the wavelength band and all matrices in it
   subroutine setup_band()
      if (matrices%waves == 'bnd') then
         call find_k()
         call calc_E_rel()
      else
         call band_no_blackbody()
      end if

   end subroutine setup_band

!****************************************************************************80
! Set the particle as Draine&Lee1984 astrosilicate
   subroutine band_astrosilicate()
      real(dp) :: lmax
      real(dp), dimension(:), allocatable :: c, diff, absdif, w, epsr, epsi
      integer :: i, n, size_param, ind

      call read_mesh()

      size_param = size(mesh%params, 1)
      if (allocated(mesh%params)) deallocate (mesh%params)
      allocate (mesh%params(size_param, matrices%bars))

! Read the astrosilicate data
      open (unit=15, file="examples/eps_Sil", status='old', &
            access='sequential', form='formatted', action='read')

      read (15, *) n
      read (15, *)
      allocate (w(n), epsr(n), epsi(n))

      do i = 1, n
         read (15, *) w(i), epsr(i), epsi(i)
         epsr(i) = epsr(i) + 1
      end do

      close (15)

! Find the closest values for dielectric constant according to
! the wavelength band
      do i = 1, matrices%bars
         ind = minloc(abs(w - 2d0*pi/mesh%ki(i)/1d-6), 1)
         mesh%params(:, i) = dcmplx(epsr(ind), epsi(ind))
      end do

   end subroutine band_astrosilicate

!****************************************************************************80
! Setup the wavelength band and all matrices inv A in it
   subroutine band_no_blackbody()
      real(dp) :: lmax
      real(dp), dimension(:), allocatable :: c, diff, absdif
      integer :: i, n

      n = matrices%bars
      allocate (c(n))
      allocate (diff(n))
      allocate (absdif(n))

! Calculate lambda_max of distribution and width of bars
      lmax = b/matrices%temp
! dlmbda = (matrices%lambda2-matrices%lambda1)/n

      if (matrices%waves == 'inv') then
         ! print*, 'inv = ', matrices%waves
         call linspace(1d0/(matrices%lambda1), 1d0/(matrices%lambda2), n, c)
         do i = 1, n
            c(i) = 1d0/c(i)
         end do
      else if (matrices%waves == 'log') then
         ! print*, 'log = ', matrices%waves
         call linspace(dlog10(matrices%lambda1), dlog10(matrices%lambda2), n, c)
         do i = 1, n
            c(i) = 10**c(i)
         end do
      else if(matrices%waves == 'isrf') then
         call linspace(1d0/(matrices%lambda1), 1d0/(matrices%lambda2), n, c)
         do i = 1, n
            c(i) = 1d0/c(i)
         end do
      else
         ! print*, matrices%waves
         call linspace(matrices%lambda1, matrices%lambda2, n, c)
      end if

! Calculate the final wavelengths and the corresponding wave numbers
      do i = 1, n
         mesh%ki(i) = 2d0*pi/c(i)
      end do
      matrices%E_rel = 1d0

      if(matrices%waves == 'isrf') then
         call calc_E_rel()
         if (2*pi/mesh%ki(1)/1d-6 < 0.2d0) matrices%E_rel(1) = maxval(matrices%E_rel)*2d0
      end if 

   end subroutine band_no_blackbody

end module setup
