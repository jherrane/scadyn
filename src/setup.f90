module setup
! Copyright (c) 2018 Joonas Herranen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common
   use io
   use mie
   use translations
   use projection

   implicit none
contains

!****************************************************************************80

   subroutine init_geometry()
      real(dp) :: ka
      integer :: i
      complex(dp) :: k

!* TETRAHEDRAL MESH *************************************************
      if (use_mie .NE. 1) then
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
         call int_points()

         if (matrices%Tmat == 0) print *, 'Initializing FFT... '
         if (allocated(mesh%nodes)) deallocate (mesh%nodes, mesh%elem_box, mesh%tetras)
         call build_grid(mesh) ! geometry
         call build_box(mesh) ! geometry
         call tetras_in_cubes(mesh) ! geometry

         if (matrices%Tmat == 0 .AND. debug == 1) then
            print *, '   Grid size            = ', (/mesh%Nx, mesh%Ny, mesh%Nz/)
            print *, '   Delta grid           = ', real(mesh%delta)
            print *, '   Number of cubes      = ', mesh%N_cubes
            print *, '   Delta cubes          = ', real(mesh%box_delta)
            print *, '   Elems. in cube (max) = ', mesh%N_tet_cube
            print *, 'Done'
         end if

         if (matrices%Tmat == 0) call construct_projectors()

! Check whether using maxval is good or not
         do i = 1, matrices%bars
            ka = mesh%ki(i)*mesh%a
            matrices%Nmaxs(i) = truncation_order(ka)
         end do
      else ! MIE *************************************************************
         call int_points_mie()
         do i = 1, matrices%bars
            ka = mesh%ki(i)*mesh%a
            matrices%Nmaxs(i) = truncation_order(ka)
         end do

         allocate (mesh%eps(1), mesh%refr(1))
         mesh%refr = dcmplx(matrices%refr_r, matrices%refr_i)
         write (*, '(2(A,F5.3))') '    Refractive index     =   ', matrices%refr_r, ' + i', matrices%refr_i
      end if
   end subroutine init_geometry

!****************************************************************************80

   subroutine construct_projectors()
      integer :: i, M_box
      real(dp) :: k_orig

      k_orig = mesh%k ! Keep the original just in case
      M_box = size(mesh%elem_box, 1)

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
         if (size(mesh%refrs, 2) > 1) then
            mesh%eps = dcmplx(real(mesh%refrs(:,i))**2-imag(mesh%refrs(:,i))**2, &
                  2d0*real(mesh%refrs(:,i))*imag(mesh%refrs(:,i)))
         else
            mesh%eps = dcmplx(real(mesh%refrs(:,1))**2-imag(mesh%refrs(:,1))**2, &
                  2d0*real(mesh%refrs(:,1))*imag(mesh%refrs(:,1)))
         end if
         
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
      type(data) :: matrices
      type(mesh_struct) :: mesh
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

   subroutine int_points()
      real(dp), dimension(:, :), allocatable :: P
      real(dp), dimension(:), allocatable :: w
      real(dp) :: d, d2
      integer :: i1

      d = 0d0
      do i1 = 1, size(mesh%node, 2)
         d2 = vlen(mesh%node(:, i1))
         if (d2 > d) then
            d = d2
         end if
      end do
      d = 1.05d0*d

      call sample_points(P, w, 20, 20)
      if (.not. allocated(mesh%P)) then
         allocate (mesh%P(size(P, 1), size(P, 2)), mesh%w(size(w)))
      end if

! These points should enclose the object
      mesh%P = P*d !* mesh%a
      mesh%w = w*d**2 !* mesh%a**2

   end subroutine int_points

!****************************************************************************80

   subroutine int_points_mie()
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
         mesh%ki(i) = matrices%ref_med*2d0*pi/c(i)
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
      real(dp), dimension(:), allocatable :: c, diff, absdif, w, epsr, epsi, &
                                             refr_r, refr_i
      integer :: i, n, size_param, ind

      call read_mesh()

      size_param = size(mesh%refrs, 1)
      if (allocated(mesh%refrs)) deallocate (mesh%refrs)
      allocate (mesh%refrs(size_param, matrices%bars))

! Read the astrosilicate data
      open (unit=15, file="examples/eps_Sil", status='old', &
            access='sequential', form='formatted', action='read')

      read (15, *) n
      read (15, *)
      allocate (w(n), epsr(n), epsi(n), refr_r(n), refr_i(n))

      do i = 1, n
         read (15, *) w(i), epsr(i), epsi(i), refr_r(i), refr_i(i)
         refr_r(i) = refr_r(i) + 1
      end do

      close (15)

! Find the closest values for dielectric constant according to
! the wavelength band
      do i = 1, matrices%bars
         ind = minloc(abs(w - 2d0*pi/mesh%ki(i)/1d-6), 1)
         mesh%refrs(:, i) = dcmplx(refr_r(ind), refr_i(ind))
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

      if (matrices%waves == 'inv') then
         call linspace(1d0/(matrices%lambda1), 1d0/(matrices%lambda2), n, c)
         do i = 1, n
            c(i) = 1d0/c(i)
         end do
      else if (matrices%waves == 'log') then
         call linspace(dlog10(matrices%lambda1), dlog10(matrices%lambda2), n, c)
         do i = 1, n
            c(i) = 10**c(i)
         end do
      else 
         call linspace(matrices%lambda1, matrices%lambda2, n, c)
      end if

! Calculate the final wavelengths and the corresponding wave numbers
      do i = 1, n
         mesh%ki(i) = matrices%ref_med*2d0*pi/c(i)
      end do
      matrices%E_rel = 1d0

      if(matrices%waves == 'isrf') then
         call calc_E_rel()
         if (2*pi/mesh%ki(1)/1d-6 < 0.2d0) matrices%E_rel(1) = maxval(matrices%E_rel)*0.5d0 
      end if 

   end subroutine band_no_blackbody


!****************************************************************************80

   subroutine diagonalize_inertia()
      integer :: i, imin, imax, imid, negs

! Test if I is "almost diagonal", which causes diasym to do funny stuff
      if (dabs(mesh%I(1, 1)*mesh%I(2, 2)*mesh%I(3, 3)) > 1d6*dabs(mesh%I(1, 2)*mesh%I(1, 3)*mesh%I(2, 3))) then
         if (use_mie .NE. 1) then
            print *, 'FALLBACK MODE IN DIAGONALIZATION OF INERTIA TENSOR'
            matrices%Ip = 0d0
            imin = minloc([mesh%I(1, 1), mesh%I(2, 2), mesh%I(3, 3)], 1)
            imax = maxloc([mesh%I(1, 1), mesh%I(2, 2), mesh%I(3, 3)], 1)
            if (imin == 1) then
               if (imax == 2) then
                  imid = 3
                  matrices%P = reshape([1d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 1d0, 0d0], [3, 3])
               end if
               if (imax == 3) then
                  imid = 2
                  matrices%P = eye(3)
               end if
            else if (imin == 2) then
               if (imax == 1) then
                  imid = 3
                  matrices%P = reshape([0d0, 1d0, 0d0, 0d0, 0d0, 1d0, 1d0, 0d0, 0d0], [3, 3])
               end if
               if (imax == 3) then
                  imid = 1
                  matrices%P = reshape([0d0, 1d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0, 1d0], [3, 3])
               end if
            else
               if (imax == 1) then
                  imid = 2
                  matrices%P = reshape([0d0, 0d0, 1d0, 0d0, 1d0, 1d0, 1d0, 0d0, 0d0], [3, 3])
               end if
               if (imax == 2) then
                  imid = 1
                  matrices%P = reshape([0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 1d0, 0d0], [3, 3])
               end if
            end if
            matrices%Ip(1) = minval([mesh%I(1, 1), mesh%I(2, 2), mesh%I(3, 3)])
            matrices%Ip(3) = maxval([mesh%I(1, 1), mesh%I(2, 2), mesh%I(3, 3)])
            matrices%Ip(2) = mesh%I(imid, imid)
         else
            matrices%Ip = [mesh%I(1, 1), mesh%I(2, 2), mesh%I(3, 3)]
         end if
      else
! P will be the rotation matrix between principal axes and laboratory
! axes, so that diag(Ip) = P'*I*P
         matrices%P = mesh%I
         call diasym(matrices%P, matrices%Ip)
      end if

! Probably the choice in DDSCAT: force eigenvectors to have mostly positive
! components
      do i = 1, 3
         negs = 0
         if (matrices%P(1, i) < 0d0) negs = negs + 1
         if (matrices%P(2, i) < 0d0) negs = negs + 1
         if (matrices%P(3, i) < 0d0) negs = negs + 1
         if (negs >= 2) matrices%P(:, i) = -matrices%P(:, i)
      end do

      matrices%I = 0d0
      matrices%I_inv = 0d0

      forall (i=1:3) matrices%I(i, i) = matrices%Ip(i)
      forall (i=1:3) matrices%I_inv(i, i) = 1d0/matrices%Ip(i)
      
      mesh%alpha = matrices%Ip/(2d0/5d0*mesh%rho*mesh%V*mesh%a**2)

   end subroutine diagonalize_inertia

!****************************************************************************80
! Computes interstellar environment values needed to solve the equations of 
! motion in the alignment problem
   subroutine interstellar_env()
      integer :: i
      real(dp) :: urad, lambdamean
      matrices%wT = (15d0/(8d0*pi*mesh%alpha(3))*k_b*matrices%Td/ &
                     (mesh%rho*mesh%a**5))**0.5
      matrices%TDG = (2d0*mesh%alpha(3)*mesh%rho*mesh%a**2)/ &
                     (5d0*matrices%Kw*matrices%B_len**2)/(4*pi/(mu))

      matrices%Tdrag = (pi*mesh%alpha(3)*mesh%rho*mesh%a)/(3*mesh%drag*matrices%nH* &
                        (2d0*pi*1.67d-27*k_b*matrices%Tgas)**0.5)
      urad = 0d0
      lambdamean = 0d0
      do i = 1,matrices%bars
         urad = urad + (matrices%E_rel(i)*matrices%E)**2/sqrt(mu/epsilon)/2d0/cc
         lambdamean = lambdamean + 1d0*pi/mesh%ki(i)
      end do 
      
      lambdamean = lambdamean/matrices%bars
      matrices%M =  urad*lambdamean*mesh%a**2*matrices%Tdrag/ &
                     (2d0*matrices%Ip(3)*matrices%wT)

      if(debug==1) then
         print*, 'The interstellar environment values'
         print*, 'lambda_mean = ', lambdamean
         print*, 'u_rad = ', urad
         print*, 'w_T = ', matrices%wT
         print*, 'T_DG (years) = ', matrices%TDG/(60*60*24*365)
         print*, 'T_drag (years) = ', matrices%Tdrag/(60*60*24*365)
         print*, 'M = ', matrices%M
      end if

   end subroutine interstellar_env

!****************************************************************************80

   subroutine mie_params()
      mesh%V = 4d0/3d0*pi*mesh%a**3
      mesh%mass = mesh%V*mesh%rho
      mesh%CM = [0d0, 0d0, 0d0]
      mesh%I = 2d0/5d0*mesh%mass*mesh%a**2*eye(3)

      matrices%P = eye(3)
   end subroutine mie_params

!****************************************************************************80
! Mass parameters for VIE-mesh
! output: mesh%V = volume of the mesh
!         mesh%mass = mass of the mesh
!         mesh%CM = coordinates of the center of mass
!         mesh%I = The moment of inertia tensor
!****************************************************************************80
   subroutine vie_params()
      real(dp) :: rho, V, totV, mass, totMass, detJ, &
                  a, b, c, ap, bp, cp
      integer  :: i1
      real(dp), dimension(3) :: p0, p1, p2, p3, COM, CM
      real(dp), dimension(4) :: x, y, z
      real(dp), dimension(3, 3) :: I, E

      rho = mesh%rho ! currently density not in tetrahedral element data

      I = 0.d0
      totV = 0.d0
      totMass = 0.d0
      CM = 0.d0

      do i1 = 1, mesh%N_tet
         ! Get vertices
         p0 = mesh%node(:, mesh%elem(1, i1))
         p1 = mesh%node(:, mesh%elem(2, i1))
         p2 = mesh%node(:, mesh%elem(3, i1))
         p3 = mesh%node(:, mesh%elem(4, i1))
         COM = [(p0(1) + p1(1) + p2(1) + p3(1))/4d0, &
                (p0(2) + p1(2) + p2(2) + p3(2))/4d0, &
                (p0(3) + p1(3) + p2(3) + p3(3))/4d0]

         V = abs(dot_product((p0 - p3), &
                             crossRR((p1 - p3), (p2 - p3))))/6d0

         x = [p0(1) - COM(1), p1(1) - COM(1), &
              p2(1) - COM(1), p3(1) - COM(1)]
         y = [p0(2) - COM(2), p1(2) - COM(2), &
              p2(2) - COM(2), p3(2) - COM(2)]
         z = [p0(3) - COM(3), p1(3) - COM(3), &
              p2(3) - COM(3), p3(3) - COM(3)]

         detJ = (x(2) - x(1))*((y(3) - y(1))*(z(4) &
                                              - z(1)) - (y(4) - y(1))*(z(3) - z(1))) - &
                (x(3) - x(1))*((y(2) - y(1))*(z(4) &
                                              - z(1)) - (y(4) - y(1))*(z(2) - z(1))) + &
                (x(4) - x(1))*((y(2) - y(1))*(z(3) &
                                              - z(1)) - (y(3) - y(1))*(z(2) - z(1)))

         a = rho*detJ*diagterm(y, z)/60d0
         b = rho*detJ*diagterm(x, z)/60d0
         c = rho*detJ*diagterm(x, y)/60d0

         ap = rho*detJ*offterm(y, z)/120d0
         bp = rho*detJ*offterm(x, z)/120d0
         cp = rho*detJ*offterm(x, y)/120d0

         E = reshape([a, -bp, -cp, -bp, b, -ap, -cp, -ap, c], &
                     [3, 3])
         mass = rho*V

         totV = totV + V
         totMass = totMass + mass
         CM = CM + mass*COM
         I = I + E + mass*(dot_product(COM, COM)*eye(3) - &
                           real_outer_product(COM, COM))

      end do

      mesh%V = totV
      mesh%mass = totMass
      mesh%CM = CM/totMass
      mesh%I = I

      mesh%node(1, :) = mesh%node(1, :) - mesh%CM(1)
      mesh%node(2, :) = mesh%node(2, :) - mesh%CM(2)
      mesh%node(3, :) = mesh%node(3, :) - mesh%CM(3)

   end subroutine vie_params

!****************************************************************************80

   function diagterm(y, z) result(a)

      real(dp)                 :: a

      real(dp), dimension(4) :: y, &
                                z

      a = (y(1)*y(1) &
           + y(1)*y(2) + y(2)*y(2) &
           + y(1)*y(3) + y(2)*y(3) + y(3)*y(3) &
           + y(1)*y(4) + y(2)*y(4) + y(3)*y(4) + y(4)*y(4) &
           + z(1)*z(1) &
           + z(1)*z(2) + z(2)*z(2) &
           + z(1)*z(3) + z(2)*z(3) + z(3)*z(3) &
           + z(1)*z(4) + z(2)*z(4) + z(3)*z(4) + z(4)*z(4))

   end function diagterm

!****************************************************************************80

   function offterm(y, z) result(ap)
      real(dp)                 :: ap
      real(dp), dimension(4) :: y, z

      ap = (2*y(1)*z(1) &
            + y(2)*z(1) &
            + y(3)*z(1) &
            + y(4)*z(1) &
            + y(1)*z(2) &
            + 2*y(2)*z(2) &
            + y(3)*z(2) &
            + y(4)*z(2) &
            + y(1)*z(3) &
            + y(2)*z(3) &
            + 2*y(3)*z(3) &
            + y(4)*z(3) &
            + y(1)*z(4) &
            + y(2)*z(4) &
            + y(3)*z(4) &
            + 2*y(4)*z(4))

   end function offterm

end module setup
