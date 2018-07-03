module mueller
   use shapebeam

   implicit none

contains

!****************************************************************************80
! Compute the mueller matrix according to the mode chosen.
   subroutine compute_mueller(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: ii
      real(dp) :: E
      CHARACTER(LEN=120) :: mueller_out
      real(dp), dimension(:, :), allocatable :: S

      if (trim(matrices%mueller_mode) == 'none') return

      mueller_out = 'out/mueller'//'-'//trim(matrices%mueller_mode)//trim(matrices%out)
      if (file_exists(mueller_out)) then
         print *, ' Mueller matrix already exists, quitting...'
         stop
      end if

! Choose wavelength
      ii = matrices%whichbar
      if (ii == 0) ii = 1

      E = matrices%E_rel(ii)*matrices%E
      mesh%k = mesh%ki(ii)

      select case (trim (matrices%mueller_mode))
      case ('ave')
         call compute_ave_mueller(matrices, mesh, 180, 1, ii, E, S)
      case ('ori')
         call compute_log_mueller(matrices, mesh, 180, 90, ii, E, S)
      case ('perf_ori')
         call compute_aligned_mueller(matrices, mesh, 180, 90, ii, E, S)
      end select

      call write_mueller(matrices, S)
   end subroutine compute_mueller

!****************************************************************************80
! Compute orientation averaged Mueller matrices for given number of theta, phi.
! The phi-dependency is lost due to averaging, so no need to give N_phi /= 1.
   subroutine compute_ave_mueller(matrices, mesh, N_theta, N_phi, ii, E, S)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, ii, halton_init, N_points, N_theta, N_phi
      real(dp) :: E, vec(3)
      real(dp), dimension(:, :), allocatable :: S
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      N_points = 720 ! Number of averaging directions

      allocate (S(N_theta*N_phi, 18))
      S = 0d0

      matrices%R = eye(3)
      halton_init = 0

      do i = 1, N_points
         vec(1) = 1d0
         vec(2) = acos(2d0*halton_seq(halton_init + i, 2) - 1d0)
         vec(3) = halton_seq(halton_init + i, 3)*2d0*pi
         matrices%khat = [dsin(vec(2))*dcos(vec(3)), dsin(vec(2))*dsin(vec(3)), dcos(vec(2))]

         matrices%khat = -matrices%khat/vlen(matrices%khat)
         matrices%R = rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0])

         S = S + update_mueller(matrices, mesh, N_theta, N_phi, ii, E, p, q, p90, q90)/N_points
         call print_bar(i, N_points)
      end do
   end subroutine compute_ave_mueller

!****************************************************************************80
! Compute the Mueller matrix of a perfectly aligned particle. The alignment
! direction is fixed so that the major axis of inertia is always in the
! +z-direction.
   subroutine compute_aligned_mueller(matrices, mesh, N_theta, N_phi, ii, E, S, J)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, ii, N_points, N_theta, N_phi
      real(dp) :: E, vec(3), phi, R0(3, 3), Qt(3, 3), omega(3), aproj(3), RR(3, 3), theta
      real(dp), dimension(:, :), allocatable :: S
      real(dp), dimension(3), optional :: J
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      N_points = 500 ! Number of points to calculate the perfect orientations

      allocate (S(N_theta*N_phi, 18))
      S = 0d0

! Hard coded orientation data
      omega = [1d0, 0d0, 0d0]
      R0 = rotate_a_to_b(matrices%P(:, 3), omega)
      Qt = matmul(R0, matrices%P)
      aproj = [Qt(1, 3), Qt(2, 3), 0d0]
      aproj = aproj/vlen(aproj)
      phi = dacos(aproj(1))

      do i = 1, N_points
         theta = dble(i - 1)*2d0*pi/(N_points - 1)
         RR = matmul(transpose(R_aa(omega, theta)), transpose(R0))

         matrices%khat = -matmul(RR, [0d0, 0d0, 1d0])
         matrices%khat = matrices%khat/vlen(matrices%khat)

         matrices%R = rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0])
         matrices%R = matmul(matrices%R,R_aa(matrices%khat, phi))

         S = S + update_mueller(matrices, mesh, N_theta, N_phi, ii, E, p, q, p90, q90)/N_points
         call print_bar(i, N_points)
      end do
   end subroutine compute_aligned_mueller

!****************************************************************************80
! Compute the Mueller matrix from the data from the dynamical simulation. The
! particle alignment state is taken from the logged orientations. Thus, the
! situation may be aligned or not.
   subroutine compute_log_mueller(matrices, mesh, N_theta, N_phi, ii, E, S)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, ii, N_points, N_theta, N_phi
      real(dp) :: E, RR(3, 3)
      real(dp), dimension(:, :), allocatable :: S
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      call read_log(matrices, mesh, 5000)
      N_points = size(matrices%RRR, 3)

      allocate (S(N_theta*N_phi, 18))
      S = 0d0

      do i = 1, N_points
         matrices%R = transpose(matmul(matrices%R_al, matrices%RRR(:, :, i)))
         matrices%khat = matmul(transpose(matrices%R), [0d0, 0d0, 1d0])
         matrices%khat = -matrices%khat/vlen(matrices%khat)
         matrices%R = transpose(rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0]))

         S = S + update_mueller(matrices, mesh, N_theta, N_phi, ii, E, p, q, p90, q90)/N_points
         call print_bar(i, N_points)
      end do
   end subroutine compute_log_mueller

!****************************************************************************80
! Updates the Mueller matrix average.
   function update_mueller(matrices, mesh, N_theta, N_phi, ii, E, p, q, p90, q90) result(S)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: N_theta, N_phi, ii
      real(dp) :: E
      real(dp), dimension(:, :), allocatable :: S
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      call rot_setup(matrices, mesh)
      call scattered_fields(matrices, mesh, E, p, q, p90, q90, ii)
      call mueller_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_theta, N_phi, S)
   end function update_mueller

! Below are WIP functionalities for the SOCpol application
!****************************************************************************80

   subroutine scattering_extinction_matrices(matrices, mesh, points)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, ind, ind2, ii, N_points, N_size, N_ia
      real(dp) :: inc_angles(2), K(4, 4), Cext_al, Csca_al, Cext_ave, Csca_ave
      complex(dp) :: Kevals(4), Kevecs(4, 4)
      real(dp), dimension(:, :), allocatable :: S_ave, S_al, SS, K_ave, K_al, KK, points
      CHARACTER(LEN=80) :: mueller_out, extinction_out

      N_points = size(points, 2) ! Number of points to calculate the perfect orientations

      mueller_out = 'mueller'
      extinction_out = 'extinction_matrix'

      inc_angles = [90d0, 180d0]
      allocate (SS(2*N_points*size(mesh%ki, 1)*size(inc_angles), 20))
      allocate (KK(2*size(mesh%ki, 1)*size(inc_angles), 20))
      SS = 0d0
      KK = 0d0
      allocate (S_al(N_points, 18), K_al(1, 18), S_ave(N_points, 18), K_ave(1, 18))

      ind = 0
      ind2 = 0
      do N_size = 1, size(mesh%ki, 1)
! Choose wavelength
         mesh%k = mesh%ki(N_size)
         do N_ia = 1, size(inc_angles, 1)            
            call mueller_ave(matrices, mesh, S_ave, K_ave, Csca_ave, Cext_ave, points, N_size)
            call mueller_align(matrices, mesh, S_al, K_al, Csca_al, Cext_al, points, N_size, inc_angles(N_ia))

            do i = 1, N_points
               ind = ind + 1
               SS(ind, 1) = N_size
               SS(ind, 2) = N_ia
               SS(ind, 3) = i
               SS(ind, 4:19) = S_ave(i, 3:18)
               SS(ind, 20) = Csca_ave
               ind = ind + 1
               SS(ind, 1) = N_size
               SS(ind, 2) = N_ia
               SS(ind, 3) = i
               SS(ind, 4:19) = S_al(i, 3:18)
               SS(ind, 20) = Csca_al

               call print_bar(ind, 2*size(mesh%ki, 1)*size(inc_angles)*N_points)
            end do
            ind2 = ind2 + 1
            KK(ind2, 1) = N_size
            KK(ind2, 2) = N_ia
            KK(ind2, 3) = 1
            KK(ind2, 4:19) = K_ave(1, 3:18)
            KK(ind2, 20) = Cext_ave
            ind2 = ind2 + 1
            KK(ind2, 1) = N_size
            KK(ind2, 2) = N_ia
            KK(ind2, 3) = 1
            KK(ind2, 4:19) = K_al(1, 3:18)
            KK(ind2, 20) = Cext_al

         end do
      end do
      call write_RT_matrix(SS, mueller_out, 1)
      call write_RT_matrix(KK, extinction_out, 2)

   end subroutine scattering_extinction_matrices

!****************************************************************************80

   subroutine mueller_ave(matrices, mesh, SS, KK, Csca_out, Cext_out, points, ii)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, ii, N_avgs, halton_init, nm, Nmax
      real(dp) :: E, vec(3), k_sph(3)
      real(dp) :: Cext, Cabs, Csca, Csca_out, Cext_out
      real(dp), dimension(:, :), allocatable :: S, SS, K, KK, points
      complex(dp), dimension(:), allocatable :: a_in, b_in
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      SS = 0d0
      KK = 0d0
      Csca = 0d0 
      Cext = 0d0 

      matrices%R = eye(3)
      N_avgs = 360 ! Number of averaging directions
      halton_init = 0
      E = matrices%E_rel(ii)*matrices%E
      Nmax = matrices%Nmaxs(ii)

      do i = 1, N_avgs
         vec(1) = 1d0
         vec(2) = acos(2d0*halton_seq(halton_init + i, 2) - 1d0)
         vec(3) = halton_seq(halton_init + i, 3)*2d0*pi
         matrices%khat = [dsin(vec(2))*dcos(vec(3)), dsin(vec(2))*dsin(vec(3)), dcos(vec(2))]

         matrices%khat = -matrices%khat/vlen(matrices%khat)
         k_sph = cart2sph(matrices%khat)
         matrices%R = rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0])

         call rot_setup(matrices, mesh)
         call incident_fields(matrices, mesh, E, a_in, b_in, ii)
         call scattered_fields(matrices, mesh, E, p, q, p90, q90, ii)
         if (allocated(S)) deallocate (S, K)
         call scattering_matrix_coeff(p, q, p90, q90, dcmplx(mesh%ki(ii)), size(points, 2), points, S)
         call extinction_matrix_coeff(p, q, p90, q90, dcmplx(mesh%ki(ii)), k_sph(2), k_sph(3), K)
         call cross_sections(p, q, a_in, b_in, dcmplx(mesh%ki(ii)), Nmax, Cext, Csca, Cabs)
         SS = SS + S/N_avgs
         KK = KK + K/N_avgs
         Csca_out = Csca_out + Csca/N_avgs
         Cext_out = Cext_out + Cext/N_avgs
      end do

   end subroutine mueller_ave

!****************************************************************************80

   subroutine mueller_align(matrices, mesh, SS, KK, Csca_out, Cext_out, points, ii, psi, xi_in)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      integer :: i, j, ii, Ntheta, Nphi, Nmax
      real(dp) :: E, omega(3), al_direction(3), theta, phi, &
                  aproj(3), k_sph(3), psi, x_B(3), n_phi(3)
      real(dp) :: Cext, Cabs, Csca, Csca_out, Cext_out, xi
      real(dp), dimension(3,3) :: Rtheta, Rphi, Qt, R0
      real(dp), optional :: xi_in
      real(dp), dimension(:, :), allocatable :: S, SS, K, KK, points
      complex(dp), dimension(:), allocatable :: a_in, b_in
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      SS = 0d0
      KK = 0d0
      Csca = 0d0 
      Cext = 0d0 

      Ntheta = 36 ! Number of averaging directions
      Nphi = 36 ! Number of precession averaging directions
      if(.NOT. present(xi_in)) then
         xi = 0d0
         Nphi = 1
      else
         xi = xi_in
      end if 

! First, set up B in the direction of psi (B still lies in the xz-plane)
      psi = matrices%B_psi
      x_B = [0d0, 1d0, 0d0]

! Rotation axis for precession averaging
      n_phi = matmul(R_aa(x_B, psi), [0d0,1d0,0d0]) ! Is actually B

      E = matrices%E_rel(ii)*matrices%E
      omega = matmul(R_aa([1d0,0d0,0d0], psi+xi),[0d0,0d0,1d0])
      R0 = rotate_a_to_b(matrices%P(:, 3), omega)
      Qt = matmul(R0, matrices%P)

      Nmax = matrices%Nmaxs(ii)
      
      do j = 1, Nphi
         phi = dble(j - 1)*2d0*pi/dble(Nphi)
         Rphi = R_aa(n_phi,phi)
         do i = 1, Ntheta
            theta = dble(i - 1)*2d0*pi/dble(Ntheta - 1)
            Rtheta = matmul(transpose(R_aa(omega, theta)), transpose(R0))
            Rtheta = matmul(transpose(Rphi),Rtheta)

            matrices%khat = -matmul(Rtheta, [0d0, 0d0, 1d0])
            matrices%khat = matrices%khat/vlen(matrices%khat)
            k_sph = cart2sph(matrices%khat)

            matrices%R = rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0])

            call rot_setup(matrices, mesh)
            call incident_fields(matrices, mesh, E, a_in, b_in, ii)
            call scattered_fields(matrices, mesh, E, p, q, p90, q90, ii)
            if (allocated(S)) deallocate (S, K)
            call scattering_matrix_coeff(p, q, p90, q90, dcmplx(mesh%ki(ii)), size(points, 2), points, S)
            call extinction_matrix_coeff(p, q, p90, q90, dcmplx(mesh%ki(ii)), k_sph(2), k_sph(3), K)
            call cross_sections(p, q, a_in, b_in, dcmplx(mesh%ki(ii)), Nmax, Cext, Csca, Cabs)
            SS = SS + S/Nphi/Ntheta
            KK = KK + K/Nphi/Ntheta
            Csca_out = Csca_out + Csca/Nphi/Ntheta
            Cext_out = Cext_out + Cext/Nphi/Ntheta
         end do 
      end do

   end subroutine mueller_align

end module mueller
