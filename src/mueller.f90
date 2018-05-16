module mueller
   use shapebeam

   implicit none

contains

!****************************************************************************80
! Compute the mueller matrix according to the mode chosen.
   subroutine compute_mueller()
      integer :: ii
      real(dp) :: E
      CHARACTER(LEN=80) :: mueller_out
      real(dp), dimension(:, :), allocatable :: S

      if (trim(matrices%mueller_mode) == 'none') return

      mueller_out = trim(matrices%mueller)//'-'//trim(matrices%mueller_mode)
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
         call compute_ave_mueller(180, 1, ii, E, S)
      case ('ori')
         call compute_log_mueller(180, 90, ii, E, S)
      case ('perf_ori')
         call compute_aligned_mueller(180, 90, ii, E, S)
      end select

      call write_mueller(S, mueller_out)
   end subroutine compute_mueller

!****************************************************************************80
! Compute orientation averaged Mueller matrices for given number of theta, phi.
! The phi-dependency is lost due to averaging, so no need to give N_phi /= 1.
   subroutine compute_ave_mueller(N_theta, N_phi, ii, E, S)
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
         matrices%R_fixk = transpose(rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0]))

         S = S + update_mueller(N_theta, N_phi, ii, E, p, q, p90, q90)/N_points
         call print_bar(i, N_points)
      end do
   end subroutine compute_ave_mueller

!****************************************************************************80
! Compute the Mueller matrix of a perfectly aligned particle. The alignment
! direction is fixed so that the major axis of inertia is always in the
! +z-direction.
   subroutine compute_aligned_mueller(N_theta, N_phi, ii, E, S)
      integer :: i, ii, N_points, N_theta, N_phi
      real(dp) :: E, vec(3), phi, R0(3, 3), Qt(3, 3), omega(3), aproj(3), RR(3, 3), theta
      real(dp), dimension(:, :), allocatable :: S
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

         matrices%khat = matmul(RR, [0d0, 0d0, 1d0])
         matrices%khat = -matrices%khat/vlen(matrices%khat)
         matrices%R_fixk = transpose(rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0]))
         matrices%R = R_aa(matrices%khat, phi)

         S = S + update_mueller(N_theta, N_phi, ii, E, p, q, p90, q90)/N_points
         call print_bar(i, N_points)
      end do
   end subroutine compute_aligned_mueller

!****************************************************************************80
! Compute the Mueller matrix from the data from the dynamical simulation. The
! particle alignment state is taken from the logged orientations. Thus, the
! situation may be aligned or not.
   subroutine compute_log_mueller(N_theta, N_phi, ii, E, S)
      integer :: i, ii, N_points, N_theta, N_phi
      real(dp) :: E, RR(3, 3)
      real(dp), dimension(:, :), allocatable :: S
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      call read_log(5000)
      N_points = size(matrices%RRR, 3)

      allocate (S(N_theta*N_phi, 18))
      S = 0d0

      do i = 1, N_points
         matrices%R = transpose(matmul(matrices%R_al, matrices%RRR(:, :, i)))
         matrices%khat = matmul(transpose(matrices%R), [0d0, 0d0, 1d0])
         matrices%khat = -matrices%khat/vlen(matrices%khat)
         matrices%R_fixk = transpose(rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0]))

         S = S + update_mueller(N_theta, N_phi, ii, E, p, q, p90, q90)/N_points
         call print_bar(i, N_points)
      end do
   end subroutine compute_log_mueller

!****************************************************************************80
! Updates the Mueller matrix average.
   function update_mueller(N_theta, N_phi, ii, E, p, q, p90, q90) result(S)
      integer :: N_theta, N_phi, ii
      real(dp) :: E
      real(dp), dimension(:, :), allocatable :: S
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      call rot_setup()
      call scattered_fields(E, p, q, p90, q90, ii)
      call mueller_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_theta, N_phi, S)
   end function update_mueller

! Below are WIP functionalities for the SOCpol application
!****************************************************************************80

   subroutine scattering_extinction_matrices(a_dist, points, al_direction)
      integer :: i, ind, ii, N_points, N_size, N_ia
      real(dp) :: al_direction(3), inc_angles(2), K(4, 4), Cext, Csca
      complex(dp) :: Kevals(4), Kevecs(4, 4)
      real(dp), dimension(:), allocatable :: a_dist
      real(dp), dimension(:, :), allocatable :: SS, SSS, KK, KKK, points
      CHARACTER(LEN=80) :: mueller_out, extinction_out

      N_points = size(points, 2) ! Number of points to calculate the perfect orientations

      mueller_out = 'mueller'
      extinction_out = 'extinction_matrix'

      inc_angles = [90d0, 180d0]
      allocate (SSS(N_points*size(a_dist, 1)*size(inc_angles), 20))
      allocate (KKK(size(a_dist, 1)*size(inc_angles), 20))
      SSS = 0d0
      KKK = 0d0
      allocate (SS(N_points, 18), KK(1, 18), )

      ind = 0
      do N_size = 1, size(a_dist, 1)
! Choose wavelength
         mesh%k = mesh%ki(N_size)
         do N_ia = 1, size(inc_angles, 1)
! We take every N_size as the critical size, below which nothing is aligned
            SS = 0d0
            KK = 0d0
            Csca = 0d0 
            Cext = 0d0 
            do ii = 1, size(a_dist, 1)
               if (a_dist(ii) <= a_dist(N_size)) then
                  call mueller_ave(SS, KK, Csca, Cext, points, ii)
               else
                  call mueller_align(SS, KK, Csca, Cext, points, ii, al_direction)
               end if
            end do

            SS = SS/size(a_dist, 1)
            KK = KK/size(a_dist, 1)

            do i = 1, N_points
               ind = ind + 1
               SSS(ind, 1) = N_size
               SSS(ind, 2) = N_ia
               SSS(ind, 3) = i
               SSS(ind, 4:19) = SS(i, 3:18)
               SSS(ind, 20) = Csca

               call print_bar(ind, size(a_dist, 1)*size(inc_angles)*N_points)
            end do

            KKK(ind/N_points, 1) = N_size
            KKK(ind/N_points, 2) = N_ia
            KKK(ind/N_points, 3) = 1
            KKK(ind/N_points, 4:19) = KK(1, 3:18)
            KKK(ind/N_points, 20) = Cext

         end do
      end do
      call write_RT_matrix(SSS, mueller_out, 1)
      call write_RT_matrix(KKK, extinction_out, 2)

   end subroutine scattering_extinction_matrices

!****************************************************************************80

   subroutine mueller_ave(SS, KK, Csca_out, Cext_out, points, ii)
      integer :: i, ii, N_avgs, halton_init, nm, Nmax
      real(dp) :: E, vec(3), k_sph(3)
      real(dp) :: Cext, Cabs, Csca, Csca_out, Cext_out
      real(dp), dimension(:, :), allocatable :: S, SS, K, KK, points
      complex(dp), dimension(:), allocatable :: a_in, b_in
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

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
         matrices%R_fixk = transpose(rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0]))

         call rot_setup()
         call incident_fields(E, a_in, b_in, ii)
         call scattered_fields(E, p, q, p90, q90, ii)
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

   subroutine mueller_align(SS, KK, Csca_out, Cext_out, points, ii, al_direction)
      integer :: i, ii, N_avgs, Nmax
      real(dp) :: E, omega(3), al_direction(3), theta, phi, RR(3, 3), Qt(3, 3), R0(3, 3), &
                  aproj(3), k_sph(3)
      real(dp) :: Cext, Cabs, Csca, Csca_out, Cext_out
      real(dp), dimension(:, :), allocatable :: S, SS, K, KK, points
      complex(dp), dimension(:), allocatable :: a_in, b_in
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      N_avgs = 36 ! Number of averaging directions

      E = matrices%E_rel(ii)*matrices%E
      omega = al_direction
      omega = omega/vlen(omega)
      R0 = rotate_a_to_b(matrices%P(:, 3), omega)
      Qt = matmul(R0, matrices%P)
      aproj = [Qt(1, 3), Qt(2, 3), 0d0]
      aproj = aproj/vlen(aproj)
      phi = dacos(aproj(1))
      Nmax = matrices%Nmaxs(ii)

      do i = 1, N_avgs
         theta = dble(i - 1)*2d0*pi/dble(N_avgs - 1)
         RR = matmul(transpose(R_aa(omega, theta)), transpose(R0))

         matrices%khat = matmul(RR, [0d0, 0d0, 1d0])
         matrices%khat = -matrices%khat/vlen(matrices%khat)
         k_sph = cart2sph(matrices%khat)

         matrices%R_fixk = transpose(rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0]))
         matrices%R = R_aa(matrices%khat, phi)

         call rot_setup()
         call incident_fields(E, a_in, b_in, ii)
         call scattered_fields(E, p, q, p90, q90, ii)
         if (allocated(S)) deallocate (S, K)
         call scattering_matrix_coeff(p, q, p90, q90, dcmplx(mesh%ki(ii)), size(points, 2), points, S)
         call extinction_matrix_coeff(p, q, p90, q90, dcmplx(mesh%ki(ii)), k_sph(2), k_sph(3), K)
         call cross_sections(p, q, a_in, b_in, dcmplx(mesh%ki(ii)), Nmax, Cext, Csca, Cabs)
         SS = SS + S/N_avgs
         KK = KK + K/N_avgs
         Csca_out = Csca_out + Csca/N_avgs
         Cext_out = Cext_out + Cext/N_avgs
      end do

   end subroutine mueller_align

end module mueller
