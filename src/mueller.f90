module mueller
! Copyright (c) 2018 Joonas Herranen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use shapebeam

   implicit none

contains

!****************************************************************************80
! Compute the mueller matrix according to the mode chosen.
   subroutine compute_mueller()
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
         call compute_ave_mueller(180, 1, ii, E, S)
      case ('ori')
         call compute_log_mueller(180, 90, ii, E, S)
      case ('perf_ori')
         if(matrices%xi_in>1d-6) then
            call compute_aligned_mueller(180, 1, ii, E, S, matrices%xi_in)
         else
            call compute_aligned_mueller(180, 1, ii, E, S)
         end if
      end select

      call write_mueller( S)
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
         matrices%R = rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0])

         S = S + update_mueller(N_theta, N_phi, ii, E, p, q, p90, q90)/N_points
         call print_bar(i, N_points)
      end do
   end subroutine compute_ave_mueller

!****************************************************************************80
! Compute the Mueller matrix of a perfectly aligned particle. The alignment
! direction is such that the major axis is either parallel to external B or
! precesses about it in angle xi
   subroutine compute_aligned_mueller(N_theta, N_phi, ii, E, S, xi_in)
      integer :: i, j, ii, NB, N_theta, N_phi, Nxi, ind
      real(dp) :: E, vec(3), phi, R0(3, 3), Qt(3, 3), omega(3), &
      aproj(3), RR(3, 3), theta, xi, B(3), Rxi(3,3), phi_B
      real(dp), dimension(:, :), allocatable :: S
      real(dp), optional :: xi_in
      complex(dp), dimension(:), allocatable :: p, q, p90, q90

      NB = 36 ! Number of points to calculate the perfect orientations

      allocate (S(N_theta*N_phi, 18))
      S = 0d0

      B = [sin(matrices%B_psi), 0d0, cos(matrices%B_psi)]
      Nxi = 20 ! Number of precession averaging directions
      if(.NOT. present(xi_in)) then
         xi = 0d0
         Nxi = 1
      else
         xi = xi_in
      end if 
      Rxi = R_aa([0d0,1d0,0d0], xi)
      ind = 0
      do j = 1,Nxi
         phi_B = dble(j - 1)*2d0*pi/Nxi
         R0 = matmul(Rxi,rotate_a_to_b(matrices%P(:, 3), B))
         R0 = matmul(R_aa(B,phi_B),R0)
         omega = matmul(R0, matrices%P(:,3))
         Qt = matmul(R0, matrices%P)
         aproj = [Qt(1, 3), Qt(2, 3), 0d0]
         aproj = aproj/vlen(aproj)
         phi = dacos(aproj(1))

         do i = 1, NB
            theta = dble(i - 1)*2d0*pi/NB
            RR = matmul(transpose(R_aa(omega, theta)), transpose(R0))

            matrices%khat = -matmul(RR, [0d0, 0d0, 1d0])
            matrices%khat = matrices%khat/vlen(matrices%khat)

            matrices%R = rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0])
            matrices%R = matmul(matrices%R,R_aa(matrices%khat, phi))

            S = S + update_mueller(N_theta, N_phi, ii, E, p, q, p90, q90)/(NB*Nxi)
            ind = ind +1 
            call print_bar(ind, Nxi*NB)
         end do
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
         matrices%R = transpose(rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0]))

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

end module mueller
