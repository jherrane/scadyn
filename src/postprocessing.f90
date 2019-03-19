module postprocessing
! Copyright (c) 2018 Joonas Herranen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use integrator
   use shapebeam
   use mueller

   implicit none

contains

!****************************************************************************80
! Test the  methods of calculating scattering forces. When incident field is 
! a planewave, the results should match always with the example. Numerical
! shenanigans seem to occur with other beams. PSA: Fixing it can be a hassle.
   subroutine test_methods()
      integer :: i, j, ii, range(2)
      real(dp), dimension(:, :), allocatable :: Q_fcoll, Q_tcoll
      real(dp), dimension(3) :: F, N, FF, NN, N_B, N_DG
      real(dp) :: urad

      call rot_setup()

      allocate (Q_fcoll(3, matrices%bars))
      allocate (Q_tcoll(3, matrices%bars))

      if (matrices%whichbar == 0) then
         range(1) = 1
         range(2) = matrices%bars
      else
         range(1) = matrices%whichbar
         range(2) = matrices%whichbar
         write (*, '(A, 20F6.3)') ' Chosen wavelength: ', 2d6*pi/mesh%ki(matrices%whichbar)
      end if

      N = dcmplx(0.0d0, 0.0d0)
      F = dcmplx(0.0d0, 0.0d0)
      Q_tcoll = 0d0
      Q_fcoll = 0d0

      do i = range(1), range(2)
         call forcetorque(i)
         ii = i
         if (matrices%whichbar > 0) ii = matrices%whichbar
         F = F + matrices%force
         N = N + matrices%torque
         urad = epsilon*(matrices%E_rel(ii)*matrices%E)**2/2d0
         Q_fcoll(:, ii) = matmul(matrices%R, &
            matrices%force/urad/pi/mesh%a**2)
         Q_tcoll(:, ii) = matmul(matrices%R, &
            matrices%torque)*(mesh%ki(ii))/urad/pi/mesh%a**2
      end do

      NN = matmul(matrices%R, N)
      FF = matmul(matrices%R, F)

      print *, ''
      write (*, '(A)') 'Analytical VSWF coefficient integration of MST:'

      write (*, '(A,3ES11.3,A)') ' F   = (', real(FF), ' ) N'
      write (*, '(A,3ES11.3,A)') ' N   = (', real(NN), ' ) Nm'
      write (*, '(A,3ES11.3,A)') ' a_F = (', real(FF)/mesh%mass, ' ) m/s^2'
      write (*, '(A,3ES11.3,A)') ' a_N = (', matmul(matrices%I_inv,real(NN)), ' ) rad/s^2'

      if (j == 2) then
         N_DG = DG_torque()
         N_B = barnett_torque()

         write (*, '(A,3ES11.3,A)') ' N_DG = (', real(N_DG), ' ) Nm'
         write (*, '(A,3ES11.3,A)') ' N_B = (', real(N_B), ' ) Nm'
      end if

      print *, ''

      print *, 'Efficiencies for each wavelength separately'
      do i = 1, size(Q_fcoll, 2)
         Q_fcoll(:, i) = matmul(matrices%R, Q_fcoll(:, i))
         Q_tcoll(:, i) = matmul(matrices%R, Q_tcoll(:, i))
      end do

      write (*, '(A, 80F7.1)') '  WL(nm):', 2d9*pi/mesh%ki
      write (*, '(A, 80F7.3)') '   Qfx:', Q_fcoll(1, :)
      write (*, '(A, 80F7.3)') '   Qfy:', Q_fcoll(2, :)
      write (*, '(A, 80F7.3)') '   Qfz:', Q_fcoll(3, :)

      write (*, '(A, 80F7.3)') '   Qtx:', Q_tcoll(1, :)
      write (*, '(A, 80F7.3)') '   Qty:', Q_tcoll(2, :)
      write (*, '(A, 80F7.3)') '   Qtz:', Q_tcoll(3, :)

   end subroutine test_methods

!****************************************************************************80
! Calculates the torque efficiency of a particle averaged over rotation
! about its major axis of inertia (a_3) for angles 0 to pi between
! a_3 and incident k0. Results are automatically comparable with DDSCAT.
   subroutine torque_efficiency(q_factor)
      real(dp), intent(out), optional :: q_factor
      integer :: i, Ntheta
      real(dp), dimension(3) :: Q_t
      real(dp), dimension(:), allocatable :: theta, costheta
      real(dp), dimension(:, :), allocatable :: Q_coll

      Ntheta = 60 ! Angle of rotation of a_3 about e_1 (when psi=0)

      allocate (theta(Ntheta), costheta(Ntheta), Q_coll(3, Ntheta))
      call linspace(-1d0,1d0, Ntheta, costheta)
      theta = acos(costheta)

! Torque efficiency calculations as in Lazarian2007b
      write (*, '(A)') '  Starting the calculation of beta-averaged torque efficiency:'
! Theta loop
      do i = 1, Ntheta
         Q_t = get_Qav(theta(i), 0d0)
! Flip the coordinate labels to match Lazarian2007b: In this code, k is along 
! z instead of x, 0-polarization is the y of L2007b. So, inverse that next.
         Q_coll(:, i) = [Q_t(3),Q_t(1),Q_t(2)]
         call print_bar(i, Ntheta)
      end do
      open (unit=1, file="out/Q"//trim(matrices%out), ACTION="write", STATUS="replace")
      if(present(q_factor)) then
         q_factor = maxval(abs(Q_coll(1,:)))/maxval(abs(Q_coll(2,:)))
         write (1, '(1E12.3)') q_factor
      else
         write (1, '(A)') 'cos(theta)  Q_{t,1} Q_{t,2} Q_{t,3}'
         do i = 1, Ntheta
            write (1, '(4E12.3)') dcos(theta(i)), Q_coll(:, i)
         end do
      end if 
      close (1)
   end subroutine torque_efficiency

!****************************************************************************80
! Calculate force efficiency as a function of displacement in optical tweezers. 
! The particle is considered to be in a constant orientation.
   subroutine force_map()
      integer :: i, j, k, Nxyz, ind, axis, Nav, iii
      real(dp) :: limit
      real(dp), dimension(3,3) :: R, R_beta
      real(dp), dimension(3) :: Q, a
      real(dp), dimension(:), allocatable :: x, beta
      real(dp), dimension(:, :), allocatable :: Q_coll
      character(1) :: out1

      Nxyz = 100
      Nav = 180
      limit = 4d0*pi/mesh%ki(1)
      
      allocate (x(Nxyz), Q_coll(3, Nxyz**3))
      allocate(beta(Nav))
      call linspace(0d0, 2d0*pi, Nav, beta)
      call linspace(-limit, limit, Nxyz, x)

      do iii = 1,3
         ind = 1
         axis = iii
         a = matrices%P(1:3, axis)
         write(out1,'(I1)') axis
         R = rotate_a_to_b(a,[0d0,0d0,1d0])

         Q_coll(:,:) = 0d0
         write (*, '(A)') '  Calculating force efficiency '//out1//' around beam focus:'
         do j = 1,3
            do i = 1, Nxyz
               matrices%x_CM = [0d0, 0d0, 0d0]
               matrices%x_CM(j) = x(i)
               do k = 1,Nav
                  R_beta = R_aa([0d0,0d0,1d0], beta(k))

   ! The ultimate rotation matrices for scattering event
                  matrices%R = matmul(R_beta, R)

                  call get_forces()
                  Q = matmul(matrices%R,matrices%Q_f)
                  Q_coll(:,ind) = Q_coll(:,ind) + Q
               end do
               call print_bar(ind, 3*Nxyz)
               ind = ind + 1 
            end do 
         end do
         Q_coll = Q_coll/Nav

         
         open (unit=1, file="out/Qf"//out1//trim(matrices%out), ACTION="write", STATUS="replace")
         ind = 1
         write (1, '(A)') 'x  y  z  Q_{f,1} Q_{f,2} Q_{f,3}'
         do j = 1,3
            do i = 1, Nxyz
               matrices%x_CM = [0d0, 0d0, 0d0]
               matrices%x_CM(j) = x(i)
               write (1, '(6E12.3)') matrices%x_CM, Q_coll(:, ind)
               ind = ind + 1
            end do 
         end do

         close (1)
      end do
   end subroutine force_map

end module postprocessing
