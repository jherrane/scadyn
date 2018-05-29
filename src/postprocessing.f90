module postprocessing
   use integrator
   use shapebeam
   use mueller

   implicit none

contains

!****************************************************************************80
! Test the numerical and analytical methods of calculating scattering forces. 
! When incident field is planewave, the results should match always. Numerical
! shenanigans seem to occur with other beams. PSA: Fixing it can be a hassle.
   subroutine test_methods()
      integer :: i, j, ii, range(2)
      real(dp), dimension(:, :), allocatable :: Q_fcoll, Q_tcoll
      complex(dp), dimension(3) :: F, N, FF, NN, N_B, N_DG

      matrices%R_fixk = transpose(rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0]))

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

      do j = 1, 2
         N = dcmplx(0.0d0, 0.0d0)
         F = dcmplx(0.0d0, 0.0d0)
         Q_tcoll = 0d0
         Q_fcoll = 0d0

         do i = range(1), range(2)
            select case (j)
            case (1); call forcetorque_num(i)
            case (2); call forcetorque(i)
            end select
            ii = i
            if (matrices%whichbar > 0) ii = matrices%whichbar
            F = F + matrices%force
            N = N + matrices%torque
            Q_fcoll(:, ii) = matmul(matrices%R, matrices%Q_f)
            Q_tcoll(:, ii) = matmul(matrices%R, matrices%Q_t)
         end do

         NN = matmul(matrices%R, N)
         FF = matmul(matrices%R, F)

         print *, ''
         select case (j)
         case (1); write (*, '(A)') 'Numerical integration of MST:'
         case (2); write (*, '(A)') 'Analytical VSWF coefficient integration of MST:'
         end select
         write (*, '(A,3ES11.3,A)') ' F = (', real(FF), ' ) N'
         write (*, '(A,3ES11.3,A)') ' N = (', real(NN), ' ) Nm'

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

         write (*, '(A, 80F7.1)') 'WL(nm):', 2d9*pi/mesh%ki
         write (*, '(A, 80F7.3)') '   Qfx:', Q_fcoll(1, :)
         write (*, '(A, 80F7.3)') '   Qfy:', Q_fcoll(2, :)
         write (*, '(A, 80F7.3)') '   Qfz:', Q_fcoll(3, :)

         write (*, '(A, 80F7.3)') '   Qtx:', Q_tcoll(1, :)
         write (*, '(A, 80F7.3)') '   Qty:', Q_tcoll(2, :)
         write (*, '(A, 80F7.3)') '   Qtz:', Q_tcoll(3, :)
      end do

   end subroutine test_methods

!****************************************************************************80
! Calculates the torque efficiency of a particle averaged over rotation
! about its major axis of inertia (a_3) for angles 0 to pi between
! a_3 and incident k0. Results are automatically comparable with DDSCAT.
   subroutine torque_efficiency()
      integer :: i, Ntheta
      real(dp), dimension(3) :: Q_t
      real(dp), dimension(:), allocatable :: theta
      real(dp), dimension(:, :), allocatable :: Q_coll

      Ntheta = 60 ! Angle of rotation of a_3 about e_1 (when psi=0)

      allocate (theta(Ntheta), Q_coll(3, Ntheta))
      call linspace(0d0, pi, Ntheta, theta)

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

      open (unit=1, file="out/Q.out", ACTION="write", STATUS="replace")
      write (1, '(A)') 'cos(theta)  Q_{t,1} Q_{t,2} Q_{t,3}'
      do i = 1, Ntheta
         write (1, '(4E12.3)') dcos(theta(i)), Q_coll(:, i)
      end do
      close (1)
   end subroutine torque_efficiency

!****************************************************************************80
! Calculates the torque efficiency of a particle averaged over rotation
! about its major axis of inertia (a_3) for angles 0 to pi between
! a_3 and incident k0. Results are automatically comparable with DDSCAT.
   subroutine RAT_efficiency(Nxi, Nphi, Npsi_in, FGH)
      integer :: i, j, k, l, Nxi, Nphi, Npsi, ind
      integer, optional :: Npsi_in
      real(dp), dimension(3) :: Q_t, n_phi
      real(dp), dimension(:, :), allocatable :: F_coll
      real(dp), dimension(:, :), allocatable, optional, intent(out) :: FGH
      real(dp), dimension(3,3) :: R_B, R_phi, R_xi
      real(dp), dimension(:), allocatable :: psi, xi, phi

      if(.NOT. present(Npsi_in)) then
         Npsi = 1
      else
         Npsi = Npsi_in
      end if

      allocate (xi(Nxi), phi(Nphi), psi(Npsi), F_coll(6, Nxi*Npsi))
      if(present(FGH)) allocate(FGH(4,Nxi*Npsi))
      call linspace(0d0, pi, Nxi, xi)
      call linspace(0d0, pi/2d0, Npsi, psi)
      if(.NOT. present(Npsi_in)) psi(1) = matrices%B_psi
      call linspace(0d0, pi*2d0, Nphi, phi)
      F_coll(:, :) = 0d0

! Radiative torque calculations, emulating work of Draine & Weingartner (1997), ApJ 480:633  
      ind = 0
      write (*, '(A)') '  Starting the calculation of phi-averaged radiative torques:'

! First, set up rotation axis for precession averaging in the direction of psi (B)     
      do i = 1, Npsi
         n_phi = matmul(R_aa([0d0, 1d0, 0d0], psi(i)), [0d0,0d0,1d0]) ! Is actually B
         R_B = rotate_a_to_b(matrices%P(:, 3), n_phi)

! Rotation to xi of a_3 is a combined rotation of angle psi+xi
         do j = 1, Nxi
            R_xi = matmul(R_aa([0d0, 1d0, 0d0], xi(j)), R_B)

! Find rotation of angle phi around the rotated a_3-axis and the whole rotation
            do k = 1, Nphi
               R_phi = R_aa(n_phi, phi(k))
               matrices%R = matmul(R_phi, R_xi) ! First a_3 to xi, then beta about a_3
               
               call get_forces()

               Q_t = matmul(matrices%R,matrices%Q_t)/Nphi

               F_coll(3, ind + 1) = F_coll(3, ind + 1) + F_align(Q_t, xi(j), phi(k), psi(i))
               F_coll(4, ind + 1) = F_coll(4, ind + 1) + H_align(Q_t, xi(j), phi(k), psi(i))
               F_coll(5, ind + 1) = F_coll(5, ind + 1) + G_align(Q_t, phi(k), psi(i))
            end do

            F_coll(1:2, ind + 1) = [xi(j), psi(i)]
            call print_bar(ind + 1, size(F_coll, 2))
            ind = ind + 1
         end do
      end do

      open (unit=1, file="out/F.out", ACTION="write", STATUS="replace")
      write (1, '(A)') 'xi   psi   F  H  G'
      do i = 1, size(F_coll, 2)
         write (1, '(6ES12.3)') dcos(F_coll(1:2, i)), F_coll(3:5, i)
      end do
      close (1)

      if(present(FGH)) then
         FGH(1,:) = F_coll(1,:)
         FGH(2,:) = F_coll(3,:)
         FGH(3,:) = F_coll(5,:)
         FGH(4,:) = F_coll(4,:)
      end if
   end subroutine RAT_efficiency

!****************************************************************************80
! Calculate the torque efficiency projected on the inertia axes for all 
! possible direction of radiation incidence. Similar analysis can be done
! taking into account the particle orientations or spin state. 
   subroutine stability_analysis(result)
      integer :: i, j, k, ind, Npoints, Nphi, Ntheta, Nbeta
      real(dp), dimension(3, 3) ::  Rbeta, Rtheta
      real(dp), dimension(3) :: n_beta, k0, a_1, a_2, a_3, N, maxdir, a
      real(dp) :: mx, tmp
      real(dp), dimension(:, :), allocatable :: Q_coll, vec 
      real(dp), dimension(:), allocatable :: theta, costheta, phi, beta
      real(dp), dimension(3), optional :: result
      
      call rot_setup()

      k0 = matrices%khat

      a_1 = matrices%P(1:3, 1)
      a_2 = matrices%P(1:3, 2)
      a_3 = matrices%P(1:3, 3)

      Ntheta = 90
      Nphi = 180
      allocate (theta(Ntheta), costheta(Ntheta), phi(Nphi))
      call linspace(0d0, 2d0*pi, Nphi, phi)
      call linspace(1d0, -1d0, Ntheta, costheta)
      do i = 1, Ntheta
         theta(i) = dacos(costheta(i))
      end do
      Npoints = Ntheta*Nphi
      allocate (vec(3, Npoints))
      vec = 0d0
      mx = 0d0

      vec = uniform_sphere(Npoints, theta, phi)

      open (unit=1, file="out/NdotQ.dat", ACTION="write", STATUS="replace")
      write (1, '(2(A,I0))') 'Ntheta = ', Ntheta, ' | Nphi = ', Nphi
      write (1, '(A)') '  x   y   z   N.a_1   N.a_2   N.a_3'

      do i = 1, Npoints
         ! The ultimate rotation matrices for scattering event
         matrices%R = rotate_a_to_b(k0, vec(:, i))
         call get_forces()

         a = matmul(transpose(matrices%R), a_3)
         N = matmul(transpose(matrices%P), matrices%Q_t)
         tmp = dot_product(N, a_3)
         if (tmp > mx) then
            mx = tmp
            maxdir = a
         end if
         write (1, '(6E12.3)') a(1), a(2), a(3), &
                dot_product(N, a_1), dot_product(N, a_2), tmp
         call print_bar(i + 1, Npoints)
      end do
      close (1)

      tmp = dacos(dot_product(k0, maxdir/vlen(maxdir)))*180d0/pi
      write (*, '(A,3F7.3,A,F7.3)') 'Stablest direction = (', maxdir, &
                                    ' ), angle between a_3 and k = ', tmp
      if(present(result)) result = maxdir
   end subroutine stability_analysis

!****************************************************************************80
! Mueller matrices needed in SOCpol.
   subroutine test_mueller(Ntheta, Nphi)
      integer :: i, j, ind, halton_init, Nphi, Ntheta
      real(dp) :: al_direction(3)
      real(dp), dimension(:), allocatable :: a_dist
      real(dp), dimension(:,:), allocatable :: points

      halton_init = 0
      allocate(points(2,Ntheta*Nphi))

      ind = 1
      do i = 1, Ntheta
         do j = 1, Nphi
            points(1, ind) = pi*(i - 1)/(Ntheta) + pi/Ntheta/2.0
            points(2, ind) = 2*pi*(j - 1)/Nphi
            ind = ind + 1
         end do
      end do

      a_dist = mesh%ki*mesh%a/mesh%ki(2)
      al_direction = [0d0, 0d0, 1d0]
      call scattering_extinction_matrices(a_dist, points, al_direction)

   end subroutine test_mueller

!****************************************************************************80
! Calculate alignment of stably spinning particle. The stability direction is
! determined by force analysis
   subroutine stable_particle_RAT()
      integer :: i, j, k, Nang, ind, N_points, numlines, last
      integer :: t1, t2, rate
      real(dp) :: E, RR(3, 3)
      complex(dp), dimension(:), allocatable :: p, q, p90, q90
      real(dp), dimension(3, 3) :: R_B, R_xi, R_init, RP
      real(dp), dimension(3) :: k0, E0, E90, Q_t, nphi, a_3, x_B, xstable
      real(dp) :: xi, w, phi, psi, tol
      real(dp), dimension(:, :), allocatable :: FGH
      real(dp), dimension(:), allocatable :: thetas

!       call stability_analysis(xstable)

! ! Rotation like this is not strictly needed for the next steps, but P(:,3) 
! ! gives now the internal alignment direction w.r.t. the incidence direction.
!       RP = rotate_a_to_b(matrices%P(:,3), xstable)
!       matrices%P = matmul(RP, matrices%P)
      call RAT_efficiency(60,20,FGH=FGH)
      allocate(matrices%FGH(size(FGH, 1), size(FGH,2)))
      matrices%FGH = FGH

! Start integration

      xi = pi/2d0
      w = 1d0
      ! print*, w, cos(xi)
      do i = 1,3000
         call ADE_update(w, xi)
         ! print*, w, cos(xi)
      end do
      ! open (unit=1, file="out/F.out", ACTION="write", STATUS="replace")
      ! write (1, '(A)') 'xi   psi   F  H  G'
      ! do i = 1, size(F_coll, 2)
      !    write (1, '(6ES12.3)') dcos(F_coll(1:2, i)), F_coll(3:5, i)
      ! end do
      ! close (1)
      
   end subroutine stable_particle_RAT

end module postprocessing
