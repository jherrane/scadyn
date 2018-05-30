module integrator
   use forces
   use shapebeam

   implicit none

contains

!****************************************************************************80
! Main routine
   subroutine integrate()
      integer :: t1, t2, rate

      write(*,'(A)') ' Integration in progress...'
      call system_clock(t1, rate)
      call initialize()
      if(int_mode < 2) call solve_eoms()
      call system_clock(t2, rate)
      write(*,'(2(A,g0))') '  Done in ', real(T2 - T1)/real(rate), ' seconds'

   end subroutine integrate

!****************************************************************************
! Do the mandatory setup
   subroutine initialize()

! Calculates the mass parameters for the mesh
      if (use_mie == 1) then
         call mie_params()
      else
         call vie_params()
      end if

      call diagonalize_inertia()
      call interstellar_env()

      call polarization()
      call init_values()

      call allocate_inc_wave()

      if (beam_shape == 1) call gaussian_beams()
      if (beam_shape == 2) call laguerre_gaussian_beams(p, l)
      if (beam_shape == 3) call bessel_beams()
   end subroutine initialize

!****************************************************************************80
! Calculates grain dynamics. All radiative torques are calculated in so
! called scattering frame, and all ! dynamics are integrated in principal frame.
   subroutine solve_eoms()
      integer :: i
      character(len=80) :: lg

      lg = matrices%out
      call start_log(lg)

! Iteration of optical force calculation
      do i = 0, it_max-1
         if (i >= it_stop) exit

         select case (matrices%which_int)
         case (1)
            call RK4_update()
         case (2)
            call vlv_update()
         case (3)
            call PCDM_update()
         case default
            call euler_update()
         end select

         call update_values()
         call append_log(lg, i)
         call print_bar(i+1, it_max)
      end do

   end subroutine solve_eoms

!****************************************************************************80
! Somehow, if the particle movement is almost periodic, the averaged RAT should
! be calculated over one quasiperiod. Otherwise results will be more skewed.
   subroutine compute_log_RAT()
      integer :: i, j, k, Nang, ind, N_points, numlines, last
      integer :: t1, t2, rate
      real(dp) :: E, RR(3, 3)
      complex(dp), dimension(:), allocatable :: p, q, p90, q90
      real(dp), dimension(3, 3) :: R_B, R_xi, R_init, RP
      real(dp), dimension(3) :: k0, E0, E90, Q_t, nphi, a_3, x_B
      real(dp) :: xi, phi, psi, tol
      real(dp), dimension(:, :), allocatable :: F_coll
      real(dp), dimension(:), allocatable :: thetas

      tol = 5d-2
      numlines = it_log
      if(it_log == 0) numlines = it_max
      call read_log(numlines)

      matrices%R = eye(3)
      call rot_setup()

      k0 = matrices%khat
      E0 = real(matrices%E0hat)
      E90 = real(matrices%E90hat)

      Nang = 60 ! Angle of rotation of a_3 about e_1 (when psi=0)
      allocate (thetas(Nang))
      call linspace(0d0, pi, Nang, thetas)

      allocate (F_coll(6, Nang))
      F_coll(:, :) = 0d0
      ind = 0

      write (*, '(A)') '  Starting the calculation of averaged radiative torques:'
! First, set up B in the direction of psi (B still lies in the xz-plane)
      psi = matrices%B_psi
      x_B = [0d0, 1d0, 0d0]

! Rotation axis for precession averaging
      nphi = matmul(R_aa(x_B, psi), k0) ! Is actually B

      call system_clock(t1, rate)
      do i = 1, numlines
         if(near_identity(matmul(matrices%RRR(:, :, i),&
            transpose(matrices%RRR(:, :, 1))), tol ) .AND. i>100) then
            write(*, '(A,I0)') '  Periodicity detected at step ', i 
            last = i
            exit
         end if
         RR = transpose(matrices%RRR(:, :, i))
         RP = matmul(RR, transpose(matrices%P(:,:)))
! If stable rotation about a_3, this a_3 will be almost (anti/)paraller
! to w
         a_3 = RP(3,:) 
! Set up rotation from a_3 to new B
         R_B = rotate_a_to_b(a_3, nphi)

         matrices%khat = matmul(RR, [0d0, 0d0, 1d0])
         matrices%khat = -matrices%khat/vlen(matrices%khat)
         matrices%R_fixk = transpose(rotate_a_to_b(matrices%khat, [0.d0, 0.d0, 1.d0]))
         ind = 0
! Xi loop
         do j = 1, Nang
            xi = thetas(j)
! Rotation to xi of a_3 is a combined rotation of angle psi+xi
            R_xi = matmul(R_aa(x_B, xi), R_B)
            matrices%R = matmul(RR, R_xi) 
            RP = matmul(matrices%R, transpose(matrices%P(:,:)))
! Find the current orientation of a_3 and from there the current phi
            a_3 = RP(3,:)
            a_3 = [a_3(1), a_3(2), 0d0]
            a_3 = a_3/vlen(a_3)
            phi = dacos(a_3(1))

            call rot_setup()

            Q_t = 0d0
            if (matrices%whichbar == 0) then
               do k = 1, matrices%bars
                  call forcetorque(k)
                  Q_t = Q_t + matrices%Q_t/matrices%bars
               end do
            else
               k = matrices%whichbar
               call forcetorque(k)
               Q_t = Q_t + matrices%Q_t
            end if

            Q_t = matmul(matrices%R, Q_t)

            F_coll(3, ind + 1) = F_coll(3, ind + 1) + F_align(Q_t, xi, phi, psi)/Nang
            F_coll(4, ind + 1) = F_coll(4, ind + 1) + H_align(Q_t, xi, phi, psi)/Nang
            F_coll(5, ind + 1) = F_coll(5, ind + 1) + G_align(Q_t, phi, psi)/Nang
            F_coll(1:2, ind + 1) = [xi, psi]

            ind = ind + 1
         end do
         if(i==numlines) last = i
         call print_bar(i,numlines)
      end do

      F_coll(3:5, :) = F_coll(3:5, :)/last

      call system_clock(t2, rate)
      write(*,'(2(A,g0))') 'Done in ', real(T2 - T1)/real(rate), ' seconds'

      open (unit=1, file="out/F.out", ACTION="write", STATUS="replace")
      write (1, '(A)') 'xi   psi   F  H  G'
      do i = 1, size(F_coll, 2)
         write (1, '(6ES12.3)') dcos(F_coll(1:2, i)), F_coll(3:5, i)
      end do
      close (1)
      
   end subroutine compute_log_RAT

!****************************************************************************80
! Calculates adaptive time step over a maximum angle the particle
! can rotate during a step
   subroutine adaptive_step()
      real(dp) :: max_w, test_dt

! Use an adaptive time step, integrate over rotation <= rot_max
      max_w = vlen(matrices%w)
      if (max_w < 1d-7) max_w = max_w + vlen(get_dotw(matrices%N, matrices%w, matrices%I, matrices%I_inv))
      test_dt = matrices%rot_max/max_w

      if (test_dt < matrices%dt0) then
         matrices%dt = test_dt
      else
         matrices%dt = matrices%dt0
      end if

   end subroutine adaptive_step

!****************************************************************************80
! Calculates the update coefficients for a linear vector ODE
! via 4th order Runge-Kutta algorithm
! input:  v(3) = velocity at previous timestep
!         dt = timestep
! output: dx(3) = RK-coefficient for calculating
!         x_next = x + dt*dx/6
   function RK4_3D(v, dt) result(dx)

      real(dp) :: dx(3), v(3), a(3), b(3), c(3), d(3), dt

      a = v
      b = v + dt*a/2d0
      c = v + dt*b/2d0
      d = v + dt*c

      dx = (a + 2d0*b + 2d0*c + d)*dt/6d0

   end function RK4_3D

!****************************************************************************80
! Calculates the update coefficients for a linear vector ODE
! via 2nd order Euler
! input:  v(3) = velocity at previous timestep
!         dt = timestep
! output: dx(3) = RK-coefficient for calculating
!         x_next = x + dt*dx
   function euler_3D(v, dt) result(dx)

      real(dp) :: dx(3), v(3), dt

      dx = 0.5d0*v*dt

   end function euler_3D

!****************************************************************************80

   function get_dotq(w, q) result(dq)

      real(dp), dimension(4) :: q, dq, quat
      real(dp), dimension(3) :: w_W, w
      real(dp), dimension(3, 3) :: R

      R = quat2mat(q)
      w_W = matmul(R, w)

      quat(1) = 0.d0
      quat(2:4) = w_W

      dq = (1d0/2d0)*quat_mult(quat, q)

   end function get_dotq

!****************************************************************************80

   function get_dotw(N, w, I, I_inv) result(dw)
      real(dp) :: N(3), w(3), dw(3), I(3, 3), I_inv(3, 3)

      dw = matmul(I_inv, N + crossRR(matmul(I, w), w))

   end function get_dotw

!****************************************************************************80

   function get_dxiw(xi, w) result(dxiw)
      real(dp) :: dxiw(2), xi, w, beta
      real(dp), dimension(:), allocatable :: xis, F, H

      allocate(xis(size(matrices%FGH,2)), F(size(matrices%FGH,2)), H(size(matrices%FGH,2)))
      xis = dcos(matrices%FGH(1,:))
      F = matrices%FGH(2,:)
      H = matrices%FGH(4,:)

      beta = matrices%Tdrag/matrices%TDG
      dxiw(1) = matrices%M*interp1D(F,xis,dcos(xi)) - &
                  (1+beta*(dsin(xi))**2)*w
      dxiw(2) = matrices%M*interp1D(H,xis,dcos(xi)) - &
                  beta*dsin(xi)*dcos(xi)

   end function get_dxiw

!****************************************************************************80

   subroutine euler_update()
      real(dp) :: qn(4), wn(3), vn(3), xn(3), P(3, 3)
      complex(dp), dimension(6) :: FN

      P = matrices%P

      call get_forces()
      call adaptive_step()

      wn = .5d0*matrices%dt*get_dotw(matrices%N, matrices%w, matrices%I, matrices%I_inv)
      qn = matrices%dt*get_dotq(matrices%w, matrices%q)
      vn = .5d0*matrices%dt*matrices%F/mesh%mass
      xn = .5d0*matrices%dt*matrices%v_CM

      matrices%qn = matrices%q + qn
      matrices%wn = matrices%w + wn
      matrices%vn = matrices%v_CM + vn
      matrices%xn = matrices%x_CM + xn

      matrices%qn = normalize_quat(matrices%qn)
      matrices%Rn = quat2mat(matrices%qn) ! Update orientation matrix

   end subroutine euler_update

!****************************************************************************80
! Runge-Kutta 4 update for the rotational equations of motion
   subroutine RK4_update()
      real(dp), dimension(4) :: qn
      real(dp) :: dt, k_q(4, 4), k_w(3, 4), coeffs(4), &
                  k_v(3, 4), k_x(3, 4), R(3, 3), wn(3), P(3, 3)

      R = matrices%R
      P = matrices%P
      wn = matrices%w
      qn = matrices%q

      coeffs = [1d0, 2d0, 2d0, 1d0]

      call get_forces()

      call adaptive_step()

      dt = matrices%dt

      k_v(:, 1) = dt*matrices%F/mesh%mass
      k_x(:, 1) = dt*matrices%v_CM

      k_w(:, 1) = dt*get_dotw(matrices%N, matrices%w, matrices%I, matrices%I_inv)
      k_q(:, 1) = dt*get_dotq(matrices%w, matrices%q)
      matrices%R = quat2mat(matrices%q + 0.5d0*k_q(:, 1))

      call get_forces()
      k_v(:, 2) = dt*matrices%F/mesh%mass
      k_x(:, 2) = dt*(matrices%v_CM + k_v(:, 1)*0.5d0)

      k_w(:, 2) = dt*get_dotw(matrices%N, matrices%w + 0.5d0*k_w(:, 1), matrices%I, matrices%I_inv)
      k_q(:, 2) = dt*get_dotq(matrices%w + 0.5d0*k_w(:, 1), matrices%q + 0.5d0*k_q(:, 1))
      matrices%R = quat2mat(matrices%q + 0.5d0*k_q(:, 2))

      call get_forces()
      k_v(:, 3) = dt*matrices%F/mesh%mass
      k_x(:, 3) = dt*(matrices%v_CM + k_v(:, 2)*0.5d0)

      k_w(:, 3) = dt*get_dotw(matrices%N, matrices%w + 0.5d0*k_w(:, 2), matrices%I, matrices%I_inv)
      k_q(:, 3) = dt*get_dotq(matrices%w + 0.5d0*k_w(:, 2), matrices%q + 0.5d0*k_q(:, 2))
      matrices%R = quat2mat(matrices%q + 0.5d0*k_q(:, 3))

      call get_forces()
      k_v(:, 4) = dt*matrices%F/mesh%mass
      k_x(:, 4) = dt*(matrices%v_CM + k_v(:, 3))

      k_w(:, 4) = dt*get_dotw(matrices%N, matrices%w + k_w(:, 3), matrices%I, matrices%I_inv)
      k_q(:, 4) = dt*get_dotq(matrices%w + k_w(:, 3), matrices%q + k_q(:, 3))

      matrices%qn = matrices%q + (1d0/6d0)*matmul(k_q, coeffs)
      matrices%wn = matrices%w + (1d0/6d0)*matmul(k_w, coeffs)
      matrices%vn = matrices%v_CM + (1d0/6d0)*matmul(k_v, coeffs)
      matrices%xn = matrices%x_CM + (1d0/6d0)*matmul(k_x, coeffs)

      matrices%qn = normalize_quat(matrices%qn)
      matrices%Rn = quat2mat(matrices%qn) ! Update orientation matrix
   end subroutine RK4_update

!****************************************************************************80
! Variational Lie-Verlet method for rotational integration
   subroutine vlv_update()
      integer :: maxiter, i1
      real(dp) :: Rn(3, 3), wnh(3), dt, I(3), Jw(3), Ff
      real(dp) :: PxW(3), hel(3), Jac(3, 3), Jwn(3), iterstop
      complex(dp), dimension(6) :: FN

      maxiter = 50

      I = matrices%Ip
      Jac = 0.d0
      PxW = 0.d0
      Jw = 0.d0

! First force calculation for getting a reasonable value for dt
      if (matrices%tt == 0d0 .AND. .NOT. matrices%E < 1d-7) then
         call get_forces()
      end if
      call adaptive_step()
      dt = matrices%dt
      if (matrices%E < 1d-7) then 
         matrices%N = 0d0
         matrices%F = 0d0
      end if

! Step 1) Newton solve
      wnh = matrices%w ! Body angular velocity
      Jw = matmul(matrices%I, wnh)
      iterstop = maxval(matrices%I)*1.d-12

      do i1 = 1, maxiter
         Ff = dot_product(wnh, Jw)
         PxW = 0.5*dt*crossRR(Jw, wnh)

         hel = matmul(matrices%I, matrices%w) - Jw + PxW - 0.25d0*dt**2*Ff*wnh + 0.5d0*dt*matrices%N
         if (vlen(hel) < iterstop) then
            exit
            ! else if(i1==maxiter) then
            !    print*, ' Maximum number of VLV Newton iterations exceeded...'
         end if

         Jac = -matrices%I + 0.5d0*dt*(reshape( &
                                       [0.d0, Jw(3) - I(1)*wnh(3), -Jw(2) + I(1)*wnh(2), &
                                        -Jw(3) + I(2)*wnh(3), 0.d0, Jw(1) - I(2)*wnh(1), &
                                        Jw(2) - I(3)*wnh(2), -Jw(1) + I(3)*wnh(1), 0.d0], [3, 3]) &
                                       - dt*real_outer_product(wnh, Jw)) - 0.25d0*dt*dt*Ff*eye(3)
         wnh = wnh - matmul(inv(Jac), hel)
         Jw = matmul(matrices%I, wnh)
      end do

! Step 2) Explicit configuration update
      Rn = matmul(matrices%R, cay(dt*wnh))
      matrices%qn = mat2quat(Rn)
      matrices%qn = normalize_quat(matrices%qn)
      matrices%Rn = quat2mat(matrices%qn)

      matrices%R = matrices%Rn
      if (matrices%E > 1d-7) then 
         call get_forces()
      else
         matrices%F = 0d0
         matrices%N = 0d0
      end if

! Step 3) Explicit angular velocity update
      Jwn = Jw + PxW + 0.25d0*dt**2d0*dot_product(wnh, Jw)*wnh + 0.5d0*dt*matrices%N
      matrices%wn = matmul(matrices%I_inv, Jwn)

      matrices%vn = matrices%v_CM + euler_3D(dble(matrices%F)/mesh%mass, dt)
      matrices%xn = matrices%x_CM + euler_3D(matrices%v_CM, dt)

   end subroutine vlv_update

!****************************************************************************80
! Improved predictor-corrector direct multiplication (PCDM) method of
! Seelen, L. J. H., Padding, J. T., & Kuipers, J. A. M. (2016).
! Improved quaternion-based integration scheme for
! rigid body motion. Acta Mechanica, 1-9. DOI: 10.1007/s00707-016-1670-x
   subroutine PCDM_update()
      real(dp), dimension(4) :: q, qn
      real(dp) ::  w(3), wn(3), Imat(3, 3), Imat_inv(3, 3), dwb(3), wb(3), wbn(3), dt
      complex(dp), dimension(6) :: FN

      call get_forces()
      call adaptive_step()

      dt = matrices%dt

! Inertia tensor
      Imat = matrices%I
      Imat_inv = matrices%I_inv

! Extract old time step data
      wb = matrices%w
      q = matrices%q
      w = quat_rotation(q, wb)
      dwb = matmul(Imat_inv, matrices%N - crossRR(wb, matmul(Imat, wb)))

! Predictions for w at n + 3/4 time step and q at n+1 time step
      wbn = wb + 0.25d0*dwb*dt
      wn = quat_rotation(q, wbn) ! body->lab frame
      qn = quat_mult(q_worm(wn, dt), q)
      matrices%R = quat2mat(qn)

! Prediction for w at time n + 1
      wbn = wb + 0.5d0*dt*dwb
      wn = quat_rotation(qn, wbn) ! body->lab frame

! New torque with predicted values
      call get_forces()

! Corrected new values
      matrices%dw = matmul(Imat_inv, matrices%N - crossRR(wbn, matmul(Imat, wbn)))
      matrices%wn = wb + matrices%dw*dt
      matrices%qn = quat_mult(q_worm(wn, dt), q)
      matrices%Rn = quat2mat(matrices%qn)

      matrices%vn = matrices%v_CM + euler_3D(dble(matrices%F)/mesh%mass, dt)
      matrices%xn = matrices%x_CM + euler_3D(matrices%v_CM, dt)

   end subroutine PCDM_update

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
      matrices%wT = (15d0/(8d0*pi*mesh%alpha(3))*k_b*matrices%Tp/ &
                     (mesh%rho*mesh%a**5))**0.5
      matrices%TDG = (2d0*mesh%alpha(3)*mesh%rho*mesh%a**2)/ &
                     (5d0*matrices%Kw*matrices%B_len**2)/(4*pi/(mu))
      matrices%Tdrag = (pi*mesh%alpha(3)*mesh%rho*mesh%a)/(3*mesh%drag*matrices%nH* &
                        (2d0*pi*1.67d-27*k_b*matrices%Tp)**0.5)
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
! Runge-Kutta 4 update for the alignment differential equation (ADE) pair
   subroutine ADE_update(w, xi)
      real(dp) :: dxiw(2), w, xi, dt, coeffs(4), k_w(4), k_xi(4)

      dt = 1d-5

      dxiw = get_dxiw(xi, w)
      k_w(1) = dxiw(1)
      k_xi(1) = dxiw(2)

      dxiw = get_dxiw(xi + k_xi(1)*dt/2, w + k_w(1)*dt/2)
      k_w(2) = dxiw(1)
      k_xi(2) = dxiw(2)

      dxiw = get_dxiw(xi + k_xi(2)*dt/2, w + k_w(2)*dt/2)
      k_w(3) = dxiw(1)
      k_xi(3) = dxiw(2)

      dxiw = get_dxiw(xi + k_xi(3)*dt/2, w + k_w(3)*dt/2)
      k_w(4) = dxiw(1)
      k_xi(4) = dxiw(2)

      w = w + dt*(k_w(1)+2*k_w(2)+2*k_w(3)+k_w(4))/6
      xi = xi + dt*(k_xi(1)+2*k_xi(2)+2*k_xi(3)+k_xi(4))/6

   end subroutine ADE_update

!****************************************************************************80
! Compute trajectory map
   subroutine trajectory_xi_w(Nxi, Nw, xi_grid, w_grid, dxi, dw, w_limits)
! Arguments
      integer, intent(in) :: Nw, Nxi
      real(dp), optional :: w_limits(2)
      real(dp), dimension(:,:), allocatable, intent(out) :: dxi, dw, w_grid, xi_grid
! Local variables
      integer :: i, j
      real(dp), dimension(:), allocatable :: w, xi
      real(dp) :: dxiw(2)
      
      allocate(w(Nw), xi(Nxi))
      allocate(dxi(Nxi, Nw), dw(Nxi, Nw))
      allocate(xi_grid(Nxi, Nw), w_grid(Nxi, Nw))

      if(.NOT. present(w_limits)) w_limits = [-100d0,100d0]
      call linspace(w_limits(1), w_limits(2), Nw, w)
      call linspace(0d0,pi,Nxi,xi)

      call meshgrid(xi_grid,w_grid,xi,w)
      do i=1,Nxi
         do j = 1,Nw
            dxiw = get_dxiw(xi_grid(i,j), w_grid(i,j))
            dxi(i,j) = dxiw(1)
            dw(i,j) = dxiw(2)
         end do
      end do 

   end subroutine trajectory_xi_w

!****************************************************************************80

   subroutine mie_params()
      mesh%V = 4d0/3d0*pi*mesh%a**3
      mesh%mass = mesh%V*mesh%rho
      mesh%CM = [0d0, 0d0, 0d0]
      mesh%I = 2d0/5d0*mesh%mass*mesh%a**2*eye(3)

      matrices%CM = mesh%CM ! Save real CM if it happens to be not near origin
      matrices%x_CM = mesh%CM
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
         p0 = mesh%coord(:, mesh%etopol(1, i1))
         p1 = mesh%coord(:, mesh%etopol(2, i1))
         p2 = mesh%coord(:, mesh%etopol(3, i1))
         p3 = mesh%coord(:, mesh%etopol(4, i1))
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

      mesh%coord(1, :) = mesh%coord(1, :) - mesh%CM(1)
      mesh%coord(2, :) = mesh%coord(2, :) - mesh%CM(2)
      mesh%coord(3, :) = mesh%coord(3, :) - mesh%CM(3)
      matrices%CM = mesh%CM ! Save real CM if it happens to be not near origin
      mesh%CM = dble([0d0, 0d0, 0d0])
      matrices%x_CM = mesh%CM

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

!****************************************************************************80
end module integrator
