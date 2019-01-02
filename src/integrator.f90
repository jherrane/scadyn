module integrator
! Copyright (c) 2018 Joonas Herranen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use forces
   use shapebeam
   use setup 

   implicit none

contains

!****************************************************************************80
! Main routine
   subroutine integrate()
      integer :: t1, t2, rate

      write(*,'(A)') ' Integration in progress...'
      call system_clock(t1, rate)
      call initialize()
      if(int_mode /= 1) call solve_eoms()
      call system_clock(t2, rate)
      write(*,'(2(A,g0))') '  Done in ', real(T2 - T1)/real(rate), ' seconds'

      if(int_mode < 2) then
         write(*,'(A)') ' Spin-up and dissipation calculations in progress...'
         call system_clock(t1, rate)
         call solve_dissipation()
         call system_clock(t2, rate)
         write(*,'(2(A,g0))') '  Done in ', real(T2 - T1)/real(rate), ' seconds'

         write(*,'(A)') ' Alignment calculation in progress...'
         call system_clock(t1, rate)
         call RAT_alignment()
         call system_clock(t2, rate)
         write(*,'(2(A,g0))') '  Done in ', real(T2 - T1)/real(rate), ' seconds'
      end if

      if(shortlog == 1) call alignment_log()

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

      if (beam_shape == 1) call laguerre_gaussian_beams(p, l)
      if (beam_shape == 2) call bessel_beams()

   end subroutine initialize

!****************************************************************************80
! Calculates grain dynamics. All radiative torques are calculated in so
! called scattering frame, and all dynamics are integrated in principal frame.
   subroutine solve_eoms()
      integer :: i, wi, writeless
      real(dp) :: J(3), E, F(3), N(3)

      call start_log()

      writeless = 0
      wi = 1
      if(it_max>100000) writeless = 1
! Iteration of optical force calculation
      do i = 1, it_max
         if (i > it_stop) exit

         if(beam_shape == 0) then
            call vlv_update()
         else
            if(i==1) call ot_calibrate()
            call ot_update()
            if(autoterminate .AND. abs(matrices%x_CM(3)) > 4d0*2d0*pi/mesh%ki(1)) then
               stop 'Particle flew out of the tweezers'
            end if
         end if

         call update_values()
         J = matmul(matrices%I,matrices%w)
         E = 0.5d0*dot_product(J,matrices%w)
         if(dot_product(J,J)/=0d0) then
            matrices%q_param = 2d0*matrices%Ip(3)*E/dot_product(J,J)
         end if
         
         if(int_mode < 2) call check_stability(i)
         if(shortlog == 0) then
            if(writeless==0)then
               call append_log(i)
            else
               if(mod(i,(it_max/100000))==0)then
                  call append_log(wi)
                  wi = wi +1
               end if
            end if
         end if
         call print_bar(i+1, it_max)
         
      end do

   end subroutine solve_eoms

!****************************************************************************80

   subroutine check_stability(n)
      integer :: n,i
      ! Calculate mean q-parameters for the variance calculation
      if(n<=window) then
         matrices%q_mean = matrices%q_mean + matrices%q_param/window
         matrices%q_list(n) = matrices%q_param
      end if 

! If the simulation has run for long enough and the particle spins stably, 
! flag the calculation to end (by changing the it_stop value).
      if (n >= window) then
         if (n==window)then
            do i=1,window
               matrices%q_var = matrices%q_var + &
               ((matrices%q_list(i)-matrices%q_mean)**2)/window
            end do
         end if
         if (is_aligned == 0) then
            is_aligned = alignment_state()
            if(n==it_stop)then
               write(*,'(A)') "  Finished integration, using average q-parameter..."
               write(*,'(A,F11.5)') "   q = ", matrices%q_mean
               matrices%q_param0 = matrices%q_mean
            end if 
         else if (alignment_found == 0 .AND. int_mode < 2) then
            write(*,'(A)') "  Found nearly stable q-parameter closer to unity, stopping..."
            write(*,'(A,F11.5)') "   q = ", matrices%q_mean
            it_stop = n + it_log + 1
            alignment_found = 1
            matrices%q_param0 = matrices%q_mean
         end if
      end if
   end subroutine check_stability

!****************************************************************************80

   subroutine solve_dissipation()
      integer :: i, i_h
      real(dp) :: q, J(3), I3, t, h, JJ, dq_max, HH

      I3 = matrices%Ip(3)
      J = matmul(matrices%I,matrices%w)
      q = matrices%q_mean
      dq_max = 5d-3
      h = matrices%Ip(3)/matrices%Ip(1)
      i_h = int((q-1d0)/dq_max)
      JJ = vlen(J)

      write(*,'(A)') '  Calculating spin-up...'
      call spin_up(t, HH, JJ)
      write(*,'(A,F11.3,A)') '  Total spin-up time: ', t/365/24/3600, ' y'
      matrices%tau_rad = t/365/24/3600
      write(*,'(A,ES11.3)') '  Average spin-up torque: ', HH
      t = 0d0

      open(unit=1, file='out/q'//trim(matrices%out), action='write', status='replace')
      write(1,'(A)') 'q, J, t'
      if(q-1d0<1d-3) then 
         q = matrices%Ip(2)/matrices%Ip(1)
         i_h = int((q-1d0)/dq_max)
      end if
      write(1,'(3(ES11.3))') q, JJ, t
      do i = 0, i_h
         call relax_step(q, t, dq_max, JJ)
         if(q-1d0<1d-3) exit
         write(1,'(3(ES11.3))') q, JJ, t
      end do
      close(1)
      write(*,'(A,I0,A)') '  Total dissipation time: ', int(t/365/24/3600), ' y'
      write(*,'(A,2ES11.3)') '  Final (J, w): ', abs(JJ), abs(JJ)/matrices%Ip(3)
      matrices%tau_int = t/365/24/3600
      matrices%w_thermal = abs(JJ)/matrices%Ip(3)


   end subroutine solve_dissipation

!****************************************************************************80

   subroutine spin_up(t, H, Jout)
      real(dp), intent(out) :: t, H
      real(dp), intent(inout) :: Jout
      integer :: N_i, i
      real(dp) :: I3, dt, Q_t(3), J, wT, R0(3,3), xp(2), xi, phi, &
                  psi, F(3), N(3)

      I3 = matrices%Ip(3)
      R0 = matrices%R
      H = 0d0
      t = 0d0
      J = 0d0
      N_i = 0

      ! Iteration of optical force calculation
      do i = 0, window-1
         if (i>20 .AND. near_identity(matmul(transpose(R0),matrices%R),6d-2)) exit

         call vlv_update(0)

         call update_values()
         xp = get_xiphi()
         xi = xp(1)
         phi = xp(2)
         psi = matrices%B_psi

         call get_forces(1)
         Q_t = matmul(matrices%R,matmul(matrices%P,(matrices%N)))
         H = H + dot_product(Q_t, rhat(xi, phi, psi))
         N_i = i
      end do

      H = H/dble(N_i)
      wT = sqrt(2*k_b*matrices%Tgas/I3)
      dt = I3*wT/abs(H)
      J = J + H*dt
      t = t + dt
      Jout = J

   end subroutine spin_up

!****************************************************************************80
! Calculate alignment of stably spinning particle. 
   subroutine RAT_alignment()
      integer :: i, j, k, ind, Npoints, Nw, Nxi
      real(dp) :: w1, xi1
      real(dp), dimension(:, :, :), allocatable :: path_w, path_xi
      real(dp), dimension(:, :), allocatable :: FGH, xi, w

      call RAT_efficiency(60,20,FGH=FGH)
      allocate(matrices%FGH(size(FGH, 1), size(FGH,2)))
      matrices%FGH = FGH

! Start integration

      Nw = 10
      Nxi = 15
      Npoints = 300
      allocate(path_w(Nxi,Nw,Npoints), path_xi(Nxi,Nw,Npoints))
      call grid_xi_w(Nxi, Nw, xi, w, w_limits=[-40d0,40d0])

      open(unit=1,file="out/path"//trim(matrices%out), action="write", status="replace")
      write(1,'(I0)') Nxi
      write(1,'(I0)') Nw
      write(1,'(I0)') Npoints
      ind = 1
      do i = 1,Nxi
         do j = 1, Nw
            w1 = w(j,i)
            xi1 = xi(j,i)
            path_w(i,j,1) = w1
            path_xi(i,j,1) = xi1
            write(1,'(2ES12.3)') xi1, w1
            do k = 2, Npoints
               call ADE_update(w1,xi1)
               path_w(i,j,k) = w1
               path_xi(i,j,k) = xi1
               write(1,'(2ES12.3)') xi1, w1
            end do
            call print_bar(ind,Nw*Nxi)
            ind = ind + 1
         end do
      end do

      close(1)
      
   end subroutine RAT_alignment

!****************************************************************************80
! Calculates the torque efficiency of a particle averaged over rotation
! about its major axis of inertia (a_3) for angles 0 to pi between
! a_3 and incident k0. Results are automatically comparable with DDSCAT.
   subroutine RAT_efficiency(Nxi, Nphi, Npsi_in, FGH)
      integer :: i, j, k, Nxi, Nphi, Npsi, ind
      integer, optional :: Npsi_in
      real(dp) :: F, G, H
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

! Radiative torque calculations, emulating work of Draine & Weingartner (1997),
 ! ApJ 480:633  
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

               call get_forces(1)

               Q_t = matmul(transpose(matrices%R),matrices%Q_t)/dble(Nphi)
               F = dot_product(Q_t, xihat(xi(j), phi(k), psi(i)))
               G = dot_product(Q_t, phihat(xi(j), phi(k), psi(i)))
               H = dot_product(Q_t, rhat(xi(j), phi(k), psi(i)))
               F_coll(3, ind + 1) = F_coll(3, ind + 1) + F
               F_coll(4, ind + 1) = F_coll(4, ind + 1) + H
               F_coll(5, ind + 1) = F_coll(5, ind + 1) + G
            end do

            F_coll(1:2, ind + 1) = [xi(j), psi(i)]
            call print_bar(ind + 1, size(F_coll, 2))
            ind = ind + 1
         end do
      end do
      close(1)

      open(unit=1, file="out/F"//trim(matrices%out), ACTION="write", STATUS="replace")
      write (1, '(A)') 'xi   psi   F  H  G'
      do i = 1, size(F_coll, 2)
         write (1, '(6ES12.3)') dcos(F_coll(1:2, i)), F_coll(3:5, i)
      end do
      close(1)

      if(present(FGH)) then
         FGH(1,:) = F_coll(1,:)
         FGH(2,:) = F_coll(3,:)
         FGH(3,:) = F_coll(5,:)
         FGH(4,:) = F_coll(4,:)
      end if
   end subroutine RAT_efficiency

!****************************************************************************80
! Calculates adaptive time step over a maximum angle the particle
! can rotate during a step
   subroutine adaptive_step(N, Fin)
      real(dp), dimension(3), optional :: Fin
      real(dp) :: dw(3), N(3), F(3), dx(3), lambda

      dw = get_dotw(N, matrices%w, matrices%I, matrices%I_inv)
      matrices%dt = matrices%rot_max/max(vlen(matrices%w), vlen(dw))

      if(present(Fin))then
         F = Fin
         lambda = 2d0*pi/minval(mesh%ki)
         dx = matrices%v_CM*matrices%dt + 0.5d0*F*matrices%dt**2
         matrices%dt = matrices%dt*min(max(vlen(dx)/(lambda)/matrices%dt,0.5d0),2d0)
      end if

   end subroutine adaptive_step

!****************************************************************************80
! Calculates the update coefficients for a linear vector ODE
! via 2nd order Euler
! input:  v(3) = velocity at previous timestep
!         dt = timestep
! output: dx(3) = RK-coefficient for calculating
!         x_next = x + dt*v
   function euler_3D(v, dt) result(dx)

      real(dp) :: dx(3), v(3), dt

      dx = v*dt

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

      allocate(xis(size(matrices%FGH,2)), &
         F(size(matrices%FGH,2)), H(size(matrices%FGH,2)))
      xis = dcos(matrices%FGH(1,:))
      F = matrices%FGH(2,:)
      H = matrices%FGH(4,:)

      beta = matrices%Tdrag/matrices%TDG
! dxi = MF/w -betasin(xi)cos(xi)
      dxiw(1) = matrices%M*interp1D(F,xis,dcos(xi))/w - &
                  beta*dsin(xi)*dcos(xi)
! dw = MH -(1+betasin^2(xi))w
      dxiw(2) = matrices%M*interp1D(H,xis,dcos(xi)) - &
                  (1d0+beta*(dsin(xi))**2)*w

   end function get_dxiw

!****************************************************************************80
! Variational Lie-Verlet method for rotational integration
   subroutine vlv_update(mode)
      real(dp), dimension(3) :: F, N
      integer, optional :: mode
      integer :: maxiter, i1
      real(dp) :: Rn(3, 3), wnh(3), dt, I(3), Jw(3), Ff, FFD(3)
      real(dp) :: PxW(3), hel(3), Jac(3, 3), Jwn(3), iterstop

      maxiter = 50

      I = matrices%Ip
      Jac = 0.d0
      PxW = 0.d0
      Jw = 0.d0

! First force calculation for getting a reasonable value for dt
      if (matrices%tt == 0d0) then
         call get_forces()
      end if

      N = matrices%N
      call adaptive_step(N)
      dt = matrices%dt

      if (matrices%E < 1d-7) then 
         N = 0d0
         matrices%F = 0d0
      end if

      if(present(mode)) then
         if(mode==0) N = 0d0
      end if

      FFD = (dble(matrices%F))/mesh%mass
      matrices%vn = matrices%v_CM + euler_3D(FFD, dt)*dt
      matrices%xn = matrices%x_CM + euler_3D(matrices%vn, dt)*dt

! Step 1) Newton solve
      wnh = matrices%w ! Body angular velocity
      Jw = matmul(matrices%I, wnh)
      iterstop = maxval(matrices%I)*1.d-12

      do i1 = 1, maxiter
         Ff = dot_product(wnh, Jw)
         PxW = 0.5*dt*crossRR(Jw, wnh)

         hel = matmul(matrices%I, matrices%w) - Jw + PxW - &
         0.25d0*dt**2*Ff*wnh + 0.5d0*dt*N

         if (vlen(hel) < iterstop) then
            exit
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
      N = matrices%N
      if(present(mode)) then
         if(mode==0) N = 0d0
      end if

! Step 3) Explicit angular velocity update
      Jwn = Jw + PxW + 0.25d0*dt**2d0*dot_product(wnh, Jw)*wnh + 0.5d0*dt*N
      matrices%wn = matmul(matrices%I_inv, Jwn)

   end subroutine vlv_update

!****************************************************************************80
! Simpler integration for optical tweezers. Includes numerical shenanigans to
! take into account drag, gravity and some other more realistic effects. 
! Especially the rotational rates are heavily capped, as they would not be 
! expected.
   subroutine ot_update()
      real(dp), dimension(3) :: F, N, F_drag, N_drag, F_G, F_m, F_mag
      integer :: maxiter, i1
      real(dp) :: Rn(3, 3), wnh(3), dt, I(3), Jw(3), Ff
      real(dp) :: PxW(3), hel(3), Jac(3, 3), Jwn(3), iterstop
      real(dp) :: rvec(3), D, Re_w, Re

      maxiter = 50

      I = matrices%Ip
      Jac = 0.d0
      PxW = 0.d0
      Jw = 0.d0

      call get_forces()

! Calculate Reynolds numbers for rotation (_w) and translation
      Re_w = matrices%rho_med*vlen(matrices%w)*mesh%a**2/matrices%mu
      Re = matrices%rho_med*vlen(matrices%v_CM)*mesh%a/matrices%mu

! Use the Faxen's law for torque on a slowly spinning sphere when Reynolds 
! number is low. When Re_w>>1, drags are newtonian. 
      if(Re_w < 1d0)then
         N_drag = -8d0*pi*matrices%mu*mesh%a**3*matrices%w
      else
         N_drag = -0.5d0*matrices%rho_med*mesh%a**5*matrices%w*vlen(matrices%w)
      end if    

! Use Stokes drag (Re low) or Newton drag (Re large). Re is not probably 
! very large.
      if(Re < 1d0)then
         F_drag = -6d0*pi*matrices%mu*mesh%a*matrices%v_CM  
      else
         F_drag = -0.5d0*pi*mesh%a**2*matrices%rho_med*vlen(matrices%v_CM)*matrices%v_CM
      end if

! At really low pressures, the continuum assumption may not be valid.
! Disregard viscosity and thus drag.
      if(matrices%rho_med<1d-1) then
         F_drag = 0d0
         N_drag = 0d0
      end if

      F_G = -(mesh%rho-matrices%rho_med)*mesh%V*9.81d0*[0d0,0d0,1d0]
      F_m = 0.5d0*matrices%rho_med*mesh%V*dble(matrices%F)/mesh%mass
      F_mag = matrices%rho_med*mesh%V*crossRR(matrices%w,matrices%v_CM)*mesh%mass

      ! print*, vlen(N), vlen(N_drag), vlen(matrices%w)
      ! print*, vlen(matrices%F), vlen(F_drag), vlen(F_G), vlen(F_m), vlen(F_mag)
      N = matrices%N + N_drag
      F = (matrices%F + F_G + F_drag + F_m + F_mag)
      call adaptive_step( N, F )
      dt = matrices%dt

! For brownian motion, we approximate viscosity with water's
      if(brownian) then
         D = k_b*matrices%Tgas/(3d0*pi*1d-3*mesh%a)
         rvec = rand_vec()*sqrt(6d0*D*dt)
         matrices%x_CM = matrices%x_CM + rvec
      end if

! Step update taking gravity, drag, added mass and Magnus
! forces into account
      matrices%vn = matrices%v_CM + F/mesh%mass*dt
      matrices%xn = matrices%x_CM + matrices%v_CM*dt + 0.5d0*(F/mesh%mass)*dt**2 

! Rotation step 1) Newton solve
      wnh = matrices%w ! Body angular velocity
      Jw = matmul(matrices%I, wnh)
      iterstop = maxval(matrices%I)*1.d-12

      do i1 = 1, maxiter
         Ff = dot_product(wnh, Jw)
         PxW = 0.5*dt*crossRR(Jw, wnh)

         hel = matmul(matrices%I, matrices%w) - Jw + PxW - &
         0.25d0*dt**2*Ff*wnh + 0.5d0*dt*N

         if (vlen(hel) < iterstop) then
            exit
         end if

         Jac = -matrices%I + 0.5d0*dt*(reshape( &
              [0.d0, Jw(3) - I(1)*wnh(3), -Jw(2) + I(1)*wnh(2), &
              -Jw(3) + I(2)*wnh(3), 0.d0, Jw(1) - I(2)*wnh(1), &
               Jw(2) - I(3)*wnh(2), -Jw(1) + I(3)*wnh(1), 0.d0], [3, 3]) &
               - dt*real_outer_product(wnh, Jw)) - 0.25d0*dt*dt*Ff*eye(3)
         wnh = wnh - matmul(inv(Jac), hel)
         Jw = matmul(matrices%I, wnh)
      end do

! Rotation step 2) Explicit configuration update
      Rn = matmul(matrices%R, cay(dt*wnh))
      matrices%qn = mat2quat(Rn)
      matrices%qn = normalize_quat(matrices%qn)
      matrices%Rn = quat2mat(matrices%qn)

! Rotation step 3) Explicit angular velocity update
      Jwn = Jw + PxW + 0.25d0*dt**2d0*dot_product(wnh, Jw)*wnh &
      + 0.5d0*dt*N
      matrices%wn = matmul(matrices%I_inv, Jwn)
      matrices%F = F
      matrices%N = N

   end subroutine ot_update

!****************************************************************************
! Adjust E-field at intensity maximum (of origin-centered (incl.LG0l) beams)
   subroutine ot_calibrate
      real(dp), dimension(3) :: F
      real(dp) :: Fnorm, Fg
      logical :: tested_gravity

      tested_gravity = .FALSE.

      if(seedling /= 0) brownian = .TRUE.
      Fg = (mesh%rho)*mesh%V*9.81d0

      F = matrices%x_CM
      matrices%x_CM = [0d0, 0d0, 0d0]
      if(l>0 .AND. p == 0)then
         matrices%x_CM(1) = sqrt(l*beam_w0**2/2)/mesh%ki(1)
         write(*,'(A,ES9.3,A)') '  LG0l-maximum at x = ', matrices%x_CM(1)
      end if
      do while (.NOT. tested_gravity)
         call get_forces()
         Fnorm = vlen(matrices%F)
         if (abs(Fg)>1.01*abs(Fnorm) )  then
            matrices%E = sqrt(Fg/Fnorm)*matrices%E
         else
            write(*,'(A,ES9.3,A)') '  E-field maximum adjusted to ', matrices%E, ' V/m'
            write(*,'(A,ES9.3,A)') '  Intensity at maximum is then ', matrices%E**2/(2d0*377d0), ' W/m^2'
            if(l==0 .AND. p == 0. .AND. beam_shape == 1) then
               write(*,'(A,ES9.3,A)') '  Corresponding LG00 beam power ', &
               matrices%E**2/(2d0*377d0)*(pi/2d0)*(2d0/matrices%NA/mesh%ki(1))**2, ' W'
            end if
            call get_forces()
            tested_gravity = .TRUE.
         end if
      end do
      matrices%x_CM = F

   end subroutine ot_calibrate

!****************************************************************************80
! Runge-Kutta 4 update for the alignment differential equation (ADE) pair
   subroutine ADE_update(w, xi)
      real(dp) :: dxiw(2), w, xi, dt, k_w(4), k_xi(4)

      dt = 1d-2

      dxiw = get_dxiw(xi, w)
      k_w(1) = dxiw(2)
      k_xi(1) = dxiw(1)

      dxiw = get_dxiw(xi + k_xi(1)*dt/2, w + k_w(1)*dt/2)
      k_w(2) = dxiw(2)
      k_xi(2) = dxiw(1)

      dxiw = get_dxiw(xi + k_xi(2)*dt/2, w + k_w(2)*dt/2)
      k_w(3) = dxiw(2)
      k_xi(3) = dxiw(1)

      dxiw = get_dxiw(xi + k_xi(3)*dt/2, w + k_w(3)*dt/2)
      k_w(4) = dxiw(2)
      k_xi(4) = dxiw(1)

      w = w + dt*(k_w(1)+2*k_w(2)+2*k_w(3)+k_w(4))/6
      xi = xi + dt*(k_xi(1)+2*k_xi(2)+2*k_xi(3)+k_xi(4))/6

   end subroutine ADE_update

!****************************************************************************80
! Compute trajectory map
   subroutine grid_xi_w(Nxi, Nw, xi_grid, w_grid, w_limits)
! Arguments
      integer, intent(in) :: Nw, Nxi
      real(dp), optional :: w_limits(2)
      real(dp), dimension(:,:), allocatable, intent(out) :: w_grid, xi_grid
! Local variables
      real(dp), dimension(:), allocatable :: w, xi
      
      allocate(w(Nw), xi(Nxi))
      allocate(xi_grid(Nxi, Nw), w_grid(Nxi, Nw))

      if(.NOT. present(w_limits)) w_limits = [-100d0,100d0]
      call linspace(w_limits(1), w_limits(2), Nw, w)
      call linspace(-1d0,1d0,Nxi,xi)
      xi = acos(xi)

      call meshgrid(xi_grid,w_grid,xi,w)

   end subroutine grid_xi_w

!****************************************************************************80
end module integrator
