module forces
! Copyright (c) 2018 Joonas Herranen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use T_matrix

   implicit none

contains

!****************************************************************************80
! Calculate forces and torques during single step
   subroutine get_forces(mode)
      real(dp), dimension(3) :: Q_t, Q_f
      real(dp) :: u, wl
      real(dp), dimension(3) :: F, N, N_DG, N_B
      integer :: i
      integer, optional :: mode

      call rot_setup()

      Q_t = 0d0
      Q_f = 0d0
      N = 0d0
      F = 0d0
      if(present(mode)) then
         u = 0d0
         wl = 0d0
         do i = 1, matrices%bars
            u = u + epsilon*(matrices%E*matrices%E_rel(i))**2
            wl = wl + 2d0*pi/mesh%ki(i)
         end do
         wl = wl/matrices%bars
      end if

! Iteration of a single wavelength, if multiple wavelengths are calculated, scale
! the efficiency sum to take account the relative intensities i.e. omit E_rel(i)
      if (matrices%whichbar == 0) then
         do i = 1, matrices%bars
            call forcetorque(i)
            F = F + matrices%force
            N = N + matrices%torque
            if(present(mode))then
               Q_t = Q_t + matrices%torque/mesh%a**2/(wl*u)/matrices%bars
            else
               Q_t = Q_t + matrices%torque/&
               (epsilon*(matrices%E)**2)/&
                  (0.5d0*pi/mesh%ki(i))/mesh%a**2/matrices%bars
            end if
         end do
      else
         call forcetorque(matrices%whichbar)
         F = F + matrices%force
         N = N + matrices%torque
         Q_t = Q_t + matrices%torque/(0.5d0*pi/mesh%ki(matrices%whichbar))/&
               (epsilon*(matrices%E_rel(matrices%whichbar)*matrices%E)**2)/&
               mesh%a**2/matrices%bars
         Q_f = Q_f + matrices%force/&
               (epsilon*(matrices%E_rel(matrices%whichbar)*matrices%E)**2/2d0)/&
               pi/mesh%a**2/matrices%bars
      end if

      if (calc_extra_torques == 1) then
         N_DG = DG_torque()
         N_B = barnett_torque()
         N = N + N_B + N_DG
      end if

      matrices%F = matmul(matrices%R, F) ! {x}_sca -> {x}_lab
      matrices%N = matmul(transpose(matrices%P), N) ! {N}_sca -> {N}_b
      matrices%Q_t = Q_t
      matrices%Q_f = Q_f

   end subroutine get_forces

!****************************************************************************80
! Calculate forces and torques using the analytical z-formulae
   subroutine forcetorque(i)
      complex(dp), dimension(:), allocatable :: a_in, b_in, a90, b90, &
                                                a, b, a_temp, b_temp, &
                                                p, q, p90, q90, &
                                                a2, b2, p2, q2, a290, &
                                                b290, p290, q290
      complex(dp), dimension(:, :), allocatable :: Taa, Tab, Tba, Tbb
      complex(8), dimension(:), allocatable :: rotD, rotD90
      integer, dimension(:, :), allocatable :: indD, indD90
      integer :: Nmax, i, las, nm
      real(dp) :: tx, ty, tz, fx, fy, fz, E

      E = matrices%E_rel(i)*matrices%E
      mesh%k = mesh%ki(i)

      Nmax = matrices%Nmaxs(i)
      las = (Nmax + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1
      nm = (Nmax + 1)**2 - 1

      allocate (a(nm), b(nm), a_in(nm), b_in(nm), a_temp(nm), b_temp(nm), &
         a90(nm), b90(nm), p(nm), q(nm), p90(nm), q90(nm), &
         a2(nm), b2(nm), a290(nm), b290(nm), p2(nm), q2(nm), p290(nm), q290(nm))

      rotD = matrices%rotDs(1:las, i)
      indD = matrices%indDs(1:las, :, i)
      rotD90 = matrices%rotD90s(1:las, i)
      indD90 = matrices%indD90s(1:las, :, i)

      Taa = matrices%Taai(1:nm, 1:nm, i)
      Tab = matrices%Tabi(1:nm, 1:nm, i)
      Tba = matrices%Tbai(1:nm, 1:nm, i)
      Tbb = matrices%Tbbi(1:nm, 1:nm, i)

      a_temp = E*matrices%as(1:nm, i)/sqrt(2d0*sqrt(mu/epsilon)*mesh%k**2)/2d0
      b_temp = E*matrices%bs(1:nm, i)/sqrt(2d0*sqrt(mu/epsilon)*mesh%k**2)/2d0

      if(beam_shape /= 0 .AND. vlen(matrices%x_CM) > 1d-8) then
         call translate(matrices%x_CM, Nmax, Nmax, dcmplx(mesh%k), &
            a_temp, b_temp, a_in, b_in, 0)
      else
         a_in = a_temp
         b_in = b_temp
      end if

      a = sparse_matmul(rotD, indD, a_in, nm)
      b = sparse_matmul(rotD, indD, b_in, nm)

      a90 = sparse_matmul(rotD90, indD90, a_in, nm)
      b90 = sparse_matmul(rotD90, indD90, b_in, nm)

      if (matrices%polarization == 1 .OR. beam_shape /= 0) then
         a90 = dcmplx(0d0)
         b90 = dcmplx(0d0)
      end if

      p = matmul(Taa, a) + matmul(Tab, b)
      q = matmul(Tbb, b) + matmul(Tba, a)
      p90 = matmul(Taa, a90) + matmul(Tab, b90)
      q90 = matmul(Tbb, b90) + matmul(Tba, a90)

      p = 2d0*p + a
      q = 2d0*q + b
      p90 = 2d0*p90+a90
      q90 = 2d0*q90+b90

! formula is for z-direction
      fz = F_z(Nmax, a, b, p, q) + F_z(Nmax, a90, b90, p90, q90)
      tz = T_z(Nmax, a, b, p, q) + T_z(Nmax, a90, b90, p90, q90)

! x-direction
      a2 = sparse_matmul(matrices%rotXs(1:las, i), matrices%indXs(1:las, :, i), a, nm)
      b2 = sparse_matmul(matrices%rotXs(1:las, i), matrices%indXs(1:las, :, i), b, nm)
      p2 = sparse_matmul(matrices%rotXs(1:las, i), matrices%indXs(1:las, :, i), p, nm)
      q2 = sparse_matmul(matrices%rotXs(1:las, i), matrices%indXs(1:las, :, i), q, nm)
      a290 = sparse_matmul(matrices%rotXs(1:las, i), matrices%indXs(1:las, :, i), a90, nm)
      b290 = sparse_matmul(matrices%rotXs(1:las, i), matrices%indXs(1:las, :, i), b90, nm)
      p290 = sparse_matmul(matrices%rotXs(1:las, i), matrices%indXs(1:las, :, i), p90, nm)
      q290 = sparse_matmul(matrices%rotXs(1:las, i), matrices%indXs(1:las, :, i), q90, nm)

      fx = F_z(Nmax, a2, b2, p2, q2) + F_z(Nmax, a290, b290, p290, q290)
      tx = T_z(Nmax, a2, b2, p2, q2) + T_z(Nmax, a290, b290, p290, q290)

! y-direction
      a2 = sparse_matmul(matrices%rotYs(1:las, i), matrices%indYs(1:las, :, i), a, nm)
      b2 = sparse_matmul(matrices%rotYs(1:las, i), matrices%indYs(1:las, :, i), b, nm)
      p2 = sparse_matmul(matrices%rotYs(1:las, i), matrices%indYs(1:las, :, i), p, nm)
      q2 = sparse_matmul(matrices%rotYs(1:las, i), matrices%indYs(1:las, :, i), q, nm)
      a290 = sparse_matmul(matrices%rotYs(1:las, i), matrices%indYs(1:las, :, i), a90, nm)
      b290 = sparse_matmul(matrices%rotYs(1:las, i), matrices%indYs(1:las, :, i), b90, nm)
      p290 = sparse_matmul(matrices%rotYs(1:las, i), matrices%indYs(1:las, :, i), p90, nm)
      q290 = sparse_matmul(matrices%rotYs(1:las, i), matrices%indYs(1:las, :, i), q90, nm)

      fy = F_z(Nmax, a2, b2, p2, q2) + F_z(Nmax, a290, b290, p290, q290)
      ty = T_z(Nmax, a2, b2, p2, q2) + T_z(Nmax, a290, b290, p290, q290)

! Averaging now, thus division by the number of polarization states
      if(matrices%polarization > 1) then
         matrices%force = dble([fx, fy, fz]/cc/2d0)
         matrices%torque = dble([tx, ty, tz]/(cc*mesh%k)/2d0)
      else
         matrices%force = dble([fx, fy, fz]/cc)
         matrices%torque = dble([tx, ty, tz]/(cc*mesh%k))
      end if
   end subroutine forcetorque

!****************************************************************************80

   subroutine relax_step(q, t, dq_max, J)
      real(dp), intent(inout) :: q, t
      real(dp), intent(in) :: J, dq_max
      real(dp) :: tau_n, tau_e, tau_int, dq, I1, I3, dt

      I1 = matrices%Ip(1)
      I3 = matrices%Ip(3)

! Calculate q-update
      tau_n = barnett_time(3d-11, 1.3d8, J/matrices%Ip(3), 1d-4, 1d-4)
      tau_e = barnett_time(1d-13, -1.76d11, J/matrices%Ip(3), 1d-10, 1d-6)
      tau_int = 1d0/(1d0/tau_n+1d0/tau_e)

      dq = -(q-1)*(1-q*I1/I3)/tau_int/(1-I1/I3)
      dt = abs(dq_max/dq)
      q = q + dq*dt
      t = t + dt

   end subroutine relax_step

!****************************************************************************80
! Calculates z-component of force analytically
   function F_z(Nmax, a, bb, p, qq) result(f)
      real(dp) :: f, g
      complex(dp), dimension(:), allocatable, intent(in) :: a, bb, p, qq
      complex(dp), dimension(:), allocatable :: b, q
      integer :: Nmax, n, m, ind

      ind = 0
      f = 0.d0
      b = bb*dcmplx(0d0, 1d0)
      q = qq*dcmplx(0d0, 1d0)
      do n = 1, Nmax
         do m = -n, n
            ind = ind + 1
            g = (2d0*m/(n*(n + 1d0)))*dimag(dconjg(a(ind))*b(ind) - dconjg(p(ind))*q(ind))
            if (ind + 2*n + 2 <= Nmax*(Nmax + 2)) then
               g = g &
                   - (2d0/(n + 1d0))*sqrt((n*(n + 2d0)*(n - m + 1d0)*(n + m + 1d0))/((2d0*n + 1d0)*(2d0*n + 3d0))) &
                   *dimag(a(ind)*dconjg(a(ind + 2*n + 2)) + b(ind)*dconjg(b(ind + 2*n + 2)) - &
                          p(ind)*dconjg(p(ind + 2*n + 2)) - q(ind)*dconjg(q(ind + 2*n + 2)))
            end if
            f = f + g
         end do
      end do
   end function F_z

!****************************************************************************80
! Calculates the z-component of torque analytically
   function T_z(Nmax, a, b, p, q) result(t)
      real(dp) :: t
      complex(dp), dimension(:), allocatable :: a, b, p, q
      integer :: Nmax, n, m, ind

      ind = 0
      t = 0.d0
      do n = 1, Nmax
         do m = -n, n
            ind = ind + 1
            t = t + dble(m)*(cdabs(a(ind))**2 + cdabs(b(ind))**2 - &
                             cdabs(p(ind))**2 - cdabs(q(ind))**2)
         end do
      end do

   end function T_z

!****************************************************************************80
! Calculates the z-component of torque analytically
   function barnett_torque() result(N)
      real(dp) :: a
      real(dp), dimension(3) :: N
      real(dp), dimension(3)    :: w, mu_Bar, B

! The approximate Barnett proportional constant X(0)*hbar/(g*mu_B)
      a = 10d0**(-14)
      w = matrices%w
      mu_Bar = a*mesh%V*w
      B = matrices%B

      N = crossRR(mu_Bar, B)
   end function barnett_torque

!****************************************************************************80
! Calculate the Barnett relaxation (nuclear/electron) time scale before next 
! step
   function barnett_time(chi0, g, w, T1, T2) result(tau)
      real(dp) :: tau, T1, T2, chi0, I1, I3, w, V, g

      V = mesh%V
      I1 = matrices%Ip(1)
      I3 = matrices%Ip(3)

      tau =(2d0*V*(chi0*T2)*(I3-I1))/((I1**2)*g**2)*w**2/(1+(I3*T1*T2*w**2)/(2d0*I1))
      tau = 1d0/tau

   end function barnett_time

!****************************************************************************80
! Calculates the z-component of torque analytically
   function DG_torque() result(N)
      real(dp) :: psi, xi, phi, Kw, V, tau, mag_B
      real(dp), dimension(3) :: N
      real(dp), dimension(3)    :: a3, w, B, B_perp, proj_Bperp_a3, psi_vec

      mag_B = vlen(matrices%B)
      Kw = 10d0**(-13)

      a3 = matrices%P(:, 3)
      V = mesh%V
      w = matrices%w
      tau = matrices%Ip(3)/(Kw*V*mag_B**2)
      B = matrices%B/mag_B
      B_perp = crossRR([0d0, 1d0, 0d0], B)
      proj_Bperp_a3 = a3 - dot_product(a3, B_perp)*a3
      proj_Bperp_a3 = proj_Bperp_a3/vlen(proj_Bperp_a3)

      psi = dacos(dot_product(B, matrices%khat))
      xi = dacos(dot_product(B, a3))
      phi = dacos(dot_product(B_perp, proj_Bperp_a3))

      if (psi > pi/2) then
         psi = psi - pi/2
      end if
      if (proj_Bperp_a3(2) < 0) phi = phi + pi

      psi_vec = [dcos(phi)*dcos(xi)*dcos(phi) - dsin(psi)*dsin(xi), &
                 dcos(xi)*dsin(phi), -dsin(psi)*dcos(xi)*dcos(phi) - dcos(psi)*dsin(xi)]

      N = -(dsin(xi)*dcos(xi)*psi_vec + dsin(xi)*dsin(xi)*a3)*matrices%Ip(3)*vlen(w)/tau

      if (mag_B < 1d-16) N = 0d0

   end function DG_torque

!****************************************************************************80
! Beta (rotation about a_3) averaged torque efficiency as a function of 
! capital theta and phi (as in the works of Draine and Weingartner)
   function get_Qav(theta, phi) result(Q)
      real(dp) :: Q(3), theta, phi, beta
      integer :: i, Nbeta

      Q = 0d0
      Nbeta = 20

      do i = 1,Nbeta
         beta = dble(i)*pi*2d0/Nbeta
         Q = Q + matmul(matrices%Rk,get_Q(theta,phi,beta))/Nbeta
      end do

   end function get_Qav

!****************************************************************************80
! Torque efficiency as a function of capital theta, phi and beta (as in the 
! works of Draine and Weingartner)
   function get_Q(theta, phi, beta) result(Q)
      real(dp) :: R_thta(3, 3), nbeta(3), R_beta(3, 3), Q(3), theta, beta, phi, &
                  a_3(3), R_phi(3, 3)

      a_3 = matrices%P(1:3, 3)

      R_thta = R_theta(theta)
      R_phi = R_aa([0d0, 0d0, 1d0], phi)

! Rotation axis for beta averaging for current theta
      nbeta = matmul(matmul(R_phi, R_thta), a_3) ! Beta rotation about a_3
      nbeta = nbeta/vlen(nbeta) ! Ensure unit length of axis vector

      R_beta = R_aa(nbeta, beta)

! The ultimate rotation matrices for scattering event
      matrices%R = matmul(R_beta, R_thta) ! First a_3 to theta, then beta about a_3

      call get_forces()
      Q = matrices%Q_t

   end function get_Q

!****************************************************************************80
! Calculates orientation angles in terms of alignment angles xi, psi and phi
   function get_thetaphi(xi, psi, phi) result(tp)
      real(dp) :: xi, psi, phi, tp(2)

      tp(1) = dacos(cos(psi)*cos(xi) - sin(psi)*sin(xi)*cos(phi))
      tp(2) = 2d0*datan2(sin(tp(1)) - sin(psi)*cos(xi) - &
                         cos(psi)*sin(xi)*cos(phi), sin(xi)*sin(phi))

   end function get_thetaphi

!****************************************************************************80
! Calculates alignment angles xi and phi from a given rotation matrix
   function get_xiphi() result(xp)
      real(dp) :: R(3,3), xi, phi, xp(2), B(3), a_3_perp_z(3)

      R = matrices%R
      B = matrices%B/vlen(matrices%B)
      a_3_perp_z = matrices%P(:,3)
      xi = dacos(dot_product(B, matmul(R,a_3_perp_z)))

      a_3_perp_z(3) = 0d0
      a_3_perp_z = a_3_perp_z/vlen(a_3_perp_z)

      phi = dacos(dot_product([1d0,0d0,0d0],a_3_perp_z))
      if(a_3_perp_z(2)<0d0) phi = phi + pi

      xp = [xi,phi]

   end function get_xiphi

end module forces
