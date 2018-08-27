module shapebeam
   use T_matrix
   use bessel
   use gaussquad

   implicit none
contains

!****************************************************************************80

   subroutine laguerre_gaussian_beams(p, l)
      integer :: i, p, l

      matrices%width = 2d0*pi/(minval(mesh%ki))
      if(matrices%whichbar /= 0)then
         call laguerre_gauss_farfield(matrices%whichbar, p, l)
      else
         do i = 1, matrices%bars
            call laguerre_gauss_farfield(i, p, l)
         end do
      end if

   end subroutine laguerre_gaussian_beams

!****************************************************************************80
! Axisymmetric LG-beam
   subroutine laguerre_gauss_farfield(i, p, l)
      integer :: i, nmax, total_modes, iii, jjj, ntheta, nphi, p, l, tp, info, &
      lwork, ind, paraxial_order
      integer, dimension(:), allocatable :: nn, mm, nn_old
      complex(dp) :: x, y, BCP(9)
      real(dp) :: k, w0, truncation_angle, beam_angle, wscaling, mode_input_power, &
      aperture_power_normalization, NA, invL, zz, w, expw
      real(dp), dimension(:), allocatable :: theta, phi, rw, LL, work, dr
      complex(dp), dimension(:), allocatable :: beam_envelope, Ex, Ey, &
                                                Etheta, Ephi, e_field, &
                                                a, b, a_nm, b_nm
      complex(dp), dimension(:, :), allocatable :: coefficient_matrix

      k = mesh%ki(i)
      NA = 0.9398d0
      nmax = matrices%Nmaxs(i)
      allocate (a_nm((nmax + 1)**2 - 1), b_nm((nmax + 1)**2 - 1))
      a_nm = dcmplx(0d0)
      b_nm = dcmplx(0d0)
      truncation_angle = 90d0
      x = dcmplx(1d0, 0d0)
      y = dcmplx(0d0, 1d0)

      paraxial_order = 2*p+abs(l)

      w0 = 1d0

      total_modes = nmax**2 + 2*nmax
      allocate (nn(total_modes), mm(total_modes), nn_old(total_modes))
      do iii = 1, total_modes
         nn(iii) = floor(sqrt(dble(iii)))
         mm(iii) = iii-nn(iii)**2-nn(iii)
      end do
      nn_old = nn
      nn = PACK(nn, nn_old >= abs(mm))
      mm = PACK(mm, nn_old >= abs(mm))
      ntheta = nmax + 1
      nphi = 2*(nmax+1)
      tp = ntheta*nphi
      allocate (theta(tp), phi(tp))
      call angular_grid(ntheta, nphi, theta, phi)
      allocate (e_field(2*tp), rw(tp), dr(tp), LL(tp), beam_envelope(tp), Ex(tp), Ey(tp), &
                Etheta(tp), Ephi(tp))

      beam_angle = dasin(NA)
      wscaling=1.0d0/dtan(abs(beam_angle))
      mode_input_power = 0d0
      aperture_power_normalization = 0d0
      do iii = 1, tp
         rw(iii) = 2d0*(wscaling*w0)**2*(dtan(theta(iii)))**2
         dr(iii) = (wscaling*w0)*dabs(dcos(theta(iii)))
      end do

      LL = laguerre(p, abs(l), rw)

      do iii = 1, tp
         beam_envelope(iii) = dcmplx(rw(iii)**(abs(l/2))*LL(iii)* &
                                   zexp(dcmplx(-rw(iii)/2, l*phi(iii) + pi/2d0*(p*2+abs(l)+1))))
         mode_input_power = mode_input_power + &
            2d0*pi*abs(beam_envelope(iii))**2*sqrt(rw(iii)/2)*abs(dr(iii))
         aperture_power_normalization = aperture_power_normalization + &
            2d0*pi*abs(beam_envelope(iii))**2*dsin(theta(iii))
      end do
      mode_input_power = sqrt(mode_input_power)
      aperture_power_normalization = sqrt(aperture_power_normalization)

      do iii = 1, tp
         beam_envelope(iii) = beam_envelope(iii)*mode_input_power/aperture_power_normalization
         if (theta(iii) < pi*(180-truncation_angle)/180) beam_envelope(iii) = dcmplx(0d0)
      end do

      Ex = x*beam_envelope
      Ey = y*beam_envelope
      do iii = 1, tp
         Etheta(iii) = -Ex(iii)*dcos(phi(iii)) - Ey(iii)*dsin(phi(iii))
         Ephi(iii) = -Ex(iii)*dsin(phi(iii)) + Ey(iii)*dcos(phi(iii))
      end do

      e_field = [Etheta, Ephi]

      allocate(coefficient_matrix(2*tp, 2*size(nn, 1)))

      do iii = 1, size(nn, 1)
         do jjj = 1, tp
! Get the spherical harmonics, only the derivative theta and phi components of the
! gradient are needed.
            BCP = vsh(nn(iii), mm(iii), theta(jjj), phi(jjj))
! Coefficient matrix A is the solution to A*e_field(=B) = expansion_coefficients (=x)
            coefficient_matrix(jjj, iii) = BCP(6)*dcmplx(0d0, 1d0)**(nn(iii) + 1)/ &
                                           dsqrt(dble(nn(iii))*(nn(iii) + 1))
            coefficient_matrix(tp + jjj, iii) = -BCP(5)*dcmplx(0d0, 1d0)**(nn(iii) + 1)/ &
                                                dsqrt(dble(nn(iii))*(nn(iii) + 1))
            coefficient_matrix(jjj, iii + size(nn, 1)) = BCP(5)*dcmplx(0d0, 1d0)**(nn(iii))/ &
                                                         dsqrt(dble(nn(iii))*(nn(iii) + 1))
            coefficient_matrix(tp + jjj, iii + size(nn, 1)) = BCP(6)*dcmplx(0d0, 1d0)**(nn(iii))/ &
                                                              dsqrt(dble(nn(iii))*(nn(iii) + 1))
         end do
      end do

      lwork = 64*max(2*tp, 2*size(nn, 1))
      allocate (work(lwork))

! Solve the linear problem, same as x = A\B in matlab
      call zgels('N', 2*tp, 2*size(nn, 1), 1, coefficient_matrix, 2*tp, e_field, &
                 2*tp, work, lwork, info)

! Solution is written in the e_field variable in zgels, then some book keeping
      allocate (a(size(nn, 1)), b(size(nn, 1)))
      a = e_field(1:size(nn, 1))
      b = e_field(size(nn, 1) + 1:2*size(nn, 1))

      do iii = 1, size(nn, 1)
         ind = nn(iii)*(nn(iii) + 1) + mm(iii)
         a_nm(ind) = a(iii)
         b_nm(ind) = b(iii)
         if (zabs(a(iii))**2 + zabs(b(iii))**2 < 1d-8) then
            a_nm(ind) = 0d0
            b_nm(ind) = 0d0
         end if
      end do

      matrices%as(1:(nmax + 1)**2 - 1, i) = a_nm
      matrices%bs(1:(nmax + 1)**2 - 1, i) = b_nm

   end subroutine laguerre_gauss_farfield


!****************************************************************************80

   subroutine fields_out(which, n)
      integer :: n, nn, i, which
      real(dp), allocatable :: grid(:, :)
      complex(dp) :: F(3), G(3)
      complex(dp), dimension(:, :), allocatable :: E

      nn = n*n

      allocate (E(3, nn))

      do i = 1, nn
         call calc_fields(matrices%as(:, which), matrices%bs(:, which), &
                          dcmplx(mesh%ki(which)), matrices%field_points(:, i), &
                          F, G, 0)
         E(:, i) = F
         ! E(:,i) = crossCC(F,G)/mu
      end do
      matrices%E_field = E

   end subroutine fields_out

!****************************************************************************80

   subroutine scat_fields_out(which, n)
      integer :: n, nn, i, which
      real(dp), allocatable :: grid(:, :)
      complex(dp), dimension(:), allocatable :: p, q, p90, q90
      complex(dp) :: F(3), G(3)
      complex(dp), dimension(:, :), allocatable :: E

      nn = n*n

! Ensure that directions are ok. They might already be...
      matrices%Rk = eye(3)
      allocate (E(3, nn))
      call rot_setup()
      call scattered_fields(1d0, p, q, p90, q90, which)

      do i = 1, nn
         call calc_fields(p, q, dcmplx(mesh%ki(which)), &
            matrices%field_points(:, i), F, G, 1)
         E(:, i) = F
         ! E(:,i) = crossCC(F,G)/mu
      end do

      matrices%E_field = E

   end subroutine scat_fields_out

!****************************************************************************80

   function field_grid() result(grid)
      integer :: n, nn, i, j, ind
      real(dp) :: lim
      real(dp), allocatable :: z(:), y(:), grid(:, :)

      lim = 2d0*matrices%width
      n = 100
      nn = n*n
      allocate (z(n), y(n), grid(3, nn))

      call linspace(-lim, lim, n, z)
      call linspace(-lim, lim, n, y)
      grid = 0d0
      do i = 1, n
         do j = 1, n
            ind = n*(j - 1) + i
            grid(1, ind) = z(i)
            grid(2, ind) = y(j)
         end do
      end do

      ! do i = 1,size(grid,2)
      !   write(*, '(3ES11.3)') grid(:,i)
      ! end do

   end function field_grid

!****************************************************************************80

   subroutine write_fields()
      integer :: i, n
      character(LEN=80) :: gridname, fieldname, scatfield
      i = 1
      if(matrices%whichbar /= 0) i = matrices%whichbar
      n = 100
      gridname = 'grid_xyz.h5'
      fieldname = 'E_field.h5'
      scatfield = 'E_scat.h5'

      matrices%field_points = field_grid()

      call fields_out(i, n)
      call write2file(dcmplx(matrices%field_points), gridname)
      call write2file(matrices%E_field, fieldname)

      call scat_fields_out(i, n)
      call write2file(matrices%E_field, scatfield)

   end subroutine write_fields

end module shapebeam
