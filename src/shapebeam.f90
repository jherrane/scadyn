module shapebeam
! Copyright (c) 2018 Joonas Herranen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use T_matrix
   use bessel

   implicit none
contains

!****************************************************************************80

   subroutine rotate_beam(i, axis, angle)
      integer :: i, nmax
      real(dp) :: axis(3), angle
      complex(dp), dimension(:), allocatable :: a_nm, b_nm
      complex(8), dimension(:), allocatable :: rotD
      integer, dimension(:, :), allocatable :: indD

      nmax = matrices%Nmaxs(i)
      allocate(rotD((nmax + 1)*(2*nmax + 1)*(2*nmax + 3)/3 - 1))
      allocate(indD((nmax + 1)*(2*nmax + 1)*(2*nmax + 3)/3 - 1,2))
      call sph_rotation_sparse_gen(mat2euler(R_aa(axis,angle)), nmax, rotD, indD)

      allocate(a_nm((nmax + 1)**2 - 1), b_nm((nmax + 1)**2 - 1))
      a_nm = matrices%as(1:(nmax + 1)**2 - 1, i)
      b_nm = matrices%bs(1:(nmax + 1)**2 - 1, i)

      matrices%as(1:(nmax + 1)**2 - 1, i) =  &
       sparse_matmul(rotD, indD, a_nm, (nmax + 1)**2 - 1)
      matrices%bs(1:(nmax + 1)**2 - 1, i) =  &
       sparse_matmul(rotD, indD, b_nm, (nmax + 1)**2 - 1)
   end subroutine rotate_beam

!****************************************************************************80

   subroutine laguerre_gaussian_beams(p, l)
      integer :: i, p, l, j
      real(dp) :: angle

      angle = 0d0
      j = 1
      if(matrices%whichbar/=0) j = matrices%whichbar
      matrices%width = 2d0*pi/(mesh%ki(j))
      if(matrices%whichbar /= 0)then
         call laguerre_gauss_farfield(matrices%whichbar, p, l)
         call rotate_beam(matrices%whichbar,[1d0,0d0,0d0],angle)
      else
         do i = 1, matrices%bars
            call laguerre_gauss_farfield(i, p, l)
            call rotate_beam(i,[1d0,0d0,0d0],angle)
         end do
      end if

   end subroutine laguerre_gaussian_beams

!****************************************************************************80
! Axisymmetric LG-beam
   subroutine laguerre_gauss_farfield(i, p, l)
      integer :: i, nmax, total_modes, iii, jjj, ntheta, nphi, p, l, tp, info, &
      lwork, ind, paraxial_order, row, p2, l2, naux, n, m
      integer, dimension(:), allocatable :: nn, mm, i1_out, i2_out, miv
      complex(dp) :: x, y, YY(3), BCP(9)
      real(dp) :: k, w0, truncation_angle, beam_angle, wscaling, mode_input_power, &
      aperture_power_normalization, NA, invL, zz, w, c, norm_paraxial
      real(dp), dimension(:), allocatable :: theta, phi, rw, LL, work, dr
      complex(dp), dimension(:), allocatable :: beam_envelope, Ex, Ey, &
                                                Etheta, Ephi, e_field, &
                                                a, b, a_nm, b_nm, a_nm2, b_nm2, fab
      complex(dp), dimension(:, :), allocatable :: coefficient_matrix
      real(dp), dimension(:,:), allocatable :: modeweights
      complex(8), dimension(:), allocatable :: rotD
      integer, dimension(:, :), allocatable :: indD

      k = mesh%ki(i)
! The numerical aperture can be above 1 only when refractive index of the medium
! is larger than unity.
      NA = matrices%NA
      if(NA/matrices%ref_med>1d0) then
         print*, '   Problem: Numerical aperture too high, limiting it to', matrices%ref_med
         NA = matrices%ref_med-1d-7
         matrices%NA = NA
      end if
      
      write (*, '(2(A,F5.3))') '  NA of the optical system =   ', matrices%NA 
      nmax = matrices%Nmaxs(i)
      allocate (a_nm((nmax + 1)**2 - 1), b_nm((nmax + 1)**2 - 1))

      a_nm = dcmplx(0d0)
      b_nm = dcmplx(0d0)
      truncation_angle = 90d0
      if(matrices%polarization == 1) then
         x = dcmplx(1d0, 0d0)
         y = dcmplx(0d0, 0d0)
      else if (matrices%polarization == 2) then
         x = dcmplx(1d0, 0d0)
         y = dcmplx(1d0, 0d0)
      else if (matrices%polarization == 3) then
         x = dcmplx(1d0, 0d0)
         y = dcmplx(0d0, 1d0)
      end if

      paraxial_order = 2*p+abs(l)
      allocate(i1_out(paraxial_order+1), i2_out(paraxial_order+1))
      do iii = 1, paraxial_order+1
         i2_out(iii) = -paraxial_order + 2*(iii-1)
         i1_out(iii) = floor(dble(paraxial_order-abs(i2_out(iii)))/2)
      end do
      row = (l+paraxial_order)/2+1
      p2 = i1_out(row)
      l2 = i2_out(row)
      w0 = paraxial_beam_waist(paraxial_order)
      beam_w0 = w0
      write (*, '(2(A,F7.3))') '  Beam waist               = ', w0

      total_modes = nmax**2 + 2*nmax
      allocate (nn(total_modes), mm(total_modes))

      naux = 0
      do iii = 1, total_modes
         nn(iii) = floor(sqrt(dble(iii)))
         mm(iii) = iii-nn(iii)**2-nn(iii)
         if(mm(iii)==l2+1 .OR. mm(iii)==l2-1) naux = naux + 1
      end do
      allocate(miv(naux))
      naux = 0
      do iii = 1, total_modes
         if(mm(iii)==l2+1 .OR. mm(iii)==l2-1) then
            naux = naux + 1
            miv(naux) = iii
         end if
      end do
      nn = nn(miv)
      mm = mm(miv)

      ntheta = 2*(nmax + 1)
      nphi = 3
      tp = ntheta*nphi
      allocate (theta(tp), phi(tp))
      call angular_grid(ntheta, nphi, theta, phi)
      allocate (e_field(2*tp), rw(tp), dr(tp), LL(tp), beam_envelope(tp), Ex(tp), &
                Ey(tp), Etheta(tp), Ephi(tp))
      
      beam_angle = dasin(NA/matrices%ref_med)
      write (*, '(2(A,F7.3))') '  Beam angle (deg)         = ', beam_angle*180d0/pi
      wscaling=1.0d0/dtan(abs(beam_angle))
      mode_input_power = 0d0
      aperture_power_normalization = 0d0
      do iii = 1, tp
         rw(iii) = 2d0*(wscaling*w0)**2*(dtan(theta(iii)))**2
         dr(iii) = (wscaling*w0)/dcos(theta(iii))**2
      end do

      LL = laguerre(p2, abs(l2), rw)

      norm_paraxial = sqrt(2d0*factorial(p2)/(pi*factorial(p2+abs(l2))))
      do iii = 1, tp
         beam_envelope(iii) = norm_paraxial*dcmplx(rw(iii)**(abs(l2/2))*LL(iii)* &
                              zexp(dcmplx(-rw(iii)/2, l2*phi(iii) + &
                              pi/2d0*(p2*2+abs(l2)+1))))
         mode_input_power = mode_input_power + &
            2d0*pi*abs(beam_envelope(iii))**2*sqrt(rw(iii)/2)*abs(dr(iii))
         aperture_power_normalization = aperture_power_normalization + &
            2d0*pi*abs(beam_envelope(iii))**2*dsin(theta(iii))
      end do
      mode_input_power = sqrt(mode_input_power)
      aperture_power_normalization = sqrt(aperture_power_normalization)

      do iii = 1, tp
         beam_envelope(iii) = beam_envelope(iii)*&
         mode_input_power/aperture_power_normalization
         if (theta(iii) < pi*(180-truncation_angle)/180) then
            beam_envelope(iii) = dcmplx(0d0)
         end if
      end do

      Ex = x*beam_envelope
      Ey = y*beam_envelope
      do iii = 1, tp
         Etheta(iii) = -Ex(iii)*dcos(phi(iii)) - Ey(iii)*dsin(phi(iii))
         Ephi(iii) = -Ex(iii)*dsin(phi(iii)) + Ey(iii)*dcos(phi(iii))
      end do

      e_field = [Etheta, Ephi]

      allocate(coefficient_matrix(2*tp, 2*size(nn, 1)))
      allocate(fab(2*size(nn,1)))

      do iii = 1, size(nn, 1)
         do jjj = 1, tp
            n = nn(iii)
            m = mm(iii)
! Get the spherical harmonics, only the derivative theta and phi components of the
! gradient are needed.
            YY = spharm(n, m, theta(jjj), phi(jjj))
! Coefficient matrix A is the solution to A*e_field(=B) = expansion_coefficients (=x)
            coefficient_matrix(jjj, iii) = &
            YY(3)*dcmplx(0d0, 1d0)**(nn(iii) + 1)/dsqrt(dble(nn(iii))*(nn(iii) + 1))

            coefficient_matrix(tp + jjj, iii) = &
            -YY(2)*dcmplx(0d0, 1d0)**(nn(iii) + 1)/ dsqrt(dble(nn(iii))*(nn(iii) + 1))
            
            coefficient_matrix(jjj, iii + size(nn, 1)) = &
            YY(2)*dcmplx(0d0, 1d0)**(nn(iii))/dsqrt(dble(nn(iii))*(nn(iii) + 1))

            coefficient_matrix(tp + jjj, iii + size(nn, 1)) = &
            YY(3)*dcmplx(0d0, 1d0)**(nn(iii))/dsqrt(dble(nn(iii))*(nn(iii) + 1))
         end do
      end do

! Solve the linear problem, same as x = A\B in matlab
      fab = matmul(pinv(coefficient_matrix),e_field)

! Solution is written in the e_field variable in zgels, then some book keeping
      allocate (a(size(nn, 1)), b(size(nn, 1)))
      a = fab(1:size(nn, 1))
      b = fab(size(nn, 1) + 1:2*size(nn, 1))

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

   function field_grid(n) result(grid)
      integer :: n, nn, i, j, ind
      real(dp) :: lim
      real(dp), allocatable :: z(:), y(:), grid(:, :)

      lim = 2d0*matrices%width
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

   end function field_grid

!****************************************************************************80
! Write fields to file. Normalize grid back to units of wavelength.
   subroutine write_fields()
      integer :: i, n
      character(LEN=80) :: gridname, fieldname, scatfield
      i = 1
      if(matrices%whichbar /= 0) i = matrices%whichbar
      n = 300
      gridname = 'grid_xyz.h5'
      fieldname = 'E_field.h5'
      scatfield = 'E_scat.h5'

      matrices%field_points = field_grid(n)

      call fields_out(i, n)
      call write2file(matrices%E_field, 'out/'//fieldname)

      call scat_fields_out(i, n)
      call write2file(matrices%E_field, 'out/'//scatfield)

      matrices%field_points = matrices%field_points/(2d0*pi/(mesh%ki(i)))
      call write2file(dcmplx(matrices%field_points), 'out/'//gridname)

   end subroutine write_fields

!****************************************************************************80

   function paraxial_beam_waist(paraxial_order) result(w)
      real(dp) :: w, w0, invL, expw, zz
      integer :: paraxial_order

      w = 1d0
      if(paraxial_order/=0)then
         invL = 1d0/dble(abs(paraxial_order))
         zz = exp(-(abs(paraxial_order)+2d0)*invL)
         w = -(1d0+2d0*sqrt(invL)+invL)
         w0 = -w

         do while(dabs(w-w0)>1d-5)
            w0 = w
            expw = exp(w)
            w = w0-(w0*expw+zz)/(expw+w0*expw)
         end do

         w = sqrt(-dble(abs(paraxial_order))/2d0*w)
      end if
   end function

end module shapebeam
