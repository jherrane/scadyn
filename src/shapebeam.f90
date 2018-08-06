module shapebeam
   use T_matrix
   use bessel
   use gaussquad

   implicit none
contains

!****************************************************************************80
! The Gaussian amplitude profile in localized approximation (kw0>=5) a'la
! MSTM 3.0 (Mackowski et al 2013)
   subroutine gaussian_beams()
      real(dp) :: width
      integer :: i

      width = 10d0/(maxval(mesh%ki))/sqrt(2d0)
      if(matrices%whichbar /= 0)then
         call gaussian_beam_shape(matrices%whichbar, matrices%Nmaxs(i), width)
      else
         do i = 1, matrices%bars
            call gaussian_beam_shape(i, matrices%Nmaxs(i), width)
         end do
      end if

   end subroutine gaussian_beams

!****************************************************************************80
! The Gaussian amplitude profile in localized approximation (kw0>=5) a'la
! MSTM 3.0 (Mackowski et al 2013)
   subroutine gaussian_beam_shape(i, Nmax, width)
      real(dp) :: gn, kw0, width
      integer :: n, m, ind, i, Nmax

      kw0 = mesh%ki(i)*width

      if (kw0 < 5d0) then
         write (*, '(A)') "    Problem: width of Gaussian beam at focal point is smaller than wavelength!"
         write (*, '(2(A, ES9.3))') "     Wavelength is ", 2d0*pi/mesh%ki(i), ", width is ", width
      end if

      ind = 0
      do n = 1, Nmax
         gn = dexp(-((dble(n) + .5d0)/kw0)**2)
         do m = -n, n
            ind = ind + 1
            matrices%as(ind, i) = gn*matrices%as(ind, i)
            matrices%bs(ind, i) = gn*matrices%bs(ind, i)
         end do
      end do

   end subroutine gaussian_beam_shape

!****************************************************************************80

   subroutine laguerre_gaussian_beams(p, l)
      real(dp) :: width
      integer :: i, p, l

      width = 1d0/(maxval(mesh%ki))

      if(matrices%whichbar /= 0)then
         call laguerre_gauss_farfield(matrices%whichbar, p, l, width)
      else
         do i = 1, matrices%bars
            call laguerre_gauss_farfield(i, p, l, width)
         end do
      end if

   end subroutine laguerre_gaussian_beams

!****************************************************************************80
! Axisymmetric LG-beam
   subroutine laguerre_gauss_farfield(i, p, l, w0)
      integer :: i, nmax, total_modes, iii, jjj, ntheta, nphi, p, l, tp, info, lwork, ind
      integer, dimension(:), allocatable :: nn, mm, nn_old
      complex(dp) :: x, y, BCP(9)
      real(dp) :: k, w0, truncation_angle
      real(dp), dimension(:), allocatable :: theta, phi, rw, LL, work
      complex(dp), dimension(:), allocatable :: beam_envelope, Ex, Ey, &
                                                Etheta, Ephi, e_field, expansion_coefficients, &
                                                a, b, a_nm, b_nm
      complex(dp), dimension(:, :), allocatable :: coefficient_matrix

      k = mesh%ki(i)
      nmax = matrices%Nmaxs(i)
      allocate (a_nm((nmax + 1)**2 - 1), b_nm((nmax + 1)**2 - 1))
      a_nm = dcmplx(0d0)
      b_nm = dcmplx(0d0)
      truncation_angle = 90d0
      x = dcmplx(1d0, 0d0)
      y = dcmplx(0d0, 1d0)

      
      ! total_modes = nmax**2 + 2*nmax
      total_modes = 2*nmax
      allocate (nn(total_modes), mm(total_modes), nn_old(total_modes))
      do iii = 1, total_modes
         nn(iii) = ceiling(dble(iii)/2d0)
         mm(iii) = int((-1d0)**(iii) + l)
         ! nn(iii) = floor(sqrt(dble(iii)))
         ! mm(iii) = iii-nn(iii)**2-nn(iii)
      end do
      nn_old = nn
      nn = PACK(nn, nn_old >= abs(mm))
      mm = PACK(mm, nn_old >= abs(mm))
      ntheta = nmax + 1
      nphi = 2*p+abs(l)
      nphi = nphi+3-mod(nphi,2)
      tp = ntheta*nphi
      allocate (theta(tp), phi(tp))
      call angular_grid(ntheta, nphi, theta, phi)
      allocate (e_field(2*tp), rw(tp), LL(tp), beam_envelope(tp), Ex(tp), Ey(tp), &
                Etheta(tp), Ephi(tp))

      do iii = 1, tp
         rw(iii) = 2d0*k**2*w0**2*(dtan(theta(iii)))**2
      end do

      LL = laguerre(p, abs(l), rw)

      do iii = 1, tp
         beam_envelope(iii) = dcmplx(rw(iii)**(abs(l/2))*LL(iii)* &
                                   zexp(dcmplx(-rw(iii)/2, l*phi(iii) + pi/2d0*(p*2+abs(l)+1))))
         if (theta(iii) < pi*(180-truncation_angle)/180) beam_envelope(iii) = dcmplx(0d0)
      end do

      Ex = x*beam_envelope
      Ey = y*beam_envelope
      do iii = 1, tp
         Etheta(iii) = -Ex(iii)*dcos(phi(iii)) - Ey(iii)*dsin(phi(iii))
         Ephi(iii) = -Ex(iii)*dsin(phi(iii)) + Ey(iii)*dcos(phi(iii))
      end do

      e_field = [Etheta, Ephi]

      allocate (coefficient_matrix(2*tp, 2*size(nn, 1)), expansion_coefficients(2*size(nn, 1)))

      do iii = 1, size(nn, 1)
         do jjj = 1, tp
            BCP = vsh(nn(iii), mm(iii), theta(jjj), phi(jjj))
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

      lwork = 64*min(2*tp, 2*size(nn, 1))
      allocate (work(lwork))

      call zgels('N', 2*tp, 2*size(nn, 1), 1, coefficient_matrix, 2*tp, e_field, &
                 2*tp, work, lwork, info)

      expansion_coefficients = e_field(1:2*size(nn, 1))

      allocate (a(size(nn, 1)), b(size(nn, 1)))
      a = expansion_coefficients(1:size(nn, 1))
      b = expansion_coefficients(size(nn, 1) + 1:2*size(nn, 1))

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

      lim = 3.5d0
      n = 100
      nn = n*n
      allocate (z(n), y(n), grid(3, nn))

      call linspace(-lim, lim, n, z)
      call linspace(-lim, lim, n, y)
      grid = 0d0
      do i = 1, n
         do j = 1, n
            ind = n*(j - 1) + i
            grid(1, ind) = z(i)*mesh%a
            grid(2, ind) = y(j)*mesh%a
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

!****************************************************************************80

   subroutine bsc_farfield(nn,mm,E,theta,phi,zero_rejection_level,a,b)
      integer, dimension(:), allocatable, intent(in) :: nn, mm
      real(dp), dimension(:), allocatable, intent(in) :: theta, phi
      complex(dp), dimension(:), allocatable, intent(in) :: E 
      complex(dp), dimension(:), allocatable, intent(out) :: a, b
      real(dp) :: zero_rejection_level 

   end subroutine bsc_farfield

end module shapebeam
