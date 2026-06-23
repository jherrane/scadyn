module possu
! Copyright (c) 2018 Johannes Markkanen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common
   use, intrinsic :: iso_c_binding
   use omp_lib

   implicit none

!include '/usr/local/include/fftw3.f03'
   include 'fftw3.f03'
!include 'fftw3-mpi.f03'

! Cached FFTW plans (rebuilt when grid dimensions change)
   integer*8, save :: plan_3d_fwd = 0, plan_3d_bwd = 0
   integer, save :: fft3_Nx = -1, fft3_Ny = -1, fft3_Nz = -1

   integer*8, dimension(nthreads), save :: plan_1d_x_fwd = 0
   integer*8, dimension(nthreads), save :: plan_1d_y_fwd = 0
   integer*8, dimension(nthreads), save :: plan_1d_z_fwd = 0
   integer*8, dimension(nthreads), save :: plan_1d_x_bwd = 0
   integer*8, dimension(nthreads), save :: plan_1d_y_bwd = 0
   integer*8, dimension(nthreads), save :: plan_1d_z_bwd = 0
   integer, save :: fft2_Nx = -1, fft2_Ny = -1, fft2_Nz = -1

contains

!****************************************************************************80

   subroutine free_fftw_plans()

      if (plan_3d_fwd /= 0) then
         call dfftw_destroy_plan(plan_3d_fwd)
         plan_3d_fwd = 0
      end if
      if (plan_3d_bwd /= 0) then
         call dfftw_destroy_plan(plan_3d_bwd)
         plan_3d_bwd = 0
      end if
      fft3_Nx = -1
      fft3_Ny = -1
      fft3_Nz = -1

      call free_fft2_plans_only()
   end subroutine free_fftw_plans

!****************************************************************************80

   subroutine ensure_fft3_plans(Nx, Ny, Nz)
      integer, intent(in) :: Nx, Ny, Nz
      complex(dp), dimension(:,:,:), allocatable :: dummy

      if (fft3_Nx == Nx .and. fft3_Ny == Ny .and. fft3_Nz == Nz) return

      if (plan_3d_fwd /= 0) call dfftw_destroy_plan(plan_3d_fwd)
      if (plan_3d_bwd /= 0) call dfftw_destroy_plan(plan_3d_bwd)
      plan_3d_fwd = 0
      plan_3d_bwd = 0

      allocate (dummy(Nx, Ny, Nz))
      call dfftw_plan_dft_3d(plan_3d_fwd, Nx, Ny, Nz, dummy, dummy, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_plan_dft_3d(plan_3d_bwd, Nx, Ny, Nz, dummy, dummy, FFTW_BACKWARD, FFTW_ESTIMATE)
      deallocate (dummy)

      fft3_Nx = Nx
      fft3_Ny = Ny
      fft3_Nz = Nz
   end subroutine ensure_fft3_plans

!****************************************************************************80

   subroutine ensure_fft2_plans(Nx, Ny, Nz)
      integer, intent(in) :: Nx, Ny, Nz
      complex(dp), dimension(:), allocatable :: dummy_x, dummy_y, dummy_z
      integer :: t

      if (fft2_Nx == Nx .and. fft2_Ny == Ny .and. fft2_Nz == Nz) return

      call free_fft2_plans_only()

      allocate (dummy_x(Nx), dummy_y(Ny), dummy_z(Nz))
      do t = 1, nthreads
         call dfftw_plan_dft_1d(plan_1d_x_fwd(t), Nx, dummy_x, dummy_x, FFTW_FORWARD, FFTW_ESTIMATE)
         call dfftw_plan_dft_1d(plan_1d_y_fwd(t), Ny, dummy_y, dummy_y, FFTW_FORWARD, FFTW_ESTIMATE)
         call dfftw_plan_dft_1d(plan_1d_z_fwd(t), Nz, dummy_z, dummy_z, FFTW_FORWARD, FFTW_ESTIMATE)
         call dfftw_plan_dft_1d(plan_1d_x_bwd(t), Nx, dummy_x, dummy_x, FFTW_BACKWARD, FFTW_ESTIMATE)
         call dfftw_plan_dft_1d(plan_1d_y_bwd(t), Ny, dummy_y, dummy_y, FFTW_BACKWARD, FFTW_ESTIMATE)
         call dfftw_plan_dft_1d(plan_1d_z_bwd(t), Nz, dummy_z, dummy_z, FFTW_BACKWARD, FFTW_ESTIMATE)
      end do
      deallocate (dummy_x, dummy_y, dummy_z)

      fft2_Nx = Nx
      fft2_Ny = Ny
      fft2_Nz = Nz
   end subroutine ensure_fft2_plans

!****************************************************************************80

   subroutine free_fft2_plans_only()
      integer :: t

      do t = 1, nthreads
         if (plan_1d_x_fwd(t) /= 0) call dfftw_destroy_plan(plan_1d_x_fwd(t))
         if (plan_1d_y_fwd(t) /= 0) call dfftw_destroy_plan(plan_1d_y_fwd(t))
         if (plan_1d_z_fwd(t) /= 0) call dfftw_destroy_plan(plan_1d_z_fwd(t))
         if (plan_1d_x_bwd(t) /= 0) call dfftw_destroy_plan(plan_1d_x_bwd(t))
         if (plan_1d_y_bwd(t) /= 0) call dfftw_destroy_plan(plan_1d_y_bwd(t))
         if (plan_1d_z_bwd(t) /= 0) call dfftw_destroy_plan(plan_1d_z_bwd(t))
         plan_1d_x_fwd(t) = 0
         plan_1d_y_fwd(t) = 0
         plan_1d_z_fwd(t) = 0
         plan_1d_x_bwd(t) = 0
         plan_1d_y_bwd(t) = 0
         plan_1d_z_bwd(t) = 0
      end do
      fft2_Nx = -1
      fft2_Ny = -1
      fft2_Nz = -1
   end subroutine free_fft2_plans_only

!****************************************************************************80
! 3D FFT (fftw3)
   subroutine FFT3d(A)
      complex(dp), dimension(:, :, :) :: A

      integer :: Nx, Ny, Nz

      Nx = size(A, 1)
      Ny = size(A, 2)
      Nz = size(A, 3)

      call ensure_fft3_plans(Nx, Ny, Nz)
      call dfftw_execute_dft(plan_3d_fwd, A, A)

   end subroutine FFT3d

!****************************************************************************80
! 3D FFT on the zero-padded grid (fftw3)
   subroutine FFT3d_2(A)
      complex(dp), dimension(:, :, :) :: A

      integer :: Nx, Ny, Nz, i1, i2, tid
      complex(dp) :: vec_x(size(A, 1))
      complex(dp) :: vec_y(size(A, 2))
      complex(dp) :: vec_z(size(A, 3))

      Nx = size(A, 1)
      Ny = size(A, 2)
      Nz = size(A, 3)

      call ensure_fft2_plans(Nx, Ny, Nz)

      !______________________ 1d-FFT x-direction_________________________!

!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Nz, Ny) shared(A, plan_1d_x_fwd)
!$OMP do
      do i1 = 1, Nz/2
         do i2 = 1, Ny/2
            tid = omp_get_thread_num() + 1
            vec_x = A(:, i2, i1)
            call dfftw_execute_dft(plan_1d_x_fwd(tid), vec_x, vec_x)
            A(:, i2, i1) = vec_x
         end do
      end do
!$OMP end do
!$OMP end parallel

      !______________________ 1d-FFT y-direction_________________________!

!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Nz, Nx) shared(A, plan_1d_y_fwd)
!$OMP do

      do i1 = 1, Nz/2
         do i2 = 1, Nx
            tid = omp_get_thread_num() + 1
            vec_y = A(i2, :, i1)
            call dfftw_execute_dft(plan_1d_y_fwd(tid), vec_y, vec_y)
            A(i2, :, i1) = vec_y
         end do
      end do
!$OMP end do
!$OMP end parallel

!______________________ 1d-FFT z-direction_________________________!

!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Ny, Nx) shared(A, plan_1d_z_fwd)
!$OMP do

      do i1 = 1, Ny
         do i2 = 1, Nx
            tid = omp_get_thread_num() + 1
            vec_z = A(i2, i1, :)
            call dfftw_execute_dft(plan_1d_z_fwd(tid), vec_z, vec_z)
            A(i2, i1, :) = vec_z
         end do
      end do
!$OMP end do
!$OMP end parallel

   end subroutine FFT3d_2

!****************************************************************************80
! unscaled 3D inverse FFT (fftw3)
! scaling in function arr2vec (common.f90)
   subroutine IFFT3d(A)
      complex(dp), dimension(:, :, :) :: A

      integer :: Nx, Ny, Nz

      Nx = size(A, 1)
      Ny = size(A, 2)
      Nz = size(A, 3)

      call ensure_fft3_plans(Nx, Ny, Nz)
      call dfftw_execute_dft(plan_3d_bwd, A, A)

   end subroutine IFFT3d

!****************************************************************************80
! 3D inverse FFT on the zero-padded grid (fftw3)
   subroutine IFFT3d_2(A)
      complex(dp), dimension(:, :, :) :: A

      integer :: Nx, Ny, Nz, i1, i2, tid
      complex(dp) :: vec_x(size(A, 1))
      complex(dp) :: vec_y(size(A, 2))
      complex(dp) :: vec_z(size(A, 3))

      Nx = size(A, 1)
      Ny = size(A, 2)
      Nz = size(A, 3)

      call ensure_fft2_plans(Nx, Ny, Nz)

!______________________ 1d-IFFT z-direction_________________________!

!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Ny, Nx) shared(A, plan_1d_z_bwd)
!$OMP do

      do i1 = 1, Ny
         do i2 = 1, Nx
            tid = omp_get_thread_num() + 1
            vec_z = A(i2, i1, :)
            call dfftw_execute_dft(plan_1d_z_bwd(tid), vec_z, vec_z)
            A(i2, i1, :) = vec_z
         end do
      end do
!$OMP end do
!$OMP end parallel

      !______________________ 1d-IFFT y-direction_________________________!

!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Nz, Nx) shared(A, plan_1d_y_bwd)
!$OMP do

      do i1 = 1, Nz/2
         do i2 = 1, Nx
            tid = omp_get_thread_num() + 1
            vec_y = A(i2, :, i1)
            call dfftw_execute_dft(plan_1d_y_bwd(tid), vec_y, vec_y)
            A(i2, :, i1) = vec_y
         end do
      end do
!$OMP end do
!$OMP end parallel

      !______________________ 1d-IFFT x-direction_________________________!

!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Nz, Ny) shared(A, plan_1d_x_bwd)
!$OMP do
      do i1 = 1, Nz/2
         do i2 = 1, Ny/2
            tid = omp_get_thread_num() + 1
            vec_x = A(:, i2, i1)
            call dfftw_execute_dft(plan_1d_x_bwd(tid), vec_x, vec_x)
            A(:, i2, i1) = vec_x
         end do
      end do
!$OMP end do
!$OMP end parallel

   end subroutine IFFT3d_2

!****************************************************************************80
! Inverse of a real matrix (lapack)
   function inv(A) result(Ainv)
      real(dp), dimension(:, :), intent(in) :: A
      real(dp), dimension(:, :), allocatable :: Ainv

      real(dp), dimension(:), allocatable :: work  ! work array for LAPACK
      integer, dimension(:), allocatable :: ipiv   ! pivot indices
      integer :: n, info

      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      n = size(A, 1)
      allocate (Ainv(n, n))
      allocate (work(n))
      allocate (ipiv(n))
      Ainv = A

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
   end function inv

!****************************************************************************80
! Inverse of a complex matrix (lapack)
   function Cinv(A) result(Ainv)
      complex(dp), dimension(:, :), intent(in) :: A
      complex(dp), dimension(:, :), allocatable :: Ainv

      complex(dp), dimension(:), allocatable :: work  ! work array for LAPACK
      integer, dimension(:), allocatable :: ipiv   ! pivot indices
      integer :: n, info

      ! External procedures defined in LAPACK
      external ZGETRF
      external ZGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      n = size(A, 1)
      allocate (Ainv(n, n))
      allocate (work(n))
      allocate (ipiv(n))
      Ainv = A

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call ZGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call ZGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
   end function Cinv

!****************************************************************************80
!       SVD (lapack)
   subroutine svd(A, U, S, VT, M, N)
      complex(dp), dimension(:, :), intent(in) :: A

      complex(dp) :: U(M, M), VT(N, N)
      real(dp) :: S(N), RWORK(5*N)

      complex(dp), dimension(:), allocatable :: WORK
      INTEGER LDA, LDU, M, N, LWORK, LDVT, INFO
      CHARACTER JOBU, JOBVT

      JOBU = 'A'
      JOBVT = 'A'
      LDA = M
      LDU = M
      LDVT = N

      LWORK = MAX(1, 3*MIN(M, N) + MAX(M, N), 5*MIN(M, N))

      ALLOCATE (work(lwork))

      CALL ZGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO)

   end subroutine svd

!****************************************************************************80
! Pseudoinverse
   function pinv(A) result(Ainv)
      complex(dp), dimension(:, :), intent(in) :: A
      complex(dp), dimension(:, :), allocatable :: Ainv

      complex(dp), dimension(:, :), allocatable :: U, V, SU
      real(dp), dimension(:), allocatable :: S
      integer :: M, N, i1

      M = size(A, 1)
      N = size(A, 2)

      allocate (U(M, M), V(N, N), S(N))
      allocate (SU(N, M))
      allocate (Ainv(N, M))

      call svd(A, U, S, V, M, N)

      do i1 = 1, N
         if (s(i1) > 1e-6) then
            SU(i1, :) = 1/s(i1)*conjg(U(:, i1))
         else
            SU(i1, :) = 0.0*conjg(U(:, i1))
         end if

      end do

      Ainv = matmul(transpose(conjg(V)), SU)
   end function pinv

end module possu
