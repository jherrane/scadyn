module possu
! Copyright (c) 2018 Johannes Markkanen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common
   use, intrinsic :: iso_c_binding

   implicit none

!include '/usr/local/include/fftw3.f03'
   include 'fftw3.f03'
!include 'fftw3-mpi.f03'
contains

!****************************************************************************80
! 3D FFT (fftw3)
   subroutine FFT3d(A)
      complex(dp), dimension(:, :, :) :: A

      integer*8 :: plan
      integer :: Nx, Ny, Nz, iret, t1, t2, rate

      Nx = size(A, 1)
      Ny = size(A, 2)
      Nz = size(A, 3)

      call dfftw_plan_dft_3d(plan, Nx, Ny, Nz, A, A, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(plan, A, A)
      call dfftw_destroy_plan(plan)

   end subroutine FFT3d

!****************************************************************************80
! 3D FFT-mpi routine (fftw3)
   subroutine FFT3d_2(A)
      complex(dp), dimension(:, :, :) :: A

      integer*8 :: plan
      integer :: Nx, Ny, Nz, iret, t1, t2, rate, i1, i2
      complex(dp) :: vec_x(size(A, 1))
      complex(dp) :: vec_y(size(A, 2))
      complex(dp) :: vec_z(size(A, 3))

      Nx = size(A, 1)
      Ny = size(A, 2)
      Nz = size(A, 3)

      !______________________ 1d-FFT x-direction_________________________!

      call dfftw_plan_dft_1d(plan, Nx, vec_x, vec_x, FFTW_FORWARD, FFTW_MEASURE)

!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Nz, Ny, plan)  &
!$OMP shared(A)
!$OMP do
      do i1 = 1, Nz/2
         do i2 = 1, Ny/2
            vec_x = A(:, i2, i1)
            call dfftw_execute_dft(plan, vec_x, vec_x)
            A(:, i2, i1) = vec_x
         end do
      end do
!$OMP end do
!$OMP end parallel

      call dfftw_destroy_plan(plan)

      !______________________ 1d-FFT y-direction_________________________!

      call dfftw_plan_dft_1d(plan, Ny, vec_y, vec_y, FFTW_FORWARD, FFTW_MEASURE)
!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Nz, Nx, plan)  &
!$OMP shared(A)
!$OMP do

      do i1 = 1, Nz/2
         do i2 = 1, Nx
            vec_y = A(i2, :, i1)
            call dfftw_execute_dft(plan, vec_y, vec_y)
            A(i2, :, i1) = vec_y
         end do
      end do
!$OMP end do
!$OMP end parallel

      call dfftw_destroy_plan(plan)

!______________________ 1d-FFT z-direction_________________________!

      call dfftw_plan_dft_1d(plan, Nz, vec_z, vec_z, FFTW_FORWARD, FFTW_MEASURE)
!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Ny, Nx, plan)  &
!$OMP shared(A)
!$OMP do

      do i1 = 1, Ny
         do i2 = 1, Nx
            vec_z = A(i2, i1, :)
            call dfftw_execute_dft(plan, vec_z, vec_z)
            A(i2, i1, :) = vec_z
         end do
      end do
!$OMP end do
!$OMP end parallel

      call dfftw_destroy_plan(plan)

   end subroutine FFT3d_2

!****************************************************************************80
! unscaled 3D inverse FFT (fftw3)
! scaling in function arr2vec (common.f90)
   subroutine IFFT3d(A)
      complex(dp), dimension(:, :, :) :: A

      integer*8 :: plan
      integer :: Nx, Ny, Nz, iret

      Nx = size(A, 1)
      Ny = size(A, 2)
      Nz = size(A, 3)

      call dfftw_plan_dft_3d(plan, Nx, Ny, Nz, A, A, FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(plan, A, A)
      call dfftw_destroy_plan(plan)

   end subroutine IFFT3d

!****************************************************************************80
! 3D FFT routine for the zero-padded grid (fftw3)
   subroutine IFFT3d_2(A)
      complex(dp), dimension(:, :, :) :: A

      integer*8 :: plan
      integer :: Nx, Ny, Nz, iret, t1, t2, rate, i1, i2
      complex(dp) :: vec_x(size(A, 1))
      complex(dp) :: vec_y(size(A, 2))
      complex(dp) :: vec_z(size(A, 3))

      Nx = size(A, 1)
      Ny = size(A, 2)
      Nz = size(A, 3)

!______________________ 1d-IFFT z-direction_________________________!

      call dfftw_plan_dft_1d(plan, Nz, vec_z, vec_z, FFTW_BACKWARD, FFTW_MEASURE)
!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Ny, Nx, plan)  &
!$OMP shared(A)
!$OMP do

      do i1 = 1, Ny
         do i2 = 1, Nx
            vec_z = A(i2, i1, :)
            call dfftw_execute_dft(plan, vec_z, vec_z)
            A(i2, i1, :) = vec_z
         end do
      end do
!$OMP end do
!$OMP end parallel

      call dfftw_destroy_plan(plan)

      !______________________ 1d-IFFT y-direction_________________________!

      call dfftw_plan_dft_1d(plan, Ny, vec_y, vec_y, FFTW_BACKWARD, FFTW_MEASURE)
!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Nz, Nx, plan)  &
!$OMP shared(A)
!$OMP do

      do i1 = 1, Nz/2
         do i2 = 1, Nx
            vec_y = A(i2, :, i1)
            call dfftw_execute_dft(plan, vec_y, vec_y)
            A(i2, :, i1) = vec_y
         end do
      end do
!$OMP end do
!$OMP end parallel

      call dfftw_destroy_plan(plan)

      !______________________ 1d-IFFT x-direction_________________________!

      call dfftw_plan_dft_1d(plan, Nx, vec_x, vec_x, FFTW_BACKWARD, FFTW_MEASURE)

!$OMP parallel num_threads(nthreads) default(private) &
!$OMP firstprivate(Nz, Ny, plan)  &
!$OMP shared(A)
!$OMP do
      do i1 = 1, Nz/2
         do i2 = 1, Ny/2
            vec_x = A(:, i2, i1)
            call dfftw_execute_dft(plan, vec_x, vec_x)
            A(:, i2, i1) = vec_x
         end do
      end do
!$OMP end do
!$OMP end parallel

      call dfftw_destroy_plan(plan)

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
