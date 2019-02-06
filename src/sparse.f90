module sparse
! Copyright (c) 2018 Johannes Markkanen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
!use omp_lib
   use common
   implicit none

contains

!****************************************************************************80

   function sparse_T_matmul(val, ind, x, N) result(Ax)
      double complex, dimension(:, :) :: val
      integer, dimension(:, :) :: ind
      double complex, dimension(:) :: x
      double complex, dimension(:), allocatable :: Ax
      integer :: N

      integer ::  rows, cols, i, j, row_t, t1, t2, rate

      allocate (Ax(N))

      Ax(:) = dcmplx(0.0, 0.0)
      rows = size(ind, 2)
      cols = size(ind, 1)

!!$omp parallel num_threads(nthreads) default(private) &
!!$omp firstprivate(cols, rows)  &
!!$omp shared(x, val, ind, Ax)

!!$omp do reduction(+:Ax)
      do i = 1, rows
         do j = 1, cols
            row_t = ind(j, i)
            Ax(row_t) = Ax(row_t) + val(j, i)*x(i)
         end do
      end do
!!$omp end do

!!$omp end parallel

   end function sparse_T_matmul

!****************************************************************************80

   function sparse_eps_matmul(eps, val, ind, x, Nbasis) result(Ax)
      double complex, dimension(:, :) :: val
      integer, dimension(:, :) :: ind
      double complex, dimension(:) :: x, eps
      double complex, dimension(:), allocatable :: Ax

      double complex :: eri
      integer :: N, rows, cols, i, j, tt, Nbasis

      N = size(ind, 2)
      allocate (Ax(N))
      Ax = dcmplx(0.0, 0.0)
      rows = size(ind, 2)
      cols = size(ind, 1)

!$OMP parallel num_threads(nthreads) default(none) &
!$OMP firstprivate(cols, rows, Nbasis)  &
!$OMP private(i,j,eri,tt) &
!$OMP shared(Ax,x, val, ind,eps)
!$OMP do
      do i = 1, rows
         !eri = eps((i-1)/Nbasis+1) - dcmplx(1.0, 0.0)
         do j = 1, cols
            !Ax(i) = Ax(i) + val(j,i) * eri * x(ind(j,i))
            Ax(i) = Ax(i) + val(j, i)*x(ind(j, i))
         end do
      end do
!$OMP end do
!$OMP end parallel

   end function sparse_eps_matmul

end module sparse
