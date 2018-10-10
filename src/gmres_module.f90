module gmres_module
! Copyright (c) 2018 Johannes and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common
   use possu
   use sparse
   use sparse_mat
   implicit none

contains

!****************************************************************************80

   subroutine compute_Ax2(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh

      double complex, dimension(:, :, :), allocatable :: SS
      double complex, dimension(:), allocatable :: Y_aim, Y_xx, Y_yy, Y_zz, Y2_x, Y2_y, Y2_z, y_new, tmp
      double complex, dimension(:, :), allocatable :: X_xyz
      integer :: x, y, z, i1, mm, T1, T2, rate, i2, Nbasis
      double precision :: k2, scale

      if (mesh%order == 0) then
         Nbasis = 1
      end if
      if (mesh%order == 1) then
         Nbasis = 4
      end if

      k2 = mesh%k**2
      mm = mesh%Nx*mesh%Ny*mesh%Nz

      scale = 1.0/(mm*8.0)

      allocate (X_xyz(Nbasis*mesh%N_tet, 3))

      x = size(matrices%Fg, 1)
      y = size(matrices%Fg, 2)
      z = size(matrices%Fg, 3)

      allocate (SS(x, y, z))
      allocate (tmp(mm))

! Old solution vector

      if (mesh%order == 0) then
         do i1 = 1, mesh%N_tet
            X_xyz(i1, 1) = matrices%x(3*(i1 - 1) + 1)
            X_xyz(i1, 2) = matrices%x(3*(i1 - 1) + 2)
            X_xyz(i1, 3) = matrices%x(3*(i1 - 1) + 3)
         end do
      end if

      if (mesh%order == 1) then
         do i1 = 1, mesh%N_tet
            do i2 = 1, 4
               X_xyz(4*(i1 - 1) + i2, 1) = matrices%x(12*(i1 - 1) + i2)
               X_xyz(4*(i1 - 1) + i2, 2) = matrices%x(12*(i1 - 1) + 4 + i2)
               X_xyz(4*(i1 - 1) + i2, 3) = matrices%x(12*(i1 - 1) + 8 + i2)
            end do
         end do
      end if

!--------------- Projections into the grid nodes ----------------------------!

      call vec2arr2(SS, mesh, sparse_T_matmul(matrices%S, matrices%indS, X_xyz(:, 1), mm))

      call FFT3d_2(SS)

      SS = matrices%Fg*SS

      call IFFT3d_2(SS)

      Y_xx = sparse_eps_matmul(mesh%param, matrices%S, matrices%indS, arr2vec(mesh, SS), Nbasis)

      call vec2arr2(SS, mesh, sparse_T_matmul(matrices%S, matrices%indS, X_xyz(:, 2), mm))
      call FFT3d_2(SS)
      SS = matrices%Fg*SS
      call IFFT3d_2(SS)
      Y_yy = sparse_eps_matmul(mesh%param, matrices%S, matrices%indS, arr2vec(mesh, SS), Nbasis)

      call vec2arr2(SS, mesh, sparse_T_matmul(matrices%S, matrices%indS, X_xyz(:, 3), mm))
      call FFT3d_2(SS)
      SS = matrices%Fg*SS
      call IFFT3d_2(SS)
      Y_zz = sparse_eps_matmul(mesh%param, matrices%S, matrices%indS, arr2vec(mesh, SS), Nbasis)

      call vec2arr2(SS, mesh, sparse_T_matmul(matrices%Sx, matrices%indS, X_xyz(:, 1), mm))
      call FFT3d_2(SS)
      SS = matrices%Fg*SS
      call IFFT3d_2(SS)
      Y2_x = sparse_eps_matmul(mesh%param, matrices%Sx, matrices%indS, arr2vec(mesh, SS), Nbasis)
      Y2_y = sparse_eps_matmul(mesh%param, matrices%Sy, matrices%indS, arr2vec(mesh, SS), Nbasis)
      Y2_z = sparse_eps_matmul(mesh%param, matrices%Sz, matrices%indS, arr2vec(mesh, SS), Nbasis)

      call vec2arr2(SS, mesh, sparse_T_matmul(matrices%Sy, matrices%indS, X_xyz(:, 2), mm))
      call FFT3d_2(SS)
      SS = matrices%Fg*SS
      call IFFT3d_2(SS)
      Y2_x = Y2_x + sparse_eps_matmul(mesh%param, matrices%Sx, matrices%indS, arr2vec(mesh, SS), Nbasis)
      Y2_y = Y2_y + sparse_eps_matmul(mesh%param, matrices%Sy, matrices%indS, arr2vec(mesh, SS), Nbasis)
      Y2_z = Y2_z + sparse_eps_matmul(mesh%param, matrices%Sz, matrices%indS, arr2vec(mesh, SS), Nbasis)

      call vec2arr2(SS, mesh, sparse_T_matmul(matrices%Sz, matrices%indS, X_xyz(:, 3), mm))
      call FFT3d_2(SS)
      SS = matrices%Fg*SS
      call IFFT3d_2(SS)
      Y2_x = Y2_x + sparse_eps_matmul(mesh%param, matrices%Sx, matrices%indS, arr2vec(mesh, SS), Nbasis)
      Y2_y = Y2_y + sparse_eps_matmul(mesh%param, matrices%Sy, matrices%indS, arr2vec(mesh, SS), Nbasis)
      Y2_z = Y2_z + sparse_eps_matmul(mesh%param, matrices%Sz, matrices%indS, arr2vec(mesh, SS), Nbasis)

!-------------------------- Construct final vector--------------------------------------
      allocate (Y_aim(3*Nbasis*mesh%N_tet))

      if (mesh%order == 0) then
         do i1 = 1, mesh%N_tet
            Y_aim(3*(i1 - 1) + 1) = Y2_x(i1) - mesh%k**2*Y_xx(i1)
            Y_aim(3*(i1 - 1) + 2) = Y2_y(i1) - mesh%k**2*Y_yy(i1)
            Y_aim(3*(i1 - 1) + 3) = Y2_z(i1) - mesh%k**2*Y_zz(i1)
         end do
      end if

      if (mesh%order == 1) then
         do i1 = 1, mesh%N_tet
            do i2 = 1, 4
               Y_aim(12*(i1 - 1) + i2) = Y2_x(4*(i1 - 1) + i2) - mesh%k**2*Y_xx(4*(i1 - 1) + i2)
               Y_aim(12*(i1 - 1) + 4 + i2) = Y2_y(4*(i1 - 1) + i2) - mesh%k**2*Y_yy(4*(i1 - 1) + i2)
               Y_aim(12*(i1 - 1) + 8 + i2) = Y2_z(4*(i1 - 1) + i2) - mesh%k**2*Y_zz(4*(i1 - 1) + i2)
            end do
         end do
      end if

!--------------------- Multiply with the sparse matrix---------------------------------

      if (mesh%order == 0) then
         y_new = matmul_Acorr_const(matrices%sp_mat, matrices%sp_ind, matrices%x) + Y_aim
      end if

      if (mesh%order == 1) then
         y_new = matmul_Acorr(matrices%sp_mat, matrices%sp_ind, matrices%x) + Y_aim
      end if

      matrices%Ax = Y_new

   end subroutine compute_Ax2

!****************************************************************************80

   subroutine gmres(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh

      double complex, dimension(:), allocatable :: r, x, w, cs, sn, g, y
      double complex, dimension(:, :), allocatable :: v, h

      integer :: N, max_iter, k, j, i, iter, m, ite, T1, T2, rate, dest, ierr, N_procs
      double precision :: err_tol, b_norm, error, nu, normav, normav2
      double complex :: temp(2), tmp, hr

      err_tol = mesh%tol
      max_iter = mesh%restart ! restart number
      m = mesh%maxit ! number of iterations / restarts

      N = size(matrices%rhs)
      allocate (x(N))

      matrices%x = matrices%rhs ! Initial guess
      x = matrices%x
      b_norm = dble(sqrt(dot_product(matrices%rhs, matrices%rhs)))

      allocate (y(m))

      allocate (r(N), w(N))
      allocate (v(N, m + 1))
      allocate (h(m + 1, m))
      allocate (cs(m + 1), sn(m + 1), g(m + 1))

      v(:, :) = dcmplx(0.0, 0.0)
      h(:, :) = dcmplx(0.0, 0.0)

      cs(:) = dcmplx(0.0, 0.0)
      sn(:) = dcmplx(0.0, 0.0)

      w(:) = dcmplx(0.0, 0.0)
! GMRES ITERATIONS
      ite = 0

      print *, 'Start iterating'
      do iter = 1, max_iter

         matrices%x = x

         call compute_Ax2(matrices, mesh)

         r = matrices%rhs - matrices%Ax
         v(:, 1) = r/sqrt(dot_product(r, r))
         g(:) = dcmplx(0.0, 0.0)
         g(1) = sqrt(dot_product(r, r))

         do i = 1, m
            call system_clock(T1, rate)
            matrices%x = v(:, i)

            call compute_Ax2(matrices, mesh)
            w = matrices%Ax

            normav = dble(sqrt(dot_product(w, w)))

            !_______Modified Gram-Schmidt____________________________
            do k = 1, i
               h(k, i) = dot_product(v(:, k), w)
               w = w - h(k, i)*v(:, k)
            end do

            h(i + 1, i) = sqrt(dot_product(w, w))
            normav2 = dble(h(i + 1, i))
            v(:, i + 1) = w

            !_____________Reorthogonalize?________________________________
            if (normav + 0.001*normav2 == normav) then
               do j = 1, i
                  hr = dot_product(v(:, j), v(:, i + 1))
                  h(j, i) = h(j, i) + hr
                  v(:, i + 1) = v(:, i + 1) - hr*v(:, j)
               end do
               h(i + 1, i) = sqrt(dot_product(v(:, i + 1), v(:, i + 1)))
            end if
            !______________________________________________________

            if (h(i + 1, i) .ne. 0.0) then
               v(:, i + 1) = v(:, i + 1)/h(i + 1, i)
            end if

            !_____ apply Givens rotations_________________________________
            if (i > 1) then
               do k = 1, i - 1
                  tmp = cs(k)*h(k, i) - sn(k)*h(k + 1, i)
                  h(k + 1, i) = sn(k)*h(k, i) + conjg(cs(k))*h(k + 1, i)
                  h(k, i) = tmp
               end do
            end if
            !________________________________________________
            nu = dble(sqrt(dot_product(H(i:i + 1, i), H(i:i + 1, i))))

            if (nu .ne. 0.0) then

               cs(i) = conjg(h(i, i)/nu)
               sn(i) = -h(i + 1, i)/nu
               H(i, i) = cs(i)*H(i, i) - sn(i)*H(i + 1, i); 
               H(i + 1, i) = 0.0; 
               temp(1:2) = g(i:i + 1)
               g(i) = cs(i)*temp(1) - sn(i)*temp(2)
               g(i + 1) = sn(i)*temp(1) + conjg(cs(i))*temp(2)
            end if

            error = abs(g(i + 1))/b_norm; 
            if (error < err_tol) then
               y = matmul(Cinv(H(1:i, 1:i)), g(1:i)); 
               x = x + matmul(V(:, 1:i), y)
               exit
            end if

            call system_clock(T2)
            print *, 'RE (', ite + 1, ')', '=', real(error), 'time/iter =', real(T2 - T1)/real(rate)
            ite = ite + 1

            error = abs(g(i + 1))/b_norm; 
         end do

         if (error < err_tol) then
            exit
         end if

         y = matmul(Cinv(H(1:m, 1:m)), g(1:m)); 
         x = x + matmul(V(:, 1:m), y)
         matrices%x = x

         call compute_Ax2(matrices, mesh)
         r = matrices%rhs - matrices%Ax

         if (error < err_tol) then
            exit
         end if

      end do

      matrices%x = x

      print *, 'GMRES converged in', ite, 'iterations'

   end subroutine gmres

end module gmres_module
