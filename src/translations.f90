module translations
! Copyright (c) 2018 Johannes Markkanen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common
   use sfunctions
   use integration_points
   implicit none

contains

!****************************************************************************80

   function a_mn(m, n) result(a)
      integer :: m, n
      real(dp) :: a

      if (n < abs(m)) then
         a = 0.0d0
      else
         a = sqrt(dble((n + 1 + abs(m))*(n + 1 - abs(m)))/dble((2*n + 1)*(2*n + 3)))
      end if

   end function a_mn

!****************************************************************************80

   function b_mn(m, n) result(b)
      integer :: m, n
      real(dp) :: b

      b = sqrt(dble((n - m - 1)*(n - m))/dble((2*n - 1)*(2*n + 1)))

   end function b_mn

!****************************************************************************80
! Translation matrix for the  scalar wave function along z-axis
! Output :: C^mn_nu
!
! Axial (mu = m) translation based on Chew's recursive formula
!
! Chew, W. C. "Recurrence relations for three-dimensional scalar addition theorem."
!  Journal of Electromagnetic Waves and Applications 6.1-4 (1992): 133-142
!
! see also. Nail A. Gumerov and Ramani Duraiswami,
!"Fast, Exact, and Stable Computation of Multipole Translation and
! Rotation Coefficients for the 3-D Helmholtz Equation"
! Note: incorrect sign in eq. (187)
   subroutine scalar_translator(r, k, Nmax, C, inout)
      real(dp) :: r
      complex(dp) :: k
      integer :: Nmax, inout
      integer :: NN, n, m, nu, mp
      complex(dp), dimension(3*Nmax + 5, 3*Nmax + 5, 3*Nmax + 5) ::  C
      complex(dp), dimension(3*Nmax + 6) :: sphj, sphy, sphh
      real(dp), dimension(3*Nmax + 6) :: sphjr, sphyr
      real(dp), dimension(3*Nmax + 5) ::a_mn_mat
      complex(dp) :: kr

      NN = 3*Nmax + 5
      kr = k*r

      C(:, :, :) = dcmplx(0.0, 0.0)

      call sphj2(NN, real(kr), NN, sphjr)
      sphj = dcmplx(sphjr, 0.0)
      sphh = sphankel(NN, kr)

! ____ C^00_nu______ !
      do nu = 0, NN - 1
         ! Regular
         if (inout == 0) then
            C(1, 1, nu + 1) = sqrt(2.0d0*nu + 1.0d0)*sphj(nu + 1)
            C(1, nu + 1, 1) = (-1)**nu*sqrt(2.0d0*nu + 1.0d0)*sphj(nu + 1)
         end if

         ! Outgoing
         if (inout == 1) then
            C(1, 1, nu + 1) = sqrt(2.0d0*nu + 1.0d0)*sphh(nu + 1)
            C(1, nu + 1, 1) = (-1)**nu*sqrt(2.0d0*nu + 1.0d0)*sphh(nu + 1)
         end if
      end do

!____ C^01_0 _______!
      C(1, 2, 1) = sqrt(3.0d0)*(-1.0d0*sqrt(1.0d0/3.0d0)*C(1, 1, 2))

      do mp = 1, NN
         a_mn_mat(mp) = a_mn(0, mp - 1)
      end do

! ____C^01_nu ______!

      do nu = 1, NN - 2

         C(1, 2, nu + 1) = -a_mn_mat(nu + 1)/a_mn_mat(1)*C(1, 1, nu + 2) &
                           + a_mn_mat(nu)/a_mn_mat(1)*C(1, 1, nu)

      end do

!____ C^0n_nu _____!

      do n = 1, NN - 1
         do nu = 1, NN - n - 2

            C(1, n + 2, nu + 1) = -a_mn_mat(nu + 1)/a_mn_mat(n + 1)*C(1, n + 1, nu + 2) &
                                  + a_mn_mat(nu)/a_mn_mat(n + 1)*C(1, n + 1, nu) &
                                  + a_mn_mat(n)/a_mn_mat(n + 1)*C(1, n, nu + 1)

         end do
      end do

! ___ C^mn_nu ____!

      do m = 0, NN - 2

         n = m

         do mp = 1, NN
            a_mn_mat(mp) = a_mn(m + 1, mp - 1) ! -1
         end do

         ! C(m,m,:)
         do nu = m + 1, NN - 2

            C(m + 2, n + 2, nu + 1) = (b_mn(-m - 1, nu)/b_mn(-m - 1, n + 1)*C(m + 1, n + 1, nu) &
                                       + b_mn(m, nu + 1)/b_mn(-m - 1, n + 1)*C(m + 1, n + 1, nu + 2)) !"Note sign difference to (187)

         end do

         do nu = m + 1, NN - 2

            C(m + 2, n + 3, nu + 1) = (-a_mn_mat(nu + 1)/a_mn_mat(n + 2)*C(m + 2, n + 2, nu + 2) &
                                       + a_mn_mat(nu)/a_mn_mat(n + 2)*C(m + 2, n + 2, nu)) ! nu +  1-> ind
         end do

         do n = m + 1, NN - 2
            do nu = m + 1, NN - 2 - n + m + 1

               C(m + 2, n + 2, nu + 1) = (-a_mn_mat(nu + 1)/a_mn_mat(n + 1)*C(m + 2, n + 1, nu + 2) &
                                          + a_mn_mat(nu)/a_mn_mat(n + 1)*C(m + 2, n + 1, nu) &
                                          + a_mn_mat(n)/a_mn_mat(n + 1)*C(m + 2, n, nu + 1))  ! nu +  1-> ind

            end do
         end do
      end do

   end subroutine scalar_translator

!****************************************************************************80
! Translation of the VSWF coefficient along z-axis
   subroutine vector_translator_sparse(r, k, N_in, N_out, AA, BB, indx, inout)
      real(dp) :: r
      complex(dp) :: k
      integer :: Nmax, N_in, N_out, inout

      complex(dp), dimension(N_out + 1, N_out, N_in) :: A, B, C0
      complex(dp), dimension(:) :: AA, BB
      integer, dimension(:, :) :: indx
      complex(dp), dimension(:, :, :), allocatable :: C
      integer :: m, n, nu, mu, las2, ind1, ind2, dm
      complex(dp) :: a1, a2
      real(dp) :: sc, aa1, aa2
      complex(dp), dimension(N_out + 1) :: a1_mat, a2_mat
      real(dp), dimension(N_out + 1) :: sc_mat

      A(:, :, :) = dcmplx(0.0, 0.0)
      B(:, :, :) = dcmplx(0.0, 0.0)
      C0(:, :, :) = dcmplx(0.0, 0.0)

      Nmax = maxval([N_in, N_out])
      allocate (C(3*Nmax + 5, 3*Nmax + 5, 3*Nmax + 5))
      call scalar_translator(r, k, Nmax, C, inout)

      do nu = 1, N_in
         do dm = 0, N_out
            aa1 = (nu - dm + 1.0d0)*(nu + dm + 1.0d0)/((2.0d0*nu + 1.0d0)*(2.0d0*nu + 3.0d0))
            a1_mat(dm + 1) = k*r/(nu + 1.0d0)*sqrt(dcmplx(aa1))

            aa2 = (nu - dm)*(nu + dm)/((2.0d0*nu + 1.0d0)*(2.0d0*nu - 1.0d0))
            a2_mat(dm + 1) = k*r/dble(nu)*sqrt(dcmplx(aa2))

            sc_mat(dm + 1) = sqrt(nu*(nu + 1.0d0))/sqrt(dcmplx((dm + 1)*((dm + 1) + 1.0d0)))
         end do

         do n = 1, N_out

            do m = 0, n

               a1 = a1_mat(m + 1)
               a2 = a2_mat(m + 1)
               sc = sc_mat(n)

               A(m + 1, n, nu) = sc*(C(m + 1, n + 1, nu + 1) - a1*C(m + 1, n + 1, nu + 2) - a2*C(m + 1, n + 1, nu))
               B(m + 1, n, nu) = sc*(-dcmplx(0.0, 1.0)*k*m*r/(nu*(nu + 1.0d0))*C(m + 1, n + 1, nu + 1))
               C0(m + 1, n, nu) = C(m + 1, n + 1, nu + 1)

            end do
         end do
      end do

      AA(:) = dcmplx(0.0, 0.0)
      BB(:) = dcmplx(0.0, 0.0)

      las2 = 0
      do nu = 1, N_in
         do mu = -nu, nu
            ind2 = nu + mu + nu**2
            do n = 1, N_out
               ind1 = n + mu + n**2
               m = mu

               if (abs(mu) <= nu .and. abs(mu) <= n) then
                  las2 = las2 + 1

                  indx(las2, 1) = ind1
                  indx(las2, 2) = ind2

                  if (m < 0) then
                     AA(las2) = A(abs(m) + 1, n, nu)
                     BB(las2) = B(abs(m) + 1, n, nu)
                  else
                     AA(las2) = A(abs(m) + 1, n, nu)
                     BB(las2) = -B(abs(m) + 1, n, nu)
                  end if

               end if

            end do
         end do
      end do

   end subroutine vector_translator_sparse

!****************************************************************************80
! Rotation matrix (recursive) for VSWFs
! Algorithm:  Choi, et al., J. Chem. Phys. vol. 11, No. 19, 1999.
   subroutine sph_rotation_sparse(theta, phi2, Nmax, spD, ind)
      real(dp) :: theta, phi, phi2
      integer :: Nmax
      complex(dp), dimension(:) :: spD
      integer, dimension(:, :) :: ind
      integer :: n, m, mu, ind1, ind2, mm, x, y, en, las
      real(dp) :: ct, st, cp, sp, R(3, 3)
      complex(dp) :: C(3, 3), invC(3, 3), D1(7, 7)
      complex(dp), dimension(:, :), allocatable :: locD
      complex(dp) :: a1, a2, a3
      real(dp) ::  a, b, b2, cm, dm, dmm, delta
      complex(dp), dimension((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1) :: D
      complex(dp), dimension(2*(Nmax + 2) + 1, 2*(Nmax + 2) + 1) :: DD

      phi = -1.0d0*phi2 !

      D(:, :) = dcmplx(0.0, 0.0)

      ct = cos(theta)
      st = sin(theta)
      cp = cos(phi)
      sp = sin(phi)

! YZ-rotation matrix
      R(1, :) = [ct*cp, -sp, st*cp]
      R(2, :) = [ct*sp, cp, st*sp]
      R(3, :) = [-st, dble(0.0), ct]

      a1 = 1.0d0/sqrt(2.0d0)
      a2 = dcmplx(0.0, 1.0)*(1.0d0)/sqrt(2.0d0)
      a3 = 1.0d0

      C(1, :) = [a1, dcmplx(0.0, 0.0), -a1]
      C(2, :) = [-a2, dcmplx(0.0, 0.0), -a2]
      C(3, :) = [0.0d0, 1.0d0, 0.0d0]

      invC(1, :) = [a1, a2, 0.0*a3]
      invC(2, :) = [0.0*a3, 0.0*a3, a3]
      invC(3, :) = [-a1, a2, 0.0*a3]

      D1(:, :) = dcmplx(0.0, 0.0)
      D1(3:5, 3:5) = (transpose(matmul(invC, matmul(R, C))))

      DD(:, :) = dcmplx(0.0, 0.0)
      DD(1:7, 1:7) = D1

      D(1:3, 1:3) = D1(3:5, 3:5)

      do n = 2, Nmax
         allocate (locD(2*n + 1, 2*n + 1))
         locD(:, :) = dcmplx(0.0, 0.0)
         do m = -n, n
            do mu = -n + 1, n - 1
               a = sqrt(dble(n + m)*(n - m)/(dble(n + mu)*(n - mu)))
               b = sqrt(dble(n + m)*(n + m - 1)/(2*dble(n + mu)*(n - mu)))

               mm = -m
               b2 = sqrt(dble(n + mm)*(n + mm - 1)/(2*dble(n + mu)*(n - mu)))

               x = mu + n + 1
               y = m + n + 1

               locD(x, y) = (D1(4, 4)*a*DD(x + 1, y + 1) &
                             + b*D1(4, 5)*DD(x + 1, y) &
                             + b2*D1(4, 3)*DD(x + 1, y + 2))
            end do

            mu = -n

            cm = sqrt(dble(n + m)*(n - m)/dble(n*(2*n - 1))); 
            dm = sqrt(dble(n + m)*(n + m - 1)/dble(2*n*(2*n - 1))); 
            dmm = sqrt(dble(n - m)*(n - m - 1)/dble(2*n*(2*n - 1))); 
            x = mu + n + 1; 
            y = m + n + 1; 
            locD(x, y) = (D1(3, 4)*cm*DD(x + 2, y + 1) &
                          + dm*D1(3, 5)*DD(x + 2, y) &
                          + dmm*D1(3, 3)*DD(x + 2, y + 2))

            mu = n; 
            x = mu + n + 1; 
            y = m + n + 1; 
            locD(x, y) = (D1(5, 4)*cm*DD(x, y + 1) &
                          + dm*D1(5, 5)*DD(x, y) &
                          + dmm*D1(5, 3)*DD(x, y + 2)); 
         end do
         DD(:, :) = dcmplx(0.0, 0.0)
         en = 2*(n + 2) - 1

         DD(3:en, 3:en) = locD

         ind1 = n**2; 
         ind2 = (n + 1)**2 - 1; 
         D(ind1:ind2, ind1:ind2) = locD

         deallocate (locD)
      end do

!______________Normalization___________________!

      las = 0
      do n = 1, Nmax
         do m = -n, n
            do mu = -n, n
               las = las + 1
               ind1 = n*n + n + m
               ind2 = n*n + n + mu

               if (m >= 0 .and. mu >= 0) then
                  delta = 1
               end if

               if (m >= 0 .and. mu < 0) then
                  delta = (-1)**mu
               end if

               if (m < 0 .and. mu >= 0) then
                  delta = (-1)**m
               end if

               if (m < 0 .and. mu < 0) then
                  delta = (-1)**(m + mu)
               end if

               spD(las) = delta*D(ind1, ind2)
               ind(las, 1) = ind1
               ind(las, 2) = ind2
            end do
         end do
      end do

   end subroutine sph_rotation_sparse

!****************************************************************************80
! Rotation matrix (recursive) for VSWFs
! Algorithm:  Choi, et al., J. Chem. Phys. vol. 11, No. 19, 1999.
   subroutine sph_rotation_sparse2(theta, phi2, Nmax1, Nmax2, spD1, spD2, ind_1, ind_2)
      real(dp) :: theta, phi, phi2
      integer :: Nmax, Nmax1, Nmax2
      complex(dp), dimension(:) :: spD1, spD2
      integer, dimension(:, :) :: ind_1, ind_2

      integer :: n, m, mu, ind1, ind2, mm, x, y, en, las, size_rot, m2, mu2, en2
      real(dp) :: ct, st, cp, sp, R(3, 3)
      complex(dp) :: C(3, 3), invC(3, 3), D1(7, 7)
      complex(dp), dimension(:, :), allocatable :: DD, locD
      complex(dp) :: a1, a2, a3
      real(dp) ::  a, b, b2, cm, dm, dmm, delta, nom_a, nom_b, nom_b2
      real(dp), allocatable, dimension(:) :: denom_mu

      if (Nmax1 > Nmax2) then
         Nmax = Nmax1
         size_rot = size(spD1)
      else
         Nmax = Nmax2
         size_rot = size(spD2)
      end if

      allocate (DD(2*(Nmax + 2) + 1, 2*(Nmax + 2) + 1))
      allocate (locD(2*Nmax + 1, 2*Nmax + 1))
      allocate (denom_mu(2*Nmax - 1))

      phi = -1.0d0*phi2 !

      ct = cos(theta)
      st = sin(theta)
      cp = cos(phi)
      sp = sin(phi)

! YZ-rotation matrix
      R(1, :) = [ct*cp, -sp, st*cp]
      R(2, :) = [ct*sp, cp, st*sp]
      R(3, :) = [-st, dble(0.0), ct]

      a1 = 1.0d0/sqrt(2.0d0)
      a2 = dcmplx(0.0, 1.0)*(1.0d0)/sqrt(2.0d0)
      a3 = 1.0d0

      C(1, :) = [a1, dcmplx(0.0, 0.0), -a1]
      C(2, :) = [-a2, dcmplx(0.0, 0.0), -a2]
      C(3, :) = [0.0d0, 1.0d0, 0.0d0]

      invC(1, :) = [a1, a2, 0.0*a3]
      invC(2, :) = [0.0*a3, 0.0*a3, a3]
      invC(3, :) = [-a1, a2, 0.0*a3]

      D1(:, :) = dcmplx(0.0, 0.0)
      D1(3:5, 3:5) = (transpose(matmul(invC, matmul(R, C))))

      DD = dcmplx(0.0, 0.0)
      DD(1:7, 1:7) = D1

      las = 0
!********** n = 1 *************************************!!
      n = 1
      do m2 = -n, n
         do mu2 = -n, n
            las = las + 1

            ind1 = n*n + n + m2
            ind2 = n*n + n + mu2

            if (m2 >= 0 .and. mu2 >= 0) then
               delta = 1
            end if

            if (m2 >= 0 .and. mu2 < 0) then
               delta = (-1)**mu2
            end if

            if (m2 < 0 .and. mu2 >= 0) then
               delta = (-1)**m2
            end if

            if (m2 < 0 .and. mu2 < 0) then
               delta = (-1)**(m2 + mu2)
            end if

            if (Nmax1 >= n) then
               spD1(las) = delta*D1(2 + m2 + n + 1, 2 + mu2 + n + 1)
               ind_1(las, 1) = ind1
               ind_1(las, 2) = ind2
            end if

            if (Nmax2 >= n) then
               spD2(las) = delta*D1(2 + m2 + n + 1, 2 + mu2 + n + 1)
               ind_2(las, 1) = ind1
               ind_2(las, 2) = ind2
            end if

         end do
      end do

!****************************************************

      do n = 2, Nmax
         locD(:, :) = dcmplx(0.0, 0.0)
         do mu = -n + 1, n - 1
            denom_mu(mu + n) = sqrt(dble((n + mu)*(n - mu)))
         end do

         do m = -n, n
            mm = -m
            nom_a = sqrt(dble((n + m)*(n - m)))
            nom_b = sqrt(dble((n + m)*(n + m - 1)))
            nom_b2 = sqrt(dble((n + mm)*(n + mm - 1)))

            do mu = -n + 1, n - 1
               a = nom_a/denom_mu(mu + n)
               b = nom_b/(sqrt(2.0d0)*denom_mu(mu + n))
               b2 = nom_b2/(sqrt(2.0d0)*denom_mu(mu + n))

               x = mu + n + 1
               y = m + n + 1

               locD(x, y) = D1(4, 4)*a*DD(x + 1, y + 1) &
                            + b*D1(4, 5)*DD(x + 1, y) &
                            + b2*D1(4, 3)*DD(x + 1, y + 2)
            end do

            mu = -n

            cm = sqrt(dble(n + m)*(n - m)/dble(n*(2*n - 1))); 
            dm = sqrt(dble(n + m)*(n + m - 1)/dble(2*n*(2*n - 1))); 
            dmm = sqrt(dble(n - m)*(n - m - 1)/dble(2*n*(2*n - 1))); 
            x = mu + n + 1; 
            y = m + n + 1; 
            locD(x, y) = (D1(3, 4)*cm*DD(x + 2, y + 1) &
                          + dm*D1(3, 5)*DD(x + 2, y) &
                          + dmm*D1(3, 3)*DD(x + 2, y + 2))

            mu = n; 
            x = mu + n + 1; 
            y = m + n + 1; 
            locD(x, y) = (D1(5, 4)*cm*DD(x, y + 1) &
                          + dm*D1(5, 5)*DD(x, y) &
                          + dmm*D1(5, 3)*DD(x, y + 2)); 
         end do

         DD(:, :) = dcmplx(0.0, 0.0)
         en = 2*(n + 2) - 1
         en2 = 2*n + 1
         DD(3:en, 3:en) = locD(1:en2, 1:en2)

         ind1 = n**2; 
         ind2 = (n + 1)**2 - 1; 
         do m2 = -n, n
            do mu2 = -n, n
               las = las + 1

               ind1 = n*n + n + m2
               ind2 = n*n + n + mu2

               if (m2 >= 0 .and. mu2 >= 0) then
                  delta = 1
               end if

               if (m2 >= 0 .and. mu2 < 0) then
                  delta = (-1)**mu2
               end if

               if (m2 < 0 .and. mu2 >= 0) then
                  delta = (-1)**m2

               end if

               if (m2 < 0 .and. mu2 < 0) then
                  delta = (-1)**(m2 + mu2)

               end if

               if (Nmax1 >= n) then
                  spD1(las) = delta*locD(m2 + n + 1, mu2 + n + 1)
                  ind_1(las, 1) = ind1
                  ind_1(las, 2) = ind2
               end if

               if (Nmax2 >= n) then
                  spD2(las) = delta*locD(m2 + n + 1, mu2 + n + 1)
                  ind_2(las, 1) = ind1
                  ind_2(las, 2) = ind2
               end if
            end do
         end do
      end do

   end subroutine sph_rotation_sparse2

!****************************************************************************80

   subroutine sph_rotation_sparse_gen(angles, Nmax, spD, ind)
      real(dp) :: rot(3, 3)
      real(dp) :: angles(3) !Euler angles [alpha,beta,gamma]
      integer :: Nmax, Nmax1, Nmax2
      complex(dp), dimension(:) :: spD
      integer, dimension(:, :) :: ind

      integer :: n, m, mu, ind1, ind2, mm, x, y, en, las, size_rot, m2, mu2, en2
      real(dp) :: ct, st, cp, sp, R(3, 3)
      complex(dp) :: C(3, 3), invC(3, 3), D1(7, 7)
      complex(dp), dimension(:, :), allocatable :: DD, locD
      complex(dp) :: a1, a2, a3
      real(dp) ::  a, b, b2, cm, dm, dmm, delta, nom_a, nom_b, nom_b2
      real(dp), allocatable, dimension(:) :: denom_mu

! Note Change of signs alpha and gamma in the rotation matrix
      R = rotation_matrix(-angles(1), angles(2), -angles(3))

      allocate (DD(2*(Nmax + 2) + 1, 2*(Nmax + 2) + 1))
      allocate (locD(2*Nmax + 1, 2*Nmax + 1))
      allocate (denom_mu(2*Nmax - 1))

      a1 = 1.0d0/sqrt(2.0d0)
      a2 = dcmplx(0.0, 1.0)*(1.0d0)/sqrt(2.0d0)
      a3 = 1.0d0

      C(1, :) = [a1, dcmplx(0.0, 0.0), -a1]
      C(2, :) = [-a2, dcmplx(0.0, 0.0), -a2]
      C(3, :) = [0.0d0, 1.0d0, 0.0d0]

      invC(1, :) = [a1, a2, 0.0*a3]
      invC(2, :) = [0.0*a3, 0.0*a3, a3]
      invC(3, :) = [-a1, a2, 0.0*a3]

      D1(:, :) = dcmplx(0.0, 0.0)
      D1(3:5, 3:5) = (transpose(matmul(invC, matmul(R, C))))

      DD = dcmplx(0.0, 0.0)
      DD(1:7, 1:7) = D1

      las = 0
!********** n = 1 *************************************!!
      n = 1
      do m2 = -n, n
         do mu2 = -n, n
            las = las + 1

            ind1 = n*n + n + m2
            ind2 = n*n + n + mu2

            if (m2 >= 0 .and. mu2 >= 0) then
               delta = 1
            end if

            if (m2 >= 0 .and. mu2 < 0) then
               delta = (-1)**mu2
            end if

            if (m2 < 0 .and. mu2 >= 0) then
               delta = (-1)**m2
            end if

            if (m2 < 0 .and. mu2 < 0) then
               delta = (-1)**(m2 + mu2)
            end if

            spD(las) = delta*D1(2 + m2 + n + 1, 2 + mu2 + n + 1)
            ind(las, 1) = ind1
            ind(las, 2) = ind2

         end do
      end do

!****************************************************

      do n = 2, Nmax
         locD(:, :) = dcmplx(0.0, 0.0)
         do mu = -n + 1, n - 1
            denom_mu(mu + n) = sqrt(dble((n + mu)*(n - mu)))
         end do

         do m = -n, n
            mm = -m

            nom_a = sqrt(dble((n + m)*(n - m)))
            nom_b = sqrt(dble((n + m)*(n + m - 1)))
            nom_b2 = sqrt(dble((n + mm)*(n + mm - 1)))

            do mu = -n + 1, n - 1

               a = nom_a/denom_mu(mu + n)
               b = nom_b/(sqrt(2.0d0)*denom_mu(mu + n))
               b2 = nom_b2/(sqrt(2.0d0)*denom_mu(mu + n))

               x = mu + n + 1
               y = m + n + 1

               locD(x, y) = D1(4, 4)*a*DD(x + 1, y + 1) &
                            + b*D1(4, 5)*DD(x + 1, y) &
                            + b2*D1(4, 3)*DD(x + 1, y + 2)
            end do

            mu = -n

            cm = sqrt(dble(n + m)*(n - m)/dble(n*(2*n - 1))); 
            dm = sqrt(dble(n + m)*(n + m - 1)/dble(2*n*(2*n - 1))); 
            dmm = sqrt(dble(n - m)*(n - m - 1)/dble(2*n*(2*n - 1))); 
            x = mu + n + 1; 
            y = m + n + 1; 
            locD(x, y) = (D1(3, 4)*cm*DD(x + 2, y + 1) &
                          + dm*D1(3, 5)*DD(x + 2, y) &
                          + dmm*D1(3, 3)*DD(x + 2, y + 2))

            mu = n; 
            x = mu + n + 1; 
            y = m + n + 1; 
            locD(x, y) = (D1(5, 4)*cm*DD(x, y + 1) &
                          + dm*D1(5, 5)*DD(x, y) &
                          + dmm*D1(5, 3)*DD(x, y + 2)); 
         end do

         DD(:, :) = dcmplx(0.0, 0.0)
         en = 2*(n + 2) - 1
         en2 = 2*n + 1
         DD(3:en, 3:en) = locD(1:en2, 1:en2)

         ind1 = n**2; 
         ind2 = (n + 1)**2 - 1; 
         do m2 = -n, n
            do mu2 = -n, n
               las = las + 1

               ind1 = n*n + n + m2
               ind2 = n*n + n + mu2

               if (m2 >= 0 .and. mu2 >= 0) then
                  delta = 1
               end if

               if (m2 >= 0 .and. mu2 < 0) then
                  delta = (-1)**mu2
               end if

               if (m2 < 0 .and. mu2 >= 0) then
                  delta = (-1)**m2

               end if

               if (m2 < 0 .and. mu2 < 0) then
                  delta = (-1)**(m2 + mu2)

               end if

               spD(las) = delta*locD(m2 + n + 1, mu2 + n + 1)
               ind(las, 1) = ind1
               ind(las, 2) = ind2
            end do
         end do
      end do

   end subroutine sph_rotation_sparse_gen

!****************************************************************************80
! Sparse matmul
   function sparse_matmul(val, ind, x, nm_out) result(Ax)
      complex(dp), dimension(:) :: val, x
      integer, dimension(:, :) :: ind
      integer :: nm_out, i1
      complex(dp), dimension(nm_out) :: Ax

      Ax(:) = dcmplx(0.0, 0.0)
      do i1 = 1, size(ind, 1)
         Ax(ind(i1, 1)) = Ax(ind(i1, 1)) + val(i1)*x(ind(i1, 2))
      end do

   end function sparse_matmul

!****************************************************************************80
! Sparse matmul congj transpose
   function sparse_matmul_H(val, ind, x, nm_out) result(Ax)
      complex(dp), dimension(:) :: val, x
      integer, dimension(:, :) :: ind
      integer :: nm_out, i1
      complex(dp), dimension(nm_out) :: Ax

      Ax(:) = dcmplx(0.0, 0.0)
      do i1 = 1, size(ind, 1)
         Ax(ind(i1, 2)) = Ax(ind(i1, 2)) + conjg(val(i1))*x(ind(i1, 1))
      end do

   end function sparse_matmul_H

!****************************************************************************80
! translation inout = 0 (in) inout = 1 (out)
   subroutine translate(tr, N_in, N_out, k, a_in, b_in, a_out, b_out, inout)
      real(dp) :: tr(3), dir(3)
      integer :: N_in, N_out
      complex(dp), dimension((N_in + 1)**2 - 1) :: a_in, b_in
      complex(dp), dimension((N_out + 1)**2 - 1) :: a_out, b_out
      complex(dp), dimension(:), allocatable :: A, B
      integer, dimension(:, :), allocatable :: indx
      complex(dp) :: k
      integer :: inout, m, n, nu, las, nm_out, nm_in, las2
      complex(dp), dimension((N_in + 1)*(2*N_in + 1)*(2*N_in + 3)/3 - 1) :: rot_in
      integer, dimension((N_in + 1)*(2*N_in + 1)*(2*N_in + 3)/3 - 1, 2) :: ind_rot_in
      complex(dp), dimension((N_out + 1)*(2*N_out + 1)*(2*N_out + 3)/3 - 1) :: rot_out
      integer, dimension((N_out + 1)*(2*N_out + 1)*(2*N_out + 3)/3 - 1, 2) :: ind_rot_out

      las2 = 0
      do n = 1, N_out
         do nu = 1, N_in
            do m = -n, n
               if (abs(m) <= n .and. abs(m) <= nu) then
                  las2 = las2 + 1
               end if
            end do
         end do
      end do

      allocate (A(las2), B(las2))
      allocate (indx(las2, 2))

      nm_out = (N_out + 1)**2 - 1
      nm_in = (N_in + 1)**2 - 1

      dir = cart2sph(tr) ! r, theta, phi

      if (dir(1) > 0.0d0) then
! Translation coefficients along the z-axis
         call vector_translator_sparse(dir(1), k, N_in, N_out, A, B, indx, inout)

         call sph_rotation_sparse2(dir(2), dir(3), N_in, N_out, rot_in, rot_out, ind_rot_in, ind_rot_out)

         ! Translate

         a_out = sparse_matmul_H(rot_out, ind_rot_out, &
                                 (sparse_matmul(A, indx, sparse_matmul(rot_in, ind_rot_in, a_in, nm_in), nm_out) + &
                                  sparse_matmul(B, indx, sparse_matmul(rot_in, ind_rot_in, b_in, nm_in), nm_out)), nm_out)

         b_out = sparse_matmul_H(rot_out, ind_rot_out, &
                                 (sparse_matmul(A, indx, sparse_matmul(rot_in, ind_rot_in, b_in, nm_in), nm_out) + &
                                  sparse_matmul(B, indx, sparse_matmul(rot_in, ind_rot_in, a_in, nm_in), nm_out)), nm_out)
      else
         if (size(a_in) > size(a_out)) then
            a_out = a_in(1:size(a_out))
            b_out = b_in(1:size(b_out))
         else
            a_out = dcmplx(0.0, 0.0)
            a_out(1:size(a_in)) = a_in
            b_out = dcmplx(0.0, 0.0)
            b_out(1:size(b_in)) = b_in
         end if
      end if

   end subroutine translate

end module
