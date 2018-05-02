module cluster
   use common
   use mie
   use sfunctions
   use omp_lib
   implicit none

contains

!*****************************************************************
!
! Calculates far fields from coefficients a_nm, b_nm
!
! Eth = theta component of the electric field * R / exp(ikR)
! Ephi = phi component of the electric field * R / exp(ikR)
!
! Note: This cannot evaluate fields at the z-axis
!
!*******************************************************************

   subroutine calc_farfields(a_nm, b_nm, k, r, theta, phi, F)
      complex(dp), dimension(:) :: a_nm, b_nm
      real(dp) :: Po(3)
      complex(dp) :: k, F(3), G(3), Ephi, Eth, i
      integer :: Nmax, n, m, ind, mm
      real(dp) :: r, theta, phi, vec(3), q
      complex(dp) :: kr, alpha, beta, gamma, cc
      complex(dp), dimension(:), allocatable :: sphj, sphy, sphh
      complex(dp), dimension(:), allocatable :: chf1, chd1, chf2, chd2
      complex(dp) :: P(3), B(3), C(3), Y, Y1, Y2, M_nm(3), N_nm(3)
      real(dp), dimension(:), allocatable :: L, L1, L2

      i = dcmplx(0.0, 1.0)
      Nmax = int(sqrt(dble(1 + size(a_nm)))) - 1

      ind = 0
      F(:) = dcmplx(0.0, 0.0)
      G(:) = dcmplx(0.0, 0.0)

      sphh = sphankel_a(Nmax + 1, k*r)

      do n = 1, Nmax

         !alpha = (-i)**(n+1) * exp(i*k*r)/(k*r)!sphh(n+1)
         !gamma = -(-i)**(n+2) * exp(i*k*r)/(k*r) !(n+1.0d0)/kr * sphh(n+1) - sphh(n+2)

         alpha = sphh(n + 1)
         gamma = -sphh(n + 2)

         allocate (L(n + 1), L1(n + 2), L2(n))

         call legendre2(n, cos(theta), L)
         call legendre2(n + 1, cos(theta), L1)
         call legendre2(n - 1, cos(theta), L2)

         q = (sqrt(n*(n + 1.0d0)))/((n*2d0 + 1.0d0)*sin(theta)); 
         do m = -n, n
            ind = ind + 1
            mm = abs(m)

            cc = sqrt((2d0*n + 1.0d0)*factorial(n - mm)/factorial(n + mm)/(4d0*pi)); 
            ! Unnormalized complex scalar spherical harmonics
            Y = L(mm + 1)*exp(dcmplx(0.0, m*phi)); 
            Y1 = L1(mm + 1)*exp(dcmplx(0.0, m*phi)); 
            if (mm == n) then
               Y2 = dcmplx(0.0, 0.0)
            else
               Y2 = L2(mm + 1)*exp(dcmplx(0.0, m*phi))
            end if

            ! transversal vector spherical harmonics

            Y1 = Y1*((n - mm + 1.0d0)/(n + 1.0d0))
            Y2 = Y2*(dble(n + mm)/dble(n))

            B(:) = dcmplx(0.0, 0.0)
            B(2) = Y1 - Y2
            B(3) = ((dcmplx(0.0, m*(2*n + 1.0)))/(n*(n + 1.0)))*Y

            B = B*q

            C(:) = dcmplx(0.0, 0.0)
            C(2) = B(3)
            C(3) = -B(2)

            ! Spherical vector wave functions
            M_nm = cc*alpha*C
            N_nm = cc*(gamma*B)

            F = F + a_nm(ind)*M_nm + b_nm(ind)*N_nm

         end do

         deallocate (L, L1, L2)

      end do

      Eth = F(2)
      Ephi = F(3)

      F = sph2cart_vec(theta, phi, F)

!G = G/(dcmplx(0.0,1.0) * omega * mu)

   end subroutine calc_farfields

!**********************************************************************
!
! Computes farfields and cross sections for a cluster
!
!**********************************************************************
   subroutine compute_cluster_fields(sphere, x, rhs, k, N_theta, N_phi, E_out, crs)
      type(data_struct), dimension(:) :: sphere
      complex(dp), dimension(:) :: x, rhs
      complex(dp) :: k
      complex(dp), dimension(:, :), allocatable :: E_out

      integer :: sph, a1, a2, b1, b2, Nmax, N_phi, N_theta, i1, i2, las
      real(dp) :: phi, theta, r(3), RR, cp(3), Cext, Csca, Cabs
      real(dp) :: Cext2, Csca2, Cabs2, rad, absb, rhat(3), Csca_int, crs(4)
      complex(dp) :: F(3), G(3), epsr
      complex(dp), dimension(:), allocatable :: c_in, d_in

      allocate (E_out(3, N_theta*N_phi))
      E_out(:, :) = dcmplx(0.0, 0.0)

      Cext = 0d0
      Csca = 0d0
      Cabs = 0d0
      Csca_int = 0d0

!$OMP parallel default(private) &
!$OMP firstprivate(N_phi, N_theta, k) &
!$OMP shared(Cabs, Cext, sphere,x, E_out,rhs)
!$OMP do reduction(+:E_out,Cabs,Cext)

      do sph = 1, size(sphere)

         a1 = sphere(sph)%ind_a1
         a2 = sphere(sph)%ind_a2
         b1 = sphere(sph)%ind_b1
         b2 = sphere(sph)%ind_b2
         cp = sphere(sph)%cp
         rad = sphere(sph)%r
         epsr = sphere(sph)%eps_r

         Nmax = sphere(sph)%Nmax
         RR = 1.0d6
         las = 1
         do i1 = 1, N_phi
            do i2 = 1, N_theta

               phi = 2*pi*(i1 - 1)/N_phi; 
               theta = pi*(i2 - 1)/(N_theta) + pi/N_theta/2; 
               r(1) = RR*sin(theta)*cos(phi); 
               r(2) = RR*sin(theta)*sin(phi); 
               r(3) = RR*cos(theta); 
               ! rhat = r/sqrt(dot_product(r,r))

               r = r - cp

               call calc_fields(x(a1:a2), x(b1:b2), k, r, F, G, 1)

               E_out(:, las) = E_out(:, las) + F

               !Csca_int = Csca_int + pi/2*4*pi*RR**2*dot_product(F,F) &
               !     *sin(theta)/N_phi/N_theta
               las = las + 1
            end do
         end do

         call cross_sections(x(a1:a2), x(b1:b2), rhs(a1:a2), rhs(b1:b2), &
                             k, Nmax, Cext2, Csca2, absb)

         call sphere_absorbtion2(x(a1:a2), x(b1:b2), dble(k), sphere(sph)%r, sphere(sph)%eps_r, Nmax, Cabs2)

         Cext = Cext + Cext2
         Cabs = Cabs + Cabs2

      end do
!$OMP end do
!$OMP end parallel
      Csca_int = 0d0
      las = 1
      do i1 = 1, N_phi
         do i2 = 1, N_theta
            phi = 2*pi*(i1 - 1)/N_phi; 
            theta = pi*(i2 - 1)/(N_theta) + pi/N_theta/2; 
            RR = 1.0d6
            Csca_int = Csca_int + pi/2.0*4.0*pi*RR**2* &
                       dot_product(E_out(:, las), E_out(:, las))* &
                       sin(theta)/N_phi/N_theta

            las = las + 1
         end do
      end do

      Csca = Cext - Cabs

      crs = [Cext, Csca, Cabs, Csca_int]
      print *, 'Cext', Cext
      print *, 'Csca', Csca
      print *, 'Csca_int', Csca_int
      print *, 'Cabs', Cabs

   end subroutine compute_cluster_fields

!**********************************************************************
!
! Computes farfields and cross sections for a cluster
!
!**********************************************************************
   subroutine compute_cluster_farfields(sphere, x, rhs, k, N_theta, N_phi, E_out, crs)
      type(data_struct), dimension(:) :: sphere
      complex(dp), dimension(:) :: x, rhs
      complex(dp) :: k
      complex(dp), dimension(:, :), allocatable :: E_out

      integer :: sph, a1, a2, b1, b2, Nmax, N_phi, N_theta, i1, i2, las
      real(dp) :: phi, theta, r(3), RR, cp(3), Cext, Csca, Cabs
      real(dp) :: Cext2, Csca2, Cabs2, rad, absb, rhat(3), Csca_int, crs(4)
      complex(dp) :: F(3), G(3), epsr, Eth, Ephi
      complex(dp), dimension(:), allocatable :: c_in, d_in

      allocate (E_out(3, N_theta*N_phi))
      E_out(:, :) = dcmplx(0.0, 0.0)

      Cext = 0.0
      Csca = 0.0
      Cabs = 0.0
      Csca_int = 0.0

!$OMP parallel default(private) &
!$OMP firstprivate(N_phi, N_theta, k) &
!$OMP shared(Cabs, Cext, sphere,x, E_out,rhs)
!$OMP do reduction(+:E_out,Cabs,Cext)

      do sph = 1, size(sphere)

         a1 = sphere(sph)%ind_a1
         a2 = sphere(sph)%ind_a2
         b1 = sphere(sph)%ind_b1
         b2 = sphere(sph)%ind_b2
         cp = sphere(sph)%cp
         rad = sphere(sph)%r
         epsr = sphere(sph)%eps_r

         Nmax = sphere(sph)%Nmax

         las = 1
         do i1 = 1, N_phi
            do i2 = 1, N_theta

               phi = 2*pi*(i1 - 1)/N_phi; 
               theta = pi*(i2 - 1)/(N_theta) + pi/N_theta/2; 
               RR = 1e6
               r(1) = RR*sin(theta)*cos(phi); 
               r(2) = RR*sin(theta)*sin(phi); 
               r(3) = RR*cos(theta); 
               ! rhat = r/sqrt(dot_product(r,r))

               r = r - cp

               !call calc_fields(x(a1:a2), x(b1:b2), k, r, F, G, 1)
               !print*, F / exp(dcmplx(0.0,1.0)*k*RR) * (k*RR)
               call calc_farfields(x(a1:a2), x(b1:b2), k, norm(r, cp), theta, phi, F)
               !print*, F
               !stop

               E_out(:, las) = E_out(:, las) + F/exp(dcmplx(0.0, 1.0)*k*RR)*(k*RR)
               !E_out(3,las) = E_out(3,las) + Ephi

               !Csca_int = Csca_int + pi/2*4*pi*RR**2*dot_product(F,F) &
               !     *sin(theta)/N_phi/N_theta
               las = las + 1
            end do
         end do

         call cross_sections(x(a1:a2), x(b1:b2), rhs(a1:a2), rhs(b1:b2), &
                             k, Nmax, Cext2, Csca2, absb)

         call sphere_absorbtion2(x(a1:a2), x(b1:b2), dble(k), sphere(sph)%r, sphere(sph)%eps_r, Nmax, Cabs2)

         Cext = Cext + Cext2
         Cabs = Cabs + Cabs2

      end do
!$OMP end do
!$OMP end parallel
      Csca_int = 0.0
      las = 1
      do i1 = 1, N_phi
         do i2 = 1, N_theta
            phi = 2*pi*(i1 - 1)/N_phi; 
            theta = pi*(i2 - 1)/(N_theta) + pi/N_theta/2; 
            Csca_int = Csca_int + pi/2.0*4.0*pi* &
                       dot_product(E_out(:, las), E_out(:, las))* &
                       sin(theta)/N_phi/N_theta

            las = las + 1
         end do
      end do

      Csca = Cext - Cabs

      crs = [Cext, Csca, Cabs, Csca_int]
      print *, 'Cext', Cext
      print *, 'Csca', Csca
      print *, 'Csca_int', Csca_int
      print *, 'Cabs', Cabs

   end subroutine compute_cluster_farfields

   subroutine mueller_matrix_far(sphere, x, x2, b_vec, b_vec2, k, N_theta, N_phi, Nmax, S_out, crs)
      type(data_struct), dimension(:) :: sphere
      complex(dp), dimension(:) :: x, x2
      complex(dp) :: k
      real(dp), dimension(:, :), allocatable :: S_out
      integer :: N_phi, N_theta, Nmax
      complex(dp), dimension(:) :: b_vec, b_vec2
      complex(dp), dimension(:, :), allocatable :: E_out, E_out2

      integer :: i1, i2, las
      real(dp) :: theta, phi, abcd(2, 2), r(3), RR, unit_th(3), unit_phi(3), crs(4)
      real(dp) :: crs1(4), crs2(4)
      complex(dp), dimension(N_phi*N_theta) :: f11, f12, f21, f22
      complex(dp) ::S1, S2, S3, S4
      complex(dp) :: i

      i = dcmplx(0.0, 1.0)

      allocate (S_out(N_phi*N_theta, 18))

!call inc_xy(Nmax, sphere, b_vec, b_vec2, k)
      call compute_cluster_farfields(sphere, x, b_vec, k, N_theta, N_phi, E_out, crs1)
      call compute_cluster_farfields(sphere, x2, b_vec2, k, N_theta, N_phi, E_out2, crs2)

      crs = (crs1 + crs2)/2

      las = 0
      do i1 = 1, N_phi
         do i2 = 1, N_theta

            theta = pi*(i2 - 1)/(N_theta) + pi/N_theta/2.0
            phi = 2*pi*(i1 - 1)/N_phi

            abcd(1, :) = [cos(phi), sin(phi)]
            abcd(2, :) = [sin(phi), -cos(phi)]

            unit_th = [cos(theta)*cos(phi), sin(phi)*cos(theta), -sin(theta)]
            unit_phi = [-sin(phi), cos(phi), 0.0d0]

            !f11(las+1) = (E_out(2,las+1))
            !f21(las+1) = (E_out(3,las+1))
            !f12(las+1) = (E_out2(2,las+1))
            !f22(las+1) = (E_out2(3,las+1))

            f11(las + 1) = dot_product(E_out(:, las + 1), unit_th); 
            f21(las + 1) = dot_product(E_out(:, las + 1), unit_phi); 
            f12(las + 1) = dot_product(E_out2(:, las + 1), unit_th); 
            f22(las + 1) = dot_product(E_out2(:, las + 1), unit_phi); 
            S1 = -i*(f21(las + 1)*abcd(1, 2) + f22(las + 1)*abcd(2, 2))
            S2 = -i*(f11(las + 1)*abcd(1, 1) + f12(las + 1)*abcd(2, 1))
            S3 = i*(f11(las + 1)*abcd(1, 2) + f12(las + 1)*abcd(2, 2))
            S4 = i*(f21(las + 1)*abcd(1, 1) + f22(las + 1)*abcd(2, 1))

            S_out(las + 1, 1) = phi
            S_out(las + 1, 2) = theta

            ! Mueller matrix
            S_out(las + 1, 3) = (abs(S1)**2 + abs(S2)**2 + abs(S3)**2 + abs(S4)**2)/2.0 !S11
            S_out(las + 1, 4) = (abs(S2)**2 - abs(S1)**2 + abs(S4)**2 - abs(S3)**2)/2.0 !S12
            S_out(las + 1, 5) = -real(S2*conjg(S3) + S1*conjg(S4)) !S13
            S_out(las + 1, 6) = -imag(S2*conjg(S3) - S1*conjg(S4)) !S14
            S_out(las + 1, 7) = (abs(S2)**2 - abs(S1)**2 + abs(S3)**2 - abs(S4)**2)/2.0 !S21
            S_out(las + 1, 8) = (abs(S1)**2 - abs(S3)**2 - abs(S4)**2 + abs(S2)**2)/2.0 !S22
            S_out(las + 1, 9) = real(S2*conjg(S3) - S1*conjg(S4)) !S23
            S_out(las + 1, 10) = -imag(S2*conjg(S3) + S1*conjg(S4)) !S24
            S_out(las + 1, 11) = real(S2*conjg(S4) + S1*conjg(S3)) !S31
            S_out(las + 1, 12) = -real(S2*conjg(S4) - S1*conjg(S3)) !S32
            S_out(las + 1, 13) = -real(S1*conjg(S2) + S3*conjg(S4)) !S33
            S_out(las + 1, 14) = imag(S2*conjg(S1) + S4*conjg(S3)) ! S34
            S_out(las + 1, 15) = imag(S4*conjg(S2) + S1*conjg(S3)) ! S41
            S_out(las + 1, 16) = imag(S4*conjg(S2) - S1*conjg(S3)) ! S42
            S_out(las + 1, 17) = imag(S1*conjg(S2) - S3*conjg(S4)) ! S43
            S_out(las + 1, 18) = -real(S1*conjg(S2) - S3*conjg(S4)) ! S44

            las = las + 1
         end do
      end do

   end subroutine mueller_matrix_far

   subroutine mueller_matrix(sphere, x, x2, b_vec, b_vec2, k, N_theta, N_phi, Nmax, S_out, crs)
      type(data_struct), dimension(:) :: sphere
      complex(dp), dimension(:) :: x, x2
      complex(dp) :: k
      real(dp), dimension(:, :), allocatable :: S_out
      integer :: N_phi, N_theta, Nmax
      complex(dp), dimension(:) :: b_vec, b_vec2
      complex(dp), dimension(:, :), allocatable :: E_out, E_out2

      integer :: i1, i2, las
      real(dp) :: theta, phi, abcd(2, 2), r(3), RR, unit_th(3), unit_phi(3), crs(4)
      real(dp) :: crs1(4), crs2(4)
      complex(dp), dimension(N_phi*N_theta) :: f11, f12, f21, f22
      complex(dp) ::S1, S2, S3, S4
      complex(dp) :: i

      i = dcmplx(0.0, 1.0)
      RR = 1.0d6
      allocate (S_out(N_phi*N_theta, 18))

!call inc_xy(Nmax, sphere, b_vec, b_vec2, k)
      call compute_cluster_fields(sphere, x, b_vec, k, N_theta, N_phi, E_out, crs1)
      call compute_cluster_fields(sphere, x2, b_vec2, k, N_theta, N_phi, E_out2, crs2)

      crs = (crs1 + crs2)/2

      las = 0
      do i1 = 1, N_phi
         do i2 = 1, N_theta

            theta = pi*(i2 - 1)/(N_theta) + pi/N_theta/2.0
            phi = 2*pi*(i1 - 1)/N_phi

            abcd(1, :) = [cos(phi), sin(phi)]
            abcd(2, :) = [sin(phi), -cos(phi)]

            r(1) = RR*sin(theta)*cos(phi)
            r(2) = RR*sin(theta)*sin(phi)
            r(3) = RR*cos(theta)

            unit_th = [cos(theta)*cos(phi), sin(phi)*cos(theta), -sin(theta)]
            unit_phi = [-sin(phi), cos(phi), 0.0d0]; 
            f11(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out(:, las + 1), unit_th); 
            f21(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out(:, las + 1), unit_phi); 
            f12(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out2(:, las + 1), unit_th); 
            f22(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out2(:, las + 1), unit_phi); 
            S1 = -i*(f21(las + 1)*abcd(1, 2) + f22(las + 1)*abcd(2, 2))
            S2 = -i*(f11(las + 1)*abcd(1, 1) + f12(las + 1)*abcd(2, 1))
            S3 = i*(f11(las + 1)*abcd(1, 2) + f12(las + 1)*abcd(2, 2))
            S4 = i*(f21(las + 1)*abcd(1, 1) + f22(las + 1)*abcd(2, 1))

            S_out(las + 1, 1) = phi
            S_out(las + 1, 2) = theta

            ! Mueller matrix
            S_out(las + 1, 3) = (abs(S1)**2 + abs(S2)**2 + abs(S3)**2 + abs(S4)**2)/2.0 !S11
            S_out(las + 1, 4) = (abs(S2)**2 - abs(S1)**2 + abs(S4)**2 - abs(S3)**2)/2.0 !S12
            S_out(las + 1, 5) = -real(S2*conjg(S3) + S1*conjg(S4)) !S13
            S_out(las + 1, 6) = -imag(S2*conjg(S3) - S1*conjg(S4)) !S14
            S_out(las + 1, 7) = (abs(S2)**2 - abs(S1)**2 + abs(S3)**2 - abs(S4)**2)/2.0 !S21
            S_out(las + 1, 8) = (abs(S1)**2 - abs(S3)**2 - abs(S4)**2 + abs(S2)**2)/2.0 !S22
            S_out(las + 1, 9) = real(S2*conjg(S3) - S1*conjg(S4)) !S23
            S_out(las + 1, 10) = -imag(S2*conjg(S3) + S1*conjg(S4)) !S24
            S_out(las + 1, 11) = real(S2*conjg(S4) + S1*conjg(S3)) !S31
            S_out(las + 1, 12) = -real(S2*conjg(S4) - S1*conjg(S3)) !S32
            S_out(las + 1, 13) = -real(S1*conjg(S2) + S3*conjg(S4)) !S33
            S_out(las + 1, 14) = imag(S2*conjg(S1) + S4*conjg(S3)) ! S34
            S_out(las + 1, 15) = imag(S4*conjg(S2) + S1*conjg(S3)) ! S41
            S_out(las + 1, 16) = imag(S4*conjg(S2) - S1*conjg(S3)) ! S42
            S_out(las + 1, 17) = imag(S1*conjg(S2) - S3*conjg(S4)) ! S43
            S_out(las + 1, 18) = -real(S1*conjg(S2) - S3*conjg(S4)) ! S44

            las = las + 1
         end do
      end do

   end subroutine mueller_matrix

!****************************************

end module cluster
