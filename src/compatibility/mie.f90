module mie
   use sfunctions
   use common
   use io
   use translations

   implicit none

contains

   subroutine sphere_absorbtion2(a_in, b_in, k, a, epsr, Nmax, absb)
      complex(dp), dimension(:) :: a_in, b_in
      real(dp) :: k, a, absb
      integer :: Nmax
      complex(dp) :: epsr

      complex(dp), dimension((Nmax + 1)**2 - 1) ::a_n, b_n, c_n, d_n
      complex(dp) :: dn, cn
      real(dp) ::  int, x
      integer :: n, m, las

      x = k*a

      call mie_coeff_nm(Nmax, x, sqrt(epsr), a_n, b_n, c_n, d_n)

      las = 1
      int = 0d0
      do n = 1, Nmax

         do m = -n, n

            cn = -dcmplx(1.0d0/a_n(las) + 1.0d0)
            dn = -dcmplx(1.0d0/b_n(las) + 1.0d0)

            int = int + 1/k**2d0*(dble(dn)*abs(a_in(las))**2d0 + &
               dble(cn)*abs(b_in(las))**2d0)

            las = las + 1
         end do

      end do

      absb = int

   end subroutine sphere_absorbtion2

!******************************************************************************

   subroutine cross_sections(a_out, b_out, a_in, b_in, k, Nmax, Cext, Csca, Cabs)
      real(dp) :: Cext, Cabs, Csca
      complex(dp) :: k
      complex(dp), dimension(:) :: a_out, a_in, b_out, b_in
      integer :: las, n, m, Nmax
      real(dp) :: ee

      Cext = 0d0
      Csca = 0d0
      las = 1
      do n = 1, Nmax
         do m = -n, n
            ee = dble(n)*(n + 1)/(2*n + 1)*factorial(n + m)/factorial(n - m)
            Csca = Csca + dble(1/k**2*(a_out(las)*conjg(a_out(las)) + b_out(las)*conjg(b_out(las))))
            Cext = Cext - real(1/k**2*(a_out(las)*conjg(a_in(las)) + b_out(las)*conjg(b_in(las))))

            las = las + 1
         end do
      end do

      Cabs = Cext - Csca

   end subroutine cross_sections

!******************************************************

   subroutine scattering_cross_section(a_out, b_out, k, Nmax, Csca)
      real(dp) :: Csca
      complex(dp) :: k
      complex(dp), dimension(:) :: a_out, b_out
      integer :: las, n, m, Nmax

      Csca = 0d0
      las = 1
      do n = 1, Nmax
         do m = -n, n

            Csca = Csca + dble(1/k**2*(a_out(las)*conjg(a_out(las)) + b_out(las)*conjg(b_out(las))))

            las = las + 1
         end do
      end do

   end subroutine scattering_cross_section

!******************************************************************************
! Bohren Huffman mie coefficients
   SUBROUTINE BHMIE(NSTOP, X, REFREL, a_n, b_n)

!Arguments:
      real(dp) :: X
      complex(dp) :: REFREL
! Local variables:
      integer :: N, NSTOP, NMX, NN
      real(dp) :: CHI, CHI0, CHI1, DX, EN, FN, P, PII, PSI, PSI0, PSI1,YMOD

      complex(dp) :: AN, AN1, BN, BN1, DREFRL, XI, XI1, Y
      complex(dp) :: D(NSTOP + 1), a_n(NSTOP), b_n(NSTOP)

!C*** Obtain pi:
      PII = 4.*ATAN(1.D0)
      DX = X
      DREFRL = REFREL
      Y = X*DREFRL
      YMOD = ABS(Y)
!
!*** Series expansion terminated after NSTOP terms
!    Logarithmic derivatives calculated from NMX on down

      NMX = NSTOP

!*** Logarithmic derivative D(J) calculated by downward recurrence
!    beginning with initial value (0.,0.) at J=NMX
!
      D(NMX) = (0.0d0, 0.0d0)
      NN = NMX - 1
      DO N = 1, NN
         EN = dble(NMX - N + 1)
         D(NMX - N) = (EN/Y) - (1./(D(NMX - N + 1) + EN/Y))
      end do
!
!*** Riccati-Bessel functions with real argument X
!    calculated by upward recurrence
!
      PSI0 = COS(DX)
      PSI1 = SIN(DX)
      CHI0 = -SIN(DX)
      CHI1 = COS(DX)
      XI1 = DCMPLX(PSI1, -CHI1)

      P = -1d0
      DO N = 1, NSTOP
         EN = dble(N)
         FN = (2.0d0*EN + 1.0d0)/(EN*(EN + 1.0d0))

         PSI = (2.0d0*EN - 1.0d0)*PSI1/DX - PSI0
         CHI = (2.0d0*EN - 1.0d0)*CHI1/DX - CHI0
         XI = DCMPLX(PSI, -CHI)

         IF (N .GT. 1) THEN
            AN1 = AN
            BN1 = BN
            a_n(N - 1) = -an1
            b_n(N - 1) = -bn1

         ENDIF

         AN = (D(N)/DREFRL + EN/DX)*PSI - PSI1
         AN = AN/((D(N)/DREFRL + EN/DX)*XI - XI1)
         BN = (DREFRL*D(N) + EN/DX)*PSI - PSI1
         BN = BN/((DREFRL*D(N) + EN/DX)*XI - XI1)

         PSI0 = PSI1
         PSI1 = PSI
         CHI0 = CHI1
         CHI1 = CHI
         XI1 = DCMPLX(PSI1, -CHI1)

      END DO

   end subroutine bhmie

!******************************************************************************

   subroutine mie_coeff_nm(N, x, mr, a_nm, b_nm, c_nm, d_nm)
      integer :: N
      real(dp) :: x
      complex(dp) :: mr
      complex(dp), dimension((N + 1)**2 - 1) :: a_nm, b_nm, c_nm, d_nm

      complex(dp), dimension(N) :: a_n, b_n
      integer :: las, i1, m

      call BHMIE(N + 1, x, mr, a_n, b_n)

      las = 0; 
      do i1 = 1, N
         do m = -i1, i1
            las = las + 1; 
            a_nm(las) = a_n(i1); 
            b_nm(las) = b_n(i1); 
            c_nm(las) = dcmplx(0d0) !c_n(i1);
            d_nm(las) = dcmplx(0d0) !d_n(i1);
         end do
      end do

   end subroutine mie_coeff_nm

!******************************************************************************
! The routine computes the SVWF expansion coefficients
! for a time harmonic x-polarized planewave
! propagating +z-direction  with the wave number k
   subroutine planewave(Nmax, k, a_nm, b_nm)
      integer :: Nmax
      complex(dp), dimension((Nmax + 1)**2 - 1) :: a_nm, b_nm
      complex(dp) :: k

      integer :: ind, n, m, mm
      real(dp) :: scale, C, E0, omega
      complex(dp) :: q

      E0 = 1d0
      omega = k*299792458.0
      ind = 0

      do n = 1, Nmax

         scale = sqrt(dble(n*(n + 1)))

         do m = -n, n
            ind = ind + 1
            mm = abs(m)

            C = scale*E0*sqrt(pi*(2*n + 1d0))/(n*(n + 1d0))*sqrt(factorial(n + 1)/factorial(n - 1))

            q = -(dcmplx(0.0, 1.0)**n*k)/dcmplx(omega*mu)

            if (mm == 1) then
               a_nm(ind) = dcmplx(0.0, 1.0)**dcmplx(n - 1.0)*dcmplx(C)
               b_nm(ind) = -dcmplx(0.0, 1.0)**dcmplx(n + 1.0)*dcmplx(C)

               if (m == -1) then
                  a_nm(ind) = -a_nm(ind)
               end if
            else
               a_nm(ind) = dcmplx(0.0, 0.0)
               b_nm(ind) = dcmplx(0.0, 0.0)
            end if
            a_nm(ind) = -a_nm(ind)
            b_nm(ind) = -b_nm(ind)
         end do

      end do

   end subroutine planewave

!******************************************************************************
! The routine computes the SVWF expansion coefficients
! for a time harmonic x-and y-polarized planewave
! propagating +z-direction  with the wave number k
   subroutine planewave2(Nmax, k, a_nm, b_nm, a_nm2, b_nm2)
      integer :: Nmax
      complex(dp), dimension((Nmax + 1)**2 - 1) :: a_nm, b_nm, a_nm2, b_nm2
      complex(dp) :: k

      integer :: ind, n, m, mm, las, nm_in
      real(dp) :: scale, C, E0, omega
      complex(dp) :: q

      complex(dp), dimension(:), allocatable :: rotD
      integer, dimension(:, :), allocatable :: indD

      E0 = 1
      omega = k*299792458.0
      ind = 0

      nm_in = (Nmax + 1)**2 - 1

      do n = 1, Nmax

         scale = sqrt(dble(n*(n + 1)))

         do m = -n, n
            ind = ind + 1
            mm = abs(m)

            C = scale*E0*sqrt(pi*(2*n + 1d0))/(n*(n + 1d0))*sqrt(factorial(n + 1)/factorial(n - 1))

            q = -(dcmplx(0.0, 1.0)**n*k)/dcmplx(omega*mu)

            if (mm == 1) then
               a_nm(ind) = dcmplx(0.0, 1.0)**dcmplx((n - 1.0))*dcmplx(C)
               b_nm(ind) = -dcmplx(0.0, 1.0)**dcmplx((n + 1.0))*dcmplx(C)

               if (m == -1) then
                  a_nm(ind) = -a_nm(ind)

               end if

            else
               a_nm(ind) = dcmplx(0.0, 0.0)
               b_nm(ind) = dcmplx(0.0, 0.0)
            end if
            a_nm(ind) = -a_nm(ind)
            b_nm(ind) = -b_nm(ind)
         end do
      end do

      las = 0
      do n = 1, Nmax
         las = las + (2*n + 1)**2
      end do

      allocate (rotD(las))
      allocate (indD(las, 2))

      call sph_rotation_sparse(0.0d0, -pi/2.0d0, Nmax, rotD, indD)
      a_nm2 = sparse_matmul(rotD, indD, a_nm, nm_in)
      b_nm2 = sparse_matmul(rotD, indD, b_nm, nm_in)

   end subroutine planewave2

!*****************************************************************
! Calculates fields from coefficients
!
! F = electric field
! G = Magnetic field
   subroutine calc_fields(a_nm, b_nm, k, Po, F, G, inou)
      complex(dp), dimension(:) :: a_nm, b_nm
      real(dp) :: Po(3)
      complex(dp) :: k, F(3), G(3)
      integer :: inou
      integer :: Nmax, n, m, ind, mm
      real(dp) :: r, theta, phi, vec(3), q, omega
      complex(dp) :: kr, alpha, beta, gamma, ccc
      complex(dp), dimension(:), allocatable :: sphj, sphy, sphh
      complex(dp) :: P(3), B(3), C(3), Y, Y1, Y2, M_nm(3), N_nm(3)
      real(dp), dimension(:), allocatable :: L, L1, L2

      Nmax = int(sqrt(dble(1 + size(a_nm)))) - 1

      vec = cart2sph(Po)
      r = vec(1)
      theta = vec(2)
      phi = vec(3)

      kr = k*r

      omega = real(k)*cc

      allocate (sphj(Nmax + 2), sphy(Nmax + 2), sphh(Nmax + 2))

! spherical bessel function at kr
      if (inou == 0) then
         call cspherebessel(Nmax + 1, kr, sphj, sphy)
         sphh = sphj
      else

         sphh = sphankel(Nmax + 1, kr)

      end if

      ind = 0
      F(:) = dcmplx(0.0, 0.0)
      G(:) = dcmplx(0.0, 0.0)

      do n = 1, Nmax
         alpha = sphh(n + 1)
         beta = sqrt(dcmplx(n*(n + 1)))/kr*sphh(n + 1)
         gamma = (n + 1.0d0)/kr*sphh(n + 1) - sphh(n + 2)

         allocate (L(n + 1), L1(n + 2), L2(n))

         call legendre2(n, cos(theta), L)
         call legendre2(n + 1, cos(theta), L1)
         call legendre2(n - 1, cos(theta), L2)

         q = (sqrt(n*(n + 1.0d0)))/((n*2d0 + 1.0d0)*sin(theta)); 
         do m = -n, n
            ind = ind + 1
            mm = abs(m)

            ccc = dcmplx(sqrt((2d0*n + 1.0d0)*factorial(n - mm)/factorial(n + mm)/(4d0*pi)))

            ! Unnormalized complex scalar spherical harmonics
            Y = dcmplx(L(mm + 1)*exp(dcmplx(0.0, m*phi)))
            Y1 = dcmplx(L1(mm + 1)*exp(dcmplx(0.0, m*phi))) 
            if (mm == n) then
               Y2 = dcmplx(0.0, 0.0)
            else
               Y2 = dcmplx(L2(mm + 1)*exp(dcmplx(0.0, m*phi)))
            end if

            ! vector spherical harmonics
            P(:) = dcmplx(0.0, 0.0)
            P(1) = Y

            Y1 = Y1*((n - mm + 1.0d0)/(n + 1.0d0))

            Y2 = Y2*(dcmplx(n + mm)/dcmplx(n))

            B(:) = dcmplx(0.0, 0.0)
            B(2) = Y1 - Y2
            B(3) = ((dcmplx(0.0, m*(2*n + 1.0)))/dcmplx(n*(n + 1.0)))*Y

            B = B*q

            C(:) = dcmplx(0.0, 0.0)
            C(2) = B(3)
            C(3) = -B(2)

            ! Spherical vector wave functions
            M_nm = ccc*alpha*C
            N_nm = ccc*(beta*P + gamma*B)

            F = F + a_nm(ind)*M_nm + b_nm(ind)*N_nm
            G = G + k*(a_nm(ind)*N_nm + b_nm(ind)*M_nm)

         end do

         deallocate (L, L1, L2)

      end do

      F = sph2cart_vec(theta, phi, F)
      G = sph2cart_vec(theta, phi, G)
      G = G/(dcmplx(0.0, omega*mu))

   end subroutine calc_fields

!*****************************************************************
! Calculates fields from coefficients
! F = electric field
! G = Magnetic field
   subroutine calc_MN(MM_nm, NN_nm, Nmax, k, Po, inou)
      real(dp) :: Po(3)
      complex(dp) :: k
      integer :: inou
      integer :: Nmax, n, m, ind, mm
      real(dp) :: r, theta, phi, vec(3), q, omega
      complex(dp) :: kr, alpha, beta, gamma, cc
      complex(dp), dimension(:), allocatable :: sphj, sphy, sphh
      complex(dp) :: P(3), B(3), C(3), Y, Y1, Y2, M_nm(3), N_nm(3)
      real(dp), dimension(:), allocatable :: L, L1, L2
      complex(dp), dimension(3, (Nmax + 1)**2 - 1) :: MM_nm, NN_nm

      vec = cart2sph(Po)
      r = vec(1)
      theta = vec(2)
      phi = vec(3)

      kr = k*r

      omega = real(k)*299792458.0

      allocate (sphj(Nmax + 2), sphy(Nmax + 2), sphh(Nmax + 2))

! spherical bessel functions at kr
      if (inou == 0) then
         call cspherebessel(Nmax + 1, kr, sphj, sphy)
         sphh = sphj
      else
         sphh = sphankel(Nmax + 1, kr)
      end if

      ind = 0
      do n = 1, Nmax
         alpha = sphh(n + 1)
         beta = sqrt(dcmplx(n*(n + 1)))/kr*sphh(n + 1)
         gamma = (n + 1.0d0)/kr*sphh(n + 1) - sphh(n + 2)

         allocate (L(n + 1), L1(n + 2), L2(n))

         call legendre2(n, cos(theta), L)
         call legendre2(n + 1, cos(theta), L1)
         call legendre2(n - 1, cos(theta), L2)

         q = (sqrt(n*(n + 1.0d0)))/((n*2d0 + 1.0d0)*sin(theta)); 
         do m = -n, n
            ind = ind + 1
            mm = abs(m)

            cc = dcmplx(sqrt((2d0*n + 1.0d0)*factorial(n - mm)/factorial(n + mm)/(4d0*pi)))
            ! Unnormalized complex scalar spherical harmonics
            Y = dcmplx(L(mm + 1)*exp(dcmplx(0.0, m*phi)))
            Y1 = dcmplx(L1(mm + 1)*exp(dcmplx(0.0, m*phi))) 
            if (mm == n) then
               Y2 = dcmplx(0.0, 0.0)
            else
               Y2 = dcmplx(L2(mm + 1)*exp(dcmplx(0.0, m*phi)))
            end if

            ! vector spherical harmonics
            P(:) = dcmplx(0.0, 0.0)
            P(1) = Y

            Y1 = Y1*((n - mm + 1.0d0)/(n + 1.0d0))

            Y2 = Y2*(dcmplx(n + mm)/dcmplx(n))

            B(:) = dcmplx(0.0, 0.0)
            B(2) = Y1 - Y2
            B(3) = ((dcmplx(0.0, m*(2*n + 1.0)))/dcmplx(n*(n + 1.0)))*Y

            B = B*q

            C(:) = dcmplx(0.0, 0.0)
            C(2) = B(3)
            C(3) = -B(2)

            ! Spherical vector wave functions
            M_nm = cc*alpha*C
            N_nm = cc*(beta*P + gamma*B)

            MM_nm(:, ind) = sph2cart_vec(theta, phi, M_nm)
            NN_nm(:, ind) = sph2cart_vec(theta, phi, N_nm)
         end do

         deallocate (L, L1, L2)

      end do

   end subroutine calc_MN

!******************************************************************************

   subroutine mueller_matrix_coeff(a_nm, b_nm, a_nm2, b_nm2, k, N_theta, N_phi, S_out)

      complex(dp), dimension(:) :: a_nm, b_nm, a_nm2, b_nm2
      complex(dp) :: k
      real(dp), dimension(:, :), allocatable :: S_out
      integer :: N_phi, N_theta
      complex(dp), dimension(3) :: E_out, E_out2, H_out, H_out2

      integer :: i1, i2, las
      real(dp) :: theta, phi, abcd(2, 2), r(3), RR, unit_th(3), unit_phi(3)
      complex(dp), dimension(N_phi*N_theta) :: f11, f12, f21, f22
      complex(dp) ::S1, S2, S3, S4
      complex(dp) :: i

      i = dcmplx(0.0, 1.0)
      RR = 1.0d6
      allocate (S_out(N_phi*N_theta, 18))

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

            call calc_fields(a_nm, b_nm, k, r, E_out, H_out, 1)
            call calc_fields(a_nm2, b_nm2, k, r, E_out2, H_out2, 1)

            unit_th = [cos(theta)*cos(phi), sin(phi)*cos(theta), -sin(theta)]
            unit_phi = [-sin(phi), cos(phi), 0.0d0]; 
            f11(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out, dcmplx(unit_th))
            f21(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out, dcmplx(unit_phi))
            f12(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out2, dcmplx(unit_th))
            f22(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out2, dcmplx(unit_phi))
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

   end subroutine mueller_matrix_coeff

!******************************************************************************

   subroutine scattering_matrix_coeff(a_nm, b_nm, a_nm2, b_nm2, k, N_points, points, S_out)

      complex(dp), dimension(:) :: a_nm, b_nm, a_nm2, b_nm2
      complex(dp) :: k
      real(dp), dimension(:, :), allocatable :: S_out, points
      integer :: N_points
      complex(dp), dimension(3) :: E_out, E_out2, H_out, H_out2

      integer :: i1, las
      real(dp) :: theta, phi, abcd(2, 2), r(3), RR, unit_th(3), unit_phi(3)
      complex(dp), dimension(N_points) :: f11, f12, f21, f22
      complex(dp) ::S1, S2, S3, S4
      complex(dp) :: i

      i = dcmplx(0.0, 1.0)
      RR = 1.0d6
      allocate (S_out(N_points, 18))

      las = 0
      do i1 = 1, N_points
         theta = points(1, i1)
         phi = points(2, i1)

         abcd(1, :) = [cos(phi), sin(phi)]
         abcd(2, :) = [sin(phi), -cos(phi)]

         r(1) = RR*sin(theta)*cos(phi)
         r(2) = RR*sin(theta)*sin(phi)
         r(3) = RR*cos(theta)

         call calc_fields(a_nm, b_nm, k, r, E_out, H_out, 1)
         call calc_fields(a_nm2, b_nm2, k, r, E_out2, H_out2, 1)

         unit_th = [cos(theta)*cos(phi), sin(phi)*cos(theta), -sin(theta)]
         unit_phi = [-sin(phi), cos(phi), 0.0d0]; 
         f11(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out, dcmplx(unit_th))
         f21(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out, dcmplx(unit_phi))
         f12(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out2, dcmplx(unit_th))
         f22(las + 1) = k*RR/exp(i*k*(RR))*dot_product(E_out2, dcmplx(unit_phi))
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

   end subroutine scattering_matrix_coeff

!******************************************************************************

   subroutine extinction_matrix_coeff(a_nm, b_nm, a_nm2, b_nm2, k, theta, phi, S_out)
      complex(dp), dimension(:) :: a_nm, b_nm, a_nm2, b_nm2
      complex(dp) :: k
      real(dp), dimension(:, :), allocatable :: S_out
      complex(dp), dimension(3) :: E_out, E_out2, H_out, H_out2

      real(dp) :: theta, phi, abcd(2, 2), r(3), RR, unit_th(3), unit_phi(3)
      complex(dp) :: f11, f12, f21, f22
      complex(dp) :: S11, S12, S21, S22
      complex(dp) :: i

      i = dcmplx(0.0, 1.0)
      RR = 1.0d6

      allocate (S_out(1, 18))

      abcd(1, :) = [cos(phi), sin(phi)]
      abcd(2, :) = [sin(phi), -cos(phi)]

      r(1) = RR*sin(theta)*cos(phi)
      r(2) = RR*sin(theta)*sin(phi)
      r(3) = RR*cos(theta)

      call calc_fields(a_nm, b_nm, k, r, E_out, H_out, 1)
      call calc_fields(a_nm2, b_nm2, k, r, E_out2, H_out2, 1)

      unit_th = [cos(theta)*cos(phi), sin(phi)*cos(theta), -sin(theta)]
      unit_phi = [-sin(phi), cos(phi), 0.0d0]; 
      f11 = k*RR/exp(i*k*(RR))*dot_product(E_out, dcmplx(unit_th))
      f21 = k*RR/exp(i*k*(RR))*dot_product(E_out, dcmplx(unit_phi)) 
      f12 = k*RR/exp(i*k*(RR))*dot_product(E_out2, dcmplx(unit_th))
      f22 = k*RR/exp(i*k*(RR))*dot_product(E_out2, dcmplx(unit_phi)) 
      S22 = -i*(f21*abcd(1, 2) + f22*abcd(2, 2))
      S11 = -i*(f11*abcd(1, 1) + f12*abcd(2, 1))
      S12 = i*(f11*abcd(1, 2) + f12*abcd(2, 2))
      S21 = i*(f21*abcd(1, 1) + f22*abcd(2, 1))

      S_out(1, 1) = phi*dble(k)/(2d0*pi)
      S_out(1, 2) = theta*dble(k)/(2d0*pi)

! Mueller matrix
      S_out(1, 3) = imag(S11+S22) !K11
      S_out(1, 4) = imag(S11-S22) !K12
      S_out(1, 5) = -imag(S12+S21) !K13
      S_out(1, 6) = real(S21-S12) !K14
      S_out(1, 7) = imag(S11-S22) !K21
      S_out(1, 8) = imag(S11+S22) !K22
      S_out(1, 9) = imag(S21-S12) !K23
      S_out(1, 10) = -real(S12+S21) !K24
      S_out(1, 11) = -imag(S12+S21) !K31
      S_out(1, 12) = -imag(S21-S12) !K32
      S_out(1, 13) = imag(S11+S22) !K33
      S_out(1, 14) = real(S22-S11) ! K34
      S_out(1, 15) = real(S21-S12) ! K41
      S_out(1, 16) = real(S12+S21) ! K42
      S_out(1, 17) = -real(S22-S11) ! K43
      S_out(1, 18) = imag(S11+S22) ! K44

      S_out = S_out*2d0*pi/dble(k)

   end subroutine extinction_matrix_coeff

!******************************************************************************

   subroutine tr_T(Taa, Tbb, Tab, Tba, k, crs)

      complex(dp), dimension(:, :) :: Taa, Tbb, Tab, Tba
      integer :: nm, mu, nu
      real(dp) :: T1, T2, crs(4)
      real(dp) :: k

      nm = size(Taa, 1)

      T1 = 0d0
      T2 = 0d0

      do mu = 1, nm
         T1 = T1 + (real(Taa(mu, mu)) + real(Tbb(mu, mu)) + real(Tab(mu, mu)) + real(Tba(mu, mu)))
         do nu = 1, nm

            T2 = T2 + (abs(Taa(mu, nu))**2 + abs(Tbb(mu, nu))**2 + &
                       abs(Tab(mu, nu))**2 + abs(Tba(mu, nu))**2)
         end do
      end do

      crs(1) = -2*pi*T1/k**2
      crs(2) = 2*pi*T2/k**2
      crs(3) = (-2*pi*T1/k**2 - 2*pi*T2/k**2)
      crs(4) = T1 + T2

      print *, 'Ave. Extinction cross section:', -2*pi*T1/k**2
      print *, 'Ave. Scattering cross section:', 2*pi*T2/k**2
      print *, 'Ave. absorption cross section:', -2*pi*T1/k**2 - 2*pi*T2/k**2
      print *, 'Tr_(real(T) + T adj(T)):', T1 + T2

   end subroutine tr_T

end module mie

