module postprocessing
use integrator

contains

!******************************************************************************

subroutine test_methods( matrices, mesh )
type(data) :: matrices
type(mesh_struct) :: mesh
integer :: i, j, ii, range(2)
real(dp), dimension(:,:), allocatable :: Q_fcoll, Q_tcoll
real(dp), dimension(3) :: k0, E0, E90
complex(dp), dimension(3) :: F, N, FF, NN
real(dp) :: Qmag

! Now, to comply with the properties of DDSCAT, the scattering happens
! with respect to the particle principal axes
matrices%khat = matrices%P(:,3)
matrices%Rexp = find_Rexp( matrices )
matrices%Rexp = transpose( matrices%Rexp )

call rot_setup(matrices)

allocate(Q_fcoll(3,matrices%bars))
allocate(Q_tcoll(3,matrices%bars))

if(matrices%whichbar == 0) then
	range(1) = 1
	range(2) = matrices%bars
else
	range(1) = matrices%whichbar
	range(2) = matrices%whichbar
	write(*,'(A, 20F6.3)') ' Chosen wavelength: ', 2d6*pi/mesh%ki(matrices%whichbar)
end if

do j = 1,2
	N = dcmplx(0.0d0, 0.0d0)
	F = dcmplx(0.0d0, 0.0d0)
	Q_tcoll = 0d0
	Q_fcoll = 0d0
	
	do i = range(1),range(2)
		select case (j)
			case(1); call forcetorque_num(i, matrices, mesh)
			case(2); call forcetorque(i, matrices, mesh)
		end select
		ii = i
		if(matrices%whichbar > 0) ii = matrices%whichbar
		F = F + matrices%force
		N = N + matrices%torque
		Q_fcoll(:,ii) = matmul(matrices%R,matrices%Q_f)
		Q_tcoll(:,ii) = matmul(matrices%R,matrices%Q_t)
	end do

	NN =  matmul(matrices%R, N)
	FF =  matmul(matrices%R, F)

	print*, ''
	select case (j)
		case(1); write(*,'(A)') 'Numerical integration of MST:'
		case(2); write(*,'(A)') 'Analytical VSWF coefficient integration of MST:'
	end select
	write(*,'(A,3ES11.3,A)') ' F = (', real(FF), ' ) N'
	write(*,'(A,3ES11.3,A)') ' N = (', real(NN), ' ) Nm'

	print*, ''
	print*, 'Efficiencies for each wavelength separately'
	write(*,'(A, 80F7.1)') 'WL(nm):', 2d9*pi/mesh%ki
	write(*,'(A, 80F7.3)') '   Qfx:', Q_fcoll(1,:)
	write(*,'(A, 80F7.3)') '   Qfy:', Q_fcoll(2,:)
	write(*,'(A, 80F7.3)') '   Qfz:', Q_fcoll(3,:)

	write(*,'(A, 80F7.3)') '   Qtx:', Q_tcoll(1,:)
	write(*,'(A, 80F7.3)') '   Qty:', Q_tcoll(2,:)
	write(*,'(A, 80F7.3)') '   Qtz:', Q_tcoll(3,:)
end do

!call print_mat(matrices%P, 'P')
!write(*,'(A, 3ES12.4)')  'Ip   =', matrices%Ip
!write(*,'(A, 3F8.4)') 'E0   =', real(matrices%E0hat)
!write(*,'(A, 3F8.4)') 'E90  =', real(matrices%E90hat)
!write(*,'(A, 3F8.4)') 'khat =', matrices%khat

!open(unit=2, file="out/Q.dat", ACTION="write", STATUS="replace")
!write(2,'(A)') "      ka        Qf_k       Qt_1       Qt_2       Qt_3       |Qt|"
!do i = range(1),range(2)
!	ii = i
!	if(matrices%whichbar > 0) ii = matrices%whichbar
!	Qmag = dsqrt(Q_tcoll(1,ii)**2 + Q_tcoll(2,ii)**2 + Q_tcoll(3,ii)**2)
!	write(2,'(9F11.7)') mesh%ki(ii)*mesh%a, &
!	dot_product(Q_fcoll(:,ii),matrices%khat), &
!	dot_product(Q_tcoll(:,ii),matrices%khat), &
!	dot_product(Q_tcoll(:,ii),-matrices%P(:,2)), &
!	dot_product(Q_tcoll(:,ii),-matrices%P(:,1)),  Qmag
!end do
!close(2)

end subroutine test_methods

!******************************************************************************
! Calculates the torque efficiency of a particle averaged over rotation
! about its major axis of inertia (a_3) for angles 0 to pi between 
! a_3 and incident k0. Results are automatically comparable with DDSCAT.
subroutine torque_efficiency( matrices, mesh )
type(data) :: matrices
type(mesh_struct) :: mesh
integer :: i, j, k, Nang, Bang
real(dp), dimension(3, 3) ::  R_beta, R_thta
real(dp), dimension(3) :: k0, E0, E90, Q_t, nbeta, a_3
real(dp) :: theta, beta
real(dp), dimension(:,:), allocatable :: Q_coll

call rot_setup(matrices)

k0 = matrices%khat
E0 = real(matrices%E0hat)
E90 = real(matrices%E90hat)

a_3 = matrices%P(1:3,3)
!~ write(*,'(A,3F7.3)') '  a_3  =', a_3

Nang = 180 ! Angle of rotation of a_3 about e_1 (when psi=0)
Bang = 360 ! Angle of averaging over beta

allocate(Q_coll(3,Nang))
Q_coll(:,:) = 0d0

write(*,'(A)') '  Starting the calculation of beta-averaged torque efficiency:'
! Theta loop
do i=0,Nang-1
	theta = dble(i)*pi/180d0
	R_thta = R_theta(matrices, theta)

	! Rotation axis for beta averaging for current theta
	nbeta = matmul(R_thta,a_3) ! Beta rotation about a_3
	nbeta = nbeta/vlen(nbeta) ! Ensure unit length of axis vector

	! Beta averaging loop
	do j = 0,Bang-1
		Q_t = 0d0

		! Find rotation of angle beta around the rotated a_3-axis
		beta = dble(j)*pi/Bang*2d0
		R_beta = R_aa(nbeta,beta)

		! The ultimate rotation matrices for scattering event
		matrices%R = matmul(R_beta,R_thta) ! First a_3 to theta, then beta about a_3
		call rot_setup(matrices)
		
		if (matrices%whichbar == 0) then
			do k = 1,matrices%bars
				call forcetorque(k, matrices, mesh)
				Q_t = Q_t + matrices%Q_t
			end do
		else
			k = matrices%whichbar
			call forcetorque(k, matrices, mesh)
			Q_t = Q_t + matrices%Q_t
		end if
		
		! Flip the coordinate labels to match Lazarian2007b
		Q_t = matmul(matrices%Rkt,Q_t)
		Q_t = [dot_product(Q_t,k0),dot_product(Q_t,E0),dot_product(Q_t,E90)]
		Q_coll(:,i+1) = Q_coll(:,i+1) + Q_t
	end do
	call print_bar(i+1,Nang)

end do
Q_coll = Q_coll/Bang
if(matrices%whichbar==0) Q_coll = Q_coll/matrices%bars

open(unit=1, file="out/Q.out", ACTION="write", STATUS="replace")
write(1,'(A)') 'cos(theta)  Q_{t,1} Q_{t,2} Q_{t,3}'
do i = 0,Nang-1
theta = dble(i)*pi/180d0
write(1,'(4E12.3)') dcos(theta), Q_coll(:,i+1)
end do
close(1)

end subroutine torque_efficiency

!******************************************************************************
! Compute the mueller matrix. Bad code ahead. 
subroutine compute_mueller(matrices,mesh)
type(data) :: matrices
type(mesh_struct) :: mesh
integer :: i
complex(dp), dimension(:) , allocatable :: p, q, p90, q90
integer :: Nmax, ii, las, nm, N_points, halton_init, N_avgs, N_theta, N_phi
real(dp) :: E, vec(3), phi, R0(3,3), Qt(3,3), omega(3), aproj(3), RR(3,3), theta
real(dp), dimension(:,:), allocatable :: SS, SSS
CHARACTER(LEN=80) :: mueller_out
CHARACTER(LEN=8) :: mode

N_points = 500 ! Number of points to calculate the perfect orientations
N_theta = 180  ! Theta range in Mueller matrices
N_avgs = 1     ! Used in randomly oriented (RO) averages
N_phi = 90     ! Number of phi-angles between 0 and 2pi for non-RO cases

mueller_out = trim(matrices%mueller)//'-'//trim(matrices%mueller_mode)
if(file_exists(mueller_out)) then
	print*, ' Mueller matrix already exists, quitting...'
	stop
end if

select case (trim(matrices%mueller_mode))
	case('ave') 
	! Note: Averaging could be done by summing phi-variations, when rotations
	! of the particle are handled in systematic fashion, which may be more
	! computationally efficient.
		matrices%R = eye(3)
		N_avgs = 720 ! Number of averaging directions
		N_points = 1 ! N_avgs replaces this
		N_phi = 1    ! Averaging occurs, so phi-dependency is lost
	case('ori')
		call read_log(matrices,mesh)
end select

allocate(SSS(N_theta*N_phi,18))
SSS = 0d0
halton_init = 0

! Choose wavelength
ii = matrices%whichbar
if(ii==0) ii = 1

E = matrices%E_rel(ii)*matrices%E
mesh%k = mesh%ki(ii)

Nmax = matrices%Nmaxs(ii)
las = (Nmax  + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1
nm = (Nmax + 1)**2 - 1

select case (trim(matrices%mueller_mode))
	case('perf_ori')
		! Hard coded orientation data
		omega = [ 1d0, 0d0, 0d0 ]
		omega = omega/vlen(omega) 
		R0 = rotate_a_to_b(matrices%P(:,3),omega)
		Qt = matmul(R0,matrices%P)
		aproj = [Qt(1,3),Qt(2,3),0d0]
		aproj = aproj/vlen(aproj)
	case('ori'); N_points = size(matrices%RRR,3)
end select

do i2 = 1, N_points
	do i1 = 1, N_avgs
		if(trim(matrices%mueller_mode)=='ave') then
			vec(1) = 1d0
			vec(2) = acos(2*halton_seq(halton_init+i1, 2)-1)
			vec(3) = halton_seq(halton_init+i1, 3)*2*pi
			matrices%khat = [dsin(vec(2))*dcos(vec(3)), dsin(vec(2))*dsin(vec(3)), dcos(vec(2))]
		else if(trim(matrices%mueller_mode)=='ori') then
			RR = transpose(matmul(matrices%R_al,matrices%RRR(:,:,i2)))
			matrices%khat = matmul(RR,[0d0,0d0,1d0])
		else if(trim(matrices%mueller_mode)=='perf_ori') then
			theta = dble(i2-1)*2*pi/(N_points-1)
			RR = matmul(transpose(R_aa(omega,theta)),transpose(R0))
			if(N_avgs == 1) then
				phi = dacos(aproj(1))
			else
				phi = dble(i1-1)*pi/2d0/(N_avgs-1)
			end if
			matrices%khat = matmul(RR,[0d0,0d0,1d0])
		end if
		
		matrices%khat = -matrices%khat/vlen(matrices%khat)
		matrices%Rexp = find_Rexp( matrices ) 
		matrices%Rexp = transpose( matrices%Rexp )
		
		if(trim(matrices%mueller_mode)=='perf_ori') matrices%R = R_aa(matrices%khat, phi)

		call rot_setup(matrices)
		call mueller_fields(matrices,E,p,q,p90,q90,ii)
		if(allocated(SS)) deallocate(SS)
		call mueller_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_theta, N_phi, Nmax, SS)
		SSS = SSS + SS
		if(trim(matrices%mueller_mode)=='ave') call print_bar(i1,N_avgs)
	end do
	call print_bar(i2,N_points)
end do

SSS = SSS/N_points/N_avgs
call write_mueller(SSS,mueller_out)

end subroutine compute_mueller

!*******************************************************************************

subroutine mueller_fields(matrices,E,p,q,p90,q90,ii)
type(data) :: matrices
real(dp) :: E
complex(dp), dimension(:) , allocatable :: a_in, b_in, a90, b90, &
a, b, p, q, p90, q90, ptemp, qtemp, p90temp, q90temp
complex(dp), dimension(:, :), allocatable :: Taa, Tab, Tba, Tbb
complex(8), dimension(:), allocatable :: rotD, rotD90, rbak
integer, dimension(:, :), allocatable :: indD, indD90, ibak
integer :: nm, las, ii, Nmax

Nmax = matrices%Nmaxs(ii)
las = (Nmax  + 1)*(2*Nmax + 1)*(2*Nmax + 3)/3 - 1
nm = (Nmax + 1)**2 - 1

allocate(rbak(las), ibak(las,2))
allocate(a(nm),b(nm),a90(nm),b90(nm),ptemp(nm),qtemp(nm),p90temp(nm),q90temp(nm))
allocate(rotD(las), indD(las,2), rotD90(las), indD90(las,2))
allocate(Taa(nm,nm),Tab(nm,nm),Tba(nm,nm),Tbb(nm,nm))

rotD = matrices%rotDs(1:las, ii)
indD = matrices%indDs(1:las, :, ii)
rotD90 = matrices%rotD90s(1:las, ii)
indD90 = matrices%indD90s(1:las, :, ii)

Taa = matrices%Taai(1:nm, 1:nm, ii)
Tab = matrices%Tabi(1:nm, 1:nm, ii)
Tba = matrices%Tbai(1:nm, 1:nm, ii)
Tbb = matrices%Tbbi(1:nm, 1:nm, ii)

a_in = E*matrices%as(1:nm, ii)
b_in = E*matrices%bs(1:nm, ii)

a = sparse_matmul(rotD, indD, a_in, nm)
b = sparse_matmul(rotD, indD, b_in, nm)
a90 = sparse_matmul(rotD90, indD90, a_in, nm)
b90 = sparse_matmul(rotD90, indD90, b_in, nm)

ptemp = matmul(Taa, a) + matmul(Tab, b)
qtemp = matmul(Tbb, b) + matmul(Tba, a)
p90temp = matmul(Taa, a90) + matmul(Tab, b90)
q90temp = matmul(Tbb, b90) + matmul(Tba, a90)

call sph_rotation_sparse_gen2(matrices%Rkt, Nmax, rbak, ibak)

p = sparse_matmul(rbak, ibak, ptemp, nm)
q = sparse_matmul(rbak, ibak, qtemp, nm)
p90 = sparse_matmul(rbak, ibak, p90temp, nm)
q90 = sparse_matmul(rbak, ibak, q90temp, nm)

end subroutine mueller_fields

!*******************************************************************************
subroutine stability_analysis(matrices,mesh)
type(data) :: matrices
type(mesh_struct) :: mesh
integer :: i, j, k, l, ind, Npoints, Nang, Bang, N_J, Nphi, Ntheta
real(dp), dimension(3, 3) ::  R_beta, R_thta
real(dp), dimension(3) :: Q_t, nbeta
real(dp) :: thta, beta
real(dp), dimension(:,:), allocatable :: Q_coll, vec
real(dp) :: res1, res2, res3, mx
real(dp), dimension(3) :: k0, a_1, a_2, a_3, knew, NN, maxdir, anew
complex(dp), dimension(3) :: N
real(dp), dimension(:), allocatable :: theta, costheta, phi 

call rot_setup(matrices)

k0 = matrices%khat

a_1 = matrices%P(1:3,1)
a_2 = matrices%P(1:3,2)
a_3 = matrices%P(1:3,3)
!~ write(*,'(A,3F7.3)') '  a_3  =', a_3

Ntheta = 90
Nphi = 180
allocate(theta(Ntheta),costheta(Ntheta),phi(Nphi))
call linspace(0d0,2d0*pi,Nphi,phi)
call linspace(1d0,-1d0,Ntheta,costheta)
do i = 1,Ntheta
	theta(i) = dacos(costheta(i))
end do
Npoints = Ntheta*Nphi
allocate(vec(3,Npoints))
vec = 0d0
mx = 0d0

!~ vec = fibonacci_sphere(Npoints,1)
vec = uniform_sphere(Npoints,theta,phi)

open(unit=1, file="out/NdotQ.dat", ACTION="write", STATUS="replace")
write(1,'(2(A,I0))') 'Ntheta = ', Ntheta, ' | Nphi = ', Nphi
write(1,'(A)') '  x   y   z   N.a_1   N.a_2   N.a_3'

do i = 1,Npoints
	N = dcmplx(0d0)
	knew = vec(:,i)       
	! The ultimate rotation matrices for scattering event
	matrices%R = rotate_a_to_b(k0,knew)
	call rot_setup(matrices)
	
	if (matrices%whichbar == 0) then
		do k = 1,matrices%bars
			call forcetorque(k, matrices, mesh)
			N = N + matrices%Q_t
		end do
	else
		k = matrices%whichbar
		call forcetorque(k, matrices, mesh)
		N = N + matrices%Q_t
	end if
	anew = matmul(transpose(matrices%R),a_3)
	NN = matmul(transpose(matrices%P) , real(N))
	res1 = dot_product(NN,a_1)
	res2 = dot_product(NN,a_2)
	res3 = dot_product(NN,a_3)
	if (res3 > mx) then
		mx = res3
		maxdir = anew
	end if
	write(1,'(6E12.3)') anew(1), anew(2), anew(3), res1, res2, res3
	call print_bar(i+1,Npoints)
end do

close(1)
res = dacos(dot_product([0d0,0d0,1d0],maxdir/vlen(maxdir)))*180d0/pi
write(*,'(A,3F7.3,A,F7.3)') 'Stablest direction = (', maxdir, ' ), angle between a_3 and k = ', res

N_J  = 5
Nang = 180 ! Angle of rotation of a_3 about e_1 (when psi=0)
Bang = 360 ! Angle of averaging over beta

allocate(Q_coll(3,Nang))
Q_coll(:,:) = 0d0
ind = 0

open(unit=1, file="out/N.out", ACTION="write", STATUS="replace")
write(1,'(A)') 'cos(theta)  N_k   N_E0    N_E90'

! Theta loop
do i=0,Nang-1
	thta = dble(i)*pi/180d0
	R_thta = R_theta(matrices, thta)
	N = 0d0

	! Rotation axis for beta averaging for current theta
	nbeta = matmul(R_thta,a_3) ! Beta rotation about a_3
	nbeta = nbeta/vlen(nbeta) ! Ensure unit length of axis vector

	! Beta averaging loop
	do j = 0,Bang-1
		Q_t = 0d0

		! Find rotation of angle beta around the rotated a_3-axis
		beta = dble(j)*pi/Bang*2d0
		R_beta = R_aa(nbeta,beta)

		! The ultimate rotation matrices for scattering event
		matrices%R = matmul(R_beta,R_thta) ! First a_3 to theta, then beta about a_3
		call rot_setup(matrices)
		if (matrices%whichbar == 0) then
			do k = 1,matrices%bars
				call forcetorque(k, matrices, mesh)
				N = N + matrices%torque
			end do
		else
			k = matrices%whichbar
			call forcetorque(k, matrices, mesh)
			N = N + matrices%torque
		end if
		! Flip the coordinate labels to match Lazarian2007b
		NN = matmul(matrices%Rkt,real(N))
		NN = [dot_product(NN,matrices%khat),&
		dot_product(NN,real(matrices%E0hat)),dot_product(NN,real(matrices%E90hat))]
		Q_coll(:,i+1) = Q_coll(:,i+1) + NN
	end do
!~ 	NN = matmul(transpose(matrices%R) , real(N))
	ind = ind + 1
	call print_bar(ind,Nang)

	write(1,'(4E14.5)') dcos(thta), Q_coll(:,i+1)
	
end do
close(1)

end subroutine stability_analysis

end module postprocessing
