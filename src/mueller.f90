module mueller
use shapebeam

implicit none

contains

!******************************************************************************
! Compute the mueller matrix. Bad code ahead. 
subroutine compute_mueller(matrices,mesh)
type(data) :: matrices
type(mesh_struct) :: mesh
integer :: i, iii1, iii2
complex(dp), dimension(:) , allocatable :: p, q, p90, q90
integer :: Nmax, ii, las, nm, N_points, halton_init, N_avgs, N_theta, N_phi, everyn
real(dp) :: E, vec(3), phi, R0(3,3), Qt(3,3), omega(3), aproj(3), RR(3,3), theta
real(dp), dimension(:,:), allocatable :: SS, SSS
CHARACTER(LEN=80) :: mueller_out
CHARACTER(LEN=8) :: mode

N_points = 500 ! Number of points to calculate the perfect orientations
N_theta = 180  ! Theta range in Mueller matrices
N_avgs = 1     ! Used in randomly oriented (RO) averages
N_phi = 90     ! Number of phi-angles between 0 and 2pi for non-RO cases
everyn = 1

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
      everyn = 1
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

do iii2 = 1, N_points, everyn
   do iii1 = 1, N_avgs
      if(trim(matrices%mueller_mode)=='ave') then
         vec(1) = 1d0
         vec(2) = acos(2*halton_seq(halton_init+iii1, 2)-1)
         vec(3) = halton_seq(halton_init+iii1, 3)*2*pi
         matrices%khat = [dsin(vec(2))*dcos(vec(3)), dsin(vec(2))*dsin(vec(3)), dcos(vec(2))]
      else if(trim(matrices%mueller_mode)=='ori') then
         RR = transpose(matmul(matrices%R_al,matrices%RRR(:,:,iii2)))
         matrices%khat = matmul(RR,[0d0,0d0,1d0])
      else if(trim(matrices%mueller_mode)=='perf_ori') then
         theta = dble(iii2-1)*2*pi/(N_points-1)
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
      call scattered_fields(matrices,E,p,q,p90,q90,ii)
      if(allocated(SS)) deallocate(SS)
      call mueller_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_theta, N_phi, Nmax, SS)
      SSS = SSS + SS
      if(trim(matrices%mueller_mode)=='ave') call print_bar(iii1,N_avgs)
   end do
   call print_bar(iii2,N_points)
end do

SSS = everyn*SSS/N_points/N_avgs
call write_mueller(SSS,mueller_out)

end subroutine compute_mueller

!******************************************************************************

subroutine test_mueller(matrices, mesh)
type(data) :: matrices
type(mesh_struct) :: mesh
integer :: i, halton_init, N_points
real(dp) :: al_direction(3), vec(3)
real(dp), dimension(:), allocatable :: a_dist
real(dp), dimension(:,:), allocatable :: points

halton_init = 0
N_points = 18*9
allocate(points(2,N_points))

do i = 1, N_points
   points(1,i) = acos(2*halton_seq(halton_init+i, 2)-1)
   points(2,i) = halton_seq(halton_init+i, 3)*2*pi
end do

a_dist = mesh%ki*mesh%a/mesh%ki(2)
al_direction = [0d0,0d0,1d0]
call scattering_extinction_matrices(matrices, mesh, a_dist, points, al_direction)

end subroutine test_mueller

!******************************************************************************

subroutine scattering_extinction_matrices(matrices, mesh, a_dist, points, al_direction)
type(data) :: matrices
type(mesh_struct) :: mesh
integer :: i, ind, tot_ind
complex(dp), dimension(:) , allocatable :: p, q, p90, q90
integer :: Nmax, ii, las, nm, N_points, halton_init, N_avgs, N_theta, N_phi, N_size, N_ia
real(dp) :: E, vec(3), phi, R0(3,3), Qt(3,3), omega(3), aproj(3), RR(3,3), theta, & 
al_direction(3), inc_angles(2)
real(dp), dimension(:), allocatable :: a_dist
real(dp), dimension(:,:), allocatable :: SS, SSS, points
CHARACTER(LEN=80) :: mueller_out
CHARACTER(LEN=8) :: mode

N_points = size(points,2) ! Number of points to calculate the perfect orientations
N_avgs = 1     ! Used in randomly oriented (RO) averages

mueller_out = trim(matrices%mueller)
if(file_exists(mueller_out)) then
   print*, ' Mueller matrix already exists, quitting...'
   stop
end if

inc_angles = [90, 180]
allocate(SSS(N_points*size(a_dist,1)*size(inc_angles),19))
SSS = 0d0
allocate(SS(N_points,18))

ind = 1
do N_size = 1, size(a_dist,1)
   ! Choose wavelength
   mesh%k = mesh%ki(N_size)

   Nmax = matrices%Nmaxs(N_size)
   do N_ia = 1,size(inc_angles,1)
      ! We take every N_size as the critical size, below which nothing is aligned
      SS = 0d0
      do ii = 1, size(a_dist,1)
         if(a_dist(ii)<=a_dist(N_size)) then
            SS = SS + mueller_ave(matrices, mesh, points, ii)/size(a_dist,1)
         else
            SS = SS + mueller_align(matrices, mesh, points, ii, al_direction)/size(a_dist,1)
         end if
      end do

      do i = 1, N_points
         SSS(ind, 1) = N_size
         SSS(ind, 2) = N_ia
         SSS(ind, 3) = i
         SSS(ind, 4:19) = SS(i, 3:18)
         ind = ind + 1
         call print_bar(ind,size(a_dist,1)*size(inc_angles)*N_points)
      end do
   end do
end do
call write_mueller2(SSS,mueller_out)

end subroutine scattering_extinction_matrices

!******************************************************************************

function mueller_ave(matrices, mesh, points, ii) result(SS)
type(data) :: matrices
type(mesh_struct) :: mesh
integer :: i, ii, N_avgs, halton_init
real(dp) :: E, vec(3)
real(dp), dimension(:,:), allocatable :: S, SS, points
complex(dp), dimension(:) , allocatable :: p, q, p90, q90

matrices%R = eye(3)
N_avgs = 720 ! Number of averaging directions
halton_init = 0
E = matrices%E_rel(ii)*matrices%E
allocate(SS(size(points,2),18))

do i = 1, N_avgs
   vec(1) = 1d0
   vec(2) = acos(2*halton_seq(halton_init+i, 2)-1)
   vec(3) = halton_seq(halton_init+i, 3)*2*pi
   matrices%khat = [dsin(vec(2))*dcos(vec(3)), dsin(vec(2))*dsin(vec(3)), dcos(vec(2))]

   matrices%khat = -matrices%khat/vlen(matrices%khat)
   matrices%Rexp = find_Rexp( matrices ) 
   matrices%Rexp = transpose( matrices%Rexp )

   call rot_setup(matrices)
   call scattered_fields(matrices,E,p,q,p90,q90,ii)
   if(allocated(S)) deallocate(S)
   call scattering_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_avgs, points, S)
   SS = SS + S
end do

SS = SS/N_avgs

end function mueller_ave

!******************************************************************************

function mueller_align(matrices, mesh, points, ii, al_direction) result(SS)
type(data) :: matrices
type(mesh_struct) :: mesh
integer :: i, ii, N_avgs, halton_init
real(dp) :: E, omega(3), al_direction(3), theta, phi, RR(3,3), Qt(3,3), R0(3,3), aproj(3)
real(dp), dimension(:,:), allocatable :: S, SS, points
complex(dp), dimension(:) , allocatable :: p, q, p90, q90

matrices%R = eye(3)
N_avgs = 720 ! Number of averaging directions

E = matrices%E_rel(ii)*matrices%E
omega = al_direction
omega = omega/vlen(omega) 
R0 = rotate_a_to_b(matrices%P(:,3),omega)
Qt = matmul(R0,matrices%P)
aproj = [Qt(1,3),Qt(2,3),0d0]
aproj = aproj/vlen(aproj)
phi = dacos(aproj(1))

allocate(SS(size(points,2),18))

do i = 1, N_avgs
   theta = dble(i-1)*2*pi/(N_avgs-1)
   RR = matmul(transpose(R_aa(omega,theta)),transpose(R0))

   matrices%khat = matmul(RR,[0d0,0d0,1d0])

   matrices%khat = -matrices%khat/vlen(matrices%khat)
   matrices%Rexp = find_Rexp( matrices ) 
   matrices%Rexp = transpose( matrices%Rexp )
   matrices%R = R_aa(matrices%khat, phi)

   call rot_setup(matrices)
   call scattered_fields(matrices,E,p,q,p90,q90,ii)
   if(allocated(S)) deallocate(S)
   call scattering_matrix_coeff(p, q, p90, q90, dcmplx(mesh%k), N_avgs, points, S)
   SS = SS + S
end do

SS = SS/N_avgs

end function mueller_align

end module mueller