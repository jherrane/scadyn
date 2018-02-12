module T_matrix
use transformation_matrices
use precorrection
use gmres_module
use setup
use common

implicit none 

contains

!******************************************************************************

subroutine calc_T(matrices, mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: i, ii, nm, sz

sz = size(mesh%ki)
if(matrices%singleT == 1) sz = 1

write(*,'(3(A,I0))') ' Construct matrices for ', sz, ' wavelengths...'
if (use_mie == 1) then
 call mie_T_matrix(matrices, mesh)
else
 do i = 1,sz
  if(size(mesh%params,2)>1) mesh%param = mesh%params(:,i)
  ii = i
  if(matrices%singleT == 1) ii = matrices%whichbar
  write(*,'(3(A,I0))') ' Step ', ii, '/', sz, ''
  nm = (matrices%Nmaxs(ii)+1)**2-1 
  mesh%k = mesh%ki(ii)

  call allocate_T(ii, matrices)

  call update_projections(matrices,ii)

  if(allocated(matrices%Fg)) deallocate(matrices%Fg)
  if(allocated(matrices%sp_mat)) deallocate(matrices%sp_mat, matrices%sp_ind)
  call build_G(matrices, mesh)

  print*,' Constructing the sparse part of the system matrix'
  if(mesh%order == 0) then
   call compute_near_zone_interactions_const(matrices,mesh)
  end if
  if(mesh%order == 1) then 
  call compute_near_zone_interactions_lin(matrices,mesh)
  end if

  print*,' Compute T-matrix...'
  call compute_T_matrix(matrices, mesh, matrices%Nmaxs(ii), matrices%Taa, &
  matrices%Tab, matrices%Tba, matrices%Tbb)

  matrices%Taai(1:nm,1:nm,ii) = matrices%Taa
  matrices%Tabi(1:nm,1:nm,ii) = matrices%Tab
  matrices%Tbai(1:nm,1:nm,ii) = matrices%Tba
  matrices%Tbbi(1:nm,1:nm,ii) = matrices%Tbb
 end do
end if

end subroutine calc_T

!******************************************************************************

subroutine allocate_T(i, matrices)
type (data) :: matrices
integer :: Nmax, i

if(allocated(matrices%Taa))then
 deallocate(matrices%Taa,matrices%Tba,matrices%Tab,matrices%Tbb)
end if

Nmax = matrices%Nmaxs(i)
allocate(matrices%Taa((Nmax+1)**2-1, (Nmax+1)**2-1))
allocate(matrices%Tab((Nmax+1)**2-1, (Nmax+1)**2-1))
allocate(matrices%Tba((Nmax+1)**2-1, (Nmax+1)**2-1))
allocate(matrices%Tbb((Nmax+1)**2-1, (Nmax+1)**2-1))

end subroutine allocate_T

!******************************************************************************

subroutine allocate_Ti(matrices)
type (data) :: matrices
integer :: Nmax

Nmax = maxval(matrices%Nmaxs)

allocate(matrices%Taai((Nmax+1)**2-1, (Nmax+1)**2-1,matrices%bars))
allocate(matrices%Tabi((Nmax+1)**2-1, (Nmax+1)**2-1,matrices%bars))
allocate(matrices%Tbai((Nmax+1)**2-1, (Nmax+1)**2-1,matrices%bars))
allocate(matrices%Tbbi((Nmax+1)**2-1, (Nmax+1)**2-1,matrices%bars))

matrices%Taai(:,:,:) = dcmplx(0.0,0.0)
matrices%Tabi(:,:,:) = dcmplx(0.0,0.0)
matrices%Tbai(:,:,:) = dcmplx(0.0,0.0)
matrices%Tbbi(:,:,:) = dcmplx(0.0,0.0)

end subroutine allocate_Ti


!******************************************************************************

subroutine mie_T_matrix(matrices, mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: i, j, ci, nm, Nmax
real(dp) :: k, ka
complex(dp) :: eps_r, m, r1
complex(dp), dimension(:), allocatable :: a_n, b_n, j0, j1, h0, j0d, j1d, h0d

do i = 1,size(mesh%ki)
 if(allocated(a_n)) deallocate(a_n,b_n,j0,j1,h0,j0d,j1d,h0d)
 
! write(*,'(3(A,I0))') ' Step ', i, '/', size(mesh%ki), ''
 
 k = mesh%ki(i)
 ka = dcmplx(k*mesh%a)
 Nmax = matrices%Nmaxs(i) 
 nm = (Nmax+1)**2-1 

 allocate(a_n(Nmax), b_n(Nmax), j0(Nmax), j1(Nmax), h0(Nmax), j0d(Nmax), &
 j1d(Nmax), h0d(Nmax))
 m = matrices%refr
 r1 = ka*m

 call sbesseljd(Nmax,dcmplx(ka),j0,j0d)
 call sbesseljd(Nmax,r1,j1,j1d)
 call shankeljd(Nmax,dcmplx(ka),h0,h0d)
 
 do j = 1,Nmax
  b_n(j) = -( j1d(j)*j0(j) - m*j0d(j)*j1(j) ) / ( j1d(j)*h0(j) - m*h0d(j)*j1(j) )
  a_n(j) = -( j0d(j)*j1(j) - m*j1d(j)*j0(j) ) / ( h0d(j)*j1(j) - m*j1d(j)*h0(j) )
 end do

 do j = 1,nm
  ci = floor(sqrt(dble(j)))
  matrices%Taai(j,j,i) = a_n(ci)
  matrices%Tbbi(j,j,i) = b_n(ci)
 end do
end do

end subroutine mie_T_matrix

!******************************************************************************

subroutine ori_ave_T(matrices,mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: i, ii, nm, sz, n, m, Nmax
complex(dp), dimension(:,:), allocatable :: Taa, Tab, Tba, Tbb, &
Taa_av, Tab_av, Tba_av, Tbb_av
complex(dp), dimension(:), allocatable :: tnaa, tnab, tnba, tnbb

sz = size(mesh%ki)
do i = 1,sz
 Nmax = matrices%Nmaxs(i)
 nm = (Nmax+1)**2-1 
 mesh%k = mesh%ki(i)

 if(allocated(Taa)) deallocate(Taa,Tab,Tba,Tbb,Taa_av,Tab_av,Tba_av,&
 Tbb_av,tnaa,tnab,tnba,tnbb)

 allocate(Taa(nm,nm),Tab(nm,nm),Tba(nm,nm),Tbb(nm,nm),&
 Taa_av(nm,nm),Tab_av(nm,nm),Tba_av(nm,nm),Tbb_av(nm,nm))
 allocate(tnaa(Nmax),tnab(Nmax),tnba(Nmax),tnbb(Nmax))

 Taa = matrices%Taai(1:nm,1:nm,i)
 Tab = matrices%Tabi(1:nm,1:nm,i) 
 Tba = matrices%Tbai(1:nm,1:nm,i) 
 Tbb = matrices%Tbbi(1:nm,1:nm,i) 

 tnaa = dcmplx(0d0)
 tnab = dcmplx(0d0)
 tnba = dcmplx(0d0)
 tnbb = dcmplx(0d0)
 ii = 1
 
 do n = 1,Nmax
  do m = -n,n
   tnaa(n) = tnaa(n) + (Taa(ii,ii))/(2*n+1)
   tnab(n) = tnab(n) + (Tab(ii,ii))/(2*n+1)
   tnba(n) = tnba(n) + (Tba(ii,ii))/(2*n+1)
   tnbb(n) = tnbb(n) + (Tbb(ii,ii))/(2*n+1)
   ii = ii + 1
  end do
 end do

 Taa_av = dcmplx(0d0)
 Tab_av = dcmplx(0d0)
 Tba_av = dcmplx(0d0)
 Tbb_av = dcmplx(0d0)

 ii = 1
 do n = 1,Nmax
  do m = -n,n
   Taa_av(ii,ii) = tnaa(n)
   Tab_av(ii,ii) = tnab(n)
   Tba_av(ii,ii) = tnba(n)
   Tbb_av(ii,ii) = tnbb(n)
   ii = ii + 1
  end do
 end do

 matrices%Taai(1:nm,1:nm,i) = Taa_av
 matrices%Tabi(1:nm,1:nm,i) = Tab_av
 matrices%Tbai(1:nm,1:nm,i) = Tba_av
 matrices%Tbbi(1:nm,1:nm,i) = Tbb_av
end do

end subroutine ori_ave_T

!******************************************************************************

subroutine compute_T_matrix(matrices, mesh, Nmax, Taa, Tab, Tba, Tbb)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: k
integer :: Nmax, nm

complex(dp) :: mat(3*mesh%N_tet,2*((Nmax+1)**2-1))
complex(dp) :: T_mat(2*((Nmax+1)**2-1), 2*((Nmax+1)**2-1))
complex(dp), dimension((Nmax+1)**2-1,(Nmax+1)**2-1) :: Taa, Tab, Tba, Tbb
complex(dp) :: sc

k = mesh%k

! Compute transformations
call vswf2constant(mesh, dcmplx(k), Nmax, mat)

do nm = 1,size(mat,2)
 matrices%rhs = mat(:,nm)
 call gmres(matrices,mesh)
 T_mat(:,nm) = matmul(transpose(conjg(mat)),matrices%x)
 call print_bar(nm,size(mat,2))
end do

nm = (Nmax+1)**2-1
sc = dcmplx(0.0d0,k**3.0d0)
Taa = T_mat(1:nm,1:nm) * sc
Tab = T_mat(1:nm, nm+1:2*nm) * sc
Tba = T_mat(nm+1:2*nm,1:nm) * sc
Tbb = T_mat(nm+1:2*nm,nm+1:2*nm) *sc

end subroutine compute_T_matrix

!*******************************************************************************

subroutine scattered_fields(matrices,E,p,q,p90,q90,ii)
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

if(.not. allocated(p)) allocate(p(nm),q(nm),p90(nm),q90(nm))

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

end subroutine scattered_fields

!******************************************************************************

function vsh(n, m, thet, ph) result(BCP)
integer :: n, m,  mm
real(dp) :: r, thet, ph,theta, phi, q
complex(dp) :: BCP(9), P(3), B(3), C(3), Y, Y1, Y2
real(dp), dimension(:), allocatable :: L, L1, L2

theta = max(1d-7, thet)
phi = max(1d-7, ph)

allocate(L(n+1),L1(n+2),L2(n))

call legendre2(n,cos(theta),L)  
call legendre2(n+1,cos(theta),L1)
call legendre2(n-1,cos(theta),L2)

q=(sqrt(n*(n+1.0d0)))/((n*2d0+1.0d0)*sin(theta));

mm = abs(m)

! Unnormalized complex scalar spherical harmonics
Y=L(mm+1)*exp(dcmplx(0.0, m*phi));
Y1=L1(mm+1)*exp(dcmplx(0.0, m*phi));

if(mm == n) then
   Y2 = dcmplx(0.0,0.0)
else 
   Y2 = L2(mm+1)*exp(dcmplx(0.0, m*phi)) 
end if

! vector spherical harmonics
P(:) = dcmplx(0.0,0.0)
P(1) = Y

Y1=Y1*((n-mm+1.0d0)/(n+1.0d0))

Y2=Y2*(dble(n+mm)/dble(n))

B(:) = dcmplx(0.0,0.0)
B(2) = Y1-Y2
B(3)=((dcmplx(0.0, m*(2*n+1.0)))/(n*(n+1.0)))*Y

B = B*q

C(:) = dcmplx(0.0,0.0)
C(2) = B(3) 
C(3) = -B(2)

BCP = [B, C, P]

end function vsh

end module T_matrix
