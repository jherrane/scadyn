module common
implicit none

! Contains:
!  PARAMETERS
!  TYPES
!  LOWFUNCTIONS
!  LOWSUBROUTINES

! PARAMETERS ******************************************************************
!******************************************************************************

integer, parameter   :: dp        = selected_real_kind(15, 307)
integer, parameter   :: nthreads  = 24

real(dp), parameter  :: pi        = 3.141592653589793d0
real(dp), parameter  :: epsilon   = 8.854187817d-12
real(dp), parameter  :: mu        = pi*4d-7
real(dp), parameter  :: mu_B      = 9.2700968d-24
real(dp), parameter  :: hbar      = 1.054572d-34
real(dp), parameter  :: k_b       = 1.38064852d-23
real(dp), parameter  :: cc        = 2.99792458d8
real(dp), parameter  :: b         = 2.8977729d-3

complex(dp), parameter              :: i1  = dcmplx(0d0,1d0)
real(dp), dimension(3), parameter   :: e_1 = [1d0, 0d0, 0d0]
real(dp), dimension(3), parameter   :: e_2 = [0d0, 1d0, 0d0]
real(dp), dimension(3), parameter   :: e_3 = [0d0, 0d0, 1d0]

integer :: seedling           = 0
integer :: use_mie            = 0
integer :: run_test           = 0
integer :: int_mode           = 0
integer :: beam_shape         = 0
integer :: calc_extra_torques = 0
integer :: p                  = 0
integer :: l                  = 0

! TYPES ***********************************************************************
!******************************************************************************

type mesh_struct
   complex(dp), dimension( :, : ), allocatable :: params
   complex(dp), dimension( : ), allocatable :: param

   real(dp), dimension( :, : ), allocatable  ::  coord, nodes, P
   real(dp), dimension( : ), allocatable ::  ki, w, radius
   real(dp), dimension( 3 )  ::  min_coord
   real(dp)  ::  delta, k, box_delta, tol, grid_size, rho, a, V, mass, CM(3), &
   I(3,3), maxrad

   integer, dimension(:,:), allocatable :: etopol, etopol_box, tetras, etopol_edges, edges
   integer :: Nx, Ny, Nz, N_cubes, N_tet, N_tet_cube, M_ex, Nx_cube, Ny_cube, &
   Nz_cube, M1_loc, M2_loc, N1_loc, N2_loc, restart, maxit, near_zone, order, &
   N_node, is_mesh

   character(LEN=80) :: meshname = 'shape.h5'
   character(LEN=80) :: projector = 'pfft'
end type mesh_struct

type data
   complex, dimension(:,:,:), allocatable :: sp_mat
   complex(dp), dimension(:,:,:), allocatable :: Fg, listS, listSx, listSy, &
   listSz, Taai, Tbbi, Tabi, Tbai, Ai
   complex(dp), dimension(:,:), allocatable :: S, Sx, Sy, Sz, E_field, &
   Taa, Tbb, Tab, Tba, S_loc, Sx_loc, Sy_loc, Sz_loc, as, bs
   complex(dp), dimension(:), allocatable ::  x, Ax, rhs, rcs
   complex(dp), dimension(3) :: E0, E90, force, torque, E0hat, E90hat
   complex( 8 ), dimension( :, : ), allocatable :: rotDs, rotD90s, rotYs, rotXs

   real(dp), dimension(:,:,:), allocatable :: RRR
   real(dp), dimension(:,:), allocatable :: field_points, T, www
   real(dp), dimension(:), allocatable :: E_rel

   real(dp), dimension(3,1000) :: x_buf, v_buf, w_buf, J_buf, N_buf, F_buf
   real(dp), dimension(3,3,1000) :: R_buf
   real(dp), dimension(1,1000) :: t_buf

   real(dp), dimension(4) :: q, qn
   real(dp), dimension(3,3) :: R, Rn, Rkt, P, I, I_inv, R_init, R90_init, Rexp, R_al
   real(dp), dimension(3) :: khat, w, x_CM, v_CM, N, wn, xn, vn, J, F, Ip, CM,&
   dw, k_orig, E0_orig, E90_orig, Q_t, Q_f, B
   real(dp) ::  khi_0, rot_max, lambda1, lambda2, temp, dt0, dt, tt, &
   E, refr, refi, tol_m, M1, M3, B_psi

   integer, dimension(:,:,:), allocatable :: listindS, indDs, indD90s, indXs, indYs
   integer, dimension(:,:), allocatable :: indS, T_ind, indS_loc, sp_ind
   integer, dimension(:), allocatable :: S_tet_loc, Nmaxs
   integer :: Nmax, polarization, bars, Tmat, it_max, which_int, &
   is_aggr, whichbar, buffer, it_log, N_mean, is_aligned, alignment_found, it_stop, &
   singleT

   character(len=38) :: waves = 'band'
   character(len=38) :: paramsfile = 'params.in'
   character(LEN=80) :: tname = 'T.h5'
   character(LEN=80) :: out = 'log'
   character(LEN=80) :: mueller = 'mueller'
   character(len=8)  :: mueller_mode = 'none'

   real(dp)          :: B_len = 0d0
end type data

type data_struct
   integer :: Nmax, ind_a1, ind_a2, ind_b1, ind_b2
   integer :: Tmat_ind, ifT
   real(dp) :: euler_angles(3)
   real(dp) :: cp(3), r
   complex(dp) :: eps_r
   complex(dp), dimension(:,:), allocatable :: A, B
end type data_struct


contains

! LOWFUNCTIONS ****************************************************************
!******************************************************************************

function crossRR(a,b) result(c)
real(dp), intent(in) :: a(3), b(3)
real(dp) :: c(3)
c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = -a(1)*b(3) + a(3)*b(1)
c(3) = a(1)*b(2) - a(2)*b(1)

end function crossRR

!******************************************************************************

function crossRC(a,b) result(c)
complex(dp), intent(in) :: b(3)
real(dp), intent(in) :: a(3)
complex(dp) :: c(3)
c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = -a(1)*b(3) + a(3)*b(1)
c(3) = a(1)*b(2) - a(2)*b(1)

end function crossRC

!******************************************************************************
function crossCC(a,b) result(c)
complex(dp), intent(in) :: a(3), b(3)
complex(dp) :: c(3)
c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = -a(1)*b(3) + a(3)*b(1)
c(3) = a(1)*b(2) - a(2)*b(1)

end function crossCC

!******************************************************************************

function norm(r, rp) result(c)
implicit none
real(dp), dimension(3), intent(in) :: r, rp
real(dp) :: c
c = sqrt((r(1)-rp(1))**2 + (r(2)-rp(2))**2 + (r(3)-rp(3))**2)
!norm=sqrt(dot_product(r-rp,r-rp))
end function norm

!******************************************************************************

function vlen(r) result(c)
implicit none
real(dp), dimension(3), intent(in) :: r
real(dp) :: c
c = sqrt((r(1))**2 + (r(2))**2 + (r(3))**2)
!norm=sqrt(dot_product(r-rp,r-rp))
end function vlen

!******************************************************************************

function tet_area(coord) result(c)
implicit none
real(dp), dimension(3,4), intent(in) :: coord
real(dp) :: c, tri1, tri2, tri3, tri4

tri1 = vlen(crossRR(coord(:,2)-coord(:,1),coord(:,3)-coord(:,1))) / 2d0
tri2 = vlen(crossRR(coord(:,3)-coord(:,4),coord(:,2)-coord(:,4))) / 2d0
tri3 = vlen(crossRR(coord(:,4)-coord(:,3),coord(:,1)-coord(:,3))) / 2d0
tri4 = vlen(crossRR(coord(:,1)-coord(:,2),coord(:,4)-coord(:,2))) / 2d0

c = tri1+tri2+tri3+tri4
!norm=sqrt(dot_product(r-rp,r-rp))
end function tet_area

!******************************************************************************

function tetra_volume(tet_coord) result(vol)
real(dp), intent(in) :: tet_coord(:,:)
real(dp) :: p41(3), p21(3), p31(3), cr(3)
real(dp) :: dt, vol
p41 = tet_coord(:,4) - tet_coord(:,1)
p21 = tet_coord(:,2) - tet_coord(:,1)
p31 = tet_coord(:,3) - tet_coord(:,1)

cr = crossRR(p21,p31)
dt = dot_product(p41,cr)
vol = sqrt(dt*dt)/6d0

end function tetra_volume

!******************************************************************************

function Gr(rf,rp,k) result(Green)
real(dp), intent(in) :: rf(3), rp(3)
real(dp), intent(in) :: k
real(dp) :: R
complex(dp) :: Green

R=norm(rf,rp)
Green = cdexp(dcmplx(0.0d0,k*R)) / dcmplx((4d0*pi*R))
end function Gr

!******************************************************************************

function grad_G(rf,rp,k) result(gradG)
real(dp), intent(in) :: rf(3), rp(3), k
complex(dp) :: gradG(3)
real(dp) :: R

R=norm(rf,rp)
gradG = cdexp(dcmplx(0.0d0,k*R)) / dcmplx((4d0*pi*R**3)) * dcmplx(-1.0d0,k*R) *&
 dcmplx(rf-rp)

end function grad_G

!******************************************************************************

function tetra_face(tet_coord, face) result(tri_coord)
real(dp), intent(in) :: tet_coord(3,4)
integer, intent(in) :: face
real(dp) :: tri_coord(3,3)

if(face == 1) then
 tri_coord(:,1) = tet_coord(:,2)
 tri_coord(:,2) = tet_coord(:,3)
 tri_coord(:,3) = tet_coord(:,4)
end if
if(face == 2) then
 tri_coord(:,1) = tet_coord(:,1)
 tri_coord(:,2) = tet_coord(:,4)
 tri_coord(:,3) = tet_coord(:,3)
end if
if(face == 3) then
 tri_coord(:,1) = tet_coord(:,1)
 tri_coord(:,2) = tet_coord(:,2)
 tri_coord(:,3) = tet_coord(:,4)
end if
if(face == 4) then
 tri_coord(:,1) = tet_coord(:,1)
 tri_coord(:,2) = tet_coord(:,3)
 tri_coord(:,3) = tet_coord(:,2)
end if

end function tetra_face

!******************************************************************************

function tri_n_vectors(tri_coord) result(nvec)
real(dp), intent(in) :: tri_coord(3,3)
real(dp) :: nvec(3)
real(dp) :: a1(3), a2(3), n(3), len_n

a1 = tri_coord(:,2) - tri_coord(:,1)
a2 = tri_coord(:,3) - tri_coord(:,1)

n=crossRR(a1,a2)
len_n = sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))

nvec = n/len_n

end function tri_n_vectors

!******************************************************************************

function tetra_n_vectors(tet_coord) result(n_vectors)
real(dp), intent(in) :: tet_coord(3,4)
real(dp) :: n_vectors(3,4), tri_coord(3,3), n_vec(3), test(3), tt
integer :: n

do n = 1,4
 tri_coord = tetra_face(tet_coord,n)
 n_vec = tri_n_vectors(tri_coord)

 test = tri_coord(:,1) - tet_coord(:,n)
 tt = dot_product(n_vec,test)
 if(tt<0) then
 print*,'Something is wrong'
 end if

 n_vectors(:,n) = n_vec

end do

end function tetra_n_vectors

!******************************************************************************

function get_tetra_vol(mesh) result(volume)
type(mesh_struct)             ::  mesh
real(dp)                      ::  volume, V, totV
integer                       ::  i1
real( dp ), dimension( 3 )    ::  p0, p1, p2, p3

totV             =  0.d0

do i1 = 1, mesh%N_tet
 ! Get vertices
 p0              =  mesh%coord( :, mesh%etopol( 1, i1 ) )
 p1              =  mesh%coord( :, mesh%etopol( 2, i1 ) )
 p2              =  mesh%coord( :, mesh%etopol( 3, i1 ) )
 p3              =  mesh%coord( :, mesh%etopol( 4, i1 ) )

 V               =  abs( dot_product( ( p0 - p3 ), &
                    crossRR( ( p1 - p3 ), ( p2 - p3 ) ) ) ) / 6d0

 totV            =  totV + V
end do

volume           =  totV

end function get_tetra_vol

!******************************************************************************

function vec2arr(mesh,vec) result(arr)
type (mesh_struct) :: mesh
complex(dp), dimension(:), intent(in) :: vec
complex(dp), dimension(:,:,:), allocatable :: arr
integer :: las, m1, l1, k1

las = 1

allocate(arr(2*mesh%Nx, 2*mesh%Ny, 2*mesh%Nz))

arr(:,:,:) = dcmplx(0.0d0, 0.0d0)

do m1 = 1, mesh%Nz
 do l1 = 1, mesh%Ny
  do k1 = 1, mesh%Nx
   arr(k1, l1, m1) = vec(las)
   las = las + 1
  end do
 end do
end do

end function vec2arr

!******************************************************************************

function arr2vec(mesh,arr) result(vec)
type (mesh_struct) :: mesh
complex(dp), dimension(:), allocatable :: vec
complex(dp), dimension(:,:,:), intent(in) :: arr
integer :: las, m1, l1, k1

las = 1
allocate(vec(mesh%Nx*mesh%Ny*mesh%Nz))

do m1 = 1, mesh%Nz
 do l1 = 1, mesh%Ny
  do k1 = 1, mesh%Nx
   vec(las) = arr(k1,l1,m1)
   las = las + 1
  end do
 end do
end do

!scale IFFT
vec = vec / (8d0*mesh%Nx*mesh%Ny*mesh%Nz)

end function arr2vec

!******************************************************************************

function calc_det3(mat) result(det)
real(dp), dimension(3,3), intent(in) :: mat
real(dp) :: det

det = mat(1,1) * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
-  mat(1,2) * (mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
+  mat(1,3) * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))

end function calc_det3

!******************************************************************************

subroutine gradshape(rot,dN,T_coord)
real(dp), dimension(3,4), intent(in) :: T_coord
real(dp), dimension(3,6) :: rot
real(dp), dimension(3,4) :: dN
real(dp) :: J(3,3), detJ, J0(3,3)

J(:,1) = T_coord(:,2) - T_coord(:,1)
J(:,2) = T_coord(:,3) - T_coord(:,1)
J(:,3) = T_coord(:,4) - T_coord(:,1)

detJ = calc_det3(J)

if(detJ <= 0) then
 print*,'Something is wrong (detJ <= 0)'
 stop
end if

J0(1,1) = J(2,2)*J(3,3) - J(2,3)*J(3,2);
J0(2,1) = -(J(1,2)*J(3,3) - J(1,3)*J(3,2));
J0(3,1) = (J(1,2)*J(2,3) - J(1,3)*J(2,2));

J0(1,2) = -(J(2,1)*J(3,3) - J(2,3)*J(3,1));
J0(2,2) = (J(1,1)*J(3,3) - J(1,3)*J(3,1));
J0(3,2) = -(J(1,1)*J(2,3) - J(1,3)*J(2,1));

J0(1,3) = J(2,1)*J(3,2) - J(2,2)*J(3,1);
J0(2,3) = -(J(1,1)*J(3,2) - J(1,2)*J(3,1));
J0(3,3) = (J(1,1)*J(2,2) - J(1,2)*J(2,1));

dN(:,2) = J0(:,1) / detJ
dN(:,3) = J0(:,2) / detJ
dN(:,4) = J0(:,3) / detJ
dN(:,1) = - dN(:,2) - dN(:,3) - dN(:,4)

rot(:,1) = 2d0*crossRR(dN(:,1),dN(:,2))
rot(:,2) = 2d0*crossRR(dN(:,1),dN(:,3))
rot(:,3) = 2d0*crossRR(dN(:,1),dN(:,4))
rot(:,4) = 2d0*crossRR(dN(:,2),dN(:,3))
rot(:,5) = 2d0*crossRR(dN(:,3),dN(:,4))
rot(:,6) = 2d0*crossRR(dN(:,2),dN(:,4))

end subroutine gradshape

!******************************************************************************

function sph_unit_vectors(theta,phi) result(vec)
real(dp), intent(in) :: theta, phi
real(dp) :: vec(3,3)

vec(:,1) = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)] ! r_hat
vec(:,2) = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)] ! theta_hat
vec(:,3) = [-sin(phi), cos(phi), dble(0.0)] ! phi_hat

end function sph_unit_vectors

!******************************************************************************

function real_outer_product(a,b) result(vec)
real(dp), dimension(:) :: a, b
real(dp), dimension(:,:), allocatable :: vec
integer :: i1, Na, Nb

Na = size(a)
Nb = size(b)

allocate(vec(Na,Nb))

do i1 = 1,Na
 vec(i1,:) = b*a(i1)
end do

end function real_outer_product

!******************************************************************************

function rotation(x, axis, angle) result(xp)
real(dp) :: x(3), angle, xp(3)
integer :: axis

if(axis == 1) then
 xp(1) = x(1)
 xp(2) = cos(angle)*x(2) - sin(angle)*x(3)
 xp(3) = sin(angle)*x(2) + cos(angle)*x(3)
end if

if(axis == 2) then
 xp(1) = cos(angle)*x(1) + sin(angle)*x(3)
 xp(2) = x(2)
 xp(3) = -sin(angle)*x(1) + cos(angle)*x(3)
end if

if(axis == 3) then
 xp(1) = cos(angle)*x(1) - sin(angle)*x(2)
 xp(2) = sin(angle)*x(1) + cos(angle)*x(2)
 xp(3) = x(3)
end if

end function rotation

!******************************************************************************
! Transform Tait-Bryan XYZ angles to a rotation matrix
function rotation_matrix(a,b,g) result(rot)
real(dp) :: a, b, g
real(dp) :: rot(3,3)

rot(1,1) = cos(b)*cos(g);
rot(1,2) = cos(a)*sin(g) + sin(a)*sin(b)*cos(g);
rot(1,3) = sin(a)*sin(g) - cos(a)*sin(b)*cos(g);

rot(2,1) = -cos(b)*sin(g);
rot(2,2) = cos(a)*cos(g) - sin(a)*sin(b)*sin(g);
rot(2,3) = sin(a)*cos(g) + cos(a)*sin(b)*sin(g);

rot(3,1) = sin(b);
rot(3,2) = -sin(a)*cos(b);
rot(3,3) = cos(a)*cos(b);

!rot = transpose(rot) ! Important transform to column-major order!

end function rotation_matrix

!******************************************************************************

function R_aa(axis,angle) result(rot)
real(dp) :: axis(3), angle, rot(3,3)

rot = dcos(angle)*eye(3) + dsin(angle)*cross_skew(axis) + (1d0-dcos(angle))*&
real_outer_product(axis, axis)

end function R_aa

!******************************************************************************

function halton_seq(index, base) result(num)
integer :: base, index, i
real(dp) :: num, f
num = 0.d0

f = 1.0d0/dble(base)

i = index

do while(i>0)
 num = num + f * dble(mod(i,base))
 i = floor(dble(i)/dble(base))
 f = f / dble(base)
end do

end function halton_seq

!******************************************************************************

function neighbour_tet(T_nodes, B_nodes) result(num)
integer :: T_nodes(4), B_nodes(4)
integer :: i, j, num
num = 0
do i=1,4
 do j=1,4
  if(T_nodes(i) == B_nodes(j)) then
   num = 1
  end if
 end do
end do

end function neighbour_tet

!******************************************************************************

function factorial(a) result(c)
integer :: i1, a
real(dp) :: c, b

b=1.0d0
do i1 = 1,a
 b = b * i1
end do
c=dble(b)

end function factorial

!******************************************************************************

function cart2sph(x) result(vec)
real(dp), intent(in) :: x(3)
real(dp) :: vec(3)

vec(1) = sqrt(x(1)**2 + x(2)**2 + x(3)**2) ! r
vec(2) = acos(x(3)/vec(1)) ! theta
vec(3) = atan2(x(2),x(1))

end function cart2sph

!******************************************************************************

function cart2sph2(x) result(vec)
real(dp), intent(in) :: x(3)
real(dp) :: vec(3)

vec(1) = sqrt(x(1)**2 + x(2)**2 + x(3)**2) ! r
vec(2) = atan2(sqrt(x(1)**2+x(2)**2),x(3))
vec(3) = atan2(x(2),x(1))

end function cart2sph2

!******************************************************************************

function sph2cart(r,theta,phi) result(x)
real(dp), intent(in) :: r, theta, phi
real(dp) :: x(3)

x(1) = r*sin(theta)*cos(phi)
x(2) = r*sin(theta)*sin(phi)
x(3) = r*cos(theta)

end function sph2cart

!******************************************************************************

function sph2cart_vec(theta, phi, vec) result(vec2)
complex(dp), intent(in) :: vec(3)
real(dp), intent(in) :: theta, phi
complex(dp) :: vec2(3)
real(dp) :: H(3,3)

H(1,:) = [sin(theta)*cos(phi), cos(theta)*cos(phi), -sin(phi)]
H(2,:) = [sin(theta)*sin(phi), cos(theta)*sin(phi), cos(phi)]
H(3,:) = [cos(theta), -sin(theta), dble(0.0)];

vec2 = matmul(dcmplx(H),vec)

end function sph2cart_vec

!******************************************************************************

function circ2cart_vec(vec) result(vec2)
complex(dp), intent(in) :: vec(3)
complex(dp) :: vec2(3)
complex(dp) :: H(3,3)

H(1,:) = [dcmplx(1d0)/sqrt(2d0),-i1/sqrt(2d0),dcmplx(0d0)]
H(2,:) = [dcmplx(0d0),dcmplx(0d0),dcmplx(1d0)]
H(3,:) = [dcmplx(1d0)/sqrt(2d0),i1/sqrt(2d0),dcmplx(0d0)];

vec2 = matmul(H,vec)

end function circ2cart_vec

!******************************************************************************

function binomial(n,k) result(c)
integer :: n, k, i1
real(dp) :: c

c = 1.0d0
do i1 = 1,k
 c = c * (n + 1.0d0 - i1)/i1
end do

end function binomial

!******************************************************************************

function rotation_angles(x) result(vec)
real(dp), intent(in) :: x(3)
real(dp) :: vec(3)

vec(1) = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
vec(3) = atan2(x(2),x(1))     ! phi
vec(2) = atan2(x(1),x(3))

end function rotation_angles

!******************************************************************************

function truncation_order(ka) result(Nmax)
real(dp) :: ka,lim
integer :: Nmax

lim = 4d0
if(ka > 1) then
 Nmax = floor(ka + lim* (ka)**(1.0d0/3.0d0))
else
 Nmax = floor(lim)
end if

end function truncation_order

!******************************************************************************

function truncation_order2(ka) result(Nmax)
real(dp) :: ka
integer :: Nmax

if(ka > 1) then
 Nmax = floor(ka + 3.0d0* (ka)**(1.0d0/3.0d0))
else
 Nmax = 4
end if

end function truncation_order2

!******************************************************************************

function cross_skew(v) result(Vx)
real(dp) :: v(3), Vx(3,3)

Vx = reshape([0.d0, v(3), -v(2), -v(3), 0.d0, v(1), &
v(2), -v(1), 0.d0],[3,3])

end function cross_skew

!******************************************************************************

function eye(dim) result(I)
real(dp), dimension(:,:), allocatable :: I
integer :: i1, dim

allocate(I(dim, dim))
I = 0.d0
forall(i1 = 1:dim) I(i1,i1) = 1.d0

end function eye

!******************************************************************************
!    Discover Euler angle vector from 3x3 matrix

!    Uses the conventions above.

!    Parameters
!    ----------
!    M : array-like, shape (3,3)
!    cy_thresh : None or scalar, optional
!       threshold below which to give up on straightforward arctan for
!       estimating x rotation.  If None (default), estimate from
!       precision of input.

!    Returns
!    -------
!    z : scalar
!    y : scalar
!    x : scalar
!       Rotations in radians around z, y, x axes, respectively

!    Notes
!    -----
!    If there was no numerical error, the routine could be derived using
!    Sympy expression for z then y then x rotation matrix, which is::

!      [                       cos(y)*cos(z),                       -cos(y)*sin(z),         sin(y)],
!      [cos(x)*sin(z) + cos(z)*sin(x)*sin(y), cos(x)*cos(z) - sin(x)*sin(y)*sin(z), -cos(y)*sin(x)],
!      [sin(x)*sin(z) - cos(x)*cos(z)*sin(y), cos(z)*sin(x) + cos(x)*sin(y)*sin(z),  cos(x)*cos(y)]

!    with the obvious derivations for z, y, and x

!       z = atan2(-r12, r11)
!       y = atan2(sin(r13),cy)
!       x = atan2(-r23, r33)

!    Problems arise when cos(y) is close to zero, because both of::

!       z = atan2(cos(y)*sin(z), cos(y)*cos(z))
!       x = atan2(cos(y)*sin(x), cos(x)*cos(y))

!    will be close to atan2(0, 0), and highly unstable.

!    The ``cy`` fix for numerical instability below is from: *Graphics
!    Gems IV*, Paul Heckbert (editor), Academic Press, 1994, ISBN:
!    0123361559.  Specifically it comes from EulerAngles.c by Ken
!    Shoemake, and deals with the case where cos(y) is close to zero:

!    See: http://www.graphicsgems.org/

!    The code appears to be licensed (from the website) as "can be used
!    without restrictions".

function mat2euler(R) result(angles)
real(dp) :: R(3,3), angles(3), tol, cy

tol = 1.d-11

cy = sqrt(R(3,3)*R(3,3)+R(2,3)*R(2,3))
if(cy > tol) then
 angles(3) = atan2(-R(1,2),R(1,1)) ! atan( cos(y)*sin(z)/ (cos(y)*cos(z)) )
 angles(2) = atan2(R(1,3),cy) ! atan( sin(y)/cy )
 angles(1) = atan2( -R(2,3),R(3,3) ) ! atan( cos(y)*sin(x)/(cos(x)*cos(y)) )
else ! cos(y) (close to) zero, so x -> 0.0 (see above)
 ! so r21 -> sin(z), r22 -> cos(z) and
 angles(3) = atan2(R(2,1),R(2,2))
 angles(2) = atan2(R(1,3),cy)
 angles(1) = 0.d0
end if

end function mat2euler

!******************************************************************************

function rotate_a_to_b(a,b) result(R)
real(dp), dimension(3), intent(in) :: a, b
real(dp), dimension(3,3) :: R
real(dp), dimension(3) :: ax, xx, ei
real(dp) :: tolerance, theta, ax_norm, s, c, t, x, y, z

tolerance = 1.d-6

theta = dot_product(a,b)
if(theta < 0.d0) then
 theta = dacos(max(theta,-1.d0))
else
 theta = dacos(min(theta,1.d0))
end if

ax = crossRR(a, b)
ax_norm = vlen(ax)
if(theta < tolerance) then
 R = eye(3)
 return
else if(pi-theta < tolerance)then
 ei = 0.d0
 ei(minloc(abs(a))) = 1.d0
 xx = crossRR(a,ei)
 ax = xx/vlen(xx)
! R = eye(3)
! R(3,3) = -1.d0
 !return
else
 ax = ax/ax_norm
end if

s = sin(theta)
c = cos(theta)
t = 1d0 - c

x = ax(1)
y = ax(2)
z = ax(3)

R = reshape([t*x*x + c, t*x*y + s*z, t*x*z - s*y, t*x*y - s*z, &
t*y*y + c, t*y*z + s*x, t*x*z + s*y, t*y*z - s*x, t*z*z + c],[3,3])

end function rotate_a_to_b

!******************************************************************************

function find_Rexp(matrices) result(R)
type (data), intent(in) :: matrices
real(dp), dimension(3,3) :: R
real(dp), dimension(3) :: k, knew

k = matrices%khat
knew = dble((/0.d0,0.d0,1.d0/))

R = rotate_a_to_b(k,knew)

end function find_Rexp

!******************************************************************************

function find_Rexp2(matrices) result(R)
type (data), intent(in) :: matrices
real(dp), dimension(3,3) :: I, R, Vx
real(dp), dimension(3) :: b, a, v
real(dp) :: tolerance, s, c

tolerance = 1.d-6

b = matrices%khat
a = dble((/0.d0,0.d0,1.d0/))
I = eye(3)

c = dot_product(a,b)
v = crossRR(a,b)
s = vlen(v)

Vx = reshape([0.d0,v(3),-v(2),-v(3),0.d0,v(1),v(2),-v(1),0.d0],[3,3])

if(c < tolerance) then
 R = -eye(3)
 return
else
 R = I + Vx + matmul(Vx,Vx)*(1d0-c)/s**2
end if

end function find_Rexp2

!******************************************************************************
! Rodrigues' rotation formula
function rodr( x ) result( expx )
real( dp ) :: x(3), expx(3,3), theta, xhat(3,3)

theta = vlen(x)
xhat = hat(x)
expx = eye(3) + sin(theta)*xhat/theta + (1d0-cos(theta))*matmul(xhat,xhat)/(theta**2)

end function rodr

!******************************************************************************
! Random rotation matrix
function rand_rot() result ( R )
real( dp ) :: R(3,3), theta, phi, vec(3), u(3)
integer :: i, size
integer, allocatable :: seed(:)

! Overly complicated random number hassle... Something may be awry somewhere.
if(seedling==0) then
 call system_clock(i)
else
 i = seedling
end if
call random_seed(size=size)
allocate(seed(size))
seed=i+37*[(i,i=0,size-1)]
call random_seed(put=seed)
call random_number(u)
deallocate(seed)

! After this everything is nice and clean
theta = dacos(2d0*u(1)-1)
phi = 2d0*pi*u(2)
vec = sph2cart(1d0,theta,phi)

R = R_aa(vec,2d0*phi*u(3))

end function rand_rot

!******************************************************************************
! Dericative of the Rodrigues' rotation formula
function drodr( x ) result( dexpx )
real( dp ) :: x(3), dexpx(3,3), theta, xhat(3,3)

theta = vlen(x)
xhat = hat(x)
dexpx = eye(3) +  (1d0-cos(theta))*xhat/theta**2 &
 + (theta-sin(theta))*matmul(xhat,xhat)/(theta**3)

end function drodr

!******************************************************************************
! Cayley formula
function cay( x ) result( cayx )
real( dp ) :: x(3), cayx(3,3), theta, xhat(3,3)

theta = vlen(x)
xhat = hat(x)
cayx = eye(3) + 4d0*xhat/(4d0+theta**2) + 2d0*matmul(xhat,xhat)/(4+theta**2)

end function cay

!******************************************************************************
! Dericative of the Cayley formula
function dcay( x ) result( dcayx )
real( dp ) :: x(3), dcayx(3,3), theta, xhat(3,3)

theta = vlen(x)
xhat = hat(x)
dcayx =  4d0*eye(3)/(4d0+theta**2) + 2d0*xhat/(4d0+theta**2)

end function dcay

!******************************************************************************
! Hatmap R^3 --> so(3)
function hat( v ) result( Vx )
real( dp )                 :: v(3),Vx(3,3)

Vx = reshape([0.d0, v(3), -v(2), -v(3), 0.d0, v(1), &
v(2), -v(1), 0.d0],[3,3])

end function hat

!******************************************************************************

function quat_mult(q1,q2) result(q3)
real(dp), dimension(4) :: q1, q2, q3

q3(1) = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4)
q3(2:4) = q1(1)*q2(2:4) + q2(1)*q1(2:4) + crossRR(q1(2:4),q2(2:4))

end function quat_mult

!******************************************************************************
! Quaternion-vector multiplication
function quat_rotation(q1,v) result(vout)
real(dp), dimension(4) :: q1, q2, q3
real(dp), dimension(3) :: v,vout

q2(1) = 0d0
q2(2:4) = v
q3 = quat_mult(q1,q2)
q3 = quat_mult(q3,quat_conj(q1))
vout = q3(2:4)

end function quat_rotation

!******************************************************************************

function quat_norm(q) result(norm)
real(dp), dimension(4) :: q
real(dp) :: norm

norm = sqrt(dot_product(q,q))

end function quat_norm

!******************************************************************************

function quat_conj(q) result(qc)
real(dp), dimension(4) :: q, qc

qc(1) = q(1)
qc(2:4) = -q(2:4)

end function quat_conj

!******************************************************************************

function quat_inv(q) result(qinv)
real(dp), dimension(4) :: q, qinv

qinv(1) = q(1)/quat_norm(q)
qinv(2:4) = -q(2:4)/quat_norm(q)

end function

!******************************************************************************

function normalize_quat(q) result(qhat)
real(dp), dimension(4):: q, qhat

qhat = q/quat_norm(q)

end function normalize_quat

!******************************************************************************

function quat2mat(q) result(M)
real(dp), dimension(4) :: q
real(dp), dimension(3,3) :: M
real(dp) :: w, x, y, z
w = q(1)
x = q(2)
y = q(3)
z = q(4)

M(1,1) = 1d0 - 2.d0*(y**2 + z**2)
M(1,2) = 2d0*(x*y - w*z)
M(1,3) = 2d0*(x*z + w*y)
M(2,1) = 2d0*(x*y + w*z)
M(2,2) = 1d0 - 2.d0*(x**2 + z**2)
M(2,3) = 2d0*(y*z - w*x)
M(3,1) = 2d0*(x*z - w*y)
M(3,2) = 2d0*(y*z + w*x)
M(3,3) = 1d0 - 2.d0*(x**2 + y**2)

end function quat2mat

!******************************************************************************

function mat2quat(M) result(q)
real(dp) :: M(3,3), q(4), q0, q1, q2, q3

q0 = (1d0+M(1,1)+M(2,2)+M(3,3))/4d0
q1 = (1d0+M(1,1)-M(2,2)-M(3,3))/4d0
q2 = (1d0-M(1,1)+M(2,2)-M(3,3))/4d0
q3 = (1d0-M(1,1)-M(2,2)+M(3,3))/4d0

if(q0<0d0) q0 = 0d0
if(q1<0d0) q1 = 0d0
if(q2<0d0) q2 = 0d0
if(q3<0d0) q3 = 0d0

q0 = sqrt(q0)
q1 = sqrt(q1)
q2 = sqrt(q2)
q3 = sqrt(q3)

if (q0 >= q1 .AND. q0 >= q2 .AND. q0 >= q3) then
 q0 = q0*1d0
 q1 = q1*sign(1.d0,M(3,2)-M(2,3))
 q2 = q2*sign(1.d0,M(1,3)-M(3,1))
 q3 = q3*sign(1.d0,M(2,1)-M(1,2))
else if ( q1 >= q0 .AND. q1 >= q2 .AND. q1 >= q3 ) then
 q0 = q0*sign(1.d0,M(3,2)-M(2,3))
 q1 = q1*1d0
 q2 = q2*sign(1.d0,M(2,1)+M(1,2))
 q3 = q3*sign(1.d0,M(1,3)+M(3,1))
else if ( q2 >= q0 .AND. q2 >= q1 .AND. q2 >= q3 ) then
 q0 = q0*sign(1.d0,M(1,3)-M(3,1))
 q1 = q1*sign(1.d0,M(2,1)+M(1,2))
 q2 = q2*1d0
 q3 = q3*sign(1.d0,M(3,2)+M(2,3))
else if ( q3 >= q0 .AND. q3 >= q1 .AND. q3 >= q2 ) then
 q0 = q0*sign(1.d0,M(2,1)-M(1,2))
 q1 = q1*sign(1.d0,M(3,1)+M(1,3))
 q2 = q2*sign(1.d0,M(3,2)+M(2,3))
 q3 = q3*1d0
end if

 q = [q0, q1, q2, q3]
 q = normalize_quat(q)

end function mat2quat

!******************************************************************************

function rand_sphere() result(vec)
integer, parameter :: seed = 86456
real(dp) :: vec(3), x, y, z, t, r(2)

call random_number(r)

z = 2d0*r(1) - 1
t = 2d0*pi*r(2)
x = dsqrt(1d0-z**2)*dcos(t)
y = dsqrt(1d0-z**2)*dsin(t)
vec = [x,y,z]

end function rand_sphere

!******************************************************************************

function fibonacci_sphere(n,randomize) result(vec)
integer :: n, randomize, i
integer :: size
integer, allocatable :: seed(:)
real(dp) :: vec(3,n), x, y, z, t, r, rnd, offset, incr


rnd = 1d0
if(randomize == 1) then
  ! Overly complicated random number hassle... Something may be awry somewhere.
  call system_clock(i)
  call random_seed(size=size)
  allocate(seed(size))
  seed=i+37*[(i,i=0,size-1)]
  call random_seed(put=seed)
  call random_number(rnd)
  deallocate(seed)
end if

offset = 2d0/n
incr = pi*(3d0-dsqrt(5d0))

do i = 1,n
  y = ((dble(i)*offset)-1d0)+0.5d0*offset
  if(y>1d0) y = 1d0
  r = dsqrt(1d0-y**2)
  t = dmod(i+rnd,dble(n))*incr
  x = r*dcos(t)
  z = r*dsin(t)
  vec(1,i) = x
  vec(2,i) = y
  vec(3,i) = z
end do

end function fibonacci_sphere

!******************************************************************************

function uniform_sphere(n,theta,phi) result(vec)
integer :: n, i, j
integer :: size
integer, allocatable :: seed(:)
real(dp), allocatable :: theta(:),phi(:)
real(dp) :: vec(3,n), x, y, z, t, r, rnd, offset, incr

do j = 1,size(phi,1)
	do i = 1,size(theta,1)
    x = dsin(theta(i))*dcos(phi(j))
    y = dsin(theta(i))*dsin(phi(j))
    z = dcos(theta(i))
    vec(1,i+(j-1)*size(theta,1)) = x
    vec(2,i+(j-1)*size(theta,1)) = y
    vec(3,i+(j-1)*size(theta,1)) = z
	end do
end do

end function uniform_sphere

!******************************************************************************
! Function to calculate a leapfrog update of an quaternion when corresponding
! angular velocity is w. As in [Seelen2016]
function q_worm(w,dt) result(qw)
real(dp), dimension(4) :: qw
real(dp) :: w(3),dt, wlen, sw, cw, wdt2

wlen = vlen(w)
wdt2 = wlen*dt/2
if(wdt2<1.d-10)then
 sw = 0d0
 cw = 1d0
 wlen = 1d0
else
 cw = cos(wdt2)
 sw = sin(wdt2)
end if

qw(1) = cw
qw(2:4) = (sw/wlen)*w
qw = normalize_quat(qw)
end function q_worm

!******************************************************************************
! Rolling mean for purposes when data contains no high peaks. Any peaks will
! affect the mean for a long time, which is not always desirable.
function rolling_mean(N, mean, new_sample) result(avg)
real(dp) :: mean, new_sample, avg
integer :: N

avg = mean - mean/N + new_sample/N

end function rolling_mean

!******************************************************************************
! Tests whether a particle is aligned internally, i.e. one principal axis is
! directed in the direction of angular velocity
function alignment_state(matrices) result(ans)
type(data) :: matrices
real(dp) :: aw, cw, w(3), wnorm
integer :: ans

ans = 0

wnorm = vlen(matrices%w)

if(wnorm>0d0)then
   w = matrices%w/vlen(matrices%w)
else
   w = matrices%w
end if

aw = dabs(dot_product(w,[1d0,0d0,0d0]))
cw = dabs(dot_product(w,[0d0,0d0,1d0]))

matrices%M1 = rolling_mean(matrices%N_mean, matrices%M1, aw)
matrices%M3 = rolling_mean(matrices%N_mean, matrices%M3, cw)

if(1d0-matrices%M1 < matrices%tol_m .OR. 1d0-matrices%M3 < matrices%tol_m) ans = 1

end function alignment_state

!******************************************************************************
function file_exists(fname) result(exists)
logical :: exists
CHARACTER(LEN=80) :: fname

inquire(file=trim(fname), exist=exists)

end function file_exists

! LOWSUBROUTINES **************************************************************
!******************************************************************************

subroutine linmap_tet(P,W, tet_coord, P0, W0)
real(dp), intent(in) :: tet_coord(:,:)
real(dp), intent(in) :: P0(:,:)
real(dp),  intent(in) :: W0(:)

real(dp) :: P(3,size(W0)), W(size(W0))
integer :: n, i1
real(dp) :: N1, N2, N3, N4
n=size(W0)

do i1 = 1,n
 N1 = 1 - P0(1,i1) - P0(2,i1) - P0(3,i1)
 N2 = P0(1,i1)
 N3 = P0(2,i1)
 N4 = P0(3,i1)

 P(:,i1) = tet_coord(:,1) * N1 + tet_coord(:,2) * N2 + tet_coord(:,3) * N3 + tet_coord(:,4) * N4

end do

W=W0 * 6d0 * tetra_volume(tet_coord)

end subroutine linmap_tet

!******************************************************************************

subroutine linmap_tri(P, W, tri_coord, P0, W0)
implicit none
real(dp), intent(in) :: tri_coord(3,3)
real(dp), intent(in) :: P0(:,:)
real(dp), intent(in) :: W0(:)

real(dp) :: P(3,size(W0)), W(size(W0))

real(dp), dimension(3) :: a, b, ww
real(dp) :: ww1
integer i1, n

n=size(W0)

a=tri_coord(:,2)-tri_coord(:,1);
b=tri_coord(:,3)-tri_coord(:,1);

ww(1) = (a(1)*a(1) +  a(2)*a(2) + a(3)*a(3)) * b(1)
ww(2) = (a(1)*a(1) +  a(2)*a(2) + a(3)*a(3)) * b(2)
ww(3) = (a(1)*a(1) +  a(2)*a(2) + a(3)*a(3)) * b(3)

ww1=sqrt((ww(1)*b(1)+ww(2)*b(2)+ww(3)*b(3))-((a(1)*b(1)+a(2)*b(2)+a(3)*b(3))*(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))));

do i1 = 1,n
 P(:,i1) = tri_coord(:,1) + a * P0(1,i1) + b * P0(2,i1)
end do

W = W0 * ww1

end subroutine linmap_tri

!******************************************************************************

subroutine vec2arr2(arr,mesh,vec)
type (mesh_struct), intent(in) :: mesh
complex(dp), dimension(:), intent(in) :: vec
complex(dp), dimension(:,:,:) :: arr
integer :: las, m1, l1, k1, Nx, Ny, Nz

Nx = mesh%Nx
Ny = mesh%Ny
Nz = mesh%Nz

arr = dcmplx(0.0d0, 0.0d0)

!!$omp parallel num_threads(nthreads) default(none) &
!!$omp firstprivate(Nx, Ny, Nz)  &
!!$omp private(m1,l1,k1,las) &
!!$omp shared(arr, vec)
!!$omp do
!!do m1 = 1, 2*Nz
!!   do l1 = 1, 2*Ny
!!      do k1 = 1, 2*Nx
!!         las = k1 + (l1-1)*Nx + (m1-1)*Ny*Nx
!!         arr(k1, l1, m1) = dcmplx(0.0, 0.0)
!         !las = las + 1
!!      end do
!!   end do
!!end do
!!$omp end do
!!$omp end parallel


!$omp parallel num_threads(nthreads) default(none) &
!$omp firstprivate(Nx, Ny, Nz)  &
!$omp private(m1,l1,k1,las) &
!$omp shared(arr, vec)
!$omp do
do m1 = 1, Nz
 do l1 = 1, Ny
  do k1 = 1, Nx
   las = k1 + (l1-1)*Nx + (m1-1)*Ny*Nx
   arr(k1, l1, m1) = vec(las)
   !las = las + 1
  end do
 end do
end do
!$omp end do
!$omp end parallel


end subroutine vec2arr2

!******************************************************************************

subroutine diasym(a,eig)
!**************************************************************
! Calls the LAPACK diagonalization subroutine DSYEV
! input:  a(3,3) = real symmetric matrix to be diagonalized
! output: a(3,3) = orthonormal eigenvectors of a
!         eig(3) = eigenvalues of a in descending order
!**************************************************************
integer lwork, info
real(dp), intent(inout) :: a(3,3)
real(dp), intent(out) :: eig(3)
real(dp) :: work(3*(3+3/2))

lwork=3*(3+3/2)
call dsyev('V','U',3,a,3,eig,work,lwork,info)
if(info>0)then
	print*, "Something wrong with diagonalization! Regards, dsyev (Lapack)"
end if

end subroutine diasym

!******************************************************************************

subroutine linspace(d1,d2,n,grid)

implicit none

integer, intent(in) :: n
real(dp), intent(in) :: d1, d2
real(dp), DIMENSION(n), intent(out) :: grid

integer :: indxi

grid(1) = d1
do indxi= 0,n-2
 grid(indxi+1) = d1+(DBLE(indxi)*(d2-d1))/DBLE(n-1)
end do
grid(n) = d2

end subroutine linspace

!******************************************************************************

subroutine angular_grid(ntheta, nphi, theta, phi)

implicit none

integer, intent(in) :: ntheta, nphi
real(dp), DIMENSION(ntheta*nphi), intent(out) :: theta, phi
real(dp) :: theta_temp(ntheta), phi_temp(nphi)
integer :: i, j, ii

call linspace(pi/ntheta,pi,ntheta,theta_temp)
theta_temp = theta_temp - 0.5d0*pi/ntheta
call linspace(2*pi/nphi,2*pi,nphi,phi_temp)
phi_temp = phi_temp - 2d0*pi/nphi

ii = 0
do j = 1,nphi
  do i = 1,ntheta
   ii = ii + 1
   theta(ii) = theta_temp(i)
   phi(ii) = phi_temp(j)
  end do
end do 

end subroutine angular_grid

!*******************************************************************************

function choose (n, k) result (res)
implicit none
integer, intent (in) :: n
integer, intent (in) :: k
integer :: res

 res = factorial (n) / (factorial (k) * factorial (n - k))
 
end function choose

!*******************************************************************************

function R_theta(matrices, theta) result(R)
type(data) :: matrices
real(dp), dimension(3, 3) ::  R_init, R_thta, R, R_pol
real(dp), dimension(3) :: a_3, e_3, e_2, k, ninit, a_2, a_1
real(dp) :: theta, theta0
integer :: Nang

Nang = 180

k = matrices%k_orig
e_2 = matrices%E0_orig
e_3 = matrices%E90_orig
a_1 = matrices%P(1:3,1)
a_2 = matrices%P(1:3,2)
a_3 = matrices%P(1:3,3)

! Rotation of a_3 to k
theta0 = dacos(dot_product(k, a_3))
if(theta0>1d-6)then
 ninit = crossRR(k,a_3)
 ninit = ninit/vlen(ninit)
 R_init = R_aa(ninit,-theta0)
else
 R_init = eye(3)
!	print*, 'check this out 1'
end if

! Rotate polarizations to align with a_1, a_2
theta0 = dacos(dot_product(e_3,matmul(R_init,a_1)))
if(theta0>1d-6)then
 ninit = crossRR(e_3,matmul(R_init,a_1))
 ninit = ninit/vlen(ninit)
 R_pol = R_aa(ninit,-theta0)
else
 R_pol = eye(3)
!	print*, 'check that out 2'
end if
R_init = matmul(R_pol,R_init)

if(dot_product(matmul(R_init,a_3),k)<1d0-1d-6)then
 print*, ' WARNING (R_theta): Did not succeed to rotate a_3 onto k'
end if

! Rotation of a_3 from k to current theta away from k
R_thta = R_aa(e_3, theta)
R = matmul(R_thta,R_init)

if(dacos(dot_product(matmul(R,a_3), k))-theta>1d-6 ) then
 print*, 'WARNING (R_theta): Rotation of a_3 to given theta failed'
end if

end function R_theta

!******************************************************************************

end module common

