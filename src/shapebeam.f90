module shapebeam
use T_matrix

implicit none
contains

!*******************************************************************************
! The Gaussian amplitude profile in localized approximation (kw0>=5) a'la
! MSTM 3.0 (Mackowski et al 2013)
subroutine gaussian_beams(matrices,mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: width
integer :: i 

width = 5d0/(minval(mesh%ki))

do i = 1,matrices%bars
    call gaussian_beam_shape(matrices,mesh,i,matrices%Nmaxs(i),width)
end do

end subroutine gaussian_beams

!*******************************************************************************
! The Gaussian amplitude profile in localized approximation (kw0>=5) a'la
! MSTM 3.0 (Mackowski et al 2013)
subroutine gaussian_beam_shape(matrices, mesh, i, Nmax, width)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: C, gn, kw0, width, maxgn
integer :: n, m, ind, i, Nmax

kw0 = mesh%ki(i)*width

if(kw0<5d0) then
    write(*,'(A)'), "    Problem: width of Gaussian beam at focal point is smaller than wavelength!"
    write(*,'(2(A, ES9.3))'), "     Wavelength is ", 2d0*pi/mesh%ki(i), ", width is ", width
end if

ind = 0
maxgn = 0d0
do n = 1,Nmax
    gn = dexp(-((dble(n)+.5d0)/kw0)**2)
    do m = -n,n
        ind = ind + 1
        matrices%as(ind,i) = gn*matrices%as(ind,i)
        matrices%bs(ind,i) = gn*matrices%bs(ind,i)
    end do 
end do


end subroutine gaussian_beam_shape

!*******************************************************************************

subroutine laguerre_gauss_farfield(nmax, p, l, w0, PP, x, y, trunc, offset)
integer :: nmax, axisymmetry, total_modes, i
integer, dimension(:), allocatable :: nn, mm
real(dp) :: p, k, w0, PP, l, x, y, trunc, offset(3), zero_rejection_level, &
medium_refractive_index, beam_wavelength0, beam_power, speed_in_medium, &
epsilon0, kappa, xcomponent, ycomponent, radial_mode, azimuthal_mode, &
truncation_angle
    
zero_rejection_level = 1e-8

medium_refractive_index = 1
! beam_wavelength0 = 1064e-9
beam_wavelength0 = 1
beam_power = 1e-3
speed_in_medium = cc/medium_refractive_index
epsilon0 = 8.854187817e-12
kappa = medium_refractive_index**2

k = 2*pi * medium_refractive_index / beam_wavelength0
xcomponent = x
ycomponent = y
radial_mode = p
azimuthal_mode = l

if(trunc>1d-7)then
    truncation_angle = trunc
else
    truncation_angle = 90
end if

total_modes = nmax**2 + 2*nmax
allocate(nn(total_modes), mm(total_modes))
do i = 1,total_modes
    nn(i) = floor(sqrt(dble(i)));
    mm(i) = i - nn(i)**2 - nn(i)
end do

end subroutine laguerre_gauss_farfield

!*******************************************************************************

subroutine fields_out(matrices,mesh,which,n)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: n,nn, i, j, ind, which
real(dp) :: lim
real(dp), allocatable :: z(:), y(:), grid(:,:)
complex(dp) :: F(3), G(3) 
complex(dp), dimension(:,:), allocatable :: E  

nn = n*n

grid = field_grid(matrices,mesh)

allocate(E(3,nn))

matrices%field_points = grid
do i = 1,nn
    call calc_fields(matrices%as(:,which),matrices%bs(:,which),&
    dcmplx(mesh%ki(which)), grid(:,i),F,G,0)
    E(:,i) = F 
end do
matrices%E_field = E

end subroutine fields_out

!*******************************************************************************

subroutine scat_fields_out(matrices, mesh, which, n)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: n, nn, i, j, ind, which
real(dp) :: lim
real(dp), allocatable :: grid(:,:)
complex(dp), dimension(:) , allocatable :: p, q, p90, q90
complex(dp) :: F(3), G(3) 
complex(dp), dimension(:,:), allocatable :: E  

nn = n*n

grid = field_grid(matrices,mesh)

! Ensure that directions are ok. They might already be...
matrices%Rkt = eye(3)
allocate(E(3,nn))
call rot_setup(matrices)
call scattered_fields(matrices,1d0,p,q,p90,q90,which)

matrices%field_points = grid
do i = 1,nn
    call calc_fields(p,q,dcmplx(mesh%ki(which)), grid(:,i),F,G,1)
    E(:,i) = F 
end do

matrices%E_field = E


end subroutine scat_fields_out

!*******************************************************************************

function field_grid(matrices,mesh) result(grid)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: n, nn, i, j, ind, which
real(dp) :: lim
real(dp), allocatable :: z(:), y(:), grid(:,:)
complex(dp) :: F(3), G(3) 
complex(dp), dimension(:,:), allocatable :: E  

lim = 3.5d0
n = 100
nn = n*n
allocate(z(n),y(n),grid(3,nn))

call linspace(-lim,lim,n,z)
call linspace(-lim,lim,n,y)
grid = 0d0
do i = 1,n
    do j = 1,n
        ind = n*(j-1)+i 
        grid(3,ind) = z(i)*mesh%a 
        grid(2,ind) = y(j)*mesh%a
    end do
end do

end function field_grid

!*******************************************************************************

subroutine write_fields(matrices,mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: width
integer :: i, n
character(LEN=80) :: gridname, fieldname, scatfield
i = 1
n = 100
gridname = 'grid_xyz.h5'
fieldname = 'E_field.h5'
scatfield = 'E_scat.h5'

call fields_out(matrices,mesh,i,n)
call write2file(dcmplx(matrices%field_points),gridname)
call write2file(matrices%E_field,fieldname)

call scat_fields_out(matrices,mesh,i,n)
call write2file(matrices%E_field,scatfield)

end subroutine write_fields

end module shapebeam