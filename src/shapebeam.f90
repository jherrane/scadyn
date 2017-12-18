module shapebeam
use T_matrix

implicit none
contains

! The Gaussian amplitude profile in localized approximation (kw0<0.2) a'la
! MSTM 3.0 (Mackowski et al 2013)
subroutine gaussian_beams(matrices,mesh)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: width
integer :: i 

width = 1d0/(5.1d0*maxval(mesh%ki))

do i = 1,matrices%bars
    call gaussian_beam_shape(matrices,mesh,i,matrices%Nmaxs(i),width)
end do

end subroutine gaussian_beams

! The Gaussian amplitude profile in localized approximation (kw0<0.2) a'la
! MSTM 3.0 (Mackowski et al 2013)
subroutine gaussian_beam_shape(matrices, mesh, i, Nmax, width)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: C, gn, kw0, width, maxgn
integer :: n, m, ind, i, Nmax

kw0 = mesh%ki(i)*width
if(kw0>5d0) then
    print*, "Problem: width of Gaussian beam at focal point too large for the wavelength!"
    print*, "Wavelength is ", 2d0*pi/mesh%ki(i)
end if

ind = 0
maxgn = 0d0
do n = 1,Nmax
    gn = dexp(-((dble(n)+.5d0)/kw0)**2)
    if(gn>maxgn) maxgn = gn

    do m = -n,n
        ind = ind + 1
        matrices%as(ind,i) = gn*matrices%as(ind,i)
        matrices%bs(ind,i) = gn*matrices%bs(ind,i)
    end do 
end do

matrices%as(:,i) = matrices%as(:,i)/maxgn
matrices%bs(:,i) = matrices%bs(:,i)/maxgn

end subroutine gaussian_beam_shape


subroutine fields_out(matrices,mesh,which)
type (mesh_struct) :: mesh
type (data) :: matrices
integer :: n, nn, i, j, ind, which
real(dp) :: lim
real(dp), allocatable :: z(:), y(:), grid(:,:)
complex(dp) :: F(3), G(3) 
complex(dp), dimension(:,:), allocatable :: E  

lim = 3d0
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

allocate(E(3,nn))

matrices%field_points = grid
do i = 1,nn
    call calc_fields(matrices%as(:,which),matrices%bs(:,which),&
    dcmplx(mesh%ki(which)), grid(:,i),F,G,0)
    E(:,i) = F 
end do
matrices%E_field = E

end subroutine fields_out

end module shapebeam