module shapebeam
use T_matrix

implicit none
contains

! The Gaussian amplitude profile in localized approximation (kw0<0.2) a'la
! MSTM 3.0 (Mackowski et al 2013)
subroutine gaussian_beams(matrices,mesh,width)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: width
integer :: i 

do i = 1,matrices%bars
    call gaussian_beam_shape(matrices,mesh,i,matrices%Nmaxs(i),width)
end do

end subroutine gaussian_beams

! The Gaussian amplitude profile in localized approximation (kw0<0.2) a'la
! MSTM 3.0 (Mackowski et al 2013)
subroutine gaussian_beam_shape(matrices, mesh, i, Nmax, width)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: C, gn, kw0, width
integer :: n, m, ind, i, Nmax

kw0 = mesh%ki(i)*width
if(kw0>0.2d0) then
    print*, "Problem: width of Gaussian beam at focal point too large for the wavelength!"
    print*, "Wavelength is ", 2d0*pi/mesh%ki(i)
end if

ind = 0
do n = 1,Nmax
    gn = dexp(-((dble(n)+.5d0)/kw0)**2)
    do m = -n,n
        ind = ind + 1
        matrices%as(ind,i) = gn*matrices%as(ind,i)
        matrices%bs(ind,i) = gn*matrices%bs(ind,i)
    end do 
end do

end subroutine gaussian_beam_shape


end module shapebeam