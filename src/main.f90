program main
use common
use io
use integrator
use postprocessing

implicit none
type (data) :: matrices
type (mesh_struct) :: mesh

call splash()

call check_paramsfile(matrices,mesh)
call read_params(matrices, mesh)
call read_arguments(matrices, mesh)
call setup_band(matrices, mesh)

call init_geometry(matrices,mesh)

call allocate_Ti(matrices)

if (matrices%Tmat == 1 .AND. file_exists(matrices%tname)) then
   call read_T(matrices, mesh)
   call fix_band(matrices, mesh)
else
   write(*,'(A,1ES8.2)') '  a = ', mesh%a
   if(matrices%singleT == 1) then
      if(.NOT. file_exists(matrices%tname)) call T_empty(matrices, mesh)
   end if

   call calc_T( matrices, mesh )
   call write_T( matrices, mesh )
end if

if(run_test == 0) then
   call integrate(matrices,mesh) 
else
   call tests(matrices,mesh)
end if

if(trim(matrices%mueller_mode) /= 'none') call compute_mueller(matrices,mesh)

contains

!******************************************************************************

subroutine tests(matrices,mesh)
type(data) :: matrices
type(mesh_struct) :: mesh

call polarization(matrices)
call allocate_inc_wave(matrices, mesh)

if (use_mie==1)then
   call mie_params(matrices,mesh)
else
   if(mesh%is_mesh == 1) then
      call vie_params(matrices,mesh)
   else
      call aggr_params(matrices,mesh)
   end if
end if

matrices%x_CM = mesh%CM
call diagonalize_inertia(matrices, mesh)
call init_values(matrices, mesh)

if(beam_shape == 1) call gaussian_beams(matrices,mesh)
if(beam_shape == 2) call laguerre_gaussian_beams(matrices,mesh,p,l)
if(beam_shape == 3) call bessel_beams(matrices,mesh)

if(run_test == 1) call test_methods(matrices, mesh)
if(run_test == 2) call torque_efficiency(matrices, mesh)
if(run_test == 3) call stability_analysis(matrices, mesh)
if(run_test == 4) call write_fields(matrices,mesh)
if(run_test == 5) call test_mueller(matrices,mesh)

end subroutine tests

end program main
