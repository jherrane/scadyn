program main
   use common
   use io
   use integrator
   use postprocessing

   implicit none
   type(data) :: matrices
   type(mesh_struct) :: mesh

   call splash('v0.5')

   call check_paramsfile(matrices)
   call read_params(matrices, mesh)
   call read_arguments(matrices, mesh)
   call setup_band(matrices, mesh)

   call init_geometry(matrices, mesh)

   call allocate_Ti(matrices, mesh)

   if (matrices%Tmat == 1 .AND. file_exists(matrices%tname)) then
      call read_T(matrices, mesh)
      call fix_band(matrices, mesh)
      write (*, '(A, 20F6.3)') ' Wavelengths in um: ', 2d6*pi/mesh%ki
   else
      if (matrices%singleT == 1) then
         call T_empty(matrices, mesh)
      end if
      write (*, '(A, 20F6.3)') ' Wavelengths in um: ', 2d6*pi/mesh%ki
      call calc_T(matrices, mesh)
      call write_T(matrices, mesh)
   end if

   if (run_test == 0) then
      call integrate(matrices, mesh)
   else
      call tests(matrices, mesh)
   end if

   call compute_mueller(matrices, mesh)

contains

!****************************************************************************80

   subroutine tests(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      real(dp) :: q
      call polarization(matrices, mesh)
      call allocate_inc_wave(matrices, mesh)

      if (use_mie == 1) then
         call mie_params(matrices, mesh)
      else if (mesh%is_mesh == 1) then
         call vie_params(matrices, mesh)
      end if

      matrices%x_CM = mesh%CM
      call diagonalize_inertia(matrices, mesh)
      call interstellar_env(matrices, mesh)
      call init_values(matrices, mesh)

      if (beam_shape == 1) call gaussian_beams(matrices, mesh)
      if (beam_shape == 2) call laguerre_gaussian_beams(matrices, mesh, p, l)
      if (beam_shape == 3) call bessel_beams(matrices, mesh)

      if (run_test == 1) call test_methods(matrices, mesh)
      if (run_test == 2) then
         ! call torque_efficiency(matrices, mesh, q)
         ! call RAT_efficiency(matrices, mesh, 60, 20, Npsi_in = 4)
         call RAT_alignment(matrices, mesh)
      end if 
      if (run_test == 3) call write_fields(matrices, mesh)
      if (run_test == 4) call test_mueller(matrices, mesh, 9, 18)

   end subroutine tests

end program main
