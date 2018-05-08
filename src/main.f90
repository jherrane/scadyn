program main
   use common
   use io
   use integrator
   use postprocessing

   implicit none

   call splash()

   call check_paramsfile()
   call read_params()
   call read_arguments()
   call setup_band()

   call init_geometry()

   call allocate_Ti()

   if (matrices%Tmat == 1 .AND. file_exists(matrices%tname)) then
      call read_T()
      call fix_band()
   else
      write (*, '(A,1ES8.2)') '  a = ', mesh%a
      if (matrices%singleT == 1) then
         if (.NOT. file_exists(matrices%tname)) call T_empty()
      end if

      call calc_T()
      call write_T()
   end if

   if (run_test == 0) then
      call integrate()
   else
      call tests()
   end if

   if (trim(matrices%mueller_mode) /= 'none') call compute_mueller()

contains

!****************************************************************************80

   subroutine tests()
      call polarization()
      call allocate_inc_wave()

      if (mesh%is_mesh == 1) then
         call vie_params()
      else if (use_mie == 1) then
         call mie_params()
      else
         call aggr_params()
      end if

      matrices%x_CM = mesh%CM
      call diagonalize_inertia()
      call init_values()

      if (beam_shape == 1) call gaussian_beams()
      if (beam_shape == 2) call laguerre_gaussian_beams(p, l)
      if (beam_shape == 3) call bessel_beams()

      if (run_test == 1) call test_methods()
      if (run_test == 2) call torque_efficiency()
      if (run_test == 3) call stability_analysis()
      if (run_test == 4) call write_fields()
      if (run_test == 5) call test_mueller()

   end subroutine tests

end program main
