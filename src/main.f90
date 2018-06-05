program main
   use common
   use io
   use integrator
   use postprocessing

   implicit none

   call splash('v0.4')

   call check_paramsfile()
   call read_params()
   call read_arguments()
   call setup_band()

   call init_geometry()

   call allocate_Ti()

   if (matrices%Tmat == 1 .AND. file_exists(matrices%tname)) then
      call read_T()
      call fix_band()
      write (*, '(A, 20F6.3)') ' Wavelengths in um: ', 2d6*pi/mesh%ki
   else
      if (matrices%singleT == 1) then
         call T_empty()
      end if
      write (*, '(A, 20F6.3)') ' Wavelengths in um: ', 2d6*pi/mesh%ki
      call calc_T()
      call write_T()
   end if

   if (run_test == 0) then
      call integrate()
      if (int_mode >= 1 .AND. int_mode < 3) then
         call compute_log_RAT()
      else if(int_mode ==3 ) then
         call stable_particle_RAT()
      end if
   else
      call tests()
   end if

   call compute_mueller()

contains

!****************************************************************************80

   subroutine tests()
      call polarization()
      call allocate_inc_wave()

      if (use_mie == 1) then
         call mie_params()
      else if (mesh%is_mesh == 1) then
         call vie_params()
      end if

      matrices%x_CM = mesh%CM
      call diagonalize_inertia()
      call init_values()

      if (beam_shape == 1) call gaussian_beams()
      if (beam_shape == 2) call laguerre_gaussian_beams(p, l)
      if (beam_shape == 3) call bessel_beams()

      if (run_test == 1) call test_methods()
      if (run_test == 2) then
         call torque_efficiency()
         call RAT_efficiency(60, 20, Npsi_in = 4)
      end if 
      if (run_test == 3) call stability_analysis()
      if (run_test == 4) call write_fields()
      if (run_test == 5) call test_mueller(9, 18)

   end subroutine tests

end program main
