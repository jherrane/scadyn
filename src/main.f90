program main
   use common
   use io
   use integrator
   use postprocessing

   implicit none

   call splash('v0.5')

   call check_paramsfile()
   call read_params()
   call read_arguments()
   call setup_band()

   call init_geometry()

   call allocate_Ti()

   if (matrices%Tmat == 1) then
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
   else
      call tests()
   end if

   call compute_mueller()

contains

!****************************************************************************80

   subroutine tests()
      real(dp) :: q
      if (use_mie == 1) then
         call mie_params()
      else
         call vie_params()
      end if

      call diagonalize_inertia()
      call interstellar_env()

      call polarization()
      call init_values()

      call allocate_inc_wave()

      if (beam_shape == 1) call laguerre_gaussian_beams(p, l)
      if (beam_shape == 2) call bessel_beams()
      
      if (run_test == 1) call test_methods()
      if (run_test == 2) then
         call torque_efficiency(q)
         call RAT_efficiency(60, 20, Npsi_in = 4)
         ! call RAT_alignment()
      end if 
      if (run_test == 3) call write_fields()

   end subroutine tests

end program main
