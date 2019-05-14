program main
! Copyright (c) 2018 Joonas Herranen and University of Helsinki
! All rights reserved.
! The MIT License is applied to this software, see LICENSE
   use common
   use io
   use integrator
   use postprocessing

   implicit none
   real(dp) :: tick, tock
   character(80) :: tstr

   call splash('v0.9', tick)

   call check_paramsfile()
   call read_params()
   call read_arguments()
   call setup_band()

   call init_geometry()

   call allocate_Ti()

   if (matrices%Tmat == 1) then
      if (use_mie /= 1) then
         call read_T()
         call fix_band()
      else
         call mie_T_matrix()
      end if
      write (*, '(A, 20F6.3)') ' Wavelengths in medium (um): ', 2d6*pi/mesh%ki
   else 
      if(use_mie /= 1) call T_empty()
      write (*, '(A, 20F6.3)') ' Wavelengths in medium (um): ', 2d6*pi/mesh%ki
      call calc_T()
   end if

   if (run_test == 0) then
      call integrate()
   else
      call tests()
   end if

   call compute_mueller()

   call curr_time(tock)
   tstr = time2str(tock-tick)
   write (*, '(2A)') ' Total execution time: ', trim(tstr)

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
         ! call torque_efficiency()
         call RAT_efficiency(90, 30, Npsi_in = 50)
         ! call RAT_alignment()
      end if 
      if (run_test == 3) call write_fields()
      if (run_test == 4) call force_map()

   end subroutine tests

end program main
