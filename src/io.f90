module io
   use common
   use h5io

   implicit none

contains
! MY ROUTINES
! CLEARTEXT READ
! CLEARTEXT WRITE

! MY ROUTINES *****************************************************************
!****************************************************************************80

   subroutine splash(v)
      character(4) :: v ! Version
      print *, '******************************************************************************'
      print *, '**                                                                          **'
      print *, '**                         JVIE T-Matrix Dynamics ',v,'                      **'
      print *, '**                                                                          **'
      print *, '******************************************************************************'
      print *, ''
      call curr_time()

   end subroutine splash

!****************************************************************************80

   subroutine curr_time()
      character(8)  :: date, fmt, fmt2
      character(10) :: time
      character(5)  :: zone
      character(3)  :: hh, mm, ss, ms
      integer, dimension(8) :: values

      fmt = '(I2.2)'
      fmt2 = '(I3.3)'
      call date_and_time(date, time, zone, values)
      write (hh, fmt) values(5)
      write (mm, fmt) values(6)
      write (ss, fmt) values(7)
      write (ms, fmt2) values(8)
      print'(3(A, I0),8A)', ' Time: ', values(3), '.', values(2), '.', values(1), ' -- ', &
         trim(hh), ':', trim(mm), ':', trim(ss), '.', trim(ms)

   end subroutine curr_time

!****************************************************************************80

   subroutine print_bar(i1, Nmax)
      integer :: i1, Nmax, k
      real(dp) :: r
      character(len=1) :: bar, back

      back = char(8)
      bar = '='
      r = 100d0/Nmax
! print the percentage and the bar without line change, then reset cursor
      if (floor(i1*r) - floor((i1 - 1)*r) > 0) then
         write (6, '(2x,1i3,1a1,2x,1a1,256a1)', advance='no') &
            ceiling(r*i1), '%', '|', (bar, k=1, 50*i1/Nmax)
         write (6, '(256a1)', advance='no') (back, k=1, (50*i1/Nmax) + 9)
         if (i1 == Nmax) then
            write (*, *) ''
         end if
      end if

   end subroutine print_bar

!****************************************************************************80

   subroutine print_mat(mat, matname)
      real(dp), dimension(:, :) :: mat
      character(len=*) :: matname
      integer :: i, m

      write (*, '(2A)') trim(matname), ' = '
      m = size(mat, 1)
      do i = 1, m
         write (*, '(10F8.4)') mat(i, :)
      end do

   end subroutine print_mat

! CLEARTEXT READ **************************************************************
!****************************************************************************80
   subroutine check_paramsfile()
! Subroutine to read input parameters
      integer ::  i
      character(len=80) :: arg

      do i = 1, command_argument_count(), 2
         call get_command_argument(i, arg)

         select case (arg)
         case ('-p', '--paramsfile')
            call get_command_argument(i + 1, matrices%paramsfile)
            write (*, '(2A)') ' Paramsfile: ', trim(matrices%paramsfile)
         end select
      end do
   end subroutine check_paramsfile

!****************************************************************************80
! Subroutine to read input parameters
   subroutine read_arguments()
      integer :: i
      character(len=80) :: arg_name, arg

      matrices%singleT = 0

      do i = 1, command_argument_count(), 2
         call get_command_argument(i, arg_name)

         select case (arg_name)

         case ('-d', '--debug')
            call get_command_argument(i, arg)
            debug = 1
         case ('-m', '--mesh')
            call get_command_argument(i + 1, arg)
            mesh%meshname = arg
            write (*, '(2A)') ' Mesh: ', trim(mesh%meshname)
         case ('-T', '--Tmat')
            call get_command_argument(i + 1, arg)
            matrices%tname = arg
            write (*, '(2A)') ' T-matrix: ', trim(matrices%tname)
         case ('-o', '--out')
            call get_command_argument(i + 1, arg)
            matrices%out = arg
            write (*, '(2A)') ' Log: ', 'out/log'//trim(matrices%out)
         case ('--mueller_mode')
            call get_command_argument(i + 1, arg)
            matrices%mueller_mode = trim(arg)
            write (*, '(2A)') ' Mueller mode: ', trim(matrices%mueller_mode)
         case ('-p', '--paramsfile') ! Skip, as already read in check_paramsfile
         case ('--refr')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%refr
         case ('--refi')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%refi
         case ('-w', '--wavelen')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%whichbar
         case ('-S', '--singleT')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%singleT
            matrices%Tmat = 0
         case ('-s', '--seed')
            call get_command_argument(i + 1, arg)
            read (arg, *) seedling
            matrices%R = rand_rot()
         case ('--Mie')
            call get_command_argument(i + 1, arg)
            read (arg, *) use_mie
         case ('-B')
            call get_command_argument(i + 1, arg)
            read (arg, *) calc_extra_torques
         case ('-x', '--xi')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%xi_in
         case ('-r', '--relax')
            call get_command_argument(i + 1, arg)
            read (arg, *) relax

         case ('-h', '--help')
            write (*, '(A)') ' Commands        Value       Description'
            write (*, '(A)') '---------------------------------------'
            write (*, '(A)') ' -d --debug                  Print more detailed info'
            write (*, '(A)') ' -m --mesh       mesh.h5     Mesh geometry'
            write (*, '(A)') ' -T --Tmat       T.h5        T-matrix file'
            write (*, '(A)') ' -o --out                    Output file identifier'
            write (*, '(A)') ' -p --paramsfile params.in   Read input parameters from file'
            write (*, '(A)') '    --refr       0.0         Real part of refractive index'
            write (*, '(A)') '    --refi       0.0         Imaginary part of refractive index'
            write (*, '(A)') ' -w --wavelen    0           Choose wavelength from the T-matrix'
            write (*, '(A)') ' -S --singleT    0           Calculate only one T-matrix, of wb'
            write (*, '(A)') ' -s --seed       0           RNG seed'
            write (*, '(A)') '    --Mie        1           Use the Mie sphere'
            write (*, '(A)') ' -B              1           Use external magnetic field'
            write (*, '(A)') ' -x --xi         0.0         Precession angle about B'
            write (*, '(A)') ' -r --relax      1           Use internal relaxation'


            stop
         case default
            print'(a,a,/)', 'Unrecognized command-line option: ', arg_name
            stop
         end select
      end do

      if (matrices%singleT == 1 .AND. matrices%whichbar == 0) then
         write (*, '(A)') ' Warning: Single T-matrix chosen but all bars used. Now whichbar = 1'
         matrices%whichbar = 1
      end if

   end subroutine read_arguments

!****************************************************************************80
! Subroutine to read a input file for nifty usage
   subroutine read_params()
! Input related variables
      character(len=150) :: buffer, label
      real(dp) :: temp, tempii
      integer :: pos
      integer, parameter :: fh = 15
      integer :: ios = 0
      integer :: line = 0

! Control file variables
      real(dp) :: khat(3), lambda1, lambda2
      integer :: whichbar
      character(len=4) :: R0 = 'i'
      character(len=8) :: mode

      open (fh, file=matrices%paramsfile)

! ios is negative if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive if an error was
! detected.  ios is zero otherwise.

      do while (ios == 0)
         read (fh, '(A)', iostat=ios) buffer
         if (ios == 0) then
            line = line + 1

            ! Find the first instance of whitespace.  Split label and data.
            pos = scan(buffer, '    ')
            label = buffer(1:pos)
            buffer = buffer(pos + 1:)

            select case (label)
! Parameters
            case ('beam_shape')
               read (buffer, *, iostat=ios) beam_shape
            case ('pl')
               read (buffer, *, iostat=ios) p, l
            case ('int_mode')
               read (buffer, *, iostat=ios) int_mode
            case ('shortlog')
               read (buffer, *, iostat=ios) shortlog
            case ('test_forces')
               read (buffer, *, iostat=ios) run_test
            case ('it_max')
               read (buffer, *, iostat=ios) it_max
            case ('window')
               read (buffer, *, iostat=ios) window
            case ('it_log')
               read (buffer, *, iostat=ios) it_log
! Matrices
            case ('E')
               read (buffer, *, iostat=ios) matrices%E
            case ('w0')
               read (buffer, *, iostat=ios) matrices%w
            case ('pos')
               read (buffer, *, iostat=ios) matrices%x_CM
            case ('R0')
               read (buffer, *, iostat=ios) R0
            case ('dt')
               read (buffer, *, iostat=ios) matrices%dt0
            case ('refr')
               read (buffer, *, iostat=ios) temp
               if (temp > 1d-7) matrices%refr = temp
            case ('refi')
               read (buffer, *, iostat=ios) tempii
            case ('tol_m')
               read (buffer, *, iostat=ios) matrices%tol_m
            case ('rot_max')
               read (buffer, *, iostat=ios) matrices%rot_max
            case ('khat')
               read (buffer, *, iostat=ios) khat
               matrices%khat = khat/vlen(khat)
            case ('polarization')
               read (buffer, *, iostat=ios) matrices%polarization
            case ('bars')
               read (buffer, *, iostat=ios) matrices%bars
            case ('lambda1')
               read (buffer, *, iostat=ios) lambda1
               matrices%lambda1 = lambda1*1.d-9
            case ('lambda2')
               read (buffer, *, iostat=ios) lambda2
               matrices%lambda2 = lambda2*1.d-9
            case ('T')
               read (buffer, *, iostat=ios) matrices%temp
            case ('Tmat')
               read (buffer, *, iostat=ios) matrices%Tmat
            case ('whichbar')
               read (buffer, *, iostat=ios) whichbar
               if (matrices%whichbar == 0) matrices%whichbar = whichbar
            case ('mueller_mode')
               read (buffer, *, iostat=ios) mode
               if (trim(mode) /= 'none') matrices%mueller_mode = mode
            case ('waves')
               read (buffer, *, iostat=ios) matrices%waves
            case ('B')
               read (buffer, *, iostat=ios) matrices%B_len
               if (matrices%B_len > 1d-14) calc_extra_torques = 1
            case ('B_psi')
               read (buffer, *, iostat=ios) matrices%B_psi
               matrices%B_psi = matrices%B_psi*pi/180
            case ('Td')
               read (buffer, *, iostat=ios) matrices%Td
            case ('Tgas')
               read (buffer, *, iostat=ios) matrices%Tgas
            case ('nH')
               read (buffer, *, iostat=ios) matrices%nH
            case ('Kw')
               read (buffer, *, iostat=ios) matrices%Kw
! Mesh
            case ('rho')
               read (buffer, *, iostat=ios) mesh%rho
            case ('a')
               read (buffer, *, iostat=ios) mesh%a
            case ('tol')
               read (buffer, *, iostat=ios) mesh%tol
            case ('maxit')
               read (buffer, *, iostat=ios) mesh%maxit
            case ('restart')
               read (buffer, *, iostat=ios) mesh%restart
            case ('order')
               read (buffer, *, iostat=ios) mesh%M_ex
            case ('cell_size')
               read (buffer, *, iostat=ios) mesh%grid_size
            case ('near_zone')
               read (buffer, *, iostat=ios) mesh%near_zone
            case ('expansion_order')
               read (buffer, *, iostat=ios) mesh%order
            case ('drag')
               read (buffer, *, iostat=ios) mesh%drag
               mesh%drag = mesh%drag*1d6
            case default
               !print *, 'Skipping invalid label at line', line

            end select
         end if
      end do

      close (fh)

      if (R0 == 'r') then
         matrices%R = rand_rot()
      else if (R0 == 't') then
         matrices%R = reshape(dble([0, 1, 0, 0, 0, 1, 1, 0, 0]), [3, 3])
      else
         matrices%R = eye(3)
      end if

      matrices%x_CM = matrices%x_CM*matrices%lambda2

      if (temp > 1d-7) matrices%refi = tempii
      it_stop = it_max
      if(it_log > it_max) it_log = 0
      matrices%B = matrices%B_len*(matmul(R_aa([0d0, 1d0, 0d0], matrices%B_psi), [0d0, 0d0, 1d0]))

      allocate (mesh%ki(matrices%bars))
      allocate (matrices%E_rel(matrices%bars))
      allocate (matrices%Nmaxs(matrices%bars))

   end subroutine read_params

!****************************************************************************80

   function get_last_line_no(fle) result(lineno)
      character(len=80) :: fle
      integer :: lineno
      integer :: ios = 0
      integer :: nbline = 0

      open (11, file=trim(fle))
      do while (.true.)
         read (11, *, iostat=ios) ! ios should have been declared as an integer
         if (ios > 0) then
            stop'problem somewhere (Called from io:get_last_line_no)'
         else if (ios < 0) then ! end of file is reached
            exit
         else
            nbline = nbline + 1
         end if
      end do
      close (11)
      lineno = nbline

   end function get_last_line_no

!****************************************************************************80

   subroutine read_log(no)
      integer, parameter :: fh = 15
      integer :: line

! Control file variables
      real(dp) :: wlast(3), t_tresh, t1, t2, tf, w_av(3), Rab(3, 3)
      integer :: firstlineno, lastlineno, no, numlines, i, n1, n2
      real(dp) :: x(3), v(3), w(3), J(3), F(3), N(3), t, R(9)

      numlines = no
      lastlineno = get_last_line_no(matrices%out)
      if(numlines>=lastlineno-23) numlines = lastlineno-23
      firstlineno = lastlineno-numlines
      open (fh, file='out/log'//trim(matrices%out))
      do i = 1, lastlineno-1
         read (fh, *)
      end do

      read (fh, *) line, x, v, w, J, N, F, t, R

      wlast = w/dsqrt(w(1)**2 + w(2)**2 + w(3)**2)

      close (fh)

      t_tresh = 4d0*pi/dsqrt(w(1)**2 + w(2)**2 + w(3)**2)
      tf = t
      t1 = t - t_tresh
      t2 = t - 10*t_tresh

      w_av = 0d0
      open (fh, file=trim(matrices%out))
      do i = 1, lastlineno
         if (i >= firstlineno) then
            read (fh, *) line, x, v, w, J, N, F, t, R
            w_av = w_av + w
            if (t > t2) then
               ! print*, line
               t2 = tf
            end if
            if (t > t1) then
               ! print*, line
               t1 = tf
            end if
            if (i == firstlineno) then
               ! print*, line
               n1 = line
            end if
         else
            read (fh, *)
         end if
      end do
      n2 = line
      w_av = w_av/(lastlineno - firstlineno + 1)
      ! print*, 'w_av = ', w_av
      close (fh)
      Rab = rotate_a_to_b(w_av, [1d0, 0d0, w_av(3)])
      ! call print_mat(Rab,'R_ab')
      matrices%R_al = Rab
      allocate (matrices%RRR(3, 3, n2 - n1 + 1))
      allocate (matrices%www(3, n2 - n1 + 1))

      open (fh, file=trim(matrices%out))
      do i = 1, lastlineno
         if (i > firstlineno) then
            read (fh, *) line, x, v, w, J, N, F, t, R
            matrices%RRR(:, :, i - firstlineno) = reshape(R, [3, 3])
            matrices%www(:, i - firstlineno) = w
         else
            read (fh, *)
         end if
      end do
      close (fh)

   end subroutine read_log

!****************************************************************************80

   subroutine read_mueller(A, fname)
      real(dp), allocatable, intent(out) :: A(:, :)
      character(len=80), intent(in) :: fname
      integer, dimension(2) :: dims ! Dataset dimensions

      integer     :: i

      dims = (/181, 18/)
      allocate (A(dims(1) - 1, dims(2)))

      open (unit=1, file=trim(fname))
      read (1, *)
      do i = 2, dims(1)
         read (1, '(2F8.4,16ES14.4)') A(i, :)
      end do

      close (1)

   end subroutine read_mueller

! CLEARTEXT WRITE *************************************************************
!****************************************************************************80

   subroutine write_array(A, fname)
      real(dp), intent(in) :: A(:, :)
      character(len=8), intent(in) :: fname
      integer(HSIZE_T), dimension(2) :: dims ! Dataset dimensions
      integer     :: i

      dims = [int(size(A, 1), 8), int(size(A, 2), 8)]

      open (unit=1, file=trim(fname), ACTION="write", STATUS="replace")

      do i = 1, int(dims(1))
         write (1, '(1024ES14.6)') A(i, :)
      end do

      close (1)

   end subroutine write_array

!****************************************************************************80

   subroutine write_mueller(A)
      real(dp), intent(in) :: A(:, :)
      character(len=120) :: fname
      integer(HSIZE_T), dimension(2) :: dims ! Dataset dimensions
      integer     :: i

      fname = 'out/mueller-'//trim(matrices%mueller_mode)//trim(matrices%out)
      dims = [int(size(A, 1), 8), int(size(A, 2), 8)]

      open (unit=1, file=trim(fname), ACTION="write", STATUS="replace")
      write (1, '(18A)') '   phi ', '    theta ', '      S11  ', '         S12  ', &
         '         S13  ', '         S14  ', '         S21  ', '         S22  ', '         S23  ', &
         '         S24  ', '         S31  ', '         S32  ', '         S33  ', '         S34  ', &
         '         S41  ', '         S42  ', '         S43  ', '         S44  '

      do i = 1, int(dims(1))
         write (1, '(2F8.4,16ES14.6)') A(i, :)
      end do

      close (1)

   end subroutine write_mueller

!****************************************************************************80

   subroutine init_values()
      matrices%k_orig = matrices%khat
      matrices%E0_orig = real(matrices%E0hat)
      matrices%E90_orig = real(matrices%E90hat)

      matrices%q_mean = 0d0

      matrices%xn = matrices%x_CM ! Input in lab frame (if needed)
      matrices%vn = matrices%v_CM ! Input in lab frame
      matrices%wn = matrices%w
      matrices%dt = matrices%dt0
      matrices%tt = 0d0
      matrices%q = mat2quat(matrices%R)
      matrices%N = dble([0d0, 0d0, 0d0])
      matrices%F = dble([0d0, 0d0, 0d0])
      matrices%J = matmul(matmul(matrices%P, matrices%R), matmul(matrices%I, matrices%w)) ! In lab frame!
      matrices%R90_init = reshape(dble([0d0, 1d0, 0d0, -1d0, 0d0, 0d0, 0d0, 0d0, 1d0]), [3, 3])
   end subroutine init_values

!****************************************************************************80

   subroutine update_values()

      matrices%x_CM = matrices%xn ! In lab frame!
      matrices%v_CM = matrices%vn ! In lab frame!
      matrices%w = matrices%wn ! In body frame!
      matrices%R = matrices%Rn
      matrices%q = matrices%qn
      matrices%tt = matrices%tt + matrices%dt
      matrices%J = matmul(matmul(matrices%P, matrices%R), &
                          matmul(matrices%I, matrices%w))
   end subroutine update_values

!****************************************************************************80

   subroutine start_log()
      integer :: i
      character(len=120) :: fname
      fname = 'out/log' // trim(matrices%out)

      matrices%buffer = 0
      matrices%q_mean = 0d0
      matrices%q_var = 0d0
      allocate(matrices%q_list(window))
      matrices%q_list = 0d0

      open (unit=1, file=fname, ACTION="write", STATUS="replace")

      write (1, '(A, A)') 'meshname =   ', trim(mesh%meshname)
      write (1, '(A, 3f7.3)') 'k_hat    = ', matrices%khat
      write (1, '(A, 3ES11.3)') 'dt       = ', matrices%dt
      write (1, '(A,I0)') 'Nmax     =   ', it_max
      write (1, '(A, 3ES11.3)') 'wB0      = ', matrices%w
      write (1, '(A, 9f7.3)') 'R0       = ', matrices%R
      write (1, '(A, 3ES11.3)') 'rho      = ', mesh%rho
      write (1, '(A, 3ES11.3)') 'CM       = ', matrices%x_CM
      write (1, '(A, 3ES11.3)') 'a        = ', mesh%a
      write (1, '(A, 80ES11.3)') 'ki       = ', mesh%ki
      write (1, '(A, 2f8.2)') 'range nm = ', matrices%lambda1*1d9, matrices%lambda2*1d9
      write (1, '(A, 3ES11.3)') 'mass     = ', mesh%mass
      write (1, '(A, 3ES11.3)') 'V        = ', mesh%V
      write (1, '(A, 3ES11.3)') 'rot_max  = ', matrices%rot_max
      write (1, '(A)') 'Ip       = '
      write (1, '(9ES11.3)') matrices%Ip(:)
      write (1, '(A)') 'Q        = '
      do i = 1, 3
         write (1, '(3f7.3)') matrices%P(i, :)
      end do
      if(shortlog==0) write (1, '(A)') ' Log: n | t | x(1:3) | w(1:3) | v(1:3) | N(1:3) | F(1:3) | R(1:9) '
      write (1, '(A)') ' '

      close (1)

   end subroutine start_log

!****************************************************************************80

   subroutine append_log(n)
      integer :: n, i, md, ind
      character(len=120) :: fname, fmt

      fname = 'out/log'//trim(matrices%out)
      fmt = '(I0, A, ES16.8, A, 5(3ES16.8, A), 9ES16.8)'
      md = mod(n, 1000)

! If the modulo if zero, we are at the final place of the buffer. Otherwise, 
! just save to buffer position given by the modulo.
      ind = md
      if (md == 0) ind = 1000
      matrices%x_buf(:, ind) = matrices%x_CM
      matrices%w_buf(:, ind) = matrices%w
      matrices%v_buf(:, ind) = matrices%v_CM
      matrices%N_buf(:, ind) = matrices%N
      matrices%F_buf(:, ind) = matrices%F
      matrices%t_buf(:, ind) = matrices%tt
      matrices%R_buf(:, :, ind) = matrices%R

      if (md == 0 .OR. n == it_stop) then
         open (unit=1, file=fname, action="write", position="append", STATUS="old")
         if (n > it_stop - it_log .OR. it_log == 0) then
            do i = 1, 1000
               if (i + matrices%buffer <= it_stop) then
                  write (1, fmt) i + matrices%buffer,  ' |',&
                     matrices%t_buf(:, i),  ' |',&
                     matrices%x_buf(:, i),  ' |',&
                     matrices%w_buf(:, i),  ' |',&
                     matrices%v_buf(:, i),  ' |',&
                     matrices%N_buf(:, i),  ' |',&
                     matrices%F_buf(:, i),  ' |',&
                     matrices%R_buf(:, :, i)
               end if
            end do
         end if
         close (1)
         matrices%buffer = n
      end if

   end subroutine append_log

!****************************************************************************80

   subroutine alignment_log()
      integer :: n, i, md, ind
      character(len=120) :: fname, fmt

      fname = 'out/log'//trim(matrices%out)
      fmt = '(A, 3ES16.8)'
      open (unit=1, file=fname, action="write", position="append", STATUS="old")

      write (1, fmt) 'J ', matrices%J
      write (1, fmt) 'q0 ', matrices%q_param0
      write (1, fmt) 'tau_rad ', matrices%tau_rad
      write (1, fmt) 'tau_int ', matrices%tau_int
      write (1, fmt) 'w_thermal ', matrices%w_thermal

      close (1)

   end subroutine alignment_log

end module
