module io
   use common
   use hdf5

   implicit none

contains
! MY ROUTINES
! CLEARTEXT READ
! CLEARTEXT WRITE
! HDF5 READ
! HDF5 WRITE

! MY ROUTINES *****************************************************************
!******************************************************************************

   subroutine splash()
      print *, '*************************************************************'
      print *, '**                                                         **'
      print *, '**             JVIE T-Matrix Dynamics v. 0.1               **'
      print *, '**                                                         **'
      print *, '*************************************************************'
      call curr_time()

   end subroutine splash

!******************************************************************************

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

!******************************************************************************

   subroutine print_bar(i1, Nmax)
      integer :: i1, Nmax, k
      character(len=1) :: bar, back
      real(dp) :: r

      back = char(8)
      bar = '='
      r = 100d0*dble(i1)/dble(Nmax)
! print the percentage and the bar without line change, then reset cursor
      if (r - dble(floor(r)) < 1d-7 .OR. i1 == 1) then
         write (6, '(2x,1i3,1a1,2x,1a1,256a1)', advance='no') &
            100*i1/Nmax, '%', '|', (bar, k=1, 50*i1/Nmax)
         write (6, '(256a1)', advance='no') (back, k=1, (50*i1/Nmax) + 9)
      end if

   end subroutine print_bar

!******************************************************************************

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

!******************************************************************************

   subroutine print_cvec(mat, matname)
      complex(dp), dimension(:) :: mat
      character(len=*) :: matname
      integer :: i, m

      write (*, '(2A)') trim(matname), ' = '
      m = size(mat, 1)
      do i = 1, m
         if (mat(i) .NE. dcmplx(0d0, 0d0)) write (*, '(F8.4,A,F8.4,A)') real(mat(i)), &
            ' + ', imag(mat(i)), 'i'
      end do

   end subroutine print_cvec

! CLEARTEXT READ **************************************************************
!******************************************************************************
   subroutine check_paramsfile(matrices)
! Subroutine to read input parameters
      type(data), intent(inout):: matrices
      integer ::  i
      character(len=80) :: arg

      do i = 1, command_argument_count(), 2
         call get_command_argument(i, arg)

         select case (arg)
         case ('-paramsfile')
            call get_command_argument(i + 1, matrices%paramsfile)
         end select
      end do
   end subroutine check_paramsfile

!******************************************************************************

   subroutine read_arguments(matrices, mesh)
! Subroutine to read input parameters
      type(data), intent(inout):: matrices
      type(mesh_struct), intent(inout) :: mesh
      integer :: i
      character(len=80) :: arg_name, arg

      matrices%singleT = 0

      do i = 1, command_argument_count(), 2
         call get_command_argument(i, arg_name)

         select case (arg_name)

         case ('-mesh')
            call get_command_argument(i + 1, arg)
            mesh%meshname = arg
            write (*, '(2A)') ' Mesh: ', trim(mesh%meshname)
         case ('-T')
            call get_command_argument(i + 1, arg)
            matrices%tname = arg
            write (*, '(2A)') ' T-matrix: ', trim(matrices%tname)
         case ('-log')
            call get_command_argument(i + 1, arg)
            matrices%out = arg
            write (*, '(2A)') ' Log: ', trim(matrices%out)
         case ('-mueller')
            call get_command_argument(i + 1, arg)
            matrices%mueller = arg
            write (*, '(2A)') ' Mueller: ', trim(matrices%mueller)
         case ('-mueller_mode')
            call get_command_argument(i + 1, arg)
            matrices%mueller_mode = trim(arg)
            write (*, '(2A)') ' Mueller mode: ', trim(matrices%mueller_mode)
         case ('-paramsfile')
            call get_command_argument(i + 1, arg)
            matrices%paramsfile = trim(arg)
         case ('-refr')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%refr
         case ('-refi')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%refi
         case ('-wb')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%whichbar
         case ('-singleT')
            call get_command_argument(i + 1, arg)
            read (arg, *) matrices%singleT
         case ('-seed')
            call get_command_argument(i + 1, arg)
            read (arg, *) seedling
            matrices%R = rand_rot()
         case ('-mie')
            call get_command_argument(i + 1, arg)
            read (arg, *) use_mie
         case ('-B')
            call get_command_argument(i + 1, arg)
            read (arg, *) calc_extra_torques

         case ('-help')
            print *, 'Command line parameters'
            print *, '-mesh mesh.h5          "Mesh geometry"'
            print *, '-T T.h5                "T-matrix file"'
            print *, '-log out/log           "Log to file"'
            print *, '-mueller out/mueller   "Mueller matrix to file"'
            print *, '-paramsfile params.in  "Read input parameters from file"'
            print *, '-refr 0.0              "Real part of refractive index"'
            print *, '-refi 0.0              "Imaginary part of refractive index"'
            print *, '-wb 0                  "Choose wavelength from the T-matrix"'
            print *, '-singleT 0             "Calculate only one T-matrix, of wb"'
            print *, '-mie 1                 "Use the Mie sphere"'
            print *, '-B 1                   "Use external magnetic field"'

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

!******************************************************************************

   subroutine read_params(matrices, mesh)
! Subroutine to read a input file for nifty usage
      type(data), intent(inout):: matrices
      type(mesh_struct), intent(inout) :: mesh
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
      character(len=4) :: R0
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
            case ('rho')
               read (buffer, *, iostat=ios) mesh%rho
            case ('a')
               read (buffer, *, iostat=ios) mesh%a
            case ('E')
               read (buffer, *, iostat=ios) matrices%E
            case ('v0')
               read (buffer, *, iostat=ios) matrices%v_CM
            case ('w0')
               read (buffer, *, iostat=ios) matrices%w
            case ('R0')
               read (buffer, *, iostat=ios) R0
               if (R0 == 'r') then
                  matrices%R = rand_rot()
               else if (R0 == 't') then
                  matrices%R = reshape(dble([0, 1, 0, 0, 0, 1, 1, 0, 0]), [3, 3])
               else
                  matrices%R = eye(3)
               end if
            case ('khi_0')
               read (buffer, *, iostat=ios) matrices%khi_0
            case ('dt')
               read (buffer, *, iostat=ios) matrices%dt0
            case ('it_max')
               read (buffer, *, iostat=ios) matrices%it_max
            case ('it_log')
               read (buffer, *, iostat=ios) matrices%it_log
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
            case ('Tmat')
               read (buffer, *, iostat=ios) matrices%Tmat
            case ('choose_integrator')
               read (buffer, *, iostat=ios) matrices%which_int
            case ('whichbar')
               read (buffer, *, iostat=ios) whichbar
               if (matrices%whichbar == 0) matrices%whichbar = whichbar
            case ('test_forces')
               read (buffer, *, iostat=ios) run_test
            case ('is_aggr')
               read (buffer, *, iostat=ios) matrices%is_aggr
            case ('mueller_mode')
               read (buffer, *, iostat=ios) mode
               if (trim(mode) /= 'none') matrices%mueller_mode = mode
            case ('waves')
               read (buffer, *, iostat=ios) matrices%waves
            case ('refr')
               read (buffer, *, iostat=ios) temp
               if (temp > 1d-7) matrices%refr = temp
            case ('refi')
               read (buffer, *, iostat=ios) tempii
            case ('N_mean')
               read (buffer, *, iostat=ios) matrices%N_mean
            case ('tol_m')
               read (buffer, *, iostat=ios) matrices%tol_m
            case ('B')
               read (buffer, *, iostat=ios) matrices%B_len
               if (matrices%B_len > 1d-14) calc_extra_torques = 1
            case ('B_psi')
               read (buffer, *, iostat=ios) matrices%B_psi
               matrices%B_psi = matrices%B_psi*pi/180
            case ('beam_shape')
               read (buffer, *, iostat=ios) beam_shape
            case ('pl')
               read (buffer, *, iostat=ios) p, l
            case default
               !print *, 'Skipping invalid label at line', line

            end select
         end if
      end do

      close (fh)

      if (temp > 1d-7) matrices%refi = tempii
      matrices%it_stop = matrices%it_max
      mesh%is_mesh = 1 - matrices%is_aggr
      matrices%B = matrices%B_len*(matmul(R_aa([0d0, 1d0, 0d0], matrices%B_psi), [0d0, 0d0, 1d0]))
      allocate (mesh%ki(matrices%bars))
      allocate (matrices%E_rel(matrices%bars))
      allocate (matrices%Nmaxs(matrices%bars))

   end subroutine read_params

!******************************************************************************

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

!******************************************************************************

   subroutine read_log(matrices)
      type(data), intent(inout):: matrices
      integer, parameter :: fh = 15
      integer :: line

! Control file variables
      real(dp) :: wlast(3), t_tresh, t1, t2, tf, w_av(3), Rab(3, 3)
      integer :: firstlineno, lastlineno, i, n1, n2
      character(len=140) :: dummy
      real(dp) :: x(3), v(3), w(3), J(3), F(3), N(3), t, R(9)

      firstlineno = 5023
      lastlineno = get_last_line_no(matrices%out)

      open (fh, file=trim(matrices%out))
      do i = 1, lastlineno - 1
         read (fh, *)
      end do

      read (fh, *) line, dummy, x, dummy, v, dummy, w, dummy, J, dummy, F, dummy, N, dummy, &
         t, dummy, R
!write(*,'(I0,3ES11.3)') line, w
!call print_mat(reshape(R,[3,3]), 'R')

      wlast = w/dsqrt(w(1)**2 + w(2)**2 + w(3)**2)
!write(*,'(3ES11.3)') wlast

      close (fh)

      t_tresh = 4d0*pi/dsqrt(w(1)**2 + w(2)**2 + w(3)**2)
      tf = t
      t1 = t - t_tresh
      t2 = t - 10*t_tresh
!print*, t1, t2

      w_av = 0d0
      open (fh, file=trim(matrices%out))
      do i = 1, lastlineno
         if (i >= firstlineno) then
            read (fh, *) line, dummy, x, dummy, v, dummy, w, dummy, J, dummy, F, dummy, N, dummy, &
               t, dummy, R
            w_av = w_av + w
            if (t > t2) then
!     print*, line
               t2 = tf
            end if
            if (t > t1) then
!     print*, line
               t1 = tf
            end if
            if (i == firstlineno) then
!     print*, line
               n1 = line
            end if
         else
            read (fh, *)
         end if
      end do
      n2 = line
      w_av = w_av/(lastlineno - firstlineno + 1)
!print*, 'w_av = ', w_av
      close (fh)
      Rab = rotate_a_to_b(w_av, [1d0, 0d0, w_av(3)])
!call print_mat(Rab,'R_ab')
      matrices%R_al = Rab
      allocate (matrices%RRR(3, 3, n2 - n1 + 1))
      allocate (matrices%www(3, n2 - n1 + 1))

      open (fh, file=trim(matrices%out))
      do i = 1, lastlineno
         if (i >= firstlineno) then
            read (fh, *) line, dummy, x, dummy, v, dummy, w, dummy, J, dummy, F, dummy, N, dummy, &
               t, dummy, R
            matrices%RRR(:, :, i - firstlineno + 1) = reshape(R, [3, 3])
            matrices%www(:, i - firstlineno + 1) = w
         else
            read (fh, *)
         end if
      end do

      close (fh)

   end subroutine read_log

!******************************************************************************

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
!******************************************************************************

   subroutine write_mueller(A, fname)
      real(dp), intent(in) :: A(:, :)
      character(len=80), intent(in) :: fname
      integer(HSIZE_T), dimension(2) :: dims ! Dataset dimensions
      integer     :: i

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

!******************************************************************************

   subroutine write_RT_matrix(A, fname, type)
      real(dp), intent(in) :: A(:, :)
      character(len=80), intent(in) :: fname
      integer(HSIZE_T), dimension(2) :: dims ! Dataset dimensions

      integer     :: i, type

      dims = [int(size(A, 1), 8), int(size(A, 2), 8)]
      open (unit=1, file=trim(fname), ACTION="write", STATUS="replace")

      if (type == 1) then
         write (1, '(19A)') ' N_size ', '     N_ia ', '    N_pts ', '      S11  ', '         S12  ', &
            '         S13  ', '         S14  ', '         S21  ', '         S22  ', '         S23  ', &
            '         S24  ', '         S31  ', '         S32  ', '         S33  ', '         S34  ', &
            '         S41  ', '         S42  ', '         S43  ', '         S44  '
      else if (type == 2) then
         write (1, '(19A)') ' N_size ', '     N_ia ', '    N_pts ', '      K11  ', '         K12  ', &
            '         K13  ', '         K14  ', '         K21  ', '         K22  ', '         K23  ', &
            '         K24  ', '         K31  ', '         K32  ', '         K33  ', '         K34  ', &
            '         K41  ', '         K42  ', '         K43  ', '         K44  '
      end if

      do i = 1, int(dims(1))
         write (1, '(3I8,16ES14.6)') int(A(i, 1:3)), A(i, 4:19)
      end do

      close (1)

   end subroutine write_RT_matrix

!******************************************************************************

   subroutine init_values(matrices)
      type(data) :: matrices

      matrices%k_orig = matrices%khat
      matrices%E0_orig = real(matrices%E0hat)
      matrices%E90_orig = real(matrices%E90hat)

      matrices%M1 = 0d0
      matrices%M3 = 0d0
      matrices%is_aligned = 0
      matrices%alignment_found = 0

      matrices%xn = matrices%x_CM ! Input in lab frame (if needed)
      matrices%vn = matrices%v_CM ! Input in lab frame
      matrices%wn = matrices%w
      matrices%dt = matrices%dt0
      matrices%tt = 0d0
      matrices%q = mat2quat(matrices%R)
      matrices%N = dble([0d0, 0d0, 0d0])
      matrices%F = dble([0d0, 0d0, 0d0])
      matrices%J = matmul(matmul(matrices%P, matrices%R), matmul(matrices%I, matrices%w)) ! In lab frame!

! R_expansion, for planewave expansion rotation from (0, 0, 1) position to khat
      matrices%Rexp = find_Rexp(matrices)
      matrices%Rexp = transpose(matrices%Rexp)
      matrices%R90_init = reshape(dble([0d0, 1d0, 0d0, -1d0, 0d0, 0d0, 0d0, 0d0, 1d0]), [3, 3])

   end subroutine init_values

!******************************************************************************

   subroutine update_values(matrices)
      type(data) :: matrices

      matrices%x_CM = matrices%xn ! In lab frame!
      matrices%v_CM = matrices%vn ! In lab frame!
      matrices%w = matrices%wn ! In body frame!
      matrices%R = matrices%Rn
      matrices%q = matrices%qn
      matrices%tt = matrices%tt + matrices%dt
      matrices%J = matmul(matmul(matrices%P, matrices%R), &
                          matmul(matrices%I, matrices%w))
   end subroutine update_values

!******************************************************************************

   subroutine start_log(fname, matrices, mesh)
      type(data) :: matrices
      type(mesh_struct), intent(in) :: mesh
      integer :: i1
      character(len=80) :: fname

      matrices%buffer = 0

!if(file_exists(fname)) then
! print*, "Log exists, shutting down... (Remember to disable this later!)"
! stop
!end if

      open (unit=1, file=fname, ACTION="write", STATUS="replace")

      write (1, '(A, A)') 'meshname =   ', mesh%meshname
      write (1, '(A, 3f7.3)') 'k_hat    = ', matrices%khat
      write (1, '(A, 3ES11.3)') 'dt       = ', matrices%dt
      write (1, '(A,I20)') 'Nmax     = ', matrices%it_max
      write (1, '(A, 3ES11.3)') 'w0       = ', matmul(matmul(matrices%P, matrices%R), &
                                                      matrices%w)
      write (1, '(A, 9f7.3)') 'R0       = ', matrices%R
      write (1, '(A, 3ES11.3)') 'rho      = ', mesh%rho
      write (1, '(A, 3ES11.3)') 'CM       = ', matrices%CM
      write (1, '(A, 3ES11.3)') 'a        = ', mesh%a
      write (1, '(A, 80ES11.3)') 'ki       = ', mesh%ki
      write (1, '(A, 2f8.2)') 'range nm = ', matrices%lambda1*1d9, matrices%lambda2*1d9
      write (1, '(A, 3ES11.3)') 'mass     = ', mesh%mass
      write (1, '(A, 3ES11.3)') 'V        = ', mesh%V
      write (1, '(A, 3ES11.3)') 'rot_max  = ', matrices%rot_max
      write (1, '(A)') 'Ip       = '
      write (1, '(9ES11.3)') matrices%Ip(:)
      write (1, '(A)') 'Q        = '
      do i1 = 1, 3
         write (1, '(9f7.3)') matrices%P(i1, :)
      end do
      write (1, '(A)') ' Log: n | x | v | w | J | N | F | t | R '
      write (1, '(A)') ' '

      close (1)

   end subroutine start_log

!******************************************************************************

   subroutine append_log(fname, n, matrices)
      type(data):: matrices
      integer :: n, i, md
      character(len=80) :: fname, fmt

      fmt = '(I0, 6(A, 3ES11.3), A, 1ES16.8, A, 9f7.3)'
      md = mod(n, 1000)

      if (n > 100000) then
         if (matrices%is_aligned == 0) then
            matrices%is_aligned = alignment_state(matrices)
         else if (matrices%alignment_found == 0) then
            print *, " ******************* HEY, IT'S ALIGNED! ********************"
            matrices%it_stop = n + matrices%it_log
            matrices%alignment_found = 1
         end if
      end if

      if (md == 0) md = 1000
      matrices%x_buf(:, md) = matrices%x_CM
      matrices%v_buf(:, md) = matrices%v_CM
      matrices%w_buf(:, md) = matrices%w
      matrices%J_buf(:, md) = matrices%J
      matrices%N_buf(:, md) = matrices%N
      matrices%F_buf(:, md) = matrices%F
      matrices%t_buf(:, md) = matrices%tt
      matrices%R_buf(:, :, md) = matrices%R

      md = mod(n, 1000)
      if (md == 0 .OR. n == matrices%it_stop) then
         open (unit=1, file=fname, action="write", position="append", STATUS="old")
         if (n > matrices%it_stop - matrices%it_log .OR. matrices%it_log == 0) then
            do i = 1, 1000
               if (i + matrices%buffer <= matrices%it_stop) then
                  write (1, fmt) i + matrices%buffer, ' |', &
                     matrices%x_buf(:, i), ' |', &
                     matrices%v_buf(:, i), ' |', &
                     matmul(matmul(matrices%P, matrices%R_buf(:, :, i)), matrices%w_buf(:, i)), ' |', &
                     matrices%J_buf(:, i), ' |', &
                     matrices%N_buf(:, i), ' |', &
                     matrices%F_buf(:, i), ' |', &
                     matrices%t_buf(:, i), ' |', &
                     matrices%R_buf(:, :, i)
               end if
            end do
         end if
         close (1)
         matrices%buffer = n
      end if

   end subroutine append_log

! HDF5 READ *******************************************************************
!******************************************************************************

   subroutine read_mesh(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh
      character(len=80) :: file

!character(len=8), PARAMETER :: file = "mesh.h5"
      character(len=5), PARAMETER :: coord_dataset = "coord"
      character(len=6), PARAMETER :: etopol_dataset = "etopol"
      character(len=7), PARAMETER :: param_r_dataset = "param_r"
      character(len=7), PARAMETER :: param_i_dataset = "param_i"

      integer(HID_T) :: file_id
      integer(HID_T) :: coord_dataset_id, etopol_dataset_id, param_r_dataset_id, param_i_dataset_id
      integer(HID_T) :: coord_dataspace_id, etopol_dataspace_id, param_r_dataspace_id, param_i_dataspace_id
      integer(HSIZE_T), dimension(2) :: dims_out, coord_dims, etopol_dims, param_r_dims, param_i_dims

      integer :: error

      real(dp), dimension(:, :), allocatable :: coord
      integer, dimension(:, :), allocatable :: etopol
      real(dp), dimension(:, :), allocatable :: param_r, param_i
      real(dp) :: vol, a_eff

      file = mesh%meshname
      call h5open_f(error)
      call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error)

!******************************************************************************

      call h5dopen_f(file_id, coord_dataset, coord_dataset_id, error)
      call h5dget_space_f(coord_dataset_id, coord_dataspace_id, error)
      call H5sget_simple_extent_dims_f(coord_dataspace_id, dims_out, coord_dims, error)
      allocate (coord(coord_dims(1), coord_dims(2)))
      call h5dread_f(coord_dataset_id, H5T_NATIVE_DOUBLE, coord, coord_dims, error)
      call h5dclose_f(coord_dataset_id, error)

!******************************************************************************

      call h5dopen_f(file_id, etopol_dataset, etopol_dataset_id, error)
      call h5dget_space_f(etopol_dataset_id, etopol_dataspace_id, error)
      call H5sget_simple_extent_dims_f(etopol_dataspace_id, dims_out, etopol_dims, error)
      allocate (etopol(etopol_dims(1), etopol_dims(2)))
      call h5dread_f(etopol_dataset_id, H5T_NATIVE_integer, etopol, etopol_dims, error)
      call h5dclose_f(etopol_dataset_id, error)

!******************************************************************************

      call h5dopen_f(file_id, param_r_dataset, param_r_dataset_id, error)
      call h5dget_space_f(param_r_dataset_id, param_r_dataspace_id, error)
      call H5sget_simple_extent_dims_f(param_r_dataspace_id, dims_out, param_r_dims, error)
      if (param_r_dims(1) < int(size(etopol, 2), 8)) then
         allocate (param_r(size(etopol, 2), param_r_dims(1)))
      else
         allocate (param_r(param_r_dims(1), 1))
      end if
      call h5dread_f(param_r_dataset_id, H5T_NATIVE_DOUBLE, param_r, param_r_dims, error)
      call h5dclose_f(param_r_dataset_id, error)

!******************************************************************************

      call h5dopen_f(file_id, param_i_dataset, param_i_dataset_id, error)
      call h5dget_space_f(param_i_dataset_id, param_i_dataspace_id, error)
      call H5sget_simple_extent_dims_f(param_i_dataspace_id, dims_out, param_i_dims, error)
      if (param_i_dims(1) < int(size(etopol, 2), 8)) then
         allocate (param_i(size(etopol, 2), param_i_dims(1)))
      else
         allocate (param_i(param_i_dims(1), 1))
      end if
      call h5dread_f(param_i_dataset_id, H5T_NATIVE_DOUBLE, param_i, param_i_dims, error)
      call h5dclose_f(param_i_dataset_id, error)

!******************************************************************************

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      if (allocated(mesh%coord)) deallocate (mesh%coord, mesh%etopol, mesh%param, mesh%params)
      allocate (mesh%coord(3, size(coord, 2)))
      mesh%coord = coord
      print *, '   Number of nodes      =', size(coord, 2)
      mesh%N_node = size(coord, 2)

      allocate (mesh%etopol(4, size(etopol, 2)))
      mesh%etopol = etopol
      print *, '   Number of elements   =', size(etopol, 2)
      mesh%N_tet = size(etopol, 2)

      allocate (mesh%params(size(param_r, 1), size(param_r, 2)))
      allocate (mesh%param(size(param_r, 1)))
      mesh%params = dcmplx(param_r, param_i)
      if (matrices%refr > 1d-7) then
         mesh%param = dcmplx(matrices%refr**2 - matrices%refi**2, &
                             2d0*matrices%refr*matrices%refi)
         write (*, '(2(A,F5.3))') '    Refractive index     =   ', matrices%refr, ' + i', matrices%refi
         write (*, '(2(A,F5.3))') '    Dielectric constant  =   ', matrices%refr**2 - matrices%refi**2, &
            ' + i', 2d0*matrices%refr*matrices%refi
      else
         mesh%param = dcmplx(param_r(:, 1), param_i(:, 1))
         if (maxval(param_r(:, 1)) - minval(param_r(:, 1)) > 1d-7 .OR. maxval(param_i(:, 1)) - minval(param_i(:, 1)) > 1d-7) then
         else
            write (*, '(2(A,F5.3))') '    Dielectric constant  =   ', param_r(1, 1), ' + i', param_i(1, 1)
         end if
      end if

      vol = get_tetra_vol(mesh)
      a_eff = (3d0*vol/4d0/pi)**(1d0/3d0)

      mesh%coord = mesh%coord*mesh%a/a_eff ! Scale coordinates to correspond the real grain
! With correction to scaling (effective radius is defined by volumes)

   end subroutine read_mesh

!******************************************************************************

   subroutine read_aggr(mesh)
      type(mesh_struct) :: mesh
      character(len=16), PARAMETER :: dataset1 = "coord" ! Dataset name
      character(len=16), PARAMETER :: dataset2 = "radius" ! Dataset name

      integer(HID_T) :: file_id, dataset1_id, dataset2_id, dataspace1_id, dataspace2_id
      integer(HSIZE_T), dimension(2) :: coord_dims, dims_out
      integer :: error

      real(dp), dimension(:, :), allocatable :: coord
      real(dp), dimension(:), allocatable :: radius

      call h5open_f(error) ! initialize interface
      call h5fopen_f(mesh%meshname, H5F_ACC_RDWR_F, file_id, error) ! open file

!* COORD **********************************************************************

      call h5dopen_f(file_id, dataset1, dataset1_id, error)
      call h5dget_space_f(dataset1_id, dataspace1_id, error)
      call H5sget_simple_extent_dims_f(dataspace1_id, dims_out, coord_dims, error)
      allocate (coord(coord_dims(1), coord_dims(2)))
      call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, coord, coord_dims, error)
      call h5dclose_f(dataset1_id, error)

!* RADIUS *********************************************************************

      call h5dopen_f(file_id, dataset2, dataset2_id, error)
      call h5dget_space_f(dataset2_id, dataspace2_id, error)
      call H5sget_simple_extent_dims_f(dataspace2_id, dims_out, coord_dims, error)
      allocate (radius(coord_dims(2)))
      call h5dread_f(dataset2_id, H5T_NATIVE_DOUBLE, radius, coord_dims, error)
      call h5dclose_f(dataset2_id, error)

!******************************************************************************

      call h5fclose_f(file_id, error) ! close file
      call h5close_f(error) ! close inteface

      allocate (mesh%coord(size(coord, 1), size(coord, 2)))
      allocate (mesh%radius(size(radius)))

      mesh%coord = coord
      mesh%radius = radius

   end subroutine read_aggr

!******************************************************************************

   subroutine read_T(matrices, mesh)
!character(len=28), intent(in) :: fname
      type(data) :: matrices
      type(mesh_struct) :: mesh

      character(len=80) :: file ! File name
      character(len=16), PARAMETER :: dataset1 = "Taai_r"
      character(len=16), PARAMETER :: dataset2 = "Taai_i"
      character(len=16), PARAMETER :: dataset3 = "Tabi_r"
      character(len=16), PARAMETER :: dataset4 = "Tabi_i"
      character(len=16), PARAMETER :: dataset5 = "Tbai_r"
      character(len=16), PARAMETER :: dataset6 = "Tbai_i"
      character(len=16), PARAMETER :: dataset7 = "Tbbi_r"
      character(len=16), PARAMETER :: dataset8 = "Tbbi_i"
      character(len=16), PARAMETER :: dataset9 = "T-ref-a"
      character(len=16), PARAMETER :: dataset10 = "T-wavlens"

      integer(HID_T) :: file_id
      integer(HID_T) :: dataset1_id, dataset2_id, dataset3_id, dataset4_id
      integer(HID_T) :: dataset5_id, dataset6_id, dataset7_id, dataset8_id
      integer(HID_T) :: dataset9_id, dataset10_id
      integer(HID_T) :: dataspace_id

      integer(HSIZE_T), dimension(3) :: dims_out, dims
      integer(HSIZE_T), dimension(1) :: dimswl, dimsinfo, dims1

      integer :: error, i

      real(dp), dimension(3) :: ref_a
      real(dp), dimension(:), allocatable :: wls, param_r, param_i
      real(dp), dimension(:, :, :), allocatable :: Taai_r, Taai_i
      real(dp), dimension(:, :, :), allocatable :: Tabi_r, Tabi_i
      real(dp), dimension(:, :, :), allocatable :: Tbai_r, Tbai_i
      real(dp), dimension(:, :, :), allocatable :: Tbbi_r, Tbbi_i
      complex(dp), dimension(:, :, :), allocatable :: Taai, Tabi, Tbai, Tbbi
      complex(dp) :: ref
      LOGICAL :: exists

      file = matrices%tname

      call h5open_f(error)
      call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error)

!******************************************************************************

      call h5dopen_f(file_id, dataset1, dataset1_id, error)
      call h5dget_space_f(dataset1_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      allocate (Taai_r(dims(1), dims(2), dims(3)))
      call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, Taai_r, dims, error)
      call h5dclose_f(dataset1_id, error)

      call h5dopen_f(file_id, dataset2, dataset2_id, error)
      call h5dget_space_f(dataset2_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      allocate (Taai_i(dims(1), dims(2), dims(3)))
      call h5dread_f(dataset2_id, H5T_NATIVE_DOUBLE, Taai_i, dims, error)
      call h5dclose_f(dataset2_id, error)

      allocate (Taai(size(Taai_r, 1), size(Taai_r, 1), dims(3)))
      Taai = dcmplx(Taai_r, Taai_i)

!******************************************************************************

      call h5dopen_f(file_id, dataset3, dataset3_id, error)
      call h5dget_space_f(dataset3_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      allocate (Tabi_r(dims(1), dims(2), dims(3)))
      call h5dread_f(dataset3_id, H5T_NATIVE_DOUBLE, Tabi_r, dims, error)
      call h5dclose_f(dataset3_id, error)

      call h5dopen_f(file_id, dataset4, dataset4_id, error)
      call h5dget_space_f(dataset4_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      allocate (Tabi_i(dims(1), dims(2), dims(3)))
      call h5dread_f(dataset4_id, H5T_NATIVE_DOUBLE, Tabi_i, dims, error)
      call h5dclose_f(dataset4_id, error)

      allocate (Tabi(size(Tabi_r, 1), size(Tabi_r, 1), dims(3)))
      Tabi = dcmplx(Tabi_r, Tabi_i)

!******************************************************************************

      call h5dopen_f(file_id, dataset5, dataset5_id, error)
      call h5dget_space_f(dataset5_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      allocate (Tbai_r(dims(1), dims(2), dims(3)))
      call h5dread_f(dataset5_id, H5T_NATIVE_DOUBLE, Tbai_r, dims, error)
      call h5dclose_f(dataset5_id, error)

      call h5dopen_f(file_id, dataset6, dataset6_id, error)
      call h5dget_space_f(dataset6_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      allocate (Tbai_i(dims(1), dims(2), dims(3)))
      call h5dread_f(dataset6_id, H5T_NATIVE_DOUBLE, Tbai_i, dims, error)
      call h5dclose_f(dataset6_id, error)

      allocate (Tbai(size(Tbai_r, 1), size(Tbai_r, 1), dims(3)))
      Tbai = dcmplx(Tbai_r, Tbai_i)

!******************************************************************************

      call h5dopen_f(file_id, dataset7, dataset7_id, error)
      call h5dget_space_f(dataset7_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      allocate (Tbbi_r(dims(1), dims(2), dims(3)))
      call h5dread_f(dataset7_id, H5T_NATIVE_DOUBLE, Tbbi_r, dims, error)
      call h5dclose_f(dataset7_id, error)

      call h5dopen_f(file_id, dataset8, dataset8_id, error)
      call h5dget_space_f(dataset8_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      allocate (Tbbi_i(dims(1), dims(2), dims(3)))
      call h5dread_f(dataset8_id, H5T_NATIVE_DOUBLE, Tbbi_i, dims, error)
      call h5dclose_f(dataset8_id, error)

      allocate (Tbbi(size(Tbbi_r, 1), size(Tbbi_r, 1), dims(3)))
      Tbbi = dcmplx(Tbbi_r, Tbbi_i)

!******************************************************************************

      call h5lexists_f(file_id, dataset9, exists, error)
      if (exists) then
         call h5dopen_f(file_id, dataset9, dataset9_id, error)
         call h5dget_space_f(dataset9_id, dataspace_id, error)
         call H5sget_simple_extent_dims_f(dataspace_id, dimsinfo, dims1, error)

         call h5dread_f(dataset9_id, H5T_NATIVE_DOUBLE, ref_a, dims1, error)
         call h5dclose_f(dataset9_id, error)

         mesh%a = ref_a(3)
         allocate (param_r(size(mesh%param)), param_i(size(mesh%param)))
         param_r = real(mesh%param)
         param_i = imag(mesh%param)
         if (maxval(param_r) - minval(param_r) > 1d-7 .OR. maxval(param_i) - minval(param_i) > 1d-7) then
            write (*, '(A)') '    Note: Geometry is inhomogeneous.'
         else
            mesh%param = dcmplx(ref_a(1), ref_a(2))
            ref = sqrt(dcmplx(ref_a(1), ref_a(2)))
            matrices%refr = real(ref)
            matrices%refi = imag(ref)
         end if

         !****************************************************************************

         call h5dopen_f(file_id, dataset10, dataset10_id, error)
         call h5dget_space_f(dataset10_id, dataspace_id, error)
         call H5sget_simple_extent_dims_f(dataspace_id, dimswl, dims1, error)
         allocate (wls(dims1(1)))
         call h5dread_f(dataset10_id, H5T_NATIVE_DOUBLE, wls, dims1, error)
         call h5dclose_f(dataset10_id, error)

         mesh%ki = 2d0*pi/wls
         do i = 1, size(wls, 1)
!   num = 0
!   do j = 1,size(Taai,1)
!     if(real(Taai(j,1,i))/=0d0) num = num + 1
!   end do
!   matrices%Nmaxs(i) = int(dsqrt(real(num)+1d0)-1d0)
            matrices%Nmaxs(i) = truncation_order(mesh%ki(i)*mesh%a)
         end do
         !****************************************************************************
      end if

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      if (allocated(matrices%Taai)) deallocate (matrices%Taai, matrices%Tabi, &
                                                matrices%Tbai, matrices%Tbbi)

      allocate (matrices%Taai(size(Taai, 1), size(Taai, 2), size(Taai, 3)), &
                matrices%Tabi(size(Tabi, 1), size(Tabi, 2), size(Tabi, 3)), &
                matrices%Tbai(size(Tbai, 1), size(Tbai, 2), size(Tbai, 3)), &
                matrices%Tbbi(size(Tbbi, 1), size(Tbbi, 2), size(Tbbi, 3)))

      matrices%Taai = Taai
      matrices%Tabi = Tabi
      matrices%Tbai = Tbai
      matrices%Tbbi = Tbbi

   end subroutine read_T

! HDF5 WRITE ******************************************************************
!******************************************************************************
   subroutine write_T(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh

      if (matrices%singleT == 1) then
         call singleT_write2file(matrices)
      else
         call T_write2file(matrices, mesh)
      end if

   end subroutine write_T

!******************************************************************************

   subroutine T_write2file(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh

      character(len=80) :: fname
      character(len=80) :: filename

      character(len=6), PARAMETER :: dsetname1 = "Taai_r"
      character(len=6), PARAMETER :: dsetname2 = "Taai_i"
      character(len=6), PARAMETER :: dsetname3 = "Tabi_r"
      character(len=6), PARAMETER :: dsetname4 = "Tabi_i"
      character(len=6), PARAMETER :: dsetname5 = "Tbai_r"
      character(len=6), PARAMETER :: dsetname6 = "Tbai_i"
      character(len=6), PARAMETER :: dsetname7 = "Tbbi_r"
      character(len=6), PARAMETER :: dsetname8 = "Tbbi_i"
      character(len=16), PARAMETER :: dsetname9 = "T-ref-a"
      character(len=16), PARAMETER :: dsetname10 = "T-wavlens"

      integer(HID_T) :: file_id
      integer(HID_T) :: dset_id1
      integer(HID_T) :: dset_id2
      integer(HID_T) :: dset_id3
      integer(HID_T) :: dset_id4
      integer(HID_T) :: dset_id5
      integer(HID_T) :: dset_id6
      integer(HID_T) :: dset_id7
      integer(HID_T) :: dset_id8
      integer(HID_T) :: dset_id9
      integer(HID_T) :: dset_id10

      integer(HID_T) :: dspace_id, dspace_id2, dspace_id3

      integer(HSIZE_T), dimension(3) :: dims
      integer(HSIZE_T), dimension(1) :: dimswl, dimsinfo
      integer     ::    rank = 3
      integer     ::   error ! Error flag

      fname = matrices%tname

      filename = fname
      dims = int8((/size(matrices%Taai, 1), size(matrices%Taai, 2), size(matrices%Taai, 3)/))
      dimswl = int8((/matrices%bars/))
      dimsinfo = int8((/3/))

      CALL h5open_f(error)
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      CALL h5screate_simple_f(1, dimsinfo, dspace_id2, error)
      CALL h5screate_simple_f(1, dimswl, dspace_id3, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(matrices%Taai), dims, error)
      CALL h5dclose_f(dset_id1, error)
      CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(matrices%Taai), dims, error)
      CALL h5dclose_f(dset_id2, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id3, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, real(matrices%Tabi), dims, error)
      CALL h5dclose_f(dset_id3, error)
      CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id4, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, imag(matrices%Tabi), dims, error)
      CALL h5dclose_f(dset_id4, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id5, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, real(matrices%Tbai), dims, error)
      CALL h5dclose_f(dset_id5, error)
      CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id6, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, imag(matrices%Tbai), dims, error)
      CALL h5dclose_f(dset_id6, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname7, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id7, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, real(matrices%Tbbi), dims, error)
      CALL h5dclose_f(dset_id7, error)
      CALL h5dcreate_f(file_id, dsetname8, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id8, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, imag(matrices%Tbbi), dims, error)
      CALL h5dclose_f(dset_id8, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname9, H5T_NATIVE_DOUBLE, dspace_id2, &
                       dset_id9, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id9, H5T_NATIVE_DOUBLE, [real(mesh%param(1)), &
                                                    imag(mesh%param(1)), mesh%a], dimsinfo, error)
      CALL h5dclose_f(dset_id9, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname10, H5T_NATIVE_DOUBLE, dspace_id3, &
                       dset_id10, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id10, H5T_NATIVE_DOUBLE, 2d0*pi/mesh%ki, dimswl, error)
      CALL h5dclose_f(dset_id10, error)

!******************************************************************

      CALL h5sclose_f(dspace_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   end subroutine T_write2file

!******************************************************************************

   subroutine T_empty(matrices, mesh)
      type(data) :: matrices
      type(mesh_struct) :: mesh

      character(len=80) :: fname
      character(len=80) :: filename

      character(len=6), PARAMETER :: dsetname1 = "Taai_r"
      character(len=6), PARAMETER :: dsetname2 = "Taai_i"
      character(len=6), PARAMETER :: dsetname3 = "Tabi_r"
      character(len=6), PARAMETER :: dsetname4 = "Tabi_i"
      character(len=6), PARAMETER :: dsetname5 = "Tbai_r"
      character(len=6), PARAMETER :: dsetname6 = "Tbai_i"
      character(len=6), PARAMETER :: dsetname7 = "Tbbi_r"
      character(len=6), PARAMETER :: dsetname8 = "Tbbi_i"
      character(len=16), PARAMETER :: dsetname9 = "T-ref-a"
      character(len=16), PARAMETER :: dsetname10 = "T-wavlens"

      integer(HID_T) :: file_id
      integer(HID_T) :: dset_id1
      integer(HID_T) :: dset_id2
      integer(HID_T) :: dset_id3
      integer(HID_T) :: dset_id4
      integer(HID_T) :: dset_id5
      integer(HID_T) :: dset_id6
      integer(HID_T) :: dset_id7
      integer(HID_T) :: dset_id8
      integer(HID_T) :: dset_id9
      integer(HID_T) :: dset_id10

      integer(HID_T) :: dspace_id, dspace_id2, dspace_id3

      integer(HSIZE_T), dimension(3) :: dims
      integer(HSIZE_T), dimension(1) :: dimswl, dimsinfo
      integer     ::    rank = 3
      integer     ::   error ! Error flag

      complex(dp), dimension(:, :, :), allocatable    ::   emptyT

      fname = matrices%tname

      filename = fname
      dims = int8((/size(matrices%Taai, 1), size(matrices%Taai, 2), size(matrices%Taai, 3)/))
      dimswl = int8((/matrices%bars/))
      dimsinfo = int8((/3/))
      allocate (emptyT(dims(1), dims(2), dims(3)))
      emptyT = dcmplx(0d0, 0d0)

      CALL h5open_f(error)
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      CALL h5screate_simple_f(1, dimsinfo, dspace_id2, error)
      CALL h5screate_simple_f(1, dimswl, dspace_id3, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id1, error)
      CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id2, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id3, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id3, error)
      CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id4, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id4, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id5, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id5, error)
      CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id6, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id6, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname7, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id7, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id7, error)
      CALL h5dcreate_f(file_id, dsetname8, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id8, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id8, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname9, H5T_NATIVE_DOUBLE, dspace_id2, &
                       dset_id9, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id9, H5T_NATIVE_DOUBLE, [real(mesh%param(1)), &
                                                    imag(mesh%param(1)), mesh%a], dimsinfo, error)
      CALL h5dclose_f(dset_id9, error)

!******************************************************************

      CALL h5dcreate_f(file_id, dsetname10, H5T_NATIVE_DOUBLE, dspace_id3, &
                       dset_id10, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id10, H5T_NATIVE_DOUBLE, 2d0*pi/mesh%ki, dimswl, error)
      CALL h5dclose_f(dset_id10, error)

!******************************************************************

      CALL h5sclose_f(dspace_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   end subroutine T_empty

!*******************************************************************************

   subroutine singleT_write2file(matrices)
      type(data) :: matrices

      character(len=80) :: filename

      character(len=6), PARAMETER :: dsetname1 = "Taai_r"
      character(len=6), PARAMETER :: dsetname2 = "Taai_i"
      character(len=6), PARAMETER :: dsetname3 = "Tabi_r"
      character(len=6), PARAMETER :: dsetname4 = "Tabi_i"
      character(len=6), PARAMETER :: dsetname5 = "Tbai_r"
      character(len=6), PARAMETER :: dsetname6 = "Tbai_i"
      character(len=6), PARAMETER :: dsetname7 = "Tbbi_r"
      character(len=6), PARAMETER :: dsetname8 = "Tbbi_i"

      integer(HID_T) :: file_id
      integer(HID_T) :: dset_id1
      integer(HID_T) :: dset_id2
      integer(HID_T) :: dset_id3
      integer(HID_T) :: dset_id4
      integer(HID_T) :: dset_id5
      integer(HID_T) :: dset_id6
      integer(HID_T) :: dset_id7
      integer(HID_T) :: dset_id8

      integer(HID_T) :: dspace_id

      integer(HSIZE_T), dimension(3) :: dims, dims_out
      real(dp), dimension(:, :, :), allocatable :: Ti_r, Ti_i
      integer     ::   error ! Error flag

      filename = matrices%tname
      dims = int8((/size(matrices%Taai, 1), size(matrices%Taai, 2), size(matrices%Taai, 3)/))

      call h5open_f(error)
      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

!******************************************************************

      call h5dopen_f(file_id, dsetname1, dset_id1, error)
      call h5dget_space_f(dset_id1, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      allocate (Ti_r(dims(1), dims(2), dims(3)))
      call h5dread_f(dset_id1, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(:, :, matrices%whichbar) = real(matrices%Taai(:, :, matrices%whichbar))
      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id1, error)

      call h5dopen_f(file_id, dsetname2, dset_id2, error)
      call h5dget_space_f(dset_id2, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      allocate (Ti_i(dims(1), dims(2), dims(3)))
      call h5dread_f(dset_id2, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(:, :, matrices%whichbar) = imag(matrices%Taai(:, :, matrices%whichbar))
      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id2, error)

!******************************************************************

      call h5dopen_f(file_id, dsetname3, dset_id3, error)
      call h5dget_space_f(dset_id3, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id3, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(:, :, matrices%whichbar) = real(matrices%Tabi(:, :, matrices%whichbar))
      CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id3, error)

      call h5dopen_f(file_id, dsetname4, dset_id4, error)
      call h5dget_space_f(dset_id4, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id4, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(:, :, matrices%whichbar) = imag(matrices%Tabi(:, :, matrices%whichbar))
      CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id4, error)

!******************************************************************

      call h5dopen_f(file_id, dsetname5, dset_id5, error)
      call h5dget_space_f(dset_id5, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id5, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(:, :, matrices%whichbar) = real(matrices%Tbai(:, :, matrices%whichbar))
      CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id5, error)

      call h5dopen_f(file_id, dsetname6, dset_id6, error)
      call h5dget_space_f(dset_id6, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id6, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(:, :, matrices%whichbar) = imag(matrices%Tbai(:, :, matrices%whichbar))
      CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id6, error)

!******************************************************************

      call h5dopen_f(file_id, dsetname7, dset_id7, error)
      call h5dget_space_f(dset_id7, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id7, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(:, :, matrices%whichbar) = real(matrices%Tbbi(:, :, matrices%whichbar))
      CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id7, error)

      call h5dopen_f(file_id, dsetname8, dset_id8, error)
      call h5dget_space_f(dset_id8, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id8, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(:, :, matrices%whichbar) = imag(matrices%Tbbi(:, :, matrices%whichbar))
      CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id8, error)

!******************************************************************

      CALL h5sclose_f(dspace_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   end subroutine singleT_write2file

!*******************************************************************************

   subroutine write2file(A, fname)
      complex(dp), intent(in) :: A(:, :)
      character(len=80):: fname
      character(len=80) :: filename
!character(len=7), PARAMETER :: filename = "A.h5" ! File name
      character(len=3), PARAMETER :: dsetname1 = "A_r" ! Dataset name
      character(len=3), PARAMETER :: dsetname2 = "A_i" ! Dataset name

      integer(HID_T) :: file_id ! File identifier
      integer(HID_T) :: dset_id1 ! Dataset identifier
      integer(HID_T) :: dset_id2 ! Dataset identifier
      integer(HID_T) :: dspace_id ! Dataspace identifier

      integer(HSIZE_T), dimension(2) :: dims ! Dataset dimensions
      integer     ::    rank = 2 ! Dataset rank

      integer     ::   error ! Error flag
!integer     :: i, j

      filename = fname
      dims = [int(size(A, 1), 8), int(size(A, 2), 8)]

      CALL h5open_f(error)

!
! Create a new file using default properties.
!
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
      CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
      CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(A), dims, error)

      CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(A), dims, error)

!
! End access to the dataset and release resources used by it.
!
      CALL h5dclose_f(dset_id1, error)
      CALL h5dclose_f(dset_id2, error)
!
! Terminate access to the data space.
!
      CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
      CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
      CALL h5close_f(error)

   end subroutine write2file

end module
