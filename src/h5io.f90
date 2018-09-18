module h5io
   use common
   use hdf5

   implicit none
   contains
   ! HDF5 READ
   ! HDF5 WRITE

   ! HDF5 READ *******************************************************************
!****************************************************************************80

   subroutine read_mesh()
      character(len=80) :: file

      ! character(len=8), PARAMETER :: file = "mesh.h5"
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

!****************************************************************************80

      call h5dopen_f(file_id, coord_dataset, coord_dataset_id, error)
      call h5dget_space_f(coord_dataset_id, coord_dataspace_id, error)
      call H5sget_simple_extent_dims_f(coord_dataspace_id, dims_out, coord_dims, error)
      allocate (coord(coord_dims(1), coord_dims(2)))
      call h5dread_f(coord_dataset_id, H5T_NATIVE_DOUBLE, coord, coord_dims, error)
      call h5dclose_f(coord_dataset_id, error)

!****************************************************************************80

      call h5dopen_f(file_id, etopol_dataset, etopol_dataset_id, error)
      call h5dget_space_f(etopol_dataset_id, etopol_dataspace_id, error)
      call H5sget_simple_extent_dims_f(etopol_dataspace_id, dims_out, etopol_dims, error)
      allocate (etopol(etopol_dims(1), etopol_dims(2)))
      call h5dread_f(etopol_dataset_id, H5T_NATIVE_integer, etopol, etopol_dims, error)
      call h5dclose_f(etopol_dataset_id, error)

!****************************************************************************80

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

!****************************************************************************80

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

!****************************************************************************80

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      if (allocated(mesh%coord)) deallocate (mesh%coord, mesh%etopol, mesh%param, mesh%params)
      allocate (mesh%coord(3, size(coord, 2)))
      mesh%coord = coord
      print *, '   Number of nodes                =', size(coord, 2)
      mesh%N_node = size(coord, 2)

      allocate (mesh%etopol(4, size(etopol, 2)))
      mesh%etopol = etopol
      print *, '   Number of elements             =', size(etopol, 2)
      mesh%N_tet = size(etopol, 2)

      allocate (mesh%params(size(param_r, 1), size(param_r, 2)))
      allocate (mesh%param(size(param_r, 1)))
      mesh%params = dcmplx(param_r, param_i)
      write (*, '(2(A,F5.3))') '    Refractive index of medium     =   ', matrices%ref_med
      if (matrices%refr > 1d-7) then
         mesh%param = dcmplx(matrices%refr**2 - matrices%refi**2, &
                             2d0*matrices%refr*matrices%refi)
         write (*, '(2(A,F5.3))') '    Refractive index               =   '&
         , matrices%refr, ' + i', matrices%refi
         write (*, '(2(A,F5.3))') '    Dielectric constant            =   '&
         , matrices%refr**2 - matrices%refi**2, &
            ' + i', 2d0*matrices%refr*matrices%refi
      else
         mesh%param = dcmplx(param_r(:, 1), param_i(:, 1))
         if (maxval(param_r(:, 1)) - minval(param_r(:, 1)) > 1d-7 &
            .OR. maxval(param_i(:, 1)) - minval(param_i(:, 1)) > 1d-7) then
         else
            write (*, '(2(A,F5.3))') '    Dielectric constant  =   '&
            , param_r(1, 1), ' + i', param_i(1, 1)
         end if
      end if

      vol = get_tetra_vol()
      a_eff = (3d0*vol/4d0/pi)**(1d0/3d0)

      mesh%coord = mesh%coord*mesh%a/a_eff ! Scale coordinates to correspond the real grain
! With correction to scaling (effective radius is defined by volumes)

   end subroutine read_mesh

!****************************************************************************80

   subroutine read_T()
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

      integer(HSIZE_T), dimension(1) :: dims_out, dims
      integer(HSIZE_T), dimension(1) :: dimswl, dimsinfo, dims1

      integer :: error, i, nm, ind

      real(dp), dimension(3) :: ref_a
      real(dp), dimension(:), allocatable :: wls, param_r, param_i
      real(dp), dimension(:), allocatable :: Taai_r, Taai_i
      real(dp), dimension(:), allocatable :: Tabi_r, Tabi_i
      real(dp), dimension(:), allocatable :: Tbai_r, Tbai_i
      real(dp), dimension(:), allocatable :: Tbbi_r, Tbbi_i
      complex(dp), dimension(:), allocatable :: Taai, Tabi, Tbai, Tbbi
      complex(dp) :: ref
      LOGICAL :: exists

      file = matrices%tname

      allocate(Taai_i(T_size), Taai_r(T_size), Tabi_r(T_size), Tabi_i(T_size),&
         Tbai_r(T_size), Tbai_i(T_size), Tbbi_r(T_size), Tbbi_i(T_size), &
         Taai(T_size), Tabi(T_size), Tbai(T_size), Tbbi(T_size))

      call h5open_f(error)
      call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error)

!****************************************************************************80

      call h5dopen_f(file_id, dataset1, dataset1_id, error)
      call h5dget_space_f(dataset1_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, Taai_r, dims, error)
      call h5dclose_f(dataset1_id, error)

      call h5dopen_f(file_id, dataset2, dataset2_id, error)
      call h5dget_space_f(dataset2_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      call h5dread_f(dataset2_id, H5T_NATIVE_DOUBLE, Taai_i, dims, error)
      call h5dclose_f(dataset2_id, error)

      Taai = dcmplx(Taai_r, Taai_i)

!****************************************************************************80

      call h5dopen_f(file_id, dataset3, dataset3_id, error)
      call h5dget_space_f(dataset3_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      call h5dread_f(dataset3_id, H5T_NATIVE_DOUBLE, Tabi_r, dims, error)
      call h5dclose_f(dataset3_id, error)

      call h5dopen_f(file_id, dataset4, dataset4_id, error)
      call h5dget_space_f(dataset4_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      call h5dread_f(dataset4_id, H5T_NATIVE_DOUBLE, Tabi_i, dims, error)
      call h5dclose_f(dataset4_id, error)

      Tabi = dcmplx(Tabi_r, Tabi_i)

!****************************************************************************80

      call h5dopen_f(file_id, dataset5, dataset5_id, error)
      call h5dget_space_f(dataset5_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      call h5dread_f(dataset5_id, H5T_NATIVE_DOUBLE, Tbai_r, dims, error)
      call h5dclose_f(dataset5_id, error)

      call h5dopen_f(file_id, dataset6, dataset6_id, error)
      call h5dget_space_f(dataset6_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      call h5dread_f(dataset6_id, H5T_NATIVE_DOUBLE, Tbai_i, dims, error)
      call h5dclose_f(dataset6_id, error)

      Tbai = dcmplx(Tbai_r, Tbai_i)

!****************************************************************************80

      call h5dopen_f(file_id, dataset7, dataset7_id, error)
      call h5dget_space_f(dataset7_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      call h5dread_f(dataset7_id, H5T_NATIVE_DOUBLE, Tbbi_r, dims, error)
      call h5dclose_f(dataset7_id, error)

      call h5dopen_f(file_id, dataset8, dataset8_id, error)
      call h5dget_space_f(dataset8_id, dataspace_id, error)
      call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
      call h5dread_f(dataset8_id, H5T_NATIVE_DOUBLE, Tbbi_i, dims, error)
      call h5dclose_f(dataset8_id, error)

      Tbbi = dcmplx(Tbbi_r, Tbbi_i)

!****************************************************************************80

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

!****************************************************************************80

         call h5dopen_f(file_id, dataset10, dataset10_id, error)
         call h5dget_space_f(dataset10_id, dataspace_id, error)
         call H5sget_simple_extent_dims_f(dataspace_id, dimswl, dims1, error)
         allocate (wls(dims1(1)))
         call h5dread_f(dataset10_id, H5T_NATIVE_DOUBLE, wls, dims1, error)
         call h5dclose_f(dataset10_id, error)

         mesh%ki = 2d0*pi/wls
         do i = 1, size(wls, 1)
            ! num = 0
            ! do j = 1,size(Taai,1)
            !    if(real(Taai(j,1,i))/=0d0) num = num + 1
            ! end do
            ! matrices%Nmaxs(i) = int(dsqrt(real(num)+1d0)-1d0)
            matrices%Nmaxs(i) = truncation_order(mesh%ki(i)*&
               (dble(maxval([mesh%Nx, mesh%Ny, mesh%Nz]))* &
                             mesh%delta)/2.0d0)
         end do

!****************************************************************************80
      end if

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      matrices%Taai = dcmplx(0d0)
      matrices%Tabi = dcmplx(0d0)
      matrices%Tbai = dcmplx(0d0)
      matrices%Tbbi = dcmplx(0d0)

      ind = 1
      do i = 1, matrices%bars
         nm = (matrices%Nmaxs(i)+1)**2-1
         matrices%Taai(1:nm,1:nm,i) = reshape(Taai(ind:(ind-1)+nm**2), [nm,nm])
         matrices%Tbai(1:nm,1:nm,i) = reshape(Tbai(ind:(ind-1)+nm**2), [nm,nm])
         matrices%Tabi(1:nm,1:nm,i) = reshape(Tabi(ind:(ind-1)+nm**2), [nm,nm])
         matrices%Tbbi(1:nm,1:nm,i) = reshape(Tbbi(ind:(ind-1)+nm**2), [nm,nm])
         ind = ind + nm**2
      end do

   end subroutine read_T

! HDF5 WRITE ******************************************************************
!****************************************************************************80
   subroutine write_T()
      if (matrices%singleT == 1) then
         call singleT_write2file()
      else
         call T_write2file()
      end if

   end subroutine write_T

!****************************************************************************80

   subroutine T_write2file()
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

      integer(HSIZE_T), dimension(1) :: dims, dimswl, dimsinfo
      integer     ::    rank = 1
      integer     ::   error ! Error flag
      integer     ::    nm, ind, i

      real(dp), dimension(T_size) :: Taai_r, Taai_i, Tabi_r, Tabi_i, &
                                     Tbai_r, Tbai_i, Tbbi_r, Tbbi_i

      fname = matrices%tname

      filename = fname
      dims = int8((/T_size/))
      dimswl = int8((/matrices%bars/))
      dimsinfo = int8((/3/))

      ind = 1
      do i = 1, matrices%bars
         nm = (matrices%Nmaxs(i) + 1)**2 - 1
         Taai_r(ind:(ind-1)+nm**2) = reshape(real(matrices%Taai(1:nm,1:nm,i)),[nm**2])
         Tbai_r(ind:(ind-1)+nm**2) = reshape(real(matrices%Tbai(1:nm,1:nm,i)),[nm**2])
         Tabi_r(ind:(ind-1)+nm**2) = reshape(real(matrices%Tabi(1:nm,1:nm,i)),[nm**2])
         Tbbi_r(ind:(ind-1)+nm**2) = reshape(real(matrices%Tbbi(1:nm,1:nm,i)),[nm**2])
         Taai_i(ind:(ind-1)+nm**2) = reshape(imag(matrices%Taai(1:nm,1:nm,i)),[nm**2])
         Tbai_i(ind:(ind-1)+nm**2) = reshape(imag(matrices%Tbai(1:nm,1:nm,i)),[nm**2])
         Tabi_i(ind:(ind-1)+nm**2) = reshape(imag(matrices%Tabi(1:nm,1:nm,i)),[nm**2])
         Tbbi_i(ind:(ind-1)+nm**2) = reshape(imag(matrices%Tbbi(1:nm,1:nm,i)),[nm**2])
         ind = ind + nm**2 
      end do

      CALL h5open_f(error)
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      CALL h5screate_simple_f(1, dimsinfo, dspace_id2, error)
      CALL h5screate_simple_f(1, dimswl, dspace_id3, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, Taai_r, dims, error)
      CALL h5dclose_f(dset_id1, error)
      CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, Taai_i, dims, error)
      CALL h5dclose_f(dset_id2, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id3, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, Tabi_r, dims, error)
      CALL h5dclose_f(dset_id3, error)
      CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id4, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, Tabi_i, dims, error)
      CALL h5dclose_f(dset_id4, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id5, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, Tbai_r, dims, error)
      CALL h5dclose_f(dset_id5, error)
      CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id6, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, Tbai_i, dims, error)
      CALL h5dclose_f(dset_id6, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname7, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id7, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, Tbbi_r, dims, error)
      CALL h5dclose_f(dset_id7, error)
      CALL h5dcreate_f(file_id, dsetname8, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id8, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, Tbbi_i, dims, error)
      CALL h5dclose_f(dset_id8, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname9, H5T_NATIVE_DOUBLE, dspace_id2, &
                       dset_id9, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id9, H5T_NATIVE_DOUBLE, [real(mesh%param(1)), &
                                                    imag(mesh%param(1)), mesh%a], dimsinfo, error)
      CALL h5dclose_f(dset_id9, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname10, H5T_NATIVE_DOUBLE, dspace_id3, &
                       dset_id10, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id10, H5T_NATIVE_DOUBLE, 2d0*pi/mesh%ki, dimswl, error)
      CALL h5dclose_f(dset_id10, error)

!****************************************************************************80

      CALL h5sclose_f(dspace_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   end subroutine T_write2file

!****************************************************************************80

   subroutine T_empty()
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

      integer(HSIZE_T), dimension(1) :: dims, dimswl, dimsinfo
      integer     ::    rank = 1
      integer     ::   error ! Error flag

      complex(dp), dimension(:), allocatable    ::   emptyT

      fname = matrices%tname
      if(file_exists(fname)) return

      filename = fname
      dims = int8((/T_size/))
      dimswl = int8((/matrices%bars/))
      dimsinfo = int8((/3/))
      allocate (emptyT(T_size))
      emptyT = dcmplx(0d0, 0d0)

      CALL h5open_f(error)
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      CALL h5screate_simple_f(1, dimsinfo, dspace_id2, error)
      CALL h5screate_simple_f(1, dimswl, dspace_id3, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id1, error)
      CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id2, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id3, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id3, error)
      CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id4, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id4, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id5, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id5, error)
      CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id6, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id6, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname7, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id7, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id7, error)
      CALL h5dcreate_f(file_id, dsetname8, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id8, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id8, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname9, H5T_NATIVE_DOUBLE, dspace_id2, &
                       dset_id9, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id9, H5T_NATIVE_DOUBLE, [real(mesh%param(1)), &
                                                    imag(mesh%param(1)), mesh%a], dimsinfo, error)
      CALL h5dclose_f(dset_id9, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname10, H5T_NATIVE_DOUBLE, dspace_id3, &
                       dset_id10, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id10, H5T_NATIVE_DOUBLE, 2d0*pi/mesh%ki, dimswl, error)
      CALL h5dclose_f(dset_id10, error)

!****************************************************************************80

      CALL h5sclose_f(dspace_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   end subroutine T_empty

!****************************************************************************80

   subroutine singleT_write2file()
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

      integer(HSIZE_T), dimension(1) :: dims, dims_out
      real(dp), dimension(:), allocatable :: Ti_r, Ti_i
      integer     ::   error ! Error flag
      integer :: i, ind1, ind2, nm

      filename = matrices%tname
      dims = int8((/T_size/))
      allocate(Ti_r(T_size), Ti_i(T_size))

      ind1 = 1
      do i = 1, matrices%bars
         nm = (matrices%Nmaxs(i)+1)**2-1
         if(i==matrices%whichbar)then
            ind2 = ind1+nm**2-1
            exit
         else
            ind1 = ind1+nm**2
         end if
      end do

      call h5open_f(error)
      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

!****************************************************************************80

      call h5dopen_f(file_id, dsetname1, dset_id1, error)
      call h5dget_space_f(dset_id1, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id1, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(ind1:ind2) = reshape(real(matrices%Taai(1:nm, 1:nm, &
         matrices%whichbar)), [nm**2])
      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id1, error)

      call h5dopen_f(file_id, dsetname2, dset_id2, error)
      call h5dget_space_f(dset_id2, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id2, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(ind1:ind2) = reshape(imag(matrices%Taai(1:nm, 1:nm, &
         matrices%whichbar)), [nm**2])
      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id2, error)

!****************************************************************************80

      call h5dopen_f(file_id, dsetname3, dset_id3, error)
      call h5dget_space_f(dset_id3, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id3, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(ind1:ind2) = reshape(real(matrices%Tabi(1:nm, 1:nm, &
         matrices%whichbar)), [nm**2])
      CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id3, error)

      call h5dopen_f(file_id, dsetname4, dset_id4, error)
      call h5dget_space_f(dset_id4, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id4, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(ind1:ind2) = reshape(imag(matrices%Tabi(1:nm, 1:nm, &
         matrices%whichbar)), [nm**2])
      CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id4, error)

!****************************************************************************80

      call h5dopen_f(file_id, dsetname5, dset_id5, error)
      call h5dget_space_f(dset_id5, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id5, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(ind1:ind2) = reshape(real(matrices%Tbai(1:nm, 1:nm, &
         matrices%whichbar)), [nm**2])
      CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id5, error)

      call h5dopen_f(file_id, dsetname6, dset_id6, error)
      call h5dget_space_f(dset_id6, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id6, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(ind1:ind2) = reshape(imag(matrices%Tbai(1:nm, 1:nm, &
         matrices%whichbar)), [nm**2])
      CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id6, error)

!****************************************************************************80

      call h5dopen_f(file_id, dsetname7, dset_id7, error)
      call h5dget_space_f(dset_id7, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id7, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(ind1:ind2) = reshape(real(matrices%Tbbi(1:nm, 1:nm, &
         matrices%whichbar)), [nm**2])
      CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id7, error)

      call h5dopen_f(file_id, dsetname8, dset_id8, error)
      call h5dget_space_f(dset_id8, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id8, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(ind1:ind2) = reshape(imag(matrices%Tbbi(1:nm, 1:nm, &
         matrices%whichbar)), [nm**2])
      CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id8, error)

!****************************************************************************80

      CALL h5sclose_f(dspace_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   end subroutine singleT_write2file

!****************************************************************************80

   subroutine write2file(A, fname)
      complex(dp), intent(in) :: A(:, :)
      character(len=80):: fname
      character(len=80) :: filename
      character(len=3), PARAMETER :: dsetname1 = "A_r" ! Dataset name
      character(len=3), PARAMETER :: dsetname2 = "A_i" ! Dataset name

      integer(HID_T) :: file_id ! File identifier
      integer(HID_T) :: dset_id1 ! Dataset identifier
      integer(HID_T) :: dset_id2 ! Dataset identifier
      integer(HID_T) :: dspace_id ! Dataspace identifier

      integer(HSIZE_T), dimension(2) :: dims ! Dataset dimensions
      integer     ::    rank = 2 ! Dataset rank

      integer     ::   error ! Error flag

      filename = fname
      dims = [int(size(A, 1), 8), int(size(A, 2), 8)]

      CALL h5open_f(error)

! Create a new file using default properties.
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
! Create the dataspace.
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
! Create and write dataset using default properties.
      CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(A), dims, error)
      CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(A), dims, error)
! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id1, error)
      CALL h5dclose_f(dset_id2, error)
! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, error)
! Close the file.
      CALL h5fclose_f(file_id, error)
! Close FORTRAN interface.
      CALL h5close_f(error)

   end subroutine write2file

end module