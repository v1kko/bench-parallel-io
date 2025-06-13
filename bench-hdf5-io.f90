program bench_hdf5_io
  use hdf5
  use mpi_f08
  use mpi_f08, only: MPI_REAL_RP => MPI_DOUBLE_PRECISION
  use, intrinsic :: iso_fortran_env, only: rp => real64
  implicit none
  type :: arr_ptr
    real(rp), pointer, contiguous , dimension(:,:,:) :: arr
  end type arr_ptr
  !
  ! input domain parameters
  !
  integer , parameter,       dimension(3) :: n    = [1024,1024,1024]
  !
  ! file names
  !
  character(len=*), parameter             :: file_i = '/scratch-shared/vazizi/fld.bin', &
                                             file_o = file_i
  !
  ! local problem sizes
  !
  integer, dimension(3)                   :: nn,lo,hi
  !
  ! MPI stuff
  !
  integer                                 :: myid,nproc,dims(3),coords(3)
  type(MPI_COMM)                          :: comm_cart
  !
  ! computational variables
  !
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p
  real(rp)                                :: time
  integer                                 :: istep
  !
  ! for parallel I/O w/ one file per task
  !
  integer    :: iunit,rlen
  integer(8) :: isize
  character(len=5) :: cmyid
  !
  ! variables for timing
  !
  real(rp) :: tw,twav,twmin,twmax
  ! 
  integer :: ierr
  !
  ! initialize MPI
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  !
  call h5open_f(ierr)
  !
  ! create processor grid
  !
  dims(:) = 0
  coords(:) = 0
  call MPI_DIMS_CREATE(nproc,2,dims(2:3),ierr)
  !
  call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims(2:3),[.false.,.false.],.true.,comm_cart,ierr)
  call MPI_CART_COORDS(comm_cart,myid,2,coords(2:3),ierr)
  dims(1) = 1
  !
  ! decompose the domain, 2D, pencil in firxt axis
  !
  call distribute_grid(n,dims,coords,[1,1,1],nn,lo,hi)

  !
  ! allocate input and output arrays
  !
  allocate(u(0:nn(1)+1,0:nn(2)+1,0:nn(3)+1), &
           v(0:nn(1)+1,0:nn(2)+1,0:nn(3)+1), &
           w(0:nn(1)+1,0:nn(2)+1,0:nn(3)+1), &
           p(0:nn(1)+1,0:nn(2)+1,0:nn(3)+1))
  !
  ! create some random data
  !
  call random_number(u)
  call random_number(v)
  call random_number(w)
  call random_number(p)
  istep    = 420
  time     = 42._rp

  ! remove existing checkpoint
  call unlink(file_i)

  !
  ! write w/ MPI I/O an time it
  !
  tw = MPI_WTIME()
  call load('w',file_o,MPI_COMM_WORLD,myid,n,[1,1,1],lo,hi,u,v,w,p,time,istep)
  tw = MPI_WTIME()-tw
  call MPI_ALLREDUCE(tw,twmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(tw,twmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(tw,twav ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(myid == 0) print*, 'WRITING (HDF5 I/O): Avrg, min & max elapsed time: '
  if(myid == 0) print*, twav/(1._rp*nproc),twmin,twmax
  if(myid == 0) print*,''
  !
  ! read w/ MPI I/O and time it
  !
  tw = MPI_WTIME()
  call load('r',file_i,MPI_COMM_WORLD,myid,n,[1,1,1],lo,hi,u,v,w,p,time,istep)
  tw = MPI_WTIME()-tw
  if(myid == 0) print*, 'Loaded field at time = ', time, 'step = ',istep,'.'
  call MPI_ALLREDUCE(tw,twmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(tw,twmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(tw,twav ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(myid == 0) print*, 'READING (HDF5 I/O): Avrg, min & max elapsed time: '
  if(myid == 0) print*, twav/(1._rp*nproc),twmin,twmax
  if(myid == 0) print*,''
  !
  call MPI_FINALIZE(ierr)
contains
  subroutine distribute_grid(ng,dims,coords,lo_g,n,lo,hi)
    implicit none
    integer, intent(in ), dimension(3) :: ng,dims,coords,lo_g
    integer, intent(out), dimension(3) :: n,lo,hi
    n(:) = ng(:)/dims(:)
    where(coords(:)+1 <= mod(ng(:),dims(:))) n(:) = n(:) + 1
    lo(:) = lo_g(:)   + (coords(:)  )*n(:)
    hi(:) = lo_g(:)-1 + (coords(:)+1)*n(:)
    where(coords(:)+1 >  mod(ng(:),dims(:)))
      lo(:) = lo(:) +    mod(ng(:),dims(:))
      hi(:) = hi(:) +    mod(ng(:),dims(:))
    end where
  end subroutine distribute_grid
  !
  subroutine load(io,filename,comm,myid,ng,nh,lo,hi,u,v,w,p,time,istep)
    !
    ! reads/writes a restart file
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    type(MPI_COMM)  , intent(in) :: comm
    integer         , intent(in) :: myid
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    real(rp), target, intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: u,v,w,p
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    type(arr_ptr) ::   io_vars(4)
    character(len=10) :: c_io_vars(4)
    real(rp), dimension(2) :: fldinfo
    type(MPI_FILE) :: fh
    integer :: nreals_myid
    
    io_vars(1)%arr => u
    io_vars(2)%arr => v
    io_vars(3)%arr => w
    io_vars(4)%arr => p
    c_io_vars(1) = "u"
    c_io_vars(2) = "v"
    c_io_vars(3) = "w"
    c_io_vars(4) = "p"

    !
    select case(io)
    case('r')
      call io_field_hdf5('r',filename,c_io_vars,MPI_COMM_WORLD,ng,nh,lo,hi,io_vars,time,istep,4)
    case('w')
      call io_field_hdf5('w',filename,c_io_vars,MPI_COMM_WORLD,ng,nh,lo,hi,io_vars,time,istep,4)
    end select
  end subroutine load
  !
  subroutine io_field_hdf5(io,filename,c_io_vars,comm,ng,nh,lo,hi,io_vars,time,istep, nvar)
    use iso_fortran_env
    use hdf5
    !
    ! collective single field data I/O using HDF5
    !
    ! written with the help of Josh Romero,
    ! with the field data I/O inspired from the AFiD code
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    type(mpi_comm)  , intent(in) :: comm
    integer         , intent(in) :: nvar
    integer         , intent(in), dimension(3)   :: ng,nh,lo,hi
    ! Must statically define kind real64 for hdf5 routine,
    ! either duplicate for other precisions (and add to interface) or do some Macro wizardry
    !real(real64), intent(inout),contiguous, dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    type(arr_ptr)    ,intent(in) ::   io_vars(:)
    character(len=10),intent(in) :: c_io_vars(:)
    real(real64), intent(inout) :: time
    integer(int32) , intent(inout) :: istep

    ! adjustable parameters
    logical :: chunk_checkpoint = .true.
    integer :: compression_level = 0
    integer :: chunk_size = 16
    integer :: ipencil_axis = 1

    real(real64), pointer, dimension(:,:,:) :: var
    integer , dimension(3) :: n
    integer , dimension(3) :: subsizes, starts
    !
    ! HDF5 variables
    !
    integer :: ndims, ierr, i
    integer(HID_T) :: file_id, group_id
    integer(HID_T) :: dspace_id, dset_id, attr_id
    integer(HID_T) :: filespace
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    !
    integer(HSIZE_T) :: dims(3), chunk(3)
    !
    integer(HID_T) :: xfer_pid, file_pid, dset_pid
    integer(HSIZE_T) , dimension(3) :: data_count
    integer(HSSIZE_T), dimension(3) :: data_offset
    integer(HSSIZE_T), dimension(3) :: halo_offset
    logical :: file_exists, group_exists, dset_exists
    character(len=20):: name, varname
    integer(hsize_t) :: data_size

    !
    n(:)        = hi(:)-lo(:)+1
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    !
    ndims = 3
    dims(:) = ng(:)
    data_count(:) = subsizes(:)
    data_offset(:) = starts(:)
    halo_offset(:) = nh(:)+1
    !

    ! Common operations
    call h5open_f(ierr)


    call h5pcreate_f(h5p_dataset_xfer_f, xfer_pid, ierr)
    call h5pset_dxpl_mpio_f(xfer_pid, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5pcreate_f(h5p_file_access_f, file_pid, ierr)
    call h5pset_fapl_mpio_f(file_pid, comm, MPI_INFO_NULL, ierr)
    call h5pcreate_f(h5p_dataset_create_f, dset_pid, ierr)

    if(chunk_checkpoint) then
      ! Sets the other dimensions of the chunk, 1 make a chunk of every independent pencil
      chunk = chunk_size
      !Change chunking axis by editing the following line
      chunk(ipencil_axis) = ng(ipencil_axis)
      !Turn chunks on/off by commenting out the following line
      call h5pset_chunk_f(dset_pid, ndims, chunk , ierr)

      !Turn compression on/off by toggeling the following line, or change the level (1 least/fast, 9 most/slow)
      call h5pset_deflate_f(dset_pid, compression_level, ierr) 
    endif

    call h5screate_simple_f(ndims, data_count+2*nh(:), memspace, ierr) 
    call h5screate_simple_f(ndims, dims, slabspace, ierr)

    call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, ierr)
    call h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, halo_offset, data_count, ierr)

    select case(io)
    case('r')
      inquire(file=filename,exist=file_exists)
      if (file_exists) then 
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr, access_prp=file_pid)
      else
        error stop "Checkpoint file "//filename//"  does not exist"
      endif

      ! find latest group
      call h5lget_name_by_idx_f(file_id, ".", H5_INDEX_NAME_F, H5_ITER_DEC_F, 0_HSIZE_T, name, ierr, data_size)

      call h5aopen_by_name_f(file_id, trim(name), 'time', attr_id, ierr)
      call h5aread_f(attr_id, H5T_IEEE_F64LE, time, (/1_HSIZE_T/), ierr)
      call h5aclose_f(attr_id, ierr)

      call h5aopen_by_name_f(file_id, trim(name), 'istep', attr_id, ierr)
      call h5aread_f(attr_id, H5T_STD_I32LE, istep, (/1_HSIZE_T/), ierr)
      call h5aclose_f(attr_id, ierr)

      do i = 1, nvar
        varname = c_io_vars(i)
        var(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) => io_vars(i)%arr
        call h5dopen_f(file_id, trim(name)//'/'//varname, dset_id, ierr)
        call h5dread_f(dset_id,H5T_IEEE_F64LE,var,dims,ierr, &
                       file_space_id=slabspace,mem_space_id=memspace,xfer_prp=xfer_pid)
        call h5dclose_f(dset_id, ierr)
      end do

      call h5fclose_f(file_id, ierr)

    case('w')
      inquire(file=filename,exist=file_exists)
      if (file_exists) then 
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, access_prp=file_pid)
      else
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, access_prp=file_pid)
      endif

      write(name,'(I8.8)') istep
      name = 'istep_'//trim(adjustl(name))
      call h5lexists_f(file_id, name, group_exists, ierr)
      if (group_exists) then
        call h5gopen_f(file_id, name, group_id, ierr)
      else
        call h5gcreate_f(file_id, name, group_id, ierr)

        call h5screate_simple_f(1, (/1_HSIZE_T/), dspace_id, ierr)
        call h5acreate_f(group_id, "time", H5T_IEEE_F64LE, dspace_id, attr_id, ierr)
        call h5awrite_f(attr_id, H5T_IEEE_F64LE, time, (/1_HSIZE_T/), ierr)
        call h5aclose_f(attr_id, ierr)

        call h5acreate_f(group_id, "istep", H5T_STD_I32LE, dspace_id, attr_id, ierr)
        call h5awrite_f(attr_id, H5T_STD_I32LE, istep, (/1_HSIZE_T/), ierr)
        call h5aclose_f(attr_id, ierr)
        call h5sclose_f(dspace_id, ierr)
      end if

      call h5screate_simple_f(ndims, dims, dspace_id, ierr)

      do i = 1, nvar
        varname = c_io_vars(i)
        call h5lexists_f(group_id, varname, dset_exists, ierr)
        if (dset_exists) then
          call h5dopen_f(group_id, varname, dset_id, ierr)
        else
          call h5dcreate_f(group_id, varname, H5T_IEEE_F64LE, dspace_id, dset_id, ierr, &
                           dcpl_id=dset_pid)
        end if

        var(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) => io_vars(i)%arr
        call h5dwrite_f(dset_id, H5T_IEEE_F64LE, var, dims, ierr, &
                        file_space_id=slabspace, mem_space_id=memspace, xfer_prp=xfer_pid)
        call h5dclose_f(dset_id, ierr)
      end do
      call h5sclose_f(dspace_id, ierr)
      call h5gclose_f(group_id, ierr)
      call h5fclose_f(file_id, ierr)
    end select
10 continue
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(slabspace, ierr)
    call h5pclose_f(file_pid, ierr)
    call h5pclose_f(xfer_pid, ierr)
    call h5close_f(ierr)
  end subroutine io_field_hdf5
end program bench_hdf5_io
