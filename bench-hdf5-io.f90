program bench_hdf5_io
  use hdf5
  use mpi_f08
  use mpi_f08, only: MPI_REAL_RP => MPI_DOUBLE_PRECISION
  use, intrinsic :: iso_fortran_env, only: rp => real64
  implicit none
  !
  ! input domain parameters
  !
  integer , parameter,       dimension(3) :: n    = [256,256,256]
  !
  ! file names
  !
  character(len=*), parameter             :: file_i = 'data/fld.bin', &
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
  call MPI_DIMS_CREATE(nproc,3,dims(:),ierr)
  !
  call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,[.false.,.false.,.false.],.true.,comm_cart,ierr)
  call MPI_CART_COORDS(comm_cart,myid,3,coords,ierr)
  !
  ! decompose the domain
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
  u(:,:,:) = 1._rp*myid
  v(:,:,:) = 2._rp*myid
  w(:,:,:) = 3._rp*myid
  p(:,:,:) = 4._rp*myid
  istep    = 420
  time     = 42._rp
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
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: u,v,w,p
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    real(rp), dimension(2) :: fldinfo
    type(MPI_FILE) :: fh
    integer :: nreals_myid
    !
    select case(io)
    case('r')
      call io_field('r',"u.h5",ng,nh,lo,hi,u)
      call io_field('r',"v.h5",ng,nh,lo,hi,v)
      call io_field('r',"w.h5",ng,nh,lo,hi,w)
      call io_field('r',"p.h5",ng,nh,lo,hi,p)
    case('w')
      call io_field('w',"u.h5",ng,nh,lo,hi,u)
      call io_field('w',"v.h5",ng,nh,lo,hi,v)
      call io_field('w',"w.h5",ng,nh,lo,hi,w)
      call io_field('w',"p.h5",ng,nh,lo,hi,p)
    end select
  end subroutine load
  !
  subroutine io_field(io,filename,ng,nh,lo,hi,var)
    implicit none
    character(len=1), intent(in)                 :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in), dimension(3)   :: ng,nh,lo,hi
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    integer , dimension(3) :: n
    integer , dimension(3) :: sizes,subsizes,starts
    !
    ! HDF5 variables
    !
    integer :: ndims, ierr
    integer(HID_T) :: file_id
    integer(HID_T) :: filespace
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    !
    integer(HID_T) :: dset
    !
    integer(HSIZE_T) :: dims(3)
    !
    integer(HID_T) :: plist_id
    integer(HSIZE_T), dimension(3) :: data_count  
    integer(HSSIZE_T), dimension(3) :: data_offset 
    integer(HSSIZE_T), dimension(3) :: halo_offset 
    type(MPI_INFO) :: info = MPI_INFO_NULL
    !
    n(:)        = hi(:)-lo(:)+1
    sizes(:)    = ng(:)
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    !
    ndims = 3
    dims(:) = ng(:)
    data_count(:) = subsizes(:)
    data_offset(:) = starts(:)
    halo_offset(:) = nh(:)
    !
    select case(io)
    case('r')
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD%MPI_VAL,info%MPI_VAL,ierr)
      call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,ierr,access_prp=plist_id)
      call h5pclose_f(plist_id,ierr)
      call h5dopen_f(file_id,'var',dset,ierr)
      call h5screate_simple_f(ndims,data_count+2*nh(:),memspace,ierr) 
      call h5dget_space_f(dset,slabspace,ierr)
      call h5sselect_hyperslab_f (slabspace,H5S_SELECT_SET_F,data_offset,data_count,ierr)
      call h5sselect_hyperslab_f (memspace,H5S_SELECT_SET_F,halo_offset,data_count,ierr)
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr) 
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      call h5dread_f(dset,H5T_NATIVE_DOUBLE,var,dims,ierr,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
      call h5pclose_f(plist_id,ierr)
      call h5dclose_f(dset,ierr)
      call h5sclose_f(memspace,ierr)
      call h5fclose_f(file_id,ierr)
    case('w')
      call h5screate_simple_f(ndims,dims,filespace,ierr)
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD%MPI_VAL,info%MPI_VAL,ierr)
      call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr,access_prp=plist_id)
      call h5pclose_f(plist_id,ierr)
      call h5dcreate_f(file_id,'var',H5T_NATIVE_DOUBLE,filespace,dset,ierr)
      call h5screate_simple_f(ndims,data_count+2*nh(:),memspace,ierr) 
      call h5dget_space_f(dset,slabspace,ierr)
      call h5sselect_hyperslab_f (slabspace,H5S_SELECT_SET_F,data_offset,data_count,ierr)
      call h5sselect_hyperslab_f (memspace,H5S_SELECT_SET_F,halo_offset,data_count,ierr)
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr) 
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,var,dims,ierr,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
      call h5pclose_f(plist_id,ierr)
      call h5dclose_f(dset,ierr)
      call h5sclose_f(memspace,ierr)
      call h5fclose_f(file_id,ierr)
    end select
  end subroutine io_field
end program bench_hdf5_io
