program bench_mpi_io
  use mpi_f08
  use mpi_f08, only: MPI_REAL_RP => MPI_DOUBLE_PRECISION
  use, intrinsic :: iso_fortran_env, only: rp => real64
  implicit none
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
  integer                                 :: myid,nproc,dims(3),coords(3),ierr
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
  ! initialize MPI
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
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
  call random_number(u)
  call random_number(v)
  call random_number(w)
  call random_number(p)
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
  if(myid == 0) print*, 'WRITING (MPI I/O): Avrg, min & max elapsed time: '
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
  if(myid == 0) print*, 'READING (MPI I/O): Avrg, min & max elapsed time: '
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
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    integer :: ierr
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(comm,filename, &
           MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      good = (product(int(ng(:),MPI_OFFSET_KIND))*4+2)*int(storage_size(1._rp)/8)
      if(filesize /= good) then
        if(myid == 0) print*, ''
        if(myid == 0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid == 0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call MPI_FINALIZE()
        error stop
      end if
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call io_field('r',fh,ng,nh,lo,hi,disp,u)
      call io_field('r',fh,ng,nh,lo,hi,disp,v)
      call io_field('r',fh,ng,nh,lo,hi,disp,w)
      call io_field('r',fh,ng,nh,lo,hi,disp,p)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
      nreals_myid = 0
      if(myid == 0) nreals_myid = 2
      call MPI_FILE_READ(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
      call MPI_BCAST(fldinfo,2,MPI_REAL_RP,0,comm,ierr)
      time  =      fldinfo(1)
      istep = nint(fldinfo(2))
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(comm,filename, &
           MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)
      disp = 0_MPI_OFFSET_KIND
      call io_field('w',fh,ng,nh,lo,hi,disp,u)
      call io_field('w',fh,ng,nh,lo,hi,disp,v)
      call io_field('w',fh,ng,nh,lo,hi,disp,w)
      call io_field('w',fh,ng,nh,lo,hi,disp,p)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
      fldinfo = [time,1._rp*istep]
      nreals_myid = 0
      if(myid == 0) nreals_myid = 2
      call MPI_FILE_WRITE(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
    end select
  end subroutine load
  !
  subroutine io_field(io,fh,ng,nh,lo,hi,disp,var)
    implicit none
    character(len=1), intent(in)                 :: io
    type(MPI_FILE)  , intent(in)                 :: fh
    integer         , intent(in), dimension(3)   :: ng,nh,lo,hi
    integer(kind=MPI_OFFSET_KIND), intent(inout) :: disp
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    integer , dimension(3) :: n
    integer , dimension(3) :: sizes,subsizes,starts
    type(MPI_DATATYPE) :: type_glob,type_loc
    integer :: ierr
    !
    n(:)        = hi(:)-lo(:)+1
    sizes(:)    = ng(:)
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_glob,ierr)
    call MPI_TYPE_COMMIT(type_glob,ierr)
    sizes(:)    = n(:) + 2*nh(:)
    subsizes(:) = n(:)
    starts(:)   = 0 + nh(:)
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_loc ,ierr)
    call MPI_TYPE_COMMIT(type_loc,ierr)
    select case(io)
    case('r')
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_READ_ALL(fh,var,1,type_loc,MPI_STATUS_IGNORE,ierr)
    case('w')
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_WRITE_ALL(fh,var,1,type_loc,MPI_STATUS_IGNORE,ierr)
    end select
    disp = disp+product(int(ng(:),MPI_OFFSET_KIND))*(storage_size(1._rp)/8)
    call MPI_TYPE_FREE(type_glob,ierr)
    call MPI_TYPE_FREE(type_loc ,ierr)
  end subroutine io_field
end program bench_mpi_io
