program bench_naive_io
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
  u(:,:,:) = 1._rp*myid
  v(:,:,:) = 2._rp*myid
  w(:,:,:) = 3._rp*myid
  p(:,:,:) = 4._rp*myid
  istep    = 420
  time     = 42._rp
  !
  ! write one file per task and time it
  !
  inquire(iolength=rlen) 1._rp
  isize = product(int(nn(:),8))*rlen
  write(cmyid,'(i5.5)') myid
  tw = MPI_WTIME()
  open(newunit=iunit,file=file_o//'-'//cmyid,access='direct',recl=isize)
    write(iunit,rec=1) u(1:nn(1),1:nn(2),1:nn(3))
    write(iunit,rec=2) v(1:nn(1),1:nn(2),1:nn(3))
    write(iunit,rec=3) w(1:nn(1),1:nn(2),1:nn(3))
    write(iunit,rec=4) p(1:nn(1),1:nn(2),1:nn(3))
  close(iunit)
  tw = MPI_WTIME()-tw
  call MPI_ALLREDUCE(tw,twmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(tw,twmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(tw,twav ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(myid == 0) print*, 'WRITING (one file per task): Avrg, min & max elapsed time: '
  if(myid == 0) print*, twav/(1._rp*nproc),twmin,twmax
  if(myid == 0) print*,''
  !
  ! read one file per task and time it
  !
  tw = MPI_WTIME()
  open(newunit=iunit,file=file_i//'-'//cmyid,access='direct',recl=isize)
    read(iunit,rec=1) u(1:nn(1),1:nn(2),1:nn(3))
    read(iunit,rec=2) v(1:nn(1),1:nn(2),1:nn(3))
    read(iunit,rec=3) w(1:nn(1),1:nn(2),1:nn(3))
    read(iunit,rec=4) p(1:nn(1),1:nn(2),1:nn(3))
  close(iunit)
  tw = MPI_WTIME()-tw
  call MPI_ALLREDUCE(tw,twmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(tw,twmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(tw,twav ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(myid == 0) print*, 'READING (one file per task): Avrg, min & max elapsed time: '
  if(myid == 0) print*, twav/(1._rp*nproc),twmin,twmax
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
end program bench_naive_io
