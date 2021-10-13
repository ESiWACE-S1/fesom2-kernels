program exchange_nod_prog
  use wallclock_mod
  use io_mod
  use mpi
  use distributed_mod
  implicit none

  integer :: n, nu1, nl1, nz, myDim_nod2D, n_it, nod2D
  integer, dimension(:), allocatable :: ulevels_nod2D
  integer, dimension(:), allocatable :: nlevels_nod2D
  integer, parameter :: MAX_LEVELS=50
  integer, parameter :: MAX_ITERATIONS=200
  integer :: nn

  integer :: myDim_elem2D
  integer :: elem, nl, nl12, nu12
  integer, dimension(:), allocatable :: ulevels
  integer, dimension(:), allocatable :: nlevels
  integer, dimension(3) :: enodes
  integer, dimension(:,:), allocatable :: elem2D_nodes
  double precision :: bignumber = huge(1.0)

  integer, dimension(:), allocatable :: nod_in_elem2D_num
  integer, dimension(:,:), allocatable :: nod_in_elem2D

  integer :: nl2, nu2, el(2), edge
  integer :: myDim_edge2D
  integer, dimension(:,:), allocatable :: edges
  integer, dimension(:,:), allocatable :: edge_tri
  integer, dimension(:), allocatable :: ulevels_edge
  integer, dimension(:), allocatable :: nlevels_edge
  double precision, dimension(:,:), allocatable :: fct_plus
  double precision, dimension(:,:), allocatable :: fct_minus

  integer, dimension(:), allocatable :: n_parts
  integer, dimension(:), allocatable :: parts_id
  character(1024) :: file_name
  
  real(8)::t1,delta, delta_orig

  call MPI_Init(MPIerr)
  call MPI_Comm_size(MPI_COMM_WORLD, npes, MPIerr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mype, MPIerr)

  write(file_name, '(i8)') npes
  file_name='dist_'//trim(adjustl(file_name))//'_rpart.out'
  call load_part(file_name, n_parts, parts_id)

  if(mype == 0) write(*,*) parts_id
  nod2D = read_nod2D_elems(0, nod_in_elem2D_num, nod_in_elem2D)

!!  myDim_nod2D = read_nod2D_levels(mype, ulevels_nod2D, nlevels_nod2D)
!!
!!  myDim_elem2D = read_elem2D_nodes(mype, elem2D_nodes)
!!  myDim_elem2D = read_elem2D_levels(mype, ulevels, nlevels)
!!
!!
!!  myDim_edge2D = read_edges_nodes(mype, edges)
!!  myDim_edge2D = read_edges_levels(mype, ulevels_edge, nlevels_edge)
!!  myDim_edge2D = read_edges_tri(mype, edge_tri)
!!
!!  allocate(fct_plus(MAX_LEVELS, myDim_nod2D))
!!  allocate(fct_minus(MAX_LEVELS, myDim_nod2D))
!!
!!  nl=maxval(nlevels(:))
!!
!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  t1=wallclock()
!!  do n_it=1, MAX_ITERATIONS
!!
!!    !!$acc update self(fct_plus,fct_minus)
!!    ! fct_minus and fct_plus must be known to neighbouring PE
!!    call exchange_nod3D(fct_plus, fct_minus)
!!
!!    !call exchange_nod3D_begin(fct_plus, fct_minus)
!!    !call exchange_nod3D_end  ! fct_plus, fct_minus
!!    !!$acc update device(fct_plus,fct_minus)
!!  end do
!!  delta=wallclock()-t1
!!  delta_orig=delta
!!  write(*,*) "done"
!!  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
!!
  call MPI_Finalize(MPIerr)


end program exchange_nod_prog
