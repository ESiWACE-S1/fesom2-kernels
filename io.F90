module io_mod
  implicit none
  private

  integer, parameter :: MAX_LEVELS=50
  integer, parameter :: MAX_PATH=1024
  integer, parameter :: MAX_NOD_IN_ELEM=10

  public :: get_myDim_nod2D
  public :: get_myDim_elem2D
  
  public :: read_nod2D_levels
  public :: read_nod2D_elems

  public :: read_elem2D_nodes
  public :: read_elem2D_levels

  public :: read_edges_nodes
  public :: read_edges_levels
  public :: read_edges_tri

contains
  function get_myDim_nod2D(mype) result(myDim_nod2D)
    implicit none
    integer, intent(in) :: mype
    integer :: myDim_nod2D
    character(MAX_PATH) :: file_name
    integer :: fileID = 10

    write(file_name, '(i8)') mype
    file_name='loops_nod2D_levels_'//trim(adjustl(file_name))//'.dat'
    open(fileID, file=file_name, status='old', action='read')
    read(fileID, '(i8)') myDim_nod2D
    close(fileID)
  end function get_myDim_nod2D

  function get_myDim_elem2D(mype) result(myDim_elem2D)
    implicit none
    integer, intent(in) :: mype
    integer :: myDim_elem2D
    character(MAX_PATH) :: file_name
    integer :: fileID = 10

    write(file_name, '(i8)') mype
    file_name='loops_elem2D_levels_'//trim(adjustl(file_name))//'.dat'
    open(fileID, file=file_name, status='old', action='read')
    read(fileID, '(i8)') myDim_elem2D
    close(fileID)
  end function get_myDim_elem2D

 
  function read_nod2D_levels(mype, ulevels_nod2D, nlevels_nod2D) result(myDim_nod2D)
    implicit none
    integer, intent(in) :: mype
    integer, intent(inout), dimension(:), allocatable :: ulevels_nod2D
    integer, intent(inout), dimension(:), allocatable :: nlevels_nod2D
    integer :: myDim_nod2D
    character(MAX_PATH) :: file_name
    integer :: fileID = 10, n, nn, nz

    myDim_nod2D = 0

    write(file_name, '(i8)') mype
    file_name='loops_nod2D_levels_'//trim(adjustl(file_name))//'.dat'
    open(fileID, file=file_name, status='old', action='read')

    read(fileID, *) myDim_nod2D

    print*, "loading",myDim_nod2D, " node levels"

    allocate(ulevels_nod2D(myDim_nod2D))
    allocate(nlevels_nod2D(myDim_nod2D))

    do n=1,myDim_nod2D !+ edim_nod2d
       ! TODO check nz<=MAX_LEVELS
       read(fileID, *) nn, nz, ulevels_nod2D(n), nlevels_nod2D(n)
       nlevels_nod2D(n) = nlevels_nod2D(n)+1
    end do

    close(fileID)

  end function read_nod2D_levels

  function read_nod2D_elems(mype, nod_in_elem2D_num, nod_in_elem2D) result(myDim_nod2D)
    implicit none
    integer, intent(in) :: mype
    integer, intent(inout), dimension(:), allocatable :: nod_in_elem2D_num
    integer, intent(inout), dimension(:,:), allocatable :: nod_in_elem2D
    integer :: myDim_nod2D
    character(MAX_PATH) :: file_name, line
    integer :: fileID = 11, n, nn

    myDim_nod2D = 0
    write(file_name, '(i8)') mype
    file_name='loops_nod2D_elems_'//trim(adjustl(file_name))//'.dat'
    open(fileID, file=file_name, status='old', action='read')
    read(fileID, '(i8)') myDim_nod2D

    print*, "loading",myDim_nod2D, " node elems"


    allocate(nod_in_elem2D_num(myDim_nod2D))
    allocate(nod_in_elem2D(MAX_NOD_IN_ELEM, myDim_nod2D))

    do n=1,myDim_nod2D !+ edim_nod2d
       !read(fileID,*) nn, nod_in_elem2D_num(n), nod_in_elem2D(1:mesh%nod_in_elem2D_num(n),n)
       read(fileID, '(a)') line
       read(line,*) nn, nod_in_elem2D_num(n)
       if(nod_in_elem2D_num(n) > MAX_NOD_IN_ELEM) then
          write(*,*) 'ERROR: please increase MAX_NOD_IN_ELEM to', nod_in_elem2D_num(n)
          stop
       end if
       read(line,*) nn, nod_in_elem2D_num(n), nod_in_elem2D(1:nod_in_elem2D_num(n), n)
    end do
    close(fileID)

  end function read_nod2D_elems

  function read_elem2D_nodes(mype, elem2D_nodes) result(myDim_elem2D)
    implicit none
    integer, intent(in) :: mype
    integer, intent(inout), dimension(:,:), allocatable :: elem2D_nodes
    integer :: myDim_elem2D
    character(MAX_PATH) :: file_name
    integer :: fileID = 12, n, nn

    myDim_elem2D = 0
    write(file_name, '(i8)') mype
    file_name='loops_elem2D_nodes_'//trim(adjustl(file_name))//'.dat'
    open(fileID, file=file_name, status='old', action='read')

    read(fileID, *) myDim_elem2D
    write(*,*) "loading",myDim_elem2D, " element nodes"

    allocate(elem2D_nodes(3, myDim_elem2D))

    do n=1,myDim_elem2D !+ edim_nod2d
       ! TODO check nz<=MAX_LEVELS
       read(fileID, *) nn, elem2D_nodes(:, n)
    end do

    close(fileID)
  end function read_elem2D_nodes


  function read_elem2D_levels(mype, ulevels, nlevels) result(myDim_elem2D)
    implicit none
    integer, intent(in) :: mype
    integer, intent(inout), dimension(:), allocatable :: ulevels
    integer, intent(inout), dimension(:), allocatable :: nlevels
    integer :: myDim_elem2D
    character(MAX_PATH) :: file_name
    integer :: fileID = 13, n, nn, nz

    myDim_elem2D = 0

    write(file_name, '(i8)') mype
    file_name='loops_elem2D_levels_'//trim(adjustl(file_name))//'.dat'
    open(fileID, file=file_name, status='old', action='read')
    read(fileID, '(i8)') myDim_elem2D

    write(*,*) "loading",myDim_elem2D, " element levels"

    allocate(ulevels(myDim_elem2D))
    allocate(nlevels(myDim_elem2D))

    do n=1, myDim_elem2D
       read(fileID, *) nn, nz, ulevels(n), nlevels(n)
       nlevels(n) = nlevels(n) + 1
    end do
    close(fileID)

  end function read_elem2D_levels

  function read_edges_nodes(mype, edges_nodes) result(myDim_edges)
    implicit none
    integer, intent(in) :: mype
    integer, intent(inout), dimension(:,:), allocatable :: edges_nodes
    integer :: myDim_edges
    character(MAX_PATH) :: file_name
    integer :: fileID = 14, n, nn

    myDim_edges = 0
    write(file_name, '(i8)') mype
    file_name='loops_edges_nodes_'//trim(adjustl(file_name))//'.dat'
    open(fileID, file=file_name, status='old', action='read')

    read(fileID, *) myDim_edges
    write(*,*) "loading",myDim_edges, " edge nodes"

    allocate(edges_nodes(2, myDim_edges))

    do n=1,myDim_edges !+ edim_nod2d
       ! TODO check nz<=MAX_LEVELS
       read(fileID, *) nn, edges_nodes(:, n)
    end do

    close(fileID)
  end function read_edges_nodes


  function read_edges_levels(mype, ulevels, nlevels) result(myDim_edges)
    implicit none
    integer, intent(in) :: mype
    integer, intent(inout), dimension(:), allocatable :: ulevels
    integer, intent(inout), dimension(:), allocatable :: nlevels
    integer :: myDim_edges
    character(MAX_PATH) :: file_name
    integer :: fileID = 15, n, nn, nz

    myDim_edges = 0

    write(file_name, '(i8)') mype
    file_name='loops_edges_levels_'//trim(adjustl(file_name))//'.dat'
    open(fileID, file=file_name, status='old', action='read')
    read(fileID, '(i8)') myDim_edges

    write(*,*) "loading",myDim_edges, " edge levels"

    allocate(ulevels(myDim_edges))
    allocate(nlevels(myDim_edges))

    do n=1, myDim_edges
       read(fileID, *) nn, nz, ulevels(n), nlevels(n)
    end do
    close(fileID)

  end function read_edges_levels

  function read_edges_tri(mype, edges_tri) result(myDim_edges)
    implicit none
    integer, intent(in) :: mype
    integer, intent(inout), dimension(:,:), allocatable :: edges_tri
    integer :: myDim_edges
    character(MAX_PATH) :: file_name
    integer :: fileID = 15, n, nn

    myDim_edges = 0
    write(file_name, '(i8)') mype
    file_name='loops_edges_edge_tri_'//trim(adjustl(file_name))//'.dat'
    open(fileID, file=file_name, status='old', action='read')

    read(fileID, *) myDim_edges
    write(*,*) "loading",myDim_edges, " edge tri"

    allocate(edges_tri(2, myDim_edges))

    do n=1,myDim_edges !+ edim_nod2d
       read(fileID, *) nn, edges_tri(:, n)
    end do

    close(fileID)
  end function read_edges_tri

end module io_mod
