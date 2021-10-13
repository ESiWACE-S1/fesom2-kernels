module distributed_mod
  use mpi
  implicit none
  private

  public :: exchange_nod3D
  public :: exchange_nod3D_begin
  public :: exchange_nod3D_end

  integer, parameter :: MAX_NEIGHBOR_PARTITIONS=32
  type com_struct
     integer                                       :: rPEnum ! the number of PE I receive info from
     integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: rPE    ! their list
     integer, dimension(MAX_NEIGHBOR_PARTITIONS+1) :: rptr   ! allocatables to the list of nodes
     integer, dimension(:), allocatable            :: rlist  ! the list of nodes
     integer                                       :: sPEnum ! send part
     integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: sPE
     integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: sptr
     integer, dimension(:), allocatable            :: slist
     integer, dimension(:), allocatable            :: req    ! request for MPI_Wait
     integer                                       :: nreq   ! number of requests for MPI_Wait
                                                             ! (to combine halo exchange of several fields)
  end type com_struct
  public :: com_struct

  type(com_struct)   :: com_nod2D
  public :: com_nod2D

  ! Nodal fields (2D; 2D integer; 3D with nl-1 or nl levels, one, two, or three values)
  integer, allocatable       :: s_mpitype_nod3D(:,:,:), r_mpitype_nod3D(:,:,:)

  ! general MPI part
  integer,public     :: MPIERR
  integer,public     :: npes
  integer,public     :: mype
  integer            :: maxPEnum=100
  integer, allocatable, dimension(:)  :: part

  public :: load_part

  contains
subroutine load_part(filename, n_parts, parts_id)
  implicit none
  character(len=*), intent(in) :: filename
  integer, dimension(:), allocatable, intent(inout) :: n_parts
  integer, dimension(:), allocatable, intent(inout) :: parts_id

  integer :: nparts, n, nn
  integer :: fileID = 18
  integer, dimension(:), allocatable :: parts

  open(fileID, file=filename, status='old', action='read')
  read(fileID, *) nparts
  allocate(n_parts(0:nparts))
  read(fileID, *) n_parts(1:nparts)
  n_parts(0) = 0
  do n=1, nparts
   n_parts(n) = n_parts(n) + n_parts(n-1)
  end do
  allocate(parts(n_parts(nparts)))
  read(fileID, *) parts(1:n_parts(nparts))
  close(fileID)

  write(*,*) 'nparts', nparts
  write(*,*) n_parts(1:nparts)

  allocate(parts_id(n_parts(npes)))
  parts_id = -1
  do n=1, npes
    do nn=n_parts(n-1)+1, n_parts(n)
      parts_id(parts(nn)) = n
    end do
  end do

  deallocate(parts)

end subroutine load_part

! ========================================================================
subroutine exchange_nod3D(nod1_array3D, nod2_array3D)
  implicit none
  double precision, intent(inout)  :: nod1_array3D(:,:)
  double precision, intent(inout)  :: nod2_array3D(:,:)

  call exchange_nod3D_begin(nod1_array3D, nod2_array3D)
  call exchange_nod3D_end
 
end subroutine exchange_nod3D

! ========================================================================
subroutine exchange_nod3D_begin(nod1_array3D,nod2_array3D)
implicit none

double precision, intent(inout) :: nod1_array3D(:,:) 
double precision, intent(inout) :: nod2_array3D(:,:) 
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
 
 integer  :: n, sn, rn
 integer  :: nz, nl1, nl2

    sn=com_nod2D%sPEnum
    rn=com_nod2D%rPEnum

    nl1 = ubound(nod1_array3D,1)

    if ((nl1<ubound(r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(r_mpitype_nod3D, 2))) then
       if (mype==0) then
          print *,'Subroutine exchange_nod3D not implemented for',nl1,'layers.'
          print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
       endif
       call MPI_Abort( MPI_COMM_WORLD, 1, MPIerr )
    endif

    nl2 = ubound(nod2_array3D,1)
    if ((nl2<ubound(r_mpitype_nod3D, 2)-1) .or. (nl2>ubound(r_mpitype_nod3D, 2))) then
       if (mype==0) then
          print *,'Subroutine exchange_nod3D not implemented for',nl2,'layers.'
          print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
       endif
       call MPI_Abort( MPI_COMM_WORLD, 1 , MPIerr)
    endif

    do n=1,rn    
      call MPI_Irecv(nod1_array3D, 1, r_mpitype_nod3D(n,nl1,1), com_nod2D%rPE(n), &
           com_nod2D%rPE(n),      MPI_COMM_WORLD, com_nod2D%req(2*n-1), MPIerr)  

      call MPI_Irecv(nod2_array3D, 1, r_mpitype_nod3D(n,nl2,1), com_nod2D%rPE(n), &
           com_nod2D%rPE(n)+npes, MPI_COMM_WORLD, com_nod2D%req(2*n  ), MPIerr) 
    end do

    do n=1, sn
      call MPI_Isend(nod1_array3D, 1, s_mpitype_nod3D(n,nl1,1), com_nod2D%sPE(n), &
           mype,      MPI_COMM_WORLD, com_nod2D%req(2*rn+2*n-1), MPIerr)

      call MPI_Isend(nod2_array3D, 1, s_mpitype_nod3D(n,nl2,1), com_nod2D%sPE(n), &
           mype+npes, MPI_COMM_WORLD, com_nod2D%req(2*rn+2*n), MPIerr)
    end do

    com_nod2D%nreq = 2*(rn+sn)

end subroutine exchange_nod3D_begin
!=======================================
! AND WAITING
!=======================================
 
subroutine exchange_nod3D_end
  implicit none
  call MPI_Waitall(com_nod2D%nreq, com_nod2D%req, MPI_STATUSES_IGNORE, MPIerr)
end subroutine exchange_nod3D_end

!=======================================================================
!!! communication_nodn(mesh)
!!subroutine init_com_nod2D(com_nod2D)
!!  implicit none
!!  type(com_struct), intent(in), target :: com_nod2D
!!  integer                  :: n, np, prank, el, r_count, s_count, q, i, j, nod, k, l
!!  integer                  :: num_send(0:npes-1), num_recv(0:npes-1), nd_count
!!  integer, allocatable     :: recv_from_pe(:), send_to_pes(:,:)
!!  logical                  :: max_laendereck_too_small=.false.
!!  integer                  :: IERR
!!  ! Assume we have 2D partitioning vector in part. Find communication rules
!!  ! Reduce allocation: find all neighboring PE
!!
!!  nd_count = count(part(1:nod2d) == mype)
!!  allocate(recv_from_pe(nod2d), send_to_pes(MAX_LAENDERECK,nd_count), &
!!       myList_nod2D(nd_count), STAT=IERR)
!!  if (IERR /= 0) then
!!     write (*,*) 'Could not allocate arrays in communication_nodn'
!!     stop
!!  endif
!!
!!  nd_count = 0
!!  do n=1,nod2D
!!     ! Checks if element el has nodes that belong to different partitions
!!     if (part(n) == mype) then
!!        nd_count = nd_count+1
!!        myList_nod2D(nd_count)=n
!!     endif
!!  end do
!!  myDim_nod2D=nd_count
!!
!!  num_send(0:npes-1) = 0
!!  num_recv(0:npes-1) = 0
!!  recv_from_pe(1:nod2d) = -1
!!  send_to_pes(1:MAX_LAENDERECK,1:nd_count) = -1
!!
!!  ! For the local nodes, run through the patch and collect all nodes in the patch
!!  ! (It would be simpler to re-use the adjacency matrix, but it is not a global
!!  !  variable... and I am lazy and want to keep the data structure simple)
!!
!!  do l=1,nd_count
!!     n = myList_nod2D(l)
!!     do i = 1, nod_in_elem2D_num(n)
!!        ! Over all elements el that the node n is part of
!!        el = nod_in_elem2D(i,n)
!!
!!        ! Checks, if elements are quads or triangles
!!        q = 4    ! quads as default
!!        if (elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3  ! triangle
!!
!!        do j = 1, q
!!           ! Over all nodes in every element el
!!           nod = elem2D_nodes(j,el)
!!           
!!           ! Checks, if node j is not in another partitionen
!!           if (part(nod) /= mype) then
!!              ! Checks, if not already considered to be received from this
!!              ! node from partition part(nod)
!!              if (recv_from_pe(nod) == -1) then  ! nod already collected to be received?
!!                 ! We have to receive this node from part(nod)
!!                 ! Add plus one node to the total number of
!!                 num_recv(part(nod)) = num_recv(part(nod)) + 1
!!                 ! ???
!!                 recv_from_pe(nod) = part(nod)   ! recv_from_pe(recv_count) = nod  ! no new information, just handy
!!              endif
!!              ! Checks, if all possible connected partition
!!              ! And we have to send n to part(nod). Do we know this already?
!!              do k=1,MAX_LAENDERECK    !???
!!                 if (send_to_pes(k,l) == part(nod)) then
!!                    exit  ! already collected
!!                 elseif (send_to_pes(k,l) == -1) then
!!                    send_to_pes(k,l) = part(nod)
!!                    num_send(part(nod)) = num_send(part(nod)) + 1
!!                    exit
!!                 elseif (k== MAX_LAENDERECK) then
!!                    max_laendereck_too_small = .true.  ! Problem
!!                 endif
!!              enddo
!!           endif
!!        enddo
!!     enddo
!!  enddo
!!
!!  if (max_laendereck_too_small) then
!!     print *,'Increase MAX_LAENDERECK in gen_modules_partitioning.F90 and recompile'
!!     stop
!!  endif
!!
!!  ! Now, build the send and recv communication data structure
!!  ! To how many PE needs to be send and which PEs
!!!!$  com_nod2D%rPEnum = count(num_recv(0:npes-1) > 0)
!!!!$  com_nod2D%sPEnum = count(num_send(0:npes-1) > 0)
!!
!!!!$  if (com_nod2D%rPEnum > MAX_NEIGHBOR_PARTITIONS .or.  &
!!!!$       com_nod2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
!!!!$     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
!!!!$     stop
!!!!$  endif
!!!!$  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
!!!!$  allocate(com_nod2D%sPE(com_nod2D%sPEnum))
!!
!!  r_count = 0
!!  s_count = 0
!!  com_nod2D%rptr(1) = 1
!!  com_nod2D%sptr(1) = 1
!!
!!  do np = 0, npes-1
!!     if(num_recv(np) /= 0) then
!!        r_count = r_count+1
!!        com_nod2D%rPE(r_count) = np
!!        com_nod2D%rptr(r_count+1) =  com_nod2D%rptr(r_count)+ num_recv(np)
!!     end if
!!     if(num_send(np) /= 0) then
!!        s_count = s_count+1
!!        com_nod2D%sPE(s_count) = np
!!        com_nod2D%sptr(s_count+1) =  com_nod2D%sptr(s_count)+ num_send(np)
!!     end if
!!  enddo
!!  com_nod2D%rPEnum = r_count
!!  com_nod2D%sPEnum = s_count
!!  if (com_nod2D%rPEnum > MAX_NEIGHBOR_PARTITIONS .or.  &
!!       com_nod2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
!!     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
!!     stop
!!  endif
!!
!!  ! Counts the number of node for each partition PE mype has to send/receive
!!  ! In ascending order of PE number
!!!!$  r_count = 0
!!!!$  s_count = 0
!!!!$  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1)) 
!!!!$  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1))
!!
!!!!$  com_nod2D%rptr(1) = 1
!!!!$  com_nod2D%sptr(1) = 1
!!!!$
!!!!$  do r_count = 1, com_nod2D%rPEnum
!!!!$     np = com_nod2D%rPE(r_count)
!!!!$     com_nod2D%rptr(r_count+1) =  com_nod2D%rptr(r_count)+ num_recv(np)
!!!!$  enddo
!!!!$  do s_count = 1, com_nod2D%sPEnum
!!!!$     np = com_nod2D%sPE(s_count)
!!!!$     com_nod2D%sptr(s_count+1) =  com_nod2D%sptr(s_count)+ num_send(np)
!!!!$  enddo
!!
!!  ! Lists themselves
!!
!!  r_count = 0
!!  eDim_nod2D=com_nod2D%rptr(com_nod2D%rPEnum+1)-1   
!!  allocate(com_nod2D%rlist(eDim_nod2D), &
!!       com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1), STAT=IERR) 
!!  if (IERR /= 0) then
!!     write (*,*) 'Could not allocate arrays in communication_nodn'
!!     stop
!!  endif
!!
!!  do np = 1,com_nod2D%rPEnum
!!     prank = com_nod2D%rPE(np)
!!     do n = 1, nod2D
!!        if (recv_from_pe(n) == prank) then
!!           r_count = r_count+1
!!           com_nod2D%rlist(r_count) = n
!!        end if
!!     end do
!!  end do
!!
!!  s_count = 0
!!!!$  allocate(com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1)) 
!!  do np = 1,com_nod2D%sPEnum
!!     prank = com_nod2D%sPE(np)
!!     do l = 1, nd_count
!!        n = myList_nod2D(l)
!!        if(any(send_to_pes(:,l) == prank)) then 
!!           s_count = s_count+1
!!           com_nod2D%slist(s_count) = n
!!        end if
!!     end do
!!  end do
!!
!!  ! Summary of this piece: mype receives
!!  ! information on external 2D nodes from
!!  ! comm_nod2D%rPEnum external PEs
!!  ! Their ranks (numbers) are in array
!!  ! comm_nod2D%rPE(:)
!!  ! Pointers to external node numbers are in
!!  ! comm_nod2D%rptr(:)
!!  ! The node numbers are in 
!!  ! comm_nod2D%list(:)
!!  ! Putting everything into structure takes many operations, but
!!  ! without the structure we will need to many names and arrays
!!  ! Do not forget that we need also send part.
!!
!!  ! mype sends its data to
!!  ! comm_nod2D%sPEnum external PEs
!!  ! Their ranks (numbers) are in array
!!  ! comm_nod2D%sPE(:)
!!  ! Pointers to external node numbers are in
!!  ! comm_nod2D%sptr(:)
!!  ! The node numbers are in 
!!  ! comm_nod2D%list(:)
!!
!!  deallocate(recv_from_pe, send_to_pes)
!!end subroutine init_com_nod2D


end module distributed_mod
