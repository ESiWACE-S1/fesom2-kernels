program oce_adv_tra_fct
  use wallclock_mod
  use io_mod
  implicit none

  integer :: n, nu1, nl1, nz, myDim_nod2D, n_it
  integer, dimension(:), allocatable :: ulevels_nod2D
  integer, dimension(:), allocatable :: nlevels_nod2D
  double precision, dimension(:,:), allocatable :: fct_ttf_min
  double precision, dimension(:,:), allocatable :: fct_ttf_max
  double precision, dimension(:,:), allocatable :: LO
  double precision, dimension(:,:), allocatable :: ttf
  integer, parameter :: MAX_LEVELS=50
  integer, parameter :: MAX_ITERATIONS=200
  integer :: mype, nn

  integer :: myDim_elem2D
  integer :: elem, nl, nl12, nu12
  integer, dimension(:), allocatable :: ulevels
  integer, dimension(:), allocatable :: nlevels
  integer, dimension(3) :: enodes
  integer, dimension(:,:), allocatable :: elem2D_nodes
  double precision, dimension(:,:,:), allocatable :: UV_rhs
  double precision :: bignumber = huge(1.0)

  integer, dimension(:), allocatable :: nod_in_elem2D_num
  integer, dimension(:,:), allocatable :: nod_in_elem2D
  double precision, dimension(:), allocatable :: tvert_min
  double precision, dimension(:), allocatable :: tvert_max

  integer :: nl2, nu2, el(2), edge
  integer :: myDim_edge2D
  integer, dimension(:,:), allocatable :: edges
  integer, dimension(:,:), allocatable :: edge_tri
  integer, dimension(:), allocatable :: ulevels_edge
  integer, dimension(:), allocatable :: nlevels_edge
  double precision, dimension(:,:), allocatable :: fct_plus
  double precision, dimension(:,:), allocatable :: fct_minus
  double precision, dimension(:,:), allocatable :: adf_h, adf_v

  double precision :: dt, flux, flux_eps
  double precision, dimension(:,:), allocatable :: areasvol
  
  real(8)::t1,delta, delta_orig

  mype = 0
  myDim_nod2D = read_nod2D_levels(mype, ulevels_nod2D, nlevels_nod2D)

  allocate(fct_ttf_min(MAX_LEVELS, myDim_nod2D))
  allocate(fct_ttf_max(MAX_LEVELS, myDim_nod2D))
  allocate(LO(MAX_LEVELS, myDim_nod2D))
  allocate(ttf(MAX_LEVELS, myDim_nod2D))
  allocate(areasvol(MAX_LEVELS, myDim_nod2D))

  myDim_elem2D = read_elem2D_nodes(mype, elem2D_nodes)
  myDim_elem2D = read_elem2D_levels(mype, ulevels, nlevels)

  allocate(UV_rhs(2,MAX_LEVELS, myDim_elem2D))

  myDim_nod2D = read_nod2D_elems(mype, nod_in_elem2D_num, nod_in_elem2D)

  allocate(tvert_min(MAX_LEVELS))
  allocate(tvert_max(MAX_LEVELS))

  myDim_edge2D = read_edges_nodes(mype, edges)
  myDim_edge2D = read_edges_levels(mype, ulevels_edge, nlevels_edge)
  myDim_edge2D = read_edges_tri(mype, edge_tri)

  allocate(fct_plus(MAX_LEVELS, myDim_nod2D))
  allocate(fct_minus(MAX_LEVELS, myDim_nod2D))
  allocate(adf_h(MAX_LEVELS, myDim_edge2D))
  allocate(adf_v(MAX_LEVELS, myDim_edge2D))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  t1=wallclock()
  do n_it=1, MAX_ITERATIONS
     ! loop a1
     do n=1,myDim_nod2D !+ edim_nod2d
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        !$acc loop vector
        do nz=nu1, nl1-1
           fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
           fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
        end do
     end do

     ! loop a2
     do elem=1, myDim_elem2D
        enodes=elem2D_nodes(:,elem)
        nu1 = ulevels(elem)
        nl1 = nlevels(elem)
        !$acc loop vector
        do nz=nu1, nl1-1
           UV_rhs(1,nz,elem)=maxval(fct_ttf_max(nz,enodes))
           UV_rhs(2,nz,elem)=minval(fct_ttf_min(nz,enodes))
        end do
        if (nl1<=nl-1) then
           !$acc loop vector
           do nz=nl1,nl-1
              UV_rhs(1,nz,elem)=-bignumber
              UV_rhs(2,nz,elem)= bignumber
           end do
        endif
     end do ! --> do elem=1, myDim_elem2D

     ! loop a1 vlimit=1
     do n=1, myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)

        !___________________________________________________________________
        !$acc loop vector
        do nz=nu1,nl1-1
           ! max,min horizontal bound in cluster around node n in every
           ! vertical layer
           ! nod_in_elem2D     --> elem indices of which node n is surrounded
           ! nod_in_elem2D_num --> max number of surrounded elem
           tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
           tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
        end do

        !___________________________________________________________________
        ! calc max,min increment of surface layer with respect to low order
        ! solution
        fct_ttf_max(nu1,n)=tvert_max(nu1)-LO(nu1,n)
        fct_ttf_min(nu1,n)=tvert_min(nu1)-LO(nu1,n)

        ! calc max,min increment from nz-1:nz+1 with respect to low order
        ! solution at layer nz
        !$acc loop vector
        do nz=nu1+1,nl1-2
           fct_ttf_max(nz,n)=maxval(tvert_max(nz-1:nz+1))-LO(nz,n)
           fct_ttf_min(nz,n)=minval(tvert_min(nz-1:nz+1))-LO(nz,n)
        end do
        ! calc max,min increment of bottom layer -1 with respect to low order
        ! solution
        nz=nl1-1
        fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
        fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)
     end do

     ! loop b1
     do n=1, myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        !$acc loop vector
        do nz=nu1,nl1-1
           fct_plus(nz,n)=0.
           fct_minus(nz,n)=0.
        end do
     end do
     ! vertical
     do n=1, myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        do nz=nu1,nl1-1
           fct_plus(nz,n) =fct_plus(nz,n) +(max(0.0d0,adf_v(nz,n))+max(0.0d0,-adf_v(nz+1,n)))
           fct_minus(nz,n)=fct_minus(nz,n)+(min(0.0d0,adf_v(nz,n))+min(0.0d0,-adf_v(nz+1,n)))
        end do
     end do
     ! horizontal
     do edge=1, myDim_edge2D
        enodes(1:2)=edges(:,edge)
        el=edge_tri(:,edge)
        nl1=nlevels(el(1))-1
        nu1=ulevels(el(1))
        nl2=0
        nu2=0
        if(el(2)>0) then
           nl2=nlevels(el(2))-1
           nu2=ulevels(el(2))
        end if

        nl12 = max(nl1,nl2)
        nu12 = nu1
        if (nu2>0) nu12 = min(nu1,nu2)

        !$acc loop vector
        do nz=nu12, nl12
           !$acc atomic update
           fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0d0, adf_h(nz,edge))
           !$acc atomic update
           fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0d0, adf_h(nz,edge))
           !$acc atomic update
           fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0d0,-adf_h(nz,edge))
           !$acc atomic update
           fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0d0,-adf_h(nz,edge))
        end do
     end do


     ! loop b2
     do n=1,myDim_nod2D
        nu1=ulevels_nod2D(n)
        nl1=nlevels_nod2D(n)
        !$acc loop vector&
        !$acc& private(flux)
        do nz=nu1,nl1-1
           flux=fct_plus(nz,n)*dt/areasvol(nz,n)+flux_eps
           fct_plus(nz,n)=min(1.0d0,fct_ttf_max(nz,n)/flux)
           flux=fct_minus(nz,n)*dt/areasvol(nz,n)-flux_eps
           fct_minus(nz,n)=min(1.0d0,fct_ttf_min(nz,n)/flux)
        end do
     end do

     ! loop b3 (similar to a1)

  end do
  delta=wallclock()-t1
  delta_orig=delta
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta



end program oce_adv_tra_fct
