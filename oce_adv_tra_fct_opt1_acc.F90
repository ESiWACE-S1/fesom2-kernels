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
  double precision, dimension(:,:), allocatable :: U_rhs
  double precision, dimension(:,:), allocatable :: V_rhs
  double precision :: bignumber = huge(1.0)

  integer, dimension(:), allocatable :: nod_in_elem2D_num
  integer, dimension(:,:), allocatable :: nod_in_elem2D
  double precision, dimension(:), allocatable :: tvert_min
  double precision, dimension(:), allocatable :: tvert_max

  integer :: nl2, nu2, edge, el(2)
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
  
  double precision :: a, an, ae
  real(8)::t1,delta, delta_orig

  mype = 0
  myDim_nod2D = read_nod2D_levels(mype, ulevels_nod2D, nlevels_nod2D)
  !$acc enter data copyin(ulevels_nod2D, nlevels_nod2D)

  allocate(fct_ttf_min(MAX_LEVELS, myDim_nod2D))
  !$acc enter data create(fct_ttf_min)
  allocate(fct_ttf_max(MAX_LEVELS, myDim_nod2D))
  !$acc enter data create(fct_ttf_max)
  allocate(LO(MAX_LEVELS, myDim_nod2D))
  !$acc enter data create(LO)
  allocate(ttf(MAX_LEVELS, myDim_nod2D))
  !$acc enter data create(ttf)
  allocate(areasvol(MAX_LEVELS, myDim_nod2D))
  !$acc enter data create(areasvol)

  myDim_elem2D = read_elem2D_nodes(mype, elem2D_nodes)
  !$acc enter data copyin(elem2D_nodes)
  myDim_elem2D = read_elem2D_levels(mype, ulevels, nlevels)
  !$acc enter data copyin(ulevels, nlevels)

  allocate(U_rhs(MAX_LEVELS, myDim_elem2D))
  !$acc enter data create(U_rhs)
  allocate(V_rhs(MAX_LEVELS, myDim_elem2D))
  !$acc enter data create(V_rhs)

  myDim_nod2D = read_nod2D_elems(mype, nod_in_elem2D_num, nod_in_elem2D)
  !$acc enter data copyin(nod_in_elem2D_num, nod_in_elem2D)

  allocate(tvert_min(MAX_LEVELS))
  !$acc enter data create(tvert_min)
  allocate(tvert_max(MAX_LEVELS))
  !$acc enter data create(tvert_max)

  myDim_edge2D = read_edges_nodes(mype, edges)
  !$acc enter data copyin(edges)
  myDim_edge2D = read_edges_levels(mype, ulevels_edge, nlevels_edge)
  !$acc enter data copyin(ulevels_edge, nlevels_edge)
  myDim_edge2D = read_edges_tri(mype, edge_tri)
  !$acc enter data copyin(edge_tri)

  allocate(fct_plus(MAX_LEVELS, myDim_nod2D))
  !$acc enter data create(fct_plus)
  allocate(fct_minus(MAX_LEVELS, myDim_nod2D))
  !$acc enter data create(fct_minus)
  allocate(adf_h(MAX_LEVELS, myDim_edge2D))
  !$acc enter data create(adf_h)
  allocate(adf_v(MAX_LEVELS, myDim_edge2D))
  !$acc enter data create(adf_v)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  t1=wallclock()
  do n_it=1, MAX_ITERATIONS
     ! loop a1
     !$acc parallel loop present(ulevels_nod2D, nlevels_nod2D, LO, ttf, fct_ttf_max, fct_ttf_min)
     do n=1,myDim_nod2D !+ edim_nod2d
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        do nz=nu1, nl1-1
           fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
           fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
        end do
     end do


     ! loop a2
     !$acc parallel loop present(elem2D_nodes, ulevels, nlevels, U_rhs, V_rhs, fct_ttf_max, fct_ttf_min) &
     !$acc & private(enodes)
     do elem=1, myDim_elem2D
        enodes=elem2D_nodes(:,elem)
        nu1 = ulevels(elem)
        nl1 = nlevels(elem)
        do nz=nu1,nl1-1
           U_rhs(nz,elem)=maxval(fct_ttf_max(nz,enodes))
           V_rhs(nz,elem)=minval(fct_ttf_min(nz,enodes))
        end do
        do nz=nl1,nl-1
           U_rhs(nz,elem)=-bignumber
           V_rhs(nz,elem)= bignumber
        end do
     end do ! --> do elem=1, myDim_elem2D

     ! loop a1 vlimit=1
     !$acc parallel loop present(ulevels_nod2D, nlevels_nod2D, tvert_max, tvert_min, U_rhs, V_rhs &
     !$acc                     &, fct_ttf_max, fct_ttf_min)
     do n=1, myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        do nz=nu1,nl1-1
           tvert_max(nz) = huge(1.0d0)
           tvert_min(nz) = -huge(1.0d0)
        end do
        !$acc loop vector
        do nn=1, nod_in_elem2D_num(n)
           elem = nod_in_elem2D(nn,n)
           do nz=nu1,nl1-1
              tvert_max(nz)= max(tvert_max(nz), U_rhs(nz,elem))
              tvert_min(nz)= min(tvert_min(nz), V_rhs(nz,elem))
           end do
        end do

        !___________________________________________________________________
        ! calc max,min increment of surface layer with respect to low order
        ! solution
        fct_ttf_max(nu1,n)=tvert_max(nu1)-LO(nu1,n)
        fct_ttf_min(nu1,n)=tvert_min(nu1)-LO(nu1,n)

        ! calc max,min increment from nz-1:nz+1 with respect to low order
        ! solution at layer nz
        do nz=nu1+1,nl1-2
           fct_ttf_max(nz,n)=max(tvert_max(nz-1),tvert_max(nz), tvert_max(nz+1))-LO(nz,n)
           fct_ttf_min(nz,n)=min(tvert_min(nz-1),tvert_min(nz), tvert_min(nz+1))-LO(nz,n)
        end do
        ! calc max,min increment of bottom layer -1 with respect to low order
        ! solution
        nz=nl1-1
        fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
        fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)
     end do

     ! loop b1
     !$acc parallel loop present(ulevels_nod2D, nlevels_nod2D, fct_plus, fct_minus)
     do n=1, myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        do nz=nu1,nl1
           fct_plus(nz,n)=0.
           fct_minus(nz,n)=0.
        end do
     end do
     ! vertical
     !$acc parallel loop present(ulevels_nod2D, nlevels_nod2D, fct_plus, fct_minus, adf_v)
     do n=1, myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        do nz=nu1,nl1
           a = adf_v(nz,n)
           an= adf_v(nz+1,n)
           fct_plus(nz,n) =fct_plus(nz,n) +(max(0.0d0,a)+max(0.0d0,-an))
           fct_minus(nz,n)=fct_minus(nz,n)+(min(0.0d0,a)+min(0.0d0,-an))
        end do
     end do
     ! horizontal
     !$acc parallel loop present(nlevels_edge,ulevels_edge,edges,fct_plus,fct_minus,adf_h)&
     !$acc& private(enodes)
     do edge=1, myDim_edge2D
        enodes(1:2)=edges(:,edge)
        nu1 = ulevels_edge(edge)
        nl1 = nlevels_edge(edge)
        do nz=nu1,nl1
           a = adf_h(nz,edge)
           !$acc atomic update
           fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0d0, a)
           !$acc atomic update
           fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0d0, a)
           !$acc atomic update
           fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0d0,-a)
           !$acc atomic update
           fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0d0,-a)
        end do
     end do


     ! loop b2
     !$acc parallel loop present(ulevels_nod2D, nlevels_nod2D, fct_plus, fct_minus)
     do n=1,myDim_nod2D
        nu1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        do nz=nu1,nl1
           flux=fct_plus(nz,n)*dt/areasvol(nz,n)+flux_eps
           fct_plus(nz,n)=min(1.0d0,fct_ttf_max(nz,n)/flux)
           flux=fct_minus(nz,n)*dt/areasvol(nz,n)-flux_eps
           fct_minus(nz,n)=min(1.0d0,fct_ttf_min(nz,n)/flux)
        end do
     end do

     ! loop b3 (similar to a1)
     !$acc parallel loop present(ulevels_nod2D, nlevels_nod2D, adf_v, adf_h, fct_plus, fct_minus)
     do n=1, myDim_nod2D
        nu1=ulevels_nod2D(n)
        nl1=nlevels_nod2D(n)

        !_______________________________________________________________________
        nz=nu1
        ae=1.0d0
        flux=adf_v(nz,n)
        if(flux>=0.0d0) then
           ae=min(ae,fct_plus(nz,n))
        else
           ae=min(ae,fct_minus(nz,n))
        end if
        adf_v(nz,n)=ae*adf_v(nz,n) 

        !_______________________________________________________________________
        do nz=nu1+1,nl1-1
           ae=1.0d0
           flux=adf_v(nz,n)
           if(flux>=0.d0) then
              ae=min(ae,fct_minus(nz-1,n))
              ae=min(ae,fct_plus(nz,n))
           else
              ae=min(ae,fct_plus(nz-1,n))
              ae=min(ae,fct_minus(nz,n))
           end if
           adf_v(nz,n)=ae*adf_v(nz,n)
        end do
        ! the bottom flux is always zero
     end do

     !Horizontal

     !$acc parallel loop present(nlevels_edge,ulevels_edge,edges,fct_plus,fct_minus,adf_h)&
     !$acc& private(enodes)
     do edge=1, myDim_edge2D
        enodes(1:2)=edges(:,edge)
        nu1 = ulevels_edge(edge)
        nl1 = nlevels_edge(edge)
        do nz=nu1,nl1
           ae=1.0d0
           flux=adf_h(nz,edge)

           if(flux>=0.d0) then
              ae=min(ae,fct_plus(nz,enodes(1)))
              ae=min(ae,fct_minus(nz,enodes(2)))
           else
              ae=min(ae,fct_minus(nz,enodes(1)))
              ae=min(ae,fct_plus(nz,enodes(2)))
           endif

           adf_h(nz,edge)=ae*adf_h(nz,edge)
        end do
     end do


  end do
  delta=wallclock()-t1

  delta_orig=delta
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta



end program oce_adv_tra_fct
