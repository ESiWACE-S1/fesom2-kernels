program oce_adv_tra_fct_loop_b1_horizontal
  use wallclock_mod
  use io_mod
        implicit none

        integer :: n, nu12, nl12, nz, myDim_edge2D, n_it, edge, myDim_nod2D, myDim_elem2D
        integer, dimension(:), allocatable :: ulevels
        integer, dimension(:), allocatable :: nlevels
        integer, dimension(:), allocatable :: ulevels_elem2D
        integer, dimension(:), allocatable :: nlevels_elem2D
        integer, dimension(:,:), allocatable :: edges
        integer, dimension(:,:), allocatable :: edge_tri
        integer, dimension(2) :: enodes
        integer :: nl1, nu1, nl2, nu2, el(2)
        double precision, dimension(:,:), allocatable :: fct_plus
        double precision, dimension(:,:), allocatable :: fct_minus
        double precision, dimension(:,:), allocatable :: adf_h
        integer, parameter :: MAX_LEVELS=50
        integer, parameter :: MAX_ITERATIONS=1000
        integer :: mype, nn
        ! LBS
        integer, dimension(:), allocatable :: ulevels_nod2D
        integer, dimension(:), allocatable :: nlevels_nod2D
        integer, dimension(:,:), allocatable :: nodenz_lbs
        integer :: n_elemnz, nlevels_sum

        double precision :: a
        !https://stackoverflow.com/a/6880672
        real(8)::t1,delta, delta_orig

        mype = 0
        myDim_edge2D = read_edges_nodes(mype, edges)
        myDim_edge2D = read_edges_levels(mype, ulevels, nlevels)
        myDim_edge2D = read_edges_tri(mype, edge_tri)

        myDim_nod2D = get_myDim_nod2D(mype)

        myDim_elem2D = read_elem2D_levels(mype, ulevels_elem2D, nlevels_elem2D)

        allocate(fct_plus(MAX_LEVELS, myDim_nod2D))
        allocate(fct_minus(MAX_LEVELS, myDim_nod2D))
        allocate(adf_h(MAX_LEVELS, myDim_edge2D))

        write(*,*) "COMPUTE: fct_minus/fct_plus(nz,enodes)+=max(0,adf_h(nz,edge))"
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           !$acc parallel loop gang present(nlevels,ulevels,edges,fct_plus,fct_minus,adf_h)&
           !$acc& private(enodes,nl12,nu12)&
#ifdef with_acc_vector_length
           !$acc& vector_length(z_vector_length)&
#endif
#ifdef with_acc_async
           !$acc& async(stream_hor_adv_tra)&
#endif
           !$acc
           do edge=1, mydim_edge2d
              enodes(1:2)=edges(:,edge)
              el=edge_tri(:,edge)
              nu1=ulevels_elem2D(el(1)); nl1=nlevels_elem2D(el(1))-1
              nl2=0             ; nu2=0
              if(el(2)>0) then
                 nu2=ulevels_elem2D(el(2))
                 nl2=nlevels_elem2D(el(2))-1
              end if
              nu12 = nu1        ; nl12 = max(nl1,nl2)
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
        end do
        delta=wallclock()-t1
        delta_orig = delta
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "PRECOMPUTED LEVELS: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           !$acc parallel loop gang present(nlevels,ulevels,edges,fct_plus,fct_minus,adf_h)&
           !$acc& private(enodes,nl12,nu12)&
#ifdef with_acc_vector_length
           !$acc& vector_length(z_vector_length)&
#endif
#ifdef with_acc_async
           !$acc& async(stream_hor_adv_tra)&
#endif
           !$acc
           do edge=1, mydim_edge2d
              enodes(1:2)=edges(:,edge)

              nl12 = nlevels(edge)
              nu12 = ulevels(edge)

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
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "PRECOMPUTED LEVELS+store: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           do edge=1, mydim_edge2d
              enodes(1:2)=edges(:,edge)

              !$acc loop vector
              do nz=ulevels(edge), nlevels(edge)
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
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        myDim_nod2D = read_nod2D_levels(mype, ulevels_nod2D, nlevels_nod2D)
        allocate(nodenz_lbs(MAX_LEVELS, myDim_nod2D))

        nlevels_sum = 0
        nodenz_lbs=-1
        nn = 0
        do n=1,myDim_nod2D
           nu1 = ulevels_nod2D(n)
           nl1 = nlevels_nod2D(n)
           do nz=nu1, nl1-1
              nn = nn + 1
              nodenz_lbs(nz,n) = nn
           end do
           nlevels_sum = nlevels_sum + nl1-1-nu1+1
        end do

end program oce_adv_tra_fct_loop_b1_horizontal
