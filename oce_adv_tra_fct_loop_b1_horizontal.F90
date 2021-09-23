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

        double precision :: a, s_
        integer :: e_
        integer :: nl

        integer :: n_node2edges
        integer, dimension(:), allocatable :: csr_node2edges !(0:myDim_nod2D)
        integer, dimension(:), allocatable :: csr_node2edges_ID
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

              do nz=nu12, nl12
                 fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0d0, adf_h(nz,edge))
                 fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0d0, adf_h(nz,edge))
                 fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0d0,-adf_h(nz,edge))
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
           do edge=1, mydim_edge2d
              enodes(1:2)=edges(:,edge)

              nl12 = nlevels(edge)
              nu12 = ulevels(edge)

              do nz=nu12, nl12
                 fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0d0, adf_h(nz,edge))
                 fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0d0, adf_h(nz,edge))
                 fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0d0,-adf_h(nz,edge))
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

              do nz=ulevels(edge), nlevels(edge)
                 a = adf_h(nz,edge)
                 fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0d0, a)
                 fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0d0, a)
                 fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0d0,-a)
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

#ifdef _OPENACC
#if 1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "ACC: PRECOMPUTED LEVELS+store: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !$acc enter data copyin(nlevels,ulevels,edges,fct_plus,fct_minus,adf_h)
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           !$acc parallel loop present(nlevels,ulevels,edges,fct_plus,fct_minus,adf_h)&
           !$acc& private(enodes)
           do edge=1, mydim_edge2d
              enodes(1:2)=edges(:,edge)

              nl12 = nlevels(edge)
              nu12 = ulevels(edge)

              do nz=nu12, nl12
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
        !$acc exit data delete(nlevels,ulevels,edges,adf_h)
        !$acc exit data copyout(fct_plus,fct_minus)
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
#endif

#if 1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "ACC: PRECOMPUTED node2edges+store: iterating over",MAX_ITERATIONS, " iterations..."
        n_node2edges=0
        allocate(csr_node2edges(0:myDim_nod2D))
        csr_node2edges(:) = 0
        do edge=1, mydim_edge2d
           enodes(1:2)=edges(:,edge)

           csr_node2edges(enodes(1)) = csr_node2edges(enodes(1)) + 1
           csr_node2edges(enodes(2)) = csr_node2edges(enodes(2)) + 1
        end do
        do n=1, myDim_nod2D
           csr_node2edges(n) = csr_node2edges(n-1) + csr_node2edges(n)
        end do
        allocate(csr_node2edges_ID(csr_node2edges(myDim_nod2D)))
        do n=myDim_nod2D, 1, -1
           csr_node2edges(n) = csr_node2edges(n-1)
        end do
        do edge=1, mydim_edge2d
           enodes(1:2)=edges(:,edge)

           csr_node2edges(enodes(1)) = csr_node2edges(enodes(1)) + 1
           csr_node2edges(enodes(2)) = csr_node2edges(enodes(2)) + 1
           csr_node2edges_ID(csr_node2edges(enodes(1))) =  edge
           csr_node2edges_ID(csr_node2edges(enodes(2))) = -edge
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !$acc enter data copyin(nlevels,ulevels,fct_plus,fct_minus,adf_h)
        !$acc enter data copyin(csr_node2edges, csr_node2edges_ID)
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           !$acc parallel loop present(nlevels,ulevels,fct_plus,fct_minus,adf_h, csr_node2edges, csr_node2edges_ID)
           do n=1, myDim_nod2D
              do nn=csr_node2edges(n-1)+1, csr_node2edges(n)
                 edge = csr_node2edges_ID(nn)
                 e_ = abs(edge)
                 s_ = edge/e_
                 !if(edge > 0) then
                    nl12 = nlevels(e_)
                    nu12 = ulevels(e_)
                    do nz=nu12, nl12
                       a = s_*adf_h(nz,e_)
                       fct_plus (nz,n)=fct_plus (nz,n) + max(0.0d0, a)
                       fct_minus(nz,n)=fct_minus(nz,n) + min(0.0d0, a)
                    end do
                 !else
                 !   edge = -edge
                 !   nl12 = nlevels(edge)
                 !   nu12 = ulevels(edge)
                 !   do nz=nu12, nl12
                 !      a = adf_h(nz,edge)
                 !      fct_plus (nz,n)=fct_plus (nz,n) + max(0.0d0, -a)
                 !      fct_minus(nz,n)=fct_minus(nz,n) + min(0.0d0, -a)
                 !   end do
                 !end if
              end do
           end do
        end do
        delta=wallclock()-t1
        !$acc exit data delete(nlevels,ulevels,edges,adf_h)
        !$acc exit data copyout(fct_plus,fct_minus)
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
#endif
#endif

#if 1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "ACC: collapse: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nl = maxval(nlevels(:))
        !$acc enter data copyin(nlevels,ulevels,edges,fct_plus,fct_minus,adf_h)
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           !$acc parallel loop present(nlevels,ulevels,edges,fct_plus,fct_minus,adf_h)&
           !$acc& private(enodes)&
           !$acc& collapse(2)
           do edge=1, mydim_edge2d
              do nz=1, nl
                 enodes(1:2)=edges(:,edge)

                 nl12 = nlevels(edge)
                 nu12 = ulevels(edge)

                 !do nz=nu12, nl12
                 if(nu12 <= nz .and. nz <= nl12) then
                    a = adf_h(nz,edge)
                    !$acc atomic update
                    fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0d0, a)
                    !$acc atomic update
                    fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0d0, a)
                    !$acc atomic update
                    fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0d0,-a)
                    !$acc atomic update
                    fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0d0,-a)
                 end if
              end do
           end do
        end do
        delta=wallclock()-t1
        !$acc exit data delete(nlevels,ulevels,edges,adf_h)
        !$acc exit data copyout(fct_plus,fct_minus)
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
#endif

end program oce_adv_tra_fct_loop_b1_horizontal
