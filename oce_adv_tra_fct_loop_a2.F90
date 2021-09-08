program oce_adv_tra_fct_loop_a2
  use wallclock_mod
  use io_mod
        implicit none

        integer :: n, nu1, nl1, nl, nz, myDim_elem2D, n_it, elem
        integer, dimension(:), allocatable :: ulevels
        integer, dimension(:), allocatable :: nlevels
        integer, dimension(3) :: enodes
        integer, dimension(:,:), allocatable :: elem2D_nodes
        double precision, dimension(:,:), allocatable :: fct_ttf_min
        double precision, dimension(:,:), allocatable :: fct_ttf_max
        double precision, dimension(:,:,:), allocatable :: UV_rhs
        double precision :: bignumber = huge(1.0)
        integer, parameter :: MAX_LEVELS=50
        integer, parameter :: MAX_ITERATIONS=1000
        integer :: mype, nn
        !
        double precision, dimension(:,:), allocatable :: U_rhs, V_rhs
        ! LBS
        integer, dimension(:), allocatable :: ulevels_nod2D
        integer, dimension(:), allocatable :: nlevels_nod2D
        integer :: myDim_nod2D, n_n1, n_n2, n_n3
        integer :: n_elemnz, nlevels_sum
        integer, dimension(:), allocatable :: csr_elemnz2nodenz !(0:n_elemnz)
        integer, dimension(:), allocatable :: csr_elemnz2nodenz_ID
        integer, dimension(:,:), allocatable :: nodenz_lbs
        double precision, dimension(:), allocatable :: U_rhs_lbs, V_rhs_lbs
        double precision, dimension(:), allocatable :: fct_ttf_min_lbs
        double precision, dimension(:), allocatable :: fct_ttf_max_lbs

        !https://stackoverflow.com/a/6880672
        real(8)::t1,delta, delta_orig

        mype = 0
        myDim_elem2D = read_elem2D_nodes(mype, elem2D_nodes)
        myDim_elem2D = read_elem2D_levels(mype, ulevels, nlevels)
        

        allocate(fct_ttf_min(MAX_LEVELS, myDim_elem2D))
        allocate(fct_ttf_max(MAX_LEVELS, myDim_elem2D))
        allocate(UV_rhs(2,MAX_LEVELS, myDim_elem2D))

        nl = maxval(nlevels(:)-ulevels(:)+1)

        write(*,*) "COMPUTE: UV_rhs(nz,elem)=maxval/minval(fct_ttf_max(nz,enodes))/bignumber"
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "DEFAULT iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !$acc parallel loop gang present(UV_rhs,nl,ulevels,nlevels,elem2D_nodes,fct_ttf_min,fct_ttf_max) &
                !$acc& private(enodes,nu1,nl1)&
#ifdef WITH_ACC_VECTOR_LENGTH
                !$acc& vector_length(z_vector_length)&
#endif
                !$acc
                do elem=1, myDim_elem2D
                        enodes=elem2D_nodes(:,elem)
                        nu1 = ulevels(elem)
                        nl1 = nlevels(elem)
                        !$acc loop vector
                        do nz=nu1, nl1-1
                                UV_rhs(1,nz,elem)=maxval(fct_ttf_max(nz,enodes))
                                UV_rhs(2,nz,elem)=minval(fct_ttf_min(nz,enodes))
                        end do
                        ! USELESS if (nl1<=nl-1) then
                                !$acc loop vector
                                do nz=nl1,nl-1
                                        UV_rhs(1,nz,elem)=-bignumber
                                        UV_rhs(2,nz,elem)= bignumber
                                end do
                        !endif
                end do ! --> do elem=1, myDim_elem2D
        end do
        delta=wallclock()-t1
        delta_orig = delta
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "MAXVAL->MAX iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !$acc parallel loop gang present(UV_rhs,nl,ulevels,nlevels,elem2D_nodes,fct_ttf_min,fct_ttf_max) &
                !$acc& private(enodes,nu1,nl1)&
#ifdef WITH_ACC_VECTOR_LENGTH
                !$acc& vector_length(z_vector_length)&
#endif
                !$acc
                do elem=1, myDim_elem2D
                        enodes=elem2D_nodes(:,elem)
                        nu1 = ulevels(elem)
                        nl1 = nlevels(elem)
                        !$acc loop vector
                        do nz=nu1, nl1-1
                                UV_rhs(1,nz,elem)=max(fct_ttf_max(nz,enodes(1)), fct_ttf_max(nz,enodes(2)), fct_ttf_max(nz,enodes(3)))
                                UV_rhs(2,nz,elem)=min(fct_ttf_min(nz,enodes(1)), fct_ttf_min(nz,enodes(2)), fct_ttf_min(nz,enodes(3)))
                        end do
                        !$acc loop vector
                        do nz=nl1,nl-1
                           UV_rhs(1,nz,elem)=-bignumber
                           UV_rhs(2,nz,elem)= bignumber
                        end do
                end do ! --> do elem=1, myDim_elem2D
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "U_rhs,V_rhs iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(U_rhs(MAX_LEVELS, myDim_elem2D))
        allocate(V_rhs(MAX_LEVELS, myDim_elem2D))

        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !$acc parallel loop gang present(UV_rhs,nl,ulevels,nlevels,elem2D_nodes,fct_ttf_min,fct_ttf_max) &
                !$acc& private(enodes,nu1,nl1)&
#ifdef WITH_ACC_VECTOR_LENGTH
                !$acc& vector_length(z_vector_length)&
#endif
                !$acc
                do elem=1, myDim_elem2D
                        enodes=elem2D_nodes(:,elem)
                        nu1 = ulevels(elem)
                        nl1 = nlevels(elem)
                        !$acc loop vector
                        do nz=nu1, nl1-1
                                U_rhs(nz,elem)=maxval(fct_ttf_max(nz,enodes))
                                V_rhs(nz,elem)=minval(fct_ttf_min(nz,enodes))
                                !U_rhs(nz,elem)=max(fct_ttf_max(nz,enodes(1)), fct_ttf_max(nz,enodes(2)), fct_ttf_max(nz,enodes(3)))
                                !V_rhs(nz,elem)=min(fct_ttf_min(nz,enodes(1)), fct_ttf_min(nz,enodes(2)), fct_ttf_min(nz,enodes(3)))
                        end do
                        !$acc loop vector
                        do nz=nl1,nl-1
                           U_rhs(nz,elem)=-bignumber
                           V_rhs(nz,elem)= bignumber
                        end do
                end do ! --> do elem=1, myDim_elem2D
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "ONLY local NL iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !$acc parallel loop gang present(UV_rhs,nl,ulevels,nlevels,elem2D_nodes,fct_ttf_min,fct_ttf_max) &
                !$acc& private(enodes,nu1,nl1)&
#ifdef WITH_ACC_VECTOR_LENGTH
                !$acc& vector_length(z_vector_length)&
#endif
                !$acc
                do elem=1, myDim_elem2D
                        enodes=elem2D_nodes(:,elem)
                        nu1 = ulevels(elem)
                        nl1 = nlevels(elem)
                        !$acc loop vector
                        do nz=nu1, nl1-1
                                U_rhs(nz,elem)=maxval(fct_ttf_max(nz,enodes))
                                V_rhs(nz,elem)=minval(fct_ttf_min(nz,enodes))
                                !U_rhs(nz,elem)=max(fct_ttf_max(nz,enodes(1)), fct_ttf_max(nz,enodes(2)), fct_ttf_max(nz,enodes(3)))
                                !V_rhs(nz,elem)=min(fct_ttf_min(nz,enodes(1)), fct_ttf_min(nz,enodes(2)), fct_ttf_min(nz,enodes(3)))
                        end do
                end do ! --> do elem=1, myDim_elem2D
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
        write(*,*) 'nlevels_sum', nlevels_sum, (real((nl-1+1)*myDim_elem2D-nlevels_sum))/real(nl*myDim_elem2D)
        write(*,*) 'fill nodenz_lbs with', nn, 'values'
        n_elemnz = 0
        do n=1,myDim_elem2D
           nu1 = ulevels(n)
           nl1 = nlevels(n)
           n_elemnz = n_elemnz + nl1-1-nu1+1
        enddo
        write(*,*) 'n_elemnz', n_elemnz
        allocate(csr_elemnz2nodenz(0:n_elemnz))
        nn = 0
        csr_elemnz2nodenz(nn) = 0
        do n=1,myDim_elem2D
           enodes=elem2D_nodes(:,n)
           nu1 = ulevels(n)
           nl1 = nlevels(n)
           do nz=nu1, nl1-1
              nn = nn + 1
              csr_elemnz2nodenz(nn) = csr_elemnz2nodenz(nn-1) + 3 ! 3 nodes at this level
           end do
           !do nz=nl1,nl-1
           !   nn = nn + 1
           !   csr_elemnz2nodenz(nn) = csr_elemnz2nodenz(nn-1) + 1 ! out of levels: set node id to -1
           !end do
        enddo
        write(*,*) 'csr_elemnz2nodenz', nn, n_elemnz, csr_elemnz2nodenz(n_elemnz)
        allocate(csr_elemnz2nodenz_ID(csr_elemnz2nodenz(n_elemnz)))
        nn = 0
        do n=1,myDim_elem2D
           enodes=elem2D_nodes(:,n)
           nu1 = ulevels(n)
           nl1 = nlevels(n)
           do nz=nu1, nl1-1
              nn = nn + 1; csr_elemnz2nodenz_ID(nn) = nodenz_lbs(nz, enodes(1))
              nn = nn + 1; csr_elemnz2nodenz_ID(nn) = nodenz_lbs(nz, enodes(2))
              nn = nn + 1; csr_elemnz2nodenz_ID(nn) = nodenz_lbs(nz, enodes(3))
           end do
           do nz=nl1,nl-1
              nn = nn + 1; csr_elemnz2nodenz_ID(nn) = -1
           end do
        enddo

        allocate(fct_ttf_max_lbs(-1:nlevels_sum))
        allocate(fct_ttf_min_lbs(-1:nlevels_sum))
        allocate(U_rhs_lbs(csr_elemnz2nodenz(n_elemnz)))
        allocate(V_rhs_lbs(csr_elemnz2nodenz(n_elemnz)))

        fct_ttf_max_lbs(-1) = -huge(1.0)
        fct_ttf_min_lbs(-1) = huge(1.0)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "LBS_TRANSFORM: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           !elemnz2nodenz
           nn=0
           do n=1, n_elemnz
              !do nn=csr_elemnz2nodenz(n-1)+1,csr_elemnz2nodenz(n) ! => 3
              nn=nn+1
              ! get the 3 nodes of the element at current level
              n_n1=csr_elemnz2nodenz_ID(nn)
              if(n_n1 == -1) then
                 U_rhs_lbs(n)=-bignumber
                 V_rhs_lbs(n)= bignumber
              else
                 n_n2=csr_elemnz2nodenz_ID(nn+2)
                 n_n3=csr_elemnz2nodenz_ID(nn+3)
                 nn = nn + 2

                 U_rhs_lbs(n)=max(fct_ttf_max_lbs(n_n1), fct_ttf_max_lbs(n_n2), fct_ttf_max_lbs(n_n3))
                 V_rhs_lbs(n)=min(fct_ttf_min_lbs(n_n1), fct_ttf_min_lbs(n_n2), fct_ttf_min_lbs(n_n3))
              endif
           end do
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "LBS on nodes/NZ iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           do elem=1, myDim_elem2D
              enodes=elem2D_nodes(:,elem)
              nu1 = ulevels(elem)
              nl1 = nlevels(elem)
              do nz=nu1, nl1-1
                 U_rhs(nz,elem)=maxval(fct_ttf_max_lbs(nodenz_lbs(nz,enodes)))
                 V_rhs(nz,elem)=minval(fct_ttf_min_lbs(nodenz_lbs(nz,enodes)))
              end do
              !$acc loop vector
              do nz=nl1,nl-1
                 U_rhs(nz,elem)=-bignumber
                 V_rhs(nz,elem)= bignumber
              end do
           end do ! --> do elem=1, myDim_elem2D
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') " timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        deallocate(ulevels)
        deallocate(nlevels)

end program oce_adv_tra_fct_loop_a2
