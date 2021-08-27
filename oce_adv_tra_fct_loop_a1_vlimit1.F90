program oce_adv_tra_fct_loop_a1_vlimit1
  use wallclock_mod
  use io_mod
        implicit none

        integer :: n, nu1, nl1, nz, myDim_nod2D, n_it, myDim_elem2D, elem, j
        integer, dimension(:), allocatable :: ulevels_nod2D
        integer, dimension(:), allocatable :: nlevels_nod2D
        integer, dimension(:), allocatable :: nod_in_elem2D_num
        integer, dimension(:,:), allocatable :: nod_in_elem2D
        double precision, dimension(:), allocatable :: tvert_min
        double precision, dimension(:), allocatable :: tvert_max
        double precision, dimension(:,:,:), allocatable :: UV_rhs
        double precision, dimension(:,:), allocatable :: ttf
        double precision, dimension(:,:), allocatable :: fct_ttf_min
        double precision, dimension(:,:), allocatable :: fct_ttf_max
        double precision, dimension(:,:), allocatable :: LO
        integer, parameter :: MAX_PATH=1024
        integer, parameter :: MAX_LEVELS=50
        integer, parameter :: MAX_NOD_IN_ELEM=10
        integer, parameter :: MAX_ITERATIONS=1000
        character(MAX_PATH) :: file_name, line
        integer :: mype, fileID, nn
        ! LBS_TRANSFORM
        integer, dimension(:), allocatable :: csr_node2levels ! 0:myDim_nod2D
        integer, dimension(:), allocatable :: csr_node2levels_ID ! level indices
        integer, dimension(:), allocatable :: csr_node2levels_nodID ! node indices
        integer :: nlevels_sum
        double precision, dimension(:), allocatable :: tvert_min_lbs
        double precision, dimension(:), allocatable :: tvert_max_lbs
        double precision, dimension(:), allocatable :: fct_ttf_min_lbs
        double precision, dimension(:), allocatable :: fct_ttf_max_lbs
        double precision, dimension(:), allocatable :: LO_lbs
        double precision, dimension(:), allocatable :: ttf_lbs
        double precision, dimension(:,:,:), allocatable :: UV_rhs_


        !https://stackoverflow.com/a/6880672
        real(8)::t1,delta, delta_orig

        mype = 0
        myDim_nod2D = read_nod2D_levels(mype, ulevels_nod2D, nlevels_nod2D)

        write(file_name, '(i8)') mype
        file_name='loops_nod2D_elems_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name)
        read(fileID, '(i8)') myDim_nod2D

        allocate(nod_in_elem2D_num(myDim_nod2D))
        allocate(nod_in_elem2D(MAX_NOD_IN_ELEM, myDim_nod2D))

        do n=1,myDim_nod2D !+ edim_nod2d
                !read(fileID,*) nn, nod_in_elem2D_num(n), nod_in_elem2D(1:mesh%nod_in_elem2D_num(n),n)
                read(fileID, '(a)') line
                read(line,*) nn, nod_in_elem2D_num(n)
                read(line,*) nn, nod_in_elem2D_num(n), nod_in_elem2D(1:nod_in_elem2D_num(n), n)
        end do
        close(fileID)

        write(file_name, '(i8)') mype
        file_name='loops_elem2D_nodes_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name)
        read(fileID, '(i8)') myDim_elem2D
        close(fileID)

        allocate(tvert_min(MAX_LEVELS))
        allocate(tvert_max(MAX_LEVELS))
        allocate(UV_rhs(2,MAX_LEVELS, myDim_elem2D))
        allocate(fct_ttf_min(MAX_LEVELS, myDim_nod2D))
        allocate(fct_ttf_max(MAX_LEVELS, myDim_nod2D))
        allocate(LO(MAX_LEVELS, myDim_nod2D))


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !$acc parallel loop gang present(ulevels_nod2D,nlevels_nod2D,nod_in_elem2D,nod_in_elem2D_num,UV_rhs&
                !$acc& ) private(tvert_min,tvert_max)&
                !$acc& private(nu1,nl1)&
#ifdef WITH_ACC_VECTOR_LENGTH
                !$acc& vector_length(z_vector_length)&
#endif
                !$acc
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

        end do
        delta=wallclock()-t1
        delta_orig=delta
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "3 LOOPS: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !$acc parallel loop gang present(ulevels_nod2D,nlevels_nod2D,nod_in_elem2D,nod_in_elem2D_num,UV_rhs&
                !$acc& ) private(tvert_min,tvert_max)&
                !$acc& private(nu1,nl1)&
#ifdef WITH_ACC_VECTOR_LENGTH
                !$acc& vector_length(z_vector_length)&
#endif
                !$acc
                do n=1, myDim_nod2D
                    nu1 = ulevels_nod2D(n)
                    nl1 = nlevels_nod2D(n)
                    do nz=nu1,nl1-1
                       tvert_max(nz) = huge(1.0)
                       tvert_min(nz) = -huge(1.0)
                    end do
                    !$acc loop vector
                    do nn=1, nod_in_elem2D_num(n)
                       elem = nod_in_elem2D(nn,n)
                       do nz=nu1,nl1-1
                          tvert_max(nz)= max(tvert_max(nz), UV_rhs(1,nz,elem))
                          tvert_min(nz)= min(tvert_min(nz), UV_rhs(2,nz,elem))
                       end do
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

        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "3 LOOPS+max(nz-1,nz,nz+1): iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !$acc parallel loop gang present(ulevels_nod2D,nlevels_nod2D,nod_in_elem2D,nod_in_elem2D_num,UV_rhs&
                !$acc& ) private(tvert_min,tvert_max)&
                !$acc& private(nu1,nl1)&
#ifdef WITH_ACC_VECTOR_LENGTH
                !$acc& vector_length(z_vector_length)&
#endif
                !$acc
                do n=1, myDim_nod2D
                    nu1 = ulevels_nod2D(n)
                    nl1 = nlevels_nod2D(n)
                    do nz=nu1,nl1-1
                       tvert_max(nz) = huge(1.0)
                       tvert_min(nz) = -huge(1.0)
                    end do
                    !$acc loop vector
                    do nn=1, nod_in_elem2D_num(n)
                       elem = nod_in_elem2D(nn,n)
                       do nz=nu1,nl1-1
                          tvert_max(nz)= max(tvert_max(nz), UV_rhs(1,nz,elem))
                          tvert_min(nz)= min(tvert_min(nz), UV_rhs(2,nz,elem))
                       end do
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
                        fct_ttf_max(nz,n)=max(tvert_max(nz-1),tvert_max(nz), tvert_max(nz+1))-LO(nz,n)
                        fct_ttf_min(nz,n)=min(tvert_min(nz-1),tvert_min(nz), tvert_min(nz+1))-LO(nz,n)
                    end do
                    ! calc max,min increment of bottom layer -1 with respect to low order
                    ! solution
                    nz=nl1-1
                    fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
                    fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)
                end do

        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(csr_node2levels(0:myDim_nod2D))
        csr_node2levels(0) = 0
        nlevels_sum=0
        do n=1,myDim_nod2D
           nu1 = ulevels_nod2D(n)
           nl1 = nlevels_nod2D(n)
           csr_node2levels(n) = csr_node2levels(n-1) + nl1-1-nu1+1
           nlevels_sum = nlevels_sum + nl1-1-nu1+1
        enddo
        write(*,*) nlevels_sum, csr_node2levels(myDim_nod2D)
        allocate(csr_node2levels_ID(0:nlevels_sum+1))
        allocate(csr_node2levels_nodID(0:nlevels_sum+1))
        csr_node2levels_ID(0) = -1
        csr_node2levels_ID(nlevels_sum+1) = -1
        csr_node2levels_nodID(0) = -1
        csr_node2levels_nodID(nlevels_sum+1) = -1
        nn = 0
        do n=1,myDim_nod2D
           nu1 = ulevels_nod2D(n)
           nl1 = nlevels_nod2D(n)
           do nz=nu1, nl1-1
              nn = nn + 1
              csr_node2levels_ID(nn) = nz
              csr_node2levels_nodID(nn) = n
           end do
        end do
        write(*,*) nlevels_sum, csr_node2levels(myDim_nod2D), nn
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "LBS_TRANSFORM linear access: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(tvert_min_lbs(nlevels_sum))
        allocate(tvert_max_lbs(nlevels_sum))
        allocate(fct_ttf_min_lbs(nlevels_sum))
        allocate(fct_ttf_max_lbs(nlevels_sum))
        allocate(ttf_lbs(nlevels_sum))
        allocate(LO_lbs(nlevels_sum))
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           do nn=1,csr_node2levels(myDim_nod2D)
              tvert_max_lbs(nn) = huge(1.0)
              tvert_min_lbs(nn) = -huge(1.0)
           end do
           do nn=1,csr_node2levels(myDim_nod2D)
              nz = csr_node2levels_ID(nn)
              n = csr_node2levels_nodID(nn)
              do j=1, nod_in_elem2D_num(n)
                 elem = nod_in_elem2D(j,n)
                 tvert_max_lbs(nn)= max(tvert_max_lbs(nn), UV_rhs(1,nz,elem))
                 tvert_min_lbs(nn)= min(tvert_min_lbs(nn), UV_rhs(2,nz,elem))
              end do
           end do
           !!!___________________________________________________________________
           !!! calc max,min increment of surface layer with respect to low order
           !!! solution
           do nn=1,csr_node2levels(myDim_nod2D)
              !!! calc max,min increment from nz-1:nz+1 with respect to low order
              !!! solution at layer nz
              if(csr_node2levels_ID(nn-1) == csr_node2levels_ID(nn+1)) then
                 fct_ttf_max_lbs(nn) = max(tvert_max_lbs(nn-1),tvert_max_lbs(nn),tvert_max_lbs(nn+1)) - LO_lbs(nn)
                 fct_ttf_min_lbs(nn) = min(tvert_min_lbs(nn-1),tvert_min_lbs(nn),tvert_min_lbs(nn+1)) - LO_lbs(nn)
              else
                 fct_ttf_max_lbs(nn) = tvert_max_lbs(nn) - LO_lbs(nn)
                 fct_ttf_min_lbs(nn) = tvert_min_lbs(nn) - LO_lbs(nn)
              end if
           end do
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "LBS_TRANSFORM linear access+permute nz/elem: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(UV_rhs_(2, myDim_nod2D, MAX_LEVELS))
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           do nn=1,csr_node2levels(myDim_nod2D)
              tvert_max_lbs(nn) = huge(1.0)
              tvert_min_lbs(nn) = -huge(1.0)
           end do
           do n=1, myDim_nod2D
              nn = csr_node2levels(n-1)
              do j=1, nod_in_elem2D_num(n)
                 elem = nod_in_elem2D(j,n)
                 !do nn=csr_node2levels(n-1)+1, csr_node2levels(n)
                 !   nz = csr_node2levels_ID(nn)
                 ! TODO nz>1 when cavity
                 !do nz=1,csr_node2levels(n)-csr_node2levels(n-1)
                 do nz=csr_node2levels_ID(nn+1), csr_node2levels_ID(csr_node2levels(n))
                    tvert_max_lbs(nn+nz)= max(tvert_max_lbs(nn+nz), UV_rhs(1,nz,elem))
                    tvert_min_lbs(nn+nz)= min(tvert_min_lbs(nn+nz), UV_rhs(2,nz,elem))
                 end do
              end do
           end do
           !!!___________________________________________________________________
           !!! calc max,min increment of surface layer with respect to low order
           !!! solution
           do nn=1,csr_node2levels(myDim_nod2D)
              !!! calc max,min increment from nz-1:nz+1 with respect to low order
              !!! solution at layer nz
              if(csr_node2levels_ID(nn-1) == csr_node2levels_ID(nn+1)) then
                 fct_ttf_max_lbs(nn) = max(tvert_max_lbs(nn-1),tvert_max_lbs(nn),tvert_max_lbs(nn+1)) - LO_lbs(nn)
                 fct_ttf_min_lbs(nn) = min(tvert_min_lbs(nn-1),tvert_min_lbs(nn),tvert_min_lbs(nn+1)) - LO_lbs(nn)
              else
                 fct_ttf_max_lbs(nn) = tvert_max_lbs(nn) - LO_lbs(nn)
                 fct_ttf_min_lbs(nn) = tvert_min_lbs(nn) - LO_lbs(nn)
              end if
           end do
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta



        deallocate(ulevels_nod2D)
        deallocate(nlevels_nod2D)
end program oce_adv_tra_fct_loop_a1_vlimit1


