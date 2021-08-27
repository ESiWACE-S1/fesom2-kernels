program oce_adv_tra_fct_loop_a1
  use wallclock_mod
        implicit none

        integer :: n, nu1, nl1, nz, myDim_nod2D, n_it
        integer, dimension(:), allocatable :: ulevels_nod2D
        integer, dimension(:), allocatable :: nlevels_nod2D
        double precision, dimension(:,:), allocatable :: fct_ttf_min
        double precision, dimension(:,:), allocatable :: fct_ttf_max
        double precision, dimension(:,:), allocatable :: LO
        double precision, dimension(:,:), allocatable :: ttf
        integer, parameter :: MAX_PATH=1024
        integer, parameter :: MAX_LEVELS=50
        integer, parameter :: MAX_ITERATIONS=1000
        character(MAX_PATH) :: file_name
        integer :: mype, fileID, nn
        ! LBS_TRANSFORM
        integer, dimension(:), allocatable :: csr_node2levels ! 0:myDim_nod2D
        integer, dimension(:), allocatable :: csr_node2levels_ID ! level indices
        integer, dimension(:), allocatable :: csr_node2levels_nodID ! node indices
        integer :: nlevels_sum
        double precision, dimension(:), allocatable :: fct_ttf_min_lbs
        double precision, dimension(:), allocatable :: fct_ttf_max_lbs
        double precision, dimension(:), allocatable :: LO_lbs
        double precision, dimension(:), allocatable :: ttf_lbs

        !https://stackoverflow.com/a/6880672
        real(8)::t1,delta, delta_orig

        mype = 0
        myDim_nod2D = read_nod2D_levels(mype, ulevels_nod2D, nlevels_nod2D)

        allocate(fct_ttf_min(MAX_LEVELS, myDim_nod2D))
        allocate(fct_ttf_max(MAX_LEVELS, myDim_nod2D))
        allocate(LO(MAX_LEVELS, myDim_nod2D))
        allocate(ttf(MAX_LEVELS, myDim_nod2D))


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !$acc parallel loop gang present(fct_ttf_max,fct_ttf_min,LO,ttf,nlevels_nod2D,ulevels_nod2D)&
                !$acc& private(nl1,nu1)&
#ifdef WITH_ACC_VECTOR_LENGTH
                !$acc& vector_length(z_vector_length)&
#endif
                !$acc
                do n=1,myDim_nod2D
                        nu1 = ulevels_nod2D(n)
                        nl1 = nlevels_nod2D(n)
                        !$acc loop vector
                        do nz=nu1, nl1-1
                                fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
                                fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
                        end do
                end do
        end do
        delta=wallclock()-t1
        delta_orig=delta
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "LOOP_SEQ: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !$acc parallel loop gang vector present(fct_ttf_max,fct_ttf_min,LO,ttf,nlevels_nod2D,ulevels_nod2D)&
                !$acc& private(nl1,nu1)&
                !$acc
                do n=1,myDim_nod2D
                        nu1 = ulevels_nod2D(n)
                        nl1 = nlevels_nod2D(n)
                        !$acc loop seq
                        do nz=nu1, nl1-1
                                fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
                                fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
                        end do
                end do
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "LBS_TRANSFORM: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(csr_node2levels(0:myDim_nod2D))
        csr_node2levels(0) = 0
        do n=1,myDim_nod2D
           nu1 = ulevels_nod2D(n)
           nl1 = nlevels_nod2D(n)
           csr_node2levels(n) = csr_node2levels(n-1) + nl1-1-nu1+1
           nlevels_sum = nlevels_sum + nl1-1-nu1+1
        enddo
        write(*,*) nlevels_sum, csr_node2levels(myDim_nod2D)
        allocate(csr_node2levels_ID(nlevels_sum))
        allocate(csr_node2levels_nodID(nlevels_sum))
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
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
        !        !$acc parallel loop gang vector present(fct_ttf_max,fct_ttf_min,LO,ttf,nlevels_nod2D,ulevels_nod2D)&
        !        !$acc& private(nl1,nu1)&
        !        !$acc
           do nn=1,csr_node2levels(myDim_nod2D)
              n = csr_node2levels_nodID(nn)
              nz = csr_node2levels_ID(nn)
              !                nu1 = ulevels_nod2D(n)
              !                nl1 = nlevels_nod2D(n)
              !                !$acc loop seq
              !                do nz=nu1, nl1-1
              fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
              fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
              !                end do
           end do
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "LBS_TRANSFORM linear access: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(fct_ttf_min_lbs(nlevels_sum))
        allocate(fct_ttf_max_lbs(nlevels_sum))
        allocate(ttf_lbs(nlevels_sum))
        allocate(LO_lbs(nlevels_sum))
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           do nn=1,csr_node2levels(myDim_nod2D)
              fct_ttf_max_lbs(nn)=max(LO_lbs(nn), ttf_lbs(nn))
              fct_ttf_min_lbs(nn)=min(LO_lbs(nn), ttf_lbs(nn))
           end do
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta


        deallocate(ulevels_nod2D)
        deallocate(nlevels_nod2D)
end program oce_adv_tra_fct_loop_a1


