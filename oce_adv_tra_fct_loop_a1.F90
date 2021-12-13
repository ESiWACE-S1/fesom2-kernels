program oce_adv_tra_fct_loop_a1
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
        integer, parameter :: MAX_ITERATIONS=1000
        integer :: mype, nn
        ! LBS_TRANSFORM
        integer, dimension(:), allocatable :: csr_node2levels ! 0:myDim_nod2D
        integer, dimension(:), allocatable :: csr_node2levels_ID ! level indices
        integer, dimension(:), allocatable :: csr_node2levels_nodID ! node indices
        integer :: nlevels_sum
        double precision, dimension(:), allocatable :: fct_ttf_min_lbs
        double precision, dimension(:), allocatable :: fct_ttf_max_lbs
        double precision, dimension(:), allocatable :: LO_lbs
        double precision, dimension(:), allocatable :: ttf_lbs

        integer :: nl
        !https://stackoverflow.com/a/6880672
        real(8)::t1,delta, delta_orig
        integer :: number_nnz, number_nnzmax

        mype = 0
        myDim_nod2D = read_nod2D_levels(mype, ulevels_nod2D, nlevels_nod2D)

        allocate(fct_ttf_min(MAX_LEVELS, myDim_nod2D))
        allocate(fct_ttf_max(MAX_LEVELS, myDim_nod2D))
        allocate(LO(MAX_LEVELS, myDim_nod2D))
        allocate(ttf(MAX_LEVELS, myDim_nod2D))

        write(*,*) "COMPUTE: fct_ttf_max/min(nz,n) = max(LO(nz,n), ttf(nz,n)"
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
                !number_nnz = 0
                do n=1,myDim_nod2D
                        nu1 = ulevels_nod2D(n)
                        nl1 = nlevels_nod2D(n)
                        do nz=nu1, nl1-1
                                !number_nnz = number_nnz+1
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
        write(*,*) "LBS_TRANSFORM: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(csr_node2levels(0:myDim_nod2D))
        nlevels_sum = 0
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

        allocate(ttf_lbs(nlevels_sum))
        allocate(LO_lbs(nlevels_sum))

        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           do n=1,myDim_nod2D
              do nn=csr_node2levels(n-1)+1,csr_node2levels(n)
                 nz = csr_node2levels_ID(nn)
                 fct_ttf_max(nz,n)=max(LO_lbs(nn), ttf_lbs(nn))
                 fct_ttf_min(nz,n)=min(LO_lbs(nn), ttf_lbs(nn))
              end do
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

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "LBS_TRANSFORM linear access+PUSHBACK: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           do nn=1,csr_node2levels(myDim_nod2D)
              fct_ttf_max_lbs(nn)=max(LO_lbs(nn), ttf_lbs(nn))
              fct_ttf_min_lbs(nn)=min(LO_lbs(nn), ttf_lbs(nn))
           end do
           nn = 0
           do n=1,myDim_nod2D
              nu1 = ulevels_nod2D(n)
              nl1 = nlevels_nod2D(n)
              do nz=nu1, nl1-1
                 nn = nn + 1
                 fct_ttf_max(nz, n) = fct_ttf_max_lbs(nn)
                 fct_ttf_min(nz, n) = fct_ttf_min_lbs(nn)
              end do
           end do
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

#ifdef _OPENACC
#if 1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "ACC: iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !$acc enter data copyin(ulevels_nod2D, ulevels_nod2D, LO, ttf)
        !$acc enter data create(fct_ttf_max, fct_ttf_min)
        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           !$acc parallel loop present(ulevels_nod2D, ulevels_nod2D, LO, ttf, fct_ttf_max, fct_ttf_min)
           do n=1,myDim_nod2D
              nu1 = ulevels_nod2D(n)
              nl1 = nlevels_nod2D(n)
              do nz=nu1, nl1-1
                 fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
                 fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
              end do
           end do
        end do
        delta=wallclock()-t1
        !$acc exit data delete(ulevels_nod2D, ulevels_nod2D, LO, ttf)
        !$acc exit data copyout(fct_ttf_max, fct_ttf_min)
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
#endif
#if 1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) "ACC: collapse iterating over",MAX_ITERATIONS, " iterations..."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nl = maxval(nlevels_nod2D(:))

        t1=wallclock()
        do n_it=1, MAX_ITERATIONS
           !number_nnzmax = 0
           do n=1,myDim_nod2D
              do nz=1, nl
                 !number_nnzmax = number_nnzmax+1
                    fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
                    fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
              end do
           end do
           fct_ttf_min(1) = 0
           fct_ttf_max(1) = 0
        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
#endif
#endif

        write(*,*) 'nnz', number_nnz, number_nnzmax, real(number_nnzmax)/real(number_nnz)

        deallocate(ulevels_nod2D)
        deallocate(nlevels_nod2D)
end program oce_adv_tra_fct_loop_a1


