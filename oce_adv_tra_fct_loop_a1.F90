program oce_adv_tra_fct_loop_a1
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

        !https://stackoverflow.com/a/6880672
        real(8)::t1,delta
        iNTEGER :: clock_start,clock_end,clock_rate
        REAL(KIND=8) :: elapsed_time

        CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate

        mype = 0
        write(file_name, '(i8)') mype
        file_name='loops_nod2D_levels_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name)

        read(fileID, *) myDim_nod2D
        write(*,*) "loading",myDim_nod2D, " nodes"

        allocate(ulevels_nod2D(myDim_nod2D))
        allocate(nlevels_nod2D(myDim_nod2D))

        do n=1,myDim_nod2D !+ edim_nod2d
                ! TODO check nz<=MAX_LEVELS
                read(fileID, *) nn, nz, ulevels_nod2D(n), nlevels_nod2D(n)
                nlevels_nod2D(n) = nlevels_nod2D(n)+1
        end do

        close(fileID)

        allocate(fct_ttf_min(MAX_LEVELS, myDim_nod2D))
        allocate(fct_ttf_max(MAX_LEVELS, myDim_nod2D))
        allocate(LO(MAX_LEVELS, myDim_nod2D))
        allocate(ttf(MAX_LEVELS, myDim_nod2D))


        write(*,*) "iterating over",MAX_ITERATIONS, " iterations..."
        t1=secnds(0.0)
        CALL SYSTEM_CLOCK(COUNT=clock_start) ! Start timing
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
        delta=secnds(t1)
        CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing
        write(*,*) "done"
        ! Calculate the elapsed time in seconds:
        elapsed_time=REAL((clock_end-clock_start)/clock_rate,8)

        write(*,*) "timing", delta, delta/real(MAX_ITERATIONS)
        write(*,*) "elapsed", elapsed_time/real(MAX_ITERATIONS)

        deallocate(ulevels_nod2D)
        deallocate(nlevels_nod2D)
end program oce_adv_tra_fct_loop_a1


