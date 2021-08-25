program oce_adv_tra_fct_loop_a1_vlimit1
  use wallclock_mod
        implicit none

        integer :: n, nu1, nl1, nz, myDim_nod2D, n_it, myDim_elem2D
        integer, dimension(:), allocatable :: ulevels_nod2D
        integer, dimension(:), allocatable :: nlevels_nod2D
        integer, dimension(:), allocatable :: nod_in_elem2D_num
        integer, dimension(:,:), allocatable :: nod_in_elem2D
        double precision, dimension(:), allocatable :: tvert_min
        double precision, dimension(:), allocatable :: tvert_max
        double precision, dimension(:,:,:), allocatable :: UV_rhs
        double precision, dimension(:,:), allocatable :: ttf
        integer, parameter :: MAX_PATH=1024
        integer, parameter :: MAX_LEVELS=50
        integer, parameter :: MAX_NOD_IN_ELEM=10
        integer, parameter :: MAX_ITERATIONS=1000
        character(MAX_PATH) :: file_name, line
        integer :: mype, fileID, nn

        !https://stackoverflow.com/a/6880672
        real(8)::t1,delta

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


        write(*,*) "iterating over",MAX_ITERATIONS, " iterations..."
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
                    !!!___________________________________________________________________
                    !!! calc max,min increment of surface layer with respect to low order
                    !!! solution
                    !!fct_ttf_max(nu1,n)=tvert_max(nu1)-LO(nu1,n)
                    !!fct_ttf_min(nu1,n)=tvert_min(nu1)-LO(nu1,n)

                    !!! calc max,min increment from nz-1:nz+1 with respect to low order
                    !!! solution at layer nz
                    !!!$acc loop vector
                    !!do nz=nu1+1,nl1-2
                    !!    fct_ttf_max(nz,n)=maxval(tvert_max(nz-1:nz+1))-LO(nz,n)
                    !!    fct_ttf_min(nz,n)=minval(tvert_min(nz-1:nz+1))-LO(nz,n)
                    !!end do
                    !!! calc max,min increment of bottom layer -1 with respect to low order
                    !!! solution
                    !!nz=nl1-1
                    !!fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
                    !!fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)
                end do

        end do
        delta=wallclock()-t1
        write(*,*) "done"
        write(*,*) "timing", delta, delta/real(MAX_ITERATIONS)

        deallocate(ulevels_nod2D)
        deallocate(nlevels_nod2D)
end program oce_adv_tra_fct_loop_a1_vlimit1


