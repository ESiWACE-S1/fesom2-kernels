program oce_adv_tra_fct_loop_a2
  use wallclock_mod
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
        integer, parameter :: MAX_PATH=1024
        integer, parameter :: MAX_LEVELS=50
        integer, parameter :: MAX_ITERATIONS=1000
        character(MAX_PATH) :: file_name
        integer :: mype, fileID, nn

        !https://stackoverflow.com/a/6880672
        real(8)::t1,delta, delta_orig

        mype = 0
        write(file_name, '(i8)') mype
        file_name='loops_elem2D_nodes_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name)

        read(fileID, *) myDim_elem2D
        write(*,*) "loading",myDim_elem2D, " elements"

        allocate(elem2D_nodes(3, myDim_elem2D))

        do n=1,myDim_elem2D !+ edim_nod2d
                ! TODO check nz<=MAX_LEVELS
                read(fileID, *) nn, elem2D_nodes(:, n)
        end do

        close(fileID)

        write(file_name, '(i8)') mype
        file_name='loops_elem2D_levels_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name)
        read(fileID, '(i8)') myDim_elem2D

        allocate(ulevels(myDim_elem2D))
        allocate(nlevels(myDim_elem2D))

        do n=1, myDim_elem2D
                read(fileID, *) nn, nz, ulevels(n), nlevels(n)
                nlevels(n) = nlevels(n) + 1
        end do
        close(fileID)


        allocate(fct_ttf_min(MAX_LEVELS, myDim_elem2D))
        allocate(fct_ttf_max(MAX_LEVELS, myDim_elem2D))
        allocate(UV_rhs(2,MAX_LEVELS, myDim_elem2D))

        nl = maxval(nlevels(:)-ulevels(:)+1)

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
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

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
        write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

        deallocate(ulevels)
        deallocate(nlevels)

end program oce_adv_tra_fct_loop_a2
