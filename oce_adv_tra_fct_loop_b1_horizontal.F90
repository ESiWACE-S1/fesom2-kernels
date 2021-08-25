program oce_adv_tra_fct_loop_b1_horizontal
  use wallclock_mod
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
        integer, parameter :: MAX_PATH=1024
        integer, parameter :: MAX_LEVELS=50
        integer, parameter :: MAX_ITERATIONS=1000
        character(MAX_PATH) :: file_name
        integer :: mype, fileID, nn

        !https://stackoverflow.com/a/6880672
        real(8)::t1,delta

        mype = 0

        write(file_name, '(i8)') mype
        file_name='loops_edges_levels_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name, status='old')
        read(fileID, '(i8)') myDim_edge2D
        write(*,*) "loading",myDim_edge2D, " edges"

        allocate(ulevels(myDim_edge2D))
        allocate(nlevels(myDim_edge2D))

        do n=1, myDim_edge2D
                read(fileID, *) nn, nz, ulevels(n), nlevels(n)
        end do
        close(fileID)

        write(file_name, '(i8)') mype
        file_name='loops_edges_nodes_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name, status='old')
        read(fileID, *) myDim_edge2D

        allocate(edges(2,myDim_edge2D))

        do n=1, myDim_edge2D
                read(fileID, *) nn, edges(:,n) ! 2 nodes per edge
        end do
        close(fileID)


        write(file_name, '(i8)') mype
        file_name='loops_nod2D_levels_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name, status='old')

        read(fileID, *) myDim_nod2D
        close(fileID)

        write(file_name, '(i8)') mype
        file_name='loops_edges_edge_tri_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name, status='old')
        read(fileID, *) myDim_edge2D

        allocate(edge_tri(2,myDim_edge2D))

        do n=1, myDim_edge2D
                read(fileID, *) nn, edge_tri(:,n)
        end do
        close(fileID)

        write(file_name, '(i8)') mype
        file_name='loops_elem2D_levels_'//trim(adjustl(file_name))//'.dat'
        open(fileID, file=file_name, status='old')
        read(fileID, '(i8)') myDim_elem2D

        allocate(ulevels_elem2D(myDim_elem2D))
        allocate(nlevels_elem2D(myDim_elem2D))
        write(*,*) "loading",myDim_elem2D, " elements"

        do n=1, myDim_elem2D
                read(fileID, *) nn, nz, ulevels_elem2D(n), nlevels_elem2D(n)
                nlevels_elem2D(n) = nlevels_elem2D(n) + 1
        end do
        close(fileID)


        allocate(fct_plus(MAX_LEVELS, myDim_nod2D))
        allocate(fct_minus(MAX_LEVELS, myDim_nod2D))
        allocate(adf_h(MAX_LEVELS, myDim_edge2D))

        write(*,*) "iterating over",max_iterations, " iterations..."
        t1=wallclock()
        do n_it=1, max_iterations
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
                        el=edge_tri(:,n)
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
        write(*,*) "done"
        write(*,*) "timing", delta, delta/real(max_iterations)

        write(*,*) "iterating over",max_iterations, " iterations..."
        t1=wallclock()
        do n_it=1, max_iterations
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
        write(*,*) "timing", delta, delta/real(max_iterations)


end program oce_adv_tra_fct_loop_b1_horizontal
