program adv_tra_vert_impl
  use wallclock_mod
  use io_mod
  implicit none

  integer :: n, nu1, nl1, nz, myDim_nod2D, n_it
  integer, dimension(:), allocatable :: ulevels_nod2D
  integer, dimension(:), allocatable :: nlevels_nod2D
  integer                            :: nzmax, nzmin, tr_num
    double precision                   :: m, zinv, dt, dt_inv, dz
    double precision                   :: c1, v_adv, v_adv_p
    double precision, dimension(:,:), allocatable :: ttf, W
    double precision, dimension(:), allocatable :: a, b, c, tr
    double precision, dimension(:), allocatable :: cp, tp
  integer, parameter :: MAX_LEVELS=50
  integer, parameter :: MAX_ITERATIONS=200
  integer :: mype, nn

  double precision, dimension(:,:), allocatable :: areasvol
  double precision, dimension(:,:), allocatable :: area
  double precision, dimension(:,:), allocatable :: hnode_new
  real(8)::t1,delta, delta_orig

    dt=1e-3
    dt_inv=1.0d0/dt

  mype = 0
  myDim_nod2D = read_nod2D_levels(mype, ulevels_nod2D, nlevels_nod2D)

  allocate(ttf(MAX_LEVELS, myDim_nod2D))
  allocate(W  (MAX_LEVELS, myDim_nod2D))
  allocate(areasvol(MAX_LEVELS, myDim_nod2D))
  allocate(area(MAX_LEVELS, myDim_nod2D))
  allocate(hnode_new(MAX_LEVELS, myDim_nod2D))

  allocate(a (MAX_LEVELS))
  allocate(b (MAX_LEVELS))
  allocate(c (MAX_LEVELS))
  allocate(tr (MAX_LEVELS))
  allocate(cp (MAX_LEVELS))
  allocate(tp (MAX_LEVELS))

  nzmax=maxval(nlevels_nod2D(:))
    !___________________________________________________________________________
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  t1=wallclock()
  do n_it=1, MAX_ITERATIONS
    ! loop over local nodes
    !$acc parallel loop gang present(W,ttf,nlevels_nod2D,ulevels_nod2D,area,hnode_new,areasvol)&
    !$acc& default(none) &
    !$acc& private(nzmax,nzmin,nz,n,a,b,c,tr,cp,tp,zinv,dz,c1,v_adv)
    do n=1,myDim_nod2D
        ! max. number of levels at node n
        nzmax=nlevels_nod2D(n)
        ! upper surface index, in case of cavity !=1
        nzmin=ulevels_nod2D(n)

        !_______________________________________________________________________
        ! Regular part of coefficients: --> surface layer
        nz=nzmin

        ! 1/dz(nz)
        zinv=1.0d0*dt    ! no .../(zbar(1)-zbar(2)) because of  ALE

        a(nz)=0.0d0
        v_adv=zinv*area(nz  ,n)/areasvol(nz,n)
        b(nz)= hnode_new(nz,n)+W(nz, n)*v_adv

        v_adv=zinv*area(nz+1,n)/areasvol(nz,n)
        b(nz)= b(nz)-min(0.0d0, W(nz+1, n))*v_adv
        c(nz)=-max(0.0d0, W(nz+1, n))*v_adv

        !_______________________________________________________________________
        ! Regular part of coefficients: --> 2nd...nl-2 layer
        !$acc loop vector
        !private(v_adv)
        do nz=nzmin+1, nzmax-2
            ! update from the vertical advection
            v_adv=zinv*area(nz  ,n)/areasvol(nz,n)
            a(nz)=min(0.0d0, W(nz, n))*v_adv
            b(nz)=hnode_new(nz,n)+max(0.0d0, W(nz, n))*v_adv

            v_adv=zinv*area(nz+1,n)/areasvol(nz,n)
            b(nz)=b(nz)-min(0.0d0, W(nz+1, n))*v_adv
            c(nz)=     -max(0.0d0, W(nz+1, n))*v_adv
        end do ! --> do nz=2, nzmax-2

        !_______________________________________________________________________
        ! Regular part of coefficients: --> nl-1 layer
        nz=nzmax-1
        ! update from the vertical advection
        v_adv=zinv*area(nz  ,n)/areasvol(nz,n)
        a(nz)=                min(0.0d0, W(nz, n))*v_adv
        b(nz)=hnode_new(nz,n)+max(0.0d0, W(nz, n))*v_adv
        c(nz)=0.0d0

        !_______________________________________________________________________
        nz=nzmin
        dz=hnode_new(nz,n) ! It would be (zbar(nz)-zbar(nz+1)) if not ALE
        tr(nz)=-(b(nz)-dz)*ttf(nz,n)-c(nz)*ttf(nz+1,n)

        !$acc loop vector
        do nz=nzmin+1,nzmax-2
            tr(nz)=-a(nz)*ttf(nz-1,n)-(b(nz)-hnode_new(nz,n))*ttf(nz,n)-c(nz)*ttf(nz+1,n)
        end do
        nz=nzmax-1
        dz=hnode_new(nz,n)
        tr(nz)=-a(nz)*ttf(nz-1,n)-(b(nz)-dz)*ttf(nz,n)

        !_______________________________________________________________________
        nz = nzmin
        cp(nz) = c(nz)/b(nz)
        tp(nz) = tr(nz)/b(nz)

        ! solve for vectors c-prime and t, s-prime
        !ERROR $acc loop vector private(m)
        do nz = nzmin+1,nzmax-1
            m = b(nz)-cp(nz-1)*a(nz)
            cp(nz) = c(nz)/m
            tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
        end do

        !_______________________________________________________________________
        ! start with back substitution
        tr(nzmax-1) = tp(nzmax-1)

        ! solve for x from the vectors c-prime and d-prime
        !ERROR $acc loop vector
        do nz = nzmax-2, nzmin, -1
            tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
        end do

        !_______________________________________________________________________
        ! update tracer
        !$acc loop vector
        do nz=nzmin,nzmax-1
            ttf(nz,n)=ttf(nz,n)+tr(nz)
        end do
    end do ! --> do n=1,myDim_nod2D

  end do
  delta=wallclock()-t1
  delta_orig=delta
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta



end program adv_tra_vert_impl
