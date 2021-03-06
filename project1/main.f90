program main

    use prec
    use riemann_solver
    use sample

    implicit none

    real(p2)              :: WL(3), WR(3), WS(3), W(3), rhosl, rhosr
    real(p2)              :: gam, t1, dx, S, x0
    integer               :: ix, nx, caseid
    real(p2), allocatable :: x(:), Wout(:,:)

    ! Primitive variables WL, WR, WS are the left, right and star vectors
    ! rhosl and rhosr are the densities on the left and right of contact surface
    caseid = 3
    nx = 1000
    dx = 1.0d0/nx
    allocate(x(nx))
    allocate(Wout(3,nx))

    ! Calculated at the middle point of each segment
    do ix = 1, nx
        x(ix) = -0.5d0 + dx * (ix - 0.5d0)
    end do
    ! Initialize variables give test case ID
    call initialize(caseid, gam, t1, x0, WL, WR)
    ! Solve for p*, u*, rho*left and rho*right
    call solve_riemann(WL, WR, gam, WS, rhosl, rhosr)

    do ix = 1, nx
        ! Evaluate x/t offsetted by the initial discontnuity location
        write(*,*) 'x(ix)', x(ix)
        S = (x(ix) - x0) / t1
        call sampling(WL, WR, gam, S, WS, rhosl, rhosr, W)
        Wout(:,ix) = W
    end do
    write(*,*) 'CS Location', WS(2) * t1 + x0

    open(unit=40,file='sol.dat',form='FORMATTED')
    write(40,'(4(A14))') 'x', 'rho', 'u', 'p'
    do ix = 1, nx
        write(40,'(4(E14.5))') x(ix), Wout(:,ix)
    end do
    close(40)

    contains
        subroutine initialize(caseid, gam, t1, x0, WL, WR)
            use prec
            implicit none
            integer, intent(in)   :: caseid
            real(p2), intent(out) :: gam, t1, x0, WL(3), WR(3)

            select case(caseid)
                case(1) ! Problem 1
                    gam   =  5.0d0/3.0d0
                    t1    =  0.04d0
                    x0    =  0.0d0
                    WL(1) =  8.0d0
                    WL(2) =  0.0d0
                    WL(3) =  480.0d0
                    WR(1) =  1.0d0
                    WR(2) =  0.0d0
                    WR(3) =  1.0d0
                case(2) ! Problem 2
                    gam   =  1.4d0
                    t1    =  0.15d0
                    x0    =  0.0d0
                    WL(1) =  1.0d0
                    WL(2) =  -2.0d0
                    WL(3) =  0.40d0
                    WR(1) =  1.0d0
                    WR(2) =  2.0d0
                    WR(3) =  0.40d0
                case(3) ! Problem 3
                    gam   =  1.4d0
                    t1    =  0.012d0
                    x0    =  0.30d0
                    WL(1) =  1.0d0
                    WL(2) =  -19.59745d0
                    WL(3) =  1000.0d0
                    WR(1) =  1.0d0
                    WR(2) =  -19.59745d0
                    WR(3) =  0.010d0
                case(4) ! Test 1
                    gam   =  1.4d0
                    t1    =  0.035d0
                    x0    =  -0.1d0
                    WL(1) =  5.99924d0
                    WL(2) =  19.5975d0
                    WL(3) =  460.894d0
                    WR(1) =  5.99242d0
                    WR(2) = -6.19633d0
                    WR(3) =  46.0950d0
                case(5) ! Test 2
                    gam   =  1.4d0
                    t1    =  0.2d0
                    x0    =  -0.2d0
                    WL(1) =  1.00000d0
                    WL(2) =  0.75000d0
                    WL(3) =  1.00000d0
                    WR(1) =  0.12500d0
                    WR(2) =  0.00000d0
                    WR(3) =  0.10000d0
                case(6) ! Test 3
                    gam   =  1.4d0
                    t1    =  0.035d0
                    x0    =  0.0d0
                    WL(1) =  1.00000d0
                    WL(2) =  0.00000d0
                    WL(3) =  0.01000d0
                    WR(1) =  1.00000d0
                    WR(2) =  0.00000d0
                    WR(3) =  100.000d0
                case default
                    write(*,*) 'Invalid test case ID'
            end select

            return
        end subroutine initialize


end program main
