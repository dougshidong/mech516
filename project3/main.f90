program main

    use prec
    use glob
    use riemann_solver
    use sample

    implicit none

    real(p2)              :: WL(3), WR(3), WS(3), W(3), rhosl, rhosr
    real(p2)              :: t1, dx, S, x0
    real(p2), allocatable :: x(:), Wout(:,:)
    integer               :: ix, caseid, nx
    real(p2)              :: a, b, t0

    ! Primitive variables WL, WR, WS are the left, right and star vectors
    ! rhosl and rhosr are the densities on the left and right of contact surface
    nx = 100
    a  = -50.0d0
    b  = 50.0d0
    dx = (b - a)/nx
    allocate(x(nx))
    allocate(Wout(3,nx))

    caseid = 7
    t0 = 0.0d0

    ! Calculated at the middle point of each segment
    do ix = 1, nx
        x(ix) = a + dx * (ix - 1) + dx/2.0d0
    end do
    ! Initialize variables give test case ID
    call initialize2(caseid, t1, x0, WL, WR)
    ! Solve for p*, u*, rho*left and rho*right
    call solve_riemann(WL, WR, gam, WS, rhosl, rhosr)

    do ix = 1, nx
        ! Evaluate x/t offsetted by the initial discontnuity location
        S = (x(ix) - x0) / t1
        call sampling(WL, WR, gam, S, WS, rhosl, rhosr, W)
        Wout(:,ix) = W
    end do
    write(*,*) 'CS Location', WS(2) * t1 + x0

    open(unit=40,file='exact100.dat',form='FORMATTED')
    write(40,'(4(A14))') 'x', 'rho', 'u', 'p'
    do ix = 1, nx
        write(40,'(4(E23.14))') x(ix), Wout(:,ix)
    end do
    close(40)

    contains
        subroutine initialize2(caseid, t1, x0, WL, WR)
            use prec
            implicit none
            integer, intent(in)   :: caseid
            real(p2), intent(out) :: t1, x0, WL(3), WR(3)

            select case(caseid)
                case(7) ! Project 2
                    t1    =  25.00d0
                    x0    =  0.0d0
                    WL(1) =  1.00000d0
                    WL(2) =  0.00000d0
                    WL(3) =  2.00000d0
                    WR(1) =  1.00000d0
                    WR(2) =  0.00000d0
                    WR(3) =  1.00000d0
                case default
                    write(*,*) 'Invalid test case ID'
            end select

            return
        end subroutine initialize2


end program main
