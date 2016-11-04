program linear_advection

implicit none

integer :: ix, nx, nt, it, flux_scheme, caseid
double precision, allocatable :: x(:), u(:), f(:)
double precision, parameter :: xi = -10.0d0, xf = 240.0d0, tf = 100.0d0, a = 2.0d0
double precision, parameter :: pi = 4.d0*atan(1.d0)
double precision :: CFL, dt, dx
integer :: ilist(4),j,k
double precision :: dlist(3)
character(len=40) :: outn

CFL = 0.8d0   ! CFL Number
nx  = 1000
flux_scheme = 2
caseid = 2

! Exercise 1
dlist = (/ 0.8, 1.02, 2.0/) ! List of CFL numbers
do j = 1,3
    do k = 1, 5 ! 
        CFL = dlist(j)
        write(outn,'(a,i1,a,i3,a)') 'p1cfl',j-1,'t',k-1,'.dat'
        call initialize
        do it = 1, k
            if(it.eq.nt) dt = tf - (nt-1)*dt
            call eval_flux
            do ix = 1, nx
                u(ix) = u(ix) - dt / dx * (f(ix+1) - f(ix))
            end do
        end do
        call output
        deallocate(x,u,f)
    end do
end do

! Exercise 2 and 3
CFL = 0.8
ilist = (/ 1000, 2000, 4000, 8000 /)
do flux_scheme = 1, 2
    do caseid = 1, 2
        do j = 1, 4
            nx = ilist(j)
            write(outn,'(a,i1,a,i1,a,i4.4,a)') 'p2case',caseid,'f',flux_scheme,'nx',nx,'.dat'
            write(*,*) outn
            call initialize
            do it = 1, nt
                if(it.eq.nt) dt = tf - (nt-1)*dt
                call eval_flux
                do ix = 1, nx
                    u(ix) = u(ix) - dt / dx * (f(ix+1) - f(ix))
                end do
            end do
            call output
            deallocate(x,u,f)
        end do
    end do
end do


contains
    subroutine initialize
    implicit none
    integer :: i
    allocate(x(nx), u(nx), f(nx+1))

    dx = (xf - xi) / nx
    do i = 1, nx
        x(i) = xi + dx/2.0d0 + (i - 1) * dx
    end do

    u = 0.5d0
    do i = 1, nx
        if(x(i).ge.0.0d0 .and. x(i).le.20.0d0) then
            if(caseid.eq.1) then
                u(i) = 0.5d0 + 0.075d0*x(i)
            else if(caseid.eq.2) then
                u(i) = 1.0d0 - 0.5d0*cos(pi*x(i)/10.0d0)
            end if
        end if
    end do

    dt = CFL * dx / a ! Time-step required to respect CFL
    nt = floor(tf/dt) + 1
    end subroutine initialize

    subroutine eval_flux
    implicit none
    integer :: i
    double precision :: ap, an
    select case(flux_scheme)
        case(1) ! CIR
            ap = 0.5d0 * (a + abs(a))
            an = 0.5d0 * (a - abs(a))
            do i = 2, nx
                f(i) = ap * u(i-1) + an * u(i)
            end do
            f(1)  = a * u(1)
            f(nx + 1) = a * u(nx)
        case(2) ! Lax-Wendroff
            do i = 2, nx
                f(i) = a * 0.5d0 * ((u(i) + u(i-1)) - a * dt/dx * (u(i) - u(i-1)))
            end do
            f(1)  = a * u(1)
            f(nx + 1) = a * u(nx)
        case default
            write(*,*) 'Invalid flux scheme'
    end select

    end subroutine eval_flux

    subroutine output
    implicit none
    open(unit=40,file=outn,form='FORMATTED')
    do ix = 1, nx
        write(40,'(2(E23.14))') x(ix), u(ix)
    end do
    close(40)
    end subroutine output

end program linear_advection



