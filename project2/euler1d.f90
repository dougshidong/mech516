program euler1d

    use prec
    use glob
    use riemann_solver
    use eval
    use flux

    implicit none
    real(p2) :: a = -50.0d0, b = 50.0d0, dx, diaphragm = 0.0d0
    real(p2) :: CFL, t, tf = 25.0d0, max_uc, dt
    real(p2), allocatable :: x(:), U(:,:), W(:,:), c(:)
    real(p2) :: WL_ini(3), WR_ini(3) ! Initial Left and Right values
    integer  :: ix, i
    character(len=40) :: outn
    outn = 'sol.dat'

    CFL = 0.9
    nx = 100
    flux_scheme = 1
    ! Initialize variables give test case ID
    WL_ini(1) = 1.0d0
    WL_ini(2) = 0.0d0
    WL_ini(3) = 2.0d0
    WR_ini(1) = 1.0d0
    WR_ini(2) = 0.0d0
    WR_ini(3) = 1.0d0
    do flux_scheme = 1, 3
        call initialize
        t = 0.0d0
        do while(t.lt.tf)
            max_uc = maxval(abs(W(2,:) + c(:)))
            dt = CFL * dx / max_uc
            if((t + dt).gt.tf) then
                dt = tf - t
                t  = tf
            else
                t = t + dt
            end if
            write(*,*) 'Time ', t, ' of', tf
            call eval_flux(U, W, dt, dx)
            do ix = 1, nx
                U(:,ix) = U(:,ix) - dt / dx * (fluxes(:,ix+1) - fluxes(:,ix))
            end do
            call update_W_c
        end do
        call output
        deallocate(x, U, W, c, fluxes)
    end do

    contains
        subroutine initialize
            implicit none

            integer :: i

            select case(flux_scheme)
                case(1)
                    outn = 'Godunov.dat'
                case(2)
                    outn = 'HLLC.dat'
                case(3)
                    outn = 'MacCormack.dat'
                case default
            end select

            call initialize_flux
            dx = (b - a)/nx
            call allocate_ini
            ! Calculated at the middle point of each segment
            do i = 1, nx
                x(i) = a + dx * (i - 1) + dx / 2.0d0
                if(x(i).le.diaphragm) then
                    W(:,i) = WL_ini
                else
                    W(:,i) = WR_ini
                end if
            end do
            print *, x
            call update_U_c

            return
        end subroutine initialize

        subroutine allocate_ini
            allocate(x(nx), c(nx))
            allocate(U(3,nx),W(3,nx))
        end subroutine

        subroutine update_U_c
            do i = 1, nx
                U(:,i) = evalU(W(:,i), gam)
                c(i)   = evalc(W(:,i), gam)
            end do
        end subroutine

        subroutine update_W_c
            do i = 1, nx
                W(:,i) = evalW(U(:,i), gam)
                c(i)   = evalc(W(:,i), gam)
            end do
        end subroutine

        subroutine output
        implicit none
        open(unit=40,file=outn,form='FORMATTED')
        do ix = 1, nx
            write(40,'(4(E23.14))') x(ix), W(1,ix), W(2,ix), W(3,ix)
        end do
        close(40)
        end subroutine output

end program euler1d
