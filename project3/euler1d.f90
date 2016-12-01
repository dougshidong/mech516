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
    integer  :: ix, iw
    character(len=40) :: outn
    outn = 'sol.dat'

    CFL = 1.0
    nw = 100
    nhalo = 2
    ntot = nw + 2*nhalo
    flux_scheme = 1
    ! Initialize variables give test case ID
    WL_ini(1) = 1.0d0
    WL_ini(2) = 0.0d0
    WL_ini(3) = 2.0d0
    WR_ini(1) = 1.0d0
    WR_ini(2) = 0.0d0
    WR_ini(3) = 1.0d0
    flux_scheme = 4
    do ilim = 1, 5
        call initialize
        t = 0.0d0
        do while(t.lt.tf)
            max_uc = maxval(abs(W(2,:) + c(:)))
            dt = CFL * dx / max_uc
            dt = 0.5
            if((t + dt).gt.tf) then
                dt = tf - t
                t  = tf
            else
                t = t + dt
            end if
            call eval_flux(U, W, dt, dx)
            write(*,*) 'Time ', t, ' of', tf, fluxes(3,42)
            do ix = 1, nw
                iw = ix + nhalo
                U(:,iw) = U(:,iw) - dt / dx * (fluxes(:,iw+1) - fluxes(:,iw))
            end do
            call update_W_c(U, W, c, gam)
            call update_wall(W)
        end do
        call output
        deallocate(x, U, W, c, fluxes, slopes)
    end do

    contains
        subroutine initialize
            implicit none

            integer :: i

            select case(ilim)
                case(1)
                    outn = 'ZERO.dat'
                case(2)
                    outn = 'AVG.dat'
                case(3)
                    outn = 'MINMOD.dat'
                case(4)
                    outn = 'SUPERBEE.dat'
                case(5)
                    outn = 'VANLEERNS.dat'
                case default
            end select

            call initialize_flux
            dx = (b - a)/nw
            call allocate_ini
            ! Calculated at the middle point of each segment
            do i = 1, nw + 2*nhalo
                x(i) = a + dx * (i-nhalo - 1) + dx / 2.0d0
                if(x(i).le.diaphragm) then
                    W(:,i) = WL_ini
                else
                    W(:,i) = WR_ini
                end if
            end do
            call update_U_c(U, W, c, gam)

            return
        end subroutine initialize

        subroutine allocate_ini
            allocate(x(ntot), c(ntot))
            allocate(U(3,ntot),W(3,ntot))
        end subroutine

        subroutine output
        implicit none
        open(unit=40,file=outn,form='FORMATTED')
        do ix = 1, nw
            iw = ix + nhalo
            write(40,'(4(E23.14))') x(iw), W(1,iw), W(2,iw), W(3,iw)
        end do
        close(40)
        end subroutine output

end program euler1d
