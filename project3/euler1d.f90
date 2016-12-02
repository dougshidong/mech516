program euler1d

    use prec
    use glob
    use riemann_solver
    use eval
    use flux

    implicit none
    real(p2) :: a = 0.0d0, b = 1.50d0, dx, diaphragm = 0.0d0
    real(p2) :: CFL, t, tf = 0.01d0, max_uc, dt, ta, tb
    real(p2), allocatable :: x(:), U(:,:), W(:,:), c(:)
    real(p2) :: W1(3), W2(3) ! Initial Left and Right values
    integer  :: ix, iw, itime, nt
    character(len=40) :: outn
    outn = 'sol.dat'

    CFL = 0.5
    nw = 150
    nhalo = 2
    ntot = nw + 2*nhalo
    ! Initial conditions
    W1(1) = 1.0d0
    W1(2) = 0.0d0
    W1(3) = 1.0d0
    W2(1) = 100.0d0
    W2(2) = 0.0d0
    W2(3) = 1000.0d0
    flux_scheme = 4
    ilim = 5

    ix = 2
    select case(ix) ! Output solutions at different time intervals
        case(1)
            ta =-0.01d0
            tb = 0.15d0
        case(2)
            ta = 0.14d0
            tb = 0.30d0
        case(3)
            ta = 0.29d0
            tb = 0.45d0
    end select
    nt = 16 ! # of time intervals
    do itime = 1, nt
        tf = ta + (tb - ta) / nt * itime
        call initialize
        write(outn,'(A4,I1,A4)') 'time',itime,'.dat'
        if(itime.ge.10) write(outn,'(A4,I2,A4)') 'time',itime,'.dat'
        t = 0.0d0
        do while(t.lt.tf)
            if(t.eq.0.0d0) ilim = 1
            if(t.ne.0.0d0) ilim = 4
            max_uc = maxval(abs(W(2,:) + c(:)))
            dt = CFL * dx / max_uc
            if((t + dt).gt.tf) then
                dt = tf - t
                t  = tf
            else
                t = t + dt
            end if
            call eval_flux(U, W, dt, dx)
            !write(*,*) 'Time ', t, ' of', tf, fluxes(3,42)
            do ix = 1, nw
                iw = ix + nhalo
                U(:,iw) = U(:,iw) - dt / dx * (fluxes(:,iw+1) - fluxes(:,iw))
            end do
            call update_W_c(U, W, c, gam)
            call update_wall(W)
        end do
        open(unit=40,file=outn,form='FORMATTED')
        do ix = 1, nw
            iw = ix + nhalo
            write(40,'(5(E23.14))') x(iw), W(1,iw), W(2,iw), W(3,iw), tf
        end do
        close(40)
        !call output
        deallocate(x, U, W, c, fluxes, slopes)
    end do

    ! Pressure wrt to time computation for nw = 150, 300, 600, 1200, 2400, 4800
    nt = 6
    nw = 75
    CFL = 0.5
    do itime = 1, nt
        nw = nw * 2
        ntot = nw + 2*nhalo
        tf = 0.45d0
        call initialize
        write(outn,'(A4,I1,A4)') 'size',itime,'.dat'
        if(itime.ge.10) write(outn,'(A4,I2,A4)') 'size',itime,'.dat'
        open(unit=40,file=outn,form='FORMATTED')
        t = 0.0d0
        do while(t.lt.tf)
            if(t.eq.0.0d0) ilim = 1
            if(t.ne.0.0d0) ilim = 4
            max_uc = maxval(abs(W(2,:) + c(:)))
            dt = CFL * dx / max_uc
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
            do ix = 1, ntot
              if(x(ix)-1.30d0.ge.0.0d0) then
                write(*,'(3(E23.14))') t, W(3,ix), x(ix)
                write(40,'(3(E23.14))') t, W(3,ix), x(ix)
                exit
              end if
            end do
        end do
        close(40)
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
                if(x(i).ge.-0.1d0 .and. x(i).le.0.1d0) then
                    W(:,i) = W2
                else
                    W(:,i) = W1
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
        !do ix = 1, nw
        !    iw = ix + nhalo
        !    write(40,'(4(E23.14))') x(iw), W(1,iw), W(2,iw), W(3,iw)
        !end do
        do ix = 1, ntot
            iw = ix
            write(40,'(4(E23.14))') x(iw), W(1,iw), W(2,iw), W(3,iw)
        end do
        close(40)
        end subroutine output

end program euler1d
