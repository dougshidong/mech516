module flux
use prec
use glob
use eval
use sample
use riemann_solver
integer :: nf, flux_scheme, ilim, MC = 1
real(p2), allocatable :: fluxes(:,:)
real(p2), allocatable :: slopes(:,:), Wvtemp(:,:), Uvtemp(:,:)

real(p2) :: WS(3), WL(3), WR(3), Wtemp(3), Utemp(3), S, rhosl, rhosr
real(p2) :: SL, SR, rl, ul, pl, rr, ur, pr, Sstar, pstar

contains

    subroutine eval_flux(U, W, dt, dx)
    implicit none
    integer :: i, ix, iw
    real(p2) :: U(3,ntot), W(3,ntot), dt, dx
    real(p2) :: dWL(3), dWR(3), k
    select case(flux_scheme)
        case(1) ! Godunov
            do i = 2, nf-1
                WL = W(:,i-1)
                WR = W(:,i)
                call solve_riemann(WL, WR, gam, WS, rhosl, rhosr)
                S = 0.0d0 ! (x(ix) - x0) / t1. Change if moving CS
                call sampling(WL, WR, gam, S, WS, rhosl, rhosr, Wtemp)
                fluxes(:,i) = evalFW(Wtemp, gam)
            end do
        case(2) ! HLLC
            do i = 2, nf-1
                SL = W(2,i-1) - evalc(W(:, i-1), gam)
                SR = W(2,i) + evalc(W(:,i), gam)
                S = 0.0d0 ! (x(ix) - x0) / t1. Change if moving CS
                if(SL.ge.S) then ! Left of left wave
                    fluxes(:,i) = evalF(U(:,i-1), gam)
                else if(SR.le.S) then ! Right of right wave
                    fluxes(:,i) = evalF(U(:,i), gam)
                else
                    rl = W(1,i-1)
                    ul = W(2,i-1)
                    pl = W(3,i-1)
                    rr = W(1,i)
                    ur = W(2,i)
                    pr = W(3,i)
                    Sstar = (pr-pl + rl*ul*(SL-ul) - rr*ur*(SR-ur))
                    Sstar = Sstar / (rl*(SL-ul) - rr*(SR-ur))
                    if(Sstar.ge.S) then ! Left of CS
                        pstar = pl + rl*(SL-ul)*(Sstar-ul)
                        Utemp = rl * (SL - ul) / (SL - Sstar)
                        !Utemp(1) = Utemp(1) * 1
                        Utemp(2) = Utemp(2) * Sstar
                        Utemp(3) = Utemp(3) * (U(3,i-1) / U(1,i-1) &
                            + (Sstar-ul) * (Sstar + pl / rl / (SL - ul)))
                        fluxes(:,i) = evalF(U(:,i-1), gam) + SL*(Utemp - U(:,i-1))
                    else ! Right of CS
                        pstar = pr + rr*(SR-ur)*(Sstar-ur)
                        Utemp = rr * (SR - ur) / (SR - Sstar)
                        !Utemp(1) = Utemp(1) * 1
                        Utemp(2) = Utemp(2) * Sstar
                        Utemp(3) = Utemp(3) * (U(3,i) / U(1,i) &
                            + (Sstar-ur) * (Sstar + pr / rr / (SR - ur)))
                        fluxes(:,i) = evalF(U(:,i), gam) + SR*(Utemp - U(:,i))
                    end if
                end if
            end do
        case(3) ! MacCormack
            if(MC.eq.1) then ! Alternate between MC1 and MC2
                MC = 2
                do i = 2, nf-1
                    fluxes(:,i) = evalF(U(:,i), gam) + &
                                  evalF(  U(:,i-1) - (dt / dx) &
                                          * ( evalF(U(:,i),gam)-evalF(U(:,i-1),gam) ) &
                                        , gam  )
                end do
            else
                MC = 1
                do i = 2, nf-1
                    fluxes(:,i) = evalF(U(:,i-1), gam) + &
                                  evalF(  U(:,i) - (dt / dx) &
                                          * ( evalF(U(:,i),gam)-evalF(U(:,i-1),gam) ) &
                                        , gam  )
                end do
            end if
            fluxes = fluxes * 0.5d0
        case(4) ! MUSCL-Hancock
            allocate(Wvtemp(3,ntot),Uvtemp(3,ntot))
            Wvtemp = W
            Uvtemp = U
            k = 0.0
            do i = 2, ntot-1
                dWL(:) = W(:,i) - W(:,i-1)
                dWR(:) = W(:,i+1) - W(:,i)
                call evalSlope(slopes(:,i), dWL, dWR, k)
            end do

            dWL(:) = 0.0d0
            dWR(:) = W(:,2) - W(:,1)
            call evalSlope(slopes(:,1), dWL, dWR, k)

            dWL(:) = W(:,ntot) - W(:,ntot-1)
            dWR(:) = 0.0d0
            call evalSlope(slopes(:,ntot), dWL, dWR, k)

            do ix = 1, nw
                iw = ix + nhalo
                Uvtemp(:,iw) = Uvtemp(:,iw) - dt / dx * &
                               (evalFW(W(:,iw) + 0.5d0*slopes(:,iw), gam) &
                              - evalFW(W(:,iw) - 0.5d0*slopes(:,iw), gam))
            end do

            call update_W(Uvtemp, Wvtemp, gam)
            call update_wall(Wvtemp)

            do i = 2, nf-1
                WL = W(:,i-1) + Wvtemp(:,i-1) + slopes(:,i-1)
                WR = W(:,i) + Wvtemp(:,i) - slopes(:,i)
                WL = 0.5d0*WL
                WR = 0.5d0*WR
                call solve_riemann(WL, WR, gam, WS, rhosl, rhosr)
                S = 0.0d0 ! (x(ix) - x0) / t1. Change if moving CS
                call sampling(WL, WR, gam, S, WS, rhosl, rhosr, Wtemp)
                fluxes(:,i) = evalFW(Wtemp, gam)
            end do
            deallocate(Wvtemp,Uvtemp)

        case default
            write(*,*) 'Invalid flux scheme'
    end select
    ! Boundary fluxes
    !fluxes(:,1) = fluxes(:,2)
    !fluxes(:,nf) = fluxes(:,nf-1)
    return
    end subroutine eval_flux

    subroutine initialize_flux
        implicit none
        nf = ntot + 1
        allocate(fluxes(3,nf))
        allocate(slopes(3,ntot))
        fluxes = 0
        return
    end subroutine initialize_flux

    subroutine evalSlope(slope, dWL, dWR, k)
        implicit none
        real(p2) :: slope(3), dWL(3), dWR(3), k, one = 1.0d0, beta
        integer  :: i
        ! k = -1  -  right slope
        ! k =  0  -  average slope
        ! k =  1  -  left slope
        select case(ilim) ! Limiter Type
            case(1) ! 0 slope
                slope = 0.0
            case(2) ! 0 slope
                slope = (1.0d0+k)*dWL + (1.0d0-k)*dWR
            case(3) ! MINMOD
                do i = 1, 3
                    slope(i) = sign(one, dWL(i)) * &
                               max( 0.0, min( abs(dWL(i)), sign(one,dWL(i))*dWR(i) ) )
                end do
            case(4) ! SUPERBEE
                beta = 2
                do i = 1, 3
                    slope(i) = 0.5d0*( sign(one, dWL(i)) + sign(one, dWR(i)) ) * &
                               max( min( beta * abs(dWL(i)), abs(dWR(i)) ), &
                                    min( abs(dWL(i)), beta * abs(dWR(i)) ) )
                end do
            case(5) ! Van Leer
                do i = 1, 3
                    slope(i) = 0.5d0*( sign(one, dWL(i)) + sign(one, dWR(i)) ) * &
                               2.0d0*min( abs(dWL(i)), &
                                          abs(dWR(i)), &
                                          abs(dWL(i) + dWR(i)) * 0.25d0 )
                end do
            case default
                write(*,*) 'Invalid limiter'
        end select
        slope = slope * 0.5d0

    end subroutine

end module flux
