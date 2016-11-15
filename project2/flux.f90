module flux
use prec
use glob
use eval
use sample
use riemann_solver
integer :: nf, flux_scheme, MC = 1
real(p2), allocatable :: fluxes(:,:)

real(p2) :: WS(3), WL(3), WR(3), Wtemp(3), Utemp(3), S, rhosl, rhosr
real(p2) :: SL, SR, rl, ul, pl, rr, ur, pr, Sstar, pstar

contains

    subroutine eval_flux(U, W, dt, dx)
    implicit none
    integer :: i
    real(p2) :: U(3,nx), W(3,nx), dt, dx
    select case(flux_scheme)
        case(1) ! Godunov
            do i = 2, nf-1
                WL = W(:,i-1)
                WR = W(:,i)
                call solve_riemann(WL, WR, gam, WS, rhosl, rhosr)
                S = 0.0d0 ! (x(ix) - x0) / t1. Change if moving CS
                call sampling(WL, WR, gam, S, WS, rhosl, rhosr, Wtemp)
                Utemp = evalU(Wtemp, gam)
                fluxes(:,i) = evalF(Utemp, gam)
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
        case default
            write(*,*) 'Invalid flux scheme'
    end select
    ! Boundary fluxes
    fluxes(:,1) = fluxes(:,2)
    fluxes(:,nf) = fluxes(:,nf-1)
    return
    end subroutine eval_flux

    subroutine initialize_flux
        implicit none
        nf = nx + 1
        allocate(fluxes(3,nf))
        fluxes = 0
        return
    end subroutine initialize_flux

end module flux
