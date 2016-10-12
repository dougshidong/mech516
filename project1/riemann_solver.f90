module riemann_solver
use prec
use eval

contains

subroutine solve_riemann(WL, WR, gam, WS, rhosl, rhosr)
    implicit none
    real(p2), intent(in)  :: WL(3), WR(3)
    real(p2), intent(in)  :: gam
    real(p2), intent(out) :: WS(3), rhosl, rhosr
    real(p2)              :: prs
    real(p2)              :: f, df, pnext, fl, fr
    real(p2)              :: resi
    real(p2), parameter   :: eps = 1e-12

    ! Initialize p* with acoustic approximation
    pnext = acoustic_approx(WL, WR, gam)
    resi  = 1
    ! Solve for p* iteratively using Newton-Raphson
    do while(resi.gt.eps)
        prs = pnext
        f   = flfrdu(prs, WL, WR, gam)
        df  = dflfrdu(prs, WL, WR, gam)
        pnext = prs - f/df
        resi = abs(pnext - prs) / prs

        write(*,'(2(A,E23.15))') 'prs', prs, 'resi', resi
    end do

    prs = pnext
    WS(3) = prs

    ! Solve for u*, rhosl, and rhosr given p*

    ! Left
    if(prs > WL(3)) then    ! Left shock
        fl =       fshock(prs, WL, gam)
        rhosl = evalrho_hugoniot(WL(1), prs/WL(3), gam)
    else                    ! Left rarefaction
        fl = frarefaction(prs, WL, gam)
        rhosl = evalrho_isentropic(WL(1), prs/WL(3), gam)
    end if
    ! Right
    if(prs > WR(3)) then    ! Right shock
        fr =       fshock(prs, WR, gam)
        rhosr = evalrho_hugoniot(WR(1), prs/WR(3), gam)
    else                    ! Right rarefaction
        fr = frarefaction(prs, WR, gam)
        rhosr = evalrho_isentropic(WR(1), prs/WR(3), gam)
    end if

    WS(2) = 0.5d0*(WL(2) + WR(2) + fr - fl)

    return
end subroutine solve_riemann
!***************************************************************    
function acoustic_approx(WL, WR, gam) result(p0)
    implicit none
    real(p2), intent(in)  :: WL(3), WR(3), gam
    real(p2)              :: p0
    real(p2)              :: rr,ur,pr,rl,ul,pl,cr,cl

    rr = WR(1)
    ur = WR(2)
    pr = WR(3)
    rl = WL(1)
    ul = WL(2)
    pl = WL(3)

    cr = evalc(WR, gam)
    cl = evalc(WL, gam)

    p0 = pl*rr*cr + pr*rl*cl + (ul-cr)*(rl*cl)*(rr*cr)
    p0 = p0 / (rl*cl + rr*cr)

    if(p0.le.0) p0 = 1e-3

    return
end function acoustic_approx
!***************************************************************    
function flfrdu(p, WL, WR, gam) result(f)
    implicit none
    real(p2), intent(in)  :: p, WL(3), WR(3), gam
    real(p2)              :: f
    real(p2)              :: fl, fr
    if(p > WL(3)) fl =       fshock(p, WL, gam)
    if(p < WL(3)) fl = frarefaction(p, WL, gam)
    if(p > WR(3)) fr =       fshock(p, WR, gam)
    if(p < WR(3)) fr = frarefaction(p, WR, gam)
    f = fl + fr + WR(2) - WL(2)
    return
end function flfrdu
!***************************************************************    
function dflfrdu(p, WL, WR, gam) result(df)
    implicit none
    real(p2), intent(in)  :: p, WL(3), WR(3), gam
    real(p2)              :: df
    real(p2)              :: dfl, dfr
    if(p > WL(3)) dfl =       dfshock(p, WL, gam)
    if(p < WL(3)) dfl = dfrarefaction(p, WL, gam)
    if(p > WR(3)) dfr =       dfshock(p, WR, gam)
    if(p < WR(3)) dfr = dfrarefaction(p, WR, gam)
    df = dfl + dfr
    return
end function dflfrdu
!***************************************************************    
function fshock(pin, WK, gam) result(fs)
    implicit none
    real(p2), intent(in)  :: pin, WK(3), gam
    real(p2)              :: fs

    fs = 2.0d0 / ((gam + 1.0d0) * WK(1)) &
         / (pin + (gam-1.0d0)/(gam+1.0d0) * WK(3))
    fs = sqrt(fs)
    fs = (pin - WK(3)) * fs
    return
end function fshock
!***************************************************************    
function frarefaction(pin, WK, gam) result(fr)
    implicit none
    real(p2), intent(in)  :: pin, WK(3), gam
    real(p2)              :: fr
    real(p2)              :: cknown

    cknown = evalc(WK, gam)
    fr = (2.0d0*cknown / (gam-1.0d0)) &
         * ((pin/WK(3))**((gam-1.0d0)/(2.0d0*gam)) - 1.0d0)

    return
end function frarefaction
!***************************************************************    
function dfshock(pin, WK, gam) result(dfs)
    implicit none
    real(p2), intent(in)  :: pin, WK(3), gam
    real(p2)              :: dfs

    dfs = 2.0d0 / ((gam + 1.0d0) * WK(1)) &
          / (pin + (gam-1.0d0)/(gam+1.0d0) * WK(3))
    if(dfs.ge.0) then
        dfs = sqrt(dfs)
    else
        write(*,*) 'Negative fshock sqrt', dfs
        stop
    end if
    dfs = (1.0d0 - (pin - WK(3)) &
          / (2.0d0*(pin + (gam-1.0d0)/(gam+1.0d0) * WK(3)))) &
          * dfs

    return
end function dfshock
!***************************************************************    
function dfrarefaction(pin, WK, gam) result(dfr)
    implicit none
    real(p2), intent(in)  :: pin, WK(3), gam
    real(p2)              :: dfr
    real(p2)              :: cknown

    cknown = evalc(WK, gam)

    dfr = (1.0d0 / (WK(1)*cknown)) &
          * (pin/WK(3))**(-(gam+1.0d0)/(2.0d0*gam))

    return
end function dfrarefaction
!***************************************************************    

end module riemann_solver
