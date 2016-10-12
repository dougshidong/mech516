module eval
use prec
contains
function evalc(W, gam) result(c)
    ! Evaluate speed of sound
    implicit none
    real(p2) :: W(3), gam, c
    c = gam * W(3) / W(1)
    if(c.ge.0.0d0) then
        c = sqrt(c)
    else
        write(*,*) 'Negative sqrt c', c
    end if
    return
end function evalc

function evalms(pratio, gam) result(ms)
    ! Evaluate Mach Number
    implicit none
    real(p2) :: pratio, gam, ms
    ms = ( (gam+1.0d0) * pratio + (gam-1.0d0) ) / (2.0d0*gam)
    if(ms.ge.0.0d0) then
        ms = sqrt(ms)
    else
        write(*,*) 'Negative sqrt ms', ms
    end if
    return
end function evalms

function evalrho_hugoniot(rhok, pratio, gam) result(rho)
    ! Evaluate density across a shock based on Hugoniot relation
    implicit none
    real(p2) :: rhok, pratio, gam, rho
    real(p2) :: gamr
    gamr = (gam-1.0d0)/(gam+1.0d0)
    rho = rhok * (pratio + gamr) / (gamr*pratio + 1.0d0)
    return
end function evalrho_hugoniot

function evalrho_isentropic(rhok, pratio, gam) result(rho)
    ! Evaluate density across a rarefaction based on isentopic relation
    implicit none
    real(p2) :: rhok, pratio, gam, rho
    rho = pratio**(1.0d0/gam) * rhok
    return
end function evalrho_isentropic

function evalc_isentropic(ck, pratio, gam) result(c)
    ! Evaluate speed of sound across rarefaction based on isentropic relation
    implicit none
    real(p2) :: ck, pratio, gam, c
    c = pratio**((gam-1.0d0)/(2.0d0*gam)) * ck
    return
end function evalc_isentropic

subroutine evalrarefaction(WK, S, gam, direction, W)
    ! Evaluate primitive variables inside expansion fan
    implicit none
    real(p2) :: WK(3), gam, W(3), ctemp, S
    integer  :: direction
    real(p2) :: twogp1, twogm1, gm1gp1, big

    ctemp = evalc(WK, gam)
    if(direction.eq.1) ctemp = -ctemp
    
    twogp1 = 2.0 / (gam+1.0)
    twogm1 = 2.0 / (gam-1.0)
    gm1gp1 = (gam-1.0) / (gam+1.0)
    big = (twogp1 + gm1gp1 / ctemp * (WK(2) - S))

    W(1) = WK(1) * big**twogm1
    W(2) = twogp1 * (ctemp + WK(2) / twogm1 + S)
    W(3) = WK(3) * big**(twogm1*gam)
    return
end subroutine evalrarefaction

end module eval
