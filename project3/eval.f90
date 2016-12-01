module eval
use prec
use glob
contains

function evalU(W, gam) result(U)
    ! Evaluate speed of sound
    implicit none
    real(p2) :: W(3), gam, U(3)
    U(1) = W(1)
    U(2) = W(1) * W(2)
    U(3) = W(3) / (gam - 1.0d0) + W(1) * W(2) * W(2) / 2.0d0
    return
end function evalU

function evalW(U, gam) result(W)
    ! Evaluate speed of sound
    implicit none
    real(p2) :: U(3), gam, W(3)
    W(1) = U(1)
    W(2) = U(2) / U(1)
    W(3) = evalp(U, gam)
    return
end function evalW

function evalF(U, gam) result(F)
    ! Evaluate convective flux vector
    implicit none
    real(p2) :: U(3), gam, F(3), p
    p = evalp(U, gam)
    F(1) = U(2)
    F(2) = p + U(2) * U(2) / U(1)
    F(3) = (U(3) + p) * U(2) / U(1)
    return
end function evalF

function evalFW(W, gam) result(F)
    ! Evaluate convective flux vector
    implicit none
    real(p2) :: W(3), gam, F(3)
    F(1) = W(1)*W(2)
    F(2) = F(1) * W(2) + W(3)
    F(3) = (W(3) / (gam - 1.0d0) + F(1) * W(2) * 0.5d0 + W(3)) * W(2)
    return
end function evalFW

function evalp(U, gam) result(p)
    ! Evaluate pressure
    implicit none
    real(p2) :: U(3), gam, p
    p = (U(3) - U(2) * U(2) / (2.0d0*U(1))) * (gam - 1.0d0)
    return
end function evalp

function evalc(W, gam) result(c)
    ! Evaluate speed of sound
    implicit none
    real(p2) :: W(3), gam, c
    c = gam * W(3) / W(1)
    if(c.ge.0.0d0) then
        c = sqrt(c)
    else
        write(*,*) 'Negative sqrt ms', c
        stop
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
        stop
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

subroutine update_wall(W)
    implicit none
    real(p2) :: W(3,ntot)
    integer  :: i
    do i = 1, nhalo
        W(:,i) = W(:,2*nhalo - i + 1) ! Copy density and pressure
        W(2,i) = -W(2,i) ! Reflect velocity

        W(:,ntot-i+1) = W(:,ntot - 2*nhalo + i) ! Copy density and pressure
        W(2,ntot-i+1) = -W(2,ntot-i+1) ! Reflect velocity
    end do
end subroutine

subroutine update_U_c(U, W, c, gam)
    implicit none
    real(p2) :: U(3,ntot), W(3,ntot), c(ntot), gam
    integer  :: i
    do i = 1, ntot
        U(:,i) = evalU(W(:,i), gam)
        c(i)   = evalc(W(:,i), gam)
    end do
end subroutine

subroutine update_W_c(U, W, c, gam)
    implicit none
    real(p2) :: U(3,ntot), W(3,ntot), c(ntot), gam
    integer  :: i
    do i = 1, ntot
        W(:,i) = evalW(U(:,i), gam)
        c(i)   = evalc(W(:,i), gam)
    end do
end subroutine

subroutine update_W(U, W, gam)
    implicit none
    real(p2) :: U(3,ntot), W(3,ntot), gam
    integer  :: i
    do i = 1, ntot
        W(:,i) = evalW(U(:,i), gam)
    end do
end subroutine

end module eval
