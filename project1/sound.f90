module eval
contains
function evalc(w, gam) result(c)
    use prec
    implicit none
    real(p2) :: w(3), c
    c = gam * w(3) / w(1)
    if(c.ge.0.0d0) then
        c = sqrt(c)
    else
        write(*,*) 'Negative sqrt c', c
    end if
    return
end function evalc

function evalms(pratio, gam) result(ms)
    use prec
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

function evalrhohugoniot(rhok, pratio, gam) result(rho)
    use prec
    implicit none
    real(p2) :: rhok, pratio, gam, rho
    real(p2) :: gamr
    gamr = (gam-1.0d0)/(gam+1.0d0)
    rho = rhok (pratio + gamr) / (gamr*pratio + 1.0d0)
    return
end function evalrhohugoniot

function evalrhoisentropic(rhok, pratio, gam) result(rho)
    use prec
    implicit none
    real(p2) :: rhok, pratio, gam, rho
    real(p2) :: gamr
    gamr = (gam-1.0d0)/(gam+1.0d0)
    rho = rhok (pratio + gamr) / (gamr*pratio + 1.0d0)
    return
end function evalrhohugoniot

nd module eval
