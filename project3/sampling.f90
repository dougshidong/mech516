module sample
use prec
use eval
contains
subroutine sampling(WL, WR, gam, S, WS, rhols, rhors, W)
    ! Given WL, WR, gam and S
    ! Return W at S
    ! Will first determine the region of the flow to properly assign W
    implicit none
    real(p2), intent(in)  :: WL(3), WR(3), WS(3), rhols, rhors
    real(p2), intent(in)  :: gam, S
    real(p2), intent(out) :: W(3)
    real(p2)              :: ck, cs, ms, SS, SH, ST, pratio
    integer               :: printt

    W = 0
    printt = 0
    if(S.lt.WS(2)) then ! Left from CS

        pratio = WS(3)/WL(3)
        ck = evalc(WL, gam)

        if(pratio.gt.1.0d0) then ! Left shock

            ms = evalms(pratio, gam) 
            SS = WL(2) - ck * ms
            if(S.lt.SS) then ! Left-side of left shock
                if(printt.eq.1) write(*,*) 'left of left shock'
                W = WL
            else ! Right-side of left shock
                if(printt.eq.1) write(*,*) 'right of left shock'
                W = WS
                W(1) = rhols
            end if

        else ! Left Rarefaction

            SH = WL(2) - ck
            if(S.lt.SH) then ! Left-side of left rarefaction
                if(printt.eq.1) write(*,*) 'left of left rare'
                W = WL
            else
                cs = evalc_isentropic(ck, pratio, gam)
                ST = WS(2) - cs
                if(S.gt.ST) then ! Right-side of left rarefaction
                    if(printt.eq.1) write(*,*) 'right of left rare'
                    W = WS
                    W(1) = rhols
                else ! Inside left rarefaction
                    if(printt.eq.1) write(*,*) 'inside of left rare'
                    call evalrarefaction(WL, S, gam, 0, W)
                end if
            end if

        end if

    else ! Right from CS
        pratio = WS(3)/WR(3)
        ck = evalc(WR, gam)

        if(pratio.gt.1.0d0) then ! Right shock

            ms = evalms(pratio, gam) 
            SS = WR(2) + ck * ms
            if(S.gt.SS) then ! Right-side of right shock
                W = WR
                if(printt.eq.1) write(*,*) 'right of right shock'
            else ! Left-side of Right shock
                W = WS
                W(1) = rhors
                if(printt.eq.1) write(*,*) 'left of right shock'
            end if

        else ! Right Rarefaction

            SH = WR(2) + ck
            if(S.gt.SH) then ! Right-side of right rarefaction
                W = WR
                if(printt.eq.1) write(*,*) 'right of right rare'
            else
                cs = evalc_isentropic(ck, pratio, gam)
                ST = WS(2) + cs
                if(S.lt.ST) then ! Left-side of right rarefaction
                    if(printt.eq.1) write(*,*) 'left of right rare'
                    W = WS
                    W(1) = rhors
                else ! Inside right rarefaction
                    if(printt.eq.1) write(*,*) 'inside of right rare'
                    call evalrarefaction(WR, S, gam, 1, W)
                end if
            end if

        end if

    end if
end subroutine sampling
end module sample
