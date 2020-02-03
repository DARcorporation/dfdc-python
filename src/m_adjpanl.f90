module m_adjpanl
    implicit none
contains
    !*==ADJPANL.f90  processed by SPAG 7.25DB at 08:51 on  3 Feb 2020
    !=========================================================================
    ! DFDC (Ducted Fan Design Code) is an aerodynamic and aeroacoustic design
    ! and analysis tool for aircraft with propulsors in ducted fan
    ! configurations.
    !
    ! This software was developed under the auspices and sponsorship of the
    ! Tactical Technology Office (TTO) of the Defense Advanced Research
    ! Projects Agency (DARPA).
    !
    ! Copyright (c) 2004, 2005, Booz Allen Hamilton Inc., All Rights Reserved
    !
    ! This program is free software; you can redistribute it and/or modify it
    ! under the terms of the GNU General Public License as published by the
    ! Free Software Foundation; either version 2 of the License, or (at your
    ! option) any later version.
    !
    ! This program is distributed in the hope that it will be useful, but
    ! WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    ! General Public License for more details.
    !
    ! You should have received a copy of the GNU General Public License along
    ! with this program; if not, write to the Free Software Foundation, Inc.,
    ! 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
    !
    ! Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (drela@mit.edu)
    ! Program Management: Brad Tousley, Paul Eremenko (eremenko@alum.mit.edu)
    !
    !=========================================================================
    !

    subroutine adjpanl
        use i_dfdc
        use m_sgutil2, only : sgshft, sgrenum, isgfind, sgcopy
        use m_geom, only : itpset, xypspl
        use m_spline, only : seval, sinvrt
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: itx = 600
        !
        ! Local variables
        !
        integer :: i, iel, ip, ip1, ip2, ipmov, iprcb, iprcb2, &
                & iprdw, iprdw2, ite, n, ncb, nduct, nmov, &
                & nnew, nold, npnew, nr, nt1, nt1new, nt2, &
                & nt2new
        logical :: lovrmax
        real :: sbeg, send, sle1, sle2, sold, stcb, stdw, ste, &
                & xle1, xle2, xlemax, xte1, xte2, xtemin, ycb, &
                & ydw, yle1, yle2, yte1, yte2
        real, dimension(itx) :: ss, st1, st2, xt1, xt2, xts1, &
                & xts2, yt1, yt2, yts1, yts2
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------
        !     Readjusts panel nodes for CB and duct to give
        !     spacing compatible with a vortex wake grid between
        !     the elements.
        !-----------------------------------------------------------------
        !
        !
        !---- Flag for spacing in overlap area to TE (T uses max # of points from CB
        !     or duct for overlap area, F simply copies spacing from opposite body)
        lovrmax = .true.
        !
        !---- Copy CB panel geometry to temporary array
        iel = 1
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        n = ip2 - ip1 + 1
        if (n>itx) then
            write (*, *) 'ADJPANL: ITX too small for N ', itx, n
            stop
        endif
        do i = 1, n
            ip = ip1 + i - 1
            xt1(i) = xp(ip)
            yt1(i) = yp(ip)
            st1(i) = sp(ip)
            xts1(i) = xps(ip)
            yts1(i) = yps(ip)
        enddo
        nt1 = n
        !
        xle1 = xple(iel)
        yle1 = yple(iel)
        sle1 = sple(iel)
        !---- Use TE point for CB foil
        xte1 = xp(ip1)
        yte1 = yp(ip1)
        !
        !---- Copy duct panel geometry to temporary array
        iel = 2
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        n = ip2 - ip1 + 1
        if (n>itx) then
            write (*, *) 'ADJPANL: ITX too small for N ', itx, n
            stop
        endif
        do i = 1, n
            ip = ip1 + i - 1
            xt2(i) = xp(ip)
            yt2(i) = yp(ip)
            st2(i) = sp(ip)
            xts2(i) = xps(ip)
            yts2(i) = yps(ip)
        enddo
        nt2 = n
        !
        xle2 = xple(iel)
        yle2 = yple(iel)
        sle2 = sple(iel)
        !---- Use inner TE point for duct foil
        xte2 = xp(ip2)
        yte2 = yp(ip2)
        !
        !---- Axial region for rotor(s) defined by overlap of LE/TE's
        xlemax = max(xle1, xle2)
        xtemin = min(xte1, xte2)
        !
        !cc      GO TO 100
        !
        !
        !
        !
        !-----------------------------------------------------------------
        !-----------------------------------------------------------------
        !---- start process with upstream rotor
        !
        nr = 1
        !---- Check rotor station for sanity
        if (xdisk(nr)<xlemax .or. xdisk(nr)>xtemin) then
            write (*, *) '*Rotor station out of bounds ', xdisk(nr)
            write (*, *) '  XLE1 =', xle1, '  XTE1 =', xte1
            write (*, *) '  XLE2 =', xle2, '  XTE2 =', xte2
            !c        CALL ASKR('Enter rotor disk X location^',XDISK(NR))
            !---- default rotor is halfway along passage
            xdisk(nr) = 0.5 * (xlemax + xtemin)
            write (*, *) '*Rotor station set to ', xdisk(nr)
        endif
        if (ldbg) write (*, *) ' '
        !
        !-----------------------------------------------------------------
        !---- find location of rotor on CB wall
        stcb = st1(1)
        call sinvrt(stcb, xdisk(nr), xt1, xts1, st1, nt1)
        ycb = seval(stcb, yt1, yts1, st1, nt1)
        rhub(nr) = ycb
        if (ldbg) write (*, *) 'Rhub,S on CB @ ', rhub(nr), stcb
        !---- Find nearest point on CB corresponding to rotor
        iprcb = isgfind(stcb, st1, nt1)
        iprotcb(nr) = iprcb
        if (ldbg) write (*, *) 'Nearest CB point to rotor @ ', iprcb, &
                & st1(iprcb)
        !---- Readjust point spacing to put grid point at rotor intersection
        sbeg = st1(1)
        send = st1(nt1)
        sold = st1(iprcb)
        call sgshft(sbeg, send, sold, stcb, st1, ss, nt1)
        !---- Generate new X,Y coordinates
        call sgsplupd(xt1, xts1, yt1, yts1, st1, nt1, xt1, xts1, yt1, yts1, ss, nt1)
        call sgcopy(ss, st1, nt1)
        !
        !
        !-----------------------------------------------------------------
        !---- find location of rotor on duct wall
        stdw = st2(nt2)
        call sinvrt(stdw, xdisk(nr), xt2, xts2, st2, nt2)
        ydw = seval(stdw, yt2, yts2, st2, nt2)
        rtip(nr) = ydw
        if (ldbg) write (*, *) 'Rtip,S on duct @ ', rtip(nr), stdw
        !---- Find nearest point on duct corresponding to rotor
        iprdw = isgfind(stdw, st2, nt2)
        iprotdw(nr) = iprdw
        if (ldbg) write (*, *) 'Nearest Duct point to rotor @ ', iprdw, &
                & st2(iprdw)
        !---- Readjust point spacing to put grid point at rotor intersection
        sbeg = sle2
        send = st2(nt2)
        sold = st2(iprdw)
        call sgshft(sbeg, send, sold, stdw, st2, ss, nt2)
        !---- Generate new X,Y coordinates
        call sgsplupd(xt2, xts2, yt2, yts2, st2, nt2, xt2, xts2, yt2, yts2, ss, nt2)
        call sgcopy(ss, st2, nt2)
        !
        !c      GO TO 100
        !
        !-----------------------------------------------------------------
        !---- Match spacing in overlap area downstream of rotor (to TE's)
        !     to define I direction slipstream grid
        !
        !---- Duct TE upstream of CB TE
        if (xte2<=xte1) then
            !
            if (ldbg) write (*, *) 'XTEduct <= XTECB ', xte2, xte1
            !---- Find point on CB corresponding to duct TE
            ste = st1(1)
            call sinvrt(ste, xte2, xt1, xts1, st1, nt1)
            if (ldbg) write (*, *) 'CB point @ Xduct TE @ ', ste
            !---- Find nearest point on CB corresponding to duct TE
            ite = 0
            ite = isgfind(ste, st1, nt1)
            if (ldbg) write (*, *) 'Nearest CB point to duct TE @ ', &
                    & ite, ste, st1(ite)
            !---- If this point can be moved (not TE point!!), align to duct TE
            if (ite>1) then
                !---- Readjust CB point spacing to put grid point at duct TE X location
                sbeg = st1(1)
                send = st1(iprcb)
                sold = st1(ite)
                call sgshft(sbeg, send, sold, ste, st1, ss, nt1)
                !---- Update CB spline definition
                call sgsplupd(xt1, xts1, yt1, yts1, st1, nt1, xt1, xts1, yt1, yts1, &
                        & ss, nt1)
                call sgcopy(ss, st1, nt1)
            endif
            !
            !---- Respace overlap zone from rotor to CB TE point
            nduct = nt2 - iprdw + 1
            ncb = iprcb - ite + 1
            !---- Check # points in overlap zone, use densest paneling from CB or duct
            if (ncb<nduct .or. (.not.lovrmax)) then
                !---- Respace CB points from rotor to duct TE point with # from duct
                nnew = nduct
                nold = ncb
                call sgcopy(st1, ss, nt1)
                call sgrenum(ss(ite), nold, nnew)
                call sgcopy(st1(iprcb), ss(ite + nnew - 1), nt1 - iprcb + 1)
                nt1new = nt1 + nnew - nold
                !---- Generate new X,Y coordinates
                call sgsplupd(xt1, xts1, yt1, yts1, st1, nt1, xt1, xts1, yt1, yts1, &
                        & ss, nt1new)
                call sgcopy(ss, st1, nt1new)
                nt1 = nt1new
                iprcb = iprotcb(nr) + nnew - nold
                iprotcb(nr) = iprcb
                !
            elseif (ncb>nduct) then
                !---- Respace duct points from rotor to duct TE point with # from CB
                nnew = ncb
                nold = nduct
                call sgcopy(st2, ss, nt2)
                call sgrenum(ss(iprdw), nold, nnew)
                nt2new = nt2 + nnew - nold
                !---- Generate new X,Y coordinates
                call sgsplupd(xt2, xts2, yt2, yts2, st2, nt2, xt2, xts2, yt2, yts2, &
                        & ss, nt2new)
                call sgcopy(ss, st2, nt2new)
                nt2 = nt2new
            endif
            !
            !
            !---- CB TE upstream of duct TE
        elseif (xte2>xte1) then
            !
            if (ldbg) write (*, *) 'XTEduct > XTECB ', xte2, xte1
            !---- Find point on duct corresponding to CB TE
            ste = st2(nt2)
            call sinvrt(ste, xte1, xt2, xts2, st2, nt2)
            !---- Find nearest point on duct corresponding to CB TE
            ite = isgfind(ste, st2, nt2)
            if (ldbg) write (*, *) 'Nearest Duct point to CB TE @ ', &
                    & ite, ste, st2(ite)
            !---- If this point can be moved (not TE point!!), align to CB TE
            if (ite<nt2) then
                !---- Readjust duct point spacing to put grid point at CB TE X location
                sbeg = st2(iprdw)
                send = st2(nt2)
                sold = st2(ite)
                call sgshft(sbeg, send, sold, ste, st2, ss, nt2)
                call sgsplupd(xt2, xts2, yt2, yts2, st2, nt2, xt2, xts2, yt2, yts2, &
                        & ss, nt2)
                call sgcopy(ss, st2, nt2)
            endif
            !
            !---- Respace overlap zone from rotor to CB TE point
            nduct = ite - iprdw + 1
            ncb = iprcb
            !---- Check # points in overlap zone, use densest paneling from CB or duct
            if (ncb>nduct .or. (.not.lovrmax)) then
                !---- Respace duct points from rotor to CB TE point with # from CB
                nnew = ncb
                nold = nduct
                call sgcopy(st2, ss, nt2)
                call sgrenum(ss(iprdw), nold, nnew)
                call sgcopy(st2(ite), ss(iprdw + nnew - 1), nt2 - ite + 1)
                nt2new = nt2 + nnew - nold
                !---- Generate new X,Y coordinates
                call sgsplupd(xt2, xts2, yt2, yts2, st2, nt2, xt2, xts2, yt2, yts2, &
                        & ss, nt2new)
                call sgcopy(ss, st2, nt2new)
                nt2 = nt2new
            elseif (ncb<nduct) then
                !---- Respace CB points from rotor to CB TE point with # from duct
                nnew = nduct
                nold = ncb
                call sgcopy(st1, ss, nt1)
                call sgrenum(ss(1), nold, nnew)
                call sgcopy(st1(nold + 1), ss(nnew + 1), nt1 - nold)
                nt1new = nt1 + nnew - nold
                !---- Generate new X,Y coordinates
                call sgsplupd(xt1, xts1, yt1, yts1, st1, nt1, xt1, xts1, yt1, yts1, &
                        & ss, nt1new)
                call sgcopy(ss, st1, nt1new)
                nt1 = nt1new
                iprcb = iprotcb(nr) + nnew - nold
                iprotcb(nr) = iprcb
            endif
            !
        endif
        !
        !      GO TO 100
        !
        !-----------------------------------------------------------------
        !-----------------------------------------------------------------
        !---- Process downstream rotor(s), fiddling grid to put move points
        !     to match downstream rotor line
        do nr = 2, nrotor
            !
            !---- Check rotor station for sanity
            if (ldbg) write (*, *) ' '
            if (xdisk(nr)<xlemax .or. xdisk(nr)>xtemin) then
                write (*, *) '*Rotor#N station out of bounds ', xdisk(nr)
                write (*, *) '  XLE1 =', xle1, '  XTE1 =', xte1
                write (*, *) '  XLE2 =', xle2, '  XTE2 =', xte2
                !c        CALL ASKR('Enter rotor disk X location^',XDISK)
                !---- default rotor #2 is halfway along passage
                xdisk(nr) = 0.5 * (xdisk(1) + xtemin)
                write (*, *) '*Rotor#N station set to ', xdisk(nr)
            endif
            !
            !-----------------------------------------------------------------
            !---- find location of rotor on CB wall
            stcb = st1(1)
            call sinvrt(stcb, xdisk(nr), xt1, xts1, st1, nt1)
            ycb = seval(stcb, yt1, yts1, st1, nt1)
            rhub(nr) = ycb
            if (ldbg) write (*, *) 'Rhub,S #N on CB @ ', rhub(nr), stcb
            !---- Find nearest point on CB corresponding to rotor
            iprcb2 = isgfind(stcb, st1, nt1)
            iprotcb(nr) = iprcb2
            if (ldbg) write (*, *) 'Nearest CB point to rotor @ ', &
                    & iprcb2, st1(iprcb2)
            !---- Readjust point spacing to put grid point at rotor intersection
            sbeg = st1(1)
            send = st1(iprcb)
            sold = st1(iprcb2)
            call sgshft(sbeg, send, sold, stcb, st1, ss, nt1)
            !---- Generate new X,Y coordinates
            call sgsplupd(xt1, xts1, yt1, yts1, st1, nt1, xt1, xts1, yt1, yts1, ss, &
                    & nt1)
            call sgcopy(ss, st1, nt1)
            !
            !-----------------------------------------------------------------
            !---- find location of rotor on duct wall
            stdw = st2(iprdw)
            call sinvrt(stdw, xdisk(nr), xt2, xts2, st2, nt2)
            ydw = seval(stdw, yt2, yts2, st2, nt2)
            rtip(nr) = ydw
            if (ldbg) write (*, *) 'Rtip,S #N on duct @ ', rtip(nr), &
                    & stdw
            !---- index of rotor#2 from rotor#1 must match offset from rotor#1 on CB
            iprdw2 = iprdw + (iprcb - iprcb2)
            iprotdw(nr) = iprdw2
            if (ldbg) write (*, *) 'Nearest Duct point to rotor @ ', &
                    & iprdw2, st2(iprdw2)
            !---- Readjust point spacing to put grid point at rotor intersection
            sbeg = st2(iprdw)
            send = st2(nt2)
            sold = st2(iprdw2)
            call sgshft(sbeg, send, sold, stdw, st2, ss, nt2)
            !---- Generate new X,Y coordinates
            call sgsplupd(xt2, xts2, yt2, yts2, st2, nt2, xt2, xts2, yt2, yts2, ss, &
                    & nt2)
            call sgcopy(ss, st2, nt2)
            !
        enddo ! loop over NROTOR=>2
        !
        !---- Paneling (points) modified to match on CB and duct and
        !     to match rotor lines
        !
        !-----------------------------------------------------------------
        !---- Before installing new panel geometry check IEL>2 to save old
        !     geometry
        if (nel>=3) then
            npnew = nt1 + nt2
            nmov = npnew - ipfrst(3) + 1
            if (nmov>0) then
                !---- Move existing elements up to give room for modified elements 1 and 2
                do iel = nel, 3, -1
                    ip1 = ipfrst(iel)
                    ip2 = iplast(iel)
                    do ip = ip2, ip1, -1
                        ipmov = ip + nmov
                        xp(ipmov) = xp(ip)
                        yp(ipmov) = yp(ip)
                    enddo
                    ipfrst(iel) = ip1 + nmov
                    iplast(iel) = ip2 + nmov
                enddo
            endif
        endif
        !
        !-----------------------------------------------------------------
        !---- Put adjusted body wall points back into paneled point arrays
        ip = 0
        !---- CB element
        iel = 1
        ip1 = ip + 1
        do i = 1, nt1
            ip = ip + 1
            xp(ip) = xt1(i)
            yp(ip) = yt1(i)
        enddo
        ip2 = ip
        ipfrst(iel) = ip1
        iplast(iel) = ip2
        i = ip1
        n = ip2 - ip1 + 1
        !---- spline and set other stuff for element IEL
        call xypspl(iel)
        !---- set indices of TE panel, if any
        call itpset(iel)
        !---- update pointers to CB wall point at rotor line(s)
        do nr = 1, nrotor
            iprotcb(nr) = ip1 - 1 + iprotcb(nr)
        enddo
        !
        !---- Put new duct wall points back into paneled point arrays
        iel = 2
        ip1 = ip + 1
        do i = 1, nt2
            ip = ip + 1
            xp(ip) = xt2(i)
            yp(ip) = yt2(i)
        enddo
        ip2 = ip
        ipfrst(iel) = ip1
        iplast(iel) = ip2
        i = ip1
        n = ip2 - ip1 + 1
        !---- spline and set other stuff for element IEL
        call xypspl(iel)
        !---- set indices of TE panel, if any
        call itpset(iel)
        !---- update pointer to duct wall point at rotor line
        do nr = 1, nrotor
            iprotdw(nr) = ip1 - 1 + iprotdw(nr)
            if (ldbg) write (*, *) 'NR IPROTCB IPROTDW  ', nr, &
                    & iprotcb(nr), iprotdw(nr)
        enddo
        !
        !---- Shift old points in any remaining defined elements to match new ordering
        if (nel>=3) then
            npnew = nt1 + nt2
            nmov = npnew + 1 - ipfrst(3)
            if (nmov<0) then
                !---- Move existing elements down to give room for modified elements 1 and 2
                do iel = 3, nel
                    ip1 = ipfrst(iel)
                    ip2 = iplast(iel)
                    do ip = ip1, ip2
                        ipmov = ip + nmov
                        xp(ipmov) = xp(ip)
                        yp(ipmov) = yp(ip)
                    enddo
                    ipfrst(iel) = ip1 + nmov
                    iplast(iel) = ip2 + nmov
                enddo
            endif
        endif
        !---- Reset total # of paneled points
        nptot = iplast(nel)
        !
        !cc        write(14,99) i,xt1(i),xx(i),xdisk(nr)
        99   format (i5, 5(1x, f12.6))
        !
        !---- Set flag for rotor defined (at least geometrically)
        !c      LROTOR = .TRUE.
        !
        !---- invalidate any existing solution
        lncvp = .false.
        lqaic = .false.
        lqgic = .false.
        lqcnt = .false.
        lsysp = .false.
        lgsys = .false.
        lgamu = .false.
        lgama = .false.
        lsigp = .false.
        lsigm = .false.
        !
    end subroutine adjpanl
    !*==SGSPLUPD.f90  processed by SPAG 7.25DB at 08:51 on  3 Feb 2020
    ! ADJPANL


    subroutine sgsplupd(x1, xs1, y1, ys1, s1, n1, x2, xs2, y2, ys2, s2, n2)
        use m_spline, only : segspl, seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: ix = 600
        !
        ! Dummy arguments
        !
        integer :: n1, n2
        real, dimension(*) :: s1, s2, x1, x2, xs1, xs2, y1, y2, &
                & ys1, ys2
        !
        ! Local variables
        !
        integer :: i
        real, dimension(ix) :: xx, yy
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------
        !     Generates points and spline definitions for spline curve
        !     X1,XS1,Y1,YS1,S1,N1 (with N1 points at S1) to new spline
        !     with spacing S2 with N2 points.
        !     New curve with spline and points returned in X2,XS2,Y2,YS2,S2,N2
        !-----------------------------------------------------------------
        !
        !
        !
        if (n1>ix .or. n2>ix) then
            write (*, *) 'SGSPLUPD number of input points exceeds IX ', &
                    & n1, n2
            stop
        endif
        !
        !---- Generate new X,Y coordinates at S2 locations
        do i = 1, n2
            xx(i) = seval(s2(i), x1, xs1, s1, n1)
            yy(i) = seval(s2(i), y1, ys1, s1, n1)
        enddo
        !
        !---- Replace points in output arrays for redefined curve
        do i = 1, n2
            x2(i) = xx(i)
            y2(i) = yy(i)
        enddo
        !---- Respline new points
        call segspl(x2, xs2, s2, n2)
        call segspl(y2, ys2, s2, n2)
        !
    end subroutine sgsplupd
end module m_adjpanl
