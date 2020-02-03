module m_pnsubs
    implicit none
contains
    !*==PANDEF.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020






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

    subroutine pandef(ieldef)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ieldef
        !
        ! Local variables
        !
        real :: anpan, sbmag, sbmax, sbmin
        integer :: iel, iel1, iel2, irpn, npanmax, npanmin
        real, dimension(nex) :: sbtot
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Sets default paneling parameters
        !-----------------------------------------------
        !
        if (ieldef==0) then
            iel1 = 1
            iel2 = nbel
        else
            iel1 = ieldef
            iel2 = ieldef
        endif
        !
        sbmin = 0.0
        sbmax = 0.0
        do iel = iel1, iel2
            sbtot(iel) = sb(iblast(iel)) - sb(ibfrst(iel))
            if (netype(iel)==0) then
                sbmin = min(sbmin, sbtot(iel))
                sbmax = max(sbmax, sbtot(iel))
            endif
        enddo
        !---- use mean element arclength for scaling
        sbmag = 0.5 * (sbmin + sbmax)
        !
        anpan = float(npandef)
        !---- set min, max # panels per element
        npanmin = npandef / 2
        npanmax = 2 * npandef
        !
        do iel = iel1, iel2
            cvex(iel) = 1.0
            smof(iel) = 1.0
            fsle(iel) = 0.6
            fste(iel) = 0.6
            !
            if (sbtot(iel)==0.0) then
                npan(iel) = 1
            elseif (netype(iel)==0) then
                !-------- solid airfoil element
                npan(iel) = int(anpan * (sbtot(iel) / sbmag)**fpandef)
                npan(iel) = max(npanmin, npan(iel))
                npan(iel) = min(npanmax, npan(iel))
                cvex(iel) = 1.0
                fsle(iel) = 0.6
                fste(iel) = 0.6
                !
            elseif (netype(iel)==1 .or. netype(iel)==2 .or. netype(iel)   &
                    & ==5) then
                !-------- wake sheet, wake line or source line element
                npan(iel) = int(anpan / 4.0 * (sbtot(iel) / sbmag)**fpandef)
                !---- not less than 2 panels for these...
                npan(iel) = max(2, npan(iel))
                cvex(iel) = 1.0
                fsle(iel) = 0.6
                fste(iel) = 5.0
            else
                !-------- point or ring element (these won't really be used)
                npan(iel) = 1
                cvex(iel) = 1.0
                fsle(iel) = 1.0
                fste(iel) = 1.0
            endif
            !
            !------ refinement-station parameters
            do irpn = 1, nrpnx
                crrat(irpn, iel) = 1.0
                srpn1(irpn, iel) = 0.0
                srpn2(irpn, iel) = 0.0
            enddo
        enddo
        !
        !---- Set flag for spacing defined, unset user-defined flag
        lrspcdef = .true.
        lrspcusr = .false.
        !
    end subroutine pandef
    !*==PANCOP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PANDEF


    subroutine pancop
        use i_dfdc
        use m_geom, only : itpset, xypspl
        use m_geutil, only : minmax
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: dummy
        integer :: i, ib, iel, ip, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------
        !     Copies XB,YB arrays directly into panel arrays.
        !------------------------------------------------------------
        !
        if (nbtot>ipx) stop 'PANCOP: Array overflow on IPX.'
        !
        !---- total number of panel nodes
        nptot = nbtot
        nel = nbel
        !c      write(60,*) 'PANCOP '
        !c      write(60,*) 'NBEL,NBTOT ',NBEL,NBTOT
        !
        !---- panel node locations
        do ib = 1, nbtot
            ip = ib
            !c        write(60,99) ib,xb(ip),yb(ip),xp(ip),yp(ip)
            xp(ip) = xb(ib)
            yp(ip) = yb(ib)
        enddo
        99   format (i5, 6f12.6)
        !
        !---- reference point location for each element
        do iel = 1, nel
            xprefe(iel) = xbrefe(iel)
            yprefe(iel) = ybrefe(iel)
            !
            xpcent(iel) = xbcen2dt(iel)
            xpcent(iel) = xbcen2dt(iel)
        enddo
        !
        !---- set new stuff for each element
        do iel = 1, nel
            !c       write(60,*) IEL,IPFRST(IEL),IPLAST(IEL),IBFRST(IEL),IBLAST(IEL)
            ipfrst(iel) = ibfrst(iel)
            iplast(iel) = iblast(iel)
            if (iplast(iel)>ipx) stop 'PANCOP: Array overflow on IPX'
            !
            !------ set prescribed panel number to actual new panel number
            npan(iel) = iplast(iel) - ipfrst(iel) + 1
            !
            !------ spline and set other stuff for element IEL
            call xypspl(iel)
            !
            !------ set indices of TE panel, if any
            call itpset(iel)
            !
            i = ipfrst(iel)
            n = iplast(iel) - ipfrst(iel) + 1
            call minmax(n, xp(i), xpmine(iel), xpmaxe(iel))
            call minmax(n, yp(i), ypmine(iel), ypmaxe(iel))
        enddo
        !
        !---- set overall min,max over all elements
        call minmax(nel, xpmine, xpmin, dummy)
        call minmax(nel, ypmine, ypmin, dummy)
        call minmax(nel, xpmaxe, dummy, xpmax)
        call minmax(nel, ypmaxe, dummy, ypmax)
        !
        !---- we now have a new geometry... invalidate any existing solution
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
        lgsame = .true.
        !
    end subroutine pancop
    !*==PANGEN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PANCOP

    subroutine pangen(lquery)
        use i_dfdc
        use m_geom, only : itpset, xypspl
        use m_spline, only : seval
        use m_sgutil, only : sgcurv
        use m_geutil, only : minmax
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: lquery
        !
        ! Local variables
        !
        real :: dummy, sbi
        integer :: i, ib, ib1, ib2, iel, ig, ip, ip1, ip2, &
                & ipdel, kel, kp, kp1, kp2, kpd, n, nbpts, ng, &
                & nppts
        logical, save :: lcplot
        real, dimension(ipx) :: sg
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------------
        !     Generates panel nodes X,Y from buffer-airfoil XB,YB arrays.
        !     If LQUERY=T, then user's interaction is requested.
        !------------------------------------------------------------------
        !
        !
        data lcplot/.false./
        !
        if (nbel==0) return
        nel = nbel
        !
        !---- process each element
        do iel = 1, nbel
            !
            ib1 = ibfrst(iel)
            ib2 = iblast(iel)
            nbpts = ib2 - ib1 + 1
            ip1 = ipfrst(iel)
            ip2 = iplast(iel)
            nppts = ip2 - ip1 + 1
            if (ldbg) then
                write (*, *) 'PANGEN IEL ', iel
                write (*, *) 'IB1,IB2 ', ib1, ib2
                write (*, *) 'IP1,IP2 ', ip1, ip2
            endif
            !
            if (nppts<=1 .and. nbpts<=1) then
                !------- don't try to panel point or ring element
                ng = 1
                sg(1) = 0.
                !
            else
                !
                if (.not.(lquery)) then
                    !--------- no initial paneling for SGCURV
                    sg(1) = 0.
                    !-------- set current spacing so SGCURV has something to plot
                elseif (nppts>1) then
                    !--------- set SG from existing paneled geometry, if defined
                    ng = nppts
                    do ip = ip1, ip2
                        ig = ip - ip1 + 1
                        sg(ig) = (sp(ip) - sp(ip1)) / (sp(ip2) - sp(ip1))
                    enddo
                elseif (nbpts>1) then
                    !--------- or, set SG from buffer geometry
                    ng = nbpts
                    do ib = ib1, ib2
                        ig = ib - ib1 + 1
                        sg(ig) = (sb(ib) - sb(ib1)) / (sb(ib2) - sb(ib1))
                    enddo
                    if (npan(iel)<=0) npan(iel) = nbpts
                else
                    write (*, *)
                    write (*, *) '  No current paneling to modify'
                    write (*, *) '  Input case and/or Execute PANE command'
                    return
                    !
                endif
                !
                !         WRITE(*,1150) IEL
                ! 1150    FORMAT(
                !     & /' ============================================================='
                !     &//'  Element', I4)
                !
                !------- surface element... generate spacing array SG
                !
                call sgcurv(lquery, lcplot, nbpts, xb(ib1), xbs(ib1), yb(ib1), &
                        & ybs(ib1), sb(ib1), ipx, ng, sg, npan(iel), cvex(iel), &
                        & smof(iel), fsle(iel), fste(iel), nrpnx, srpn1(1, iel)&
                        &, srpn2(1, iel), crrat(1, iel))
            endif
            !
            if (lquery) then
                !---- Set flag for spacing defined, unset user-defined flag
                lrspcdef = .true.
                lrspcusr = .true.
            endif
            !
            !------ starting and ending indices for this element IEL
            if (iel==1) then
                ip1 = 1
                ip2 = ng
            else
                ip1 = iplast(iel - 1) + 1
                ip2 = iplast(iel - 1) + ng
            endif
            if (ip2>ipx) stop 'PANGEN: Array overflow on IPX'
            !
            if (lquery) then
                ipdel = ip2 - iplast(iel)
                if (nptot + ipdel>ipx) stop 'PANGEN: Array overflow on IPX'
                !
                !-------- move current points in subsequent elements to make/take-up room
                if (ipdel>0) then
                    kp1 = nptot
                    kp2 = iplast(iel) + 1
                    kpd = -1
                else
                    kp1 = iplast(iel) + 1
                    kp2 = nptot
                    kpd = +1
                endif
                do kp = kp1, kp2, kpd
                    sp(kp + ipdel) = sp(kp)
                    xp(kp + ipdel) = xp(kp)
                    yp(kp + ipdel) = yp(kp)
                    xps(kp + ipdel) = xps(kp)
                    yps(kp + ipdel) = yps(kp)
                enddo
                do kel = iel + 1, nel
                    ipfrst(kel) = ipfrst(kel) + ipdel
                    iplast(kel) = iplast(kel) + ipdel
                enddo
            endif
            !
            !------ set panel nodes at fractional arc length locations SG
            do ip = ip1, ip2
                ig = ip - ip1 + 1
                sbi = sb(ib1) + (sb(ib2) - sb(ib1)) * sg(ig)
                !
                xp(ip) = seval(sbi, xb(ib1), xbs(ib1), sb(ib1), nbpts)
                yp(ip) = seval(sbi, yb(ib1), ybs(ib1), sb(ib1), nbpts)
            enddo
            ipfrst(iel) = ip1
            iplast(iel) = ip2
            !
            !------ centroid location
            xpcent(iel) = xbcen2dt(iel)
            ypcent(iel) = ybcen2dt(iel)

            !------ spline and set other stuff for element IEL
            call xypspl(iel)
            !
            !------ set indices of TE panel, if any
            call itpset(iel)
            !
            i = ipfrst(iel)
            n = iplast(iel) - ipfrst(iel) + 1
            call minmax(n, xp(i), xpmine(iel), xpmaxe(iel))
            call minmax(n, yp(i), ypmine(iel), ypmaxe(iel))
        enddo
        !
        !---- set total number of panel vertices
        nptot = iplast(nel)
        !
        !---- set overall min,max over all elements
        call minmax(nel, xpmine, xpmin, dummy)
        call minmax(nel, ypmine, ypmin, dummy)
        call minmax(nel, xpmaxe, dummy, xpmax)
        call minmax(nel, ypmaxe, dummy, ypmax)
        !
        !---- default reference point location for each element
        do iel = 1, nel
            xprefe(iel) = 0.0
            yprefe(iel) = 0.0
        enddo
        !
        !---- we now have a new geometry... invalidate any existing solution
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
        lqspec = .false.
        lgsame = .false.
        !
    end subroutine pangen
    !*==PANWRT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PANGEN



    subroutine panwrt
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: iel, irpn, lu
        character(80) :: pfndef
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------
        !     Writes paneling-parameter file
        !-----------------------------------------
        !
        1000 format (a)
        !
        lu = 8
        20   do
            !
            !---- use existing paneling-parameter filename (if any) as the default
            write (*, 1200) pfile(1:60)
            1200    format (/' Enter output filename:  ', a)
            !
            read (*, 1000) pfndef
            if (pfndef(1:1)/=' ') then
                pfile = pfndef
            elseif (pfile(1:1)==' ') then
                cycle
            endif
            exit
        enddo
        !
        open (lu, file = pfile, status = 'UNKNOWN', err = 20)
        rewind (lu)
        !
        write (lu, *) nel, nrpnx
        do iel = 1, nel
            write (lu, *) npan(iel)
            write (lu, *) cvex(iel), smof(iel), fsle(iel), fste(iel)
            do irpn = 1, nrpnx
                write (lu, *) srpn1(irpn, iel), srpn2(irpn, iel), &
                        & crrat(irpn, iel)
            enddo
        enddo
        close (lu)
        !
    end subroutine panwrt
    !*==PANGET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PANWRT



    subroutine panget(fname, error)
        use i_dfdc
        use m_userio, only : asks, strip
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: error
        character(*) :: fname
        !
        ! Local variables
        !
        integer :: iel, irpn, lu, nel0, nfname, nrpn0
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------
        !     Reads paneling-parameter file FNAME
        !-----------------------------------------
        !
        lu = 8
        !
        if (fname(1:1)==' ') call asks(&
                &'Enter paneling-parameter filename^'&
                &, fname)
        !
        !---- strip off leading blanks and get number of characters
        call strip(fname, nfname)
        !
        write (*, *)
        !
        open (lu, file = fname, status = 'OLD', err = 90)
        rewind (lu)
        !
        read (lu, *, err = 80) nel0, nrpn0
        nel0 = min(nel0, nex)
        nrpn0 = min(nrpn0, nrpnx)
        !
        do iel = 1, nel0
            read (lu, *, err = 80) npan(iel)
            read (lu, *, err = 80) cvex(iel), smof(iel), fsle(iel), &
                    & fste(iel)
            do irpn = 1, nrpn0
                read (lu, *, err = 80) srpn1(irpn, iel), srpn2(irpn, iel), &
                        & crrat(irpn, iel)
            enddo
        enddo
        close (lu)
        write (*, *) 'Paneling parameters read from file ', &
                & fname(1:nfname)
        lrspcdef = .true.
        lrspcusr = .true.
        error = .false.
        return
        !
        80   close (lu)
        write (*, *) 'READ error on panel-parameter file ', &
                & fname(1:nfname)
        error = .true.
        return
        !
        !      WRITE(*,*) 'OPEN error on panel-parameter file ', FNAME(1:NFNAME)
        90   error = .true.
    end subroutine panget
end module m_pnsubs
