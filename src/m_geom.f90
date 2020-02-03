module m_geom
    implicit none
contains
    !*==XYPSPL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! BAKSUB
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

    subroutine xypspl(iel)
        use i_dfdc
        use m_spline, only : scalc, seval, segspl
        use m_geutil, only : lefind
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iel
        !
        ! Local variables
        !
        integer :: i, ip1, ip2, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------
        !     Splines panel-node coordinates
        !     and sets other element-related stuff.
        !---------------------------------------------
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        !
        i = ip1
        n = ip2 - ip1 + 1
        if (n>1) then
            call scalc(xp(i), yp(i), sp(i), n)
            call segspl(xp(i), xps(i), sp(i), n)
            call segspl(yp(i), yps(i), sp(i), n)
        endif
        !
        if (lxbod(iel)) then
            !----- axisymmetric body on axis
            if (xp(ip1)<xp(ip2)) then
                sple(iel) = sp(ip1)
                xple(iel) = xp(ip1)
                yple(iel) = yp(ip1)
                xpte(iel) = xp(ip2)
                ypte(iel) = yp(ip2)
            else
                sple(iel) = sp(ip2)
                xple(iel) = xp(ip2)
                yple(iel) = yp(ip2)
                xpte(iel) = xp(ip1)
                ypte(iel) = yp(ip1)
            endif
            !
        elseif (lbody(iel) .and. (lv1zr(iel) .or. lv2zr(iel))) then
            !----- axisymmetric body not closed on axis
            if (xp(ip1)<xp(ip2)) then
                sple(iel) = sp(ip1)
                xple(iel) = xp(ip1)
                yple(iel) = yp(ip1)
                xpte(iel) = xp(ip2)
                ypte(iel) = yp(ip2)
            else
                sple(iel) = sp(ip2)
                xple(iel) = xp(ip2)
                yple(iel) = yp(ip2)
                xpte(iel) = xp(ip1)
                ypte(iel) = yp(ip1)
            endif
            !
        elseif (lbody(iel) .and. .not.lxbod(iel)) then
            !----- body off axis
            call lefind(sple(iel), xp(i), xps(i), yp(i), yps(i), sp(i), n)
            xple(iel) = seval(sple(iel), xp(i), xps(i), sp(i), n)
            yple(iel) = seval(sple(iel), yp(i), yps(i), sp(i), n)
            xpte(iel) = 0.5 * (xp(ip1) + xp(ip2))
            ypte(iel) = 0.5 * (yp(ip1) + yp(ip2))
            !
        else
            !----- zero-thickness body
            sple(iel) = sp(ip1)
            xple(iel) = xp(ip1)
            yple(iel) = yp(ip1)
            xpte(iel) = xp(ip2)
            ypte(iel) = yp(ip2)
        endif
        !
        !---- set location for plotted element index
        if (netype(iel)==0 .or. netype(iel)==1 .or. netype(iel)==2) then
            !----- surface or axis line... set plot location at skin centroid
            xelnum(iel) = xpcent(iel)
            yelnum(iel) = ypcent(iel)
        else
            !----- point,ring,source line or vortex wake element...
            !      just set plot location at the point itself
            xelnum(iel) = xp(ip1)
            yelnum(iel) = yp(ip1)
        endif
        !
    end subroutine xypspl
    !*==CVPGEN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! XYPSPL




    subroutine cvpgen
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: dsp, dxp, dyp
        integer :: ic, iel, ip
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Sets control-point locations on paneled airfoil.
        !     Sets locations of point sources.
        !---------------------------------------------------------
        !
        if (ldbg) write (*, *) 'Setting control points for elements'
        !
        do iel = 1, nex
            icfrst(iel) = 0
            iclast(iel) = 0
        enddo
        !
        do ic = 1, icx
            ipco(ic) = 0
            ipcp(ic) = 0
        enddo
        !
        !
        !---- initialize counters for pointer accumulation
        ic = 0
        !
        do iel = 1, nel
            if (netype(iel)==3 .or. netype(iel)==4) then
                !------- ring or point singularity... no control points
                ip = ipfrst(iel)
                !
                !------- set limits so as to skip over IC do-loops
                icfrst(iel) = 0
                iclast(iel) = -1
                !
                cycle
            endif
            !
            !------ first control point in element IEL
            icfrst(iel) = ic + 1
            !
            !------ go over all panels on element IEL
            do ip = ipfrst(iel), iplast(iel) - 1
                ic = ic + 1
                if (ic>ncx) stop 'CVPGEN: Array overflow on NCX'
                !
                dxp = xp(ip + 1) - xp(ip)
                dyp = yp(ip + 1) - yp(ip)
                dsp = sqrt(dxp**2 + dyp**2)
                !
                if (dsp==0.0) then
                    !--------- zero-length panel gets dummy control point
                    ictype(ic) = -1
                    !
                else
                    !--------- ordinary panel
                    ictype(ic) = 0
                    if (iel2ir(iel)==1) then
                        ictype(ic) = 0
                    elseif (iel2ir(iel)>1 .and. iel2ir(iel)<nrp) then
                        !c             ICTYPE(IC) = -2
                    elseif (iel2ir(iel)==nrp) then
                        ictype(ic) = 0
                    endif
                    !
                endif
                !
                !-------- panel nodes defining sheet strength at IC
                ipco(ic) = ip
                ipcp(ic) = ip + 1
            enddo
            !
            !------ last control point in element IEL
            iclast(iel) = ic
            !
        enddo
        !
        !---- total number of control points
        nctot = ic
        !
        !---- compute normal vectors at vortex nodes, compute geometric quantities
        do iel = 1, nel
            !------ don't do point singularities
            if (netype(iel)==0 .or. netype(iel)==1 .or. netype(iel)       &
                    & ==2 .or. netype(iel)==5 .or. netype(iel)==6 .or.          &
                    & netype(iel)==7) then
                call xycset(iel)
                call anpset(iel)
            endif
        enddo
        !
        !---- invalidate any existing solution
        lqaic = .false.
        lsysp = .false.
        lgsys = .false.
        lgamu = .false.
        lgama = .false.
        !
        lncvp = .true.
        !
    end subroutine cvpgen
    !*==XYCSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CVPGEN



    subroutine xycset(iel)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iel
        !
        ! Local variables
        !
        real, save :: dcfrac
        real :: dsca, dsp, dspinv, dxp, dyp, stan, xtan, ytan
        integer :: ic, ic1, ic2, ice, ip, ip1, ip2, ipo, ipp
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Sets control points and their unit normal vectors
        !     for element IEL.
        !----------------------------------------------------------------
        !
        !c    DATA DCFRAC / 0.1 /
        data dcfrac/0.05/
        !c    DATA DCFRAC / 0.02 /
        !
        do ic = icfrst(iel), iclast(iel)
            ipo = ipco(ic)
            ipp = ipcp(ic)
            !
            dxp = xp(ipp) - xp(ipo)
            dyp = yp(ipp) - yp(ipo)
            dsp = sqrt(dxp**2 + dyp**2)
            !
            if (dsp==0.0) then
                !------- zero-length panel
                dspinv = 0.0
            else
                !------- ordinary panel
                dspinv = 1.0 / dsp
            endif
            !
            !------ assign control point at center of panel
            xc(ic) = xp(ipo) + 0.5 * dxp
            yc(ic) = yp(ipo) + 0.5 * dyp
            dsc(ic) = dsp
            dsc_dxy(1, ic) = dxp * dspinv
            dsc_dxy(2, ic) = dyp * dspinv
            !
            !------ normal vector at control point
            anc(1, ic) = dyp * dspinv
            anc(2, ic) = -dxp * dspinv
            !
            !------ geometry sensitivities of normal vector for design/optimization stuff
            anc_dxy(1, 1, ic) = anc(1, ic) * anc(2, ic) * dspinv
            anc_dxy(1, 2, ic) = anc(2, ic)**2 * dspinv
            anc_dxy(2, 1, ic) = -anc(1, ic)**2 * dspinv
            anc_dxy(2, 2, ic) = -anc(2, ic) * anc(1, ic) * dspinv
        enddo
        !
        !
        !------ set up extra QNDOF control points for each body element
        if (.not.lbody(iel)) then
            ice = ncx + iel
            xc(ice) = 0.
            yc(ice) = 0.
            !
        elseif (lxbod(iel)) then
            ic = iclast(iel)
            ip = iplast(iel)
            !
            ice = ncx + iel
            xc(ice) = xp(ip) - dcfrac * dsc(ic)
            yc(ice) = 0.
            anc(1, ice) = 1.0
            anc(2, ice) = 0.
            !
        else
            ic1 = icfrst(iel)
            ic2 = iclast(iel)
            !
            ip1 = ipfrst(iel)
            ip2 = iplast(iel)
            !
            xtan = xp(ip1) - xp(ip1 + 1) + xp(ip2) - xp(ip2 - 1)
            ytan = yp(ip1) - yp(ip1 + 1) + yp(ip2) - yp(ip2 - 1)
            stan = sqrt(xtan**2 + ytan**2)
            !
            dsca = 0.5 * (dsc(ic1) + dsc(ic2))
            !
            ice = ncx + iel
            xc(ice) = 0.5 * (xp(ip1) + xp(ip2)) - dcfrac * dsca * xtan / stan
            yc(ice) = 0.5 * (yp(ip1) + yp(ip2)) - dcfrac * dsca * ytan / stan
            anc(1, ice) = xtan / stan
            anc(2, ice) = ytan / stan
            !
        endif
        !
    end subroutine xycset
    !*==ANPSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! XYCSET



    subroutine anpset(iel)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iel
        !
        ! Local variables
        !
        real :: ds, dsq, dsqm, dsqp, dx, dxm, dxp, dx_dxm, &
                & dx_dxp, dx_dym, dx_dyp, dy, dym, dyp, dy_dxm, &
                & dy_dxp, dy_dym, dy_dyp, tmp1, tmp2
        integer :: ip, ipm, ipp
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Sets unit normal vectors on vortex node locations
        !     for element IEL.
        !----------------------------------------------------------------
        !
        !---- don't process single-point elements
        if (ipfrst(iel)==iplast(iel)) return
        !
        do ip = ipfrst(iel), iplast(iel)
            !
            !------ examine the two panels IPM..IP, IP..IPP adjoining vortex node IP
            ipm = max(ip - 1, ipfrst(iel))
            ipp = min(ip + 1, iplast(iel))
            !
            dxm = xp(ip) - xp(ipm)
            dym = yp(ip) - yp(ipm)
            !
            dxp = xp(ipp) - xp(ip)
            dyp = yp(ipp) - yp(ip)
            !
            dsqm = dxm * dxm + dym * dym
            dsqp = dxp * dxp + dyp * dyp
            !
            if (dsqm==0.0) then
                !------- first point on element, or panel IPM..IP has zero length
                dx = dxp
                dy = dyp
                dx_dxm = 0.
                dx_dym = 0.
                dx_dxp = 1.0
                dx_dyp = 0.
                dy_dxm = 0.
                dy_dym = 0.
                dy_dxp = 0.
                dy_dyp = 1.0
            elseif (dsqp==0.0) then
                !------- last point on element, or panel IP..IPP has zero length
                dx = dxm
                dy = dym
                dx_dxm = 1.0
                dx_dym = 0.
                dx_dxp = 0.
                dx_dyp = 0.
                dy_dxm = 0.
                dy_dym = 1.0
                dy_dxp = 0.
                dy_dyp = 0.
            else
                !------- usual interior point... normal vector is weighted average
                dx = dxm * dsqp + dxp * dsqm
                dy = dym * dsqp + dyp * dsqm
                dx_dxm = dsqp + dxp * 2.0 * dxm
                dx_dym = +dxp * 2.0 * dym
                dx_dxp = dsqm + dxm * 2.0 * dxp
                dx_dyp = +dxm * 2.0 * dyp
                dy_dxm = +dyp * 2.0 * dxm
                dy_dym = dsqp + dyp * 2.0 * dym
                dy_dxp = +dym * 2.0 * dxp
                dy_dyp = dsqm + dym * 2.0 * dyp
            endif
            !
            dsq = dx * dx + dy * dy
            if (dsq==0.0) then
                !------- internal error
                write (*, *) '? ANPSET: zero avg panel length. IEL IP =', &
                        & iel, ip
                return
            endif
            !
            !------ set unit normal ANP
            ds = sqrt(dsq)
            anp(1, ip) = dy / ds
            anp(2, ip) = -dx / ds
            !
            tmp1 = -anp(1, ip) / dsq
            anp_xym(1, 1, ip) = -dy_dxm / ds - tmp1 * (dx * dx_dxm + dy * dy_dxm)
            anp_xym(1, 2, ip) = -dy_dym / ds - tmp1 * (dx * dx_dym + dy * dy_dym)
            !
            anp_xyo(1, 1, ip) = dy_dxm / ds + tmp1 * (dx * dx_dxm + dy * dy_dxm)       &
                    & - dy_dxp / ds - tmp1 * (dx * dx_dxp + dy * dy_dxp)
            anp_xyo(1, 2, ip) = dy_dym / ds + tmp1 * (dx * dx_dym + dy * dy_dym)       &
                    & - dy_dyp / ds - tmp1 * (dx * dx_dyp + dy * dy_dyp)
            !
            anp_xyp(1, 1, ip) = dy_dxp / ds + tmp1 * (dx * dx_dxp + dy * dy_dxp)
            anp_xyp(1, 2, ip) = dy_dyp / ds + tmp1 * (dx * dx_dyp + dy * dy_dyp)
            !
            tmp2 = -anp(2, ip) / dsq
            anp_xym(2, 1, ip) = dx_dxm / ds - tmp2 * (dx * dx_dxm + dy * dy_dxm)
            anp_xym(2, 2, ip) = dx_dym / ds - tmp2 * (dx * dx_dym + dy * dy_dym)
            !
            anp_xyo(2, 1, ip) = -dx_dxm / ds + tmp2 * (dx * dx_dxm + dy * dy_dxm)      &
                    & + dx_dxp / ds - tmp2 * (dx * dx_dxp + dy * dy_dxp)
            anp_xyo(2, 2, ip) = -dx_dym / ds + tmp2 * (dx * dx_dym + dy * dy_dym)      &
                    & + dx_dyp / ds - tmp2 * (dx * dx_dyp + dy * dy_dyp)
            !
            anp_xyp(2, 1, ip) = -dx_dxp / ds + tmp2 * (dx * dx_dxp + dy * dy_dxp)
            anp_xyp(2, 2, ip) = -dx_dyp / ds + tmp2 * (dx * dx_dyp + dy * dy_dyp)
        enddo
        !
    end subroutine anpset
    !*==ITPSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ANPSET


    subroutine itpset(iel)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iel
        !
        ! Local variables
        !
        integer :: ip1, ip2
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------
        !     Sets TE panel endpoint indices
        !--------------------------------------
        !
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        !
        !---- default case... no TE panel
        ipte1(iel) = -999
        ipte2(iel) = -999
        !
        !---- no TE panel on closed axis body or non-solid body
        if (lxbod(iel) .or. netype(iel)/=0) return
        !
        !---- TE panel must be chosen...
        if (ltpan(iel)) then
            !----- usual TE panel closing off body
            ipte1(iel) = ip1
            ipte2(iel) = ip2
            if (yp(ip1)==0.0 .and. yp(ip2)>0.0) then
                ipte1(iel) = 0
                ipte2(iel) = ip2
            elseif (yp(ip2)==0.0 .and. yp(ip1)>0.0) then
                ipte1(iel) = ip1
                ipte2(iel) = 0
            endif
        endif
        !
    end subroutine itpset
end module m_geom
