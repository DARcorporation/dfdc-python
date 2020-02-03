module m_wakesubs
    implicit none
contains
    !*==WAKERESET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020

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

    subroutine wakereset
        use i_dfdc
        use m_dfdcsubs, only : updrotwak
        use m_inigrd, only : updgrd
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: dx, dx0, dy, dy0
        integer :: iel, ielo, ip, ip1, ip1o, ip2o, ipm, ipte
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Realigns the upper and lower wakes bounding the duct
        !     vortex wake grid with the current flow.
        !
        !     A new wake grid system is defined by updating and re-relaxing
        !     the grid.  All wake elements are re-initialized with the new
        !     geometry.
        !-------------------------------------------------------------
        !
        !---- CB wake
        iel = ir2iel(1)
        call wakmov(iel, isplot(iel))
        !
        !---- duct wake (check for tip gap...)
        if (tgap<=0.0) then
            iel = ir2iel(nrp)
            call wakmov(iel, isplot(iel))
        else
            !---- if tip gap move next inner wake (the one with circulation)
            iel = ir2iel(nrp - 1)
            call wakmov(iel, isplot(iel))
            !---- move outer wake streamline to match moved streamline
            ip1 = ipfrst(iel)
            ielo = ir2iel(nrp)
            ip1o = ipfrst(ielo)
            ip2o = iplast(ielo)
            ipte = ip1 + igtedw - 1
            dx0 = xp(ip1o) - xp(ipte)
            dy0 = yp(ip1o) - yp(ipte)
            !c        write(*,*) 'WAKERESET first dx0,dy0 ',dx0,dy0
            do ip = ip1o, ip2o
                ipm = ipte + (ip - ip1o)
                dx = xp(ip) - xp(ipm)
                dy = yp(ip) - yp(ipm)
                !c          write(*,*) 'WAKERESET ip,dx,dy ',ip,dx,dy
                yp(ip) = yp(ip) - (dy - dy0)
            enddo
        endif
        !
        !---- change grid boundaries, relax grid and update vortex wakes
        call updgrd
        call updrotwak
        !
        !---- invalidate existing solution
        lncvp = .false.
        lqaic = .false.
        lqgic = .false.
        lqcnt = .false.
        lgsys = .false.
        lgamu = .false.
        lgama = .false.
        !
    end subroutine wakereset
    !*==WAKMOV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! WAKERESET




    subroutine wakmov(iela, isq)
        use i_dfdc
        use m_geom, only : xycset, xypspl, anpset
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iela, isq
        !
        ! Local variables
        !
        integer :: iel, ip, ip1, ip2, ir
        real :: xp2
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Moves existing wake (element IELA) by
        !     aligning panels with current velocities.
        !     Panel side flag ISQ indicates which velocity to use.
        !       ISQ=0 uses mean surface velocities QC
        !       ISQ=1 uses right surface velocities QCR
        !       ISQ=2 uses left  surface velocities QCL
        !
        !     Panel points are moved preserving existing panel lengths
        !     while aligning the panels with the current velocity
        !     vectors on the panel.
        !-------------------------------------------------------------
        !
        if (netype(iela)==0) then
            !----- this element is a body panel
            write (*, *) 'WAKMOV: Element', iela, &
                    &' is body, cannot be moved!'
            return
        endif
        !
        iel = iela
        !
        !---- Save X location for end of wake
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        xp2 = xp(ip2)
        !
        call wakmovr(iel, xp, yp, isq)
        !
        !---- Reset last point to initial (unchanged) XWAKE location
        xp(ip2) = xp2
        !
        !---- respline node powgrints
        call xypspl(iel)
        !---- reset control points and normal vectors
        call xycset(iel)
        call anpset(iel)
        !
        !      WRITE(*,*) 'Generating dq/d(gam,sig) influence of moved wake...'
        !      CALL QAIC1(.FALSE., 1,NCTOT, IEL,IEL)
        !      DO KEL = 1, NEL
        !        ICE = NCX + KEL
        !        IF(LBODY(KEL)) CALL QAIC1(.FALSE.,ICE,ICE,IEL,IEL)
        !      ENDDO
        !
        !      IP1 = IPFRST(IEL)
        !      IP2 = IPLAST(IEL)
        !      CALL GSYS(.FALSE., .FALSE., .TRUE., IP1,IP2, 1,0)
        !
        !---- free existing unit wake strength vectors (if any), mark vectors invalid
        do ip = ip1, ip2
            luset(iugam(ip)) = .false.
            luset(iusig(ip)) = .false.
            iugam(ip) = 0
            iusig(ip) = 0
        enddo
        do ir = 1, nrp
            luset(iuvwk(ir)) = .false.
            iuvwk(ir) = 0
        enddo
        !
    end subroutine wakmov
    !*==WAKMOVR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! WAKMOV



    subroutine wakmovr(iel, xpnew, ypnew, isq)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iel, isq
        real, dimension(*) :: xpnew, ypnew
        !
        ! Local variables
        !
        real, dimension(ipx) :: dsp
        real :: dxp, dyp, qmag, qstag, qx, qy, un, vn
        integer :: ic, ic1, ic2, ip, ip1, ip2, ipo, ipp
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------
        !---- Routine to move points on element to new positions
        !     to align with velocity vector on side ISQ (0/1/2)
        !----------------------------------------------------------
        !
        !
        !---- speed below this is treated as stagnation
        qstag = 0.0001 * qref
        !
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        !
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        !
        !---- calculate velocity and its Jacobians at all wake points
        !cc      CALL QAIC1(.FALSE., IC1,IC2, 1,NEL)
        !---- set velocities on both sides of panels at all control points
        !cc      CALL QQCALC(NCTOT,ANC, IPCO,IPCP, GAM,SIG, QC, QCL,QCR )
        !
        do ip = ip1, ip2 - 1
            dxp = xp(ip + 1) - xp(ip)
            dyp = yp(ip + 1) - yp(ip)
            dsp(ip) = sqrt(dxp**2 + dyp**2)
        enddo
        !
        !---- march down wake setting wake points
        ip = ipco(ic1)
        xpnew(ip) = xp(ip)
        ypnew(ip) = yp(ip)
        !
        do ic = ic1, ic2
            ipo = ipco(ic)
            ipp = ipcp(ic)
            !
            if (isq==0) then
                qx = qc(1, ic)
                qy = qc(2, ic)
            elseif (isq==1) then
                qx = qcr(1, ic)
                qy = qcr(2, ic)
            elseif (isq==2) then
                qx = qcl(1, ic)
                qy = qcl(2, ic)
            else
                qx = qc(1, ic)
                qy = qc(2, ic)
            endif
            !
            qmag = sqrt(qx**2 + qy**2)
            !------- stagnation... don't try to move points
            if (qmag<qstag) return
            !
            un = qx / qmag
            vn = qy / qmag
            xpnew(ipp) = xpnew(ipo) + un * dsp(ipo)
            ypnew(ipp) = ypnew(ipo) + vn * dsp(ipo)
            !
        enddo
        !
    end subroutine wakmovr
end module m_wakesubs
