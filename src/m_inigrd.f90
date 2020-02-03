module m_inigrd
    implicit none
contains
    !*==INIGRD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! UVGRD
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


    subroutine inigrd
        use i_dfdc
        use m_userio, only : askr
        use m_xutils, only : setexp
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: amass, amass1, dmass1, ds, ds1, dx, dxw, dy, &
                & frac, frc, smax, xcb, xdw, xl, xlast, xtecb, &
                & xtedw, xu, ycb, ydw, yl, ylast, yu, yy
        integer :: i, iel, ilwr, ip, ip1cb, ip1dw, ip2cb, ip2dw, &
                & ipc, ipcb, ipd, ipdw, ipup, ite, iupr, j, k, &
                & nr
        real, dimension(ipx) :: s
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Sets up vortex wake grid system from rotor and wake
        !     geometry
        !------------------------------------------------------
        !
        if (ldbg) then
            write (*, *) ' '
            write (*, *) 'INIGRD setup grid for slipstream'
            !
            write (lundbg, *) ' '
            write (lundbg, *) 'INIGRD setup grid for slipstream'
        endif
        !
        !---- Find most-upstream rotor point on duct wall and set as start of grid
        nrupstrm = 1
        ipup = iprotdw(nrupstrm)
        do nr = 2, nrotor
            if (iprotdw(nr)<ipup) then
                nrupstrm = nr
                ipup = iprotdw(nr)
            endif
        enddo
        !---- Setup rotor at upstream end of grid
        do nr = 1, nrotor
            igrotor(nr) = 1 + iprotdw(nr) - ipup
        enddo
        !
        nr = nrupstrm
        !---- Indices of rotor line on CB and duct define upstream grid boundary
        !---- find axial location of rotor and TE's on cb and duct walls
        iel = 1
        ip1cb = ipfrst(iel)
        ip2cb = iplast(iel)
        ipcb = iprotcb(nr)
        if (ipcb==0) then
            write (*, *) 'Error locating rotor on CB wall'
            stop
        endif
        if (ldbg) write (lundbg, *) 'Rotor on CB @ ', ipcb
        xcb = xp(ipcb)
        ycb = yp(ipcb)
        xtecb = xp(ip1cb)
        !c    IF(LDBG) WRITE(*,*) 'XR, XTECB ',XCB,XTECB
        !
        iel = 2
        ip1dw = ipfrst(iel)
        ip2dw = iplast(iel)
        ipdw = iprotdw(nr)
        if (ipdw==0) then
            write (*, *) 'Error locating rotor on Duct wall'
            stop
        endif
        if (ldbg) write (lundbg, *) 'Rotor on Duct wall @ ', ipdw
        xdw = xp(ipdw)
        ydw = yp(ipdw)
        xtedw = xp(ip2dw)
        !c    IF(LDBG) WRITE(*,*) 'XR XTEDW ',XDW,XTEDW
        !
        !
        !---- Define lower streamline (J=1)
        i = igrotor(nr) - 1
        j = 1
        do ip = ipcb, ip1cb, -1
            i = i + 1
            xg(i, j) = xp(ip)
            yg(i, j) = yp(ip)
            ip2ig(ip) = i
        enddo
        xlast = xg(i, j)
        ylast = yg(i, j)
        igtecb = i
        !c      write(*,*) 'IGTECB,II defined as ',IGTECB,II
        !c    IF(LDBG) WRITE(*,*) 'L IGTECB,XYLAST ',IGTECB,XLAST,YLAST
        !
        !---- Add points at X values of duct wall if it extends downstream of CB
        ipd = ipdw + igtecb - 1
        if (ipd<ip2dw) then
            do ip = ipd + 1, ip2dw
                i = i + 1
                xg(i, j) = xp(ip)
                yg(i, j) = ylast
            enddo
            xlast = xg(i, j)
            ylast = yg(i, j)
            ite = i
            !c     IF(LDBG) WRITE(*,*) 'L ITE,XYLAST ',ITE,XLAST,YLAST
        endif
        !
        !---- Check downstream wake end location
        if (ldbg) write (lundbg, *) 'INIGRD XWAKE ', xwake
        if (xdwklen<0.0 .or. xdwklen>5.0) then
            write (*, *) 'Rotor wake length out of bounds ', xdwklen
            write (*, *) '  CB   TE at X,Y =', xbte(1), ybte(1)
            write (*, *) '  Duct TE at X,Y =', xbte(2), ybte(2)
            write (*, *) '  Rotor Diameter =', 2.0 * rtip(1)
            call askr('Enter length of wake from TE dXlen/D^', xdwklen)
        endif
        xwake = max(xbte(1), xbte(2)) + xdwklen * 2.0 * rtip(1)
        !
        !---- Add points downstream to end of wake
        dx = xg(i, j) - xg(i - 1, j)
        dy = yg(i, j) - yg(i - 1, j)
        ds = sqrt(dx * dx + dy * dy)
        !
        smax = xwake - xlast
        dxw = smax / float(nwake - 1)
        ds1 = 0.2 * dxw
        call setexp(s, ds1, smax, nwake)
        do k = 2, nwake
            i = i + 1
            xg(i, j) = xlast + s(k)
            yg(i, j) = ylast
        enddo
        ilwr = i
        !
        !
        !---- Define upper streamline (J=JJ=NRP from rotor line points)
        jj = nrp
        !
        i = igrotor(nr) - 1
        j = jj
        do ip = ipdw, ip2dw
            i = i + 1
            xg(i, j) = xp(ip)
            yg(i, j) = yp(ip)
            ip2ig(ip) = i
        enddo
        xlast = xg(i, j)
        ylast = yg(i, j)
        igtedw = i
        !c      write(*,*) 'IGTEDW,II defined as ',IGTEDW,II
        !c    IF(LDBG) WRITE(*,*) 'U IGTEDW,XYLAST ',IGTEDW,XLAST,YLAST
        !
        !---- Add points at X values of CB wall if it extends downstream of duct
        ipc = ipcb - igtedw + 1
        if (ipc>ip1cb) then
            do ip = ipc - 1, ip1cb, -1
                i = i + 1
                xg(i, j) = xp(ip)
                yg(i, j) = ylast
            enddo
            xlast = xg(i, j)
            ylast = yg(i, j)
            ite = i
            !c     IF(LDBG) WRITE(*,*) 'U ITE,XYLAST ',ITE,XLAST,YLAST
        endif
        !
        !---- Add points downstream to end of wake
        dx = xg(i, j) - xg(i - 1, j)
        dy = yg(i, j) - yg(i - 1, j)
        ds = sqrt(dx * dx + dy * dy)
        !
        smax = xwake - xlast
        dxw = smax / float(nwake - 1)
        ds1 = 0.2 * dxw
        call setexp(s, ds1, smax, nwake)
        do k = 2, nwake
            i = i + 1
            xg(i, j) = xlast + s(k)
            yg(i, j) = ylast
        enddo
        iupr = i
        !
        if (ilwr/=iupr) then
            write (*, *) 'INIGRD: incompatible # points on J=1 and J=JJ', &
                    & ilwr, iupr
            write (*, *) 'IP1CB,IP2CB,IPCB ', ip1cb, ip2cb, ipcb
            write (*, *) 'IP1DW,IP2DW,IPDW ', ip1dw, ip2dw, ipdw
            stop
        endif
        !---- Set total # of streamwise points
        ii = iupr
        !
        !
        !---- Add points at inlet and outlet boundaries
        !---- Inlet points are simply the rotor line points XR,YR
        i = 1
        xgmin = xrp(1, nr)
        do j = 1, jj
            xg(i, j) = xrp(j, nr)
            yg(i, j) = yrp(j, nr)
            xgmin = min(xgmin, xg(i, j))
        enddo
        !
        !---- Sweep downstream setting intermediate J lines preserving relative
        !     massflow between upper and lower walls
        amass1 = pi * (yg(1, jj)**2 - yg(1, 1)**2)
        do i = 2, ii
            xl = xg(i, 1)
            yl = yg(i, 1)
            xu = xg(i, jj)
            yu = yg(i, jj)
            amass = pi * (yu**2 - yl**2)
            !---- Create radial spacing by preserving relative area
            do j = 2, jj - 1
                dmass1 = pi * (yg(1, j)**2 - yg(1, j - 1)**2)
                frc = dmass1 / amass1
                yy = sqrt(frc * amass / pi + yg(i, j - 1)**2)
                frac = (yy - yl) / (yu - yl)
                xg(i, j) = xl + frac * (xu - xl)
                yg(i, j) = yy
            enddo
        enddo
        !
        !---- Done, initial grid defined!
        if (ldbg) then
            write (*, *) 'INIGRD grid defined with II,JJ ', ii, jj
            write (lundbg, *) 'INIGRD grid defined with II,JJ ', ii, jj
        endif
        !
        !     Add initial grid smoothing relaxation?
        !---- Grid smoothing by elliptic solver
        call rlxgrd
        !
        !---- invalidate existing solution to incorporate new vortex wakes
        lncvp = .false.
        lqaic = .false.
        lqgic = .false.
        lqcnt = .false.
        lgsys = .false.
        lgamu = .false.
        lgama = .false.
        !
    end subroutine inigrd
    !*==UPDGRD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! INIGRD


    subroutine updgrd
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: i, iel, ip, ip1, ip2, j
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Updates grid upper and lower boundaries with
        !     points from current upper and lower vortex wakes
        !------------------------------------------------------
        !
        if (ldbg) then
            write (*, *) ' '
            write (*, *) 'Updating grid system boundaries'
        endif
        !
        !---- Redefine lower streamline (J=1) with IR=1 wake points
        j = 1
        iel = ir2iel(1)
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        do ip = ip1, ip2
            i = igtecb + (ip - ip1)
            xg(i, j) = xp(ip)
            yg(i, j) = yp(ip)
        enddo
        if (ldbg) write (lundbg, *)                                      &
                &'Updated lower boundary from element '&
                &, iel
        !
        !---- Redefine upper streamline (J=JJ) with IR=NRP wake points
        j = jj
        iel = ir2iel(nrp)
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        do ip = ip1, ip2
            i = igtedw + (ip - ip1)
            xg(i, j) = xp(ip)
            yg(i, j) = yp(ip)
        enddo
        if (ldbg) write (lundbg, *)                                      &
                &'Updated upper boundary from element '&
                &, iel
        !
        !---- Grid boundaries redefined!
        if (ldbg) write (*, *) 'UPDGRD grid boundaries updated'
        !
        !---- Grid smoothing by elliptic solver
        call rlxgrd
        !
    end subroutine updgrd
    !*==RLXGRD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! UPDGRD




    subroutine rlxgrd
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: i, itmaxs, j, kbcbot, kbcinl, kbcout, kbctop
        real :: toler
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Runs elliptic grid solution to smooth grid
        !------------------------------------------------------
        !
        if (ii<=0 .or. jj<=0) then
            write (*, *) 'RLXGRD:  no grid defined II,JJ ', ii, jj
            return
        endif
        !
        !---- Set up streamfunction array
        ypos(1) = 0.0
        do j = 2, jj
            ypos(j) = yg(1, j)**2 - yg(1, j - 1)**2 + ypos(j - 1)
        enddo
        !---- Set up streamwise spacing array
        do i = 1, ii
            xpos(i) = xg(i, 1)
        enddo
        !
        !---- Run grid solver
        itmaxs = -400
        toler = 1.0e-9
        kbcinl = 0 ! inlet plane Dirichlet, fix nodes
        kbcout = 1 ! outlet plane Neumann, move nodes on boundary
        kbcbot = 0 ! bottom streamline Dirichlet, fix nodes
        kbctop = 0 ! top streamline Dirichlet, fix nodes
        !      WRITE(*,*)
        !      WRITE(*,*) 'Grid relaxation iteration'
        call axell(ix, ii, jj, xg, yg, xpos, ypos, itmaxs, toler, kbcinl, kbcout, &
                & kbcbot, kbctop)
        !
    end subroutine rlxgrd
    !*==SETGRDFLW.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! RLXGRD

    !     &  XG(IX,JX),    YG(IX,JX),
    !     &  QG(IX,JX),    QXG(IX,JX),   QYG(IX,JX),  QTG(IX,JX),
    !     &  RG(IX,JX),    PG(IX,JX),    POG(IX,JX),
    !     &  BGAMG(IX,JX), DSG(IX,JX),   DHG(IX,JX),
    !     &  XPOS(IX),     YPOS(JX)


    subroutine setgrdflw
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: dbgamg, ddhg, ddsg, sigc, vmc
        integer :: i, ic, ic1, ic2, iel, ig, ip1, ip2, ir, j, &
                & n, nnc
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Sets grid flow data from circulation and entropy on
        !     rotor lines
        !---------------------------------------------------------
        !
        !---- Clear grid data for accumulation over sources
        call clrgrdflw
        !
        !---- Set B*GAM circulation, delta entropy, and delta enthalpy on gridC
        !     from rotor and drag sources

        do n = 1, nrotor
            !
            iel = ielrotor(n)
            ic1 = icfrst(iel)
            ic2 = iclast(iel)
            !---- Set values at rotor I station due to blade row
            ig = igrotor(n)
            !
            nnc = ic2 - ic1 + 1
            if (nnc/=jj - 1) then
                write (*, *) 'Rotor source element-grid mismatch ', nnc, &
                        & jj - 1
                stop
            endif
            !
            do j = 1, jj - 1
                ir = j
                !---- Circulation and enthalpy
                dbgamg = bgam(ir, n)
                ddhg = omega(n) * bgam(ir, n) * pi2i
                !---- Entropy from drag sources on rotor
                ic = ic1 + ir - 1
                ip1 = ipco(ic)
                ip2 = ipcp(ic)
                vmc = sqrt(qc(1, ic)**2 + qc(2, ic)**2)
                sigc = 0.5 * (sig(ip1) + sig(ip2))
                ddsg = vmc * sigc
                !
                !---- Convect values downstream from rotor
                do i = ig, ii - 1
                    bgamg(i, j) = bgamg(i, j) + dbgamg
                    dhg(i, j) = dhg(i, j) + ddhg
                    dsg(i, j) = dsg(i, j) + ddsg
                    !c        RG(I,J)    = RG(1,J)
                enddo
                !
            enddo
        enddo
        !
    end subroutine setgrdflw
    !*==CLRGRDFLW.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETGRDFLW


    subroutine clrgrdflw
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: i, j
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Clears (initializes to zero) grid flow BGAMG,DHG,DSG
        !     This is normally done before accumulating disk
        !     contributions
        !---------------------------------------------------------
        !
        !---- Clear grid data for accumulation of circulation, enthalpy, entropy
        do j = 1, jj
            do i = 1, ii
                bgamg(i, j) = 0.0
                dhg(i, j) = 0.0
                dsg(i, j) = 0.0
                rg(i, j) = rho
            enddo
        enddo
        !
    end subroutine clrgrdflw
    !*==ROTBG2GRD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CLRGRDFLW


    subroutine rotbg2grd(n)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        !
        ! Local variables
        !
        real :: dbgamg, ddhg
        integer :: i, ig, ir, j
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Updates grid flow (circulation) from BGAM on rotor N
        !     rotor lines
        !---------------------------------------------------------
        !
        !---- Set B*GAM circulation, delta enthalpy on grid from blade/rotors
        !---- Set values at grid IG station due to blade row
        ig = igrotor(n)
        !
        do j = 1, jj - 1
            ir = j
            !---- Circulation and enthalpy
            dbgamg = bgam(ir, n)
            ddhg = omega(n) * bgam(ir, n) * pi2i
            !
            !---- Convect values downstream from rotor
            do i = ig, ii - 1
                bgamg(i, j) = bgamg(i, j) + dbgamg
                dhg(i, j) = dhg(i, j) + ddhg
            enddo
            !
        enddo
        !
    end subroutine rotbg2grd
    !*==GETGRDFLW.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ROTBG2GRD



    subroutine getgrdflw
        use i_dfdc
        use m_grdutils, only : uvgrdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: i, io, ip, j, jo, jp
        real :: xcg, xoo, xop, xpo, xpp, ycg, yoo, yop, ypo, &
                & ypp
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Sets grid initial flow data from rotor line conditions
        !---------------------------------------------------------
        !
        call uvgrdc(ix, ii, jj, xg, yg, xpos, ypos, qxg, qyg)
        !---- Print grid velocities
        write (44, *) 'I,J, XC,YC, U, V'
        do i = 1, ii - 1
            io = i
            ip = i + 1
            do j = 1, jj - 1
                jo = j
                jp = j + 1
                xoo = xg(io, jo)
                xpo = xg(ip, jo)
                xop = xg(io, jp)
                xpp = xg(ip, jp)
                yoo = yg(io, jo)
                ypo = yg(ip, jo)
                yop = yg(io, jp)
                ypp = yg(ip, jp)
                xcg = 0.25 * (xoo + xpo + xop + xpp)
                ycg = 0.25 * (yoo + ypo + yop + ypp)
                write (44, 100) i, j, xcg, ycg, qxg(i, j), qyg(i, j)
            enddo
        enddo
        100  format (2i5, 5(1x, f11.6))
        !
    end subroutine getgrdflw
    !*==INLSFCN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GETGRDFLW



    subroutine inlsfcn
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: j, nr
        real :: vx, yscl
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Calculates streamfunction (YPOS) at grid inlet using
        !     velocities at rotor and current grid points
        !---------------------------------------------------------
        !
        nr = nrupstrm
        !---- Set up streamfunction array using grid and rotor line velocity
        ypos(1) = 0.0
        do j = 1, jj
            vx = vabs(1, j - 1, nr)
            ypos(j) = 0.5 * vx * (yg(1, j)**2 - yg(1, j - 1)**2) + ypos(j - 1)
            !c        YPOS(J) = 0.5*VX*YG(1,J)**2
        enddo
        yscl = 1.0 / ypos(jj)
        !c      DO J = 1, JJ
        !c        write(*,*) 'j ypos ',j,ypos(j)
        !c        YPOS(J) = YSCL*YPOS(J)
        !c      END DO
        !
    end subroutine inlsfcn
    !*==RLXGRD2.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! INLSFCN


    subroutine rlxgrd2
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: i, itmaxs, kbcbot, kbcinl, kbcout, kbctop
        real :: toler
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Runs elliptic grid solution to smooth grid
        !------------------------------------------------------
        !
        if (ii<=0 .or. jj<=0) then
            write (*, *) 'RLXGRD:  no grid defined II,JJ ', ii, jj
            return
        endif
        !
        !---- Set up streamfunction array
        call inlsfcn
        !
        !---- Set up streamwise spacing array
        do i = 1, ii
            xpos(i) = xg(i, 1)
        enddo
        !
        !---- Run grid solver
        itmaxs = -400
        toler = 1.0e-9
        kbcinl = 0 ! inlet plane Dirichlet, fix nodes
        kbcout = 1 ! outlet plane Neumann, move nodes on boundary
        kbcbot = 0 ! bottom streamline Dirichlet, fix nodes
        kbctop = 0 ! top streamline Dirichlet, fix nodes
        call axell(ix, ii, jj, xg, yg, xpos, ypos, itmaxs, toler, kbcinl, kbcout, &
                & kbcbot, kbctop)
        !
    end subroutine rlxgrd2
    !*==AXELL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! RLXGRD2


    subroutine axell(ix, ii, jj, x, y, xpos, ypos, itmaxs, toler, kbcinl, &
            & kbcout, kbcbot, kbctop)
        use m_spline, only : deval, spline, seval, scalc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: idim = 500
        !
        ! Dummy arguments
        !
        integer :: ii, itmaxs, ix, jj, kbcbot, kbcinl, kbcout, &
                & kbctop
        real :: toler
        real, dimension(ix, *) :: x, y
        real, dimension(*) :: xpos, ypos
        !
        ! Local variables
        !
        real :: a, ad1, ad2, ainv, aja, alf, b, bet, cetm, &
                & cetp, cxim, cxip, detav, detl, detm, detp, detq, &
                & dmax, ds, dset1, dset2, dset3, dx, dxde1, dxde2, &
                & dxdet, dxdx1, dxdx2, dxdxi, dxiav, dxil, dxim, &
                & dxip, dxiq, dxy, dy, dyde1, dyde2, dydet, dydx1, &
                & dydx2, dydxi, gam, rez, rlx, rlx1, rlx2, rlx3, &
                & srlx, xdse, xlo, xmm, xmo, xmp, xol, xom, xoo, &
                & xop, xoq, xpm, xpo, xpp, xqo, xs, xsn, yaa, &
                & yam, yap, ydse, ylo, ymm, ymo, ymp, yol, yom, &
                & yoo, yop, yoq, ypm, ypo
        real, dimension(idim) :: c, sb, sbi, st, sti, xb, xbs, &
                & xt, xts, yb, ybs, yt, yts
        real, dimension(2, idim) :: d
        integer :: iback, ifin, il, im, io, ip, ipass, iq, &
                & itmax, jl, jm, jo, jp, jq
        real :: ypp, yqo, ys, ysn, z_s, z_xoo, z_yoo
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Axisymmetric elliptic grid smoother.
        !     Uses SLOR implicit along i direction.
        !
        !   Input:
        !     IX        first grid-array dimension
        !     II,JJ     grid i,j size
        !     X(i,j)    grid z coordinates
        !     Y(i,j)    grid r coordinates
        !     XPOS(i)   grid s values
        !     YPOS(j)   grid e values  (axisymmetric streamfunction)
        !                ( YPOS(j) = (j-1)**2 gives uniform spacing in r )
        !     ITMAXS    max number of SLOR smoothing passes
        !               if ITMAXS < 0, overrelaxation is NOT used
        !     TOLER     convergence tolerance (movement distance)
        !     KBCINL    1 = Neumann, 0 = Dirichlet   BC applied at i=1
        !     KBCOUT    1 = Neumann, 0 = Dirichlet   BC applied at i=II
        !     KBCBOT    1 = Neumann, 0 = Dirichlet   BC applied at j=1
        !     KBCTOP    1 = Neumann, 0 = Dirichlet   BC applied at j=JJ
        !
        !   Output:
        !     X(i,j)    smoothed-grid z coordinates
        !     Y(i,j)    smoothed-grid r coordinates
        !
        !.............................................................
        !
        !   Solves Thompson's grid-generation equations with source term
        !   which makes e(x,r) satisfy the axisymmetric streamfunction
        !   equation.  Hence,  e = constant  lines are streamlines.
        !
        !       s_zz + s_rr = 0
        !       e_zz + e_rr = e_r / r
        !
        !   These equations for s(z,r), e(z,r) are inverted to
        !
        !       a z_ss  -  2b z_se  +  c z_ee  =  (J/r) z_e z_s
        !
        !       a r_ss  -  2b r_se  +  c r_ee  =  (J/r) r_e z_s
        !
        !   where
        !                    2        2                    2        2
        !         a  =  (z_e)  + (r_e)           c  = (z_s)  + (r_s)
        !
        !         b  =  z_e z_s  +  r_e r_s
        !
        !         J  =  z_s r_e  -  z_e r_s
        !
        !
        !   The equivalent form which is actually differenced is
        !
        !       a z_ss  -  2b z_se  +  c/r (r z_e)_e  -  b/r r_s z_e  =  0
        !
        !       a r_ss  -  2b r_se  +  c/r (r r_e)_e  -  b/r r_s r_e  =  0
        !
        !-------------------------------------------------------------
        !
        !
        if (ii>idim) stop 'ELLIP: Array overflow.  Increase IDIM.'
        !
        itmax = iabs(itmaxs)
        !
        !---- convergence tolerance for each phase
        dset1 = 1.0e-1
        dset2 = 5.0e-3
        dset3 = 5.0e-7
        !
        if (itmaxs>0) then
            !----- over-relaxation parameter for each phase
            rlx1 = 1.00      !          DMAX > DSET1
            rlx2 = 1.10      !  DSET1 > DMAX > DSET2
            rlx3 = 1.40      !  DSET2 > DMAX > DSET3
            !CC    STOP              !  DSET3 > DMAX
        else
            rlx1 = 1.0
            rlx2 = 1.0
            rlx3 = 1.0
        endif
        !
        !---- spline bottom & top wall shapes
        jo = 1
        do io = 1, ii
            xb(io) = x(io, jo)
            yb(io) = y(io, jo)
        enddo
        !
        dx = abs(xb(ii) - xb(1))
        dy = abs(yb(ii) - yb(1))
        dxy = max(xb(1), yb(1), dx, dy)
        !
        call scalc(xb, yb, sb, ii)
        call spline(xb, xbs, sb, ii)
        call spline(yb, ybs, sb, ii)
        do io = 1, ii
            sbi(io) = sb(io)
        enddo
        !
        jo = jj
        do io = 1, ii
            xt(io) = x(io, jo)
            yt(io) = y(io, jo)
        enddo
        call scalc(xt, yt, st, ii)
        call spline(xt, xts, st, ii)
        call spline(yt, yts, st, ii)
        do io = 1, ii
            sti(io) = st(io)
        enddo
        !
        rlx = rlx1
        !
        do ipass = 1, itmax
            !
            dmax = 0.
            !
            jo = 1
            if (kbcbot==1 .and. jj>2) then
                !-------- relax to get dX/dj = 0 on lower streamline
                do io = 2, ii - 1
                    ip = io + 1
                    im = io - 1
                    !
                    jp = jo + 1
                    jq = jo + 2
                    !
                    xoo = x(io, jo)
                    xop = x(io, jp)
                    xoq = x(io, jq)
                    !
                    yoo = y(io, jo)
                    yop = y(io, jp)
                    yoq = y(io, jq)
                    !
                    xs = deval(sbi(io), xb, xbs, sb, ii)
                    ys = deval(sbi(io), yb, ybs, sb, ii)
                    !
                    detp = ypos(jp) - ypos(jo)
                    detq = ypos(jq) - ypos(jp)
                    !
                    dxde1 = (xop - xoo) / detp
                    dxde2 = (xoq - xop) / detq
                    !
                    dyde1 = (yop - yoo) / detp
                    dyde2 = (yoq - yop) / detq
                    !
                    rez = xs * (dxde1 - detp * (dxde2 - dxde1) / (detp + detq))          &
                            & + ys * (dyde1 - detp * (dyde2 - dyde1) / (detp + detq))
                    !
                    z_xoo = xs * (-1.0 / detp - 1.0 / (detp + detq))
                    z_yoo = ys * (-1.0 / detp - 1.0 / (detp + detq))
                    !
                    z_s = z_xoo * xs + z_yoo * ys
                    !
                    ds = -rez / z_s
                    !
                    sbi(io) = sbi(io) + ds
                    !
                    x(io, jo) = seval(sbi(io), xb, xbs, sb, ii)
                    y(io, jo) = seval(sbi(io), yb, ybs, sb, ii)
                enddo
            endif
            !
            jo = jj
            if (kbctop==1 .and. jj>2) then
                !-------- relax to get dX/dj = 0 on upper streamline
                do io = 2, ii - 1
                    ip = io + 1
                    im = io - 1
                    !
                    jm = jo - 1
                    jl = jo - 2
                    !
                    xoo = x(io, jo)
                    xom = x(io, jm)
                    xol = x(io, jl)
                    !
                    yoo = y(io, jo)
                    yom = y(io, jm)
                    yol = y(io, jl)
                    !
                    xs = deval(sti(io), xt, xts, st, ii)
                    ys = deval(sti(io), yt, yts, st, ii)
                    !
                    detm = ypos(jo) - ypos(jm)
                    detl = ypos(jm) - ypos(jl)
                    !
                    dxde1 = (xoo - xom) / detm
                    dxde2 = (xom - xol) / detl
                    !
                    dyde1 = (yoo - yom) / detm
                    dyde2 = (yom - yol) / detl
                    !
                    rez = xs * (dxde1 + detm * (dxde1 - dxde2) / (detm + detl))          &
                            & + ys * (dyde1 + detm * (dyde1 - dyde2) / (detm + detl))
                    !
                    z_xoo = xs * (1.0 / detm + 1.0 / (detm + detl))
                    z_yoo = ys * (1.0 / detm + 1.0 / (detm + detl))
                    !
                    z_s = z_xoo * xs + z_yoo * ys
                    !
                    ds = -rez / z_s
                    !
                    sti(io) = sti(io) + ds
                    !
                    x(io, jo) = seval(sti(io), xt, xts, st, ii)
                    y(io, jo) = seval(sti(io), yt, yts, st, ii)
                enddo
            endif
            !
            !
            !------ go over all interior streamlines
            do jo = 2, jj - 1
                jm = jo - 1
                jp = jo + 1
                !
                if (kbcinl==1) then
                    !---------- Neumann BC is specified on i=1 plane:  relax node positions
                    !
                    io = 1
                    ip = io + 1
                    iq = io + 2
                    !
                    xoo = x(io, jo)
                    xpo = x(ip, jo)
                    xqo = x(iq, jo)
                    !
                    yoo = y(io, jo)
                    ypo = y(ip, jo)
                    yqo = y(iq, jo)
                    !
                    xs = x(io, jp) - x(io, jm)
                    ys = y(io, jp) - y(io, jm)
                    !
                    xsn = xs / sqrt(xs**2 + ys**2)
                    ysn = ys / sqrt(xs**2 + ys**2)
                    !
                    dxip = xpos(ip) - xpos(io)
                    dxiq = xpos(iq) - xpos(ip)
                    !
                    dxdx1 = (xpo - xoo) / dxip
                    dxdx2 = (xqo - xpo) / dxiq
                    !
                    dydx1 = (ypo - yoo) / dxip
                    dydx2 = (yqo - ypo) / dxiq
                    !
                    !---------- 2nd-order 3-point difference for tangential velocity
                    rez = xsn * (dxdx1 - dxip * (dxdx2 - dxdx1) / (dxip + dxiq))         &
                            & + ysn * (dydx1 - dxip * (dydx2 - dydx1) / (dxip + dxiq))
                    !
                    z_xoo = xsn * (-1.0 / dxip - 1.0 / (dxip + dxiq))
                    z_yoo = ysn * (-1.0 / dxip - 1.0 / (dxip + dxiq))
                    !
                    z_s = z_xoo * xsn + z_yoo * ysn
                    !
                    srlx = 1.00 * rlx
                    ds = -srlx * rez / z_s
                    !
                    x(io, jo) = x(io, jo) + ds * xsn
                    y(io, jo) = y(io, jo) + ds * ysn
                endif
                !
                if (kbcout==1) then
                    !---------- Neumann BC is specified on i=II plane:  relax node positions
                    !
                    io = ii
                    im = io - 1
                    il = io - 2
                    !
                    xoo = x(io, jo)
                    xmo = x(im, jo)
                    xlo = x(il, jo)
                    !
                    yoo = y(io, jo)
                    ymo = y(im, jo)
                    ylo = y(il, jo)
                    !
                    xs = x(io, jp) - x(io, jm)
                    ys = y(io, jp) - y(io, jm)
                    !
                    xsn = xs / sqrt(xs**2 + ys**2)
                    ysn = ys / sqrt(xs**2 + ys**2)
                    !
                    dxim = xpos(io) - xpos(im)
                    dxil = xpos(im) - xpos(il)
                    !
                    dxdx1 = (xoo - xmo) / dxim
                    dxdx2 = (xmo - xlo) / dxil
                    !
                    dydx1 = (yoo - ymo) / dxim
                    dydx2 = (ymo - ylo) / dxil
                    !
                    !---------- 2nd-order 3-point difference for tangential velocity
                    rez = xsn * (dxdx1 + dxim * (dxdx1 - dxdx2) / (dxim + dxil))         &
                            & + ysn * (dydx1 + dxim * (dydx1 - dydx2) / (dxim + dxil))
                    !
                    z_xoo = xsn * (1.0 / dxim + 1.0 / (dxim + dxil))
                    z_yoo = ysn * (1.0 / dxim + 1.0 / (dxim + dxil))
                    !
                    z_s = z_xoo * xsn + z_yoo * ysn
                    !
                    srlx = 1.00 * rlx
                    ds = -srlx * rez / z_s
                    !
                    x(io, jo) = x(io, jo) + ds * xsn
                    y(io, jo) = y(io, jo) + ds * ysn
                endif

                !-------- relax all points on this streamline by SLOR
                do io = 2, ii - 1
                    im = io - 1
                    ip = io + 1
                    !
                    xmm = x(im, jm)
                    xom = x(io, jm)
                    xpm = x(ip, jm)
                    xmo = x(im, jo)
                    xoo = x(io, jo)
                    xpo = x(ip, jo)
                    xmp = x(im, jp)
                    xop = x(io, jp)
                    xpp = x(ip, jp)
                    ymm = y(im, jm)
                    yom = y(io, jm)
                    ypm = y(ip, jm)
                    ymo = y(im, jo)
                    yoo = y(io, jo)
                    ypo = y(ip, jo)
                    ymp = y(im, jp)
                    yop = y(io, jp)
                    ypp = y(ip, jp)
                    !
                    yap = 0.5 * (yoo + yop)
                    yam = 0.5 * (yoo + yom)
                    yaa = 0.5 * (yap + yam)
                    !
                    dxim = xpos(io) - xpos(im)
                    dxip = xpos(ip) - xpos(io)
                    dxiav = 0.5 * (dxim + dxip)
                    !
                    detm = ypos(jo) - ypos(jm)
                    detp = ypos(jp) - ypos(jo)
                    detav = 0.5 * (detm + detp)
                    !
                    !----------------
                    !            DXDETP = (XOP - XOO)/DETP * YAP/YAA
                    !            DYDETP = (YOP - YOO)/DETP * YAP/YAA
                    !            ALFP = DXDETP**2 + DYDETP**2
                    !C
                    !            DXDETM = (XOO - XOM)/DETM * YAM/YAA
                    !            DYDETM = (YOO - YOM)/DETM * YAM/YAA
                    !            ALFM = DXDETM**2 + DYDETM**2
                    !C
                    !            DXDET = 0.5*(DXDETP + DXDETM)
                    !            DYDET = 0.5*(DYDETP + DYDETM)
                    !            ALF = 0.5*(ALFP + ALFM)
                    !C
                    !C
                    !            DXDXIP = (XPO - XOO)/DXIP
                    !            DYDXIP = (YPO - YOO)/DXIP
                    !            GAMP = DXDXIP**2 + DYDXIP**2
                    !C
                    !            DXDXIM = (XOO - XMO)/DXIM
                    !            DYDXIM = (YOO - YMO)/DXIM
                    !            GAMM = DXDXIM**2 + DYDXIM**2
                    !C
                    !            DXDXI = 0.5*(DXDXIP + DXDXIM)
                    !            DYDXI = 0.5*(DYDXIP + DYDXIM)
                    !            GAM = 0.5*(GAMP + GAMM)
                    !C
                    !C
                    !            BET = DXDET*DXDXI + DYDET*DYDXI
                    !
                    !----------------
                    !            DXDET = 0.5*( (XOP - XOO)/DETP * YAP/YAA
                    !     &                  + (XOO - XOM)/DETM * YAM/YAA )
                    !            DYDET = 0.5*( (YOP - YOO)/DETP * YAP/YAA
                    !     &                  + (YOO - YOM)/DETM * YAM/YAA )
                    !            DXDXI = 0.5*( (XPO - XOO)/DXIP
                    !     &                  + (XOO - XMO)/DXIM )
                    !            DYDXI = 0.5*( (YPO - YOO)/DXIP
                    !     &                  + (YOO - YMO)/DXIM )

                    !
                    dxdet = 0.5 * (xop - xom) / detav
                    dydet = 0.5 * (yop - yom) / detav
                    dxdxi = 0.5 * (xpo - xmo) / dxiav
                    dydxi = 0.5 * (ypo - ymo) / dxiav
                    !
                    !            DXDET2 = 0.5*( (XOP - XOO) * YAP/YAA
                    !     &                   + (XOO - XOM) * YAM/YAA ) / DETAV
                    !            DYDET2 = 0.5*( (YOP - YOO) * YAP/YAA
                    !     &                   + (YOO - YOM) * YAM/YAA ) / DETAV
                    !            ALF = DXDET2**2 + DYDET2**2
                    !            GAM = DXDXI**2 + DYDXI**2
                    !            BET = DXDET2*DXDXI + DYDET2*DYDXI
                    !            AJA = DYDET2*DXDXI - DXDET2*DYDXI
                    !
                    alf = dxdet**2 + dydet**2
                    gam = dxdxi**2 + dydxi**2
                    bet = dxdet * dxdxi + dydet * dydxi
                    aja = dydet * dxdxi - dxdet * dydxi


                    !----------------
                    !
                    !            CXIM = 1.0/(DXIM*DXIAV)
                    !            CXIP = 1.0/(DXIP*DXIAV)
                    !            CETM = 1.0/(DETM*DETAV) * YAM/YAA
                    !            CETP = 1.0/(DETP*DETAV) * YAP/YAA
                    !
                    cxim = detm * detp / (dxim * dxiav)
                    cxip = detm * detp / (dxip * dxiav)
                    cetm = detp / detav * yam / yaa
                    cetp = detm / detav * yap / yaa
                    !
                    b = -alf * cxim
                    a = alf * (cxim + cxip) + gam * (cetm + cetp)
                    c(io) = -alf * cxip
                    if (io==2) b = 0.0
                    !
                    xdse = detm * detp * (xpp - xmp - xpm + xmm) / (4.0 * dxiav * detav)
                    ydse = detm * detp * (ypp - ymp - ypm + ymm) / (4.0 * dxiav * detav)
                    !
                    d(1, io) = alf * ((xpo - xoo) * cxip - (xoo - xmo) * cxim)            &
                            & - 2.0 * bet * xdse + &
                            & gam * ((xop - xoo) * cetp - (xoo - xom) * cetm)            &
                            & - bet * dydxi * dxdet * detm * detp / yaa
                    ! !!! (XPP-XMP-XPM+XMM) / (4.0*DXIAV*DETAV)&
                    !
                    d(2, io) = alf * ((ypo - yoo) * cxip - (yoo - ymo) * cxim)            &
                            & - 2.0 * bet * ydse + &
                            & gam * ((yop - yoo) * cetp - (yoo - yom) * cetm)            &
                            & - bet * dydxi * dydet * detm * detp / yaa
                    ! !!! (YPP-YMP-YPM+YMM) / (4.0*DXIAV*DETAV)&
                    !
                    ainv = 1.0 / (a - b * c(im))
                    c(io) = c(io) * ainv
                    d(1, io) = (d(1, io) - b * d(1, im)) * ainv
                    d(2, io) = (d(2, io) - b * d(2, im)) * ainv
                    !
                enddo
                !
                d(1, ii) = 0.
                d(2, ii) = 0.
                !
                ifin = ii - 1
                do iback = 2, ifin
                    io = ii - iback + 1
                    ip = io + 1
                    d(1, io) = d(1, io) - c(io) * d(1, ip)
                    d(2, io) = d(2, io) - c(io) * d(2, ip)
                    !
                    x(io, jo) = x(io, jo) + rlx * d(1, io)
                    y(io, jo) = y(io, jo) + rlx * d(2, io)
                    ad1 = abs(d(1, io))
                    ad2 = abs(d(2, io))
                    dmax = max(dmax, ad1, ad2)
                enddo
                !
            enddo
            !
            !        IF(MOD(IPASS,10).EQ.0) THEN
            !          WRITE(*,*) IPASS, '  Dmax = ', DMAX, RLX
            !        ENDIF
            !
            if (dmax<toler * dxy) return
            !
            rlx = rlx1
            if (dmax<dset1 * dxy) rlx = rlx2
            if (dmax<dset2 * dxy) rlx = rlx3
            if (dmax<dset3 * dxy) return
            !
        enddo
        !
    end subroutine axell
end module m_inigrd
