module m_dfdcsubs
    implicit none
contains
    !*==DFINIT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
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
    !     Version 070-ES1
    !     Philip Carter, Esotec Developments, February 2009
    !     philip (at) esotec (dot) org
    !
    !     Changes from 0.70:
    !
    !     DFINIT
    !       SVERSION (line 66)
    !       LCOORD   (line 424)
    !       LLOAD    (line 436)
    !       IDEVPS = 4   ! Color PostScript
    !
    !     DFLOAD
    !       File version displayed in exponential format on load.
    !       Reporting formats cleaned up.
    !
    !     DFSAVE
    !       File version written in exponential format (allows ES versions).
    !       ANAME (geometry name) correctly saved.
    !       RPM2 written.
    !
    !     Version 070-ES2
    !     ROTINIT   LLOFT condition imposed on tip gap paneling (for lofting)
    !
    !=========================================================================
    !
    subroutine dfinit(ldebug)
        use i_dfdc
        use m_atmo, only : atmo
        use m_aero, only : putaero
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: ldebug
        !
        ! Local variables
        !
        real :: a0, cdmin, cldmin, clmax, clmin, cmcon, dcdcl2, &
                & dcdcl2s, dclda, dclda_stall, dcl_stall, mcrit, &
                & reref, rexp, rpm, toc, xisect
        character(10), save :: digits
        integer :: i, ib, ic, iel, ip, ipact, ir, irpn, k1, &
                & k2, kp, l, n, npol
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------
        !     Sets all initial defaults
        !--------------------------------------
        data digits/'0123456789'/
        !
        !---- Set version (in case this code is used by DLL w/o main DFDC routine)
        !
        !      VERSION  =  0.70e02          ! This gets saved to file
        !      SVERSION = '070-ES2'         ! A more flexible format
        !
        !---------------------------------------------------------------
        !---- Debug output
        !cc      LDBG = .FALSE.
        ldbg = ldebug
        !---- Logical unit for debugging printout
        lundbg = 40
        !---------------------------------------------------------------
        !
        !---- set digit string array for constructing filename
        do kp = 1, 99
            k2 = kp / 10 + 1
            k1 = kp - 10 * (kp / 10) + 1
            pnum(kp) = digits(k2:k2) // digits(k1:k1)
        enddo
        !
        !---------------------------------------------------------------
        !
        do ip = 1, ipx
            gam(ip) = 0.
            sig(ip) = 0.
            gth(ip) = 0.
            gamvsp(ip) = 0.
            sigvsp(ip) = 0.
            vmavg(ip) = 0.
            izero(ip) = 0
            ip2ig(ip) = 0
        enddo
        do ic = 1, icx
            ic2ig(ic) = 0
        enddo
        do iel = 1, nex
            gamset(iel) = 0.
            sigset(iel) = 0.
        enddo
        do ib = 1, ibx
            one(ib) = 1.0
            zer(ib) = 0.
        enddo
        !
        !---- no panel nodes, control points
        nbtot = 0
        nptot = 0
        nctot = 0
        do iel = 1, nex
            ibfrst(iel) = 0
            iblast(iel) = -1
            ipfrst(iel) = 0
            iplast(iel) = -1
            icfrst(iel) = 0
            iclast(iel) = -1
            nbtype(iel) = -999
            !
            iel2ir(iel) = 0
        enddo
        !
        !---------------------------------------------------------------
        !---- no flow yet
        qinf = 0.0
        !---- start off with incompressible flow
        mach = 0.
        mach1 = 0.
        !---- acceleration due to gravity for scaling centrifugal blade tension (m/s^2)
        gee = 9.81
        !---- default unity reference velocity
        qref = 1.0
        !---- setup for SL conditions, Std temperature
        alth = 0.0
        deltat = 0.0
        !! sea level atmosphere parameters
        !     RHO =  1.226      ! fluid density         kg/m**3
        !     RMU =  1.78E-05   ! dynamic viscosity     kg/m-s
        !     VSO =  340.0      ! speed of sound        m/s
        call atmo(alth, deltat, vso, rho, rmu)
        !---------------------------------------------------------------
        !
        !---- no drag objects defined
        ndobj = 0
        ieldrgobj(1) = 0
        ldrgobj = .false.
        !
        !---------------------------------------------------------------
        !---- no rotor yet
        nrotor = 0
        !---- no blade defined yet
        lbldef = .false.
        !
        do n = 1, nrx
            !---- disk type set to undefined
            irtype(n) = 0
            !---- no defining points for rotor or actuator disk
            nrdef(n) = 0
            irtypdef(n) = 0
            !---- no elements defined for rotor source line
            ielrotor(n) = 0
            !
            atip(n) = 0.0
            adisk(n) = 0.0
            xdisk(n) = 999.  ! default rotor disk location
            nrbld(n) = 5
            omega(n) = 0.0
        enddo
        !---- default blade is rotor with RPM=1000
        rpm = 1000.0
        omega(1) = pi * rpm / 30.0
        !
        tgap = 0.0  ! no tip gap on rotor
        tgapzl = 0.0  ! tip gap of zero loss (default is zero)
        !
        nrp = 0
        nrc = 0
        !
        nrsta = 11   ! default # of radial stations
        xdwklen = 1.5
        ! default wake length in duct diameters
        xwake = -999.
        nwake = 0
        !
        tduct = 0.0
        ttot = 0.0
        tvis = 0.0
        qtot = 0.0
        qvis = 0.0
        ptot = 0.0
        pvis = 0.0
        !
        !
        do ir = 1, irx
            chdes(ir) = 1.0
            cldes(ir) = 0.8
            clpos(ir) = cldes(ir)
            clneg(ir) = -cldes(ir)
            ir2iel(ir) = 0
            !
            do n = 1, nrx
                bgam(ir, n) = 0.0
                do l = 1, 3
                    vind(l, ir, n) = 0.0
                    vabs(l, ir, n) = 0.0
                    vrel(l, ir, n) = 0.0
                enddo
                chr(ir, n) = chdes(ir)
                clr(ir, n) = cldes(ir)
            enddo
        enddo
        !
        !---------------------------------------------------------------
        !--- Setup aerodynamic data for blade elements
        !--- Default aero properties
        !
        a0 = 0.       ! zero lift angle of attack   radians
        dclda = 6.0     ! lift curve slope            /radian
        clmax = 1.5    ! stall Cl
        clmin = -0.5    ! negative stall Cl
        dcl_stall = 0.1
        ! CL increment from incipient to total stall
        dclda_stall = 0.1
        ! stalled lift curve slope    /radian
        cmcon = -0.1   ! section Cm  (for pitch-axis moments)
        cdmin = 0.013 ! minimum Cd
        cldmin = 0.5   ! Cl at minimum Cd
        dcdcl2 = 0.03  ! d(Cd)/d(Cl**2)
        dcdcl2s = 0.0   ! d(Cd)/d(Cl**2) secondary (annulus) drag (dead)
        reref = 200000. ! Reynolds Number at which Cd values apply
        rexp = -0.4    ! Exponent for Re scaling of Cd:  Cd ~ Re**exponent
        mcrit = 0.8     ! Critical Mach number
        toc = 0.1     ! section thickness ratio (dead)
        !
        xpaxis = 0.35   ! x/c location of pitch axis for blade plots, lofting
        !
        !--- Install default data into aero section #1 for each disk
        do n = 1, nrx
            naero(n) = 1
            xisect = 0.0
            call putaero(n, naero(n), xisect, a0, clmax, clmin, dclda, &
                    & dclda_stall, dcl_stall, cdmin, cldmin, dcdcl2, dcdcl2s, &
                    & cmcon, mcrit, toc, reref, rexp)
            !--- set pointers from blade sections to aero sections
            do i = 1, irx
                iaero(i, n) = 1
            enddo
        enddo
        !
        !---- no BL solution yet
        lvisc = .false.
        ampmax = 9.0
        uren = 0.0
        do iel = 1, nex
            cxvis(iel) = 0.0
            sfix(1, iel) = 1.0
            sfix(2, iel) = 1.0
            strn(1, iel) = 1.0
            strn(2, iel) = 1.0
            ssep(1, iel) = 1.0
            ssep(2, iel) = 1.0
            icblbeg(1, iel) = 0
            icblbeg(2, iel) = 0
            icblend(1, iel) = 0
            icblend(2, iel) = 0
        enddo
        !---- no inflow velocities defined
        ninfl = 0
        !
        !---------------------------------------------------------------
        !---- no grid dimensions yet
        ii = 0
        jj = 0
        !
        !---- by default input geometry is respaced
        lrspc = .true.
        !---------------------------------------------------------------
        !---- default number of panels per element, panel-number exponent
        npandef = 61
        fpandef = 0.7
        !---- set panel parameters, no default spacing or user-defined spacing
        lrspcusr = .false.
        lrspcdef = .false.
        do iel = 1, nex
            npan(iel) = -999
            cvex(iel) = 1.0
            smof(iel) = 1.0
            fsle(iel) = 0.6
            fste(iel) = 0.6
            !------ refinement-station parameters
            do irpn = 1, nrpnx
                crrat(irpn, iel) = 1.0
                srpn1(irpn, iel) = 0.0
                srpn2(irpn, iel) = 0.0
            enddo
        enddo
        !
        !---- assume no case-paneling parameter file exists
        pfile = ' '
        !
        !---------------------------------------------------------------
        !---- maximum TE gap (relative to perimeter) to deem airfoil as being closed
        dtebod = 0.1
        !---- minimum TE gap (relative to perimeter) to give closed airfoil a TE panel
        dtepan = 0.00001
        !
        !---------------------------------------------------------------
        !---- no geometry yet
        nel = 0
        nbel = 0
        name = ' '
        !
        !---- no accumulated point sequence yet
        npol = 0
        lpacc = .false.
        !
        ipact = 0
        !
        npseq = 0
        !
        !---------------------------------------------------------------
        !---- Convergence flag for iterative solutions, no converged solution yet
        lconv = .false.
        !---- Iterative solver parameters
        rlxsolv = 0.4
        epssolv = 0.0002
        itrmaxsolv = 50
        !---- VMAVG is not initialized
        lvmav = .false.
        vavginit = qinf + 1.0
        !
        !---- vortex wakes will not be automatically realigned
        lwrlx = .false.
        !
        !---------------------------------------------------------------
        !---- don't normalize airfoil upon input
        lnorm = .false.
        !
        !---- no geometry modes defined yet
        nmod = 0
        lxymov = .false.
        !
        !---------------------------------------------------------------
        !---- no valid control-point pointers yet
        lncvp = .false.
        !
        !---- no valid Aero,Geometry influence matrices yet
        lqaic = .false.
        lqgic = .false.
        lvaic = .false.
        lvgic = .false.
        !
        !---- no valid control-point velocities yet
        lqcnt = .false.
        !
        !---- no valid system pointers yet
        lsysp = .false.
        !
        !---- no factored system matrix for gamma yet
        lgsys = .false.
        lqsys = .false.
        !
        !---- no gamma solution yet, no unit solutions
        lgama = .false.
        lgamu = .false.
        !
        !---- no source-influence pointers or matrix
        lsigp = .false.
        lsigm = .false.
        !
        !---- don't plot reference cp or force data
        lpref = .false.
        lfref = .false.
        !
        !---- default filename prefix
        prefix = ' '
        !
        cpfile = ' '
        frfile = ' '
        !---------------------------------------------------------------
        !---- no target elements yet for QDES or GDES
        ielqdes = 0
        ielgdes = 1
        !
        !---- enforce geometry slope-matching at grafted-geometry endpoints
        lgslop = .true.
        !
        !---- no geometric symmetry assumed
        lgsymm = .false.
        !
        !---- same geometry between panel and buffer
        lgsame = .false.
        !
        !---- plot home positions
        lhompl = .true.
        !
        !---- no specified speed distributions yet
        nqsp = 0
        nseg = 0
        lqspec = .false.
        lqsppl = .false.
        lgsppl = .false.
        !
        !---- enforce slope-matching and gamma curvature at inverse-segment endpoints
        lqslop = .true.
        lqcurv = .false.
        !
        !---- inverse-segment endpoints on edge of a surface airfoil are free to move
        lgnfix = .false.
        !
        !---- do not display viscous Q(s)
        lqsvis = .false.
        !
        !---- do not reverse s or Q direction in Q(s) plots
        lqsrev = .false.
        lqqrev = .false.
        !
        !---- maximum underrelaxation factor
        rlxmax = 1.0
        !
        !---------------------------------------------------------------
        !---- show point-strength and wake-strength values
        lpsho = .true.
        lwsho = .true.
        !
        !---- plot point-singularities, panel vertices, element numbers
        lpplt = .true.
        lvplt = .true.
        leplt = .true.
        !
        !---- do not plot symbols on Cp or Q line plots
        lsplt = .false.
        !
        !---- plot displacement surfaces, don't use streamline integration to generate
        ldint = .false.
        !
        !---- plot grid on Cp or Q plots
        lpgrid = .true.
        !
        !---- plot stagnation point on Cp or Q plots
        lpstag = .true.
        !
        !---- plot grid on aero, airfoil geometry plots, ticks, parameters
        lagrid = .true.
        lggrid = .true.
        lgtick = .true.
        lgparm = .true.
        !
        lgeoplx = .false.
        ! for xfoil plot routines
        !
        !---- do not print neg. rpm Beta and Alfa output in local coordinates
        lcoord = .false.
        !
        !-----a case file has not been loaded
        lload = .false.
        !
        !
    end subroutine dfinit
    !*==GENGEOM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! DFINIT


    subroutine gengeom
        use i_dfdc
        use m_geom, only : cvpgen
        use m_pnsubs, only : pancop, pangen, pandef
        use m_inigrd, only : inigrd
        use m_adjpanl, only : adjpanl
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        logical :: lquery
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !     Generates case geometry from buffer airfoils
        !     and does the necessary geometry processing.
        !----------------------------------------------------
        !
        !---- process current buffer geometry to identify element types and
        !     set element parameters
        !      write(80,*) 'before geproc'
        !      CALL XYBDMP
        !      CALL GEPROC
        !      IF(LDBG) THEN
        !        CALL GEPRINT(LUNDBG)
        !      ENDIF
        !
        !========================================================================
        !---- Check for respacing flag to redistribute points on buffer geometry
        if (lrspc) then
            !
            !---- If paneling parameters are unset then reset default paneling
            !     parameters for buffer geometry elements
            !      write(*,*) 'gengeom  npan',NPAN(1)
            !          WRITE(*,*) 'Setting default respacing parameters'
            if (.not.lrspcdef) call pandef(0)
            !
            !        WRITE(*,*) 'Repaneling buffer geometry'
            lquery = .false.
            call pangen(lquery)
            !
        else
            !        WRITE(*,*) 'Copying buffer geometry into panel arrays'
            call pancop
        endif
        !
        !----------------------------------------------------------------
        !---- At this point the input airfoils are in the panel geometry
        lgsame = .true.
        !
        !---- Adjust paneling on airfoils for rotor and wake grid
        call adjpanl
        !      write(80,*) 'after ADJPANL PDMP'
        !      CALL XYPDMP(80)
        !
        !---- Reset # of elements to buffer geometry, remaining elements are added
        !     to this: drag elements, rotor line sources, then vortex wakes
        nel = nbel
        !
        !---- Set up drag area element(s)
        call setdrgobj
        !
        !---- Set up rotor points
        call rotorinit
        !
        !---- Set up rotor wake grid from repaneled airfoils
        call inigrd
        !
        !---- Set up rotor wake elements from wake grid
        call setrotwak
        !
        !---- Sets up control points and pointers to panels
        call cvpgen
        !
        !---- Initialize rotor/wake pointers
        call rotpinit
        !
        !---- List all defined elements
        if (ldbg) then
            call pelist(6)
            call pelist(lundbg)
        endif
        call wakebox
        !
        !----------------------------------------------------------------
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
        nqsp = 0
        nseg = 0
        !
        write (*, *)
        write (*, *) 'Geometry loaded'
        !
    end subroutine gengeom
    !*==PELIST.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GENGEOM


    subroutine pelist(lu)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: lu
        !
        ! Local variables
        !
        integer :: i, iel, ip1, ip2, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Lists elements defined in panel geometry
        !------------------------------------------------------
        !
        !---- print out element info
        do iel = 1, nel
            write (lu, 1001)
            write (lu, 1005) iel
            !
            ip1 = ipfrst(iel)
            ip2 = iplast(iel)
            i = ip1
            n = ip2 - ip1 + 1
            !
            write (lu, 1010) n
            !
            if (netype(iel)==0) then
                if (.not.(lbody(iel))) then
                    write (lu, 2000) 'Open surface'
                elseif (lxbod(iel)) then
                    write (lu, 2000) 'Closed body on axis'
                else
                    write (lu, 2000) 'Closed body'
                endif
                !
            elseif (netype(iel)==1) then
                write (lu, 2000) 'Prescribed wake surface'
                !
            elseif (netype(iel)==2) then
                write (lu, 2000) 'Prescribed line singularity on axis'
                !
            elseif (netype(iel)==3) then
                write (lu, 2000) 'Prescribed ring singularity'
                !
            elseif (netype(iel)==4) then
                write (lu, 2000) 'Prescribed point singularity'
                !
            elseif (netype(iel)==5) then
                write (lu, 2000) 'Source-only line singularity'
                !
            elseif (netype(iel)==6) then
                write (lu, 2000) 'Rotor source-only line singularity'
                !
            elseif (netype(iel)==7) then
                write (lu, 2000) 'Rotor vortex wake line'
                !
            endif
            !
            if (ltpan(iel)) write (lu, 2000) 'TE-gap panel will be used'
            !
            if (iewake(iel)/=0) write (lu, 2100)                          &
                    & 'Prescribed-wake element', &
                    & iewake(iel), '  is attached'
            !
            if (netype(iel)==0 .or. netype(iel)==1 .or. netype(iel)==2)  &
                    & then
                write (lu, 1050) xple(iel), yple(iel), xpte(iel), &
                        & ypte(iel)
            else
                write (lu, 1060) 'First point at', xp(ip1), yp(ip1)
                write (lu, 1060) 'Last  point at', xp(ip2), yp(ip2)
            endif
            !
        enddo
        !
        write (lu, 1001)
        !
        return
        !...............................................................
        1001 format (/' -----------------------------------------------')
        1005 format (' Element', i3, ' ...')
        1010 format ('   Number of input coordinate points:', i4)
        1020 format ('     delta(x) =', f10.5, '   Scale =', f10.5, &
                &/'     delta(y) =', f10.5, '   Angle =', f10.5, ' deg')
        1050 format ('    LE x,y  =', 2f10.5/'    TE x,y  =', 2f10.5)
        1060 format (' ', a, ' x,y  =', 2f10.5)
        2000 format ('   ', a)
        2100 format ('   ', a, i3, a)
    end subroutine pelist
    !*==GEPROC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PELIST



    subroutine geproc
        use i_dfdc
        use m_geutil, only : minmax
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: db11sq, db12sq, db1tsq, db21sq, db22sq, db2tsq, &
                & dsbmin, dsbsq, dsbtol, dummy
        integer :: i, ib, ib1, ib2, iel, kb1, kb2, kel, n
        logical :: lyzero
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !     Processes buffer airfoil geometry to set
        !     element type, various limits and pointers.
        !----------------------------------------------------
        !
        !---- process all elements
        do iel = 1, nbel
            call elproc(iel)
            !
            i = ibfrst(iel)
            n = iblast(iel) - ibfrst(iel) + 1
            call minmax(n, xb(i), xbmine(iel), xbmaxe(iel))
            call minmax(n, yb(i), ybmine(iel), ybmaxe(iel))
        enddo
        !
        !---- set overall min,max over all elements
        call minmax(nbel, xbmine, xbmin, dummy)
        call minmax(nbel, ybmine, ybmin, dummy)
        call minmax(nbel, xbmaxe, dummy, xbmax)
        call minmax(nbel, ybmaxe, dummy, ybmax)
        !
        do iel = 1, nbel
            !------ default reference point for each element
            xbrefe(iel) = 0.
            ybrefe(iel) = 0.
            !------ clear element movement accumulators
            dxbsum(iel) = 0.
            dybsum(iel) = 0.
            agbsum(iel) = 0.
            xfbsum(iel) = 1.0
            yfbsum(iel) = 1.0
            !------ first assume this element doesn't have a wake
            iewake(iel) = 0
        enddo
        !
        !---- set element-type indicators and wake pointers
        do iel = 1, nbel
            ib1 = ibfrst(iel)
            ib2 = iblast(iel)
            !
            !------ is this element entirely on the y=0 line?
            lyzero = .true.
            do ib = ib1, ib2
                if (yb(ib)/=0.0) lyzero = .false.
            enddo
            !
            if (ib1==ib2) then
                !------- only one point defines element...
                if (lyzero) then
                    !-------- axisymm: point singularity on axis
                    netype(iel) = 4
                else
                    !-------- axisymm: ring singularity
                    netype(iel) = 3
                endif
                !
            elseif (lyzero) then
                !------- line element on axis
                netype(iel) = 2
                !
                !------- check for source line element defined in buffer geometry (old scheme)
            elseif (nbtype(iel)==5) then
                netype(iel) = nbtype(iel)
                !
            else
                !------- assume this element is a solid surface
                netype(iel) = 0
                !
                !------- Now check for wake elements
                !------- a wake element is attached to TE of some element before it
                do kel = 1, iel - 1
                    !--------- check only solid boundarys (NETYPE=0) for attached wakes
                    if (netype(kel)==0) then
                        kb1 = ibfrst(kel)
                        kb2 = iblast(kel)
                        db11sq = (xb(ib1) - xb(kb1))**2 + (yb(ib1) - yb(kb1))**2
                        db21sq = (xb(ib2) - xb(kb1))**2 + (yb(ib2) - yb(kb1))**2
                        db12sq = (xb(ib1) - xb(kb2))**2 + (yb(ib1) - yb(kb2))**2
                        db22sq = (xb(ib2) - xb(kb2))**2 + (yb(ib2) - yb(kb2))**2
                        db1tsq = (xb(ib1) - xbte(kel))**2 + (yb(ib1) - ybte(kel)) &
                                & **2
                        db2tsq = (xb(ib2) - xbte(kel))**2 + (yb(ib2) - ybte(kel)) &
                                & **2
                        !
                        dsbsq = min(db11sq, db12sq, db21sq, db22sq, db1tsq, db2tsq)
                        dsbmin = sqrt(dsbsq)
                        !
                        !----------- this element IEL is KEL's wake if it starts at KEL's TE point
                        dsbtol = 0.0001 * abs(sb(kb2) - sb(kb1))
                        if (dsbmin<dsbtol) then
                            netype(iel) = 1
                            iewake(kel) = iel
                            exit
                        endif
                    endif
                enddo
                !
                !
            endif
            !
        enddo
        !
        !---- set wake-termination box outline
        !c      CALL WAKEBOX
        !
    end subroutine geproc
    !*==GEPRINT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GEPROC


    subroutine geprint(lu)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: lu
        !
        ! Local variables
        !
        integer :: i, ib1, ib2, iel, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !     Processes buffer airfoil geometry to set
        !     various limits and pointers.
        !----------------------------------------------------
        !
        !---- print out element info
        do iel = 1, nbel
            write (lu, 1001)
            write (lu, 1005) iel
            !
            ib1 = ibfrst(iel)
            ib2 = iblast(iel)
            i = ib1
            n = ib2 - ib1 + 1
            !
            write (lu, 1010) n
            !
            if (area2da(iel)>=0.0) then
                write (lu, 2000) 'Counterclockwise ordering'
            else
                write (lu, 2000) 'Clockwise ordering'
            endif
            !
            if (netype(iel)==0) then
                if (.not.(lbody(iel))) then
                    write (lu, 2000) 'Open surface'
                elseif (lxbod(iel)) then
                    write (lu, 2000) 'Closed body on axis'
                else
                    write (lu, 2000) 'Closed body'
                endif
                !
            elseif (netype(iel)==1) then
                write (lu, 2000) 'Prescribed wake surface'
                !
            elseif (netype(iel)==2) then
                write (lu, 2000) 'Prescribed line singularity on axis'
                !
            elseif (netype(iel)==3) then
                write (lu, 2000) 'Prescribed ring singularity'
                !
            elseif (netype(iel)==4) then
                write (lu, 2000) 'Prescribed point singularity'
                !
            endif
            !
            if (ltpan(iel)) write (lu, 2000) 'TE-gap panel will be used'
            !
            if (iewake(iel)/=0) write (lu, 2100)                          &
                    & 'Prescribed-wake element', &
                    & iewake(iel), '  is attached'
            !
            if (netype(iel)==0 .or. netype(iel)==1 .or. netype(iel)==2)  &
                    & then
                write (lu, 1050) xble(iel), yble(iel), xbte(iel), &
                        & ybte(iel)
            else
                write (lu, 1060) xb(ib1), yb(ib1)
            endif
            !
        enddo
        !
        write (lu, 1001)
        !
        return
        !
        !...............................................................
        1001 format (/' -----------------------------------------------')
        1005 format (' Element', i3, ' ...')
        1010 format ('   Number of input coordinate points:', i4)
        1020 format ('     delta(x) =', f10.5, '   Scale =', f10.5, &
                &/'     delta(y) =', f10.5, '   Angle =', f10.5, ' deg')
        1050 format ('    LE x,y  =', 2f10.5/'    TE x,y  =', 2f10.5)
        1060 format ('       x,y  =', 2f10.5)
        2000 format ('   ', a)
        2100 format ('   ', a, i3, a)
    end subroutine geprint
    !*==ELPROC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GEPRINT



    subroutine elproc(iel)
        use i_dfdc
        use m_spline, only : segspl, scalc, seval, sinvrt
        use m_geutil, only : lefind, axcalc, aecalc
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
        real :: dsbte, perim, sblen, ss, xtmp, ytmp
        integer :: i, ib, ib1, ib2, ibcut, n
        logical :: lcut, ltmp
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Processes buffer geometry element IEL.
        !------------------------------------------------------
        !
        ib1 = ibfrst(iel)
        ib2 = iblast(iel)
        !
        i = ib1
        n = ib2 - ib1 + 1
        call scalc(xb(i), yb(i), sb(i), n)
        !
        !----- check for body that crosses the axis (defined by upper and lower sides)
        lcut = .false.
        do ib = ib1, ib2 - 1
            if (yb(ib)>0.0 .and. yb(ib + 1)<=0.0) then
                lcut = .true.
                ibcut = ib + 1
                exit
            endif
        enddo
        !
        if (lcut) then
            !----- element crosses axis, cut it at the axis and use only upper half
            if (ldbg) then
                write (*, *) 'Cutting element ', iel, ' at X axis '
                write (lundbg, *) 'Cutting element ', iel, ' at X axis '
            endif
            call segspl(xb(i), xbs(i), sb(i), n)
            call segspl(yb(i), ybs(i), sb(i), n)
            ss = sb(ibcut)
            call sinvrt(ss, 0.0, yb(ib1), ybs(ib1), sb(ib1), n)
            xb(ibcut) = seval(ss, xb(ib1), xbs(ib1), sb(ib1), n)
            yb(ibcut) = 0.0
            iblast(iel) = ibcut
            !
            !----- recalculate indices for element and arc length
            ib1 = ibfrst(iel)
            ib2 = iblast(iel)
            i = ib1
            n = ib2 - ib1 + 1
            call scalc(xb(i), yb(i), sb(i), n)
            !c       DO IB = IB1, IB2
            !c         write(50,*) IB, XB(IB), YB(IB)
            !c       END DO
        endif
        !
        !
        if (n==1) then
            !---- single point element (doublet/source point)
            lv1zr(iel) = .false.
            lv2zr(iel) = .false.
            lxbod(iel) = .false.
            lbody(iel) = .false.
            ltpan(iel) = .false.
            lqnzr(iel) = .false.
            !
        else
            !---- multiple point element (doublet/source line or vortex/source panels)
            sblen = abs(sb(ib2) - sb(ib1))
            !
            !---- are first,last points effectively on the axis?
            lv1zr(iel) = abs(yb(ib1))<sblen * dtepan
            lv2zr(iel) = abs(yb(ib2))<sblen * dtepan
            !
            if (lv1zr(iel)) then
                !---- TE gap is really a TE disk at point IB2
                dsbte = abs(yb(ib2))
            else
                !---- usual TE gap
                dsbte = sqrt((xb(ib1) - xb(ib2))**2 + (yb(ib1) - yb(ib2))**2)
            endif
            !
            !---- axisymmetric body closed on the axis ?
            lxbod(iel) = lv1zr(iel) .and. lv2zr(iel)
            !
            !---- finite-thickness body ?
            lbody(iel) = lv1zr(iel) .or. lv2zr(iel) .or.                   &
                    & dsbte<abs(sb(ib2) - sb(ib1)) * dtebod
            !
            !---- add TE-gap panel ?
            ltpan(iel) = lbody(iel) .and. dsbte>abs(sb(ib2) - sb(ib1)) * dtepan
            !
            !---- explicitly zero out the body normal velocity QNDOF?
            !c     LQNZR(IEL) = .TRUE.
            lqnzr(iel) = .false.
            !
        endif
        !
        !---- 2D buffer geometry data
        !---- set 2D area, centroid, section moduli - for area
        call aecalc(n, xb(i), yb(i), one, 0, perim, area2da(iel), xbcen2da(iel), &
                & ybcen2da(iel), eixx2da(iel), eiyy2da(iel), eixy2da(iel))
        !---- set 2D area, centroid, section moduli - per unit skin thickness
        call aecalc(n, xb(i), yb(i), one, 2, perim, area2dt(iel), xbcen2dt(iel), &
                & ybcen2dt(iel), eixx2dt(iel), eiyy2dt(iel), eixy2dt(iel))
        !
        !---- Axisymmetric buffer geometry data
        !---- set axisymmetric surface area, volume, centroid, radii of gyration
        call axcalc(n, xb(i), yb(i), one, 0, volumv(iel), asurfv(iel), &
                & xbcenv(iel), ybcenv(iel), rgyrxv(iel), rgyryv(iel))
        !
        !---- assume element will be used and written in input order by default
        lrevel(iel) = .false.
        !
        !----- area is negative (clockwise order)... mark for reversing the order
        if (lbody(iel) .and. area2da(iel)<0.0) lrevel(iel) = .true.
        !
        if (lrevel(iel)) then
            !---- reverse the order of the points
            do ib = ib1, (ib2 + ib1) / 2
                xtmp = xb(ib2 - ib + ib1)
                ytmp = yb(ib2 - ib + ib1)
                xb(ib2 - ib + ib1) = xb(ib)
                yb(ib2 - ib + ib1) = yb(ib)
                xb(ib) = xtmp
                yb(ib) = ytmp
            enddo
            ltmp = lv1zr(iel)
            lv1zr(iel) = lv2zr(iel)
            lv2zr(iel) = ltmp
        endif
        !
        call scalc(xb(i), yb(i), sb(i), n)
        call segspl(xb(i), xbs(i), sb(i), n)
        call segspl(yb(i), ybs(i), sb(i), n)
        !
        if (lbody(iel) .and. .not.lxbod(iel)) then
            !----- body off axis
            if (lv2zr(iel)) then
                sble(iel) = sb(ib2)
                xble(iel) = xb(ib2)
                yble(iel) = yb(ib2)
                xbte(iel) = xb(ib1)
                ybte(iel) = yb(ib1)
            else
                call lefind(sble(iel), xb(i), xbs(i), yb(i), ybs(i), sb(i), n)
                xble(iel) = seval(sble(iel), xb(i), xbs(i), sb(i), n)
                yble(iel) = seval(sble(iel), yb(i), ybs(i), sb(i), n)
                xbte(iel) = 0.5 * (xb(ib1) + xb(ib2))
                ybte(iel) = 0.5 * (yb(ib1) + yb(ib2))
            endif
            !
            !----- other body... set LE,TE from leftmost and rightmost points
        elseif (xb(ib1)<xb(ib2)) then
            sble(iel) = sb(ib1)
            xble(iel) = xb(ib1)
            yble(iel) = yb(ib1)
            xbte(iel) = xb(ib2)
            ybte(iel) = yb(ib2)
        else
            sble(iel) = sb(ib2)
            xble(iel) = xb(ib2)
            yble(iel) = yb(ib2)
            xbte(iel) = xb(ib1)
            ybte(iel) = yb(ib1)
            !
        endif
        !
        !---- set Cp-side flag
        if (lbody(iel)) then
            !----- outer surface only (to right side of +s CCW direction)
            isplot(iel) = +1
        else
            !----- default... plot both sides
            isplot(iel) = 0
        endif
        !
    end subroutine elproc
    !*==WAKEBOX.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ELPROC



    subroutine wakebox
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: dxbox, dybox
        integer :: iel
        logical :: lwbset
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------
        !     Sets limits of "wake box".
        !     A wake is terminated when it leaves this box.
        !---------------------------------------------------
        !
        xwbox(1) = 1.0e23
        xwbox(2) = -1.0e23
        ywbox(1) = 1.0e23
        ywbox(2) = -1.0e23
        lwbset = .false.
        do iel = 1, nel
            if (netype(iel)==0 .or. netype(iel)==7) then
                !------- box defined only by solid or vortex wake elements
                xwbox(1) = min(xwbox(1), xpmine(iel))
                xwbox(2) = max(xwbox(2), xpmaxe(iel))
                ywbox(1) = min(ywbox(1), ypmine(iel))
                ywbox(2) = max(ywbox(2), ypmaxe(iel))
                lwbset = .true.
            endif
        enddo
        if (.not.lwbset) then
            !----- use default if no solid elements exist
            xwbox(1) = xbmin
            xwbox(2) = xbmax
            ywbox(1) = ybmin
            ywbox(2) = ybmax
        endif
        !
        !c      DBOX = MAX( XWBOX(2)-XWBOX(1) , YWBOX(2)-YWBOX(1) )
        !c      XWBOX(1) = XWBOX(1) - DBOX
        !c      XWBOX(2) = XWBOX(2) + DBOX
        !c      YWBOX(1) = YWBOX(1) - DBOX
        !c      YWBOX(2) = YWBOX(2) + DBOX
        !
        dxbox = xwbox(2) - xwbox(1)
        dybox = ywbox(2) - ywbox(1)
        xwbox(1) = xwbox(1) - 0.5 * dxbox
        xwbox(2) = xwbox(2) - 0.02 * dxbox
        ywbox(1) = 0.0
        ywbox(2) = ywbox(2) + dybox
        !
    end subroutine wakebox
    !*==SETDRGOBJ.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! WAKEBOX


    subroutine setdrgobj
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: nd
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Sets up drag area elements into paneled geometry
        !------------------------------------------------------
        !
        if (ndobj<1) return
        !
        nptot = iplast(nel)
        !
        !---- Use rotor points to define paneled element for rotor source line
        do nd = 1, ndobj
            call addelem5(xddef(1, nd), yddef(1, nd), nddef(nd))
            ieldrgobj(nd) = nel
        enddo
        ldrgobj = .true.
        !
    end subroutine setdrgobj
    !*==ROTORINIT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine rotorinit
        use i_dfdc
        use m_aero, only : setiaero
        use m_spline, only : segspl, seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real, dimension(irx) :: betar0, bgam0, chr0, cldes0, t1, &
                & t1s, t2, t2s, t3, t3s, t4, t4s, &
                & yrc0
        real :: cld1, cld2, cln1, cln2, clp1, clp2, dr, dx, dy, &
                & tg, tgapzlt, xcb, xdw, ycb, ydw
        integer :: i, ic, ipcb, ipdw, ir, n, nrc0, nrcsav
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Sets up rotor geometry from panel geometry
        !------------------------------------------------------
        !
        if (nrotor==0) then
            write (*, *) 'No actuator disk/blades defined'
            return
        endif
        !
        nrcsav = nrc
        if (nrcsav>0) then
            n = 1  ! assume at least one disk defined
            do ir = 1, nrcsav
                yrc0(ir) = yrc(ir, n)
                cldes0(ir) = cldes(ir)
            enddo
        endif
        !
        do n = 1, nrotor
            !
            !---- get radial location of rotor on CB and duct wall
            if (iprotcb(n)==0) then
                write (*, *) 'Error locating disk on CB wall'
                stop
            endif
            ipcb = iprotcb(n)
            xcb = xp(ipcb)
            ycb = yp(ipcb)
            rhub(n) = ycb
            if (ldbg) then
                write (lundbg, *) 'Disk on CB @IP= ', ipcb
                write (lundbg, *) 'Xhub,Rhub on CB ', xcb, rhub(n)
            endif
            !
            if (iprotdw(n)==0) then
                write (*, *) 'Error locating rotor on Duct wall'
                stop
            endif
            ipdw = iprotdw(n)
            xdw = xp(ipdw)
            ydw = yp(ipdw)
            rtip(n) = ydw
            if (ldbg) then
                write (lundbg, *) 'Disk on Duct @IP= ', ipdw
                write (lundbg, *) 'Xtip,Rtip on Duct wall ', xdw, rtip(n)
            endif
            !
            !---- Update swept disk area
            adisk(n) = pi * (rtip(n)**2 - rhub(n)**2)
            atip(n) = pi * (rtip(n)**2)
            if (ldbg) write (*, *) 'Disk # ', n, ' swept area ', &
                    & adisk(n)
            !
            !
            if (irtype(n)/=0) then ! disk defined
                !---- Save old blade and circulation distribution for interpolation to new pts
                nrc0 = nrcsav
                do ir = 1, nrc0
                    yrc0(ir) = yrc(ir, n)
                    bgam0(ir) = bgam(ir, n)
                    chr0(ir) = chr(ir, n)
                    betar0(ir) = betar(ir, n)
                    !          IF(N.EQ.1) CLDES0(IR) = CLDES(IR)
                enddo
            endif
            !
            !---- Set discrete points on rotor line using linear interpolation hub-tip
            nrp = nrsta
            dx = (xdw - xcb) / float(nrp - 1)
            dy = (ydw - ycb) / float(nrp - 1)
            do ir = 1, nrp
                xrp(ir, n) = xcb + float(ir - 1) * dx
                yrp(ir, n) = ycb + float(ir - 1) * dy
                !c        write(*,*) 'rotorpts ',i,xrp(i,n),yrp(i,n)
            enddo
            !
            !
            ! *** In what follows, condition on LLOFT is to allow paneling to geometric
            !     tip gap for lofting, while changing nothing for OPER.
            !
            !---- If effective tip gap is specified add two points to last interval
            !     at tip gap from duct wall, shift original rotor points inwards
            !
            if (lloft) then
                tgapzlt = 0.0
            else
                tgapzlt = tgapzl
            endif
            !
            tg = max(0.0, tgap - tgapzlt)
            if (tg>0.0) then
                dx = xdw - xcb
                dr = rtip(n) - rhub(n)
                do ir = 1, nrp
                    yrp(ir, n) = yrp(ir, n) - 2.0 * tg * (yrp(ir, n) - rhub(n)) / dr
                    xrp(ir, n) = xrp(ir, n) - dx * (rtip(n) - yrp(ir, n)) / dr
                    !c          write(*,*) 'rotorpts ',i,xrp(i,n),yrp(i,n)
                enddo
                !--- Add two points, one at RTIP-TGAP, one at shroud wall at RTIP
                nrp = nrp + 1
                yrp(nrp, n) = rtip(n) - tg
                xrp(nrp, n) = xdw - dx * (rtip(n) - yrp(nrp, n)) / dr
                nrp = nrp + 1
                xrp(nrp, n) = xdw
                yrp(nrp, n) = ydw
            endif
            !
            !---- Rotor interval centers and initial circulation
            nrc = nrp - 1
            do ir = 1, nrc
                xrc(ir, n) = 0.5 * (xrp(ir, n) + xrp(ir + 1, n))
                yrc(ir, n) = 0.5 * (yrp(ir, n) + yrp(ir + 1, n))
                bgam(ir, n) = 0.0
                !c        write(*,*) 'rotorcpts ',ir,xrc(ir,n),yrc(ir,n)
            enddo
            !
            !=====================================================================
            !
            if (irtype(n)==0) then ! disk undefined
                !
                !---- If actuator disk has been input interpolate BGAM to rotor center points
                if (irtypdef(n)==1 .and. nrdef(n)>=2) then
                    do i = 1, nrdef(n)
                        t1(i) = bgamdef(i, n)
                    enddo
                    !---- Spline rotor definition arrays
                    call segspl(t1, t1s, yrdef(1, n), nrdef(n))
                    !---- Set BGAM for points on rotor line
                    do ir = 1, nrc
                        bgam(ir, n) = seval(yrc(ir, n), t1, t1s, yrdef(1, n), &
                                & nrdef(n))
                    enddo
                    if (tgap>0.0) bgam(nrc, n) = 0.0
                    irtype(n) = 1
                    lbldef = .false.
                    lvmav = .false.
                endif
                !---- If a rotor has been input (Y,CH,BETA) interpolate CH,BETA to rotor
                !     center points
                if (irtypdef(n)==2 .and. nrdef(n)>=2) then
                    do i = 1, nrdef(n)
                        t1(i) = chrdef(i, n)
                        t2(i) = betadef(i, n)
                    enddo
                    !---- Spline rotor definition arrays
                    call segspl(t1, t1s, yrdef(1, n), nrdef(n))
                    call segspl(t2, t2s, yrdef(1, n), nrdef(n))
                    !---- Set CH,BETA for points on rotor line
                    do ir = 1, nrc
                        chr(ir, n) = seval(yrc(ir, n), t1, t1s, yrdef(1, n), nrdef(n)&
                                &)
                        betar(ir, n) = seval(yrc(ir, n), t2, t2s, yrdef(1, n), &
                                & nrdef(n))
                    enddo
                    irtype(n) = 2
                    lbldef = .true.
                    lvmav = .false.
                    call setiaero
                endif
                !
            else
                !
                !---- Interpolate old rotor/circulation to new rotor points
                do i = 1, nrc0
                    t1(i) = chr0(i)
                    t2(i) = betar0(i)
                    t3(i) = bgam0(i)
                enddo
                !---- Spline rotor definition arrays
                call segspl(t1, t1s, yrc0, nrc0)
                call segspl(t2, t2s, yrc0, nrc0)
                call segspl(t3, t3s, yrc0, nrc0)
                !---- Set CH,BETA,BGAM for points on new rotor line
                do ir = 1, nrc
                    chr(ir, n) = seval(yrc(ir, n), t1, t1s, yrc0, nrc0)
                    betar(ir, n) = seval(yrc(ir, n), t2, t2s, yrc0, nrc0)
                    bgam(ir, n) = seval(yrc(ir, n), t3, t3s, yrc0, nrc0)
                enddo
                !cc        IF(TGAP.GT.0.0) BGAM(NRC,N) = 0.0
                call setiaero
                !
            endif
            !
        enddo
        ! end NROTOR loop
        !
        !---- Move design CLs and BB data to new radial stations
        !
        if (nrcsav>0) then
            clp1 = clpos(1)
            clp2 = clpos(nrcsav)
            cln1 = clneg(1)
            cln2 = clneg(nrcsav)
            cld1 = cldes(1)
            cld2 = cldes(nrcsav)
            !
            do ir = 1, nrc
                clpos(ir) = clp1 + float(ir - 1) / float(nrc - 1) * (clp2 - clp1)
                clneg(ir) = cln1 + float(ir - 1) / float(nrc - 1) * (cln2 - cln1)
                cldes(ir) = cld1 + float(ir - 1) / float(nrc - 1) * (cld2 - cld1)
            enddo
            !
            do n = 1, nrotor
                if (lbbloft(n)) then
                    do ic = 1, nrcsav
                        t4(ic) = bbvfac(ic, n)
                    enddo
                    call segspl(t4, t4s, yrc0, nrcsav)
                    do ic = 1, nrc
                        bbvfac(ic, n) = seval(yrc(ic, n), t4, t4s, yrc0, nrcsav)
                    enddo
                endif
            enddo
        endif
        !
        !        DO I = 1, NRCSAV
        !          T4(I) = CLPOS(I)
        !        END DO
        !        CALL SEGSPL(T4,T4S,YRC0,NRCSAV)
        !        DO I = 1, NRC
        !          CLPOS(I) = SEVAL(YRC(I,1),T4,T4S,YRC0,NRCSAV)
        !        END DO
        !
        !        DO I = 1, NRCSAV
        !          T4(I) = CLNEG(I)
        !        END DO
        !        CALL SEGSPL(T4,T4S,YRC0,NRCSAV)
        !        DO I = 1, NRC
        !          CLNEG(I) = SEVAL(YRC(I,2),T4,T4S,YRC0,NRCSAV)
        !        END DO
        !
        !        DO I = 1, NRCSAV
        !          T4(I) = CLDES(I)
        !        END DO
        !        CALL SEGSPL(T4,T4S,YRC0,NRCSAV)
        !        DO I = 1, NRC
        !          CLDES(I) = SEVAL(YRC(I,1),T4,T4S,YRC0,NRCSAV)
        !        END DO
        !      ENDIF
        !
    end subroutine rotorinit
    !*==SETROTWAK.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ROTORINIT



    subroutine setrotwak
        use i_dfdc
        use m_spline, only : segspl, seval
        use m_aero, only : setiaero
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: iel, ig, ig1, ip, ip1, ip2, ir, n, nr
        real, dimension(irx) :: t1, t1s, t2, t2s, t3, t3s, y1
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Sets up rotor source line and wake elements from wake grid
        !------------------------------------------------------
        !
        !---- Add rotor source element after elements now in panel geometry
        !     (foils, drag area objects)
        !c      write(*,*) 'SETROTWAK entered NEL ',NEL
        !
        do nr = 1, nrotor
            !
            !---- Define paneled source line element for rotor lines from grid points
            ig = igrotor(nr)
            do ir = 1, nrp
                xrp(ir, nr) = xg(ig, ir)
                yrp(ir, nr) = yg(ig, ir)
            enddo
            call addelem6(xrp(1, nr), yrp(1, nr), nrp)
            ielrotor(nr) = nel
            !
            !---- Interpolate old rotor params/circulation to new rotor points
            nrc = nrp - 1
            do ir = 1, nrc
                y1(ir) = yrc(ir, nr)
                t1(ir) = chr(ir, nr)
                t2(ir) = betar(ir, nr)
                t3(ir) = bgam(ir, nr)
            enddo
            !---- Redefine rotor center points
            do ir = 1, nrc
                xrc(ir, nr) = 0.5 * (xrp(ir, nr) + xrp(ir + 1, nr))
                yrc(ir, nr) = 0.5 * (yrp(ir, nr) + yrp(ir + 1, nr))
            enddo
            !---- Spline rotor definition arrays
            call segspl(t1, t1s, y1, nrc)
            call segspl(t2, t2s, y1, nrc)
            call segspl(t3, t3s, y1, nrc)
            !---- Set CH,BETA,BGAM for points on new rotor line
            do ir = 1, nrc
                chr(ir, nr) = seval(yrc(ir, nr), t1, t1s, y1, nrc)
                betar(ir, nr) = seval(yrc(ir, nr), t2, t2s, y1, nrc)
                bgam(ir, nr) = seval(yrc(ir, nr), t3, t3s, y1, nrc)
            enddo
            !cc        IF(TGAP.GT.0.0) BGAM(NRC,NR) = 0.0
            call setiaero
            !
        enddo
        !
        !---- Add vortex wake elements to panel geometry
        !
        !----- set up wake element from centerbody TE (element 1)
        ir = 1
        iel = 1
        ig1 = igtecb
        n = ii - ig1 + 1
        call addwake(xg(ig1, ir), yg(ig1, ir), n)
        !---- save wake element index for this rotor point
        ir2iel(ir) = nel
        iel2ir(nel) = ir
        !---- plot Cp only on interior flow side
        isplot(nel) = 2
        iewake(iel) = nel
        !---- set up pointers from wake points to grid IG streamwise index
        ip1 = ipfrst(nel)
        ip2 = iplast(nel)
        do ip = ip1, ip2
            ip2ig(ip) = ig1 + ip - ip1
        enddo
        !
        !----- set up intermediate rotor wakes (from rotor line)
        do ir = 2, nrp - 1
            ig1 = 1
            n = ii
            call addwake(xg(ig1, ir), yg(ig1, ir), n)
            !---- save wake element index for this rotor point
            ir2iel(ir) = nel
            iel2ir(nel) = ir
            !---- plot Cp only on interior flow side
            isplot(nel) = -1
            !---- set up pointers from wake points to grid IG streamwise index
            ip1 = ipfrst(nel)
            ip2 = iplast(nel)
            do ip = ip1, ip2
                ip2ig(ip) = ig1 + ip - ip1
            enddo
        enddo

        !----- and wake from duct TE (element 2)
        ir = nrp
        iel = 2
        ig1 = igtedw
        n = ii - ig1 + 1
        call addwake(xg(ig1, ir), yg(ig1, ir), n)
        !---- save wake element index for this rotor point
        ir2iel(ir) = nel
        iel2ir(nel) = ir
        !---- plot Cp only on interior flow side
        isplot(nel) = 1
        iewake(iel) = nel
        !---- set up pointers from wake points to grid IG streamwise index
        ip1 = ipfrst(nel)
        ip2 = iplast(nel)
        do ip = ip1, ip2
            ip2ig(ip) = ig1 + ip - ip1
        enddo
        !
        !----------------------------------------------------------------
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
    end subroutine setrotwak
    !*==ADDWAKE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETROTWAK


    subroutine addwake(x, y, n)
        use i_dfdc
        use m_geom, only : xypspl
        use m_geutil, only : minmax
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(*) :: x, y
        !
        ! Local variables
        !
        integer :: i, iel, ip, ip1, ip2
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Add wake specified by N points X,Y to panel geometry
        !---------------------------------------------------------
        !
        !---- add wake to panel geometry
        ip1 = iplast(nel) + 1
        ip2 = ip1 - 1 + n
        iel = nel + 1
        if (ldbg) write (lundbg, 10) iel, n, ip1, ip2
        !---- add wake points
        do ip = ip1, ip2
            i = ip - ip1 + 1
            xp(ip) = x(i)
            yp(ip) = y(i)
            if (ldbg) write (lundbg, *) ip, xp(ip), yp(ip)
        enddo
        !---- centroid location
        xpcent(iel) = 0.5 * (xp(ip1) + xp(ip2))
        ypcent(iel) = 0.5 * (yp(ip1) + yp(ip2))
        !---- panel limits
        call minmax(n, xp(ip1), xpmine(iel), xpmaxe(iel))
        call minmax(n, yp(ip1), ypmine(iel), ypmaxe(iel))
        !
        !---- setup pointers and panel data
        ipfrst(iel) = ip1
        iplast(iel) = ip2
        netype(iel) = 7
        ltpan(iel) = .true.
        !      LTPAN(IEL) = .FALSE.
        ipte2(iel) = ip2
        ipte1(iel) = 0
        call xypspl(iel)
        !      CALL XYCSET(IEL)
        !
        nel = iel
        nptot = ip2
        !
        if (ldbg) write (*, 10) iel, n, ip1, ip2
        !
        10   format ('Vortex wake element ', i4, ' added with #pts ', i5, ' IP1 ', &
                & i5, ' IP2 ', i5)
        !
    end subroutine addwake
    !*==UPDROTWAK.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ADDWAKE



    subroutine updrotwak
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: i, i1, iel, ip, ip1, ip2, ir, n, nip
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Updates up rotor wake element geometry from wake grid
        !     Assumes that points have moved but indices remain the same
        !------------------------------------------------------
        !
        !c      write(*,*) 'UPDROTWAK entered NEL ',NEL
        !
        !----- update wake element from centerbody TE (element 1)
        ir = 1
        iel = ir2iel(ir)
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        i1 = igtecb
        n = ii - i1 + 1
        nip = ip2 - ip1 + 1
        if (n/=nip) then
            write (*, *) 'UPDROTWAK wake grid mismatch on IR ', ir, n, &
                    & nip
            stop
        endif
        i = i1
        do ip = ip1, ip2
            xp(ip) = xg(i, ir)
            yp(ip) = yg(i, ir)
            i = i + 1
        enddo
        !
        !----- set up intermediate rotor wakes (from rotor line)
        do ir = 2, nrp - 1
            iel = ir2iel(ir)
            ip1 = ipfrst(iel)
            ip2 = iplast(iel)
            i1 = 1
            n = ii
            nip = ip2 - ip1 + 1
            if (n/=nip) then
                write (*, *) 'UPDROTWAK wake grid mismatch on IR ', ir, n, &
                        & nip
                stop
            endif
            i = i1
            do ip = ip1, ip2
                xp(ip) = xg(i, ir)
                yp(ip) = yg(i, ir)
                i = i + 1
            enddo
        enddo

        !----- and wake from duct TE (element 2)
        ir = nrp
        iel = ir2iel(ir)
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        i1 = igtedw
        n = ii - i1 + 1
        nip = ip2 - ip1 + 1
        if (n/=nip) then
            write (*, *) 'UPDROTWAK wake grid mismatch on IR ', ir, n, &
                    & nip
            stop
        endif
        i = i1
        do ip = ip1, ip2
            xp(ip) = xg(i, ir)
            yp(ip) = yg(i, ir)
            i = i + 1
        enddo
        !
        !----------------------------------------------------------------
        !---- invalidate any existing solution
        lncvp = .false.
        lqaic = .false.
        lqgic = .false.
        lqcnt = .false.
        lgsys = .false.
        lgamu = .false.
        lgama = .false.
        !
    end subroutine updrotwak
    !*==ROTPINIT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! UPDROTWAK


    subroutine rotpinit
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: ic, ic1, ic2, iel, ig1, ig2, ip, ip1, ip2, &
                & ipcb, ipdw, ir, n, nr
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Sets up pointers for panel points on rotor wakes
        !     and pointers from panel centers to streamline grid
        !------------------------------------------------------
        !
        !---- Clear pointers to streamline rows from body and wake points
        !     Used for additional gamma from vortex wakes
        do ip = 1, ipx
            ip2ir(ip) = 0
        enddo
        !
        !---- Clear pointers to grid centers from body and wake centers
        do ic = 1, icx
            ic2ig(ic) = 0
        enddo
        !
        nr = 1
        !---- Pointers for vortex wake on centerbody
        ir = 1
        iel = 1
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        ipcb = iprotcb(nr)
        !c    IF(LDBG) WRITE(22,*) IR, IEL, IP1, IP2, IPCB
        do ip = ipcb, ip1, -1
            ip2ir(ip) = ir
            !c      IF(LDBG) WRITE(22,*) IP,IR
        enddo
        !---- set pointers from panel center to grid center
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        do ic = ic1, ic2
            ig1 = ip2ig(ipco(ic))
            ig2 = ip2ig(ipcp(ic))
            ic2ig(ic) = min(ig1, ig2)
            !c        IF(IC2IG(IC).NE.0) WRITE(22,*) IC,IC2IG(IC)
        enddo
        !
        !---- Set pointers for vortex wake elements
        do ir = 1, nrp
            iel = ir2iel(ir)
            if (iel<=0) then
                write (*, *) 'IR2IEL out of range in ROTPINIT ', iel
                stop
            endif
            ip1 = ipfrst(iel)
            ip2 = iplast(iel)
            do ip = ip1, ip2
                ip2ir(ip) = ir
            enddo
            !---- set pointers from panel center to grid center
            ic1 = icfrst(iel)
            ic2 = iclast(iel)
            do ic = ic1, ic2
                ig1 = ip2ig(ipco(ic))
                ig2 = ip2ig(ipcp(ic))
                ic2ig(ic) = min(ig1, ig2)
                !c          IF(IC2IG(IC).NE.0) WRITE(22,*) IC,IC2IG(IC)
            enddo
        enddo
        !
        !---- Pointers for vortex wake on duct
        ir = nrp
        iel = 2
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        n = ip2 - ip1 + 1
        ipdw = iprotdw(nr)
        !c    IF(LDBG) WRITE(22,*) IR, IEL, IP1, IP2, IPDW
        do ip = ipdw, ip2
            ip2ir(ip) = ir
            !c      IF(LDBG) WRITE(22,*) IP,IR
        enddo
        !---- set pointers from panel center to grid center
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        do ic = ic1, ic2
            ig1 = ip2ig(ipco(ic))
            ig2 = ip2ig(ipcp(ic))
            ic2ig(ic) = min(ig1, ig2)
            !c        IF(IC2IG(IC).NE.0) WRITE(22,*) IC,IC2IG(IC)
        enddo
        !
    end subroutine rotpinit
    !*==ADDELEM5_2PT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ROTPINIT



    subroutine addelem5_2pt(x1, y1, x2, y2, npts)
        use i_dfdc
        use m_geom, only : xypspl
        use m_geutil, only : minmax
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: npts
        real :: x1, x2, y1, y2
        !
        ! Local variables
        !
        real :: dx, dy, fnpt1
        integer :: i, iel, ip, ip1, ip2, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Add line source element (NETYPE=5) to paneled geometry
        !     Specified by X1,Y1 start and X2,Y2 end coordinates and
        !     # of points.  Intermediate ponts are linearly interpolated
        !     between endpoints.
        !------------------------------------------------------
        !
        if (npts<=1) then
            write (*, *) 'Error in ADDELEM12 NPTS<2'
            stop
        endif
        !
        ip = iplast(nel)
        iel = nel
        !
        ip1 = ip + 1
        ip2 = ip + npts
        !
        dx = x2 - x1
        dy = y2 - y1
        fnpt1 = float(npts - 1)
        do i = 1, npts
            ip = ip + 1
            xp(ip) = x1 + dx * float(i - 1) / fnpt1
            yp(ip) = y1 + dy * float(i - 1) / fnpt1
        enddo
        !---- centroid location
        xpcent(iel) = 0.5 * (xp(ip1) + xp(ip2))
        ypcent(iel) = 0.5 * (yp(ip1) + yp(ip2))
        !---- panel limits
        call minmax(n, xp(ip1), xpmine(iel), xpmaxe(iel))
        call minmax(n, yp(ip1), ypmine(iel), ypmaxe(iel))
        !
        !--- add element
        iel = iel + 1
        ipfrst(iel) = ip1
        iplast(iel) = ip2
        netype(iel) = 5
        ltpan(iel) = .false.
        call xypspl(iel)
        !
        if (ldbg) then
            write (lundbg, *) ' '
            write (lundbg, *) 'Adding line source element ', iel
            write (lundbg, 100) 'Start ', ip1, x1, y1
            write (lundbg, 100) 'End   ', ip2, x2, y2
            !
            write (*, 10) iel, npts, ip1, ip2
        endif
        !
        10   format ('Line source element ', i4, ' added with #pts ', i5, ' IP1 ', &
                & i5, ' IP2 ', i5)
        100  format (a, i5, 4(1x, g13.5))
        !
        nel = iel
        nptot = ip2
        !
    end subroutine addelem5_2pt
    !*==ADDELEM5.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ADDELEM5_2PT



    subroutine addelem5(x, y, n)
        use i_dfdc
        use m_geom, only : xypspl
        use m_geutil, only : minmax
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(*) :: x, y
        !
        ! Local variables
        !
        integer :: i, iel, ip, ip1, ip2
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Add line source element (NETYPE=5) to paneled geometry
        !     Specified by X,Y coordinates and # of points
        !------------------------------------------------------
        !
        if (n<=1) then
            write (*, *) 'Error in ADDELEM5 NPTS<2'
            stop
        endif
        !
        !---- add sourceline points to panel geometry
        ip1 = iplast(nel) + 1
        ip2 = ip1 - 1 + n
        iel = nel + 1
        if (ldbg) write (lundbg, 10) iel, n, ip1, ip2
        !
        !---- add points to panel arrays
        do ip = ip1, ip2
            i = ip - ip1 + 1
            xp(ip) = x(i)
            yp(ip) = y(i)
            if (ldbg) write (lundbg, *) ip, xp(ip), yp(ip)
        enddo
        !---- centroid location
        xpcent(iel) = 0.5 * (xp(ip1) + xp(ip2))
        ypcent(iel) = 0.5 * (yp(ip1) + yp(ip2))
        !---- panel limits
        call minmax(n, xp(ip1), xpmine(iel), xpmaxe(iel))
        call minmax(n, yp(ip1), ypmine(iel), ypmaxe(iel))
        !
        !--- add element
        ipfrst(iel) = ip1
        iplast(iel) = ip2
        netype(iel) = 5
        ltpan(iel) = .false.
        call xypspl(iel)
        !
        10   format (/'Line source element ', i4, ' added with #pts ', i5, ' IP1 ', &
                & i5, ' IP2 ', i5)
        100  format (a, i5, 4(1x, g13.5))
        !
        nel = iel
        nptot = ip2
        !
    end subroutine addelem5
    !*==ADDELEM6.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ADDELEM5



    subroutine addelem6(x, y, n)
        use i_dfdc
        use m_geom, only : xypspl
        use m_geutil, only : minmax
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(*) :: x, y
        !
        ! Local variables
        !
        integer :: i, iel, ip, ip1, ip2
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Add element for rotor source line to panel geometry.
        !     Specified by N points X,Y.
        !     Defined as TYPE=6 element
        !---------------------------------------------------------
        !
        !---- add wake to panel geometry
        ip1 = iplast(nel) + 1
        ip2 = ip1 - 1 + n
        iel = nel + 1
        if (ldbg) write (lundbg, 10) iel, n, ip1, ip2
        !
        !---- add points to panel arrays
        do ip = ip1, ip2
            i = ip - ip1 + 1
            xp(ip) = x(i)
            yp(ip) = y(i)
            if (ldbg) write (lundbg, *) ip, xp(ip), yp(ip)
        enddo
        !---- centroid location
        xpcent(iel) = 0.5 * (xp(ip1) + xp(ip2))
        ypcent(iel) = 0.5 * (yp(ip1) + yp(ip2))
        !---- panel limits
        call minmax(n, xp(ip1), xpmine(iel), xpmaxe(iel))
        call minmax(n, yp(ip1), ypmine(iel), ypmaxe(iel))
        !
        !---- setup pointers and panel data
        ipfrst(iel) = ip1
        iplast(iel) = ip2
        netype(iel) = 6
        ltpan(iel) = .false.
        call xypspl(iel)
        !      CALL XYCSET(IEL)
        !
        nel = iel
        nptot = ip2
        !
        10   format ('Rotor source line element ', i4, ' added with #pts ', i5, &
                &' IP1 ', i5, ' IP2 ', i5)
        !
    end subroutine addelem6
    !*==DFLOAD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ADDELEM6




    subroutine dfload(fnamin, ferror)
        use i_dfdc
        use m_userio, only : getflt, rdline, strip
        use m_pnsubs, only : panget
        use m_airio, only : areadnr
        use m_atmo, only : atmo
        use m_aero, only : putaero
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: itx = 300
        !
        ! Dummy arguments
        !
        logical :: ferror
        character(*) :: fnamin
        !
        ! Local variables
        !
        real :: a0, a0deg, cdmin, cldmin, clmax, clmin, cmcon, &
                & dcdcl2, dcdcl2s, dclda, dclda_stall, dcl_stall, &
                & filevers, mcrit, reref, rexp, rpm1, rpm2, toc, &
                & xisect
        logical :: error, lopen
        character(128) :: fname, line
        integer :: i, ib, ibnext, icnt, id, iel, iftype, ir, &
                & irpn, it, iv, k, lu, n, nbin, nd1, nf, &
                & ninput, nr, nrpnin, nrpnrd, ntel
        integer, dimension(nex) :: nt
        real, dimension(10) :: rinput
        real, dimension(itx, nex) :: xt, yt
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Reads previously saved duct case data and geometry
        !     from file FNAMIN in DFDC version 0.5 format.  This
        !     format separates data into groups by keywords.
        !     Once read the geometry is processed from buffer to
        !     paneled geometry.
        !---------------------------------------------------------
        !
        !---- local arrays for calling AREAD
        !
        !
        !---- reset some flags for new dataset
        ferror = .false.
        nrotor = 0
        lbldef = .false.
        lconv = .false.
        nrp = 0
        nrc = 0
        do n = 1, nrx
            irtype(n) = 0
            nrdef(n) = 0
            irtypdef(n) = 0
        enddo
        !
        lrspcdef = .false.
        lrspcusr = .false.
        !
        ndobj = 0
        !
        lopen = .false.
        lu = 1
        fname = fnamin
        call strip(fname, nf)
        if (fname/=' ') then
            open (lu, file = fname, status = 'OLD', err = 98)
            lopen = .true.
        else
            return
        endif
        !
        icnt = 0
        !
        !---- Read header line and version from DFDC input file
        call rdline(lu, line, icnt)
        if (line=='END' .or. line=='ERR') goto 210
        !
        if (line(1:4)/='DFDC') then
            write (*, *) 'Not DFDC input file...may be duct geometry file'
            close (lu)
            !c        CALL LOADG(FNAME)
            return
        endif
        iv = index(line, 'Version')
        read (line(iv + 8:len(line)), *) filevers
        write (*, 1005) filevers
        !
        1000 format (a)
        1005 format (/, ' Reading file from DFDC Version ', e8.2)
        1010 format (' File  ', a, ' not found'/)
        1020 format (' File  ', a, &
                &' has incompatible format'/' Loading not completed'/)
        1025 format (' Case Name: ', a)
        !
        !
        !--- Get Case title from line #2
        call rdline(lu, line, icnt)
        if (line=='END' .or. line=='ERR') goto 210
        name = line
        call strip(name, nname)
        !
        write (*, 1025) name
        do
            !
            !-----------------------------------------------------
            !--- Main loop to find keywords and start read actions
            !    for data groups
            call rdline(lu, line, icnt)
            if (line=='END') then
                !
                !------------------------------------------------------------------------
                nbel = ntel
                !
                !---- Process airfoil inputs for CB and duct
                !---- put element-indexed coordinates XT,YT into linear bufffer airfoil
                !     geometry arrays XB,YB
                ibnext = 1
                do iel = 1, nbel
                    !------ first point in element IEL
                    ibfrst(iel) = ibnext
                    !------ initalize accumulation counter
                    ib = ibfrst(iel) - 1
                    !------ go over all points on element IEL
                    if (ldbg) write (lundbg, *) 'Airfoil element'
                    do it = 1, nt(iel)
                        ib = ib + 1
                        if (ib>ibx) stop 'LOAD: Array overflow on IBX'
                        xb(ib) = xt(it, iel)
                        yb(ib) = yt(it, iel)
                        if (ldbg) write (lundbg, *) it, xt(it, iel), &
                                & yt(it, iel)
                    enddo
                    !---- set buffer airfoils to NBTYPE=0
                    nbtype(iel) = 0
                    iblast(iel) = ib
                    ibnext = ib + 1
                enddo
                !---- set total number of points in buffer geometry
                nbtot = iblast(nbel)
                !
                !
                !=======================================================================
                !---- If panel spacing data is not in case file try reading explicitly
                !     from xxx.pan paneling input file
                if (.not.lrspcdef) then
                    !---- get file prefix (if it has a .xxx suffix)
                    k = index(fname, '.') - 1
                    if (k<=0) then
                        !----- no "." prefix/suffix separator
                        !      just tack on ".pan" to full filename
                        call strip(fname, nf)
                        k = nf
                    endif
                    prefix = fname(1:k)
                    nprefix = k
                    pfile = prefix(1:nprefix) // '.pan'
                    call panget(pfile, error)
                    !
                    if (error) then
                        write (*, 2100) pfile(1:k + 4)
                    else
                        write (*, 2110) pfile(1:k + 4)
                    endif
                endif
                !
                2100       format (' No repaneling parameters from file ', a)
                2110       format (' Repaneling parameters from file ', a)

                !
                !=======================================================================
                !---- process current buffer geometry to identify element types and
                !     set element parameters
                call geproc
                if (ldbg) call geprint(lundbg)
                !
                !---- process buffer foil geometry and set up paneled geometry and grid
                !c    CALL GENGEOM
                !
                !---- take note of new case
                lload = .true.
                do i = 1, nrotor
                    lbbloft(i) = .false.
                enddo
                goto 300
            elseif (line=='ERR') then
                goto 210
                !
                !
                !---- OPER data
            elseif (index(line, 'OPER')/=0) then
                call rdline(lu, line, icnt)
                !c        READ(LINE,*,ERR=210) QINF,QREF,RPM
                ninput = 4
                call getflt(line, rinput, ninput, error)
                if (ninput>=4) then
                    qinf = rinput(1)
                    qref = rinput(2)
                    rpm1 = rinput(3)
                    rpm2 = rinput(4)
                elseif (ninput>=3) then
                    qinf = rinput(1)
                    qref = rinput(2)
                    rpm1 = rinput(3)
                    rpm2 = 0.0
                elseif (ninput>=2) then
                    qinf = rinput(1)
                    qref = rinput(2)
                    rpm1 = 0.0
                    rpm2 = 0.0
                endif
                !
                call rdline(lu, line, icnt)
                read (line, *, err = 210) rho, vso, rmu, alth
                call rdline(lu, line, icnt)
                read (line, *, err = 210) xdwklen, nwake
                call rdline(lu, line, icnt)
                read (line, *, err = 210) lwrlx
                call rdline(lu, line, icnt)
                if (index(line, 'ENDOPER')/=0) then
                    write (*, 1030) qinf, qref, rpm1, rpm2, rho, vso, &
                            & rmu, alth, xdwklen, nwake
                else
                    write (*, *) 'No ENDOPER to OPER section'
                    stop
                endif
                !
                1030       format (' OPER data read from file', /, '  Qinf   = ', f9.3, 5x, &
                        &'Qref   = ', f9.3, /, '  RPM1   = ', f9.1, 5x, &
                        & 'RPM2   = ', f9.1, /, '  Rho    = ', f9.5, 5x, &
                        & 'VSound = ', f9.3, /, '  rMU    = ', e9.3, 5x, &
                        & 'Alt    = ', f9.5, /, '  Xdwake = ', f9.5, 5x, &
                        & 'Nwake  = ', i6)
                !
                !---- Input RPM's set speed of disks
                omega(1) = pi * rpm1 / 30.0
                omega(2) = pi * rpm2 / 30.0
                if (rho==0.0 .or. vso==0.0 .or. rmu==0.0) then
                    deltat = 0.0
                    call atmo(alth, deltat, vso, rho, rmu)
                endif
                !
                !
                !---- AERO data
            elseif (index(line, 'AERO')/=0) then
                !--- Read aero section definitions
                nr = nrotor + 1
                call rdline(lu, line, icnt)
                read (line, *, err = 210) naero(nr)
                do n = 1, naero(nr)
                    call rdline(lu, line, icnt)
                    read (line, *, err = 210) xisect
                    call rdline(lu, line, icnt)
                    read (line, *, err = 210) a0deg, dclda, clmax, clmin
                    call rdline(lu, line, icnt)
                    read (line, *, err = 210) dclda_stall, dcl_stall, cmcon, &
                            & mcrit
                    call rdline(lu, line, icnt)
                    read (line, *, err = 210) cdmin, cldmin, dcdcl2
                    call rdline(lu, line, icnt)
                    !c          READ(LINE,*,ERR=210) REREF,REXP
                    ninput = 4
                    call getflt(line, rinput, ninput, error)
                    if (ninput>=4) then
                        reref = rinput(1)
                        rexp = rinput(2)
                        toc = rinput(3)
                        dcdcl2s = rinput(4)
                    elseif (ninput>=3) then
                        reref = rinput(1)

                        rexp = rinput(2)
                        toc = rinput(3)
                        dcdcl2s = 0.020
                    elseif (ninput>=2) then
                        reref = rinput(1)
                        rexp = rinput(2)
                        toc = 0.1
                        dcdcl2s = 0.020
                    endif
                    a0 = a0deg * dtr
                    call putaero(nr, n, xisect, a0, clmax, clmin, dclda, &
                            & dclda_stall, dcl_stall, cdmin, cldmin, dcdcl2, &
                            & dcdcl2s, cmcon, mcrit, toc, reref, rexp)
                enddo
                call rdline(lu, line, icnt)
                if (index(line, 'ENDAERO')/=0) then
                    write (*, 2010) nr, naero(nr)
                else
                    write (*, *) 'No ENDAERO to AERO section'
                    stop
                endif
                !
                2010       format (/, ' AERO data read for Disk', i2, ': ', i3, ' sections')
                !
                !
                !---- ROTOR data
            elseif (index(line, 'ROTOR')/=0) then
                nr = nrotor + 1
                call rdline(lu, line, icnt)
                read (line, *, err = 210) xdisk(nr), nrbld(nr), nrsta
                call rdline(lu, line, icnt)
                read (line, *, err = 210) nrdef(nr)
                !--- Get rotor definition data
                do ir = 1, nrdef(nr)
                    call rdline(lu, line, icnt)
                    read (line, *, end = 210, err = 210) yrdef(ir, nr), &
                            & chrdef(ir, nr), betadef(ir, nr)
                enddo
                !
                call rdline(lu, line, icnt)
                if (index(line, 'ENDROTOR')/=0) then
                    write (*, 2020) xdisk(nr), nrbld(nr), nrsta
                    if (nrdef(nr)>=2) then
                        write (*, 2030) nrdef(nr)
                        do ir = 1, nrdef(nr)
                            write (*, 15) yrdef(ir, nr), chrdef(ir, nr), &
                                    & betadef(ir, nr)
                            betadef(ir, nr) = betadef(ir, nr) * dtr
                        enddo
                        15               format (1x, f12.6, 2x, f12.6, 4x, f12.6)
                        irtypdef(nr) = 2
                        nrotor = nr
                    else
                        write (*, *) 'Rotor blade defined by too few stations'
                    endif
                else
                    write (*, *) 'No ENDROTOR found'
                    stop
                endif
                !
                2020       format (/, ' ROTOR data read from file', /, '  XDisk: ', f5.3, &
                        &'   Blades:', i2, '   Rotor pnts:', i3)
                2030       format ('  Rotor blade defined with', i3, ' points', /, &
                        &'         R            CH            Beta')
                !
                !
                !---- ACTDISK data
            elseif (index(line, 'ACTDISK')/=0) then
                nr = nrotor + 1
                call rdline(lu, line, icnt)
                read (line, *, err = 210) xdisk(nr), nrsta
                call rdline(lu, line, icnt)
                read (line, *, err = 210) nrdef(nr)
                !--- Get actuator disk definition data
                do ir = 1, nrdef(nr)
                    call rdline(lu, line, icnt)
                    read (line, *, end = 210, err = 210) yrdef(ir, nr), &
                            & bgamdef(ir, nr)
                enddo
                !
                call rdline(lu, line, icnt)
                if (index(line, 'ENDACTDISK')/=0) then
                    write (*, 2040) xdisk(nr), nrsta
                    if (nrdef(nr)>=2) then
                        write (*, 2050) nrdef(nr)
                        do ir = 1, nrdef(nr)
                            write (*, 15) yrdef(ir, nr), bgamdef(ir, nr)
                        enddo
                        irtypdef(nr) = 1
                        nrotor = nr
                    else
                        write (*, *)                                           &
                                &'Actuator disk defined by too few stations'
                    endif
                else
                    write (*, *) 'No ENDACTDISK found'
                    stop
                endif
                !
                2040       format (/, ' ACTDISK data read from file', /, '  XDisk: ', f5.3, &
                        &'   Rotor points:', i3)
                2050       format ('  Actuator Disk defined with', i3, ' points', /, &
                        &'         R           BGAM')
                !
                !
                !---- DRAGOBJ data
            elseif (index(line, 'DRAGOBJ')/=0) then
                nd1 = ndobj + 1
                if (nd1>ndrgx) then
                    write (*, *)                                              &
                            &'Number of drag objects exceeds dimension NDGX'
                else
                    call rdline(lu, line, icnt)
                    read (line, *, err = 210) nddef(nd1)
                    !--- Get drag object definition data
                    do id = 1, nddef(nd1)
                        call rdline(lu, line, icnt)
                        read (line, *, end = 210, err = 210) xddef(id, nd1), &
                                & yddef(id, nd1), cdadef(id, nd1)
                    enddo
                    !
                    call rdline(lu, line, icnt)
                    if (index(line, 'ENDDRAGOBJ')/=0) then
                        ndobj = nd1
                        write (*, 2060) nd1
                        if (nddef(nd1)>=2) then
                            write (*, 2070) nddef(nd1)
                            do id = 1, nddef(nd1)
                                write (*, 15) xddef(id, nd1), yddef(id, nd1), &
                                        & cdadef(id, nd1)
                            enddo
                        else
                            write (*, *)                                        &
                                    &'Drag object defined by too few stations'
                        endif
                    else
                        write (*, *) 'No ENDDRAGOBJ found'
                        stop
                    endif
                    ! ENDDRAGOBJ check
                    !
                endif
                ! NDRGX check
                !
                2060       format (/, ' DRAGOBJ ', i2, '  read from file')
                2070       format ('  Drag object defined with', i3, ' points', /, &
                        &'         X             R              CDA')
                !
                !
                !---- GEOM data (CB and duct coordinates in a single XFOIL multi-element file)
            elseif (index(line, 'GEOM')/=0) then
                !
                !---- Read the combined CB and duct airfoil file
                !
                call areadnr(lu, itx, nex, xt, yt, nt, ntel, aname, ispars, iftype)
                if (iftype==0 .or. ntel/=2) then
                    !----- read error occurred for two elements
                    write (*, 2080) ntel
                    stop
                endif
                !
                2080       format (' Error reading GEOM data', /, &
                        & ' Coordinates read for', i2, ' foils')
                !
                !---- User-defined respacing specs
                !---- PANE data (respacing data for surface points
            elseif (index(line, 'PANE')/=0) then
                call rdline(lu, line, icnt)
                read (line, *, err = 210) nbin, nrpnin
                nrpnrd = min(nrpnin, nrpnx)
                do iel = 1, nbin
                    call rdline(lu, line, icnt)
                    read (line, *, err = 210) npan(iel)
                    call rdline(lu, line, icnt)
                    read (line, *, err = 210) cvex(iel), smof(iel), fsle(iel), &
                            & fste(iel)
                    do irpn = 1, nrpnrd
                        call rdline(lu, line, icnt)
                        read (line, *, err = 210) srpn1(irpn, iel), &
                                & srpn2(irpn, iel), &
                                & crrat(irpn, iel)
                    enddo
                enddo
                call rdline(lu, line, icnt)
                if (index(line, 'ENDPANE')/=0) then
                    write (*, 2090) nbin
                    lrspcdef = .true.
                    lrspcusr = .true.
                endif
                !
            endif
            !
            2090    format (' PANELING for', i2, ' elements read from file')
        enddo
        !...............................................................
        98   write (*, 1050) fname(1:nf)
        ferror = .true.
        goto 300
        !
        write (*, 1100) fname(1:nf)
        ferror = .true.
        goto 300
        !
        210  write (*, 1150) icnt, fname(1:nf)
        ferror = .true.
        !
        1050 format (/' File OPEN error:  ', a)
        1100 format (/' File READ error:  ', a)
        1150 format (/' File READ error on line ', i3, ':  ', a)
        !
        !---- close file
        300  if (lopen) close (lu)
    end subroutine dfload
    !*==DFSAVE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! DFLOAD



    subroutine dfsave(fnamein)
        use i_dfdc
        use m_userio, only : asks, strip
        use m_aero, only : getaero
        use m_airio, only : awrite
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: fnamein
        !
        ! Local variables
        !
        real :: a0, a0deg, betadeg, cdmin, cldmin, clmax, clmin, &
                & cmcon, dcdcl2, dcdcl2s, dclda, dclda_stall, &
                & dcl_stall, mcrit, reref, rexp, rpm1, rpm2, toc, &
                & xisect
        character(1) :: ans
        character(80) :: fblnk
        character(128) :: fname
        integer :: i, iel, iftype, irpn, lu, n, nelsav, nf, nr
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------------------------
        !     Save rotor and operating state in DFDC Version 0.5 format
        !     This format saves data in groups delimited by keywords
        !     to separate data
        !--------------------------------------------------------------------------
        !
        lu = 2
        fname = fnamein
        call strip(fname, nf)
        !
        if (fname(1:1)==' ') call asks('Enter filename^', fname)
        open (lu, file = fname, status = 'OLD', err = 5)
        write (*, *)
        write (*, *) 'Output file exists.  Overwrite?  Y'
        read (*, 1000) ans
        if (index('Nn', ans)==0) goto 6
        !
        close (lu)
        write (*, *) 'Current duct case not saved.'
        return
        !
        5    open (lu, file = fname, status = 'NEW', err = 90)
        6    rewind (lu)
        !
        !
        !--- Version header and case name
        if (name==' ') name = 'Saved ducted fan'
        write (lu, 1100) version, name
        !
        !--- OPER data
        !
        if (nrotor<2) omega(2) = 0.
        !
        !--- Velocity, reference velocity and RPM
        write (lu, 1102)
        rpm1 = 30.0 * omega(1) / pi
        rpm2 = 30.0 * omega(2) / pi
        write (lu, 1200) qinf, qref, rpm1, rpm2
        !--- Altitude and atmospheric data
        write (lu, 1103)
        write (lu, 1200) rho, vso, rmu, alth
        !--- XDwake, #wake points
        write (lu, 1104)
        write (lu, 1202) xdwklen, nwake
        !--- Wake relaxation flag
        write (lu, 1105)
        write (lu, 1204) lwrlx
        write (lu, 1106)
        !
        !--- Save data for each disk (act disk or bladed)
        do nr = 1, nrotor
            !
            !--- AERO data
            !--- Save aero data for defined aero sections
            write (lu, 1108) naero(nr)
            do n = 1, naero(nr)
                call getaero(nr, n, xisect, a0, clmax, clmin, dclda, dclda_stall, &
                        & dcl_stall, cdmin, cldmin, dcdcl2, dcdcl2s, cmcon, &
                        & mcrit, toc, reref, rexp)
                write (lu, 1200) xisect
                a0deg = a0 * 180.0 / pi
                write (lu, 1109)
                write (lu, 1200) a0deg, dclda, clmax, clmin
                write (lu, 1110)
                write (lu, 1200) dclda_stall, dcl_stall, cmcon, mcrit
                write (lu, 1111)
                write (lu, 1200) cdmin, cldmin, dcdcl2
                write (lu, 1112)
                write (lu, 1200) reref, rexp, toc, dcdcl2s
            enddo
            write (lu, 1113)
            !
            !--- Save rotor blade if defined
            if (irtype(nr)==2) then
                write (lu, 1130)
                !--- Rotor axial location, #blades, #radial stations
                write (lu, 1202) xdisk(nr), nrbld(nr), nrsta
                !--- Save blade definition with chord,twist and body velocity
                write (lu, 1132) nrc
                do i = 1, nrc
                    betadeg = betar(i, nr) / dtr
                    write (lu, 1200) yrc(i, nr), chr(i, nr), betadeg
                enddo
                write (lu, 1134)
            endif
            !
            !--- Save actuator disk if defined
            if (irtype(nr)==1) then
                write (lu, 1122)
                !--- Rotor axial location, #radial stations
                write (lu, 1202) xdisk(nr), nrsta
                !--- Save actuator disk circulation
                write (lu, 1124) nrc
                do i = 1, nrc
                    write (lu, 1200) yrc(i, nr), bgam(i, nr)
                enddo
                write (lu, 1126)
            endif
            !
        enddo
        !
        !--- Save drag object(s) if defined
        if (ndobj>0) then
            do n = 1, ndobj
                if (nddef(n)>2) then
                    write (lu, 1136) nddef(n)
                    !--- Save drag object points and CDA values
                    do i = 1, nddef(n)
                        write (lu, 1200) xddef(i, n), yddef(i, n), cdadef(i, n)
                    enddo
                    write (lu, 1137)
                endif
            enddo
        endif
        !
        !--- CB and duct geometry
        write (lu, 1140)
        !---- save only 2 elements in buffer geometry into opened file
        nelsav = 2
        iftype = 2
        fblnk = ' '
        call awrite(fblnk, lu, nelsav, ipfrst, iplast, xp, yp, aname, ispars, &
                & iftype)
        write (lu, 1142)
        !
        !--- Save user-defined paneling specs
        if (lrspcusr) then
            write (lu, 1150) nbel, nrpnx
            do iel = 1, nbel
                write (lu, 1152) npan(iel), cvex(iel), smof(iel), &
                        & fsle(iel), fste(iel)
                write (lu, 1154)
                do irpn = 1, nrpnx
                    write (lu, 1200) srpn1(irpn, iel), srpn2(irpn, iel), &
                            & crrat(irpn, iel)
                enddo
            enddo
            write (lu, 1156)
        endif
        !
        close (lu)
        return
        !
        90   write (*, *) 'Bad filename.'
        write (*, *) 'Current duct case not saved.'
        return
        !
        !...................................................................
        1000 format (a)
        1100 format ('DFDC Version ', e8.2, ' '/a32)
        !
        1102 format (/'OPER', /'!   Vinf         Vref         RPM1         RPM2'&
                &)
        1103 format ('!   Rho          Vso          Rmu          Alt')
        1104 format ('!  XDwake             Nwake')
        1105 format ('!        Lwkrlx')
        1106 format ('ENDOPER')
        !
        1108 format (/'AERO', /'!  #sections'/1(1x, i5), /'!  Xisection')
        1109 format ('!  A0deg        dCLdA        CLmax        CLmin')
        1110 format ('! dCLdAstall   dCLstall      Cmconst      Mcrit')
        1111 format ('!   CDmin      CLCDmin       dCDdCL^2')
        1112 format ('!   REref        REexp        TOC         dCDdCL^2')
        1113 format ('ENDAERO')
        !
        1122 format (/'ACTDISK', /'!  Xdisk       NRsta')
        1124 format ('!  #stations'/1(1x, i5), /'!           r         BGam')
        1126 format ('ENDACTDISK')
        !
        1130 format (/'ROTOR', /'!  Xdisk               Nblds       NRsta')
        1132 format ('!  #stations'/1(1x, i5), &
                &/'!     r        Chord         Beta')
        1134 format ('ENDROTOR')
        !
        1136 format (/'DRAGOBJ', /'!  #pts'/1(1x, i5), &
                &/'!     x            r            CDA')
        1137 format ('ENDDRAGOBJ')
        !
        1140 format (/'GEOM')
        1142 format ('ENDGEOM')
        !
        1150 format (/'PANELING', /'!  #elements   #refinement zones'/(2(1x, i5))&
                &)
        1152 format ('!  #panel nodes'/1(1x, i5)                                &
                &/'!  curv_expon  curv_smooth   dsL/dsAvg    dsR/dsAvg', &
                & /4(1x, g12.5))
        1154 format ('!  s1/smax      s2/smax     ds/dsAvg')
        1156 format ('ENDPANELING')
        !
        1200 format (5(1x, g12.5))
        1202 format (1(1x, g12.5), 2(8x, i5))
        1204 format (12x, l1)
        !
        !x123456789012x123456789012x123456789012x123456789012x123456789012
        !!         Rho          Vso          Rmu           Alt')
        !
    end subroutine dfsave
end module m_dfdcsubs
