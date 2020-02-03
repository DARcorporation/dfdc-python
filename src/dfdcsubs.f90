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
    SUBROUTINE DFINIT(LDEBUG)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        LOGICAL :: LDEBUG
        !
        ! Local variables
        !
        REAL :: A0, CDMIN, CLDMIN, CLMAX, CLMIN, CMCON, DCDCL2, &
                & DCDCL2S, DCLDA, DCLDA_STALL, DCL_STALL, MCRIT, &
                & REREF, REXP, RPM, TOC, XISECT
        CHARACTER(10), SAVE :: DIGITS
        INTEGER :: I, IB, IC, IEL, IP, IPACT, IR, IRPN, K1, &
                & K2, KP, L, N, NPOL
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------
        !     Sets all initial defaults
        !--------------------------------------
        DATA DIGITS/'0123456789'/
        !
        !---- Set version (in case this code is used by DLL w/o main DFDC routine)
        !
        !      VERSION  =  0.70e02          ! This gets saved to file
        !      SVERSION = '070-ES2'         ! A more flexible format
        !
        !---------------------------------------------------------------
        !---- Debug output
        !cc      LDBG = .FALSE.
        LDBG = LDEBUG
        !---- Logical unit for debugging printout
        LUNDBG = 40
        !---------------------------------------------------------------
        !
        !---- set digit string array for constructing filename
        DO KP = 1, 99
            K2 = KP / 10 + 1
            K1 = KP - 10 * (KP / 10) + 1
            PNUM(KP) = DIGITS(K2:K2) // DIGITS(K1:K1)
        ENDDO
        !
        !---------------------------------------------------------------
        !
        DO IP = 1, IPX
            GAM(IP) = 0.
            SIG(IP) = 0.
            GTH(IP) = 0.
            GAMVSP(IP) = 0.
            SIGVSP(IP) = 0.
            VMAVG(IP) = 0.
            IZERO(IP) = 0
            IP2IG(IP) = 0
        ENDDO
        DO IC = 1, ICX
            IC2IG(IC) = 0
        ENDDO
        DO IEL = 1, NEX
            GAMSET(IEL) = 0.
            SIGSET(IEL) = 0.
        ENDDO
        DO IB = 1, IBX
            ONE(IB) = 1.0
            ZER(IB) = 0.
        ENDDO
        !
        !---- no panel nodes, control points
        NBTOT = 0
        NPTOT = 0
        NCTOT = 0
        DO IEL = 1, NEX
            IBFRST(IEL) = 0
            IBLAST(IEL) = -1
            IPFRST(IEL) = 0
            IPLAST(IEL) = -1
            ICFRST(IEL) = 0
            ICLAST(IEL) = -1
            NBTYPE(IEL) = -999
            !
            IEL2IR(IEL) = 0
        ENDDO
        !
        !---------------------------------------------------------------
        !---- no flow yet
        QINF = 0.0
        !---- start off with incompressible flow
        MACH = 0.
        MACH1 = 0.
        !---- acceleration due to gravity for scaling centrifugal blade tension (m/s^2)
        GEE = 9.81
        !---- default unity reference velocity
        QREF = 1.0
        !---- setup for SL conditions, Std temperature
        ALTH = 0.0
        DELTAT = 0.0
        !! sea level atmosphere parameters
        !     RHO =  1.226      ! fluid density         kg/m**3
        !     RMU =  1.78E-05   ! dynamic viscosity     kg/m-s
        !     VSO =  340.0      ! speed of sound        m/s
        CALL ATMO(ALTH, DELTAT, VSO, RHO, RMU)
        !---------------------------------------------------------------
        !
        !---- no drag objects defined
        NDOBJ = 0
        IELDRGOBJ(1) = 0
        LDRGOBJ = .FALSE.
        !
        !---------------------------------------------------------------
        !---- no rotor yet
        NROTOR = 0
        !---- no blade defined yet
        LBLDEF = .FALSE.
        !
        DO N = 1, NRX
            !---- disk type set to undefined
            IRTYPE(N) = 0
            !---- no defining points for rotor or actuator disk
            NRDEF(N) = 0
            IRTYPDEF(N) = 0
            !---- no elements defined for rotor source line
            IELROTOR(N) = 0
            !
            ATIP(N) = 0.0
            ADISK(N) = 0.0
            XDISK(N) = 999.  ! default rotor disk location
            NRBLD(N) = 5
            OMEGA(N) = 0.0
        ENDDO
        !---- default blade is rotor with RPM=1000
        RPM = 1000.0
        OMEGA(1) = PI * RPM / 30.0
        !
        TGAP = 0.0  ! no tip gap on rotor
        TGAPZL = 0.0  ! tip gap of zero loss (default is zero)
        !
        NRP = 0
        NRC = 0
        !
        NRSTA = 11   ! default # of radial stations
        XDWKLEN = 1.5
        ! default wake length in duct diameters
        XWAKE = -999.
        NWAKE = 0
        !
        TDUCT = 0.0
        TTOT = 0.0
        TVIS = 0.0
        QTOT = 0.0
        QVIS = 0.0
        PTOT = 0.0
        PVIS = 0.0
        !
        !
        DO IR = 1, IRX
            CHDES(IR) = 1.0
            CLDES(IR) = 0.8
            CLPOS(IR) = CLDES(IR)
            CLNEG(IR) = -CLDES(IR)
            IR2IEL(IR) = 0
            !
            DO N = 1, NRX
                BGAM(IR, N) = 0.0
                DO L = 1, 3
                    VIND(L, IR, N) = 0.0
                    VABS(L, IR, N) = 0.0
                    VREL(L, IR, N) = 0.0
                ENDDO
                CHR(IR, N) = CHDES(IR)
                CLR(IR, N) = CLDES(IR)
            ENDDO
        ENDDO
        !
        !---------------------------------------------------------------
        !--- Setup aerodynamic data for blade elements
        !--- Default aero properties
        !
        A0 = 0.       ! zero lift angle of attack   radians
        DCLDA = 6.0     ! lift curve slope            /radian
        CLMAX = 1.5    ! stall Cl
        CLMIN = -0.5    ! negative stall Cl
        DCL_STALL = 0.1
        ! CL increment from incipient to total stall
        DCLDA_STALL = 0.1
        ! stalled lift curve slope    /radian
        CMCON = -0.1   ! section Cm  (for pitch-axis moments)
        CDMIN = 0.013 ! minimum Cd
        CLDMIN = 0.5   ! Cl at minimum Cd
        DCDCL2 = 0.03  ! d(Cd)/d(Cl**2)
        DCDCL2S = 0.0   ! d(Cd)/d(Cl**2) secondary (annulus) drag (dead)
        REREF = 200000. ! Reynolds Number at which Cd values apply
        REXP = -0.4    ! Exponent for Re scaling of Cd:  Cd ~ Re**exponent
        MCRIT = 0.8     ! Critical Mach number
        TOC = 0.1     ! section thickness ratio (dead)
        !
        XPAXIS = 0.35   ! x/c location of pitch axis for blade plots, lofting
        !
        !--- Install default data into aero section #1 for each disk
        DO N = 1, NRX
            NAERO(N) = 1
            XISECT = 0.0
            CALL PUTAERO(N, NAERO(N), XISECT, A0, CLMAX, CLMIN, DCLDA, &
                    & DCLDA_STALL, DCL_STALL, CDMIN, CLDMIN, DCDCL2, DCDCL2S, &
                    & CMCON, MCRIT, TOC, REREF, REXP)
            !--- set pointers from blade sections to aero sections
            DO I = 1, IRX
                IAERO(I, N) = 1
            ENDDO
        ENDDO
        !
        !---- no BL solution yet
        LVISC = .FALSE.
        AMPMAX = 9.0
        UREN = 0.0
        DO IEL = 1, NEX
            CXVIS(IEL) = 0.0
            SFIX(1, IEL) = 1.0
            SFIX(2, IEL) = 1.0
            STRN(1, IEL) = 1.0
            STRN(2, IEL) = 1.0
            SSEP(1, IEL) = 1.0
            SSEP(2, IEL) = 1.0
            ICBLBEG(1, IEL) = 0
            ICBLBEG(2, IEL) = 0
            ICBLEND(1, IEL) = 0
            ICBLEND(2, IEL) = 0
        ENDDO
        !---- no inflow velocities defined
        NINFL = 0
        !
        !---------------------------------------------------------------
        !---- no grid dimensions yet
        II = 0
        JJ = 0
        !
        !---- by default input geometry is respaced
        LRSPC = .TRUE.
        !---------------------------------------------------------------
        !---- default number of panels per element, panel-number exponent
        NPANDEF = 61
        FPANDEF = 0.7
        !---- set panel parameters, no default spacing or user-defined spacing
        LRSPCUSR = .FALSE.
        LRSPCDEF = .FALSE.
        DO IEL = 1, NEX
            NPAN(IEL) = -999
            CVEX(IEL) = 1.0
            SMOF(IEL) = 1.0
            FSLE(IEL) = 0.6
            FSTE(IEL) = 0.6
            !------ refinement-station parameters
            DO IRPN = 1, NRPNX
                CRRAT(IRPN, IEL) = 1.0
                SRPN1(IRPN, IEL) = 0.0
                SRPN2(IRPN, IEL) = 0.0
            ENDDO
        ENDDO
        !
        !---- assume no case-paneling parameter file exists
        PFILE = ' '
        !
        !---------------------------------------------------------------
        !---- maximum TE gap (relative to perimeter) to deem airfoil as being closed
        DTEBOD = 0.1
        !---- minimum TE gap (relative to perimeter) to give closed airfoil a TE panel
        DTEPAN = 0.00001
        !
        !---------------------------------------------------------------
        !---- no geometry yet
        NEL = 0
        NBEL = 0
        NAME = ' '
        !
        !---- no accumulated point sequence yet
        NPOL = 0
        LPACC = .FALSE.
        !
        IPACT = 0
        !
        NPSEQ = 0
        !
        !---------------------------------------------------------------
        !---- Convergence flag for iterative solutions, no converged solution yet
        LCONV = .FALSE.
        !---- Iterative solver parameters
        RLXSOLV = 0.4
        EPSSOLV = 0.0002
        ITRMAXSOLV = 50
        !---- VMAVG is not initialized
        LVMAV = .FALSE.
        VAVGINIT = QINF + 1.0
        !
        !---- vortex wakes will not be automatically realigned
        LWRLX = .FALSE.
        !
        !---------------------------------------------------------------
        !---- don't normalize airfoil upon input
        LNORM = .FALSE.
        !
        !---- no geometry modes defined yet
        NMOD = 0
        LXYMOV = .FALSE.
        !
        !---------------------------------------------------------------
        !---- no valid control-point pointers yet
        LNCVP = .FALSE.
        !
        !---- no valid Aero,Geometry influence matrices yet
        LQAIC = .FALSE.
        LQGIC = .FALSE.
        LVAIC = .FALSE.
        LVGIC = .FALSE.
        !
        !---- no valid control-point velocities yet
        LQCNT = .FALSE.
        !
        !---- no valid system pointers yet
        LSYSP = .FALSE.
        !
        !---- no factored system matrix for gamma yet
        LGSYS = .FALSE.
        LQSYS = .FALSE.
        !
        !---- no gamma solution yet, no unit solutions
        LGAMA = .FALSE.
        LGAMU = .FALSE.
        !
        !---- no source-influence pointers or matrix
        LSIGP = .FALSE.
        LSIGM = .FALSE.
        !
        !---- don't plot reference Cp or force data
        LPREF = .FALSE.
        LFREF = .FALSE.
        !
        !---- default filename prefix
        PREFIX = ' '
        !
        CPFILE = ' '
        FRFILE = ' '
        !---------------------------------------------------------------
        !---- no target elements yet for QDES or GDES
        IELQDES = 0
        IELGDES = 1
        !
        !---- enforce geometry slope-matching at grafted-geometry endpoints
        LGSLOP = .TRUE.
        !
        !---- no geometric symmetry assumed
        LGSYMM = .FALSE.
        !
        !---- same geometry between panel and buffer
        LGSAME = .FALSE.
        !
        !---- plot home positions
        LHOMPL = .TRUE.
        !
        !---- no specified speed distributions yet
        NQSP = 0
        NSEG = 0
        LQSPEC = .FALSE.
        LQSPPL = .FALSE.
        LGSPPL = .FALSE.
        !
        !---- enforce slope-matching and gamma curvature at inverse-segment endpoints
        LQSLOP = .TRUE.
        LQCURV = .FALSE.
        !
        !---- inverse-segment endpoints on edge of a surface airfoil are free to move
        LGNFIX = .FALSE.
        !
        !---- do not display viscous Q(s)
        LQSVIS = .FALSE.
        !
        !---- do not reverse s or Q direction in Q(s) plots
        LQSREV = .FALSE.
        LQQREV = .FALSE.
        !
        !---- maximum underrelaxation factor
        RLXMAX = 1.0
        !
        !---------------------------------------------------------------
        !---- show point-strength and wake-strength values
        LPSHO = .TRUE.
        LWSHO = .TRUE.
        !
        !---- plot point-singularities, panel vertices, element numbers
        LPPLT = .TRUE.
        LVPLT = .TRUE.
        LEPLT = .TRUE.
        !
        !---- do not plot symbols on Cp or Q line plots
        LSPLT = .FALSE.
        !
        !---- plot displacement surfaces, don't use streamline integration to generate
        LDPLT = .TRUE.
        LDINT = .FALSE.
        !
        !---- plot grid on Cp or Q plots
        LPGRID = .TRUE.
        !
        !---- plot stagnation point on Cp or Q plots
        LPSTAG = .TRUE.
        !
        !---- plot grid on aero, airfoil geometry plots, ticks, parameters
        LAGRID = .TRUE.
        LGGRID = .TRUE.
        LGTICK = .TRUE.
        LGPARM = .TRUE.
        !
        LGEOPLX = .FALSE.
        ! for xfoil plot routines
        !
        !---- do not print neg. rpm Beta and Alfa output in local coordinates
        LCOORD = .FALSE.
        !
        !-----a case file has not been loaded
        LLOAD = .FALSE.
        !
        !
    END SUBROUTINE DFINIT
    !*==GENGEOM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! DFINIT
    
    
    SUBROUTINE GENGEOM
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        LOGICAL :: LQUERY
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
        IF (LRSPC) THEN
            !
            !---- If paneling parameters are unset then reset default paneling
            !     parameters for buffer geometry elements
            !      write(*,*) 'gengeom  npan',NPAN(1)
            !          WRITE(*,*) 'Setting default respacing parameters'
            IF (.NOT.LRSPCDEF) CALL PANDEF(0)
            !
            !        WRITE(*,*) 'Repaneling buffer geometry'
            LQUERY = .FALSE.
            CALL PANGEN(LQUERY)
            !
        ELSE
            !        WRITE(*,*) 'Copying buffer geometry into panel arrays'
            CALL PANCOP
        ENDIF
        !
        !----------------------------------------------------------------
        !---- At this point the input airfoils are in the panel geometry
        LGSAME = .TRUE.
        !
        !---- Adjust paneling on airfoils for rotor and wake grid
        CALL ADJPANL
        !      write(80,*) 'after ADJPANL PDMP'
        !      CALL XYPDMP(80)
        !
        !---- Reset # of elements to buffer geometry, remaining elements are added
        !     to this: drag elements, rotor line sources, then vortex wakes
        NEL = NBEL
        !
        !---- Set up drag area element(s)
        CALL SETDRGOBJ
        !
        !---- Set up rotor points
        CALL ROTORINIT
        !
        !---- Set up rotor wake grid from repaneled airfoils
        CALL INIGRD
        !
        !---- Set up rotor wake elements from wake grid
        CALL SETROTWAK
        !
        !---- Sets up control points and pointers to panels
        CALL CVPGEN
        !
        !---- Initialize rotor/wake pointers
        CALL ROTPINIT
        !
        !---- List all defined elements
        IF (LDBG) THEN
            CALL PELIST(6)
            CALL PELIST(LUNDBG)
        ENDIF
        CALL WAKEBOX
        !
        !----------------------------------------------------------------
        !---- invalidate any existing solution
        LNCVP = .FALSE.
        LQAIC = .FALSE.
        LQGIC = .FALSE.
        LQCNT = .FALSE.
        LSYSP = .FALSE.
        LGSYS = .FALSE.
        LGAMU = .FALSE.
        LGAMA = .FALSE.
        LSIGP = .FALSE.
        LSIGM = .FALSE.
        !
        NQSP = 0
        NSEG = 0
        !
        WRITE (*, *)
        WRITE (*, *) 'Geometry loaded'
        !
    END SUBROUTINE GENGEOM
    !*==PELIST.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GENGEOM
    
    
    SUBROUTINE PELIST(LU)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: LU
        !
        ! Local variables
        !
        INTEGER :: I, IEL, IP1, IP2, N
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Lists elements defined in panel geometry
        !------------------------------------------------------
        !
        !---- print out element info
        DO IEL = 1, NEL
            WRITE (LU, 1001)
            WRITE (LU, 1005) IEL
            !
            IP1 = IPFRST(IEL)
            IP2 = IPLAST(IEL)
            I = IP1
            N = IP2 - IP1 + 1
            !
            WRITE (LU, 1010) N
            !
            IF (NETYPE(IEL)==0) THEN
                IF (.NOT.(LBODY(IEL))) THEN
                    WRITE (LU, 2000) 'Open surface'
                ELSEIF (LXBOD(IEL)) THEN
                    WRITE (LU, 2000) 'Closed body on axis'
                ELSE
                    WRITE (LU, 2000) 'Closed body'
                ENDIF
                !
            ELSEIF (NETYPE(IEL)==1) THEN
                WRITE (LU, 2000) 'Prescribed wake surface'
                !
            ELSEIF (NETYPE(IEL)==2) THEN
                WRITE (LU, 2000) 'Prescribed line singularity on axis'
                !
            ELSEIF (NETYPE(IEL)==3) THEN
                WRITE (LU, 2000) 'Prescribed ring singularity'
                !
            ELSEIF (NETYPE(IEL)==4) THEN
                WRITE (LU, 2000) 'Prescribed point singularity'
                !
            ELSEIF (NETYPE(IEL)==5) THEN
                WRITE (LU, 2000) 'Source-only line singularity'
                !
            ELSEIF (NETYPE(IEL)==6) THEN
                WRITE (LU, 2000) 'Rotor source-only line singularity'
                !
            ELSEIF (NETYPE(IEL)==7) THEN
                WRITE (LU, 2000) 'Rotor vortex wake line'
                !
            ENDIF
            !
            IF (LTPAN(IEL)) WRITE (LU, 2000) 'TE-gap panel will be used'
            !
            IF (IEWAKE(IEL)/=0) WRITE (LU, 2100)                          &
                    & 'Prescribed-wake element', &
                    & IEWAKE(IEL), '  is attached'
            !
            IF (NETYPE(IEL)==0 .OR. NETYPE(IEL)==1 .OR. NETYPE(IEL)==2)  &
                    & THEN
                WRITE (LU, 1050) XPLE(IEL), YPLE(IEL), XPTE(IEL), &
                        & YPTE(IEL)
            ELSE
                WRITE (LU, 1060) 'First point at', XP(IP1), YP(IP1)
                WRITE (LU, 1060) 'Last  point at', XP(IP2), YP(IP2)
            ENDIF
            !
        ENDDO
        !
        WRITE (LU, 1001)
        !
        RETURN
        !...............................................................
        1001 FORMAT (/' -----------------------------------------------')
        1005 FORMAT (' Element', I3, ' ...')
        1010 FORMAT ('   Number of input coordinate points:', I4)
        1020 FORMAT ('     delta(x) =', F10.5, '   Scale =', F10.5, &
                &/'     delta(y) =', F10.5, '   Angle =', F10.5, ' deg')
        1050 FORMAT ('    LE x,y  =', 2F10.5/'    TE x,y  =', 2F10.5)
        1060 FORMAT (' ', A, ' x,y  =', 2F10.5)
        2000 FORMAT ('   ', A)
        2100 FORMAT ('   ', A, I3, A)
    END SUBROUTINE PELIST
    !*==GEPROC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PELIST
    
    
    
    SUBROUTINE GEPROC
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        REAL :: DB11SQ, DB12SQ, DB1TSQ, DB21SQ, DB22SQ, DB2TSQ, &
                & DSBMIN, DSBSQ, DSBTOL, DUMMY
        INTEGER :: I, IB, IB1, IB2, IEL, KB1, KB2, KEL, N
        LOGICAL :: LYZERO
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !     Processes buffer airfoil geometry to set
        !     element type, various limits and pointers.
        !----------------------------------------------------
        !
        !---- process all elements
        DO IEL = 1, NBEL
            CALL ELPROC(IEL)
            !
            I = IBFRST(IEL)
            N = IBLAST(IEL) - IBFRST(IEL) + 1
            CALL MINMAX(N, XB(I), XBMINE(IEL), XBMAXE(IEL))
            CALL MINMAX(N, YB(I), YBMINE(IEL), YBMAXE(IEL))
        ENDDO
        !
        !---- set overall min,max over all elements
        CALL MINMAX(NBEL, XBMINE, XBMIN, DUMMY)
        CALL MINMAX(NBEL, YBMINE, YBMIN, DUMMY)
        CALL MINMAX(NBEL, XBMAXE, DUMMY, XBMAX)
        CALL MINMAX(NBEL, YBMAXE, DUMMY, YBMAX)
        !
        DO IEL = 1, NBEL
            !------ default reference point for each element
            XBREFE(IEL) = 0.
            YBREFE(IEL) = 0.
            !------ clear element movement accumulators
            DXBSUM(IEL) = 0.
            DYBSUM(IEL) = 0.
            AGBSUM(IEL) = 0.
            XFBSUM(IEL) = 1.0
            YFBSUM(IEL) = 1.0
            !------ first assume this element doesn't have a wake
            IEWAKE(IEL) = 0
        ENDDO
        !
        !---- set element-type indicators and wake pointers
        DO IEL = 1, NBEL
            IB1 = IBFRST(IEL)
            IB2 = IBLAST(IEL)
            !
            !------ is this element entirely on the y=0 line?
            LYZERO = .TRUE.
            DO IB = IB1, IB2
                IF (YB(IB)/=0.0) LYZERO = .FALSE.
            ENDDO
            !
            IF (IB1==IB2) THEN
                !------- only one point defines element...
                IF (LYZERO) THEN
                    !-------- axisymm: point singularity on axis
                    NETYPE(IEL) = 4
                ELSE
                    !-------- axisymm: ring singularity
                    NETYPE(IEL) = 3
                ENDIF
                !
            ELSEIF (LYZERO) THEN
                !------- line element on axis
                NETYPE(IEL) = 2
                !
                !------- check for source line element defined in buffer geometry (old scheme)
            ELSEIF (NBTYPE(IEL)==5) THEN
                NETYPE(IEL) = NBTYPE(IEL)
                !
            ELSE
                !------- assume this element is a solid surface
                NETYPE(IEL) = 0
                !
                !------- Now check for wake elements
                !------- a wake element is attached to TE of some element before it
                DO KEL = 1, IEL - 1
                    !--------- check only solid boundarys (NETYPE=0) for attached wakes
                    IF (NETYPE(KEL)==0) THEN
                        KB1 = IBFRST(KEL)
                        KB2 = IBLAST(KEL)
                        DB11SQ = (XB(IB1) - XB(KB1))**2 + (YB(IB1) - YB(KB1))**2
                        DB21SQ = (XB(IB2) - XB(KB1))**2 + (YB(IB2) - YB(KB1))**2
                        DB12SQ = (XB(IB1) - XB(KB2))**2 + (YB(IB1) - YB(KB2))**2
                        DB22SQ = (XB(IB2) - XB(KB2))**2 + (YB(IB2) - YB(KB2))**2
                        DB1TSQ = (XB(IB1) - XBTE(KEL))**2 + (YB(IB1) - YBTE(KEL)) &
                                & **2
                        DB2TSQ = (XB(IB2) - XBTE(KEL))**2 + (YB(IB2) - YBTE(KEL)) &
                                & **2
                        !
                        DSBSQ = MIN(DB11SQ, DB12SQ, DB21SQ, DB22SQ, DB1TSQ, DB2TSQ)
                        DSBMIN = SQRT(DSBSQ)
                        !
                        !----------- this element IEL is KEL's wake if it starts at KEL's TE point
                        DSBTOL = 0.0001 * ABS(SB(KB2) - SB(KB1))
                        IF (DSBMIN<DSBTOL) THEN
                            NETYPE(IEL) = 1
                            IEWAKE(KEL) = IEL
                            EXIT
                        ENDIF
                    ENDIF
                ENDDO
                !
                !
            ENDIF
            !
        ENDDO
        !
        !---- set wake-termination box outline
        !c      CALL WAKEBOX
        !
    END SUBROUTINE GEPROC
    !*==GEPRINT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GEPROC
    
    
    SUBROUTINE GEPRINT(LU)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: LU
        !
        ! Local variables
        !
        INTEGER :: I, IB1, IB2, IEL, N
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !     Processes buffer airfoil geometry to set
        !     various limits and pointers.
        !----------------------------------------------------
        !
        !---- print out element info
        DO IEL = 1, NBEL
            WRITE (LU, 1001)
            WRITE (LU, 1005) IEL
            !
            IB1 = IBFRST(IEL)
            IB2 = IBLAST(IEL)
            I = IB1
            N = IB2 - IB1 + 1
            !
            WRITE (LU, 1010) N
            !
            IF (AREA2DA(IEL)>=0.0) THEN
                WRITE (LU, 2000) 'Counterclockwise ordering'
            ELSE
                WRITE (LU, 2000) 'Clockwise ordering'
            ENDIF
            !
            IF (NETYPE(IEL)==0) THEN
                IF (.NOT.(LBODY(IEL))) THEN
                    WRITE (LU, 2000) 'Open surface'
                ELSEIF (LXBOD(IEL)) THEN
                    WRITE (LU, 2000) 'Closed body on axis'
                ELSE
                    WRITE (LU, 2000) 'Closed body'
                ENDIF
                !
            ELSEIF (NETYPE(IEL)==1) THEN
                WRITE (LU, 2000) 'Prescribed wake surface'
                !
            ELSEIF (NETYPE(IEL)==2) THEN
                WRITE (LU, 2000) 'Prescribed line singularity on axis'
                !
            ELSEIF (NETYPE(IEL)==3) THEN
                WRITE (LU, 2000) 'Prescribed ring singularity'
                !
            ELSEIF (NETYPE(IEL)==4) THEN
                WRITE (LU, 2000) 'Prescribed point singularity'
                !
            ENDIF
            !
            IF (LTPAN(IEL)) WRITE (LU, 2000) 'TE-gap panel will be used'
            !
            IF (IEWAKE(IEL)/=0) WRITE (LU, 2100)                          &
                    & 'Prescribed-wake element', &
                    & IEWAKE(IEL), '  is attached'
            !
            IF (NETYPE(IEL)==0 .OR. NETYPE(IEL)==1 .OR. NETYPE(IEL)==2)  &
                    & THEN
                WRITE (LU, 1050) XBLE(IEL), YBLE(IEL), XBTE(IEL), &
                        & YBTE(IEL)
            ELSE
                WRITE (LU, 1060) XB(IB1), YB(IB1)
            ENDIF
            !
        ENDDO
        !
        WRITE (LU, 1001)
        !
        RETURN
        !
        !...............................................................
        1001 FORMAT (/' -----------------------------------------------')
        1005 FORMAT (' Element', I3, ' ...')
        1010 FORMAT ('   Number of input coordinate points:', I4)
        1020 FORMAT ('     delta(x) =', F10.5, '   Scale =', F10.5, &
                &/'     delta(y) =', F10.5, '   Angle =', F10.5, ' deg')
        1050 FORMAT ('    LE x,y  =', 2F10.5/'    TE x,y  =', 2F10.5)
        1060 FORMAT ('       x,y  =', 2F10.5)
        2000 FORMAT ('   ', A)
        2100 FORMAT ('   ', A, I3, A)
    END SUBROUTINE GEPRINT
    !*==ELPROC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GEPRINT
    
    
    
    SUBROUTINE ELPROC(IEL)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IEL
        !
        ! Local variables
        !
        REAL :: DSBTE, PERIM, SBLEN, SS, XTMP, YTMP
        INTEGER :: I, IB, IB1, IB2, IBCUT, N
        LOGICAL :: LCUT, LTMP
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Processes buffer geometry element IEL.
        !------------------------------------------------------
        !
        IB1 = IBFRST(IEL)
        IB2 = IBLAST(IEL)
        !
        I = IB1
        N = IB2 - IB1 + 1
        CALL SCALC(XB(I), YB(I), SB(I), N)
        !
        !----- check for body that crosses the axis (defined by upper and lower sides)
        LCUT = .FALSE.
        DO IB = IB1, IB2 - 1
            IF (YB(IB)>0.0 .AND. YB(IB + 1)<=0.0) THEN
                LCUT = .TRUE.
                IBCUT = IB + 1
                EXIT
            ENDIF
        ENDDO
        !
        IF (LCUT) THEN
            !----- element crosses axis, cut it at the axis and use only upper half
            IF (LDBG) THEN
                WRITE (*, *) 'Cutting element ', IEL, ' at X axis '
                WRITE (LUNDBG, *) 'Cutting element ', IEL, ' at X axis '
            ENDIF
            CALL SEGSPL(XB(I), XBS(I), SB(I), N)
            CALL SEGSPL(YB(I), YBS(I), SB(I), N)
            SS = SB(IBCUT)
            CALL SINVRT(SS, 0.0, YB(IB1), YBS(IB1), SB(IB1), N)
            XB(IBCUT) = SEVAL(SS, XB(IB1), XBS(IB1), SB(IB1), N)
            YB(IBCUT) = 0.0
            IBLAST(IEL) = IBCUT
            !
            !----- recalculate indices for element and arc length
            IB1 = IBFRST(IEL)
            IB2 = IBLAST(IEL)
            I = IB1
            N = IB2 - IB1 + 1
            CALL SCALC(XB(I), YB(I), SB(I), N)
            !c       DO IB = IB1, IB2
            !c         write(50,*) IB, XB(IB), YB(IB)
            !c       END DO
        ENDIF
        !
        !
        IF (N==1) THEN
            !---- single point element (doublet/source point)
            LV1ZR(IEL) = .FALSE.
            LV2ZR(IEL) = .FALSE.
            LXBOD(IEL) = .FALSE.
            LBODY(IEL) = .FALSE.
            LTPAN(IEL) = .FALSE.
            LQNZR(IEL) = .FALSE.
            !
        ELSE
            !---- multiple point element (doublet/source line or vortex/source panels)
            SBLEN = ABS(SB(IB2) - SB(IB1))
            !
            !---- are first,last points effectively on the axis?
            LV1ZR(IEL) = ABS(YB(IB1))<SBLEN * DTEPAN
            LV2ZR(IEL) = ABS(YB(IB2))<SBLEN * DTEPAN
            !
            IF (LV1ZR(IEL)) THEN
                !---- TE gap is really a TE disk at point IB2
                DSBTE = ABS(YB(IB2))
            ELSE
                !---- usual TE gap
                DSBTE = SQRT((XB(IB1) - XB(IB2))**2 + (YB(IB1) - YB(IB2))**2)
            ENDIF
            !
            !---- axisymmetric body closed on the axis ?
            LXBOD(IEL) = LV1ZR(IEL) .AND. LV2ZR(IEL)
            !
            !---- finite-thickness body ?
            LBODY(IEL) = LV1ZR(IEL) .OR. LV2ZR(IEL) .OR.                   &
                    & DSBTE<ABS(SB(IB2) - SB(IB1)) * DTEBOD
            !
            !---- add TE-gap panel ?
            LTPAN(IEL) = LBODY(IEL) .AND. DSBTE>ABS(SB(IB2) - SB(IB1)) * DTEPAN
            !
            !---- explicitly zero out the body normal velocity QNDOF?
            !c     LQNZR(IEL) = .TRUE.
            LQNZR(IEL) = .FALSE.
            !
        ENDIF
        !
        !---- 2D buffer geometry data
        !---- set 2D area, centroid, section moduli - for area
        CALL AECALC(N, XB(I), YB(I), ONE, 0, PERIM, AREA2DA(IEL), XBCEN2DA(IEL), &
                & YBCEN2DA(IEL), EIXX2DA(IEL), EIYY2DA(IEL), EIXY2DA(IEL))
        !---- set 2D area, centroid, section moduli - per unit skin thickness
        CALL AECALC(N, XB(I), YB(I), ONE, 2, PERIM, AREA2DT(IEL), XBCEN2DT(IEL), &
                & YBCEN2DT(IEL), EIXX2DT(IEL), EIYY2DT(IEL), EIXY2DT(IEL))
        !
        !---- Axisymmetric buffer geometry data
        !---- set axisymmetric surface area, volume, centroid, radii of gyration
        CALL AXCALC(N, XB(I), YB(I), ONE, 0, VOLUMV(IEL), ASURFV(IEL), &
                & XBCENV(IEL), YBCENV(IEL), RGYRXV(IEL), RGYRYV(IEL))
        !
        !---- assume element will be used and written in input order by default
        LREVEL(IEL) = .FALSE.
        !
        !----- area is negative (clockwise order)... mark for reversing the order
        IF (LBODY(IEL) .AND. AREA2DA(IEL)<0.0) LREVEL(IEL) = .TRUE.
        !
        IF (LREVEL(IEL)) THEN
            !---- reverse the order of the points
            DO IB = IB1, (IB2 + IB1) / 2
                XTMP = XB(IB2 - IB + IB1)
                YTMP = YB(IB2 - IB + IB1)
                XB(IB2 - IB + IB1) = XB(IB)
                YB(IB2 - IB + IB1) = YB(IB)
                XB(IB) = XTMP
                YB(IB) = YTMP
            ENDDO
            LTMP = LV1ZR(IEL)
            LV1ZR(IEL) = LV2ZR(IEL)
            LV2ZR(IEL) = LTMP
        ENDIF
        !
        CALL SCALC(XB(I), YB(I), SB(I), N)
        CALL SEGSPL(XB(I), XBS(I), SB(I), N)
        CALL SEGSPL(YB(I), YBS(I), SB(I), N)
        !
        IF (LBODY(IEL) .AND. .NOT.LXBOD(IEL)) THEN
            !----- body off axis
            IF (LV2ZR(IEL)) THEN
                SBLE(IEL) = SB(IB2)
                XBLE(IEL) = XB(IB2)
                YBLE(IEL) = YB(IB2)
                XBTE(IEL) = XB(IB1)
                YBTE(IEL) = YB(IB1)
            ELSE
                CALL LEFIND(SBLE(IEL), XB(I), XBS(I), YB(I), YBS(I), SB(I), N)
                XBLE(IEL) = SEVAL(SBLE(IEL), XB(I), XBS(I), SB(I), N)
                YBLE(IEL) = SEVAL(SBLE(IEL), YB(I), YBS(I), SB(I), N)
                XBTE(IEL) = 0.5 * (XB(IB1) + XB(IB2))
                YBTE(IEL) = 0.5 * (YB(IB1) + YB(IB2))
            ENDIF
            !
            !----- other body... set LE,TE from leftmost and rightmost points
        ELSEIF (XB(IB1)<XB(IB2)) THEN
            SBLE(IEL) = SB(IB1)
            XBLE(IEL) = XB(IB1)
            YBLE(IEL) = YB(IB1)
            XBTE(IEL) = XB(IB2)
            YBTE(IEL) = YB(IB2)
        ELSE
            SBLE(IEL) = SB(IB2)
            XBLE(IEL) = XB(IB2)
            YBLE(IEL) = YB(IB2)
            XBTE(IEL) = XB(IB1)
            YBTE(IEL) = YB(IB1)
            !
        ENDIF
        !
        !---- set Cp-side flag
        IF (LBODY(IEL)) THEN
            !----- outer surface only (to right side of +s CCW direction)
            ISPLOT(IEL) = +1
        ELSE
            !----- default... plot both sides
            ISPLOT(IEL) = 0
        ENDIF
        !
    END SUBROUTINE ELPROC
    !*==WAKEBOX.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ELPROC
    
    
    
    SUBROUTINE WAKEBOX
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        REAL :: DXBOX, DYBOX
        INTEGER :: IEL
        LOGICAL :: LWBSET
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------
        !     Sets limits of "wake box".
        !     A wake is terminated when it leaves this box.
        !---------------------------------------------------
        !
        XWBOX(1) = 1.0E23
        XWBOX(2) = -1.0E23
        YWBOX(1) = 1.0E23
        YWBOX(2) = -1.0E23
        LWBSET = .FALSE.
        DO IEL = 1, NEL
            IF (NETYPE(IEL)==0 .OR. NETYPE(IEL)==7) THEN
                !------- box defined only by solid or vortex wake elements
                XWBOX(1) = MIN(XWBOX(1), XPMINE(IEL))
                XWBOX(2) = MAX(XWBOX(2), XPMAXE(IEL))
                YWBOX(1) = MIN(YWBOX(1), YPMINE(IEL))
                YWBOX(2) = MAX(YWBOX(2), YPMAXE(IEL))
                LWBSET = .TRUE.
            ENDIF
        ENDDO
        IF (.NOT.LWBSET) THEN
            !----- use default if no solid elements exist
            XWBOX(1) = XBMIN
            XWBOX(2) = XBMAX
            YWBOX(1) = YBMIN
            YWBOX(2) = YBMAX
        ENDIF
        !
        !c      DBOX = MAX( XWBOX(2)-XWBOX(1) , YWBOX(2)-YWBOX(1) )
        !c      XWBOX(1) = XWBOX(1) - DBOX
        !c      XWBOX(2) = XWBOX(2) + DBOX
        !c      YWBOX(1) = YWBOX(1) - DBOX
        !c      YWBOX(2) = YWBOX(2) + DBOX
        !
        DXBOX = XWBOX(2) - XWBOX(1)
        DYBOX = YWBOX(2) - YWBOX(1)
        XWBOX(1) = XWBOX(1) - 0.5 * DXBOX
        XWBOX(2) = XWBOX(2) - 0.02 * DXBOX
        YWBOX(1) = 0.0
        YWBOX(2) = YWBOX(2) + DYBOX
        !
    END SUBROUTINE WAKEBOX
    !*==SETDRGOBJ.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! WAKEBOX
    
    
    SUBROUTINE SETDRGOBJ
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: ND
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Sets up drag area elements into paneled geometry
        !------------------------------------------------------
        !
        IF (NDOBJ<1) RETURN
        !
        NPTOT = IPLAST(NEL)
        !
        !---- Use rotor points to define paneled element for rotor source line
        DO ND = 1, NDOBJ
            CALL ADDELEM5(XDDEF(1, ND), YDDEF(1, ND), NDDEF(ND))
            IELDRGOBJ(ND) = NEL
        ENDDO
        LDRGOBJ = .TRUE.
        !
    END SUBROUTINE SETDRGOBJ
    !*==ROTORINIT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE ROTORINIT
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        REAL, DIMENSION(IRX) :: BETAR0, BGAM0, CHR0, CLDES0, T1, &
                & T1S, T2, T2S, T3, T3S, T4, T4S, &
                & YRC0
        REAL :: CLD1, CLD2, CLN1, CLN2, CLP1, CLP2, DR, DX, DY, &
                & TG, TGAPZLT, XCB, XDW, YCB, YDW
        INTEGER :: I, IC, IPCB, IPDW, IR, N, NRC0, NRCSAV
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Sets up rotor geometry from panel geometry
        !------------------------------------------------------
        !
        IF (NROTOR==0) THEN
            WRITE (*, *) 'No actuator disk/blades defined'
            RETURN
        ENDIF
        !
        NRCSAV = NRC
        IF (NRCSAV>0) THEN
            N = 1  ! assume at least one disk defined
            DO IR = 1, NRCSAV
                YRC0(IR) = YRC(IR, N)
                CLDES0(IR) = CLDES(IR)
            ENDDO
        ENDIF
        !
        DO N = 1, NROTOR
            !
            !---- get radial location of rotor on CB and duct wall
            IF (IPROTCB(N)==0) THEN
                WRITE (*, *) 'Error locating disk on CB wall'
                STOP
            ENDIF
            IPCB = IPROTCB(N)
            XCB = XP(IPCB)
            YCB = YP(IPCB)
            RHUB(N) = YCB
            IF (LDBG) THEN
                WRITE (LUNDBG, *) 'Disk on CB @IP= ', IPCB
                WRITE (LUNDBG, *) 'Xhub,Rhub on CB ', XCB, RHUB(N)
            ENDIF
            !
            IF (IPROTDW(N)==0) THEN
                WRITE (*, *) 'Error locating rotor on Duct wall'
                STOP
            ENDIF
            IPDW = IPROTDW(N)
            XDW = XP(IPDW)
            YDW = YP(IPDW)
            RTIP(N) = YDW
            IF (LDBG) THEN
                WRITE (LUNDBG, *) 'Disk on Duct @IP= ', IPDW
                WRITE (LUNDBG, *) 'Xtip,Rtip on Duct wall ', XDW, RTIP(N)
            ENDIF
            !
            !---- Update swept disk area
            ADISK(N) = PI * (RTIP(N)**2 - RHUB(N)**2)
            ATIP(N) = PI * (RTIP(N)**2)
            IF (LDBG) WRITE (*, *) 'Disk # ', N, ' swept area ', &
                    & ADISK(N)
            !
            !
            IF (IRTYPE(N)/=0) THEN ! disk defined
                !---- Save old blade and circulation distribution for interpolation to new pts
                NRC0 = NRCSAV
                DO IR = 1, NRC0
                    YRC0(IR) = YRC(IR, N)
                    BGAM0(IR) = BGAM(IR, N)
                    CHR0(IR) = CHR(IR, N)
                    BETAR0(IR) = BETAR(IR, N)
                    !          IF(N.EQ.1) CLDES0(IR) = CLDES(IR)
                ENDDO
            ENDIF
            !
            !---- Set discrete points on rotor line using linear interpolation hub-tip
            NRP = NRSTA
            DX = (XDW - XCB) / FLOAT(NRP - 1)
            DY = (YDW - YCB) / FLOAT(NRP - 1)
            DO IR = 1, NRP
                XRP(IR, N) = XCB + FLOAT(IR - 1) * DX
                YRP(IR, N) = YCB + FLOAT(IR - 1) * DY
                !c        write(*,*) 'rotorpts ',i,xrp(i,n),yrp(i,n)
            ENDDO
            !
            !
            ! *** In what follows, condition on LLOFT is to allow paneling to geometric
            !     tip gap for lofting, while changing nothing for OPER.
            !
            !---- If effective tip gap is specified add two points to last interval
            !     at tip gap from duct wall, shift original rotor points inwards
            !
            IF (LLOFT) THEN
                TGAPZLT = 0.0
            ELSE
                TGAPZLT = TGAPZL
            ENDIF
            !
            TG = MAX(0.0, TGAP - TGAPZLT)
            IF (TG>0.0) THEN
                DX = XDW - XCB
                DR = RTIP(N) - RHUB(N)
                DO IR = 1, NRP
                    YRP(IR, N) = YRP(IR, N) - 2.0 * TG * (YRP(IR, N) - RHUB(N)) / DR
                    XRP(IR, N) = XRP(IR, N) - DX * (RTIP(N) - YRP(IR, N)) / DR
                    !c          write(*,*) 'rotorpts ',i,xrp(i,n),yrp(i,n)
                ENDDO
                !--- Add two points, one at RTIP-TGAP, one at shroud wall at RTIP
                NRP = NRP + 1
                YRP(NRP, N) = RTIP(N) - TG
                XRP(NRP, N) = XDW - DX * (RTIP(N) - YRP(NRP, N)) / DR
                NRP = NRP + 1
                XRP(NRP, N) = XDW
                YRP(NRP, N) = YDW
            ENDIF
            !
            !---- Rotor interval centers and initial circulation
            NRC = NRP - 1
            DO IR = 1, NRC
                XRC(IR, N) = 0.5 * (XRP(IR, N) + XRP(IR + 1, N))
                YRC(IR, N) = 0.5 * (YRP(IR, N) + YRP(IR + 1, N))
                BGAM(IR, N) = 0.0
                !c        write(*,*) 'rotorcpts ',ir,xrc(ir,n),yrc(ir,n)
            ENDDO
            !
            !=====================================================================
            !
            IF (IRTYPE(N)==0) THEN ! disk undefined
                !
                !---- If actuator disk has been input interpolate BGAM to rotor center points
                IF (IRTYPDEF(N)==1 .AND. NRDEF(N)>=2) THEN
                    DO I = 1, NRDEF(N)
                        T1(I) = BGAMDEF(I, N)
                    ENDDO
                    !---- Spline rotor definition arrays
                    CALL SEGSPL(T1, T1S, YRDEF(1, N), NRDEF(N))
                    !---- Set BGAM for points on rotor line
                    DO IR = 1, NRC
                        BGAM(IR, N) = SEVAL(YRC(IR, N), T1, T1S, YRDEF(1, N), &
                                & NRDEF(N))
                    ENDDO
                    IF (TGAP>0.0) BGAM(NRC, N) = 0.0
                    IRTYPE(N) = 1
                    LBLDEF = .FALSE.
                    LVMAV = .FALSE.
                ENDIF
                !---- If a rotor has been input (Y,CH,BETA) interpolate CH,BETA to rotor
                !     center points
                IF (IRTYPDEF(N)==2 .AND. NRDEF(N)>=2) THEN
                    DO I = 1, NRDEF(N)
                        T1(I) = CHRDEF(I, N)
                        T2(I) = BETADEF(I, N)
                    ENDDO
                    !---- Spline rotor definition arrays
                    CALL SEGSPL(T1, T1S, YRDEF(1, N), NRDEF(N))
                    CALL SEGSPL(T2, T2S, YRDEF(1, N), NRDEF(N))
                    !---- Set CH,BETA for points on rotor line
                    DO IR = 1, NRC
                        CHR(IR, N) = SEVAL(YRC(IR, N), T1, T1S, YRDEF(1, N), NRDEF(N)&
                                &)
                        BETAR(IR, N) = SEVAL(YRC(IR, N), T2, T2S, YRDEF(1, N), &
                                & NRDEF(N))
                    ENDDO
                    IRTYPE(N) = 2
                    LBLDEF = .TRUE.
                    LVMAV = .FALSE.
                    CALL SETIAERO
                ENDIF
                !
            ELSE
                !
                !---- Interpolate old rotor/circulation to new rotor points
                DO I = 1, NRC0
                    T1(I) = CHR0(I)
                    T2(I) = BETAR0(I)
                    T3(I) = BGAM0(I)
                ENDDO
                !---- Spline rotor definition arrays
                CALL SEGSPL(T1, T1S, YRC0, NRC0)
                CALL SEGSPL(T2, T2S, YRC0, NRC0)
                CALL SEGSPL(T3, T3S, YRC0, NRC0)
                !---- Set CH,BETA,BGAM for points on new rotor line
                DO IR = 1, NRC
                    CHR(IR, N) = SEVAL(YRC(IR, N), T1, T1S, YRC0, NRC0)
                    BETAR(IR, N) = SEVAL(YRC(IR, N), T2, T2S, YRC0, NRC0)
                    BGAM(IR, N) = SEVAL(YRC(IR, N), T3, T3S, YRC0, NRC0)
                ENDDO
                !cc        IF(TGAP.GT.0.0) BGAM(NRC,N) = 0.0
                CALL SETIAERO
                !
            ENDIF
            !
        ENDDO
        ! end NROTOR loop
        !
        !---- Move design CLs and BB data to new radial stations
        !
        IF (NRCSAV>0) THEN
            CLP1 = CLPOS(1)
            CLP2 = CLPOS(NRCSAV)
            CLN1 = CLNEG(1)
            CLN2 = CLNEG(NRCSAV)
            CLD1 = CLDES(1)
            CLD2 = CLDES(NRCSAV)
            !
            DO IR = 1, NRC
                CLPOS(IR) = CLP1 + FLOAT(IR - 1) / FLOAT(NRC - 1) * (CLP2 - CLP1)
                CLNEG(IR) = CLN1 + FLOAT(IR - 1) / FLOAT(NRC - 1) * (CLN2 - CLN1)
                CLDES(IR) = CLD1 + FLOAT(IR - 1) / FLOAT(NRC - 1) * (CLD2 - CLD1)
            ENDDO
            !
            DO N = 1, NROTOR
                IF (LBBLOFT(N)) THEN
                    DO IC = 1, NRCSAV
                        T4(IC) = BBVFAC(IC, N)
                    ENDDO
                    CALL SEGSPL(T4, T4S, YRC0, NRCSAV)
                    DO IC = 1, NRC
                        BBVFAC(IC, N) = SEVAL(YRC(IC, N), T4, T4S, YRC0, NRCSAV)
                    ENDDO
                ENDIF
            ENDDO
        ENDIF
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
    END SUBROUTINE ROTORINIT
    !*==SETROTWAK.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ROTORINIT
    
    
    
    SUBROUTINE SETROTWAK
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: IEL, IG, IG1, IP, IP1, IP2, IR, N, NR
        REAL, DIMENSION(IRX) :: T1, T1S, T2, T2S, T3, T3S, Y1
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
        DO NR = 1, NROTOR
            !
            !---- Define paneled source line element for rotor lines from grid points
            IG = IGROTOR(NR)
            DO IR = 1, NRP
                XRP(IR, NR) = XG(IG, IR)
                YRP(IR, NR) = YG(IG, IR)
            ENDDO
            CALL ADDELEM6(XRP(1, NR), YRP(1, NR), NRP)
            IELROTOR(NR) = NEL
            !
            !---- Interpolate old rotor params/circulation to new rotor points
            NRC = NRP - 1
            DO IR = 1, NRC
                Y1(IR) = YRC(IR, NR)
                T1(IR) = CHR(IR, NR)
                T2(IR) = BETAR(IR, NR)
                T3(IR) = BGAM(IR, NR)
            ENDDO
            !---- Redefine rotor center points
            DO IR = 1, NRC
                XRC(IR, NR) = 0.5 * (XRP(IR, NR) + XRP(IR + 1, NR))
                YRC(IR, NR) = 0.5 * (YRP(IR, NR) + YRP(IR + 1, NR))
            ENDDO
            !---- Spline rotor definition arrays
            CALL SEGSPL(T1, T1S, Y1, NRC)
            CALL SEGSPL(T2, T2S, Y1, NRC)
            CALL SEGSPL(T3, T3S, Y1, NRC)
            !---- Set CH,BETA,BGAM for points on new rotor line
            DO IR = 1, NRC
                CHR(IR, NR) = SEVAL(YRC(IR, NR), T1, T1S, Y1, NRC)
                BETAR(IR, NR) = SEVAL(YRC(IR, NR), T2, T2S, Y1, NRC)
                BGAM(IR, NR) = SEVAL(YRC(IR, NR), T3, T3S, Y1, NRC)
            ENDDO
            !cc        IF(TGAP.GT.0.0) BGAM(NRC,NR) = 0.0
            CALL SETIAERO
            !
        ENDDO
        !
        !---- Add vortex wake elements to panel geometry
        !
        !----- set up wake element from centerbody TE (element 1)
        IR = 1
        IEL = 1
        IG1 = IGTECB
        N = II - IG1 + 1
        CALL ADDWAKE(XG(IG1, IR), YG(IG1, IR), N)
        !---- save wake element index for this rotor point
        IR2IEL(IR) = NEL
        IEL2IR(NEL) = IR
        !---- plot Cp only on interior flow side
        ISPLOT(NEL) = 2
        IEWAKE(IEL) = NEL
        !---- set up pointers from wake points to grid IG streamwise index
        IP1 = IPFRST(NEL)
        IP2 = IPLAST(NEL)
        DO IP = IP1, IP2
            IP2IG(IP) = IG1 + IP - IP1
        ENDDO
        !
        !----- set up intermediate rotor wakes (from rotor line)
        DO IR = 2, NRP - 1
            IG1 = 1
            N = II
            CALL ADDWAKE(XG(IG1, IR), YG(IG1, IR), N)
            !---- save wake element index for this rotor point
            IR2IEL(IR) = NEL
            IEL2IR(NEL) = IR
            !---- plot Cp only on interior flow side
            ISPLOT(NEL) = -1
            !---- set up pointers from wake points to grid IG streamwise index
            IP1 = IPFRST(NEL)
            IP2 = IPLAST(NEL)
            DO IP = IP1, IP2
                IP2IG(IP) = IG1 + IP - IP1
            ENDDO
        ENDDO
    
        !----- and wake from duct TE (element 2)
        IR = NRP
        IEL = 2
        IG1 = IGTEDW
        N = II - IG1 + 1
        CALL ADDWAKE(XG(IG1, IR), YG(IG1, IR), N)
        !---- save wake element index for this rotor point
        IR2IEL(IR) = NEL
        IEL2IR(NEL) = IR
        !---- plot Cp only on interior flow side
        ISPLOT(NEL) = 1
        IEWAKE(IEL) = NEL
        !---- set up pointers from wake points to grid IG streamwise index
        IP1 = IPFRST(NEL)
        IP2 = IPLAST(NEL)
        DO IP = IP1, IP2
            IP2IG(IP) = IG1 + IP - IP1
        ENDDO
        !
        !----------------------------------------------------------------
        !---- invalidate any existing solution
        LNCVP = .FALSE.
        LQAIC = .FALSE.
        LQGIC = .FALSE.
        LQCNT = .FALSE.
        LSYSP = .FALSE.
        LGSYS = .FALSE.
        LGAMU = .FALSE.
        LGAMA = .FALSE.
        LSIGP = .FALSE.
        LSIGM = .FALSE.
        !
    END SUBROUTINE SETROTWAK
    !*==ADDWAKE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETROTWAK
    
    
    SUBROUTINE ADDWAKE(X, Y, N)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(*) :: X, Y
        !
        ! Local variables
        !
        INTEGER :: I, IEL, IP, IP1, IP2
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Add wake specified by N points X,Y to panel geometry
        !---------------------------------------------------------
        !
        !---- add wake to panel geometry
        IP1 = IPLAST(NEL) + 1
        IP2 = IP1 - 1 + N
        IEL = NEL + 1
        IF (LDBG) WRITE (LUNDBG, 10) IEL, N, IP1, IP2
        !---- add wake points
        DO IP = IP1, IP2
            I = IP - IP1 + 1
            XP(IP) = X(I)
            YP(IP) = Y(I)
            IF (LDBG) WRITE (LUNDBG, *) IP, XP(IP), YP(IP)
        ENDDO
        !---- centroid location
        XPCENT(IEL) = 0.5 * (XP(IP1) + XP(IP2))
        YPCENT(IEL) = 0.5 * (YP(IP1) + YP(IP2))
        !---- panel limits
        CALL MINMAX(N, XP(IP1), XPMINE(IEL), XPMAXE(IEL))
        CALL MINMAX(N, YP(IP1), YPMINE(IEL), YPMAXE(IEL))
        !
        !---- setup pointers and panel data
        IPFRST(IEL) = IP1
        IPLAST(IEL) = IP2
        NETYPE(IEL) = 7
        LTPAN(IEL) = .TRUE.
        !      LTPAN(IEL) = .FALSE.
        IPTE2(IEL) = IP2
        IPTE1(IEL) = 0
        CALL XYPSPL(IEL)
        !      CALL XYCSET(IEL)
        !
        NEL = IEL
        NPTOT = IP2
        !
        IF (LDBG) WRITE (*, 10) IEL, N, IP1, IP2
        !
        10   FORMAT ('Vortex wake element ', I4, ' added with #pts ', I5, ' IP1 ', &
                & I5, ' IP2 ', I5)
        !
    END SUBROUTINE ADDWAKE
    !*==UPDROTWAK.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ADDWAKE
    
    
    
    SUBROUTINE UPDROTWAK
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: I, I1, IEL, IP, IP1, IP2, IR, N, NIP
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
        IR = 1
        IEL = IR2IEL(IR)
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        I1 = IGTECB
        N = II - I1 + 1
        NIP = IP2 - IP1 + 1
        IF (N/=NIP) THEN
            WRITE (*, *) 'UPDROTWAK wake grid mismatch on IR ', IR, N, &
                    & NIP
            STOP
        ENDIF
        I = I1
        DO IP = IP1, IP2
            XP(IP) = XG(I, IR)
            YP(IP) = YG(I, IR)
            I = I + 1
        ENDDO
        !
        !----- set up intermediate rotor wakes (from rotor line)
        DO IR = 2, NRP - 1
            IEL = IR2IEL(IR)
            IP1 = IPFRST(IEL)
            IP2 = IPLAST(IEL)
            I1 = 1
            N = II
            NIP = IP2 - IP1 + 1
            IF (N/=NIP) THEN
                WRITE (*, *) 'UPDROTWAK wake grid mismatch on IR ', IR, N, &
                        & NIP
                STOP
            ENDIF
            I = I1
            DO IP = IP1, IP2
                XP(IP) = XG(I, IR)
                YP(IP) = YG(I, IR)
                I = I + 1
            ENDDO
        ENDDO
    
        !----- and wake from duct TE (element 2)
        IR = NRP
        IEL = IR2IEL(IR)
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        I1 = IGTEDW
        N = II - I1 + 1
        NIP = IP2 - IP1 + 1
        IF (N/=NIP) THEN
            WRITE (*, *) 'UPDROTWAK wake grid mismatch on IR ', IR, N, &
                    & NIP
            STOP
        ENDIF
        I = I1
        DO IP = IP1, IP2
            XP(IP) = XG(I, IR)
            YP(IP) = YG(I, IR)
            I = I + 1
        ENDDO
        !
        !----------------------------------------------------------------
        !---- invalidate any existing solution
        LNCVP = .FALSE.
        LQAIC = .FALSE.
        LQGIC = .FALSE.
        LQCNT = .FALSE.
        LGSYS = .FALSE.
        LGAMU = .FALSE.
        LGAMA = .FALSE.
        !
    END SUBROUTINE UPDROTWAK
    !*==ROTPINIT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! UPDROTWAK
    
    
    SUBROUTINE ROTPINIT
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: IC, IC1, IC2, IEL, IG1, IG2, IP, IP1, IP2, &
                & IPCB, IPDW, IR, N, NR
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
        DO IP = 1, IPX
            IP2IR(IP) = 0
        ENDDO
        !
        !---- Clear pointers to grid centers from body and wake centers
        DO IC = 1, ICX
            IC2IG(IC) = 0
        ENDDO
        !
        NR = 1
        !---- Pointers for vortex wake on centerbody
        IR = 1
        IEL = 1
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        IPCB = IPROTCB(NR)
        !c    IF(LDBG) WRITE(22,*) IR, IEL, IP1, IP2, IPCB
        DO IP = IPCB, IP1, -1
            IP2IR(IP) = IR
            !c      IF(LDBG) WRITE(22,*) IP,IR
        ENDDO
        !---- set pointers from panel center to grid center
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        DO IC = IC1, IC2
            IG1 = IP2IG(IPCO(IC))
            IG2 = IP2IG(IPCP(IC))
            IC2IG(IC) = MIN(IG1, IG2)
            !c        IF(IC2IG(IC).NE.0) WRITE(22,*) IC,IC2IG(IC)
        ENDDO
        !
        !---- Set pointers for vortex wake elements
        DO IR = 1, NRP
            IEL = IR2IEL(IR)
            IF (IEL<=0) THEN
                WRITE (*, *) 'IR2IEL out of range in ROTPINIT ', IEL
                STOP
            ENDIF
            IP1 = IPFRST(IEL)
            IP2 = IPLAST(IEL)
            DO IP = IP1, IP2
                IP2IR(IP) = IR
            ENDDO
            !---- set pointers from panel center to grid center
            IC1 = ICFRST(IEL)
            IC2 = ICLAST(IEL)
            DO IC = IC1, IC2
                IG1 = IP2IG(IPCO(IC))
                IG2 = IP2IG(IPCP(IC))
                IC2IG(IC) = MIN(IG1, IG2)
                !c          IF(IC2IG(IC).NE.0) WRITE(22,*) IC,IC2IG(IC)
            ENDDO
        ENDDO
        !
        !---- Pointers for vortex wake on duct
        IR = NRP
        IEL = 2
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        N = IP2 - IP1 + 1
        IPDW = IPROTDW(NR)
        !c    IF(LDBG) WRITE(22,*) IR, IEL, IP1, IP2, IPDW
        DO IP = IPDW, IP2
            IP2IR(IP) = IR
            !c      IF(LDBG) WRITE(22,*) IP,IR
        ENDDO
        !---- set pointers from panel center to grid center
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        DO IC = IC1, IC2
            IG1 = IP2IG(IPCO(IC))
            IG2 = IP2IG(IPCP(IC))
            IC2IG(IC) = MIN(IG1, IG2)
            !c        IF(IC2IG(IC).NE.0) WRITE(22,*) IC,IC2IG(IC)
        ENDDO
        !
    END SUBROUTINE ROTPINIT
    !*==ADDELEM5_2PT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ROTPINIT
    
    
    
    SUBROUTINE ADDELEM5_2PT(X1, Y1, X2, Y2, NPTS)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: NPTS
        REAL :: X1, X2, Y1, Y2
        !
        ! Local variables
        !
        REAL :: DX, DY, FNPT1
        INTEGER :: I, IEL, IP, IP1, IP2, N
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
        IF (NPTS<=1) THEN
            WRITE (*, *) 'Error in ADDELEM12 NPTS<2'
            STOP
        ENDIF
        !
        IP = IPLAST(NEL)
        IEL = NEL
        !
        IP1 = IP + 1
        IP2 = IP + NPTS
        !
        DX = X2 - X1
        DY = Y2 - Y1
        FNPT1 = FLOAT(NPTS - 1)
        DO I = 1, NPTS
            IP = IP + 1
            XP(IP) = X1 + DX * FLOAT(I - 1) / FNPT1
            YP(IP) = Y1 + DY * FLOAT(I - 1) / FNPT1
        ENDDO
        !---- centroid location
        XPCENT(IEL) = 0.5 * (XP(IP1) + XP(IP2))
        YPCENT(IEL) = 0.5 * (YP(IP1) + YP(IP2))
        !---- panel limits
        CALL MINMAX(N, XP(IP1), XPMINE(IEL), XPMAXE(IEL))
        CALL MINMAX(N, YP(IP1), YPMINE(IEL), YPMAXE(IEL))
        !
        !--- add element
        IEL = IEL + 1
        IPFRST(IEL) = IP1
        IPLAST(IEL) = IP2
        NETYPE(IEL) = 5
        LTPAN(IEL) = .FALSE.
        CALL XYPSPL(IEL)
        !
        IF (LDBG) THEN
            WRITE (LUNDBG, *) ' '
            WRITE (LUNDBG, *) 'Adding line source element ', IEL
            WRITE (LUNDBG, 100) 'Start ', IP1, X1, Y1
            WRITE (LUNDBG, 100) 'End   ', IP2, X2, Y2
            !
            WRITE (*, 10) IEL, NPTS, IP1, IP2
        ENDIF
        !
        10   FORMAT ('Line source element ', I4, ' added with #pts ', I5, ' IP1 ', &
                & I5, ' IP2 ', I5)
        100  FORMAT (A, I5, 4(1X, G13.5))
        !
        NEL = IEL
        NPTOT = IP2
        !
    END SUBROUTINE ADDELEM5_2PT
    !*==ADDELEM5.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ADDELEM5_2PT
    
    
    
    SUBROUTINE ADDELEM5(X, Y, N)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(*) :: X, Y
        !
        ! Local variables
        !
        INTEGER :: I, IEL, IP, IP1, IP2
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Add line source element (NETYPE=5) to paneled geometry
        !     Specified by X,Y coordinates and # of points
        !------------------------------------------------------
        !
        IF (N<=1) THEN
            WRITE (*, *) 'Error in ADDELEM5 NPTS<2'
            STOP
        ENDIF
        !
        !---- add sourceline points to panel geometry
        IP1 = IPLAST(NEL) + 1
        IP2 = IP1 - 1 + N
        IEL = NEL + 1
        IF (LDBG) WRITE (LUNDBG, 10) IEL, N, IP1, IP2
        !
        !---- add points to panel arrays
        DO IP = IP1, IP2
            I = IP - IP1 + 1
            XP(IP) = X(I)
            YP(IP) = Y(I)
            IF (LDBG) WRITE (LUNDBG, *) IP, XP(IP), YP(IP)
        ENDDO
        !---- centroid location
        XPCENT(IEL) = 0.5 * (XP(IP1) + XP(IP2))
        YPCENT(IEL) = 0.5 * (YP(IP1) + YP(IP2))
        !---- panel limits
        CALL MINMAX(N, XP(IP1), XPMINE(IEL), XPMAXE(IEL))
        CALL MINMAX(N, YP(IP1), YPMINE(IEL), YPMAXE(IEL))
        !
        !--- add element
        IPFRST(IEL) = IP1
        IPLAST(IEL) = IP2
        NETYPE(IEL) = 5
        LTPAN(IEL) = .FALSE.
        CALL XYPSPL(IEL)
        !
        10   FORMAT (/'Line source element ', I4, ' added with #pts ', I5, ' IP1 ', &
                & I5, ' IP2 ', I5)
        100  FORMAT (A, I5, 4(1X, G13.5))
        !
        NEL = IEL
        NPTOT = IP2
        !
    END SUBROUTINE ADDELEM5
    !*==ADDELEM6.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ADDELEM5
    
    
    
    SUBROUTINE ADDELEM6(X, Y, N)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(*) :: X, Y
        !
        ! Local variables
        !
        INTEGER :: I, IEL, IP, IP1, IP2
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
        IP1 = IPLAST(NEL) + 1
        IP2 = IP1 - 1 + N
        IEL = NEL + 1
        IF (LDBG) WRITE (LUNDBG, 10) IEL, N, IP1, IP2
        !
        !---- add points to panel arrays
        DO IP = IP1, IP2
            I = IP - IP1 + 1
            XP(IP) = X(I)
            YP(IP) = Y(I)
            IF (LDBG) WRITE (LUNDBG, *) IP, XP(IP), YP(IP)
        ENDDO
        !---- centroid location
        XPCENT(IEL) = 0.5 * (XP(IP1) + XP(IP2))
        YPCENT(IEL) = 0.5 * (YP(IP1) + YP(IP2))
        !---- panel limits
        CALL MINMAX(N, XP(IP1), XPMINE(IEL), XPMAXE(IEL))
        CALL MINMAX(N, YP(IP1), YPMINE(IEL), YPMAXE(IEL))
        !
        !---- setup pointers and panel data
        IPFRST(IEL) = IP1
        IPLAST(IEL) = IP2
        NETYPE(IEL) = 6
        LTPAN(IEL) = .FALSE.
        CALL XYPSPL(IEL)
        !      CALL XYCSET(IEL)
        !
        NEL = IEL
        NPTOT = IP2
        !
        10   FORMAT ('Rotor source line element ', I4, ' added with #pts ', I5, &
                &' IP1 ', I5, ' IP2 ', I5)
        !
    END SUBROUTINE ADDELEM6
    !*==DFLOAD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ADDELEM6
    
    
    
    
    SUBROUTINE DFLOAD(FNAMIN, FERROR)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        INTEGER, PARAMETER :: ITX = 300
        !
        ! Dummy arguments
        !
        LOGICAL :: FERROR
        CHARACTER(*) :: FNAMIN
        !
        ! Local variables
        !
        REAL :: A0, A0DEG, CDMIN, CLDMIN, CLMAX, CLMIN, CMCON, &
                & DCDCL2, DCDCL2S, DCLDA, DCLDA_STALL, DCL_STALL, &
                & FILEVERS, MCRIT, REREF, REXP, RPM1, RPM2, TOC, &
                & XISECT
        LOGICAL :: ERROR, LOPEN
        CHARACTER(128) :: FNAME, LINE
        INTEGER :: I, IB, IBNEXT, ICNT, ID, IEL, IFTYPE, IR, &
                & IRPN, IT, IV, K, LU, N, NBIN, ND1, NF, &
                & NINPUT, NR, NRPNIN, NRPNRD, NTEL
        INTEGER, DIMENSION(NEX) :: NT
        REAL, DIMENSION(10) :: RINPUT
        REAL, DIMENSION(ITX, NEX) :: XT, YT
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
        FERROR = .FALSE.
        NROTOR = 0
        LBLDEF = .FALSE.
        LCONV = .FALSE.
        NRP = 0
        NRC = 0
        DO N = 1, NRX
            IRTYPE(N) = 0
            NRDEF(N) = 0
            IRTYPDEF(N) = 0
        ENDDO
        !
        LRSPCDEF = .FALSE.
        LRSPCUSR = .FALSE.
        !
        NDOBJ = 0
        !
        LOPEN = .FALSE.
        LU = 1
        FNAME = FNAMIN
        CALL STRIP(FNAME, NF)
        IF (FNAME/=' ') THEN
            OPEN (LU, FILE = FNAME, STATUS = 'OLD', ERR = 98)
            LOPEN = .TRUE.
        ELSE
            RETURN
        ENDIF
        !
        ICNT = 0
        !
        !---- Read header line and version from DFDC input file
        CALL RDLINE(LU, LINE, ICNT)
        IF (LINE=='END' .OR. LINE=='ERR') GOTO 210
        !
        IF (LINE(1:4)/='DFDC') THEN
            WRITE (*, *) 'Not DFDC input file...may be duct geometry file'
            CLOSE (LU)
            !c        CALL LOADG(FNAME)
            RETURN
        ENDIF
        IV = INDEX(LINE, 'Version')
        READ (LINE(IV + 8:LEN(LINE)), *) FILEVERS
        WRITE (*, 1005) FILEVERS
        !
        1000 FORMAT (A)
        1005 FORMAT (/, ' Reading file from DFDC Version ', E8.2)
        1010 FORMAT (' File  ', A, ' not found'/)
        1020 FORMAT (' File  ', A, &
                &' has incompatible format'/' Loading not completed'/)
        1025 FORMAT (' Case Name: ', A)
        !
        !
        !--- Get Case title from line #2
        CALL RDLINE(LU, LINE, ICNT)
        IF (LINE=='END' .OR. LINE=='ERR') GOTO 210
        NAME = LINE
        CALL STRIP(NAME, NNAME)
        !
        WRITE (*, 1025) NAME
        DO
            !
            !-----------------------------------------------------
            !--- Main loop to find keywords and start read actions
            !    for data groups
            CALL RDLINE(LU, LINE, ICNT)
            IF (LINE=='END') THEN
                !
                !------------------------------------------------------------------------
                NBEL = NTEL
                !
                !---- Process airfoil inputs for CB and duct
                !---- put element-indexed coordinates XT,YT into linear bufffer airfoil
                !     geometry arrays XB,YB
                IBNEXT = 1
                DO IEL = 1, NBEL
                    !------ first point in element IEL
                    IBFRST(IEL) = IBNEXT
                    !------ initalize accumulation counter
                    IB = IBFRST(IEL) - 1
                    !------ go over all points on element IEL
                    IF (LDBG) WRITE (LUNDBG, *) 'Airfoil element'
                    DO IT = 1, NT(IEL)
                        IB = IB + 1
                        IF (IB>IBX) STOP 'LOAD: Array overflow on IBX'
                        XB(IB) = XT(IT, IEL)
                        YB(IB) = YT(IT, IEL)
                        IF (LDBG) WRITE (LUNDBG, *) IT, XT(IT, IEL), &
                                & YT(IT, IEL)
                    ENDDO
                    !---- set buffer airfoils to NBTYPE=0
                    NBTYPE(IEL) = 0
                    IBLAST(IEL) = IB
                    IBNEXT = IB + 1
                ENDDO
                !---- set total number of points in buffer geometry
                NBTOT = IBLAST(NBEL)
                !
                !
                !=======================================================================
                !---- If panel spacing data is not in case file try reading explicitly
                !     from xxx.pan paneling input file
                IF (.NOT.LRSPCDEF) THEN
                    !---- get file prefix (if it has a .xxx suffix)
                    K = INDEX(FNAME, '.') - 1
                    IF (K<=0) THEN
                        !----- no "." prefix/suffix separator
                        !      just tack on ".pan" to full filename
                        CALL STRIP(FNAME, NF)
                        K = NF
                    ENDIF
                    PREFIX = FNAME(1:K)
                    NPREFIX = K
                    PFILE = PREFIX(1:NPREFIX) // '.pan'
                    CALL PANGET(PFILE, ERROR)
                    !
                    IF (ERROR) THEN
                        WRITE (*, 2100) PFILE(1:K + 4)
                    ELSE
                        WRITE (*, 2110) PFILE(1:K + 4)
                    ENDIF
                ENDIF
                !
                2100       FORMAT (' No repaneling parameters from file ', A)
                2110       FORMAT (' Repaneling parameters from file ', A)
    
                !
                !=======================================================================
                !---- process current buffer geometry to identify element types and
                !     set element parameters
                CALL GEPROC
                IF (LDBG) CALL GEPRINT(LUNDBG)
                !
                !---- process buffer foil geometry and set up paneled geometry and grid
                !c    CALL GENGEOM
                !
                !---- take note of new case
                LLOAD = .TRUE.
                DO I = 1, NROTOR
                    LBBLOFT(I) = .FALSE.
                ENDDO
                GOTO 300
            ELSEIF (LINE=='ERR') THEN
                GOTO 210
                !
                !
                !---- OPER data
            ELSEIF (INDEX(LINE, 'OPER')/=0) THEN
                CALL RDLINE(LU, LINE, ICNT)
                !c        READ(LINE,*,ERR=210) QINF,QREF,RPM
                NINPUT = 4
                CALL GETFLT(LINE, RINPUT, NINPUT, ERROR)
                IF (NINPUT>=4) THEN
                    QINF = RINPUT(1)
                    QREF = RINPUT(2)
                    RPM1 = RINPUT(3)
                    RPM2 = RINPUT(4)
                ELSEIF (NINPUT>=3) THEN
                    QINF = RINPUT(1)
                    QREF = RINPUT(2)
                    RPM1 = RINPUT(3)
                    RPM2 = 0.0
                ELSEIF (NINPUT>=2) THEN
                    QINF = RINPUT(1)
                    QREF = RINPUT(2)
                    RPM1 = 0.0
                    RPM2 = 0.0
                ENDIF
                !
                CALL RDLINE(LU, LINE, ICNT)
                READ (LINE, *, ERR = 210) RHO, VSO, RMU, ALTH
                CALL RDLINE(LU, LINE, ICNT)
                READ (LINE, *, ERR = 210) XDWKLEN, NWAKE
                CALL RDLINE(LU, LINE, ICNT)
                READ (LINE, *, ERR = 210) LWRLX
                CALL RDLINE(LU, LINE, ICNT)
                IF (INDEX(LINE, 'ENDOPER')/=0) THEN
                    WRITE (*, 1030) QINF, QREF, RPM1, RPM2, RHO, VSO, &
                            & RMU, ALTH, XDWKLEN, NWAKE
                ELSE
                    WRITE (*, *) 'No ENDOPER to OPER section'
                    STOP
                ENDIF
                !
                1030       FORMAT (' OPER data read from file', /, '  Qinf   = ', F9.3, 5X, &
                        &'Qref   = ', F9.3, /, '  RPM1   = ', F9.1, 5X, &
                        & 'RPM2   = ', F9.1, /, '  Rho    = ', F9.5, 5X, &
                        & 'VSound = ', F9.3, /, '  rMU    = ', E9.3, 5X, &
                        & 'Alt    = ', F9.5, /, '  Xdwake = ', F9.5, 5X, &
                        & 'Nwake  = ', I6)
                !
                !---- Input RPM's set speed of disks
                OMEGA(1) = PI * RPM1 / 30.0
                OMEGA(2) = PI * RPM2 / 30.0
                IF (RHO==0.0 .OR. VSO==0.0 .OR. RMU==0.0) THEN
                    DELTAT = 0.0
                    CALL ATMO(ALTH, DELTAT, VSO, RHO, RMU)
                ENDIF
                !
                !
                !---- AERO data
            ELSEIF (INDEX(LINE, 'AERO')/=0) THEN
                !--- Read aero section definitions
                NR = NROTOR + 1
                CALL RDLINE(LU, LINE, ICNT)
                READ (LINE, *, ERR = 210) NAERO(NR)
                DO N = 1, NAERO(NR)
                    CALL RDLINE(LU, LINE, ICNT)
                    READ (LINE, *, ERR = 210) XISECT
                    CALL RDLINE(LU, LINE, ICNT)
                    READ (LINE, *, ERR = 210) A0DEG, DCLDA, CLMAX, CLMIN
                    CALL RDLINE(LU, LINE, ICNT)
                    READ (LINE, *, ERR = 210) DCLDA_STALL, DCL_STALL, CMCON, &
                            & MCRIT
                    CALL RDLINE(LU, LINE, ICNT)
                    READ (LINE, *, ERR = 210) CDMIN, CLDMIN, DCDCL2
                    CALL RDLINE(LU, LINE, ICNT)
                    !c          READ(LINE,*,ERR=210) REREF,REXP
                    NINPUT = 4
                    CALL GETFLT(LINE, RINPUT, NINPUT, ERROR)
                    IF (NINPUT>=4) THEN
                        REREF = RINPUT(1)
                        REXP = RINPUT(2)
                        TOC = RINPUT(3)
                        DCDCL2S = RINPUT(4)
                    ELSEIF (NINPUT>=3) THEN
                        REREF = RINPUT(1)
    
                        REXP = RINPUT(2)
                        TOC = RINPUT(3)
                        DCDCL2S = 0.020
                    ELSEIF (NINPUT>=2) THEN
                        REREF = RINPUT(1)
                        REXP = RINPUT(2)
                        TOC = 0.1
                        DCDCL2S = 0.020
                    ENDIF
                    A0 = A0DEG * DTR
                    CALL PUTAERO(NR, N, XISECT, A0, CLMAX, CLMIN, DCLDA, &
                            & DCLDA_STALL, DCL_STALL, CDMIN, CLDMIN, DCDCL2, &
                            & DCDCL2S, CMCON, MCRIT, TOC, REREF, REXP)
                ENDDO
                CALL RDLINE(LU, LINE, ICNT)
                IF (INDEX(LINE, 'ENDAERO')/=0) THEN
                    WRITE (*, 2010) NR, NAERO(NR)
                ELSE
                    WRITE (*, *) 'No ENDAERO to AERO section'
                    STOP
                ENDIF
                !
                2010       FORMAT (/, ' AERO data read for Disk', I2, ': ', I3, ' sections')
                !
                !
                !---- ROTOR data
            ELSEIF (INDEX(LINE, 'ROTOR')/=0) THEN
                NR = NROTOR + 1
                CALL RDLINE(LU, LINE, ICNT)
                READ (LINE, *, ERR = 210) XDISK(NR), NRBLD(NR), NRSTA
                CALL RDLINE(LU, LINE, ICNT)
                READ (LINE, *, ERR = 210) NRDEF(NR)
                !--- Get rotor definition data
                DO IR = 1, NRDEF(NR)
                    CALL RDLINE(LU, LINE, ICNT)
                    READ (LINE, *, END = 210, ERR = 210) YRDEF(IR, NR), &
                            & CHRDEF(IR, NR), BETADEF(IR, NR)
                ENDDO
                !
                CALL RDLINE(LU, LINE, ICNT)
                IF (INDEX(LINE, 'ENDROTOR')/=0) THEN
                    WRITE (*, 2020) XDISK(NR), NRBLD(NR), NRSTA
                    IF (NRDEF(NR)>=2) THEN
                        WRITE (*, 2030) NRDEF(NR)
                        DO IR = 1, NRDEF(NR)
                            WRITE (*, 15) YRDEF(IR, NR), CHRDEF(IR, NR), &
                                    & BETADEF(IR, NR)
                            BETADEF(IR, NR) = BETADEF(IR, NR) * DTR
                        ENDDO
                        15               FORMAT (1X, F12.6, 2X, F12.6, 4X, F12.6)
                        IRTYPDEF(NR) = 2
                        NROTOR = NR
                    ELSE
                        WRITE (*, *) 'Rotor blade defined by too few stations'
                    ENDIF
                ELSE
                    WRITE (*, *) 'No ENDROTOR found'
                    STOP
                ENDIF
                !
                2020       FORMAT (/, ' ROTOR data read from file', /, '  XDisk: ', F5.3, &
                        &'   Blades:', I2, '   Rotor pnts:', I3)
                2030       FORMAT ('  Rotor blade defined with', I3, ' points', /, &
                        &'         R            CH            Beta')
                !
                !
                !---- ACTDISK data
            ELSEIF (INDEX(LINE, 'ACTDISK')/=0) THEN
                NR = NROTOR + 1
                CALL RDLINE(LU, LINE, ICNT)
                READ (LINE, *, ERR = 210) XDISK(NR), NRSTA
                CALL RDLINE(LU, LINE, ICNT)
                READ (LINE, *, ERR = 210) NRDEF(NR)
                !--- Get actuator disk definition data
                DO IR = 1, NRDEF(NR)
                    CALL RDLINE(LU, LINE, ICNT)
                    READ (LINE, *, END = 210, ERR = 210) YRDEF(IR, NR), &
                            & BGAMDEF(IR, NR)
                ENDDO
                !
                CALL RDLINE(LU, LINE, ICNT)
                IF (INDEX(LINE, 'ENDACTDISK')/=0) THEN
                    WRITE (*, 2040) XDISK(NR), NRSTA
                    IF (NRDEF(NR)>=2) THEN
                        WRITE (*, 2050) NRDEF(NR)
                        DO IR = 1, NRDEF(NR)
                            WRITE (*, 15) YRDEF(IR, NR), BGAMDEF(IR, NR)
                        ENDDO
                        IRTYPDEF(NR) = 1
                        NROTOR = NR
                    ELSE
                        WRITE (*, *)                                           &
                                &'Actuator disk defined by too few stations'
                    ENDIF
                ELSE
                    WRITE (*, *) 'No ENDACTDISK found'
                    STOP
                ENDIF
                !
                2040       FORMAT (/, ' ACTDISK data read from file', /, '  XDisk: ', F5.3, &
                        &'   Rotor points:', I3)
                2050       FORMAT ('  Actuator Disk defined with', I3, ' points', /, &
                        &'         R           BGAM')
                !
                !
                !---- DRAGOBJ data
            ELSEIF (INDEX(LINE, 'DRAGOBJ')/=0) THEN
                ND1 = NDOBJ + 1
                IF (ND1>NDRGX) THEN
                    WRITE (*, *)                                              &
                            &'Number of drag objects exceeds dimension NDGX'
                ELSE
                    CALL RDLINE(LU, LINE, ICNT)
                    READ (LINE, *, ERR = 210) NDDEF(ND1)
                    !--- Get drag object definition data
                    DO ID = 1, NDDEF(ND1)
                        CALL RDLINE(LU, LINE, ICNT)
                        READ (LINE, *, END = 210, ERR = 210) XDDEF(ID, ND1), &
                                & YDDEF(ID, ND1), CDADEF(ID, ND1)
                    ENDDO
                    !
                    CALL RDLINE(LU, LINE, ICNT)
                    IF (INDEX(LINE, 'ENDDRAGOBJ')/=0) THEN
                        NDOBJ = ND1
                        WRITE (*, 2060) ND1
                        IF (NDDEF(ND1)>=2) THEN
                            WRITE (*, 2070) NDDEF(ND1)
                            DO ID = 1, NDDEF(ND1)
                                WRITE (*, 15) XDDEF(ID, ND1), YDDEF(ID, ND1), &
                                        & CDADEF(ID, ND1)
                            ENDDO
                        ELSE
                            WRITE (*, *)                                        &
                                    &'Drag object defined by too few stations'
                        ENDIF
                    ELSE
                        WRITE (*, *) 'No ENDDRAGOBJ found'
                        STOP
                    ENDIF
                    ! ENDDRAGOBJ check
                    !
                ENDIF
                ! NDRGX check
                !
                2060       FORMAT (/, ' DRAGOBJ ', I2, '  read from file')
                2070       FORMAT ('  Drag object defined with', I3, ' points', /, &
                        &'         X             R              CDA')
                !
                !
                !---- GEOM data (CB and duct coordinates in a single XFOIL multi-element file)
            ELSEIF (INDEX(LINE, 'GEOM')/=0) THEN
                !
                !---- Read the combined CB and duct airfoil file
                !
                CALL AREADNR(LU, ITX, NEX, XT, YT, NT, NTEL, ANAME, ISPARS, IFTYPE)
                IF (IFTYPE==0 .OR. NTEL/=2) THEN
                    !----- read error occurred for two elements
                    WRITE (*, 2080) NTEL
                    STOP
                ENDIF
                !
                2080       FORMAT (' Error reading GEOM data', /, &
                        & ' Coordinates read for', I2, ' foils')
                !
                !---- User-defined respacing specs
                !---- PANE data (respacing data for surface points
            ELSEIF (INDEX(LINE, 'PANE')/=0) THEN
                CALL RDLINE(LU, LINE, ICNT)
                READ (LINE, *, ERR = 210) NBIN, NRPNIN
                NRPNRD = MIN(NRPNIN, NRPNX)
                DO IEL = 1, NBIN
                    CALL RDLINE(LU, LINE, ICNT)
                    READ (LINE, *, ERR = 210) NPAN(IEL)
                    CALL RDLINE(LU, LINE, ICNT)
                    READ (LINE, *, ERR = 210) CVEX(IEL), SMOF(IEL), FSLE(IEL), &
                            & FSTE(IEL)
                    DO IRPN = 1, NRPNRD
                        CALL RDLINE(LU, LINE, ICNT)
                        READ (LINE, *, ERR = 210) SRPN1(IRPN, IEL), &
                                & SRPN2(IRPN, IEL), &
                                & CRRAT(IRPN, IEL)
                    ENDDO
                ENDDO
                CALL RDLINE(LU, LINE, ICNT)
                IF (INDEX(LINE, 'ENDPANE')/=0) THEN
                    WRITE (*, 2090) NBIN
                    LRSPCDEF = .TRUE.
                    LRSPCUSR = .TRUE.
                ENDIF
                !
            ENDIF
            !
            2090    FORMAT (' PANELING for', I2, ' elements read from file')
        ENDDO
        !...............................................................
        98   WRITE (*, 1050) FNAME(1:NF)
        FERROR = .TRUE.
        GOTO 300
        !
        WRITE (*, 1100) FNAME(1:NF)
        FERROR = .TRUE.
        GOTO 300
        !
        210  WRITE (*, 1150) ICNT, FNAME(1:NF)
        FERROR = .TRUE.
        !
        1050 FORMAT (/' File OPEN error:  ', A)
        1100 FORMAT (/' File READ error:  ', A)
        1150 FORMAT (/' File READ error on line ', I3, ':  ', A)
        !
        !---- close file
        300  IF (LOPEN) CLOSE (LU)
    END SUBROUTINE DFLOAD
    !*==DFSAVE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! DFLOAD
    
    
    
    SUBROUTINE DFSAVE(FNAMEIN)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        CHARACTER(*) :: FNAMEIN
        !
        ! Local variables
        !
        REAL :: A0, A0DEG, BETADEG, CDMIN, CLDMIN, CLMAX, CLMIN, &
                & CMCON, DCDCL2, DCDCL2S, DCLDA, DCLDA_STALL, &
                & DCL_STALL, MCRIT, REREF, REXP, RPM1, RPM2, TOC, &
                & XISECT
        CHARACTER(1) :: ANS
        CHARACTER(80) :: FBLNK
        CHARACTER(128) :: FNAME
        INTEGER :: I, IEL, IFTYPE, IRPN, LU, N, NELSAV, NF, NR
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------------------------
        !     Save rotor and operating state in DFDC Version 0.5 format
        !     This format saves data in groups delimited by keywords
        !     to separate data
        !--------------------------------------------------------------------------
        !
        LU = 2
        FNAME = FNAMEIN
        CALL STRIP(FNAME, NF)
        !
        IF (FNAME(1:1)==' ') CALL ASKS('Enter filename^', FNAME)
        OPEN (LU, FILE = FNAME, STATUS = 'OLD', ERR = 5)
        WRITE (*, *)
        WRITE (*, *) 'Output file exists.  Overwrite?  Y'
        READ (*, 1000) ANS
        IF (INDEX('Nn', ANS)==0) GOTO 6
        !
        CLOSE (LU)
        WRITE (*, *) 'Current duct case not saved.'
        RETURN
        !
        5    OPEN (LU, FILE = FNAME, STATUS = 'NEW', ERR = 90)
        6    REWIND (LU)
        !
        !
        !--- Version header and case name
        IF (NAME==' ') NAME = 'Saved ducted fan'
        WRITE (LU, 1100) VERSION, NAME
        !
        !--- OPER data
        !
        IF (NROTOR<2) OMEGA(2) = 0.
        !
        !--- Velocity, reference velocity and RPM
        WRITE (LU, 1102)
        RPM1 = 30.0 * OMEGA(1) / PI
        RPM2 = 30.0 * OMEGA(2) / PI
        WRITE (LU, 1200) QINF, QREF, RPM1, RPM2
        !--- Altitude and atmospheric data
        WRITE (LU, 1103)
        WRITE (LU, 1200) RHO, VSO, RMU, ALTH
        !--- XDwake, #wake points
        WRITE (LU, 1104)
        WRITE (LU, 1202) XDWKLEN, NWAKE
        !--- Wake relaxation flag
        WRITE (LU, 1105)
        WRITE (LU, 1204) LWRLX
        WRITE (LU, 1106)
        !
        !--- Save data for each disk (act disk or bladed)
        DO NR = 1, NROTOR
            !
            !--- AERO data
            !--- Save aero data for defined aero sections
            WRITE (LU, 1108) NAERO(NR)
            DO N = 1, NAERO(NR)
                CALL GETAERO(NR, N, XISECT, A0, CLMAX, CLMIN, DCLDA, DCLDA_STALL, &
                        & DCL_STALL, CDMIN, CLDMIN, DCDCL2, DCDCL2S, CMCON, &
                        & MCRIT, TOC, REREF, REXP)
                WRITE (LU, 1200) XISECT
                A0DEG = A0 * 180.0 / PI
                WRITE (LU, 1109)
                WRITE (LU, 1200) A0DEG, DCLDA, CLMAX, CLMIN
                WRITE (LU, 1110)
                WRITE (LU, 1200) DCLDA_STALL, DCL_STALL, CMCON, MCRIT
                WRITE (LU, 1111)
                WRITE (LU, 1200) CDMIN, CLDMIN, DCDCL2
                WRITE (LU, 1112)
                WRITE (LU, 1200) REREF, REXP, TOC, DCDCL2S
            ENDDO
            WRITE (LU, 1113)
            !
            !--- Save rotor blade if defined
            IF (IRTYPE(NR)==2) THEN
                WRITE (LU, 1130)
                !--- Rotor axial location, #blades, #radial stations
                WRITE (LU, 1202) XDISK(NR), NRBLD(NR), NRSTA
                !--- Save blade definition with chord,twist and body velocity
                WRITE (LU, 1132) NRC
                DO I = 1, NRC
                    BETADEG = BETAR(I, NR) / DTR
                    WRITE (LU, 1200) YRC(I, NR), CHR(I, NR), BETADEG
                ENDDO
                WRITE (LU, 1134)
            ENDIF
            !
            !--- Save actuator disk if defined
            IF (IRTYPE(NR)==1) THEN
                WRITE (LU, 1122)
                !--- Rotor axial location, #radial stations
                WRITE (LU, 1202) XDISK(NR), NRSTA
                !--- Save actuator disk circulation
                WRITE (LU, 1124) NRC
                DO I = 1, NRC
                    WRITE (LU, 1200) YRC(I, NR), BGAM(I, NR)
                ENDDO
                WRITE (LU, 1126)
            ENDIF
            !
        ENDDO
        !
        !--- Save drag object(s) if defined
        IF (NDOBJ>0) THEN
            DO N = 1, NDOBJ
                IF (NDDEF(N)>2) THEN
                    WRITE (LU, 1136) NDDEF(N)
                    !--- Save drag object points and CDA values
                    DO I = 1, NDDEF(N)
                        WRITE (LU, 1200) XDDEF(I, N), YDDEF(I, N), CDADEF(I, N)
                    ENDDO
                    WRITE (LU, 1137)
                ENDIF
            ENDDO
        ENDIF
        !
        !--- CB and duct geometry
        WRITE (LU, 1140)
        !---- save only 2 elements in buffer geometry into opened file
        NELSAV = 2
        IFTYPE = 2
        FBLNK = ' '
        CALL AWRITE(FBLNK, LU, NELSAV, IPFRST, IPLAST, XP, YP, ANAME, ISPARS, &
                & IFTYPE)
        WRITE (LU, 1142)
        !
        !--- Save user-defined paneling specs
        IF (LRSPCUSR) THEN
            WRITE (LU, 1150) NBEL, NRPNX
            DO IEL = 1, NBEL
                WRITE (LU, 1152) NPAN(IEL), CVEX(IEL), SMOF(IEL), &
                        & FSLE(IEL), FSTE(IEL)
                WRITE (LU, 1154)
                DO IRPN = 1, NRPNX
                    WRITE (LU, 1200) SRPN1(IRPN, IEL), SRPN2(IRPN, IEL), &
                            & CRRAT(IRPN, IEL)
                ENDDO
            ENDDO
            WRITE (LU, 1156)
        ENDIF
        !
        CLOSE (LU)
        RETURN
        !
        90   WRITE (*, *) 'Bad filename.'
        WRITE (*, *) 'Current duct case not saved.'
        RETURN
        !
        !...................................................................
        1000 FORMAT (A)
        1100 FORMAT ('DFDC Version ', E8.2, ' '/A32)
        !
        1102 FORMAT (/'OPER', /'!   Vinf         Vref         RPM1         RPM2'&
                &)
        1103 FORMAT ('!   Rho          Vso          Rmu          Alt')
        1104 FORMAT ('!  XDwake             Nwake')
        1105 FORMAT ('!        Lwkrlx')
        1106 FORMAT ('ENDOPER')
        !
        1108 FORMAT (/'AERO', /'!  #sections'/1(1X, I5), /'!  Xisection')
        1109 FORMAT ('!  A0deg        dCLdA        CLmax        CLmin')
        1110 FORMAT ('! dCLdAstall   dCLstall      Cmconst      Mcrit')
        1111 FORMAT ('!   CDmin      CLCDmin       dCDdCL^2')
        1112 FORMAT ('!   REref        REexp        TOC         dCDdCL^2')
        1113 FORMAT ('ENDAERO')
        !
        1122 FORMAT (/'ACTDISK', /'!  Xdisk       NRsta')
        1124 FORMAT ('!  #stations'/1(1X, I5), /'!           r         BGam')
        1126 FORMAT ('ENDACTDISK')
        !
        1130 FORMAT (/'ROTOR', /'!  Xdisk               Nblds       NRsta')
        1132 FORMAT ('!  #stations'/1(1X, I5), &
                &/'!     r        Chord         Beta')
        1134 FORMAT ('ENDROTOR')
        !
        1136 FORMAT (/'DRAGOBJ', /'!  #pts'/1(1X, I5), &
                &/'!     x            r            CDA')
        1137 FORMAT ('ENDDRAGOBJ')
        !
        1140 FORMAT (/'GEOM')
        1142 FORMAT ('ENDGEOM')
        !
        1150 FORMAT (/'PANELING', /'!  #elements   #refinement zones'/(2(1X, I5))&
                &)
        1152 FORMAT ('!  #panel nodes'/1(1X, I5)                                &
                &/'!  curv_expon  curv_smooth   dsL/dsAvg    dsR/dsAvg', &
                & /4(1X, G12.5))
        1154 FORMAT ('!  s1/smax      s2/smax     ds/dsAvg')
        1156 FORMAT ('ENDPANELING')
        !
        1200 FORMAT (5(1X, G12.5))
        1202 FORMAT (1(1X, G12.5), 2(8X, I5))
        1204 FORMAT (12X, L1)
        !
        !x123456789012x123456789012x123456789012x123456789012x123456789012
        !!         Rho          Vso          Rmu           Alt')
        !
    END SUBROUTINE DFSAVE
end module m_dfdcsubs
