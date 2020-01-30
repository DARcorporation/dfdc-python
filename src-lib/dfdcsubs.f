
C=========================================================================
C DFDC (Ducted Fan Design Code) is an aerodynamic and aeroacoustic design
C and analysis tool for aircraft with propulsors in ducted fan
C configurations.
C 
C This software was developed under the auspices and sponsorship of the
C Tactical Technology Office (TTO) of the Defense Advanced Research
C Projects Agency (DARPA).
C 
C Copyright (c) 2004, 2005, Booz Allen Hamilton Inc., All Rights Reserved
C
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version.
C 
C This program is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C 
C You should have received a copy of the GNU General Public License along
C with this program; if not, write to the Free Software Foundation, Inc.,
C 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
C
C Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (drela@mit.edu)
C Program Management: Brad Tousley, Paul Eremenko (eremenko@alum.mit.edu)
C
C=========================================================================
C
C     Version 070-ES1
C     Philip Carter, Esotec Developments, February 2009
C     philip (at) esotec (dot) org
C
C     Changes from 0.70:
C
C     DFINIT
C       SVERSION (line 66)
C       LCOORD   (line 424)
C       LLOAD    (line 436)
C       IDEVPS = 4   ! Color PostScript
C
C     DFLOAD
C       File version displayed in exponential format on load.
C       Reporting formats cleaned up.
C
C     DFSAVE
C       File version written in exponential format (allows ES versions).
C       ANAME (geometry name) correctly saved.
C       RPM2 written. 
C
C     Version 070-ES2
C     ROTINIT   LLOFT condition imposed on tip gap paneling (for lofting)
C     
C=========================================================================
C
      SUBROUTINE DFINIT(LDEBUG)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      LOGICAL LDEBUG
C--------------------------------------
C     Sets all initial defaults
C--------------------------------------
      CHARACTER*10 DIGITS
      DATA DIGITS / '0123456789' /
C
C---- Set version (in case this code is used by DLL w/o main DFDC routine)
C
C      VERSION  =  0.70e02          ! This gets saved to file
C      SVERSION = '070-ES2'         ! A more flexible format
C
C---------------------------------------------------------------
C---- Debug output
ccc      LDBG = .FALSE.
      LDBG = LDEBUG
C---- Logical unit for debugging printout
      LUNDBG = 40
C---------------------------------------------------------------
C
C---- set digit string array for constructing filename
      DO KP = 1, 99
        K2 =          KP/10   + 1
        K1 = KP - 10*(KP/10)  + 1
        PNUM(KP) = DIGITS(K2:K2) // DIGITS(K1:K1)
      ENDDO
C
C---------------------------------------------------------------
C
      DO IP = 1, IPX
        GAM(IP) = 0.
        SIG(IP) = 0.
        GTH(IP) = 0.
        GAMVSP(IP) = 0.
        SIGVSP(IP) = 0.
        VMAVG(IP)  = 0.
        IZERO(IP)  = 0
        IP2IG(IP)  = 0
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
C
C---- no panel nodes, control points
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
C
        IEL2IR(IEL) = 0
      ENDDO
C
C---------------------------------------------------------------
C---- no flow yet
      QINF = 0.0
C---- start off with incompressible flow
      MACH  = 0.
      MACH1 = 0.
C---- acceleration due to gravity for scaling centrifugal blade tension (m/s^2)
      GEE = 9.81
C---- default unity reference velocity
      QREF = 1.0
C---- setup for SL conditions, Std temperature
      ALTH   = 0.0
      DELTAT = 0.0
C! sea level atmosphere parameters
C     RHO =  1.226      ! fluid density         kg/m**3
C     RMU =  1.78E-05   ! dynamic viscosity     kg/m-s
C     VSO =  340.0      ! speed of sound        m/s
      CALL ATMO(ALTH,DELTAT,VSO,RHO,RMU) 
C---------------------------------------------------------------
C
C---- no drag objects defined
      NDOBJ  = 0
      IELDRGOBJ(1) = 0
      LDRGOBJ = .FALSE.
C
C---------------------------------------------------------------
C---- no rotor yet
      NROTOR = 0
C---- no blade defined yet
      LBLDEF = .FALSE.
C
      DO N = 1, NRX
C---- disk type set to undefined
       IRTYPE(N)   = 0
C---- no defining points for rotor or actuator disk
       NRDEF(N)    = 0
       IRTYPDEF(N) = 0
C---- no elements defined for rotor source line
       IELROTOR(N) = 0
C
       ATIP(N)  = 0.0
       ADISK(N) = 0.0
       XDISK(N) = 999.   ! default rotor disk location
       NRBLD(N) = 5
       OMEGA(N) = 0.0
      END DO
C---- default blade is rotor with RPM=1000
      RPM = 1000.0
      OMEGA(1) = PI*RPM/30.0
C
      TGAP   = 0.0    ! no tip gap on rotor
      TGAPZL = 0.0    ! tip gap of zero loss (default is zero)
C
      NRP = 0
      NRC = 0 
C
      NRSTA = 11     ! default # of radial stations
      XDWKLEN = 1.5  ! default wake length in duct diameters
      XWAKE  = -999.
      NWAKE  = 0
C
      TDUCT = 0.0
      TTOT = 0.0
      TVIS = 0.0
      QTOT = 0.0
      QVIS = 0.0
      PTOT = 0.0
      PVIS = 0.0
C
C
      DO IR = 1, IRX
        CHDES(IR) = 1.0
        CLDES(IR) = 0.8
        CLPOS(IR) =  CLDES(IR)
        CLNEG(IR) = -CLDES(IR)
        IR2IEL(IR) = 0
C
        DO N = 1, NRX
         BGAM(IR,N) = 0.0
         DO L = 1, 3
           VIND(L,IR,N) = 0.0
           VABS(L,IR,N) = 0.0
           VREL(L,IR,N) = 0.0
         END DO
         CHR(IR,N) = CHDES(IR)
         CLR(IR,N) = CLDES(IR)
       END DO
      END DO
C
C---------------------------------------------------------------
C--- Setup aerodynamic data for blade elements
C--- Default aero properties 
C
      A0    = 0.         ! zero lift angle of attack   radians
      DCLDA =  6.0       ! lift curve slope            /radian
      CLMAX =   1.5      ! stall Cl
      CLMIN =  -0.5      ! negative stall Cl
      DCL_STALL =   0.1  ! CL increment from incipient to total stall
      DCLDA_STALL = 0.1  ! stalled lift curve slope    /radian
      CMCON   = -0.1     ! section Cm  (for pitch-axis moments)
      CDMIN   =  0.013   ! minimum Cd
      CLDMIN  =  0.5     ! Cl at minimum Cd
      DCDCL2  =  0.03    ! d(Cd)/d(Cl**2)
      DCDCL2S =  0.0     ! d(Cd)/d(Cl**2) secondary (annulus) drag (dead)
      REREF =  200000.   ! Reynolds Number at which Cd values apply
      REXP  =  -0.4      ! Exponent for Re scaling of Cd:  Cd ~ Re**exponent
      MCRIT =  0.8       ! Critical Mach number
      TOC   =  0.1       ! section thickness ratio (dead)
C
      XPAXIS = 0.35     ! x/c location of pitch axis for blade plots, lofting
C
C--- Install default data into aero section #1 for each disk 
      DO N = 1, NRX
       NAERO(N) = 1
       XISECT = 0.0
       CALL PUTAERO(N,NAERO(N),XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &              CMCON,MCRIT,TOC,REREF,REXP)
C--- set pointers from blade sections to aero sections
       DO I=1, IRX
         IAERO(I,N) = 1
       END DO
      END DO
C
C---- no BL solution yet
      LVISC = .FALSE.
      AMPMAX = 9.0
      UREN   = 0.0
      DO IEL = 1, NEX
        CXVIS(IEL)  = 0.0
        SFIX(1,IEL) = 1.0
        SFIX(2,IEL) = 1.0
        STRN(1,IEL) = 1.0
        STRN(2,IEL) = 1.0
        SSEP(1,IEL) = 1.0
        SSEP(2,IEL) = 1.0
        ICBLBEG(1,IEL) = 0
        ICBLBEG(2,IEL) = 0
        ICBLEND(1,IEL) = 0
        ICBLEND(2,IEL) = 0
      END DO
C---- no inflow velocities defined
      NINFL  = 0
C
C---------------------------------------------------------------
C---- no grid dimensions yet
      II = 0
      JJ = 0
C
C---- by default input geometry is respaced
      LRSPC = .TRUE.
C---------------------------------------------------------------
C---- default number of panels per element, panel-number exponent
      NPANDEF = 61
      FPANDEF = 0.7
C---- set panel parameters, no default spacing or user-defined spacing 
      LRSPCUSR = .FALSE.
      LRSPCDEF = .FALSE.
      DO IEL = 1, NEX
        NPAN(IEL) = -999
        CVEX(IEL) = 1.0
        SMOF(IEL) = 1.0
        FSLE(IEL) = 0.6
        FSTE(IEL) = 0.6
C------ refinement-station parameters
        DO IRPN = 1, NRPNX
          CRRAT(IRPN,IEL) = 1.0
          SRPN1(IRPN,IEL) = 0.0
          SRPN2(IRPN,IEL) = 0.0
        ENDDO
      ENDDO
C
C---- assume no case-paneling parameter file exists
      PFILE = ' '
C
C---------------------------------------------------------------
C---- maximum TE gap (relative to perimeter) to deem airfoil as being closed
      DTEBOD = 0.1
C---- minimum TE gap (relative to perimeter) to give closed airfoil a TE panel
      DTEPAN = 0.00001
C
C---------------------------------------------------------------
C---- no geometry yet
      NEL = 0
      NBEL = 0
      NAME = ' '
C
C---- no accumulated point sequence yet
      NPOL = 0
      LPACC = .FALSE.
C
      IPACT = 0
C
      NPSEQ = 0
C
C---------------------------------------------------------------
C---- Convergence flag for iterative solutions, no converged solution yet
      LCONV = .FALSE.
C---- Iterative solver parameters
      RLXSOLV = 0.4
      EPSSOLV = 0.0002
      ITRMAXSOLV = 50
C---- VMAVG is not initialized
      LVMAV = .FALSE.
      VAVGINIT = QINF + 1.0
C
C---- vortex wakes will not be automatically realigned 
      LWRLX = .FALSE.
C
C---------------------------------------------------------------
C---- don't normalize airfoil upon input
      LNORM = .FALSE.
C
C---- no geometry modes defined yet
      NMOD = 0
      LXYMOV = .FALSE.
C
C---------------------------------------------------------------
C---- no valid control-point pointers yet
      LNCVP = .FALSE.
C
C---- no valid Aero,Geometry influence matrices yet
      LQAIC = .FALSE.
      LQGIC = .FALSE.
      LVAIC = .FALSE.
      LVGIC = .FALSE.
C
C---- no valid control-point velocities yet
      LQCNT = .FALSE.
C
C---- no valid system pointers yet
      LSYSP = .FALSE.
C
C---- no factored system matrix for gamma yet
      LGSYS = .FALSE.
      LQSYS = .FALSE.
C
C---- no gamma solution yet, no unit solutions
      LGAMA = .FALSE.
      LGAMU = .FALSE.
C
C---- no source-influence pointers or matrix
      LSIGP = .FALSE.
      LSIGM = .FALSE.
C
C---- don't plot reference Cp or force data
      LPREF = .FALSE.
      LFREF = .FALSE.
C
C---- default filename prefix
      PREFIX = ' '
C
      CPFILE = ' '
      FRFILE = ' '
C---------------------------------------------------------------
C---- no target elements yet for QDES or GDES
      IELQDES = 0
      IELGDES = 1
C
C---- enforce geometry slope-matching at grafted-geometry endpoints
      LGSLOP = .TRUE.
C
C---- no geometric symmetry assumed
      LGSYMM = .FALSE.
C
C---- same geometry between panel and buffer
      LGSAME = .FALSE.
C
C---- plot home positions
      LHOMPL = .TRUE.
C
C---- no specified speed distributions yet
      NQSP = 0
      NSEG = 0
      LQSPEC = .FALSE.
      LQSPPL = .FALSE.
      LGSPPL = .FALSE.
C
C---- enforce slope-matching and gamma curvature at inverse-segment endpoints
      LQSLOP = .TRUE.
      LQCURV = .FALSE.
C
C---- inverse-segment endpoints on edge of a surface airfoil are free to move
      LGNFIX = .FALSE.
C
C---- do not display viscous Q(s)
      LQSVIS = .FALSE.
C
C---- do not reverse s or Q direction in Q(s) plots
      LQSREV = .FALSE.
      LQQREV = .FALSE.
C
C---- maximum underrelaxation factor
      RLXMAX = 1.0
C
C---------------------------------------------------------------
C---- show point-strength and wake-strength values
      LPSHO = .TRUE.
      LWSHO = .TRUE.
C
C---- plot point-singularities, panel vertices, element numbers
      LPPLT = .TRUE.
      LVPLT = .TRUE.
      LEPLT = .TRUE.
C
C---- do not plot symbols on Cp or Q line plots
      LSPLT = .FALSE.
C
C---- plot displacement surfaces, don't use streamline integration to generate
      LDPLT = .TRUE.
      LDINT = .FALSE.
C
C---- plot grid on Cp or Q plots
      LPGRID = .TRUE.
C
C---- plot stagnation point on Cp or Q plots
      LPSTAG = .TRUE.
C
C---- plot grid on aero, airfoil geometry plots, ticks, parameters
      LAGRID = .TRUE.
      LGGRID = .TRUE.
      LGTICK = .TRUE.
      LGPARM = .TRUE.
C
      LGEOPLX= .FALSE.  ! for xfoil plot routines
C
C---- do not print neg. rpm Beta and Alfa output in local coordinates
      LCOORD = .FALSE.
C
C-----a case file has not been loaded
      LLOAD  = .FALSE.
C
C---------------------------------------------------------------
C---- Plotting flag
      IDEV = 1   ! X11 window only
c     IDEV = 2   ! B&W PostScript output file only (no color)
c     IDEV = 3   ! both X11 and B&W PostScript file
c     IDEV = 4   ! Color PostScript output file only 
c     IDEV = 5   ! both X11 and Color PostScript file 
C
C---- Re-plotting flag (for hardcopy)
C      IDEVPS = 2   ! B&W PostScript
      IDEVPS = 4   ! Color PostScript
C
C---- PostScript output logical unit and file specification
      IPSLU = 0  ! output to file  plot.ps   on LU 4    (default case)
c     IPSLU = ?  ! output to file  plot?.ps  on LU 10+?
C
C---- screen fraction taken up by plot window upon opening
C---- portrait window cannot be taller than desktop!
C---- sized for laptops... 
      SCRNFR  = 0.8
      PORTFAC = 1.1   ! Scales portrait windows (times SCRNFR)
      PAR = 0.6       ! plot aspect ratio (for xrotor subs)
C
C---- Default plot size in inches
C-    (Default plot window is 11.0 x 8.5)
C-   (Must not be larger than XPAGE if objects are to fit on paper page)
      SIZE = 10.0
C
C---- plot-window dimensions in inches for plot blowup calculations
C-    currently,  11.0 x 8.5  default window is hard-wired in libPlt
      XPAGE = 11.0
      YPAGE = 8.5
C
C---- page margins in inches, inset from XPAGE,YPAGE
      XMARG = 0.0
      YMARG = 0.0
C
C---- default input mode is via cursor
      LCURS = .TRUE.
C
C---- grid for spline user input plots
      LGRID = .TRUE.
C
C---- default X window orientation is landscape
      LLAND = .TRUE.
C
C---- no legend by default
      LEGEND = .FALSE.
C
C---- set flag to plot Cp on both sides, set color for each geometry element
      DO IEL = 1, NEX
        ISPLOT(IEL) = 0
CCC     ICOLEL(IEL) = 3 + MOD(IEL-1,8)
C
C------ reserve magenta for non-element plotting stuff, 
C-       and also leave out violet since it is too close to magenta
        ICOLEL(IEL) = 3 + MOD(IEL-1,6)
      ENDDO
C
C   Default-colormap color indices...
C
C     3  red
C     4  orange
C     5  yellow
C     6  green
C     7  cyan
C     8  blue
C     9  violet
C    10  magenta
C
C---------------------------------------------------------------
C---- baseline character size, plot-symbol size
      CHGT = 0.013   
      SHGT = 0.004
      CSIZE = CHGT   ! CH in xfoil code
C
C---- Cp(x) and Cp(y) plot aspect ratios
      CPXAR = 0.50
      CPYAR = 0.65
C
C---- axis limits for x-y plots
      CPLMAX =  1.0
      CPLMIN = -2.0
      CPLDEL = -0.5
C
      QPLMAX = 1.6
      QPLMIN = 0.0
      QPLDEL = 0.2
C
      XPLMAX =  1.0
      XPLMIN =  0.0
      XPLDEL =  0.1
C
      YPLMAX =  1.0
      YPLMIN =  0.0
      YPLDEL =  0.1
C
C---- default geometry box (set to 0,0,0,0 to disable boxing on plots)
      XGBOX(1) = 0.0
      YGBOX(1) = 0.0
      XGBOX(2) = 0.0
      YGBOX(2) = 0.0
C
C---- Plot type set to 0
      IPTYPE = 0
C---- scaling factors for Cp-vector and velocity-vector plots
      PVFAC = 0.2
      QVFAC = 0.1
C
C---- symbol_size/SH for each type of element
      SHF(0) = 0.20  ! solid surface
      SHF(1) = 0.25  ! wake (source+vortex specified)
      SHF(2) = 0.25  ! axis line (specified)
      SHF(3) = 0.50  ! line or ring (specified)
      SHF(4) = 0.75  ! point (specified)
      SHF(5) = 0.20  ! source line
      SHF(6) = 0.20  ! rotor source line
      SHF(7) = 0.20  ! vortex wake line
C
      SSIZEL(0) = SHGT*SHF(0)
      SSIZEL(1) = SHGT*SHF(1)
      SSIZEL(2) = SHGT*SHF(2)
      SSIZEL(3) = SHGT*SHF(3)
      SSIZEL(4) = SHGT*SHF(4)
      SSIZEL(5) = SHGT*SHF(5)
      SSIZEL(6) = SHGT*SHF(6)
      SSIZEL(7) = SHGT*SHF(7)
C
C---------------------------------------------------------------
C---- wake up Xplot11 routines
      CALL PLINITIALIZE
      CALL COLORMAPDEFAULT
C
C---- set up color spectrum
      NCOLOR = 64
      CALL COLORSPECTRUMHUES(NCOLOR,'RYGCBM')
C
C---- no active plot yet
      LPLOT = .FALSE.
C
      RETURN
      END ! DFINIT





      SUBROUTINE GENGEOM
C----------------------------------------------------
C     Generates case geometry from buffer airfoils 
C     and does the necessary geometry processing.
C----------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL ERROR, LQUERY
C
C---- process current buffer geometry to identify element types and 
C     set element parameters
c      write(80,*) 'before geproc'
c      CALL XYBDMP
c      CALL GEPROC
c      IF(LDBG) THEN
c        CALL GEPRINT(LUNDBG)
c      ENDIF
C
C========================================================================
C---- Check for respacing flag to redistribute points on buffer geometry
      IF(LRSPC) THEN 
C
C---- If paneling parameters are unset then reset default paneling 
C     parameters for buffer geometry elements
c      write(*,*) 'gengeom  npan',NPAN(1)
        IF(.NOT.LRSPCDEF) THEN
c          WRITE(*,*) 'Setting default respacing parameters'
          CALL PANDEF(0)
        ENDIF
C
c        WRITE(*,*) 'Repaneling buffer geometry'
        LQUERY = .FALSE.
        CALL PANGEN(LQUERY)
C
      ELSE
c        WRITE(*,*) 'Copying buffer geometry into panel arrays'
        CALL PANCOP
      ENDIF
C
C----------------------------------------------------------------
C---- At this point the input airfoils are in the panel geometry
      LGSAME = .TRUE.
C
C---- Adjust paneling on airfoils for rotor and wake grid
      CALL ADJPANL
c      write(80,*) 'after ADJPANL PDMP'
c      CALL XYPDMP(80)
C
C---- Reset # of elements to buffer geometry, remaining elements are added
C     to this: drag elements, rotor line sources, then vortex wakes
      NEL   = NBEL
C
C---- Set up drag area element(s)
      CALL SETDRGOBJ
C
C---- Set up rotor points
      CALL ROTORINIT
C
C---- Set up rotor wake grid from repaneled airfoils 
      CALL INIGRD
C
C---- Set up rotor wake elements from wake grid
      CALL SETROTWAK
C
C---- Sets up control points and pointers to panels
      CALL CVPGEN
C
C---- Initialize rotor/wake pointers 
      CALL ROTPINIT
C
C---- List all defined elements
      IF(LDBG) THEN
        CALL PELIST(6)
        CALL PELIST(LUNDBG)
      ENDIF
      CALL WAKEBOX
C
C----------------------------------------------------------------
C---- invalidate any existing solution
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
C
      NQSP = 0
      NSEG = 0
C
      WRITE(*,*)
      WRITE(*,*) 'Geometry loaded'
C
      RETURN
      END ! GENGEOM


      SUBROUTINE PELIST(LU)
C------------------------------------------------------
C     Lists elements defined in panel geometry
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- print out element info
      DO IEL = 1, NEL
        WRITE(LU,1001)
        WRITE(LU,1005) IEL
C
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        I = IP1
        N = IP2 - IP1 + 1
C       
        WRITE(LU,1010) N
C
        IF    (NETYPE(IEL).EQ.0) THEN
          IF(LBODY(IEL)) THEN
            IF(LXBOD(IEL)) THEN
             WRITE(LU,2000) 'Closed body on axis'
            ELSE
             WRITE(LU,2000) 'Closed body'
            ENDIF
          ELSE
            WRITE(LU,2000) 'Open surface'
          ENDIF
C
        ELSEIF(NETYPE(IEL).EQ.1) THEN
          WRITE(LU,2000) 'Prescribed wake surface'
C
        ELSEIF(NETYPE(IEL).EQ.2) THEN
          WRITE(LU,2000) 'Prescribed line singularity on axis'
C
        ELSEIF(NETYPE(IEL).EQ.3) THEN
          WRITE(LU,2000) 'Prescribed ring singularity'
C
        ELSEIF(NETYPE(IEL).EQ.4) THEN
          WRITE(LU,2000) 'Prescribed point singularity'
C
        ELSEIF(NETYPE(IEL).EQ.5) THEN
          WRITE(LU,2000) 'Source-only line singularity'
C
        ELSEIF(NETYPE(IEL).EQ.6) THEN
          WRITE(LU,2000) 'Rotor source-only line singularity'
C
        ELSEIF(NETYPE(IEL).EQ.7) THEN
          WRITE(LU,2000) 'Rotor vortex wake line'
C
        ENDIF
C
        IF(LTPAN(IEL)) THEN
          WRITE(LU,2000) 'TE-gap panel will be used'
        ENDIF
C
        IF(IEWAKE(IEL).NE.0) THEN
          WRITE(LU,2100) 'Prescribed-wake element', IEWAKE(IEL),
     &                  '  is attached'
        ENDIF
C
        IF( NETYPE(IEL).EQ.0 .OR.
     &      NETYPE(IEL).EQ.1 .OR.
     &      NETYPE(IEL).EQ.2      ) THEN
          WRITE(LU,1050) XPLE(IEL), YPLE(IEL),
     &                  XPTE(IEL), YPTE(IEL)
        ELSE
          WRITE(LU,1060) 'First point at',XP(IP1), YP(IP1)
          WRITE(LU,1060) 'Last  point at',XP(IP2), YP(IP2)
        ENDIF
C
      ENDDO
C
      WRITE(LU,1001)
C
      RETURN
C...............................................................
 1001 FORMAT(/' -----------------------------------------------')
 1005 FORMAT( ' Element', I3,' ...')
 1010 FORMAT( '   Number of input coordinate points:', I4)
 1020 FORMAT( '     delta(x) =', F10.5,'   Scale =', F10.5,
     &       /'     delta(y) =', F10.5,'   Angle =', F10.5,' deg')
 1050 FORMAT( '    LE x,y  =', 2F10.5 
     &       /'    TE x,y  =', 2F10.5 )
 1060 FORMAT( ' ',A,' x,y  =', 2F10.5 )
 2000 FORMAT( '   ', A)
 2100 FORMAT( '   ', A, I3, A)
      END ! PELIST



      SUBROUTINE GEPROC
C----------------------------------------------------
C     Processes buffer airfoil geometry to set
C     element type, various limits and pointers.
C----------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LYZERO
C
C---- process all elements
      DO IEL=1, NBEL
        CALL ELPROC(IEL)
C
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL MINMAX(N,XB(I),XBMINE(IEL),XBMAXE(IEL))
        CALL MINMAX(N,YB(I),YBMINE(IEL),YBMAXE(IEL))
      ENDDO
C
C---- set overall min,max over all elements
      CALL MINMAX(NBEL,XBMINE,XBMIN,DUMMY)
      CALL MINMAX(NBEL,YBMINE,YBMIN,DUMMY)
      CALL MINMAX(NBEL,XBMAXE,DUMMY,XBMAX)
      CALL MINMAX(NBEL,YBMAXE,DUMMY,YBMAX)
C
      DO IEL = 1, NBEL
C------ default reference point for each element
        XBREFE(IEL) = 0.
        YBREFE(IEL) = 0.
C------ clear element movement accumulators
        DXBSUM(IEL) = 0.
        DYBSUM(IEL) = 0.
        AGBSUM(IEL) = 0.
        XFBSUM(IEL) = 1.0
        YFBSUM(IEL) = 1.0
C------ first assume this element doesn't have a wake
        IEWAKE(IEL) = 0
      ENDDO
C
C---- set element-type indicators and wake pointers
      DO IEL = 1, NBEL
        IB1 = IBFRST(IEL)
        IB2 = IBLAST(IEL)
C
C------ is this element entirely on the y=0 line?
        LYZERO = .TRUE.
        DO IB = IB1, IB2
          IF(YB(IB) .NE. 0.0) LYZERO = .FALSE.
        ENDDO
C
        IF(IB1.EQ.IB2) THEN
C------- only one point defines element...
         IF(LYZERO) THEN
C-------- axisymm: point singularity on axis
          NETYPE(IEL) = 4
         ELSE
C-------- axisymm: ring singularity
          NETYPE(IEL) = 3
         ENDIF
C
        ELSEIF(LYZERO) THEN
C------- line element on axis
         NETYPE(IEL) = 2
C
C------- check for source line element defined in buffer geometry (old scheme)
        ELSEIF(NBTYPE(IEL).EQ.5) THEN
           NETYPE(IEL) = NBTYPE(IEL)
C
        ELSE
C------- assume this element is a solid surface
           NETYPE(IEL) = 0
C
C------- Now check for wake elements 
C------- a wake element is attached to TE of some element before it
         DO KEL = 1, IEL-1
C--------- check only solid boundarys (NETYPE=0) for attached wakes
           IF(NETYPE(KEL).EQ.0) THEN
             KB1 = IBFRST(KEL)
             KB2 = IBLAST(KEL)
             DB11SQ = (XB(IB1)-XB(KB1))**2 + (YB(IB1)-YB(KB1))**2
             DB21SQ = (XB(IB2)-XB(KB1))**2 + (YB(IB2)-YB(KB1))**2
             DB12SQ = (XB(IB1)-XB(KB2))**2 + (YB(IB1)-YB(KB2))**2
             DB22SQ = (XB(IB2)-XB(KB2))**2 + (YB(IB2)-YB(KB2))**2
             DB1TSQ = (XB(IB1)-XBTE(KEL))**2 + (YB(IB1)-YBTE(KEL))**2
             DB2TSQ = (XB(IB2)-XBTE(KEL))**2 + (YB(IB2)-YBTE(KEL))**2
C
             DSBSQ = MIN( DB11SQ,
     &                    DB12SQ,
     &                    DB21SQ,
     &                    DB22SQ,
     &                    DB1TSQ,
     &                    DB2TSQ)
             DSBMIN = SQRT(DSBSQ)
C
C----------- this element IEL is KEL's wake if it starts at KEL's TE point
             DSBTOL = 0.0001 * ABS(SB(KB2) - SB(KB1))
             IF(DSBMIN.LT.DSBTOL) THEN
              NETYPE(IEL) = 1
              IEWAKE(KEL) = IEL
              GO TO 25
             ENDIF
           ENDIF
         ENDDO
C
 25      CONTINUE
C
        ENDIF
C
      ENDDO
C
C---- set wake-termination box outline
cc      CALL WAKEBOX
C
      RETURN
      END ! GEPROC


      SUBROUTINE GEPRINT(LU)
C----------------------------------------------------
C     Processes buffer airfoil geometry to set
C     various limits and pointers.
C----------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- print out element info
      DO IEL = 1, NBEL
        WRITE(LU,1001)
        WRITE(LU,1005) IEL
C
        IB1 = IBFRST(IEL)
        IB2 = IBLAST(IEL)
        I = IB1
        N = IB2 - IB1 + 1
C       
        WRITE(LU,1010) N
C
        IF(AREA2DA(IEL).GE.0.0) THEN
         WRITE(LU,2000) 'Counterclockwise ordering'
        ELSE
         WRITE(LU,2000) 'Clockwise ordering'
        ENDIF
C
        IF    (NETYPE(IEL).EQ.0) THEN
          IF(LBODY(IEL)) THEN
            IF(LXBOD(IEL)) THEN
             WRITE(LU,2000) 'Closed body on axis'
            ELSE
             WRITE(LU,2000) 'Closed body'
            ENDIF
          ELSE
            WRITE(LU,2000) 'Open surface'
          ENDIF
C
        ELSEIF(NETYPE(IEL).EQ.1) THEN
          WRITE(LU,2000) 'Prescribed wake surface'
C
        ELSEIF(NETYPE(IEL).EQ.2) THEN
          WRITE(LU,2000) 'Prescribed line singularity on axis'
C
        ELSEIF(NETYPE(IEL).EQ.3) THEN
          WRITE(LU,2000) 'Prescribed ring singularity'
C
        ELSEIF(NETYPE(IEL).EQ.4) THEN
          WRITE(LU,2000) 'Prescribed point singularity'
C
        ENDIF
C
        IF(LTPAN(IEL)) THEN
          WRITE(LU,2000) 'TE-gap panel will be used'
        ENDIF
C
        IF(IEWAKE(IEL).NE.0) THEN
          WRITE(LU,2100) 'Prescribed-wake element', IEWAKE(IEL),
     &                  '  is attached'
        ENDIF
C
        IF( NETYPE(IEL).EQ.0 .OR.
     &      NETYPE(IEL).EQ.1 .OR.
     &      NETYPE(IEL).EQ.2      ) THEN
          WRITE(LU,1050) XBLE(IEL), YBLE(IEL),
     &                   XBTE(IEL), YBTE(IEL)
        ELSE
          WRITE(LU,1060) XB(IB1), YB(IB1)
        ENDIF
C
      ENDDO
C
      WRITE(LU,1001)
C
      RETURN
C
C...............................................................
 1001 FORMAT(/' -----------------------------------------------')
 1005 FORMAT( ' Element', I3,' ...')
 1010 FORMAT( '   Number of input coordinate points:', I4)
 1020 FORMAT( '     delta(x) =', F10.5,'   Scale =', F10.5,
     &       /'     delta(y) =', F10.5,'   Angle =', F10.5,' deg')
 1050 FORMAT( '    LE x,y  =', 2F10.5 
     &       /'    TE x,y  =', 2F10.5 )
 1060 FORMAT( '       x,y  =', 2F10.5 )
 2000 FORMAT( '   ', A)
 2100 FORMAT( '   ', A, I3, A)
      END ! GEPRINT



      SUBROUTINE ELPROC(IEL)
C------------------------------------------------------
C     Processes buffer geometry element IEL.
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LTMP, LCUT
C
      IB1 = IBFRST(IEL)
      IB2 = IBLAST(IEL)
C
      I = IB1
      N = IB2 - IB1 + 1
      CALL SCALC(XB(I),YB(I),SB(I),N)
C
C----- check for body that crosses the axis (defined by upper and lower sides)
      LCUT = .FALSE.
      DO IB=IB1, IB2-1
        IF(YB(IB).GT.0.0 .AND. YB(IB+1).LE.0.0) THEN
          LCUT = .TRUE.
          IBCUT = IB+1
          GO TO 50
        ENDIF
      END DO
C
 50   IF(LCUT) THEN
C----- element crosses axis, cut it at the axis and use only upper half 
       IF(LDBG) THEN
         WRITE(*,*) 'Cutting element ', IEL, ' at X axis '
         WRITE(LUNDBG,*) 'Cutting element ', IEL, ' at X axis '
       ENDIF
       CALL SEGSPL(XB(I),XBS(I),SB(I),N)
       CALL SEGSPL(YB(I),YBS(I),SB(I),N)
       SS = SB(IBCUT)
       CALL SINVRT(SS,0.0,YB(IB1),YBS(IB1),SB(IB1),N)
       XB(IBCUT) = SEVAL(SS,XB(IB1),XBS(IB1),SB(IB1),N)
       YB(IBCUT) = 0.0
       IBLAST(IEL) = IBCUT
C
C----- recalculate indices for element and arc length
       IB1 = IBFRST(IEL)
       IB2 = IBLAST(IEL)
       I = IB1
       N = IB2 - IB1 + 1
       CALL SCALC(XB(I),YB(I),SB(I),N)
cc       DO IB = IB1, IB2
cc         write(50,*) IB, XB(IB), YB(IB)
cc       END DO
      ENDIF      
C
C
      IF(N.EQ.1) THEN
C---- single point element (doublet/source point)
       LV1ZR(IEL) = .FALSE.
       LV2ZR(IEL) = .FALSE.
       LXBOD(IEL) = .FALSE.
       LBODY(IEL) = .FALSE.
       LTPAN(IEL) = .FALSE.
       LQNZR(IEL) = .FALSE.
C
      ELSE
C---- multiple point element (doublet/source line or vortex/source panels)
       SBLEN = ABS(SB(IB2)-SB(IB1))
C
C---- are first,last points effectively on the axis?
       LV1ZR(IEL) = ABS(YB(IB1)) .LT. SBLEN*DTEPAN
       LV2ZR(IEL) = ABS(YB(IB2)) .LT. SBLEN*DTEPAN
C
       IF(LV1ZR(IEL)) THEN
C---- TE gap is really a TE disk at point IB2
        DSBTE = ABS(YB(IB2))
       ELSE
C---- usual TE gap
        DSBTE = SQRT((XB(IB1)-XB(IB2))**2 + (YB(IB1)-YB(IB2))**2)
       ENDIF
C
C---- axisymmetric body closed on the axis ?
       LXBOD(IEL) = LV1ZR(IEL) .AND. LV2ZR(IEL)
C
C---- finite-thickness body ?
       LBODY(IEL) = LV1ZR(IEL) .OR. LV2ZR(IEL)
     &         .OR. DSBTE .LT. ABS(SB(IB2)-SB(IB1))*DTEBOD
C
C---- add TE-gap panel ?
       LTPAN(IEL) = LBODY(IEL)
     &        .AND. DSBTE .GT. ABS(SB(IB2)-SB(IB1))*DTEPAN
C
C---- explicitly zero out the body normal velocity QNDOF? 
cc     LQNZR(IEL) = .TRUE.
       LQNZR(IEL) = .FALSE.
C
      ENDIF
C
C---- 2D buffer geometry data
C---- set 2D area, centroid, section moduli - for area
      CALL AECALC(N,XB(I),YB(I),ONE, 0, PERIM,
     &            AREA2DA(IEL),XBCEN2DA(IEL),YBCEN2DA(IEL),
     &            EIXX2DA(IEL),EIYY2DA(IEL),EIXY2DA(IEL) )
C---- set 2D area, centroid, section moduli - per unit skin thickness
      CALL AECALC(N,XB(I),YB(I),ONE, 2, PERIM,
     &            AREA2DT(IEL),XBCEN2DT(IEL),YBCEN2DT(IEL),
     &            EIXX2DT(IEL),EIYY2DT(IEL),EIXY2DT(IEL) )
C
C---- Axisymmetric buffer geometry data
C---- set axisymmetric surface area, volume, centroid, radii of gyration
      CALL AXCALC(N,XB(I),YB(I),ONE, 0, 
     &            VOLUMV(IEL), ASURFV(IEL),XBCENV(IEL),YBCENV(IEL),
     &            RGYRXV(IEL),RGYRYV(IEL) )
C
C---- assume element will be used and written in input order by default
      LREVEL(IEL) = .FALSE.
C
      IF(LBODY(IEL) .AND. AREA2DA(IEL).LT.0.0) THEN
C----- area is negative (clockwise order)... mark for reversing the order
       LREVEL(IEL) = .TRUE.
      ENDIF
C
      IF(LREVEL(IEL)) THEN
C---- reverse the order of the points 
       DO IB=IB1, (IB2+IB1)/2
         XTMP = XB(IB2-IB+IB1)
         YTMP = YB(IB2-IB+IB1)
         XB(IB2-IB+IB1) = XB(IB)
         YB(IB2-IB+IB1) = YB(IB)
         XB(IB) = XTMP
         YB(IB) = YTMP
       ENDDO
       LTMP = LV1ZR(IEL)
       LV1ZR(IEL) = LV2ZR(IEL)
       LV2ZR(IEL) = LTMP
      ENDIF
C
      CALL SCALC(XB(I),YB(I),SB(I),N)
      CALL SEGSPL(XB(I),XBS(I),SB(I),N)
      CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
      IF(LBODY(IEL) .AND. .NOT.LXBOD(IEL)) THEN
C----- body off axis
       IF(LV2ZR(IEL)) THEN
         SBLE(IEL) = SB(IB2) 
         XBLE(IEL) = XB(IB2)
         YBLE(IEL) = YB(IB2)
         XBTE(IEL) = XB(IB1)
         YBTE(IEL) = YB(IB1)
        ELSE
         CALL LEFIND(SBLE(IEL),XB(I),XBS(I),YB(I),YBS(I),SB(I),N)
         XBLE(IEL) = SEVAL(SBLE(IEL),XB(I),XBS(I),SB(I),N)
         YBLE(IEL) = SEVAL(SBLE(IEL),YB(I),YBS(I),SB(I),N)
         XBTE(IEL) = 0.5*(XB(IB1)+XB(IB2))
         YBTE(IEL) = 0.5*(YB(IB1)+YB(IB2))
       ENDIF
C
      ELSE
C----- other body... set LE,TE from leftmost and rightmost points
       IF(XB(IB1).LT.XB(IB2)) THEN
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
       ENDIF
C
      ENDIF
C
C---- set Cp-side flag
      IF(LBODY(IEL)) THEN
C----- outer surface only (to right side of +s CCW direction)
       ISPLOT(IEL) = +1
      ELSE
C----- default... plot both sides
       ISPLOT(IEL) = 0
      ENDIF
C
      RETURN
      END ! ELPROC



      SUBROUTINE WAKEBOX
C---------------------------------------------------
C     Sets limits of "wake box".
C     A wake is terminated when it leaves this box.
C---------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LWBSET
C
      XWBOX(1) =  1.0E23
      XWBOX(2) = -1.0E23
      YWBOX(1) =  1.0E23
      YWBOX(2) = -1.0E23
      LWBSET = .FALSE.
      DO IEL = 1, NEL
        IF(NETYPE(IEL).EQ.0 .OR. NETYPE(IEL).EQ.7) THEN
C------- box defined only by solid or vortex wake elements
         XWBOX(1) = MIN(XWBOX(1),XPMINE(IEL))
         XWBOX(2) = MAX(XWBOX(2),XPMAXE(IEL))
         YWBOX(1) = MIN(YWBOX(1),YPMINE(IEL))
         YWBOX(2) = MAX(YWBOX(2),YPMAXE(IEL))
         LWBSET = .TRUE.
        ENDIF
      ENDDO
      IF(.NOT.LWBSET) THEN
C----- use default if no solid elements exist
       XWBOX(1) = XBMIN
       XWBOX(2) = XBMAX
       YWBOX(1) = YBMIN
       YWBOX(2) = YBMAX
      ENDIF
C
cc      DBOX = MAX( XWBOX(2)-XWBOX(1) , YWBOX(2)-YWBOX(1) )
cc      XWBOX(1) = XWBOX(1) - DBOX
cc      XWBOX(2) = XWBOX(2) + DBOX
cc      YWBOX(1) = YWBOX(1) - DBOX
cc      YWBOX(2) = YWBOX(2) + DBOX
C
      DXBOX = XWBOX(2)-XWBOX(1)
      DYBOX = YWBOX(2)-YWBOX(1)
      XWBOX(1) = XWBOX(1) - 0.5*DXBOX
      XWBOX(2) = XWBOX(2) - 0.02*DXBOX 
      YWBOX(1) = 0.0
      YWBOX(2) = YWBOX(2) + DYBOX
C
      RETURN
      END ! WAKEBOX


      SUBROUTINE SETDRGOBJ
C------------------------------------------------------
C     Sets up drag area elements into paneled geometry 
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(NDOBJ.LT.1) THEN
        RETURN
      ENDIF
C
      NPTOT = IPLAST(NEL)
C
C---- Use rotor points to define paneled element for rotor source line
      DO ND = 1, NDOBJ
        CALL ADDELEM5(XDDEF(1,ND),YDDEF(1,ND),NDDEF(ND))
        IELDRGOBJ(ND) = NEL
      END DO
      LDRGOBJ = .TRUE.
C
      RETURN
      END


      SUBROUTINE ROTORINIT
C------------------------------------------------------
C     Sets up rotor geometry from panel geometry 
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION YRC0(IRX),BGAM0(IRX),CHR0(IRX),BETAR0(IRX)
      DIMENSION CLDES0(IRX)
      DIMENSION T1(IRX),T1S(IRX),
     &          T2(IRX),T2S(IRX),
     &          T3(IRX),T3S(IRX),
     &          T4(IRX),T4S(IRX)
C
      IF(NROTOR.EQ.0) THEN
        WRITE(*,*) 'No actuator disk/blades defined'
        RETURN
      ENDIF
C
      NRCSAV = NRC
      IF(NRCSAV.GT.0) THEN
        N = 1   ! assume at least one disk defined
        DO IR = 1, NRCSAV
          YRC0(IR)   = YRC(IR,N)
          CLDES0(IR) = CLDES(IR)
        END DO
      ENDIF
C
      DO N = 1, NROTOR
C
C---- get radial location of rotor on CB and duct wall
      IF(IPROTCB(N).EQ.0) THEN
        WRITE(*,*) 'Error locating disk on CB wall'
        STOP
      ENDIF
      IPCB = IPROTCB(N)
      XCB = XP(IPCB)
      YCB = YP(IPCB)
      RHUB(N) = YCB
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'Disk on CB @IP= ',IPCB
        WRITE(LUNDBG,*) 'Xhub,Rhub on CB ', XCB,RHUB(N)
      ENDIF
C      
      IF(IPROTDW(N).EQ.0) THEN
        WRITE(*,*) 'Error locating rotor on Duct wall'
        STOP
      ENDIF
      IPDW = IPROTDW(N)
      XDW = XP(IPDW)
      YDW = YP(IPDW)
      RTIP(N) = YDW
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'Disk on Duct @IP= ',IPDW
        WRITE(LUNDBG,*) 'Xtip,Rtip on Duct wall ', XDW,RTIP(N)
      ENDIF
C
C---- Update swept disk area
      ADISK(N) = PI*(RTIP(N)**2 - RHUB(N)**2)
      ATIP(N)  = PI*(RTIP(N)**2)
      IF(LDBG) THEN
        WRITE(*,*) 'Disk # ',N,' swept area ',ADISK(N)
      ENDIF
C
C
      IF(IRTYPE(N).NE.0) THEN   ! disk defined
C---- Save old blade and circulation distribution for interpolation to new pts
        NRC0 = NRCSAV
        DO IR = 1, NRC0
          YRC0(IR)   = YRC(IR,N)
          BGAM0(IR)  = BGAM(IR,N)
          CHR0(IR)   = CHR(IR,N)
          BETAR0(IR) = BETAR(IR,N)
c          IF(N.EQ.1) CLDES0(IR) = CLDES(IR)
        END DO
      ENDIF
C
C---- Set discrete points on rotor line using linear interpolation hub-tip
      NRP = NRSTA
      DX = (XDW-XCB)/FLOAT(NRP-1)
      DY = (YDW-YCB)/FLOAT(NRP-1)
      DO IR = 1, NRP
        XRP(IR,N) = XCB + FLOAT(IR-1)*DX
        YRP(IR,N) = YCB + FLOAT(IR-1)*DY
cc        write(*,*) 'rotorpts ',i,xrp(i,n),yrp(i,n)
      END DO
C
C
C *** In what follows, condition on LLOFT is to allow paneling to geometric
C     tip gap for lofting, while changing nothing for OPER.
C
C---- If effective tip gap is specified add two points to last interval 
C     at tip gap from duct wall, shift original rotor points inwards
C
      IF(LLOFT) THEN
         TGAPZLT = 0.0
      ELSE
         TGAPZLT = TGAPZL
      ENDIF
C
      TG = MAX(0.0,TGAP - TGAPZLT)
      IF(TG.GT.0.0) THEN
        DX = XDW - XCB
        DR = RTIP(N) - RHUB(N)
        DO IR = 1, NRP
          YRP(IR,N) = YRP(IR,N) - 2.0*TG*(YRP(IR,N)-RHUB(N))/DR
          XRP(IR,N) = XRP(IR,N) -     DX*(RTIP(N)-YRP(IR,N))/DR 
cc          write(*,*) 'rotorpts ',i,xrp(i,n),yrp(i,n)
        END DO
C--- Add two points, one at RTIP-TGAP, one at shroud wall at RTIP
        NRP = NRP + 1
        YRP(NRP,N) = RTIP(N) - TG
        XRP(NRP,N) = XDW  - DX*(RTIP(N)-YRP(NRP,N))/DR 
        NRP = NRP + 1
        XRP(NRP,N) = XDW 
        YRP(NRP,N) = YDW
      ENDIF
C
C---- Rotor interval centers and initial circulation
      NRC = NRP - 1
      DO IR = 1, NRC
        XRC(IR,N)  = 0.5*(XRP(IR,N)+XRP(IR+1,N))
        YRC(IR,N)  = 0.5*(YRP(IR,N)+YRP(IR+1,N))
        BGAM(IR,N) = 0.0
cc        write(*,*) 'rotorcpts ',ir,xrc(ir,n),yrc(ir,n)
      END DO
C
C=====================================================================
C
      IF(IRTYPE(N).EQ.0) THEN   ! disk undefined
C
C---- If actuator disk has been input interpolate BGAM to rotor center points
        IF(IRTYPDEF(N).EQ.1 .AND. NRDEF(N).GE.2) THEN
          DO I = 1, NRDEF(N)
            T1(I) = BGAMDEF(I,N)
          END DO
C---- Spline rotor definition arrays
          CALL SEGSPL(T1,T1S,YRDEF(1,N),NRDEF(N))
C---- Set BGAM for points on rotor line
          DO IR = 1, NRC
            BGAM(IR,N) = SEVAL(YRC(IR,N),T1,T1S,YRDEF(1,N),NRDEF(N))
          END DO
          IF(TGAP.GT.0.0) BGAM(NRC,N) = 0.0
          IRTYPE(N) = 1
          LBLDEF = .FALSE.
          LVMAV = .FALSE.
        ENDIF
C---- If a rotor has been input (Y,CH,BETA) interpolate CH,BETA to rotor 
C     center points
        IF(IRTYPDEF(N).EQ.2 .AND. NRDEF(N).GE.2) THEN
          DO I = 1, NRDEF(N)
            T1(I) = CHRDEF(I,N)
            T2(I) = BETADEF(I,N)
          END DO
C---- Spline rotor definition arrays
          CALL SEGSPL(T1,T1S,YRDEF(1,N),NRDEF(N))
          CALL SEGSPL(T2,T2S,YRDEF(1,N),NRDEF(N))
C---- Set CH,BETA for points on rotor line
          DO IR = 1, NRC
            CHR(IR,N)   = SEVAL(YRC(IR,N),T1,T1S,YRDEF(1,N),NRDEF(N))
            BETAR(IR,N) = SEVAL(YRC(IR,N),T2,T2S,YRDEF(1,N),NRDEF(N))
          END DO
          IRTYPE(N) = 2
          LBLDEF = .TRUE.
          LVMAV = .FALSE.
          CALL SETIAERO
        ENDIF
C
      ELSE
C
C---- Interpolate old rotor/circulation to new rotor points
        DO I = 1, NRC0
          T1(I) = CHR0(I)
          T2(I) = BETAR0(I)
          T3(I) = BGAM0(I)
        END DO
C---- Spline rotor definition arrays
        CALL SEGSPL(T1,T1S,YRC0,NRC0)
        CALL SEGSPL(T2,T2S,YRC0,NRC0)
        CALL SEGSPL(T3,T3S,YRC0,NRC0)
C---- Set CH,BETA,BGAM for points on new rotor line
        DO IR = 1, NRC
          CHR(IR,N)   = SEVAL(YRC(IR,N),T1,T1S,YRC0,NRC0)
          BETAR(IR,N) = SEVAL(YRC(IR,N),T2,T2S,YRC0,NRC0)
          BGAM(IR,N)  = SEVAL(YRC(IR,N),T3,T3S,YRC0,NRC0)
        END DO
ccc        IF(TGAP.GT.0.0) BGAM(NRC,N) = 0.0
        CALL SETIAERO
C
      ENDIF
C
      END DO ! end NROTOR loop
C
C---- Move design CLs and BB data to new radial stations
C
      IF(NRCSAV.GT.0) THEN
        CLP1 = CLPOS(1)
        CLP2 = CLPOS(NRCSAV)
        CLN1 = CLNEG(1)
        CLN2 = CLNEG(NRCSAV)
        CLD1 = CLDES(1)
        CLD2 = CLDES(NRCSAV)
C
        DO IR = 1, NRC
          CLPOS(IR) = CLP1 + FLOAT(IR-1)/FLOAT(NRC-1) * (CLP2-CLP1)
          CLNEG(IR) = CLN1 + FLOAT(IR-1)/FLOAT(NRC-1) * (CLN2-CLN1)
          CLDES(IR) = CLD1 + FLOAT(IR-1)/FLOAT(NRC-1) * (CLD2-CLD1)
        ENDDO
C
        DO N = 1, NROTOR
          IF(LBBLOFT(N)) THEN
            DO IC = 1,NRCSAV
              T4(IC) = BBVFAC(IC,N)
            ENDDO
            CALL SEGSPL(T4,T4S,YRC0,NRCSAV)
            DO IC = 1, NRC
              BBVFAC(IC,N) = SEVAL(YRC(IC,N),T4,T4S,YRC0,NRCSAV)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
C        DO I = 1, NRCSAV
C          T4(I) = CLPOS(I)
C        END DO
C        CALL SEGSPL(T4,T4S,YRC0,NRCSAV)
C        DO I = 1, NRC
C          CLPOS(I) = SEVAL(YRC(I,1),T4,T4S,YRC0,NRCSAV)
C        END DO
C
C        DO I = 1, NRCSAV
C          T4(I) = CLNEG(I)
C        END DO
C        CALL SEGSPL(T4,T4S,YRC0,NRCSAV)
C        DO I = 1, NRC
C          CLNEG(I) = SEVAL(YRC(I,2),T4,T4S,YRC0,NRCSAV)
C        END DO
C
C        DO I = 1, NRCSAV
C          T4(I) = CLDES(I)
C        END DO
C        CALL SEGSPL(T4,T4S,YRC0,NRCSAV)
C        DO I = 1, NRC
C          CLDES(I) = SEVAL(YRC(I,1),T4,T4S,YRC0,NRCSAV)
C        END DO
C      ENDIF
C
      RETURN
      END ! ROTORINIT



      SUBROUTINE SETROTWAK
C------------------------------------------------------
C     Sets up rotor source line and wake elements from wake grid
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION Y1(IRX)
      DIMENSION T1(IRX),T1S(IRX),
     &          T2(IRX),T2S(IRX),
     &          T3(IRX),T3S(IRX)
C
C---- Add rotor source element after elements now in panel geometry 
C     (foils, drag area objects)
cc      write(*,*) 'SETROTWAK entered NEL ',NEL
C
      DO NR = 1, NROTOR
C
C---- Define paneled source line element for rotor lines from grid points
       IG = IGROTOR(NR)
       DO IR = 1, NRP
        XRP(IR,NR) = XG(IG,IR)
        YRP(IR,NR) = YG(IG,IR)
       END DO
       CALL ADDELEM6(XRP(1,NR),YRP(1,NR),NRP)
       IELROTOR(NR) = NEL
C
C---- Interpolate old rotor params/circulation to new rotor points
       NRC = NRP - 1
       DO IR = 1, NRC
         Y1(IR) = YRC(IR,NR)
         T1(IR) = CHR(IR,NR)
         T2(IR) = BETAR(IR,NR)
         T3(IR) = BGAM(IR,NR)
       END DO
C---- Redefine rotor center points
       DO IR = 1, NRC
         XRC(IR,NR)  = 0.5*(XRP(IR,NR)+XRP(IR+1,NR))
         YRC(IR,NR)  = 0.5*(YRP(IR,NR)+YRP(IR+1,NR))
       END DO 
C---- Spline rotor definition arrays
       CALL SEGSPL(T1,T1S,Y1,NRC)
       CALL SEGSPL(T2,T2S,Y1,NRC)
       CALL SEGSPL(T3,T3S,Y1,NRC)
C---- Set CH,BETA,BGAM for points on new rotor line
        DO IR = 1, NRC
          CHR(IR,NR)   = SEVAL(YRC(IR,NR),T1,T1S,Y1,NRC)
          BETAR(IR,NR) = SEVAL(YRC(IR,NR),T2,T2S,Y1,NRC)
          BGAM(IR,NR)  = SEVAL(YRC(IR,NR),T3,T3S,Y1,NRC)
        END DO
ccc        IF(TGAP.GT.0.0) BGAM(NRC,NR) = 0.0
        CALL SETIAERO
C
      END DO
C
C---- Add vortex wake elements to panel geometry 
C
C----- set up wake element from centerbody TE (element 1)
      IR = 1
      IEL = 1
      IG1 = IGTECB
      N = II - IG1 + 1
      CALL ADDWAKE(XG(IG1,IR),YG(IG1,IR),N)
C---- save wake element index for this rotor point
      IR2IEL(IR)  = NEL
      IEL2IR(NEL) = IR
C---- plot Cp only on interior flow side
      ISPLOT(NEL) = 2
      IEWAKE(IEL) = NEL
C---- set up pointers from wake points to grid IG streamwise index
      IP1 = IPFRST(NEL)
      IP2 = IPLAST(NEL)
      DO IP = IP1, IP2
        IP2IG(IP) = IG1 + IP-IP1 
      END DO
C
C----- set up intermediate rotor wakes (from rotor line)
      DO IR = 2, NRP-1
        IG1 = 1
        N = II
        CALL ADDWAKE(XG(IG1,IR),YG(IG1,IR),N)
C---- save wake element index for this rotor point
        IR2IEL(IR)  = NEL
        IEL2IR(NEL) = IR
C---- plot Cp only on interior flow side
        ISPLOT(NEL) = -1
C---- set up pointers from wake points to grid IG streamwise index
        IP1 = IPFRST(NEL)
        IP2 = IPLAST(NEL)
        DO IP = IP1, IP2
          IP2IG(IP) = IG1 + IP-IP1 
        END DO
      END DO

C----- and wake from duct TE (element 2)
      IR = NRP
      IEL = 2
      IG1 = IGTEDW
      N = II - IG1 + 1
      CALL ADDWAKE(XG(IG1,IR),YG(IG1,IR),N)
C---- save wake element index for this rotor point
      IR2IEL(IR)  = NEL
      IEL2IR(NEL) = IR
C---- plot Cp only on interior flow side
      ISPLOT(NEL) = 1
      IEWAKE(IEL) = NEL
C---- set up pointers from wake points to grid IG streamwise index
      IP1 = IPFRST(NEL)
      IP2 = IPLAST(NEL)
      DO IP = IP1, IP2
        IP2IG(IP) = IG1 + IP-IP1 
      END DO
C
C----------------------------------------------------------------
C---- invalidate any existing solution
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
C
      RETURN
      END ! SETROTWAK


      SUBROUTINE ADDWAKE(X,Y,N)
C---------------------------------------------------------
C     Add wake specified by N points X,Y to panel geometry 
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION X(*),Y(*)
C
C---- add wake to panel geometry
      IP1 = IPLAST(NEL) + 1
      IP2 = IP1 - 1 + N
      IEL = NEL + 1
      IF(LDBG) THEN
        WRITE(LUNDBG,10) IEL,N,IP1,IP2
      ENDIF
C---- add wake points 
      DO IP = IP1, IP2
        I = IP - IP1 + 1
        XP(IP) = X(I)
        YP(IP) = Y(I)
        IF(LDBG) WRITE(LUNDBG,*) IP,XP(IP),YP(IP) 
      END DO
C---- centroid location
      XPCENT(IEL) = 0.5*(XP(IP1) + XP(IP2))
      YPCENT(IEL) = 0.5*(YP(IP1) + YP(IP2))
C---- panel limits
      CALL MINMAX(N,XP(IP1),XPMINE(IEL),XPMAXE(IEL))
      CALL MINMAX(N,YP(IP1),YPMINE(IEL),YPMAXE(IEL))
C
C---- setup pointers and panel data
      IPFRST(IEL) = IP1
      IPLAST(IEL) = IP2
      NETYPE(IEL) = 7
      LTPAN(IEL) = .TRUE.
c      LTPAN(IEL) = .FALSE.
      IPTE2(IEL) = IP2
      IPTE1(IEL) = 0
      CALL XYPSPL(IEL)
c      CALL XYCSET(IEL)
C
      NEL   = IEL
      NPTOT = IP2
C
      IF(LDBG) THEN
        WRITE(*,10) IEL,N,IP1,IP2
      ENDIF
C
 10   FORMAT('Vortex wake element ',I4,' added with #pts ',I5,
     &       ' IP1 ',I5,' IP2 ',I5)
C
      RETURN
      END ! ADDWAKE



      SUBROUTINE UPDROTWAK
C------------------------------------------------------
C     Updates up rotor wake element geometry from wake grid
C     Assumes that points have moved but indices remain the same
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
cc      write(*,*) 'UPDROTWAK entered NEL ',NEL
C
C----- update wake element from centerbody TE (element 1)
      IR = 1
      IEL = IR2IEL(IR)
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      I1  = IGTECB
      N   = II - I1 + 1
      NIP = IP2 - IP1 + 1
      IF(N.NE.NIP) THEN 
        WRITE(*,*) 'UPDROTWAK wake grid mismatch on IR ',IR,N,NIP
        STOP
      ENDIF
      I = I1
      DO IP = IP1, IP2
        XP(IP) = XG(I,IR)
        YP(IP) = YG(I,IR)
        I = I + 1
      END DO
C
C----- set up intermediate rotor wakes (from rotor line)
      DO IR = 2, NRP-1
        IEL = IR2IEL(IR)
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        I1  = 1
        N   = II
        NIP = IP2 - IP1 + 1
        IF(N.NE.NIP) THEN
        WRITE(*,*) 'UPDROTWAK wake grid mismatch on IR ',IR,N,NIP
          STOP
        ENDIF
        I = I1
        DO IP = IP1, IP2
          XP(IP) = XG(I,IR)
          YP(IP) = YG(I,IR)
          I = I + 1
        END DO
      END DO

C----- and wake from duct TE (element 2)
      IR = NRP
      IEL = IR2IEL(IR)
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      I1  = IGTEDW
      N   = II - I1 + 1
      NIP = IP2 - IP1 + 1
      IF(N.NE.NIP) THEN
        WRITE(*,*) 'UPDROTWAK wake grid mismatch on IR ',IR,N,NIP
        STOP
      ENDIF
      I = I1
      DO IP = IP1, IP2
        XP(IP) = XG(I,IR)
        YP(IP) = YG(I,IR)
        I = I + 1
      END DO
C
C----------------------------------------------------------------
C---- invalidate any existing solution
      LNCVP = .FALSE.
      LQAIC = .FALSE.
      LQGIC = .FALSE.
      LQCNT = .FALSE.
      LGSYS = .FALSE.
      LGAMU = .FALSE.
      LGAMA = .FALSE.
C
      RETURN
      END ! UPDROTWAK


      SUBROUTINE ROTPINIT
C------------------------------------------------------
C     Sets up pointers for panel points on rotor wakes
C     and pointers from panel centers to streamline grid
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Clear pointers to streamline rows from body and wake points
C     Used for additional gamma from vortex wakes
      DO IP = 1, IPX
        IP2IR(IP) = 0
      END DO
C
C---- Clear pointers to grid centers from body and wake centers
      DO IC = 1, ICX
        IC2IG(IC) = 0
      END DO
C
      NR = 1
C---- Pointers for vortex wake on centerbody
      IR = 1
      IEL = 1       
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      IPCB = IPROTCB(NR)
cc    IF(LDBG) WRITE(22,*) IR, IEL, IP1, IP2, IPCB
      DO IP = IPCB, IP1, -1
        IP2IR(IP) = IR
cc      IF(LDBG) WRITE(22,*) IP,IR
      END DO
C---- set pointers from panel center to grid center
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      DO IC = IC1, IC2
        IG1 = IP2IG(IPCO(IC))
        IG2 = IP2IG(IPCP(IC))
        IC2IG(IC) = MIN(IG1,IG2)
cc        IF(IC2IG(IC).NE.0) WRITE(22,*) IC,IC2IG(IC)
      END DO
C
C---- Set pointers for vortex wake elements
      DO IR = 1, NRP
        IEL = IR2IEL(IR)       
        IF(IEL.LE.0) THEN
          WRITE(*,*) 'IR2IEL out of range in ROTPINIT ',IEL
          STOP
        ENDIF
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        DO IP = IP1, IP2
          IP2IR(IP) = IR
        END DO
C---- set pointers from panel center to grid center
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        DO IC = IC1, IC2
          IG1 = IP2IG(IPCO(IC))
          IG2 = IP2IG(IPCP(IC))
          IC2IG(IC) = MIN(IG1,IG2)
cc          IF(IC2IG(IC).NE.0) WRITE(22,*) IC,IC2IG(IC)
        END DO
      END DO
C
C---- Pointers for vortex wake on duct 
      IR = NRP
      IEL = 2       
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      N = IP2 - IP1 + 1
      IPDW = IPROTDW(NR)
cc    IF(LDBG) WRITE(22,*) IR, IEL, IP1, IP2, IPDW
      DO IP = IPDW, IP2
        IP2IR(IP) = IR
cc      IF(LDBG) WRITE(22,*) IP,IR        
      END DO
C---- set pointers from panel center to grid center
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      DO IC = IC1, IC2
        IG1 = IP2IG(IPCO(IC))
        IG2 = IP2IG(IPCP(IC))
        IC2IG(IC) = MIN(IG1,IG2)
cc        IF(IC2IG(IC).NE.0) WRITE(22,*) IC,IC2IG(IC)
      END DO
C
      RETURN
      END ! ROTPINIT



      SUBROUTINE ADDELEM5_2PT(X1,Y1,X2,Y2,NPTS)
C------------------------------------------------------
C     Add line source element (NETYPE=5) to paneled geometry 
C     Specified by X1,Y1 start and X2,Y2 end coordinates and 
C     # of points.  Intermediate ponts are linearly interpolated
C     between endpoints.
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(NPTS.LE.1) THEN 
        WRITE(*,*) 'Error in ADDELEM12 NPTS<2'
        STOP
      ENDIF
C
      IP = IPLAST(NEL)
      IEL = NEL
C
      IP1 = IP + 1
      IP2 = IP + NPTS
C
      DX = X2 - X1
      DY = Y2 - Y1
      FNPT1 = FLOAT(NPTS-1)
      DO I = 1, NPTS
        IP = IP + 1
        XP(IP) = X1 + DX*FLOAT(I-1)/FNPT1
        YP(IP) = Y1 + DY*FLOAT(I-1)/FNPT1
      END DO
C---- centroid location
      XPCENT(IEL) = 0.5*(XP(IP1) + XP(IP2))
      YPCENT(IEL) = 0.5*(YP(IP1) + YP(IP2))
C---- panel limits
      CALL MINMAX(N,XP(IP1),XPMINE(IEL),XPMAXE(IEL))
      CALL MINMAX(N,YP(IP1),YPMINE(IEL),YPMAXE(IEL))
C
C--- add element 
      IEL = IEL + 1
      IPFRST(IEL) = IP1
      IPLAST(IEL) = IP2
      NETYPE(IEL) = 5
      LTPAN(IEL) = .FALSE.
      CALL XYPSPL(IEL)
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*) ' '
        WRITE(LUNDBG,*) 'Adding line source element ',IEL
        WRITE(LUNDBG,100) 'Start ',IP1,X1,Y1
        WRITE(LUNDBG,100) 'End   ',IP2,X2,Y2
C
        WRITE(*,10) IEL,NPTS,IP1,IP2
      ENDIF
C
 10   FORMAT('Line source element ',I4,' added with #pts ',I5,
     &       ' IP1 ',I5,' IP2 ',I5)
 100  FORMAT(A,I5,4(1X,G13.5))
C
      NEL   = IEL
      NPTOT = IP2
C
      RETURN
      END ! ADDELEM5_2PT



      SUBROUTINE ADDELEM5(X,Y,N)
C------------------------------------------------------
C     Add line source element (NETYPE=5) to paneled geometry 
C     Specified by X,Y coordinates and # of points
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION X(*),Y(*)
C
      IF(N.LE.1) THEN 
        WRITE(*,*) 'Error in ADDELEM5 NPTS<2'
        STOP
      ENDIF
C
C---- add sourceline points to panel geometry
      IP1 = IPLAST(NEL) + 1
      IP2 = IP1 - 1 + N
      IEL = NEL + 1
      IF(LDBG) THEN
        WRITE(LUNDBG,10) IEL,N,IP1,IP2
      ENDIF
C
C---- add points to panel arrays 
      DO IP = IP1, IP2
        I = IP - IP1 + 1
        XP(IP) = X(I)
        YP(IP) = Y(I)
        IF(LDBG) WRITE(LUNDBG,*) IP,XP(IP),YP(IP) 
      END DO
C---- centroid location
      XPCENT(IEL) = 0.5*(XP(IP1) + XP(IP2))
      YPCENT(IEL) = 0.5*(YP(IP1) + YP(IP2))
C---- panel limits
      CALL MINMAX(N,XP(IP1),XPMINE(IEL),XPMAXE(IEL))
      CALL MINMAX(N,YP(IP1),YPMINE(IEL),YPMAXE(IEL))
C
C--- add element 
      IPFRST(IEL) = IP1
      IPLAST(IEL) = IP2
      NETYPE(IEL) = 5
      LTPAN(IEL) = .FALSE.
      CALL XYPSPL(IEL)
C
 10   FORMAT(/'Line source element ',I4,' added with #pts ',I5,
     &        ' IP1 ',I5,' IP2 ',I5)
 100  FORMAT(A,I5,4(1X,G13.5))
C
      NEL   = IEL
      NPTOT = IP2
C
      RETURN
      END ! ADDELEM5



      SUBROUTINE ADDELEM6(X,Y,N)
C---------------------------------------------------------
C     Add element for rotor source line to panel geometry.  
C     Specified by N points X,Y.
C     Defined as TYPE=6 element
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION X(*),Y(*)
C
C---- add wake to panel geometry
      IP1 = IPLAST(NEL) + 1
      IP2 = IP1 - 1 + N
      IEL = NEL + 1
      IF(LDBG) THEN
        WRITE(LUNDBG,10) IEL,N,IP1,IP2
      ENDIF
C
C---- add points to panel arrays 
      DO IP = IP1, IP2
        I = IP - IP1 + 1
        XP(IP) = X(I)
        YP(IP) = Y(I)
        IF(LDBG) WRITE(LUNDBG,*) IP,XP(IP),YP(IP) 
      END DO
C---- centroid location
      XPCENT(IEL) = 0.5*(XP(IP1) + XP(IP2))
      YPCENT(IEL) = 0.5*(YP(IP1) + YP(IP2))
C---- panel limits
      CALL MINMAX(N,XP(IP1),XPMINE(IEL),XPMAXE(IEL))
      CALL MINMAX(N,YP(IP1),YPMINE(IEL),YPMAXE(IEL))
C
C---- setup pointers and panel data
      IPFRST(IEL) = IP1
      IPLAST(IEL) = IP2
      NETYPE(IEL) = 6
      LTPAN(IEL) = .FALSE.
      CALL XYPSPL(IEL)
c      CALL XYCSET(IEL)
C
      NEL   = IEL
      NPTOT = IP2
C
 10   FORMAT('Rotor source line element ',I4,' added with #pts ',I5,
     &       ' IP1 ',I5,' IP2 ',I5)
C
      RETURN
      END ! ADDELEM6




      SUBROUTINE DFLOAD(FNAMIN,FERROR)
C---------------------------------------------------------
C     Reads previously saved duct case data and geometry 
C     from file FNAMIN in DFDC version 0.5 format.  This
C     format separates data into groups by keywords.
C     Once read the geometry is processed from buffer to 
C     paneled geometry.
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*)   FNAMIN
      CHARACTER*128   FNAME
      CHARACTER*128   IFILE, LINE, FILECB, FILEDCT
      LOGICAL LOPEN, LREPANEL
      DIMENSION RINPUT(10)
C
C---- local arrays for calling AREAD
      PARAMETER (ITX=300)
      DIMENSION XT(ITX,NEX), YT(ITX,NEX)
      DIMENSION NT(NEX)
C
      LOGICAL ERROR, FERROR
C
C---- reset some flags for new dataset
      FERROR = .FALSE.
      NROTOR = 0
      LBLDEF = .FALSE.
      LCONV  = .FALSE.
      NRP       = 0
      NRC       = 0
      DO N = 1, NRX
       IRTYPE(N) = 0
       NRDEF(N)    = 0
       IRTYPDEF(N) = 0
      END DO
C
      LRSPCDEF = .FALSE.
      LRSPCUSR = .FALSE.
C
      NDOBJ       = 0
C
      LOPEN = .FALSE.
      LU = 1
      FNAME = FNAMIN
      CALL STRIP(FNAME,NF)
      IF(FNAME.NE.' ') THEN
        OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=98)
        LOPEN = .TRUE.
       ELSE
        RETURN
      ENDIF
C
      ICNT = 0
C
C---- Read header line and version from DFDC input file
      CALL RDLINE(LU,LINE,ICNT)
      IF(LINE.EQ.'END' .OR. LINE.EQ.'ERR') GO TO 210
C
      IF(LINE(1:4).NE.'DFDC') THEN
        WRITE(*,*) 'Not DFDC input file...may be duct geometry file'
        CLOSE(LU)
cc        CALL LOADG(FNAME)
        RETURN
      ENDIF
      IV = INDEX(LINE,'Version') 
      READ(LINE(IV+8:LEN(LINE)),*) FILEVERS
      WRITE(*,1005) FILEVERS
C
 1000 FORMAT(A)
 1005 FORMAT(/,' Reading file from DFDC Version ',E8.2)
 1010 FORMAT(' File  ',A,' not found'/)
 1020 FORMAT(' File  ',A,' has incompatible format'/
     &       ' Loading not completed'/)
 1025 FORMAT(' Case Name: ',A)
C
C
C--- Get Case title from line #2
      CALL RDLINE(LU,LINE,ICNT)
      IF(LINE.EQ.'END' .OR. LINE.EQ.'ERR') GO TO 210
      NAME = LINE
      CALL STRIP(NAME,NNAME)
C
      WRITE(*,1025) NAME
C
C-----------------------------------------------------
C--- Main loop to find keywords and start read actions 
C    for data groups
 10   CALL RDLINE(LU,LINE,ICNT)
      IF(LINE.EQ.'END') THEN
        GO TO 110
      ELSEIF(LINE.EQ.'ERR') THEN 
        GO TO 210
C
C
C---- OPER data
      ELSEIF(INDEX(LINE,'OPER').NE.0) THEN
        CALL RDLINE(LU,LINE,ICNT)
cc        READ(LINE,*,ERR=210) QINF,QREF,RPM
          NINPUT = 4
          CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
          IF(NINPUT.GE.4) THEN
            QINF = RINPUT(1)            
            QREF = RINPUT(2)            
            RPM1 = RINPUT(3)            
            RPM2 = RINPUT(4)            
          ELSEIF(NINPUT.GE.3) THEN
            QINF = RINPUT(1)            
            QREF = RINPUT(2)            
            RPM1 = RINPUT(3)
            RPM2 = 0.0            
          ELSEIF(NINPUT.GE.2) THEN
            QINF = RINPUT(1)            
            QREF = RINPUT(2)            
            RPM1 = 0.0
            RPM2 = 0.0            
          ENDIF
C
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) RHO,VSO,RMU,ALTH
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) XDWKLEN,NWAKE
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) LWRLX
        CALL RDLINE(LU,LINE,ICNT)
        IF(INDEX(LINE,'ENDOPER').NE.0) THEN
          WRITE(*,1030)QINF,QREF,RPM1,RPM2,RHO,VSO,RMU,ALTH,
     &                 XDWKLEN,NWAKE
        ELSE
          WRITE(*,*) 'No ENDOPER to OPER section'
          STOP
        ENDIF
C
 1030   FORMAT( ' OPER data read from file',/,
     &          '  Qinf   = ',F9.3,5X,'Qref   = ',F9.3,/,
     &          '  RPM1   = ',F9.1,5X,'RPM2   = ',F9.1,/,
     &          '  Rho    = ',F9.5,5X,'VSound = ',F9.3,/,
     &          '  rMU    = ',E9.3,5X,'Alt    = ',F9.5,/,
     &          '  Xdwake = ',F9.5,5X,'Nwake  = ',I6)
C
C---- Input RPM's set speed of disks
        OMEGA(1) = PI*RPM1/30.0
        OMEGA(2) = PI*RPM2/30.0
        IF(RHO.EQ.0.0 .OR. VSO.EQ.0.0 .OR. RMU.EQ.0.0) THEN
           DELTAT = 0.0
           CALL ATMO(ALTH,DELTAT,VSO,RHO,RMU) 
        ENDIF
C
C
C---- AERO data
      ELSEIF(INDEX(LINE,'AERO').NE.0) THEN
C--- Read aero section definitions
        NR = NROTOR + 1
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) NAERO(NR)
        DO N = 1, NAERO(NR)
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) XISECT
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) A0DEG,DCLDA,CLMAX,CLMIN
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) DCLDA_STALL,DCL_STALL,CMCON,MCRIT
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) CDMIN,CLDMIN,DCDCL2
          CALL RDLINE(LU,LINE,ICNT)
cc          READ(LINE,*,ERR=210) REREF,REXP
          NINPUT = 4
          CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
          IF(NINPUT.GE.4) THEN
            REREF = RINPUT(1)            
            REXP  = RINPUT(2)            
            TOC   = RINPUT(3)            
            DCDCL2S = RINPUT(4)            
          ELSEIF(NINPUT.GE.3) THEN
            REREF = RINPUT(1)            

            REXP  = RINPUT(2)            
            TOC   = RINPUT(3)            
            DCDCL2S = 0.020
          ELSEIF(NINPUT.GE.2) THEN
            REREF = RINPUT(1)            
            REXP  = RINPUT(2)            
            TOC     = 0.1
            DCDCL2S = 0.020
          ENDIF
          A0 = A0DEG *DTR
          CALL PUTAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &                 DCLDA,DCLDA_STALL,DCL_STALL,
     &                 CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &                 CMCON,MCRIT,TOC,REREF,REXP)
        END DO
        CALL RDLINE(LU,LINE,ICNT)
        IF(INDEX(LINE,'ENDAERO').NE.0) THEN
          WRITE(*,2010) NR,NAERO(NR)
        ELSE
          WRITE(*,*) 'No ENDAERO to AERO section'
          STOP
        ENDIF
C
 2010 FORMAT(/,' AERO data read for Disk',I2,': ',I3,' sections')
C
C
C---- ROTOR data
      ELSEIF(INDEX(LINE,'ROTOR').NE.0) THEN
        NR = NROTOR + 1
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) XDISK(NR), NRBLD(NR), NRSTA
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) NRDEF(NR)
C--- Get rotor definition data
        DO IR = 1, NRDEF(NR)
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,END=210,ERR=210) YRDEF(IR,NR),CHRDEF(IR,NR),
     &                                 BETADEF(IR,NR)
        END DO
C
 13     CALL RDLINE(LU,LINE,ICNT)
        IF(INDEX(LINE,'ENDROTOR').NE.0) THEN
          WRITE(*,2020) XDISK(NR), NRBLD(NR), NRSTA
          IF(NRDEF(NR).GE.2) THEN
           WRITE(*,2030) NRDEF(NR)
           DO IR = 1, NRDEF(NR)
             WRITE(*,15) YRDEF(IR,NR),CHRDEF(IR,NR),BETADEF(IR,NR)
             BETADEF(IR,NR) = BETADEF(IR,NR)*DTR
           END DO
 15        FORMAT(1X,F12.6,2X,F12.6,4X,F12.6)
           IRTYPDEF(NR) = 2
           NROTOR = NR
          ELSE
           WRITE(*,*) 'Rotor blade defined by too few stations'
          ENDIF
        ELSE
          WRITE(*,*) 'No ENDROTOR found'
          STOP
        ENDIF
C
 2020 FORMAT(/,
     & ' ROTOR data read from file',/,
     & '  XDisk: ',F5.3,'   Blades:',I2,'   Rotor pnts:',I3)
 2030 FORMAT(
     & '  Rotor blade defined with',I3,' points',/,
     & '         R            CH            Beta'  )
C
C
C---- ACTDISK data
      ELSEIF(INDEX(LINE,'ACTDISK').NE.0) THEN
        NR = NROTOR + 1
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) XDISK(NR), NRSTA
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) NRDEF(NR)
C--- Get actuator disk definition data
        DO IR = 1, NRDEF(NR)
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,END=210,ERR=210) YRDEF(IR,NR),BGAMDEF(IR,NR)
        END DO   
C
 18     CALL RDLINE(LU,LINE,ICNT)
        IF(INDEX(LINE,'ENDACTDISK').NE.0) THEN
          WRITE(*,2040) XDISK(NR), NRSTA
          IF(NRDEF(NR).GE.2) THEN
           WRITE(*,2050) NRDEF(NR)
           DO IR = 1, NRDEF(NR)
             WRITE(*,15) YRDEF(IR,NR),BGAMDEF(IR,NR)
           END DO
           IRTYPDEF(NR) = 1
           NROTOR = NR
          ELSE
           WRITE(*,*) 'Actuator disk defined by too few stations'
          ENDIF
        ELSE
          WRITE(*,*) 'No ENDACTDISK found'
          STOP
        ENDIF
C
 2040 FORMAT(/,
     & ' ACTDISK data read from file',/,
     & '  XDisk: ',F5.3,'   Rotor points:',I3)
 2050 FORMAT(
     & '  Actuator Disk defined with',I3,' points',/,
     & '         R           BGAM'  )
C
C
C---- DRAGOBJ data
      ELSEIF(INDEX(LINE,'DRAGOBJ').NE.0) THEN
        ND1 = NDOBJ + 1
        IF(ND1.GT.NDRGX) THEN
          WRITE(*,*) 'Number of drag objects exceeds dimension NDGX'
        ELSE
         CALL RDLINE(LU,LINE,ICNT)
         READ(LINE,*,ERR=210) NDDEF(ND1)
C--- Get drag object definition data
         DO ID = 1, NDDEF(ND1)
           CALL RDLINE(LU,LINE,ICNT)
           READ(LINE,*,END=210,ERR=210) XDDEF(ID,ND1),YDDEF(ID,ND1),
     &                                  CDADEF(ID,ND1)
         END DO   
C
 19      CALL RDLINE(LU,LINE,ICNT)
         IF(INDEX(LINE,'ENDDRAGOBJ').NE.0) THEN
           NDOBJ = ND1
           WRITE(*,2060) ND1
           IF(NDDEF(ND1).GE.2) THEN
            WRITE(*,2070) NDDEF(ND1)
            DO ID = 1, NDDEF(ND1)
              WRITE(*,15) XDDEF(ID,ND1),YDDEF(ID,ND1),CDADEF(ID,ND1)
            END DO
           ELSE
            WRITE(*,*) 'Drag object defined by too few stations'
           ENDIF
         ELSE
           WRITE(*,*) 'No ENDDRAGOBJ found'
           STOP
         ENDIF ! ENDDRAGOBJ check
C
        ENDIF ! NDRGX check
C
 2060 FORMAT(/,' DRAGOBJ ',I2,'  read from file')
 2070 FORMAT(
     & '  Drag object defined with',I3,' points',/,
     & '         X             R              CDA')
C
C
C---- GEOM data (CB and duct coordinates in a single XFOIL multi-element file)
      ELSEIF(INDEX(LINE,'GEOM').NE.0) THEN
C
C---- Read the combined CB and duct airfoil file
C
        CALL AREADNR(LU, ITX,NEX, XT,YT,
     &               NT,NTEL,
     &               ANAME,ISPARS,IFTYPE)
        IF(IFTYPE.EQ.0 .OR. NTEL.NE.2) THEN
C----- read error occurred for two elements
         WRITE(*,2080) NTEL
         STOP
        ENDIF
C
 2080 FORMAT(
     $ ' Error reading GEOM data',/,
     $ ' Coordinates read for',I2,' foils')
C
C---- User-defined respacing specs
C---- PANE data (respacing data for surface points
      ELSEIF(INDEX(LINE,'PANE').NE.0) THEN
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) NBIN, NRPNIN
        NRPNRD = MIN(NRPNIN,NRPNX)
        DO IEL = 1, NBIN
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) NPAN(IEL)
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) CVEX(IEL), SMOF(IEL), 
     &                         FSLE(IEL), FSTE(IEL)
          DO IRPN = 1, NRPNRD
            CALL RDLINE(LU,LINE,ICNT)
            READ(LINE,*,ERR=210) SRPN1(IRPN,IEL),
     &                           SRPN2(IRPN,IEL),
     &                           CRRAT(IRPN,IEL)
          ENDDO
        ENDDO
        CALL RDLINE(LU,LINE,ICNT) 
        IF(INDEX(LINE,'ENDPANE').NE.0) THEN
          WRITE(*,2090) NBIN
          LRSPCDEF = .TRUE.
          LRSPCUSR = .TRUE.
        ENDIF
C
      ENDIF
      GO TO 10
C
 2090 FORMAT( ' PANELING for', I2, ' elements read from file' )
C
C------------------------------------------------------------------------
 110  CONTINUE
      NBEL = NTEL
C
C---- Process airfoil inputs for CB and duct
C---- put element-indexed coordinates XT,YT into linear bufffer airfoil 
C     geometry arrays XB,YB
      IBNEXT = 1
      DO IEL=1, NBEL
C------ first point in element IEL
        IBFRST(IEL) = IBNEXT
C------ initalize accumulation counter
        IB = IBFRST(IEL) - 1
C------ go over all points on element IEL
        IF(LDBG) WRITE(LUNDBG,*) 'Airfoil element'
        DO IT = 1, NT(IEL)
          IB = IB+1
          IF(IB.GT.IBX) STOP 'LOAD: Array overflow on IBX'
          XB(IB) = XT(IT,IEL)
          YB(IB) = YT(IT,IEL)
          IF(LDBG) WRITE(LUNDBG,*) IT,XT(IT,IEL),YT(IT,IEL)
        ENDDO
C---- set buffer airfoils to NBTYPE=0  
        NBTYPE(IEL) = 0
        IBLAST(IEL) = IB
        IBNEXT = IB + 1
      END DO
C---- set total number of points in buffer geometry
      NBTOT = IBLAST(NBEL)
C
C
C=======================================================================
C---- If panel spacing data is not in case file try reading explicitly
C     from xxx.pan paneling input file
      IF(.NOT.LRSPCDEF) THEN
C---- get file prefix (if it has a .xxx suffix)
        K = INDEX(FNAME,'.') - 1
        IF(K.LE.0) THEN
C----- no "." prefix/suffix separator 
C      just tack on ".pan" to full filename
         CALL STRIP(FNAME,NF)
         K = NF
        ENDIF
        PREFIX = FNAME(1:K)
        NPREFIX = K
        PFILE = PREFIX(1:NPREFIX) // '.pan'
        CALL PANGET(PFILE,ERROR)
C
        IF(ERROR) THEN
         WRITE(*,2100) PFILE(1:K+4)
        ELSE
         WRITE(*,2110) PFILE(1:K+4)
        ENDIF
      ENDIF
C
 2100 FORMAT( ' No repaneling parameters from file ',A)
 2110 FORMAT( ' Repaneling parameters from file ',A)

C
C=======================================================================
C---- process current buffer geometry to identify element types and 
C     set element parameters
      CALL GEPROC
      IF(LDBG) THEN
        CALL GEPRINT(LUNDBG)
      ENDIF
C
C---- process buffer foil geometry and set up paneled geometry and grid
cc    CALL GENGEOM
C
C---- take note of new case
      LLOAD = .TRUE.
      DO I=1,NROTOR
        LBBLOFT(I) = .FALSE.
      ENDDO
      GO TO 300
C...............................................................
   98 CONTINUE
      WRITE(*,1050) FNAME(1:NF)
      FERROR = .TRUE.
      GO TO 300
C
   99 CONTINUE
      WRITE(*,1100) FNAME(1:NF)
      FERROR = .TRUE.
      GO TO 300
C
 210  CONTINUE
      WRITE(*,1150) ICNT,FNAME(1:NF)
      FERROR = .TRUE.
C
 1050 FORMAT(/' File OPEN error:  ', A)
 1100 FORMAT(/' File READ error:  ', A)
 1150 FORMAT(/' File READ error on line ',I3,':  ', A)
C
C---- close file
 300  IF(LOPEN) CLOSE(LU)
      RETURN
      END ! DFLOAD



      SUBROUTINE DFSAVE(FNAMEIN)
C--------------------------------------------------------------------------
C     Save rotor and operating state in DFDC Version 0.5 format
C     This format saves data in groups delimited by keywords 
C     to separate data
C--------------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*) FNAMEIN
      CHARACTER ANS*1, FNAME*128, FBLNK*80
C
      LU = 2
      FNAME = FNAMEIN
      CALL STRIP(FNAME,NF)
C
      IF(FNAME(1:1) .EQ. ' ') THEN
        CALL ASKS('Enter filename^',FNAME)
      ENDIF
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=5)
      WRITE(*,*)
      WRITE(*,*) 'Output file exists.  Overwrite?  Y'
      READ (*,1000) ANS
      IF(INDEX('Nn',ANS).EQ.0) GO TO 6
C
      CLOSE(LU)
      WRITE(*,*) 'Current duct case not saved.'
      RETURN
C
 5    OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=90)
 6    REWIND(LU)
C
C      
C--- Version header and case name
      IF(NAME.EQ.' ') NAME = 'Saved ducted fan'
      WRITE(LU,1100) VERSION,NAME
C
C--- OPER data
C
      IF(NROTOR.LT.2) OMEGA(2)=0.
C
C--- Velocity, reference velocity and RPM
      WRITE(LU,1102) 
      RPM1 = 30.0*OMEGA(1)/PI
      RPM2 = 30.0*OMEGA(2)/PI
      WRITE(LU,1200) QINF,QREF,RPM1,RPM2
C--- Altitude and atmospheric data
      WRITE(LU,1103) 
      WRITE(LU,1200) RHO,VSO,RMU,ALTH
C--- XDwake, #wake points
      WRITE(LU,1104)
      WRITE(LU,1202) XDWKLEN, NWAKE
C--- Wake relaxation flag
      WRITE(LU,1105)
      WRITE(LU,1204) LWRLX
      WRITE(LU,1106)
C
C--- Save data for each disk (act disk or bladed)
      DO NR = 1, NROTOR
C
C--- AERO data
C--- Save aero data for defined aero sections
      WRITE(LU,1108) NAERO(NR)
      DO N = 1, NAERO(NR)
       CALL GETAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &              CMCON,MCRIT,TOC,REREF,REXP)
       WRITE(LU,1200) XISECT
       A0DEG = A0 *180.0/PI
       WRITE(LU,1109)
       WRITE(LU,1200) A0DEG,DCLDA,CLMAX,CLMIN
       WRITE(LU,1110) 
       WRITE(LU,1200) DCLDA_STALL,DCL_STALL,CMCON,MCRIT
       WRITE(LU,1111) 
       WRITE(LU,1200) CDMIN,CLDMIN,DCDCL2
       WRITE(LU,1112)
       WRITE(LU,1200) REREF,REXP,TOC,DCDCL2S
      END DO
      WRITE(LU,1113)
C
C--- Save rotor blade if defined
      IF(IRTYPE(NR).EQ.2) THEN
       WRITE(LU,1130)
C--- Rotor axial location, #blades, #radial stations
       WRITE(LU,1202) XDISK(NR),NRBLD(NR), NRSTA
C--- Save blade definition with chord,twist and body velocity
       WRITE(LU,1132) NRC
       DO I=1, NRC
         BETADEG = BETAR(I,NR)/DTR
         WRITE(LU,1200) YRC(I,NR),CHR(I,NR),BETADEG
       END DO
       WRITE(LU,1134)
      ENDIF
C
C--- Save actuator disk if defined
      IF(IRTYPE(NR).EQ.1) THEN
       WRITE(LU,1122)
C--- Rotor axial location, #radial stations
       WRITE(LU,1202) XDISK(NR),NRSTA
C--- Save actuator disk circulation
       WRITE(LU,1124) NRC
       DO I=1, NRC
         WRITE(LU,1200) YRC(I,NR),BGAM(I,NR)
       END DO
       WRITE(LU,1126)
      ENDIF
C
      END DO
C
C--- Save drag object(s) if defined
      IF(NDOBJ.GT.0) THEN
        DO N=1, NDOBJ
          IF(NDDEF(N).GT.2) THEN
            WRITE(LU,1136) NDDEF(N)
C--- Save drag object points and CDA values
            DO I=1, NDDEF(N)
              WRITE(LU,1200) XDDEF(I,N),YDDEF(I,N),CDADEF(I,N)
            END DO
           WRITE(LU,1137)
          ENDIF
        END DO
      ENDIF
C
C--- CB and duct geometry
      WRITE(LU,1140)
C---- save only 2 elements in buffer geometry into opened file
      NELSAV = 2
      IFTYPE = 2
      FBLNK = ' '
      CALL AWRITE(FBLNK,LU,
     &            NELSAV,IPFRST,IPLAST,XP,YP,
     &            ANAME, ISPARS, IFTYPE)
      WRITE(LU,1142)
C
C--- Save user-defined paneling specs
      IF(LRSPCUSR) THEN
       WRITE(LU,1150) NBEL, NRPNX
       DO IEL = 1, NBEL
         WRITE(LU,1152) NPAN(IEL),
     &                  CVEX(IEL), SMOF(IEL), FSLE(IEL), FSTE(IEL)
         WRITE(LU,1154)
         DO IRPN = 1, NRPNX
           WRITE(LU,1200) SRPN1(IRPN,IEL),
     &                    SRPN2(IRPN,IEL),
     &                    CRRAT(IRPN,IEL)
         ENDDO
       ENDDO
       WRITE(LU,1156)
      ENDIF
C
      CLOSE(LU)
      RETURN
C
 90   WRITE(*,*) 'Bad filename.'
      WRITE(*,*) 'Current duct case not saved.'
      RETURN
C 
C...................................................................
 1000 FORMAT(A)
 1100 FORMAT('DFDC Version ',E8.2,' '/A32)
C
 1102 FORMAT(/'OPER',
     &       /'!   Vinf         Vref         RPM1         RPM2')
 1103 FORMAT( '!   Rho          Vso          Rmu          Alt')
 1104 FORMAT( '!  XDwake             Nwake')
 1105 FORMAT( '!        Lwkrlx') 
 1106 FORMAT( 'ENDOPER') 
C
 1108 FORMAT(/'AERO',
     &       /'!  #sections'/1(1X,I5),
     &       /'!  Xisection')
 1109 FORMAT( '!  A0deg        dCLdA        CLmax        CLmin')
 1110 FORMAT( '! dCLdAstall   dCLstall      Cmconst      Mcrit')
 1111 FORMAT( '!   CDmin      CLCDmin       dCDdCL^2')
 1112 FORMAT( '!   REref        REexp        TOC         dCDdCL^2')
 1113 FORMAT( 'ENDAERO') 
C
 1122 FORMAT(/'ACTDISK',
     &       /'!  Xdisk       NRsta')
 1124 FORMAT( '!  #stations'/1(1X,I5),
     &       /'!           r         BGam')
 1126 FORMAT( 'ENDACTDISK') 
C
 1130 FORMAT(/'ROTOR',
     &       /'!  Xdisk               Nblds       NRsta')
 1132 FORMAT( '!  #stations'/1(1X,I5),
     &       /'!     r        Chord         Beta')
 1134 FORMAT( 'ENDROTOR') 
C
 1136 FORMAT(/'DRAGOBJ',
     &       /'!  #pts'/1(1X,I5),
     &       /'!     x            r            CDA')
 1137 FORMAT( 'ENDDRAGOBJ') 
C
 1140 FORMAT(/'GEOM')
 1142 FORMAT( 'ENDGEOM') 
C
 1150 FORMAT(/'PANELING',
     &       /'!  #elements   #refinement zones'/(2(1X,I5)))
 1152 FORMAT( '!  #panel nodes'/1(1X,I5)
     &       /'!  curv_expon  curv_smooth   dsL/dsAvg    dsR/dsAvg',
     &       /4(1X,G12.5))
 1154 FORMAT('!  s1/smax      s2/smax     ds/dsAvg')
 1156 FORMAT( 'ENDPANELING') 
C
 1200 FORMAT(5(1X,G12.5))
 1202 FORMAT(1(1X,G12.5),2(8X,I5))
 1204 FORMAT(12X,L1)
C
Cx123456789012x123456789012x123456789012x123456789012x123456789012
C!         Rho          Vso          Rmu           Alt')
C
      END ! DFSAVE



