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

      SUBROUTINE WAKERESET
C-------------------------------------------------------------
C     Realigns the upper and lower wakes bounding the duct
C     vortex wake grid with the current flow.  
C
C     A new wake grid system is defined by updating and re-relaxing
C     the grid.  All wake elements are re-initialized with the new
C     geometry.
C-------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- CB wake
      IEL = IR2IEL(1)
      CALL WAKMOV(IEL,ISPLOT(IEL))
C
C---- duct wake
      IEL = IR2IEL(NRP)
      CALL WAKMOV(IEL,ISPLOT(IEL))
C
C---- change grid boundaries, relax grid and update vortex wakes
      CALL UPDGRD
      CALL UPDROTWAK 
C
C---- invalidate existing solution
      LNCVP = .FALSE.
      LQAIC = .FALSE.
      LQGIC = .FALSE.
      LQCNT = .FALSE.
      LGSYS = .FALSE.
      LGAMU = .FALSE.
      LGAMA = .FALSE.
C
      RETURN
      END ! WAKERESET




      SUBROUTINE WAKMOV(IELA,ISQ)
C-------------------------------------------------------------
C     Moves existing wake (element IELA) by
C     aligning panels with current velocities.
C     Panel side flag ISQ indicates which velocity to use.
C       ISQ=0 uses mean surface velocities QC
C       ISQ=1 uses right surface velocities QCR
C       ISQ=2 uses left  surface velocities QCL
C
C     Panel points are moved preserving existing panel lengths 
C     while aligning the panels with the current velocity 
C     vectors on the panel.
C-------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(NETYPE(IELA).EQ.0) THEN
C----- this element is a body panel
       WRITE(*,*) 'WAKMOV: Element',IELA,' is body, cannot be moved!'
       RETURN
      ENDIF
C
      IEL = IELA
C
C---- Save X location for end of wake
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      XP2 = XP(IP2)
C
      CALL WAKMOVR(IEL,XP,YP,ISQ)
C
C---- Reset last point to initial (unchanged) XWAKE location
      XP(IP2) = XP2
C
C---- respline node points
      CALL XYPSPL(IEL)
C---- reset control points and normal vectors
      CALL XYCSET(IEL)
      CALL ANPSET(IEL)
C
      WRITE(*,*) 'Generating dq/d(gam,sig) influence of moved wake...'
      CALL QAIC1(.FALSE., 1,NC, IEL,IEL)
      DO KEL = 1, NEL
        ICE = NCX + KEL
        IF(LBODY(KEL)) CALL QAIC1(.FALSE.,ICE,ICE,IEL,IEL)
      ENDDO
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      CALL GSYS(.FALSE., .FALSE., .TRUE., IP1,IP2, 1,0)
C
C---- free existing unit wake strength vectors (if any), mark vectors invalid
      DO IP = IP1, IP2
        LUSET(IUGAM(IP)) = .FALSE.
        LUSET(IUSIG(IP)) = .FALSE.
        IUGAM(IP) = 0
        IUSIG(IP) = 0
      ENDDO
      DO IR = 1, NRP
        LUSET(IUVWK(IR)) = .FALSE.
        IUVWK(IR) = 0
      ENDDO
C
      RETURN
      END ! WAKMOV



      SUBROUTINE WAKMOVR(IEL,XPNEW,YPNEW,ISQ)
C----------------------------------------------------------
C---- Routine to move points on element to new positions
C     to align with velocity vector on side ISQ (0/1/2)
C----------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION XPNEW(*), YPNEW(*)
C
      DIMENSION DSP(IPX)
C
C---- speed below this is treated as stagnation
      QSTAG = 0.0001*QREF
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
C
C---- calculate velocity and its Jacobians at all wake points
ccc      CALL QAIC1(.FALSE., IC1,IC2, 1,NEL)
C---- set velocities on both sides of panels at all control points
ccc      CALL QQCALC(NC,ANC, IPCO,IPCP, GAM,SIG, QC, QCL,QCR )
C
      DO IP = IP1, IP2-1
        DXP = XP(IP+1) - XP(IP)
        DYP = YP(IP+1) - YP(IP)
        DSP(IP) = SQRT(DXP**2 + DYP**2)
      ENDDO
C
C---- march down wake setting wake points
      IP = IPCO(IC1)
      XPNEW(IP) = XP(IP)
      YPNEW(IP) = YP(IP)
C
      DO IC = IC1, IC2
        IPO = IPCO(IC)
        IPP = IPCP(IC)
C
        IF(ISQ .EQ. 0) THEN
          QX = QC(1,IC)
          QY = QC(2,IC)
        ELSEIF(ISQ .EQ. 1) THEN
          QX = QCR(1,IC)
          QY = QCR(2,IC)
        ELSEIF(ISQ .EQ. 2) THEN
          QX = QCL(1,IC)
          QY = QCL(2,IC)
        ELSE
          QX = QC(1,IC)
          QY = QC(2,IC)
        ENDIF
C
        QMAG = SQRT(QX**2 + QY**2)
        IF(QMAG.LT.QSTAG) THEN
C------- stagnation... don't try to move points
         RETURN
        ENDIF
C
        UN = QX/QMAG
        VN = QY/QMAG
        XPNEW(IPP) = XPNEW(IPO) + UN*DSP(IPO)
        YPNEW(IPP) = YPNEW(IPO) + VN*DSP(IPO)
C
      END DO
C
      RETURN
      END ! WAKMOVR




      SUBROUTINE WAKMOV0(IELA)
C-----------------------------------------------
C     Moves existing wake of element IELA by
C     aligning wake with current velocities.
C-----------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(IEWAKE(IELA).EQ.0) THEN
C----- this element has no wake
       WRITE(*,*) 'WAKINT: Element',IELA,'  has no wake!'
       RETURN
      ENDIF
C
      IEL = IEWAKE(IELA)
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
      CALL WAKMOV1(IEL,XP,YP)
C
      CALL XYCSET(IEL)
      CALL ANPSET(IEL)
C
C---- spline new wake
      I = IPFRST(IEL)
      N = IPLAST(IEL) - IPFRST(IEL) + 1
      CALL SCALC(XP(I),YP(I),SP(I),N)
      CALL SEGSPL(XP(I),XPS(I),SP(I),N)
      CALL SEGSPL(YP(I),YPS(I),SP(I),N)
C
C
      WRITE(*,*) 'Generating dq/d(gam,sig) influence of moved wake...'
      CALL QAIC1(.FALSE., 1,NC, IEL,IEL)
      DO KEL = 1, NEL
        ICE = NCX + KEL
        IF(LBODY(KEL)) CALL QAIC1(.FALSE.,ICE,ICE,IEL,IEL)
      ENDDO
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      CALL GSYS(.FALSE., .FALSE., .FALSE., IP1,IP2, 1,0)
C
C
C---- free existing unit wake strength vectors (if any), mark vectors invalid
      DO IP = IP1, IP2
        LUSET(IUGAM(IP)) = .FALSE.
        LUSET(IUSIG(IP)) = .FALSE.
        IUGAM(IP) = 0
        IUSIG(IP) = 0
      ENDDO
C
      RETURN
      END ! WAKMOV0


      SUBROUTINE WAKMOV1(IEL,XPNEW,YPNEW)
      INCLUDE 'DFDC.INC'
      DIMENSION XPNEW(*), YPNEW(*)
C
      DIMENSION DSP(IPX)
C
C---- speed below this is treated as stagnation
      QSTAG = 0.0001*QREF
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
C
C---- calculate velocity and its Jacobians at all wake points
      CALL QAIC1(.TRUE., IC1,IC2, 1,NEL)
C
      DO IP = IP1, IP2-1
        DXP = XP(IP+1) - XP(IP)
        DYP = YP(IP+1) - YP(IP)
        DSP(IP) = SQRT(DXP**2 + DYP**2)
      ENDDO
C
C---- march down wake setting wake points
      IP = IPCO(IC1)
      XPNEW(IP) = XP(IP)
      YPNEW(IP) = YP(IP)
C
      DO 50 IC = IC1, IC2
        IPO = IPCO(IC)
        IPP = IPCP(IC)
C
        QMAG = SQRT(QC(1,IC)**2 + QC(2,IC)**2)
        IF(QMAG.LT.QSTAG) THEN
C------- stagnation... don't try to move wake
         GO TO 51
        ENDIF
C
        UN = QC(1,IC)/QMAG
        VN = QC(2,IC)/QMAG
C
        XPNEW(IPP) = XPNEW(IPO) + UN*DSP(IPO)
        YPNEW(IPP) = YPNEW(IPO) + VN*DSP(IPO)
 50   CONTINUE
 51   CONTINUE
C
      RETURN
      END ! WAKMOV1


      SUBROUTINE WAKINT(XS,YS,DS,IELA)
C-----------------------------------------------
C     Trace a streamline to be a wake
C     Wake then becomes an additional element
C     added to the current list.
C-----------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- iteration limit for wake point convergence
      DATA NITERV / 5 /
C
C---- wake node position convergence tolerance
      DATA TOLER / 0.001 /
C
c      IF(IEWAKE(IELA).NE.0) THEN
C----- this element already has a wake
c       WRITE(*,*) 'WAKINT: Element',IELA,'  already has a wake!'
c       RETURN
c      ENDIF
C
      IF(NEL.GE.NEX) THEN
       WRITE(*,*) 'Cannot generate wake.  Will exceed array limit NEX.'
       RETURN
      ENDIF
C
c      write(*,*) 'Enter starting point for streamline x,y,ds'
c      read(*,*) xs,ys,ds
      if(DS.LE.0.0) then
        write(*,*) 'Enter starting delta s for streamline'
        read(*,*) ds
      endif
C
      WGEOM = 1.05
C
C---- speed below this is treated as stagnation
      QSTAG = 0.0001*QREF
C
C---- first,last points of element to receive wake
c      IP1A = IPFRST(IELA)
c      IP2A = IPLAST(IELA)
C
C---- set length of first wake panel, estimate tangent vector
c      IF(LBODY(IELA) .AND. .NOT.LXBOD(IELA)) THEN
C----- closed body (if axisymmetric geometry, body is not on the axis)
c       SDEL1 = 0.5*(SP(IP1A+1) - SP(IP1A) + SP(IP2A) - SP(IP2A-1))
c       DXT = XPS(IP2A) - XPS(IP1A)
c       DYT = YPS(IP2A) - YPS(IP1A)
C
c      ELSEIF(LXBOD(IELA)) THEN
C----- axisymmetric body on the axis
c       IC1 = ICFRST(IELA)
c       IC2 = ICLAST(IELA)
c       CALL QTCALC(1,ANC(1,IC1),IPCO(IC1),IPCP(IC1), 
c     &             GAM,QC(1,IC1), QCT1)
c       CALL QTCALC(1,ANC(1,IC2),IPCO(IC2),IPCP(IC2), 
c     &             GAM,QC(1,IC2), QCT2)
C

c       SDEL1 = SP(IP2A) - SP(IP2A-1)
c       DXT = XPS(IP2A)
c       DYT = YPS(IP2A)
C
c      ELSE
C----- zero-thickness surface
c       SDEL1 = SP(IP2A) - SP(IP2A-1)
c       DXT = XPS(IP2A)
c       DYT = YPS(IP2A)
C
c      ENDIF
c      DST = SQRT(DXT*DXT + DYT*DYT)
C
c      SPDELT = SDEL1


      SPDELT = DS
C
C
C---- element index of newly created wake
      IEL = NEL + 1
C
      IP = NP + 1
      IC = NC + 1
C
      IF(IP.GE.IPX) THEN
       WRITE(*,*) 'WAKINT: Array limit IPX exceeded.'
       RETURN
      ENDIF
      IF(IC.GE.NCX) THEN
       WRITE(*,*) 'WAKINT: Array limit NCX exceeded.'
       RETURN
      ENDIF
C
C---- set first wake point, estimate segment to second wake point
c      XP(IP) = XPTE(IELA)
c      YP(IP) = YPTE(IELA)
      XP(IP) = xs
      YP(IP) = ys
      XPDEL  = ds
      YPDEL  = 0.0
C
C---- first node and control point indices on wake
      IPFRST(IEL) = IP
      ICFRST(IEL) = IC
C
      NETYPE(IEL) = 1
C
C---- march down wake setting wake points
      NSTEP = IPX - IP
      DO 50 ISTEP = 1, NSTEP
        ICTYPE(IC) = 1
        IPCO(IC) = IP
        IPCP(IC) = IP+1
C
C------ iterate on location of next wake point IP
        DO ITERV = 1, NITERV
          XC(IC) = XP(IP) + 0.5*XPDEL
          YC(IC) = YP(IP) + 0.5*YPDEL
C
C-------- calculate velocity and its Jacobians at point IC
          CALL QAIC1(.TRUE., IC,IC, 1,NEL)
C
          QMAG = SQRT(QC(1,IC)**2 + QC(2,IC)**2)
          IF(QMAG.LT.QSTAG) THEN
C--------- stagnation... don't try to iterate on wake trajectory
           GO TO 30
          ENDIF
C
C-------- streamline curvature at point IC
          UN = QC(1,IC)/QMAG
          VN = QC(2,IC)/QMAG
          CURVC = (  UN*UN* QC_XC(2,IC)
     &             + UN*VN*(QC_YC(2,IC)-QC_XC(1,IC))
     &             - VN*VN* QC_YC(1,IC)             ) / QMAG
C
C-------- set up 2x2 Newton system for panel node IP+1 ...
C
          SPDEL = SQRT(XPDEL**2 + YPDEL**2)
C
C-------- residual 1:  panel aligns with velocity vector  dr x q = 0
          Z1 = XPDEL*QC(2,IC)
     &       - YPDEL*QC(1,IC)
          Z1_XP =  QC(2,IC) + XPDEL*QC_XC(2,IC)*0.5
     &                      - YPDEL*QC_XC(1,IC)*0.5
          Z1_YP = -QC(1,IC) + XPDEL*QC_YC(2,IC)*0.5
     &                      - YPDEL*QC_YC(1,IC)*0.5
C
C-------- residual 2: panel length is imposed   sqrt(dr.dr) - ds = 0
          Z2 = SPDEL - SPDELT
          Z2_XP = XPDEL/SPDEL
          Z2_YP = YPDEL/SPDEL
C
          DET =   Z1_XP*Z2_YP - Z1_YP*Z2_XP
          DXP = -(Z1   *Z2_YP - Z1_YP*Z2   )/DET
          DYP = -(Z1_XP*Z2    - Z1   *Z2_XP)/DET
          DSP = SQRT(DXP**2 + DYP**2)
C
          DSPTOL = SPDELT*TOLER
          DSPMAX = SPDELT*0.2
C
          IF    (DSP .LT. DSPTOL) THEN
           GO TO 30
          ELSEIF(DSP .GT. DSPMAX) THEN
           RLX = DSPMAX/DSP
          ELSE
           RLX = 1.0
          ENDIF
C
          XPDEL = XPDEL + RLX*DXP
          YPDEL = YPDEL + RLX*DYP
        ENDDO
        WRITE(*,*) 'WAKINT: Wake point iteration failed. dxy =', DXP,DYP
C
 30     CONTINUE
C------ set new wake point
        IF(IP.GE.IPX) THEN
         WRITE(*,*) 'WAKINT: Array limit IPX exceeded.'
         GO TO 90
        ENDIF
C
        XP(IP+1) = XP(IP) + XPDEL
        YP(IP+1) = YP(IP) + YPDEL
C
C------ see if wake is long enough
        IF(XP(IP+1).LT.XWBOX(1) .OR.
     &     XP(IP+1).GT.XWBOX(2) .OR.
     &     YP(IP+1).LT.YWBOX(1) .OR.
     &     YP(IP+1).GT.YWBOX(2)     ) THEN
C-------- terminate wake trajectory integration
          GO TO 90
        ENDIF
C
        IP = IP + 1
        IC = IC + 1
        IF(IP.GE.IPX) THEN
         WRITE(*,*) 'WAKINT: Array limit IPX exceeded.'
         RETURN
        ENDIF
        IF(IC.GE.NCX) THEN
         WRITE(*,*) 'WAKINT: Array limit NCX exceeded.'
         RETURN
        ENDIF

C------ set new target panel length 
        XPDELT = XPDELT * WGEOM
        YPDELT = YPDELT * WGEOM
        SPDELT = SPDELT * WGEOM
C
 50   CONTINUE
      WRITE(*,*) 'WAKINT:  Wake incomplete.  Array limit reached.'
C
 90   CONTINUE
      IPLAST(IEL) = IP+1
      ICLAST(IEL) = IC
C
      NP = IP+1
      NC = IC
      NEL = IEL
C
c      IEWAKE(IELA) = IEL
      LBODY(IEL) = .FALSE.
      LXBOD(IEL) = LXBOD(IELA)
      XPCENT(IEL) = 0.5*(XP(IPFRST(IEL)) + XP(IPLAST(IEL)))
      YPCENT(IEL) = 0.5*(YP(IPFRST(IEL)) + YP(IPLAST(IEL)))
C
      CALL XYPSPL(IEL)
C
      CALL XYCSET(IEL)
      CALL ANPSET(IEL)
C
C---- set system indices for new wake
      DO IC = ICFRST(IEL), ICLAST(IEL)
        KSYSVNC(IC) = -999
        KSYSGMG(IC) = -999
        IELDOFC(IC) = 0
      ENDDO
      DO IP = IPFRST(IEL), IPLAST(IEL)
        KSYSGAM(IP) = -999
        JSYSGAM(IP) = 0
        JSYSDXY(IP) = 0
        JAICGAM(IP) = 0
        JAICSIG(IP) = 0
      ENDDO
      KSYSKUT(IEL) = -999
      KSYSQNC(IEL) = -999
      KSYSGSS(IEL) = -999
      JSYSQNC(IEL) = 0
C
c      JAICG = NAICGAM
c      JAICS = NAICSIG
c      DO JV = IPFRST(IEL), IPLAST(IEL)
C------ wake gamma will be imposed externally (e.g. BL curvature relation)
c        JAICG = JAICG + 1
c        JAICGAM(JV) = JAICG
C
C------ wake sigma will be imposed externally (e.g. BL blowing model)
c        JAICS = JAICS + 1
c        JAICSIG(JV) = JAICS
c      ENDDO
C
c      NAICGAM = JAICG
c      NAICSIG = JAICS
C
C
      WRITE(*,*) 'Generating dq/d(gam,sig) influence of new wake...'
      CALL QAIC1(.FALSE., 1,NC, IEL,IEL)
      DO KEL = 1, NEL
        ICE = NCX + KEL
        IF(LBODY(KEL)) CALL QAIC1(.FALSE.,ICE,ICE,IEL,IEL)
      ENDDO
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      CALL GSYS(.FALSE.,.FALSE., .FALSE., IP1,IP2, 1,0)
C
      LQCNT = .TRUE.
      LQAIC = .TRUE.
      LQGIC = .FALSE.
C
C---- free existing unit wake strength vectors (if any), mark vectors invalid
      DO IP = IP1, IP2
        LUSET(IUGAM(IP)) = .FALSE.
        LUSET(IUSIG(IP)) = .FALSE.
        IUGAM(IP) = 0
        IUSIG(IP) = 0
      ENDDO
C
      RETURN
      END ! WAKINT

