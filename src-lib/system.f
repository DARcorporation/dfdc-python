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
C This version assumes delta(B*GAM)*WT/WM defined at each rotor center 

      SUBROUTINE SYSP
C-------------------------------------------
C     Sets up pointers associating....
C        matrix rows with control points,
C        matrix columns with unknowns
C-------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) WRITE(*,*) 'Generating flow-tangency system pointers'
C
C---- set what type of equation gets assigned to each control point,
C-    and which unknown vortex gets associated with each control point
C
C---- -999 indicates un-enforced control point
C----    0 indicates undefined variable
      DO IC = 1, NCTOT
        KSYSVNC(IC) = -999
        KSYSGMG(IC) = -999
        IELDOFC(IC) = 0
      ENDDO
      DO IP = 1, NPTOT
        KSYSGAM(IP) = -999
        JSYSGAM(IP) = 0
        JSYSDXY(IP) = 0
        JAICGAM(IP) = 0
        JAICSIG(IP) = 0
        JAICXYP(IP) = 0
        JAICGTH(IP) = 0
      ENDDO
      DO IEL = 1, NEL
        KSYSKUT(IEL) = -999
        KSYSQNC(IEL) = -999
        KSYSGSS(IEL) = -999
        KSYSQTT(IEL) = -999
        JSYSQNC(IEL) = 0
        JAICVWK(IEL) = 0
      ENDDO
C
C
C---- initialize system-row and system-column counters
      KSYS = 0
      JSYS = 0
C
C---- initialize gamma,sigma,xy column counters
      JAICG = 0
      JAICS = 0
      JAICZ = 0
C---- initialize wake gamma (point) column counters
      JAICW = 0
C---- initialize vortex wake row column counters
      JAICV = 0
C
C---- go over all elements
      DO IEL = 1, NEL
C
ccC------ LVeZR = f if first,last vortex node is not an unknown 
cc        LV1ZR(IEL) = YP(IPFRST(IEL)) .EQ. 0.0
cc        LV2ZR(IEL) = YP(IPLAST(IEL)) .EQ. 0.0
C
        IF(NETYPE(IEL).EQ.0) THEN
C-------- normal solid-wall airfoil... 
C
C--------------------------------------------------
C-------- set up variable indices
C
          DO JP = IPFRST(IEL), IPLAST(IEL)
C---------- wall gamma will be unknown
            JSYS = JSYS + 1
            JSYSGAM(JP) = JSYS
            JPSYS(JSYS) = JP
C---------- wall sigma will be imposed externally (e.g. BL blowing model)
            JAICS = JAICS + 1
            JAICSIG(JP) = JAICS
          ENDDO
C
          IF(LBODY(IEL)) THEN
C--------- declare normal-velocity DOF, and set element which it influences
           JSYS = JSYS + 1
           JSYSQNC(IEL) = JSYS
           JPSYS(JSYS) = 0
           NDOF = IEL
          ELSE
           NDOF = 0
          ENDIF
C
          IF(LXYMOV) THEN
           DO JP = IPFRST(IEL), IPLAST(IEL)
C----------- wall x,y movement will be imposed externally (e.g. optimization)
             JAICZ = JAICZ + 1
             JAICXYP(JP) = JAICZ
           ENDDO
          ENDIF
C
C--------------------------------------------------
C-------- set up equation indices
C
          IF(LV1ZR(IEL)) THEN
C---------- enforce zero-gamma condition at first point
            KSYS = KSYS + 1
            KSYSGAM(IPFRST(IEL)) = KSYS
          ENDIF
C
C-------- enforce flow tangency at all control points, note Qn DOF if any
          DO IC = ICFRST(IEL), ICLAST(IEL)
            KSYS = KSYS + 1
C
            IF(ICTYPE(IC).GE.0) THEN
C----------- usual control point
             KSYSVNC(IC) = KSYS
             IELDOFC(IC) = NDOF
            ELSE
C----------- zero-length panel... set gamma-continuity equation
             KSYSGMG(IC) = KSYS
            ENDIF
          ENDDO
C
          IF(LV2ZR(IEL)) THEN
C---------- enforce zero-gamma condition at last point
            KSYS = KSYS + 1
            KSYSGAM(IPLAST(IEL)) = KSYS
          ENDIF
C
          IF(.NOT.(LV1ZR(IEL).OR.LV2ZR(IEL))) THEN
C--------- enforce Kutta condition
           KSYS = KSYS + 1
           KSYSKUT(IEL) = KSYS
          ENDIF
C
          IF(LBODY(IEL) .AND. .NOT.(LV1ZR(IEL).AND.LV2ZR(IEL))) THEN
C--------- will have normal-velocity DOF... enforce some constraint on it
           KSYS = KSYS + 1
           IF(LQNZR(IEL)) THEN
C---------- direct constraint
            KSYSQNC(IEL) = KSYS
           ELSE
cccC---------- indirect constraint: set zero loading curvature at TE
ccc            KSYSGSS(IEL) = KSYS
C
C---------- indirect constraint: set zero velocity at interior point near TE
            KSYSQTT(IEL) = KSYS
           ENDIF
          ENDIF
C
C-------- check for vortex wake gamma on body element
          DO JP = IPFRST(IEL), IPLAST(IEL)
            IF(IP2IR(JP).GT.0) THEN
C---------- GTH will be imposed or solved for at this point
              JAICW = JAICW + 1
              JAICGTH(JP) = JAICW
            ENDIF
          ENDDO
C
        ELSEIF(NETYPE(IEL).EQ.1) THEN
C-------- wake element...
C---------- wake gamma will be imposed externally (e.g. BL curvature relation)
          DO JP = IPFRST(IEL), IPLAST(IEL)
            JAICG = JAICG + 1
            JAICGAM(JP) = JAICG
C---------- wake sigma will be imposed externally (e.g. BL blowing model)
cc            JAICS = JAICS + 1
cc            JAICSIG(JP) = JAICS
          ENDDO
C
          IF(LXYMOV) THEN
           DO JP = IPFRST(IEL), IPLAST(IEL)
C----------- wake x,y movement will be imposed externally (e.g. optimization)
             JAICZ = JAICZ + 1
             JAICXYP(JP) = JAICZ
           ENDDO
          ENDIF
C
        ELSEIF(NETYPE(IEL).GE.2 .AND. NETYPE(IEL).LE.4) THEN
C-------- axis line, ring or point singularities...(NETYPE=2,3,4)
C-------- sigma, doublet will be imposed externally
          DO JP = IPFRST(IEL), IPLAST(IEL)
            JAICG = JAICG + 1
            JAICGAM(JP) = JAICG
            JAICS = JAICS + 1
            JAICSIG(JP) = JAICS
          ENDDO
C
        ELSEIF(NETYPE(IEL).EQ.5) THEN
C-------- drag area source element NETYPE=5 (source RHS only)
C-------- sigma will be imposed externally (e.g. drag source)
          DO JP = IPFRST(IEL), IPLAST(IEL)
            JAICS = JAICS + 1
            JAICSIG(JP) = JAICS
          ENDDO
C
        ELSEIF(NETYPE(IEL).EQ.6) THEN
C-------- rotor source element NETYPE=6 (source RHS only)
C-------- sigma will be imposed externally (e.g. drag source)
          DO JP = IPFRST(IEL), IPLAST(IEL)
            JAICS = JAICS + 1
            JAICSIG(JP) = JAICS
          ENDDO
C
        ELSEIF(NETYPE(IEL).EQ.7) THEN
C-------- vortex wake element (RHS)
C-------- wake gamma will be imposed externally
          IF(JAICV.LT.NRP) THEN
            JAICV = JAICV + 1
            IR = IEL2IR(IEL)
            JAICVWK(IR) = JAICV
cc            write(19,*) 'iel,ir,jaicv ',iel,ir,jaicv
          ENDIF
C-------- set up points with vortex wake gamma 
          DO JP = IPFRST(IEL), IPLAST(IEL)
C---------- GTH will be imposed or solved for at this point
            JAICW = JAICW + 1
            JAICGTH(JP) = JAICW
          ENDDO
C
        ELSE
          WRITE(*,*) 'Unknown NETYPE in SYSP ', NETYPE(IEL)
C
        ENDIF
      ENDDO
C
C---- set system size
      NSYS = KSYS
C
      IF(KSYS.NE.JSYS) THEN
       WRITE(*,*) 'Error in SYSP, KSYS JSYS =', KSYS, JSYS
      ENDIF
C
      NAICGAM = JAICG
      NAICSIG = JAICS
      NAICXYP = JAICZ
      NAICVWK = JAICV
      NAICGTH = JAICW
C
      IF(LDBG) THEN
       WRITE(*,*) ' '
       WRITE(*,*) 'SYSP system setup and pointers'
       WRITE(*,*) '  NP      = ',NPTOT
       WRITE(*,*) '  NC      = ',NCTOT
       WRITE(*,*) '  NSYS    = ',NSYS
       WRITE(*,*) ' '
       WRITE(*,*) '  NAICGAM = ',NAICGAM
       WRITE(*,*) '  NAICSIG = ',NAICSIG
       WRITE(*,*) '  NAICXYP = ',NAICXYP
       WRITE(*,*) '  NAICVWK = ',NAICVWK
       WRITE(*,*) '  NAICGTH = ',NAICGTH
       WRITE(*,*) ' '
       WRITE(*,*) ' '
C
       WRITE(LUNDBG,*) 'SYSP system setup and pointers'
       WRITE(LUNDBG,*) '  NP      = ',NPTOT
       WRITE(LUNDBG,*) '  NC      = ',NCTOT
       WRITE(LUNDBG,*) '  NSYS    = ',NSYS
       WRITE(LUNDBG,*) ' '
       WRITE(LUNDBG,*) '  NAICGAM = ',NAICGAM
       WRITE(LUNDBG,*) '  NAICSIG = ',NAICSIG
       WRITE(LUNDBG,*) '  NAICXYP = ',NAICXYP
       WRITE(LUNDBG,*) '  NAICVWK = ',NAICVWK
       WRITE(LUNDBG,*) '  NAICGTH = ',NAICGTH
       WRITE(LUNDBG,*) ' '
      ENDIF
C
      LSYSP = .TRUE.
C
      RETURN
      END ! SYSP


      SUBROUTINE CLRSYS
      INCLUDE 'DFDC.INC'
C
C---- clear all system arrays
      DO K = 1, NSYS
        RES(K) = 0.
        DO J = 0, NSYS
          SYS(K,J) = 0.
        ENDDO
        AICQFF(K,1) = 0.
        DO J = 0, NAICGAM
          AICGAM(K,J) = 0.
        ENDDO
        DO J = 0, NAICSIG
          AICSIG(K,J) = 0.
        ENDDO
        DO J = 0, NAICXYP
          AICXYP(K,1,J) = 0.
          AICXYP(K,2,J) = 0.
        ENDDO
        DO J = 0, NAICVWK
          AICVWK(K,J) = 0.
        ENDDO
      ENDDO
C
      RETURN
      END ! CLRSYS


      SUBROUTINE GSYS(LSYS,LQFF,LVWK, IP1,IP2, IZ1,IZ2)
C-----------------------------------------------------------
C     Generates influence matrices of unknown gamma 
C     and freestream on flow tangency and other residuals.
C
C     Input:  LSYS     if T, set up surface-gamma Jacobian
C             LQFF     if T, set up freestream Jacobian
C             IP1,IP2  set up Jacobian for gamma,sigma of nodes IP1..IP2
C             IZ1,IZ2  set up Jacobian for x,y of nodes IZ1..IZ2
C
C     Input:  NSYS         system size
C             GAM(j)       panel-node vortex sheet strengths
C             QNDOF(e)     element normalwash DOF variables
C             ANC(.i)      control-point normal vectors
C             QC(.i)       control-point velocities
C             QC_GAM(.ij)  AIC matrix  dQC(.i)/dGAM(j) on surfaces
C             QC_SIG(.ij)  AIC matrix  dQC(.i)/dSIG(j) on surfaces
C             QC_XYP(.i.j) AIC matrix  dQC(.i)/dXP(j),dYP(j)
C             QC_GTH(.ij)  AIC matrix  dQC(.i)/dGTH(j) on surfaces
C             KSYS...      active-equation pointers
C             JSYS...      active-variable pointers
C
C     Output: linear system components
C
C        LHS system
C             RES(k)       Qn(k)  for current QC(j), QNDOF(e)
C                      or  GAM-continuity residual
C                      or  Kutta condition residual
C                      or  TE loading regularity residual
C             SYS(kj)      dRES(k)/dGAM(j) , dRES(k)/dQNDOF(e)   if LSYS=t
C
C        RHS columns  
C             AICQFF(k)    dRES(k)/dQinf                         if LQFF=t
C             AICVWK(e)    dRES(k)/dWBGAM on rotor               if LVWK=t
C
C             AICGAM(kj)   dRES(k)/dGAM(j)    j = IP1..IP2
C             AICSIG(kj)   dRES(k)/dSIG(j)    j = IP1..IP2
C
C             AICXYP(k.j)  dRES(k)/dXP(j)     j = IZ1..IZ2
C                          dRES(k)/dYP(j)     j = IZ1..IZ2
C
C-----------------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LSYS,LQFF,LVWK
C
      IF(NSYS.EQ.0) RETURN
C
      IF(LDBG) THEN
       WRITE(*,*) 'Setting up direct-problem Newton system...'
       WRITE(LUNDBG,*) 'GSYS direct problem system setup'
       WRITE(LUNDBG,*) '   LSYS,LQFF,LVWK ',LSYS,LQFF,LVWK
       WRITE(LUNDBG,*) '   IP1,IP2        ',IP1,IP2
       WRITE(LUNDBG,*) '   IZ1,IZ2        ',IZ1,IZ2
      ENDIF
C
C----------------------------------------------------------
C---- clear system arrays
      IF(LSYS) THEN
       DO K = 1, NSYS
         RES(K) = 0.
         DO J = 1, NSYS
           SYS(K,J) = 0.
         ENDDO
       ENDDO
      ENDIF
C
C---- clear RHS arrays
      IF(LQFF) THEN
       DO K = 1, NSYS
         AICQFF(K,1) = 0.
       ENDDO
      ENDIF
C
      IF(LVWK) THEN
       DO IR = 1, NRP
         J = JAICVWK(IR)
         DO K = 1, NSYS
           AICVWK(K,J) = 0.
         ENDDO
       ENDDO
      ENDIF
C
      DO K = 1, NSYS
       DO JP = IP1, IP2
         J = JAICGAM(JP)
         AICGAM(K,J) = 0.
C
         J = JAICSIG(JP)
         AICSIG(K,J) = 0.
C
         J = JAICGTH(JP)
         AICGTH(K,J) = 0.
       ENDDO
C
       DO JP = IZ1, IZ2
         J = JAICXYP(JP)
         AICXYP(K,1,J) = 0.
         AICXYP(K,2,J) = 0.
       ENDDO
      ENDDO
C
C----------------------------------------------------------
C----------------------------------------------------------
C---- set flow tangency at all control points
      DO 10 IC = 1, NCTOT
C
c        write(77,*) 'IC ',IC
c        do jp = 1, np
c          write(77,*) 'JP, QC1-2 ',JP,QC_GTH(1,JP,IC),QC_GTH(2,JP,IC)
c        end do
C------ K = system row index for control point IC
        K = KSYSVNC(IC)
C
C------ skip un-enforced control points
        IF(K.LT.0) GO TO 10
C
C----------------------------------------------------------
        IF(LSYS) THEN
C------- flow-tangency residual
         RES(K) = ANC(1,IC)*QC(1,IC)
     &          + ANC(2,IC)*QC(2,IC)
C
C
         DO JP = 1, NPTOT
C--------- dRES/dgamma_surface
C-         J = system column index of GAM(JP)  (0 if GAM is on r.h.s.)
           J = JSYSGAM(JP)
           SYS(K,J) = ANC(1,IC)*QC_GAM(1,JP,IC)
     &              + ANC(2,IC)*QC_GAM(2,JP,IC)
         ENDDO
        ENDIF
C
C----------------------------------------------------------
        IF(LQFF) THEN
C------- dRES/dQinf
         AICQFF(K,1) = ANC(1,IC)
        ENDIF
C
C----------------------------------------------------------
        IF(LVWK) THEN
cc          write(19,*) 'ic ',ic
C------- dRES/dBCirc rotor vortex wake influence
          DO JP = 1, NPTOT
            IR = IP2IR(JP)
            IF(IR.NE.0) THEN
C---- we are on a vortex wake element
C
C---- accumulate wake influences by wake rows using wake index
              J = JAICVWK(IR)
              AICVWK(K,J) = AICVWK(K,J) +
     &                     ( ANC(1,IC)*QC_GTH(1,JP,IC) +
     &                       ANC(2,IC)*QC_GTH(2,JP,IC) ) / YP(JP)
cc              write(19,*) 'jp,ir,j ',jp,ir,j
C
C---- set wake point influence individually
              J = JAICGTH(JP)
              AICGTH(K,J) = ANC(1,IC)*QC_GTH(1,JP,IC)
     &                    + ANC(2,IC)*QC_GTH(2,JP,IC)
            ENDIF
          ENDDO
        ENDIF
C
C----------------------------------------------------------
        DO JP = IP1, IP2
C-------- dRES/dgamma, dRES/dsigma
C-        J = r.h.s. column index  (0 if GAM,SIG is a l.h.s. unknown)
          J = JAICGAM(JP)
          AICGAM(K,J) = ANC(1,IC)*QC_GAM(1,JP,IC)
     &                + ANC(2,IC)*QC_GAM(2,JP,IC)
C
          J = JAICSIG(JP)
          AICSIG(K,J) = ANC(1,IC)*QC_SIG(1,JP,IC)
     &                + ANC(2,IC)*QC_SIG(2,JP,IC)
        ENDDO
C
C----------------------------------------------------------
        DO JP = IZ1, IZ2
C-------- dRES/dxy
C-        J = r.h.s. column index
          J = JAICXYP(JP)
          AICXYP(K,1,J) = ANC(1,IC)*QC_XP(1,JP,IC)
     &                  + ANC(2,IC)*QC_XP(2,JP,IC)
          AICXYP(K,2,J) = ANC(1,IC)*QC_YP(1,JP,IC)
     &                  + ANC(2,IC)*QC_YP(2,JP,IC)
        ENDDO
C
C----------------------------------------------------------
C------ local source contribution to flow-tangency residual
        IPO = IPCO(IC)
        IPP = IPCP(IC)
C
        IF(LSYS) THEN
         RES(K) = RES(K) - 0.25*(SIG(IPO) + SIG(IPP))
        ENDIF
C
        IF(IPO.GE.IP1 .AND. IPO.LE.IP2) THEN
         J = JAICSIG(IPO)
         AICSIG(K,J) = AICSIG(K,J) - 0.25
        ENDIF
        IF(IPP.GE.IP1 .AND. IPP.LE.IP2) THEN
         J = JAICSIG(IPP)
         AICSIG(K,J) = AICSIG(K,J) - 0.25
        ENDIF
C
C------ derivatives w.r.t. dn movement of the two nodes of panel IC
        RES_DX = ANC_DXY(1,1,IC)*QC(1,IC)
     &         + ANC_DXY(2,1,IC)*QC(2,IC)
        RES_DY = ANC_DXY(1,2,IC)*QC(1,IC)
     &         + ANC_DXY(2,2,IC)*QC(2,IC)
C
        IF(IPO.GE.IZ1 .AND. IPO.LE.IZ2) THEN
         J = JAICXYP(IPO)
         AICXYP(K,1,J) = AICXYP(K,1,J) - RES_DX
         AICXYP(K,2,J) = AICXYP(K,2,J) - RES_DY
        ENDIF
        IF(IPP.GE.IZ1 .AND. IPP.LE.IZ2) THEN
         J = JAICXYP(IPP)
         AICXYP(K,1,J) = AICXYP(K,1,J) + RES_DX
         AICXYP(K,2,J) = AICXYP(K,2,J) + RES_DY
        ENDIF
C
C----------------------------------------------------------
        IF(LSYS .AND. IELDOFC(IC).GT.0) THEN
C------- this control point has a normal-velocity QNDOF contributing to Q.n
         IE = IELDOFC(IC)
         RES(K) = RES(K) - QNDOF(IE)
C
C------- J = system column index of QNDOF(IE)
         J = JSYSQNC(IE)
         IF(J.LE.0) THEN
           WRITE(*, *) '? Unassigned normal-velocity DOF', IE
           WRITE(LUNDBG,*) '? Unassigned normal-velocity DOF', IE
         ENDIF
C
         SYS(K,J) = -1.0
        ENDIF
        IF(LDBG) WRITE(LUNDBG,*) K, 'V.n=Qn ',RES(K) !@@@
 10   CONTINUE
C
C---- set gamma-continuity at all zero-length panels
      DO 15 IC = 1, NCTOT
C------ K = system row index of gamma-continuity equation (if any)
        K = KSYSGMG(IC)
C
        IF(LSYS .AND. K.GT.0) THEN
         JPO = IPCO(IC)
         JPP = IPCP(IC)
         RES(K) = GAM(JPP) - GAM(JPO)
C
         SYS(K,JPO) = -1.0
         SYS(K,JPP) =  1.0
         IF(LDBG) WRITE(LUNDBG,*) K, 'dGam=0 ',RES(K) !@@@
        ENDIF
 15   CONTINUE
C
C
      DO 20 IEL = 1, NEL
C------ zero net curvature of gamma at TE (triangular loading towards TE)
        K = KSYSGSS(IEL)
        IF(LSYS .AND. K.GT.0) THEN

          JP1O = IPFRST(IEL) + 2
          JP1M = IPFRST(IEL) + 1
          JP1L = IPFRST(IEL)
          JP2O = IPLAST(IEL)
          JP2M = IPLAST(IEL) - 1
          JP2L = IPLAST(IEL) - 2
C
          DS1O = SP(JP1O) - SP(JP1M)
          DS1M = SP(JP1M) - SP(JP1L)
          DS2O = SP(JP2O) - SP(JP2M)
          DS2M = SP(JP2M) - SP(JP2L)
C
C-------- 1st derivative
c          F1O = 1.0/DS1O + 1.0/(DS1O+DS1M)
c          F1L = DS1O/((DS1O+DS1M)*DS1M)
c          F1M = -(F1O+F1L)
c          F2O = 1.0/DS2O + 1.0/(DS2O+DS2M)
c          F2L = DS2O/((DS2O+DS2M)*DS2M)
c          F2M = -(F2O+F2L)
C
C-------- 2nd derivative
          F1O =   1.0/DS1O
          F1M = -(1.0/DS1O + 1.0/DS1M)
          F1L =              1.0/DS1M
          F2O =   1.0/DS2O
          F2M = -(1.0/DS2O + 1.0/DS2M)
          F2L =              1.0/DS2M
C
          J1O = JSYSGAM(JP1O)
          J1M = JSYSGAM(JP1M)
          J1L = JSYSGAM(JP1L)
          J2O = JSYSGAM(JP2O)
          J2M = JSYSGAM(JP2M)
          J2L = JSYSGAM(JP2L)
C
          IF    (LV1ZR(IEL)) THEN
           RES(K)     = -F2O*GAM(JP2O)
     &                  -F2M*GAM(JP2M)
     &                  -F2L*GAM(JP2L)
           SYS(K,J2O) = -F2O
           SYS(K,J2M) = -F2M
           SYS(K,J2L) = -F2L
C
          ELSEIF(LV2ZR(IEL)) THEN
           RES(K)     =  F1O*GAM(JP1O)
     &                  +F1M*GAM(JP1M)
     &                  +F1L*GAM(JP1L)
           SYS(K,J1O) =  F1O
           SYS(K,J1M) =  F1M
           SYS(K,J1L) =  F1L
C
          ELSE
           RES(K)     =  F1O*GAM(JP1O)
     &                  +F1M*GAM(JP1M)
     &                  +F1L*GAM(JP1L)
     &                  -F2O*GAM(JP2O)
     &                  -F2M*GAM(JP2M)
     &                  -F2L*GAM(JP2L)
           SYS(K,J1O) =  F1O
           SYS(K,J1M) =  F1M
           SYS(K,J1L) =  F1L
           SYS(K,J2O) = -F2O
           SYS(K,J2M) = -F2M
           SYS(K,J2L) = -F2L
C
          ENDIF
          IF(LDBG) WRITE(LUNDBG,*) K, 'Qss=0 ',RES(K)  !@@@
        ENDIF
C
C------ Kutta condition
        K = KSYSKUT(IEL)
        IF(K.GT.0) THEN
          IF(LSYS) THEN
C
            IF(LBODY(IEL)) THEN
C------------ closed body
              CALL SWDOTS(IEL,DS1SW,DS2SW,DN1SW,DN2SW)
              JP1 = IPFRST(IEL)
              JP2 = IPLAST(IEL)
C
C@@@
              ds1sw = -1.0
              ds2sw =  1.0
              dn1sw = 0.
              dn2sw = 0.
C------------ basic Kutta condition
              RES(K) = DS2SW*GAM(JP2) - DS1SW*GAM(JP1)
     &               + DN2SW*SIG(JP2) - DN1SW*SIG(JP1)
ccc     &               + CIRDOT*2.0/(GAM(JP1)-GAM(JP2))
C
              J1 = JSYSGAM(JP1)
              J2 = JSYSGAM(JP2)
              SYS(K,J1) = -DS1SW
              SYS(K,J2) =  DS2SW
C
              IF(JP1.GE.IP1 .AND. JP1.LE.IP2) THEN
               J1 = JAICSIG(JP1)
               AICSIG(K,J1) = -DN1SW
              ENDIF
              IF(JP2.GE.IP1 .AND. JP2.LE.IP2) THEN
ccc               J2 = JAICGAM(JP2)   !!! bug
               J2 = JAICSIG(JP2)
               AICSIG(K,J2) =  DN2SW
              ENDIF
C
              IF(IEWAKE(IEL).NE.0) THEN
C------------- add wake panel contribution
               IELW = IEWAKE(IEL)
               JPW = IPFRST(IELW)
               RES(K) = RES(K) - GAM(JPW)
               RES(K) = RES(K) - GTH(JPW)   !!HHY
               IF(JPW.GE.IP1 .AND. JP2.LE.IP2) THEN
                JW = JAICGAM(JPW)
                AICGAM(K,JW) = -1.0
c                JW = JAICGTH(JPW)
c                AICGTH(K,JW) = -1.0
               ENDIF
              ENDIF
C
            ELSE
C------------ zero-thickness plate
              JP = IPLAST(IEL)
C
              RES(K) = GAM(JP)
C
              J = JSYSGAM(JP)
              SYS(K,J) = 1.0
C
              IF(IEWAKE(IEL).NE.0) THEN
C------------- add wake panel contribution
               IELW = IEWAKE(IEL)
               JPW = IPFRST(IELW)
               RES(K) = RES(K) - GAM(JPW)
               RES(K) = RES(K) - GTH(JPW)   !!HHY
C
               IF(JPW.GE.IP1 .AND. JP2.LE.IP2) THEN
                JW = JAICGAM(JPW)
                AICGAM(K,JW) = -1.0
c                JW = JAICGTH(JPW)
c                AICGTH(K,JW) = -1.0
               ENDIF
              ENDIF
C
            ENDIF
          ENDIF
          IF(LDBG) WRITE(LUNDBG,*) K, 'Kutta ',RES(K) !@@@
        ENDIF
C
C------ normal-velocity QNDOF constraint
        K = KSYSQNC(IEL)
        IF(LSYS .AND. K.GT.0) THEN

         J = JSYSQNC(IEL)
         RES(K) = QNDOF(IEL)
CCC  &          - QNDOFSP(IEL)
         SYS(K,J) = 1.0
         IF(LDBG) WRITE(LUNDBG,*) K, 'Qn=0 ',RES(K)  !@@@
        ENDIF
C
C------ zero internal tangential velocity constraint
        K = KSYSQTT(IEL)
        IF(K.GT.0) THEN
C
         IC = NCX + IEL
C
         IF(LSYS) THEN

ccc debug stuff to fiddle with net massflow from QNdof
ccc            write(*,*) 'enter Qtt'   !@@@
ccc            read(*,*) qttf    !@@@
               qttf = 0.0

C-------- flow-tangency residual
          RES(K) = ANC(1,IC)*QC(1,IC)
     &           + ANC(2,IC)*QC(2,IC)  -  qttf*QINF   !@@@
C
          IF(LDBG) WRITE(LUNDBG,*) K, 'Qtt=0 ',RES(K)  !@@@

          DO JP = 1, NPTOT
C---------- dRES/dgamma_surface
C-          J = system column index of GAM(JP)  (0 if GAM is on r.h.s.)
            J = JSYSGAM(JP)
            SYS(K,J) = ANC(1,IC)*QC_GAM(1,JP,IC)
     &               + ANC(2,IC)*QC_GAM(2,JP,IC)
          ENDDO
         ENDIF
C
         IF(LQFF) THEN
C-------- dRES/dQinf
          AICQFF(K,1) = ANC(1,IC)   -  qttf   !@@@
         ENDIF
C
C----------------------------------------------------------
         IF(LVWK) THEN
cc          write(19,*) 'icx ',ic
C------- dRES/dBCirc rotor vortex wake influence
          DO JP = 1, NPTOT
            IR = IP2IR(JP)
            IF(IR.NE.0) THEN
C---- we are on a vortex wake element
C
C---- accumulate wake influences by wake rows using wake index
              J = JAICVWK(IR)
              AICVWK(K,J) = AICVWK(K,J) +
     &                     ( ANC(1,IC)*QC_GTH(1,JP,IC) +
     &                       ANC(2,IC)*QC_GTH(2,JP,IC) ) / YP(JP)
cc              write(19,*) 'jp,ir,j ',jp,ir,j
C
C---- set wake point influence individually
              J = JAICGTH(JP)
              AICGTH(K,J) = ANC(1,IC)*QC_GTH(1,JP,IC)
     &                    + ANC(2,IC)*QC_GTH(2,JP,IC)
            ENDIF
          ENDDO
         ENDIF
C
C----------------------------------------------------------
         DO JP = IP1, IP2
C--------- dRES/dgamma, dRES/dsigma
C-         J = r.h.s. column index  (0 if GAM,SIG is a l.h.s. unknown)
           J = JAICGAM(JP)
           AICGAM(K,J) = ANC(1,IC)*QC_GAM(1,JP,IC)
     &                 + ANC(2,IC)*QC_GAM(2,JP,IC)
C
           J = JAICSIG(JP)
           AICSIG(K,J) = ANC(1,IC)*QC_SIG(1,JP,IC)
     &                 + ANC(2,IC)*QC_SIG(2,JP,IC)
         ENDDO
        ENDIF
C
 20   CONTINUE
C
C---- set zero vortex strengths explicitly (if selected)
      IF(LSYS) THEN
       DO 30 JP = IP1, IP2
         K = KSYSGAM(JP)
         IF(K.GT.0) THEN
          J = JSYSGAM(JP)
          RES(K) = GAM(JP)
          SYS(K,J) = 1.0
          IF(LDBG) WRITE(LUNDBG,*) K, 'gam=0 ',RES(K)  !@@@
         ENDIF
 30    CONTINUE
      ENDIF
C
      IF(LSYS) THEN
       LGSYS = .TRUE.
       LGAMU = .FALSE.
       LQGIC = .FALSE.
      ENDIF
C
      RETURN
      END ! GSYS



      SUBROUTINE SYSPQ
C-------------------------------------------
C     Modifies pointers from SYSP to solve
C     mixed-inverse problem
C-------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) THEN
       WRITE(*,*) 'Modifying system pointers for mixed-inverse...'
       WRITE(LUNDBG,*) 'SYSPQ modify system pointers for mixed inverse'
       WRITE(LUNDBG,*) '   NSYS =',NSYS,' (before augmentation)'
      ENDIF
C
C---- -999 indicates un-enforced control point
C----    0 indicates undefined variable
      DO ISEG = 0, NSEG
        DO K = 1, 4
          KSYSDNS(K,ISEG) = -999
          JSYSQSP(K,ISEG) = 0
        ENDDO
      ENDDO
C
C---- assume that no panel node lies inside an inverse segment
      DO IP = 1, NPTOT
        KSEGP(IP) = 0
      ENDDO
C
C---- initialize system-row and system-column counters to current count
      KSYS = NSYS
      JSYS = NSYS
C
      DO ISEG = 1, NSEG
        IPS1 = IPSEG1(ISEG)
        IPS2 = IPSEG2(ISEG)
C
C------ for nodes inside target segment, change variable from dgam to dn
        DO IP = IPS1, IPS2
          JSYSDXY(IP) = JSYSGAM(IP)
          JSYSGAM(IP) = 0
C
          KSEGP(IP) = ISEG
        ENDDO
C
        IEL = IELSEG(ISEG)
C
        IF(LNFIX1(ISEG) .OR. LNFIX2(ISEG)) THEN
C------- use constant term if at least one endpoint is fixed
         JSYS = JSYS+1
         JSYSQSP(1,ISEG) = JSYS
C
         KSYS = KSYS+1
         IF(LNFIX1(ISEG)) THEN
          KSYSDNS(1,ISEG) = KSYS
         ELSE
          KSYSDNS(2,ISEG) = KSYS
         ENDIF
        ENDIF
C
        IF(LNFIX1(ISEG) .AND. LNFIX2(ISEG)) THEN
C------- also use linear term if the other endpoint is also fixed
         JSYS = JSYS+1
         JSYSQSP(2,ISEG) = JSYS
C
         KSYS = KSYS+1
         IF(LNFIX2(ISEG)) THEN
          KSYSDNS(2,ISEG) = KSYS
         ELSE
          KSYSDNS(1,ISEG) = KSYS
         ENDIF
        ENDIF
C
        IF(LDFIX1(IEL)) THEN
C------- match slope on the left
         KSYS = KSYS+1
         JSYS = JSYS+1
         KSYSDNS(3,ISEG) = KSYS
         JSYSQSP(3,ISEG) = JSYS
        ENDIF
C
        IF(LDFIX2(IEL)) THEN
C------- match slope on the right
         KSYS = KSYS+1
         JSYS = JSYS+1
         KSYSDNS(4,ISEG) = KSYS
         JSYSQSP(4,ISEG) = JSYS
        ENDIF
C
      ENDDO
C
C---- new augmented-system size
      NSYS = KSYS
      IF(LDBG) WRITE(LUNDBG,*) '   NSYS =',NSYS,' augmented for inverse'
C
      IF(KSYS.NE.JSYS) THEN
       WRITE(*,*) 'Error in SYSPQ, KSYS JSYS =', KSYS, JSYS
      ENDIF
C
C---- standard direct-problem pointers no longer valid
      LSYSP = .FALSE.
C
      RETURN
      END ! SYSPQ


      SUBROUTINE QSYS
C-----------------------------------------------------------
C     Sets up Newton system for mixed-inverse problem.
C     The equation residuals are the same as the 
C     direct problem set up in GSYS.  The GAM variables
C     in the inverse segments ISEG = 1..NSEG are replaced
C     by DNP normal-movement variables.  
C
C     In addition, the system is augmented by 1,2,3, or 4
C     geometry-regularity conditions per inverse segment.
C
C     Unlike the direct problem, this system is nonlinear,
C     so the Jacobian matrix depends on the current solution.
C-----------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(NSYS.EQ.0) RETURN
C
      IF(LDBG) THEN
       WRITE(*,*) 'Setting up mixed-inverse Newton system...'
       WRITE(LUNDBG,*) 'QSYS mixed-inverse problem system setup'
       WRITE(LUNDBG,*) '   NSYS =', NSYS
      ENDIF
C
C---- clear system arrays
      DO K = 1, NSYS
        RES(K) = 0.
        DO J = 0, NSYS
          SYS(K,J) = 0.
        ENDDO
        DO J = 0, NAICGAM
          AICGAM(K,J) = 0.
        ENDDO
        DO J = 0, NAICSIG
          AICSIG(K,J) = 0.
        ENDDO
        AICQFF(K,1) = 0.
      ENDDO
C
C---- set flow tangency at all control points
      DO 10 IC = 1, NCTOT
        K = KSYSVNC(IC)
C
C------ skip un-enforced control points
        IF(K.LT.0) GO TO 10
C
C------- flow tangency at control point IC
         RES(K)   = ANC(1,IC)*QC(1,IC)
     &            + ANC(2,IC)*QC(2,IC)
C
         RES_QC1  = ANC(1,IC)
         RES_QC2  = ANC(2,IC)
         RES_ANC1 =           QC(1,IC)
         RES_ANC2 =           QC(2,IC)
         RES_UINF = ANC(1,IC)
C
C------- derivatives w.r.t. dn movement of the two nodes of panel IC
         JPO = IPCO(IC)
         JPP = IPCP(IC)
         RES_DX = RES_ANC1*ANC_DXY(1,1,IC) + RES_ANC2*ANC_DXY(2,1,IC)
         RES_DY = RES_ANC1*ANC_DXY(1,2,IC) + RES_ANC2*ANC_DXY(2,2,IC)
         RES_XPO = -RES_DX
         RES_YPO = -RES_DY
         RES_XPP =  RES_DX
         RES_YPP =  RES_DY
C
         JO = JSYSDXY(JPO)
         JP = JSYSDXY(JPP)
         SYS(K,JO) = SYS(K,JO) + RES_XPO*ANP(1,JPO) + RES_YPO*ANP(2,JPO)
         SYS(K,JP) = SYS(K,JP) + RES_XPP*ANP(1,JPP) + RES_YPP*ANP(2,JPP)
C
C------- derivatives w.r.t. vorticity and position of all nodes
         DO JP = 1, NPTOT
           RES_GAM = RES_QC1*QC_GAM(1,JP,IC)
     &             + RES_QC2*QC_GAM(2,JP,IC)
           RES_SIG = RES_QC1*QC_SIG(1,JP,IC)
     &             + RES_QC2*QC_SIG(2,JP,IC)
C
           RES_XP  = RES_QC1*QC_XP (1,JP,IC)
     &             + RES_QC2*QC_XP (2,JP,IC)
           RES_YP  = RES_QC1*QC_YP (1,JP,IC)
     &             + RES_QC2*QC_YP (2,JP,IC)
C
           CALL GXYLIN(K,JP,RES_GAM,RES_SIG,RES_XP,RES_YP)
         ENDDO
C
         IF(IELDOFC(IC).GT.0) THEN
C-------- flow tangency contribution of normal-velocity DOF
          IEL = IELDOFC(IC)
          RES(K) = RES(K) - QNDOF(IEL)
C
          J = JSYSQNC(IEL)
          SYS(K,J) = -1.0
         ENDIF
         IF(LDBG) WRITE(LUNDBG,*) K, 'V.n=Qn ',RES(K) !@@@
C
         AICQFF(K,1) = RES_UINF
 10   CONTINUE
C
C---- set gamma or dn-continuity at all zero-length panels
      DO 15 IC = 1, NCTOT
        K = KSYSGMG(IC)
C
C------ skip if there's no constraint
        IF(K.LT.0) GO TO 15
C
        JPO = IPCO(IC)
        JPP = IPCP(IC)
C
        IF(JSYSGAM(JPO).GT.0 .AND. JSYSGAM(JPP).GT.0) THEN
         RES(K) = GAM(JPP) - GAM(JPO)
         RES_GAMP =  1.0
         RES_GAMO = -1.0
         RES_SIGP = 0.
         RES_SIGO = 0.
         RES_XPP = 0.
         RES_XPO = 0.
         RES_YPP = 0.
         RES_YPO = 0.
        ELSE
CCC      RES(K) = DN(JPP) - DN(JPO)
         RES(K) = 0.
         RES_GAMP = 0.
         RES_GAMO = 0.
         RES_XPP = ANP(1,JPP)
         RES_XPO = ANP(1,JPO)
         RES_YPP = ANP(2,JPP)
         RES_YPO = ANP(2,JPO)
        ENDIF
        IF(LDBG) WRITE(LUNDBG,*) K, 'dGam=0 ',RES(K) !@@@
C
        CALL GXYLIN(K,JPO,RES_GAMO,RES_SIGO,RES_XPO,RES_YPO)
        CALL GXYLIN(K,JPP,RES_GAMP,RES_SIGP,RES_XPP,RES_YPP)
 15   CONTINUE
C
      DO 20 IEL = 1, NEL
C------ zero net curvature of gamma at TE (triangular loading towards TE)
        K = KSYSGSS(IEL)
        IF(K.GT.0) THEN

          JP1O = IPFRST(IEL) + 2
          JP1M = IPFRST(IEL) + 1
          JP1L = IPFRST(IEL)
          JP2O = IPLAST(IEL)
          JP2M = IPLAST(IEL) - 1
          JP2L = IPLAST(IEL) - 2
C
          DS1O = SP(JP1O) - SP(JP1M)
          DS1M = SP(JP1M) - SP(JP1L)
          DS2O = SP(JP2O) - SP(JP2M)
          DS2M = SP(JP2M) - SP(JP2L)
C
C-------- 1st derivative
c          F1O = 1.0/DS1O + 1.0/(DS1O+DS1M)
c          F1L = DS1O/((DS1O+DS1M)*DS1M)
c          F1M = -(F1O+F1L)
c          F2O = 1.0/DS2O + 1.0/(DS2O+DS2M)
c          F2L = DS2O/((DS2O+DS2M)*DS2M)
c          F2M = -(F2O+F2L)
C
C-------- 2nd derivative
          F1O =   1.0/DS1O
          F1M = -(1.0/DS1O + 1.0/DS1M)
          F1L =              1.0/DS1M
          F2O =   1.0/DS2O
          F2M = -(1.0/DS2O + 1.0/DS2M)
          F2L =              1.0/DS2M
C
          IF    (LV1ZR(IEL)) THEN
           RES(K) = -F2O*GAM(JP2O)
     &              -F2M*GAM(JP2M)
     &              -F2L*GAM(JP2L)
           CALL GXYLIN(K,JP2O,-F2O,0.0, 0.0,0.0)
           CALL GXYLIN(K,JP2M,-F2M,0.0, 0.0,0.0)
           CALL GXYLIN(K,JP2L,-F2L,0.0, 0.0,0.0)
C
          ELSEIF(LV2ZR(IEL)) THEN
           RES(K) =  F1O*GAM(JP1O)
     &              +F1M*GAM(JP1M)
     &              +F1L*GAM(JP1L)
           CALL GXYLIN(K,JP1O, F1O,0.0, 0.0,0.0)
           CALL GXYLIN(K,JP1M, F1M,0.0, 0.0,0.0)
           CALL GXYLIN(K,JP1L, F1L,0.0, 0.0,0.0)
C
          ELSE
           RES(K) =  F1O*GAM(JP1O)
     &              +F1M*GAM(JP1M)
     &              +F1L*GAM(JP1L)
     &              -F2O*GAM(JP2O)
     &              -F2M*GAM(JP2M)
     &              -F2L*GAM(JP2L)
           CALL GXYLIN(K,JP1O, F1O,0.0, 0.0,0.0)
           CALL GXYLIN(K,JP1M, F1M,0.0, 0.0,0.0)
           CALL GXYLIN(K,JP1L, F1L,0.0, 0.0,0.0)
           CALL GXYLIN(K,JP2O,-F2O,0.0, 0.0,0.0)
           CALL GXYLIN(K,JP2M,-F2M,0.0, 0.0,0.0)
           CALL GXYLIN(K,JP2L,-F2L,0.0, 0.0,0.0)
C
          ENDIF
          IF(LDBG) WRITE(LUNDBG,*) K, 'Qss=0 ',RES(K)  !@@@
        ENDIF
C
C------ Kutta condition
        K = KSYSKUT(IEL)
        IF(K.GT.0) THEN
C
         IF(LBODY(IEL)) THEN
          JP1 = IPFRST(IEL)
          JP2 = IPLAST(IEL)
          RES(K) = GAM(JP1) + GAM(JP2)
          CALL GXYLIN(K,JP1, 1.0,0.0, 0.0,0.0)
          CALL GXYLIN(K,JP2, 1.0,0.0, 0.0,0.0)
C
         ELSE
          JP = IPLAST(IEL)
          RES(K) = GAM(JP)
          CALL GXYLIN(K,JP, 1.0,0.0, 0.0,0.0)
C
         ENDIF
C
         IF(IEWAKE(IEL).NE.0) THEN
           IELW = IEWAKE(IEL)
           JP = IPFRST(IELW)
           RES(K) = RES(K) - GAM(JP)
c           RES(K) = RES(K) - GTH(JP)  !!! HHY
           CALL GXYLIN(K,JP,-1.0,0.0, 0.0,0.0)
         ENDIF
         IF(LDBG) WRITE(LUNDBG,*) K, 'Kutta ',RES(K) !@@@
        ENDIF
C
C------ normal-velocity DOF
        K = KSYSQNC(IEL)
        IF(K.GT.0) THEN
C
         J = JSYSQNC(IEL)
         RES(K) = QNDOF(IEL)
CCC  &          - QNDOFSP(IEL)
         SYS(K,J) = 1.0
         IF(LDBG) WRITE(LUNDBG,*) K, 'Qn=0 ',RES(K)  !@@@
        ENDIF
C
C------ zero internal tangential velocity constraint
        K = KSYSQTT(IEL)
        IF(K.GT.0) THEN
C
         IC = NCX + IEL
C
c            write(*,*) 'enter Qtan/Qinf'   !@@@
c            read(*,*) qttf    !@@@
            qttf = 0.0

C-------- flow-tangency residual
          RES(K) = ANC(1,IC)*QC(1,IC)
     &           + ANC(2,IC)*QC(2,IC)  -  qttf*QINF   !@@@
C
          IF(LDBG) WRITE(LUNDBG,*) K, 'Qtan=0 ',RES(K)  !@@@
          RES_QC1  = ANC(1,IC)
          RES_QC2  = ANC(2,IC)
          RES_ANC1 =           QC(1,IC)
          RES_ANC2 =           QC(2,IC)
          RES_UINF = ANC(1,IC)
C
cC-------- derivatives w.r.t. dn movement of the two nodes of panel IC
c          JPO = IPCO(IC)
c          JPP = IPCP(IC)
c          RES_DX = RES_ANC1*ANC_DXY(1,1,IC) + RES_ANC2*ANC_DXY(2,1,IC)
c          RES_DY = RES_ANC1*ANC_DXY(1,2,IC) + RES_ANC2*ANC_DXY(2,2,IC)
c          RES_XPO = -RES_DX
c          RES_YPO = -RES_DY
c          RES_XPP =  RES_DX
c          RES_YPP =  RES_DY
C
c          JO = JSYSDXY(JPO)
c          JP = JSYSDXY(JPP)
c          SYS(K,JO) = SYS(K,JO) + RES_XPO*ANP(1,JPO) + RES_YPO*ANP(2,JPO)
c          SYS(K,JP) = SYS(K,JP) + RES_XPP*ANP(1,JPP) + RES_YPP*ANP(2,JPP)
C
C-------- derivatives w.r.t. vorticity and position of all nodes
          DO JP = 1, NPTOT
            RES_GAM = RES_QC1*QC_GAM(1,JP,IC)
     &              + RES_QC2*QC_GAM(2,JP,IC)
            RES_SIG = RES_QC1*QC_SIG(1,JP,IC)
     &              + RES_QC2*QC_SIG(2,JP,IC)
C
            RES_XP  = RES_QC1*QC_XP (1,JP,IC)
     &              + RES_QC2*QC_XP (2,JP,IC)
            RES_YP  = RES_QC1*QC_YP (1,JP,IC)
     &              + RES_QC2*QC_YP (2,JP,IC)
C 
            CALL GXYLIN(K,JP,RES_GAM,RES_SIG,RES_XP,RES_YP)
          ENDDO
C
        ENDIF
 20   CONTINUE
C
C
C---- set zero vortex strengths explicitly
      DO 30 JP = 1, NPTOT
        K = KSYSGAM(JP)
C
C------ skip non-specified vortices
        IF(K.LT.0) GO TO 30
C
        RES(K) = GAM(JP)
        CALL GXYLIN(K,JP, 1.0,0.0, 0.0,0.0)
        IF(LDBG) WRITE(LUNDBG,*) K, 'gam=0 ',RES(K)  !@@@
 30   CONTINUE
C
C
C---- segment-fixing conditions
      DO ISEG = 1, NSEG
        K = KSYSDNS(1,ISEG)
        IF(K.GT.0) THEN

         JP = IPSEG1(ISEG)
CCC      RES(K) = DN(JP)
         RES(K) = 0.
         J = JSYSDXY(JP)
         SYS(K,J) = 1.0
         IF(LDBG) WRITE(LUNDBG,*) K, 'dn1=0 ',RES(K)  !@@@
        ENDIF
C
        K = KSYSDNS(2,ISEG)
        IF(K.GT.0) THEN

         JP = IPSEG2(ISEG)
CCC      RES(K) = DN(JP)
         RES(K) = 0.
         J = JSYSDXY(JP)
         SYS(K,J) = 1.0
         IF(LDBG) WRITE(LUNDBG,*) K, 'dn2=0 ',RES(K)  !@@@
        ENDIF
C
        K = KSYSDNS(3,ISEG)
        IF(K.GT.0) THEN

C
c         JP = IPSEG1(ISEG) + 1
cCCC      RES(K) = DN(JP)
c         RES(K) = 0.
c         J = JSYSDXY(JP)
c         SYS(K,J) = 1.0
C
         KQSP = KQTARG
         JPM = IPSEG1(ISEG) - 1
         JPO = IPSEG1(ISEG)
         JPP = IPSEG1(ISEG) + 1
         DSM = SSPEC(JPO) - SSPEC(JPM)
         DSP = SSPEC(JPP) - SSPEC(JPO)
         FRM =   1.0/DSM
         FRO = -(1.0/DSM + 1.0/DSP)
         FRP =             1.0/DSP
         RES(K) = FRM*GAM(JPM)
     &          + FRO*GAM(JPO)
     &          + FRP*GAM(JPP)
     &          - FRM*QSPEC(JPM,KQSP)
     &          - FRO*QSPEC(JPO,KQSP)
     &          - FRP*QSPEC(JPP,KQSP)
         CALL GXYLIN(K,JPM, FRM,0.0, 0.0,0.0)
         CALL GXYLIN(K,JPO, FRO,0.0, 0.0,0.0)
         CALL GXYLIN(K,JPP, FRP,0.0, 0.0,0.0)
         IF(LDBG) WRITE(LUNDBG,*) K, 'gam_ss1=0 ',RES(K)  !@@@
        ENDIF
C
        K = KSYSDNS(4,ISEG)
        IF(K.GT.0) THEN
C
C
c         JP = IPSEG2(ISEG) - 1
cCCC      RES(K) = DN(JP)
c         RES(K) = 0.
c         J = JSYSDXY(JP)
c         SYS(K,J) = 1.0
C
         KQSP = KQTARG
         JPM = IPSEG2(ISEG) - 1
         JPO = IPSEG2(ISEG)
         JPP = IPSEG2(ISEG) + 1
         DSM = SSPEC(JPO) - SSPEC(JPM)
         DSP = SSPEC(JPP) - SSPEC(JPO)
         FRM =   1.0/DSM
         FRO = -(1.0/DSM + 1.0/DSP)
         FRP =             1.0/DSP
         RES(K) = FRM*GAM(JPM)
     &          + FRO*GAM(JPO)
     &          + FRP*GAM(JPP)
     &          - FRM*QSPEC(JPM,KQSP)
     &          - FRO*QSPEC(JPO,KQSP)
     &          - FRP*QSPEC(JPP,KQSP)
         CALL GXYLIN(K,JPM, FRM,0.0, 0.0,0.0)
         CALL GXYLIN(K,JPO, FRO,0.0, 0.0,0.0)
         CALL GXYLIN(K,JPP, FRP,0.0, 0.0,0.0)
         IF(LDBG) WRITE(LUNDBG,*) K, 'gam_ss2=0 ',RES(K)  !@@@
        ENDIF
      ENDDO
C
      LQSYS = .TRUE.
C
      RETURN
      END ! QSYS




      SUBROUTINE GXYLIN(K,JP, RES_GAM,RES_SIG, RES_XP,RES_YP)
      INCLUDE 'DFDC.INC'
C
      J = JSYSGAM(JP)
      SYS(K,J) = SYS(K,J) + RES_GAM
C
      J = JAICGAM(JP)
      AICGAM(K,J) = AICGAM(K,J) + RES_GAM
C
      J = JAICSIG(JP)
      AICSIG(K,J) = AICSIG(K,J) + RES_SIG
C
      ISEG = KSEGP(JP)
      J1 = JSYSQSP(1,ISEG)
      J2 = JSYSQSP(2,ISEG)
      J3 = JSYSQSP(3,ISEG)
      J4 = JSYSQSP(4,ISEG)
      SYS(K,J1) = SYS(K,J1) + RES_GAM*FSPEC(1,JP)
      SYS(K,J2) = SYS(K,J2) + RES_GAM*FSPEC(2,JP)
      SYS(K,J3) = SYS(K,J3) + RES_GAM*FSPEC(3,JP)
      SYS(K,J4) = SYS(K,J4) + RES_GAM*FSPEC(4,JP)
C
      J = JSYSDXY(JP)
      SYS(K,J) = SYS(K,J) + RES_XP*ANP(1,JP)
     &                    + RES_YP*ANP(2,JP)
C
      RETURN
      END ! GXYLIN
 


      SUBROUTINE GUCALC(LGUQI,LVWK,IP1,IP2)
C-----------------------------------------------------
C     Sets up V.n residuals at all control points.
C     Calculates unit-freestream vortex sheet strengths 
C     GAMU by back-substituting into already-factored
C     Jacobian matrix.
C-----------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LGUQI, LVWK
C
      IF(LDBG) WRITE(*,*) 'Generating unit vorticity distributions...'
C
C---------------------------------------------------------------------
C---- assign unit-solution indices
      IF(LGUQI) THEN
       IF(IUQINF.EQ.0) THEN
C------ freestream x-velocity
        IF(LDBG) WRITE(*,*) 'Generating solution for unit Uinf'
        CALL GUMAKE(IUQINF,AICQFF(1,1))
       ENDIF
      ENDIF
C
      IF(LVWK) THEN
       DO IR = 1, NRP
        IF(IUVWK(IR).EQ.0) THEN
C------ rotor circulation, radial station IR
         IF(LDBG) WRITE(*,*) 'Generating solution for unit GAM IR ',IR
         J = JAICVWK(IR)
         CALL GUMAKE(IUVWK(IR),AICVWK(1,J))
         CALL GUVWK(IR)
        ENDIF
       END DO
      ENDIF
C
      DO IP = IP1, IP2
        J = JAICGAM(IP)
        IF(IUGAM(IP).EQ.0 .AND. J.GT.0) THEN
         IF(LDBG) WRITE(*,*) 'Generating solution for unit gamma',IP
         CALL GUMAKE(IUGAM(IP),AICGAM(1,J))
         GAMU(IP,IUGAM(IP)) = 1.0
        ENDIF
C
        J = JAICSIG(IP)
        IF(IUSIG(IP).EQ.0 .AND. J.GT.0) THEN
         IF(LDBG) WRITE(*,*) 'Generating solution for unit sigma',IP
         CALL GUMAKE(IUSIG(IP),AICSIG(1,J))
         SIGU(IP,IUSIG(IP)) = 1.0
        ENDIF
      ENDDO
C
      RETURN
      END ! GUCALC



      SUBROUTINE GUMAKE(IU,AIC)
C-------------------------------------------------------------------
C     Finds the index IU of the first available unit r.h.s. vector
C     Generates GAMU(.IU),SIGU(.IU),QNDOFU(.IU) unit solutions.
C-------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      REAL AIC(*)
C
C---- find first available vector index KU
      DO KU = 1, NU
        IF(.NOT.LUSET(KU)) GO TO 10
      ENDDO
C
C---- no free vectors found... add new one to list
      NU = NU + 1
      IF(NU.GT.NUX) THEN
       WRITE(*,*) 'GUMAKE: Array limit exceeded. Increase NUX to', NU
      ENDIF
C
      KU = MIN( NU , NUX )
C
C-----------------------------------
 10   CONTINUE
C
C---- set Jacobian vector
      DO J = 1, NSYS
        RES(J) = AIC(J)
      ENDDO
C
C---- calculate response
      IF(NSYS.GT.0) THEN
        CALL BAKSUB(NSYSX,NSYS,SYS(1,1),RES(1))
      ENDIF
      RES(0) = 0.
C
C---- first clear entire response vectors
      DO IP = 1, NPTOT
        GAMU(IP,KU) = 0.
        SIGU(IP,KU) = 0.
      ENDDO
      DO IEL = 1, NEL
        QNDOFU(IEL,KU) = 0.
      ENDDO
C
C---- store response of surface gamma
      DO IP = 1, NPTOT
        J = JSYSGAM(IP)
        GAMU(IP,KU) = -RES(J)
      ENDDO
C
C---- store response of normal-velocity DOFs
      DO IEL = 1, NEL
        J = JSYSQNC(IEL)
        QNDOFU(IEL,KU) = -RES(J)
      ENDDO
C
      LUSET(KU) = .TRUE.
C
      IU = KU
C
      RETURN
      END ! GUMAKE



      SUBROUTINE GSOLVE
C-------------------------------------------------------------------
C     Does a direct solve from current RHS (knowns)
C     Generates GAM(.),SIG(.),QNDOF(.)
C-------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION GAMK(IPX),SIGK(IPX)
C
C---- set imposed (known) gamma and sigma strengths
      DO IP = 1, NPTOT
        GAMK(IP) = GAMVSP(IP) 
        SIGK(IP) = SIGVSP(IP) 
      ENDDO
C
      DO IEL = 1, NEL
        DO IP = IPFRST(IEL), IPLAST(IEL)
          GAMK(IP) = GAMK(IP) + GAMSET(IEL)
          SIGK(IP) = SIGK(IP) + SIGSET(IEL)
        ENDDO
      ENDDO
C
C---- set Jacobian vector (RHS)
      DO K = 1, NSYS
        RES(K) = AICQFF(K,1)*QINF
        DO IP = 1, NPTOT
          J = JAICGTH(IP)
          RES(K) = RES(K) + AICGTH(K,J)*GTH(IP)
        ENDDO
        DO IP = 1, NPTOT
          J = JAICGAM(IP)
          RES(K) = RES(K) + AICGAM(K,J)*GAMK(IP)
        ENDDO
        DO IP = 1, NPTOT
          J = JAICSIG(IP)
          RES(K) = RES(K) + AICSIG(K,J)*SIGK(IP)
        ENDDO
      ENDDO
C
C---- calculate solution
      IF(NSYS.GT.0) THEN
        CALL BAKSUB(NSYSX,NSYS,SYS(1,1),RES(1))
      ELSE
        WRITE(*,*) 'GSOLVE error NSYS=',NSYS
      ENDIF
      RES(0) = 0.
C
C---- first clear unknown solution vectors
      DO IP = 1, NPTOT
        GAM(IP) = 0.
        SIG(IP) = 0.
      ENDDO
      DO IEL = 1, NEL
        QNDOF(IEL) = 0.
      ENDDO
C
C---- store solution for surface gamma + known gamma values
C---- store solution for surface sigma + known sigma values
      DO IP = 1, NPTOT
        J = JSYSGAM(IP)
        GAM(IP) = -RES(J) + GAMK(IP)
        SIG(IP) = SIG(IP) + SIGK(IP)
      ENDDO
C
C---- store solution for normal-velocity DOFs
      DO IEL = 1, NEL
        J = JSYSQNC(IEL)
        QNDOF(IEL) = -RES(J)
      ENDDO
C
C---- Set strength of TE panels (for debugging use)
      CALL SETGSTE
C
      LGAMA = .TRUE.
      LQCNT = .FALSE.
C
      RETURN
      END ! GSOLVE


      SUBROUTINE GSOLVE0
C-------------------------------------------------------------------
C     Does a direct solve from current RHS (knowns)
C     Generates GAM(.),SIG(.),QNDOF(.)
C-------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION GAMK(IPX),SIGK(IPX)
C
C---- set known gamma and sigma strengths
      DO IP = 1, NPTOT
        GAMK(IP) = GAMVSP(IP)
        SIGK(IP) = SIGVSP(IP)
      ENDDO
C
      DO IEL = 1, NEL
        DO IP = IPFRST(IEL), IPLAST(IEL)
          GAMK(IP) = GAMK(IP) + GAMSET(IEL)
          SIGK(IP) = SIGK(IP) + SIGSET(IEL)
        ENDDO
      ENDDO
C
C---- set Jacobian vector (RHS)
      DO K = 1, NSYS
        RES(K) = AICQFF(K,1)*QINF
c        DO IR = 1, NRP
c          RES(K) = RES(K) + AICVWK(K,IR)*WBGAM(IR)
c        ENDDO
        DO IP = 1, NPTOT
          J = JAICGAM(IP)
          RES(K) = RES(K) + AICGAM(K,J)*GAMK(IP)
        ENDDO
        DO IP = 1, NPTOT
          J = JAICSIG(IP)
          RES(K) = RES(K) + AICSIG(K,J)*SIGK(IP)
        ENDDO
      ENDDO
C
C---- calculate response
      IF(NSYS.GT.0) THEN
        CALL BAKSUB(NSYSX,NSYS,SYS(1,1),RES(1))
      ENDIF
      RES(0) = 0.
C
C---- first clear unknown solution vectors
      DO IP = 1, NPTOT
        GAM(IP) = 0.
        GTH(IP) = 0.
        SIG(IP) = 0.
      ENDDO
      DO IEL = 1, NEL
        QNDOF(IEL) = 0.
      ENDDO
C
C---- store solution for surface gamma + known gamma values
C---- store solution for surface sigma + known sigma values
      DO IP = 1, NPTOT
        J = JSYSGAM(IP)
        GAM(IP) = -RES(J) + GAMK(IP)
        SIG(IP) = SIG(IP) + SIGK(IP)
C---- set GTH for wake point from WBGAM for wake row
c        IF(IP2IR(IP).NE.0) THEN
c          IR = IP2IR(IP)
c          GTH(IP) = WBGAM(IR)/YP(IP)
c        ENDIF
      ENDDO
C
C---- store solution for normal-velocity DOFs
      DO IEL = 1, NEL
        J = JSYSQNC(IEL)
        QNDOF(IEL) = -RES(J)
      ENDDO
C
C---- Set strength of TE panels (for debugging use)
      CALL SETGSTE
      LGAMA = .TRUE.
      LQCNT = .FALSE.
C
      RETURN
      END ! GSOLVE0


      SUBROUTINE GUVWK(IRADD)
C-------------------------------------------------------------------
C     Creates wake gamma "mode" shape corresponding to unit
C     rotor row IRADD circulation (actually Wt/Wm *delta BGAM /2*pi)
C
C     Creates GTHU(.IU) unit vorticity distribution
C-------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
cc      WRITE(*,*) 'GAMADDVWK: add rotor wake gammas to unit soln ',IRADD
C
C---- store imposed unit rotor WBGAM circulation gamma distribution (on wake
C     and body nodes) as GTHU
      DO JP = 1, NPTOT
C---- clear unit response 
        KU = IUVWK(IRADD)
        GTHU(JP,KU) = 0.
C
        IR = IP2IR(JP)
        IF(IR.NE.0) THEN
C---- we are on a rotor wake element, use rotor index to locate influence
         IF(IR.EQ.IRADD) THEN 
           IF(IUVWK(IR).GT.0) THEN
             KU = IUVWK(IR)
             GTHU(JP,KU) = 1.0 / YP(JP)
           ELSE
             WRITE(*,*) 'GUVWK: Unassigned rotor RHS IR =',IR
           ENDIF
         ENDIF
C
        ENDIF
      ENDDO
C
      RETURN
      END ! GUVWK



      SUBROUTINE GUSUM
C-----------------------------------------------------
C     Superimposes unit solutions for current...
C       * freestream Qinf
C       * rotor circulations
C       * specified wake  gamma,sigma
C       * specified ring  Gamma,Sigma
C       * specified point Doublet,Source
C       * specified source-only elements
C-----------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION UVAL(0:NUX)
C
      IF(LDBG) WRITE(*,*) 'Superimposing unit-strength distributions'
C
      DO IU = 0, NU
        UVAL(IU) = 0.
      ENDDO
C
C---- put prescribed strengths into indexed array (inactive quantity has IU=0)
C
C----- axisymmetric case
      UVAL(IUQINF) = QINF
C
C---- Rotor wake circulation (actually Wt/Wm *delta BGAM /2*pi)
c      DO IR = 1, NRP
c        UVAL(IUVWK(IR)) = WBGAM(IR)
cc       write(18,*) ir,iuvwk(ir),wbgam(ir)
c      ENDDO
C
C---- Additional gamma specified from viscous response (from BL effects)
      DO IP = 1, NPTOT
cc        write(24,*) IP,IUGAM(IP),GAMVSP(IP)
        UVAL(IUGAM(IP)) = GAMVSP(IP)
        UVAL(IUSIG(IP)) = SIGVSP(IP)
      ENDDO
C
      DO IEL = 1, NEL
        DO IP = IPFRST(IEL), IPLAST(IEL)
          UVAL(IUGAM(IP)) = UVAL(IUGAM(IP)) + GAMSET(IEL)
          UVAL(IUSIG(IP)) = UVAL(IUSIG(IP)) + SIGSET(IEL)
        ENDDO
      ENDDO
C
C---- clear all singularities for accumulation
      DO IP = 1, NPTOT
        GAM(IP) = 0.
        SIG(IP) = 0.
        GTH(IP) = 0.
      ENDDO
      DO IEL = 1, NEL
        QNDOF(IEL) = 0.
      ENDDO
C
C---- accumulate unit contributions weighted by prescribed strengths
      DO 100 IU = 1, NU
C------ no sense accumulating zeros
        IF(UVAL(IU) .EQ. 0.0) GO TO 100
cc          write(20,*) iu,uval(iu)
C
        DO IP = 1, NPTOT
          GAM(IP) = GAM(IP) + GAMU(IP,IU)*UVAL(IU)
          SIG(IP) = SIG(IP) + SIGU(IP,IU)*UVAL(IU)
          GTH(IP) = GTH(IP) + GTHU(IP,IU)*UVAL(IU)
cc          write(20,99) ip,iu,gamu(ip,iu),GTHU(IP,IU)
        ENDDO
        DO IEL = 1, NEL
          QNDOF(IEL) = QNDOF(IEL) + QNDOFU(IEL,IU)*UVAL(IU)
        ENDDO
C
 100  CONTINUE
C
c      DO IP = 1, NPTOT
c        write(24,*) IP,GAM(IP),GAMVSP(IP)
c      ENDDO
C
C---- Set strength of TE panels (for debugging use)
      CALL SETGSTE
C
 99   format(2i5,5(1x,F12.6))
C
      LGAMA = .TRUE.
      LQCNT = .FALSE.
C
      RETURN
      END ! GUSUM



      SUBROUTINE QCSUM
C-----------------------------------------------------
C     Computes velocities at control points for
C     current panel and point-singularity strengths.
C-----------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LANS
      DIMENSION QCT(2,ICX)
C
      IF(LDBG) WRITE(*,*) 'Computing control-point velocities'
C
      DO 20 IC = 1, NCTOT
        IF    (ICTYPE(IC).EQ.-999) THEN
         WRITE(*,*) '? unassigned control point ', IC
C
        ELSEIF(ICTYPE(IC).EQ.-1) THEN
C------- zero-length panel (skip)
C
        ELSEIF(ICTYPE(IC).EQ.-2) THEN
C------- vortex wake panels w/o control points (skip)
C
        ELSE
C------- usual velocity contributions for all vortex and source strengths
         QC(1,IC) = 0.
         QC(2,IC) = 0.
C
         DO JP=1, NPTOT
           QC(1,IC) = QC(1,IC) + QC_GAM(1,JP,IC)*GAM(JP)
     &                         + QC_SIG(1,JP,IC)*SIG(JP)
           QC(2,IC) = QC(2,IC) + QC_GAM(2,JP,IC)*GAM(JP)
     &                         + QC_SIG(2,JP,IC)*SIG(JP)
C---- Add contribution from rotor wake vorticity
           QC(1,IC) = QC(1,IC) + QC_GTH(1,JP,IC)*GTH(JP)
           QC(2,IC) = QC(2,IC) + QC_GTH(2,JP,IC)*GTH(JP)
         ENDDO
C
C---- Add freestream flow contribution
         QC(1,IC) = QC(1,IC) + QINF
        ENDIF
C
 20   CONTINUE
C
      LQCNT = .TRUE.
C
      RETURN
      END ! QCSUM




      SUBROUTINE QCUSET
C-----------------------------------------------------
C     Computes velocities at control points for
C     current unit-forcing solutions
C-----------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) WRITE(*,*) 'Computing control-pt unit-forcing velocities'
C
      DO 20 IC = 1, NCTOT
        IF    (ICTYPE(IC).EQ.-999) THEN
         WRITE(*,*) '? unassigned control point ', IC
C
        ELSEIF(ICTYPE(IC).EQ.-1) THEN
C------- zero-length panel (skip)
C
        ELSEIF(ICTYPE(IC).EQ.-2) THEN
C------- vortex wake panel w/o control point (skip)
C
        ELSE
C------- usual velocity contributions for all vortex strengths
         DO IU=1, NU
           QCU(1,IC,IU) = 0.
           QCU(2,IC,IU) = 0.
C
           DO JP=1, NPTOT
cc             IF(JSYSGAM(JP).GT.0) THEN
C---- Contribution from solution unit responses
              QCU(1,IC,IU) = QCU(1,IC,IU) + QC_GAM(1,JP,IC)*GAMU(JP,IU)
              QCU(2,IC,IU) = QCU(2,IC,IU) + QC_GAM(2,JP,IC)*GAMU(JP,IU)
C---- Contribution from rotor wake vorticity unit forcing
              QCU(1,IC,IU) = QCU(1,IC,IU) + QC_GTH(1,JP,IC)*GTHU(JP,IU)
              QCU(2,IC,IU) = QCU(2,IC,IU) + QC_GTH(2,JP,IC)*GTHU(JP,IU)
cc             ENDIF
           ENDDO
         ENDDO
C
C---- Contribution from freestream
         IU = IUQINF
         QCU(1,IC,IU) = QCU(1,IC,IU) + 1.0
C
C---- Contribution from specified vortex and source panels
         DO IP = 1, NPTOT
           IU = IUGAM(IP)
           QCU(1,IC,IU) = QCU(1,IC,IU) + QC_GAM(1,IP,IC)
           QCU(2,IC,IU) = QCU(2,IC,IU) + QC_GAM(2,IP,IC)
C
           IU = IUSIG(IP)
           QCU(1,IC,IU) = QCU(1,IC,IU) + QC_SIG(1,IP,IC)
           QCU(2,IC,IU) = QCU(2,IC,IU) + QC_SIG(2,IP,IC)
         ENDDO
C
        ENDIF
 20   CONTINUE
C
      RETURN
      END ! QCUSET



      SUBROUTINE QCUSUM
C-----------------------------------------------------
C     Computes velocities at control points for
C     current panel and point-singularity strengths.
C-----------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION UVAL(0:NUX)
C
      IF(LDBG) WRITE(*,*) 'Superimposing unit-strength velocities'
C
      DO IU = 0, NU
        UVAL(IU) = 0.
      ENDDO
C
C---- put prescribed strengths into indexed array (inactive quantity has IU=0)
C
C----- freestream 
      UVAL(IUQINF) = QINF
C
C---- Rotor wake circulation (actually Wt/Wm *delta BGAM /2*pi)
c      DO IR = 1, NRP
c        UVAL(IUVWK(IR)) = WBGAM(IR)
cc       write(18,*) ir,iuvwk(ir),wbgam(ir)
c      ENDDO
C
      DO IP = 1, NPTOT
        UVAL(IUGAM(IP)) = GAMVSP(IP)
        UVAL(IUSIG(IP)) = SIGVSP(IP)
      ENDDO
C---- Specified singularities
      DO IEL = 1, NEL
        DO IP = IPFRST(IEL), IPLAST(IEL)
          UVAL(IUGAM(IP)) = UVAL(IUGAM(IP)) + GAMSET(IEL)
          UVAL(IUSIG(IP)) = UVAL(IUSIG(IP)) + SIGSET(IEL)
        ENDDO
      ENDDO
C
C---- accumulate unit contributions weighted by prescribed strengths
      DO IC = 1, NCTOT
        QC(1,IC) = 0.
        QC(2,IC) = 0.
      ENDDO
C
      DO 100 IU = 1, NU
C------ no sense accumulating zeros
        IF(UVAL(IU) .EQ. 0.0) GO TO 100
C
        DO IC = 1, NCTOT
          QC(1,IC) = QC(1,IC) + QCU(1,IC,IU)*UVAL(IU)
          QC(2,IC) = QC(2,IC) + QCU(2,IC,IU)*UVAL(IU)
        ENDDO
C
 100  CONTINUE
C
      LQCNT = .TRUE.
C
      RETURN
      END ! QCUSUM




      SUBROUTINE SETGSTE
C-----------------------------------------------------
C     Sets TE vortex and source strength 
C     from vortex and source strength on body panels 
C-----------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) WRITE(*,*) 'Setting TE GAMT,SIGT'
C
C---- Set TE panel strengths
cc      write(30,*) 'SETGSTE IEL,GAMT1-2,SIGT1-2'
C
      DO IEL = 1, NEL
        DO K = 1, 2
          GAMT(K,IEL) = 0.
          SIGT(K,IEL) = 0.
        ENDDO
C
C---- Body elements with GAM
        IF(LTPAN(IEL) .AND. NETYPE(IEL).EQ.0) THEN
          IP1 = IPTE1(IEL)
          IP2 = IPTE2(IEL)
C
          IF(IP1.NE.0) THEN
           DO K = 1, 2
             GAMT(K,IEL) = GAMT(K,IEL) + GAMT_GAM(K,1,IEL)*GAM(IP1)
     &                                 + GAMT_SIG(K,1,IEL)*SIG(IP1)
             SIGT(K,IEL) = SIGT(K,IEL) + SIGT_GAM(K,1,IEL)*GAM(IP1)
     &                                 + SIGT_SIG(K,1,IEL)*SIG(IP1)
           ENDDO
          ENDIF
          IF(IP2.NE.0) THEN
           DO K = 1, 2
             GAMT(K,IEL) = GAMT(K,IEL) + GAMT_GAM(K,2,IEL)*GAM(IP2)
     &                                 + GAMT_SIG(K,2,IEL)*SIG(IP2)
             SIGT(K,IEL) = SIGT(K,IEL) + SIGT_GAM(K,2,IEL)*GAM(IP2)
     &                                 + SIGT_SIG(K,2,IEL)*SIG(IP2)
           ENDDO
          ENDIF
        ENDIF
C
C---- vortex wake elements with GTH specified
        IF(LTPAN(IEL) .AND. NETYPE(IEL).EQ.7) THEN
          IP1 = IPTE1(IEL)
          IP2 = IPTE2(IEL)
C
          IF(IP1.NE.0) THEN
           DO K = 1, 2
             GAMT(K,IEL) = GAMT(K,IEL) + GAMT_GAM(K,1,IEL)*GTH(IP1)
             SIGT(K,IEL) = SIGT(K,IEL) + SIGT_GAM(K,1,IEL)*GTH(IP1)
           ENDDO
          ENDIF
          IF(IP2.NE.0) THEN
           DO K = 1, 2
             GAMT(K,IEL) = GAMT(K,IEL) + GAMT_GAM(K,2,IEL)*GTH(IP2)
             SIGT(K,IEL) = SIGT(K,IEL) + SIGT_GAM(K,2,IEL)*GTH(IP2)
           ENDDO
          ENDIF
        ENDIF
C
cc        write(30,99) iel,GAMT(1,IEL),GAMT(2,IEL),SIGT(1,IEL),SIGT(2,IEL)
 99     format(i5,4(1x,f12.6)) 
      ENDDO
C
      RETURN
      END ! SETGSTE



      SUBROUTINE NWDOTS(IEL,DS1NW,DS2NW,DN1NW,DN2NW)
      INCLUDE 'DFDC.INC'
      REAL ANW(2)
C
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
C
c      IELW = IEWAKE(IEL)
c      IF(IELW .NE. 0) THEN
c       ICW = ICFRST(IELW)
c       ANW(1) = ANC(1,ICW)
c       ANW(2) = ANC(2,ICW)
c      ELSE
       ANW(1) = ANC(1,IC2) - ANC(1,IC1)
       ANW(2) = ANC(2,IC2) - ANC(2,IC1)
       ANWABS = SQRT(ANW(1)**2 + ANW(2)**2)
       ANW(1) = ANW(1)/ANWABS
       ANW(2) = ANW(2)/ANWABS
c      ENDIF
C
      DS1NW = ANC(2,IC1)*ANW(1) - ANC(1,IC1)*ANW(2)
      DS2NW = ANC(2,IC2)*ANW(1) - ANC(1,IC2)*ANW(2)
C
      DN1NW = ANC(1,IC1)*ANW(1) + ANC(2,IC1)*ANW(2)
      DN2NW = ANC(1,IC2)*ANW(1) + ANC(2,IC2)*ANW(2)
C
      RETURN
      END ! NWDOTS


      SUBROUTINE SWDOTS(IEL,DS1SW,DS2SW,DN1SW,DN2SW)
      INCLUDE 'DFDC.INC'
      REAL ANW(2)
C
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
C
c      IELW = IEWAKE(IEL)
c      IF(IELW .NE. 0) THEN
c       ICW = ICFRST(IELW)
c       ANW(1) = ANC(1,ICW)
c       ANW(2) = ANC(2,ICW)
c      ELSE
       ANW(1) = ANC(1,IC2) - ANC(1,IC1)
       ANW(2) = ANC(2,IC2) - ANC(2,IC1)
       ANWABS = SQRT(ANW(1)**2 + ANW(2)**2)
       ANW(1) = ANW(1)/ANWABS
       ANW(2) = ANW(2)/ANWABS
c      ENDIF
C
      DS1SW = ANC(2,IC1)*ANW(2) + ANC(1,IC1)*ANW(1)
      DS2SW = ANC(2,IC2)*ANW(2) + ANC(1,IC2)*ANW(1)
C
      DN1SW = ANC(1,IC1)*ANW(2) - ANC(2,IC1)*ANW(1)
      DN2SW = ANC(1,IC2)*ANW(2) - ANC(2,IC2)*ANW(1)
C
      RETURN
      END ! SWDOTS


