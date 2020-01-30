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

      SUBROUTINE XYPSPL(IEL)
      INCLUDE 'DFDC.INC'
C---------------------------------------------
C     Splines panel-node coordinates 
C     and sets other element-related stuff.
C---------------------------------------------
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
      I = IP1
      N = IP2 - IP1 + 1
      IF(N.GT.1) THEN
       CALL SCALC(XP(I),YP(I),SP(I),N)
       CALL SEGSPL(XP(I),XPS(I),SP(I),N)
       CALL SEGSPL(YP(I),YPS(I),SP(I),N)
      ENDIF
C
      IF(LXBOD(IEL)) THEN
C----- axisymmetric body on axis
       IF(XP(IP1).LT.XP(IP2)) THEN
        SPLE(IEL) = SP(IP1)
        XPLE(IEL) = XP(IP1)
        YPLE(IEL) = YP(IP1)
        XPTE(IEL) = XP(IP2)
        YPTE(IEL) = YP(IP2)
       ELSE
        SPLE(IEL) = SP(IP2)
        XPLE(IEL) = XP(IP2)
        YPLE(IEL) = YP(IP2)
        XPTE(IEL) = XP(IP1)
        YPTE(IEL) = YP(IP1)
       ENDIF
C
      ELSEIF(LBODY(IEL) .AND. (LV1ZR(IEL) .OR. LV2ZR(IEL))) THEN
C----- axisymmetric body not closed on axis 
       IF(XP(IP1).LT.XP(IP2)) THEN
        SPLE(IEL) = SP(IP1)
        XPLE(IEL) = XP(IP1)
        YPLE(IEL) = YP(IP1)
        XPTE(IEL) = XP(IP2)
        YPTE(IEL) = YP(IP2)
       ELSE
        SPLE(IEL) = SP(IP2)
        XPLE(IEL) = XP(IP2)
        YPLE(IEL) = YP(IP2)
        XPTE(IEL) = XP(IP1)
        YPTE(IEL) = YP(IP1)
       ENDIF
C
      ELSEIF(LBODY(IEL) .AND. .NOT.LXBOD(IEL)) THEN
C----- body off axis
       CALL LEFIND(SPLE(IEL),XP(I),XPS(I),YP(I),YPS(I),SP(I),N)
       XPLE(IEL) = SEVAL(SPLE(IEL),XP(I),XPS(I),SP(I),N)
       YPLE(IEL) = SEVAL(SPLE(IEL),YP(I),YPS(I),SP(I),N)
       XPTE(IEL) = 0.5*(XP(IP1)+XP(IP2))
       YPTE(IEL) = 0.5*(YP(IP1)+YP(IP2))
C
      ELSE
C----- zero-thickness body
       SPLE(IEL) = SP(IP1)
       XPLE(IEL) = XP(IP1)
       YPLE(IEL) = YP(IP1)
       XPTE(IEL) = XP(IP2)
       YPTE(IEL) = YP(IP2)
      ENDIF
C
C---- set location for plotted element index
      IF(NETYPE(IEL).EQ.0 .OR.
     &   NETYPE(IEL).EQ.1 .OR.
     &   NETYPE(IEL).EQ.2     ) THEN
C----- surface or axis line... set plot location at skin centroid
       XELNUM(IEL) = XPCENT(IEL)
       YELNUM(IEL) = YPCENT(IEL)
      ELSE
C----- point,ring,source line or vortex wake element... 
C      just set plot location at the point itself
       XELNUM(IEL) = XP(IP1)
       YELNUM(IEL) = YP(IP1)
      ENDIF
C
      RETURN
      END ! XYPSPL




      SUBROUTINE CVPGEN
C---------------------------------------------------------
C     Sets control-point locations on paneled airfoil.
C     Sets locations of point sources.
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) WRITE(*,*) 'Setting control points for elements'
C
      DO IEL=1, NEX
        ICFRST(IEL) = 0
        ICLAST(IEL) = 0
      ENDDO
C
      DO IC=1, ICX
        IPCO(IC) = 0
        IPCP(IC) = 0
      ENDDO
C
C
C---- initialize counters for pointer accumulation
      IC = 0
C
      DO 100 IEL=1, NEL
        IF( NETYPE(IEL).EQ.3 .OR.
     &      NETYPE(IEL).EQ.4     ) THEN
C------- ring or point singularity... no control points
         IP = IPFRST(IEL)
C
C------- set limits so as to skip over IC do-loops
         ICFRST(IEL) = 0
         ICLAST(IEL) = -1
C
         GO TO 100
        ENDIF
C
C------ first control point in element IEL
        ICFRST(IEL) = IC+1
C
C------ go over all panels on element IEL
        DO 50 IP = IPFRST(IEL), IPLAST(IEL)-1
          IC = IC+1
          IF(IC.GT.NCX) STOP 'CVPGEN: Array overflow on NCX'
C
          DXP = XP(IP+1) - XP(IP)
          DYP = YP(IP+1) - YP(IP)
          DSP = SQRT(DXP**2 + DYP**2)
C
          IF(DSP.EQ.0.0) THEN
C--------- zero-length panel gets dummy control point
           ICTYPE(IC) = -1
C
          ELSE
C--------- ordinary panel
           ICTYPE(IC) = 0
           IF(IEL2IR(IEL).EQ.1) THEN
             ICTYPE(IC) = 0
           ELSEIF(IEL2IR(IEL).GT.1 .AND. IEL2IR(IEL).LT.NRP) THEN
cc             ICTYPE(IC) = -2
           ELSEIF(IEL2IR(IEL).EQ.NRP) THEN
             ICTYPE(IC) = 0
           ENDIF
C
          ENDIF
C
C-------- panel nodes defining sheet strength at IC
          IPCO(IC) = IP
          IPCP(IC) = IP+1
 50     CONTINUE
C
C------ last control point in element IEL
        ICLAST(IEL) = IC
C
 100  CONTINUE
C
C---- total number of control points
      NCTOT = IC
C
C---- compute normal vectors at vortex nodes, compute geometric quantities
      DO IEL=1, NEL
C------ don't do point singularities
        IF( NETYPE(IEL).EQ.0 .OR.
     &      NETYPE(IEL).EQ.1 .OR.
     &      NETYPE(IEL).EQ.2 .OR.
     &      NETYPE(IEL).EQ.5 .OR.
     &      NETYPE(IEL).EQ.6 .OR.
     &      NETYPE(IEL).EQ.7     )  THEN
         CALL XYCSET(IEL)
         CALL ANPSET(IEL)
        ENDIF
      ENDDO
C
C---- invalidate any existing solution
      LQAIC = .FALSE.
      LSYSP = .FALSE.
      LGSYS = .FALSE.
      LGAMU = .FALSE.
      LGAMA = .FALSE.
C
      LNCVP = .TRUE.
C
      RETURN
      END ! CVPGEN



      SUBROUTINE XYCSET(IEL)
C----------------------------------------------------------------
C     Sets control points and their unit normal vectors
C     for element IEL.
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
cc    DATA DCFRAC / 0.1 /
      DATA DCFRAC / 0.05 /
cc    DATA DCFRAC / 0.02 /
C
      DO IC = ICFRST(IEL), ICLAST(IEL)
        IPO = IPCO(IC)
        IPP = IPCP(IC)
C
        DXP = XP(IPP) - XP(IPO)
        DYP = YP(IPP) - YP(IPO)
        DSP = SQRT(DXP**2 + DYP**2)
C
        IF(DSP.EQ.0.0) THEN
C------- zero-length panel
         DSPINV = 0.0
        ELSE
C------- ordinary panel
         DSPINV = 1.0/DSP
        ENDIF
C
C------ assign control point at center of panel
        XC(IC) = XP(IPO) + 0.5*DXP
        YC(IC) = YP(IPO) + 0.5*DYP
        DSC(IC) = DSP
        DSC_DXY(1,IC) =  DXP * DSPINV
        DSC_DXY(2,IC) =  DYP * DSPINV
C
C------ normal vector at control point
        ANC(1,IC) =  DYP*DSPINV
        ANC(2,IC) = -DXP*DSPINV
C
C------ geometry sensitivities of normal vector for design/optimization stuff
        ANC_DXY(1,1,IC) =  ANC(1,IC)*ANC(2,IC) * DSPINV
        ANC_DXY(1,2,IC) =  ANC(2,IC)**2        * DSPINV
        ANC_DXY(2,1,IC) = -ANC(1,IC)**2        * DSPINV
        ANC_DXY(2,2,IC) = -ANC(2,IC)*ANC(1,IC) * DSPINV
      ENDDO
C
C
C------ set up extra QNDOF control points for each body element
      IF(.NOT.LBODY(IEL)) THEN
       ICE = NCX + IEL
       XC(ICE) = 0.
       YC(ICE) = 0.
C
      ELSEIF(LXBOD(IEL)) THEN
       IC = ICLAST(IEL)
       IP = IPLAST(IEL)
C
       ICE = NCX + IEL
       XC(ICE) = XP(IP) - DCFRAC*DSC(IC)
       YC(ICE) = 0.
       ANC(1,ICE) = 1.0
       ANC(2,ICE) = 0.
C
      ELSE
       IC1 = ICFRST(IEL)
       IC2 = ICLAST(IEL)
C
       IP1 = IPFRST(IEL)
       IP2 = IPLAST(IEL)
C
       XTAN = XP(IP1) - XP(IP1+1) + XP(IP2) - XP(IP2-1)
       YTAN = YP(IP1) - YP(IP1+1) + YP(IP2) - YP(IP2-1)
       STAN = SQRT(XTAN**2 + YTAN**2)
C
       DSCA = 0.5*(DSC(IC1) + DSC(IC2))
C
       ICE = NCX + IEL
       XC(ICE) = 0.5*(XP(IP1)+XP(IP2)) - DCFRAC*DSCA*XTAN/STAN
       YC(ICE) = 0.5*(YP(IP1)+YP(IP2)) - DCFRAC*DSCA*YTAN/STAN
       ANC(1,ICE) = XTAN/STAN
       ANC(2,ICE) = YTAN/STAN
C
      ENDIF
C
      RETURN
      END ! XYCSET



      SUBROUTINE ANPSET(IEL)
C----------------------------------------------------------------
C     Sets unit normal vectors on vortex node locations
C     for element IEL.
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- don't process single-point elements
      IF(IPFRST(IEL).EQ.IPLAST(IEL)) RETURN
C
      DO 10 IP = IPFRST(IEL), IPLAST(IEL)
C
C------ examine the two panels IPM..IP, IP..IPP adjoining vortex node IP
        IPM = MAX( IP-1 , IPFRST(IEL) )
        IPP = MIN( IP+1 , IPLAST(IEL) )
C
        DXM = XP(IP) - XP(IPM)
        DYM = YP(IP) - YP(IPM)
C
        DXP = XP(IPP) - XP(IP)
        DYP = YP(IPP) - YP(IP)
C
        DSQM = DXM*DXM + DYM*DYM
        DSQP = DXP*DXP + DYP*DYP
C
        IF    (DSQM.EQ.0.0) THEN
C------- first point on element, or panel IPM..IP has zero length
         DX = DXP
         DY = DYP
         DX_DXM = 0.
         DX_DYM = 0.
         DX_DXP = 1.0
         DX_DYP = 0.
         DY_DXM = 0.
         DY_DYM = 0.
         DY_DXP = 0.
         DY_DYP = 1.0
        ELSEIF(DSQP.EQ.0.0) THEN
C------- last point on element, or panel IP..IPP has zero length
         DX = DXM
         DY = DYM
         DX_DXM = 1.0
         DX_DYM = 0.
         DX_DXP = 0.
         DX_DYP = 0.
         DY_DXM = 0.
         DY_DYM = 1.0
         DY_DXP = 0.
         DY_DYP = 0.
        ELSE
C------- usual interior point... normal vector is weighted average
         DX = DXM*DSQP + DXP*DSQM
         DY = DYM*DSQP + DYP*DSQM
         DX_DXM = DSQP            + DXP*2.0*DXM
         DX_DYM =                 + DXP*2.0*DYM
         DX_DXP =            DSQM + DXM*2.0*DXP
         DX_DYP =                 + DXM*2.0*DYP
         DY_DXM =                 + DYP*2.0*DXM
         DY_DYM = DSQP            + DYP*2.0*DYM
         DY_DXP =                 + DYM*2.0*DXP
         DY_DYP =            DSQM + DYM*2.0*DYP
        ENDIF
C
        DSQ = DX*DX + DY*DY
        IF(DSQ.EQ.0.0) THEN
C------- internal error
         WRITE(*,*) '? ANPSET: zero avg panel length. IEL IP =',IEL,IP
         RETURN
        ENDIF
C
C------ set unit normal ANP
        DS = SQRT(DSQ)
        ANP(1,IP) =  DY/DS
        ANP(2,IP) = -DX/DS
C
        TMP1 = -ANP(1,IP)/DSQ
        ANP_XYM(1,1,IP) = - DY_DXM/DS - TMP1*(DX*DX_DXM + DY*DY_DXM)
        ANP_XYM(1,2,IP) = - DY_DYM/DS - TMP1*(DX*DX_DYM + DY*DY_DYM)
C
        ANP_XYO(1,1,IP) =   DY_DXM/DS + TMP1*(DX*DX_DXM + DY*DY_DXM)
     &                    - DY_DXP/DS - TMP1*(DX*DX_DXP + DY*DY_DXP)
        ANP_XYO(1,2,IP) =   DY_DYM/DS + TMP1*(DX*DX_DYM + DY*DY_DYM)
     &                    - DY_DYP/DS - TMP1*(DX*DX_DYP + DY*DY_DYP)
C
        ANP_XYP(1,1,IP) =   DY_DXP/DS + TMP1*(DX*DX_DXP + DY*DY_DXP)
        ANP_XYP(1,2,IP) =   DY_DYP/DS + TMP1*(DX*DX_DYP + DY*DY_DYP)
C
        TMP2 = -ANP(2,IP)/DSQ
        ANP_XYM(2,1,IP) =   DX_DXM/DS - TMP2*(DX*DX_DXM + DY*DY_DXM)
        ANP_XYM(2,2,IP) =   DX_DYM/DS - TMP2*(DX*DX_DYM + DY*DY_DYM)
C
        ANP_XYO(2,1,IP) = - DX_DXM/DS + TMP2*(DX*DX_DXM + DY*DY_DXM)
     &                    + DX_DXP/DS - TMP2*(DX*DX_DXP + DY*DY_DXP)
        ANP_XYO(2,2,IP) = - DX_DYM/DS + TMP2*(DX*DX_DYM + DY*DY_DYM)
     &                    + DX_DYP/DS - TMP2*(DX*DX_DYP + DY*DY_DYP)
C
        ANP_XYP(2,1,IP) = - DX_DXP/DS + TMP2*(DX*DX_DXP + DY*DY_DXP)
        ANP_XYP(2,2,IP) = - DX_DYP/DS + TMP2*(DX*DX_DYP + DY*DY_DYP)
 10   CONTINUE
C
      RETURN
      END ! ANPSET


      SUBROUTINE ITPSET(IEL)
C--------------------------------------
C     Sets TE panel endpoint indices
C--------------------------------------
      INCLUDE 'DFDC.INC'
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
C---- default case... no TE panel
      IPTE1(IEL) = -999
      IPTE2(IEL) = -999
C
C---- no TE panel on closed axis body or non-solid body
      IF(LXBOD(IEL) .OR. NETYPE(IEL).NE.0) RETURN
C
C---- TE panel must be chosen...
      IF(LTPAN(IEL)) THEN
C----- usual TE panel closing off body
       IPTE1(IEL) = IP1
       IPTE2(IEL) = IP2
       IF    (YP(IP1).EQ.0.0 .AND. YP(IP2).GT.0.0) THEN
        IPTE1(IEL) = 0
        IPTE2(IEL) = IP2
       ELSEIF(YP(IP2).EQ.0.0 .AND. YP(IP1).GT.0.0) THEN
        IPTE1(IEL) = IP1
        IPTE2(IEL) = 0
       ENDIF
      ENDIF
C
      RETURN
      END ! ITPSET
