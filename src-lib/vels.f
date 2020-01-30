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


      SUBROUTINE QCPFOR
C---------------------------------------------
C     Computes velocities, Cp, forces, etc.
C     for current singularity distribitions
C---------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION ANT(2)
C
C---- compute velocities at control points (average between two panel sides)
      IF(.NOT.LQCNT) CALL QCSUM
C
C---- set velocities on both sides of panels at all control points
      CALL QQCALC(NCTOT,ANC, IPCO,IPCP, GAM,SIG,GTH, QC, QCL,QCR )
C
C---- also set QCL,QCR sensitivities to imposed freestream, singularities, etc.
      DO IU=1, NU
        CALL QQCALC(NCTOT,ANC, IPCO,IPCP, 
     &              GAMU(1,IU),SIGU(1,IU),GTHU(1,IU), 
     &              QCU(1,1,IU), QCLU(1,1,IU),QCRU(1,1,IU) )
      ENDDO
C
C---- set velocities and Cp on both sides of panels at all control points
      CALL CPCALC(NCTOT, QCL, QINF,QREF, CPL, CPL_QCL,CPL_QINF)
      CALL CPCALC(NCTOT, QCR, QINF,QREF, CPR, CPR_QCR,CPR_QINF)
C
C---- add in contribution from rotor enthalpy and entropy losses to Cp's
C     CB and duct walls and on upper and lower wakes
      QRSQ = QREF**2
      DO IC = 1, NCTOT
        IP1 = IPCO(IC)
        IP2 = IPCP(IC)
        IR1 = IP2IR(IP1)
        IR2 = IP2IR(IP2)
        IG  = IC2IG(IC)
C---- check for control point in rotor slipstream (in grid)
        IF(IG.GT.0) THEN
          IR = IR1
          IF(IR.EQ.1) THEN
           JG = IR
           IF(IP2.LE.IPLAST(1)) THEN
C---- on CB downstream of rotor
            CIRC = BGAMG(IG,JG)
            VTT  = CIRC*PI2I/YC(IC)
            VTSQ = VTT**2
            DELTACP = 2.0*(DHG(IG,JG)-DSG(IG,JG))/QRSQ - VTSQ/QRSQ
            CPR(IC) = CPR(IC) + DELTACP
           ELSE
C---- on CB wake element
            CIRC = BGAMG(IG,JG)
            VTT  = CIRC*PI2I/YC(IC)
            VTSQ = VTT**2
            DELTACP = 2.0*(DHG(IG,JG)-DSG(IG,JG))/QRSQ - VTSQ/QRSQ
            CPL(IC) = CPL(IC) + DELTACP
           ENDIF
C
          ELSEIF(IR.EQ.NRP) THEN
           JG = IR-1
           IF(IP1.LT.IPLAST(2)) THEN
C---- on DUCT downstream of rotor
            CIRC = BGAMG(IG,JG)
            VTT  = CIRC*PI2I/YC(IC)
            VTSQ = VTT**2
            DELTACP = 2.0*(DHG(IG,JG)-DSG(IG,JG))/QRSQ - VTSQ/QRSQ
            CPR(IC) = CPR(IC) + DELTACP
           ELSE
C---- on DUCT wake element
            CIRC = BGAMG(IG,JG)
            VTT  = CIRC*PI2I/YC(IC)
            VTSQ = VTT**2
            DELTACP = 2.0*(DHG(IG,JG)-DSG(IG,JG))/QRSQ - VTSQ/QRSQ
            CPR(IC) = CPR(IC) + DELTACP
           ENDIF
          ELSE
           JG = IR-1
           CIRC = BGAMG(IG,JG)
           VTT  = CIRC*PI2I/YC(IC)
           VTSQ = VTT**2
           DELTACP = 2.0*(DHG(IG,JG)-DSG(IG,JG))/QRSQ - VTSQ/QRSQ
           CPR(IC) = CPR(IC) + DELTACP
C
           JG = IR
           CIRC = BGAMG(IG,JG)
           VTT  = CIRC*PI2I/YC(IC)
           VTSQ = VTT**2
           DELTACP = 2.0*(DHG(IG,JG)-DSG(IG,JG))/QRSQ - VTSQ/QRSQ
           CPL(IC) = CPL(IC) + DELTACP
          ENDIF
        ENDIF
      END DO
C
C---- also set CPL,CPR sensitivities
      DO IU=1, NU
        DO IC=1, NCTOT
          CPLU(IC,IU) = CPL_QCL(1,IC)*QCLU(1,IC,IU)
     &                + CPL_QCL(2,IC)*QCLU(2,IC,IU)
          CPRU(IC,IU) = CPR_QCR(1,IC)*QCRU(1,IC,IU)
     &                + CPR_QCR(2,IC)*QCRU(2,IC,IU)
        ENDDO
      ENDDO
C
C---- integrate Cp to get forces and moment, and all sensitivities
      COSA = 1.0
      SINA = 0.0
      CX(0) = 0.
      CY(0) = 0.
      CM(0) = 0.
      CD(0) = 0.
      CXVIS(0) = 0.
      DO IU=1, NU
        CXU(0,IU) = 0.
        CYU(0,IU) = 0.
        CMU(0,IU) = 0.
        CDU(0,IU) = 0.
      ENDDO
      DO IEL = 1, NEL
        I = ICFRST(IEL)
        N = ICLAST(IEL) - ICFRST(IEL) + 1
        IF(N.GE.1) THEN
         CALL FCOEFF(N,CPL(I),CPR(I),XC(I),YC(I),ANC(1,I),DSC(I), 
     &               CX(IEL),CY(IEL),CM(IEL))
C
         CD(IEL) = CX(IEL)
C
C---- Totals go in index 0     
         CX(0) = CX(0) + CX(IEL)
c         CY(0) = CY(0) + CY(IEL)
c         CM(0) = CM(0) + CM(IEL)
c         CD(0) = CD(0) + CD(IEL)
C
c        cosa = u / sqrt(u**2 + v**2)         
c        sina = v / sqrt(u**2 + v**2)
c        cosa_u =  v**2 / q**3 = sina**2/q
c        cosa_v =  -u*v / q**3 = -sina cosa/q
c        sina_u =  -v*u / q**3 = -sina cosa/q
c        sina_v =  u**2 / q**3 = cosa**2/q
c
C------- also set sensitivities due to freestream
         IU = IUQINF 
         IF(LUSET(IU)) THEN
           CALL FCOEFF(N,CPLU(I,IU),CPRU(I,IU),
     &                 XC(I),YC(I),ANC(1,I),DSC(I), 
     &                 CXU(IEL,IU),CYU(IEL,IU),CMU(IEL,IU))
           CDU(IEL,IU) = CXU(IEL,IU)
C---- Totals go in index 0     
           CXU(0,IU) = CXU(0,IU) + CXU(IEL,IU)
c           CYU(0,IU) = CYU(0,IU) + CYU(IEL,IU)
c           CMU(0,IU) = CMU(0,IU) + CMU(IEL,IU)
c           CDU(0,IU) = CDU(0,IU) + CDU(IEL,IU)
         ENDIF
C
C---- Add in viscous forces if BL has been calculated
         IF(LVISC) THEN
           CX(IEL) = CX(IEL) + CXVIS(IEL)
           CD(IEL) = CX(IEL)
C---- Add viscous forces to total forces in index 0     
           CX(0)    = CX(0)    + CXVIS(IEL)
           CXVIS(0) = CXVIS(0) + CXVIS(IEL)
         ENDIF
C
C
C---- Elements with TE panels
C     Add forces due to pressures on TE panel
         IF(LTPAN(IEL) .AND. NETYPE(IEL).EQ.0) THEN
          IPT1 = IPTE1(IEL)
          IPT2 = IPTE2(IEL)
          ICT1 = 0
          ICT2 = 0
          IF(IPT1.EQ.IPFRST(IEL)) ICT1 = ICFRST(IEL)
          IF(IPT1.EQ.IPLAST(IEL)) ICT1 = ICLAST(IEL)
          IF(IPT2.EQ.IPFRST(IEL)) ICT2 = ICFRST(IEL)
          IF(IPT2.EQ.IPLAST(IEL)) ICT2 = ICLAST(IEL)
c          write(12,*) 'iel ',iel
c          write(12,*) 'ipt1,ipt2 ',ipt1,ipt2
c          write(12,*) 'ipfrst,iplast ',ipfrst(iel),iplast(iel)
c          write(12,*) 'ict1,ict2 ',ict1,ict2
c          write(12,*) 'icfrst,iclast ',icfrst(iel),iclast(iel)
C
          IF(ICT1.NE.0 .AND. ICT2.NE.0) THEN
            CPTL = 0.5*(CPL(ICT1) + CPL(ICT2))
            CPTR = 0.5*(CPR(ICT1) + CPR(ICT2))
          ELSEIF(ICT1.EQ.0) THEN
            CPTL = CPL(ICT2)
            CPTR = CPR(ICT2)
          ELSEIF(ICT2.EQ.0) THEN
            CPTL = CPL(ICT1)
            CPTR = CPR(ICT1)
          ENDIF
          XCT = 0.5*(XPT(1,IEL) + XPT(2,IEL))
          YCT = 0.5*(YPT(1,IEL) + YPT(2,IEL))
          DXT = XPT(1,IEL) - XPT(2,IEL)
          DYT = YPT(1,IEL) - YPT(2,IEL)
          DSCT = SQRT(DXT*DXT + DYT*DYT)
          ANT(1) =  DYT/DSCT
          ANT(2) = -DXT/DSCT
c          write(12,*) 'cptl,cptr ',cptl,cptr
c          write(12,*) 'xct,yct ',xct,yct
c          write(12,*) 'dsct ',dsct
c          write(12,*) 'ant ',ant(1),ant(2)
C
          CALL FCOEFF(1,CPTL,CPTR,XCT,YCT,ANT,DSCT, 
     &                CXTE,CYTE,CMTE)
c          write(12,*) 'CXTE ',cxte
c          write(12,*) 'CYTE ',cyte
c          write(12,*) 'CMTE ',cmte
C
C---- Add in TE forces
          CX(IEL) = CX(IEL) + CXTE
c          CY(IEL) = CY(IEL) + CYTE
c          CD(IEL) = CX(IEL)
C---- Add TE forces to total forces in index 0     
          CX(0) = CX(0) + CXTE
c          CY(0) = CY(0) + CYTE
c          CM(0) = CM(0) + CMTE
c          CD(0) = CD(0) + CXTE

C------- also set sensitivities due to freestream
          IU = IUQINF 
          IF(LUSET(IU)) THEN
            IF(ICT1.NE.0 .AND. ICT2.NE.0) THEN
              CPTL = 0.5*(CPLU(ICT1,IU) + CPLU(ICT2,IU))
              CPTR = 0.5*(CPRU(ICT1,IU) + CPRU(ICT2,IU))
            ELSEIF(ICT1.EQ.0) THEN
              CPTL = CPLU(ICT2,IU)
              CPTR = CPRU(ICT2,IU)
            ELSEIF(ICT2.EQ.0) THEN
              CPTL = CPLU(ICT1,IU)
              CPTR = CPRU(ICT1,IU)
            ENDIF
            CALL FCOEFF(1,CPTL,CPTR,XCT,YCT,ANT,DSCT, 
     &                  CXTE,CYTE,CMTE)
            CXU(IEL,IU) = CXU(IEL,IU) + CXTE
c            CYU(IEL,IU) = CYU(IEL,IU) + CYTE
c            CDU(IEL,IU) = CXU(IEL,IU)
C---- Totals go in index 0     
            CXU(0,IU) = CXU(0,IU) + CXTE
c            CYU(0,IU) = CYU(0,IU) + CYTE
c            CMU(0,IU) = CMU(0,IU) + CMTE
c            CDU(0,IU) = CDU(0,IU) + CXTE
          ENDIF
         ENDIF
C
        ENDIF
      ENDDO
C
C---- Save duct + centerbody (IEL=1,2) X force (thrust)
      QUE = 0.5*RHO*QREF**2
      TDUCT = -(CX(1)+CX(2))*QUE
C
      RETURN
      END ! QCPFOR



      SUBROUTINE QQCALC(N,AN, IPO,IPP, GAM,SIG,GTH, Q, QL,QR)
C-------------------------------------------------------------
C     Sets velocity on both sides of vortex,source sheet.
C-------------------------------------------------------------
      DIMENSION AN(2,*)
      DIMENSION IPO(*), IPP(*)
      DIMENSION GAM(*), SIG(*), GTH(*)
      DIMENSION Q(2,*), QL(2,*), QR(2,*)
C
      DO I = 1, N
        IO = IPO(I)
        IP = IPP(I)
C
        GAMC = 0.5*(GAM(IO) + GAM(IP))
        SIGC = 0.5*(SIG(IO) + SIG(IP))
C---- Add vorticity contribution from additional vortex wake gamma
        GAMC = GAMC + 0.5*(GTH(IO) + GTH(IP))
C
C------ 1/2 of sheet strength adds,subtracts to left,right speed
        DQT = 0.5*GAMC
        DQN = 0.5*SIGC
C
        QR(1,I) = Q(1,I) + AN(1,I)*DQN + AN(2,I)*DQT
        QR(2,I) = Q(2,I) + AN(2,I)*DQN - AN(1,I)*DQT
C
        QL(1,I) = Q(1,I) - AN(1,I)*DQN - AN(2,I)*DQT
        QL(2,I) = Q(2,I) - AN(2,I)*DQN + AN(1,I)*DQT
      ENDDO
C
      RETURN
      END ! QQCALC



      SUBROUTINE QTCALC(N,AN,IPO,IPP, GAM,GTH, Q, QT)
C-------------------------------------------------------------
C     Sets tangential velocity component
C-------------------------------------------------------------
      DIMENSION AN(2,*)
      DIMENSION IPO(*), IPP(*)
      DIMENSION GAM(*), GTH(*)
      DIMENSION Q(2,*), QT(*)
C
      DIMENSION AT(2)
C
      DO I = 1, N
        IO = IPO(I)
        IP = IPP(I)
C
        GAMC = 0.5*(GAM(IO) + GAM(IP))
C---- Add vorticity from vortex wake gamma
        GAMC = GAMC + 0.5*(GTH(IO) + GTH(IP))
C
C------ set tangential vector along sheet (AN  points to right side)
        AT(1) = -AN(2,I)
        AT(2) =  AN(1,I)
C
C------ 1/2 of sheet strength subtracts to right-side speed (KVSIDE=+1)
        DQT = -0.5*GAMC
C
        QT(I) = Q(1,I)*AT(1) + Q(2,I)*AT(2) + DQT
      ENDDO
C
      RETURN
      END ! QTCALC



      SUBROUTINE QJAC(DX,DY,DU,DV,
     &             QJ,QJ_DX,QJ_DY,QJ_DU,QJ_DV )
C------------------------------------------------------------------
C     Computes velocity jacobian from position,velocity deltas,
C     and zero velocity div and curl conditions.
C------------------------------------------------------------------
      REAL QJ(2,2), QJ_DX(2,2),QJ_DY(2,2),QJ_DU(2,2),QJ_DV(2,2)
C
      DSQ = DX**2 + DY**2
      IF(DSQ.EQ.0.0) THEN
       DSQI = 1.0
      ELSE
       DSQI = 1.0/DSQ
      ENDIF
C
      QJ(1,1) = (DU*DX - DV*DY)*DSQI
      QJ(1,2) = (DU*DY + DV*DX)*DSQI
C
      QJ_DX(1,1) = ( DU - 2.0*QJ(1,1)*DX)*DSQI
      QJ_DY(1,1) = (-DV - 2.0*QJ(1,1)*DY)*DSQI
C
      QJ_DX(1,2) = ( DV - 2.0*QJ(1,2)*DX)*DSQI
      QJ_DY(1,2) = ( DU - 2.0*QJ(1,2)*DY)*DSQI
C
      QJ_DU(1,1) =   DX*DSQI
      QJ_DV(1,1) =  -DY*DSQI
C
      QJ_DU(1,2) =   DY*DSQI
      QJ_DV(1,2) =   DX*DSQI
C
C
      QJ(2,1) = QJ(1,2)
      QJ_DX(2,1) = QJ_DX(1,2)
      QJ_DY(2,1) = QJ_DY(1,2)
      QJ_DU(2,1) = QJ_DU(1,2)
      QJ_DV(2,1) = QJ_DV(1,2)
C
      QJ(2,2) = -QJ(1,1)
      QJ_DX(2,2) = -QJ_DX(1,1)
      QJ_DY(2,2) = -QJ_DY(1,1)
      QJ_DU(2,2) = -QJ_DU(1,1)
      QJ_DV(2,2) = -QJ_DV(1,1)
C
      RETURN
      END ! QJAC


      SUBROUTINE QCURV(U,V,QJ,
     &                 CRV,CRV_U,CRV_V,CRV_QJ)
C-----------------------------------------------------
C     Computes streamline curvature from velocity
C     and velocity Jacobian.
C-----------------------------------------------------
      REAL QJ(2,2), CRV_QJ(2,2)
C
      QSQ = U**2 + V**2
      Q = SQRT(QSQ)
C
      IF(Q.EQ.0.0) THEN
       QI = 1.0
      ELSE
       QI = 1.0/Q
      ENDIF
C
      UN = U*QI
      VN = V*QI
      UN_U =  VN*VN*QI
      UN_V = -UN*VN*QI
      VN_U = -UN*VN*QI
      VN_V =  UN*UN*QI
C
      QI_U = -UN*QI**2
      QI_V = -VN*QI**2
C
C
      CRV    = (  UN*UN*QJ(2,1) -  UN*VN*QJ(1,1)
     &          + UN*VN*QJ(2,2) -  VN*VN*QJ(1,2) ) * QI
C
C---- CRV( QJ UN(U V) VN(U V) QI(U V) )
      CRV_QJ(1,1) = -UN*VN*QI
      CRV_QJ(1,2) = -VN*VN*QI
      CRV_QJ(2,1) =  UN*UN*QI
      CRV_QJ(2,2) =  UN*VN*QI
C
      CRV_UN = ( 2.0*UN*QJ(2,1) -     VN*QJ(1,1) 
     &          +    VN*QJ(2,2)                  ) * QI
      CRV_VN = (                -     UN*QJ(1,1) 
     &          +    UN*QJ(2,2) - 2.0*VN*QJ(1,2) ) * QI
      CRV_QI =    UN*UN*QJ(2,1) -  UN*VN*QJ(1,1)
     &          + UN*VN*QJ(2,2) -  VN*VN*QJ(1,2)
C
C---- CRV( QJ U V )
      CRV_U = CRV_UN*UN_U + CRV_VN*VN_U + CRV_QI*QI_U
      CRV_V = CRV_UN*UN_V + CRV_VN*VN_V + CRV_QI*QI_V
C
      RETURN
      END ! QCURV



      SUBROUTINE CV2SET(ICM,ICO,
     &                  QC,XC,YC, ANC, GAMC,SIGC,
     &                  CV, CV_QC,CV_RC, CV_ANC,
     &                  CV_GAMC,CV_SIGC )
      REAL QC(2,*), XC(*),YC(*), ANC(2,*)
      REAL GAMC(*), SIGC(*)
      REAL CV_QC(2,-1:0), CV_RC(2,-1:0), CV_ANC(2,-1:0)
      REAL CV_GAMC(-1:0), CV_SIGC(-1:0)
C
      REAL QSM(2),QSO(2), QSA(2)
      REAL DR(2), DQ(2), QJ(2,2), QJ_DR(2,2,2), QJ_DQ(2,2,2)
      REAL CV_Q(2), CV_QJ(2,2)
      REAL CV_QSM(2),CV_QSO(2)
C
      GAMCM = GAMC(ICM)
      SIGCM = SIGC(ICM)
C
      GAMCO = GAMC(ICO)
      SIGCO = SIGC(ICO)
C
C---- 1/2 of sheet strength adds,subtracts to left,right speed
      QSGN = +0.5
C
      QSM(1) = QC(1,ICM) + QSGN*(ANC(1,ICM)*SIGCM + ANC(2,ICM)*GAMCM)
      QSM(2) = QC(2,ICM) + QSGN*(ANC(2,ICM)*SIGCM - ANC(1,ICM)*GAMCM)
C
      QSO(1) = QC(1,ICO) + QSGN*(ANC(1,ICO)*SIGCO + ANC(2,ICO)*GAMCO)
      QSO(2) = QC(2,ICO) + QSGN*(ANC(2,ICO)*SIGCO - ANC(1,ICO)*GAMCO)
C
      DQ(1) = QSO(1) - QSM(1)
      DQ(2) = QSO(2) - QSM(2)
C
      QSA(1) = (QSO(1) + QSM(1))*0.5
      QSA(2) = (QSO(2) + QSM(2))*0.5
C
      DR(1) = XC(ICO) - XC(ICM)
      DR(2) = YC(ICO) - YC(ICM)
C
C
      CALL QJAC( DR(1),DR(2),
     &           DQ(1),DQ(2),
     &           QJ,QJ_DR(1,1,1),QJ_DR(1,1,2),
     &              QJ_DQ(1,1,1),QJ_DQ(1,1,2) )
C
      CALL QCURV( QSA(1),QSA(2),QJ,
     &            CV, CV_Q(1) ,CV_Q(2) ,CV_QJ )
C
      DO K = 1, 2
        CV_DQ = CV_QJ(1,1)*QJ_DQ(1,1,K)
     &        + CV_QJ(1,2)*QJ_DQ(1,2,K)
     &        + CV_QJ(2,1)*QJ_DQ(2,1,K)
     &        + CV_QJ(2,2)*QJ_DQ(2,2,K)
        CV_DR = CV_QJ(1,1)*QJ_DR(1,1,K)
     &        + CV_QJ(1,2)*QJ_DR(1,2,K)
     &        + CV_QJ(2,1)*QJ_DR(2,1,K)
     &        + CV_QJ(2,2)*QJ_DR(2,2,K)
C
        CV_QSM(K) = -CV_DQ + 0.5*CV_Q(K)
        CV_QSO(K) =  CV_DQ + 0.5*CV_Q(K)
C
        CV_QC(K,-1) = CV_QSM(K)
        CV_QC(K, 0) = CV_QSO(K)
C
        CV_RC(K,-1) = -CV_DR
        CV_RC(K, 0) =  CV_DR
C
        CV_ANC(K,-1) = CV_QSM(K)*QSGN*SIGCM
        CV_ANC(K, 0) = CV_QSO(K)*QSGN*SIGCO
      ENDDO
C
      CV_ANC(1,-1) = CV_ANC(1,-1) - CV_QSM(2)*QSGN*GAMCM
      CV_ANC(1, 0) = CV_ANC(1, 0) - CV_QSO(2)*QSGN*GAMCO
      CV_ANC(2,-1) = CV_ANC(2,-1) + CV_QSM(1)*QSGN*GAMCM
      CV_ANC(2, 0) = CV_ANC(2, 0) + CV_QSO(1)*QSGN*GAMCO
C
      CV_SIGC(-1) = QSGN*(CV_QSM(1)*ANC(1,ICM) + CV_QSM(2)*ANC(2,ICM))
      CV_SIGC( 0) = QSGN*(CV_QSO(1)*ANC(1,ICO) + CV_QSO(2)*ANC(2,ICO))
C
      CV_GAMC(-1) = QSGN*(CV_QSM(1)*ANC(2,ICM) - CV_QSM(2)*ANC(1,ICM))
      CV_GAMC( 0) = QSGN*(CV_QSO(1)*ANC(2,ICO) - CV_QSO(2)*ANC(1,ICO))
C
      RETURN
      END ! CV2SET


      SUBROUTINE CV3SET(ICM,ICO,ICP,
     &                  QC,XC,YC, ANC, ANSGN, GAMC,SIGC,
     &                  CV, CV_QC,CV_RC, CV_ANC,
     &                  CV_GAMC,CV_SIGC )
      REAL QC(2,*), XC(*),YC(*), ANC(2,*)
      REAL GAMC(*), SIGC(*)
      REAL CV_QC(2,-1:1), CV_RC(2,-1:1), CV_ANC(2,-1:1)
      REAL CV_GAMC(-1:1), CV_SIGC(-1:1)
C
      REAL QSM(2),QSO(2),QSP(2)
      REAL DR(2), DQ(2), QJ(2,2), QJ_DR(2,2,2), QJ_DQ(2,2,2)
      REAL CV_Q(2), CV_QJ(2,2)
      REAL CV_QSM(2),CV_QSO(2),CV_QSP(2)
C
      GAMCM = GAMC(ICM)
      SIGCM = SIGC(ICM)
C
      GAMCO = GAMC(ICO)
      SIGCO = SIGC(ICO)
C
      GAMCP = GAMC(ICP)
      SIGCP = SIGC(ICP)
C
C---- 1/2 of sheet strength adds,subtracts to left,right speed
      QSGN = 0.5*ANSGN
C
      QSM(1) = QC(1,ICM) + QSGN*(ANC(1,ICM)*SIGCM + ANC(2,ICM)*GAMCM)
      QSM(2) = QC(2,ICM) + QSGN*(ANC(2,ICM)*SIGCM - ANC(1,ICM)*GAMCM)
C
      QSO(1) = QC(1,ICO) + QSGN*(ANC(1,ICO)*SIGCO + ANC(2,ICO)*GAMCO)
      QSO(2) = QC(2,ICO) + QSGN*(ANC(2,ICO)*SIGCO - ANC(1,ICO)*GAMCO)
C
      QSP(1) = QC(1,ICP) + QSGN*(ANC(1,ICP)*SIGCP + ANC(2,ICP)*GAMCP)
      QSP(2) = QC(2,ICP) + QSGN*(ANC(2,ICP)*SIGCP - ANC(1,ICP)*GAMCP)
C
      DQ(1) = QSP(1) - QSM(1)
      DQ(2) = QSP(2) - QSM(2)
C
      DR(1) = XC(ICP) - XC(ICM)
      DR(2) = YC(ICP) - YC(ICM)
C
      QXN = ANSGN*(QSO(1)*ANC(2,ICO) - QSO(2)*ANC(1,ICO))
      CSGN = SIGN( 1.0 , QXN )
C
      CALL QJAC( DR(1),DR(2),
     &           DQ(1),DQ(2),
     &           QJ,QJ_DR(1,1,1),QJ_DR(1,1,2),
     &              QJ_DQ(1,1,1),QJ_DQ(1,1,2) )
C
      CALL QCURV( QSO(1),QSO(2),QJ,
     &            CV, CV_Q(1) ,CV_Q(2) ,CV_QJ )
C
      CV           = CSGN*CV
      DO K = 1, 2
        CV_Q(K)    = CSGN*CV_Q(K)
        CV_QJ(K,1) = CSGN*CV_QJ(K,1)
        CV_QJ(K,2) = CSGN*CV_QJ(K,2)
      ENDDO
C
      DO K = 1, 2
        CV_DQ = CV_QJ(1,1)*QJ_DQ(1,1,K)
     &        + CV_QJ(1,2)*QJ_DQ(1,2,K)
     &        + CV_QJ(2,1)*QJ_DQ(2,1,K)
     &        + CV_QJ(2,2)*QJ_DQ(2,2,K)
        CV_DR = CV_QJ(1,1)*QJ_DR(1,1,K)
     &        + CV_QJ(1,2)*QJ_DR(1,2,K)
     &        + CV_QJ(2,1)*QJ_DR(2,1,K)
     &        + CV_QJ(2,2)*QJ_DR(2,2,K)
C
        CV_QSM(K) = -CV_DQ
        CV_QSO(K) =  CV_Q(K)
        CV_QSP(K) =  CV_DQ
C
        CV_QC(K,-1) = CV_QSM(K)
        CV_QC(K, 0) = CV_QSO(K)
        CV_QC(K,+1) = CV_QSP(K)
C
        CV_RC(K,-1) = -CV_DR
        CV_RC(K, 0) = 0.
        CV_RC(K,+1) =  CV_DR
C
        CV_ANC(K,-1) = CV_QSM(K)*QSGN*SIGCM
        CV_ANC(K, 0) = CV_QSO(K)*QSGN*SIGCO
        CV_ANC(K,+1) = CV_QSP(K)*QSGN*SIGCP
      ENDDO
C
      CV_ANC(1,-1) = CV_ANC(1,-1) - CV_QSM(2)*QSGN*GAMCM
      CV_ANC(1, 0) = CV_ANC(1, 0) - CV_QSO(2)*QSGN*GAMCO
      CV_ANC(1,+1) = CV_ANC(1,+1) - CV_QSP(2)*QSGN*GAMCP
      CV_ANC(2,-1) = CV_ANC(2,-1) + CV_QSM(1)*QSGN*GAMCM
      CV_ANC(2, 0) = CV_ANC(2, 0) + CV_QSO(1)*QSGN*GAMCO
      CV_ANC(2,+1) = CV_ANC(2,+1) + CV_QSP(1)*QSGN*GAMCP
C
      CV_SIGC(-1) = QSGN*(CV_QSM(1)*ANC(1,ICM) + CV_QSM(2)*ANC(2,ICM))
      CV_SIGC( 0) = QSGN*(CV_QSO(1)*ANC(1,ICO) + CV_QSO(2)*ANC(2,ICO))
      CV_SIGC(+1) = QSGN*(CV_QSP(1)*ANC(1,ICP) + CV_QSP(2)*ANC(2,ICP))
C
      CV_GAMC(-1) = QSGN*(CV_QSM(1)*ANC(2,ICM) - CV_QSM(2)*ANC(1,ICM))
      CV_GAMC( 0) = QSGN*(CV_QSO(1)*ANC(2,ICO) - CV_QSO(2)*ANC(1,ICO))
      CV_GAMC(+1) = QSGN*(CV_QSP(1)*ANC(2,ICP) - CV_QSP(2)*ANC(1,ICP))
C
      RETURN
      END ! CV3SET


      SUBROUTINE GETUV(XX,YY,US,VS)
C-----------------------------------------------
C     Get EIF velocity US,VS at point XX,YY
C-----------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- local arrays for calling QFCALC
      DIMENSION QF(2),
     &          QF_GAM(2,IPX),
     &          QF_SIG(2,IPX),
     &          QF_GTH(2,IPX),
     &          CPF_QF(2)
C
C------ we will evaluate only this point
      IP1 = 1
      IP2 = 1
      NF  = 1 
C------ assume the point is not on a panel
      IPFO = 0
      IPFP = 0
      IFTYPE = 0
C------ evaluate velocity components
      CALL QFCALC(IP1,IP2, XX,YY, IPFO,IPFP, IFTYPE,
     &            QF,QF_GAM,QF_SIG,QF_GTH)
C------ add freestream
      US = QF(1) + QINF
      VS = QF(2)
C
      RETURN
      END
           



