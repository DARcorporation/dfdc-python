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

      SUBROUTINE QAIC(LXYJAC)
      INCLUDE 'DFDC.INC'
      LOGICAL LXYJAC
C
      IF(LDBG) THEN
       IF(LXYJAC) THEN
         WRITE(*,*) 'Generating dq/d(gam,sig,x,y) Jacobian matrices...'
       ELSE
         WRITE(*,*) 'Generating dq/d(gam,sig) Jacobian matrices...'
       ENDIF
      ENDIF
C
      CALL QAIC1(LXYJAC,1,NCTOT,1,NEL)
      DO KEL = 1, NEL
        ICE = NCX + KEL
        IF(LBODY(KEL)) CALL QAIC1(LXYJAC,ICE,ICE,1,NEL)
      ENDDO
C
      LQCNT = .TRUE.
      LQAIC = .TRUE.
      LQGIC = LXYJAC
C
      RETURN
      END ! QAIC



      SUBROUTINE QAIC1(LXYJAC, IC1,IC2, IEL1,IEL2)
C---------------------------------------------------------------
C     Computes velocities and Jacobians at control points IC1...IC2
C     w.r.t aero quantities at nodes on elements IEL1..IEL2
C     
C     If LXYJAC=T, also computes Jacobians w.r.t. panel geometry
C
C     Input:  XC,YC(i)  control points
C             XP,YP(j)  panel geometry
C             GAM(j)    panel vorticity  (or line circ. , or point doublet)
C             SIG(j)    panel source     (or line source, or point source )
C             GTH(j)    vortex wake panel vorticity
C             QINF      freestream speed
C             ALFA      freestream angle
C
C     Output: QC(.i)         velocity at XC(i),YC(i)
C             QC_GAM(.ji)    dQC(.i)/dGAM(j)
C             QC_SIG(.ji)
C             QC_GTH(.ji)    dQC(.i)/dGTH(j)
C
C     Additonal output if LXYJAC=T:
C             QC_XP(.ji) 
C             QC_YP(.ji)
C             QC_XC(.i)
C             QC_YC(.i)
C
C     Note:   TE panel information and strengths are set by
C             this routine
C---------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LXYJAC
C
      DIMENSION QC_GAMTE(2,2,ICX), QC_SIGTE(2,2,ICX)
      DIMENSION QC_XPT(2,2,ICX),   QC_YPT (2,2,ICX)
C
      DIMENSION RCORE(IPX)
C
C---- maximum TE gap (relative to perimeter) to give closed airfoil a TE panel
      DATA TFRAC / 0.00001 /
C
C---- clear control-point velocities for accumulation in PANAIC and PNTAIC
      DO IC = IC1, IC2
        QC(1,IC) = 0.
        QC(2,IC) = 0.
        DO IEL = IEL1, IEL2
          DO IP = IPFRST(IEL), IPLAST(IEL)
            QC_GAM(1,IP,IC) = 0.
            QC_GAM(2,IP,IC) = 0.
            QC_SIG(1,IP,IC) = 0.
            QC_SIG(2,IP,IC) = 0.
            QC_XP (1,IP,IC) = 0.
            QC_XP (2,IP,IC) = 0.
            QC_YP (1,IP,IC) = 0.
            QC_YP (2,IP,IC) = 0.
            QC_GTH(1,IP,IC) = 0.
            QC_GTH(2,IP,IC) = 0.
          ENDDO
CHHY- Uncomment to save TE panel influences separately
cc          DO KP = 1, 2
cc            QC_GAMT(1,KP,IEL,IC) = 0.
cc            QC_GAMT(2,KP,IEL,IC) = 0.
cc            QC_SIGT(1,KP,IEL,IC) = 0.
cc            QC_SIGT(2,KP,IEL,IC) = 0.
cc          ENDDO
        ENDDO
      ENDDO
C
      KC = IC1
      NUMIC = IC2 - IC1 + 1
C
C---- Set AIC matrix for all singularity types
      DO 20 IEL = IEL1, IEL2
        IF(    NETYPE(IEL).EQ.0
     &    .OR. NETYPE(IEL).EQ.1 
     &    .OR. NETYPE(IEL).EQ.5 
     &    .OR. NETYPE(IEL).EQ.6 
     &    .OR. NETYPE(IEL).EQ.7 ) THEN
C-------- wall,wake,source-only line and vortex wake sheets
          IF(LXYJAC) THEN
C--------- linearize QC w.r.t GAM,SIG,XP,YP,XC,YC
           CALL PANAICD(IPFRST(IEL),IPLAST(IEL),XP,YP,
     &           NUMIC,XC(KC),YC(KC), IPCO(KC),IPCP(KC), ICTYPE(KC),
     &               GAM, SIG,
     &               IPX, QC(1,KC), QC_GAM(1,1,KC), QC_SIG(1,1,KC),
     &                              QC_XP (1,1,KC), QC_YP (1,1,KC),
     &                              QC_XC (1,  KC), QC_YC (1,  KC)  )
           DO IC = IC1, IC2
             IPO = IPCO(IC)
             IPP = IPCP(IC)
             IF(IPO.NE.0 .AND. IPP.NE.0) THEN
              QC_XP(1,IPO,IC) = QC_XP(1,IPO,IC) + QC_XC(1,IC)*0.5
              QC_XP(2,IPO,IC) = QC_XP(2,IPO,IC) + QC_XC(2,IC)*0.5
              QC_XP(1,IPP,IC) = QC_XP(1,IPP,IC) + QC_XC(1,IC)*0.5
              QC_XP(2,IPP,IC) = QC_XP(2,IPP,IC) + QC_XC(2,IC)*0.5
              QC_YP(1,IPO,IC) = QC_YP(1,IPO,IC) + QC_YC(1,IC)*0.5
              QC_YP(2,IPO,IC) = QC_YP(2,IPO,IC) + QC_YC(2,IC)*0.5
              QC_YP(1,IPP,IC) = QC_YP(1,IPP,IC) + QC_YC(1,IC)*0.5
              QC_YP(2,IPP,IC) = QC_YP(2,IPP,IC) + QC_YC(2,IC)*0.5
             ENDIF
           ENDDO
          ELSE
C--------- linearize QC w.r.t GAM,SIG only
           CALL PANAIC(IPFRST(IEL),IPLAST(IEL), XP,YP,
     &             NUMIC,XC(KC),YC(KC), IPCO(KC),IPCP(KC), ICTYPE(KC),
     &                 GAM, SIG,
     &                 IPX, QC(1,KC), QC_GAM(1,1,KC), QC_SIG(1,1,KC) )
          ENDIF
C
        ELSEIF(NETYPE(IEL).EQ.2) THEN
C------- line on axis
         DO IP = IPFRST(IEL), IPLAST(IEL)
           RCORE(IP) = 0.01
         ENDDO
         CALL AXLAIC(IPFRST(IEL),IPLAST(IEL), XP,YP,
     &               NUMIC,XC(KC),YC(KC), IPCO(KC),IPCP(KC), ICTYPE(KC),
     &               GAM, SIG, RCORE,
     &               IPX, QC(1,KC), QC_GAM(1,1,KC), QC_SIG(1,1,KC) )
C
        ELSEIF(NETYPE(IEL).EQ.3) THEN
C------- ring
         CALL LINAIC(IPFRST(IEL),IPLAST(IEL), XP,YP,
     &               NUMIC,XC(KC),YC(KC),
     &               GAM,SIG,
     &               IPX, QC(1,KC), QC_GAM(1,1,KC), QC_SIG(1,1,KC) )
C
        ELSEIF(NETYPE(IEL).EQ.4) THEN
C------- point
         CALL PNTAIC(IPFRST(IEL),IPLAST(IEL), XP,YP,
     &               NUMIC,XC(KC),YC(KC),
     &               GAM,SIG,
     &               IPX, QC(1,KC), QC_GAM(1,1,KC), QC_SIG(1,1,KC) )
        ENDIF
C
 20   CONTINUE
C
C---- Accumulate velocity from wake vorticity, GTH
C---- Save velocity AIC matrix (w/o TE panel influences) for wake vorticity 
      DO IEL = IEL1, IEL2
        IF(    NETYPE(IEL).EQ.0
     &    .OR. NETYPE(IEL).EQ.7 ) THEN
C
          IP1 = IPFRST(IEL)
          IP2 = IPLAST(IEL)
          DO IC = IC1, IC2
            DO IP = IP1, IP2
              QC_GTH(1,IP,IC) = QC_GAM(1,IP,IC) 
              QC_GTH(2,IP,IC) = QC_GAM(2,IP,IC) 
              QC(1,IC) = QC(1,IC) + QC_GTH(1,IP,IC)*GTH(IP)
              QC(2,IC) = QC(2,IC) + QC_GTH(2,IP,IC)*GTH(IP)
            END DO
          END DO
C
        ENDIF
      END DO
C
C---- Add contributions of TE panels...
      DO 30 IEL = IEL1, IEL2
        IF(.NOT. (LTPAN(IEL) .AND. 
     &           (NETYPE(IEL).EQ.0 .OR. NETYPE(IEL).EQ.7)) ) THEN
C------- no TE panel... zero out panel strengths
         DO K = 1, 2
           XPT(K,IEL) = 0.
           YPT(K,IEL) = 0.
           GAMT(K,IEL) = 0.
           SIGT(K,IEL) = 0.
           DO J = 1, 2
             GAMT_GAM(K,J,IEL) = 0.
             GAMT_SIG(K,J,IEL) = 0.
             SIGT_GAM(K,J,IEL) = 0.
             SIGT_SIG(K,J,IEL) = 0.
           ENDDO
         ENDDO
         GO TO 30
        ENDIF
C---- clear accumulating control-point velocities for each TE panel
        DO IC = IC1, IC2
          DO KP = 1, 2
            QC_GAMTE(1,KP,IC) = 0.
            QC_GAMTE(2,KP,IC) = 0.
            QC_SIGTE(1,KP,IC) = 0.
            QC_SIGTE(2,KP,IC) = 0.
            QC_XPT (1,KP,IC) = 0.
            QC_XPT (2,KP,IC) = 0.
            QC_YPT (1,KP,IC) = 0.
            QC_YPT (2,KP,IC) = 0.
          ENDDO
          QC_XC(1,IC) = 0.
          QC_XC(2,IC) = 0.
          QC_YC(1,IC) = 0.
          QC_YC(2,IC) = 0.
        ENDDO
C
C------ nodes on element spanned by TE panel
        IP1 = IPTE1(IEL)
        IP2 = IPTE2(IEL)
C
C------ set panel coordinates, strengths, and all Jacobians
        CALL TPANSET(IP1,IP2,GAM,SIG,XP,YP,
     &               XPT(1,IEL),XPT_XP(1,1,IEL),
     &               YPT(1,IEL),YPT_YP(1,1,IEL),
     &   GAMT(1,IEL),GAMT_GAM(1,1,IEL),GAMT_SIG(1,1,IEL),
     &               GAMT_XP (1,1,IEL),GAMT_YP (1,1,IEL),
     &               GAMT_DX (1,  IEL),GAMT_DY (1,  IEL),
     &   SIGT(1,IEL),SIGT_GAM(1,1,IEL),SIGT_SIG(1,1,IEL),
     &               SIGT_XP (1,1,IEL),SIGT_YP (1,1,IEL),
     &               SIGT_DX (1,  IEL),SIGT_DY (1,  IEL) )
C
C------ skip influence computation if TE panel has negligible length
        DXPT = XPT(2,IEL) - XPT(1,IEL)
        DYPT = YPT(2,IEL) - YPT(1,IEL)
        DSPT = SQRT(DXPT**2 + DYPT**2)
C
        SPLEN = ABS( SP(IPLAST(IEL)) - SP(IPFRST(IEL)) )
        IF(DSPT.LT.TFRAC*SPLEN) THEN
         WRITE(*,*) 'Negligible TE panel ignored for element', IEL
         GO TO 30
        ENDIF
C
C
        IF(LXYJAC) THEN
C-------- linearize QC w.r.t GAM,SIG,XP,YP,XC,YC
          CALL PANAICD(1,2,XPT(1,IEL),YPT(1,IEL),
     &               NUMIC,XC(KC),YC(KC), IZERO,IZERO, ICTYPE(KC),
     &               GAMT(1,IEL), SIGT(1,IEL),
     &               2, QC(1,KC), QC_GAMTE(1,1,KC), QC_SIGTE(1,1,KC),
     &                            QC_XPT (1,1,KC),  QC_YPT (1,1,KC),
     &                            QC_XC  (1,  KC),  QC_YC  (1,  KC) )
          DO IC = IC1, IC2
            IPO = IPCO(IC)
            IPP = IPCP(IC)
            IF(IPO.NE.0 .AND. IPP.NE.0) THEN
             QC_XP(1,IPO,IC) = QC_XP(1,IPO,IC) + QC_XC(1,IC)*0.5
             QC_XP(2,IPO,IC) = QC_XP(2,IPO,IC) + QC_XC(2,IC)*0.5
             QC_XP(1,IPP,IC) = QC_XP(1,IPP,IC) + QC_XC(1,IC)*0.5
             QC_XP(2,IPP,IC) = QC_XP(2,IPP,IC) + QC_XC(2,IC)*0.5
             QC_YP(1,IPO,IC) = QC_YP(1,IPO,IC) + QC_YC(1,IC)*0.5
             QC_YP(2,IPO,IC) = QC_YP(2,IPO,IC) + QC_YC(2,IC)*0.5
             QC_YP(1,IPP,IC) = QC_YP(1,IPP,IC) + QC_YC(1,IC)*0.5
             QC_YP(2,IPP,IC) = QC_YP(2,IPP,IC) + QC_YC(2,IC)*0.5
            ENDIF
          ENDDO
C
        ELSE
C-------- linearize QC w.r.t GAM,SIG only
          CALL PANAIC(1,2,XPT(1,IEL),YPT(1,IEL),
     &              NUMIC,XC(KC),YC(KC), IZERO,IZERO, ICTYPE(KC),
     &              GAMT(1,IEL), SIGT(1,IEL),
     &              2, QC(1,KC), QC_GAMTE(1,1,KC), QC_SIGTE(1,1,KC) )
        ENDIF
C
CHHY- Uncomment to save TE panel influences separately
C---- Save TE panel velocity influences
cc        DO K = 1, 2
cc          DO IC = IC1, IC2
cc            QC_GAMT(K,1,IEL,IC) = QC_GAMTE(K,1,IC)
cc            QC_GAMT(K,2,IEL,IC) = QC_GAMTE(K,2,IC)
cc            QC_SIGT(K,1,IEL,IC) = QC_SIGTE(K,1,IC)
cc            QC_SIGT(K,2,IEL,IC) = QC_SIGTE(K,2,IC)
cc          END DO
cc        END DO
C
        DO IC = IC1, IC2
          DO KP = 1, 2
            IF(KP.EQ.1) THEN
             IP = IP1
             JPO = IP1+1
             JPM = IP1
            ELSE
             IP = IP2
             JPO = IP2
             JPM = IP2-1
            ENDIF
            IF(IP.GT.0) THEN
             DO K = 1, 2
               QC_GAM(K,IP,IC) = QC_GAM(K,IP,IC)
     &                         + QC_GAMTE(K,1,IC)*GAMT_GAM(1,KP,IEL)
     &                         + QC_GAMTE(K,2,IC)*GAMT_GAM(2,KP,IEL)
     &                         + QC_SIGTE(K,1,IC)*SIGT_GAM(1,KP,IEL)
     &                         + QC_SIGTE(K,2,IC)*SIGT_GAM(2,KP,IEL)
               QC_SIG(K,IP,IC) = QC_SIG(K,IP,IC)
     &                         + QC_GAMTE(K,1,IC)*GAMT_SIG(1,KP,IEL)
     &                         + QC_GAMTE(K,2,IC)*GAMT_SIG(2,KP,IEL)
     &                         + QC_SIGTE(K,1,IC)*SIGT_SIG(1,KP,IEL)
     &                         + QC_SIGTE(K,2,IC)*SIGT_SIG(2,KP,IEL)
C
               IF(LXYJAC) THEN
                QC_XP(K,IP,IC) = QC_XP(K,IP,IC)
     &                         + QC_XPT (K,1,IC)* XPT_XP(1,KP,IEL)
     &                         + QC_XPT (K,2,IC)* XPT_XP(2,KP,IEL)
     &                         + QC_GAMTE(K,1,IC)*GAMT_XP(1,KP,IEL)
     &                         + QC_GAMTE(K,2,IC)*GAMT_XP(2,KP,IEL)
     &                         + QC_SIGTE(K,1,IC)*SIGT_XP(1,KP,IEL)
     &                         + QC_SIGTE(K,2,IC)*SIGT_XP(2,KP,IEL)
                QC_YP(K,IP,IC) = QC_YP(K,IP,IC)
     &                         + QC_YPT (K,1,IC)* YPT_YP(1,KP,IEL)
     &                         + QC_YPT (K,2,IC)* YPT_YP(2,KP,IEL)
     &                         + QC_GAMTE(K,1,IC)*GAMT_YP(1,KP,IEL)
     &                         + QC_GAMTE(K,2,IC)*GAMT_YP(2,KP,IEL)
     &                         + QC_SIGTE(K,1,IC)*SIGT_YP(1,KP,IEL)
     &                         + QC_SIGTE(K,2,IC)*SIGT_YP(2,KP,IEL)
                QC_XP(K,JPO,IC) = QC_XP(K,JPO,IC)
     &                         + QC_GAMTE(K,KP,IC)*GAMT_DX(KP,IEL)
     &                         + QC_SIGTE(K,KP,IC)*SIGT_DX(KP,IEL)
                QC_YP(K,JPO,IC) = QC_YP(K,JPO,IC)
     &                         + QC_GAMTE(K,KP,IC)*GAMT_DY(KP,IEL)
     &                         + QC_SIGTE(K,KP,IC)*SIGT_DY(KP,IEL)
                QC_XP(K,JPM,IC) = QC_XP(K,JPM,IC)
     &                         - QC_GAMTE(K,KP,IC)*GAMT_DX(KP,IEL)
     &                         - QC_SIGTE(K,KP,IC)*SIGT_DX(KP,IEL)
                QC_YP(K,JPM,IC) = QC_YP(K,JPM,IC)
     &                         - QC_GAMTE(K,KP,IC)*GAMT_DY(KP,IEL)
     &                         - QC_SIGTE(K,KP,IC)*SIGT_DY(KP,IEL)
               ENDIF
             ENDDO ! K loop
            ENDIF
          END DO ! KP loop
        END DO ! IC loop
C
C---- Special treatment for vortex wake last GTH point to include TE panel
        IF(LTPAN(IEL) .AND. NETYPE(IEL).EQ.7) THEN
          IP = IPLAST(IEL)
          DO IC = IC1, IC2
C---- Remove velocity w/o TE panel
            QC(1,IC) = QC(1,IC) - QC_GTH(1,IP,IC)*GTH(IP)
            QC(2,IC) = QC(2,IC) - QC_GTH(2,IP,IC)*GTH(IP)
C---- Save influence with TE panel
            QC_GTH(1,IP,IC) = QC_GAM(1,IP,IC) 
            QC_GTH(2,IP,IC) = QC_GAM(2,IP,IC) 
C---- Add back in velocity with TE panel
            QC(1,IC) = QC(1,IC) + QC_GTH(1,IP,IC)*GTH(IP)
            QC(2,IC) = QC(2,IC) + QC_GTH(2,IP,IC)*GTH(IP)
          END DO
        ENDIF
C
 30   CONTINUE
C
C---- add on freestream
      DO IC = IC1, IC2
        QC(1,IC) = QC(1,IC) + QINF
      ENDDO
C
      RETURN
      END ! QAIC1


 
      SUBROUTINE TPANSET(IPTE1,IPTE2,GAM,SIG,XP,YP,
     &                   XPT,XPT_XP,
     &                   YPT,YPT_YP,
     &                   GAMT,GAMT_GAM,GAMT_SIG,
     &                        GAMT_XP ,GAMT_YP ,
     &                        GAMT_DX ,GAMT_DY ,
     &                   SIGT,SIGT_GAM,SIGT_SIG,
     &                        SIGT_XP ,SIGT_YP ,
     &                        SIGT_DX ,SIGT_DY  )
C----------------------------------------------------------
C     Sets geometry, strengths, and associated Jacobians
C     of TE panel on element IEL
C----------------------------------------------------------
      DIMENSION GAM(*), SIG(*), XP(*), YP(*)
C
      DIMENSION XPT(2), XPT_XP(2,2),
     &          YPT(2), YPT_YP(2,2)
      DIMENSION GAMT(2), GAMT_GAM(2,2), GAMT_SIG(2,2),
     &          SIGT(2), SIGT_GAM(2,2), SIGT_SIG(2,2)
      DIMENSION GAMT_XP(2,2), GAMT_YP(2,2),GAMT_DX(2), GAMT_DY(2),
     &          SIGT_XP(2,2), SIGT_YP(2,2),SIGT_DX(2), SIGT_DY(2)
C
C---- local setup arrays
      DIMENSION ANX(2), ANX_DX(2), ANX_DY(2),
     &          ANY(2), ANY_DX(2), ANY_DY(2)
      DIMENSION ANXT_XP(2), ANXT_YP(2),
     &          ANYT_XP(2), ANYT_YP(2)
C
      DIMENSION DOT(2),DOT_XP(2,2),DOT_YP(2,2),DOT_DX(2),DOT_DY(2),
     &          CRS(2),CRS_XP(2,2),CRS_YP(2,2),CRS_DX(2),CRS_DY(2)
C
      DO K = 1, 2
        DO J = 1, 2
          XPT_XP(K,J) = 0.
          YPT_YP(K,J) = 0.
          GAMT_GAM(K,J) = 0.
          GAMT_SIG(K,J) = 0.
          SIGT_GAM(K,J) = 0.
          SIGT_SIG(K,J) = 0.
          GAMT_XP(K,J) = 0.
          GAMT_YP(K,J) = 0.
          SIGT_XP(K,J) = 0.
          SIGT_YP(K,J) = 0.
        ENDDO
        ANX_DX(K) = 0.
        ANY_DX(K) = 0.
        ANX_DY(K) = 0.
        ANY_DY(K) = 0.
      ENDDO
C
C---- geometry nodes spanned by TE panel
C-    (if node index = 0, then that TE panel node lies on axis)
      IP1 = IPTE1
      IP2 = IPTE2
C
C---- set node x,y of TE panel, and dx,dy of adjacent panels
C
      DO K = 1, 2
        IF(K.EQ.1) THEN
         IP    = IP1
         IPOPP = IP2
         KOPP = 2
         JPO = IP+1
         JPM = IP
         ANYAXI =  1.0
        ELSE
         IP    = IP2
         IPOPP = IP1
         KOPP = 1
         JPO = IP
         JPM = IP-1
         ANYAXI = -1.0
        ENDIF
C        
        IF(IP.GT.0) THEN
C------- usual TE panel
         XPT(K) = XP(IP)
         YPT(K) = YP(IP)
         XPT_XP(K,K) = 1.0
         YPT_YP(K,K) = 1.0
         DX = XP(JPO) - XP(JPM)
         DY = YP(JPO) - YP(JPM)
         DS = SQRT(DX**2 + DY**2)
         ANX(K) =  DY/DS
         ANY(K) = -DX/DS
         ANX_DX(K) = -DX*DY/DS**3
         ANX_DY(K) =  DX*DX/DS**3
         ANY_DX(K) = -DY*DY/DS**3
         ANY_DY(K) =  DY*DX/DS**3
        ELSE
C------- node IP is on axis, at same x location as node IPOPP
         XPT(K) = XP(IPOPP)
         YPT(K) = 0.
         XPT_XP(K,KOPP) = 1.0
         ANX(K) = 0.
         ANY(K) = ANYAXI
        ENDIF
      ENDDO
C
C---- now set TE panel normal vector ANXT,ANYT
      DXT = XPT(1) - XPT(2)
      DYT = YPT(1) - YPT(2)
      DSTSQ = DXT**2 + DYT**2
      IF(DSTSQ.EQ.0.0) THEN
       DSTI = 1.0
      ELSE
       DSTI = 1.0/SQRT(DSTSQ)
      ENDIF
      ANXT =  DYT*DSTI
      ANYT = -DXT*DSTI
      ANXT_DXT = -DXT*DYT*DSTI**3
      ANXT_DYT =  DXT*DXT*DSTI**3
      ANYT_DXT = -DYT*DYT*DSTI**3
      ANYT_DYT =  DYT*DXT*DSTI**3
C
      DO J = 1, 2
        ANXT_XP(J) = ANXT_DXT*(XPT_XP(1,J)-XPT_XP(2,J))
        ANXT_YP(J) = ANXT_DYT*(YPT_YP(1,J)-YPT_YP(2,J))
C
        ANYT_XP(J) = ANYT_DXT*(XPT_XP(1,J)-XPT_XP(2,J))
        ANYT_YP(J) = ANYT_DYT*(YPT_YP(1,J)-YPT_YP(2,J))
      ENDDO
C
C
C---- form dot,cross products with adjacent-panel normal vectors
      DO K = 1, 2
        DOT(K) = ANXT*ANX(K) + ANYT*ANY(K)
        CRS(K) = ANXT*ANY(K) - ANYT*ANX(K)
        DOT_DX(K) = ANXT*ANX_DX(K) + ANYT*ANY_DX(K)
        DOT_DY(K) = ANXT*ANX_DY(K) + ANYT*ANY_DY(K)
        CRS_DX(K) = ANXT*ANY_DX(K) - ANYT*ANX_DX(K)
        CRS_DY(K) = ANXT*ANY_DY(K) - ANYT*ANX_DY(K)
        DO J = 1, 2
          DOT_XP(K,J) = ANXT_XP(J)*ANX(K) + ANYT_XP(J)*ANY(K)
          DOT_YP(K,J) = ANXT_YP(J)*ANX(K) + ANYT_YP(J)*ANY(K)
          CRS_XP(K,J) = ANXT_XP(J)*ANY(K) - ANYT_XP(J)*ANX(K)
          CRS_YP(K,J) = ANXT_YP(J)*ANY(K) - ANYT_YP(J)*ANX(K)
        ENDDO
      ENDDO
C
C
      DO K = 1, 2
        IF(K.EQ.1) THEN
         IP    = IP1
         IPOPP = IP2
         KOPP = 2
        ELSE
         IP    = IP2
         IPOPP = IP1
         KOPP = 1
        ENDIF
C        
        IF(IP.GT.0) THEN
         GAMT(K) = DOT(K)*GAM(IP) - CRS(K)*SIG(IP)
         SIGT(K) = CRS(K)*GAM(IP) + DOT(K)*SIG(IP)
         GAMT_SIG(K,K) = -CRS(K)
         GAMT_GAM(K,K) =  DOT(K)
         SIGT_SIG(K,K) =  DOT(K)
         SIGT_GAM(K,K) =  CRS(K)
         GAMT_DX(K) = DOT_DX(K)*GAM(IP) - CRS_DX(K)*SIG(IP)
         GAMT_DY(K) = DOT_DY(K)*GAM(IP) - CRS_DY(K)*SIG(IP)
         SIGT_DX(K) = CRS_DX(K)*GAM(IP) + DOT_DX(K)*SIG(IP)
         SIGT_DY(K) = CRS_DY(K)*GAM(IP) + DOT_DY(K)*SIG(IP)
         DO J = 1, 2
          GAMT_XP(K,J) = DOT_XP(K,J)*GAM(IP) - CRS_XP(K,J)*SIG(IP)
          GAMT_YP(K,J) = DOT_YP(K,J)*GAM(IP) - CRS_YP(K,J)*SIG(IP)
          SIGT_XP(K,J) = CRS_XP(K,J)*GAM(IP) + DOT_XP(K,J)*SIG(IP)
          SIGT_YP(K,J) = CRS_YP(K,J)*GAM(IP) + DOT_YP(K,J)*SIG(IP)
         ENDDO
        ELSE
         GAMT(K) = 0.
         SIGT(K) = CRS(KOPP)*GAM(IPOPP) + DOT(KOPP)*SIG(IPOPP)
         SIGT_SIG(K,KOPP) =  DOT(KOPP)
         SIGT_GAM(K,KOPP) =  CRS(KOPP)
         SIGT_DX(K) = CRS_DX(K)*GAM(IPOPP) + DOT_DX(K)*SIG(IPOPP)
         SIGT_DY(K) = CRS_DY(K)*GAM(IPOPP) + DOT_DY(K)*SIG(IPOPP)
         DO J = 1, 2
          SIGT_XP(K,J) = CRS_XP(K,J)*GAM(IPOPP) + DOT_XP(K,J)*SIG(IPOPP)
          SIGT_YP(K,J) = CRS_YP(K,J)*GAM(IPOPP) + DOT_YP(K,J)*SIG(IPOPP)
         ENDDO
        ENDIF
      ENDDO
C
      RETURN
      END ! TPANSET


      SUBROUTINE QFCALC(IF1,IF2, XF,YF, IPFO,IPFP,IFTYPE,
     &                  QF,QF_GAM,QF_SIG,QF_GTH)
C---------------------------------------------------------------
C     Calculates velocity components QF(.) at field points 
C     IF1...IF2 in coordinate arrays XF,YF w.r.t aero quantities 
C     at all nodes on all elements (1..NEL)
C---------------------------------------------------------------
C     Input:  XF,YF(i)  field points
C             IPFO(.)   first corner point index (0 if not on panel)
C             IPFP(.)   second corner point index (0 if not on panel)
C             IFTYPE(.) type for field point (0 if not on panel)
C
C     Assumes:
C             GAM(j)    panel vorticity  (or line circ. , or point doublet)
C             SIG(j)    panel source     (or line source, or point source )
C             GTH(j)    vortex wake panel vorticity
C             QINF      freestream speed
C
C             GAMT(j)   TE panel vorticity 
C             SIGT(j)   TE panel source    
C
C     Output: QF(.i)         velocity at XC(i),YC(i)
C             QF_GAM(.ji)    dQC(.i)/dGAM(j)
C             QF_SIG(.ji)
C             QF_GTH(.ji)    dQC(.i)/dGTH(j)
C---------------------------------------------------------------
C    Note:    This routine assumes that TE panel strengths are set
C             before entry.
C---------------------------------------------------------------
CHHY- Use these to save TE panel influences separately
cc      SUBROUTINE QFCALC(IF1,IF2, XF,YF, IPFO,IPFP,IFTYPE,
cc     &                  QF,QF_GAM,QF_SIG, QF_GAMT,QF_SIGT)
cc      DIMENSION          QF_GAMT(2,2,NEX,*), QF_SIGT(2,2,NEX,*)
C--------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION XF(*), YF(*)
      DIMENSION IPFO(*), IPFP(*), IFTYPE(*)
      DIMENSION QF(2,*), QF_GAM (2,IPX,*),   QF_SIG (2,IPX,*)
      DIMENSION          QF_GTH (2,IPX,*)
C
      DIMENSION QF_GAMTE(2,2,IPX), QF_SIGTE(2,2,IPX)
      DIMENSION RCORE(IPX)
C
      DO K = 1, 2
        DO IF = IF1, IF2
          QF(K,IF) = 0.
          DO IP = 1, NPTOT
            QF_GAM(K,IP,IF) = 0.
            QF_SIG(K,IP,IF) = 0.
            QF_GTH(K,IP,IF) = 0.
          ENDDO
CHHY- Uncomment to save TE panel influences separately
cc          DO IEL = 1, NEL
cc            QF_GAMT(K,1,IEL,IF) = 0.
cc            QF_SIGT(K,2,IEL,IF) = 0.
cc          ENDDO
        ENDDO
      ENDDO
C
      NF = IF2 - IF1 + 1
      KF = IF1
C
C---- Set contribution of all vortex/source lines, rings, and points
      DO 20 IEL = 1, NEL
        IF(    NETYPE(IEL).EQ.0
     &    .OR. NETYPE(IEL).EQ.1
     &    .OR. NETYPE(IEL).EQ.5
     &    .OR. NETYPE(IEL).EQ.6
     &    .OR. NETYPE(IEL).EQ.7) THEN
C------- sheets
         CALL PANAIC(IPFRST(IEL),IPLAST(IEL),XP,YP,
     &               NF,XF(KF),YF(KF), IPFO(KF),IPFP(KF), IFTYPE(KF),
     &               GAM,SIG,
     &               IPX, QF(1,KF), QF_GAM(1,1,KF), QF_SIG(1,1,KF) )
C
        ELSEIF(NETYPE(IEL).EQ.2) THEN
C------- line on axis
         DO IP = IPFRST(IEL), IPLAST(IEL)
           RCORE(IP) = 0.01
         ENDDO
         CALL AXLAIC(IPFRST(IEL),IPLAST(IEL),XP,YP,
     &               NF,XF(KF),YF(KF), IPFO(KF),IPFP(KF), IFTYPE(KF),
     &               GAM, SIG, RCORE,
     &               IPX, QF(1,KF), QF_GAM(1,1,KF), QF_SIG(1,1,KF) )
C     
        ELSEIF(NETYPE(IEL).EQ.3) THEN
C------- lines, rings
         CALL LINAIC(IPFRST(IEL),IPLAST(IEL), XP,YP,
     &               NF,XF(KF),YF(KF),
     &               GAM,SIG,
     &               IPX, QF(1,KF), QF_GAM(1,1,KF), QF_SIG(1,1,KF) )
C     
        ELSEIF(NETYPE(IEL).EQ.4) THEN
C------- point
         I = IPFRST(IEL)
         N = 1
         CALL PNTAIC(IPFRST(IEL),IPLAST(IEL), XP,YP,
     &               NF,XF(KF),YF(KF),
     &               GAM,SIG,
     &               IPX, QF(1,KF), QF_GAM(1,1,KF), QF_SIG(1,1,KF) )
        ENDIF
 20   CONTINUE
C
C
C
C---- Accumulate velocity from wake vorticity, GTH
C---- Save velocity AIC matrix (w/o TE panel influences) for wake vorticity 
      DO IEL = 1, NEL
        IF(    NETYPE(IEL).EQ.0
     &    .OR. NETYPE(IEL).EQ.7 ) THEN
C
          IP1 = IPFRST(IEL)
          IP2 = IPLAST(IEL)
          DO IF = IF1, IF2
            DO IP = IP1, IP2
              QF_GTH(1,IP,IF) = QF_GAM(1,IP,IF) 
              QF_GTH(2,IP,IF) = QF_GAM(2,IP,IF) 
              QF(1,IF) = QF(1,IF) + QF_GTH(1,IP,IF)*GTH(IP)
              QF(2,IF) = QF(2,IF) + QF_GTH(2,IP,IF)*GTH(IP)
            END DO
          END DO
C
        ENDIF
      END DO

C
C---- Add contribution of all TE panels
      DO 30 IEL = 1, NEL
        IF(LTPAN(IEL) .AND. 
     &     (NETYPE(IEL).EQ.0 .OR. NETYPE(IEL).EQ.7)) THEN
C---- clear field-point velocities for each TE panel
           DO IF = IF1, IF2
             DO KP = 1, 2
               QF_GAMTE(1,KP,IF) = 0.
               QF_GAMTE(2,KP,IF) = 0.
               QF_SIGTE(1,KP,IF) = 0.
               QF_SIGTE(2,KP,IF) = 0.
             ENDDO
           ENDDO
C
          CALL PANAIC(1,2,XPT(1,IEL),YPT(1,IEL),
     &                NF,XF(KF),YF(KF), IZERO,IZERO, IFTYPE(KF),
     &                GAMT(1,IEL), SIGT(1,IEL),
     &                2, QF(1,KF), QF_GAMTE(1,1,KF), QF_SIGTE(1,1,KF)) 
C
CHHY- Uncomment to save TE panel influences separately
C---- Save TE velocity influences 
cc          DO IF = IF1, IF2
cc            DO K = 1, 2
cc              QF_GAMT(K,1,IEL,IF) = QF_GAMTE(K,1,IF)
cc              QF_GAMT(K,2,IEL,IF) = QF_GAMTE(K,2,IF)
cc              QF_SIGT(K,1,IEL,IF) = QF_SIGTE(K,1,IF)
cc              QF_SIGT(K,2,IEL,IF) = QF_SIGTE(K,2,IF)
cc            END DO
cc          END DO
C
C---- Add TE influences to panel nodes (assumes TPANSET has been called) 
          DO IF = IF1, IF2
            DO KP = 1, 2
              IF(KP.EQ.1) THEN
               IP = IPTE1(IEL)
              ELSE
               IP = IPTE2(IEL)
              ENDIF
              IF(IP.GT.0) THEN
                DO K = 1, 2
                  QF_GAM(K,IP,IF) = QF_GAM(K,IP,IF)
     &                            + QF_GAMTE(K,1,IF)*GAMT_GAM(1,KP,IEL)
     &                            + QF_GAMTE(K,2,IF)*GAMT_GAM(2,KP,IEL)
     &                            + QF_SIGTE(K,1,IF)*SIGT_GAM(1,KP,IEL)
     &                            + QF_SIGTE(K,2,IF)*SIGT_GAM(2,KP,IEL)
                  QF_SIG(K,IP,IF) = QF_SIG(K,IP,IF)
     &                            + QF_GAMTE(K,1,IF)*GAMT_SIG(1,KP,IEL)
     &                            + QF_GAMTE(K,2,IF)*GAMT_SIG(2,KP,IEL)
     &                            + QF_SIGTE(K,1,IF)*SIGT_SIG(1,KP,IEL)
     &                            + QF_SIGTE(K,2,IF)*SIGT_SIG(2,KP,IEL)
                ENDDO ! K loop
              ENDIF
            END DO ! KP loop
          END DO ! IF loop
C
        ENDIF
C
C---- Special treatment for vortex wake last GTH point to include TE panel
        IF(LTPAN(IEL) .AND. NETYPE(IEL).EQ.7) THEN
          IP = IPLAST(IEL)
          DO IF = IF1, IF2
C---- Remove velocity w/o TE panel
cc            QF(1,IF) = QF(1,IF) - QF_GTH(1,IP,IF)*GTH(IP)
cc            QF(2,IF) = QF(2,IF) - QF_GTH(2,IP,IF)*GTH(IP)
C---- Save influence with TE panel
            QF_GTH(1,IP,IF) = QF_GAM(1,IP,IF) 
            QF_GTH(2,IP,IF) = QF_GAM(2,IP,IF) 
C---- Add back in velocity with TE panel
cc            QF(1,IF) = QF(1,IF) + QF_GTH(1,IP,IF)*GTH(IP)
cc            QF(2,IF) = QF(2,IF) + QF_GTH(2,IP,IF)*GTH(IP)
          END DO
        ENDIF

 30   ENDDO
C
C---- Add freestream
c      DO IF = 1, NF
c        QF(1,IF) = QF(1,IF) + QINF
c      ENDDO
C
      RETURN
      END ! QFCALC
