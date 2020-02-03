module m_qaic
    implicit none
contains
    !*==QAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PANGET
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

    SUBROUTINE QAIC(LXYJAC)
        use i_dfdc
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        LOGICAL :: LXYJAC
        !
        ! Local variables
        !
        INTEGER :: ICE, KEL
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        IF (LDBG) THEN
            IF (LXYJAC) THEN
                WRITE (*, *)                                                 &
                        &'Generating dq/d(gam,sig,x,y) Jacobian matrices...'
            ELSE
                WRITE (*, *) 'Generating dq/d(gam,sig) Jacobian matrices...'
            ENDIF
        ENDIF
        !
        CALL QAIC1(LXYJAC, 1, NCTOT, 1, NEL)
        DO KEL = 1, NEL
            ICE = NCX + KEL
            IF (LBODY(KEL)) CALL QAIC1(LXYJAC, ICE, ICE, 1, NEL)
        ENDDO
        !
        LQCNT = .TRUE.
        LQAIC = .TRUE.
        LQGIC = LXYJAC
        !
    END SUBROUTINE QAIC
    !*==QAIC1.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QAIC



    SUBROUTINE QAIC1(LXYJAC, IC1, IC2, IEL1, IEL2)
        use i_dfdc
        use m_aic, only : pntaic, axlaic, panaic, panaicd, linaic
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IC1, IC2, IEL1, IEL2
        LOGICAL :: LXYJAC
        !
        ! Local variables
        !
        REAL :: DSPT, DXPT, DYPT, SPLEN
        INTEGER :: IC, IEL, IP, IP1, IP2, IPO, IPP, J, JPM, &
                & JPO, K, KC, KP, NUMIC
        REAL, DIMENSION(2, 2, ICX) :: QC_GAMTE, QC_SIGTE, QC_XPT, &
                & QC_YPT
        REAL, DIMENSION(IPX) :: RCORE
        REAL, SAVE :: TFRAC
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------
        !     Computes velocities and Jacobians at control points IC1...IC2
        !     w.r.t aero quantities at nodes on elements IEL1..IEL2
        !
        !     If LXYJAC=T, also computes Jacobians w.r.t. panel geometry
        !
        !     Input:  XC,YC(i)  control points
        !             XP,YP(j)  panel geometry
        !             GAM(j)    panel vorticity  (or line circ. , or point doublet)
        !             SIG(j)    panel source     (or line source, or point source )
        !             GTH(j)    vortex wake panel vorticity
        !             QINF      freestream speed
        !             ALFA      freestream angle
        !
        !     Output: QC(.i)         velocity at XC(i),YC(i)
        !             QC_GAM(.ji)    dQC(.i)/dGAM(j)
        !             QC_SIG(.ji)
        !             QC_GTH(.ji)    dQC(.i)/dGTH(j)
        !
        !     Additonal output if LXYJAC=T:
        !             QC_XP(.ji)
        !             QC_YP(.ji)
        !             QC_XC(.i)
        !             QC_YC(.i)
        !
        !     Note:   TE panel information and strengths are set by
        !             this routine
        !---------------------------------------------------------------
        !
        !
        !
        !---- maximum TE gap (relative to perimeter) to give closed airfoil a TE panel
        DATA TFRAC/0.00001/
        !
        !---- clear control-point velocities for accumulation in PANAIC and PNTAIC
        DO IC = IC1, IC2
            QC(1, IC) = 0.
            QC(2, IC) = 0.
            DO IEL = IEL1, IEL2
                DO IP = IPFRST(IEL), IPLAST(IEL)
                    QC_GAM(1, IP, IC) = 0.
                    QC_GAM(2, IP, IC) = 0.
                    QC_SIG(1, IP, IC) = 0.
                    QC_SIG(2, IP, IC) = 0.
                    QC_XP(1, IP, IC) = 0.
                    QC_XP(2, IP, IC) = 0.
                    QC_YP(1, IP, IC) = 0.
                    QC_YP(2, IP, IC) = 0.
                    QC_GTH(1, IP, IC) = 0.
                    QC_GTH(2, IP, IC) = 0.
                ENDDO
                !HHY- Uncomment to save TE panel influences separately
                !c          DO KP = 1, 2
                !c            QC_GAMT(1,KP,IEL,IC) = 0.
                !c            QC_GAMT(2,KP,IEL,IC) = 0.
                !c            QC_SIGT(1,KP,IEL,IC) = 0.
                !c            QC_SIGT(2,KP,IEL,IC) = 0.
                !c          ENDDO
            ENDDO
        ENDDO
        !
        KC = IC1
        NUMIC = IC2 - IC1 + 1
        !
        !---- Set AIC matrix for all singularity types
        DO IEL = IEL1, IEL2
            IF (NETYPE(IEL)==0 .OR. NETYPE(IEL)==1 .OR. NETYPE(IEL)       &
                    & ==5 .OR. NETYPE(IEL)==6 .OR. NETYPE(IEL)==7) THEN
                !-------- wall,wake,source-only line and vortex wake sheets
                IF (LXYJAC) THEN
                    !--------- linearize QC w.r.t GAM,SIG,XP,YP,XC,YC
                    CALL PANAICD(IPFRST(IEL), IPLAST(IEL), XP, YP, NUMIC, XC(KC), &
                            & YC(KC), IPCO(KC), IPCP(KC), ICTYPE(KC), GAM, SIG, &
                            & IPX, QC(1, KC), QC_GAM(1, 1, KC), QC_SIG(1, 1, KC), &
                            & QC_XP(1, 1, KC), QC_YP(1, 1, KC), QC_XC(1, KC), &
                            & QC_YC(1, KC))
                    DO IC = IC1, IC2
                        IPO = IPCO(IC)
                        IPP = IPCP(IC)
                        IF (IPO/=0 .AND. IPP/=0) THEN
                            QC_XP(1, IPO, IC) = QC_XP(1, IPO, IC) + QC_XC(1, IC) * 0.5
                            QC_XP(2, IPO, IC) = QC_XP(2, IPO, IC) + QC_XC(2, IC) * 0.5
                            QC_XP(1, IPP, IC) = QC_XP(1, IPP, IC) + QC_XC(1, IC) * 0.5
                            QC_XP(2, IPP, IC) = QC_XP(2, IPP, IC) + QC_XC(2, IC) * 0.5
                            QC_YP(1, IPO, IC) = QC_YP(1, IPO, IC) + QC_YC(1, IC) * 0.5
                            QC_YP(2, IPO, IC) = QC_YP(2, IPO, IC) + QC_YC(2, IC) * 0.5
                            QC_YP(1, IPP, IC) = QC_YP(1, IPP, IC) + QC_YC(1, IC) * 0.5
                            QC_YP(2, IPP, IC) = QC_YP(2, IPP, IC) + QC_YC(2, IC) * 0.5
                        ENDIF
                    ENDDO
                ELSE
                    !--------- linearize QC w.r.t GAM,SIG only
                    CALL PANAIC(IPFRST(IEL), IPLAST(IEL), XP, YP, NUMIC, XC(KC), &
                            & YC(KC), IPCO(KC), IPCP(KC), ICTYPE(KC), GAM, SIG, &
                            & IPX, QC(1, KC), QC_GAM(1, 1, KC), QC_SIG(1, 1, KC))
                ENDIF
                !
            ELSEIF (NETYPE(IEL)==2) THEN
                !------- line on axis
                DO IP = IPFRST(IEL), IPLAST(IEL)
                    RCORE(IP) = 0.01
                ENDDO
                CALL AXLAIC(IPFRST(IEL), IPLAST(IEL), XP, YP, NUMIC, XC(KC), &
                        & YC(KC), IPCO(KC), IPCP(KC), ICTYPE(KC), GAM, SIG, &
                        & RCORE, IPX, QC(1, KC), QC_GAM(1, 1, KC), QC_SIG(1, 1, KC)&
                        &)
                !
            ELSEIF (NETYPE(IEL)==3) THEN
                !------- ring
                CALL LINAIC(IPFRST(IEL), IPLAST(IEL), XP, YP, NUMIC, XC(KC), &
                        & YC(KC), GAM, SIG, IPX, QC(1, KC), QC_GAM(1, 1, KC), &
                        & QC_SIG(1, 1, KC))
                !
            ELSEIF (NETYPE(IEL)==4) THEN
                !------- point
                CALL PNTAIC(IPFRST(IEL), IPLAST(IEL), XP, YP, NUMIC, XC(KC), &
                        & YC(KC), GAM, SIG, IPX, QC(1, KC), QC_GAM(1, 1, KC), &
                        & QC_SIG(1, 1, KC))
            ENDIF
            !
        ENDDO
        !
        !---- Accumulate velocity from wake vorticity, GTH
        !---- Save velocity AIC matrix (w/o TE panel influences) for wake vorticity
        DO IEL = IEL1, IEL2
            IF (NETYPE(IEL)==0 .OR. NETYPE(IEL)==7) THEN
                !
                IP1 = IPFRST(IEL)
                IP2 = IPLAST(IEL)
                DO IC = IC1, IC2
                    DO IP = IP1, IP2
                        QC_GTH(1, IP, IC) = QC_GAM(1, IP, IC)
                        QC_GTH(2, IP, IC) = QC_GAM(2, IP, IC)
                        QC(1, IC) = QC(1, IC) + QC_GTH(1, IP, IC) * GTH(IP)
                        QC(2, IC) = QC(2, IC) + QC_GTH(2, IP, IC) * GTH(IP)
                    ENDDO
                ENDDO
                !
            ENDIF
        ENDDO
        !
        !---- Add contributions of TE panels...
        DO IEL = IEL1, IEL2
            IF (.NOT.(LTPAN(IEL) .AND. (NETYPE(IEL)==0 .OR. NETYPE(IEL)==7&
                    &))) THEN
                !------- no TE panel... zero out panel strengths
                DO K = 1, 2
                    XPT(K, IEL) = 0.
                    YPT(K, IEL) = 0.
                    GAMT(K, IEL) = 0.
                    SIGT(K, IEL) = 0.
                    DO J = 1, 2
                        GAMT_GAM(K, J, IEL) = 0.
                        GAMT_SIG(K, J, IEL) = 0.
                        SIGT_GAM(K, J, IEL) = 0.
                        SIGT_SIG(K, J, IEL) = 0.
                    ENDDO
                ENDDO
                CYCLE
            ENDIF
            !---- clear accumulating control-point velocities for each TE panel
            DO IC = IC1, IC2
                DO KP = 1, 2
                    QC_GAMTE(1, KP, IC) = 0.
                    QC_GAMTE(2, KP, IC) = 0.
                    QC_SIGTE(1, KP, IC) = 0.
                    QC_SIGTE(2, KP, IC) = 0.
                    QC_XPT(1, KP, IC) = 0.
                    QC_XPT(2, KP, IC) = 0.
                    QC_YPT(1, KP, IC) = 0.
                    QC_YPT(2, KP, IC) = 0.
                ENDDO
                QC_XC(1, IC) = 0.
                QC_XC(2, IC) = 0.
                QC_YC(1, IC) = 0.
                QC_YC(2, IC) = 0.
            ENDDO
            !
            !------ nodes on element spanned by TE panel
            IP1 = IPTE1(IEL)
            IP2 = IPTE2(IEL)
            !
            !------ set panel coordinates, strengths, and all Jacobians
            CALL TPANSET(IP1, IP2, GAM, SIG, XP, YP, XPT(1, IEL), XPT_XP(1, 1, IEL), &
                    & YPT(1, IEL), YPT_YP(1, 1, IEL), GAMT(1, IEL), &
                    & GAMT_GAM(1, 1, IEL), GAMT_SIG(1, 1, IEL), &
                    & GAMT_XP(1, 1, IEL), GAMT_YP(1, 1, IEL), GAMT_DX(1, IEL), &
                    & GAMT_DY(1, IEL), SIGT(1, IEL), SIGT_GAM(1, 1, IEL), &
                    & SIGT_SIG(1, 1, IEL), SIGT_XP(1, 1, IEL), &
                    & SIGT_YP(1, 1, IEL), SIGT_DX(1, IEL), SIGT_DY(1, IEL))
            !
            !------ skip influence computation if TE panel has negligible length
            DXPT = XPT(2, IEL) - XPT(1, IEL)
            DYPT = YPT(2, IEL) - YPT(1, IEL)
            DSPT = SQRT(DXPT**2 + DYPT**2)
            !
            SPLEN = ABS(SP(IPLAST(IEL)) - SP(IPFRST(IEL)))
            IF (DSPT<TFRAC * SPLEN) THEN
                WRITE (*, *) 'Negligible TE panel ignored for element', IEL
                CYCLE
            ENDIF
            !
            !
            IF (LXYJAC) THEN
                !-------- linearize QC w.r.t GAM,SIG,XP,YP,XC,YC
                CALL PANAICD(1, 2, XPT(1, IEL), YPT(1, IEL), NUMIC, XC(KC), YC(KC), &
                        & IZERO, IZERO, ICTYPE(KC), GAMT(1, IEL), SIGT(1, IEL), &
                        & 2, QC(1, KC), QC_GAMTE(1, 1, KC), QC_SIGTE(1, 1, KC), &
                        & QC_XPT(1, 1, KC), QC_YPT(1, 1, KC), QC_XC(1, KC), &
                        & QC_YC(1, KC))
                DO IC = IC1, IC2
                    IPO = IPCO(IC)
                    IPP = IPCP(IC)
                    IF (IPO/=0 .AND. IPP/=0) THEN
                        QC_XP(1, IPO, IC) = QC_XP(1, IPO, IC) + QC_XC(1, IC) * 0.5
                        QC_XP(2, IPO, IC) = QC_XP(2, IPO, IC) + QC_XC(2, IC) * 0.5
                        QC_XP(1, IPP, IC) = QC_XP(1, IPP, IC) + QC_XC(1, IC) * 0.5
                        QC_XP(2, IPP, IC) = QC_XP(2, IPP, IC) + QC_XC(2, IC) * 0.5
                        QC_YP(1, IPO, IC) = QC_YP(1, IPO, IC) + QC_YC(1, IC) * 0.5
                        QC_YP(2, IPO, IC) = QC_YP(2, IPO, IC) + QC_YC(2, IC) * 0.5
                        QC_YP(1, IPP, IC) = QC_YP(1, IPP, IC) + QC_YC(1, IC) * 0.5
                        QC_YP(2, IPP, IC) = QC_YP(2, IPP, IC) + QC_YC(2, IC) * 0.5
                    ENDIF
                ENDDO
                !
            ELSE
                !-------- linearize QC w.r.t GAM,SIG only
                CALL PANAIC(1, 2, XPT(1, IEL), YPT(1, IEL), NUMIC, XC(KC), YC(KC), &
                        & IZERO, IZERO, ICTYPE(KC), GAMT(1, IEL), SIGT(1, IEL), &
                        & 2, QC(1, KC), QC_GAMTE(1, 1, KC), QC_SIGTE(1, 1, KC))
            ENDIF
            !
            !HHY- Uncomment to save TE panel influences separately
            !---- Save TE panel velocity influences
            !c        DO K = 1, 2
            !c          DO IC = IC1, IC2
            !c            QC_GAMT(K,1,IEL,IC) = QC_GAMTE(K,1,IC)
            !c            QC_GAMT(K,2,IEL,IC) = QC_GAMTE(K,2,IC)
            !c            QC_SIGT(K,1,IEL,IC) = QC_SIGTE(K,1,IC)
            !c            QC_SIGT(K,2,IEL,IC) = QC_SIGTE(K,2,IC)
            !c          END DO
            !c        END DO
            !
            DO IC = IC1, IC2
                DO KP = 1, 2
                    IF (KP==1) THEN
                        IP = IP1
                        JPO = IP1 + 1
                        JPM = IP1
                    ELSE
                        IP = IP2
                        JPO = IP2
                        JPM = IP2 - 1
                    ENDIF
                    IF (IP>0) THEN
                        DO K = 1, 2
                            QC_GAM(K, IP, IC) = QC_GAM(K, IP, IC)                  &
                                    & + QC_GAMTE(K, 1, IC)               &
                                            & * GAMT_GAM(1, KP, IEL)              &
                                    & + QC_GAMTE(K, 2, IC)               &
                                            & * GAMT_GAM(2, KP, IEL)              &
                                    & + QC_SIGTE(K, 1, IC)               &
                                            & * SIGT_GAM(1, KP, IEL)              &
                                    & + QC_SIGTE(K, 2, IC)               &
                                            & * SIGT_GAM(2, KP, IEL)
                            QC_SIG(K, IP, IC) = QC_SIG(K, IP, IC)                  &
                                    & + QC_GAMTE(K, 1, IC)               &
                                            & * GAMT_SIG(1, KP, IEL)              &
                                    & + QC_GAMTE(K, 2, IC)               &
                                            & * GAMT_SIG(2, KP, IEL)              &
                                    & + QC_SIGTE(K, 1, IC)               &
                                            & * SIGT_SIG(1, KP, IEL)              &
                                    & + QC_SIGTE(K, 2, IC)               &
                                            & * SIGT_SIG(2, KP, IEL)
                            !
                            IF (LXYJAC) THEN
                                QC_XP(K, IP, IC) = QC_XP(K, IP, IC) + QC_XPT(K, 1, IC)&
                                        & * XPT_XP(1, KP, IEL) + QC_XPT(K, 2, IC)           &
                                        & * XPT_XP(2, KP, IEL) + QC_GAMTE(K, 1, IC)         &
                                        & * GAMT_XP(1, KP, IEL) + QC_GAMTE(K, 2, IC)        &
                                        & * GAMT_XP(2, KP, IEL) + QC_SIGTE(K, 1, IC)        &
                                        & * SIGT_XP(1, KP, IEL) + QC_SIGTE(K, 2, IC)        &
                                        & * SIGT_XP(2, KP, IEL)
                                QC_YP(K, IP, IC) = QC_YP(K, IP, IC) + QC_YPT(K, 1, IC)&
                                        & * YPT_YP(1, KP, IEL) + QC_YPT(K, 2, IC)           &
                                        & * YPT_YP(2, KP, IEL) + QC_GAMTE(K, 1, IC)         &
                                        & * GAMT_YP(1, KP, IEL) + QC_GAMTE(K, 2, IC)        &
                                        & * GAMT_YP(2, KP, IEL) + QC_SIGTE(K, 1, IC)        &
                                        & * SIGT_YP(1, KP, IEL) + QC_SIGTE(K, 2, IC)        &
                                        & * SIGT_YP(2, KP, IEL)
                                QC_XP(K, JPO, IC) = QC_XP(K, JPO, IC)               &
                                        & + QC_GAMTE(K, KP, IC) * GAMT_DX(KP, IEL)          &
                                        & + QC_SIGTE(K, KP, IC) * SIGT_DX(KP, IEL)
                                QC_YP(K, JPO, IC) = QC_YP(K, JPO, IC)               &
                                        & + QC_GAMTE(K, KP, IC) * GAMT_DY(KP, IEL)          &
                                        & + QC_SIGTE(K, KP, IC) * SIGT_DY(KP, IEL)
                                QC_XP(K, JPM, IC) = QC_XP(K, JPM, IC)               &
                                        & - QC_GAMTE(K, KP, IC) * GAMT_DX(KP, IEL)          &
                                        & - QC_SIGTE(K, KP, IC) * SIGT_DX(KP, IEL)
                                QC_YP(K, JPM, IC) = QC_YP(K, JPM, IC)               &
                                        & - QC_GAMTE(K, KP, IC) * GAMT_DY(KP, IEL)          &
                                        & - QC_SIGTE(K, KP, IC) * SIGT_DY(KP, IEL)
                            ENDIF
                        ENDDO   ! K loop
                    ENDIF
                ENDDO  ! KP loop
            ENDDO ! IC loop
            !
            !---- Special treatment for vortex wake last GTH point to include TE panel
            IF (LTPAN(IEL) .AND. NETYPE(IEL)==7) THEN
                IP = IPLAST(IEL)
                DO IC = IC1, IC2
                    !---- Remove velocity w/o TE panel
                    QC(1, IC) = QC(1, IC) - QC_GTH(1, IP, IC) * GTH(IP)
                    QC(2, IC) = QC(2, IC) - QC_GTH(2, IP, IC) * GTH(IP)
                    !---- Save influence with TE panel
                    QC_GTH(1, IP, IC) = QC_GAM(1, IP, IC)
                    QC_GTH(2, IP, IC) = QC_GAM(2, IP, IC)
                    !---- Add back in velocity with TE panel
                    QC(1, IC) = QC(1, IC) + QC_GTH(1, IP, IC) * GTH(IP)
                    QC(2, IC) = QC(2, IC) + QC_GTH(2, IP, IC) * GTH(IP)
                ENDDO
            ENDIF
            !
        ENDDO
        !
        !---- add on freestream
        DO IC = IC1, IC2
            QC(1, IC) = QC(1, IC) + QINF
        ENDDO
        !
    END SUBROUTINE QAIC1
    !*==TPANSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QAIC1



    SUBROUTINE TPANSET(IPTE1, IPTE2, GAM, SIG, XP, YP, XPT, XPT_XP, YPT, &
            & YPT_YP, GAMT, GAMT_GAM, GAMT_SIG, GAMT_XP, GAMT_YP, &
            & GAMT_DX, GAMT_DY, SIGT, SIGT_GAM, SIGT_SIG, SIGT_XP, &
            & SIGT_YP, SIGT_DX, SIGT_DY)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IPTE1, IPTE2
        REAL, DIMENSION(*) :: GAM, SIG, XP, YP
        REAL, DIMENSION(2) :: GAMT, GAMT_DX, GAMT_DY, SIGT, &
                & SIGT_DX, SIGT_DY, XPT, YPT
        REAL, DIMENSION(2, 2) :: GAMT_GAM, GAMT_SIG, GAMT_XP, &
                & GAMT_YP, SIGT_GAM, SIGT_SIG, &
                & SIGT_XP, SIGT_YP, XPT_XP, YPT_YP
        !
        ! Local variables
        !
        REAL, DIMENSION(2) :: ANX, ANXT_XP, ANXT_YP, ANX_DX, &
                & ANX_DY, ANY, ANYT_XP, ANYT_YP, &
                & ANY_DX, ANY_DY, CRS, CRS_DX, CRS_DY, &
                & DOT, DOT_DX, DOT_DY
        REAL :: ANXT, ANXT_DXT, ANXT_DYT, ANYAXI, ANYT, ANYT_DXT, &
                & ANYT_DYT, DS, DSTI, DSTSQ, DX, DXT, DY, DYT
        REAL, DIMENSION(2, 2) :: CRS_XP, CRS_YP, DOT_XP, DOT_YP
        INTEGER :: IP, IP1, IP2, IPOPP, J, JPM, JPO, K, KOPP
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------
        !     Sets geometry, strengths, and associated Jacobians
        !     of TE panel on element IEL
        !----------------------------------------------------------
        !
        !
        !---- local setup arrays
        !
        !
        DO K = 1, 2
            DO J = 1, 2
                XPT_XP(K, J) = 0.
                YPT_YP(K, J) = 0.
                GAMT_GAM(K, J) = 0.
                GAMT_SIG(K, J) = 0.
                SIGT_GAM(K, J) = 0.
                SIGT_SIG(K, J) = 0.
                GAMT_XP(K, J) = 0.
                GAMT_YP(K, J) = 0.
                SIGT_XP(K, J) = 0.
                SIGT_YP(K, J) = 0.
            ENDDO
            ANX_DX(K) = 0.
            ANY_DX(K) = 0.
            ANX_DY(K) = 0.
            ANY_DY(K) = 0.
        ENDDO
        !
        !---- geometry nodes spanned by TE panel
        !-    (if node index = 0, then that TE panel node lies on axis)
        IP1 = IPTE1
        IP2 = IPTE2
        !
        !---- set node x,y of TE panel, and dx,dy of adjacent panels
        !
        DO K = 1, 2
            IF (K==1) THEN
                IP = IP1
                IPOPP = IP2
                KOPP = 2
                JPO = IP + 1
                JPM = IP
                ANYAXI = 1.0
            ELSE
                IP = IP2
                IPOPP = IP1
                KOPP = 1
                JPO = IP
                JPM = IP - 1
                ANYAXI = -1.0
            ENDIF
            !
            IF (IP>0) THEN
                !------- usual TE panel
                XPT(K) = XP(IP)
                YPT(K) = YP(IP)
                XPT_XP(K, K) = 1.0
                YPT_YP(K, K) = 1.0
                DX = XP(JPO) - XP(JPM)
                DY = YP(JPO) - YP(JPM)
                DS = SQRT(DX**2 + DY**2)
                ANX(K) = DY / DS
                ANY(K) = -DX / DS
                ANX_DX(K) = -DX * DY / DS**3
                ANX_DY(K) = DX * DX / DS**3
                ANY_DX(K) = -DY * DY / DS**3
                ANY_DY(K) = DY * DX / DS**3
            ELSE
                !------- node IP is on axis, at same x location as node IPOPP
                XPT(K) = XP(IPOPP)
                YPT(K) = 0.
                XPT_XP(K, KOPP) = 1.0
                ANX(K) = 0.
                ANY(K) = ANYAXI
            ENDIF
        ENDDO
        !
        !---- now set TE panel normal vector ANXT,ANYT
        DXT = XPT(1) - XPT(2)
        DYT = YPT(1) - YPT(2)
        DSTSQ = DXT**2 + DYT**2
        IF (DSTSQ==0.0) THEN
            DSTI = 1.0
        ELSE
            DSTI = 1.0 / SQRT(DSTSQ)
        ENDIF
        ANXT = DYT * DSTI
        ANYT = -DXT * DSTI
        ANXT_DXT = -DXT * DYT * DSTI**3
        ANXT_DYT = DXT * DXT * DSTI**3
        ANYT_DXT = -DYT * DYT * DSTI**3
        ANYT_DYT = DYT * DXT * DSTI**3
        !
        DO J = 1, 2
            ANXT_XP(J) = ANXT_DXT * (XPT_XP(1, J) - XPT_XP(2, J))
            ANXT_YP(J) = ANXT_DYT * (YPT_YP(1, J) - YPT_YP(2, J))
            !
            ANYT_XP(J) = ANYT_DXT * (XPT_XP(1, J) - XPT_XP(2, J))
            ANYT_YP(J) = ANYT_DYT * (YPT_YP(1, J) - YPT_YP(2, J))
        ENDDO
        !
        !
        !---- form dot,cross products with adjacent-panel normal vectors
        DO K = 1, 2
            DOT(K) = ANXT * ANX(K) + ANYT * ANY(K)
            CRS(K) = ANXT * ANY(K) - ANYT * ANX(K)
            DOT_DX(K) = ANXT * ANX_DX(K) + ANYT * ANY_DX(K)
            DOT_DY(K) = ANXT * ANX_DY(K) + ANYT * ANY_DY(K)
            CRS_DX(K) = ANXT * ANY_DX(K) - ANYT * ANX_DX(K)
            CRS_DY(K) = ANXT * ANY_DY(K) - ANYT * ANX_DY(K)
            DO J = 1, 2
                DOT_XP(K, J) = ANXT_XP(J) * ANX(K) + ANYT_XP(J) * ANY(K)
                DOT_YP(K, J) = ANXT_YP(J) * ANX(K) + ANYT_YP(J) * ANY(K)
                CRS_XP(K, J) = ANXT_XP(J) * ANY(K) - ANYT_XP(J) * ANX(K)
                CRS_YP(K, J) = ANXT_YP(J) * ANY(K) - ANYT_YP(J) * ANX(K)
            ENDDO
        ENDDO
        !
        !
        DO K = 1, 2
            IF (K==1) THEN
                IP = IP1
                IPOPP = IP2
                KOPP = 2
            ELSE
                IP = IP2
                IPOPP = IP1
                KOPP = 1
            ENDIF
            !
            IF (IP>0) THEN
                GAMT(K) = DOT(K) * GAM(IP) - CRS(K) * SIG(IP)
                SIGT(K) = CRS(K) * GAM(IP) + DOT(K) * SIG(IP)
                GAMT_SIG(K, K) = -CRS(K)
                GAMT_GAM(K, K) = DOT(K)
                SIGT_SIG(K, K) = DOT(K)
                SIGT_GAM(K, K) = CRS(K)
                GAMT_DX(K) = DOT_DX(K) * GAM(IP) - CRS_DX(K) * SIG(IP)
                GAMT_DY(K) = DOT_DY(K) * GAM(IP) - CRS_DY(K) * SIG(IP)
                SIGT_DX(K) = CRS_DX(K) * GAM(IP) + DOT_DX(K) * SIG(IP)
                SIGT_DY(K) = CRS_DY(K) * GAM(IP) + DOT_DY(K) * SIG(IP)
                DO J = 1, 2
                    GAMT_XP(K, J) = DOT_XP(K, J) * GAM(IP) - CRS_XP(K, J) * SIG(IP)
                    GAMT_YP(K, J) = DOT_YP(K, J) * GAM(IP) - CRS_YP(K, J) * SIG(IP)
                    SIGT_XP(K, J) = CRS_XP(K, J) * GAM(IP) + DOT_XP(K, J) * SIG(IP)
                    SIGT_YP(K, J) = CRS_YP(K, J) * GAM(IP) + DOT_YP(K, J) * SIG(IP)
                ENDDO
            ELSE
                GAMT(K) = 0.
                SIGT(K) = CRS(KOPP) * GAM(IPOPP) + DOT(KOPP) * SIG(IPOPP)
                SIGT_SIG(K, KOPP) = DOT(KOPP)
                SIGT_GAM(K, KOPP) = CRS(KOPP)
                SIGT_DX(K) = CRS_DX(K) * GAM(IPOPP) + DOT_DX(K) * SIG(IPOPP)
                SIGT_DY(K) = CRS_DY(K) * GAM(IPOPP) + DOT_DY(K) * SIG(IPOPP)
                DO J = 1, 2
                    SIGT_XP(K, J) = CRS_XP(K, J) * GAM(IPOPP) + DOT_XP(K, J)      &
                            & * SIG(IPOPP)
                    SIGT_YP(K, J) = CRS_YP(K, J) * GAM(IPOPP) + DOT_YP(K, J)      &
                            & * SIG(IPOPP)
                ENDDO
            ENDIF
        ENDDO
        !
    END SUBROUTINE TPANSET
    !*==QFCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! TPANSET


    SUBROUTINE QFCALC(IF1, IF2, XF, YF, IPFO, IPFP, IFTYPE, QF, QF_GAM, QF_SIG, &
            & QF_GTH)
        use i_dfdc
        use m_aic, only : axlaic, linaic, panaic, pntaic
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IF1, IF2
        INTEGER, DIMENSION(*) :: IFTYPE, IPFO, IPFP
        REAL, DIMENSION(2, *) :: QF
        REAL, DIMENSION(2, IPX, *) :: QF_GAM, QF_GTH, QF_SIG
        REAL, DIMENSION(*) :: XF, YF
        !
        ! Local variables
        !
        INTEGER :: I, IEL, IF, IP, IP1, IP2, K, KF, KP, N, NF
        REAL, DIMENSION(2, 2, IPX) :: QF_GAMTE, QF_SIGTE
        REAL, DIMENSION(IPX) :: RCORE
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------
        !     Calculates velocity components QF(.) at field points
        !     IF1...IF2 in coordinate arrays XF,YF w.r.t aero quantities
        !     at all nodes on all elements (1..NEL)
        !---------------------------------------------------------------
        !     Input:  XF,YF(i)  field points
        !             IPFO(.)   first corner point index (0 if not on panel)
        !             IPFP(.)   second corner point index (0 if not on panel)
        !             IFTYPE(.) type for field point (0 if not on panel)
        !
        !     Assumes:
        !             GAM(j)    panel vorticity  (or line circ. , or point doublet)
        !             SIG(j)    panel source     (or line source, or point source )
        !             GTH(j)    vortex wake panel vorticity
        !             QINF      freestream speed
        !
        !             GAMT(j)   TE panel vorticity
        !             SIGT(j)   TE panel source
        !
        !     Output: QF(.i)         velocity at XC(i),YC(i)
        !             QF_GAM(.ji)    dQC(.i)/dGAM(j)
        !             QF_SIG(.ji)
        !             QF_GTH(.ji)    dQC(.i)/dGTH(j)
        !---------------------------------------------------------------
        !    Note:    This routine assumes that TE panel strengths are set
        !             before entry.
        !---------------------------------------------------------------
        !HHY- Use these to save TE panel influences separately
        !c      SUBROUTINE QFCALC(IF1,IF2, XF,YF, IPFO,IPFP,IFTYPE,
        !c     &                  QF,QF_GAM,QF_SIG, QF_GAMT,QF_SIGT)
        !c      DIMENSION          QF_GAMT(2,2,NEX,*), QF_SIGT(2,2,NEX,*)
        !--------------------------------------------------------
        !
        !
        DO K = 1, 2
            DO IF = IF1, IF2
                QF(K, IF) = 0.
                DO IP = 1, NPTOT
                    QF_GAM(K, IP, IF) = 0.
                    QF_SIG(K, IP, IF) = 0.
                    QF_GTH(K, IP, IF) = 0.
                ENDDO
                !HHY- Uncomment to save TE panel influences separately
                !c          DO IEL = 1, NEL
                !c            QF_GAMT(K,1,IEL,IF) = 0.
                !c            QF_SIGT(K,2,IEL,IF) = 0.
                !c          ENDDO
            ENDDO
        ENDDO
        !
        NF = IF2 - IF1 + 1
        KF = IF1
        !
        !---- Set contribution of all vortex/source lines, rings, and points
        DO IEL = 1, NEL
            IF (NETYPE(IEL)==0 .OR. NETYPE(IEL)==1 .OR. NETYPE(IEL)       &
                    & ==5 .OR. NETYPE(IEL)==6 .OR. NETYPE(IEL)==7) THEN
                !------- sheets
                CALL PANAIC(IPFRST(IEL), IPLAST(IEL), XP, YP, NF, XF(KF), YF(KF), &
                        & IPFO(KF), IPFP(KF), IFTYPE(KF), GAM, SIG, IPX, &
                        & QF(1, KF), QF_GAM(1, 1, KF), QF_SIG(1, 1, KF))
                !
            ELSEIF (NETYPE(IEL)==2) THEN
                !------- line on axis
                DO IP = IPFRST(IEL), IPLAST(IEL)
                    RCORE(IP) = 0.01
                ENDDO
                CALL AXLAIC(IPFRST(IEL), IPLAST(IEL), XP, YP, NF, XF(KF), YF(KF), &
                        & IPFO(KF), IPFP(KF), IFTYPE(KF), GAM, SIG, RCORE, IPX, &
                        & QF(1, KF), QF_GAM(1, 1, KF), QF_SIG(1, 1, KF))
                !
            ELSEIF (NETYPE(IEL)==3) THEN
                !------- lines, rings
                CALL LINAIC(IPFRST(IEL), IPLAST(IEL), XP, YP, NF, XF(KF), YF(KF), &
                        & GAM, SIG, IPX, QF(1, KF), QF_GAM(1, 1, KF), &
                        & QF_SIG(1, 1, KF))
                !
            ELSEIF (NETYPE(IEL)==4) THEN
                !------- point
                I = IPFRST(IEL)
                N = 1
                CALL PNTAIC(IPFRST(IEL), IPLAST(IEL), XP, YP, NF, XF(KF), YF(KF), &
                        & GAM, SIG, IPX, QF(1, KF), QF_GAM(1, 1, KF), &
                        & QF_SIG(1, 1, KF))
            ENDIF
        ENDDO
        !
        !
        !
        !---- Accumulate velocity from wake vorticity, GTH
        !---- Save velocity AIC matrix (w/o TE panel influences) for wake vorticity
        DO IEL = 1, NEL
            IF (NETYPE(IEL)==0 .OR. NETYPE(IEL)==7) THEN
                !
                IP1 = IPFRST(IEL)
                IP2 = IPLAST(IEL)
                DO IF = IF1, IF2
                    DO IP = IP1, IP2
                        QF_GTH(1, IP, IF) = QF_GAM(1, IP, IF)
                        QF_GTH(2, IP, IF) = QF_GAM(2, IP, IF)
                        QF(1, IF) = QF(1, IF) + QF_GTH(1, IP, IF) * GTH(IP)
                        QF(2, IF) = QF(2, IF) + QF_GTH(2, IP, IF) * GTH(IP)
                    ENDDO
                ENDDO
                !
            ENDIF
        ENDDO

        !
        !---- Add contribution of all TE panels
        DO IEL = 1, NEL
            IF (LTPAN(IEL) .AND. (NETYPE(IEL)==0 .OR. NETYPE(IEL)==7))   &
                    & THEN
                !---- clear field-point velocities for each TE panel
                DO IF = IF1, IF2
                    DO KP = 1, 2
                        QF_GAMTE(1, KP, IF) = 0.
                        QF_GAMTE(2, KP, IF) = 0.
                        QF_SIGTE(1, KP, IF) = 0.
                        QF_SIGTE(2, KP, IF) = 0.
                    ENDDO
                ENDDO
                !
                CALL PANAIC(1, 2, XPT(1, IEL), YPT(1, IEL), NF, XF(KF), YF(KF), &
                        & IZERO, IZERO, IFTYPE(KF), GAMT(1, IEL), SIGT(1, IEL), &
                        & 2, QF(1, KF), QF_GAMTE(1, 1, KF), QF_SIGTE(1, 1, KF))
                !
                !HHY- Uncomment to save TE panel influences separately
                !---- Save TE velocity influences
                !c          DO IF = IF1, IF2
                !c            DO K = 1, 2
                !c              QF_GAMT(K,1,IEL,IF) = QF_GAMTE(K,1,IF)
                !c              QF_GAMT(K,2,IEL,IF) = QF_GAMTE(K,2,IF)
                !c              QF_SIGT(K,1,IEL,IF) = QF_SIGTE(K,1,IF)
                !c              QF_SIGT(K,2,IEL,IF) = QF_SIGTE(K,2,IF)
                !c            END DO
                !c          END DO
                !
                !---- Add TE influences to panel nodes (assumes TPANSET has been called)
                DO IF = IF1, IF2
                    DO KP = 1, 2
                        IF (KP==1) THEN
                            IP = IPTE1(IEL)
                        ELSE
                            IP = IPTE2(IEL)
                        ENDIF
                        IF (IP>0) THEN
                            DO K = 1, 2
                                QF_GAM(K, IP, IF) = QF_GAM(K, IP, IF)               &
                                        & + QF_GAMTE(K, 1, IF) * GAMT_GAM(1, KP, IEL)        &
                                        & + QF_GAMTE(K, 2, IF) * GAMT_GAM(2, KP, IEL)        &
                                        & + QF_SIGTE(K, 1, IF) * SIGT_GAM(1, KP, IEL)        &
                                        & + QF_SIGTE(K, 2, IF) * SIGT_GAM(2, KP, IEL)
                                QF_SIG(K, IP, IF) = QF_SIG(K, IP, IF)               &
                                        & + QF_GAMTE(K, 1, IF) * GAMT_SIG(1, KP, IEL)        &
                                        & + QF_GAMTE(K, 2, IF) * GAMT_SIG(2, KP, IEL)        &
                                        & + QF_SIGTE(K, 1, IF) * SIGT_SIG(1, KP, IEL)        &
                                        & + QF_SIGTE(K, 2, IF) * SIGT_SIG(2, KP, IEL)
                            ENDDO    ! K loop
                        ENDIF
                    ENDDO   ! KP loop
                ENDDO  ! IF loop
                !
            ENDIF
            !
            !---- Special treatment for vortex wake last GTH point to include TE panel
            IF (LTPAN(IEL) .AND. NETYPE(IEL)==7) THEN
                IP = IPLAST(IEL)
                DO IF = IF1, IF2
                    !---- Remove velocity w/o TE panel
                    !c            QF(1,IF) = QF(1,IF) - QF_GTH(1,IP,IF)*GTH(IP)
                    !c            QF(2,IF) = QF(2,IF) - QF_GTH(2,IP,IF)*GTH(IP)
                    !---- Save influence with TE panel
                    QF_GTH(1, IP, IF) = QF_GAM(1, IP, IF)
                    QF_GTH(2, IP, IF) = QF_GAM(2, IP, IF)
                    !---- Add back in velocity with TE panel
                    !c            QF(1,IF) = QF(1,IF) + QF_GTH(1,IP,IF)*GTH(IP)
                    !c            QF(2,IF) = QF(2,IF) + QF_GTH(2,IP,IF)*GTH(IP)
                ENDDO
            ENDIF

        ENDDO
        !
        !---- Add freestream
        !      DO IF = 1, NF
        !        QF(1,IF) = QF(1,IF) + QINF
        !      ENDDO
        !
    END SUBROUTINE QFCALC
end module m_qaic
