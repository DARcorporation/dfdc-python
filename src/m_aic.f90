module m_aic
    implicit none
contains
    !*==PANAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
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

    SUBROUTINE PANAIC(JFRST, JLAST, XP, YP, NF, XF, YF, JPAN1, JPAN2, IFTYPE, &
            & GAM, SIG, JDIM, QF, QF_GAM, QF_SIG)
        use m_lamp, only : lamp, lampc
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: JDIM, JFRST, JLAST, NF
        REAL, DIMENSION(*) :: GAM, SIG, XF, XP, YF, YP
        INTEGER, DIMENSION(*) :: IFTYPE, JPAN1, JPAN2
        REAL, DIMENSION(2, *) :: QF
        REAL, DIMENSION(2, JDIM, *) :: QF_GAM, QF_SIG
        !
        ! Local variables
        !
        REAL :: ANWX, ANWY, ASWX, ASWY, UG1, UG2, US1, US2, &
                & VG1, VG2, VS1, VS2
        INTEGER :: I, J, J1, J2
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Computes field-point velocities due to vortex/source panels.
        !     Computes the velocities' Jacobians w.r.t. the panel strengths.
        !
        !  Input:
        !     JFRST     first panel node
        !     JLAST     last panel node
        !     XP,YP(j)  panel node x,y locations
        !
        !     NF        number of field points
        !     XF,YF(i)  field-point x locations  i = 1..NF
        !     JPAN1(i)  j indices of the two nodes defining the panel
        !     JPAN2(i)      which contains field point i at its midpoint
        !     IFTYPE(i) dummy field point if IFTYPE < 0  (skip)
        !
        !     GAM(j)   sheet vortex strength at X,Y(j)  (+ clockwise)
        !     SIG(j)   sheet source strength at X,Y(j)
        !
        !     JDIM     Second dimension of QF_XXX arrays
        !
        !  Output:
        !     QF(1:2,i)     x,y-velocities at XF,YF(i)
        !     QF_GAM(.j.)   dQF(..)/dGAM(j)
        !     QF_SIG(.j.)   dQF(..)/dSIG(j)
        !----------------------------------------------------------------------
        !
        !---- unit vectors tangent and normal to wall, if any
        ASWX = 1.
        ASWY = 0.
        !
        ANWX = -ASWY
        ANWY = ASWX
        !
        DO I = 1, NF
            IF (IFTYPE(I)>=0) THEN
                !
                DO J = JFRST, JLAST - 1
                    J1 = J
                    J2 = J + 1
                    !
                    IF (J1==JPAN1(I) .AND. J2==JPAN2(I)) THEN
                        !--------- field point is at the center of the J..J+1 panel
                        CALL LAMPC(XP(J1), YP(J1), XP(J2), YP(J2), UG1, VG1, UG2, &
                                & VG2, US1, VS1, US2, VS2)
                    ELSE
                        !--------- field point is somewhere else
                        CALL LAMP(XP(J1), YP(J1), XP(J2), YP(J2), XF(I), YF(I), UG1, &
                                & VG1, UG2, VG2, US1, VS1, US2, VS2)
                    ENDIF
                    !
                    !-------- accumulate velocity components and their sensitivities
                    !         Note:  GAM is defined positive clockwise
                    !                PANE,LAMP assume positive counterclockwise
                    QF(1, I) = QF(1, I) - UG1 * GAM(J1) - UG2 * GAM(J2)            &
                            & + US1 * SIG(J1) + US2 * SIG(J2)
                    !
                    QF(2, I) = QF(2, I) - VG1 * GAM(J1) - VG2 * GAM(J2)            &
                            & + VS1 * SIG(J1) + VS2 * SIG(J2)
                    !
                    QF_GAM(1, J1, I) = QF_GAM(1, J1, I) - UG1
                    QF_GAM(1, J2, I) = QF_GAM(1, J2, I) - UG2
                    QF_SIG(1, J1, I) = QF_SIG(1, J1, I) + US1
                    QF_SIG(1, J2, I) = QF_SIG(1, J2, I) + US2
                    !
                    QF_GAM(2, J1, I) = QF_GAM(2, J1, I) - VG1
                    QF_GAM(2, J2, I) = QF_GAM(2, J2, I) - VG2
                    QF_SIG(2, J1, I) = QF_SIG(2, J1, I) + VS1
                    QF_SIG(2, J2, I) = QF_SIG(2, J2, I) + VS2
                ENDDO
            ENDIF
        ENDDO
        !
    END SUBROUTINE PANAIC
    !*==LINAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    SUBROUTINE LINAIC(JFRST, JLAST, XP, YP, NF, XF, YF, VOR, SOU, JDIM, QF, &
            & QF_VOR, QF_SOU)
        use m_lamp, only : ring
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: JDIM, JFRST, JLAST, NF
        REAL, DIMENSION(2, *) :: QF
        REAL, DIMENSION(2, JDIM, *) :: QF_SOU, QF_VOR
        REAL, DIMENSION(*) :: SOU, VOR, XF, XP, YF, YP
        !
        ! Local variables
        !
        REAL :: ANWX, ANWY, ASWX, ASWY, UG, US, VG, VS
        INTEGER :: I, J
        REAL, SAVE :: PI
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------
        !     Generates line-singularity AIC matrices.
        !
        !  Input:
        !     JFRST
        !     JLAST
        !     XP,YP(j)  Singularity x,y locations       j = JFRST..JLAST
        !
        !     NF        Number of field points
        !     XF,YF(i)  Field-point x,y locations  i = 1..NF
        !
        !     VOR      Vortex line (or ring) strength (+ clockwise)
        !     SOU      Source line (or ring) strength
        !
        !     JDIM     Second dimension of VIC array
        !
        !  Output:
        !     QF(1:2,i)     x,y-velocities at XF,YF(i)
        !     QF_VOR(.j.)   dQF(..)/dVOR(j)
        !     QF_SOU(.j.)   dQF(..)/dSOU(j)
        !--------------------------------------------------------
        DATA PI/3.1415926535897932384/
        !
        !---- unit vectors tangent and normal to wall, if any
        ASWX = 0.
        ASWY = 1.
        ANWX = 1.
        ANWY = 0.
        !
        DO I = 1, NF
            DO J = JFRST, JLAST
                CALL RING(XP(J), YP(J), XF(I), YF(I), UG, VG, US, VS)
                !
                !-------- accumulate velocity components and their sensitivities
                !-        Note:  GAM is defined positive clockwise
                !-               PANE,LAMP assume positive counterclockwise
                QF(1, I) = QF(1, I) - UG * VOR(J) + US * SOU(J)
                QF(2, I) = QF(2, I) - VG * VOR(J) + VS * SOU(J)
                !
                QF_VOR(1, J, I) = QF_VOR(1, J, I) - UG
                QF_SOU(1, J, I) = QF_SOU(1, J, I) + US
                !
                QF_VOR(2, J, I) = QF_VOR(2, J, I) - VG
                QF_SOU(2, J, I) = QF_SOU(2, J, I) + VS
            ENDDO
        ENDDO
        !
    END SUBROUTINE LINAIC
    !*==AXLAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! LINAIC



    SUBROUTINE AXLAIC(JFRST, JLAST, XP, YP, NF, XF, YF, JPAN1, JPAN2, IFTYPE, &
            & DBL, SRC, RCORE, JDIM, QF, QF_DBL, QF_SRC)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: JDIM, JFRST, JLAST, NF
        REAL, DIMENSION(*) :: DBL, RCORE, SRC, XF, XP, YF, YP
        INTEGER, DIMENSION(*) :: IFTYPE, JPAN1, JPAN2
        REAL, DIMENSION(2, *) :: QF
        REAL, DIMENSION(2, JDIM, *) :: QF_DBL, QF_SRC
        !
        ! Local variables
        !
        REAL :: DX, EPS, EPS1, EPS2, R1, R2, RCORE1, RCORE2, &
                & RMX1, RMX2, RSQ1, RSQ2, UG1, UG2, US1, US2, &
                & VG1, VG2, VS1, VS2, X1, X2, Y
        INTEGER :: I, J, J1, J2
        REAL, SAVE :: PI
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------
        !     Generates axis-line singularity AIC matrices.
        !
        !  Input:
        !     JFRST
        !     JLAST
        !     XP,YP(j)  Singularity x,r locations       j = JFRST, JLAST
        !
        !     NF        Number of field points
        !     XF,YF(i)  Field-point x,r locations  i = 1..NF
        !
        !     DBL(j)    Line-doublet strength (positive when flow is to -x on axis)
        !     SRC(j)    Line-source  strength
        !
        !     JDIM      Second dimension of VIC array
        !
        !  Output:
        !     QF(1:2,i)     x,r-velocities at XF,YF(i)
        !     QF_DBL(.j.)   dQF(..)/dDBL(j)
        !     QF_SRC(.j.)   dQF(..)/dSRC(j)
        !
        !     Special case:  if YF(i)=0, then QF(i) is actually q*r
        !--------------------------------------------------------
        DATA PI/3.1415926535897932384/
        !
        DO I = 1, NF
            IF (IFTYPE(I)>=0) THEN
                !
                DO J = JFRST, JLAST - 1
                    J1 = J
                    J2 = J + 1
                    !
                    DX = XP(J2) - XP(J1)
                    !
                    X1 = XF(I) - XP(J1)
                    X2 = XF(I) - XP(J2)
                    Y = YF(I)
                    !
                    RCORE1 = RCORE(J1)
                    RCORE2 = RCORE(J2)
                    !
                    !cc       IF(J1.EQ.JPAN1(I) .AND. J2.EQ.JPAN2(I)) THEN
                    !--------- field point is at the center of the J..J+1 line segment
                    RSQ1 = X1 * X1 + Y * Y + 0.25 * RCORE1**2
                    RSQ2 = X2 * X2 + Y * Y + 0.25 * RCORE2**2
                    !          ELSE
                    !C--------- field point is somewhere else
                    !           RSQ1 = X1*X1 + Y*Y
                    !           RSQ2 = X2*X2 + Y*Y
                    !          ENDIF
                    !
                    R1 = SQRT(RSQ1)
                    R2 = SQRT(RSQ2)
                    !
                    EPS1 = Y * Y / RSQ1
                    EPS2 = Y * Y / RSQ2
                    !
                    EPS = MIN(EPS1, EPS2)
                    !
                    !-------- compute r-x using exact expression,
                    !-          or 4th-order Taylor series about Y=0 to avoid roundoff with sqrt()
                    !          IF(EPS .GT. 0.02) THEN
                    RMX1 = R1 - X1
                    RMX2 = R2 - X2
                    !          ELSE
                    !           RSG1 = SIGN(R1,X1)
                    !           RSG2 = SIGN(R2,X2)
                    !           RMX1 = (R1 - RSG1) + RSG1*EPS1*(0.5 + 0.125*EPS1)
                    !           RMX2 = (R2 - RSG2) + RSG2*EPS2*(0.5 + 0.125*EPS2)
                    !          ENDIF
                    !
                    !          IF(MIN(EPS1,EPS2) .LT. 0.0001
                    !     &       .AND. RSG1.GT.0.0
                    !     &       .AND. RSG2.GT.0.0) THEN
                    !           RMX1 = R2
                    !           RMX2 = R1
                    !          ENDIF
                    !C
                    US1 = -(1.0 / R1 + LOG(RMX1 / RMX2) / DX) / (4.0 * PI)
                    US2 = (1.0 / R2 + LOG(RMX1 / RMX2) / DX) / (4.0 * PI)
                    VS1 = Y * ((1.0 / R1 - 1.0 / DX) / RMX1 + (1.0 / DX) / RMX2) / (4.0 * PI)
                    VS2 = Y * ((1.0 / DX) / RMX1 - (1.0 / R2 + 1.0 / DX) / RMX2) / (4.0 * PI)
                    !
                    !-------- no doublet line currently implemented
                    UG1 = 0.
                    UG2 = 0.
                    VG1 = 0.
                    VG2 = 0.
                    !
                    !-------- accumulate velocity components and their sensitivities
                    QF(1, I) = QF(1, I) - UG1 * DBL(J1) - UG2 * DBL(J2)            &
                            & + US1 * SRC(J1) + US2 * SRC(J2)
                    !
                    QF(2, I) = QF(2, I) - VG1 * DBL(J1) - VG2 * DBL(J2)            &
                            & + VS1 * SRC(J1) + VS2 * SRC(J2)
                    !
                    QF_DBL(1, J1, I) = QF_DBL(1, J1, I) - UG1
                    QF_DBL(1, J2, I) = QF_DBL(1, J2, I) - UG2
                    QF_SRC(1, J1, I) = QF_SRC(1, J1, I) + US1
                    QF_SRC(1, J2, I) = QF_SRC(1, J2, I) + US2
                    !
                    QF_DBL(2, J1, I) = QF_DBL(2, J1, I) - VG1
                    QF_DBL(2, J2, I) = QF_DBL(2, J2, I) - VG2
                    QF_SRC(2, J1, I) = QF_SRC(2, J1, I) + VS1
                    QF_SRC(2, J2, I) = QF_SRC(2, J2, I) + VS2
                ENDDO
            ENDIF
        ENDDO
        !
    END SUBROUTINE AXLAIC
    !*==PNTAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! AXLAIC



    SUBROUTINE PNTAIC(JFRST, JLAST, XP, YP, NF, XF, YF, DBL, SRC, JDIM, QF, &
            & QF_DBL, QF_SRC)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: JDIM, JFRST, JLAST, NF
        REAL, DIMENSION(*) :: DBL, SRC, XF, XP, YF, YP
        REAL, DIMENSION(2, *) :: QF
        REAL, DIMENSION(2, JDIM, *) :: QF_DBL, QF_SRC
        !
        ! Local variables
        !
        REAL :: DX, DY, R32, RSQ, UG, US, VG, VS
        INTEGER :: I, J
        REAL, SAVE :: PI
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------
        !     Generates point-source and point-doublet AIC matrices
        !     (axisymmetric case only)
        !
        !  Input:
        !     JFRST
        !     JLAST
        !     XP,YP(j)  Singularity x,r locations       j = JFRST..JLAST
        !
        !     NF        Number of field points
        !     XF,YF(i)  Field-point x,r locations  i = 1..NF
        !
        !     DBL       Point-doublet strength (positive when flow is to -x on axis)
        !     SRC       Point-source  strength
        !
        !     JDIM      Second dimension of VIC array
        !
        !  Output:
        !     QF(1:2,i)     x,r-velocities at XF,YF(i)
        !     QF_DBL(.j.)   dQF(..)/dDBL(j)
        !     QF_SRC(.j.)   dQF(..)/dSRC(j)
        !--------------------------------------------------------
        DATA PI/3.1415926535897932384/
        !
        DO I = 1, NF
            DO J = JFRST, JLAST
                !-------- point singularity on axis... DBL is x-doublet, SRC is point source
                DX = XF(I) - XP(J)
                DY = YF(I)
                RSQ = DX * DX + DY * DY
                R32 = SQRT(RSQ)**3
                UG = (2.0 * DX * DX - DY * DY) / (4.0 * PI * R32 * RSQ)
                VG = (3.0 * DX * DY) / (4.0 * PI * R32 * RSQ)
                US = DX / (4.0 * PI * R32)
                VS = DY / (4.0 * PI * R32)
                !
                QF(1, I) = QF(1, I) - UG * DBL(J) + US * SRC(J)
                QF(2, I) = QF(2, I) - VG * DBL(J) + VS * SRC(J)
                !
                QF_DBL(1, J, I) = QF_DBL(1, J, I) - UG
                QF_DBL(2, J, I) = QF_DBL(2, J, I) - VG
                QF_SRC(1, J, I) = QF_SRC(1, J, I) + US
                QF_SRC(2, J, I) = QF_SRC(2, J, I) + VS
            ENDDO
        ENDDO
        !
    END SUBROUTINE PNTAIC
    !*==PANAICD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PNTAIC



    SUBROUTINE PANAICD(JFRST, JLAST, XP, YP, NF, XF, YF, JPAN1, JPAN2, IFTYPE, &
            & GAM, SIG, JDIM, QF, QF_GAM, QF_SIG, QF_XP, QF_YP, &
            & QF_XF, QF_YF)
        use m_lamp, only : dlampc, dlamp
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: JDIM, JFRST, JLAST, NF
        REAL, DIMENSION(*) :: GAM, SIG, XF, XP, YF, YP
        INTEGER, DIMENSION(*) :: IFTYPE, JPAN1, JPAN2
        REAL, DIMENSION(2, *) :: QF, QF_XF, QF_YF
        REAL, DIMENSION(2, JDIM, *) :: QF_GAM, QF_SIG, QF_XP, QF_YP
        !
        ! Local variables
        !
        REAL :: ANWX, ANWY, ASWX, ASWY
        INTEGER :: I, J, J1, J2, KQ, KR
        REAL, DIMENSION(2) :: QG1, QG2, QS1, QS2, R1, R2, RF
        REAL, DIMENSION(2, 2) :: QG1_R1, QG1_R2, QG1_RF, QG2_R1, &
                & QG2_R2, QG2_RF, QS1_R1, QS1_R2, &
                & QS1_RF, QS2_R1, QS2_R2, QS2_RF
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Sames as PANAIC, but also returns Jacobians w.r.t. the geometry
        !----------------------------------------------------------------------
        !
        !---- local arrays for calling DLAMP
        !
        !---- unit vectors tangent and normal to wall, if any
        ASWX = 1.
        ASWY = 0.
        !
        ANWX = -ASWY
        ANWY = ASWX
        !
        DO I = 1, NF
            IF (IFTYPE(I)>=0) THEN
                !
                RF(1) = XF(I)
                RF(2) = YF(I)
                !
                DO J = JFRST, JLAST - 1
                    J1 = J
                    J2 = J + 1
                    !
                    IF (J1==JPAN1(I) .AND. J2==JPAN2(I)) THEN
                        !--------- field point is at the center of the J..J+1 panel
                        R1(1) = XP(J1)
                        R1(2) = YP(J1)
                        R2(1) = XP(J2)
                        R2(2) = YP(J2)
                        CALL DLAMPC(R1, R2, QG1, QG1_R1, QG1_R2, QG2, QG2_R1, QG2_R2, &
                                & QS1, QS1_R1, QS1_R2, QS2, QS2_R1, QS2_R2)
                        DO KQ = 1, 2
                            DO KR = 1, 2
                                QG1_RF(KQ, KR) = 0.
                                QG2_RF(KQ, KR) = 0.
                                QS1_RF(KQ, KR) = 0.
                                QS2_RF(KQ, KR) = 0.
                            ENDDO
                        ENDDO
                    ELSE
                        !--------- field point is somewhere else
                        R1(1) = XP(J1)
                        R1(2) = YP(J1)
                        R2(1) = XP(J2)
                        R2(2) = YP(J2)
                        CALL DLAMP(R1, R2, RF, QG1, QG1_R1, QG1_R2, QG1_RF, QG2, &
                                & QG2_R1, QG2_R2, QG2_RF, QS1, QS1_R1, QS1_R2, &
                                & QS1_RF, QS2, QS2_R1, QS2_R2, QS2_RF)
                    ENDIF
                    !
                    !
                    !-------- accumulate velocity components and their sensitivities
                    !-        Note:  GAM is defined positive clockwise
                    !-               PANE,LAMP assume positive counterclockwise
                    DO KQ = 1, 2
                        QF(KQ, I) = QF(KQ, I) - QG1(KQ) * GAM(J1) - QG2(KQ)       &
                                & * GAM(J2) + QS1(KQ) * SIG(J1) + QS2(KQ)       &
                                & * SIG(J2)
                        !
                        QF_GAM(KQ, J1, I) = QF_GAM(KQ, J1, I) - QG1(KQ)
                        QF_GAM(KQ, J2, I) = QF_GAM(KQ, J2, I) - QG2(KQ)
                        QF_SIG(KQ, J1, I) = QF_SIG(KQ, J1, I) + QS1(KQ)
                        QF_SIG(KQ, J2, I) = QF_SIG(KQ, J2, I) + QS2(KQ)
                        !
                        QF_XP(KQ, J1, I) = QF_XP(KQ, J1, I) - QG1_R1(KQ, 1) * GAM(J1)&
                                & - QG2_R1(KQ, 1) * GAM(J2) + QS1_R1(KQ, 1)&
                                & * SIG(J1) + QS2_R1(KQ, 1) * SIG(J2)
                        QF_XP(KQ, J2, I) = QF_XP(KQ, J2, I) - QG1_R2(KQ, 1) * GAM(J1)&
                                & - QG2_R2(KQ, 1) * GAM(J2) + QS1_R2(KQ, 1)&
                                & * SIG(J1) + QS2_R2(KQ, 1) * SIG(J2)
                        QF_XF(KQ, I) = QF_XF(KQ, I) - QG1_RF(KQ, 1) * GAM(J1)      &
                                & - QG2_RF(KQ, 1) * GAM(J2) + QS1_RF(KQ, 1)   &
                                & * SIG(J1) + QS2_RF(KQ, 1) * SIG(J2)
                        !
                        QF_YP(KQ, J1, I) = QF_YP(KQ, J1, I) - QG1_R1(KQ, 2) * GAM(J1)&
                                & - QG2_R1(KQ, 2) * GAM(J2) + QS1_R1(KQ, 2)&
                                & * SIG(J1) + QS2_R1(KQ, 2) * SIG(J2)
                        QF_YP(KQ, J2, I) = QF_YP(KQ, J2, I) - QG1_R2(KQ, 2) * GAM(J1)&
                                & - QG2_R2(KQ, 2) * GAM(J2) + QS1_R2(KQ, 2)&
                                & * SIG(J1) + QS2_R2(KQ, 2) * SIG(J2)
                        QF_YF(KQ, I) = QF_YF(KQ, I) - QG1_RF(KQ, 2) * GAM(J1)      &
                                & - QG2_RF(KQ, 2) * GAM(J2) + QS1_RF(KQ, 2)   &
                                & * SIG(J1) + QS2_RF(KQ, 2) * SIG(J2)
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
        !
    END SUBROUTINE PANAICD
    !*==LINAICD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PANAICD



    SUBROUTINE LINAICD(JFRST, JLAST, XP, YP, NF, XF, YF, VOR, SOU, JDIM, QF, &
            & QF_VOR, QF_SOU, QF_XP, QF_YP, QF_XF, QF_YF)
        use m_lamp, only : dring
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: JDIM, JFRST, JLAST, NF
        REAL, DIMENSION(2, *) :: QF, QF_XF, QF_YF
        REAL, DIMENSION(2, JDIM, *) :: QF_SOU, QF_VOR, QF_XP, QF_YP
        REAL, DIMENSION(*) :: SOU, VOR, XF, XP, YF, YP
        !
        ! Local variables
        !
        REAL :: ANWX, ANWY, ASWX, ASWY, UG, UG_XF, UG_XP, UG_YF, &
                & UG_YP, US, US_XF, US_XP, US_YF, US_YP, VG, &
                & VG_XF, VG_XP, VG_YF, VG_YP, VS, VS_XF, VS_XP, &
                & VS_YF, VS_YP
        INTEGER :: I, J
        REAL, SAVE :: PI
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Sames as LINAIC, but also returns Jacobians w.r.t. the geometry
        !----------------------------------------------------------------------
        DATA PI/3.1415926535897932384/
        !
        !---- unit vectors tangent and normal to wall, if any
        ASWX = 0.
        ASWY = 1.
        ANWX = 1.
        ANWY = 0.
        !
        DO I = 1, NF
            DO J = JFRST, JLAST
                CALL DRING(XP(J), YP(J), XF(I), YF(I), UG, UG_XP, UG_YP, UG_XF, &
                        & UG_YF, VG, VG_XP, VG_YP, VG_XF, VG_YF, US, US_XP, US_YP, &
                        & US_XF, US_YF, VS, VS_XP, VS_YP, VS_XF, VS_YF)
                !
                !-------- accumulate velocity components and their sensitivities
                !-        Note:  GAM is defined positive clockwise
                !-               PANE,LAMP assume positive counterclockwise
                QF(1, I) = QF(1, I) - UG * VOR(J) + US * SOU(J)
                QF(2, I) = QF(2, I) - VG * VOR(J) + VS * SOU(J)
                !
                QF_VOR(1, I, J) = QF_VOR(1, I, J) - UG
                QF_SOU(1, I, J) = QF_SOU(1, I, J) + US
                !
                QF_VOR(2, I, J) = QF_VOR(2, I, J) - VG
                QF_SOU(2, I, J) = QF_SOU(2, I, J) + VS
                !
                QF_XP(1, I, J) = QF_XP(1, I, J) - UG_XP * VOR(J) + US_XP * SOU(J)
                QF_YP(1, I, J) = QF_YP(1, I, J) - UG_YP * VOR(J) + US_YP * SOU(J)
                QF_XF(1, I) = QF_XF(1, I) - UG_XF * VOR(J) + US_XF * SOU(J)
                QF_YF(1, I) = QF_YF(1, I) - UG_YF * VOR(J) + US_YF * SOU(J)
                !
                QF_XP(2, I, J) = QF_XP(2, I, J) - VG_XP * VOR(J) + VS_XP * SOU(J)
                QF_YP(2, I, J) = QF_YP(2, I, J) - VG_YP * VOR(J) + VS_YP * SOU(J)
                QF_XF(2, I) = QF_XF(2, I) - VG_XF * VOR(J) + VS_XF * SOU(J)
                QF_YF(2, I) = QF_YF(2, I) - VG_YF * VOR(J) + VS_YF * SOU(J)
                !
            ENDDO
        ENDDO
        !
    END SUBROUTINE LINAICD
    !*==PNTAICD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! LINAICD



    SUBROUTINE PNTAICD(JFRST, JLAST, XP, YP, NF, XF, YF, DBL, SRC, JDIM, QF, &
            & QF_DBL, QF_SRC, QF_XP, QF_YP, QF_XF, QF_YF)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: JDIM, JFRST, JLAST, NF
        REAL, DIMENSION(*) :: DBL, SRC, XF, XP, YF, YP
        REAL, DIMENSION(2, *) :: QF, QF_XF, QF_YF
        REAL, DIMENSION(2, JDIM, *) :: QF_DBL, QF_SRC, QF_XP, QF_YP
        !
        ! Local variables
        !
        REAL :: DX, DY, R32, R52, RSQ, UG, UG_XF, UG_XP, UG_YF, &
                & UG_YP, US, US_XF, US_XP, US_YF, US_YP, VG, &
                & VG_XF, VG_XP, VG_YF, VG_YP, VS, VS_XF, VS_XP, &
                & VS_YF, VS_YP
        INTEGER :: I, J
        REAL, SAVE :: PI
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Sames as PNTAIC, but also returns Jacobians w.r.t. the geometry
        !----------------------------------------------------------------------
        DATA PI/3.1415926535897932384/
        !
        DO I = 1, NF
            DO J = JFRST, JLAST
                !-------- point singularity on axis... DBL is x-doublet, SRC is point source
                DX = XF(I) - XP(J)
                DY = YF(I)
                !
                RSQ = DX * DX + DY * DY
                R32 = SQRT(RSQ)**3
                R52 = R32 * RSQ
                !
                UG = (2.0 * DX * DX - DY * DY) / (4.0 * PI * R52)
                VG = (3.0 * DX * DY) / (4.0 * PI * R52)
                US = DX / (4.0 * PI * R32)
                VS = DY / (4.0 * PI * R32)
                !
                UG_XF = 4.0 * DX / (4.0 * PI * R52) - 5.0 * DX * UG / RSQ
                UG_YF = -2.0 * DY / (4.0 * PI * R52) - 5.0 * DY * UG / RSQ
                VG_XF = 3.0 * DY / (4.0 * PI * R52) - 5.0 * DX * VG / RSQ
                VG_YF = 3.0 * DX / (4.0 * PI * R52) - 5.0 * DY * VG / RSQ
                !
                US_XF = 1.0 / (4.0 * PI * R32) - 3.0 * DX * US / RSQ
                US_YF = -3.0 * DY * US / RSQ
                VS_XF = -3.0 * DX * VS / RSQ
                VS_YF = 1.0 / (4.0 * PI * R32) - 3.0 * DY * VS / RSQ
                !
                UG_XP = -UG_XF
                UG_YP = 0.
                VG_XP = -VG_XF
                VG_YP = 0.
                !
                US_XP = -US_XF
                US_YP = 0.
                VS_XP = -VS_XF
                VS_YP = 0.
                !
                !
                QF(1, I) = QF(1, I) - UG * DBL(J) + US * SRC(J)
                QF(2, I) = QF(2, I) - VG * DBL(J) + VS * SRC(J)
                !
                QF_DBL(1, I, J) = QF_DBL(1, I, J) - UG
                QF_DBL(2, I, J) = QF_DBL(2, I, J) - VG
                QF_SRC(1, I, J) = QF_SRC(1, I, J) + US
                QF_SRC(2, I, J) = QF_SRC(2, I, J) + VS
                !
                QF_XP(1, I, J) = QF_XP(1, I, J) - UG_XP * DBL(J) + US_XP * SRC(J)
                QF_YP(1, I, J) = QF_YP(1, I, J) - UG_YP * DBL(J) + US_YP * SRC(J)
                QF_XF(1, I) = QF_XF(1, I) - UG_XF * DBL(J) + US_XF * SRC(J)
                QF_YF(1, I) = QF_YF(1, I) - UG_YF * DBL(J) + US_YF * SRC(J)
                !
                QF_XP(2, I, J) = QF_XP(2, I, J) - VG_XP * DBL(J) + VS_XP * SRC(J)
                QF_YP(2, I, J) = QF_YP(2, I, J) - VG_YP * DBL(J) + VS_YP * SRC(J)
                QF_XF(2, I) = QF_XF(2, I) - VG_XF * DBL(J) + VS_XF * SRC(J)
                QF_YF(2, I) = QF_YF(2, I) - VG_YF * DBL(J) + VS_YF * SRC(J)
                !
            ENDDO
        ENDDO
        !
    END SUBROUTINE PNTAICD
end module m_aic
