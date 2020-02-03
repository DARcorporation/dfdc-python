module m_sgutil
    implicit none
contains
    !*==SGCURV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020

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

    SUBROUTINE SGCURV(LQUERY, LCPLOT, NB, XB, XPB, YB, YPB, SB, NGX, NG, SG, &
            & NPTS, CVEX, SMOF, FSL, FSR, NREF, SREF1, SREF2, CRRAT)
        use m_spline, only : trisol, segspl, seval, curv
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        INTEGER, PARAMETER :: IX = 1601, NREFX = 20
        !
        ! Dummy arguments
        !
        REAL :: CVEX, FSL, FSR, SMOF
        LOGICAL :: LCPLOT, LQUERY
        INTEGER :: NB, NG, NGX, NPTS, NREF
        REAL, DIMENSION(NREF) :: CRRAT, SREF1, SREF2
        REAL, DIMENSION(NB) :: SB, XB, XPB, YB, YPB
        REAL, DIMENSION(NGX) :: SG
        !
        ! Local variables
        !
        REAL, DIMENSION(IX) :: AA, BB, CC, CV, CV0, CVINT, CVS, &
                & S, SNEW
        REAL :: CURVK, CURVL, CVASQ, CVAVG, CVLE, CVM, CVMAX, &
                & CVO, CVTE, DLE, DSA, DSAVG, DSB, DSLE, DSM, &
                & DSMAX, DSMIN, DSP, DSTE, DTE, FRAC, GCV, SBTOT, &
                & SGI, SM, SO, STOT
        REAL, DIMENSION(NREFX) :: CVR
        INTEGER :: I, IB, IG, IPASS, IREF, K, KK, KMULT, N
        INTEGER, SAVE :: NPASS
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------------
        !     Sets spacing array SG based on curvature and a few parameters
        !
        !  Input:  LQUERY   if T, requests user input for parameters before paneling
        !          LCPLOT   if T, spacing parameters are plotted
        !          NB       number of coordinates
        !          XB(i)    x coordinates
        !          XPB(i)   dXB/dSB spline derivative
        !          YB(i)    y coordinates
        !          YPB(i)   dYB/dSB spline derivative
        !          SB(i)    spline parameter (typically arc length)
        !          NGX      max number of points which can be generated
        !
        !  Input, also Output if interactively modified:
        !          NG       number of points actually generated
        !          SG(k)    fractional spacing array  0.0 .. 1.0
        !          NPTS     number of points to be generated
        !          CVEX     curvature-exponent parameter
        !          SMOF     smoothing-length factor
        !          FSL      fractional spacing at left  endpoint
        !          FSR      fractional spacing at right endpoint
        !          NREF     number of local-refinement regions
        !          SREF1(r) first point of refinement region
        !          SREF2(r) last  point of refinement region
        !          CRRAT(r) fractional spacing over refinement region
        !
        !
        !---------------------------------------------------------------------
        !
        !
        !
        !
        !
        !---- max number of point-setting passes
        DATA NPASS/1/
        !
        IF (NB>IX) STOP 'SGCURV: Local array overflow on IX'
        IF (NREF>NREFX) STOP 'SGCURV: Local array overflow on NREFX'
        !
        !
        SBTOT = SB(NB) - SB(1)
        IF (SBTOT<=0.0) THEN
            WRITE (*, *) '? SGCURV: Zero element arc length'
            RETURN
        ENDIF
        !
        !---- if no interaction or no spacing is present, go compute it first
        IF (.NOT.LQUERY .OR. NG==0) THEN
        ENDIF
        !================================================================
        !==== Set surface panel node distribution based on curvature
        !
        NG = NPTS
        !
        IF (NG<2) THEN
            WRITE (*, *) '? SGCURV: Must specify at least two panel nodes'
            RETURN
        ENDIF
        !
        KMULT = MIN(8, (IX - 1) / NB)
        !
        !---- set up working coordinate array
        I = 1
        IB = 1
        S(I) = SB(IB)
        DO IB = 1, NB - 1
            DSB = SB(IB + 1) - SB(IB)
            !
            IF (DSB>0.0) THEN
                KK = KMULT
            ELSE
                KK = 1
            ENDIF
            !
            DO K = 1, KK
                FRAC = FLOAT(K) / FLOAT(KK)
                I = I + 1
                IF (I>IX) STOP 'SGCURV:  Array overflow on IX'
                S(I) = SB(IB + 1) * FRAC + SB(IB) * (1.0 - FRAC)
            ENDDO
        ENDDO
        !
        !---- number of points in working array
        N = I
        !
        STOT = S(N) - S(1)
        DO IPASS = 1, NPASS
            !
            DSAVG = STOT / FLOAT(NG - 1)
            !
            DSLE = FSL * DSAVG
            DSTE = FSR * DSAVG
            !
            !---- set up curvature array (nondimensionalized with perimeter)
            DO I = 1, N
                CV(I) = CURV(S(I), XB, XPB, YB, YPB, SB, NB) * STOT
                CV(I) = ABS(CV(I)) + 1.0
            ENDDO
            !
            !---- reset curvature at corner nodes from adjacent nodes
            DO I = 2, N - 2
                IF (S(I)==S(I + 1)) THEN
                    CV(I) = 0.5 * (CV(I - 1) + CV(I + 2))
                    CV(I + 1) = CV(I)
                ENDIF
            ENDDO
            !
            !---- raise curvature to power sqrt(CVEX), limiting to prevent overflows
            DO I = 1, N
                GCV = MIN(0.5 * CVEX * LOG(CV(I)), 12.0)
                CV(I) = EXP(GCV)
            ENDDO
            !
            !---- set max and average curvature
            CVMAX = 0.01
            CVAVG = 0.
            DO I = 1, N - 1
                CVASQ = 0.5 * (CV(I)**2 + CV(I + 1)**2)
                CVMAX = MAX(CVMAX, SQRT(CVASQ))
                CVAVG = CVAVG + CVASQ * (S(I + 1) - S(I))
            ENDDO
            CVAVG = SQRT(CVAVG / (S(N) - S(1)))
            !
            CVAVG = MAX(CVAVG, 1.0)
            CVMAX = MAX(CVMAX, 1.0)
            !
            !---- set artificial curvature at ends to get approximate ds/c there
            CVLE = CVAVG * DSAVG / DSLE
            CVTE = CVAVG * DSAVG / DSTE
            CV(1) = CVLE
            CV(N) = CVTE
            !
            DO I = 1, N
                CV0(I) = CV(I)
            ENDDO
            !
            !---- set curvature smoothing length
            CURVL = 0.008 * STOT
            !
            !---- set up implicit system for smoothed curvature array CV
            CURVK = CURVL**2 * SMOF
            AA(1) = 1.0
            CC(1) = 0.
            DO I = 2, N - 1
                DSM = S(I) - S(I - 1)
                DSP = S(I + 1) - S(I)
                DSA = 0.5 * (DSM + DSP)
                IF (DSM==0.0 .OR. DSP==0.0) THEN
                    BB(I) = 0.0
                    AA(I) = 1.0
                    CC(I) = 0.0
                ELSE
                    BB(I) = -CURVK / (DSM * DSA)
                    AA(I) = CURVK / (DSM * DSA) + CURVK / (DSP * DSA) + 1.0
                    CC(I) = -CURVK / (DSP * DSA)
                ENDIF
            ENDDO
            AA(N) = 1.0
            BB(N) = 0.
            !
            !---- set artificial curvature at the bunching points
            DO IREF = 1, NREF
                IF (CRRAT(IREF)>0.0) THEN
                    CVR(IREF) = CVAVG / MAX(CRRAT(IREF), 0.0001)
                ELSE
                    CVR(IREF) = 0.0
                ENDIF
                !
                DO I = 2, N - 1
                    IF (CRRAT(IREF)>0.0) THEN
                        !--------- set modified curvature if point is in refinement area
                        SGI = (S(I) - S(1)) / STOT
                        !c            BB(I) = 0.
                        !c            AA(I) = 1.0
                        !c            CC(I) = 0.
                        IF (SGI>SREF1(IREF) .AND. SGI<SREF2(IREF)) CV(I)    &
                                & = CVR(IREF)
                    ENDIF
                ENDDO
            ENDDO
            !
            !---- calculate smoothed curvature array
            CALL TRISOL(AA, BB, CC, CV, N)
            !
            !---- spline curvature array
            CALL SEGSPL(CV, CVS, S, N)
            !
            !---- Integrate exponentiated curvature with arc length
            CVINT(1) = 0.
            DO I = 2, N
                SM = S(I - 1)
                SO = S(I)
                CVM = SEVAL(SM, CV, CVS, S, N)
                CVO = SEVAL(SO, CV, CVS, S, N)
                CVM = MAX(CVM, 1.0)
                CVO = MAX(CVO, 1.0)
                !
                CVINT(I) = CVINT(I - 1) + 0.5 * (CVO + CVM) * (SO - SM)
            ENDDO
            !
            DO I = 1, N
                CVINT(I) = CVINT(I) + S(I)
            ENDDO
            !
            !---- Calculate normalized surface spacing distribution arrays
            CALL SETCRV(SNEW, NG, S, CVINT, N)
            !
            DO IG = 1, NG
                SG(IG) = (SNEW(IG) - SNEW(1)) / (SNEW(NG) - SNEW(1))
            ENDDO
            !
            !---- calculate actual obtained grid spacings
            DLE = ABS(SNEW(2) - SNEW(1))
            DTE = ABS(SNEW(NG) - SNEW(NG - 1))
            DSMIN = DLE
            DSMAX = 0.0
            DO IG = 2, NG
                DSMIN = MIN(DSMIN, ABS(SNEW(IG) - SNEW(IG - 1)))
                DSMAX = MAX(DSMAX, ABS(SNEW(IG) - SNEW(IG - 1)))
            ENDDO
            !
        ENDDO
    END SUBROUTINE SGCURV
    !*==SETCRV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SGCURV



    SUBROUTINE SETCRV(SNEW, NPTS, SB, HH, N)
        use m_spline, only : seval, splina
        !------------------------------------------------------------------
        !     Sets spacing array SNEW(i) based on surface curvature.
        !
        !     Input: NPTS    number of points defining surface
        !            SB(.)   surface arc length
        !            HH(.)   integrated spacing function
        !            N       number of spacing array points to be generated
        !
        !     Output: SNEW(.)  spacing array to be generated
        !             LFUDGL   T if curvature was augmented al left endpoint
        !             LFUDGT   T if curvature was augmented al left endpoint
        !------------------------------------------------------------------
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        INTEGER, PARAMETER :: NMAX = 1601
        !
        ! Dummy arguments
        !
        INTEGER :: N, NPTS
        REAL, DIMENSION(N) :: HH, SB
        REAL, DIMENSION(NPTS) :: SNEW
        !
        ! Local variables
        !
        REAL :: DH, HHI, RN
        INTEGER :: I
        REAL, DIMENSION(NMAX) :: SBP
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !
        IF (N>NMAX) STOP 'SETCRV: array overflow'
        !
        RN = FLOAT(NPTS - 1)
        !
        !
        !---- spline  HH(i)  =  h(s)  =  s + c(s)
        CALL SPLINA(HH, SBP, SB, N)
        !
        !---- invert derivative from dh/ds to ds/dh
        DO I = 1, N
            SBP(I) = 1.0 / SBP(I)
        ENDDO
        !
        !---- now,  SB(i),SBP(i)  =  s(h), s'(h)
        !
        !
        DH = HH(N) / RN
        !
        SNEW(1) = 0.0
        DO I = 2, NPTS - 1
            HHI = FLOAT(I - 1) * DH
            SNEW(I) = SEVAL(HHI, SB, SBP, HH, N)
        ENDDO
        SNEW(NPTS) = SB(N)
        !
    END SUBROUTINE SETCRV
end module m_sgutil
