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
!    1-D Spline Package.
!    Interpolates a function x(s) over an s interval
!
!                                              Mark Drela
!                                              1985
!
!    Usage:
!
!C---- fill S(i), X(i) arrays
!      S(i) = ...
!      X(i) = ...
!
!C---- calculate spline coefficients for x(s) from discrete values
!C-    (can use SPLIND or SPLINA instead)
!      CALL SPLINE(X,XS,S,N)
!
!C---- evaluate splined x(s) and/or its derivatives
!C-    at any number of s points SS
!      XX   = SEVAL(SS,X,XS,S,N)
!      XXS  = DEVAL(SS,X,XS,S,N)
!      XXSS = D2VAL(SS,X,XS,S,N)
!
!C---- alternative to calling SEVAL,DEVAL,D2VAL separately
!C-    (slightly more efficient if all three quantities are needed)
!      CALL SEVALL(SS,X,XS,S,N, XX,XXS,XXSS)
!
!


SUBROUTINE SPLINE(X, XS, S, N)
    DIMENSION X(N), XS(N), S(N)
    PARAMETER (NMAX = 1601)
    DIMENSION A(NMAX), B(NMAX), C(NMAX)
    !-------------------------------------------------------
    !     Calculates spline coefficients for X(S).          |
    !     Natural end conditions are used (zero 3rd         |
    !      derivative over first, last intervals).          |
    !                                                       |
    !     To evaluate the spline at some value of S,        |
    !     use SEVAL and/or DEVAL.                           |
    !                                                       |
    !     S        independent variable array (input)       |
    !     X        dependent variable array   (input)       |
    !     XS       dX/dS array                (calculated)  |
    !     N        number of points           (input)       |
    !                                                       |
    !-------------------------------------------------------
    IF(N.GT.NMAX) STOP 'SPLINE: array overflow, increase NMAX'
    !
    IF(N.EQ.1) THEN
        XS(1) = 0.
        RETURN
    ENDIF
    !
    DO 1 I = 2, N - 1
        DSM = S(I) - S(I - 1)
        DSP = S(I + 1) - S(I)
        B(I) = DSP
        A(I) = 2.0 * (DSM + DSP)
        C(I) = DSM
        XS(I) = 3.0 * ((X(I + 1) - X(I)) * DSM / DSP + (X(I) - X(I - 1)) * DSP / DSM)
    1 CONTINUE
    !
    !---- set zero 3rd derivative end conditions
    A(1) = 1.0
    C(1) = 1.0
    XS(1) = 2.0 * (X(2) - X(1)) / (S(2) - S(1))
    !
    B(N) = 1.0
    A(N) = 1.0
    XS(N) = 2.0 * (X(N) - X(N - 1)) / (S(N) - S(N - 1))
    !
    IF(N.EQ.2) THEN
        !----- if only two points are present, specify zero 2nd derivative instead
        !-     (straight line interpolation will result)
        B(N) = 1.0
        A(N) = 2.0
        XS(N) = 3.0 * (X(N) - X(N - 1)) / (S(N) - S(N - 1))
    ENDIF
    !
    !---- solve for derivative array XS
    CALL TRISOL(A, B, C, XS, N)
    !
    RETURN
END
! SPLINE



SUBROUTINE SPLIND(X, XS, S, N, XS1, XS2)
    DIMENSION X(N), XS(N), S(N)
    PARAMETER (NMAX = 1601)
    DIMENSION A(NMAX), B(NMAX), C(NMAX)
    !-------------------------------------------------------
    !     Calculates spline coefficients for X(S).          |
    !     Same as SPLINE, but also allows specified-slope   |
    !     or zero-curvature end conditions to be imposed.   |
    !                                                       |
    !     To evaluate the spline at some value of S,        |
    !     use SEVAL and/or DEVAL.                           |
    !                                                       |
    !     S        independent variable array (input)       |
    !     X        dependent variable array   (input)       |
    !     XS       dX/dS array                (calculated)  |
    !     N        number of points           (input)       |
    !     XS1,XS2  endpoint derivatives       (input)       |
    !              If =  999.0, then usual zero second      |
    !              derivative end condition(s) are used     |
    !              If = -999.0, then zero third             |
    !              derivative end condition(s) are used     |
    !                                                       |
    !     Note: specifying both XS1,XS2 = -999.0            |
    !           is equivalent to using SPLINE.              |
    !                                                       |
    !-------------------------------------------------------
    IF(N.GT.NMAX) STOP 'SPLIND: array overflow, increase NMAX'
    !
    IF(N.EQ.1) THEN
        XS(1) = 0.
        RETURN
    ENDIF
    !
    DO 1 I = 2, N - 1
        DSM = S(I) - S(I - 1)
        DSP = S(I + 1) - S(I)
        B(I) = DSP
        A(I) = 2.0 * (DSM + DSP)
        C(I) = DSM
        XS(I) = 3.0 * ((X(I + 1) - X(I)) * DSM / DSP + (X(I) - X(I - 1)) * DSP / DSM)
    1 CONTINUE
    !
    IF(XS1.EQ.999.0) THEN
        !----- set zero second derivative end condition
        A(1) = 2.0
        C(1) = 1.0
        XS(1) = 3.0 * (X(2) - X(1)) / (S(2) - S(1))
    ELSE IF(XS1.EQ.-999.0) THEN
        !----- set zero third derivative end condition
        A(1) = 1.0
        C(1) = 1.0
        XS(1) = 2.0 * (X(2) - X(1)) / (S(2) - S(1))
    ELSE
        !----- set specified first derivative end condition
        A(1) = 1.0
        C(1) = 0.
        XS(1) = XS1
    ENDIF
    !
    IF(XS2.EQ.999.0) THEN
        B(N) = 1.0
        A(N) = 2.0
        XS(N) = 3.0 * (X(N) - X(N - 1)) / (S(N) - S(N - 1))
    ELSE IF(XS2.EQ.-999.0) THEN
        B(N) = 1.0
        A(N) = 1.0
        XS(N) = 2.0 * (X(N) - X(N - 1)) / (S(N) - S(N - 1))
    ELSE
        A(N) = 1.0
        B(N) = 0.
        XS(N) = XS2
    ENDIF
    !
    IF(N.EQ.2 .AND. XS1.EQ.-999.0 .AND. XS2.EQ.-999.0) THEN
        !----- cannot have zero 3rd derivative at both endpoints of a single interval
        B(N) = 1.0
        A(N) = 2.0
        XS(N) = 3.0 * (X(N) - X(N - 1)) / (S(N) - S(N - 1))
    ENDIF
    !
    !---- solve for derivative array XS
    CALL TRISOL(A, B, C, XS, N)
    !
    RETURN
END
! SPLIND


SUBROUTINE SPLINA(X, XS, S, N)
    DIMENSION X(N), XS(N), S(N)
    LOGICAL LEND
    !-------------------------------------------------------
    !     Calculates spline coefficients for X(S) by a      |
    !     simple averaging of adjacent segment slopes.      |
    !                                                       |
    !     Interpolated X(S) is less likely to oscillate     |
    !     than with SPLINE, but does not have continuity    |
    !     in curvature.                                     |
    !                                                       |
    !     To evaluate the spline at some value of S,        |
    !     use SEVAL and/or DEVAL.                           |
    !                                                       |
    !     S        independent variable array (input)       |
    !     X        dependent variable array   (input)       |
    !     XS       dX/dS array                (calculated)  |
    !     N        number of points           (input)       |
    !                                                       |
    !-------------------------------------------------------
    !
    IF(N.EQ.1) THEN
        XS(1) = 0.
        RETURN
    ENDIF
    !
    LEND = .TRUE.
    DO 1 I = 1, N - 1
        DS = S(I + 1) - S(I)
        IF (DS.EQ.0.) THEN
            XS(I) = XS1
            LEND = .TRUE.
        ELSE
            DX = X(I + 1) - X(I)
            XS2 = DX / DS
            IF (LEND) THEN
                XS(I) = XS2
                LEND = .FALSE.
            ELSE
                XS(I) = 0.5 * (XS1 + XS2)
            ENDIF
        ENDIF
        XS1 = XS2
    1 CONTINUE
    XS(N) = XS1
    !
    RETURN
END
! SPLINA


SUBROUTINE TRISOL(A, B, C, D, KK)
    DIMENSION A(KK), B(KK), C(KK), D(KK)
    !-----------------------------------------
    !     Solves KK long, tri-diagonal system |
    !                                         |
    !             A C          D              |
    !             B A C        D              |
    !               B A .      .              |
    !                 . . C    .              |
    !                   B A    D              |
    !                                         |
    !     The righthand side D is replaced by |
    !     the solution.  A, C are destroyed.  |
    !-----------------------------------------
    !
    DO 1 K = 2, KK
        KM = K - 1
        C(KM) = C(KM) / A(KM)
        D(KM) = D(KM) / A(KM)
        A(K) = A(K) - B(K) * C(KM)
        D(K) = D(K) - B(K) * D(KM)
    1 CONTINUE
    !
    D(KK) = D(KK) / A(KK)
    !
    DO 2 K = KK - 1, 1, -1
        D(K) = D(K) - C(K) * D(K + 1)
    2 CONTINUE
    !
    RETURN
END
! TRISOL


FUNCTION GEVAL(SS, X, XS, S, N)
    DIMENSION X(N), XS(N), S(N)
    !--------------------------------------------------
    !     Calculates int( X(SS) ) dS                   |
    !     XS array must have been calculated by SPLINE |
    !--------------------------------------------------
    IF(N.EQ.1) THEN
        GEVAL = X(1) * (SS - S(1))
        RETURN
    ENDIF
    !
    ILOW = 1
    I = N
    !
    10 IF(I - ILOW .LE. 1) GO TO 11
    !
    IMID = (I + ILOW) / 2
    IF(SS .LT. S(IMID)) THEN
        I = IMID
    ELSE
        ILOW = IMID
    ENDIF
    GO TO 10
    !
    11 CONTINUE
    !
    !---- first integrate up to I-1 point
    GEVAL = 0.
    DO K = 2, I - 1
        DS = S(K) - S(K - 1)
        !
        !------ Int X(t) dt  for t = 0..1
        DGEV = 0.5 * (X(K) + X(K - 1)) + (XS(K - 1) - XS(K)) * DS / 12.0
        !
        GEVAL = GEVAL + DGEV * DS
    ENDDO
    !
    !---- now integrate up to SS value in I-1..I interval
    DS = S(I) - S(I - 1)
    T = (SS - S(I - 1)) / DS
    CX1 = DS * XS(I - 1) - X(I) + X(I - 1)
    CX2 = DS * XS(I) - X(I) + X(I - 1)
    !
    DGEV = 0.5 * T * T * X(I)&
            + (T - 0.5 * T * T) * X(I - 1)&
            + (6.0 - 8.0 * T + 3.0 * T * T) * T * T * CX1 / 12.0&
            + (- 4.0 * T + 3.0 * T * T) * T * T * CX2 / 12.0
    !
    GEVAL = GEVAL + DGEV * DS
    !
    RETURN
END
! GEVAL


FUNCTION SEVAL(SS, X, XS, S, N)
    DIMENSION X(N), XS(N), S(N)
    !--------------------------------------------------
    !     Calculates X(SS)                             |
    !     XS array must have been calculated by SPLINE |
    !--------------------------------------------------
    IF(N.EQ.1) THEN
        SEVAL = X(1) + XS(1) * (SS - S(1))
        RETURN
    ENDIF
    !
    ILOW = 1
    I = N
    !
    10 IF(I - ILOW .LE. 1) GO TO 11
    !
    IMID = (I + ILOW) / 2
    IF(SS .LT. S(IMID)) THEN
        I = IMID
    ELSE
        ILOW = IMID
    ENDIF
    GO TO 10
    !
    11 DS = S(I) - S(I - 1)
    T = (SS - S(I - 1)) / DS
    CX1 = DS * XS(I - 1) - X(I) + X(I - 1)
    CX2 = DS * XS(I) - X(I) + X(I - 1)
    SEVAL = T * X(I) + (1.0 - T) * X(I - 1) + (T - T * T) * ((1.0 - T) * CX1 - T * CX2)
    RETURN
END
! SEVAL


FUNCTION DEVAL(SS, X, XS, S, N)
    DIMENSION X(N), XS(N), S(N)
    !--------------------------------------------------
    !     Calculates dX/dS(SS)                         |
    !     XS array must have been calculated by SPLINE |
    !--------------------------------------------------
    IF(N.EQ.1) THEN
        DEVAL = XS(1)
        RETURN
    ENDIF
    !
    ILOW = 1
    I = N
    !
    10 IF(I - ILOW .LE. 1) GO TO 11
    !
    IMID = (I + ILOW) / 2
    IF(SS .LT. S(IMID)) THEN
        I = IMID
    ELSE
        ILOW = IMID
    ENDIF
    GO TO 10
    !
    11 DS = S(I) - S(I - 1)
    T = (SS - S(I - 1)) / DS
    CX1 = DS * XS(I - 1) - X(I) + X(I - 1)
    CX2 = DS * XS(I) - X(I) + X(I - 1)
    DEVAL = X(I) - X(I - 1) + (1. - 4.0 * T + 3.0 * T * T) * CX1 + T * (3.0 * T - 2.) * CX2
    DEVAL = DEVAL / DS
    RETURN
END
! DEVAL


FUNCTION D2VAL(SS, X, XS, S, N)
    DIMENSION X(N), XS(N), S(N)
    !--------------------------------------------------
    !     Calculates d2X/dS2(SS)                       |
    !     XS array must have been calculated by SPLINE |
    !--------------------------------------------------
    IF(N.EQ.1) THEN
        D2VAL = 0.0
        RETURN
    ENDIF
    !
    ILOW = 1
    I = N
    !
    10 IF(I - ILOW .LE. 1) GO TO 11
    !
    IMID = (I + ILOW) / 2
    IF(SS .LT. S(IMID)) THEN
        I = IMID
    ELSE
        ILOW = IMID
    ENDIF
    GO TO 10
    !
    11 DS = S(I) - S(I - 1)
    T = (SS - S(I - 1)) / DS
    CX1 = DS * XS(I - 1) - X(I) + X(I - 1)
    CX2 = DS * XS(I) - X(I) + X(I - 1)
    D2VAL = (6. * T - 4.) * CX1 + (6. * T - 2.0) * CX2
    D2VAL = D2VAL / DS**2
    RETURN
END
! D2VAL


SUBROUTINE SEVALL(SS, X, XS, S, N, &
        XX, XXS, XXSS)
    DIMENSION X(N), XS(N), S(N)
    !--------------------------------------------------
    !     Calculates all spline derivatives.           |
    !     (Combines SEVAL, DEVAL, D2VAL)               |
    !     XS array must have been calculated by SPLINE |
    !--------------------------------------------------
    IF(N.EQ.1) THEN
        XX = X(1) + XS(1) * (SS - S(1))
        XXS = XS(1)
        XXSS = 0.
        RETURN
    ENDIF
    !
    ILOW = 1
    I = N
    !
    10 IF(I - ILOW .LE. 1) GO TO 11
    !
    IMID = (I + ILOW) / 2
    IF(SS .LT. S(IMID)) THEN
        I = IMID
    ELSE
        ILOW = IMID
    ENDIF
    GO TO 10
    !
    11 DS = S(I) - S(I - 1)
    T = (SS - S(I - 1)) / DS
    !
    F0 = X(I - 1)
    F1 = DS * XS(I - 1)
    F2 = -DS * (2.0 * XS(I - 1) + XS(I)) + 3.0 * (X(I) - X(I - 1))
    F3 = DS * (XS(I - 1) + XS(I)) - 2.0 * (X(I) - X(I - 1))
    !
    XX = F0 + T * (F1 + T * (F2 + T * F3))
    XXS = F1 + T * (2.0 * F2 + T * 3.0 * F3)
    XXSS = 2.0 * F2 + T * 6.0 * F3
    !
    XXS = XXS / DS
    XXSS = XXSS / DS**2
    !
    RETURN
END
! SEVALL



SUBROUTINE SEVLIN(SS, X, S, N, XX, XXS)
    DIMENSION X(N), S(N)
    !------------------------------------------------------------
    !     Calculates X(SS) and dX/ds(SS) using piecewise-linear  |
    !     interpolation. This is intended for intepolating very  |
    !     noisy data for which a cubic spline is inappropriate.  |
    !------------------------------------------------------------
    IF(N.EQ.1) THEN
        XX = X(1)
        XXS = 0.
        RETURN
    ENDIF
    !
    ILOW = 1
    I = N
    !
    10 IF(I - ILOW .LE. 1) GO TO 11
    !
    IMID = (I + ILOW) / 2
    IF(SS .LT. S(IMID)) THEN
        I = IMID
    ELSE
        ILOW = IMID
    ENDIF
    GO TO 10
    !
    11 DS = S(I) - S(I - 1)
    T = (SS - S(I - 1)) / DS
    XX = T * X(I) + (1.0 - T) * X(I - 1)
    XXS = (X(I) - X(I - 1)) / DS
    !
    RETURN
END
! SEVLIN



FUNCTION CURV(SS, X, XS, Y, YS, S, N)
    DIMENSION X(N), XS(N), Y(N), YS(N), S(N)
    !-----------------------------------------------
    !     Calculates curvature of splined 2-D curve |
    !     at S = SS                                 |
    !                                               |
    !     S        arc length array of curve        |
    !     X, Y     coordinate arrays of curve       |
    !     XS,YS    derivative arrays                |
    !              (calculated earlier by SPLINE)   |
    !-----------------------------------------------
    IF(N.EQ.1) THEN
        CURV = 0.
        RETURN
    ENDIF
    !
    ILOW = 1
    I = N
    !
    10 IF(I - ILOW .LE. 1) GO TO 11
    !
    IMID = (I + ILOW) / 2
    IF(SS .LT. S(IMID)) THEN
        I = IMID
    ELSE
        ILOW = IMID
    ENDIF
    GO TO 10
    !
    11 DS = S(I) - S(I - 1)
    T = (SS - S(I - 1)) / DS
    !
    F1 = DS * XS(I - 1)
    F2 = -DS * (2.0 * XS(I - 1) + XS(I)) + 3.0 * (X(I) - X(I - 1))
    F3 = DS * (XS(I - 1) + XS(I)) - 2.0 * (X(I) - X(I - 1))
    !
    XD = F1 + T * (2.0 * F2 + T * 3.0 * F3)
    XDD = 2.0 * F2 + T * 6.0 * F3
    !
    !
    G1 = DS * YS(I - 1)
    G2 = -DS * (2.0 * YS(I - 1) + YS(I)) + 3.0 * (Y(I) - Y(I - 1))
    G3 = DS * (YS(I - 1) + YS(I)) - 2.0 * (Y(I) - Y(I - 1))
    !
    YD = G1 + T * (2.0 * G2 + T * 3.0 * G3)
    YDD = 2.0 * G2 + T * 6.0 * G3
    !
    !
    CURV = (XD * YDD - YD * XDD) / SQRT((XD * XD + YD * YD)**3)
    !
    RETURN
END
! CURV


FUNCTION CURVS(SS, X, XS, Y, YS, S, N)
    DIMENSION X(N), XS(N), Y(N), YS(N), S(N)
    !-----------------------------------------------
    !     Calculates curvature derivative of        |
    !     splined 2-D curve at S = SS               |
    !                                               |
    !     S        arc length array of curve        |
    !     X, Y     coordinate arrays of curve       |
    !     XS,YS    derivative arrays                |
    !              (calculated earlier by SPLINE)   |
    !-----------------------------------------------
    IF(N.EQ.1) THEN
        CURVS = 0.
        RETURN
    ENDIF
    !
    ILOW = 1
    I = N
    !
    10 IF(I - ILOW .LE. 1) GO TO 11
    !
    IMID = (I + ILOW) / 2
    IF(SS .LT. S(IMID)) THEN
        I = IMID
    ELSE
        ILOW = IMID
    ENDIF
    GO TO 10
    !
    11 DS = S(I) - S(I - 1)
    T = (SS - S(I - 1)) / DS
    !
    CX1 = DS * XS(I - 1) - X(I) + X(I - 1)
    CX2 = DS * XS(I) - X(I) + X(I - 1)
    XD = X(I) - X(I - 1) + (1.0 - 4.0 * T + 3.0 * T * T) * CX1 + T * (3.0 * T - 2.0) * CX2
    XDD = (6.0 * T - 4.0) * CX1 + (6.0 * T - 2.0) * CX2
    XDDD = 6.0 * CX1 + 6.0 * CX2
    !
    CY1 = DS * YS(I - 1) - Y(I) + Y(I - 1)
    CY2 = DS * YS(I) - Y(I) + Y(I - 1)
    YD = Y(I) - Y(I - 1) + (1.0 - 4.0 * T + 3.0 * T * T) * CY1 + T * (3.0 * T - 2.0) * CY2
    YDD = (6.0 * T - 4.0) * CY1 + (6.0 * T - 2.0) * CY2
    YDDD = 6.0 * CY1 + 6.0 * CY2
    !

    F1 = DS * XS(I - 1)
    F2 = -DS * (2.0 * XS(I - 1) + XS(I)) + 3.0 * (X(I) - X(I - 1))
    F3 = DS * (XS(I - 1) + XS(I)) - 2.0 * (X(I) - X(I - 1))
    !
    XD = F1 + T * (2.0 * F2 + T * 3.0 * F3)
    XDD = 2.0 * F2 + T * 6.0 * F3
    XDDD = 6.0 * F3
    !
    !
    G1 = DS * YS(I - 1)
    G2 = -DS * (2.0 * YS(I - 1) + YS(I)) + 3.0 * (Y(I) - Y(I - 1))
    G3 = DS * (YS(I - 1) + YS(I)) - 2.0 * (Y(I) - Y(I - 1))
    !
    YD = G1 + T * (2.0 * G2 + T * 3.0 * G3)
    YDD = 2.0 * G2 + T * 6.0 * G3
    YDDD = 6.0 * G3
    !
    SQRTB = SQRT(XD * XD + YD * YD)
    BOT = SQRTB**3
    DBOTDT = 3.0 * SQRTB * (XD * XDD + YD * YDD)
    !
    TOP = XD * YDD - YD * XDD
    DTOPDT = XD * YDDD - YD * XDDD
    !
    CURVS = (DTOPDT * BOT - DBOTDT * TOP) / BOT**2 / DS
    !
    RETURN
END
! CURVS


SUBROUTINE SINVRT(SI, XI, X, XS, S, N)
    DIMENSION X(N), XS(N), S(N)
    !----------------------------------------------------
    !     Calculates the "inverse" spline function S(X). |
    !     Since S(X) can be multi-valued or not defined, |
    !      this is not a "black-box" routine.  The call- |
    !      ing program must pass via SI a sufficiently   |
    !      good initial guess for S(XI).                 |
    !                                                    |
    !     XI      specified X value       (input)        |
    !     SI      calculated S(XI) value  (input,output) |
    !     X,XS,S  usual spline arrays     (input)        |
    !                                                    |
    !----------------------------------------------------
    !
    DO 10 ITER = 1, 10
        CALL SEVALL(SI, X, XS, S, N, XX, XXS, XXSS)
        IF(XXS.EQ.0.0) GO TO 11
        !
        DS = (XI - XX) / XXS
        SI = SI + DS
        IF(ABS(DS / (S(N) - S(1))) .LT. 1.0E-5) RETURN
    10   CONTINUE
    11   WRITE(*, *) 'SINVRT: spline inversion failed.  Continuing...'
    RETURN
    !
END
! SINVRT


SUBROUTINE SCALC(X, Y, S, N)
    DIMENSION X(N), Y(N), S(N)
    !----------------------------------------
    !     Calculates the arc length array S  |
    !     for a 2-D array of points (X,Y).   |
    !----------------------------------------
    !
    S(1) = 0.
    DO 10 I = 2, N
        S(I) = S(I - 1) + SQRT((X(I) - X(I - 1))**2 + (Y(I) - Y(I - 1))**2)
    10 CONTINUE
    !
    RETURN
END
! SCALC


SUBROUTINE SEGSPL(X, XS, S, N)
    DIMENSION X(N), XS(N), S(N)
    !-----------------------------------------------
    !     Splines X(S) array just like SPLINE,      |
    !     but allows derivative discontinuities     |
    !     at segment joints.  Segment joints are    |
    !     defined by identical successive S values. |
    !-----------------------------------------------
    IF(N.EQ.1) THEN
        XS(1) = 0.
        RETURN
    ENDIF
    !
    XX = 1.0 / (S(2) - S(1))
    IF(S(1).EQ.S(2)) THEN
        XX = 1.0 / (S(2) - S(1))
        !c        STOP 'SEGSPL:  First input point duplicated'
    ENDIF
    IF(S(N).EQ.S(N - 1)) THEN
        XX = 1.0 / (S(N) - S(N - 1))
        !c        STOP 'SEGSPL:  Last  input point duplicated'
    ENDIF
    !
    ISEG0 = 1
    DO 10 ISEG = 2, N - 2
        IF(S(ISEG).EQ.S(ISEG + 1)) THEN
            NSEG = ISEG - ISEG0 + 1
            CALL SPLINE(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG)
            ISEG0 = ISEG + 1
        ENDIF
    10 CONTINUE
    !
    NSEG = N - ISEG0 + 1
    CALL SPLINE(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG)
    !
    RETURN
END
! SEGSPL


SUBROUTINE SEGSPD(X, XS, S, N, XS1, XS2)
    DIMENSION X(N), XS(N), S(N)
    !-----------------------------------------------
    !     Splines X(S) array just like SPLIND,      |
    !     but allows derivative discontinuities     |
    !     at segment joints.  Segment joints are    |
    !     defined by identical successive S values. |
    !-----------------------------------------------
    IF(N.EQ.1) THEN
        XS(1) = 0.
        RETURN
    ENDIF
    !
    IF(S(1).EQ.S(2)) STOP 'SEGSPD:  First input point duplicated'
    IF(S(N).EQ.S(N - 1)) STOP 'SEGSPD:  Last  input point duplicated'
    !
    ISEG0 = 1
    DO 10 ISEG = 2, N - 2
        IF(S(ISEG).EQ.S(ISEG + 1)) THEN
            NSEG = ISEG - ISEG0 + 1
            CALL SPLIND(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG, XS1, XS2)
            ISEG0 = ISEG + 1
        ENDIF
    10 CONTINUE
    !
    NSEG = N - ISEG0 + 1
    CALL SPLIND(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG, XS1, XS2)
    !
    RETURN
END
! SEGSPD



SUBROUTINE INTERS(OK, SS1, SS2, &
        X1, XS1, Y1, YS1, S1, N1, &
        X2, XS2, Y2, YS2, S2, N2)
    LOGICAL OK
    DIMENSION X1(N1), XS1(N1), Y1(N1), YS1(N1), S1(N1)
    DIMENSION X2(N2), XS2(N2), Y2(N2), YS2(N2), S2(N2)
    !-------------------------------------------------------
    !     Finds spline coordinate values SS1, SS2 at the
    !     intersection of two space curves (X1,Y1), (X2,Y2).
    !-------------------------------------------------------
    LOGICAL CLIP1, CLIP2
    DATA EPS / 1.0E-5 /
    !
    OK = .TRUE.
    !cc      SS1 = S1(1)
    !cc      SS2 = S2(1)
    RS1 = 1.0E12
    RS2 = 1.0E12
    DS1 = 0.0
    DS2 = 0.0
    !
    DO 1000 ITER = 1, 12
        !
        RLX = 1.0
        SS1OLD = SS1
        SS2OLD = SS2
        RS1OLD = ABS(RS1)
        RS2OLD = ABS(RS2)
        !
        DO 10 IRLX = 1, 16
            !
            CLIP1 = .FALSE.
            CLIP2 = .FALSE.
            SS1 = SS1OLD + RLX * DS1
            SS2 = SS2OLD + RLX * DS2
            !
            IF(SS1.LT.S1(1) .OR. SS1.GT.S1(N1)) THEN
                CLIP1 = .TRUE.
                SS1 = MAX(SS1, S1(1))
                SS1 = MIN(SS1, S1(N1))
            ENDIF
            IF(SS2.LT.S2(1) .OR. SS2.GT.S2(N2)) THEN
                CLIP2 = .TRUE.
                SS2 = MAX(SS2, S2(1))
                SS2 = MIN(SS2, S2(N2))
            ENDIF
            !
            XX1 = SEVAL(SS1, X1, XS1, S1, N1)
            XX2 = SEVAL(SS2, X2, XS2, S2, N2)
            YY1 = SEVAL(SS1, Y1, YS1, S1, N1)
            YY2 = SEVAL(SS2, Y2, YS2, S2, N2)
            !
            RS1 = XX1 - XX2
            RS2 = YY1 - YY2
            !
            IF(ABS(RS1).LT.RS1OLD .AND.&
                    ABS(RS2).LT.RS2OLD) GO TO 11
            !
            RLX = 0.5 * RLX
            !
        10     CONTINUE
        WRITE(*, *) 'INTERS: Under-relaxation loop failed.'
        11     CONTINUE
        !
        A11 = DEVAL(SS1, X1, XS1, S1, N1)
        A12 = -DEVAL(SS2, X2, XS2, S2, N2)
        A21 = DEVAL(SS1, Y1, YS1, S1, N1)
        A22 = -DEVAL(SS2, Y2, YS2, S2, N2)
        !
        DET = A11 * A22 - A12 * A21
        DS1 = -(RS1 * A22 - A12 * RS2) / DET
        DS2 = -(A11 * RS2 - RS1 * A21) / DET
        !
        IF(ABS(DS1) .LT. EPS * (S1(N1) - S1(1)) .AND.&
                ABS(DS2) .LT. EPS * (S2(N2) - S2(1))) RETURN
        !
    1000 CONTINUE
    WRITE(*, *) 'INTERS: Convergence failed. Res =', RS1, RS2
    IF(CLIP1)&
            WRITE(*, *)'        S1 clip:', S1(1), S1(N1), SS1, DS1
    IF(CLIP2)&
            WRITE(*, *)'        S2 clip:', S2(1), S2(N2), SS2, DS2
    OK = .FALSE.
    !
    RETURN
END
! INTERS