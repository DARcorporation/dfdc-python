module m_spline
    implicit none
contains
    !*==SPLINE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GAMSOL
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
        INTEGER :: N
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        REAL, DIMENSION(NMAX) :: A, B, C
        REAL :: DSM, DSP
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
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
        IF (N>NMAX) STOP 'SPLINE: array overflow, increase NMAX'
        !
        IF (N==1) THEN
            XS(1) = 0.
            RETURN
        ENDIF
        !
        DO I = 2, N - 1
            DSM = S(I) - S(I - 1)
            DSP = S(I + 1) - S(I)
            B(I) = DSP
            A(I) = 2.0 * (DSM + DSP)
            C(I) = DSM
            XS(I) = 3.0 * ((X(I + 1) - X(I)) * DSM / DSP + (X(I) - X(I - 1)) * DSP / DSM)
        ENDDO
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
        IF (N==2) THEN
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
    END SUBROUTINE SPLINE
    !*==SPLIND.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SPLINE
    
    
    
    SUBROUTINE SPLIND(X, XS, S, N, XS1, XS2)
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
        INTEGER :: N
        REAL :: XS1, XS2
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        REAL, DIMENSION(NMAX) :: A, B, C
        REAL :: DSM, DSP
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
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
        IF (N>NMAX) STOP 'SPLIND: array overflow, increase NMAX'
        !
        IF (N==1) THEN
            XS(1) = 0.
            RETURN
        ENDIF
        !
        DO I = 2, N - 1
            DSM = S(I) - S(I - 1)
            DSP = S(I + 1) - S(I)
            B(I) = DSP
            A(I) = 2.0 * (DSM + DSP)
            C(I) = DSM
            XS(I) = 3.0 * ((X(I + 1) - X(I)) * DSM / DSP + (X(I) - X(I - 1)) * DSP / DSM)
        ENDDO
        !
        IF (XS1==999.0) THEN
            !----- set zero second derivative end condition
            A(1) = 2.0
            C(1) = 1.0
            XS(1) = 3.0 * (X(2) - X(1)) / (S(2) - S(1))
        ELSEIF (XS1==-999.0) THEN
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
        IF (XS2==999.0) THEN
            B(N) = 1.0
            A(N) = 2.0
            XS(N) = 3.0 * (X(N) - X(N - 1)) / (S(N) - S(N - 1))
        ELSEIF (XS2==-999.0) THEN
            B(N) = 1.0
            A(N) = 1.0
            XS(N) = 2.0 * (X(N) - X(N - 1)) / (S(N) - S(N - 1))
        ELSE
            A(N) = 1.0
            B(N) = 0.
            XS(N) = XS2
        ENDIF
        !
        IF (N==2 .AND. XS1==-999.0 .AND. XS2==-999.0) THEN
            !----- cannot have zero 3rd derivative at both endpoints of a single interval
            B(N) = 1.0
            A(N) = 2.0
            XS(N) = 3.0 * (X(N) - X(N - 1)) / (S(N) - S(N - 1))
        ENDIF
        !
        !---- solve for derivative array XS
        CALL TRISOL(A, B, C, XS, N)
        !
    END SUBROUTINE SPLIND
    !*==SPLINA.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SPLIND
    
    
    SUBROUTINE SPLINA(X, XS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        REAL :: DS, DX, XS1, XS2
        INTEGER :: I
        LOGICAL :: LEND
        !
        !*** End of declarations rewritten by SPAG
        !
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
        IF (N==1) THEN
            XS(1) = 0.
            RETURN
        ENDIF
        !
        LEND = .TRUE.
        DO I = 1, N - 1
            DS = S(I + 1) - S(I)
            IF (DS==0.) THEN
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
        ENDDO
        XS(N) = XS1
        !
    END SUBROUTINE SPLINA
    !*==TRISOL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SPLINA
    
    
    SUBROUTINE TRISOL(A, B, C, D, KK)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: KK
        REAL, DIMENSION(KK) :: A, B, C, D
        !
        ! Local variables
        !
        INTEGER :: K, KM
        !
        !*** End of declarations rewritten by SPAG
        !
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
        DO K = 2, KK
            KM = K - 1
            C(KM) = C(KM) / A(KM)
            D(KM) = D(KM) / A(KM)
            A(K) = A(K) - B(K) * C(KM)
            D(K) = D(K) - B(K) * D(KM)
        ENDDO
        !
        D(KK) = D(KK) / A(KK)
        !
        DO K = KK - 1, 1, -1
            D(K) = D(K) - C(K) * D(K + 1)
        ENDDO
        !
    END SUBROUTINE TRISOL
    !*==GEVAL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! TRISOL
    
    
    REAL FUNCTION GEVAL(SS, X, XS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SS
        REAL :: GEVAL
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        REAL :: CX1, CX2, DGEV, DS, T
        INTEGER :: I, ILOW, IMID, K
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates int( X(SS) ) dS                   |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        IF (N==1) THEN
            GEVAL = X(1) * (SS - S(1))
            RETURN
        ENDIF
        !
        ILOW = 1
        I = N
        !
        DO WHILE (I - ILOW>1)
            !
            IMID = (I + ILOW) / 2
            IF (SS<S(IMID)) THEN
                I = IMID
            ELSE
                ILOW = IMID
            ENDIF
        ENDDO
        !
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
        DGEV = 0.5 * T * T * X(I) + (T - 0.5 * T * T) * X(I - 1) + (6.0 - 8.0 * T + 3.0 * T * T)    &
                & * T * T * CX1 / 12.0 + (-4.0 * T + 3.0 * T * T) * T * T * CX2 / 12.0
        !
        GEVAL = GEVAL + DGEV * DS
        !
    END FUNCTION GEVAL
    !*==SEVAL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GEVAL
    
    
    REAL FUNCTION SEVAL(SS, X, XS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SS
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        REAL :: CX1, CX2, DS, T
        INTEGER :: I, ILOW, IMID
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates X(SS)                             |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        IF (N==1) THEN
            SEVAL = X(1) + XS(1) * (SS - S(1))
            RETURN
        ENDIF
        !
        ILOW = 1
        I = N
        !
        DO WHILE (I - ILOW>1)
            !
            IMID = (I + ILOW) / 2
            IF (SS<S(IMID)) THEN
                I = IMID
            ELSE
                ILOW = IMID
            ENDIF
        ENDDO
        !
        DS = S(I) - S(I - 1)
        T = (SS - S(I - 1)) / DS
        CX1 = DS * XS(I - 1) - X(I) + X(I - 1)
        CX2 = DS * XS(I) - X(I) + X(I - 1)
        SEVAL = T * X(I) + (1.0 - T) * X(I - 1) + (T - T * T) * ((1.0 - T) * CX1 - T * CX2)
    END FUNCTION SEVAL
    !*==DEVAL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEVAL
    
    
    REAL FUNCTION DEVAL(SS, X, XS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SS
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        REAL :: CX1, CX2, DS, T
        INTEGER :: I, ILOW, IMID
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates dX/dS(SS)                         |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        IF (N==1) THEN
            DEVAL = XS(1)
            RETURN
        ENDIF
        !
        ILOW = 1
        I = N
        !
        DO WHILE (I - ILOW>1)
            !
            IMID = (I + ILOW) / 2
            IF (SS<S(IMID)) THEN
                I = IMID
            ELSE
                ILOW = IMID
            ENDIF
        ENDDO
        !
        DS = S(I) - S(I - 1)
        T = (SS - S(I - 1)) / DS
        CX1 = DS * XS(I - 1) - X(I) + X(I - 1)
        CX2 = DS * XS(I) - X(I) + X(I - 1)
        DEVAL = X(I) - X(I - 1) + (1. - 4.0 * T + 3.0 * T * T) * CX1 + T * (3.0 * T - 2.) * CX2
        DEVAL = DEVAL / DS
    END FUNCTION DEVAL
    !*==D2VAL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! DEVAL
    
    
    REAL FUNCTION D2VAL(SS, X, XS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SS
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        REAL :: CX1, CX2, DS, T
        INTEGER :: I, ILOW, IMID
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates d2X/dS2(SS)                       |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        IF (N==1) THEN
            D2VAL = 0.0
            RETURN
        ENDIF
        !
        ILOW = 1
        I = N
        !
        DO WHILE (I - ILOW>1)
            !
            IMID = (I + ILOW) / 2
            IF (SS<S(IMID)) THEN
                I = IMID
            ELSE
                ILOW = IMID
            ENDIF
        ENDDO
        !
        DS = S(I) - S(I - 1)
        T = (SS - S(I - 1)) / DS
        CX1 = DS * XS(I - 1) - X(I) + X(I - 1)
        CX2 = DS * XS(I) - X(I) + X(I - 1)
        D2VAL = (6. * T - 4.) * CX1 + (6. * T - 2.0) * CX2
        D2VAL = D2VAL / DS**2
    END FUNCTION D2VAL
    !*==SEVALL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! D2VAL
    
    
    SUBROUTINE SEVALL(SS, X, XS, S, N, XX, XXS, XXSS)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SS, XX, XXS, XXSS
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        REAL :: DS, F0, F1, F2, F3, T
        INTEGER :: I, ILOW, IMID
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates all spline derivatives.           |
        !     (Combines SEVAL, DEVAL, D2VAL)               |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        IF (N==1) THEN
            XX = X(1) + XS(1) * (SS - S(1))
            XXS = XS(1)
            XXSS = 0.
            RETURN
        ENDIF
        !
        ILOW = 1
        I = N
        !
        DO WHILE (I - ILOW>1)
            !
            IMID = (I + ILOW) / 2
            IF (SS<S(IMID)) THEN
                I = IMID
            ELSE
                ILOW = IMID
            ENDIF
        ENDDO
        !
        DS = S(I) - S(I - 1)
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
    END SUBROUTINE SEVALL
    !*==SEVLIN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEVALL
    
    
    
    SUBROUTINE SEVLIN(SS, X, S, N, XX, XXS)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SS, XX, XXS
        REAL, DIMENSION(N) :: S, X
        !
        ! Local variables
        !
        REAL :: DS, T
        INTEGER :: I, ILOW, IMID
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------
        !     Calculates X(SS) and dX/ds(SS) using piecewise-linear  |
        !     interpolation. This is intended for intepolating very  |
        !     noisy data for which a cubic spline is inappropriate.  |
        !------------------------------------------------------------
        IF (N==1) THEN
            XX = X(1)
            XXS = 0.
            RETURN
        ENDIF
        !
        ILOW = 1
        I = N
        !
        DO WHILE (I - ILOW>1)
            !
            IMID = (I + ILOW) / 2
            IF (SS<S(IMID)) THEN
                I = IMID
            ELSE
                ILOW = IMID
            ENDIF
        ENDDO
        !
        DS = S(I) - S(I - 1)
        T = (SS - S(I - 1)) / DS
        XX = T * X(I) + (1.0 - T) * X(I - 1)
        XXS = (X(I) - X(I - 1)) / DS
        !
    END SUBROUTINE SEVLIN
    !*==CURV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEVLIN
    
    
    
    REAL FUNCTION CURV(SS, X, XS, Y, YS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SS
        REAL, DIMENSION(N) :: S, X, XS, Y, YS
        !
        ! Local variables
        !
        REAL :: DS, F1, F2, F3, G1, G2, G3, T, XD, XDD, YD, &
                & YDD
        INTEGER :: I, ILOW, IMID
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Calculates curvature of splined 2-D curve |
        !     at S = SS                                 |
        !                                               |
        !     S        arc length array of curve        |
        !     X, Y     coordinate arrays of curve       |
        !     XS,YS    derivative arrays                |
        !              (calculated earlier by SPLINE)   |
        !-----------------------------------------------
        IF (N==1) THEN
            CURV = 0.
            RETURN
        ENDIF
        !
        ILOW = 1
        I = N
        !
        DO WHILE (I - ILOW>1)
            !
            IMID = (I + ILOW) / 2
            IF (SS<S(IMID)) THEN
                I = IMID
            ELSE
                ILOW = IMID
            ENDIF
        ENDDO
        !
        DS = S(I) - S(I - 1)
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
    END FUNCTION CURV
    !*==CURVS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CURV
    
    
    REAL FUNCTION CURVS(SS, X, XS, Y, YS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SS
        REAL, DIMENSION(N) :: S, X, XS, Y, YS
        !
        ! Local variables
        !
        REAL :: BOT, CX1, CX2, CY1, CY2, DBOTDT, DS, DTOPDT, &
                & F1, F2, F3, G1, G2, G3, SQRTB, T, TOP, XD, &
                & XDD, XDDD, YD, YDD, YDDD
        INTEGER :: I, ILOW, IMID
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Calculates curvature derivative of        |
        !     splined 2-D curve at S = SS               |
        !                                               |
        !     S        arc length array of curve        |
        !     X, Y     coordinate arrays of curve       |
        !     XS,YS    derivative arrays                |
        !              (calculated earlier by SPLINE)   |
        !-----------------------------------------------
        IF (N==1) THEN
            CURVS = 0.
            RETURN
        ENDIF
        !
        ILOW = 1
        I = N
        !
        DO WHILE (I - ILOW>1)
            !
            IMID = (I + ILOW) / 2
            IF (SS<S(IMID)) THEN
                I = IMID
            ELSE
                ILOW = IMID
            ENDIF
        ENDDO
        !
        DS = S(I) - S(I - 1)
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
    END FUNCTION CURVS
    !*==SINVRT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CURVS
    
    
    SUBROUTINE SINVRT(SI, XI, X, XS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SI, XI
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        REAL :: DS, XX, XXS, XXSS
        INTEGER :: ITER
        !
        !*** End of declarations rewritten by SPAG
        !
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
        DO ITER = 1, 10
            CALL SEVALL(SI, X, XS, S, N, XX, XXS, XXSS)
            IF (XXS==0.0) EXIT
            !
            DS = (XI - XX) / XXS
            SI = SI + DS
            IF (ABS(DS / (S(N) - S(1)))<1.0E-5) RETURN
        ENDDO
        WRITE (*, *) 'SINVRT: spline inversion failed.  Continuing...'
        !
    END SUBROUTINE SINVRT
    !*==SCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SINVRT
    
    
    SUBROUTINE SCALC(X, Y, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(N) :: S, X, Y
        !
        ! Local variables
        !
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------
        !     Calculates the arc length array S  |
        !     for a 2-D array of points (X,Y).   |
        !----------------------------------------
        !
        S(1) = 0.
        DO I = 2, N
            S(I) = S(I - 1) + SQRT((X(I) - X(I - 1))**2 + (Y(I) - Y(I - 1))**2)
        ENDDO
        !
    END SUBROUTINE SCALC
    !*==SEGSPL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SCALC
    
    
    SUBROUTINE SEGSPL(X, XS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        INTEGER :: ISEG, ISEG0, NSEG
        REAL :: XX
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Splines X(S) array just like SPLINE,      |
        !     but allows derivative discontinuities     |
        !     at segment joints.  Segment joints are    |
        !     defined by identical successive S values. |
        !-----------------------------------------------
        IF (N==1) THEN
            XS(1) = 0.
            RETURN
        ENDIF
        !
        XX = 1.0 / (S(2) - S(1))
        !c        STOP 'SEGSPL:  First input point duplicated'
        IF (S(1)==S(2)) XX = 1.0 / (S(2) - S(1))
        !c        STOP 'SEGSPL:  Last  input point duplicated'
        IF (S(N)==S(N - 1)) XX = 1.0 / (S(N) - S(N - 1))
        !
        ISEG0 = 1
        DO ISEG = 2, N - 2
            IF (S(ISEG)==S(ISEG + 1)) THEN
                NSEG = ISEG - ISEG0 + 1
                CALL SPLINE(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG)
                ISEG0 = ISEG + 1
            ENDIF
        ENDDO
        !
        NSEG = N - ISEG0 + 1
        CALL SPLINE(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG)
        !
    END SUBROUTINE SEGSPL
    !*==SEGSPD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEGSPL
    
    
    SUBROUTINE SEGSPD(X, XS, S, N, XS1, XS2)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: XS1, XS2
        REAL, DIMENSION(N) :: S, X, XS
        !
        ! Local variables
        !
        INTEGER :: ISEG, ISEG0, NSEG
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Splines X(S) array just like SPLIND,      |
        !     but allows derivative discontinuities     |
        !     at segment joints.  Segment joints are    |
        !     defined by identical successive S values. |
        !-----------------------------------------------
        IF (N==1) THEN
            XS(1) = 0.
            RETURN
        ENDIF
        !
        IF (S(1)==S(2)) STOP 'SEGSPD:  First input point duplicated'
        IF (S(N)==S(N - 1)) STOP 'SEGSPD:  Last  input point duplicated'
        !
        ISEG0 = 1
        DO ISEG = 2, N - 2
            IF (S(ISEG)==S(ISEG + 1)) THEN
                NSEG = ISEG - ISEG0 + 1
                CALL SPLIND(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG, XS1, XS2)
                ISEG0 = ISEG + 1
            ENDIF
        ENDDO
        !
        NSEG = N - ISEG0 + 1
        CALL SPLIND(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG, XS1, XS2)
        !
    END SUBROUTINE SEGSPD
    !*==INTERS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEGSPD
    
    
    
    SUBROUTINE INTERS(OK, SS1, SS2, X1, XS1, Y1, YS1, S1, N1, X2, XS2, Y2, YS2, S2, &
            & N2)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N1, N2
        LOGICAL :: OK
        REAL :: SS1, SS2
        REAL, DIMENSION(N1) :: S1, X1, XS1, Y1, YS1
        REAL, DIMENSION(N2) :: S2, X2, XS2, Y2, YS2
        !
        ! Local variables
        !
        REAL :: A11, A12, A21, A22, DET, DS1, DS2, RLX, RS1, &
                & RS1OLD, RS2, RS2OLD, SS1OLD, SS2OLD, XX1, XX2, &
                & YY1, YY2
        LOGICAL :: CLIP1, CLIP2
        REAL, SAVE :: EPS
        INTEGER :: IRLX, ITER
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Finds spline coordinate values SS1, SS2 at the
        !     intersection of two space curves (X1,Y1), (X2,Y2).
        !-------------------------------------------------------
        DATA EPS/1.0E-5/
        !
        OK = .TRUE.
        !cc      SS1 = S1(1)
        !cc      SS2 = S2(1)
        RS1 = 1.0E12
        RS2 = 1.0E12
        DS1 = 0.0
        DS2 = 0.0
        !
        DO ITER = 1, 12
            !
            RLX = 1.0
            SS1OLD = SS1
            SS2OLD = SS2
            RS1OLD = ABS(RS1)
            RS2OLD = ABS(RS2)
            !
            DO IRLX = 1, 16
                !
                CLIP1 = .FALSE.
                CLIP2 = .FALSE.
                SS1 = SS1OLD + RLX * DS1
                SS2 = SS2OLD + RLX * DS2
                !
                IF (SS1<S1(1) .OR. SS1>S1(N1)) THEN
                    CLIP1 = .TRUE.
                    SS1 = MAX(SS1, S1(1))
                    SS1 = MIN(SS1, S1(N1))
                ENDIF
                IF (SS2<S2(1) .OR. SS2>S2(N2)) THEN
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
                IF (ABS(RS1)<RS1OLD .AND. ABS(RS2)<RS2OLD) GOTO 11
                !
                RLX = 0.5 * RLX
                !
            ENDDO
            WRITE (*, *) 'INTERS: Under-relaxation loop failed.'
            !
            11      A11 = DEVAL(SS1, X1, XS1, S1, N1)
            A12 = -DEVAL(SS2, X2, XS2, S2, N2)
            A21 = DEVAL(SS1, Y1, YS1, S1, N1)
            A22 = -DEVAL(SS2, Y2, YS2, S2, N2)
            !
            DET = A11 * A22 - A12 * A21
            DS1 = -(RS1 * A22 - A12 * RS2) / DET
            DS2 = -(A11 * RS2 - RS1 * A21) / DET
            !
            IF (ABS(DS1)<EPS * (S1(N1) - S1(1)) .AND. ABS(DS2)                &
                    & <EPS * (S2(N2) - S2(1))) RETURN
            !
        ENDDO
        WRITE (*, *) 'INTERS: Convergence failed. Res =', RS1, RS2
        IF (CLIP1) WRITE (*, *) '        S1 clip:', S1(1), S1(N1), &
                & SS1, DS1
        IF (CLIP2) WRITE (*, *) '        S2 clip:', S2(1), S2(N2), &
                & SS2, DS2
        OK = .FALSE.
        !
    END SUBROUTINE INTERS
end module m_spline
