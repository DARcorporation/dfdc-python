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
c
c    1-D Spline Package.
c    Interpolates a function x(s) over an s interval
c
c                                              Mark Drela
c                                              1985
c
c    Usage:
c
cC---- fill S(i), X(i) arrays
c      S(i) = ...
c      X(i) = ...
c
cC---- calculate spline coefficients for x(s) from discrete values
cC-    (can use SPLIND or SPLINA instead)
c      CALL SPLINE(X,XS,S,N)
c
cC---- evaluate splined x(s) and/or its derivatives
cC-    at any number of s points SS
c      XX   = SEVAL(SS,X,XS,S,N)
c      XXS  = DEVAL(SS,X,XS,S,N)
c      XXSS = D2VAL(SS,X,XS,S,N)
c
cC---- alternative to calling SEVAL,DEVAL,D2VAL separately
cC-    (slightly more efficient if all three quantities are needed)
c      CALL SEVALL(SS,X,XS,S,N, XX,XXS,XXSS)
c
c


      SUBROUTINE SPLINE(X,XS,S,N)
      DIMENSION X(N),XS(N),S(N)
      PARAMETER (NMAX=1601)
      DIMENSION A(NMAX),B(NMAX),C(NMAX)
C-------------------------------------------------------
C     Calculates spline coefficients for X(S).          |
C     Natural end conditions are used (zero 3rd         |
C      derivative over first, last intervals).          |
C                                                       |
C     To evaluate the spline at some value of S,        |
C     use SEVAL and/or DEVAL.                           |
C                                                       |
C     S        independent variable array (input)       |
C     X        dependent variable array   (input)       |
C     XS       dX/dS array                (calculated)  |
C     N        number of points           (input)       |
C                                                       |
C-------------------------------------------------------
      IF(N.GT.NMAX) STOP 'SPLINE: array overflow, increase NMAX'
C     
      IF(N.EQ.1) THEN
       XS(1) = 0.
       RETURN
      ENDIF
C
      DO 1 I=2, N-1
        DSM = S(I) - S(I-1)
        DSP = S(I+1) - S(I)
        B(I) = DSP
        A(I) = 2.0*(DSM+DSP)
        C(I) = DSM
        XS(I) = 3.0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
    1 CONTINUE
C
C---- set zero 3rd derivative end conditions
      A(1) = 1.0
      C(1) = 1.0
      XS(1) = 2.0*(X(2)-X(1)) / (S(2)-S(1))
C
      B(N) = 1.0
      A(N) = 1.0
      XS(N) = 2.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
C
      IF(N.EQ.2) THEN
C----- if only two points are present, specify zero 2nd derivative instead
C-     (straight line interpolation will result)
       B(N) = 1.0
       A(N) = 2.0
       XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ENDIF
C
C---- solve for derivative array XS
      CALL TRISOL(A,B,C,XS,N)
C
      RETURN
      END ! SPLINE



      SUBROUTINE SPLIND(X,XS,S,N,XS1,XS2)
      DIMENSION X(N),XS(N),S(N)
      PARAMETER (NMAX=1601)
      DIMENSION A(NMAX),B(NMAX),C(NMAX)
C-------------------------------------------------------
C     Calculates spline coefficients for X(S).          |
C     Same as SPLINE, but also allows specified-slope   |
C     or zero-curvature end conditions to be imposed.   |
C                                                       |
C     To evaluate the spline at some value of S,        |
C     use SEVAL and/or DEVAL.                           |
C                                                       |
C     S        independent variable array (input)       |
C     X        dependent variable array   (input)       |
C     XS       dX/dS array                (calculated)  |
C     N        number of points           (input)       |
C     XS1,XS2  endpoint derivatives       (input)       |
C              If =  999.0, then usual zero second      |
C              derivative end condition(s) are used     |
C              If = -999.0, then zero third             |
C              derivative end condition(s) are used     |
C                                                       |
C     Note: specifying both XS1,XS2 = -999.0            |
C           is equivalent to using SPLINE.              |
C                                                       |
C-------------------------------------------------------
      IF(N.GT.NMAX) STOP 'SPLIND: array overflow, increase NMAX'
C     
      IF(N.EQ.1) THEN
       XS(1) = 0.
       RETURN
      ENDIF
C
      DO 1 I=2, N-1
        DSM = S(I) - S(I-1)
        DSP = S(I+1) - S(I)
        B(I) = DSP
        A(I) = 2.0*(DSM+DSP)
        C(I) = DSM
        XS(I) = 3.0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
    1 CONTINUE
C
      IF(XS1.EQ.999.0) THEN
C----- set zero second derivative end condition
       A(1) = 2.0
       C(1) = 1.0
       XS(1) = 3.0*(X(2)-X(1)) / (S(2)-S(1))
      ELSE IF(XS1.EQ.-999.0) THEN
C----- set zero third derivative end condition
       A(1) = 1.0
       C(1) = 1.0
       XS(1) = 2.0*(X(2)-X(1)) / (S(2)-S(1))
      ELSE
C----- set specified first derivative end condition
       A(1) = 1.0
       C(1) = 0.
       XS(1) = XS1
      ENDIF
C
      IF(XS2.EQ.999.0) THEN
       B(N) = 1.0
       A(N) = 2.0
       XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ELSE IF(XS2.EQ.-999.0) THEN
       B(N) = 1.0
       A(N) = 1.0
       XS(N) = 2.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ELSE
       A(N) = 1.0
       B(N) = 0.
       XS(N) = XS2
      ENDIF
C
      IF(N.EQ.2 .AND. XS1.EQ.-999.0 .AND. XS2.EQ.-999.0) THEN
C----- cannot have zero 3rd derivative at both endpoints of a single interval
       B(N) = 1.0
       A(N) = 2.0
       XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ENDIF
C
C---- solve for derivative array XS
      CALL TRISOL(A,B,C,XS,N)
C
      RETURN
      END ! SPLIND


      SUBROUTINE SPLINA(X,XS,S,N)
      DIMENSION X(N),XS(N),S(N)
      LOGICAL LEND
C-------------------------------------------------------
C     Calculates spline coefficients for X(S) by a      |
C     simple averaging of adjacent segment slopes.      |
C                                                       |
C     Interpolated X(S) is less likely to oscillate     |
C     than with SPLINE, but does not have continuity    |
C     in curvature.                                     |
C                                                       |
C     To evaluate the spline at some value of S,        |
C     use SEVAL and/or DEVAL.                           |
C                                                       |
C     S        independent variable array (input)       |
C     X        dependent variable array   (input)       |
C     XS       dX/dS array                (calculated)  |
C     N        number of points           (input)       |
C                                                       |
C-------------------------------------------------------
C     
      IF(N.EQ.1) THEN
       XS(1) = 0.
       RETURN
      ENDIF
C
      LEND = .TRUE.
      DO 1 I=1, N-1
        DS = S(I+1)-S(I)
        IF (DS.EQ.0.) THEN
          XS(I) = XS1
          LEND = .TRUE.
         ELSE
          DX = X(I+1)-X(I)
          XS2 = DX / DS
          IF (LEND) THEN
            XS(I) = XS2
            LEND = .FALSE.
           ELSE
            XS(I) = 0.5*(XS1 + XS2)
          ENDIF
        ENDIF
        XS1 = XS2
    1 CONTINUE
      XS(N) = XS1
C
      RETURN
      END ! SPLINA


      SUBROUTINE TRISOL(A,B,C,D,KK)
      DIMENSION A(KK),B(KK),C(KK),D(KK)
C-----------------------------------------
C     Solves KK long, tri-diagonal system |
C                                         |
C             A C          D              |
C             B A C        D              |
C               B A .      .              |
C                 . . C    .              |
C                   B A    D              |
C                                         |
C     The righthand side D is replaced by |
C     the solution.  A, C are destroyed.  |
C-----------------------------------------
C
      DO 1 K=2, KK
        KM = K-1
        C(KM) = C(KM) / A(KM)
        D(KM) = D(KM) / A(KM)
        A(K) = A(K) - B(K)*C(KM)
        D(K) = D(K) - B(K)*D(KM)
    1 CONTINUE
C
      D(KK) = D(KK)/A(KK)
C
      DO 2 K=KK-1, 1, -1
        D(K) = D(K) - C(K)*D(K+1)
    2 CONTINUE
C
      RETURN
      END ! TRISOL


      FUNCTION GEVAL(SS,X,XS,S,N)
      DIMENSION X(N),XS(N),S(N)
C--------------------------------------------------
C     Calculates int( X(SS) ) dS                   |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      IF(N.EQ.1) THEN
       GEVAL = X(1)*(SS-S(1))
       RETURN
      ENDIF
C
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 CONTINUE
C
C---- first integrate up to I-1 point
      GEVAL = 0.
      DO K = 2, I-1
        DS = S(K) - S(K-1)
C
C------ Int X(t) dt  for t = 0..1
        DGEV = 0.5*(X(K) + X(K-1)) + (XS(K-1) - XS(K))*DS/12.0
C
        GEVAL = GEVAL + DGEV*DS
      ENDDO
C
C---- now integrate up to SS value in I-1..I interval
      DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
C
      DGEV =      0.5*T*T *X(I)
     &     + (T - 0.5*T*T)*X(I-1)
     &     + (6.0 - 8.0*T + 3.0*T*T)*T*T*CX1/12.0
     &     + (    - 4.0*T + 3.0*T*T)*T*T*CX2/12.0
C
      GEVAL = GEVAL + DGEV*DS
C
      RETURN
      END ! GEVAL


      FUNCTION SEVAL(SS,X,XS,S,N)
      DIMENSION X(N),XS(N),S(N)
C--------------------------------------------------
C     Calculates X(SS)                             |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      IF(N.EQ.1) THEN
       SEVAL = X(1) + XS(1)*(SS-S(1))
       RETURN
      ENDIF
C
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      SEVAL = T*X(I) + (1.0-T)*X(I-1) + (T-T*T)*((1.0-T)*CX1 - T*CX2)
      RETURN
      END ! SEVAL


      FUNCTION DEVAL(SS,X,XS,S,N)
      DIMENSION X(N),XS(N),S(N)
C--------------------------------------------------
C     Calculates dX/dS(SS)                         |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      IF(N.EQ.1) THEN
       DEVAL = XS(1)
       RETURN
      ENDIF
C
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      DEVAL = X(I) - X(I-1) + (1.-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.)*CX2
      DEVAL = DEVAL/DS
      RETURN
      END ! DEVAL


      FUNCTION D2VAL(SS,X,XS,S,N)
      DIMENSION X(N),XS(N),S(N)
C--------------------------------------------------
C     Calculates d2X/dS2(SS)                       |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      IF(N.EQ.1) THEN
       D2VAL = 0.0
       RETURN
      ENDIF
C
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      D2VAL = (6.*T-4.)*CX1 + (6.*T-2.0)*CX2
      D2VAL = D2VAL/DS**2
      RETURN
      END ! D2VAL


      SUBROUTINE SEVALL(SS,X,XS,S,N,
     &                  XX, XXS, XXSS )
      DIMENSION X(N),XS(N),S(N)
C--------------------------------------------------
C     Calculates all spline derivatives.           |
C     (Combines SEVAL, DEVAL, D2VAL)               |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      IF(N.EQ.1) THEN
       XX = X(1) + XS(1)*(SS-S(1))
       XXS = XS(1)
       XXSS = 0.
       RETURN
      ENDIF
C
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
C
      F0 = X(I-1)
      F1 = DS*XS(I-1)
      F2 = -DS*(2.0*XS(I-1) + XS(I)) + 3.0*(X(I) - X(I-1))
      F3 =  DS*(    XS(I-1) + XS(I)) - 2.0*(X(I) - X(I-1))
C
      XX = F0 + T*(F1 + T*(    F2 + T*    F3))
      XXS =        F1 + T*(2.0*F2 + T*3.0*F3)
      XXSS =               2.0*F2 + T*6.0*F3
C
      XXS = XXS/DS
      XXSS = XXSS/DS**2
C
      RETURN
      END ! SEVALL



      SUBROUTINE SEVLIN(SS,X,S,N, XX,XXS)
      DIMENSION X(N),S(N)
C------------------------------------------------------------
C     Calculates X(SS) and dX/ds(SS) using piecewise-linear  |
C     interpolation. This is intended for intepolating very  |
C     noisy data for which a cubic spline is inappropriate.  |
C------------------------------------------------------------
      IF(N.EQ.1) THEN
       XX = X(1)
       XXS = 0.
       RETURN
      ENDIF
C
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      XX = T*X(I) + (1.0-T)*X(I-1)
      XXS =  (X(I) - X(I-1))/DS
C
      RETURN
      END ! SEVLIN



      FUNCTION CURV(SS,X,XS,Y,YS,S,N)
      DIMENSION X(N), XS(N), Y(N), YS(N), S(N)
C-----------------------------------------------
C     Calculates curvature of splined 2-D curve |
C     at S = SS                                 |
C                                               |
C     S        arc length array of curve        |
C     X, Y     coordinate arrays of curve       |
C     XS,YS    derivative arrays                |
C              (calculated earlier by SPLINE)   |
C-----------------------------------------------
      IF(N.EQ.1) THEN
       CURV = 0.
       RETURN
      ENDIF
C     
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
C
      F1 = DS*XS(I-1)
      F2 = -DS*(2.0*XS(I-1) + XS(I)) + 3.0*(X(I) - X(I-1))
      F3 =  DS*(    XS(I-1) + XS(I)) - 2.0*(X(I) - X(I-1))
C
      XD = F1 + T*(2.0*F2 + T*3.0*F3)
      XDD =        2.0*F2 + T*6.0*F3
C
C
      G1 = DS*YS(I-1)
      G2 = -DS*(2.0*YS(I-1) + YS(I)) + 3.0*(Y(I) - Y(I-1))
      G3 =  DS*(    YS(I-1) + YS(I)) - 2.0*(Y(I) - Y(I-1))
C
      YD = G1 + T*(2.0*G2 + T*3.0*G3)
      YDD =        2.0*G2 + T*6.0*G3
C
C
      CURV = (XD*YDD - YD*XDD) / SQRT((XD*XD + YD*YD)**3)
C
      RETURN
      END ! CURV


      FUNCTION CURVS(SS,X,XS,Y,YS,S,N)
      DIMENSION X(N), XS(N), Y(N), YS(N), S(N)
C-----------------------------------------------
C     Calculates curvature derivative of        |
C     splined 2-D curve at S = SS               |
C                                               |
C     S        arc length array of curve        |
C     X, Y     coordinate arrays of curve       |
C     XS,YS    derivative arrays                |
C              (calculated earlier by SPLINE)   |
C-----------------------------------------------
      IF(N.EQ.1) THEN
       CURVS = 0.
       RETURN
      ENDIF
C     
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
C
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      XD = X(I) - X(I-1) + (1.0-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.0)*CX2
      XDD = (6.0*T-4.0)*CX1 + (6.0*T-2.0)*CX2
      XDDD = 6.0*CX1 + 6.0*CX2
C
      CY1 = DS*YS(I-1) - Y(I) + Y(I-1)
      CY2 = DS*YS(I)   - Y(I) + Y(I-1)
      YD = Y(I) - Y(I-1) + (1.0-4.0*T+3.0*T*T)*CY1 + T*(3.0*T-2.0)*CY2
      YDD = (6.0*T-4.0)*CY1 + (6.0*T-2.0)*CY2
      YDDD = 6.0*CY1 + 6.0*CY2
C


      F1 = DS*XS(I-1)
      F2 = -DS*(2.0*XS(I-1) + XS(I)) + 3.0*(X(I) - X(I-1))
      F3 =  DS*(    XS(I-1) + XS(I)) - 2.0*(X(I) - X(I-1))
C
      XD = F1 + T*(2.0*F2 + T*3.0*F3)
      XDD =        2.0*F2 + T*6.0*F3
      XDDD =                  6.0*F3
C
C
      G1 = DS*YS(I-1)
      G2 = -DS*(2.0*YS(I-1) + YS(I)) + 3.0*(Y(I) - Y(I-1))
      G3 =  DS*(    YS(I-1) + YS(I)) - 2.0*(Y(I) - Y(I-1))
C
      YD = G1 + T*(2.0*G2 + T*3.0*G3)
      YDD =        2.0*G2 + T*6.0*G3
      YDDD =                  6.0*G3
C
      SQRTB = SQRT(XD*XD + YD*YD)
      BOT = SQRTB**3
      DBOTDT = 3.0*SQRTB*(XD*XDD + YD*YDD)
C
      TOP = XD*YDD - YD*XDD      
      DTOPDT = XD*YDDD - YD*XDDD
C
      CURVS = (DTOPDT*BOT - DBOTDT*TOP) / BOT**2  / DS
C
      RETURN
      END ! CURVS


      SUBROUTINE SINVRT(SI,XI,X,XS,S,N)
      DIMENSION X(N),XS(N),S(N)
C----------------------------------------------------
C     Calculates the "inverse" spline function S(X). |
C     Since S(X) can be multi-valued or not defined, |
C      this is not a "black-box" routine.  The call- |
C      ing program must pass via SI a sufficiently   |
C      good initial guess for S(XI).                 |
C                                                    |
C     XI      specified X value       (input)        |
C     SI      calculated S(XI) value  (input,output) |
C     X,XS,S  usual spline arrays     (input)        |
C                                                    |
C----------------------------------------------------
C
      DO 10 ITER=1, 10
        CALL SEVALL(SI,X,XS,S,N, XX,XXS,XXSS)
        IF(XXS.EQ.0.0) GO TO 11
C
        DS = (XI-XX)/XXS
        SI = SI + DS
        IF(ABS(DS/(S(N)-S(1))) .LT. 1.0E-5) RETURN
 10   CONTINUE
 11   WRITE(*,*) 'SINVRT: spline inversion failed.  Continuing...'
      RETURN
C
      END ! SINVRT


      SUBROUTINE SCALC(X,Y,S,N)
      DIMENSION X(N),Y(N),S(N)
C----------------------------------------
C     Calculates the arc length array S  |
C     for a 2-D array of points (X,Y).   |
C----------------------------------------
C
      S(1) = 0.
      DO 10 I=2, N
        S(I) = S(I-1) + SQRT((X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2)
   10 CONTINUE
C
      RETURN
      END ! SCALC


      SUBROUTINE SEGSPL(X,XS,S,N)
      DIMENSION X(N),XS(N),S(N)
C-----------------------------------------------
C     Splines X(S) array just like SPLINE,      |
C     but allows derivative discontinuities     |
C     at segment joints.  Segment joints are    |
C     defined by identical successive S values. |
C-----------------------------------------------
      IF(N.EQ.1) THEN
       XS(1) = 0.
       RETURN
      ENDIF
C
      XX = 1.0/(S(2)-S(1))
      IF(S(1).EQ.S(2)  ) THEN
        XX = 1.0/(S(2)-S(1))
cc        STOP 'SEGSPL:  First input point duplicated'
      ENDIF
      IF(S(N).EQ.S(N-1)) THEN
        XX = 1.0/(S(N)-S(N-1))
cc        STOP 'SEGSPL:  Last  input point duplicated'
      ENDIF
C
      ISEG0 = 1
      DO 10 ISEG=2, N-2
        IF(S(ISEG).EQ.S(ISEG+1)) THEN
         NSEG = ISEG - ISEG0 + 1
         CALL SPLINE(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG)
         ISEG0 = ISEG+1
        ENDIF
   10 CONTINUE
C
      NSEG = N - ISEG0 + 1
      CALL SPLINE(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG)
C
      RETURN
      END ! SEGSPL


      SUBROUTINE SEGSPD(X,XS,S,N,XS1,XS2)
      DIMENSION X(N),XS(N),S(N)
C-----------------------------------------------
C     Splines X(S) array just like SPLIND,      |
C     but allows derivative discontinuities     |
C     at segment joints.  Segment joints are    |
C     defined by identical successive S values. |
C-----------------------------------------------
      IF(N.EQ.1) THEN
       XS(1) = 0.
       RETURN
      ENDIF
C
      IF(S(1).EQ.S(2)  ) STOP 'SEGSPD:  First input point duplicated'
      IF(S(N).EQ.S(N-1)) STOP 'SEGSPD:  Last  input point duplicated'
C
      ISEG0 = 1
      DO 10 ISEG=2, N-2
        IF(S(ISEG).EQ.S(ISEG+1)) THEN
         NSEG = ISEG - ISEG0 + 1
         CALL SPLIND(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG,XS1,XS2)
         ISEG0 = ISEG+1
        ENDIF
   10 CONTINUE
C
      NSEG = N - ISEG0 + 1
      CALL SPLIND(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG,XS1,XS2)
C
      RETURN
      END ! SEGSPD



      SUBROUTINE INTERS(OK,SS1,SS2,
     &                  X1,XS1,Y1,YS1,S1,N1,
     &                  X2,XS2,Y2,YS2,S2,N2 )
      LOGICAL OK
      DIMENSION X1(N1),XS1(N1),Y1(N1),YS1(N1),S1(N1)
      DIMENSION X2(N2),XS2(N2),Y2(N2),YS2(N2),S2(N2)
C-------------------------------------------------------
C     Finds spline coordinate values SS1, SS2 at the
C     intersection of two space curves (X1,Y1), (X2,Y2).
C-------------------------------------------------------
      LOGICAL CLIP1, CLIP2
      DATA EPS / 1.0E-5 /
C
      OK = .TRUE.
ccc      SS1 = S1(1)
ccc      SS2 = S2(1)
      RS1 = 1.0E12
      RS2 = 1.0E12
      DS1 = 0.0
      DS2 = 0.0
C
      DO 1000 ITER=1, 12
C
        RLX = 1.0
        SS1OLD = SS1
        SS2OLD = SS2
        RS1OLD = ABS(RS1)
        RS2OLD = ABS(RS2)
C
        DO 10 IRLX=1, 16
C
          CLIP1 = .FALSE.
          CLIP2 = .FALSE.
          SS1 = SS1OLD + RLX*DS1
          SS2 = SS2OLD + RLX*DS2
C          
          IF(SS1.LT.S1(1) .OR. SS1.GT.S1(N1)) THEN
           CLIP1 = .TRUE.
           SS1 = MAX(SS1,S1(1 ))
           SS1 = MIN(SS1,S1(N1))
          ENDIF
          IF(SS2.LT.S2(1) .OR. SS2.GT.S2(N2)) THEN
           CLIP2 = .TRUE.
           SS2 = MAX(SS2,S2(1 ))
           SS2 = MIN(SS2,S2(N2))
          ENDIF
C
          XX1 = SEVAL(SS1,X1,XS1,S1,N1)
          XX2 = SEVAL(SS2,X2,XS2,S2,N2)
          YY1 = SEVAL(SS1,Y1,YS1,S1,N1)
          YY2 = SEVAL(SS2,Y2,YS2,S2,N2)
C
          RS1 = XX1 - XX2
          RS2 = YY1 - YY2
C
          IF(ABS(RS1).LT.RS1OLD .AND.
     &       ABS(RS2).LT.RS2OLD     ) GO TO 11
C
          RLX = 0.5*RLX
C
 10     CONTINUE
        WRITE(*,*) 'INTERS: Under-relaxation loop failed.'
 11     CONTINUE
C
        A11 =  DEVAL(SS1,X1,XS1,S1,N1)
        A12 = -DEVAL(SS2,X2,XS2,S2,N2)
        A21 =  DEVAL(SS1,Y1,YS1,S1,N1)
        A22 = -DEVAL(SS2,Y2,YS2,S2,N2)
C
        DET =   A11*A22 - A12*A21
        DS1 = -(RS1*A22 - A12*RS2)/DET
        DS2 = -(A11*RS2 - RS1*A21)/DET
C
        IF(ABS(DS1) .LT. EPS*(S1(N1)-S1(1)) .AND.
     &     ABS(DS2) .LT. EPS*(S2(N2)-S2(1))      ) RETURN
C
 1000 CONTINUE
      WRITE(*,*) 'INTERS: Convergence failed. Res =', RS1, RS2
      IF(CLIP1)
     & WRITE(*,*)'        S1 clip:', S1(1), S1(N1), SS1, DS1
      IF(CLIP2)
     & WRITE(*,*)'        S2 clip:', S2(1), S2(N2), SS2, DS2
      OK = .FALSE.
C
      RETURN
      END ! INTERS
