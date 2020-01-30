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

      SUBROUTINE SETEXP(S,DS1,SMAX,NN)
C........................................................
C     Sets geometrically stretched array S:
C
C       S(i+1) - S(i)  =  r * [S(i) - S(i-1)]
C
C       S     (output)  array to be set  
C       DS1   (input)   first S increment:  S(2) - S(1)
C       SMAX  (input)   final S value:      S(NN)
C       NN    (input)   number of points
C........................................................
      REAL S(NN)
C
      SIGMA = SMAX/DS1
      NEX = NN-1
      RNEX = FLOAT(NEX)
      RNI = 1.0/RNEX
C
C---- solve quadratic for initial geometric ratio guess
      AAA = RNEX*(RNEX-1.0)*(RNEX-2.0) / 6.0
      BBB = RNEX*(RNEX-1.0) / 2.0
      CCC = RNEX - SIGMA
C
      DISC = BBB**2 - 4.0*AAA*CCC
      DISC = MAX( 0.0 , DISC )
C
      IF(NEX.LE.1) THEN
       STOP 'SETEXP: Cannot fill array.  N too small.'
      ELSE IF(NEX.EQ.2) THEN
       RATIO = -CCC/BBB  +  1.0
      ELSE
       RATIO = (-BBB + SQRT(DISC))/(2.0*AAA)  +  1.0
      ENDIF
C
      IF(RATIO.EQ.1.0) GO TO 11
C
C---- Newton iteration for actual geometric ratio
      DO 1 ITER=1, 100
        SIGMAN = (RATIO**NEX - 1.0) / (RATIO - 1.0)
        RES = SIGMAN**RNI - SIGMA**RNI
        DRESDR = RNI*SIGMAN**RNI
     &         * (RNEX*RATIO**(NEX-1) - SIGMAN) / (RATIO**NEX - 1.0)
C
        DRATIO = -RES/DRESDR
        RATIO = RATIO + DRATIO
C
        IF(ABS(DRATIO) .LT. 1.0E-5) GO TO 11
C
    1 CONTINUE
      WRITE(*,*) 'SETEXP: Convergence failed.  Continuing anyway ...'
C
C---- set up stretched array using converged geometric ratio
   11 S(1) = 0.0
      DS = DS1
      DO 2 N=2, NN
        S(N) = S(N-1) + DS
        DS = DS*RATIO
    2 CONTINUE
C
      RETURN
      END


      SUBROUTINE SETEX2(S,DS1,DSN,SMAX,N)
C==========================================================
C     Sets array S stretched so that a prescribed spacing is 
C     obtained at each end.  The interior spacing is a blend 
C     of two geometric stretchings "shot" from each end.
C
C       S     (output)  array to be set  
C       DS1   (input)   approximate first S increment:  S(2) - S(1)
C       DSN   (input)   approximate last  S increment:  S(N) - S(N-1)
C       SMAX  (input)   final S value:      S(N)
C       N     (input)   number of points
C==========================================================
C
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION S(N)
C
      PARAMETER (NDIM=500)
      DIMENSION S1(NDIM), SN(NDIM)
C
C---- with bigger FEND, the actual end increments will get closer 
C-    to DS1 & DSN, but the interior spacing might get screwy.
      DATA FEND / 2.0 /
C
      IF(N.GT.NDIM) STOP 'SETEX2:  Array overflow.'
C
C---- calculate spacing arrays each having the prescribed end increment
      CALL SETEXP(S1,DS1,SMAX,N)
      CALL SETEXP(SN,DSN,SMAX,N)
C
cc      S(1) = 0.
cc      S(N) = SMAX
C
C---- blend spacing arrays with power-function weights
      DO 10 I=1, N
        IN = N-I+1
        SS1 = S1(I)
        SSN = SMAX - SN(IN)
C
C------ power function of integer index
        WT1 = FLOAT(N-I)**FEND
        WTN = FLOAT(I-1)**FEND
C
C------ power function of coordinate
CCC     WT1 = (1.0 - SSN/SMAX)**FEND
CCC     WTN = (      SS1/SMAX)**FEND
C
        S(I) = (SS1*WT1 + SSN*WTN) / (WT1 + WTN)
   10 CONTINUE
C
C---- check for monotonicity
      SGN = SIGN( 1.0 , S(N)-S(1) )
      DO 20 I=2, N
        IF(SGN*S(I) .LE. SGN*S(I-1)) THEN
         WRITE(*,*) 'SETEX2: Warning. Returned array not monotonic.'
         RETURN
        ENDIF
   20 CONTINUE
C
      RETURN
      END ! SETEX2



      FUNCTION ATANC(Y,X,THOLD)
      IMPLICIT REAL (A-H,M,O-Z)
C---------------------------------------------------------------
C     ATAN2 function with branch cut checking.
C
C     Increments position angle of point X,Y from some previous
C     value THOLD due to a change in position, ensuring that the
C     position change does not cross the ATAN2 branch cut
C     (which is in the -x direction).  For example:
C
C       ATANC( -1.0 , -1.0 , 0.75*pi )  returns  1.25*pi , whereas
C       ATAN2( -1.0 , -1.0 )            returns  -.75*pi .
C
C     Typically, ATANC is used to fill an array of angles:
C
C        THETA(1) = ATAN2( Y(1) , X(1) )
C        DO i=2, N
C          THETA(i) = ATANC( Y(i) , X(i) , THETA(i-1) )
C        END DO
C
C     This will prevent the angle array THETA(i) from jumping by 
C     +/- 2 pi when the path X(i),Y(i) crosses the negative x axis.
C
C     Input:
C       X,Y     point position coordinates
C       THOLD   position angle of nearby point
C
C     Output:
C       ATANC   position angle of X,Y
C---------------------------------------------------------------
      DATA  PI /3.1415926535897932384/
      DATA TPI /6.2831853071795864769/
C
C---- set new position angle, ignoring branch cut in ATAN2 function for now
      THNEW = ATAN2( Y , X )
      DTHET = THNEW - THOLD
C
C---- angle change cannot exceed +/- pi, so get rid of any multiples of 2 pi 
      DTCORR = DTHET - TPI*INT( (DTHET + SIGN(PI,DTHET))/TPI )
C
C---- set correct new angle
      ATANC = THOLD + DTCORR
C
      RETURN
      END ! ATANC




      SUBROUTINE HSORT(N,A,INDX)
      DIMENSION A(N)
      DIMENSION INDX(N)
C--------------------------------------
C     Heapsort algorithm.
C     Returns INDX(.) such that
C
C       A(INDX(i)) < A(INDX(i+1))
C
C     Stolen from Numerical Recipes.
C--------------------------------------
C
      DO I = 1, N
        INDX(I) = I
      ENDDO
C
      IF(N.LE.1) RETURN
C
      L = N/2 + 1
      IR = N
C
 10   CONTINUE
      IF(L.GT.1) THEN
        L = L-1
        INDXT = INDX(L)
        Q = A(INDXT)
      ELSE
        INDXT = INDX(IR)
        Q = A(INDXT)
        INDX(IR) = INDX(1)
C
        IR = IR - 1
        IF(IR.EQ.1) THEN
          INDX(1) = INDXT
          RETURN
        ENDIF
      ENDIF
C
      I = L
      J = L+L
C
 20   IF(J.LE.IR) THEN
        IF(J.LT.IR) THEN
          IF(A(INDX(J)) .LT. A(INDX(J+1))) J = J+1
        ENDIF
        IF(Q .LT. A(INDX(J))) THEN
          INDX(I) = INDX(J)
C
          I = J
          J = J+J
        ELSE
          J = IR+1
        ENDIF
        GO TO 20
      ENDIF
C
      INDX(I) = INDXT
      GO TO 10
      END


 
