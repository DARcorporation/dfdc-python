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

      FUNCTION ISGFIND(SS,S,N)
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION S(N)
C
C---- Find closest value in S1 array to SS and return index
C     Assumes array S in ascending value order and uses binary search to
C     locate interval with S(I-1)<=SS<S(I)
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
C---- Check S(I-1) and S(I) for closest value to SS
 11   ISGFIND = I
      IF(I.GT.1) THEN
        IF(ABS(SS-S(I)).GT.ABS(SS-S(I-1))) THEN
         ISGFIND = I-1        
        ENDIF
      ENDIF
C
      RETURN
      END


      SUBROUTINE SGCOPY(S1,S2,N)
C---- Copy N elements from S1 to S2
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION S1(N), S2(N)
C
      DO I=1, N
        S2(I) = S1(I)
      END DO
C
      RETURN
      END


      SUBROUTINE SGCOPF(S1,S2,N,SOFF1,SWT1,SOFF2,SWT2, FOFF,
     &                  F1,FX1,X1,N1, F2,FX2,X2,N2)
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION S1(N), S2(N)
      DIMENSION F1(N1),FX1(N1),X1(N1)
      DIMENSION F2(N2),FX2(N2),X2(N2)
C
C---- reset array S1 so that  f1(s1) = f2(s2) + foff ,
C-    with  s1 = S1*SWT1 + SOFF1  and  s2 = S2*SWT2 + SOFF2
C
      DELS1 = S1(N) - S1(1)
      DO 20 I=2, N-1
        FEV2 = SEVAL((S2(I)*SWT2+SOFF2),F2,FX2,X2,N2) + FOFF
C
        DO ITER=1, 10
          FEV1 = SEVAL((S1(I)*SWT1+SOFF1),F1,FX1,X1,N1)
          DEV1 = DEVAL((S1(I)*SWT1+SOFF1),F1,FX1,X1,N1)
          RES = FEV1 - FEV2
          RES_S1 = DEV1 * SWT1
          DS1 = -RES/RES_S1
          S1(I) = S1(I) + DS1
          IF(ABS(DS1/DELS1) .LT. 1.0E-5) GO TO 20
        END DO
        WRITE(*,*) 'SGCOPF: Convergence failed.  dS/Smax =', DS1/DELS1
C
   20 CONTINUE
C
      RETURN
      END


      SUBROUTINE SGAVG(S1,S2,N,C1)
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION S1(N), S2(N)
C
C---- average arrays S1 S2,
C-    preserving S1 spacing at right point and S2 spacing at left endpoint
      DO I=1, N
        F1 = FLOAT(I-1) * C1
        F2 = FLOAT(N-I)
        S1(I) = (S1(I)*F1 + S2(I)*F2) / (F1 + F2)
        S2(I) = S1(I)
      END DO
C
      RETURN
      END


      SUBROUTINE SGAVG1(S1,S2,N)
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION S1(N), S2(N)
C
C---- impose spacing of array S1 on S2 at right endpoint
      DO I=1, N
        F1 = FLOAT(I-1)
        F2 = FLOAT(N-I)
        S2(I) = (S1(I)*F1 + S2(I)*F2) / (F1 + F2)
      END DO
C
      RETURN
      END


      SUBROUTINE SGAVG2(S1,S2,N)
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION S1(N), S2(N)
C
C---- impose spacing of array S2 on S1 at left endpoint
      DO I=1, N
        F1 = FLOAT(I-1)
        F2 = FLOAT(N-I)
        S1(I) = (S1(I)*F1 + S2(I)*F2) / (F1 + F2)
      END DO
C
      RETURN
      END


      SUBROUTINE SGSHFT(SBEG,SEND,SOLD,SNEW,S1,S2,N)
C
C---- Shift S1 values between limits SBEG and SEND such that
C     value at SOLD is shifted to SNEW in S2.  S1 values
C     between SBEG and SOLD are scaled to correspond to 
C     the range SBEG to SNEW. S1 values between SOLD and
C     SEND are scaled to the range SNEW to SEND.
C
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION S1(N), S2(N)
C
      IF(SBEG.GE.SEND) THEN
        WRITE(*,*) 'SGSHFT: Bad input range SBEG,SEND ',SBEG,SEND
        STOP
      ENDIF
C
      IF(SBEG.GT.SOLD .OR. SEND.LT.SOLD) THEN
        WRITE(*,*) 'SGSHFT: Bad input range SBEG/SOLD/SEND ',
     &             SBEG,SOLD,SEND
        STOP
      ENDIF
C
      IF(SBEG.GT.SNEW .OR. SEND.LT.SNEW) THEN
        WRITE(*,*) 'SGSHFT: Bad input range SBEG/SNEW/SEND ',
     &             SBEG,SNEW,SEND
        STOP
      ENDIF
C
C---- Shift S1 values in SBEG-SEND range to new values such
C     that SOLD becomes SNEW (linearly interpolated from both
C     ends of range).
      FAC1 = (SNEW-SBEG)/(SOLD-SBEG)
      FAC2 = (SEND-SNEW)/(SEND-SOLD)
      DO I=1, N
        IF(S1(I).LE.SBEG) THEN
          S2(I) = S1(I)
        ELSEIF(S1(I).GT.SBEG .AND. S1(I).LE.SOLD) THEN
          S2(I) = SBEG + (S1(I)-SBEG)*FAC1      
        ELSEIF(S1(I).GT.SOLD .AND. S1(I).LT.SEND) THEN
          S2(I) = SEND - (SEND-S1(I))*FAC2      
        ELSE 
          S2(I) = S1(I)
        ENDIF
      END DO
C
      RETURN
      END



      SUBROUTINE SGRENUM(S,N1,N2)
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION S(2000)
C----------------------------------------------------------------------
C     Interpolates input array S(1:N1) into new number of points 1:N2.
C     The interpolation is a spline in the array index.
C     The first and last elements in S will remain the same.
C----------------------------------------------------------------------
      PARAMETER (NMAX=2000)
      DIMENSION X(NMAX), XI(NMAX), C(NMAX)
C
      IF(N1.GT.NMAX) STOP 'NEWNUM:  Array overflow'
C
C---- save old input array
      DO I=1, N1
        X(I) = S(I)
      END DO
C
C---- spline X(i)  (set up and solve tridiagonal system on the fly)
      C(1) = 0.5
      XI(1) = 1.5*(X(2)-X(1))
      DO I=2, N1-1
        C(I) = 1.0 / (4.0 - C(I-1))
        XI(I) = (3.0*(X(I+1) - X(I-1)) - XI(I-1)) * C(I)
      END DO
      I = N1
      XI(I) = (3.0*(X(I)-X(I-1)) - XI(I-1)) / (2.0 - C(I-1))
C
      DO I=N1-1, 1, -1
        XI(I) = XI(I) - C(I)*XI(I+1)
      END DO
C
C---- evaluate s(i) spline at new points
      DO J=1, N2
        FRAC = FLOAT(J-1)/FLOAT(N2-1)
        RI = 1.0 + FRAC*FLOAT(N1-1)
        I  = MIN( INT(RI) , N1-1 )
C
        T = RI - FLOAT(I)
        CX1 = XI(I)   - X(I+1) + X(I)
        CX2 = XI(I+1) - X(I+1) + X(I)
        S(J) = T*X(I+1) + (1.0-T)*X(I) + (T-T*T)*((1.0-T)*CX1 - T*CX2)
      END DO
C
C---- make sure new endpoints are exactly the same as old ones
      S(1)  = X(1)
      S(N2) = X(N1)
C
      RETURN
      END ! SGRENUM


