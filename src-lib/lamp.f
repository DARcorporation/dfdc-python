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

      SUBROUTINE LAMP( X1, R1, X2, R2,  XF, RF, 
     &                 UG1, VG1, UG2, VG2, 
     &                 US1, VS1, US2, VS2 )
C-----------------------------------------------------------------------
C
C     Computes the velocities at XF,RF induced by a vortex + source 
C     sheet "lampshade"  extending from X1,R1 to X2,R2.  
C
C      Vortex sheet density = gamma  (positive counterclockwise)
C      Source sheet density = sigma
C
C     Both densities are assumed to be linear in meridional 
C     arc length over the panel.
C
C
C  Input:
C  ------
C     X1,R1   x,r of one   lampshade edge
C     X2,R2   x,r of other lampshade edge
C     XF,RF   x,r of the field point
C
C  Output: 
C  -------
C     UG1,VG1   x,r velocities at XF,RF for unit gamma at X1,R1
C     UG2,VG2   x,r velocities at XF,RF for unit gamma at X2,R2
C     US1,VS1   x,r velocities at XF,RF for unit sigma at X1,R1
C     US2,VS2   x,r velocities at XF,RF for unit sigma at X2,R2
C
C     Total x,r velocities for given endpoint sheet densities are:
C
C       U = UG1*gamma1 + UG2*gamma2 + US1*sigma1 + US2*sigma2
C       V = VG1*gamma1 + VG2*gamma2 + VS1*sigma1 + VS2*sigma2
C
C-----------------------------------------------------------------------
C
C     Integrates point vortex/source ring velocities over 
C     the lampshade using a Romberg sequence applied to the 
C     simple midpoint rule:
C
C        |       *       |  ->  I1  -   I21  -   I321  -   I4321
C                                   /        /         /     .
C        |   *       *   |  ->  I2  -   I32  -   I432        .
C                                   /        /     .       8th order
C        | *   *   *   * |  ->  I3  -   I43        .       
C                                   /    .       6th order
C        |* * * * * * * *|  ->  I4       .        
C                .               .      4th order
C                .               .      
C                               2nd order
C               etc             
C
C  The first column I1,I2... are the 2nd-order integral approximations
C  computed with the stock midpoint rule.  The subsequent columns
C  are the Richardson extrapolations to higher-order accuracy.
C
C
C  Algorithm for four Romberg stages:
C
C    I21 = (4*I2 - I1) / 3           | extrapolation of 2nd-order results
C    I32 = (4*I3 - I2) / 3           |
C    I43 = (4*I4 - I3) / 3           |
C 
C    I321 = (16*I32 - I21) / 15      | extrapolation of 4th-order results
C    I432 = (16*I43 - I32) / 15      |
C
C    I4321 = (64*I432 - I321) / 63   | extrapolation of 6th-order results
C
C
C  Quantities stored in the NROMX arrays after each IROM stage:
C
C    IROM  =   1      2       3       4     
C            ----   ----    -----   ------
C    U(1)  =  I1     I21     I321    I4321  
C    U(2)  =         I2      I32     I432   
C    U(3)  =                 I3      I43    
C    U(4)  =                         I4      
C
C-----------------------------------------------------------------------
C
C---- NROMX = max number of Romberg stages
C-      This limits the integration resolution and cost.
C-      The influence of this panel on a field point 
C-      which is closer than  O(panel_length)/2**NROMX 
C-      will not be accurately represented.
      PARAMETER (NROMX=10)
C
      DIMENSION UG1I(NROMX), VG1I(NROMX),
     &          UG2I(NROMX), VG2I(NROMX),
     &          US1I(NROMX), VS1I(NROMX),
     &          US2I(NROMX), VS2I(NROMX)
C
      PARAMETER(PI=3.14159265358979)
C
C
C---- Romberg convergence tolerance
      DATA ROMTOL / 1.0E-6 /
ccc      DATA ROMTOL / 1.0E-12 /
C
C---- reference length for convergence tolerance 
      REFL = 0.5*(R1 + R2)
C
C---- evaluate integrals over  0..t..1  on increasingly fine grids
      DO 100 IROM=1, NROMX
        NT = 2**IROM / 2
C
        UG1I(IROM) = 0.
        VG1I(IROM) = 0.
        UG2I(IROM) = 0.
        VG2I(IROM) = 0.
        US1I(IROM) = 0.
        VS1I(IROM) = 0.
        US2I(IROM) = 0.
        VS2I(IROM) = 0.
C
C------ visit the midpoints of each of the NT intervals
        DO 10 IT=1, NT
          T = (FLOAT(IT) - 0.5) / FLOAT(NT)
          TB = 1.0 - T
C
          DT = 1.0/FLOAT(NT)
C
          XT = X1*TB + X2*T
          RT = R1*TB + R2*T
C
C-------- get induced velocities for vortex,source ring at XT,RT
          CALL RING( XT, RT, XF, RF, UGT, VGT, UST, VST )
C
C-------- accumulate the separate unit-gamma, unit-sigma integrals
          UG1I(IROM) = UG1I(IROM) + DT*UGT*TB
          VG1I(IROM) = VG1I(IROM) + DT*VGT*TB
          UG2I(IROM) = UG2I(IROM) + DT*UGT*T
          VG2I(IROM) = VG2I(IROM) + DT*VGT*T
          US1I(IROM) = US1I(IROM) + DT*UST*TB
          VS1I(IROM) = VS1I(IROM) + DT*VST*TB
          US2I(IROM) = US2I(IROM) + DT*UST*T
          VS2I(IROM) = VS2I(IROM) + DT*VST*T
 10     CONTINUE
C
C------ Romberg sequence using all previous grid results
        DO KROM = IROM, 2, -1
C-------- weight needed to cancel lowest-order error terms in KROM level
          W = 2.0 ** (2*(IROM-KROM+1))
C
C-------- put Richardson extrapolation for KROM level into KROM-1 level
          UG1I(KROM-1) = (W*UG1I(KROM) - UG1I(KROM-1)) / (W-1.0)
          VG1I(KROM-1) = (W*VG1I(KROM) - VG1I(KROM-1)) / (W-1.0)
          UG2I(KROM-1) = (W*UG2I(KROM) - UG2I(KROM-1)) / (W-1.0)
          VG2I(KROM-1) = (W*VG2I(KROM) - VG2I(KROM-1)) / (W-1.0)
          US1I(KROM-1) = (W*US1I(KROM) - US1I(KROM-1)) / (W-1.0)
          VS1I(KROM-1) = (W*VS1I(KROM) - VS1I(KROM-1)) / (W-1.0)
          US2I(KROM-1) = (W*US2I(KROM) - US2I(KROM-1)) / (W-1.0)
          VS2I(KROM-1) = (W*VS2I(KROM) - VS2I(KROM-1)) / (W-1.0)
        ENDDO
C
        IF(IROM.GT.1) THEN
C------- compare the best-current and best-previous integrals
         ERRUG1 = UG1I(1) - UG1I(2)
         ERRVG1 = VG1I(1) - VG1I(2)
         ERRUG2 = UG2I(1) - UG2I(2)
         ERRVG2 = VG2I(1) - VG2I(2)
         ERRUS1 = US1I(1) - US1I(2)
         ERRVS1 = VS1I(1) - VS1I(2)
         ERRUS2 = US2I(1) - US2I(2)
         ERRVS2 = VS2I(1) - VS2I(2)
C
         ERR = MAX(
     &     ABS(ERRUG1),
     &     ABS(ERRVG1),
     &     ABS(ERRUG2),
     &     ABS(ERRVG2),
     &     ABS(ERRUS1),
     &     ABS(ERRVS1),
     &     ABS(ERRUS2),
     &     ABS(ERRVS2) )
         
c         write(13,1200) 
c     &     ABS(ERRUG1),
c     &     ABS(ERRVG1),
c     &     ABS(ERRUG2),
c     &     ABS(ERRVG2),
c     &     ABS(ERRUS1),
c     &     ABS(ERRVS1),
c     &     ABS(ERRUS2),
c     &     ABS(ERRVS2)
c 1200    format(8e16.9)

         IF(ERR*REFL .LT. ROMTOL) GO TO 101
        ENDIF
 100  CONTINUE
      WRITE(*,*) 'LAMP: Romberg convergence failed.  Error =', ERR
C
 101  CONTINUE

ccc      write(*,*) IROM, ERR
C
C---- return best final results
      DELSQ = (X1-X2)**2 + (R1-R2)**2
      DELS = SQRT(DELSQ)
C
      UG1 = UG1I(1)*DELS
      VG1 = VG1I(1)*DELS
      UG2 = UG2I(1)*DELS
      VG2 = VG2I(1)*DELS
      US1 = US1I(1)*DELS
      VS1 = VS1I(1)*DELS
      US2 = US2I(1)*DELS
      VS2 = VS2I(1)*DELS
C
      RETURN
      END ! LAMP



      SUBROUTINE LAMPC( X1, R1, X2, R2,
     &                  UG1, VG1, UG2, VG2, 
     &                  US1, VS1, US2, VS2 )
C---------------------------------------------------------
C     Same as LAMP, but the field point is assumed 
C     to be at the lampshade-panel midpoint.
C
C     The 1/r and log(r) singularities in the integrands 
C     are removed from the numerical integration,
C     and are computed analytically.
C
C     The induced velocites returned by this routine 
C     are the average of the two values on each side 
C     of the sheet.  The sheet jumps are not included.
C---------------------------------------------------------
C
C---- max number of Romberg integration stages
      PARAMETER (NROMX=9)
C
      DIMENSION UG1I(NROMX), VG1I(NROMX),
     &          UG2I(NROMX), VG2I(NROMX),
     &          US1I(NROMX), VS1I(NROMX),
     &          US2I(NROMX), VS2I(NROMX)
C
      PARAMETER(PI=3.14159265358979)
C
C
C---- Romberg convergence tolerance (actual error may be much less than this)
      DATA ROMTOL / 1.0E-6 /
ccc      DATA ROMTOL / 1.0E-12 /
C
C---- lampshade meridional length**2
      DELSQ = (X1-X2)**2 + (R1-R2)**2
C
C---- reference length for convergence tolerance 
      REFL = 0.5*(R1 + R2)
C
C
C---- field point is assumed to be at midpoint
      XF = 0.5*(X1 + X2)
      RF = 0.5*(R1 + R2)
C
C---- evaluate integrals on increasingly fine grids 
C-     (start with two intervals to avoid landing right on the midpoint)
      DO 100 IROM=1, NROMX
        NT = 2**IROM
C
        UG1I(IROM) = 0.
        VG1I(IROM) = 0.
        UG2I(IROM) = 0.
        VG2I(IROM) = 0.
        US1I(IROM) = 0.
        VS1I(IROM) = 0.
        US2I(IROM) = 0.
        VS2I(IROM) = 0.
C
C------ visit the midpoints of each of the NT intervals
        DO 10 IT=1, NT
          T = (FLOAT(IT) - 0.5) / FLOAT(NT)
          TB = 1.0 - T
C
          DT = 1.0/FLOAT(NT)
C
          XT = X1*TB + X2*T
          RT = R1*TB + R2*T
C
          CALL RING( XT, RT, XF, RF, UGT, VGT, UST, VST )
C
C-------- singular parts of velocities in the limit  XT,RT -> XF,RF
          DSQ = (XT-XF)**2 + (RT-RF)**2
          UGA =  (RT-RF)/(4.0*PI*DSQ)
     &        - 0.5*LOG(DSQ/(64.0*RF**2)) / (8.0*PI*RF)
          VGA = -(XT-XF)/(4.0*PI*DSQ)
          USA = -(XT-XF)/(4.0*PI*DSQ)
          VSA = -(RT-RF)/(4.0*PI*DSQ)
     &        - 0.5*LOG(DSQ/      RF**2 ) / (8.0*PI*RF)
C
C-------- accumulate integrals, with singular parts (at t=0.5) removed
          UG1I(IROM) = UG1I(IROM) + DT*(UGT*TB - UGA)
          VG1I(IROM) = VG1I(IROM) + DT*(VGT*TB - VGA)
          UG2I(IROM) = UG2I(IROM) + DT*(UGT*T  - UGA)
          VG2I(IROM) = VG2I(IROM) + DT*(VGT*T  - VGA)
          US1I(IROM) = US1I(IROM) + DT*(UST*TB - USA)
          VS1I(IROM) = VS1I(IROM) + DT*(VST*TB - VSA)
          US2I(IROM) = US2I(IROM) + DT*(UST*T  - USA)
          VS2I(IROM) = VS2I(IROM) + DT*(VST*T  - VSA)
 10     CONTINUE
C
C------ Romberg sequence using all previous grid results
        DO KROM = IROM, 2, -1
C-------- weight needed to cancel lowest-order error terms in KROM level
          W = 2.0**(2*(IROM-KROM+1))
          WM1 = W - 1.0
C
C-------- put Richardson extrapolation for KROM level into KROM-1 level
          UG1I(KROM-1) = (W*UG1I(KROM) - UG1I(KROM-1)) / WM1
          VG1I(KROM-1) = (W*VG1I(KROM) - VG1I(KROM-1)) / WM1
          UG2I(KROM-1) = (W*UG2I(KROM) - UG2I(KROM-1)) / WM1
          VG2I(KROM-1) = (W*VG2I(KROM) - VG2I(KROM-1)) / WM1
          US1I(KROM-1) = (W*US1I(KROM) - US1I(KROM-1)) / WM1
          VS1I(KROM-1) = (W*VS1I(KROM) - VS1I(KROM-1)) / WM1
          US2I(KROM-1) = (W*US2I(KROM) - US2I(KROM-1)) / WM1
          VS2I(KROM-1) = (W*VS2I(KROM) - VS2I(KROM-1)) / WM1
        ENDDO
C
        IF(IROM.GT.1) THEN
C------- compare the best-current and best-previous integrals
         ERRUG1 = UG1I(1) - UG1I(2)
         ERRVG1 = VG1I(1) - VG1I(2)
         ERRUG2 = UG2I(1) - UG2I(2)
         ERRVG2 = VG2I(1) - VG2I(2)
         ERRUS1 = US1I(1) - US1I(2)
         ERRVS1 = VS1I(1) - VS1I(2)
         ERRUS2 = US2I(1) - US2I(2)
         ERRVS2 = VS2I(1) - VS2I(2)
C
         ERR = MAX(
     &     ABS(ERRUG1),
     &     ABS(ERRVG1),
     &     ABS(ERRUG2),
     &     ABS(ERRVG2),
     &     ABS(ERRUS1),
     &     ABS(ERRVS1),
     &     ABS(ERRUS2),
     &     ABS(ERRVS2) )
C
         IF(ERR*REFL .LT. ROMTOL) GO TO 101
        ENDIF
 100  CONTINUE
      WRITE(*,*) 'LAMPC: Romberg convergence failed.  Error =', ERR
C
 101  CONTINUE
C
      DELSQ = (X1-X2)**2 + (R1-R2)**2
      DELS = SQRT(DELSQ)
C
C---- analytically-integrated singular parts which were removed
      UGAI = (1.0 + LOG(16.0*RF/DELS)) / (4.0*PI*RF)
      VGAI = 0.
      USAI = 0.
      VSAI = (1.0 + LOG( 2.0*RF/DELS)) / (4.0*PI*RF)
C
C---- return final results, with removed parts added back on
      UG1 = (UG1I(1) + UGAI*0.5)*DELS
      VG1 = (VG1I(1) + VGAI*0.5)*DELS
      UG2 = (UG2I(1) + UGAI*0.5)*DELS
      VG2 = (VG2I(1) + VGAI*0.5)*DELS
      US1 = (US1I(1) + USAI*0.5)*DELS
      VS1 = (VS1I(1) + VSAI*0.5)*DELS
      US2 = (US2I(1) + USAI*0.5)*DELS
      VS2 = (VS2I(1) + VSAI*0.5)*DELS
C
      RETURN
      END ! LAMPC



      SUBROUTINE GLAMP( R1,     R2,     RF, 
     &         QG1, QG1_RF, 
     &         QG2, QG2_RF, 
     &         QS1, QS1_RF, 
     &         QS2, QS2_RF )
C-----------------------------------------------------------------------
C     Same as LAMP, but also returns velocity gradient
C-----------------------------------------------------------------------
      DIMENSION R1(2), R2(2), RF(2)
      DIMENSION QG1(2), QG1_RF(2,2),
     &          QG2(2), QG2_RF(2,2),
     &          QS1(2), QS1_RF(2,2),
     &          QS2(2), QS2_RF(2,2)
C
      PARAMETER (NROMX=10)
C
      DIMENSION QG1I(2,NROMX),
     &          QG2I(2,NROMX),
     &          QS1I(2,NROMX),
     &          QS2I(2,NROMX)
      DIMENSION 
     &  QG1I_RF(2,2,NROMX),
     &  QG2I_RF(2,2,NROMX),
     &  QS1I_RF(2,2,NROMX),
     &  QS2I_RF(2,2,NROMX) 
C
      DIMENSION RT(2)
      DIMENSION QGT(2), QGT_RT(2,2), QGT_RF(2,2),
     &          QST(2), QST_RT(2,2), QST_RF(2,2)
C
      PARAMETER(PI=3.14159265358979)
C
C---- Romberg convergence tolerance
      DATA ROMTOL / 1.0E-6 /
ccc   DATA ROMTOL / 1.0E-12 /
C
C---- reference length for convergence tolerance 
      REFL = 0.5*(R1(2) + R2(2))
C
C---- evaluate integrals over  0..t..1  on increasingly fine grids
      DO 100 IROM=1, NROMX
        NT = 2**IROM / 2
C
        DO K = 1, 2
          QG1I(K,IROM) = 0.
          QG2I(K,IROM) = 0.
          QS1I(K,IROM) = 0.
          QS2I(K,IROM) = 0.
          DO J = 1, 2
            QG1I_RF(K,J,IROM) = 0.
            QG2I_RF(K,J,IROM) = 0.
            QS1I_RF(K,J,IROM) = 0.
            QS2I_RF(K,J,IROM) = 0.
          ENDDO
        ENDDO
C
C------ visit the midpoints of each of the NT intervals
        DO 10 IT=1, NT
          T = (FLOAT(IT) - 0.5) / FLOAT(NT)
          TB = 1.0 - T
C
          DT = 1.0/FLOAT(NT)
C
          RT(1) = R1(1)*TB + R2(1)*T
          RT(2) = R1(2)*TB + R2(2)*T
C
C-------- get induced velocities for vortex,source ring at XT,RT
          CALL DRING( RT(  1),    RT(  2),    RF(  1),    RF(  2),
     &     QGT(1),QGT_RT(1,1),QGT_RT(1,2),QGT_RF(1,1),QGT_RF(1,2), 
     &     QGT(2),QGT_RT(2,1),QGT_RT(2,2),QGT_RF(2,1),QGT_RF(2,2),
     &     QST(1),QST_RT(1,1),QST_RT(1,2),QST_RF(1,1),QST_RF(1,2),
     &     QST(2),QST_RT(2,1),QST_RT(2,2),QST_RF(2,1),QST_RF(2,2) )
C
C-------- accumulate the separate unit-gamma, unit-sigma integrals
          DO K = 1, 2
           QG1I(K,IROM) = QG1I(K,IROM) + DT*QGT(K)*TB
           QG2I(K,IROM) = QG2I(K,IROM) + DT*QGT(K)*T
           QS1I(K,IROM) = QS1I(K,IROM) + DT*QST(K)*TB
           QS2I(K,IROM) = QS2I(K,IROM) + DT*QST(K)*T
           DO J = 1, 2
            QG1I_RF(K,J,IROM) = QG1I_RF(K,J,IROM) + DT*QGT_RF(K,J)*TB
            QG2I_RF(K,J,IROM) = QG2I_RF(K,J,IROM) + DT*QGT_RF(K,J)*T     
            QS1I_RF(K,J,IROM) = QS1I_RF(K,J,IROM) + DT*QST_RF(K,J)*TB
            QS2I_RF(K,J,IROM) = QS2I_RF(K,J,IROM) + DT*QST_RF(K,J)*T     
           ENDDO
          ENDDO
 10     CONTINUE
C
C------ Romberg sequence using all previous grid results
        DO KROM = IROM, 2, -1
C-------- weight needed to cancel lowest-order error terms in KROM level
          W = 2.0 ** (2*(IROM-KROM+1))
          WM1 = W - 1.0
C
C-------- put Richardson extrapolation for KROM level into KROM-1 level
          DO K = 1, 2
            QG1I(K,KROM-1) = (W*QG1I(K,KROM) - QG1I(K,KROM-1)) / WM1
            QG2I(K,KROM-1) = (W*QG2I(K,KROM) - QG2I(K,KROM-1)) / WM1
            QS1I(K,KROM-1) = (W*QS1I(K,KROM) - QS1I(K,KROM-1)) / WM1
            QS2I(K,KROM-1) = (W*QS2I(K,KROM) - QS2I(K,KROM-1)) / WM1
            DO J = 1, 2
              QG1I_RF(K,J,KROM-1) = (W*QG1I_RF(K,J,KROM)
     &                               - QG1I_RF(K,J,KROM-1)) / WM1
              QG2I_RF(K,J,KROM-1) = (W*QG2I_RF(K,J,KROM)
     &                               - QG2I_RF(K,J,KROM-1)) / WM1
              QS1I_RF(K,J,KROM-1) = (W*QS1I_RF(K,J,KROM)
     &                               - QS1I_RF(K,J,KROM-1)) / WM1
              QS2I_RF(K,J,KROM-1) = (W*QS2I_RF(K,J,KROM)
     &                               - QS2I_RF(K,J,KROM-1)) / WM1
            ENDDO
          ENDDO
        ENDDO
C
        IF(IROM.GT.1) THEN
C------- compare the best-current and best-previous integrals
         ERRUG1 = QG1I(1,1) - QG1I(1,2)
         ERRVG1 = QG1I(2,1) - QG1I(2,2)
         ERRUG2 = QG2I(1,1) - QG2I(1,2)
         ERRVG2 = QG2I(2,1) - QG2I(2,2)
         ERRUS1 = QS1I(1,1) - QS1I(1,2)
         ERRVS1 = QS1I(2,1) - QS1I(2,2)
         ERRUS2 = QS2I(1,1) - QS2I(1,2)
         ERRVS2 = QS2I(2,1) - QS2I(2,2)
C
         ERR = MAX(
     &     ABS(ERRUG1),
     &     ABS(ERRVG1),
     &     ABS(ERRUG2),
     &     ABS(ERRVG2),
     &     ABS(ERRUS1),
     &     ABS(ERRVS1),
     &     ABS(ERRUS2),
     &     ABS(ERRVS2) )
         
c         write(13,1200) 
c     &     ABS(ERRUG1),
c     &     ABS(ERRVG1),
c     &     ABS(ERRUG2),
c     &     ABS(ERRVG2),
c     &     ABS(ERRUS1),
c     &     ABS(ERRVS1),
c     &     ABS(ERRUS2),
c     &     ABS(ERRVS2)
c 1200    format(8e16.9)

         IF(ERR*REFL .LT. ROMTOL) GO TO 101
        ENDIF
 100  CONTINUE
      WRITE(*,*) 'GLAMP: Romberg convergence failed.  Error =', ERR
C
 101  CONTINUE

ccc      write(*,*) IROM, ERR
C
C---- return best final results
      DELSQ = (R1(1)-R2(1))**2 + (R1(2)-R2(2))**2
      DELS  = SQRT(DELSQ)
C
      DO K = 1, 2
        QG1(K) = QG1I(K,1)*DELS
        QG2(K) = QG2I(K,1)*DELS
        QS1(K) = QS1I(K,1)*DELS
        QS2(K) = QS2I(K,1)*DELS
        DO J = 1, 2
          QG1_RF(K,J) = QG1I_RF(K,J,1)*DELS
          QG2_RF(K,J) = QG2I_RF(K,J,1)*DELS
          QS1_RF(K,J) = QS1I_RF(K,J,1)*DELS
          QS2_RF(K,J) = QS2I_RF(K,J,1)*DELS
        ENDDO
      ENDDO
C
      RETURN
      END ! GLAMP


      SUBROUTINE GLAMPC( R1, R2,
     &         QG1, QG1_RF, 
     &         QG2, QG2_RF, 
     &         QS1, QS1_RF, 
     &         QS2, QS2_RF )
C-----------------------------------------------------------------------
C     Same as LAMPC, but also returns velocity gradient
C-----------------------------------------------------------------------
      DIMENSION R1(2), R2(2)
      DIMENSION QG1(2), QG1_RF(2,2),
     &          QG2(2), QG2_RF(2,2),
     &          QS1(2), QS1_RF(2,2),
     &          QS2(2), QS2_RF(2,2)
C
      PARAMETER (NROMX=10)
C
      DIMENSION QG1I(2,NROMX),
     &          QG2I(2,NROMX),
     &          QS1I(2,NROMX),
     &          QS2I(2,NROMX)
      DIMENSION 
     &  QG1I_RF(2,2,NROMX),
     &  QG2I_RF(2,2,NROMX),
     &  QS1I_RF(2,2,NROMX),
     &  QS2I_RF(2,2,NROMX) 
C
      DIMENSION RF(2)
      DIMENSION RT(2)
      DIMENSION QGT(2), QGT_RT(2,2), QGT_RF(2,2),
     &          QST(2), QST_RT(2,2), QST_RF(2,2)
      DIMENSION QGA(2), QGA_RT(2,2), QGA_RF(2,2),
     &          QSA(2), QSA_RT(2,2), QSA_RF(2,2)
      DIMENSION DSQ_RF(2)
C
      DIMENSION QGAI(2), QGAI_RF(2,2),
     &          QSAI(2), QSAI_RF(2,2)
C
      PARAMETER(PI=3.14159265358979)
C
C---- Romberg convergence tolerance
      DATA ROMTOL / 1.0E-6 /
ccc   DATA ROMTOL / 1.0E-12 /
C
C---- reference length for convergence tolerance 
      REFL = 0.5*(R1(2) + R2(2))
C
      RF(1) = 0.5*(R1(1) + R2(1))
      RF(2) = 0.5*(R1(2) + R2(2))
C
C---- evaluate integrals over  0..t..1  on increasingly fine grids
      DO 100 IROM=1, NROMX
        NT = 2**IROM
C
        DO K = 1, 2
          QG1I(K,IROM) = 0.
          QG2I(K,IROM) = 0.
          QS1I(K,IROM) = 0.
          QS2I(K,IROM) = 0.
          DO J = 1, 2
            QG1I_RF(K,J,IROM) = 0.
            QG2I_RF(K,J,IROM) = 0.
            QS1I_RF(K,J,IROM) = 0.
            QS2I_RF(K,J,IROM) = 0.
          ENDDO
        ENDDO
C
C------ visit the midpoints of each of the NT intervals
        DO 10 IT=1, NT
          T = (FLOAT(IT) - 0.5) / FLOAT(NT)
          TB = 1.0 - T
C
          DT = 1.0/FLOAT(NT)
C
          RT(1) = R1(1)*TB + R2(1)*T
          RT(2) = R1(2)*TB + R2(2)*T
C
C-------- get induced velocities for vortex,source ring at XT,RT
          CALL DRING( RT(  1),    RT(  2),    RF(  1),    RF(  2),
     &     QGT(1),QGT_RT(1,1),QGT_RT(1,2),QGT_RF(1,1),QGT_RF(1,2), 
     &     QGT(2),QGT_RT(2,1),QGT_RT(2,2),QGT_RF(2,1),QGT_RF(2,2),
     &     QST(1),QST_RT(1,1),QST_RT(1,2),QST_RF(1,1),QST_RF(1,2),
     &     QST(2),QST_RT(2,1),QST_RT(2,2),QST_RF(2,1),QST_RF(2,2) )
C
C-------- singular parts of velocities in the limit  XT,RT -> XF,RF
          DSQ = (RT(1)-RF(1))**2 + (RT(2)-RF(2))**2
          DSQ_RF(1) = -2.0*(RT(1)-RF(1))
          DSQ_RF(2) = -2.0*(RT(2)-RF(2))
C
C
          PID4 = 4.0*PI*DSQ
          PIR8 = 8.0*PI*RF(2)
          DRSQ = DSQ/RF(2)**2
C
          DR1 = (RT(1)-RF(1))/PID4
          DR2 = (RT(2)-RF(2))/PID4
C
          QGA(1) =  DR2 - 0.5*LOG(DRSQ/64.0) / PIR8
          QGA(2) = -DR1
          QSA(1) = -DR1
          QSA(2) = -DR2 - 0.5*LOG(DRSQ     ) / PIR8
C
          QGA_RF(1,1) = (-DR2 - 0.5/PIR8)/DSQ * DSQ_RF(1)
          QGA_RF(1,2) = (-DR2 - 0.5/PIR8)/DSQ * DSQ_RF(2) - 1.0/PID4
     &             + (1.0 + 0.5*LOG(DRSQ/64.0)) / (RF(2)*PIR8)
C
          QGA_RF(2,1) =   DR1            /DSQ * DSQ_RF(1) + 1.0/PID4
          QGA_RF(2,2) =   DR1            /DSQ * DSQ_RF(2)
C
          QSA_RF(1,1) =   DR1            /DSQ * DSQ_RF(1) + 1.0/PID4
          QSA_RF(1,2) =   DR1            /DSQ * DSQ_RF(2)
C
          QSA_RF(2,1) = ( DR2 - 0.5/PIR8)/DSQ * DSQ_RF(1)
          QSA_RF(2,2) = ( DR2 - 0.5/PIR8)/DSQ * DSQ_RF(2) + 1.0/PID4
     &             + (1.0 + 0.5*LOG(DRSQ     )) / (RF(2)*PIR8)
C
C-------- accumulate integrals, with singular parts (at t=0.5) removed
          DO K = 1, 2
           QG1I(K,IROM) = QG1I(K,IROM) + DT*(QGT(K)*TB - QGA(K))
           QG2I(K,IROM) = QG2I(K,IROM) + DT*(QGT(K)*T  - QGA(K))
           QS1I(K,IROM) = QS1I(K,IROM) + DT*(QST(K)*TB - QSA(K))
           QS2I(K,IROM) = QS2I(K,IROM) + DT*(QST(K)*T  - QSA(K))
           DO J = 1, 2
            QG1I_RF(K,J,IROM) =
     &      QG1I_RF(K,J,IROM) + DT*(QGT_RF(K,J)*TB - QGA_RF(K,J))
            QG2I_RF(K,J,IROM) =
     &      QG2I_RF(K,J,IROM) + DT*(QGT_RF(K,J)*T  - QGA_RF(K,J))
            QS1I_RF(K,J,IROM) =
     &      QS1I_RF(K,J,IROM) + DT*(QST_RF(K,J)*TB - QSA_RF(K,J))
            QS2I_RF(K,J,IROM) =
     &      QS2I_RF(K,J,IROM) + DT*(QST_RF(K,J)*T  - QSA_RF(K,J))        
           ENDDO
          ENDDO
 10     CONTINUE
C
C------ Romberg sequence using all previous grid results
        DO KROM = IROM, 2, -1
C-------- weight needed to cancel lowest-order error terms in KROM level
          W = 2.0 ** (2*(IROM-KROM+1))
          WM1 = W - 1.0
C
C-------- put Richardson extrapolation for KROM level into KROM-1 level
          DO K = 1, 2
            QG1I(K,KROM-1) = (W*QG1I(K,KROM) - QG1I(K,KROM-1)) / WM1
            QG2I(K,KROM-1) = (W*QG2I(K,KROM) - QG2I(K,KROM-1)) / WM1
            QS1I(K,KROM-1) = (W*QS1I(K,KROM) - QS1I(K,KROM-1)) / WM1
            QS2I(K,KROM-1) = (W*QS2I(K,KROM) - QS2I(K,KROM-1)) / WM1
            DO J = 1, 2
              QG1I_RF(K,J,KROM-1) = (W*QG1I_RF(K,J,KROM)
     &                               - QG1I_RF(K,J,KROM-1)) / WM1
              QG2I_RF(K,J,KROM-1) = (W*QG2I_RF(K,J,KROM)
     &                               - QG2I_RF(K,J,KROM-1)) / WM1
              QS1I_RF(K,J,KROM-1) = (W*QS1I_RF(K,J,KROM)
     &                               - QS1I_RF(K,J,KROM-1)) / WM1
              QS2I_RF(K,J,KROM-1) = (W*QS2I_RF(K,J,KROM)
     &                               - QS2I_RF(K,J,KROM-1)) / WM1
            ENDDO
          ENDDO
        ENDDO
C
        IF(IROM.GT.1) THEN
C------- compare the best-current and best-previous integrals
         ERRUG1 = QG1I(1,1) - QG1I(1,2)
         ERRVG1 = QG1I(2,1) - QG1I(2,2)
         ERRUG2 = QG2I(1,1) - QG2I(1,2)
         ERRVG2 = QG2I(2,1) - QG2I(2,2)
         ERRUS1 = QS1I(1,1) - QS1I(1,2)
         ERRVS1 = QS1I(2,1) - QS1I(2,2)
         ERRUS2 = QS2I(1,1) - QS2I(1,2)
         ERRVS2 = QS2I(2,1) - QS2I(2,2)
C
         ERR = MAX(
     &     ABS(ERRUG1),
     &     ABS(ERRVG1),
     &     ABS(ERRUG2),
     &     ABS(ERRVG2),
     &     ABS(ERRUS1),
     &     ABS(ERRVS1),
     &     ABS(ERRUS2),
     &     ABS(ERRVS2) )
         
c         write(13,1200) 
c     &     ABS(ERRUG1),
c     &     ABS(ERRVG1),
c     &     ABS(ERRUG2),
c     &     ABS(ERRVG2),
c     &     ABS(ERRUS1),
c     &     ABS(ERRVS1),
c     &     ABS(ERRUS2),
c     &     ABS(ERRVS2)
c 1200    format(8e16.9)

         IF(ERR*REFL .LT. ROMTOL) GO TO 101
        ENDIF
 100  CONTINUE
      WRITE(*,*) 'GLAMPC: Romberg convergence failed.  Error =', ERR
C
 101  CONTINUE

ccc      write(*,*) IROM, ERR
C
C---- return best final results
      DELSQ = (R1(1)-R2(1))**2 + (R1(2)-R2(2))**2
      DELS  = SQRT(DELSQ)
C
C---- analytically-integrated singular parts which were removed
      QGAI(1) = (1.0 + LOG(16.0*RF(2)/DELS)) / (8.0*PI*RF(2))
      QGAI(2) = 0.
      QSAI(1) = 0.
      QSAI(2) = (1.0 + LOG( 2.0*RF(2)/DELS)) / (8.0*PI*RF(2))
C
      QGAI_RF(1,1) = 0.
      QGAI_RF(2,1) = 0.
      QSAI_RF(1,1) = 0.
      QSAI_RF(2,1) = 0.
C
      QGAI_RF(1,2) = (1.0/RF(2)) / (8.0*PI*RF(2)) - QGAI(1)/RF(2)
      QGAI_RF(2,2) = 0.
      QSAI_RF(1,2) = 0.
      QSAI_RF(2,2) = (1.0/RF(2)) / (8.0*PI*RF(2)) - QSAI(2)/RF(2)
C
C
      DO K = 1, 2
        QG1(K) = (QG1I(K,1) + QGAI(K))*DELS
        QG2(K) = (QG2I(K,1) + QGAI(K))*DELS
        QS1(K) = (QS1I(K,1) + QSAI(K))*DELS
        QS2(K) = (QS2I(K,1) + QSAI(K))*DELS
        DO J = 1, 2
          QG1_RF(K,J) = (QG1I_RF(K,J,1) + QGAI_RF(K,J))*DELS
          QG2_RF(K,J) = (QG2I_RF(K,J,1) + QGAI_RF(K,J))*DELS
          QS1_RF(K,J) = (QS1I_RF(K,J,1) + QSAI_RF(K,J))*DELS
          QS2_RF(K,J) = (QS2I_RF(K,J,1) + QSAI_RF(K,J))*DELS
        ENDDO
      ENDDO
C
      RETURN
      END ! GLAMPC



      SUBROUTINE DLAMP( R1,     R2,     RF, 
     &         QG1, QG1_R1, QG1_R2, QG1_RF, 
     &         QG2, QG2_R1, QG2_R2, QG2_RF, 
     &         QS1, QS1_R1, QS1_R2, QS1_RF, 
     &         QS2, QS2_R1, QS2_R2, QS2_RF )
C-----------------------------------------------------------------------
C     Same as LAMP, but also returns derivatives.
C-----------------------------------------------------------------------
      DIMENSION R1(2), R2(2), RF(2)
      DIMENSION QG1(2), QG1_R1(2,2), QG1_R2(2,2), QG1_RF(2,2),
     &          QG2(2), QG2_R1(2,2), QG2_R2(2,2), QG2_RF(2,2),
     &          QS1(2), QS1_R1(2,2), QS1_R2(2,2), QS1_RF(2,2),
     &          QS2(2), QS2_R1(2,2), QS2_R2(2,2), QS2_RF(2,2)
C
      PARAMETER (NROMX=10)
C
      DIMENSION QG1I(2,NROMX),
     &          QG2I(2,NROMX),
     &          QS1I(2,NROMX),
     &          QS2I(2,NROMX)
      DIMENSION 
     &  QG1I_R1(2,2,NROMX),QG1I_R2(2,2,NROMX),QG1I_RF(2,2,NROMX),
     &  QG2I_R1(2,2,NROMX),QG2I_R2(2,2,NROMX),QG2I_RF(2,2,NROMX),
     &  QS1I_R1(2,2,NROMX),QS1I_R2(2,2,NROMX),QS1I_RF(2,2,NROMX),
     &  QS2I_R1(2,2,NROMX),QS2I_R2(2,2,NROMX),QS2I_RF(2,2,NROMX) 
C
      DIMENSION RT(2)
      DIMENSION QGT(2), QGT_RT(2,2), QGT_RF(2,2),
     &          QST(2), QST_RT(2,2), QST_RF(2,2)
      DIMENSION DELS_R1(2), DELS_R2(2)
C
      PARAMETER(PI=3.14159265358979)
C
C---- Romberg convergence tolerance
      DATA ROMTOL / 1.0E-6 /
ccc   DATA ROMTOL / 1.0E-12 /
C
C---- reference length for convergence tolerance 
      REFL = 0.5*(R1(2) + R2(2))
C
C---- evaluate integrals over  0..t..1  on increasingly fine grids
      DO 100 IROM=1, NROMX
        NT = 2**IROM / 2
C
        DO K = 1, 2
          QG1I(K,IROM) = 0.
          QG2I(K,IROM) = 0.
          QS1I(K,IROM) = 0.
          QS2I(K,IROM) = 0.
          DO J = 1, 2
            QG1I_R1(K,J,IROM) = 0.
            QG2I_R1(K,J,IROM) = 0.
            QS1I_R1(K,J,IROM) = 0.
            QS2I_R1(K,J,IROM) = 0.
            QG1I_R2(K,J,IROM) = 0.
            QG2I_R2(K,J,IROM) = 0.
            QS1I_R2(K,J,IROM) = 0.
            QS2I_R2(K,J,IROM) = 0.
            QG1I_RF(K,J,IROM) = 0.
            QG2I_RF(K,J,IROM) = 0.
            QS1I_RF(K,J,IROM) = 0.
            QS2I_RF(K,J,IROM) = 0.
          ENDDO
        ENDDO
C
C------ visit the midpoints of each of the NT intervals
        DO 10 IT=1, NT
          T = (FLOAT(IT) - 0.5) / FLOAT(NT)
          TB = 1.0 - T
C
          DT = 1.0/FLOAT(NT)
C
          RT(1) = R1(1)*TB + R2(1)*T
          RT(2) = R1(2)*TB + R2(2)*T
C
C-------- get induced velocities for vortex,source ring at XT,RT
          CALL DRING( RT(  1),    RT(  2),    RF(  1),    RF(  2),
     &     QGT(1),QGT_RT(1,1),QGT_RT(1,2),QGT_RF(1,1),QGT_RF(1,2), 
     &     QGT(2),QGT_RT(2,1),QGT_RT(2,2),QGT_RF(2,1),QGT_RF(2,2),
     &     QST(1),QST_RT(1,1),QST_RT(1,2),QST_RF(1,1),QST_RF(1,2),
     &     QST(2),QST_RT(2,1),QST_RT(2,2),QST_RF(2,1),QST_RF(2,2) )
C
C-------- accumulate the separate unit-gamma, unit-sigma integrals
          DO K = 1, 2
           QG1I(K,IROM) = QG1I(K,IROM) + DT*QGT(K)*TB
           QG2I(K,IROM) = QG2I(K,IROM) + DT*QGT(K)*T
           QS1I(K,IROM) = QS1I(K,IROM) + DT*QST(K)*TB
           QS2I(K,IROM) = QS2I(K,IROM) + DT*QST(K)*T
           DO J = 1, 2
            QG1I_R1(K,J,IROM) = QG1I_R1(K,J,IROM) + DT*QGT_RT(K,J)*TB*TB
            QG2I_R1(K,J,IROM) = QG2I_R1(K,J,IROM) + DT*QGT_RT(K,J)*T *TB 
            QS1I_R1(K,J,IROM) = QS1I_R1(K,J,IROM) + DT*QST_RT(K,J)*TB*TB
            QS2I_R1(K,J,IROM) = QS2I_R1(K,J,IROM) + DT*QST_RT(K,J)*T *TB 
C        
            QG1I_R2(K,J,IROM) = QG1I_R2(K,J,IROM) + DT*QGT_RT(K,J)*TB*T
            QG2I_R2(K,J,IROM) = QG2I_R2(K,J,IROM) + DT*QGT_RT(K,J)*T *T  
            QS1I_R2(K,J,IROM) = QS1I_R2(K,J,IROM) + DT*QST_RT(K,J)*TB*T
            QS2I_R2(K,J,IROM) = QS2I_R2(K,J,IROM) + DT*QST_RT(K,J)*T *T  
C        
            QG1I_RF(K,J,IROM) = QG1I_RF(K,J,IROM) + DT*QGT_RF(K,J)*TB
            QG2I_RF(K,J,IROM) = QG2I_RF(K,J,IROM) + DT*QGT_RF(K,J)*T     
            QS1I_RF(K,J,IROM) = QS1I_RF(K,J,IROM) + DT*QST_RF(K,J)*TB
            QS2I_RF(K,J,IROM) = QS2I_RF(K,J,IROM) + DT*QST_RF(K,J)*T     
           ENDDO
          ENDDO
 10     CONTINUE
C
C------ Romberg sequence using all previous grid results
        DO KROM = IROM, 2, -1
C-------- weight needed to cancel lowest-order error terms in KROM level
          W = 2.0 ** (2*(IROM-KROM+1))
          WM1 = W - 1.0
C
C-------- put Richardson extrapolation for KROM level into KROM-1 level
          DO K = 1, 2
            QG1I(K,KROM-1) = (W*QG1I(K,KROM) - QG1I(K,KROM-1)) / WM1
            QG2I(K,KROM-1) = (W*QG2I(K,KROM) - QG2I(K,KROM-1)) / WM1
            QS1I(K,KROM-1) = (W*QS1I(K,KROM) - QS1I(K,KROM-1)) / WM1
            QS2I(K,KROM-1) = (W*QS2I(K,KROM) - QS2I(K,KROM-1)) / WM1
            DO J = 1, 2
              QG1I_R1(K,J,KROM-1) = (W*QG1I_R1(K,J,KROM)
     &                               - QG1I_R1(K,J,KROM-1)) / WM1
              QG2I_R1(K,J,KROM-1) = (W*QG2I_R1(K,J,KROM)
     &                               - QG2I_R1(K,J,KROM-1)) / WM1
              QS1I_R1(K,J,KROM-1) = (W*QS1I_R1(K,J,KROM)
     &                               - QS1I_R1(K,J,KROM-1)) / WM1
              QS2I_R1(K,J,KROM-1) = (W*QS2I_R1(K,J,KROM)
     &                               - QS2I_R1(K,J,KROM-1)) / WM1
C
              QG1I_R2(K,J,KROM-1) = (W*QG1I_R2(K,J,KROM)
     &                               - QG1I_R2(K,J,KROM-1)) / WM1
              QG2I_R2(K,J,KROM-1) = (W*QG2I_R2(K,J,KROM)
     &                               - QG2I_R2(K,J,KROM-1)) / WM1
              QS1I_R2(K,J,KROM-1) = (W*QS1I_R2(K,J,KROM)
     &                               - QS1I_R2(K,J,KROM-1)) / WM1
              QS2I_R2(K,J,KROM-1) = (W*QS2I_R2(K,J,KROM)
     &                               - QS2I_R2(K,J,KROM-1)) / WM1
C
              QG1I_RF(K,J,KROM-1) = (W*QG1I_RF(K,J,KROM)
     &                               - QG1I_RF(K,J,KROM-1)) / WM1
              QG2I_RF(K,J,KROM-1) = (W*QG2I_RF(K,J,KROM)
     &                               - QG2I_RF(K,J,KROM-1)) / WM1
              QS1I_RF(K,J,KROM-1) = (W*QS1I_RF(K,J,KROM)
     &                               - QS1I_RF(K,J,KROM-1)) / WM1
              QS2I_RF(K,J,KROM-1) = (W*QS2I_RF(K,J,KROM)
     &                               - QS2I_RF(K,J,KROM-1)) / WM1
            ENDDO
          ENDDO
        ENDDO
C
        IF(IROM.GT.1) THEN
C------- compare the best-current and best-previous integrals
         ERRUG1 = QG1I(1,1) - QG1I(1,2)
         ERRVG1 = QG1I(2,1) - QG1I(2,2)
         ERRUG2 = QG2I(1,1) - QG2I(1,2)
         ERRVG2 = QG2I(2,1) - QG2I(2,2)
         ERRUS1 = QS1I(1,1) - QS1I(1,2)
         ERRVS1 = QS1I(2,1) - QS1I(2,2)
         ERRUS2 = QS2I(1,1) - QS2I(1,2)
         ERRVS2 = QS2I(2,1) - QS2I(2,2)
C
         ERR = MAX(
     &     ABS(ERRUG1),
     &     ABS(ERRVG1),
     &     ABS(ERRUG2),
     &     ABS(ERRVG2),
     &     ABS(ERRUS1),
     &     ABS(ERRVS1),
     &     ABS(ERRUS2),
     &     ABS(ERRVS2) )
         
c         write(13,1200) 
c     &     ABS(ERRUG1),
c     &     ABS(ERRVG1),
c     &     ABS(ERRUG2),
c     &     ABS(ERRVG2),
c     &     ABS(ERRUS1),
c     &     ABS(ERRVS1),
c     &     ABS(ERRUS2),
c     &     ABS(ERRVS2)
c 1200    format(8e16.9)

         IF(ERR*REFL .LT. ROMTOL) GO TO 101
        ENDIF
 100  CONTINUE
      WRITE(*,*) 'DLAMP: Romberg convergence failed.  Error =', ERR
C
 101  CONTINUE

ccc      write(*,*) IROM, ERR
C
C---- return best final results
      DELSQ = (R1(1)-R2(1))**2 + (R1(2)-R2(2))**2
      DELS  = SQRT(DELSQ)
      DELS_R1(1) =  (R1(1)-R2(1))/DELS
      DELS_R2(1) = -(R1(1)-R2(1))/DELS
      DELS_R1(2) =  (R1(2)-R2(2))/DELS
      DELS_R2(2) = -(R1(2)-R2(2))/DELS
C
      DO K = 1, 2
        QG1(K) = QG1I(K,1)*DELS
        QG2(K) = QG2I(K,1)*DELS
        QS1(K) = QS1I(K,1)*DELS
        QS2(K) = QS2I(K,1)*DELS
        DO J = 1, 2
          QG1_R1(K,J) = QG1I_R1(K,J,1)*DELS + QG1I(K,1)*DELS_R1(J)
          QG2_R1(K,J) = QG2I_R1(K,J,1)*DELS + QG2I(K,1)*DELS_R1(J)
          QS1_R1(K,J) = QS1I_R1(K,J,1)*DELS + QS1I(K,1)*DELS_R1(J)
          QS2_R1(K,J) = QS2I_R1(K,J,1)*DELS + QS2I(K,1)*DELS_R1(J)
C
          QG1_R2(K,J) = QG1I_R2(K,J,1)*DELS + QG1I(K,1)*DELS_R2(J)
          QG2_R2(K,J) = QG2I_R2(K,J,1)*DELS + QG2I(K,1)*DELS_R2(J)
          QS1_R2(K,J) = QS1I_R2(K,J,1)*DELS + QS1I(K,1)*DELS_R2(J)
          QS2_R2(K,J) = QS2I_R2(K,J,1)*DELS + QS2I(K,1)*DELS_R2(J)
C
          QG1_RF(K,J) = QG1I_RF(K,J,1)*DELS
          QG2_RF(K,J) = QG2I_RF(K,J,1)*DELS
          QS1_RF(K,J) = QS1I_RF(K,J,1)*DELS
          QS2_RF(K,J) = QS2I_RF(K,J,1)*DELS
        ENDDO
      ENDDO
C
      RETURN
      END ! DLAMP



      SUBROUTINE DLAMPC( R1,     R2,
     &          QG1, QG1_R1, QG1_R2,
     &          QG2, QG2_R1, QG2_R2,
     &          QS1, QS1_R1, QS1_R2,
     &          QS2, QS2_R1, QS2_R2 )
C-----------------------------------------------------------------------
C     Same as LAMPC, but also returns derivatives.
C-----------------------------------------------------------------------
      DIMENSION R1(2), R2(2)
      DIMENSION QG1(2), QG1_R1(2,2), QG1_R2(2,2),
     &          QG2(2), QG2_R1(2,2), QG2_R2(2,2),
     &          QS1(2), QS1_R1(2,2), QS1_R2(2,2),
     &          QS2(2), QS2_R1(2,2), QS2_R2(2,2)
C
      PARAMETER (NROMX=10)
C
      DIMENSION QG1I(2,NROMX),
     &          QG2I(2,NROMX),
     &          QS1I(2,NROMX),
     &          QS2I(2,NROMX)
      DIMENSION 
     &  QG1I_R1(2,2,NROMX),QG1I_R2(2,2,NROMX),QG1I_RF(2,2,NROMX),
     &  QG2I_R1(2,2,NROMX),QG2I_R2(2,2,NROMX),QG2I_RF(2,2,NROMX),
     &  QS1I_R1(2,2,NROMX),QS1I_R2(2,2,NROMX),QS1I_RF(2,2,NROMX),
     &  QS2I_R1(2,2,NROMX),QS2I_R2(2,2,NROMX),QS2I_RF(2,2,NROMX) 
C
      DIMENSION RF(2)
      DIMENSION RT(2)
      DIMENSION QGT(2), QGT_RT(2,2), QGT_RF(2,2),
     &          QST(2), QST_RT(2,2), QST_RF(2,2)
      DIMENSION QGA(2), QGA_RT(2,2), QGA_RF(2,2),
     &          QSA(2), QSA_RT(2,2), QSA_RF(2,2)
      DIMENSION DSQ_RT(2), DSQ_RF(2), DELS_R1(2), DELS_R2(2)
C
      DIMENSION QGAI(2), QGAI_R1(2,2), QGAI_R2(2,2), QGAI_RF(2,2),
     &          QSAI(2), QSAI_R1(2,2), QSAI_R2(2,2), QSAI_RF(2,2)
C
      PARAMETER(PI=3.14159265358979)
C
C---- Romberg convergence tolerance
      DATA ROMTOL / 1.0E-6 /
ccc   DATA ROMTOL / 1.0E-12 /
C
C---- reference length for convergence tolerance 
      REFL = 0.5*(R1(2) + R2(2))
C
      RF(1) = 0.5*(R1(1) + R2(1))
      RF(2) = 0.5*(R1(2) + R2(2))
C
C---- evaluate integrals over  0..t..1  on increasingly fine grids
      DO 100 IROM=1, NROMX
        NT = 2**IROM
C
        DO K = 1, 2
          QG1I(K,IROM) = 0.
          QG2I(K,IROM) = 0.
          QS1I(K,IROM) = 0.
          QS2I(K,IROM) = 0.
          DO J = 1, 2
            QG1I_R1(K,J,IROM) = 0.
            QG2I_R1(K,J,IROM) = 0.
            QS1I_R1(K,J,IROM) = 0.
            QS2I_R1(K,J,IROM) = 0.
            QG1I_R2(K,J,IROM) = 0.
            QG2I_R2(K,J,IROM) = 0.
            QS1I_R2(K,J,IROM) = 0.
            QS2I_R2(K,J,IROM) = 0.
            QG1I_RF(K,J,IROM) = 0.
            QG2I_RF(K,J,IROM) = 0.
            QS1I_RF(K,J,IROM) = 0.
            QS2I_RF(K,J,IROM) = 0.
          ENDDO
        ENDDO
C
C------ visit the midpoints of each of the NT intervals
        DO 10 IT=1, NT
          T = (FLOAT(IT) - 0.5) / FLOAT(NT)
          TB = 1.0 - T
C
          DT = 1.0/FLOAT(NT)
C
          RT(1) = R1(1)*TB + R2(1)*T
          RT(2) = R1(2)*TB + R2(2)*T
C
C-------- get induced velocities for vortex,source ring at XT,RT
          CALL DRING( RT(  1),    RT(  2),    RF(  1),    RF(  2),
     &     QGT(1),QGT_RT(1,1),QGT_RT(1,2),QGT_RF(1,1),QGT_RF(1,2), 
     &     QGT(2),QGT_RT(2,1),QGT_RT(2,2),QGT_RF(2,1),QGT_RF(2,2),
     &     QST(1),QST_RT(1,1),QST_RT(1,2),QST_RF(1,1),QST_RF(1,2),
     &     QST(2),QST_RT(2,1),QST_RT(2,2),QST_RF(2,1),QST_RF(2,2) )
C
C-------- singular parts of velocities in the limit  XT,RT -> XF,RF
          DSQ = (RT(1)-RF(1))**2 + (RT(2)-RF(2))**2
          DSQ_RT(1) =  2.0*(RT(1)-RF(1))
          DSQ_RT(2) =  2.0*(RT(2)-RF(2))
          DSQ_RF(1) = -2.0*(RT(1)-RF(1))
          DSQ_RF(2) = -2.0*(RT(2)-RF(2))
C
C
          PID4 = 4.0*PI*DSQ
          PIR8 = 8.0*PI*RF(2)
          DRSQ = DSQ/RF(2)**2
C
          DR1 = (RT(1)-RF(1))/PID4
          DR2 = (RT(2)-RF(2))/PID4
C
          QGA(1) =  DR2 - 0.5*LOG(DRSQ/64.0) / PIR8
          QGA(2) = -DR1
          QSA(1) = -DR1
          QSA(2) = -DR2 - 0.5*LOG(DRSQ     ) / PIR8
C
          QGA_RT(1,1) = (-DR2 - 0.5/PIR8)/DSQ * DSQ_RT(1)
          QGA_RT(1,2) = (-DR2 - 0.5/PIR8)/DSQ * DSQ_RT(2) + 1.0/PID4
          QGA_RF(1,1) = (-DR2 - 0.5/PIR8)/DSQ * DSQ_RF(1)
          QGA_RF(1,2) = (-DR2 - 0.5/PIR8)/DSQ * DSQ_RF(2) - 1.0/PID4
     &             + (1.0 + 0.5*LOG(DRSQ/64.0)) / (RF(2)*PIR8)
C
          QGA_RT(2,1) =   DR1            /DSQ * DSQ_RT(1) - 1.0/PID4
          QGA_RT(2,2) =   DR1            /DSQ * DSQ_RT(2)
          QGA_RF(2,1) =   DR1            /DSQ * DSQ_RF(1) + 1.0/PID4
          QGA_RF(2,2) =   DR1            /DSQ * DSQ_RF(2)
C
          QSA_RT(1,1) =   DR1            /DSQ * DSQ_RT(1) - 1.0/PID4
          QSA_RT(1,2) =   DR1            /DSQ * DSQ_RT(2)
          QSA_RF(1,1) =   DR1            /DSQ * DSQ_RF(1) + 1.0/PID4
          QSA_RF(1,2) =   DR1            /DSQ * DSQ_RF(2)
C
          QSA_RT(2,1) = ( DR2 - 0.5/PIR8)/DSQ * DSQ_RT(1)
          QSA_RT(2,2) = ( DR2 - 0.5/PIR8)/DSQ * DSQ_RT(2) - 1.0/PID4
          QSA_RF(2,1) = ( DR2 - 0.5/PIR8)/DSQ * DSQ_RF(1)
          QSA_RF(2,2) = ( DR2 - 0.5/PIR8)/DSQ * DSQ_RF(2) + 1.0/PID4
     &             + (1.0 + 0.5*LOG(DRSQ     )) / (RF(2)*PIR8)
C
C-------- accumulate integrals, with singular parts (at t=0.5) removed
          DO K = 1, 2
           QG1I(K,IROM) = QG1I(K,IROM) + DT*(QGT(K)*TB - QGA(K))
           QG2I(K,IROM) = QG2I(K,IROM) + DT*(QGT(K)*T  - QGA(K))
           QS1I(K,IROM) = QS1I(K,IROM) + DT*(QST(K)*TB - QSA(K))
           QS2I(K,IROM) = QS2I(K,IROM) + DT*(QST(K)*T  - QSA(K))
           DO J = 1, 2
            QG1I_R1(K,J,IROM) =
     &      QG1I_R1(K,J,IROM) + DT*(QGT_RT(K,J)*TB - QGA_RT(K,J))*TB
            QG2I_R1(K,J,IROM) =
     &      QG2I_R1(K,J,IROM) + DT*(QGT_RT(K,J)*T  - QGA_RT(K,J))*TB 
            QS1I_R1(K,J,IROM) =
     &      QS1I_R1(K,J,IROM) + DT*(QST_RT(K,J)*TB - QSA_RT(K,J))*TB
            QS2I_R1(K,J,IROM) =
     &      QS2I_R1(K,J,IROM) + DT*(QST_RT(K,J)*T  - QSA_RT(K,J))*TB 
C        
            QG1I_R2(K,J,IROM) =
     &      QG1I_R2(K,J,IROM) + DT*(QGT_RT(K,J)*TB - QGA_RT(K,J))*T
            QG2I_R2(K,J,IROM) =
     &      QG2I_R2(K,J,IROM) + DT*(QGT_RT(K,J)*T  - QGA_RT(K,J))*T
            QS1I_R2(K,J,IROM) =
     &      QS1I_R2(K,J,IROM) + DT*(QST_RT(K,J)*TB - QSA_RT(K,J))*T
            QS2I_R2(K,J,IROM) =
     &      QS2I_R2(K,J,IROM) + DT*(QST_RT(K,J)*T  - QSA_RT(K,J))*T
C        
            QG1I_RF(K,J,IROM) =
     &      QG1I_RF(K,J,IROM) + DT*(QGT_RF(K,J)*TB - QGA_RF(K,J))
            QG2I_RF(K,J,IROM) =
     &      QG2I_RF(K,J,IROM) + DT*(QGT_RF(K,J)*T  - QGA_RF(K,J))
            QS1I_RF(K,J,IROM) =
     &      QS1I_RF(K,J,IROM) + DT*(QST_RF(K,J)*TB - QSA_RF(K,J))
            QS2I_RF(K,J,IROM) =
     &      QS2I_RF(K,J,IROM) + DT*(QST_RF(K,J)*T  - QSA_RF(K,J))        
           ENDDO
          ENDDO
 10     CONTINUE
C
C------ Romberg sequence using all previous grid results
        DO KROM = IROM, 2, -1
C-------- weight needed to cancel lowest-order error terms in KROM level
          W = 2.0 ** (2*(IROM-KROM+1))
          WM1 = W - 1.0
C
C-------- put Richardson extrapolation for KROM level into KROM-1 level
          DO K = 1, 2
            QG1I(K,KROM-1) = (W*QG1I(K,KROM) - QG1I(K,KROM-1)) / WM1
            QG2I(K,KROM-1) = (W*QG2I(K,KROM) - QG2I(K,KROM-1)) / WM1
            QS1I(K,KROM-1) = (W*QS1I(K,KROM) - QS1I(K,KROM-1)) / WM1
            QS2I(K,KROM-1) = (W*QS2I(K,KROM) - QS2I(K,KROM-1)) / WM1
            DO J = 1, 2
              QG1I_R1(K,J,KROM-1) = (W*QG1I_R1(K,J,KROM)
     &                               - QG1I_R1(K,J,KROM-1)) / WM1
              QG2I_R1(K,J,KROM-1) = (W*QG2I_R1(K,J,KROM)
     &                               - QG2I_R1(K,J,KROM-1)) / WM1
              QS1I_R1(K,J,KROM-1) = (W*QS1I_R1(K,J,KROM)
     &                               - QS1I_R1(K,J,KROM-1)) / WM1
              QS2I_R1(K,J,KROM-1) = (W*QS2I_R1(K,J,KROM)
     &                               - QS2I_R1(K,J,KROM-1)) / WM1
C
              QG1I_R2(K,J,KROM-1) = (W*QG1I_R2(K,J,KROM)
     &                               - QG1I_R2(K,J,KROM-1)) / WM1
              QG2I_R2(K,J,KROM-1) = (W*QG2I_R2(K,J,KROM)
     &                               - QG2I_R2(K,J,KROM-1)) / WM1
              QS1I_R2(K,J,KROM-1) = (W*QS1I_R2(K,J,KROM)
     &                               - QS1I_R2(K,J,KROM-1)) / WM1
              QS2I_R2(K,J,KROM-1) = (W*QS2I_R2(K,J,KROM)
     &                               - QS2I_R2(K,J,KROM-1)) / WM1
C
              QG1I_RF(K,J,KROM-1) = (W*QG1I_RF(K,J,KROM)
     &                               - QG1I_RF(K,J,KROM-1)) / WM1
              QG2I_RF(K,J,KROM-1) = (W*QG2I_RF(K,J,KROM)
     &                               - QG2I_RF(K,J,KROM-1)) / WM1
              QS1I_RF(K,J,KROM-1) = (W*QS1I_RF(K,J,KROM)
     &                               - QS1I_RF(K,J,KROM-1)) / WM1
              QS2I_RF(K,J,KROM-1) = (W*QS2I_RF(K,J,KROM)
     &                               - QS2I_RF(K,J,KROM-1)) / WM1
            ENDDO
          ENDDO
        ENDDO
C
        IF(IROM.GT.1) THEN
C------- compare the best-current and best-previous integrals
         ERRUG1 = QG1I(1,1) - QG1I(1,2)
         ERRVG1 = QG1I(2,1) - QG1I(2,2)
         ERRUG2 = QG2I(1,1) - QG2I(1,2)
         ERRVG2 = QG2I(2,1) - QG2I(2,2)
         ERRUS1 = QS1I(1,1) - QS1I(1,2)
         ERRVS1 = QS1I(2,1) - QS1I(2,2)
         ERRUS2 = QS2I(1,1) - QS2I(1,2)
         ERRVS2 = QS2I(2,1) - QS2I(2,2)
C
         ERR = MAX(
     &     ABS(ERRUG1),
     &     ABS(ERRVG1),
     &     ABS(ERRUG2),
     &     ABS(ERRVG2),
     &     ABS(ERRUS1),
     &     ABS(ERRVS1),
     &     ABS(ERRUS2),
     &     ABS(ERRVS2) )
         
c         write(13,1200) 
c     &     ABS(ERRUG1),
c     &     ABS(ERRVG1),
c     &     ABS(ERRUG2),
c     &     ABS(ERRVG2),
c     &     ABS(ERRUS1),
c     &     ABS(ERRVS1),
c     &     ABS(ERRUS2),
c     &     ABS(ERRVS2)
c 1200    format(8e16.9)

         IF(ERR*REFL .LT. ROMTOL) GO TO 101
        ENDIF
 100  CONTINUE
      WRITE(*,*) 'DLAMPC: Romberg convergence failed.  Error =', ERR
C
 101  CONTINUE

ccc      write(*,*) IROM, ERR
C
C---- return best final results
      DELSQ = (R1(1)-R2(1))**2 + (R1(2)-R2(2))**2
      DELS  = SQRT(DELSQ)
      DELS_R1(1) =  (R1(1)-R2(1))/DELS
      DELS_R2(1) = -(R1(1)-R2(1))/DELS
      DELS_R1(2) =  (R1(2)-R2(2))/DELS
      DELS_R2(2) = -(R1(2)-R2(2))/DELS
C
C---- analytically-integrated singular parts which were removed
      QGAI(1) = (1.0 + LOG(16.0*RF(2)/DELS)) / (8.0*PI*RF(2))
      QGAI(2) = 0.
      QSAI(1) = 0.
      QSAI(2) = (1.0 + LOG( 2.0*RF(2)/DELS)) / (8.0*PI*RF(2))
      DO J = 1, 2
        QGAI_R1(1,J) = (-DELS_R1(J)/DELS) / (8.0*PI*RF(2))
        QGAI_R1(2,J) = 0.
        QSAI_R1(1,J) = 0.
        QSAI_R1(2,J) = (-DELS_R1(J)/DELS) / (8.0*PI*RF(2))
C
        QGAI_R2(1,J) = (-DELS_R2(J)/DELS) / (8.0*PI*RF(2))
        QGAI_R2(2,J) = 0.
        QSAI_R2(1,J) = 0.
        QSAI_R2(2,J) = (-DELS_R2(J)/DELS) / (8.0*PI*RF(2))
      ENDDO
C
      QGAI_RF(1,1) = 0.
      QGAI_RF(2,1) = 0.
      QSAI_RF(1,1) = 0.
      QSAI_RF(2,1) = 0.
C
      QGAI_RF(1,2) = (1.0/RF(2)) / (8.0*PI*RF(2)) - QGAI(1)/RF(2)
      QGAI_RF(2,2) = 0.
      QSAI_RF(1,2) = 0.
      QSAI_RF(2,2) = (1.0/RF(2)) / (8.0*PI*RF(2)) - QSAI(2)/RF(2)
C
C
      DO K = 1, 2
        QG1(K) = (QG1I(K,1) + QGAI(K))*DELS
        QG2(K) = (QG2I(K,1) + QGAI(K))*DELS
        QS1(K) = (QS1I(K,1) + QSAI(K))*DELS
        QS2(K) = (QS2I(K,1) + QSAI(K))*DELS
        DO J = 1, 2
          QG1_R1(K,J) = (QG1I_R1(K,J,1) + QGAI_R1(K,J))*DELS
     &                + (QG1I   (K,  1) + QGAI   (K  ))*DELS_R1(J)
     &            + 0.5*(QG1I_RF(K,J,1) + QGAI_RF(K,J))*DELS
          QG2_R1(K,J) = (QG2I_R1(K,J,1) + QGAI_R1(K,J))*DELS
     &                + (QG2I   (K,  1) + QGAI   (K  ))*DELS_R1(J)
     &            + 0.5*(QG2I_RF(K,J,1) + QGAI_RF(K,J))*DELS
          QS1_R1(K,J) = (QS1I_R1(K,J,1) + QSAI_R1(K,J))*DELS
     &                + (QS1I   (K,  1) + QSAI   (K  ))*DELS_R1(J)
     &            + 0.5*(QS1I_RF(K,J,1) + QSAI_RF(K,J))*DELS
          QS2_R1(K,J) = (QS2I_R1(K,J,1) + QSAI_R1(K,J))*DELS
     &                + (QS2I   (K,  1) + QSAI   (K  ))*DELS_R1(J)
     &            + 0.5*(QS2I_RF(K,J,1) + QSAI_RF(K,J))*DELS
C
          QG1_R2(K,J) = (QG1I_R2(K,J,1) + QGAI_R2(K,J))*DELS
     &                + (QG1I   (K,  1) + QGAI   (K  ))*DELS_R2(J)
     &            + 0.5*(QG1I_RF(K,J,1) + QGAI_RF(K,J))*DELS
          QG2_R2(K,J) = (QG2I_R2(K,J,1) + QGAI_R2(K,J))*DELS
     &                + (QG2I   (K,  1) + QGAI   (K  ))*DELS_R2(J)
     &            + 0.5*(QG2I_RF(K,J,1) + QGAI_RF(K,J))*DELS
          QS1_R2(K,J) = (QS1I_R2(K,J,1) + QSAI_R2(K,J))*DELS
     &                + (QS1I   (K,  1) + QSAI   (K  ))*DELS_R2(J)
     &            + 0.5*(QS1I_RF(K,J,1) + QSAI_RF(K,J))*DELS
          QS2_R2(K,J) = (QS2I_R2(K,J,1) + QSAI_R2(K,J))*DELS
     &                + (QS2I   (K,  1) + QSAI   (K  ))*DELS_R2(J)
     &            + 0.5*(QS2I_RF(K,J,1) + QSAI_RF(K,J))*DELS
        ENDDO
      ENDDO
C
      RETURN
      END ! DLAMPC




      SUBROUTINE RING( XV, RV, XF, RF, UX, UR, SX, SR )
C-----------------------------------------------------------------------
C     Computes the velocities induced by a vortex/source ring 
C     located at XV with radius RV with unit circulation 
C     and unit source/length density.
C
C     Adapted from routines provided by J. Kerwin.
C-----------------------------------------------------------------------
C  Input:
C     XV,RV  x,r of the ring
C     XF,RF  x,r of the field point
C
C  Output:
C     UX,UR  velocity at XF,RF for unit circulation (+ counterclockwise)
C     SX,SR  velocity at XF,RF for unit source/perimeter
C-----------------------------------------------------------------------
      PARAMETER(PI=3.14159265358979)
C
      IF(RV.LE.0.0) THEN
C------ zero-radius ring
        UX = 0.
        UR = 0.
        SX = 0.
        SR = 0.
        RETURN
      ENDIF
C
C---- this fails if R=1 and X=0  (on the ring itself)
      R = RF/RV
      X = (XV-XF)/RV
C
      IF(R.EQ.1.0 .AND. X.EQ.0.0) THEN
       UX = 0.
       UR = 0.
       SX = 0.
       SR = 0.
       RETURN
      ENDIF
C
      IF (RF .EQ. 0.0) THEN
C----- Control Point on the axis
       UX = 1.0 / SQRT(1.0+X**2)**3 / (2.0*RV)
       UR = 0.0
       SX =  -X / SQRT(1.0+X**2)**3 / (2.0*RV)
       SR = 0.0
C
      ELSE
C----- Control Point not on X-axis
       XRP = X**2 + (1.0+R)**2
       XRM = X**2 + (1.0-R)**2
C
       SRP = SQRT(XRP)
C
       AK = XRM/XRP
       CALL ELLEK(AK,ELE,ELK)
C
       F = 2.0/XRM
C
       UX = ( 1.0/ SRP   )*(ELK - ELE*(1.0 + F*(R-1.0)))/(2.0*PI*RV)
       UR = (   X/(SRP*R))*(ELK - ELE*(1.0 + F* R     ))/(2.0*PI*RV)
C
       SX =  (  X/ SRP   )*(    - ELE*       F         )/(2.0*PI*RV)
       SR =  (1.0/(SRP*R))*(ELK - ELE*(1.0 + F*(R-R*R)))/(2.0*PI*RV)
C
ccC----- streamfunction due to vortex
cc       PSI = ((1.0 - 2.0*XRP*R)*ELK - ELE)*RV / (2.0*PI*SQRT(XRP))
      ENDIF
C
      RETURN
      END




      SUBROUTINE DRING(        XV,    RV,    XF,    RF, 
     &                  UX, UX_XV, UX_RV, UX_XF, UX_RF,
     &                  UR, UR_XV, UR_RV, UR_XF, UR_RF,
     &                  SX, SX_XV, SX_RV, SX_XF, SX_RF,
     &                  SR, SR_XV, SR_RV, SR_XF, SR_RF )
C-----------------------------------------------------------------------
C     Same as RING, but also returns AIC derivatives w.r.t. geometry
C-----------------------------------------------------------------------
      PARAMETER(PI=3.14159265358979)
C
      UX = 0.
      UR = 0.
      SX = 0.
      SR = 0.
C
      UX_X = 0.
      UR_X = 0.
      SX_X = 0.
      SR_X = 0.
C
      UX_R = 0.
      UR_R = 0.
      SX_R = 0.
      SR_R = 0.
C
      RVI = 0.
C
      IF(RV.LE.0.0) THEN
C----- zero-radius ring
       GO TO 90
      ENDIF
C
      RVI = 1.0/RV
C
C---- this fails if R=1 and X=0  (on the ring itself)
      R =     RF *RVI
      X = (XV-XF)*RVI
C
      IF(R.EQ.1.0 .AND. X.EQ.0.0) THEN
       GO TO 90
      ENDIF
C
      IF (RF .EQ. 0.0) THEN
C----- Control Point on the axis
       UX = 1.0 / SQRT(1.0+X**2)**3 * PI  ! / (2.0*PI*RV)
       SX =  -X / SQRT(1.0+X**2)**3 * PI  ! / (2.0*PI*RV)
C
       UX_X = -3.0*X*UX/(1.0+X**2)
       SX_X = -3.0*X*UX/(1.0+X**2) - UX
C
      ELSE
C----- Control Point not on X-axis
       XRP = X**2 + (1.0+R)**2
       XRM = X**2 + (1.0-R)**2
C
       XRP_X =  2.0*X
       XRM_X =  2.0*X
       XRP_R =  2.0*(1.0+R)
       XRM_R = -2.0*(1.0-R)
C
       SRP = SQRT(XRP)
       SRP_X = (0.5/SRP)*XRP_X
       SRP_R = (0.5/SRP)*XRP_R
C
       AK   =  XRM/XRP
       AK_X = (XRM_X - AK*XRP_X)/XRP
       AK_R = (XRM_R - AK*XRP_R)/XRP
C
       CALL DELLEK(AK, ELE, DELE,
     &                 ELK, DELK )
       ELE_X = DELE*AK_X
       ELE_R = DELE*AK_R
       ELK_X = DELK*AK_X
       ELK_R = DELK*AK_R
C
       F = 2.0/XRM
       F_X = (-F/XRM)*XRM_X
       F_R = (-F/XRM)*XRM_R
C
       UX   = ( 1.0/ SRP   )*(ELK   - ELE  *(1.0 + F  *(R-1.0)))
       UX_X = ( 1.0/ SRP   )*(ELK_X - ELE_X*(1.0 + F  *(R-1.0))
     &                              - ELE  *(      F_X*(R-1.0)))
     &      -   (UX/ SRP)*SRP_X
       UX_R = ( 1.0/ SRP   )*(ELK_R - ELE_R*(1.0 + F  *(R-1.0))
     &                              - ELE  *(      F_R*(R-1.0))
     &                              - ELE  *       F           )
     &      -   (UX/SRP)*SRP_R
C
       UR   = (   X/(SRP*R))*(ELK   - ELE  *(1.0 + F  * R     ))
       UR_X = (   X/(SRP*R))*(ELK_X - ELE_X*(1.0 + F  * R     )
     &                              - ELE  *(      F_X* R     ))
     &      + ( 1.0/(SRP*R))*(ELK   - ELE  *(1.0 + F  * R     ))
     &      -   (UR/SRP)*SRP_X
       UR_R = (   X/(SRP*R))*(ELK_R - ELE_R*(1.0 + F  * R     )
     &                              - ELE  *(      F_R* R     
     &                                           + F          ))
     &      -   (UR/SRP)*SRP_R
     &      -    UR/R
C
       SX   = (   X/ SRP   )*(      - ELE  *       F           )
       SX_X = (   X/ SRP   )*(      - ELE_X*       F
     &                              - ELE  *       F_X         )
     &      + ( 1.0/ SRP   )*(      - ELE  *       F           )
     &      -   (SX/SRP)*SRP_X
       SX_R = (   X/ SRP   )*(      - ELE_R*       F 
     &                              - ELE  *       F_R         )
     &      -   (SX/SRP)*SRP_R
C
       SR   = ( 1.0/(SRP*R))*(ELK   - ELE  *(1.0 + F  *(R-R*R)))
       SR_X = ( 1.0/(SRP*R))*(ELK_X - ELE_X*(1.0 + F  *(R-R*R))
     &                              - ELE  *(      F_X*(R-R*R)))
     &      -   (SR/ SRP)*SRP_X
       SR_R = ( 1.0/(SRP*R))*(ELK_R - ELE_R*(1.0 + F  *(R-R*R))
     &                              - ELE  *(      F_R*(R-R*R)
     &                                           + F  *(1.0-2.0*R)))
     &      -   (SR/SRP)*SRP_R
     &      -    SR/R
C
      ENDIF
C
 90   CONTINUE
cc    X = (XV-XF)/RV
      X_XV =  RVI
      X_RV = -X*RVI
      X_XF = -RVI
C
cc    R = RF/RV
      R_RV = -R*RVI
      R_RF =  RVI
C
      HRIP = RVI/(2.0*PI)
C
      UX    = HRIP*UX
      UR    = HRIP*UR
      SX    = HRIP*SX
      SR    = HRIP*SR
C
      UX_XV = HRIP* UX_X*X_XV
      UR_XV = HRIP* UR_X*X_XV
      SX_XV = HRIP* SX_X*X_XV
      SR_XV = HRIP* SR_X*X_XV
      UX_RV = HRIP*(UX_X*X_RV + UX_R*R_RV) - UX*RVI
      UR_RV = HRIP*(UR_X*X_RV + UR_R*R_RV) - UR*RVI
      SX_RV = HRIP*(SX_X*X_RV + SX_R*R_RV) - SX*RVI
      SR_RV = HRIP*(SR_X*X_RV + SR_R*R_RV) - SR*RVI
      UX_XF = HRIP* UX_X*X_XF
      UR_XF = HRIP* UR_X*X_XF
      SX_XF = HRIP* SX_X*X_XF
      SR_XF = HRIP* SR_X*X_XF
      UX_RF = HRIP*             UX_R*R_RF
      UR_RF = HRIP*             UR_R*R_RF
      SX_RF = HRIP*             SX_R*R_RF
      SR_RF = HRIP*             SR_R*R_RF
C
      RETURN
      END




      SUBROUTINE ELLEK(AK,ELE,ELK)
C-----------------------------------------------------------------------
C              ELLIPTIC FUNCTIONS ROUTINE
C
C     Adapted from routines provided by J. Kerwin.
C-----------------------------------------------------------------------
C   Input
C     AK     elliptic-integral argument
C
C   Output
C     ELE    complete elliptic integral of the second kind
C     ELK    complete elliptic integral of the first  kind
C_______________________________________________________________________
C
      ALK = -LOG(AK)
C
      ELE = 1.00000000000
     &    +(0.44325141463
     &    +(0.06260601220
     &    +(0.04757383546
     &    + 0.01736506451*AK)*AK)*AK)*AK
     &  +( (0.24998368310
     &    +(0.09200180037
     &    +(0.04069697526
     &    + 0.00526449639*AK)*AK)*AK)*AK )*ALK
C          
      ELK = 1.38629436112
     &    +(0.09666344259
     &    +(0.03590092383
     &    +(0.03742563713
     &    + 0.01451196212*AK)*AK)*AK)*AK
     &  +(  0.50000000000
     &    +(0.12498593597
     &    +(0.06880248576
     &    +(0.03328355346
     &    + 0.00441787012*AK)*AK)*AK)*AK )*ALK
C
      RETURN
      END      



      SUBROUTINE DELLEK(AK, ELE, ELE_AK,
     &                      ELK, ELK_AK )
C-----------------------------------------------------------------------
C              ELLIPTIC FUNCTIONS + DERIVATIVE ROUTINE
C
C     Adapted from routines provided by J. Kerwin.
C-----------------------------------------------------------------------
C   Input
C     AK     elliptic-integral argument
C
C   Output
C     ELE     complete elliptic integral of the second kind
C     ELK     complete elliptic integral of the first  kind
C     ELE_AK  d(ELE)/d(AK)
C     ELK_AK  d(ELK)/d(AK)
C_______________________________________________________________________
C
      ALK = -LOG(AK)
C
      ELE = 1.00000000000
     &    +(0.44325141463
     &    +(0.06260601220
     &    +(0.04757383546
     &    + 0.01736506451*AK)*AK)*AK)*AK
     &  +( (0.24998368310
     &    +(0.09200180037
     &    +(0.04069697526
     &    + 0.00526449639*AK)*AK)*AK)*AK )*ALK
C          
      ELK = 1.38629436112
     &    +(0.09666344259
     &    +(0.03590092383
     &    +(0.03742563713
     &    + 0.01451196212*AK)*AK)*AK)*AK
     &  +(  0.50000000000
     &    +(0.12498593597
     &    +(0.06880248576
     &    +(0.03328355346
     &    + 0.00441787012*AK)*AK)*AK)*AK )*ALK
C
      ELE_AK =
     &      0.19326773153
     &    +(0.03321022403
     &    +(0.10202453112
     &    + 0.06419576165*AK)*AK)*AK
     &  +(  0.24998368310
     &    +(0.18400360074
     &    +(0.12209092578
     &    + 0.02105798556*AK)*AK)*AK )*ALK
C           
      ELK_AK =
     &     -0.028322493380
     &    +(0.002999361900
     &    +(0.078993357930
     &    + 0.053629978360*AK)*AK)*AK
     &    - 0.50000000000/AK
     &    +(0.12498593597
     &    +(0.13760497152
     &    +(0.09985066038
     &    + 0.01767148048*AK)*AK)*AK )*ALK
C
      RETURN
      END      







