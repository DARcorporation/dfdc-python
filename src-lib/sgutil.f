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

      SUBROUTINE SGCURV(LQUERY,LCPLOT,NB,XB,XPB,YB,YPB,SB,
     &                  NGX,NG,SG,
     &                  NPTS,CVEX,SMOF,FSL,FSR,
     &                  NREF,SREF1,SREF2,CRRAT)
      IMPLICIT REAL (A-H,M,O-Z)
      LOGICAL LQUERY, LCPLOT
      DIMENSION XB(NB),XPB(NB), 
     &          YB(NB),YPB(NB),SB(NB)
      DIMENSION SG(NGX)
      DIMENSION SREF1(NREF),SREF2(NREF),CRRAT(NREF)
C---------------------------------------------------------------------
C     Sets spacing array SG based on curvature and a few parameters
C
C  Input:  LQUERY   if T, requests user input for parameters before paneling
C          LCPLOT   if T, spacing parameters are plotted
C          NB       number of coordinates
C          XB(i)    x coordinates
C          XPB(i)   dXB/dSB spline derivative
C          YB(i)    y coordinates
C          YPB(i)   dYB/dSB spline derivative
C          SB(i)    spline parameter (typically arc length)
C          NGX      max number of points which can be generated
C
C  Input, also Output if interactively modified:
C          NG       number of points actually generated
C          SG(k)    fractional spacing array  0.0 .. 1.0
C          NPTS     number of points to be generated
C          CVEX     curvature-exponent parameter
C          SMOF     smoothing-length factor
C          FSL      fractional spacing at left  endpoint
C          FSR      fractional spacing at right endpoint
C          NREF     number of local-refinement regions
C          SREF1(r) first point of refinement region
C          SREF2(r) last  point of refinement region
C          CRRAT(r) fractional spacing over refinement region
C
C
C---------------------------------------------------------------------
      PARAMETER (IX=1601)
      DIMENSION S(IX), CV(IX), CVS(IX), CVINT(IX), SNEW(IX)
      DIMENSION AA(IX), BB(IX), CC(IX)
      DIMENSION CV0(IX)
C
      PARAMETER (NREFX=20)
      DIMENSION CVR(NREFX)
C
      CHARACTER*4 COMAND
      CHARACTER*128 COMARG
      CHARACTER*1 ANS
C
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
C
      CHARACTER*4 DUMMY
C
C---- max number of point-setting passes
      DATA NPASS / 1 /
C
      IF(NB  .GT.IX   ) STOP 'SGCURV: Local array overflow on IX'
      IF(NREF.GT.NREFX) STOP 'SGCURV: Local array overflow on NREFX'
C
C
      SBTOT = SB(NB) - SB(1)
      IF(SBTOT.LE.0.0) THEN
       WRITE(*,*) '? SGCURV: Zero element arc length'
       RETURN
      ENDIF
C
C---- if no interaction or no spacing is present, go compute it first
      IF(.NOT.LQUERY .OR. NG.EQ.0) GO TO 100
C
C=====================================================================
C==== begin user-interaction loop
C
 10   CONTINUE
C
C---- plot current spacing distribution
      CALL PLNODE(NB,XB,XPB,YB,YPB,SB, SG,NG,
     &               XSIZ,YSIZ,SREF1,SREF2,NREF)
C
      WRITE(*,2010) NPTS, FSL,FSR, CVEX, SMOF
      DO IREF = 1, NREF
        IF(CRRAT(IREF) .NE. 0.0  .AND.  SREF1(IREF) .LT. SREF2(IREF))
     &      WRITE(*,2015) IREF ,SREF1(IREF),SREF2(IREF),CRRAT(IREF)
      ENDDO
C      
 2010 FORMAT(
     &    '  __________________________________________________'
     &  //'   Z oom'
     &   /'   U nzoom'
     &  //'   N i   number of panel nodes         ', I5
     &   /'   D rr  left,right end-spacing ratios ', 2F7.3
     &   /'   C r   curvature scaling exponent    ', F10.4
     &   /'   S r   curvature smoothing factor    ', F10.4
     &   /'   # rrr local refinement parameters...')
 2015 FORMAT(
     &    '  ',
     &    I1,'                                   ', 2F8.4, F8.3)
C
 15   CALL ASKC('  Enter  key, parameters  (<Return> if spacing OK)^',
     &           COMAND,COMARG)
C
      ANS = COMAND(1:1)
C
C---------------------------------------------
C---- first test for simple commands
      IF (ANS.EQ.' ') THEN
        CALL CLRZOOM
        RETURN
      ELSEIF(INDEX('Zz',ANS).NE.0) THEN
        CALL USETZOOM(.TRUE.,.TRUE.)
        CALL REPLOT(IDEV)
        GO TO 10
      ELSEIF(INDEX('Uu',ANS).NE.0) THEN
        CALL CLRZOOM
        CALL REPLOT(IDEV)
        GO TO 10
      ELSEIF(INDEX('Pp',ANS).NE.0) THEN
        LCPLOT = .NOT.LCPLOT
        GO TO 10
      ENDIF
C
C---- extract numerical arguments, if any
      DO I=1, 20
        IINPUT(I) = 0
        RINPUT(I) = 0.0
      ENDDO
      NINPUT = 20
      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
      NINPUT = 20
      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C
      IF (INDEX('Nn',ANS).NE.0) GO TO 30
      IF (INDEX('Dd',ANS).NE.0) GO TO 40
      IF (INDEX('Cc',ANS).NE.0) GO TO 50
      IF (INDEX('Ss',ANS).NE.0) GO TO 60
      IF (INDEX('123456789',ANS).NE.0) GO TO 90
C
C---- no valid command... ask again
      GO TO 10
C
C---------------------------------------------
 30   CONTINUE
      IF(NINPUT.LT.1) THEN
 31    IINPUT(1) = NPTS
       WRITE(*,2030) IINPUT(1)
 2030  FORMAT(' Enter number of points:',I5)
       CALL READI(1,IINPUT,ERROR)
       IF(ERROR) GO TO 31
      ENDIF
C
      NPTS = IINPUT(1)
      IF(NPTS.LT.2) THEN
       WRITE(*,*) 'Must have at least 2 points'
       NINPUT = 0
       GO TO 30
      ENDIF
C
      IF(NPTS.GT.NGX) THEN
       NPTS = NGX
       WRITE(*,*) 'Number of points limited to array limit of', NGX
      ENDIF
      GO TO 100
C
C---------------------------------------------
 40   CONTINUE
      IF(NINPUT.LT.2) THEN
 41    RINPUT(1) = FSL
       RINPUT(2) = FSR
       WRITE(*,2040) RINPUT(1), RINPUT(2)
 2040  FORMAT(' Enter new  dsL/dsAvg, dsR/dsAvg:',2F9.4)
       CALL READR(2,RINPUT,ERROR)
       IF(ERROR) GO TO 41
      ENDIF
C
      FSL = RINPUT(1)
      FSR = RINPUT(2)
C
      IF(FSL.LE.0.0 .OR. FSR.LE.0.0) THEN
       NINPUT = 0
       GO TO 40
      ENDIF
      GO TO 100
C
C---------------------------------------------
 50   CONTINUE
      IF(NINPUT.LT.1) THEN
 51    RINPUT(1) = CVEX
       WRITE(*,2050) RINPUT(1)
 2050  FORMAT(' Enter new curvature exponent (0 to ~5):',F9.4)
       CALL READR(1,RINPUT,ERROR)
       IF(ERROR) GO TO 51
      ENDIF
C
      CVEX = RINPUT(1)
      IF(CVEX .LT. 0.0) THEN
       NINPUT = 0
       GO TO 50
      ENDIF
C
      CVEX = MAX(CVEX , 0.0001)
      GO TO 100
C
C---------------------------------------------
 60   CONTINUE
      IF(NINPUT.LT.1) THEN
 61    RINPUT(1) = SMOF
       WRITE(*,2060) RINPUT(1)
 2060  FORMAT(' Enter new curvature smoothing factor (0 to ~5):',F9.4)
       CALL READR(1,RINPUT,ERROR)
       IF(ERROR) GO TO 61
      ENDIF
C
      SMOF = RINPUT(1)
      IF(SMOF .LT. 0.0) THEN
       NINPUT = 0
       GO TO 60
      ENDIF
C
      GO TO 100
C
C---------------------------------------------
 90   CONTINUE
      IREF = IINPUT(1)
      IF(IREF.GT.NREF) THEN
       WRITE(*,*) 'Number of local station refinements limited to',NREF
       GO TO 10
      ENDIF
C
      IF(NINPUT.LT.3) THEN
 91    RINPUT(2) = SREF1(IREF)
       RINPUT(3) = SREF2(IREF)
       RINPUT(4) = CRRAT(IREF)
       WRITE(*,2090) IREF, RINPUT(2), RINPUT(3), RINPUT(4)
 2090  FORMAT(' Enter new station-', I1,
     &  '  s1/smax, s2/smax, ds/dsAvg:', 2F8.4,F8.3)
       CALL READR(3,RINPUT(2),ERROR)
       IF(ERROR) GO TO 91
      ENDIF
C
      SREF1(IREF) = ABS(RINPUT(2))
      SREF2(IREF) = ABS(RINPUT(3))
      CRRAT(IREF) = ABS(RINPUT(4))
      IF(SREF1(IREF).GT.SREF2(IREF)) THEN
       SREF1(IREF) = RINPUT(3)
       SREF2(IREF) = RINPUT(2)
      ENDIF
      GO TO 100
C
C================================================================
C==== Set surface panel node distribution based on curvature 
 100  CONTINUE
C
      NG = NPTS
C
      IF(NG.LT.2) THEN
       WRITE(*,*) '? SGCURV: Must specify at least two panel nodes'
       RETURN
      ENDIF
C
      KMULT = MIN( 8 , (IX-1)/NB )
C
C---- set up working coordinate array
      I = 1
      IB = 1
      S(I) = SB(IB)
      DO IB = 1, NB-1
        DSB = SB(IB+1) - SB(IB)
C
        IF(DSB.GT.0.0) THEN
         KK = KMULT
        ELSE
         KK = 1
        ENDIF
C
        DO K = 1, KK
          FRAC = FLOAT(K)/FLOAT(KK)
          I = I+1
          IF(I.GT.IX) STOP 'SGCURV:  Array overflow on IX'
          S(I) = SB(IB+1)*FRAC + SB(IB)*(1.0-FRAC)
        ENDDO
      ENDDO
C
C---- number of points in working array
      N = I
C
      STOT = S(N) - S(1)
      DO 200 IPASS = 1, NPASS
C
      DSAVG = STOT / FLOAT(NG-1)
C
      DSLE = FSL * DSAVG
      DSTE = FSR * DSAVG
C
C---- set up curvature array (nondimensionalized with perimeter)
      DO I = 1, N
        CV(I) = CURV(S(I),XB,XPB,YB,YPB,SB,NB) * STOT
        CV(I) = ABS(CV(I)) + 1.0
      ENDDO
C
C---- reset curvature at corner nodes from adjacent nodes
      DO I = 2, N-2
        IF(S(I) .EQ. S(I+1)) THEN
         CV(I) = 0.5*(CV(I-1) + CV(I+2))
         CV(I+1) = CV(I)
        ENDIF
      ENDDO
C
C---- raise curvature to power sqrt(CVEX), limiting to prevent overflows
      DO I = 1, N
        GCV = MIN( 0.5*CVEX*LOG(CV(I)) , 12.0 )
        CV(I) = EXP(GCV)
      ENDDO
C
C---- set max and average curvature 
      CVMAX = 0.01
      CVAVG = 0.
      DO I = 1, N-1
        CVASQ = 0.5*(CV(I)**2 + CV(I+1)**2)
        CVMAX = MAX( CVMAX , SQRT(CVASQ) )
        CVAVG = CVAVG + CVASQ*(S(I+1)-S(I))
      ENDDO
      CVAVG = SQRT(CVAVG/(S(N)-S(1)))
C
      CVAVG = MAX( CVAVG , 1.0 )
      CVMAX = MAX( CVMAX , 1.0 )
C
C---- set artificial curvature at ends to get approximate ds/c there
      CVLE = CVAVG * DSAVG/DSLE
      CVTE = CVAVG * DSAVG/DSTE
      CV(1) = CVLE
      CV(N) = CVTE
C
      DO I=1, N
        CV0(I) = CV(I)
      ENDDO
C
C---- set curvature smoothing length
      CURVL = 0.008*STOT
C
C---- set up implicit system for smoothed curvature array CV
      CURVK = CURVL**2 * SMOF
      AA(1) = 1.0
      CC(1) = 0.
      DO I=2, N-1
        DSM = S(I) - S(I-1)
        DSP = S(I+1) - S(I)
        DSA = 0.5*(DSM+DSP)
        IF(DSM.EQ.0.0 .OR. DSP.EQ.0.0) THEN
         BB(I) = 0.0
         AA(I) = 1.0
         CC(I) = 0.0
        ELSE
         BB(I) = - CURVK/(DSM*DSA)
         AA(I) =   CURVK/(DSM*DSA) + CURVK/(DSP*DSA)  +  1.0
         CC(I) =                   - CURVK/(DSP*DSA)
        ENDIF
      ENDDO
      AA(N) = 1.0
      BB(N) = 0.
C
C---- set artificial curvature at the bunching points
      DO IREF = 1, NREF
        IF(CRRAT(IREF) .GT. 0.0) THEN
         CVR(IREF) = CVAVG / MAX( CRRAT(IREF) , 0.0001 )
        ELSE
         CVR(IREF) = 0.0
        ENDIF
C
        DO I = 2, N-1
          IF(CRRAT(IREF).GT.0.0) THEN
C--------- set modified curvature if point is in refinement area
           SGI = (S(I)-S(1))/STOT
           IF(SGI.GT.SREF1(IREF) .AND. SGI.LT.SREF2(IREF)) THEN
cc            BB(I) = 0.
cc            AA(I) = 1.0
cc            CC(I) = 0.
            CV(I) = CVR(IREF)
           ENDIF
          ENDIF
        ENDDO
      ENDDO
C
C---- calculate smoothed curvature array
      CALL TRISOL(AA,BB,CC,CV,N)
C
C---- spline curvature array
      CALL SEGSPL(CV,CVS,S,N)
C
C---- Integrate exponentiated curvature with arc length
      CVINT(1) = 0.
      DO I = 2, N
        SM = S(I-1)
        SO = S(I)
        CVM = SEVAL(SM,CV,CVS,S,N)
        CVO = SEVAL(SO,CV,CVS,S,N)
        CVM = MAX( CVM , 1.0 )
        CVO = MAX( CVO , 1.0 )
C
        CVINT(I) = CVINT(I-1) + 0.5*(CVO+CVM)*(SO-SM)
      END DO
C
      DO I = 1, N
        CVINT(I) = CVINT(I) + S(I)
      END DO
C
C---- plot spacing parameter
      IF(LCPLOT) THEN
       CALL PLTCRV(S,CV0,CV,CVINT,N,CVAVG)
       CALL ASKS('Hit <cr>^',DUMMY)
      ENDIF
C
C---- Calculate normalized surface spacing distribution arrays
      CALL SETCRV(SNEW,NG,S,CVINT,N)
C
      DO IG = 1, NG
        SG(IG) = (SNEW(IG)-SNEW(1))/(SNEW(NG)-SNEW(1))
      END DO
C
C---- calculate actual obtained grid spacings
      DLE = ABS( SNEW(2)  - SNEW(1)    )
      DTE = ABS( SNEW(NG) - SNEW(NG-1) )
      DSMIN = DLE
      DSMAX = 0.0
      DO IG = 2, NG
        DSMIN = MIN( DSMIN , ABS(SNEW(IG)-SNEW(IG-1)) )
        DSMAX = MAX( DSMAX , ABS(SNEW(IG)-SNEW(IG-1)) )
      END DO
C
 200  CONTINUE
C
C-----------------------------------------------------------------
C---- display what parameters were used, and what we got
C
c      WRITE(*,3020) NG, FSL, FSR, CVEX, SMOF
c      DO IREF = 1, NREF
c        IF(CRRAT(IREF) .NE. 0.0  .AND.  SREF1(IREF) .LT. SREF2(IREF))
c     &      WRITE(*,3040) IREF ,SREF1(IREF),SREF2(IREF),CRRAT(IREF)
c      ENDDO
C
c 3020 FORMAT(
c     &  /'  ........................................................'
c     & //'  Number of panel nodes =', I5
c     &  /'  Prescribed panel node spacing parameters ...'
c     &  /'     dsL/dsAvg =',F7.3, '      dsR/dsAvg =',F7.3,
c     &  /'    curvature scaling exponent =',F8.4,
c     &  /'    curvature smoothing factor =',F8.4 )
c 3040 FORMAT(
c     &   '    station',I2,' refinement between s/smax =',2F8.4,
c     &   ' ,    ds/dsAvg =',F8.3)
C
c      WRITE(*,3100)
c      WRITE(*,3110) DLE/DSAVG,   DTE/DSAVG,
c     &            DSMIN/DSAVG, DSMAX/DSAVG
C
c 3100 FORMAT(/'  Actual resulting surface spacings ...')
c 3110 FORMAT( '     dsL/dsAvg =',F8.4,'     dsR/dsAvg =',F8.4,
c     &       /'   dsmin/dsAvg =',F8.4,'   dsmax/dsAvg =',F8.4 )
C
C
      IF(LQUERY) THEN
C----- go back to user-interaction menu
       GO TO 10
C
      ELSE
C----- we're done
       RETURN
C
      ENDIF
C
      END ! SGCURV



      SUBROUTINE SETCRV(SNEW,NPTS,SB,HH,N)
C------------------------------------------------------------------
C     Sets spacing array SNEW(i) based on surface curvature.
C
C     Input: NPTS    number of points defining surface
C            SB(.)   surface arc length
C            HH(.)   integrated spacing function
C            N       number of spacing array points to be generated
C
C     Output: SNEW(.)  spacing array to be generated
C             LFUDGL   T if curvature was augmented al left endpoint
C             LFUDGT   T if curvature was augmented al left endpoint
C------------------------------------------------------------------
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION SNEW(NPTS), SB(N), HH(N)
C
      PARAMETER (NMAX=1601)
      DIMENSION SBP(NMAX)
C
      IF(N.GT.NMAX) STOP 'SETCRV: array overflow'
C
      RN = FLOAT(NPTS-1)
C
C
C---- spline  HH(i)  =  h(s)  =  s + c(s)
      CALL SPLINA(HH,SBP,SB,N)
C
C---- invert derivative from dh/ds to ds/dh
      DO I=1, N
        SBP(I) = 1.0/SBP(I)
      ENDDO
C
C---- now,  SB(i),SBP(i)  =  s(h), s'(h)
C
C
      DH = HH(N)/RN
C
      SNEW(1) = 0.0
      DO I = 2, NPTS-1
        HHI = FLOAT(I-1)*DH
        SNEW(I) = SEVAL(HHI,SB,SBP,HH,N)
      ENDDO
      SNEW(NPTS) = SB(N)
C
      RETURN
      END ! SETCRV



      SUBROUTINE PLTCRV(S,CV0,CV,CVINT,N,CVAVG)
C---------------------------------------------------
C     Plots the contour, showing the node spacing
C---------------------------------------------------
C
      DIMENSION  S(N),CV0(N),CV(N),CVINT(N)
C
      I = 1
      CMAX = MAX( CV0(I) , CV(I) )
      DMAX = CVINT(I)
      DO I = 1, N
        CMAX = MAX( CMAX , CV0(I) , CV(I) )
        DMAX = MAX( DMAX , CVINT(I) )
      ENDDO
C
      CALL PLTINI2
      CALL PLOT(0.05,0.05,-3)
C
      XSIZ = 0.9
      YSIZ = 0.9 * 8.5/11.0
C
      SWT = XSIZ/(S(N)-S(1))
      CWT = YSIZ/CMAX
      DWT = YSIZ/DMAX
C
      CALL PLOT(0.0,0.0,3)
      CALL PLOT(1.0,0.0,2)
C
      CALL GETCOLOR(ICOL0)
C
      CALL NEWCOLORNAME('blue')
      IPEN = 3
      DO I = 1, N
        CALL PLOT(S(I)*SWT,CV0(I)*CWT,IPEN)
        IPEN = 2
      ENDDO
C
      CALL NEWCOLORNAME('cyan')
      IPEN = 3
      DO I=2, N
        CALL PLOT(S(I)*SWT,CV(I)*CWT,IPEN)
        IPEN = 2
      ENDDO
C
      CALL NEWCOLORNAME('orange')
      IPEN = 3
      DO I=2, N
        CALL PLOT(S(I)*SWT,CVINT(I)*DWT,IPEN)
        IPEN = 2
      ENDDO
C
      CALL NEWCOLORNAME('violet')
      CALL PLOT(S(1)*SWT,CVAVG*CWT,3)
      CALL PLOT(S(N)*SWT,CVAVG*CWT,2)
C
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
C
      RETURN
      END



      SUBROUTINE PLNODE(NB,XB,XPB,YB,YPB,SB,SG,NPTS,
     &                  XSIZ,YSIZ,SREF1,SREF2,NREF)
C
C---- Plot the airfoil, showing the node spacing
C
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION XB(NB),XPB(NB),
     &          YB(NB),YPB(NB), SB(NB)
      DIMENSION SG(NPTS)
      DIMENSION SREF1(NREF), SREF2(NREF)
C
      INCLUDE 'PLOT.INC'
      DATA XADD, YADD / 0.04, 0.08 /
C
      CALL PLTINI2
C
      XSIZ = XWIND/SIZE - 2.0*XADD
      YSIZ = YWIND/SIZE - 2.0*YADD
C
      XMIN = XB(1)
      YMIN = YB(1)
      XMAX = XB(1)
      YMAX = YB(1)
      DO I = 1, NB
        XMIN = MIN(XB(I),XMIN)
        XMAX = MAX(XB(I),XMAX)
        YMIN = MIN(YB(I),YMIN)
        YMAX = MAX(YB(I),YMAX)
      ENDDO
      SCL = 1.0 / MAX( (XMAX-XMIN)/XSIZ , (YMAX-YMIN)/YSIZ )
      XOFF = XMIN - XADD/SCL
      YOFF = YMIN - YADD/SCL
C
      CHRD = MAX( XMAX-XMIN , YMAX-YMIN )
      DSN = XSIZ * 0.01*CHRD
C
      SBTOT = SB(NB) - SB(1)
C
      CALL GETCOLOR(ICOL0)
C
C---- plot airfoil
      DO 302 IG=1, NPTS
        SBI = SB(1) + SG(IG)*SBTOT
        XX =  SEVAL(SBI,XB,XPB,SB,NB)
        YY =  SEVAL(SBI,YB,YPB,SB,NB)
        XP = SCL*(XX-XOFF)
        YP = SCL*(YY-YOFF)
C     
        IF(IG.EQ.1) THEN
         CALL PLOT(XP,YP,3)
         GO TO 302
        ENDIF
C
        DO IREF = 1, NREF
          IF(SREF1(IREF).GT.SG(IG-1) .AND. 
     &       SREF1(IREF).LE.SG(IG  )      ) THEN
           SBI = SB(1) + SREF1(IREF)*SBTOT
           XX =  SEVAL(SBI,XB,XPB,SB,NB)
           YY =  SEVAL(SBI,YB,YPB,SB,NB)
           XPR = SCL*(XX-XOFF)
           YPR = SCL*(YY-YOFF)
           CALL PLOT(XPR,YPR,2)
           CALL NEWCOLORNAME('magenta')
          ENDIF
C     
          IF(SREF2(IREF).GT.SG(IG-1) .AND. 
     &       SREF2(IREF).LE.SG(IG  )      ) THEN
           SBI = SB(1) + SREF2(IREF)*SBTOT
           XX =  SEVAL(SBI,XB,XPB,SB,NB)
           YY =  SEVAL(SBI,YB,YPB,SB,NB)
           XPR = SCL*(XX-XOFF)
           YPR = SCL*(YY-YOFF)
           CALL PLOT(XPR,YPR,2)
           CALL NEWCOLOR(ICOL0)
          ENDIF
        ENDDO
C     
        CALL PLOT(XP,YP,2)
 302  CONTINUE
C     
C     
C-----plot normal tick marks outside airfoil
      CALL NEWCOLORNAME('orange')
      DO 304 IG=1, NPTS
        SBI = SB(1) + SG(IG )*SBTOT
        XX =  SEVAL(SBI,XB,XPB,SB,NB)
        YY =  SEVAL(SBI,YB,YPB,SB,NB)
        XN =  DEVAL(SBI,YB,YPB,SB,NB)
        YN = -DEVAL(SBI,XB,XPB,SB,NB)
        XP = SCL*(XX-XOFF)
        YP = SCL*(YY-YOFF)
C     
        CALL PLOT(XP+SCL*DSN*XN,YP+SCL*DSN*YN,3)
        CALL PLOT(XP,YP,2)
 304  CONTINUE
C     
C-----fractional arc length tick marks inside element
      CALL NEWCOLORNAME('cyan')
      DO 306 IT=1, 9
        SFRAC = FLOAT(IT)/10.0
        SBI = SB(1) + SFRAC*SBTOT
        XX =  SEVAL(SBI,XB,XPB,SB,NB)
        YY =  SEVAL(SBI,YB,YPB,SB,NB)
        XN =  DEVAL(SBI,YB,YPB,SB,NB)
        YN = -DEVAL(SBI,XB,XPB,SB,NB)
        XP = SCL*(XX-XOFF)
        YP = SCL*(YY-YOFF)
        SCLF = -0.6*SCL
        IF(MOD(IT,5).EQ.0) SCLF = -1.2*SCL
        CALL PLOT(XP,YP,3)
        CALL PLOT(XP+SCLF*DSN*XN,YP+SCLF*DSN*YN,2)
 306  CONTINUE
C     
c      DO 308 IT=1, 9
c        SCLF = -0.15*SCL
c        IF(MOD(IT,5).EQ.0) SCLF = -0.30*SCL
cC
c        SFRAC = 0.1 * FLOAT(IT)/10.0
c        SBI = SB(1) + SFRAC*SBTOT
c        XX =  SEVAL(SBI,XB,XPB,SB,NB)
c        YY =  SEVAL(SBI,YB,YPB,SB,NB)
c        XN =  DEVAL(SBI,YB,YPB,SB,NB)
c        YN = -DEVAL(SBI,XB,XPB,SB,NB)
c        XP = SCL*(XX-XOFF)
c        YP = SCL*(YY-YOFF)
c        CALL PLOT(XP,YP,3)
c        CALL PLOT(XP+SCLF*DSN*XN,YP+SCLF*DSN*YN,2)
cC
c        SFRAC = 0.1 * FLOAT(IT)/10.0  +  0.9
c        SBI = SB(1) + SFRAC*SBTOT
c        XX =  SEVAL(SBI,XB,XPB,SB,NB)
c        YY =  SEVAL(SBI,YB,YPB,SB,NB)
c        XN =  DEVAL(SBI,YB,YPB,SB,NB)
c        YN = -DEVAL(SBI,XB,XPB,SB,NB)
c        XP = SCL*(XX-XOFF)
c        YP = SCL*(YY-YOFF)
c        CALL PLOT(XP,YP,3)
c        CALL PLOT(XP+SCLF*DSN*XN,YP+SCLF*DSN*YN,2)
c 308  CONTINUE
C     
      CALL NEWCOLOR(ICOL0)
C
      CALL PLFLUSH
      CALL PLEND
C
      RETURN
      END ! PLNODE

