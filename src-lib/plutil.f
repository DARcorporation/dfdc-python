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
C
C v070-ES2 
C PORTFAC scales portrait oriented windows
C
C=========================================================================

      SUBROUTINE PLTINI2
      INCLUDE 'PLOT.INC'
      LOGICAL LPLT,LLND
C
C---- terminate old plot if any
      IF(LPLOT) CALL PLEND
C
C---- initialize new plot
      IF(LLAND) THEN
        SIGNFR =  SCRNFR
      ELSE
        SIGNFR = -SCRNFR
      ENDIF
      CALL PLOPEN(SIGNFR,IPSLU,IDEV)
      LPLOT = .TRUE.
C
C---- set X-window size in inches (might have been resized by user)
      CALL GETWINSIZE(XWIND,YWIND)
C
C---- draw plot page outline offset by margins
      CALL NEWPEN(5)
      IF(XMARG .GT. 0.0) THEN
        CALL PLOTABS(      XMARG,      YMARG,3)
        CALL PLOTABS(      XMARG,YPAGE-YMARG,2)
        CALL PLOTABS(XPAGE-XMARG,      YMARG,3)
        CALL PLOTABS(XPAGE-XMARG,YPAGE-YMARG,2)
      ENDIF
      IF(YMARG .GT. 0.0) THEN
        CALL PLOTABS(      XMARG,      YMARG,3)
        CALL PLOTABS(XPAGE-XMARG,      YMARG,2)
        CALL PLOTABS(      XMARG,YPAGE-YMARG,3)
        CALL PLOTABS(XPAGE-XMARG,YPAGE-YMARG,2)
      ENDIF
      CALL NEWPEN(1)
C
      CALL PLOTABS(XMARG,YMARG,-3)
      CALL NEWCLIPABS( XMARG, XPAGE-XMARG, YMARG, YPAGE-YMARG )
C
      CALL NEWFACTOR(SIZE)
C
      RETURN
      END ! PLTINI2


      SUBROUTINE PLTINI(SCRN,IPSL,IDV,PSIZE,LPLT,LLND)
      INCLUDE 'PLOT.INC'
      LOGICAL LPLT,LLND
C
C---- terminate old plot if any
      IF(LPLOT) CALL PLEND
C
C
C---- portrait window scaling is applied here (v070-ES2)
C---- initialize new plot
      IF(LLND) THEN
        SIGNFR =  SCRNFR
      ELSE
        SIGNFR = -SCRNFR * PORTFAC
      ENDIF
      CALL PLOPEN(SIGNFR,IPSLU,IDEV)
      LPLOT = .TRUE.
      LPLT  = .TRUE.
C
C---- set X-window size in inches (might have been resized by user)
      CALL GETWINSIZE(XWIND,YWIND)
C
C---- draw plot page outline offset by margins
      CALL NEWPEN(5)
      IF(XMARG .GT. 0.0) THEN
        CALL PLOTABS(      XMARG,      YMARG,3)
        CALL PLOTABS(      XMARG,YPAGE-YMARG,2)
        CALL PLOTABS(XPAGE-XMARG,      YMARG,3)
        CALL PLOTABS(XPAGE-XMARG,YPAGE-YMARG,2)
      ENDIF
      IF(YMARG .GT. 0.0) THEN
        CALL PLOTABS(      XMARG,      YMARG,3)
        CALL PLOTABS(XPAGE-XMARG,      YMARG,2)
        CALL PLOTABS(      XMARG,YPAGE-YMARG,3)
        CALL PLOTABS(XPAGE-XMARG,YPAGE-YMARG,2)
      ENDIF
      CALL NEWPEN(1)
C
      CALL PLOTABS(XMARG,YMARG,-3)
      CALL CLRCLIP
c      CALL NEWCLIPABS( XMARG, XPAGE-XMARG, YMARG, YPAGE-YMARG )
C
      CALL NEWFACTOR(PSIZE)
C
      RETURN
      END ! PLTINI



      SUBROUTINE PLSUBS(XC,YC,CHX,STRING,ANGLE,NC,PLFONT)
C----------------------------------------------------------------
C     Plots character string as a subscript with font routine PLFONT.
C
C      XC,YC  = user coordinates of character to be subscripted 
C      CHX    = character width (user coordinates)
C      STRING = subscript character string to plot with NC characters
C      ANGLE  = angle of character (radians, positive is righthanded rotation)
C      NC     = number of subscript characters to plot
C               if NC<0 the length of the string is determined automatically 
C----------------------------------------------------------------
      CHARACTER*(*) STRING
      EXTERNAL PLFONT
      DATA  PI /3.1415926535897932384/
C
C---- subscript character reduction factor, and x,y-shift/chx
      DATA CHFAC, CHDX, CHDY / 0.7, 0.9, -0.4 /
C
      SINA = SIN(ANGLE*PI/180.0)
      COSA = COS(ANGLE*PI/180.0)
C
      XX = XC
      YY = YC
C
      IF (XC.EQ.999. .OR. YC.EQ.999.) THEN
        CALL GETLASTXY(XCHR,YCHR)
        IF(XC.EQ.999.) XX = XCHR
        IF(YC.EQ.999.) YY = YCHR
      ENDIF
C
      X = XX + CHX*(CHDX*COSA - CHDY*SINA)
      Y = YY + CHX*(CHDX*SINA + CHDY*COSA)
      CALL PLFONT(X,Y,CHX*CHFAC,STRING,ANGLE,NC)
C
      RETURN
      END



      SUBROUTINE PLSUPS(XC,YC,CHX,STRING,ANGLE,NC,PLFONT)
C----------------------------------------------------------------
C     Plots character string as a superscript with font routine PLFONT.
C
C      XC,YC  = user coordinates of character to be superscripted
C      CHX    = character width (user coordinates)
C      STRING = superscript character string to plot with NC characters
C      ANGLE  = angle of character (radians, positive is righthanded rotation)
C      NC     = number of superscript characters to plot
C               if NC<0 the length of the string is determined automatically 
C----------------------------------------------------------------
      CHARACTER*(*) STRING
      EXTERNAL PLFONT
      DATA  PI /3.1415926535897932384/
C
C---- superscript character reduction factor, and x,y-shift/chx
      DATA CHFAC, CHDX, CHDY / 0.7, 0.95, 0.7 /
C
      SINA = SIN(ANGLE*PI/180.0)
      COSA = COS(ANGLE*PI/180.0)
C
      XX = XC
      YY = YC
C
      IF (XC.EQ.999. .OR. YC.EQ.999.) THEN
        CALL GETLASTXY(XCHR,YCHR)
        IF(XC.EQ.999.) XX = XCHR
        IF(YC.EQ.999.) YY = YCHR
      ENDIF
C
      X = XX + CHX*(CHDX*COSA - CHDY*SINA)
      Y = YY + CHX*(CHDX*SINA + CHDY*COSA)
      CALL PLFONT(X,Y,CHX*CHFAC,STRING,ANGLE,NC)
C
      RETURN
      END



      SUBROUTINE PLNSIG(X,Y,CH,RNUM,ANG,NSIG)
C---------------------------------------------------
C     Same as PLNUMB, but sets number of digits 
C     after decimal point to give at least NSIG 
C     significant digits in plotted number RNUM.
C---------------------------------------------------
C
      IF(RNUM.EQ.0.0) THEN
       NDIG = 1
C
      ELSE
       GNUM = LOG10( ABS(RNUM) )
C
C----- number of insignificant zeros after decimal point
       NZER = INT(-GNUM+100.0) - 100
C
C----- set required number of digits after decimal point
       NDIG = MAX( 0 , NZER+NSIG )
C
      ENDIF
C
      CALL PLNUMB(X,Y,CH,RNUM,ANG,NDIG)
C
      RETURN
      END



      SUBROUTINE SPPLOT(N,X,XS,S,SOFF,XOFF,SWT,XWT,KK)
      DIMENSION X(N),XS(N),S(N)
C------------------------------------------------
C     Plots X(S) using spline evaluation,
C     with KK sub-intervals per spline interval.
C------------------------------------------------
C
      SPLT(SS) = (SS-SOFF)*SWT
      XPLT(XX) = (XX-XOFF)*XWT
C
      I = 1
      CALL PLOT(SPLT(S(I)),XPLT(X(I)),3)
      DO I=2, N
        DO K=1, KK-1
          SK = S(I-1) + (S(I)-S(I-1))*FLOAT(K)/FLOAT(KK)
          XK = SEVAL(SK,X,XS,S,N)
          CALL PLOT(SPLT(SK),XPLT(XK),2)
        ENDDO
        CALL PLOT(SPLT(S(I)),XPLT(X(I)),2)
      ENDDO
C
      RETURN
      END



      SUBROUTINE XYDASH(N,X,Y,XOFF,XFAC,YOFF,YFAC,FRAC)
      DIMENSION X(*),Y(*)
C
      GAP = 0.5*(1.0 - FRAC)
      DO I = 1, N-1
        DX = X(I+1) - X(I)
        DY = Y(I+1) - Y(I)
C
        X1 = X(I)   + GAP*DX
        Y1 = Y(I)   + GAP*DY
        X2 = X(I+1) - GAP*DX
        Y2 = Y(I+1) - GAP*DY
C
        CALL PLOT((X1-XOFF)*XFAC,(Y1-YOFF)*YFAC,3)
        CALL PLOT((X2-XOFF)*XFAC,(Y2-YOFF)*YFAC,2)
      ENDDO
C
      RETURN
      END ! XYDASH



      SUBROUTINE SETLIM(N,Y,ANNMIN,ANNMAX,DANN,NANN)
      DIMENSION Y(N)
C--------------------------------------------
C     Sets axis plotting limits for Y array,
C     with special treatment for zero array,
C     and separate checking for positive
C     and negative axis sides.
C--------------------------------------------
C
      YMIN = 0.0
      YMAX = 0.0
      DO I=1, N
        YMIN = MIN( YMIN , Y(I) )
        YMAX = MAX( YMAX , Y(I) )
      ENDDO
C
      IF(YMAX .EQ. YMIN) THEN
       YMIN = -1.0
       YMAX =  1.0
       DANN = 0.5
       NANN = 2
      ELSEIF(ABS(YMAX) .GT. ABS(YMIN)) THEN
       CALL SCALIT2(1,YMAX,0.0,SF,NANN)
       DANN = 1.0 / (SF*FLOAT(NANN))
       IF(YMIN/YMAX .GT. 0.9) DANN = 2.0*DANN
      ELSE
       CALL SCALIT2(1,YMIN,0.0,SF,NANN)
       DANN = 1.0 / (SF*FLOAT(NANN))
       IF(YMAX/YMIN .GT. 0.9) DANN = 2.0*DANN
      ENDIF       
C
      NANMAX = INT(YMAX/DANN + 1000.99) - 1000
      NANMIN = INT(YMIN/DANN + 1000.01) - 1000
C
      NANN = NANMAX - NANMIN
      IF(NANN .GT. 6) THEN
        DANN = 2.0*DANN
        NANMAX = INT(YMAX/DANN + 1000.99) - 1000
        NANMIN = INT(YMIN/DANN + 1000.01) - 1000
      ENDIF
C
      ANNMAX = DANN * FLOAT(NANMAX)
      ANNMIN = DANN * FLOAT(NANMIN)
C
      RETURN
      END



      SUBROUTINE SCALIT2(N,Y,YOFF,YSF,NANN)
      DIMENSION Y(*)
C.............................................................
C
C     Determines scaling factor for the offset Y array so that
C     YSF*(Ymax-YOFF) will be O(1.0), but less than 1.0.
C
C     ANN = 1.0/YSF is therefore a "nice" plot axis max annotation.
C
C     Y(1:N)   array whose scaling factor is to be determined
C     YOFF     offset of Y array  (Y-YOFF is actually scaled)
C     YSF      Y scaling factor
C     NANN     recommended number of Y annotation intervals
C.............................................................
C
      AGH = ALOG10(1.5)
      AG2 = ALOG10(2.0)
      AG3 = ALOG10(3.0)
      AG4 = ALOG10(4.0)
      AG5 = ALOG10(5.0)
      AG8 = ALOG10(8.0)
C
      YMAX = ABS(Y(1)-YOFF)
      DO I=2, N
        YMAX = MAX( YMAX , ABS(Y(I)-YOFF) )
      ENDDO
C
      IF(YMAX.EQ.0.0) THEN
        YSF = 1.0
        NANN = 5
        RETURN
      ENDIF
C
      YLOG = ALOG10(YMAX) - 0.001
C
C---- find log of nearest power of 10 above YMAX
      YLOG1 = AINT(YLOG+100.0) - 99.0
 
C---- find log of nearest 1.5x(power of 10) above YMAX
      YLOGH = YLOG1 + AGH
      IF(YLOGH-1.0.GT.YLOG) YLOGH = YLOGH - 1.0
C
C---- find log of nearest 2x(power of 10) above YMAX
      YLOG2 = YLOG1 + AG2
      IF(YLOG2-1.0.GT.YLOG) YLOG2 = YLOG2 - 1.0
C
C---- find log of nearest 3x(power of 10) above YMAX
      YLOG3 = YLOG1 + AG3
      IF(YLOG3-1.0.GT.YLOG) YLOG3 = YLOG3 - 1.0
C
C---- find log of nearest 4x(power of 10) above YMAX
      YLOG4 = YLOG1 + AG4
      IF(YLOG4-1.0.GT.YLOG) YLOG4 = YLOG4 - 1.0
C
C---- find log of nearest 5x(power of 10) above YMAX
      YLOG5 = YLOG1 + AG5
      IF(YLOG5-1.0.GT.YLOG) YLOG5 = YLOG5 - 1.0
C
C---- find log of nearest 8x(power of 10) above YMAX
      YLOG8 = YLOG1 + AG8
      IF(YLOG8-1.0.GT.YLOG) YLOG8 = YLOG8 - 1.0
C
C---- find log of smallest upper bound
      GMIN = MIN( YLOG1 , YLOGH, YLOG2 , YLOG3, YLOG4 , YLOG5 , YLOG8 )
C
      NANN = 5
      IF (GMIN.EQ.YLOGH .OR. GMIN.EQ.YLOG3) NANN = 3
      IF (GMIN.EQ.YLOG2 .OR. GMIN.EQ.YLOG4 .OR. GMIN.EQ.YLOG8) NANN = 4
C
C---- set scaling factor
      YSF = 10.0**(-GMIN)
C
      RETURN
      END ! SCALIT2



      SUBROUTINE SCALIT(N,Y,YOFF,YSF)
      DIMENSION Y(*)
C-------------------------------------------------------------
C     Y(1:N)  array whose scaling factor is to be determined
C     YOFF    offset of Y array  (Y-YOFF is actually scaled)
C     YSF     Y scaling factor
C-------------------------------------------------------------
C
      AG2 = LOG10(2.0)
      AG5 = LOG10(5.0)
C
      YMAX = ABS(Y(1) - YOFF)
      DO I=2, N
        YMAX = MAX( YMAX , ABS(Y(I)-YOFF) )
      ENDDO
C
      IF(YMAX.EQ.0.0) THEN
        YSF = 1.0
        NANN = 5
        RETURN
      ENDIF
C
      YLOG = LOG10(YMAX)
C
C---- find log of nearest power of 10 above YMAX
      YLOG1 = AINT(YLOG+100.0) - 99.0
 
C---- find log of nearest 2x(power of 10) above YMAX
      YLOG2 = YLOG1 + AG2
      IF(YLOG2-1.0.GT.YLOG) YLOG2 = YLOG2 - 1.0
C
C---- find log of nearest 5x(power of 10) above YMAX
      YLOG5 = YLOG1 + AG5
      IF(YLOG5-1.0.GT.YLOG) YLOG5 = YLOG5 - 1.0
C
C---- find log of smallest upper bound
      GMIN = MIN( YLOG1 , YLOG2 , YLOG5 )
C
C---- set scaling factor
      YSF = 10.0**(-GMIN)
C
      RETURN
      END


      SUBROUTINE AXISADJ2(xmin,xmax,xspan,deltax,ntics)
C...Make scaled axes with engineering increments between tics
C
C   Input:    xmin, xmax - input range for which scaled axis is desired
C
C   Output:   xmin, xmax - adjusted range for scaled axis
C             xspan      - adjusted span of scaled axis
C             deltax     - increment to be used for scaled axis
C             ntics      - number of axis tics (each deltax long)
C
      real    xmin,xmax,xspan,deltax,xinc,xinctbl(4)
      integer ntics,i
      data    xinctbl / 0.1, 0.2, 0.5, 1.0 /
c
      xspan1 = xmax-xmin
      if (xspan1.eq.0.) xspan1 = 1.
c
      xpon = ifix(log10(xspan1))
      xspan = xspan1 / 10.**xpon
c
      do i = 1, 4
        xinc = xinctbl(i)
        ntics = 1 + ifix(xspan/xinc)
        if (ntics.LE.12) go to 1
      end do
c
   1  deltax = xinc*10.**xpon
      xmin = deltax*  ifloor(xmin/deltax)
      xmax = deltax*iceiling(xmax/deltax)
      xspan = xmax - xmin
      return
      end



      SUBROUTINE OFFGET(XOFF,YOFF,XSF,YSF,XWIND,YWIND,LSAME,LCURS,SH)
      LOGICAL LSAME, LCURS
      CHARACTER*1 KCHAR
C---------------------------------------------------
C     Sets new blowup parameters from cursor input.
C---------------------------------------------------
C
C---- get current color
      CALL GETCOLOR(ICOL0)
C
C---- set new crosshair color
      CALL NEWCOLORNAME('MAGENTA')
C
C
      IF(LCURS) THEN
       WRITE(*,*)
       WRITE(*,*) 'Mark off corners of blowup area'
       WRITE(*,*) '(2 identical points default to current area)'
C
       CALL GETCURSORXY(XX1,YY1,KCHAR)
       CALL PLSYMB(XX1,YY1,SH,3,0.0,0)
       WRITE(*,*) 'x,y =', XX1/XSF+XOFF, YY1/YSF+YOFF
C
       CALL GETCURSORXY(XX2,YY2,KCHAR)
       CALL PLSYMB(XX2,YY2,SH,3,0.0,0)
       WRITE(*,*) 'x,y =', XX2/XSF+XOFF, YY2/YSF+YOFF
C
      ELSE
       WRITE(*,*)
       WRITE(*,*) 'Enter x,y coordinates of blowup area corners'
       WRITE(*,*) '(2 identical points default to current area)'
       WRITE(*,*)
    1  WRITE(*,*) 'Point 1: '
       READ(*,*,ERR=1) XX1, YY1
    2  WRITE(*,*) 'Point 2: '
       READ(*,*,ERR=2) XX2, YY2
C
      ENDIF
C
C---- restore to initial color
      CALL NEWCOLOR(ICOL0)
C
      IF(XX1.EQ.XX2 .AND. YY1.EQ.YY2) RETURN
C
C---- note negative scale factors
      XSGN = SIGN(1.0,XSF)
      YSGN = SIGN(1.0,YSF)
C
C---- center and dimensions of selected box in current plot units
      XCEN = 0.5*(XX1+XX2)/XSF + XOFF
      YCEN = 0.5*(YY1+YY2)/YSF + YOFF
      XDIF = ABS((XX2 - XX1)/XSF)
      YDIF = ABS((YY2 - YY1)/YSF)
C
      IF(XDIF.EQ.0.0) XDIF = 1.0E-5
      IF(YDIF.EQ.0.0) YDIF = 1.0E-5
C
C---- set new scales,offsets
      XOFF = MIN(XX1,XX2)/XSF + XOFF
      YOFF = MIN(YY1,YY2)/YSF + YOFF
      XSF = XSGN*XWIND/XDIF
      YSF = YSGN*YWIND/YDIF
C
      IF(LSAME) THEN
C----- set equal x,y scales
       SF = MIN( XSGN*XSF , YSGN*YSF )
       XSF = XSGN*SF
       YSF = YSGN*SF
C
C----- re-center the blowup
       XOFF = XCEN - 0.5*XDIF
       YOFF = YCEN - 0.5*YDIF
      ENDIF
C
      RETURN
      END ! OFFGET



      SUBROUTINE DIST(XOFF,YOFF,XSF,YSF,SH)
C--------------------------------------------------
C     Displays distance between two cursor points.
C--------------------------------------------------
      CHARACTER*1 KCHAR
C
      WRITE(*,*)
      WRITE(*,*) 'Click mouse or hit a key on each point'
      WRITE(*,*)
      CALL GETCURSORXY(XX1,YY1,KCHAR)
      CALL PLSYMB(XX1,YY1,SH,3,0.0,0)
      CALL PLFLUSH
      XX1 = XX1/XSF + XOFF
      YY1 = YY1/YSF + YOFF
      WRITE(*,1010) XX1,YY1
C
      CALL GETCURSORXY(XX2,YY2,KCHAR)
      CALL PLSYMB(XX2,YY2,SH,3,0.0,0)
      CALL PLFLUSH
      XX2 = XX2/XSF + XOFF
      YY2 = YY2/YSF + YOFF
      WRITE(*,1020) XX2,YY2
C
      DX = XX2 - XX1
      DY = YY2 - YY1
      DS = SQRT(DX*DX + DY*DY)
      WRITE(*,1050) DX, DY, DS
C
 1010 FORMAT(' x1 =', F10.6, '    y1 =', F10.6)
 1020 FORMAT(' x2 =', F10.6, '    y2 =', F10.6)
 1050 FORMAT(' dx =', F10.6, '    dy =', F10.6,'    ds =', F10.6)
C
      RETURN
      END




      SUBROUTINE PLTGSP(X,XS,Y,YS,S,N, XOFF,XSF,YOFF,YSF, KK)
      DIMENSION X(*), XS(*), Y(*), YS(*), S(*)
C----------------------------------------
C     Plots passed-in X(S),Y(S) contour
C     with KK spline sub-intervals
C----------------------------------------
      XMOD(XTMP) = XSF * (XTMP - XOFF)
      YMOD(YTMP) = YSF * (YTMP - YOFF)
C
      I = 1
      CALL PLOT(XMOD(X(I)),YMOD(Y(I)),3)
C
      DO I = 2, N
        DS = S(I) - S(I-1)
        DO K = 1, KK-1
          SK = S(I-1) + DS*FLOAT(K)/FLOAT(KK)
          XK = SEVAL(SK,X,XS,S,N)
          YK = SEVAL(SK,Y,YS,S,N)
          CALL PLOT(XMOD(XK),YMOD(YK),2)
        ENDDO
C
        CALL PLOT(XMOD(X(I)),YMOD(Y(I)),2)
      ENDDO
C
      RETURN
      END


 

      SUBROUTINE ARROW(X,Y,DX,DY)
      CALL PLOT(X,Y,3)
      CALL PLOT(X+DX,Y+DY,2)
      X1 = X + 0.85*DX + 0.02*DY
      Y1 = Y + 0.85*DY - 0.02*DX
      X2 = X + 0.85*DX - 0.02*DY
      Y2 = Y + 0.85*DY + 0.02*DX
      CALL PLOT(X1,Y1,2)
      CALL PLOT(X2,Y2,2)
      CALL PLOT(X+DX,Y+DY,2)
      RETURN
      END
 
 
      SUBROUTINE YDASH(X1,X2,Y)
      CALL NEWPEN(1)
      DX = (X2-X1)/50.0
      DO I = 1, 51
        X = X1 + DX*FLOAT(I-1)
        CALL PLOT(X-0.08*DX,Y,3)
        CALL PLOT(X+0.08*DX,Y,2)
      ENDDO
      RETURN
      END


      SUBROUTINE BOX(X,Y,IWID)
      DIMENSION X(2), Y(2)
C
      CALL NEWPEN(IWID)
      CALL PLOT(X(1),Y(1),3)
      CALL PLOT(X(2),Y(1),2)
      CALL PLOT(X(2),Y(2),2)
      CALL PLOT(X(1),Y(2),2)
      CALL PLOT(X(1),Y(1),2)
      RETURN
      END


      FUNCTION LABORT(XC,YC)
      LOGICAL LABORT
      COMMON /COM_ABORT/ XABORT(2), YABORT(2)
C
C---- return T if location XC,YC falls within abort window
C
      LABORT = XC .GE. XABORT(1) .AND. 
     &         XC .LE. XABORT(2) .AND.
     &         YC .GE. YABORT(1) .AND.
     &         YC .LE. YABORT(2)
C
      RETURN
      END


      SUBROUTINE PLCIRC(X0,Y0,RAD, NSEG)
      DTR = ATAN(1.0) / 45.0
C
      IF(NSEG.EQ.0) THEN
C------ draw solid circle
        CALL PLOT(X0+RAD,Y0,3)
        DO I=1, 360
          T = FLOAT(I) * DTR
          CALL PLOT(X0+RAD*COS(T),Y0+RAD*SIN(T),2)
        ENDDO
      ELSE
C------ draw dashed circle in NSEG segments
        KSEG = 360/NSEG
        IGAP = MAX( KSEG/5 , 1 )
        DO ISEG=1, NSEG
          I1 = (ISEG-1)*KSEG  + IGAP
          I2 =  ISEG   *KSEG  - IGAP
          IPEN = 3
          DO I=I1, I2
            T = FLOAT(I) * DTR
            CALL PLOT(X0+RAD*COS(T),Y0+RAD*SIN(T),IPEN)
            IPEN = 2
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END



      SUBROUTINE PABORT
      COMMON /COM_ABORT/ XABORT(2), YABORT(2)
C
      CALL GETWINSIZE(XWIND,YWIND)
      CALL GETORIGIN(XORG,YORG)
      CALL GETFACTORS(XSCALE,YSCALE)
C
C---- save and disable current clipping
      CALL GETCLIPABS(XMIN,XMAX,YMIN,YMAX)
      CALL CLRCLIP
C
      CALL GETCOLOR(ICOL0)
      CALL NEWCOLORNAME('red')
C
C---- set abort window in lower right corner
      XABORT(1) = (XWIND - 1.0 - XORG)/XSCALE
      XABORT(2) = (XWIND - 0.1 - XORG)/XSCALE
      YABORT(1) = (        0.1 - YORG)/YSCALE
      YABORT(2) = (        0.5 - YORG)/YSCALE
C
C---- plot abort window
      CALL PLOT(XABORT(1),YABORT(1),3)
      CALL PLOT(XABORT(2),YABORT(1),2)
      CALL PLOT(XABORT(2),YABORT(2),2)
      CALL PLOT(XABORT(1),YABORT(2),2)
      CALL PLOT(XABORT(1),YABORT(1),2)
C
      CHA = MIN( (XABORT(2)-XABORT(1))/8.0 , (YABORT(2)-YABORT(1))/1.5 )
      XCA = 0.5*(XABORT(2)+XABORT(1)) - 2.5*CHA
      YCA = 0.5*(YABORT(2)+YABORT(1)) - 0.5*CHA
      CALL PLCHAR(XCA,YCA,CHA,'ABORT',0.0,5)
C
C---- restore color and clipping
      CALL NEWCOLOR(ICOL0)
      CALL NEWCLIPABS(XMIN,XMAX,YMIN,YMAX)
C
      RETURN
      END
C
C---------------------------------------------------------------------
