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


      SUBROUTINE MODIFY(IX,IFRST,ILAST,NSIDE,NLINE,
     &                  X,Y,YD, LBLEND,
     &                  IMOD1,IMOD2,ISMOD,ILMOD,
     &                  XMOD,YMOD, XWIN,YWIN,
     &                  XOFF,YOFF,XSF,YSF, LMODPL)
      DIMENSION IFRST(NSIDE), ILAST(NSIDE)
      DIMENSION X(IX), Y(IX,NLINE), YD(IX,NLINE)
      DIMENSION XMOD(2),YMOD(2), XWIN(2),YWIN(2)
      LOGICAL LBLEND, LMODPL
C--------------------------------------------------------------------------
C     Allows user to modify functions Y1(X),Y2(X)... via cursor input.
C
C     Cursor-specified Xu,Yu values are sorted by Xu and splined.
C     The resulting spline function Yu(X) is interrogated at input 
C     X(i) points to obtain the modified Y(i,L) values.
C
C  Input:   IX          first dimension of X,Y arrays
C           IFRST(s)    first i index in segment s
C           ILAST(s)    last  i index in segment s
C           NSIDE       number of X segments :   s = 1..NSIDE
C           NLINE       number of Y functions:   l = 1..NLINE
C           X(i)        X values
C           Y(i,l)      Y values
C           YD(i,l)     spline derivative array dY/dX (used only if LSLOPE=T)
C           LBLEND      if T, blends input Yu(Xu) with Y(X) at input endpoints
C           XMOD(2)     x-limits of box for cursor input
C           YMOD(2)     y-limits of box for cursor input
C           XOFF        plot offsets,scales used to plot Y(X)
C           YOFF         "
C           XSF          "
C           YSF          "
C           LMODPL      if T, plot modified Y(i,l) points
C
C  Output:  Y(i,l)      modified Y values
C           IMOD1       first i index of modified Y(i,l) values
C           IMOD2       last  i index of modified Y(i,l) values
C           ISMOD       index s of segment containing IMOD1,IMOD2
C           ILMOD       index l of Y(i,l) function which was modified
C--------------------------------------------------------------------------
C
C---- local arrays for accumulating user-specified points
      PARAMETER (NUX=100)
      DIMENSION XU(NUX), YU(NUX), YUD(NUX)
      DIMENSION IUSORT(NUX)
      LOGICAL LDONE
C
      LOGICAL LGUI
      CHARACTER*1 CHKEY
C
      DATA SH /0.010/
C
      CALL GETCOLOR(ICOL0)
      CALL GETPEN(IPEN0)
C
      KDONE  = 1
      KERASE = 2
      KABORT = 3
      KINSIDE = 4
C
      XDWIN = XWIN(2) - XWIN(1)
      YDWIN = YWIN(2) - YWIN(1)
C
      CALL NEWPEN(5)
C
      X1 = XWIN(1) + 0.71*XDWIN
      X2 = XWIN(1) + 0.79*XDWIN
      Y1 = YWIN(1) + 0.01*YDWIN
      Y2 = YWIN(1) + 0.05*YDWIN
      CALL GUIBOX(KABORT, X1,X2,Y1,Y2, 'RED'   , ' Abort ')
C
      X1 = XWIN(1) + 0.81*XDWIN
      X2 = XWIN(1) + 0.89*XDWIN
      Y1 = YWIN(1) + 0.01*YDWIN
      Y2 = YWIN(1) + 0.05*YDWIN
      CALL GUIBOX(KERASE, X1,X2,Y1,Y2, 'YELLOW', ' Erase ')
C
      X1 = XWIN(1) + 0.91*XDWIN
      X2 = XWIN(1) + 0.99*XDWIN
      Y1 = YWIN(1) + 0.01*YDWIN
      Y2 = YWIN(1) + 0.05*YDWIN
      CALL GUIBOX(KDONE , X1,X2,Y1,Y2, 'GREEN', ' Done ')
C
      X1 = XMOD(1)
      X2 = XMOD(2)
      Y1 = YMOD(1)
      Y2 = YMOD(2)
      CALL GUIBOX(KINSIDE, X1,X2,Y1,Y2, 'ORANGE' , ' ' )
C
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      CALL NEWPEN(IPEN0)
C
      WRITE(*,*)
      WRITE(*,*) 'Click on new values...'
      WRITE(*,*)
C
      NUBEG = 1
C
 5    CONTINUE
      CALL NEWCOLORNAME('MAGENTA')
      DO NU = NUBEG, NUX
C
C------ fetch x-y point coordinates from user
        CALL GETCURSORXY(XU(NU),YU(NU),CHKEY)
C
        IF    (LGUI(KABORT,XU(NU),YU(NU))
     &         .OR. INDEX('Aa',CHKEY).GT.0) THEN
C------- return with no changes
         GO TO 90
C
        ELSEIF(LGUI(KERASE,XU(NU),YU(NU))
     &         .OR. INDEX('Ee',CHKEY).GT.0) THEN
         IF(NU.LE.1) THEN
          WRITE(*,*) 'No more points to clear'
          NUBEG = 1
         ELSE
C-------- clear previous point, overplot it white to clear it from screen
          NUBEG = NU - 1
          CALL NEWCOLORNAME('WHITE')
          CALL PLSYMB(XU(NUBEG),YU(NUBEG),SH,3,0.0,0)
         ENDIF
C
C------- keep accepting points starting from NUBEG
         GO TO 5
C
        ELSEIF(LGUI(KDONE,XU(NU),YU(NU))
     &         .OR. INDEX('Dd',CHKEY).GT.0) THEN
C------- go process inputs
         GO TO 10
C
        ELSEIF(LGUI(KINSIDE,XU(NU),YU(NU))) THEN
C------- normal click inside modify-window: plot small cross at input point
         CALL PLSYMB(XU(NU),YU(NU),SH,3,0.0,0)
C
        ELSE
C------- must be somewhere outside
         GO TO 10
C
        ENDIF
      ENDDO
      WRITE(*,*) 'MODIFY: User-input array limit NUX reached'
C
C---- pick up here when finished with input
 10   CONTINUE
      IF(INDEX('Dd',CHKEY).GT.0) THEN
C----- last point was entered with a "D" ...  add it to list
       CALL PLSYMB(XU(NU),YU(NU),SH,3,0.0,0)
      ELSE
C----- discard last point
       NU = NU-1
      ENDIF
C
      CALL PLFLUSH
C
      IF(NU.LT.2) THEN
       WRITE(*,*)
       WRITE(*,*) 'Need at least 2 points'
       GO TO 90
      ENDIF
C
C---- set first-specified point
      XUSP1 = XU(1)
      YUSP1 = YU(1)
C
C---- undo plot offsets and scales
      DO IU = 1, NU
        XU(IU) = XU(IU)/XSF + XOFF
        YU(IU) = YU(IU)/YSF + YOFF
      ENDDO
C
C---- sort XU,YU points in XU (use spline array YUD as temporary storage)
      CALL HSORT(NU,XU,IUSORT)
C
      DO KSORT = 1, NU
        IU = IUSORT(KSORT)
        YUD(KSORT) = XU(IU)
      ENDDO
      DO IU = 1, NU
        XU(IU) = YUD(IU)
      ENDDO
C
      DO KSORT = 1, NU
        IU = IUSORT(KSORT)
        YUD(KSORT) = YU(IU)
      ENDDO
      DO IU = 1, NU
        YU(IU) = YUD(IU)
      ENDDO
C
C---- remove doubled endpoints and tripled interior points
      DO IPASS = 1, 12345
        LDONE = .TRUE.
        IU = 2
        IF(XU(IU).EQ.XU(IU-1)) THEN
         LDONE = .FALSE.
         IUREM = IU
        ENDIF
        DO IU = 3, NU
          IF( XU(IU).EQ.XU(IU-1) .AND.
     &        XU(IU).EQ.XU(IU-2)      ) THEN
           LDONE = .FALSE.
           IUREM = IU
          ENDIF
        ENDDO
        IU = NU
        IF(XU(IU).EQ.XU(IU-1)) THEN
         LDONE = .FALSE.
         IUREM = IU
        ENDIF
C
        IF(LDONE) THEN
         GO TO 30
        ELSE
         DO IU = IUREM, NU-1
           XU(IU) = XU(IU+1)
           YU(IU) = YU(IU+1)
         ENDDO
         NU = NU - 1
        ENDIF
      ENDDO
C
C---- pick up here when no more points to be removed
 30   CONTINUE
      IF(NU.LT.2) THEN
       WRITE(*,*)
       WRITE(*,*) 'Need at least 2 points'
       GO TO 90
      ENDIF
C
C
C---- find which X,Y input point is closest to first-specified point
      ISMOD = 1
      ILMOD = 1
C
C---- go over all surface points
      DSQMIN = 1.0E24
      DO IL = 1, NLINE
        DO IS = 1, NSIDE
          DO I = IFRST(IS), ILAST(IS)
C---------- convert input arrays to plot coordinates
            XUI = (X(I   )-XOFF)*XSF
            YUI = (Y(I,IL)-YOFF)*YSF
            DSQ = (XUI-XUSP1)**2 + (YUI-YUSP1)**2
C
            IF(DSQ .LT. DSQMIN) THEN
C------------ this point is the closest so far... note its indices
              DSQMIN = DSQ
              ISMOD = IS
              ILMOD = IL
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
C---- set side and function to be modified
      IS = ISMOD
      IL = ILMOD
C
      IF(LBLEND) THEN
C----- reset Y and dY/dX at first and last points of modified interval
       X1 = X(IFRST(IS))
       X2 = X(ILAST(IS))
       I = IFRST(IS)
       N = ILAST(IS) - IFRST(IS) + 1
C
       IU = 1
       IF(XU(IU).GE.X1 .AND. XU(IU).LE.X2) THEN
C------ set function and derivative at left endpoint
        YU(IU) = SEVAL(XU(IU),Y(I,IL),YD(I,IL),X(I),N)
        YD1    = DEVAL(XU(IU),Y(I,IL),YD(I,IL),X(I),N)
       ELSE
        YD1 = -999.0
       ENDIF
C
       IU = NU
       IF(XU(IU).GE.X1 .AND. XU(IU).LE.X2) THEN
        YU(IU) = SEVAL(XU(IU),Y(I,IL),YD(I,IL),X(I),N)
        YD2    = DEVAL(XU(IU),Y(I,IL),YD(I,IL),X(I),N)
       ELSE
        YD2 = -999.0
       ENDIF
C
      ELSE
C----- use natural spline end conditions (zero 3rd derivative)
       YD1 = -999.0
       YD2 = -999.0
C
      ENDIF
C
C---- spline input function values
      CALL SPLIND(YU,YUD,XU,NU,YD1,YD2)
C
C
C---- go over all points on modified segment
      IMOD1 = IFRST(IS)
      DO I = IFRST(IS), ILAST(IS)
        XI = X(I)
C
        IF    (XI .LT. XU( 1)) THEN
C------- current point is before modified interval...try next point
         IMOD1 = I
        ELSEIF(XI .LE. XU(NU)) THEN
C------- stuff new point into Vspec array and plot it
         Y(I,IL) = SEVAL(XI,YU,YUD,XU,NU)
        ELSE
C------- went past modified interval...finish up
         IMOD2 = I
         GO TO 50
        ENDIF
      ENDDO
      IMOD2 = ILAST(IS)
 50   CONTINUE
C
C
      IF(LMODPL) THEN
C----- plot modified function over modified interval
       CALL NEWCOLORNAME('MAGENTA')
       IPEN = 3
       DO I = IMOD1, IMOD2
         XP = (X(I   )-XOFF)*XSF
         YP = (Y(I,IL)-YOFF)*YSF
         CALL PLOT(XP,YP,IPEN)
         IPEN = 2
       ENDDO
       CALL PLFLUSH
      ENDIF
C
C---- return normally
      CALL NEWCOLOR(ICOL0)
      RETURN
C
C-------------------------------------------------
 90   CONTINUE
      WRITE(*,*) 'No changes made'
      IMOD1 = IFRST(1)
      IMOD2 = IFRST(1) - 1
      ISMOD = 1
      ILMOD = 1
      CALL NEWCOLOR(ICOL0)
      RETURN
C
      END ! MODIFY



      SUBROUTINE MODIXY(IX,IFRST,ILAST,NSIDE,
     &                  X,Y,XD,YD,S, LBLEND,
     &                  IMOD1,IMOD2,ISMOD,
     &                  XMOD,YMOD, XWIN,YWIN,
     &                  XOFF,YOFF,XSF,YSF, LMODPL)
      DIMENSION IFRST(NSIDE), ILAST(NSIDE)
      DIMENSION X(IX),Y(IX), XD(IX),YD(IX), S(IX)
      DIMENSION XMOD(2),YMOD(2), XWIN(2),YWIN(2)
      LOGICAL LBLEND, LMODPL
C--------------------------------------------------------------------------
C     Allows user to modify contours X(S),Y(S) via cursor input.
C
C     Cursor-specified Xu,Yu values are splined in Su.
C     The resulting spline functions Xu(Su),Yu(Su) are interrogated 
C     at input S(i) points to obtain the modified X(i),Y(i) values.
C
C  Input:   IX          first dimension of X,Y arrays
C           IFRST(s)    first i index in segment s
C           ILAST(s)    last  i index in segment s
C           NSIDE       number of X segments :   s = 1..NSIDE
C           X(i)        X values
C           Y(i)        Y values
C           XD(i)       spline derivative array dX/dS (used only if LSLOPE=T)
C           YD(i)       spline derivative array dY/dS (used only if LSLOPE=T)
C           S(i)        S values
C           LBLEND      if T, blends input Yu(Xu) with Y(X) at input endpoints
C           XMOD(2)     x-limits of box for cursor input
C           YMOD(2)     y-limits of box for cursor input
C           XOFF        plot offsets,scales used to plot X(S),Y(S)
C           YOFF         "
C           XSF          "
C           YSF          "
C           LMODPL      if T, plot modified X(i),Y(i) points
C
C  Output:  X(i)        modified X values
C           Y(i)        modified Y values
C           IMOD1       first i index of modified X(i),Y(i) values
C           IMOD2       last  i index of modified X(i),Y(i) values
C           ISMOD       index s of segment containing IMOD1,IMOD2
C--------------------------------------------------------------------------
C
C---- local arrays for accumulating user-specified points
      PARAMETER (NUX=100)
      DIMENSION XU(NUX), YU(NUX), XUD(NUX), YUD(NUX), SU(NUX)
      LOGICAL LDONE
C
      LOGICAL LGUI
      CHARACTER*1 CHKEY
C
      DATA SH /0.010/
C
      CALL GETCOLOR(ICOL0)
      CALL GETPEN(IPEN0)
C
      KDONE  = 1
      KERASE = 2
      KABORT = 3
      KINSIDE = 4
C
      XDWIN = XWIN(2) - XWIN(1)
      YDWIN = YWIN(2) - YWIN(1)
C
      CALL NEWPEN(5)
C
      X1 = XWIN(1) + 0.71*XDWIN
      X2 = XWIN(1) + 0.79*XDWIN
      Y1 = YWIN(1) + 0.01*YDWIN
      Y2 = YWIN(1) + 0.05*YDWIN
      CALL GUIBOX(KABORT, X1,X2,Y1,Y2, 'RED'   , ' Abort ')
C
      X1 = XWIN(1) + 0.81*XDWIN
      X2 = XWIN(1) + 0.89*XDWIN
      Y1 = YWIN(1) + 0.01*YDWIN
      Y2 = YWIN(1) + 0.05*YDWIN
      CALL GUIBOX(KERASE, X1,X2,Y1,Y2, 'YELLOW', ' Erase ')
C
      X1 = XWIN(1) + 0.91*XDWIN
      X2 = XWIN(1) + 0.99*XDWIN
      Y1 = YWIN(1) + 0.01*YDWIN
      Y2 = YWIN(1) + 0.05*YDWIN
      CALL GUIBOX(KDONE , X1,X2,Y1,Y2, 'GREEN', ' Done ')
C
      X1 = XMOD(1)
      X2 = XMOD(2)
      Y1 = YMOD(1)
      Y2 = YMOD(2)
      CALL GUIBOX(KINSIDE, X1,X2,Y1,Y2, 'ORANGE' , ' ' )
C
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      CALL NEWPEN(IPEN0)
C
      WRITE(*,*)
      WRITE(*,*) 'Click on new values...'
      WRITE(*,*)
C
      NUBEG = 1
C
 5    CONTINUE
      CALL NEWCOLORNAME('MAGENTA')
      DO NU = NUBEG, NUX
C
C------ fetch x-y point coordinates from user
        CALL GETCURSORXY(XU(NU),YU(NU),CHKEY)
C
        IF    (LGUI(KABORT,XU(NU),YU(NU))
     &         .OR. INDEX('Aa',CHKEY).GT.0) THEN
C------- return with no changes
         GO TO 90
C
        ELSEIF(LGUI(KERASE,XU(NU),YU(NU))
     &         .OR. INDEX('Ee',CHKEY).GT.0) THEN
         IF(NU.LE.1) THEN
          WRITE(*,*) 'No more points to clear'
          NUBEG = 1
         ELSE
C-------- clear previous point, overplot it white to clear it from screen
          NUBEG = NU - 1
          CALL NEWCOLORNAME('WHITE')
          CALL PLSYMB(XU(NUBEG),YU(NUBEG),SH,3,0.0,0)
         ENDIF
C
C------- keep accepting points starting from NUBEG
         GO TO 5
C
        ELSEIF(LGUI(KDONE,XU(NU),YU(NU))
     &         .OR. INDEX('Dd',CHKEY).GT.0) THEN
C------- go process inputs
         GO TO 10
C
        ELSEIF(LGUI(KINSIDE,XU(NU),YU(NU))) THEN
C------- normal click inside modify-window: plot small cross at input point
         CALL PLSYMB(XU(NU),YU(NU),SH,3,0.0,0)
C
        ELSE
C------- must be somewhere outside
         GO TO 10
C
        ENDIF
      ENDDO
      WRITE(*,*) 'MODIXY: User-input array limit NUX reached'
C
C---- pick up here when finished with input
 10   CONTINUE
      IF(INDEX('Dd',CHKEY).GT.0) THEN
C----- last point was entered with a "D" ...  add it to list
       CALL PLSYMB(XU(NU),YU(NU),SH,3,0.0,0)
      ELSE
C----- discard last point
       NU = NU-1
      ENDIF
C
      CALL PLFLUSH
C
      IF(NU.LT.2) THEN
       WRITE(*,*)
       WRITE(*,*) 'Need at least 2 points'
       GO TO 90
      ENDIF
C
C---- set first- and last-specified point
      XUSP1 = XU(1)
      YUSP1 = YU(1)
C
      XUSP2 = XU(NU)
      YUSP2 = YU(NU)
C
C---- undo plot offsets and scales
      DO IU = 1, NU
        XU(IU) = XU(IU)/XSF + XOFF
        YU(IU) = YU(IU)/YSF + YOFF
      ENDDO
C
C---- remove doubled endpoints and tripled interior points
      DO IPASS = 1, 12345
        LDONE = .TRUE.
        IU = 2
        IF(XU(IU).EQ.XU(IU-1)) THEN
         LDONE = .FALSE.
         IUREM = IU
        ENDIF
        DO IU = 3, NU
          IF( XU(IU).EQ.XU(IU-1) .AND.
     &        XU(IU).EQ.XU(IU-2)      ) THEN
           LDONE = .FALSE.
           IUREM = IU
          ENDIF
        ENDDO
        IU = NU
        IF(XU(IU).EQ.XU(IU-1)) THEN
         LDONE = .FALSE.
         IUREM = IU
        ENDIF
C
        IF(LDONE) THEN
         GO TO 30
        ELSE
         DO IU = IUREM, NU-1
           XU(IU) = XU(IU+1)
           YU(IU) = YU(IU+1)
         ENDDO
         NU = NU - 1
        ENDIF
      ENDDO
C
C---- pick up here when no more points to be removed
 30   CONTINUE
      IF(NU.LT.2) THEN
       WRITE(*,*)
       WRITE(*,*) 'Need at least 2 points'
       GO TO 90
      ENDIF
C
C
C---- find which X,Y input point is closest to first-specified point
      ISMOD = 1
      IMOD1 = IFRST(ISMOD)
      XUI = (X(IMOD1)-XOFF)*XSF
      YUI = (Y(IMOD1)-YOFF)*YSF
      DSQMIN = (XUI-XUSP1)**2 + (YUI-YUSP1)**2
      DO IS = 1, NSIDE
        DO I = IFRST(IS), ILAST(IS)
C-------- convert input arrays to plot coordinates
          XUI = (X(I)-XOFF)*XSF
          YUI = (Y(I)-YOFF)*YSF
          DSQ = (XUI-XUSP1)**2 + (YUI-YUSP1)**2
C
          IF(DSQ .LT. DSQMIN) THEN
C---------- this point is the closest so far... note its indices
            DSQMIN = DSQ
            ISMOD = IS
            IMOD1 = I
          ENDIF
        ENDDO
      ENDDO
C
C---- set side and function to be modified
      IS = ISMOD
C
C
C---- find which X,Y input point is closest to last-specified point,
C-    but check only element IS
      IMOD2 = IFRST(IS)
      XUI = (X(IMOD2)-XOFF)*XSF
      YUI = (Y(IMOD2)-YOFF)*YSF
      DSQMIN = (XUI-XUSP2)**2 + (YUI-YUSP2)**2
        DO I = IFRST(IS), ILAST(IS)
C-------- convert input arrays to plot coordinates
          XUI = (X(I)-XOFF)*XSF
          YUI = (Y(I)-YOFF)*YSF
          DSQ = (XUI-XUSP2)**2 + (YUI-YUSP2)**2
C
          IF(DSQ .LT. DSQMIN) THEN
C---------- this point is the closest so far... note its indices
            DSQMIN = DSQ
            IMOD2 = I
          ENDIF
        ENDDO
C
      IF    (IMOD1.EQ.IMOD2) THEN
       WRITE(*,*)
       WRITE(*,*) 'Graft endpoints must be distinct'
       GO TO 90
      ELSEIF(IMOD1.GT.IMOD2) THEN
C----- reverse the input-point ordering to get increasing S values
       DO IU = 1, NU/2
         XTMP = XU(IU)
         YTMP = YU(IU)
         XU(IU) = XU(NU-IU+1)
         YU(IU) = YU(NU-IU+1)
         XU(NU-IU+1) = XTMP
         YU(NU-IU+1) = YTMP
       ENDDO
       ITMP = IMOD1
       IMOD1 = IMOD2
       IMOD2 = ITMP
      ENDIF
C
C
C---- reset X,Y and dX/dS,dY/dS at first and last points of modified interval
      IU = 1
      IF(LBLEND .OR. IMOD1.NE.IFRST(IS)) THEN
C----- reset 1st input point to match contour, except if non-blended endpoint
       XU(IU) = X(IMOD1)
       YU(IU) = Y(IMOD1)
      ENDIF
      IF(LBLEND .AND. IMOD1.NE.IFRST(IS)) THEN
C----- match derivatives to current contour, except at the endpoints
       XUD1 = XD(IMOD1)
       YUD1 = YD(IMOD1)
      ELSE
C----- do not constrain 1st derivatives (set zero 3rd derivative instead)
       XUD1 = -999.0
       YUD1 = -999.0
      ENDIF
C
      IU = NU
      IF(LBLEND .OR. IMOD2.NE.ILAST(IS)) THEN
C----- reset 1st input point to match contour, except if non-blended endpoint
       XU(IU) = X(IMOD2)
       YU(IU) = Y(IMOD2)
      ENDIF
      IF(LBLEND .AND. IMOD2.NE.IFRST(IS)) THEN
C----- match derivatives to current contour
       XUD2 = XD(IMOD2)
       YUD2 = YD(IMOD2)
      ELSE
C----- do not constrain 1st derivatives (set zero 3rd derivative instead)
       XUD2 = -999.0
       YUD2 = -999.0
      ENDIF
C
C---- set spline parameter
      CALL SCALC(XU,YU,SU,NU)
C
C---- shift and rescale spline parameter SU to match current S
      SU1 = SU(1)
      SU2 = SU(NU)
      DO IU = 1, NU
        SFRAC = (SU(IU)-SU1)/(SU2-SU1)
        SU(IU) = S(IMOD1)*(1.0-SFRAC) + S(IMOD2)*SFRAC
      ENDDO
C
C---- spline input function values
      CALL SPLIND(XU,XUD,SU,NU,XUD1,XUD2)
      CALL SPLIND(YU,YUD,SU,NU,YUD1,YUD2)
C
C
C---- go over all points on modified segment
      DO I = IMOD1, IMOD2
        SI = S(I)
        X(I) = SEVAL(SI,XU,XUD,SU,NU)
        Y(I) = SEVAL(SI,YU,YUD,SU,NU)
      ENDDO
C
      IF(LMODPL) THEN
C----- plot modified function over modified interval
       CALL NEWCOLORNAME('MAGENTA')
       IPEN = 3
       DO I = IMOD1, IMOD2
         XP = (X(I)-XOFF)*XSF
         YP = (Y(I)-YOFF)*YSF
         CALL PLOT(XP,YP,IPEN)
         IPEN = 2
       ENDDO
       CALL PLFLUSH
      ENDIF
C
C---- return normally
      CALL NEWCOLOR(ICOL0)
      RETURN
C
C-------------------------------------------------
 90   CONTINUE
      WRITE(*,*) 'No changes made'
      IMOD1 = IFRST(1)
      IMOD2 = IFRST(1) - 1
      ISMOD = 1
      CALL NEWCOLOR(ICOL0)
      RETURN
C
      END ! MODIXY
