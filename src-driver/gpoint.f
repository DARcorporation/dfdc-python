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

      SUBROUTINE PNTMOD(IEL,LPLTEL)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      LOGICAL LPLTEL(*)
C
      DIMENSION XBOX(2), YBOX(2)
      DIMENSION RINPUT(2)
      LOGICAL LGUI
      CHARACTER*1 ANS, KCHAR
C
 1000 FORMAT(A)
      CALL GETCOLOR(ICOL0)
C
 10   WRITE(*,1010)
 1010 FORMAT(/'     -------------------------------'
     &       /'      C orner set           Z oom'
     &       /'      A dd    points        U nzoom'
     &       /'      D elete points'
     &       /'      M ove   points'
     &       /'      N ew-element input')
      WRITE(*,1100)
 1100 FORMAT(/'     Option: ',$)
      READ(*,1000) ANS
C
C--------------------------------------------------------------------
      IF(ANS.EQ.' ') THEN
       RETURN
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Zz',ANS).NE.0) THEN
       IF(LPLOT) THEN
        CALL USETZOOM(.TRUE.,.TRUE.)
        CALL REPLOT(IDEV)
       ENDIF
       GO TO 10
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Uu',ANS).NE.0) THEN
       IF(LPLOT) THEN
        CALL CLRZOOM
        CALL REPLOT(IDEV)
       ENDIF
       GO TO 10
C
      ELSEIF(INDEX('CADMNcadmn',ANS).EQ.0) THEN
       WRITE(*,*) 'Unrecognized command'
       GO TO 10
C
      ENDIF
C
      IF(INDEX('Nn',ANS).NE.0) THEN
C------ prepare stuff for new element
C
        IF(NBEL.EQ.NEX) THEN
         WRITE(*,*) 'Element arrays will overflow.  No action taken.'
         GO TO 10
        ENDIF
C
C------ initialize new-point index accumulator
        IF(NBEL.EQ.0) THEN
         IBADD = 0
        ELSE
         IBADD = IBLAST(NBEL)
        ENDIF
C
        NBEL = NBEL+1
        IBFRST(NBEL) = IBADD + 1
        LPLTEL(NBEL) = .TRUE.
C
        CALL CLRHOM(NBEL)
        XBREFE(NBEL) = 0.
        YBREFE(NBEL) = 0.
        NETYPE(NBEL) = 0
C
C------ set initial element limits to current plot limits
        XBMINE(NBEL) = XMIN
        YBMINE(NBEL) = YMIN
        XBMAXE(NBEL) = XMAX
        YBMAXE(NBEL) = YMAX
C
        WRITE(*,1600) NBEL
 1600   FORMAT(/' Click on new points for element', I4, ' ...')
      ENDIF
C
C...........................................................
C.... start cursor input loop
C
 15   CONTINUE
C
C---- box in which modification cursor inputs are recognized
      CALL GETZOOMABS(XZOFF,YZOFF,XZFAC,YZFAC)
      XBOX(1) = MIN(0.0  ,       1.0001*XMARG)/XZFAC - XZOFF
      XBOX(2) = MIN(XWIND, XPAGE-1.0001*XMARG)/XZFAC - XZOFF
      YBOX(1) = MAX(0.0  ,       1.0001*YMARG)/YZFAC - YZOFF
      YBOX(2) = MAX(YWIND, YPAGE-1.0001*YMARG)/YZFAC - YZOFF
C
      KDONE = 1
      X1 = XBOX(1) + 0.91*(XBOX(2)-XBOX(1))
      X2 = XBOX(1) + 0.99*(XBOX(2)-XBOX(1))
      Y1 = YBOX(1) + 0.01*(YBOX(2)-YBOX(1))
      Y2 = YBOX(1) + 0.05*(YBOX(2)-YBOX(1))
      CALL GUIBOX(KDONE, X1,X2,Y1,Y2, 'GREEN', ' Done ')
C
c      KERASE = 2
c      X1 = XBOX(1) + 0.81*(XBOX(2)-XBOX(1))
c      X2 = XBOX(1) + 0.89*(XBOX(2)-XBOX(1))
c      Y1 = YBOX(1) + 0.01*(YBOX(2)-YBOX(1))
c      Y2 = YBOX(1) + 0.05*(YBOX(2)-YBOX(1))
c      CALL GUIBOX(KDONE, X1,X2,Y1,Y2, 'YELLOW', ' Erase ')
cC
c      KABORT = 3
c      X1 = XBOX(1) + 0.71*(XBOX(2)-XBOX(1))
c      X2 = XBOX(1) + 0.79*(XBOX(2)-XBOX(1))
c      Y1 = YBOX(1) + 0.01*(YBOX(2)-YBOX(1))
c      Y2 = YBOX(1) + 0.05*(YBOX(2)-YBOX(1))
c      CALL GUIBOX(KDONE, X1,X2,Y1,Y2, 'GREEN', ' Abort ')
C
      IF(INDEX('CANcan',ANS).NE.0) THEN
       IF(NB.EQ.IBX) THEN
        WRITE(*,*)
     &   'Buffer geometry arrays will overflow.  No action taken.'
        GO TO 10
       ENDIF
      ENDIF
C
C--------------------------------------------------------------------
      IF    (INDEX('Cc',ANS).NE.0) THEN
        CALL POINTF(XB,YB,NBEL,IBFRST,IBLAST,
     &            XBOFF,XBSF,YBOFF,YBSF,
     &            IELC,IBC,XBC,YBC,KDONE)
        IF(IBC.EQ.0) GO TO 10
C
        IF(IBC.EQ.IBFRST(IELC) .OR.
     &     IBC.EQ.IBLAST(IELC)) THEN
         WRITE(*,*) 'Cannot double end point.  No action taken.'
         GO TO 10
        ENDIF
C
C------ add doubled point
        DO IB = NB, IBC, -1
          XB(IB+1) = XB(IB)
          YB(IB+1) = YB(IB)
          SB(IB+1) = SB(IB)
          XBS(IB+1) = XBS(IB)
          YBS(IB+1) = YBS(IB)
        ENDDO
        NB = NB+1
C
        IBLAST(IELC) = IBLAST(IELC) + 1
        DO IEL = IELC+1, NBEL
          IBFRST(IEL) = IBFRST(IEL) + 1
          IBLAST(IEL) = IBLAST(IEL) + 1
        ENDDO
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Aa',ANS).NE.0) THEN
C------ determine interval  IBC-1...IBC  which is to contain added point
        KK = 20
        CALL POINTG(XB,XBS,YB,YBS,SB,NBEL,IBFRST,IBLAST,
     &            XBOFF,XBSF,YBOFF,YBSF,
     &            IELC,IBC,XBC,YBC,KDONE,KK)
        IF(IBC.EQ.0) THEN
         CALL ELPROC(IELC)
         GO TO 10
        ENDIF
C
C------ make room for new point
        NADD = 1
        CALL BSHIFT(IELC,IBC,NADD)
C
C------ set new point
        XB(IBC) = XBC
        YB(IBC) = YBC
        WRITE(*,2030) IBC, IELC, XBC, YBC
 2030   FORMAT(' New point',I4,', element',I3,':    x y  =', 2G12.4)
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Dd',ANS).NE.0) THEN
        CALL POINTF(XB,YB,NBEL,IBFRST,IBLAST,
     &            XBOFF,XBSF,YBOFF,YBSF,
     &            IELC,IBC,XBC,YBC,KDONE)
        IF(IBC.EQ.0) THEN
         CALL ELPROC(IELC)
         GO TO 10
        ENDIF
C
C------ remove point
        NADD = -1
        CALL BSHIFT(IELC,IBC,NADD)
C
        WRITE(*,2040) IBC, IELC, XBC, YBC
 2040   FORMAT(' Deleted point',I4,', element',I3,':    x y  =', 2G12.4)
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Mm',ANS).NE.0) THEN
        CALL POINTF(XB,YB,NBEL,IBFRST,IBLAST,
     &            XBOFF,XBSF,YBOFF,YBSF,
     &            IELC,IBC,XBC,YBC,KDONE)
        IF(IBC.EQ.0) THEN
         CALL ELPROC(IELC)
         GO TO 10
        ENDIF
C
        SHB = SHGT*XBSF
        CALL PLSYMB((XB(IBC)-XBOFF)*XBSF,
     &              (YB(IBC)-YBOFF)*YBSF,SHB,1,0.0,0)
C
        WRITE(*,2050) IBC, XB(IBC), YB(IBC)
 2050   FORMAT(' Move point',I4,' :  x =',F10.6,'   y =',F10.6,
     &         '   to cursor click ...')
C
        CALL GETCURSORXY(XCRS,YCRS,KCHAR)
        IF(LGUI(KDONE,XCRS,YCRS) .OR. INDEX('Dd',KCHAR).NE.0) GO TO 10
C
C------ go from screen to user's coordinates x,y
        XB(IBC) = XCRS/XBSF + XBOFF
        YB(IBC) = YCRS/YBSF + YBOFF
C
        IF(XB(IBC).EQ.XBC .AND. YB(IBC).EQ.YBC) THEN
          RINPUT(1) = XBC
          RINPUT(2) = YBC
          CALL ASKRN('Enter new x,y coordinates^',RINPUT,2)
          XB(IBC) = RINPUT(1)
          YB(IBC) = RINPUT(2)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Nn',ANS).NE.0) THEN
C
        CALL GETCURSORXY(XCRS,YCRS,KCHAR)
        IF(LGUI(KDONE,XCRS,YCRS) .OR. INDEX('Dd',KCHAR).NE.0) THEN
         CALL PANDEF(NBEL)
         CALL ELPROC(NBEL)
         GO TO 10
        ENDIF
C
        CALL NEWCOLOR(ICOLEL(NBEL))
        SHB = SHGT*XBSF
        CALL PLSYMB(XCRS,YCRS,SHB,1,0.0,0)
C
C------ go from screen to user's coordinates x,y
        IBADD = IBADD + 1
        XB(IBADD) = XCRS/XBSF + XBOFF
        YB(IBADD) = YCRS/YBSF + YBOFF
        IBLAST(NBEL) = IBADD
        NB = IBADD
C
        WRITE(*,2070) IBADD, XB(IBADD), YB(IBADD)
 2070   FORMAT(' i x y =', I4, 2X, 2G12.5)
C
        I = IBFRST(NBEL)
        N = IBLAST(NBEL) - IBFRST(NBEL) + 1
        CALL MINMAX(N,XB(I),XBMIN1,XBMAX1)
        CALL MINMAX(N,YB(I),YBMIN1,YBMAX1)
        XBMINE(NBEL) = MIN(XBMINE(NBEL),XBMIN1)
        XBMAXE(NBEL) = MAX(XBMAXE(NBEL),XBMAX1)
        YBMINE(NBEL) = MIN(YBMINE(NBEL),YBMIN1)
        YBMAXE(NBEL) = MAX(YBMAXE(NBEL),YBMAX1)
C
C--------------------------------------------------------------------
      ENDIF
C
C---- spline modified or new geometry
      IEL = NBEL
      I = IBFRST(IEL)
      N = IBLAST(IEL) - IBFRST(IEL) + 1
      CALL SCALC(XB(I),YB(I),SB(I),N)
      CALL SEGSPL(XB(I),XBS(I),SB(I),N)
      CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
C---- replot modified geometry
      CALL GPLINI(.FALSE.,NBEL,LPLTEL,
     &            XBMINE,XBMAXE,YBMINE,YBMAXE,
     &            XBOFF,XBSF,YBOFF,YBSF)
      CALL GBPLOT(LPLTEL)

c      I = IBFRST(IELC)
c      N = IBLAST(IELC) - IBFRST(IELC) + 1
c      NTGSPL = 20
cC
c      CALL NEWCOLORNAME('MAGENTA')
c      CALL NEWPEN(2)
c      CALL PLTGSP(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, 
c     &            XBOFF,XBSF,YBOFF,YBSF, NTGSPL)
cC
c      CALL NEWCOLOR(ICOL0)

      GO TO 15
C
      END ! PNTMOD



      SUBROUTINE BSHIFT(IELC,IBC,NADD)
C------------------------------------------------------------
C     Adds or removes NADD spaces in buffer geometry array,
C     starting at point IBC on element IELC.
C     Note: BSHIFT alters buffer coordinates XB,YB and spline data SB,XBS,YBS
C           the pointer IBFRST, IBLAST are updated only for IELC+1 to NBEL
C           updated # buffer points NB,
C           does not alter IBC (usually IBFRST(IELC))
C------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(NB+NADD .GT. IBX) THEN
       WRITE(*,*)
     &   'BSHIFT: Buffer arrays will overflow.  No action taken.'
       RETURN
      ENDIF
C
      IF(IELC.LT.1 .OR. IELC.GT.NBEL .OR. NADD.EQ.0) RETURN
C
      IF(NADD.GT.0) THEN
       IB1 = NB
       IB2 = IBC
       IBDEL = -1
      ELSE
       IB1 = IBC+1
       IB2 = NB
       IBDEL = +1
      ENDIF
C
      DO IB = IB1, IB2, IBDEL
        SB(IB+NADD) = SB(IB)
        XB(IB+NADD) = XB(IB)
        YB(IB+NADD) = YB(IB)
        XBS(IB+NADD) = XBS(IB)
        YBS(IB+NADD) = YBS(IB)
      ENDDO
      NB = NB + NADD
C
      IBLAST(IELC) = IBLAST(IELC) + NADD
      DO IEL = NBEL, IELC+1, -1
        IBFRST(IEL) = IBFRST(IEL) + NADD
        IBLAST(IEL) = IBLAST(IEL) + NADD
      ENDDO
C
      RETURN
      END



      SUBROUTINE POINTF(X,Y,NEL,IFRST,ILAST,
     &                  XOFF,XSF,YOFF,YSF,
     &                  IELC,IC,XC,YC,KDONE)
      DIMENSION X(*),Y(*)
      DIMENSION IFRST(*), ILAST(*)
C--------------------------------------------------------
C     Finds the node IC nearest to cursor location XC,YC.
C     Returns IC=0 if XC,YC falls inside GUI(KABORT) box.
C--------------------------------------------------------
      LOGICAL LGUI
      CHARACTER*1 KCHAR
C
CCC      XMOD(XTMP) = XSF * (XTMP - XOFF)
CCC      YMOD(YTMP) = YSF * (YTMP - YOFF)
C
      WRITE(*,*)
      WRITE(*,*) 'Specify point with cursor'
C
C---- read geometry point coordinates
      CALL GETCURSORXY(XCRS,YCRS,KCHAR)
C
      IF(LGUI(KDONE,XCRS,YCRS) .OR. INDEX('Dd',KCHAR).NE.0) THEN
       IC = 0
       RETURN
      ENDIF
C
C---- go from screen to internal coordinates X,Y
      XC = XCRS/XSF + XOFF
      YC = YCRS/YSF + YOFF
C
C---- find closest node
      IEL = 1
      I = IFRST(IEL)
C
      IELC = IEL
      IC = I
      DSQMIN = (X(IC) - XC)**2 + (Y(IC) - YC)**2
C
      DO IEL = 1, NEL
        DO 10 I = IFRST(IEL), ILAST(IEL)
          DSQ = (X(I) - XC)**2 + (Y(I) - YC)**2
          IF(DSQ .LT. DSQMIN) THEN
            DSQMIN = DSQ
            IELC = IEL
            IC = I
          ENDIF
 10     CONTINUE
      ENDDO
C
      RETURN
      END ! POINTF



      SUBROUTINE POINTG(X,XS,Y,YS,S,NEL,IFRST,ILAST,
     &                  XOFF,XSF,YOFF,YSF,
     &                  IELC,IC,XC,YC,KDONE,KK)
      DIMENSION X(*),XS(*),Y(*),YS(*),S(*)
      DIMENSION IFRST(*), ILAST(*)
C--------------------------------------------------------
C     Finds the interval IC-1..IC with spline nearest
C     to cursor location XC,YC.  The spline is sampled
C     at KK sub-intervals. 
C     Returns IC=0 if XC,YC falls inside GUI(KABORT) box.
C--------------------------------------------------------
      LOGICAL LGUI
      CHARACTER*1 KCHAR
CCC      XMOD(XTMP) = XSF * (XTMP - XOFF)
CCC      YMOD(YTMP) = YSF * (YTMP - YOFF)
C
      WRITE(*,*)
      WRITE(*,*) 'Specify point with cursor'
C
C---- read geometry point coordinates
      CALL GETCURSORXY(XCRS,YCRS,KCHAR)
C
      IF(LGUI(KDONE,XCRS,YCRS) .OR. INDEX('Dd',KCHAR).NE.0) THEN
       IC = 0
       RETURN
      ENDIF
C
C---- go from screen to internal coordinates X,Y
      XC = XCRS/XSF + XOFF
      YC = YCRS/YSF + YOFF
C
C---- find closest spline node
      IEL = 1
      I = IFRST(IEL)
C
      IELC = IEL
      IC = I + 1
      KC = 0
      DSQMIN = (X(I) - XC)**2 + (Y(I) - YC)**2
C
      DO IEL = 1, NEL
        DO 10 I = IFRST(IEL)+1, ILAST(IEL)
          IM = I-1
          DS = S(I) - S(IM)
C
C-------- skip zero-width spline interval
          IF(DS .EQ. 0.0) GO TO 10
C
C-------- search sub-interval points
          DO K = 1, KK
            SK = S(IM) + DS*FLOAT(K)/FLOAT(KK)
            XK = SEVAL(SK,X(IM),XS(IM),S(IM),2)
            YK = SEVAL(SK,Y(IM),YS(IM),S(IM),2)
            DSQ = (XK - XC)**2 + (YK - YC)**2
            IF(DSQ .LT. DSQMIN) THEN
              DSQMIN = DSQ
              IELC = IEL
              IC = I
              KC = K
            ENDIF
          ENDDO
 10     CONTINUE
      ENDDO
C
      IF(KC.EQ.KK .AND. IC.LT.ILAST(IELC)) THEN
C------ spline node is the nearest point -- see on which side we are
        DOTP = (X(IC)-XC)*XS(IC) + (Y(IC)-YC)*YS(IC)
        IF(DOTP .LT. 0.0) IC = IC + 1
      ENDIF
C
      RETURN
      END ! POINTG

