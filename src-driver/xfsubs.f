
C***********************************************************************
C    XFOIL subroutines for ESLOFT
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************
C
C    This source file consists of subroutines borrowed from XFOIL 6.97.
C    Some have been modified to work with ESLOFT. Others are unchanged.
C    Philip Carter, Esotec Developments, 4 April 2009
C    philip (at) esotec (dot) org
C
C***********************************************************************



      SUBROUTINE XPLTAIR(XX,XXP,YY,YYP,SS,NN,XOFF,XSF,YOFF,YSF,ICOLOR)
C------------------------------------------------------------------
C     Plots passed-in airfoil
C     Borrowed from xfoil sources (PLTAIR in xfoil)
C------------------------------------------------------------------
      DIMENSION XX(NN), XXP(NN), YY(NN), YYP(NN), SS(NN)
C
      XMOD(XTMP) = XSF * (XTMP - XOFF)
      YMOD(YTMP) = YSF * (YTMP - YOFF)
C
      NT = 20
ccc      NT = 50
C
      CALL GETCOLOR(ICOL0)
      CALL NEWCOLOR(ICOLOR)
C
      DO 60 I=2, NN
        DS = SS(I) - SS(I-1)
        CALL PLOT(XMOD(XX(I-1)),YMOD(YY(I-1)),3)
C
C------ subdivide current panel into NT segments for smoother airfoil plot
        DO 610 IT=1, NT
          ST = SS(I-1) + DS*FLOAT(IT)/FLOAT(NT)
          XT = SEVAL(ST,XX,XXP,SS,NN)
          YT = SEVAL(ST,YY,YYP,SS,NN)
          CALL PLOT(XMOD(XT),YMOD(YT),2)
  610   CONTINUE
   60 CONTINUE
C
      CALL NEWCOLOR(ICOL0)
C
      CALL PLFLUSH
C
      RETURN
      END  ! XPLTAIR
C
C--------------------------------------------------------------------




      SUBROUTINE PLTINIX
C-------------------------------------------------------------------
C    ** Modified from xfoil PLTINI**
C-------------------------------------------------------------------
      INCLUDE 'PLOT.INC'
C
C---- terminate old plot if any
      IF(LPLOT) CALL PLEND
C
C---- initialize new plot
      IF(LLAND) THEN
        SIGNFR =  SCRNFR
      ELSE
        SIGNFR = -SCRNFR*PORTFAC
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
      END ! PLTINIX
C
C----------------------------------------------------------------------





      SUBROUTINE GOFINIX(PLOTARX)
C----------------------------------------------------------
C     Sets initial airfoil scaling and offset parameters 
C     ** Modified from xfoil sources **  
C----------------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      CS = CSIZE
C
C---- get airfoil bounding box
C---- special treatment for trans plots
C
      IF(IXTYPE.EQ.3) THEN
         IF(OUTFAC.EQ.MTI) THEN
            TRANFAC = MTI
         ELSE
            TRANFAC = 100.0
         ENDIF
C
         XBMINX = TRXX(1,1)
         YBMINX = TRYY(1,1)
         XBMAXX = TRXX(1,1)
         YBMAXX = TRYY(1,1)
         DO I=1,NXPLOT
            IAF=IXPLOT(I)
            DO J=1,NPP
               XTEMP = TRXX(J,IAF)*TRANFAC
               YTEMP = TRYY(J,IAF)*TRANFAC
               XBMINX= MIN(XBMINX,XTEMP)
               XBMAXX= MAX(XBMAXX,XTEMP)
               YBMINX= MIN(YBMINX,YTEMP)
               YBMAXX= MAX(YBMAXX,YTEMP)
            ENDDO
         ENDDO
C
         XBMINX= XBMINX * 1.25
         XBMAXX= XBMAXX * 1.25
C
C---- normalized airfoils
C
      ELSE
         XBMINX = XBB(1)
         YBMINX = YBB(1)
         XBMAXX = XBB(1)
         YBMAXX = YBB(1)
         DO I=1, NBB
           XBMINX = MIN(XBMINX,XBB(I))
           YBMINX = MIN(YBMINX,YBB(I))
           XBMAXX = MAX(XBMAXX,XBB(I))
           YBMAXX = MAX(YBMAXX,YBB(I))
         ENDDO
         YBMINX = YBMINX * 1.02
         YBMAXX = YBMAXX * 1.02
      ENDIF
C
C
C---- set camber and thickness distributions
C      CALL GETCAM(XCM,YCM,NCM,XTK,YTK,NTK,
C     &             XBB,XBBP,YBB,YBBP,SBB,NBB )
C
      XRANGE = XBMAXX - XBMINX
      YRANGE = YBMAXX - YBMINX
C
C---- set x,y scaling factors needed for O(1) size plot with "nice" limits
      CALL SCALIT(1,0.95*XRANGE,0.0,XSFX)
      CALL SCALIT(1,0.95*YRANGE,0.0,YSFX)
C  
C---- grid increment as a fraction of a nice upper bound on delta x
cc      DXYGX = 0.1 / XSFX
C
      DXYGX = 0.1 / MIN(XSFX,YSFX)
C  
C---- set "nice" grid limits as integer multiples of DXYGX
c      XGMAX = DXYGX*(INT(XBMAX/DXYGX+1000.05) - 999)
c      XGMIN = DXYGX*(INT(XBMIN/DXYGX-1000.05) + 999)
c      YGMAX = DXYGX*(INT(YBMAX/DXYGX+1000.25) - 999)
c      YGMIN = DXYGX*(INT(YBMIN/DXYGX-1000.25) + 999)
C  
C---- set "nice" grid limits as integer multiples of DXYGX
C
      XGMAXX = DXYGX*(INT(XBMAXX/DXYGX+1001.01) - 1000)
      XGMINX = DXYGX*(INT(XBMINX/DXYGX-1001.01) + 1000)
      YGMAXX = DXYGX*(INT(YBMAXX/DXYGX+1001.01) - 1000)
      YGMINX = DXYGX*(INT(YBMINX/DXYGX-1001.01) + 1000)
C
C---- set minimum scaling factor to fit airfoil or grid
C
      IF(LGGRIDX) THEN
        XRANGE = XGMAXX - XGMINX
        YRANGE = YGMAXX - YGMINX
      ELSE
        XRANGE = XBMAXX - XBMINX
        YRANGE = YBMAXX - YBMINX
      ENDIF
C
      RANGE = MAX(XRANGE,YRANGE)
C
      SF = MIN( 1.0/XRANGE , PLOTARX/YRANGE )
      XSFX = SF
      YSFX = SF
      CHGX = 0.75*CS * RANGE*SF
C--- HHY 4/24/01 keep the character size from getting too low

      CHGX = MAX(CHGX,0.0075)
C
      IF(LGGRIDX) THEN
C------ set offsets to position grid, with space for numerical axis annotations
        XOFFX = XGMINX - 0.05*RANGE - 3.0*CHGX/SF
        YOFFX = YGMINX - 0.05*RANGE - 2.0*CHGX/SF
      ELSE
C------ set offsets to position airfoil
        XOFFX = XBMINX - 0.05*RANGE
        YOFFX = YBMINX - 0.05*RANGE
      ENDIF
C
      RETURN
      END   ! GOFINIX
C
C------------------------------------------------------------------------




      SUBROUTINE PLOTGX(NDSK,ICOLOR)
C--------------------------------------------------------------
C     Plots buffer airfoil with ticked chord line or grid
C--------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      CHARACTER*80 NAMXFP,XHEAD
      CHARACTER*2 TRANLEG
C
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
      XMOD(XTMP) = XSFX * (XTMP - XOFFX)
      YMOD(YTMP) = YSFX * (YTMP - YOFFX)
C
C---- node tick mark size and corner symbol size
C
      DTICK = GTICKX*(SBB(NBB)-SBB(1))
      SSH   = DTICK * 3.0
      CS    = CSIZE * 0.7
      CSL   = CS    * 1.2
      CSD   = CS    * 1.0
C
      CALL NCALC(XBB,YBB,SBB,NBB,WX1,WX2)
C
      CALL GETCOLOR(ICOL0)
      CALL NEWCOLOR(ICOLOR)
C
      IF(LGGRIDX) THEN
        CALL GRDAIR(XGMINX,XGMAXX,YGMINX,YGMAXX,DXYGX,DXYGX,CHGX,
     &              .TRUE.,.TRUE.,XOFFX,XSFX,YOFFX,YSFX, LMASK2)
        XL0 = XMOD(XGMINX)
        YL0 = YMOD(YGMAXX) + 2.0*CS
      ELSE
C------ plot chord line and tick marks every 10% chord
C        CALL NEWPEN(1)
C        CALL PLOT(XMOD(0.0),YMOD(0.0),3)
C        CALL PLOT(XMOD(1.0),YMOD(0.0),2)
C        DO ITICK=1, 10
C          XPLT = FLOAT(ITICK)/10.0
C          CALL PLOT(XMOD(XPLT),YMOD(0.003),3)
C          CALL PLOT(XMOD(XPLT),YMOD(-.003),2)
C        ENDDO
C
        XL0 = XMOD(XBMINX)
        YL0 = YMOD(YBMAXX) + 2.0*CS
      ENDIF
C      IF(LPLCAM)  YL0 = YSFX*(YCMAX-DYOFFC-YOFFX) + 2.0*CS
C
      CALL PLFLUSH
C
      CALL NEWPEN(2)
      CALL XPLTAIR(XBB,XBBP,YBB,YBBP,SBB,NBB,
     &             XOFFX,XSFX,YOFFX,YSFX,ICOLOR)
C
      IF(LGTICKX) THEN
C----- draw tiny tick mark normal to airfoil surface at each panel node
       DO I=2, NBB-1
         CALL PLOT(XMOD(XBB(I)             ),
     &             YMOD(YBB(I)             ),3)
         CALL PLOT(XMOD(XBB(I)-DTICK*WX1(I)),
     &             YMOD(YBB(I)-DTICK*WX2(I)),2)
       ENDDO
      ENDIF
c
cC---- plot symbol at nose
c      CALL NSFIND(STLE,XB,XBP,YB,YBP,SB,NB)
c      XT = SEVAL(STLE,XB,XBP,SB,NB)
c      YT = SEVAL(STLE,YB,YBP,SB,NB)
c      CALL PLSYMB(XMOD(XT),YMOD(YT),0.005*XSFX,5,0.0,0)
c
C---- put symbol at any doubled point
C      DO I=1, NBB-1
C        IF(SBB(I) .EQ. SBB(I+1))
C     &     CALL PLSYMB(XMOD(XBB(I)),YMOD(YBB(I)),SSH,5,0.0,0)
C      ENDDO
C
C      IF(LPLCAM) THEN
C       CALL PLTCAM(' ')
C      ENDIF
C
      IF(LGPARMX) THEN        
        CALL NEWPEN(2)
        CALL GPARPLX(XL0,YL0,CS,.TRUE.,NAMEX,
     &              CHORDBX,AREABX,RADBLEX,ANGBTEX,
     &              EI11BA,EI22BA,APX1BA,APX2BA,
     &              EI11BT,EI22BT,APX1BT,APX2BT,
     &              THICKBX,CAMBRBX,TETHX)
      ENDIF
C
C---- plot header
C
      IAF=1
      CALL GETLNAMES(LONAME,NROTOR,NDSK,IAF,1,NAMXFP)
C
      CALL STRIP(SVERSION,NEV)
      NEVP=NEV+7
C
      YTIT = 0.93*YWIND
      XTIT = XMOD(XGMINX)
      TITAD= 0.5*CSL
      XVLEG= XMOD(XGMAXX)-CSD*FLOAT(NEVP)
      YVLEG= YMOD(YGMAXX)-CSD*2.1
C
      CALL PLOTABS(XTIT,YTIT,3)
      CALL NEWPEN(2)
      CALL NEWCOLORNAME('green')
C
C      CALL PLOT(XWIND*0.85,0.0,2)
C
      IF(IXTYPE.EQ.1)THEN
         XHEAD= 'Parent Airfoils: '
         CALL PLCHAR(XTIT,0.8,CSL,XHEAD,0.0,17)
      ELSEIF(IXTYPE.EQ.2)THEN
         XHEAD= 'Blended Sections: '
         CALL PLCHAR(XTIT,0.8,CSL,XHEAD,0.0,18)
      ELSEIF(IXTYPE.EQ.3) THEN
         XHEAD= 'Transformed Sections: '
         CALL PLCHAR(XTIT,0.8,CSL,XHEAD,0.0,22)
      ENDIF
C
      CALL PLCHAR(999.,0.8,CSL,NAMXFP,0.0,-1)
C
      CALL NEWCOLORNAME('black')
      YCM= YMOD(YGMAXX- DXYGX* 1.5)
      XCM= XMOD(XGMINX)-CS* 3.0
C
      IF(IXTYPE.EQ.3) THEN
         IF(OUTFAC.EQ.1.0.OR.OUTFAC.EQ.100.0.OR.
     &      OUTFAC.EQ.1000.0) THEN
            TRANLEG = 'cm'
         ELSE
            TRANLEG = 'in'
         ENDIF
         CALL PLCHAR(XCM,YCM,CSL,TRANLEG,0.0,2)
      ENDIF
C
      CALL NEWCOLORNAME('cyan')
      CALL PLCHAR(XVLEG,YVLEG,CSD,'DFDC v',0.0,6)
      CALL PLCHAR(999. ,YVLEG,CSD,SVERSION,0.0,-1)
C
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
C
      LGEOPLX = .TRUE.
      NOVER = 0
C
      RETURN
      END   ! PLOTGX
C
C---------------------------------------------------------------------



      SUBROUTINE OVERX(ICOLOR)
C----------------------------------------------------
C     Overlays plot of airfoil from coordinate file.
C----------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      XMOD(XTMP) = XSFX * (XTMP - XOFFX)
      YMOD(YTMP) = YSFX * (YTMP - YOFFX)
C
      DTICK = GTICKX*(SBB(NBB)-SBB(1))
      SSH   = DTICK * 3.0
      CS    = CSIZE
C
      CALL NCALC(XBB,YBB,SBB,NBB,WX1,WX2)
C
C      IF(.NOT.LPLOT) THEN
C       CALL PLTINIX
ccc       CALL PLOT(0.05,0.30,-3)
C      ENDIF
C
      CALL NEWPEN(2)
      CALL XPLTAIR(XBB,XBBP,YBB,YBBP,SBB,NBB,
     &             XOFFX,XSFX,YOFFX,YSFX,ICOLOR)
C
      CALL GETCOLOR(ICOL0)
      CALL NEWCOLOR(ICOLOR)
C
      IF(LGTICKX) THEN
       DO I=2, NBB-1
         CALL PLOT(XMOD(XBB(I)             ),
     &             YMOD(YBB(I)             ),3)
         CALL PLOT(XMOD(XBB(I)-DTICK*WX1(I)),
     &             YMOD(YBB(I)-DTICK*WX2(I)),2)
       ENDDO
      ENDIF
C
      CALL NEWCOLOR(ICOL0)
      CALL PLNEWPX(ICOLOR)
C
      RETURN
      END  ! OVERX
C
C----------------------------------------------------------------------



      SUBROUTINE PLNEWPX(ICOLOR)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      LOGICAL LLABEL
C
      XMOD(XTMP) = XSFX * (XTMP - XOFFX)
      YMOD(YTMP) = YSFX * (YTMP - YOFFX)
C
      CH = CSIZE * 0.7
C
      CALL GETCOLOR(ICOL0)
      CALL NEWCOLOR(ICOLOR)
      CALL NEWPEN(2)
C
      NOVER = NOVER + 1
C
      IF(NOVER.GT.21) THEN
         IPARLEV=2
         NOVERX = NOVER-22
      ELSEIF(NOVER.GT.10) THEN
         IPARLEV=1
         NOVERX = NOVER-11
      ELSE
         IPARLEV= 0
         NOVERX = NOVER
      ENDIF
C
      IF(NOVERX.EQ.0) THEN
         LLABEL= .TRUE.
      ELSE
         LLABEL= .FALSE.
      ENDIF
C
      IF(LGGRIDX) THEN
       XL0 = XMOD(XGMINX) +  9.0 *CH*FLOAT(NOVERX)
       YL0 = YMOD(YGMAXX) +  2.0 *CH + 15.0*CH*IPARLEV
      ELSE
       XL0 = XMOD(XBMINX) +  9.0 *CH*FLOAT(NOVERX)
       YL0 = YMOD(YBMAXX) +  2.0 *CH + 15.0*CH*IPARLEV
      ENDIF
      
C      IF(LPLCAM)  YL0 = YSFX*(YCMAX-YOFFX-DYOFFC) + 2.0*CH
C
      IF(LGPARMX) THEN
        CALL GPARPLX(XL0,YL0,CH,LLABEL,NAMEX,
     &              CHORDBX,AREABX,RADBLEX,ANGBTEX,
     &              EI11BA,EI22BA,APX1BA,APX2BA,
     &              EI11BT,EI22BT,APX1BT,APX2BT,
     &              THICKBX,CAMBRBX,TETHX)
      ENDIF
C
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
C
      RETURN
      END   ! PLNEWPX
C
C---------------------------------------------------------------



      SUBROUTINE GOFSETX
C----------------------------------------------------------
C     Sets grid-overlay parameters
C----------------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
C---- airfoil extent
      XBMINX = XBB(1)
      YBMINX = YBB(1)
      XBMAXX = XBB(1)
      YBMAXX = YBB(1)
      DO I=1, NBB
        XBMINX = MIN(XBMINX,XBB(I))
        YBMINX = MIN(YBMINX,YBB(I))
        XBMAXX = MAX(XBMAXX,XBB(I))
        YBMAXX = MAX(YBMAXX,YBB(I))
      ENDDO
C
      RANGE = MAX( (XWIND/SIZE)/XSFX , (YWIND/SIZE)/YSFX )
C
C---- set bounding-box corner locations in user coordinates
      XG1 = XOFFX + 0.1*RANGE + 4.0*CHGX/XSFX
      YG1 = YOFFX + 0.1*RANGE + 2.0*CHGX/YSFX
      XG2 = XOFFX - 0.1*RANGE + (XWIND/SIZE)/XSFX
      YG2 = YOFFX - 0.1*RANGE + (YWIND/SIZE)/YSFX
C
C---- crunch down onto airfoil limits
      XG1 = MAX(XG1,XBMINX)
      XG2 = MIN(XG2,XBMAXX)
      YG1 = MAX(YG1,YBMINX)
      YG2 = MIN(YG2,YBMAXX)
C
C---- set x,y scaling factors needed for O(1) size plot with "nice" limits
      CALL SCALIT(1,0.95*(XG2-XG1),0.0,GXSF)
      CALL SCALIT(1,0.95*(YG2-YG1),0.0,GYSF)
C  
      GSF = GXSF
ccc   GSF = MIN(GXSF,GYSF)
C
C---- grid increment as a fraction of a nice upper bound on delta x
      DXYGX = 0.1 / GSF
C  
C---- set "nice" grid limits as integer multiples of DXYGX
      XGMAXX = DXYGX*(INT(XG2/DXYGX+1001.01) - 1000)
      XGMINX = DXYGX*(INT(XG1/DXYGX-1001.01) + 1000)
      YGMAXX = DXYGX*(INT(YG2/DXYGX+1001.01) - 1000)
      YGMINX = DXYGX*(INT(YG1/DXYGX-1001.01) + 1000)
C
      RETURN
      END   ! GOFSETX
C
C--------------------------------------------------------------------------




      SUBROUTINE OFFGETX
C---------------------------------------------------
C     Sets new blowup parameters from cursor input.
C     ** Modified from xfoil sources **
C---------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      LOGICAL LSAMEX,LCURSX
      CHARACTER*1 KCHAR

C---- crosshair "+" symbol size
      DATA SH / 2.0 /
C
C---- these settings for future reference...
      LCURSX= .TRUE.
      LSAMEX= .TRUE.
C
C---- get current color
      CALL GETCOLOR(ICOL0)
C
C---- set new crosshair color
      CALL NEWCOLORNAME('red')
C
      XWS = XWIND/SIZE
      YWS = YWIND/SIZE
C
      IF(LCURSX) THEN
C
       WRITE(*,*)
       WRITE(*,*) 'Mark off corners of blowup area'
       WRITE(*,*) '(2 identical points default to current area)'
C
       CALL GETCURSORXY(XX1,YY1,KCHAR)
       CALL PLSYMB(XX1,YY1,SH,3,0.0,0)
       CALL PLFLUSH
       WRITE(*,*) 'x,y =', XX1/XSFX+XOFFX, YY1/YSFX+YOFFX
C
       CALL GETCURSORXY(XX2,YY2,KCHAR)
       CALL PLSYMB(XX2,YY2,SH,3,0.0,0)
       CALL PLFLUSH
       WRITE(*,*) 'x,y =', XX2/XSFX+XOFFX, YY2/YSFX+YOFFX
C
      ELSE
C
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
C
      CALL NEWCOLOR(ICOL0)
C
      IF(XX1.EQ.XX2 .AND. YY1.EQ.YY2) RETURN
C
C
      XCEN = 0.5*(XX1+XX2)/XSFX + XOFFX
      YCEN = 0.5*(YY1+YY2)/YSFX + YOFFX
      XDIF = ABS(XX2 - XX1)/XSFX
      YDIF = ABS(YY2 - YY1)/YSFX
C
      IF(XDIF.EQ.0.0) XDIF = 1.0E-5
      IF(YDIF.EQ.0.0) YDIF = 1.0E-5
C
      XOFFX = MIN(XX1,XX2)/XSFX + XOFFX
      YOFFX = MIN(YY1,YY2)/YSFX + YOFFX
      XSFX  = XWS/XDIF
      YSFX  = YWS/YDIF
C
      IF(LSAMEX) THEN
C------ set equal x,y scales
        SF = MIN( XSFX , YSFX )
        XSFX = SF
        YSFX = SF
C
C------ re-center the blowup
        XOFFX = XCEN - 0.5*XDIF
        YOFFX = YCEN - 0.5*YDIF
      ENDIF
C
      RETURN
      END   ! OFFGETX
C
C
C=====================================================================
C
C---- XFOIL standalone subroutines from here on, thank goodness...



      SUBROUTINE GEOPARX(X,XP,Y,YP,S,N, T,
     &             SLE,CHORD,AREA,RADLE,ANGTE,
     &             EI11A,EI22A,APX1A,APX2A,
     &             EI11T,EI22T,APX1T,APX2T,
     &             THICK,CAMBR,XTHICK,XCAMBR)
      DIMENSION X(*), XP(*), Y(*), YP(*), S(*), T(*)
C
      PARAMETER (IBX=600)
      DIMENSION
     &     XCAM(2*IBX), YCAM(2*IBX), YCAMP(2*IBX),
     &     XTHK(2*IBX), YTHK(2*IBX), YTHKP(2*IBX)
C------------------------------------------------------
C     Sets geometric parameters for airfoil shape
C------------------------------------------------------
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
C
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
C
      CHSQ = (XTE-XLE)**2 + (YTE-YLE)**2
      CHORD = SQRT(CHSQ)
C
      CURVLE = CURV(SLE,X,XP,Y,YP,S,N)
C
      RADLE = 0.0
      IF(ABS(CURVLE) .GT. 0.001*(S(N)-S(1))) RADLE = 1.0 / CURVLE
C
      ANG1 = ATAN2( -YP(1) , -XP(1) )
      ANG2 = ATANC(  YP(N) ,  XP(N) , ANG1 )
      ANGTE = ANG2 - ANG1
C
      DO I=1, N
        T(I) = 1.0
      ENDDO
C
      CALL AECALCX(N,X,Y,T, 1, 
     &            AREA,XCENA,YCENA,EI11A,EI22A,APX1A,APX2A)
C
      CALL AECALCX(N,X,Y,T, 2, 
     &            SLEN,XCENT,YCENT,EI11T,EI22T,APX1T,APX2T)
C
C--- Old, approximate thickness,camber routine (on discrete points only)
C      CALL TCCALC(X,XP,Y,YP,S,N, THICK,XTHICK, CAMBR,XCAMBR )
C
C--- More accurate thickness and camber estimates
cc      CALL GETCAM(XCAM,YCAM,NCAM,XTHK,YTHK,NTHK,
cc     &            X,XP,Y,YP,S,N )
cc      CALL GETMAX(XCAM,YCAM,YCAMP,NCAM,XCAMBR,CAMBR)
cc      CALL GETMAX(XTHK,YTHK,YTHKP,NTHK,XTHICK,THICK)
cc      THICK = 2.0*THICK
C
C      WRITE(*,1000) THICK,XTHICK,CAMBR,XCAMBR
C
      CALL TCSING(X,XP,Y,YP,S,N,
     &            THICK,XTHICK,CAMBR,XCAMBR)
C
C
      RETURN
      END ! GEOPARX
C
C----------------------------------------------------------------------



      SUBROUTINE GPARPLX(X0,Y0,CH, LABEL, NAME,
     &                   CHORD,AREA,RADLE,ANGTE,
     &                   EI11A,EI22A,APX1A,APX2A,
     &                   EI11T,EI22T,APX1T,APX2T,
     &                   THICK,CAMBR,TET)
C------------------------------------------------------------------
C     Modified from xfoil sources
C-----------------------------------------------------------------
      LOGICAL LABEL
      EXTERNAL PLCHAR
      CHARACTER NAME*(*)
C
      RTD = 45.0/ATAN(1.0)
C
C      XSPACE = 30.0*CH
      YSPACE =  2.0*CH
C
      X = X0
      Y = Y0
C
      IF(LABEL) THEN
       CALL PLCHAR(X,Y,CH,'    = ',0.0, 6)
       CALL PLMATH(X,Y,CH,'Oq    ',0.0, 7)
       CALL PLSUBS(X+CH,Y,CH,'TE',0.0, 2, PLCHAR)
      ENDIF
      CALL PLNUMB(X+7.0*CH,Y,CH,ANGTE*RTD        ,0.0, 2)
      CALL PLMATH(999.,Y,CH,'"'              ,0.0, 1)
      Y = Y + YSPACE
C
      IF(LABEL) THEN
       CALL PLCHAR(X,Y,CH,'t   = ',0.0, 6)
       CALL PLSUBS(X,Y,CH,'TE',0.0, 2, PLCHAR)
      ENDIF
      CALL PLNUMB(X+7.0*CH,Y,CH,TET,0.0, 5)
      Y = Y + YSPACE
C
      IF(LABEL) THEN
       CALL PLCHAR(X,Y,CH,'r   = ',0.0, 6)
       CALL PLSUBS(X,Y,CH,'LE',0.0, 2, PLCHAR)
      ENDIF
      CALL PLNUMB(X+7.0*CH,Y,CH,RADLE,0.0, 5)
      Y = Y + YSPACE
C
      IF(LABEL) THEN
       CALL PLCHAR(X,Y,CH,'camb= ',0.0, 6)
      ENDIF
      CALL PLNUMB(X+7.0*CH,Y,CH,CAMBR,0.0, 5)
      Y = Y + YSPACE
C
      IF(LABEL) THEN
       CALL PLCHAR(X,Y,CH,'area= ',0.0, 6)
      ENDIF
      CALL PLNUMB(X+7.0*CH,Y,CH, AREA,0.0, 5)
      Y = Y + YSPACE
C
      IF(LABEL) THEN
       CALL PLCHAR(X,Y,CH,'t/c = ',0.0, 6)
      ENDIF
      CALL PLNUMB(X+7.0*CH,Y,CH,THICK,0.0, 5)
      Y = Y + YSPACE
C
C
c      X = X0  +  XSPACE
c      Y = Y0
cC
c      Y = Y + YSPACE
cC
c      CALL PLMATH(X,Y,1.4*CH,'I',0.0,1)
c      CALL PLMATH(X,Y,CH,'       2     ',0.0,-1)
c      CALL PLCHAR(X,Y,CH,' (y-y ) ds = ',0.0,-1)
c      CALL PLNUMB(999.,Y,CH, 1000.0*EI11T,0.0,4)
c      CALL PLMATH(999.,Y,CH,'#'   ,0.0,1)
c      CALL PLCHAR(999.,Y,CH, '10' ,0.0,2)
c      CALL PLMATH(999.,Y,CH,   '3',0.0,1)
c      CALL PLSUBS(X+4.0*CH,Y,CH,'o',0.0,1,PLCHAR)
c      Y = Y + YSPACE
cC
c      CALL PLMATH(X,Y,1.4*CH,'I',0.0,1)
c      CALL PLMATH(X,Y,CH,'       2     ',0.0,-1)
c      CALL PLCHAR(X,Y,CH,' (y-y ) dA = ',0.0,-1)
c      CALL PLNUMB(999.,Y,CH, 1000.0*EI11A,0.0,4)
c      CALL PLMATH(999.,Y,CH,'#'   ,0.0,1)
c      CALL PLCHAR(999.,Y,CH, '10' ,0.0,2)
c      CALL PLMATH(999.,Y,CH,   '3',0.0,1)
c      CALL PLSUBS(X+4.0*CH,Y,CH,'o',0.0,1,PLCHAR)
c      Y = Y + YSPACE
cC
c      CALL PLMATH(X,Y,CH,'             ',0.0,-1)
c      CALL PLCHAR(X,Y,CH,'      area = ',0.0,-1)
c      CALL PLNUMB(999.,Y,CH, AREA,0.0, 5)
c      Y = Y + YSPACE
C
C--- Plot airfoil name over data list
      CALL PLCHAR(X+7.0*CH,Y,CH,NAME,0.0, 12)
C
      RETURN
      END ! GPARPL
C
C---------------------------------------------------------------------




      SUBROUTINE AECALCX(N,X,Y,T, ITYPE, 
     &                  AREA,XCEN,YCEN,EI11,EI22,APX1,APX2)
      DIMENSION X(*),Y(*),T(*)
C---------------------------------------------------------------
C     Calculates geometric properties of shape X,Y
C
C     Input:
C       N      number of points
C       X(.)   shape coordinate point arrays
C       Y(.)
C       T(.)   skin-thickness array, used only if ITYPE = 2
C       ITYPE  = 1 ...   integration is over whole area  dx dy
C              = 2 ...   integration is over skin  area   t ds
C
C     Output:
C       XCEN,YCEN  centroid location
C       EI11,EI22  principal moments of inertia
C       APX1,APX2  principal-axis angles
C---------------------------------------------------------------
      DATA PI / 3.141592653589793238 /
C
      SINT  = 0.0
      AINT  = 0.0
      XINT  = 0.0
      YINT  = 0.0
      XXINT = 0.0
      XYINT = 0.0
      YYINT = 0.0
C
      DO 10 IO = 1, N
        IF(IO.EQ.N) THEN
          IP = 1
        ELSE
          IP = IO + 1
        ENDIF
C
        DX =  X(IO) - X(IP)
        DY =  Y(IO) - Y(IP)
        XA = (X(IO) + X(IP))*0.50
        YA = (Y(IO) + Y(IP))*0.50
        TA = (T(IO) + T(IP))*0.50
C
        DS = SQRT(DX*DX + DY*DY)
        SINT = SINT + DS

        IF(ITYPE.EQ.1) THEN
C-------- integrate over airfoil cross-section
          DA = YA*DX
          AINT  = AINT  +       DA
          XINT  = XINT  + XA   *DA
          YINT  = YINT  + YA   *DA/2.0
          XXINT = XXINT + XA*XA*DA
          XYINT = XYINT + XA*YA*DA/2.0
          YYINT = YYINT + YA*YA*DA/3.0
        ELSE
C-------- integrate over skin thickness
          DA = TA*DS
          AINT  = AINT  +       DA
          XINT  = XINT  + XA   *DA
          YINT  = YINT  + YA   *DA
          XXINT = XXINT + XA*XA*DA
          XYINT = XYINT + XA*YA*DA
          YYINT = YYINT + YA*YA*DA
        ENDIF
C
 10   CONTINUE
C
      AREA = AINT
C
      IF(AINT .EQ. 0.0) THEN
        XCEN  = 0.0
        YCEN  = 0.0
        EI11  = 0.0
        EI22  = 0.0
        APX1 = 0.0
        APX2 = ATAN2(1.0,0.0)
        RETURN
      ENDIF
C
C
C---- calculate centroid location
      XCEN = XINT/AINT
      YCEN = YINT/AINT
C
C---- calculate inertias
      EIXX = YYINT - YCEN*YCEN*AINT
      EIXY = XYINT - XCEN*YCEN*AINT
      EIYY = XXINT - XCEN*XCEN*AINT
C
C---- set principal-axis inertias, EI11 is closest to "up-down" bending inertia
      EISQ  = 0.25*(EIXX - EIYY)**2  + EIXY**2
      SGN = SIGN( 1.0 , EIYY-EIXX )
      EI11 = 0.5*(EIXX + EIYY) - SGN*SQRT(EISQ)
      EI22 = 0.5*(EIXX + EIYY) + SGN*SQRT(EISQ)
C
      IF(EI11.EQ.0.0 .OR. EI22.EQ.0.0) THEN
C----- vanishing section stiffness
       APX1 = 0.0
       APX2 = ATAN2(1.0,0.0)
C
      ELSEIF(EISQ/(EI11*EI22) .LT. (0.001*SINT)**4) THEN
C----- rotationally-invariant section (circle, square, etc.)
       APX1 = 0.0
       APX2 = ATAN2(1.0,0.0)
C
      ELSE
C----- normal airfoil section
       C1 = EIXY
       S1 = EIXX-EI11
C
       C2 = EIXY
       S2 = EIXX-EI22
C
       IF(ABS(S1).GT.ABS(S2)) THEN
         APX1 = ATAN2(S1,C1)
         APX2 = APX1 + 0.5*PI
       ELSE
         APX2 = ATAN2(S2,C2)
         APX1 = APX2 - 0.5*PI
       ENDIF

       IF(APX1.LT.-0.5*PI) APX1 = APX1 + PI
       IF(APX1.GT.+0.5*PI) APX1 = APX1 - PI
       IF(APX2.LT.-0.5*PI) APX2 = APX2 + PI
       IF(APX2.GT.+0.5*PI) APX2 = APX2 - PI
C
      ENDIF
C
      RETURN
      END ! AECALCX
C
C-------------------------------------------------------------------------




      SUBROUTINE GETCAM (XCM,YCM,NCM,XTK,YTK,NTK,
     &                   X,XP,Y,YP,S,N )
C------------------------------------------------------
C     Finds camber and thickness 
C     distribution for input airfoil 
C------------------------------------------------------
      REAL XCM(*), YCM(*)
      REAL XTK(*), YTK(*)
      REAL X(*),XP(*),Y(*),YP(*),S(*)
C
      CALL XLFIND(SL,X,XP,Y,YP,S,N)
      XL = SEVAL(SL,X,XP,S,N)
      YL = SEVAL(SL,Y,YP,S,N)
C
C---- go over each point, finding opposite points, getting camber and thickness
      DO 10 I=1, N
C------ coordinates of point on the opposite side with the same x value
        CALL SOPPS(SOPP, S(I), X,XP,Y,YP,S,N,SL)
        XOPP = SEVAL(SOPP,X,XP,S,N)
        YOPP = SEVAL(SOPP,Y,YP,S,N)
C
C------ get camber and thickness
        XCM(I) = 0.5*(X(I)+XOPP)
        YCM(I) = 0.5*(Y(I)+YOPP)
        XTK(I) = 0.5*(X(I)+XOPP)
        YTK(I) = 0.5*(Y(I)-YOPP)
        YTK(I) = ABS(YTK(I))
c        if (XOPP.gt.0.9) then
c         write(*,*) 'cm i,x,y ',i,xcm(i),ycm(i)
c         write(*,*) 'tk i,x,y ',i,xtk(i),ytk(i)
c        endif
   10 CONTINUE
C
C---- Tolerance for nominally identical points
      TOL = 1.0E-5 * (S(N)-S(1))
ccc      TOL = 1.0E-3 * (S(N)-S(1))    ! Bad bug -- was losing x=1.0 point
C
C---- Sort the camber points
      NCM = N+1
      XCM(N+1) = XL
      YCM(N+1) = YL
      CALL SORTOL(TOL,NCM,XCM,YCM)
C
C--- Reorigin camber from LE so camberlines start at Y=0  4/24/01 HHY 
C    policy now to generate camber independent of Y-offsets 
      YOF = YCM(1)
      DO I = 1, NCM
        YCM(I) = YCM(I) - YOF
      END DO
C
C---- Sort the thickness points
      NTK = N+1
      XTK(N+1) = XL
      YTK(N+1) = 0.0
      CALL SORTOL(TOL,NTK,XTK,YTK)
C
      RETURN
      END ! GETCAM
C
C-----------------------------------------------------------------------



      SUBROUTINE SOPPS(SOPP, SI, X,XP,Y,YP,S,N, SLE)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C--------------------------------------------------
C     Calculates arc length SOPP of point 
C     which is opposite of point SI, on the 
C     other side of the airfoil baseline
C--------------------------------------------------
C
C---- reference length for testing convergence
      SLEN = S(N) - S(1)
C
C---- set chordline vector
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD = SQRT((XTE-XLE)**2 + (YTE-YLE)**2)
      DXC = (XTE-XLE) / CHORD
      DYC = (YTE-YLE) / CHORD
C
      IF(SI.LT.SLE) THEN
       IN = 1
       INOPP = N
      ELSE
       IN = N
       INOPP = 1
      ENDIF
      SFRAC = (SI-SLE)/(S(IN)-SLE)
      SOPP = SLE + SFRAC*(S(INOPP)-SLE)
C     
      IF(ABS(SFRAC) .LE. 1.0E-5) THEN
       SOPP = SLE
       RETURN
      ENDIF
C
C---- XBAR = x coordinate in chord-line axes
      XI  = SEVAL(SI , X,XP,S,N)
      YI  = SEVAL(SI , Y,YP,S,N)
      XLE = SEVAL(SLE, X,XP,S,N)
      YLE = SEVAL(SLE, Y,YP,S,N)
      XBAR = (XI-XLE)*DXC + (YI-YLE)*DYC
C
C---- converge on exact opposite point with same XBAR value
      DO 300 ITER=1, 12
        XOPP  = SEVAL(SOPP,X,XP,S,N)
        YOPP  = SEVAL(SOPP,Y,YP,S,N)
        XOPPD = DEVAL(SOPP,X,XP,S,N)
        YOPPD = DEVAL(SOPP,Y,YP,S,N)
C
        RES  = (XOPP -XLE)*DXC + (YOPP -YLE)*DYC - XBAR
        RESD =  XOPPD     *DXC +  YOPPD     *DYC
C
        IF(ABS(RES)/SLEN .LT. 1.0E-5) GO TO 305
        IF(RESD .EQ. 0.0) GO TO 303
C
        DSOPP = -RES/RESD
        SOPP = SOPP + DSOPP
C
        IF(ABS(DSOPP)/SLEN .LT. 1.0E-5) GO TO 305
 300  CONTINUE
 303  WRITE(*,*)
     &      'SOPPS: Opposite-point location failed. Continuing...'
      SOPP = SLE + SFRAC*(S(INOPP)-SLE)
C
 305  CONTINUE
      RETURN
      END ! SOPPS
C
C------------------------------------------------------------------------




      SUBROUTINE TCCALC(X,XP,Y,YP,S,N, 
     &                  THICK,XTHICK, CAMBR,XCAMBR )
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C---------------------------------------------------------------
C     Calculates max thickness and camber at airfoil points
C
C     Note: this routine does not find the maximum camber or 
C           thickness exactly as it only looks at discrete points
C
C     Input:
C       N      number of points
C       X(.)   shape coordinate point arrays
C       Y(.)
C
C     Output:
C       THICK  max thickness
C       CAMBR  max camber
C---------------------------------------------------------------
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD = SQRT((XTE-XLE)**2 + (YTE-YLE)**2)
C
C---- set unit chord-line vector
      DXC = (XTE-XLE) / CHORD
      DYC = (YTE-YLE) / CHORD
C
      THICK = 0.
      XTHICK = 0.
      CAMBR = 0.
      XCAMBR = 0.
C
C---- go over each point, finding the y-thickness and camber
      DO 30 I=1, N
        XBAR = (X(I)-XLE)*DXC + (Y(I)-YLE)*DYC
        YBAR = (Y(I)-YLE)*DXC - (X(I)-XLE)*DYC
C
C------ set point on the opposite side with the same chord x value
        CALL SOPPS(SOPP, S(I), X,XP,Y,YP,S,N, SLE)
        XOPP = SEVAL(SOPP,X,XP,S,N)
        YOPP = SEVAL(SOPP,Y,YP,S,N)
C
        YBAROP = (YOPP-YLE)*DXC - (XOPP-XLE)*DYC
C
        YC = 0.5*(YBAR+YBAROP)
        YT =  ABS(YBAR-YBAROP)
C
        IF(ABS(YC) .GT. ABS(CAMBR)) THEN
         CAMBR = YC
         XCAMBR = XOPP
        ENDIF
        IF(ABS(YT) .GT. ABS(THICK)) THEN
         THICK = YT
         XTHICK = XOPP
        ENDIF
   30 CONTINUE
C
      RETURN
      END ! TCCALC
C
C------------------------------------------------------------------



      SUBROUTINE GRDAIR(XGMIN,XGMAX, YGMIN,YGMAX,DXGN,DYGN,CHG,
     &                  LXAXIS,LYAXIS,
     &                  XOFF,XSF,YOFF,YSF, LMASK)
      LOGICAL LXAXIS,LYAXIS
C----------------------------------------
C     Plots grid with axes.
C     Intended for airfoil plot.
C----------------------------------------
C
      XMOD(XTMP) = XSF * (XTMP - XOFF)
      YMOD(YTMP) = YSF * (YTMP - YOFF)
C
      CALL NEWPEN(1)
C
C---- plot outline
      CALL PLOT(XMOD(XGMIN),YMOD(YGMIN),3)
      CALL PLOT(XMOD(XGMAX),YMOD(YGMIN),2)
      CALL PLOT(XMOD(XGMAX),YMOD(YGMAX),2)
      CALL PLOT(XMOD(XGMIN),YMOD(YGMAX),2)
      CALL PLOT(XMOD(XGMIN),YMOD(YGMIN),2)
C
      IF(LXAXIS)
     &  CALL XAXIS(XMOD(XGMIN),YMOD(YGMIN),(XGMAX-XGMIN)*XSF,
     &             DXGN*XSF, XGMIN,DXGN,CHG,-2)
      IF(LYAXIS)
     &  CALL YAXIS(XMOD(XGMIN),YMOD(YGMIN),(YGMAX-YGMIN)*YSF,
     &             DYGN*YSF, YGMIN,DYGN,CHG,-2)
C
C---- fine grid
      NXG = INT((XGMAX-XGMIN)/DXGN + 0.1)
      NYG = INT((YGMAX-YGMIN)/DYGN + 0.1)
      NXG = MAX(1,NXG)
      NYG = MAX(1,NYG)
C
      X0 = XMOD(XGMIN)
      Y0 = YMOD(YGMIN)
      DXG = (XMOD(XGMAX)-X0)/NXG
      DYG = (YMOD(YGMAX)-Y0)/NYG
C
      CALL NEWCOLORNAME('cyan')
      CALL PLGRID(X0,Y0,NXG,DXG,NYG,DYG, LMASK)
      CALL NEWCOLORNAME('black')
C
      RETURN
      END ! GRDAIR
C
C--------------------------------------------------------------------




      SUBROUTINE XLFIND(SLE,X,XP,Y,YP,S,N)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C------------------------------------------------------
C     Locates leftmost (minimum x) point location SLE
C
C     The defining condition is
C         
C      X' = 0     at  S = SLE
C
C     i.e. the surface tangent is vertical
C------------------------------------------------------
C
      DSLEN = S(N) - S(1)
C
C---- convergence tolerance
      DSEPS = (S(N)-S(1)) * 1.0E-5
C
C---- get first guess for SLE
      DO 10 I=3, N-2
        DX = X(I+1) - X(I)
        IF(DX .GT. 0.0) GO TO 11
   10 CONTINUE
C
   11 SLE = S(I)
C
C---- check for sharp LE case
      IF(S(I) .EQ. S(I-1)) THEN
ccc        WRITE(*,*) 'Sharp LE found at ',I,SLE
        RETURN
      ENDIF
C
C---- Newton iteration to get exact SLE value
      DO 20 ITER=1, 50
        DXDS = DEVAL(SLE,X,XP,S,N)
        DXDD = D2VAL(SLE,X,XP,S,N)
C
C------ drive DXDS to zero
        RES  = DXDS
        RESS = DXDD
C
C------ Newton delta for SLE 
        DSLE = -RES/RESS
C
        DSLE = MAX( DSLE , -0.01*ABS(DSLEN) )
        DSLE = MIN( DSLE ,  0.01*ABS(DSLEN) )
        SLE = SLE + DSLE
        IF(ABS(DSLE) .LT. DSEPS) RETURN
   20 CONTINUE
      WRITE(*,*) 'XLFIND:  Left point not found.  Continuing...'
      SLE = S(I)
      RETURN
      END ! XLFIND
C
C--------------------------------------------------------------------


 
      SUBROUTINE SORTOL(TOL,KK,S,W)
      DIMENSION S(KK), W(KK)
      LOGICAL DONE
C
C---- sort arrays
      DO IPASS=1, 1234
        DONE = .TRUE.
        DO N=1, KK-1
          NP = N+1
          IF(S(NP).LT.S(N)) THEN
           TEMP = S(NP)
           S(NP) = S(N)
           S(N) = TEMP
           TEMP = W(NP)
           W(NP) = W(N)
           W(N) = TEMP
           DONE = .FALSE.
          ENDIF
        END DO
        IF(DONE) GO TO 10
      END DO
      WRITE(*,*) 'Sort failed'
C
C---- search for near-duplicate pairs and eliminate extra points
C---- Modified 4/24/01 HHY to check list until ALL duplicates removed
C     This cures a bug for sharp LE foils where there were 3 LE points in
C     camber, thickness lists from GETCAM.
C
 10   KKS = KK
      DONE = .TRUE.
      DO 20 K=1, KKS
        IF(K.GE.KK) GO TO 20
        DSQ = (S(K)-S(K+1))**2 + (W(K)-W(K+1))**2
        IF(DSQ.GE.TOL*TOL) GO TO 20
C------- eliminate extra point pairs
ccc         write(*,*) 'extra on point ',k,kks
         KK = KK-1
         DO KT=K+1, KK
           S(KT) = S(KT+1)
           W(KT) = W(KT+1)
         END DO
         DONE = .FALSE.
   20 CONTINUE
      IF(.NOT.DONE) GO TO 10
C
      RETURN
      END  ! SORTOL
C
C-------------------------------------------------------------------------


 
      SUBROUTINE NORMX(NPX,NAFX,IAF,XT,YT,NT)
C----------------------------------------------------------------
C     Scales coordinates to get unit chord
C     Modified from xfoil NORM to handle multi-dim airfoil arrays
C----------------------------------------------------------------
      DIMENSION XT(NPX,NAFX),YT(NPX,NAFX)
      DIMENSION X(NPX),XP(NPX),Y(NPX),YP(NPX),S(NPX)
C
      N=NT
      DO I=1,N
         X(I)= XT(I,IAF)
         Y(I)= YT(I,IAF)
      ENDDO
C
      CALL SCALC (X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)
C
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
C
      XMAX = 0.5*(X(1) + X(N))
      XMIN = SEVAL(SLE,X,XP,S,N)
      YMIN = SEVAL(SLE,Y,YP,S,N)
C
      FUDGE = 1.0/(XMAX-XMIN)
      DO 40 I=1, N
        X(I) = (X(I)-XMIN)*FUDGE
        Y(I) = (Y(I)-YMIN)*FUDGE
C        S(I) = S(I)*FUDGE
   40 CONTINUE
C
      DO I=1,N
         XT(I,IAF)= X(I)
         YT(I,IAF)= Y(I)
      ENDDO
C
      WRITE(*,*)'Airfoil normalized'
C
      RETURN
      END  ! NORMX
C
C--------------------------------------------------------------------


 
      SUBROUTINE ROTATE(X,Y,N,ALFA)
      DIMENSION X(N), Y(N)
C
      SA = SIN(ALFA)
      CA = COS(ALFA)
CCC      XOFF = 0.25*(1.0-CA)
CCC      YOFF = 0.25*SA
      XOFF = 0.
      YOFF = 0.
      DO 8 I=1, N
        XT = X(I)
        YT = Y(I)
        X(I) = CA*XT + SA*YT + XOFF
        Y(I) = CA*YT - SA*XT + YOFF
    8 CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------------



      SUBROUTINE GETMAX(X,Y,YP,N,XMAX,YMAX)
      REAL X(*), Y(*), YP(*)
C------------------------------------------------
C     Calculates camber or thickness highpoint 
C     and x position
C------------------------------------------------
C
      XLEN = X(N) - X(1)
      XTOL = XLEN * 1.0E-5
C
      CALL SEGSPL(Y,YP,X,N)
C
C---- get approx max point and rough interval size
      YMAX0 = Y(1)
      XMAX0 = X(1)
      DO 5 I = 2, N
        IF (ABS(Y(I)).GT.ABS(YMAX0)) THEN
          YMAX0 = Y(I)
          XMAX0 = 0.5*(X(I-1) + X(I))
          DDX = 0.5*ABS(X(I+1) - X(I-1))
        ENDIF
 5    CONTINUE
      XMAX = XMAX0
C
C---- do a Newton loop to refine estimate
      DO 10 ITER=1, 10
        YMAX  = SEVAL(XMAX,Y,YP,X,N)
        RES   = DEVAL(XMAX,Y,YP,X,N)
        RESP  = D2VAL(XMAX,Y,YP,X,N)
        IF (ABS(XLEN*RESP) .LT. 1.0E-6) GO TO 20
          DX = -RES/RESP
          DX = SIGN( MIN(0.5*DDX,ABS(DX)) , DX)
          XMAX = XMAX + DX
          IF(ABS(DX) .LT. XTOL) GO TO 20
   10 CONTINUE
      WRITE(*,*)
     &  'GETMAX: Newton iteration for max camber/thickness failed.'
      YMAX = YMAX0
      XMAX = XMAX0
C
 20   RETURN
      END ! GETMAX
C
C----------------------------------------------------------------------



      SUBROUTINE SINVRTER(SI,XI,X,XS,S,N,ERROR)
      DIMENSION X(N),XS(N),S(N)
      LOGICAL ERROR
C----------------------------------------------------
C     Modified from SINVRT with error reporting
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
      ERROR= .FALSE.
      DO 10 ITER=1, 10
        CALL SEVALL(SI,X,XS,S,N, XX,XXS,XXSS)
        IF(XXS.EQ.0.0) GO TO 11
C
        DS = (XI-XX)/XXS
        SI = SI + DS
        IF(ABS(DS/(S(N)-S(1))) .LT. 1.0E-5) RETURN
 10   CONTINUE
C
 11   ERROR = .TRUE.
      WRITE(*,*) 'SINVRT: spline inversion failed.  Continuing...'
      RETURN
C
      END ! SINVRTER
C
C------------------------------------------------------------------
