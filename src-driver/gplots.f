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

      SUBROUTINE GPLINI(LOFSET,NEL,LPLTEL,
     &                  XMINE,XMAXE,YMINE,YMAXE,
     &                  XOFF,XSF,YOFF,YSF)
C-----------------------------------------------------
C     Sets up for geometry plot.
C     If LOFSET=T, initializes plot offsets,scales.
C     Plots geometry axes and grid.
C-----------------------------------------------------
      LOGICAL LOFSET, LPLTEL(*)
      DIMENSION XMINE(*),XMAXE(*),YMINE(*),YMAXE(*)
C
      EXTERNAL PLCHAR
      LOGICAL LBLANK
C
      INCLUDE 'PLOT.INC'
      INCLUDE 'MASKS.INC'
C
C---- Statement functions used to offset and scale all plots with blowups
      XMOD(XTMP) = XSF * (XTMP - XOFF)
      YMOD(YTMP) = YSF * (YTMP - YOFF)
C
C---- character width
      CHG = 0.5*CHGT
C
C---- additional left,right,bottom,top margins for annotations, etc.
      XADD1 = 0.02 + 3.0*CHG
      XADD2 = 0.04
      YADD1 = 0.04 + 2.0*CHG
      YADD2 = 0.04
C
C---- start new plot
      CALL PLTINI2
C
C---- re-origin to additional left,bottom margins
      CALL PLOT(XADD1,YADD1,-3)
C
C---- available plotting box
      XBOX = MIN(XWIND-XMARG , XPAGE-2.0*XMARG)/SIZE - XADD1 - XADD2
      YBOX = MIN(YWIND-YMARG , YPAGE-2.0*YMARG)/SIZE - YADD1 - YADD2
C
C---- set current plot limits as the min,max of all plotted elements
      XMIN =  1.0E24
      XMAX = -1.0E24
      YMIN =  1.0E24
      YMAX = -1.0E24
      LBLANK = .TRUE.
      DO IEL = 1, NEL
        IF(LPLTEL(IEL)) THEN
         XMIN = MIN(XMIN,XMINE(IEL))
         XMAX = MAX(XMAX,XMAXE(IEL))
         YMIN = MIN(YMIN,YMINE(IEL))
         YMAX = MAX(YMAX,YMAXE(IEL))
         LBLANK = .FALSE.
        ENDIF
      ENDDO
C
      IF(LBLANK) THEN
       XMIN = 0.0
       XMAX = 1.0
       YMIN = 0.0
       YMAX = 0.5
      ENDIF
C
C---- adjust them towards "nice" limits
      AXMIN = 0.4*XMIN
      AXMAX = 0.4*XMAX
      AYMIN = 0.4*YMIN
      AYMAX = 0.4*YMAX
      CALL AXISADJ2(AXMIN,AXMAX,XSPAN,XDEL,NTICS)
      CALL AXISADJ2(AYMIN,AYMAX,YSPAN,YDEL,NTICS)
C
C---- equalize annotation delta and reset limits to integral multiples of delta
      DEL = MAX(XDEL,YDEL)
      XMIN = DEL*AINT(XMIN/DEL + MIN(0.0,SIGN(0.999,XMIN)))
      XMAX = DEL*AINT(XMAX/DEL + MAX(0.0,SIGN(0.999,XMAX)))
      YMIN = DEL*AINT(YMIN/DEL + MIN(0.0,SIGN(0.999,YMIN)))
      YMAX = DEL*AINT(YMAX/DEL + MAX(0.0,SIGN(0.999,YMAX)))
C
C---- set deltas, and make sure there's at least one delta between limits
      XDEL = DEL
      YDEL = DEL
      XMAX = MAX( XMIN+XDEL , XMAX )
      YMAX = MAX( YMIN+YDEL , YMAX )
C
c      write(*,*) 'X  min,max', xmin, xmax
c      write(*,*) 'Y  min,max', ymin, ymax

      IF(LOFSET) THEN
C------ initialize offsets and scales (no blowup)
        XDIF = XMAX - XMIN
        YDIF = YMAX - YMIN
        IF(XDIF.EQ.0.0) XDIF = 1.0
        IF(YDIF.EQ.0.0) YDIF = 1.0
        XGFAC = 1.0/XDIF
        YGFAC = 1.0/YDIF
C
C------ set x,y scale factors
        XSF = XGFAC*XBOX
        YSF = YGFAC*YBOX
C
C------ equalize x,y scales to prevent distortion
        SSF = MIN(XSF,YSF)
        XSF = SSF
        YSF = SSF
C
C------ set offsets
        XOFF = XMIN
        YOFF = YMIN
      ENDIF
C
C-------------------------------------------------------------
      if(.false.) then

C---- current plot window limits
      XWMIN = XOFF
      XWMAX = XOFF + XBOX/XSF
      YWMIN = YOFF
      YWMAX = YOFF + YBOX/YSF
C
c      write(*,*) 'XW min,max', xwmin, xwmax
c      write(*,*) 'YW min,max', ywmin, ywmax

C---- scrunch plot limits down onto current plot window
      XMIN = MAX( XMIN , XWMIN )
      XMAX = MIN( XMAX , XWMAX )
      YMIN = MAX( YMIN , YWMIN )
      YMAX = MIN( YMAX , YWMAX )
      XMAX = MAX( XMIN+XDEL , XMAX )
      YMAX = MAX( YMIN+YDEL , YMAX )
C
c      write(*,*) 'X  min,max', xmin, xmax
c      write(*,*) 'Y  min,max', ymin, ymax

C---- set new "nice" limits after scrunching
      AXMIN = 0.4*XMIN
      AXMAX = 0.4*XMAX
      AYMIN = 0.4*YMIN
      AYMAX = 0.4*YMAX
      CALL AXISADJ2(AXMIN,AXMAX,XSPAN,XDEL,NTICS)
      CALL AXISADJ2(AYMIN,AYMAX,YSPAN,YDEL,NTICS)
C
C---- equalize annotation delta and reset limits to integral multiples of delta
      DEL = MAX(XDEL,YDEL)
      XMIN = DEL*AINT(XMIN/DEL + MIN(0.0,SIGN(0.999,XMIN)))
      XMAX = DEL*AINT(XMAX/DEL + MAX(0.0,SIGN(0.999,XMAX)))
      YMIN = DEL*AINT(YMIN/DEL + MIN(0.0,SIGN(0.999,YMIN)))
      YMAX = DEL*AINT(YMAX/DEL + MAX(0.0,SIGN(0.999,YMAX)))
C
      XDEL = DEL
      YDEL = DEL
C
C---- set new offsets to put lower-left corner at origin
      XOFF = XMIN
      YOFF = YMIN
C
      endif
C-------------------------------------------------------------
C
C---- set plot limits in offset/scaled coordinates
      X0 = XMOD(XMIN)
      X1 = XMOD(XMAX)
      Y0 = YMOD(YMIN)
      Y1 = YMOD(YMAX)
C
C---- annotation deltas in offset/scaled coordinates
      XD = XMOD(XMIN+XDEL) - XMOD(XMIN)
      YD = YMOD(YMIN+YDEL) - YMOD(YMIN)
C
      CALL NEWPEN(1)
      CALL XAXIS(X0,Y0,X1-X0,XD,XMIN,XDEL, CHG,-2)
      CALL YAXIS(X0,Y0,Y1-Y0,YD,YMIN,YDEL, CHG,-2)
C
      NXG = 2   * INT( (XMAX-XMIN)/XDEL + 0.01 )
      DXG = 0.5 * XD
C
      NYG = 2   * INT( (YMAX-YMIN)/YDEL + 0.01 )
      DYG = 0.5 * YD
C
      CALL PLGRID(X0,Y0,NXG,DXG,NYG,DYG,LMASK2)
C
C----- plot the axis for axisymmetric geometry
      CALL NEWPEN(2)
      CALL NEWPAT(LMASK1)
      YAX = YMOD(0.0)
      CALL PLOT(X0,YAX,3)
      CALL PLOT(X1,YAX,2)
      CALL NEWPAT(LMASK0)
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! GPLINI



      SUBROUTINE GVPLOT(LPLTEL,LSEGPL)
C------------------------------------------------
C     Plots current geometry XP,YP
C------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      LOGICAL LPLTEL(*), LSEGPL
C
      INTEGER KCOLEL(NEX)
C
C---- Statement functions used to offset and scale all plots with blowups
      XMOD(XTMP) = XPSF * (XTMP - XPOFF)
      YMOD(YTMP) = YPSF * (YTMP - YPOFF)
C
C---- character size for element numbers
      CHN = 0.5*CHGT
C
C---- number of spline sub-intervals
      NTGSPL = 20
C
      CALL GETCOLOR(ICOL0)
C
      DO IEL = 1, NEL
        KCOLEL(IEL) = ICOLEL(IEL)
      ENDDO
C
 5    CONTINUE
C
C---- plot individual elements
      DO 10 IEL = 1, NEL
        IF(.NOT.LPLTEL(IEL)) GO TO 10
C
        I = IPFRST(IEL)
        N = IPLAST(IEL) - IPFRST(IEL) + 1
C
C------ set color for this element
        CALL NEWCOLOR(KCOLEL(IEL))
C
C------ plot airfoil contour using spline sub-intervals
        CALL NEWPEN(2)
        CALL PLTGSP(XP(I),XPS(I),YP(I),YPS(I),SP(I),N, 
     &              XPOFF,XPSF,YPOFF,YPSF, NTGSPL)
C
C------ plot point singularities, or panel nodes if requested
        IF(N.EQ.1 .OR. LVPLT) THEN
         CALL NEWPEN(2)
         SSIZ = SSIZEL(NETYPE(IEL))
         DO I = IPFRST(IEL), IPLAST(IEL)
           CALL PLSYMB(XMOD(XP(I)),YMOD(YP(I)),SSIZ,1,0.0,0)
         ENDDO
        ENDIF
C
C------ plot element numbers
        IF(LEPLT) THEN
         CALL NEWPEN(2)
         XNUM = XMOD(XELNUM(IEL)) + SSIZEL(NETYPE(IEL)) + 0.25*CHN
         YNUM = YMOD(YELNUM(IEL)) + SSIZEL(NETYPE(IEL)) + 0.25*CHN
         REL = FLOAT(IEL)
         CALL PLNUMB(XNUM,YNUM,CHN,REL,0.0,-1)
        ENDIF
 10   CONTINUE
C
      IF(LSEGPL) THEN
C----- overplot lines inside segments with black lines
       CALL NEWCOLOR(ICOL0)
       CALL NEWPEN(2)
       DO 20 IEL = 1, NEL
         IF(.NOT.LPLTEL(IEL)) GO TO 20
C
         DO ISEG = 1, NSEG
           IF(IPSEG1(ISEG).GE.IPFRST(IEL) .AND.
     &        IPSEG2(ISEG).LE.IPLAST(IEL)       ) THEN
C
            I = IPSEG1(ISEG)
            N = IPSEG2(ISEG) - IPSEG1(ISEG) + 1
C
            CALL PLTGSP(XP(I),XPS(I),YP(I),YPS(I),SP(I),N, 
     &                  XPOFF,XPSF,YPOFF,YPSF, NTGSPL)
           ENDIF
         ENDDO
 20    CONTINUE
      ENDIF
C
      CALL NEWCOLOR(ICOL0)
C
      CALL PLFLUSH
      LGSPPL = .TRUE.
      LQSPPL = .FALSE.
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! GVPLOT



      SUBROUTINE GBPLOT(LPLTEL)
C------------------------------------------------
C     Plots current buffer geometry XB,YB
C------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      LOGICAL LPLTEL(*)
C
      INCLUDE 'MASKS.INC'
C
C---- Statement functions used to offset and scale all plots with blowups
      XMOD(XTMP) = XBSF * (XTMP - XBOFF)
      YMOD(YTMP) = YBSF * (YTMP - YBOFF)
C
C---- character size for element numbers
      CHN = 0.5*CHGT
C
C---- size for reference-point symbol
      SHR = 0.5*CHGT
C
C---- number of spline sub-intervals
      NTGSPL = 20
C
      CALL GETPAT(IPAT0)
      CALL GETCOLOR(ICOL0)
C
C---- plot individual elements
      DO 10 IEL = 1, NBEL
        IF(.NOT.LPLTEL(IEL)) GO TO 10
C
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
C
C------ set color for this element
        CALL NEWCOLOR(ICOLEL(IEL))
C
        IF(LHOMPL) THEN
C------- plot home position
         CALL NEWPEN(1)
         CALL NEWPAT(LMASK3)
         CALL PLHOME(IEL)
         CALL NEWPAT(IPAT0)
        ENDIF
C
C------ plot airfoil contour using spline sub-intervals
        CALL NEWPEN(2)
        CALL PLTGSP(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, 
     &              XBOFF,XBSF,YBOFF,YBSF, NTGSPL)
C
        IF(N.EQ.1 .OR. LVPLT) THEN
C------- plot point singularities, or panel nodes if requested
         CALL NEWPEN(2)
         SSIZ = SSIZEL(NETYPE(IEL))
         DO I = IBFRST(IEL), IBLAST(IEL)
           CALL PLSYMB(XMOD(XB(I)),YMOD(YB(I)),SSIZ,1,0.0,0)
         ENDDO
        ENDIF
C
        IF(LEPLT) THEN
C------- plot element number
         CALL NEWPEN(2)
         XNUM = XMOD(XELNUM(IEL)) + SSIZEL(NETYPE(IEL)) + 0.25*CHN
         YNUM = YMOD(YELNUM(IEL)) + SSIZEL(NETYPE(IEL)) + 0.25*CHN
         REL = FLOAT(IEL)
         CALL PLNUMB(XNUM,YNUM,CHN,REL,0.0,-1)
        ENDIF
C
c        IF(LRPLT) THEN
        IF(IEL.EQ.IELGDES) THEN
C------- plot refernce point
         CALL NEWPEN(2)
         XREF = XMOD(XBREFE(IEL))
         YREF = YMOD(YBREFE(IEL))
         REL = FLOAT(IEL)
         CALL PLSYMB(XREF,YREF,    SHR,1,0.0,0)
         CALL PLSYMB(XREF,YREF,1.5*SHR,3,0.0,0)
        ENDIF
c        ENDIF
C
 10   CONTINUE
C
      CALL NEWCOLOR(ICOL0)
C
      CALL PLFLUSH
      LGPLOT = .TRUE.
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! GBPLOT



      SUBROUTINE PLHOME(IEL)
      INCLUDE 'DFDC.INC'
C
      PARAMETER (KK = 2)
      PARAMETER (ITX = KK*IBX)
      DIMENSION XT(ITX), YT(ITX)
C
C---- Statement functions used to offset and scale all plots with blowups
      XMOD(XTMP) = XBSF * (XTMP - XBOFF)
      YMOD(YTMP) = YBSF * (YTMP - YBOFF)
C
C
      IB1 = IBFRST(IEL)
      IB2 = IBLAST(IEL)
C
      I = IB1
      N = IB2 - IB1 + 1
C
C---- accumulate rotation
      SA = SIN(-AGBSUM(IEL)*PI/180.0)
      CA = COS(-AGBSUM(IEL)*PI/180.0)
C
C---- set home geometry
      IT = 1
      IB = IB1
      XK = XB(IB) - DXBSUM(IEL)
      YK = YB(IB) - DYBSUM(IEL)
      XT(IT) = XMOD( (CA*XK + SA*YK)/XFBSUM(IEL) )
      YT(IT) = YMOD( (CA*YK - SA*XK)/YFBSUM(IEL) )
C
      DO IB = IB1, IB2-1
        DSB = SB(IB+1) - SB(IB)
C
        DO K = 1, KK
          IT = IT+1
          SK = SB(IB) + DSB*FLOAT(K)/FLOAT(KK)
          XK = SEVAL(SK,XB(I),XBS(I),SB(I),N) - DXBSUM(IEL)
          YK = SEVAL(SK,YB(I),YBS(I),SB(I),N) - DYBSUM(IEL)
          XT(IT) = XMOD( (CA*XK + SA*YK)/XFBSUM(IEL) )
          YT(IT) = YMOD( (CA*YK - SA*XK)/YFBSUM(IEL) )
        ENDDO
      ENDDO
      NT = IT
C
      CALL POLYLINE(XT,YT,NT,0)
C
cc      CALL XYLINE(NT,XT,YT,0.0,1.0,0.0,1.0,4)
C
c      I = IB1
c      N = IB2 - IB1 + 1
c      CALL SCALC(XT(I),YT(I),ST(I),N)
c      CALL SEGSPL(XT(I),XTS(I),ST(I),N)
c      CALL SEGSPL(YT(I),YTS(I),ST(I),N)
c      CALL PLTGSP(XT(I),XTS(I),YT(I),YTS(I),ST(I),N, 
c     &            XBOFF,XBSF,YBOFF,YBSF, 2)
C
c      IB = IB1
c      CALL PLOT(XMOD(XT(IB)),YMOD(YT(IB)),3)
c      DO IB = IB1+1, IB2
c        CALL PLOT(XMOD(XT(IB)),YMOD(YT(IB)),2)
c      ENDDO
C
      RETURN
      END ! PLHOME


