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


      SUBROUTINE QSPLIM(IEL)
C---------------------------------------------------
C     Sets up Qspec(s) plot limits for element IEL
C---------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION QTMP(2), STMP(2)
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
      I = IP1
      N = IP2 - IP1 + 1
C
C---- set limits for QSGAM, and all QSPEC distributions
      CALL SETLIM(N,QSGAM(I),QTMP(1),QTMP(2),QDEL,NANN)
      DO KQSP = 1, NQSP
        CALL SETLIM(N,QSPEC(I,KQSP),QMIN,QMAX,QDEL,NANN)
        QTMP(1) = MIN(QTMP(1),QMIN)
        QTMP(2) = MAX(QTMP(2),QMAX)
      ENDDO
      CALL SETLIM(2,QTMP,QSPMIN(IEL),QSPMAX(IEL),QSPDEL(IEL),NANN)
      IF(NANN.LE.4) QSPDEL(IEL) = 0.5*QSPDEL(IEL)
C
C---- set limits for SSPEC
      STMP(1) = SSPEC(IP1)
      STMP(2) = SSPEC(IP2)
      CALL SETLIM(2,STMP,SSPMIN(IEL),SSPMAX(IEL),SSPDEL(IEL),NANN)
C
c      write(*,*) 'QSPLIM(IEL) ',IEL
c      write(*,*) 'QSPMIN      ',QSPMIN(IEL)
c      write(*,*) 'QSPMAX      ',QSPMAX(IEL)
c      write(*,*) 'QSPDEL      ',QSPDEL(IEL)
c      write(*,*) 'SSPMIN      ',SSPMIN(IEL)
c      write(*,*) 'SSPMAX      ',SSPMAX(IEL)
c      write(*,*) 'SSPDEL      ',SSPDEL(IEL)
C
      RETURN
      END



      SUBROUTINE QPLINI(LOFSET,NEL,IELPLT, LXREV,LYREV,
     &                  XMINE,XMAXE,YMINE,YMAXE,
     &                  XOFF,XSF,YOFF,YSF)
C-----------------------------------------------------
C     Plots axes and grid overlay for Qspec(s) plot.
C
C     If NGRID=0, s spacing is SDEL
C     If NGRID>0, s spacing is implicitly determined
C         from NGX equal increments of Fgrid(s)
C-----------------------------------------------------
      LOGICAL LOFSET
      LOGICAL LXREV, LYREV
      DIMENSION XMINE(*),XMAXE(*),YMINE(*),YMAXE(*)
C
      INCLUDE 'PLOT.INC'
      INCLUDE 'MASKS.INC'
C
C---- Statement functions used to offset and scale all plots with blowups
      XMOD(XTMP) = XSF * (XTMP - XOFF)
      YMOD(YTMP) = YSF * (YTMP - YOFF)
C
C---- character widths
      CHA = 0.6*CHGT
      CHL = 0.8*CHGT
      CHG = 1.1*CHGT
C
C---- additional left,right,bottom,top margins for annotations, etc.
      XADD1 = 0.03 + 4.0*CHA
      XADD2 = 0.05
      YADD1 = 0.03 + 2.0*CHA
      YADD2 = 0.05
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
C---- current plot limits
      IEL = IELPLT
C
      XMIN = XMINE(IEL)
      XMAX = XMAXE(IEL)
C
      YMIN = YMINE(IEL)
      YMAX = YMAXE(IEL)
C
C---- adjust them towards "nice" limits
      AXMIN = 0.5*XMIN
      AXMAX = 0.5*XMAX
      AYMIN = 0.5*YMIN
      AYMAX = 0.5*YMAX
      CALL AXISADJ2(AXMIN,AXMAX,XSPAN,XDEL,NTICS)
      CALL AXISADJ2(AYMIN,AYMAX,YSPAN,YDEL,NTICS)
C
C---- reset limits to integral multiples of each delta
      XMIN = XDEL*AINT(XMIN/XDEL + MIN(0.0,SIGN(0.999,XMIN)))
      XMAX = XDEL*AINT(XMAX/XDEL + MAX(0.0,SIGN(0.999,XMAX)))
      YMIN = YDEL*AINT(YMIN/YDEL + MIN(0.0,SIGN(0.999,YMIN)))
      YMAX = YDEL*AINT(YMAX/YDEL + MAX(0.0,SIGN(0.999,YMAX)))
C
      IF(LXREV) THEN
C----- plot x in reverse direction
       TMP = XMIN
       XMIN = XMAX
       XMAX = TMP
       XDEL = -XDEL
      ENDIF
C
      IF(LYREV) THEN
C----- plot y in reverse direction
       TMP = YMIN
       YMIN = YMAX
       YMAX = TMP
       YDEL = -YDEL
      ENDIF
C
      IF(LOFSET) THEN
C------ initialize offsets and scales (no blowup)
        XDIF = XMAX - XMIN
        YDIF = YMAX - YMIN
        IF(XDIF.EQ.0.0) XDIF = 1.0
        IF(YDIF.EQ.0.0) YDIF = 1.0
        XFAC = 1.0/XDIF
        YFAC = 1.0/YDIF
C
        XSF = XFAC*XBOX
        YSF = YFAC*YBOX
C
        XOFF = XMIN
        YOFF = YMIN
      ENDIF
C
C-------------------------------------------------------------
C---- current plot window limits
      XWMIN = XOFF
      XWMAX = XOFF + XBOX/XSF
      YWMIN = YOFF
      YWMAX = YOFF + YBOX/YSF
C
C---- put the window and plot limits in the usual (non-reversed direction)
      IF(LXREV) THEN
       TMP = XWMIN
       XWMIN = XWMAX
       XWMAX = TMP
C
       TMP = XMIN
       XMIN = XMAX
       XMAX = TMP
      ENDIF
C
      IF(LYREV) THEN
       TMP = YWMIN
       YWMIN = YWMAX
       YWMAX = TMP
C
       TMP = YMIN
       YMIN = YMAX
       YMAX = TMP
      ENDIF
C
C---- scrunch plot limits down onto current plot window
      XMIN = MAX( XMIN , XWMIN )
      XMAX = MIN( XMAX , XWMAX )
      YMIN = MAX( YMIN , YWMIN )
      YMAX = MIN( YMAX , YWMAX )
C
C---- set new "nice" limits after scrunching
      AXMIN = 0.5*XMIN
      AXMAX = 0.5*XMAX
      AYMIN = 0.5*YMIN
      AYMAX = 0.5*YMAX
      CALL AXISADJ2(AXMIN,AXMAX,XSPAN,XDEL,NTICS)
      CALL AXISADJ2(AYMIN,AYMAX,YSPAN,YDEL,NTICS)
C
      XMIN = XDEL*AINT(XMIN/XDEL + MIN(0.0,SIGN(0.999,XMIN)))
      XMAX = XDEL*AINT(XMAX/XDEL + MAX(0.0,SIGN(0.999,XMAX)))
      YMIN = YDEL*AINT(YMIN/YDEL + MIN(0.0,SIGN(0.999,YMIN)))
      YMAX = YDEL*AINT(YMAX/YDEL + MAX(0.0,SIGN(0.999,YMAX)))
C
C---- back to reverse again
      IF(LXREV) THEN
       TMP = XMIN
       XMIN = XMAX
       XMAX = TMP
       XDEL = -XDEL
      ENDIF
C
      IF(LYREV) THEN
       TMP = YMIN
       YMIN = YMAX
       YMAX = TMP
       YDEL = -YDEL
      ENDIF
C
C---- set new offsets to ensure lower-left corner is at origin
      XOFF = XMIN
      YOFF = YMIN
C-------------------------------------------------------------
C
C---- set plot limits in offset/scaled coordinates
      X0 = XMOD(XMIN)
      X1 = XMOD(XMAX)
      Y0 = YMOD(YMIN)
      Y1 = YMOD(YMAX)
C
C---- annotation deltas in offset/scaled coordinates
      XD = ABS( XMOD(XMIN+XDEL) - XMOD(XMIN) )
      YD = ABS( YMOD(YMIN+YDEL) - YMOD(YMIN) )
C
C---- plot bottom and top lines of plot box
      CALL NEWPEN(1)
      CALL PLOT(X0,Y0,3)
      CALL PLOT(X1,Y0,2)
      CALL PLOT(X0,Y1,3)
      CALL PLOT(X1,Y1,2)
C
      CALL NEWPEN(2)
C
C---- plot x-axis at bottom of box
      CALL XAXIS(X0,Y0,X1-X0,XD,XMIN,XDEL, CHA,-2)
C
C---- plot dense dotted line on y=0 if it is inside current y limits
      YZ = YMOD(0.0)
      IF(YZ.GT.Y0 .AND. YZ.LT.Y1) THEN
       CALL GETPAT(IPAT0)
       CALL NEWPAT(LMASK1)
       CALL PLOT(X0,YZ,3)
       CALL PLOT(X1,YZ,2)
       CALL NEWPAT(IPAT0)
      ENDIF
C
C---- plot y-axes on left and right sides of plot box
      CALL YAXIS(X0,Y0,Y1-Y0,YD,YMIN,YDEL, CHA,-2)
      CALL YAXIS(X1,Y0,Y1-Y0,YD,YMIN,YDEL,-CHA,-3)
C
C
      CALL NEWPEN(1)
C
      NXG = 2   * INT( (XMAX-XMIN)/XDEL + 0.01 )
      DXG = 0.5 * XD
C
      NYG = 2   * INT( (YMAX-YMIN)/YDEL + 0.01 )
      DYG = 0.5 * YD
C
      CALL PLGRID(X0,Y0,NXG,DXG,NYG,DYG,LMASK2)
C
C
cC---- plot sonic lines if within range
c      IF( QDASH.LE.YMAX) THEN
c       CALL YDASH(XQMOD(XMIN),XQMOD(XMAX),YQMOD( QDASH))
c      ENDIF
c      IF(-QDASH.GE.YMIN) THEN
c       CALL YDASH(XQMOD(XMIN),XQMOD(XMAX),YQMOD(-QDASH))
c      ENDIF
C
C
      CALL NEWPEN(3)
C
      XPLT = X1 - 1.5*XD - 0.5*CHL
      YPLT = Y0          - 2.0*CHL
      CALL PLMATH(XPLT,YPLT,CHL,'^',0.0,1)
      CALL PLCHAR(XPLT,YPLT,CHL,'s',0.0,1)
C
      XPLT = X0          - 2.5*CHG
      YPLT = Y1 - 1.5*YD - 0.5*CHG
      CALL PLMATH(XPLT,YPLT,CHG,'g',0.0,1)
C
      RETURN
      END ! QPLINI



      SUBROUTINE QSPLOT
C------------------------------------------------
C     Plots Gamma(s) and Qspec(s) distributions.
C------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
C---- Statement functions used to offset and scale all plots with blowups
      XQMOD(XTMP) = XQSF * (XTMP - XQOFF)
      YQMOD(YTMP) = YQSF * (YTMP - YQOFF)
C
C---- statement function for compressible Karman-Tsien velocity
cc      QCOMP(G) = G*(1.0-TKLAM) / (1.0 - TKLAM*(G/QINF)**2)
      QCOMP(G) = G
C
C---- plot-symbol size (grow it weakly with blowup factor)
      SHT = SHGT * SQRT(ABS(XQSF))
C
      CALL GETCOLOR(ICOL0)
C
      IF(LQSLOP) THEN
       NTQSPL = 8
      ELSE
       NTQSPL = 1
      ENDIF
C
C
      IEL = IELQDES
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
      CALL NEWCOLOR(ICOLEL(IEL))
      DO IP = IP1, IP2
C------ turn off color at the start of any segment
        DO ISEG = 1, NSEG
          IF(IP.EQ.IPSEG1(ISEG)) CALL NEWCOLOR(ICOL0)
        ENDDO
C
        XPLT =       SSPEC(IP)
        YPLT = QCOMP(QSGAM(IP))
        CALL PLSYMB(XQMOD(XPLT),YQMOD(YPLT),SHT,3,0.,0)
C
C------ turn on element color at the end of any segment
        DO ISEG = 1, NSEG
          IF(IP.EQ.IPSEG2(ISEG)) CALL NEWCOLOR(ICOLEL(IEL))
        ENDDO
      ENDDO
C
C---- plot individual Qspec lines
      CALL NEWCOLOR(ICOLEL(IEL))
      DO KQSP = 1, NQSP
        CALL QSPPLT(IP1,IP2,KQSP,NTQSPL)
      ENDDO
C
C---- overplot lines inside active segments without color
      CALL NEWCOLOR(ICOL0)
      DO ISEG = 1, NSEG
        DO KQSP = 1, NQSP
          CALL QSPPLT(IPSEG1(ISEG),IPSEG2(ISEG),KQSP,NTQSPL)
        ENDDO
      ENDDO
C
c      IF(LQSVIS) THEN
c       CALL NEWCOLORNAME('MAGENTA')
c       DO IP = IP1+1, IP2
c         DSP = SP(IP) - SP(IP-1)
c         DQV = QCOMP(QVIS(IP)) - QCOMP(QVIS(IP-1))
c         SP1 = (S(IP-1) + 0.25*DSP)/SP(NV)
c         SP2 = (S(IP)   - 0.25*DSP)/SP(NV)
c         QV1 = QCOMP(QVIS(IP-1)) + 0.25*DQV
c         QV2 = QCOMP(QVIS(IP)  ) - 0.25*DQV
c         CALL PLOT(XQMOD(SP1),YQMOD(QV1),3)
c         CALL PLOT(XQMOD(SP2),YQMOD(QV2),2)
c       ENDDO
c       CALL NEWCOLOR(ICOL0)
c      ENDIF
C
      CALL PLFLUSH
      LQSPPL = .TRUE.
      LGSPPL = .FALSE.
C
C---- plot tick marks bounding all segments on target distribution
      KQSP = KQTARG
      CALL NEWCOLORNAME('MAGENTA')
      DO ISEG = 1, NSEG
        IPS1 = IPSEG1(ISEG)
        IPS2 = IPSEG2(ISEG)
        XPLT1 =       SSPEC(IPS1)
        XPLT2 =       SSPEC(IPS2)
        YPLT1 = QCOMP(QSPEC(IPS1,KQSP))
        YPLT2 = QCOMP(QSPEC(IPS2,KQSP))
        CALL PLOT(XQMOD(XPLT1),YQMOD(YPLT1)-0.03,3)
        CALL PLOT(XQMOD(XPLT1),YQMOD(YPLT1)+0.03,2)
        CALL PLOT(XQMOD(XPLT2),YQMOD(YPLT2)-0.03,3)
        CALL PLOT(XQMOD(XPLT2),YQMOD(YPLT2)+0.03,2)
      ENDDO
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! QSPLOT



      SUBROUTINE QSPPLT(IPPL1,IPPL2,KQSP,NT)
C----------------------------------------------------
C     Plots Qspec(s) via spline, 
C     between the two specified indices.
C     NT is the number of spline sub-intervals
C----------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Statement functions used to offset and scale all plots with blowups
      XQMOD(XTMP) = XQSF * (XTMP - XQOFF)
      YQMOD(YTMP) = YQSF * (YTMP - YQOFF)
C
C---- statement function for compressible Karman-Tsien velocity
ccc   QCOMP(G) = G*(1.0-TKLAM) / (1.0 - TKLAM*(G/QINF)**2)
      QCOMP(G) = G
C
      I = IPPL1
      N = IPPL2 - IPPL1 + 1
C
      DO IP = IPPL1+1, IPPL2
        DS = SSPEC(IP) - SSPEC(IP-1)
C
        XPLT =       SSPEC(IP-1)
        YPLT = QCOMP(QSPEC(IP-1,KQSP))
        CALL PLOT(XQMOD(XPLT),YQMOD(YPLT),3)
C
        DO IT=1, NT
          SSPT = SSPEC(IP-1) + DS*FLOAT(IT)/FLOAT(NT)
          QSPT = SEVAL(SSPT,QSPEC(I,KQSP),QSPECS(I,KQSP),SSPEC(I),N)
          XPLT =       SSPT
          YPLT = QCOMP(QSPT)
          CALL PLOT(XQMOD(XPLT),YQMOD(YPLT),2)
        ENDDO
      ENDDO
C
      RETURN
      END ! QSPPLT


 
      SUBROUTINE IQSGET(IEL,ISEG)
C------------------------------------------------------------
C     Sets target segment endpoint indices 
C     from cursor input on q(s) plot
C------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION IPSNEW(2)
      CHARACTER*1 KCHAR
C
C---- Statement functions used to offset and scale all plots with blowups
      XQMOD(XTMP) = XQSF * (XTMP - XQOFF)
      YQMOD(YTMP) = YQSF * (YTMP - YQOFF)
C
C---- statement function for compressible Karman-Tsien velocity
ccc   QCOMP(G) = G*(1.0-TKLAM) / (1.0 - TKLAM*(G/QINF)**2)
      QCOMP(G) = G
C
C
      CALL GETCOLOR(ICOL0)
C
      IPSNEW(1) = 0
      IPSNEW(2) = 0
C
      WRITE(*,1200) ISEG
 1200 FORMAT(/' Mark segment',I4,'  endpoints with cursor...'/)
C
      DO 10 IE = 1, 2
C------ get cursor location from user
        CALL GETCURSORXY(XE,YE,KCHAR)
        DSQMIN = 1.0E24
        IPSNEW(IE) = IPFRST(IEL)
        KQMIN = 1
C
C------ search all Qspec lines only for first selected point
        IF(IE.EQ.1) THEN
         KQSP1 = 1
         KQSPN = NQSP
        ELSE
         KQSP1 = KQTARG
         KQSPN = KQTARG
        ENDIF
C
C------ find plot point closest to cursor point
        DO KQSP = KQSP1, KQSPN
          DO IP = IPFRST(IEL), IPLAST(IEL)
            GCOMP = QCOMP(QSPEC(IP,KQSP))
            XPNT = XQMOD(SSPEC(IP))
            YPNT = YQMOD(GCOMP)
            DSQ = (XE - XPNT)**2 + (YE - YPNT)**2
            IF(DSQ.LE.DSQMIN) THEN
             DSQMIN = DSQ
             IPSNEW(IE) = IP
             KQMIN = KQSP
            ENDIF
          ENDDO
        ENDDO
C
C------ nearest point to first clicked point sets target line
        IF(IE.EQ.1) KQTARG = KQMIN
C
        CALL NEWCOLORNAME('MAGENTA')
        IP = IPSNEW(IE)
        QSCOMP = QCOMP(QSPEC(IP,KQTARG))
        CALL PLOT(XQMOD(SSPEC(IP)),YQMOD(QSCOMP)-0.03,3)
        CALL PLOT(XQMOD(SSPEC(IP)),YQMOD(QSCOMP)+0.03,2)
        CALL NEWCOLOR(ICOL0)
        CALL PLFLUSH
   10 CONTINUE
C
      IF(IPSNEW(1).EQ.IPSNEW(2)) THEN
       WRITE(*,*) '**  Endpoints must be distinct  **'
       WRITE(*,*) '**  NEW SEGMENT NOT MARKED OFF  **'
       RETURN
      ENDIF
C
      IPSEG1(ISEG) = MIN(IPSNEW(1),IPSNEW(2))
      IPSEG2(ISEG) = MAX(IPSNEW(1),IPSNEW(2))
      IELSEG(ISEG) = IEL
C
      RETURN
      END ! IQSGET

 
      SUBROUTINE IGSGET(LPLTEL,ISEG)
C------------------------------------------------------------
C     Sets target segment endpoint indices 
C     from cursor input on geometry plot.
C------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      LOGICAL LPLTEL(*)
C
      DIMENSION IPSNEW(2)
      CHARACTER*1 KCHAR
C
C---- Statement functions used to offset and scale all plots with blowups
      XMOD(XTMP) = XPSF * (XTMP - XPOFF)
      YMOD(YTMP) = YPSF * (YTMP - YPOFF)
C
C---- plot-symbol size
      SSIZ = 0.5*CHGT
C
      CALL GETCOLOR(ICOL0)
C
      IPSNEW(1) = 0
      IPSNEW(2) = 0
C
      WRITE(*,1200) ISEG
 1200 FORMAT(/' Mark segment',I4,'  endpoints with cursor...'/)
C
      DO 10 IE = 1, 2
C------ get cursor location from user
        CALL GETCURSORXY(XE,YE,KCHAR)
        DSQMIN = 1.0E24
        IELMIN = 1
        IPSNEW(IE) = IPFRST(IELMIN)
C
C------ search all elements lines only for first selected point
        IF(IE.EQ.1) THEN
         KEL1 = 1
         KEL2 = NEL
        ELSE
         KEL1 = IELMIN
         KEL2 = IELMIN
        ENDIF
C
C------ find plot point closest to cursor point
        DO 8 KEL = KEL1, KEL2
          IF(.NOT.LPLTEL(KEL)) GO TO 8
C
          IF(NETYPE(KEL).EQ.0 .OR.
     &       NETYPE(KEL).EQ.1 .OR.
     &       NETYPE(KEL).EQ.2     ) THEN
           DO IP = IPFRST(KEL), IPLAST(KEL)
             XPNT = XMOD(XP(IP))
             YPNT = YMOD(YP(IP))
             DSQ = (XE - XPNT)**2 + (YE - YPNT)**2
             IF(DSQ.LE.DSQMIN) THEN
              DSQMIN = DSQ
              IPSNEW(IE) = IP
              IELMIN = KEL
             ENDIF
           ENDDO
          ENDIF
 8      CONTINUE
C
        CALL NEWCOLORNAME('MAGENTA')
        IP = IPSNEW(IE)
        XPNT = XMOD(XP(IP))
        YPNT = YMOD(YP(IP))
        CALL PLSYMB(XPNT,YPNT,SSIZ,1,0.0,0)
        CALL NEWCOLOR(ICOL0)
        CALL PLFLUSH
        CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
 10   CONTINUE
C
      IF(IPSNEW(1).EQ.IPSNEW(2)) THEN
       WRITE(*,*) '**  Endpoints must be distinct  **'
       WRITE(*,*) '**  NEW SEGMENT NOT MARKED OFF  **'
       RETURN
      ENDIF
C
      IPSEG1(ISEG) = MIN(IPSNEW(1),IPSNEW(2))
      IPSEG2(ISEG) = MAX(IPSNEW(1),IPSNEW(2))
      IELSEG(ISEG) = IELMIN
C
      RETURN
      END ! IGSGET
