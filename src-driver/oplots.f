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
C     Version 070-ES1
C     Philip Carter, Esotec Developments, February 2009
C     philip (at) esotec (dot) org
C
C     Changes from 0.70:
C
C     CPINI, QPINI, XYPINI:
C     VERSION argument changed to SVERSION (string).
C
C     CPAXES, QPAXES, XYAXES:
C     CHARACTER*(*) SVERSION statement added
C     PLCHAR(SVERSION) replaces PLNUMB(VERSION)
C
C     Version 070-ES3  September 2013
C     Blades plotted to rotor points (GETROTG, 4 places)
C  
C======================================================================== 
      SUBROUTINE CPPLOT(LPLTEL,ICPXY)
C-----------------------------------------
C     Plots Cp vs x or y, integrated forces, 
C     parameters, and reference data.
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      CHARACTER*80 FNAME, LINE
      PARAMETER (NINPX=4)
      REAL AINP(NINPX)
      LOGICAL ERROR
      PARAMETER (NREFX=500)
      REAL XREF(NREFX),YREF(NREFX)
C
      LOGICAL LPLTEL(NEX)
C
      CALL CPPINI(LPLTEL,ICPXY)
      CALL CPLINE(LPLTEL,ICPXY,.FALSE.)
C
C---- plot various coefficients
      IF(ICPXY.EQ.1) THEN
       XPLT = ( XPLMIN + 0.65*(XPLMAX-XPLMIN))*XPLFAC
       YPLT = -CPLMIN*CPLFAC - 0.4*CHGT
      ELSE
       XPLT = (-CPLMIN - 0.4*(CPLMAX-CPLMIN))*CPLFAC
       YPLT =  YPLMAX*YPLFAC + 15.0*CHGT
      ENDIF
C
      IF(LPREF) THEN
       IF(CPFILE(1:1).NE.' ') THEN
C------ offer existing default
        WRITE(*,1100) CPFILE
 1100   FORMAT(/' Enter filename:  ', A)
        READ(*,1000) FNAME
 1000   FORMAT(A)
        CALL STRIP(FNAME,NFN)
        IF(NFN.EQ.0) FNAME = CPFILE
       ELSE
C------ just ask for filename
        CALL ASKS('Enter filename^',FNAME)
       ENDIF
C
       LU = 9
       OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=20)
       K = 0
       DO ILINE = 1, 12345
 15      READ(LU,1000,END=20) LINE
         IF(INDEX('#!',LINE(1:1)) .NE. 0) GO TO 15
C
         NINP = NINPX
         CALL GETFLT(LINE,AINP,NINP,ERROR)
         IF(ERROR) GO TO 15
         IF(NINP.GE.2) THEN
          K = K+1
          XREF(K) = AINP(1)
          YREF(K) = AINP(2)
          IF(K.GE.NREFX) GO TO 20
         ELSE
          GO TO 15
         ENDIF
       ENDDO
 20    CONTINUE
       CLOSE(LU)
       NREF = K
C
       CALL NEWPEN(2)
       SHL = 0.006
       ISYML = 5
       IF(ICPXY.EQ.1) THEN
        CALL XYSYMB(NREF,XREF,YREF,0.0,XPLFAC,0.0,-CPLFAC,SHL,ISYML)
       ELSE
        CALL XYSYMB(NREF,YREF,XREF,0.0,-CPLFAC,0.0,YPLFAC,SHL,ISYML)
       ENDIF
C
C----- set new default filename
       CPFILE = FNAME
       GO TO 50
C
 40    WRITE(*,*) 'File OPEN error'
 50    CONTINUE
      ENDIF
C
      RPM = 30.0*OMEGA(1)/PI
      CALL COEFPL(XPLT,YPLT,CHGT,
     &            NAME,NNAME,
     &            QREF,MACH,RHO,
     &            QINF,RPM,
     &            TTOT, TDUCT,
     &            NEL,NETYPE,ICOLEL,GAMSET,SIGSET,
     &            CX(0),CY(0),CD(0),CM(0),FOM)
C
      CALL PLFLUSH
      IPTYPE = JCPLT
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! CPPLOT
 

      SUBROUTINE CPPINI(LPLTEL,ICPXY)
C-----------------------------------------
C     Sets up for Cp(x) or Cp(y) plot
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      LOGICAL LPLTEL(NEX)
C
C---- additional left and top margins for plot
      DATA XADD, YADD / 0.05 , 0.05 /
C
C---- character size for axes annotation
      CHA = 0.85*CHGT
C
C---- additional left and bottom margins for axis annotations
      XANN = 3.0*CHA
      YANN = 2.0*CHA
C
C---- determine airfoil box size and location (ignore wakes)
      XMIN =  1.0E30
      XMAX = -1.0E30
      YMIN =  1.0E30
      YMAX = -1.0E30
      DO IEL = 1, NEL
        IF(LPLTEL(IEL)) THEN
         IF(NETYPE(IEL).EQ.0 .OR.
     &      NETYPE(IEL).EQ.1 .OR.
     &      NETYPE(IEL).EQ.5 .OR.
     &      NETYPE(IEL).EQ.6 .OR.
     &      NETYPE(IEL).EQ.7      ) THEN
          I = IPFRST(IEL)
          N = IPLAST(IEL) - IPFRST(IEL) + 1
          CALL MINMAX(N,XP(I),XEMIN,XEMAX)
          CALL MINMAX(N,YP(I),YEMIN,YEMAX)
          XMIN = MIN(XMIN,XEMIN)
          XMAX = MAX(XMAX,XEMAX)
          YMIN = MIN(YMIN,YEMIN)
          YMAX = MAX(YMAX,YEMAX)
         ENDIF
        ENDIF
      ENDDO
C
C---- if valid geometry box defined, use those X,Y limits for plots
      IF(XGBOX(1).NE.XGBOX(2) .AND. YGBOX(1).NE.YGBOX(2)) THEN
        XMIN = XGBOX(1)
        XMAX = XGBOX(2)
        YMIN = YGBOX(1)
        YMAX = YGBOX(2)
      ENDIF
C
C---- fudge limits slightly so airfoil isn't up against axes
      FUDGE = MAX( 0.02*(XMAX-XMIN) , 0.02*(YMAX-YMIN) )
      XPLMIN = XMIN - FUDGE
      XPLMAX = XMAX + FUDGE
      YPLMIN = YMIN - FUDGE
      YPLMAX = YMAX + FUDGE
C
C---- axisymmetric case, go to axis for reference
      YPLMIN = 0.
C
C---- set rounded limits and annotation deltas
      CALL AXISADJ2(XPLMIN,XPLMAX,XSPAN,XPLDEL,NTICS)
      CALL AXISADJ2(YPLMIN,YPLMAX,YSPAN,YPLDEL,NTICS)
C
C---- make sure at least one annotation interval is present
      XPLMAX = MAX(XPLMAX,XPLMIN+XPLDEL)
      YPLMAX = MAX(YPLMAX,YPLMIN+YPLDEL)
C
C---- start plot
      CALL PLTINI2
      CALL GETCOLOR(ICOL0)
C
C---- plotting box aspect ratio
      XBOX = MAX( XPAGE - 2.0*XMARG , 0.0001 )
      YBOX = MAX( YPAGE - 2.0*YMARG , 0.0001 )
      BOXPAR = YBOX/XBOX
C
C---- x,y geometry scaling factors
      IF(ICPXY.EQ.1) THEN
       XPLFAC =            (1.0-2.0*XADD-XANN)/(XPLMAX-XPLMIN)
       YPLFAC =     (0.5*BOXPAR-2.0*YADD-YANN)/(YPLMAX-YPLMIN)
      ELSE
       XPLFAC =            (0.85-2.0*XADD-XANN)/(XPLMAX-XPLMIN)
       YPLFAC =     (0.5*BOXPAR-2.0*YADD-YANN)/(YPLMAX-YPLMIN)
      ENDIF
C
c      XPLFAC = MAX(XPLFAC,YPLFAC)
      YPLFAC = XPLFAC
C
      CPXAR = 0.5*BOXPAR
      CPYAR = 0.25
C
C---- Cp scaling factor
      IF(ICPXY.EQ.1) THEN
       CPLFAC =  CPXAR/(CPLMAX-CPLMIN)
      ELSE
       CPLFAC = -CPYAR/(CPLMAX-CPLMIN)
      ENDIF
C
C---- set origin, airfoil scales, offsets
      IF(ICPXY.EQ.1) THEN
C----- Cp vs x  plot
       XORG =          XADD - XPLMIN*XPLFAC + XANN
       YORG = BOXPAR - YADD + CPLMIN*CPLFAC
C
       FACA =  XPLFAC
       XOFA = 0.
       YOFA = -YMIN + (-YORG + YADD)/FACA
C
      ELSE
C----- Cp vs y  plot
c       XORG = 1.0    - XADD + CPLMIN*CPLFAC
       XORG = XADD +  CPLMIN*CPLFAC
       YORG = YADD -  YPLMIN*XPLFAC + YANN
C
       FACA =  XPLFAC
       XOFA = -XMIN + (-CPLMAX*CPLFAC+XADD)/FACA
c       XOFA = -XMIN + (-XORG + XADD)/FACA
       YOFA = 0.
C
      ENDIF
C
      CALL PLOT( XORG , YORG , -3)
C
C---- plot Cp(x) or Cp(y) axes
      CALL CPAXES(ICPXY,LPGRID, CHA,
     &            XPLMIN,XPLMAX,XPLDEL,XPLFAC,
     &            YPLMIN,YPLMAX,YPLDEL,YPLFAC,
     &            CPLMIN,CPLMAX,CPLDEL,CPLFAC,
     &            XOFA,YOFA,FACA,
     &            'DFDC',SVERSION)
C
C---- plot wake grid
      CALL PLTGRID(SSIZEL,
     &             XG,YG,IX,II,JJ,
     &             XOFA,YOFA,FACA)
C
C---- plot blade or disk outlines
      CALL GETCOLOR(ICOL0)
      CALL NEWPEN(3)
      DO N = 1, NROTOR  
       IF(IRTYPE(N).EQ.1) THEN
        IEL = IELROTOR(N)
        CALL NEWCOLOR(ICOLEL(IEL))
        CALL PLTRBLD(SSIZEL,XPAXIS,
     &               XRC(1,N),YRC(1,N),ZER,ZER,NRC,
     &               XOFA,YOFA,FACA)
       ELSEIF(IRTYPE(N).EQ.2) THEN
        IEL = IELROTOR(N)
        CALL NEWCOLOR(ICOLEL(IEL))
C
C---- get blade geometry at rotor points (esotec)
        CALL GETROTG(N)
        CALL PLTRBLD(SSIZEL,XPAXIS,
     &               XRP(1,N),YRP(1,N),CHRP(1,N),BETARP(1,N),NRP,
     &               XOFA,YOFA,FACA)
C
C        CALL PLTRBLD(SSIZEL,XPAXIS,
C     &               XRC(1,N),YRC(1,N),CHR(1,N),BETAR(1,N),NRC,
C     &               XOFA,YOFA,FACA)
       ENDIF
      END DO
      CALL NEWCOLOR(ICOL0)
C
      CALL PLTAIR(LVPLT,LEPLT,
     &            NEL,IPFRST,IPLAST,XP,YP, LTPAN,
     &            NETYPE,SSIZEL,LPLTEL,ICOLEL, 
     &            XELNUM,YELNUM,
     &            XOFA,YOFA,FACA )
C
C---- plot stagnation point
      IF(LPSTAG) THEN
        IEL = 2
        XPLT = (XSTG(IEL)+XOFA)*FACA
        YPLT = (YSTG(IEL)+YOFA)*FACA
        SSIZ = .008
        CALL PLSYMB(XPLT,YPLT,SSIZ,1,0.0,0)
      ENDIF
C
c      IF(LWALL) THEN
c       CALL NEWCOLORNAME('brown')
c       CALL PLTWAL(IGEOM,
c     &             XWALL,YWALL,AWALL,
c     &             XOFA,YOFA,FACA )
c       CALL NEWCOLOR(ICOL0)
c      ENDIF
C
      RETURN
      END ! CPPINI


      SUBROUTINE CPAXES(ICPXY,LGRID, CSIZ,
     &                  XPLMIN,XPLMAX,XPLDEL,XPLFAC,
     &                  YPLMIN,YPLMAX,YPLDEL,YPLFAC,
     &                  CPLMIN,CPLMAX,CPLDEL,CPLFAC,
     &                  XOFA,YOFA,FACA,
     &                  CODE,SVERSION)
C---------------------------------------------------------------
C     Plots axes and airfoil for Cp vs x  or  Cp vs y plot
C---------------------------------------------------------------
      LOGICAL LGRID
      CHARACTER*(*) CODE
      CHARACTER*(*) SVERSION
C
      EXTERNAL PLCHAR
      INCLUDE 'MASKS.INC'
C
      CALL GETCOLOR(ICOL0)
      CALL GETPAT  (IPAT0)
C
      IF(ICPXY.EQ.1) THEN
       X0 =  XPLMIN*XPLFAC
       X1 =  XPLMAX*XPLFAC
       XD =  XPLDEL*XPLFAC
       XA =  XPLMIN*XPLFAC
C
       Y0 = -CPLMAX*CPLFAC
       Y1 = -CPLMIN*CPLFAC
       YD = -CPLDEL*CPLFAC
       YA = MAX( -CPLMAX*CPLFAC , 0.0 )
C
       XAMIN = XPLMIN
       XADEL = XPLDEL
       YAMIN = CPLMAX
       YADEL = CPLDEL
C
       NDIGX = -3
       NDIGY = -2
C
       IGXFAC = 5
       IGYFAC = 5
      ELSE
       X0 = -CPLMIN*CPLFAC
       X1 = -CPLMAX*CPLFAC
       XD = CPLDEL*CPLFAC
       XA = MAX(-CPLMIN*CPLFAC , 0.0 )
C
       Y0 =  YPLMIN*YPLFAC
       Y1 =  YPLMAX*YPLFAC
       YD =  YPLDEL*YPLFAC
       YA =  YPLMIN*YPLFAC
C
       XAMIN = CPLMIN
       XADEL = -CPLDEL
       YAMIN = YPLMIN
       YADEL = YPLDEL
C
       NDIGX = -2
       NDIGY = -3
C
       IGXFAC = 1
       IGYFAC = 2
      ENDIF
C
      IF(LGRID) THEN
       NXG = INT((X1-X0)/XD + 0.01) * IGXFAC
       NYG = INT((Y1-Y0)/YD + 0.01) * IGYFAC
       DXG = XD / FLOAT(IGXFAC)
       DYG = YD / FLOAT(IGYFAC)
       CALL NEWPEN(1)
       CALL PLGRID(X0,Y0, NXG,DXG, NYG,DYG, LMASK2 )
       CHA = 0.75*CSIZ
       IF(ICPXY.EQ.1) THEN
        CALL XAXIS(X0, Y0, X1-X0, XD, XAMIN,XADEL,CHA,-2)
        CALL YAXIS(XA, Y0, Y1-Y0, YD, YAMIN,YADEL,CHA, NDIGY)
        CALL XAXIS(X0, YA, X1-X0, XD, XAMIN,XADEL,CHA, NDIGX)
       ELSE
        CALL XAXIS(X0, YA, X1-X0, XD, XAMIN,XADEL,CHA, NDIGX)
        CALL YAXIS(0.0,Y0, Y1-Y0, YD, YAMIN,YADEL,CHA,-2)
       ENDIF
      ENDIF
C
C---- plot axes
      CALL NEWPEN(2)
C
      CALL NEWPEN(3)
      CHL = 1.30*CSIZ
      CHI = 0.65*CSIZ
      LENC = LEN(CODE)
C
      IF(ICPXY.EQ.1) THEN
        XLAB = X0 - 3.5*CHL
        YLAB = (AINT(0.7*(Y0+Y1)/YD) + 0.5)*YD - 0.45*CHL
        CALL PLCHAR(XLAB,YLAB,CHL,'C',0.0,1)
        CALL PLSUBS(XLAB,YLAB,CHL,'p',0.0,1,PLCHAR)
        IF(LGRID) THEN
         XLAB = X1 - 0.5*XD - 0.5*CHL
         YLAB = Y0 - 1.9*CHL
         CALL PLCHAR(XLAB,YLAB,CHL,'x',0.0,1)
        ENDIF
        XLI = X0 + 1.5*CHI
        YLI = Y1 - 0.3*CHI
C
      ELSE
        XLAB = (AINT(0.7*(X0+X1)/XD) + 0.5)*XD - 0.9*CHL
        YLAB = Y0 - 3.1*CHL
        CALL PLCHAR(XLAB,YLAB,CHL,'C',0.0,1)
        CALL PLSUBS(XLAB,YLAB,CHL,'p',0.0,1,PLCHAR)
        IF(LGRID) THEN
         XLAB = 0.0 - 3.5*CHL
         YLAB = Y1 - 0.5*YD - 0.5*CHL
         CALL PLCHAR(XLAB,YLAB,CHL,'y',0.0,1)
        ENDIF
        XLI = X0 + 1.5*CHI
        YLI = Y0 + 4.5*CHI
C
      ENDIF
C
C---- plot code identifier
      CALL NEWPEN(2)
      CALL PLCHAR(XLI        ,YLI-1.5*CHI,CHI,CODE   ,0.0,LENC)
      CALL PLCHAR(XLI        ,YLI-3.5*CHI,CHI,'V'    ,0.0,1)
      CALL PLCHAR(XLI+1.5*CHI,YLI-3.5*CHI,CHI,SVERSION,0.0,9)
C
C----- axis for axisymmetric cases
      CALL NEWPEN(1)
      CALL NEWPAT(LMASK1)
      CALL PLOT((XPLMIN+XOFA)*FACA-0.02,YOFA*FACA,3)
      CALL PLOT((XPLMAX+XOFA)*FACA+0.02,YOFA*FACA,2)
      CALL NEWPAT(IPAT0)
C
      CALL NEWCOLOR(ICOL0)
C
      RETURN
      END ! CPAXES


      SUBROUTINE CPLINE(LPLTEL,ICPXY,LPLINV)
C-----------------------------------------
C     Plots Cp vs x curves
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      LOGICAL LPLTEL(NEX)
      LOGICAL LPLINV
C
C---- size and type of line-plot symbol
      SHL = 0.003
      ISYML = 3
C
      CALL GETCOLOR(ICOL0)
C
      DO 10 IEL=1, NEL 
        IF(.NOT.LPLTEL(IEL)) GO TO 10
        IF(ICFRST(IEL).LE.0) GO TO 10
C
        CALL NEWCOLOR(ICOLEL(IEL))
C
        I = ICFRST(IEL)
        N = ICLAST(IEL) - ICFRST(IEL) + 1
C
        IF(ISPLOT(IEL).EQ. 0 .OR.
     &     ISPLOT(IEL).EQ. 1      ) THEN
C-------- plot Cp on right side
          IF(ICPXY.EQ.1) THEN
           CALL NEWPEN(2)
           CALL XYLINE(N,XC(I),CPR(I),0.0,XPLFAC,0.0,-CPLFAC,1)
           IF(LSPLT) THEN
            CALL XYSYMB(N,XC(I),CPR(I),0.0,XPLFAC,0.0,-CPLFAC,SHL,ISYML)
           ENDIF
C
          ELSE
           CALL NEWPEN(2)
c           CALL XYLINE(N,CPR(I),YC(I),CPLMAX,-CPLFAC,0.0,YPLFAC,1)
           CALL XYLINE(N,CPR(I),YC(I),0.0,-CPLFAC,0.0,YPLFAC,1)
           IF(LSPLT) THEN
            CALL XYSYMB(N,CPR(I),YC(I),0.0,-CPLFAC,0.0,YPLFAC,SHL,ISYML)
           ENDIF
C
          ENDIF
        ENDIF

        IF(ISPLOT(IEL).EQ. 0 .OR.
     &     ISPLOT(IEL).EQ. 2      ) THEN
C-------- plot Cp on left side
          IF(ICPXY.EQ.1) THEN
           CALL NEWPEN(2)
           CALL XYLINE(N,XC(I),CPL(I),0.0,XPLFAC,0.0,-CPLFAC,1)
           IF(LSPLT) THEN
            CALL XYSYMB(N,XC(I),CPL(I),0.0,XPLFAC,0.0,-CPLFAC,SHL,ISYML)
           ENDIF
C
          ELSE
           CALL NEWPEN(2)
c           CALL XYLINE(N,CPL(I),YC(I),CPLMAX,-CPLFAC,0.0,YPLFAC,1)
           CALL XYLINE(N,CPL(I),YC(I),0.0,-CPLFAC,0.0,YPLFAC,1)
           IF(LSPLT) THEN
            CALL XYSYMB(N,CPL(I),YC(I),0.0,-CPLFAC,0.0,YPLFAC,SHL,ISYML)
           ENDIF
C
          ENDIF
C
        ENDIF
C
 10   CONTINUE
C
      CALL NEWCOLOR(ICOL0)
      RETURN
      END ! CPLINE



      SUBROUTINE QPLOT(LPLTEL,IQPXY)
C-----------------------------------------
C     Plots Q vs x or y, integrated forces, 
C     parameters, and reference data.
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      CHARACTER*80 FNAME, LINE
      PARAMETER (NINPX=4)
      REAL AINP(NINPX)
      LOGICAL ERROR
      PARAMETER (NREFX=500)
      REAL XREF(NREFX),YREF(NREFX)
C
      LOGICAL LPLTEL(NEX)
C
      CALL QPINI(LPLTEL,IQPXY)
      CALL QPLINE(LPLTEL,IQPXY,.FALSE.)
C
C---- plot various coefficients
      IF(IQPXY.EQ.1) THEN
       XPLT = ( XPLMIN + 0.65*(XPLMAX-XPLMIN))*XPLFAC
       YPLT = QPLMAX*QPLFAC - 0.4*CHGT
      ELSE
       XPLT = (QPLMIN + 0.3*(QPLMAX-QPLMIN))*QPLFAC
       YPLT =  YPLMAX*YPLFAC + 15.0*CHGT
      ENDIF
C
      RPM = 30.0*OMEGA(1)/PI
      CALL COEFPL(XPLT,YPLT,CHGT,
     &            NAME,NNAME,
     &            QREF,MACH,RHO,
     &            QINF,RPM,
     &            TTOT, TDUCT,
     &            NEL,NETYPE,ICOLEL,GAMSET,SIGSET,
     &            CX(0),CY(0),CD(0),CM(0),FOM)
C
      CALL PLFLUSH
      IPTYPE = JQPLT
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! QPLOT


      SUBROUTINE QPINI(LPLTEL,IQPXY)
C-----------------------------------------
C     Sets up for Q(x) or Q(y) plot
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      LOGICAL LPLTEL(NEX)
C
C---- additional left and top margins for plot
      DATA XADD, YADD / 0.05 , 0.05 /
C
C---- character size for axes annotation
      CHA = 0.85*CHGT
C
C---- additional left and bottom margins for axis annotations
      XANN = 3.0*CHA
      YANN = 2.0*CHA
C
C---- determine airfoil box size and location (ignore wakes)
      XMIN =  1.0E30
      XMAX = -1.0E30
      YMIN =  1.0E30
      YMAX = -1.0E30
      DO IEL = 1, NEL
        IF(LPLTEL(IEL)) THEN
         IF(NETYPE(IEL).EQ.0 .OR.
     &      NETYPE(IEL).EQ.1 .OR.
     &      NETYPE(IEL).EQ.5 .OR.
     &      NETYPE(IEL).EQ.6 .OR.
     &      NETYPE(IEL).EQ.7      ) THEN
          I = IPFRST(IEL)
          N = IPLAST(IEL) - IPFRST(IEL) + 1
          CALL MINMAX(N,XP(I),XEMIN,XEMAX)
          CALL MINMAX(N,YP(I),YEMIN,YEMAX)
          XMIN = MIN(XMIN,XEMIN)
          XMAX = MAX(XMAX,XEMAX)
          YMIN = MIN(YMIN,YEMIN)
          YMAX = MAX(YMAX,YEMAX)
         ENDIF
        ENDIF
      ENDDO
C
C---- if valid geometry box defined, use those X,Y limits for plots
      IF(XGBOX(1).NE.XGBOX(2) .AND. YGBOX(1).NE.YGBOX(2)) THEN
        XMIN = XGBOX(1)
        XMAX = XGBOX(2)
        YMIN = YGBOX(1)
        YMAX = YGBOX(2)
      ENDIF
C
C---- fudge limits slightly so airfoil isn't up against axes
      FUDGE = MAX( 0.02*(XMAX-XMIN) , 0.02*(YMAX-YMIN) )
      XPLMIN = XMIN - FUDGE
      XPLMAX = XMAX + FUDGE
      YPLMIN = YMIN - FUDGE
      YPLMAX = YMAX + FUDGE
C
C---- axisymmetric case, go to axis for reference
      YPLMIN = 0.
C
C---- set rounded limits and annotation deltas
      CALL AXISADJ2(XPLMIN,XPLMAX,XSPAN,XPLDEL,NTICS)
      CALL AXISADJ2(YPLMIN,YPLMAX,YSPAN,YPLDEL,NTICS)
C
C---- make sure at least one annotation interval is present
      XPLMAX = MAX(XPLMAX,XPLMIN+XPLDEL)
      YPLMAX = MAX(YPLMAX,YPLMIN+YPLDEL)
C
C---- start plot
      CALL PLTINI2
      CALL GETCOLOR(ICOL0)
C
C---- plotting box aspect ratio
      XBOX = MAX( XPAGE - 2.0*XMARG , 0.0001 )
      YBOX = MAX( YPAGE - 2.0*YMARG , 0.0001 )
      BOXPAR = YBOX/XBOX
C
C---- x,y geometry scaling factors
      IF(IQPXY.EQ.1) THEN
       XPLFAC =            (1.0-2.0*XADD-XANN)/(XPLMAX-XPLMIN)
       YPLFAC =     (0.5*BOXPAR-2.0*YADD-YANN)/(YPLMAX-YPLMIN)
      ELSE
       XPLFAC =            (0.85-2.0*XADD-XANN)/(XPLMAX-XPLMIN)
       YPLFAC =     (0.5*BOXPAR-2.0*YADD-YANN)/(YPLMAX-YPLMIN)
      ENDIF
C
c      XPLFAC = MAX(XPLFAC,YPLFAC)
      YPLFAC = XPLFAC
C
      CPXAR = 0.5*BOXPAR
      CPYAR = 0.25
C
C---- Q scaling factor
      IF(IQPXY.EQ.1) THEN
       QPLFAC = CPXAR/(QPLMAX-QPLMIN)
      ELSE
       QPLFAC = CPYAR/(QPLMAX-QPLMIN)
      ENDIF
C
C---- set origin, airfoil scales, offsets
      IF(IQPXY.EQ.1) THEN
C----- Q vs x  plot
       XORG =          XADD - XPLMIN*XPLFAC + XANN
       YORG = BOXPAR - YADD - (QPLMAX-QPLMIN)*QPLFAC
       FACA =  XPLFAC
       XOFA = 0.
       YOFA = -YMIN + (-YORG + YADD)/FACA
C
      ELSE
C----- Q vs y  plot
       XORG = XADD                  + XANN
       YORG = YADD -  YPLMIN*XPLFAC + YANN
C
       FACA =  XPLFAC
       XOFA = -XMIN + (QPLMAX*QPLFAC+XADD)/FACA
c       XOFA = -XMIN + (-XORG + XADD)/FACA
       YOFA = 0.
C
      ENDIF
C
      CALL PLOT( XORG , YORG , -3)
C
C---- plot Q(x) or Q(y) axes
      CALL QPAXES(IQPXY,LPGRID, CHA,
     &            XPLMIN,XPLMAX,XPLDEL,XPLFAC,
     &            YPLMIN,YPLMAX,YPLDEL,YPLFAC,
     &            QPLMIN,QPLMAX,QPLDEL,QPLFAC,
     &            XOFA,YOFA,FACA,
     &            'DFDC',SVERSION)
C
C---- plot wake grid
      CALL PLTGRID(SSIZEL,
     &             XG,YG,IX,II,JJ,
     &             XOFA,YOFA,FACA)
C
C---- plot blade or disk outlines
      CALL GETCOLOR(ICOL0)
      CALL NEWPEN(3)
      DO N = 1, NROTOR  
       IF(IRTYPE(N).EQ.1) THEN
        IEL = IELROTOR(N)
        CALL NEWCOLOR(ICOLEL(IEL))
        CALL PLTRBLD(SSIZEL,XPAXIS,
     &               XRC(1,N),YRC(1,N),ZER,ZER,NRC,
     &               XOFA,YOFA,FACA)
       ELSEIF(IRTYPE(N).EQ.2) THEN
        IEL = IELROTOR(N)
        CALL NEWCOLOR(ICOLEL(IEL))
C
C---- get blade geometry at rotor points (esotec)
        CALL GETROTG(N)
        CALL PLTRBLD(SSIZEL,XPAXIS,
     &               XRP(1,N),YRP(1,N),CHRP(1,N),BETARP(1,N),NRP,
     &               XOFA,YOFA,FACA)
C
C        CALL PLTRBLD(SSIZEL,XPAXIS,
C     &               XRC(1,N),YRC(1,N),CHR(1,N),BETAR(1,N),NRC,
C     &               XOFA,YOFA,FACA)
       ENDIF
      END DO
      CALL NEWCOLOR(ICOL0)
C
      CALL PLTAIR(LVPLT,LEPLT,
     &            NEL,IPFRST,IPLAST,XP,YP, LTPAN,
     &            NETYPE,SSIZEL,LPLTEL,ICOLEL, 
     &            XELNUM,YELNUM,
     &            XOFA,YOFA,FACA )
C
C
C---- plot stagnation point
      IF(LPSTAG) THEN
        IEL = 2
        XPLT = (XSTG(IEL)+XOFA)*FACA
        YPLT = (YSTG(IEL)+YOFA)*FACA
        SSIZ = .008
        CALL PLSYMB(XPLT,YPLT,SSIZ,1,0.0,0)
      ENDIF
C
c      IF(LWALL) THEN
c       CALL NEWCOLORNAME('brown')
c       CALL PLTWAL(IGEOM,
c     &             XWALL,YWALL,AWALL,
c     &             XOFA,YOFA,FACA )
c       CALL NEWCOLOR(ICOL0)
c      ENDIF
C
      RETURN
      END ! QPINI


      SUBROUTINE QPAXES(IQPXY,LGRID, CSIZ,
     &                  XPLMIN,XPLMAX,XPLDEL,XPLFAC,
     &                  YPLMIN,YPLMAX,YPLDEL,YPLFAC,
     &                  QPLMIN,QPLMAX,QPLDEL,QPLFAC,
     &                  XOFA,YOFA,FACA,
     &                  CODE,SVERSION)
C---------------------------------------------------------------
C     Plots axes and airfoil for Q vs x  or  Q vs y plot
C---------------------------------------------------------------
      LOGICAL LGRID
      CHARACTER*(*) CODE
      CHARACTER*(*) SVERSION
C
      EXTERNAL PLCHAR
      INCLUDE 'MASKS.INC'
C
      CALL GETCOLOR(ICOL0)
      CALL GETPAT  (IPAT0)
C
      IF(IQPXY.EQ.1) THEN
       X0 =  XPLMIN*XPLFAC
       X1 =  XPLMAX*XPLFAC
       XD =  XPLDEL*XPLFAC
       XA =  XPLMIN*XPLFAC
C
       Y0 = QPLMIN*QPLFAC
       Y1 = QPLMAX*QPLFAC
       YD = QPLDEL*QPLFAC
       YA = MAX( QPLMIN*QPLFAC , 0.0 )
C
       XAMIN = XPLMIN
       XADEL = XPLDEL
       YAMIN = QPLMIN
       YADEL = QPLDEL
C
       NDIGX = -3
       NDIGY = -2
C
       IGXFAC = 2
       IGYFAC = 4
      ELSE
       X0 = QPLMIN*QPLFAC
       X1 = QPLMAX*QPLFAC
       XD = QPLDEL*QPLFAC
       XA = MAX( QPLMIN*QPLFAC , 0.0 )
C
       Y0 =  YPLMIN*YPLFAC
       Y1 =  YPLMAX*YPLFAC
       YD =  YPLDEL*YPLFAC
       YA =  YPLMIN*YPLFAC
C
       XAMIN = QPLMIN
       XADEL = QPLDEL
       YAMIN = YPLMIN
       YADEL = YPLDEL
C
       NDIGX = -2
       NDIGY = -3
C
       IGXFAC = 2
       IGYFAC = 2
      ENDIF
C
      IF(LGRID) THEN
       NXG = INT((X1-X0)/XD + 0.01) * IGXFAC
       NYG = INT((Y1-Y0)/YD + 0.01) * IGYFAC
       DXG = XD / FLOAT(IGXFAC)
       DYG = YD / FLOAT(IGYFAC)
       CALL NEWPEN(1)
       CALL PLGRID(X0,Y0, NXG,DXG, NYG,DYG, LMASK2 )
       CHA = 0.75*CSIZ
       IF(IQPXY.EQ.1) THEN
        CALL XAXIS(X0,Y0, X1-X0, XD, XAMIN,XADEL,CHA,-2)
       ELSE
        CALL YAXIS(X0,Y0, Y1-Y0, YD, YAMIN,YADEL,CHA,-2)
       ENDIF
      ENDIF
C
C---- plot axes
      CALL NEWPEN(2)
      CALL XAXIS(X0,YA, X1-X0, XD, XAMIN,XADEL,CHA, NDIGX)
      CALL YAXIS(XA,Y0, Y1-Y0, YD, YAMIN,YADEL,CHA, NDIGY)
C
      CALL NEWPEN(3)
c      CHL = 1.30*CSIZ
      CHL = 1.10*CSIZ
      CHI = 0.65*CSIZ
      LENC = LEN(CODE)
C
      IF(IQPXY.EQ.1) THEN
        XLAB = X0 - 6.5*CHL
        YLAB = (AINT(0.7*(Y0+Y1)/YD) + 0.5)*YD - 0.45*CHL
        CALL PLCHAR(XLAB,YLAB,CHL,'Q/Qref',0.0,6)
        IF(LGRID) THEN
         XLAB = X1 - 0.5*XD - 0.5*CHL
         YLAB = Y0 - 1.9*CHL
         CALL PLCHAR(XLAB,YLAB,CHL,'x',0.0,1)
        ENDIF
        XLI = X0 + 1.5*CHI
        YLI = Y1 - 0.3*CHI
C
      ELSE
        XLAB = (AINT(0.7*(X0+X1)/XD) + 0.5)*XD - 3.0*CHL
        YLAB = Y0 - 3.5*CHL
        CALL PLCHAR(XLAB,YLAB,CHL,'Q/Qref',0.0,6)
        IF(LGRID) THEN
         XLAB = X0 - 3.5*CHL
         YLAB = Y1 - 0.5*YD - 0.5*CHL
         CALL PLCHAR(XLAB,YLAB,CHL,'y',0.0,1)
        ENDIF
C        XLI = X1 -     CHI*FLOAT(MAX(LENC,6))
C        YLI = Y1 - 0.3*CHI
C
        XLI = X1 - 9.0*CHI
        YLI = Y0 + 4.5*CHI
C
      ENDIF
C
C---- plot code identifier
      CALL NEWPEN(2)
      CALL PLCHAR(XLI        ,YLI-1.5*CHI,CHI,CODE   ,0.0,LENC)
      CALL PLCHAR(XLI        ,YLI-3.5*CHI,CHI,'V'    ,0.0,1)
      CALL PLCHAR(XLI+1.5*CHI,YLI-3.5*CHI,CHI,SVERSION,0.0,9)
C
C----- axis for axisymmetric cases
      CALL NEWPEN(1)
      CALL NEWPAT(LMASK1)
      CALL PLOT((XPLMIN+XOFA)*FACA-0.02,YOFA*FACA,3)
      CALL PLOT((XPLMAX+XOFA)*FACA+0.02,YOFA*FACA,2)
      CALL NEWPAT(IPAT0)
C
      CALL NEWCOLOR(ICOL0)
C
      RETURN
      END ! QPAXES


      SUBROUTINE QPLINE(LPLTEL,IQPXY,LPLINV)
C-----------------------------------------
C     Plots Q vs x curves
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      DIMENSION QM(IPX)
C
      LOGICAL LPLTEL(NEX)
      LOGICAL LPLINV
C
C---- size and type of line-plot symbol
      SHL = 0.003
      ISYML = 3
C
      CALL GETCOLOR(ICOL0)
C
      DO 10 IEL=1, NEL 
        IF(.NOT.LPLTEL(IEL)) GO TO 10
        IF(ICFRST(IEL).LE.0) GO TO 10
C
        CALL NEWCOLOR(ICOLEL(IEL))
C
        I = ICFRST(IEL)
        N = ICLAST(IEL) - ICFRST(IEL) + 1
C
        IF(ISPLOT(IEL).EQ. 0 .OR.
     &     ISPLOT(IEL).EQ. 1      ) THEN
C-------- get velocity magnitudes QM
         IC1 = ICFRST(IEL) 
         IC2 = ICLAST(IEL) 
         DO IC = IC1, IC2 
            QM(IC) = SQRT(QCR(1,IC)**2 + QCR(2,IC)**2)/QREF
          END DO
C-------- plot Q on right side
          IF(IQPXY.EQ.1) THEN
           CALL NEWPEN(2)
           CALL XYLINE(N,XC(I),QM(I),0.0,XPLFAC,0.0,QPLFAC,1)
           IF(LSPLT) THEN
            CALL XYSYMB(N,XC(I),QM(I),0.0,XPLFAC,0.0,QPLFAC,SHL,ISYML)
           ENDIF
C
          ELSE
           CALL NEWPEN(2)
           CALL XYLINE(N,QM(I),YC(I),0.0,QPLFAC,0.0,YPLFAC,1)
           IF(LSPLT) THEN
            CALL XYSYMB(N,QM(I),YC(I),0.0,QPLFAC,0.0,YPLFAC,SHL,ISYML)
           ENDIF
C
          ENDIF
        ENDIF

        IF(ISPLOT(IEL).EQ. 0 .OR.
     &     ISPLOT(IEL).EQ. 2      ) THEN
C-------- get velocity magnitudes QM
         IC1 = ICFRST(IEL) 
         IC2 = ICLAST(IEL) 
          DO IC = IC1, IC2 
            QM(IC) = SQRT(QCL(1,IC)**2 + QCL(2,IC)**2)/QREF
          END DO
C-------- plot Cp on left side
          IF(IQPXY.EQ.1) THEN
           CALL NEWPEN(2)
           CALL XYLINE(N,XC(I),QM(I),0.0,XPLFAC,0.0,QPLFAC,1)
           IF(LSPLT) THEN
            CALL XYSYMB(N,XC(I),QM(I),0.0,XPLFAC,0.0,QPLFAC,SHL,ISYML)
           ENDIF
C
          ELSE
           CALL NEWPEN(2)
           CALL XYLINE(N,QM(I),YC(I),0.0,QPLFAC,0.0,YPLFAC,1)
           IF(LSPLT) THEN
            CALL XYSYMB(N,QM(I),YC(I),0.0,QPLFAC,0.0,YPLFAC,SHL,ISYML)
           ENDIF
C
          ENDIF
C
        ENDIF
C
 10   CONTINUE
C
      CALL NEWCOLOR(ICOL0)
      RETURN
      END ! QPLINE



      SUBROUTINE XYPLOT(LPLTEL      )
C-----------------------------------------
C     Plots x,y geometry with reference grid
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      LOGICAL LPLTEL(NEX)
      CALL XYPINI(LPLTEL)
      CALL PLFLUSH
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
cc      IPTYPE = XYPLT
      RETURN
      END

     

      SUBROUTINE XYPINI(LPLTEL      )
C-----------------------------------------
C     Sets up for geometry plot
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      LOGICAL LPLTEL(NEX)
C
C---- additional left and top margins for plot
      DATA XADD, YADD / 0.05 , 0.05 /
C
C---- character size for axes annotation
      CHA = 0.85*CHGT
C
C---- additional left and bottom margins for axis annotations
      XANN = 3.0*CHA
      YANN = 2.0*CHA
C
C---- determine airfoil box size and location (ignore wakes)
      XMIN =  1.0E30
      XMAX = -1.0E30
      YMIN =  1.0E30
      YMAX = -1.0E30
      DO IEL = 1, NEL
        IF(LPLTEL(IEL)) THEN
         IF(NETYPE(IEL).EQ.0 .OR.
     &      NETYPE(IEL).EQ.1 .OR.
     &      NETYPE(IEL).EQ.5 .OR.
     &      NETYPE(IEL).EQ.6 .OR.
     &      NETYPE(IEL).EQ.7      ) THEN
          I = IPFRST(IEL)
          N = IPLAST(IEL) - IPFRST(IEL) + 1
          CALL MINMAX(N,XP(I),XEMIN,XEMAX)
          CALL MINMAX(N,YP(I),YEMIN,YEMAX)
          XMIN = MIN(XMIN,XEMIN)
          XMAX = MAX(XMAX,XEMAX)
          YMIN = MIN(YMIN,YEMIN)
          YMAX = MAX(YMAX,YEMAX)
         ENDIF
        ENDIF
      ENDDO
C
C---- if valid geometry box defined, use those X,Y limits for plots
      IF(XGBOX(1).NE.XGBOX(2) .AND. YGBOX(1).NE.YGBOX(2)) THEN
        XMIN = XGBOX(1)
        XMAX = XGBOX(2)
        YMIN = YGBOX(1)
        YMAX = YGBOX(2)
      ENDIF
C
C---- fudge limits slightly so airfoil isn't up against axes
      FUDGE = MAX( 0.02*(XMAX-XMIN) , 0.02*(YMAX-YMIN) )
      XPLMIN = XMIN - FUDGE
      XPLMAX = XMAX + FUDGE
      YPLMIN = YMIN - FUDGE
      YPLMAX = YMAX + FUDGE
C
C---- axisymmetric case, go to axis for reference
      YPLMIN = 0.
C
C---- set rounded limits and annotation deltas
      CALL AXISADJ2(XPLMIN,XPLMAX,XSPAN,XPLDEL,NTICS)
      CALL AXISADJ2(YPLMIN,YPLMAX,YSPAN,YPLDEL,NTICS)
C
C---- make sure at least one annotation interval is present
      XPLMAX = MAX(XPLMAX,XPLMIN+XPLDEL)
      YPLMAX = MAX(YPLMAX,YPLMIN+YPLDEL)
C
C---- start plot
      CALL PLTINI2
      CALL GETCOLOR(ICOL0)
C
C---- plotting box aspect ratio
      XBOX = MAX( XPAGE - 2.0*XMARG , 0.0001 )
      YBOX = MAX( YPAGE - 2.0*YMARG , 0.0001 )
      BOXPAR = YBOX/XBOX
C
C---- x,y geometry scaling factors
      XPLFAC =    (1.0-2.0*XADD-XANN)/(XPLMAX-XPLMIN)
      YPLFAC = (BOXPAR-2.0*YADD-YANN)/(YPLMAX-YPLMIN)
c      XYPLFAC = MAX(XPLFAC,YPLFAC)
C
C---- set origin, airfoil scales, offsets
C----- geometry plot
      XORG = XADD - XPLMIN*XPLFAC + XANN
      YORG = YADD + 0.4*BOXPAR - 0.5*(YPLMAX-YPLMIN)*XPLFAC + YANN
C
      FACA =  XPLFAC
      XOFA =  0.
      YOFA = -YMIN + (YADD)/FACA
C
      CALL PLOT( XORG , YORG , -3)
C
C---- plot x,y axes
      CALL XYAXES(LPGRID, CHA,
     &            XPLMIN,XPLMAX,XPLDEL,XPLFAC,
     &            YPLMIN,YPLMAX,YPLDEL,YPLFAC,
     &            XOFA,YOFA,FACA,
     &            'DFDC',SVERSION)
C
C---- plot wake grid
      CALL PLTGRID(SSIZEL,
     &             XG,YG,IX,II,JJ,
     &             XOFA,YOFA,FACA)
C
C---- plot blade or disk outlines
      CALL GETCOLOR(ICOL0)
      CALL NEWPEN(3)
      DO N = 1, NROTOR  
       IF(IRTYPE(N).EQ.1) THEN
        IEL = IELROTOR(N)
        CALL NEWCOLOR(ICOLEL(IEL))
        CALL PLTRBLD(SSIZEL,XPAXIS,
     &               XRC(1,N),YRC(1,N),ZER,ZER,NRC,
     &               XOFA,YOFA,FACA)
       ELSEIF(IRTYPE(N).EQ.2) THEN
        IEL = IELROTOR(N)
        CALL NEWCOLOR(ICOLEL(IEL))
C
C---- get blade geometry at rotor points (esotec)
        CALL GETROTG(N)
        CALL PLTRBLD(SSIZEL,XPAXIS,
     &               XRP(1,N),YRP(1,N),CHRP(1,N),BETARP(1,N),NRP,
     &               XOFA,YOFA,FACA)
C
C        CALL PLTRBLD(SSIZEL,XPAXIS,
C     &               XRC(1,N),YRC(1,N),CHR(1,N),BETAR(1,N),NRC,
C     &               XOFA,YOFA,FACA)
       ENDIF
      END DO
      CALL NEWCOLOR(ICOL0)
C
      CALL PLTAIR(LVPLT,LEPLT,
     &            NEL,IPFRST,IPLAST,XP,YP, LTPAN,
     &            NETYPE,SSIZEL,LPLTEL,ICOLEL, 
     &            XELNUM,YELNUM,
     &            XOFA,YOFA,FACA )
C
      RETURN
      END ! XYPINI


      SUBROUTINE XYAXES(LGRID, CSIZ,
     &                  XPLMIN,XPLMAX,XPLDEL,XPLFAC,
     &                  YPLMIN,YPLMAX,YPLDEL,YPLFAC,
     &                  XOFA,YOFA,FACA,
     &                  CODE,SVERSION)
C---------------------------------------------------------------
C     Plots axes and airfoil for Cp vs x  or  Cp vs y plot
C---------------------------------------------------------------
      LOGICAL LGRID
      CHARACTER*(*) CODE
      CHARACTER*(*) SVERSION
C
      EXTERNAL PLCHAR
      INCLUDE 'MASKS.INC'
C
      CALL GETCOLOR(ICOL0)
      CALL GETPAT  (IPAT0)
C
      X0 =  (XPLMIN+XOFA)*FACA
      X1 =  (XPLMAX+XOFA)*FACA
      XD =   XPLDEL      *FACA
      XA =  (XPLMIN+XOFA)*FACA
C     
      Y0 =  (YPLMIN+YOFA)*FACA
      Y1 =  (YPLMAX+YOFA)*FACA
      YD =   YPLDEL      *FACA
      YA =  (YPLMIN+YOFA)*FACA
C     
      XAMIN = XPLMIN
      XADEL = XPLDEL
      YAMIN = YPLMIN
      YADEL = YPLDEL
C     
      NDIGX = -3
      NDIGY = -2
C     
      IGXFAC = 1
      IGYFAC = 1
C
      IF(LGRID) THEN
       NXG = INT((X1-X0)/XD + 0.01) * IGXFAC
       NYG = INT((Y1-Y0)/YD + 0.01) * IGYFAC
       DXG = XD / FLOAT(IGXFAC)
       DYG = YD / FLOAT(IGYFAC)
       CALL NEWPEN(1)
       CALL PLGRID(X0,Y0, NXG,DXG, NYG,DYG, LMASK2 )
       CHA = 0.75*CSIZ
c        CALL XAXIS(X0,Y0, X1-X0, XD, XAMIN,XADEL,CHA,-2)
c        CALL YAXIS(X0,Y0, Y1-Y0, YD, YAMIN,YADEL,CHA,-2)
      ENDIF
C
C---- plot axes
      CALL NEWPEN(2)
      CALL XAXIS(X0,YA, X1-X0, XD, XAMIN,XADEL,CHA, NDIGY)
      CALL YAXIS(XA,Y0, Y1-Y0, YD, YAMIN,YADEL,CHA, NDIGY)
C
      CALL NEWPEN(3)
      CHL = 1.30*CSIZ
      CHI = 0.65*CSIZ
      LENC = LEN(CODE)
C     
      IF(LGRID) THEN
       XLAB = X1 - 0.5*XD - 0.5*CHL
       YLAB = Y0 - 1.9*CHL
       CALL PLCHAR(XLAB,YLAB,CHL,'x',0.0,1)
       XLAB = X0 - 3.5*CHL
       YLAB = Y1 - 0.5*YD - 0.5*CHL
       CALL PLCHAR(XLAB,YLAB,CHL,'y',0.0,1)
      ENDIF
C    
      XLI = X1 -     CHI*FLOAT(MAX(LENC,6))
      YLI = Y1 - 1.0*CHI
C
C---- plot code identifier
      CALL NEWPEN(2)
      CALL PLCHAR(XLI        ,YLI-1.5*CHI,CHI,CODE   ,0.0,LENC)
      CALL PLCHAR(XLI        ,YLI-3.5*CHI,CHI,'V'    ,0.0,1)
      CALL PLCHAR(XLI+1.5*CHI,YLI-3.5*CHI,CHI,SVERSION,0.0,9)
C
C----- axis for axisymmetric cases
      CALL NEWPEN(1)
      CALL NEWPAT(LMASK1)
      CALL PLOT((XPLMIN+XOFA)*FACA-0.02,YOFA*FACA,3)
      CALL PLOT((XPLMAX+XOFA)*FACA+0.02,YOFA*FACA,2)
      CALL NEWPAT(IPAT0)
C
      CALL NEWCOLOR(ICOL0)
C
      RETURN
      END ! XYAXES



      SUBROUTINE PLTAIR(LVPLT,LEPLT,
     &                  NEL,IFRST,ILAST,X,Y, LTPAN,
     &                  NETYPE,SSIZEL,LPLTEL,ICOLEL,
     &                  XELNUM,YELNUM,
     &                  XOFA,YOFA,FACA)
C----------------------------------------------
C     Plots axes and airfoil for Cp vs x plot
C----------------------------------------------
      LOGICAL LVPLT, LEPLT
      DIMENSION IFRST(*), ILAST(*)
      DIMENSION X(*),Y(*)
      LOGICAL LTPAN(*)
C
      DIMENSION NETYPE(*)
      DIMENSION SSIZEL(0:7)
      LOGICAL LPLTEL(*)
      DIMENSION ICOLEL(*)
      DIMENSION XELNUM(*), YELNUM(*)
C
      INCLUDE 'MASKS.INC'
C
      CALL GETCOLOR(ICOL0)
      CALL GETPAT  (IPAT0)
C
C---- plot all elements
      YMAX = 0.
      DO 10 IEL=1, NEL
C------ skip this element if it is to be invisible
        IF(.NOT.LPLTEL(IEL)) GO TO 10
C
        CALL NEWCOLOR(ICOLEL(IEL))
C
        I1 = IFRST(IEL)
        I2 = ILAST(IEL)
C
        I = I1
        N = I2 - I1 + 1
C------- solid-surface element
        IF    (NETYPE(IEL).EQ.0) THEN
         CALL NEWPEN(3)
         CALL XYLINE(N,X(I),Y(I),-XOFA,FACA,-YOFA,FACA,1)
C------- wake element or axis line, source-only or vortex wake
        ELSEIF(NETYPE(IEL).EQ.1
     &    .OR. NETYPE(IEL).EQ.2
     &    .OR. NETYPE(IEL).EQ.5
     &    .OR. NETYPE(IEL).EQ.6
     &    .OR. NETYPE(IEL).EQ.7) THEN
         CALL NEWPEN(1)
         CALL NEWPAT(LMASK1)
         CALL XYLINE(N,X(I),Y(I),-XOFA,FACA,-YOFA,FACA,1)
cc       CALL XYDASH(N,X(I),Y(I),-XOFA,FACA,-YOFA,FACA,0.65)
         CALL NEWPAT(IPAT0)
        ENDIF
C
        IF(LTPAN(IEL)) THEN
C------- plot TE-gap panel
         CALL NEWPEN(1)
         CALL NEWPAT(LMASK1)
C
C------- assume normal TE panel connecting surface endpoints
         XTE1 = X(I1)
         YTE1 = Y(I1)
         XTE2 = X(I2)
         YTE2 = Y(I2)
C
C-------- vortex wake with closure use downstream point
          IF(NETYPE(IEL).EQ.7) THEN
            XTE1 = X(I2)
            YTE1 = Y(I2)
            XTE2 = X(I2)
            YTE2 = 0.0
C-------- axisymmetric... check if TE panel is actually a disk on axis
          ELSEIF(YTE1.EQ.0.0) THEN
C--------- point 1 is on axis... set TE disk from endpoint 2 down to axis
           XTE1 = XTE2
          ELSEIF(YTE2.EQ.0.0) THEN
C--------- point 2 is on axis... set TE disk from endpoint 1 down to axis
           XTE2 = XTE1
          ENDIF
         ENDIF
         CALL PLOT((XTE1+XOFA)*FACA,(YTE1+YOFA)*FACA,3)
         CALL PLOT((XTE2+XOFA)*FACA,(YTE2+YOFA)*FACA,2)
         CALL NEWPAT(IPAT0)
C
C------ plot point singularities, or panel nodes if requested
        IF(N.EQ.1 .OR. LVPLT) THEN
         CALL NEWPEN(2)
         SSIZ = SSIZEL(NETYPE(IEL))
         DO I = IFRST(IEL), ILAST(IEL)
           CALL PLSYMB((X(I)+XOFA)*FACA,(Y(I)+YOFA)*FACA,SSIZ,1,0.0,0)
         ENDDO
        ENDIF
C
C------ plot element numbers
        IF(LEPLT) THEN
         CALL NEWPEN(2)
         CSN = 0.006
         XNUM = (XELNUM(IEL)+XOFA)*FACA + SSIZEL(NETYPE(IEL))
         YNUM = (YELNUM(IEL)+YOFA)*FACA + SSIZEL(NETYPE(IEL))
         REL = FLOAT(IEL)
         CALL PLNUMB(XNUM,YNUM,CSN,REL,0.0,-1)
        ENDIF
C
 10   CONTINUE
C
      CALL NEWCOLOR(ICOL0)
      CALL NEWPAT  (IPAT0)
C
      RETURN
      END ! PLTAIR


      SUBROUTINE PLTRBLD(SSIZEL,XPAX,
     &                   XR,YR,CHR,BETA,NR,
     &                   XOFA,YOFA,FACA)
C----------------------------------------------
C     Plots blade outline on geometry plot 
C----------------------------------------------
      DIMENSION XR(NR),YR(NR),CHR(NR),BETA(NR)
      DIMENSION SSIZEL(0:7)
C
      INCLUDE 'MASKS.INC'
C
      IF(NR.LE.2) THEN
        RETURN
      ENDIF
C
      CALL GETPAT(IPAT0)
C
C------ plot rotor line
      CALL NEWPEN(2)
C---- Plot LE line
      SINB = SIN(BETA(1))
      COSB = COS(BETA(1))
      DX = CHR(1)*SINB
      DY = CHR(1)*COSB
      X = XR(1) - XPAX*DX
      Y = YR(1)
      CALL PLOT((X+XOFA)*FACA,(Y+YOFA)*FACA,3)
      ITYP = 1
      SSIZ = SSIZEL(ITYP)
      DO I = 2, NR
       SINB = SIN(BETA(I))
       COSB = COS(BETA(I))
       DX = CHR(I)*SINB
       DY = CHR(I)*COSB
       X = XR(I) - XPAX*DX
       Y = YR(I)
       CALL PLOT((X+XOFA)*FACA,(Y+YOFA)*FACA,2)
cc       CALL PLSYMB((X+XOFA)*FACA,(Y+YOFA)*FACA,SSIZ,1,0.0,0)
      ENDDO
C---- Plot TE line
      SINB = SIN(BETA(1))
      COSB = COS(BETA(1))
      DX = CHR(1)*SINB
      DY = CHR(1)*COSB
      X = XR(1) + (1.0-XPAX)*DX
      Y = YR(1)
      CALL PLOT((X+XOFA)*FACA,(Y+YOFA)*FACA,3)
      ITYP = 1
      SSIZ = SSIZEL(ITYP)
      DO I = 2, NR
       SINB = SIN(BETA(I))
       COSB = COS(BETA(I))
       DX = CHR(I)*SINB
       DY = CHR(I)*COSB
       X = XR(I) + (1.0-XPAX)*DX
       Y = YR(I)
       CALL PLOT((X+XOFA)*FACA,(Y+YOFA)*FACA,2)
cc       CALL PLSYMB((X+XOFA)*FACA,(Y+YOFA)*FACA,SSIZ,1,0.0,0)
      ENDDO
C
      RETURN
      END ! PLTRBLD


      SUBROUTINE PLTADISK(SSIZEL,
     &                    X1,Y1,X2,Y2,
     &                    XOFA,YOFA,FACA)
C----------------------------------------------
C     Plots actuator disk line
C----------------------------------------------
      DIMENSION SSIZEL(0:7)
C
cc      INCLUDE 'DFDC.INC'
      INCLUDE 'MASKS.INC'
C
      CALL GETPAT(IPAT0)
C------ plot disk line
      CALL NEWPEN(2)
      CALL PLOT((X1+XOFA)*FACA,(Y1+YOFA)*FACA,3)
      CALL PLOT((X2+XOFA)*FACA,(Y2+YOFA)*FACA,2)
C
      RETURN
      END ! PLTADISK






      SUBROUTINE AIRLIM(N,X,Y,XMIN,XMAX,YMIN,YMAX)
      DIMENSION X(N),Y(N)
C-----------------------------------------
C     Sets airfoil width and thickness 
C     for airfoil plot space allocation.
C-----------------------------------------
C
      XMIN = X(1)
      XMAX = X(1)
      YMIN = Y(1)
      YMAX = Y(1)
      DO I=1, N
        XMIN = MIN(XMIN,X(I))
        XMAX = MAX(XMAX,X(I))
        YMIN = MIN(YMIN,Y(I))
        YMAX = MAX(YMAX,Y(I))
      ENDDO
      AIRDX = XMAX - XMIN
      AIRDY = YMAX - YMIN
C
      IF(AIRDX.EQ.0.0 .AND.
     &   AIRDY.EQ.0.0      ) RETURN
C
C---- round up to nearest 10% of max dimension
      AIRDIM = MAX( AIRDX, AIRDY )
      AIRDX = 0.05*AIRDIM * AINT(AIRDX/(0.05*AIRDIM) + 1.2)
      AIRDY = 0.05*AIRDIM * AINT(AIRDY/(0.05*AIRDIM) + 1.2)
C
      XAVG = 0.5*(XMAX+XMIN)
      YAVG = 0.5*(YMAX+YMIN)
C
      XMIN = XAVG - 0.5*AIRDX
      XMAX = XAVG + 0.5*AIRDX
      YMIN = YAVG - 0.5*AIRDY
      YMAX = YAVG + 0.5*AIRDY
C
C---- fudge y-space again to 25% of plot width
      DDY = MIN( AIRDY , 0.25*AIRDX ) - AIRDY
C
C---- fudge y limits to match fudged y space, keeping average y the same
      YMIN = YMIN - 0.5*DDY
      YMAX = YMAX + 0.5*DDY
C
      RETURN
      END ! AIRLIM



      SUBROUTINE PLTWAL(XWALL,YWALL,AWALL,
     &                  XOFA,YOFA,FACA )
      DATA NT / 50 /
C
      SINA = -1.0
      COSA =  0.0
C
      S1 = -1.0/FACA
      S2 =  1.0/FACA
      X1 = XWALL + S1*COSA
      X2 = XWALL + S2*COSA
      Y1 = YWALL + S1*SINA
      Y2 = YWALL + S2*SINA
C
      CALL NEWPEN(5)
      CALL PLOT((X1+XOFA)*FACA,(Y1+YOFA)*FACA,3)
      CALL PLOT((X2+XOFA)*FACA,(Y2+YOFA)*FACA,2)
C
      CALL NEWPEN(2)
      DS = (S2-S1)/FLOAT(NT)
      DTICK = 0.2*DS
      DO IT = 1, NT
        ST = S1 + DS*(FLOAT(IT)-0.5)
        XT = XWALL + ST*COSA
        YT = YWALL + ST*SINA
        CALL PLOT((XT+XOFA)*FACA,(YT+YOFA)*FACA,3)
C
        XTICK = XT - DTICK*COSA + DTICK*SINA
        YTICK = YT - DTICK*COSA - DTICK*SINA
        CALL PLOT((XTICK+XOFA)*FACA,(YTICK+YOFA)*FACA,2)
      ENDDO
C
      RETURN
      END ! PLTWAL




      SUBROUTINE COEFPL(XL0,YL0,CSIZ,
     &            NAME,NNAME,
     &            QREF,AMACH,RHO,
     &            QINF,RPM,
     &            TROTOR, TDUCT,
     &            NEL,NETYPE,ICOL,GAM,SIG,
     &            CX,CY,CD,CM,FOM)
C------------------------------------------------------------------
C     Plots force coefficients for single-point  Cp vs x  plot.
C
C     XL,YL   upper-left corner of label block, 
C             returned as location of lower-left corner
C
C------------------------------------------------------------------
      CHARACTER*(*) NAME
C
      DIMENSION NETYPE(*), ICOL(*)
      DIMENSION GAM(*),    SIG(*)
C
      CHARACTER*1 CINDEX
      EXTERNAL PLCHAR, PLMATH
C
cc      CHN = 1.20*CSIZ
cc      CHC = 1.00*CSIZ
cc      CHS = 0.85*CSIZ
C
      CSIZ1 = CSIZ*0.8
      CHN = 1.20*CSIZ1
      CHC = 1.00*CSIZ1
      CHS = 0.85*CSIZ1
C
      XL = XL0
      YL = YL0
      YSPACE = 2.2*CHC
C
      CALL GETCOLOR(ICOL0)
C
      CALL NEWPEN(3)
      XPLT = XL +  6.0*CHC - 0.5*FLOAT(NNAME)*CHN
      YL = YL - CHN
      CALL PLCHAR(XPLT,YL,CHN,NAME,0.0,NNAME)
C
      YL = YL - 0.2*CSIZ
      CALL NEWPEN(2)
C
      YL = YL - YSPACE
      CALL PLCHAR(XL        ,YL,CHC,'q   = ',0.0,6)
      CALL PLMATH(XL+1.0*CHC,YL,CHC, '&'    ,0.0,1)
      CALL PLNUMB(XL+5.0*CHC,YL,CHC, QINF   ,0.0,4)
C
      YL = YL - YSPACE
      CALL PLCHAR(XL        ,YL,CHC,'q   = ',0.0,6)
      CALL PLCHAR(XL+1.0*CHC,YL,CHC, 'ref'  ,0.0,3)
      CALL PLNUMB(XL+5.0*CHC,YL,CHC, QREF   ,0.0,4)
C
      YL = YL - YSPACE
      CALL PLCHAR(XL        ,YL,CHC,'RPM = ',0.0,6)
      CALL PLNUMB(XL+6.0*CHC,YL,CHC, RPM    ,0.0,0)
C
      TTOTAL = TROTOR + TDUCT
C
      YL = YL - YSPACE
      CALL PLCHAR(XL        ,YL,CHC,'Ttotal=',0.0,7)
      CALL PLNUMB(XL+7.0*CHC,YL,CHC, TTOTAL ,0.0,2)
C
      YL = YL - YSPACE
      CALL PLCHAR(XL        ,YL,CHC,'Trotor=',0.0,7)
      CALL PLNUMB(XL+7.0*CHC,YL,CHC, TROTOR ,0.0,2)
C
      YL = YL - YSPACE
      CALL PLCHAR(XL        ,YL,CHC,'Tduct =',0.0,7)
      CALL PLNUMB(XL+7.0*CHC,YL,CHC, TDUCT,0.0,2)
C
      YL = YL - YSPACE
      CALL PLCHAR(XL        ,YL,CHC,'  FOM =',0.0,7)
      CALL PLNUMB(XL+7.0*CHC,YL,CHC, FOM,0.0,3)
C
C----- axisymmetric
c       YL = YL - YSPACE
c       CALL PLCHAR(XL        ,YL,CHC,'C  = ',0.0,5)
c       CALL PLSUBS(XL        ,YL,CHC, 'D'   ,0.0,1,PLCHAR)
c       CALL PLNUMB(XL+5.0*CHC,YL,CHC, CD    ,0.0,5)
C
C
      XL = XL0 - 1.0*CHC
      YL = YL  - 0.5*CHC
C
      DO 20 IEL = 1, NEL
        IF(NETYPE(IEL).EQ.0) GO TO 20
C
        IF(GAM(IEL).NE.0.0 .OR. SIG(IEL).NE.0.0) THEN
         YL = YL - YSPACE
         CALL NEWCOLOR(ICOL(IEL))
         CINDEX = CHAR(ICHAR('0')+IEL)
C
         IF    (NETYPE(IEL).EQ.1) THEN
          CALL PLMATH(XL,YL,CHC,'g  s  = ',0.0,8)
         ELSEIF(NETYPE(IEL).EQ.2) THEN
          CALL PLCHAR(XL,YL,CHC,'d  s  = ',0.0,8)
         ELSEIF(NETYPE(IEL).EQ.3) THEN
          CALL PLMATH(XL,YL,CHC,'G  S  = ',0.0,8)
         ELSEIF(NETYPE(IEL).EQ.4) THEN
          CALL PLCHAR(XL,YL,CHC,'D  S  = ',0.0,8)
         ELSEIF(NETYPE(IEL).EQ.5) THEN
          CALL PLCHAR(XL,YL,CHC,'D  S  = ',0.0,8)
         ELSEIF(NETYPE(IEL).EQ.6) THEN
          CALL PLCHAR(XL,YL,CHC,'D  S  = ',0.0,8)
         ENDIF
C
         CALL PLSUBS(XL+0.0*CHC,YL,CHC, CINDEX ,0.0,1,PLCHAR)
         CALL PLSUBS(XL+3.0*CHC,YL,CHC, CINDEX ,0.0,1,PLCHAR)
         CALL PLNUMB(XL+8.0*CHC,YL,CHC, GAM(IEL),0.0,4)
         CALL PLCHAR(999.      ,YL,CHC,'  '     ,0.0,2)
         CALL PLNUMB(999.      ,YL,CHC, SIG(IEL),0.0,4)
        ENDIF
 20   CONTINUE
C
      CALL NEWCOLOR(ICOL0)
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END



      SUBROUTINE PQVEC(LVEC,ITYPE)
C-------------------------------------------------------
C     Plots airfoil with 
C       ITYPE=1:  normal pressure force vectors
C       ITYPE=2:  velocity vectors
C       ITYPE=3:  normal velocity vectors
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      LOGICAL LVEC(NEX)
C
      INCLUDE 'MASKS.INC'
C
C---- additional left and top margins for plot
      DATA XADD, YADD / 0.06 , 0.04 /
C
      CALL PLTINI2
C
C---- available plotting box
      XBOX = MIN( XWIND - XMARG , XPAGE - 2.0*XMARG ) / SIZE
      YBOX = MIN( YWIND - YMARG , YPAGE - 2.0*YMARG ) / SIZE
      XBOX = MAX( XBOX , 0.05 )
      YBOX = MAX( YBOX , 0.05 )
C
C---- set geometric limits
      XMIN =  1.0E20
      XMAX = -1.0E20
      YMIN =  1.0E20
      YMAX = -1.0E20
      DO IEL=1, NEL
        IF(LVEC(IEL)) THEN
         DO IC=ICFRST(IEL), ICLAST(IEL)
           XMIN = MIN(XMIN,XC(IC))
           XMAX = MAX(XMAX,XC(IC))
           YMIN = MIN(YMIN,YC(IC))
           YMAX = MAX(YMAX,YC(IC))
         ENDDO
        ENDIF
      ENDDO
C
C---- set pressure vector scale VSF
      XRANGE = MAX(1.0E-9, XMAX-XMIN)
      YRANGE = MAX(1.0E-9, YMAX-YMIN)
C
      IF(ITYPE.EQ.1) THEN
       VFAC = PVFAC / MIN( XBOX/XRANGE , YBOX/YRANGE ) / (CPLMAX-CPLMIN)
      ELSE
       VFAC = QVFAC / MIN( XBOX/XRANGE , YBOX/YRANGE ) / (QPLMAX-QPLMIN)
      ENDIF
C
C---- set limits again, including pressure vectors
      DO IEL = 1, NEL
        IF(LVEC(IEL)) THEN
         DO IC = ICFRST(IEL), ICLAST(IEL)
           IF    (ITYPE.EQ.1) THEN
            DXR =  ABS(CPR(IC))*VFAC*ANC(1,IC)
            DYR =  ABS(CPR(IC))*VFAC*ANC(2,IC)
            DXL = -ABS(CPL(IC))*VFAC*ANC(1,IC)
            DYL = -ABS(CPL(IC))*VFAC*ANC(2,IC)
           ELSEIF(ITYPE.EQ.2) THEN
            DXR = QCR(1,IC)/QREF *VFAC
            DYR = QCR(2,IC)/QREF *VFAC
            DXL = QCL(1,IC)/QREF *VFAC
            DYL = QCL(2,IC)/QREF *VFAC
           ELSEIF(ITYPE.EQ.3) THEN
            QNR = QCR(1,IC)*ANC(1,IC) + QCR(2,IC)*ANC(2,IC)
            QNL = QCL(1,IC)*ANC(1,IC) + QCL(2,IC)*ANC(2,IC)
            DXR =  ABS(QNR)/QREF *ANC(1,IC)*VFAC
            DYR =  ABS(QNR)/QREF *ANC(2,IC)*VFAC
            DXL = -ABS(QNL)/QREF *ANC(1,IC)*VFAC
            DYL = -ABS(QNL)/QREF *ANC(2,IC)*VFAC
           ENDIF
C
           IF(ISPLOT(IEL).EQ. 0 .OR.
     &        ISPLOT(IEL).EQ. 1      ) THEN
            XMIN = MIN( XMIN, XC(IC)+DXR )
            XMAX = MAX( XMAX, XC(IC)+DXR )
            YMIN = MIN( YMIN, YC(IC)+DYR )
            YMAX = MAX( YMAX, YC(IC)+DYR )
           ENDIF
           IF(ISPLOT(IEL).EQ. 0 .OR.
     &        ISPLOT(IEL).EQ. 2      ) THEN
            XMIN = MIN( XMIN, XC(IC)+DXL )
            XMAX = MAX( XMAX, XC(IC)+DXL )
            YMIN = MIN( YMIN, YC(IC)+DYL )
            YMAX = MAX( YMAX, YC(IC)+DYL )
           ENDIF
         ENDDO
        ENDIF
      ENDDO
C
C---- set scale, offsets, to center airfoil+vectors in plot area
      XRANGE = MAX(1.0E-9, XMAX-XMIN)
      YRANGE = MAX(1.0E-9, YMAX-YMIN)
      XBOXA = XBOX*(1.0-2.0*XADD)
      YBOXA = YBOX*(1.0-2.0*YADD)
      GSF = MIN( XBOXA/XRANGE , 
     &           YBOXA/YRANGE  )
      XOFG = XMIN - 0.5*(XBOXA-GSF*XRANGE)/GSF - XADD/GSF
      YOFG = YMIN - 0.5*(YBOXA-GSF*YRANGE)/GSF - YADD/GSF
C
      VSF = VFAC*GSF
C
C===================================================================
C
      CALL GETCOLOR(ICOL0)
      CALL GETPAT(IPAT0)
C
C----- plot axis
      CALL NEWPEN(1)
      CALL NEWPAT(LMASK1)
      AXL = XMIN - 0.02*XRANGE
      AXR = XMAX + 0.02*XRANGE
      CALL PLOT((AXL-XOFG)*GSF,(0.0-YOFG)*GSF,3)
      CALL PLOT((AXR-XOFG)*GSF,(0.0-YOFG)*GSF,2)
      CALL NEWPAT(IPAT0)
C
C---- plot unit-length vector
      CHL = 0.8*CHGT
      CHS = 0.5*CHGT
      CALL NEWPEN(2)
      CALL ARROW(XADD,YADD,VSF,0.0)
      XPLT = XADD + VSF + CHL
      YPLT = YADD + 0.5*CHL
      IF(ITYPE.EQ.1) THEN
       CALL PLCHAR(XPLT        ,YPLT        ,CHL,'C =1'   ,0.0,4)
       CALL PLCHAR(XPLT+0.9*CHL,YPLT-0.3*CHL,CHS, 'p'     ,0.0,1)
      ELSE
       CALL PLCHAR(XPLT        ,YPLT        ,CHL,'Q/Q  =1',0.0,7)
       CALL PLCHAR(XPLT+2.9*CHL,YPLT-0.3*CHL,CHS,    'ref',0.0,3)
      ENDIF
C
C
C---- plot wake grid
      CALL PLTGRID(SSIZEL,
     &             XG,YG,IX,II,JJ,
     &             -XOFG,-YOFG,GSF)
C
C---- plot blade or disk outlines
      CALL GETCOLOR(ICOL0)
      CALL NEWPEN(3)
      DO N = 1, NROTOR  
       IF(IRTYPE(N).EQ.1) THEN
        IEL = IELROTOR(N)
        CALL NEWCOLOR(ICOLEL(IEL))
        CALL PLTRBLD(SSIZEL,XPAXIS,
     &               XRC(1,N),YRC(1,N),ZER,ZER,NRC,
     &               -XOFG,-YOFG,GSF)
ccc     &               XOFA,YOFA,FACA)
       ELSEIF(IRTYPE(N).EQ.2) THEN
        IEL = IELROTOR(N)
        CALL NEWCOLOR(ICOLEL(IEL))
C
C---- get blade geometry at rotor points (esotec)
        CALL GETROTG(N)
        CALL PLTRBLD(SSIZEL,XPAXIS,
     &               XRP(1,N),YRP(1,N),CHRP(1,N),BETARP(1,N),NRP,
     &               -XOFG,-YOFG,GSF)
C
C        CALL PLTRBLD(SSIZEL,XPAXIS,
C     &               XRC(1,N),YRC(1,N),CHR(1,N),BETAR(1,N),NRC,
C     &               -XOFG,-YOFG,GSF)
       ENDIF
      END DO
      CALL NEWCOLOR(ICOL0)
C
C
      DO 50 IEL = 1, NEL
        IF(.NOT.LVEC(IEL)) GO TO 50
C
        CALL NEWPEN(3)
C
        IF( NETYPE(IEL).EQ.3 .OR.
     &      NETYPE(IEL).EQ.4   ) THEN
cc        IF(IPFRST(IEL).EQ.IPLAST(IEL)) THEN
C------- point source... plot with symbol
         CALL NEWCOLORNAME('MAGENTA')
         IP = IPFRST(IEL)
         SSIZ = SSIZEL(NETYPE(IEL))
         CALL PLSYMB((XP(IP)-XOFG)*GSF,(YP(IP)-YOFG)*GSF,SSIZ,1,0.0,0)
         CALL NEWCOLOR(ICOL0)
        ELSE
C------- usual paneled contour
         IF( NETYPE(IEL).EQ.1 .OR.
     &       NETYPE(IEL).EQ.2 .OR.
     &       NETYPE(IEL).EQ.5 .OR.
     &       NETYPE(IEL).EQ.6 .OR.
     &       NETYPE(IEL).EQ.7     ) THEN
          CALL NEWPAT(LMASK1)
         ENDIF
C
         CALL NEWCOLOR(ICOL0)
         IP = IPFRST(IEL)
         CALL PLOT((XP(IP)-XOFG)*GSF,(YP(IP)-YOFG)*GSF,3)
         DO IP = IPFRST(IEL)+1, IPLAST(IEL)
           CALL PLOT((XP(IP)-XOFG)*GSF,(YP(IP)-YOFG)*GSF,2)
         ENDDO
         CALL NEWPAT(IPAT0)
C
         IF(LVPLT) THEN
          SSIZ = SSIZEL(NETYPE(IEL))
          DO IP = IPFRST(IEL), IPLAST(IEL)
           CALL PLSYMB((XP(IP)-XOFG)*GSF,(YP(IP)-YOFG)*GSF,SSIZ,1,0.0,0)
          ENDDO
         ENDIF
        ENDIF
C
C------ plot pressure vectors
        CALL NEWPEN(2)
        DO IC = ICFRST(IEL), ICLAST(IEL)
          IF(ISPLOT(IEL).EQ. 0 .OR.
     &       ISPLOT(IEL).EQ. 1      ) THEN
           CALL NEWCOLORNAME('CYAN')
           IF    (ITYPE.EQ.1) THEN
            DXR = -CPR(IC)*VSF*ANC(1,IC)
            DYR = -CPR(IC)*VSF*ANC(2,IC)
           ELSEIF(ITYPE.EQ.2) THEN
            DXR = QCR(1,IC)/QREF *VSF
            DYR = QCR(2,IC)/QREF *VSF
           ELSEIF(ITYPE.EQ.3) THEN
            QNR = QCR(1,IC)*ANC(1,IC) + QCR(2,IC)*ANC(2,IC)
            DXR = QNR/QREF *ANC(1,IC)*VSF
            DYR = QNR/QREF *ANC(2,IC)*VSF
           ENDIF
           XX = (XC(IC)-XOFG)*GSF
           YY = (YC(IC)-YOFG)*GSF
           IF((ITYPE.EQ.1 .AND. CPR(IC).GE.0.0) .OR.
     &        (ITYPE.EQ.3 .AND. QNR    .LE.0.0)     ) THEN
            CALL ARROW(XX-DXR,YY-DYR,DXR,DYR)
           ELSE
            CALL ARROW(XX    ,YY    ,DXR,DYR)
           ENDIF
          ENDIF
C
          IF(ISPLOT(IEL).EQ. 0 .OR.
     &       ISPLOT(IEL).EQ. 2      ) THEN
           CALL NEWCOLORNAME('YELLOW')
           IF    (ITYPE.EQ.1) THEN
            DXL =  CPL(IC)*VSF*ANC(1,IC)
            DYL =  CPL(IC)*VSF*ANC(2,IC)
           ELSEIF(ITYPE.EQ.2) THEN
            DXL = QCL(1,IC)/QREF *VSF
            DYL = QCL(2,IC)/QREF *VSF
           ELSEIF(ITYPE.EQ.3) THEN
            QNL = QCL(1,IC)*ANC(1,IC) + QCL(2,IC)*ANC(2,IC)
            DXL = QNL/QREF *ANC(1,IC)*VSF
            DYL = QNL/QREF *ANC(2,IC)*VSF
           ENDIF
           XX = (XC(IC)-XOFG)*GSF
           YY = (YC(IC)-YOFG)*GSF
           IF((ITYPE.EQ.1 .AND. CPL(IC).GE.0.0) .OR.
     &        (ITYPE.EQ.3 .AND. QNL    .GE.0.0)     ) THEN
            CALL ARROW(XX-DXL,YY-DYL,DXL,DYL)
           ELSE
            CALL ARROW(XX    ,YY    ,DXL,DYL)
           ENDIF
          ENDIF
        ENDDO
C
 50   CONTINUE
C
C---- plot stagnation point
      IF(LPSTAG) THEN
        CALL NEWCOLORNAME('red')
        IEL = 2
        XPLT = (XSTG(IEL)-XOFG)*GSF
        YPLT = (YSTG(IEL)-YOFG)*GSF
        SSIZ = 2.0*SSIZEL(NETYPE(IEL))
        CALL PLSYMB(XPLT,YPLT,SSIZ,1,0.0,0)
      ENDIF
C
c      IF(LWALL) THEN
c       CALL NEWCOLORNAME('brown')
c       CALL PLTWAL(IGEOM,
c     &             XWALL,YWALL,AWALL,
c     &             -XOFG,-YOFG,GSF )
c      ENDIF
C
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
      IPTYPE = JAPLT
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! PQVEC

c--- extra code bits...
c      IF(IPROTCB2.GT.0) CALL PLTROTR2(SSIZEL,
c     &             XP(IPROTCB2),YP(IPROTCB2),XP(IPROTDW2),YP(IPROTDW2),
c     &             XOFA,YOFA,FACA)







