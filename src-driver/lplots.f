C**********************************************************************
C    Plot routines for ESLOFT
C    Copyright (C) April 2009, Philip Carter, Esotec Developments.
C    philip (at) esotec (dot) org
C    Acknowledgements to Hal Youngren.
C
C    First released with DFDC v070-ES2.
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


      SUBROUTINE AFPLOT(NDSK)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C------------------------------------------------------
C     Plots esloft airfoil set
C------------------------------------------------------
      EXTERNAL PLCHAR,PLMATH
      DIMENSION XLIN(4), YLIN(4)
      DIMENSION W1(NAFX), W2(NAFX), W3(NAFX), 
     &          W4(NAFX), W5(NAFX), W6(NAFX) 
C
      DIMENSION XZ(NPX),YZ(NPX),XZP(NPX),YZP(NPX),SZ(NPX)
      DIMENSION THUP(NAFX),THDN(NAFX)
C
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
C---- set blade disk to plot
      N = NDSK
C
C---- character size for axis numbers, labels
      CS  = CSIZE * 0.75
      CSL = CS    * 1.2
C
C---- get extremities of airfoils (thanks to Dan Aingie)
C
      DO I=1,NAF
        THUP(I) = MAXVAL(YYE(:,I))
        THDN(I) = ABS(MINVAL(YYE(:,I)))
      ENDDO
C
      CALL PLOTABS(1.0,1.0,-3)
      CALL NEWPEN(2)
C
      SCL  = 0.65
      XSF  = SCL
      YSF  = SCL
      XOFF = 0.0
      YOFF = -SCL*THDN(1)
      XOFFN= 0.31 * SCL
      SPACEMX= 0.12   ! max space between airfoils
      YTOP   = 1.35   ! max height of top airfoil
C
C---- figure out spacing
C
      THTOT= 0.0
      DO I=1,NAF
         THTOT= THTOT + THUP(I)+THDN(I)
      ENDDO
C
      SPACET  = (YTOP-THTOT)/FLOAT(NAF-1)
      AFSPACE = MIN(SPACEMX,SPACET)
C
C----Initialize splines and plot
C
      NZ=NPP
C
      DO 200 I=1,NAF
        DO J=1,NZ
          XZ(J)=XXE(J,I)
          YZ(J)=YYE(J,I)
        ENDDO
C
        IF(I.NE.1) THEN
          YOFF = YOFF - (THUP(I-1)+THDN(I)+AFSPACE)
        ENDIF
        YOFFN = -SCL*(YOFF+(THDN(I)-THUP(I))/2.0) - 1.0*CS
C         
        CALL SCALC (XZ,YZ, SZ,NZ)
        CALL SEGSPL(XZ,XZP,SZ,NZ)
        CALL SEGSPL(YZ,YZP,SZ,NZ)
C
       CALL XPLTAIR(XZ,XZP,YZ,YZP,SZ,NZ,XOFF,XSF,YOFF,YSF,7)
       CALL NEWPEN(2)
       CALL NEWCOLOR(1)
       CALL PLCHAR(XOFFN,YOFFN,CS,ENAME(I,N),0.0,-1)
C
 200  CONTINUE
C
      CALL STRIP(SVERSION,NV)
      XL1 = 0.0
      XL2 = 0.7 - CS*(6.0 + FLOAT(NV))
      YL  = 0.0
      YTIT= YPAGE * 1.22
      XTIT= XPAGE * 0.06

c      XL1=0.0
c      XL2=0.0
c      YL =0.0
c      YTIT= 0.94*YWIND
c      XTIT= 0.8
C
      CALL PLOTABS(XTIT,YTIT,-3)
      CALL NEWPEN(2)
      CALL NEWCOLOR(1)
      CALL PLCHAR(XL1,YL,CSL,'ESLOFT Parents: ',0.0,16)
      CALL PLCHAR(999.,YL,CSL,SNAME,0.0,-1)
C
c      CALL PLOTABS(0.78*XWIND,YTIT,-3)
      CALL NEWCOLORNAME('cyan')
      CALL PLCHAR(XL2,YL,CS,'DFDC v',0.0,6)
      CALL PLCHAR(999.,YL,CS,SVERSION,0.0,-1)
C
      CALL NEWCOLOR(1)
      CALL PLOTABS(XTIT,YTIT-0.08,-3)
      CALL PLOT (0.0,0.0,3)
      CALL PLOT (0.7,0.0,2)
C
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! AFPLOT
C
C------------------------------------------------------------------




      SUBROUTINE LOFTPLT(NR)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C---------------------------------------------------
C     Plots loft stations and thickness distributions
C---------------------------------------------------
      CHARACTER*80 PNAME
      CHARACTER*30 TLEG
      DIMENSION Z1(NLSX), Z2(NLSX), ZY(NLSX), ZS(NLSX), 
     &          Z1P(NLSX),Z2P(NLSX),ZYP(NLSX),
     &          THMM(NLSX),TOCP(NLSX)
      LOGICAL ERROR
C
      EXTERNAL PLCHAR,PLMATH
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
C
C---- set number of loft stations
      NL = NLOFT2
C
C----------------------------------------------------------------
C
C---- plot aspect ratio
      PLPAR = 0.67
C
C---- character size for axis numbers, labels
      CS  = 1.15*CSIZE !  t,t/c,R axes
      CSL = CS*1.1    !  title and legends
      CSS = CS*1.0    !  blade stuff
      CSV = CS*0.95    !  version
C
C---- blade plot axis
      BLAXIS = 1.13
C
C---------------------------------------------------------------
C---- set x axis parameters
C
      RMAX = YLOFT(NL)
      RCM  = RMAX*100.0
C
      IF(RCM.LE.10.0) THEN
         DLAN = 1.0
      ELSEIF(RCM.LE.20.0) THEN
         DLAN = 2.0
      ELSEIF(RCM.LE.50.0) THEN
         DLAN = 5.0
      ELSEIF(RCM.LE.100.0) THEN
         DLAN = 10.0
      ELSEIF(RCM.LE.150.0) THEN
         DLAN = 15.0
      ELSEIF(RCM.LE.200.0) THEN
         DLAN = 20.0
      ELSEIF(RCM.LE.300.0) THEN
         DLAN = 30.0
      ELSEIF(RCM.LE.400.0) THEN
         DLAN = 40.0
      ELSEIF(RCM.LE.500.0) THEN
         DLAN = 50.0
      ELSEIF(RCM.LE.750.0) THEN
         DLAN = 75.0
      ELSEIF(RCM.LE.1000.0) THEN
         DLAN = 100.0
      ELSE
         DLAN = 200.0
      ENDIF
C
C      XLMAX = INT(RMAX+DLAN-0.001)
C      DLXAN = DLAN/(XLMAX*RMAX)
C
      NXINT = INT(0.99*RCM/DLAN)+1
      DLX = 1.0/NXINT
      XLSCAL = 100.0/(NXINT*DLAN)
      XGSCAL = XLSCAL/100.0
C
      RMAXIS = NXINT*DLAN/100.  !  for setting legend locations
C
C---- load dummies for thickness_mm and t/c_% and set up axes
C
      DO I=1,NL
         THMM(I)=THLOFT(I)*1000.0
         TOCP(I)=TCLOFT(I)*100.0
      ENDDO
C
      CALL MINMAX(NL,THMM,THMIN,THMAX)
      CALL MINMAX(NL,TOCP,TCMIN,TCMAX)
C
      CALL AXISADJ(THMIN,THMAX,THSPAN,DTH,NTH)
      CALL AXISADJ(TCMIN,TCMAX,TCSPAN,DTC,NTOC)
C
C---- plot x axis and grid
C
      CALL GETCOLOR(ICOL0)
      CALL NEWCOLORNAME('black')
      CALL NEWPEN(2)
      CALL XAXIS(0.0,0.0, 1.0,DLX, 0.0,DLAN, CS,-1)
C
c      CALL NEWPEN(2)
      CALL PLCHAR(0.5-1.5*CSL,-3.5*CSL,CSL,'r cm',0.0,4)
C
      CALL NEWPEN(1)
      NLXG = NXINT
      NLYG = (NTOC-1)
      DLY= 1.0/NLYG
C
      CALL NEWCOLORNAME('cyan')
      CALL PLGRID(0.0,0.0, NLXG,DLX, NLYG,DLY*PLPAR, LMASK2 )
      CALL PLFLUSH
C
C
C---- find legend locations
C
      XFAC1 = 0.4
      XFAC2 = 0.7
C
      DO I=1,NL
         FAC = (YLOFT(I)-YLOFT(1))/(YLOFT(NL)-YLOFT(1))
         IF(FAC.LE.XFAC1) IL1 = I
         IF(FAC.LE.XFAC2) IL2 = I
      ENDDO
C
      YTLEG = PLPAR + 0.8*CS
      XTLEG = 1.0 - 25.5*CS
C
      IF    (ITTYPE.EQ.1) THEN
         TLEG = '   Linear thickness/chord'
      ELSEIF(ITTYPE.EQ.2) THEN
         TLEG = 'Parabolic thickness/chord'
      ELSEIF(ITTYPE.EQ.3) THEN
         TLEG = '  Splined thickness/chord'
      ELSEIF(ITTYPE.EQ.4) THEN
         TLEG = 'Linear thickness'
      ELSEIF(ITTYPE.EQ.5) THEN
         TLEG = 'Parabolic thickness'
      ELSEIF(ITTYPE.EQ.6) THEN
         TLEG = 'Splined thickness'
      ENDIF
C
C---- plot thickness vs R
C
      CALL NEWCOLORNAME('orange')
      THOFF = THMIN
      THFAC = PLPAR/THSPAN
      YAX   = PLPAR
      DYAX  = YAX/FLOAT(NTH-1)
C  
      CALL NEWPEN(2)
      CALL YAXIS(0.0,0.0,YAX,DYAX,THMIN,DTH, CS,1)
C
      CALL NEWPEN(4)
      CALL XYLINE(NL,YLOFT,THMM,0.0,XLSCAL,THOFF,THFAC,1)
C
      CALL NEWPEN(2)
      XPLT = 0.0 - 4.5*CSL
      YPLT = THFAC*(THMAX-1.5*DTH-THOFF) - 0.5*CSL
      CALL PLCHAR(XPLT,YPLT,CSL,'t mm',0.0,4)
C
      XPLT =  YLOFT(IL1)/RMAXIS
      YPLT =   0.4*CSL + THFAC*(THMM(IL1)-THOFF)
      CALL PLCHAR(XPLT,YPLT,CSL,'t',0.0,1)
C
      IF(ITTYPE.GE.4) THEN
        CALL PLCHAR(0.0,YTLEG,CS,TLEG,0.0,20)
      ENDIF
C
      CALL PLFLUSH
C
C
C---- plot thickness/chord vs R
C
      CALL NEWCOLORNAME('green')
      TCOFF = TCMIN
      TCFAC = PLPAR/TCSPAN
      YAX   = PLPAR
      DYAX  = YAX/FLOAT(NTOC-1)
C
      CALL NEWPEN(2)
      CALL YAXIS(1.0,0.0,YAX,DYAX,TCMIN,DTC,-CS,1)
C
      CALL NEWPEN(4)
      CALL XYLINE(NL,YLOFT,TOCP,0.0,XLSCAL,TCOFF,TCFAC,1)
C
      CALL NEWPEN(2)
      XPLT = 1.0 + 1.0*CSL
      YPLT = TCFAC*(TCMAX-1.5*DTC-TCOFF) - 0.5*CSL
      CALL PLCHAR(XPLT,YPLT,CSL,'t/c %',0.0,5)
C
      XPLT =  YLOFT(IL2)/RMAXIS
      YPLT =   0.4*CSL + TCFAC*(TOCP(IL2)-TCOFF)
      CALL PLCHAR(XPLT,YPLT,CSL,'t/c',0.0,3) 
C
      IF(ITTYPE.LE.3) THEN
        CALL PLCHAR(XTLEG,YTLEG,CS,TLEG,0.0,25)
      ENDIF
C
      CALL PLFLUSH
C
C
C---- plot untwisted blade view
C
      DO I=1,NL
        CDAXIS= AXHUB-(AXHUB-AXTIP)*
     &      (YLOFT(I)-YLOFT(IHUB))/(YLOFT(ITIP)-YLOFT(IHUB))
C
        Z1(I) =  CDAXIS     * CDLOFT(I) * XLSCAL
        Z2(I) = (CDAXIS-1.0)* CDLOFT(I) * XLSCAL
        ZY(I) =  YLOFT(I)   * XLSCAL
      ENDDO
C
      CALL PLOT(0.0,BLAXIS,-3)  ! BLAXIS locates blade, set in header
C
      CALL NEWPEN(1)
      CALL NEWCOLORNAME('blue')
C
C----- plot radial axis
C
      XTIPAX = RMAX * XLSCAL
C
C      CALL PLOT(0.0,0.0,3)
      CALL PLOT(XTIPAX,0.0,2)
C
C----- plot axis tick mark
C
      CALL PLOT(0.0,-.005,3)
      CALL PLOT(0.0,0.005,2)
C
C----- plot blade shape
C
C      CALL NEWPEN(4)
C      CALL XYLINE(NL,WY,W1,0.0,1.0,0.0,1.0,1)
C      CALL XYLINE(NL,WY,W2,0.0,1.0,0.0,1.0,1)
C
C---- plot blade LE shape
C
      CALL NEWPEN(2)
      NT = 20   ! plot lines per loft interval
C
      CALL SCALC (ZY,Z1,ZS,NL)
      CALL SEGSPL(ZY,ZYP,ZS,NL)
      CALL SEGSPL(Z1,Z1P,ZS,NL)
C
      DO 60 I=2, NL
        DS = ZS(I) - ZS(I-1)
        CALL PLOT(ZY(I-1),Z1(I-1),3)
        DO IT=1, NT
          ST = ZS(I-1) + DS*FLOAT(IT)/FLOAT(NT)
          XT = SEVAL(ST,ZY,ZYP,ZS,NL)
          YT = SEVAL(ST,Z1,Z1P,ZS,NL)
          CALL PLOT(XT,YT,2)
        ENDDO
   60 CONTINUE
C
C---- plot blade TE shape
C
      CALL SCALC (ZY,Z2,ZS,NL)
      CALL SEGSPL(ZY,ZYP,ZS,NL)
      CALL SEGSPL(Z2,Z2P,ZS,NL)
C
      DO 70 I=2, NL
        DS = ZS(I) - ZS(I-1)
        CALL PLOT(ZY(I-1),Z2(I-1),3)
        DO IT=1, NT
          ST = ZS(I-1) + DS*FLOAT(IT)/FLOAT(NT)
          XT = SEVAL(ST,ZY,ZYP,ZS,NL)
          YT = SEVAL(ST,Z2,Z2P,ZS,NL)
          CALL PLOT(XT,YT,2)
        ENDDO
   70 CONTINUE
C
C
C---- plot root and tip lines
C
      CALL PLOT(ZY(1),Z1(1),3)
      CALL PLOT(ZY(1),Z2(1),2)
C      CALL PLOT(YLOFT(1),CHRP(1,N)*(XPAXIS-1.0),2)
C
      CALL PLOT(ZY(NL),Z1(NL),3)
      CALL PLOT(ZY(NL),Z2(NL),2)
C
C---- plot blade label
C
      YBTIT= -2.1*CSS
      CALL PLCHAR(0.0,YBTIT,CSS,'Untwisted',0.0,9)      
      CALL PLFLUSH
C
C---- plot loft stations
C
      XSTN = MIN(Z2(1),Z2(NL)) - 0.04
C
      CALL GETPAT(IOLDPAT)
      CALL NEWPAT(LMASK2)
      CALL NEWPEN(1)
      CALL NEWCOLORNAME('black')
C
      DO I=1,NL
         CALL PLOT(ZY(I),0.0,3)
         CALL PLOT(ZY(I),XSTN,2)
      ENDDO
      CALL NEWPAT(IOLDPAT)
C
C---- plot station numbers
C
      STNUM = -2.0*CSS
      XNUM  = XSTN +STNUM
      XHT   = XNUM +STNUM
      CALL NEWPEN(2)
C
      DO I=1,NL
         IF(I.LE.9) THEN
            YNUM = ZY(I)-0.5*CSS
         ELSE
            YNUM = ZY(I)-CSS
         ENDIF
C
         RI=FLOAT(I)
         CALL PLNUMB(YNUM,XNUM,CSS,RI,0.0,-1)
C
         YHT = ZY(I)-0.5*CSS        
         IF(I.EQ.IHUB) THEN
            CALL PLCHAR(YHT,XHT,CSS,'H',0.0,1)
         ELSEIF(I.EQ.ITIP) THEN
            CALL PLCHAR(YHT,XHT,CSS,'T',0.0,1)
         ENDIF
      ENDDO
C
      CALL PLCHAR(0.0,XNUM,CSS,'Stations',0.0,8)
      CALL PLFLUSH
C
C---- plot parent locations
C
      CALL NEWCOLORNAME('green')
      XAF = MAX(Z1(1),Z1(NL)) + 0.04
      XANUM = XAF + 1.0*CSS
      CALL NEWPAT(LMASK2)
      CALL NEWPEN(2)
c
      DO J=1,NAF
        K=NTC(J)
        DO I=1,K
          CALL NEWPEN(1)
          WAF = TCLOC(J,I) * XLSCAL 
          CALL PLOT(WAF,0.0,3)
          CALL PLOT(WAF,XAF,2)
          YNUM = WAF-0.5*CSS
          RJ=FLOAT(J)
          CALL NEWPEN(2)
          CALL PLNUMB(YNUM,XANUM,CSS,RJ,0.0,-1)
        ENDDO
      ENDDO
C
      CALL NEWPAT(IOLDPAT)
      CALL PLCHAR(0.0,XANUM,CSS,'Parents',0.0,7)
C
C---- plot header
C
      CALL STRIP(SVERSION,NV)
      ZV = 6.0 + FLOAT(NV)
C
      YLIN  = 0.27
      YTIT  = YLIN + 0.8*CSL
      XTIT = 1.0 - ZV*CSV
C
      ISTN=1
      CALL GETLNAMES(LONAME,NROTOR,NR,ISTN,1,PNAME)
C
      CALL NEWPEN(2)
      CALL NEWCOLORNAME('black')
      CALL PLOT(0.0,YLIN,3)
      CALL PLOT(1.0,YLIN,2)
C
      CALL PLCHAR(0.0,  YTIT,CSL,'Lofted Blade: ',0.0,14)
      CALL PLCHAR(999., YTIT,CSL, PNAME          ,0.0,-1)
      CALL NEWCOLORNAME('cyan')
      CALL PLCHAR(XTIT, YTIT,CSV,'DFDC v' ,0.0,6)
      CALL PLCHAR(999., YTIT,CSV, SVERSION       ,0.0,-1)
C
C
C---- all done
C
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      XYOFF(1) = 0.
      XYOFF(2) = 0.
      XYFAC(1) = 1.0
      XYFAC(2) = PLPAR
C
      RETURN
      END !  LOFTPLT
C
C--------------------------------------------------------------------




      SUBROUTINE XFOILPLOT(NDSK,INPLOT)
C-----------------------------------------------------------------------
C     Control center for xfoil plotting routines
C     Parameters taken from storage for multi-plot blow, reset, replot
C
C     INPLOT = 1  no blowup
C     INPLOT = 2  blowup
C     INPLOT = 3  replot current plot at current scale and origin
C
C     IXTYPE = 1  parent airfoils       |
C     IXTYPE = 2  blended airfoils      | global variable
C     IXTYPE = 3  transformed airfoils  |
C-----------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      PARAMETER (NCOLX=45)
      DIMENSION ICOLMAP(NCOLX)
      N = NDSK
C
      DATA ICOLMAP/
     &     1,7,9,4,6,8,3,5,10,
     &     1,7,9,4,6,8,3,5,10,
     &     1,7,9,4,6,8,3,5,10,
     &     1,7,9,4,6,8,3,5,10,
     &     1,7,9,4,6,8,3,5,10/
C
      IF(IXTYPE.EQ.1) THEN      !  plotting parent airfoils
         LGGRIDX = .TRUE.
         LGTICKX = .TRUE.
         LGPARMX = .TRUE.
         IF(INPLOT.EQ.2) LGPARMX = .FALSE.
         PLOTARX =  0.55
      ELSEIF(IXTYPE.EQ.2) THEN  ! plotting blended sections
         LGGRIDX = .TRUE.
         LGTICKX = .TRUE.
         LGPARMX = .TRUE.
         IF(INPLOT.EQ.2) LGPARMX = .FALSE.
         PLOTARX =  0.55
      ELSEIF(IXTYPE.EQ.3) THEN  ! plotting transformed sections
         LGGRIDX = .TRUE.
         LGTICKX = .TRUE.
         LGPARMX = .FALSE.
         PLOTARX =  0.70
         IF(OUTFAC.EQ.1.0.OR.OUTFAC.EQ.100.0.OR.
     &      OUTFAC.EQ.1000.0) THEN
            TRANFAC = 100.0
         ELSE
            TRANFAC = MTI
         ENDIF
      ELSE
         WRITE(*,*)
         WRITE(*,*)'Plot type error. Aborted.'
         RETURN
      ENDIF      
C
C---- now the stuff that changes each plot
C
      NBB = NPP
      DO I=1,NXPLOT
         IAF = IXPLOT(I)
C
         IF(IXTYPE.EQ.1) THEN      !  parent airfoils
            NAMEX= ENAME(IAF,N)
            DO J=1,NBB
              XBB(J)=XXE(J,IAF)
              YBB(J)=YYE(J,IAF)
            ENDDO
         ELSEIF(IXTYPE.EQ.2) THEN  !  blended sections
            CALL GETLNAMES(LONAME,NROTOR,NDSK,IAF,6,NAMEX)
            DO J=1,NBB
              XBB(J)=BLXX(J,IAF)
              YBB(J)=BLYY(J,IAF)
            ENDDO
         ELSE                      !  transformed sections
            CALL GETLNAMES(LONAME,NROTOR,NDSK,IAF,6,NAMEX)
            DO J=1,NBB
              XBB(J)=TRXX(J,IAF)*TRANFAC
              YBB(J)=TRYY(J,IAF)*TRANFAC
            ENDDO
         ENDIF
C
         CALL SCALC (XBB,YBB, SBB,NBB)
         CALL SEGSPL(XBB,XBBP,SBB,NBB)
         CALL SEGSPL(YBB,YBBP,SBB,NBB)
C
C---- load geometry parameters
C
         IF(IXTYPE.EQ.1) THEN      !  parent airfoils
            AREABX  = PAREA(IAF,N)
            RADBLEX = PLERAD(IAF,N)
            ANGBTEX = PANGTE(IAF,N)
            THICKBX = TOCE(IAF,N)
            CAMBRBX = CAMLO(IAF,N)
            TETHX   = TETH(IAF,N)
         ELSEIF(IXTYPE.EQ.2) THEN  !  blended sections
            AREABX  = BLENDATA(1,IAF)
            RADBLEX = BLENDATA(2,IAF)
            ANGBTEX = BLENDATA(3,IAF)
            THICKBX = BLENDATA(4,IAF)
            CAMBRBX = BLENDATA(6,IAF)
            TETHX   = BLENDATA(8,IAF)
         ENDIF
C
         ICOLOR = ICOLMAP(I)
C
C---- data is ready, now to direct it
C
         IF(I.EQ.1) THEN
           IF(INPLOT.EQ.1) THEN      ! plot at normal scale and origin
             CALL PLTINIX
             CALL GOFINIX(PLOTARX)
             CALL PLOTGX(NDSK,ICOLOR)
           ELSEIF(INPLOT.EQ.2) THEN  ! plot at specified scale and origin
             CALL GOFSETX
             CALL PLTINIX
             CALL PLOTGX(NDSK,ICOLOR)
           ELSEIF(INPLOT.EQ.3) THEN  ! replot at current scale and origin
             CALL PLTINIX
             CALL PLOTGX(NDSK,ICOLOR)
           ELSE
             WRITE(*,*)
             WRITE(*,*)'XFOILPLOT: plot index outside range'
             RETURN             
           ENDIF
         ELSE
           CALL OVERX(ICOLOR)        ! all subsequent plots
         ENDIF
C
      ENDDO
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END   ! XFOILPLOT
C
C----------------------------------------------------------------------




      SUBROUTINE TCMULT(NAFX,NPX,IAF,XX,YY,NP,
     &                  TYMAX,TXMAX, CYMAX,CXMAX)
C-----------------------------------------------------------
C     Reports buffer airfoil thickness and camber
C     splined from 2D airfoil arrays
C     based on xfoil TCBUF
C-----------------------------------------------------------
      PARAMETER (IBX= 2000)  ! keep this at least 8x NPX 
      DIMENSION XCM(IBX),YCM(IBX),XTK(IBX),YTK(IBX),
     &          YCMP(IBX),YTKP(IBX)
      DIMENSION XX(NPX,NAFX),YY(NPX,NAFX)
      DIMENSION XB(NPX),XBP(NPX),YB(NPX),YBP(NPX),SB(NPX)
C
      NB= NP
      DO I=1,NB
        XB(I) = XX(I,IAF)
        YB(I) = YY(I,IAF)
      ENDDO
C
      CALL SCALC (XB,YB, SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
C--- find the current buffer airfoil camber and thickness
C
      CALL GETCAM(XCM,YCM,NCM,XTK,YTK,NTK,
     &            XB,XBP,YB,YBP,SB,NB )
      CALL GETMAX(XCM,YCM,YCMP,NCM,CXMAX,CYMAX)
      CALL GETMAX(XTK,YTK,YTKP,NTK,TXMAX,TYMAT)
C
      TYMAX= 2.0*TYMAT
C
      WRITE(*,1000) TYMAX,TXMAX, CYMAX,CXMAX
 1000 FORMAT( ' Max thickness = ',F8.4,'  at x = ',F7.3,
     &       /' Max camber    = ',F8.4,'  at x = ',F7.3)
C
      RETURN
      END ! TCMULT
C
C------------------------------------------------------------------------





      SUBROUTINE TCSING (XB,XBP,YB,YBP,SB,NB,
     &                   TYMAX,TXMAX,CYMAX,CXMAX)
C-----------------------------------------------------------
C     Reports splined airfoil thickness and camber
C     for single airfoil
C     based on xfoil TCBUF
C------------------------------------------------------------
      PARAMETER (IBX= 2000)  ! keep this at least 8x NPX
      DIMENSION XCM(IBX),YCM(IBX),XTK(IBX),YTK(IBX),
     &          YCMP(IBX),YTKP(IBX)
      DIMENSION XB(*),XBP(*),YB(*),YBP(*),SB(*)
C
C--- find the current buffer airfoil camber and thickness
C
      CALL GETCAM(XCM,YCM,NCM,XTK,YTK,NTK,
     &            XB,XBP,YB,YBP,SB,NB )
      CALL GETMAX(XCM,YCM,YCMP,NCM,CXMAX,CYMAX)
      CALL GETMAX(XTK,YTK,YTKP,NTK,TXMAX,TYMAT)
C
      TYMAX= 2.0*TYMAT
C
C      WRITE(*,1000) TYMAX,TXMAX, CYMAX,CXMAX
C 1000 FORMAT( ' Max thickness = ',F8.4,'  at x = ',F7.3,
C     &       /' Max camber    = ',F8.4,'  at x = ',F7.3)
C
      RETURN
      END ! TCSING
C
C------------------------------------------------------------------------




      SUBROUTINE PLOTLOFDATA(NR)
C----Plots loft data on blade disk # NR
C
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      LOGICAL LERROR, LXRELIMIT, LYRELIMIT
      LOGICAL LPLT_DATA, LEGEND2D
      CHARACTER*4 OPT
      CHARACTER*9 XVAR, YVAR
      CHARACTER DATALBL*20, LEGENDDATA*80
      CHARACTER*80 LINE, FNAME*128
      CHARACTER*80 TITLE1, TITLE2, PNAME
C
      DIMENSION RINPUT(10), XLIM(2), YLIM(2)
C
C--- Temporary data for plotting 
      PARAMETER (IXP=100)
      DIMENSION  NPTDATA(IXP)
      DIMENSION  INFODATA(2,IXP)
      DIMENSION  XDATA(IXP,10), YDATA(IXP,10)
      DIMENSION  LEGENDDATA(10)
C
      PAR = 0.80
      LUTEMP = 3
      LXRELIMIT = .TRUE.
      LYRELIMIT = .TRUE.
      XMINSPEC = 999.
      XMAXSPEC = 999.
      YMINSPEC = 999.
      YMAXSPEC = 999.
C      DATALBL = PLTTYPE
C
 1    FORMAT(A)
      NDATA = 0
C
C---- first plot
C
      ILPLOT =1
      NCASESX=1
      XVAR   = 'r cm'
      YVAR   = 't mm'
C
      IAF=1
      CALL GETLNAMES(LONAME,NROTOR,NR,IAF,1,PNAME)
C
      GO TO 100
C
C--- Make the plot
 5    IF(LXRELIMIT) THEN
        XLIM(1) = -999.
        XLIM(2) = -999.
      ENDIF
      IF(LYRELIMIT) THEN
        YLIM(1) = -999.
        YLIM(2) = -999.
      ENDIF
      IF(NDATA.GT.0) THEN 
       CALL PLOTXY2(IXP,
     &              NDATA,NPTDATA,
     &              INFODATA,XDATA,YDATA,
     &              XMINSPEC,XMAXSPEC,YMINSPEC,YMAXSPEC,
     &              XVAR,YVAR,XLIM,YLIM,
     &              LEGEND2D,LEGENDDATA,SVERSION,
     &              TITLE1, PNAME, LPLOT,
     &              PAR,CSIZE,XMARG,YMARG,NCOLOR,SCRNFR,IDEV)
       LPLOT = .TRUE.
       CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
      ENDIF
C  
C---- Display options 
 10   WRITE(*,*) ' '
      WRITE(*,20)
 20   FORMAT(/' Lofted Section Data Plots'
     &       /' ------------------------------------------'
     &       /'   1  blade thickness'
     &       /'   2  max thickness/chord'
     &       /'   3  max thickness/chord x/c'
     &       /'   4  max camber'
     &       /'   5  max camber x/c'
     &       /'   6  section area'
     &       /'   7  leading edge radius'
     &       /'   8  trailing edge thickness'
     &       /'   9  trailing edge angle'
     &       /'   10 zero-lift alpha - aero and esloft'
     &       /'   11 beta - dfdc and esloft'
     &       /'   12 chord'
     &       /' ------------------------------------------'
     &       /'   A  toggle abscissa cm/inches',
     &       /'   B  toggle local/global beta coordinates',
     &       /'   L  limits for plot',
     &       /'   Z  zoom plot with cursor',
     &       /'   R  reset plot limits',
     &       /'   AN annotate plot',
     &       /'   H  hardcopy plot',
     &       /'   W  write plot data to file',
C
     &      //' Select option (or <return>):  ', $)
C
      READ(*,1,ERR=10) OPT
      CALL LC2UC(OPT)
C
      IF(OPT.EQ.' ') THEN
        RETURN
C
C***************************************************
C---- Change Abscissa 
       ELSE IF(OPT.EQ.'A') THEN
          IF(XVAR.EQ.'r cm') THEN
             XVAR = 'r in'
          ELSE
             XVAR = 'r cm'
          ENDIF
C
        XMINSPEC = 999.
        XMAXSPEC = 999.
        LXRELIMIT = .TRUE.
C
        GO TO 100
C
C***************************************************
C---- toggle local/global coords
C
       ELSE IF(OPT.EQ.'B') THEN
        LCOORD = .NOT.LCOORD
        IF(LCOORD) THEN
           WRITE(*,*)
           WRITE(*,*)'Local beta angles selected'
        ELSE
           WRITE(*,*)
           WRITE(*,*)'Global beta angles selected'
        ENDIF
C
        IF(OMEGA(NR).GT.0.0) THEN
           WRITE(*,*)'Applies only to neg-BGAM disks'
           GO TO 10
        ENDIF
C
        ILPLOT= 11
        YVAR  = 'beta deg'
        NCASESX=2
        GO TO 100
C
C***************************************************
C---- Make hardcopy
       ELSE IF(OPT.EQ.'H') THEN
        CALL REPLOT(IDEVPS)
        GO TO 10
C
C***************************************************
C---- Write plot data to file
C
       ELSE IF(OPT.EQ.'W') THEN
        IF(NDATA.GT.0) THEN
          CALL ASKS('Enter plot data save filename^',FNAME)
          OPEN(LUTEMP,FILE=FNAME,STATUS='UNKNOWN',ERR=10)
C
          WRITE(LUTEMP,52) TITLE1,PNAME
C
          DO IP= 1,NDATA
            WRITE(LUTEMP,53) LEGENDDATA(IP),XVAR,YVAR
            NPT = NPTDATA(IP)
            DO I = 1,NPT
              WRITE(LUTEMP,51) XDATA(I,IP),YDATA(I,IP) 
            END DO
          ENDDO
C
          CLOSE(LUTEMP)
        ENDIF
C
        GO TO 10
C
 51    FORMAT(2(G13.6,2X))
 52    FORMAT(1X,A40,/,1X,A40)
 53    FORMAT(/,1X,A40,/,1X,A10,5X,A10)

C***************************************************
C---- Annotate plot
       ELSE IF(OPT.EQ.'AN') THEN
        CALL ANNOT(CSIZE)
        GO TO 10
C
C***************************************************
C---- Zoom plot limits
       ELSE IF(OPT.EQ.'Z') THEN
        CALL OFFSET2D(XLIM,YLIM)
        LXRELIMIT = .FALSE.
        LYRELIMIT = .FALSE.
        GO TO 5
C
C***************************************************
C---- Reset plot limits
       ELSE IF(OPT.EQ.'R') THEN
        LXRELIMIT = .TRUE.
        LYRELIMIT = .TRUE.
        GO TO 5
C
C***************************************************
C---- Set plot limits
       ELSE IF(OPT.EQ.'L') THEN
        WRITE(*,*) 'Limits for plot:'
        WRITE(*,32) XLIM
 32     FORMAT('Xlimits  xmin: ',G12.6,' xmax: ',G12.6)
        READ(*,1) LINE
        IF(LINE.NE.' ') THEN 
          READ(LINE,*,ERR=10) XLIM(1), XLIM(2)
          LXRELIMIT = .FALSE.
        ENDIF
        WRITE(*,34) YLIM
 34     FORMAT('Ylimits  ymin: ',G12.6,' ymax: ',G12.6)
        READ(*,1) LINE
        IF(LINE.NE.' ') THEN 
          READ(LINE,*,ERR=10) YLIM(1), YLIM(2)
          LYRELIMIT = .FALSE.
        ENDIF
        GO TO 5
C
C
C***************************************************
C---- Change Ordinate and plot
C
         ELSE IF(OPT.EQ.'1') THEN
          ILPLOT= 1
          YVAR  = 't mm'
          NCASESX=1
         ELSE IF(OPT.EQ.'2') THEN
          ILPLOT= 2
          YVAR  = 't/c %'
          NCASESX=1
         ELSE IF(OPT.EQ.'3') THEN
          ILPLOT= 3
          YVAR  = 'x/c'
          NCASESX=1
         ELSE IF(OPT.EQ.'4') THEN
          ILPLOT= 4
          YVAR  = 'camb %'
          NCASESX=1
         ELSE IF(OPT.EQ.'5') THEN
          ILPLOT= 5
          YVAR  = 'x/c'
          NCASESX=1
         ELSE IF(OPT.EQ.'6') THEN
          ILPLOT= 6
          YVAR  = 'area scm'
          NCASESX=1
         ELSE IF(OPT.EQ.'7') THEN
          ILPLOT= 7
          YVAR  = 'r_le mm'
          NCASESX=1
         ELSE IF(OPT.EQ.'8') THEN
          ILPLOT= 8
          YVAR  = 't_te mm'
          NCASESX=1
         ELSE IF(OPT.EQ.'9') THEN
          ILPLOT= 9
          YVAR  = 'a_te deg'
          NCASESX=1
         ELSE IF(OPT.EQ.'10') THEN
          ILPLOT= 10
          YVAR  = 'A0 deg'
          NCASESX=2
         ELSE IF(OPT.EQ.'11') THEN
          ILPLOT= 11
          YVAR  = 'beta deg'
          NCASESX=2
         ELSE IF(OPT.EQ.'12') THEN
          ILPLOT= 12
          YVAR  = 'chord cm'
          NCASESX=1
C
         ELSE
            WRITE(*,*)'Command not recognized'
            GO TO 10
        ENDIF
        LYRELIMIT = .TRUE.
        YMINSPEC = 999.
        YMAXSPEC = 999.
C
C
C--- Put selected data into arrays
C
 100  CONTINUE
      TITLE1 = 'ESLOFT Section Data'
      NDT = 0
      NPT = NLOFT2
      DO NS= 1, NCASESX
        NDT = NDT + 1
C
        IFLW = 1
        INFODATA(1,NDT) = IFLW
C
C---- set blade disk to plot
C        N = NR
C        BLDS = FLOAT(NRBLD(N))
C        RPM = 30.0*OMEGA(NR)/PI      
C
        DO I = 1, NPT
C--- Install X data
          IF(XVAR.EQ.'R in') THEN
            XDATA(I,NDT) = YLOFT(I)* MTI
          ELSE
            XDATA(I,NDT) = YLOFT(I)* 100.0
          ENDIF
C
C--- Install Y data-------------------------------------------------
C---- thickness
C
          IF(ILPLOT.EQ.1) THEN
            YDATA(I,NDT) = BLENDATA(4,I)*CDLOFT(I)*1000.0
            LEGENDDATA(NDT)= 'Blade thickness'
C
C---- thickness/chord
C
          ELSEIF(ILPLOT.EQ.2) THEN
            YDATA(I,NDT) = BLENDATA(4,I)*100.0
            LEGENDDATA(NDT)= 'Max thickness/chord'
C
C---- thickness/chord x/c
C
          ELSEIF(ILPLOT.EQ.3) THEN
            YDATA(I,NDT) = BLENDATA(5,I)
            LEGENDDATA(NDT)= 'x/c max_thickness/chord'
C
C---- camber
C
          ELSEIF(ILPLOT.EQ.4) THEN
            YDATA(I,NDT) = BLENDATA(6,I) * 100.0
            LEGENDDATA(NDT)= 'Max camber'
C
C---- camber x/c
C
          ELSEIF(ILPLOT.EQ.5) THEN
            YDATA(I,NDT) = BLENDATA(7,I)
            LEGENDDATA(NDT)= 'x/c max_camber'
C
C---- area
C
          ELSEIF(ILPLOT.EQ.6) THEN
            YDATA(I,NDT)= BLENDATA(1,I)*CDLOFT(I)*CDLOFT(I)*1.0E4
            LEGENDDATA(NDT)= 'Section area'
C
C---- rad le
C
          ELSEIF(ILPLOT.EQ.7) THEN
            YDATA(I,NDT) = BLENDATA(2,I)*CDLOFT(I)*1000.0
            LEGENDDATA(NDT)= 'Leading edge radius'
C
C---- thick te
C
          ELSEIF(ILPLOT.EQ.8) THEN
            YDATA(I,NDT) = BLENDATA(8,I)*CDLOFT(I)*1000.0
            LEGENDDATA(NDT)= 'Trailing edge thickness'
C
C---- angle te
C
          ELSEIF(ILPLOT.EQ.9) THEN
            YDATA(I,NDT) = BLENDATA(3,I)/DTR
            LEGENDDATA(NDT)= 'Trailing edge angle'
C
C---- zero-lift alpha raw and corrected
C
          ELSEIF(ILPLOT.EQ.10) THEN
             IF(NS.EQ.1) THEN
                ATEMP=AZLOFT(I)
                LEGENDDATA(NDT)= 'Zero-lift alpha ESLOFT'
             ELSE
                ATEMP=AZDFDC(I)
                LEGENDDATA(NDT)= 'Zero-lift alpha AERO'
             ENDIF
C
            YDATA(I,NDT) = ATEMP
C
C--- beta raw and corrected in global or local coords
C
          ELSEIF(ILPLOT.EQ.11) THEN
C
             IF(NS.EQ.1) THEN
                BTEMP=BECORR(I)
                LEGENDDATA(NDT)= 'Beta ESLOFT'
             ELSE
                BTEMP=BELOFT(I)
                LEGENDDATA(NDT)= 'Beta DFDC'
             ENDIF
C
             IF(OMEGA(NR).LE.0.0) THEN
                IF(LCOORD) THEN
                   YDATA(I,NDT) = (PI-BTEMP)/DTR
                   TITLE1 = 'ESLOFT Section Data (local coords)'
                ELSE
                   YDATA(I,NDT) = BTEMP/DTR
                   TITLE1 = 'ESLOFT Section Data (global coords)'
                ENDIF
             ELSE
                YDATA(I,NDT) = BTEMP/DTR                
             ENDIF
C
C---- chord
C
          ELSEIF(ILPLOT.EQ.12) THEN
            YDATA(I,NDT) = CDLOFT(I)*100.
            LEGENDDATA(NDT)= 'Blade chord'
C
C
          ENDIF
        END DO
        NPTDATA(NDT) = NPT
C
      END DO
      NDATA = NDT
      GO TO 5
C    
      END     !  PLOTLOFDATA
C
C------------------------------------------------------------------------


