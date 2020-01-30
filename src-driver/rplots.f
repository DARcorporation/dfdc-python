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
C
C=========================================================================
C
C     Version 070-ES1
C     Philip Carter, Esotec Developments, February 2009
C     philip (at) esotec (dot) org
C
C     Changes from 0.70:
C
C     Disk index appended to plot titles for multi-disk cases (GETPNAME).
C
C     Plot 12:
C       Fixed titling bug (TITLE1).
C       Alpha plot added.
C       Beta and Alpha plots can use local or global coordinates (BANG).
C       Cleaned up legend.
C
C     Version 070-ES1a  4 March 2009
C
C     Subroutine GETROTG
C     Splines chord and beta to panel boundaries (CHRP, BETARP)
C
C     BLDPLT, ROTRPLT
C     Blade geometry plots show true geometry
C
C     BLDPLT, ROTRPLT, BLDCLPLT, ROTVPLT, PLOTBLDDATA
C     Numeric output fixed and updated, correct chords and betas
C     Tip gap properly handled (last point omitted)     
C
C=========================================================================


      SUBROUTINE BLDPLT(NR,VIEW)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      EXTERNAL PLCHAR,PLMATH
C
      CHARACTER*(*) VIEW
      CHARACTER*80  PNAME
C------------------------------------------
C     Plots blade planform and projections
C------------------------------------------
      DIMENSION XLIN(4), YLIN(4)
      DIMENSION W1(IRX), W2(IRX), W3(IRX), 
     &          W4(IRX), W5(IRX), W6(IRX) 
      DIMENSION XI(IRX)
C
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
C---- character size for axis numbers, labels
      CS  = 1.1*CSIZE
      CSL = 1.4*CSIZE
C
      IF(NROTOR.LE.1) THEN
         PNAME=NAME
      ELSE
         CALL GETPNAME(PNAME,NAME,NR)
      ENDIF
C
C---- case title
      CALL NEWPEN(2)
      XT = 0.5
      YT = -0.030
      CALL PLCHAR(XT,YT,CSL,PNAME,0.0,-1)
      CALL PLOT(0.0,0.15,-3)
C
C---- vertical space/R between projections
      YSPACE = 0.375
C
C---- set blade disk to plot
      N = NR
C
C---- modified input - geometery taken at panel boundaries
C----------------------------------------------------------------
C---- get physical blade geometry
      CALL GETROTG(N)
C
C---- mean radius
      YMNR = SQRT(0.5*ADISK(N)/PI + RHUB(N)**2)
C
C---- blade solidity
      CALL SPLINE(CHRP(1,N),W1,YRP,NRP)
      CMNR = SEVAL(YMNR,CHRP(1,N),W1,YRP,NRP)
      SIGMR = FLOAT(NRBLD(N))*CMNR/(2.0*PI*YMNR)
C
C---- blade angles and twist
      BHUB = BETARP(1,N)
      BTIP = BETARP(NRP,N)
C      B2   = BETAR(2,N)
C      BEXT = BHUB - YRC(1,N)*(BHUB-B2)/(YRC(1,N)-YRC(2,N))
      BTWIST = BHUB - BTIP
      CALL SPLINE(BETARP(1,N),W2,YRP,NRP)
      BETAMR = SEVAL(YMNR,BETARP(1,N),W2,YRP,NRP)
C----------------------------------------------------------------
C
C      CALL SPLINE(CHR(1,N),W1,YRC,NRC)
C      CMNR = SEVAL(YMNR,CHR(1,N),W1,YRC,NRC)
C      SIGMR = FLOAT(NRBLD(N))*CMNR/(2.0*PI*YMNR)
C
C---- blade angles and twist
C      BHUB = BETAR(1,N)
C      BTIP = BETAR(NRC,N)
C      B2   = BETAR(2,N)
C      BEXT = BHUB - YRC(1,N)*(BHUB-B2)/(YRC(1,N)-YRC(2,N))
C      BTWIST = BHUB - BTIP
C      CALL SPLINE(BETAR(1,N),W2,YRC,NRC)
C      BETAMR = SEVAL(YMNR,BETAR(1,N),W2,YRC,NRC)
C------------------------
C
cc      CALL GETWINSIZE(XWIND,YWIND)
cc      write(*,*) 'xwind,ywind ',xwind,ywind
cc      CALL GETFACTORS(XSCALE,YSCALE)
cc      write(*,*) 'xscale,yscale ',xscale,yscale
C
      IVIEW = INDEX(VIEW,'E') + INDEX(VIEW,'e')
C
      IF(IVIEW.NE.0) THEN
       DXEND = 1.02
       CALL PLOT(DXEND,0.1*YSPACE,-3)
C
       CALL GETCOLOR(ICOL0)
       CALL NEWPEN(2)
       XL1 = 0.0
       XL2 = XL1 + 16.0*CS
       XL3 = XL2 + 16.0*CS
       XL4 = XL3 + 16.0*CS
       XL5 = XL4 + 16.0*CS
C
       YL = 0.0
       CALL PLCHAR(XL1 ,YL,CS,'blds = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, FLOAT(NRBLD(N)),0.0,-1)
       CALL PLSUBS(XL2 ,YL,0.8*CS,'mean' ,0.0,4,PLCHAR)
       CALL PLMATH(XL2 ,YL,CS,'s    = '  ,0.0,7)
       CALL PLNUMB(999.,YL,CS, SIGMR     ,0.0,4)
C
       YL = YL - 2.5*CS
       CALL PLCHAR(XL1 ,YL,CS,'Rtip = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, RTIP(N)    ,0.0,4)
       CALL PLCHAR(XL2 ,YL,CS,'Ctip = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, CHRP(NRP,N) ,0.0,4)
       CALL PLSUBS(XL3 ,YL,CS,'tip'       ,0.0,3,PLCHAR)
       CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, BETARP(NRP,N)/DTR,0.0,2)
C
       YL = YL - 2.5*CS
       CALL PLCHAR(XL1 ,YL,CS,'Rhub = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, RHUB(N)    ,0.0,4)
       CALL PLCHAR(XL2 ,YL,CS,'Chub = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, CHRP(1,N)   ,0.0,4)
       CALL PLSUBS(XL3 ,YL,CS,'hub'       ,0.0,3,PLCHAR)
       CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, BETARP(1,N)/DTR ,0.0,2)
C
       YL = YL - 2.5*CS
       CALL PLCHAR(XL1 ,YL,CS,'Rmean= '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, YMNR       ,0.0,4)
       CALL PLCHAR(XL2 ,YL,CS,'Cmean= '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, CMNR       ,0.0,4)
       CALL PLSUBS(XL3 ,YL,CS,'mean'      ,0.0,4,PLCHAR)
       CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, BETAMR/DTR ,0.0,2)
C
       YL = YL - 2.5*CS
       CALL PLCHAR(XL1 ,YL,CS,'Adisk= '  ,0.0,7)
       CALL PLNUMB(999.,YL,CS, ADISK(N)  ,0.0,4)
       CALL PLCHAR(XL2 ,YL,CS,'Atip = '  ,0.0,7)
       CALL PLNUMB(999.,YL,CS, ATIP(N)   ,0.0,4)
       CALL PLSUBS(XL3 ,YL,CS,'twist'     ,0.0,5,PLCHAR)
       CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, BTWIST/DTR ,0.0,2)
C
       CALL PLOT(-DXEND,-0.5*YSPACE,-3)
       CALL NEWCOLOR(ICOL0)
      ENDIF
C
C
      XIHUB = YRP(1,N)/RTIP(N)
      CHUMAX = 0.
      CHLMAX = 0.
C
      IF(TGAP.GT.0.0 .AND. OMEGA(N).NE.0.0) THEN
         NRL = NRP-1
      ELSE
         NRL = NRP
      ENDIF
C
      DO 10 I=1, NRL
        SINB = SIN(BETARP(I,N))
        COSB = COS(BETARP(I,N))
C
        XI(I) =  YRP(I,N)/RTIP(N)
C------ projection chords in front and behind radial axis
        W1(I) =  XPAXIS     *SINB*CHRP(I,N)/RTIP(N)
        W2(I) = (XPAXIS-1.0)*SINB*CHRP(I,N)/RTIP(N)
C
        W3(I) =  XPAXIS     *COSB*CHRP(I,N)/RTIP(N)
        W4(I) = (XPAXIS-1.0)*COSB*CHRP(I,N)/RTIP(N)
C
        W5(I) =  XPAXIS     *CHRP(I,N)/RTIP(N)
        W6(I) = (XPAXIS-1.0)*CHRP(I,N)/RTIP(N)
C
        CHUMAX = MAX(CHUMAX,     CHRP(I,N)/RTIP(N))
        CHLMAX = MAX(CHLMAX,SINB*CHRP(I,N)/RTIP(N))
   10 CONTINUE
C
C---- add vertical space for wide chords
ccc      YSPACE = MAX( YSPACE , 1.05*CHUMAX )
C
C---- plot untwisted view if selected
      IVIEW = INDEX(VIEW,'U') + INDEX(VIEW,'u')
      IF(IVIEW.GT.0) THEN
        CALL PLOT(0.0,(FLOAT(IVIEW)-0.5)*YSPACE,-3)
        CALL NEWPEN(1)
        CALL GETCOLOR(ICOL0)
        CALL NEWCOLORNAME('yellow')
C
C------ plot radial axis
        CALL PLOT(0.0,0.0,3)
        CALL PLOT(1.0,0.0,2)
C
C------ plot radial tick marks
        DO NT=0, 5
          XT = 0.2*FLOAT(NT)
          CALL PLOT(XT,-.005,3)
          CALL PLOT(XT,0.005,2)
        ENDDO
C
C------ plot blade shape
        CALL NEWPEN(4)
        CALL XYLINE(NRL,XI,W5,0.0,1.0,0.0,1.0,1)
        CALL XYLINE(NRL,XI,W6,0.0,1.0,0.0,1.0,1)
C
C------ plot hub line
        CALL NEWPEN(2)
        CALL PLOT(XIHUB,CHUMAX * XPAXIS     ,3)
        CALL PLOT(XIHUB,CHUMAX *(XPAXIS-1.0),2)
C        CALL PLOT(XIHUB,CHRP(1,N)*(XPAXIS-1.0),2)
C
C---- label for plot view
        XT = XI(1) + 1.0*CS
        YT = -0.35*YSPACE
        CALL PLCHAR(XT,YT,CS,'Untwisted',0.0,-1)
C
        CALL PLOT(0.0,-(FLOAT(IVIEW)-0.5)*YSPACE,-3)
        CALL NEWCOLOR(ICOL0)
      ENDIF
C
C
C---- plot lateral view if selected
      IVIEW = INDEX(VIEW,'L') + INDEX(VIEW,'l')
      IF(IVIEW.GT.0) THEN
        CALL PLOT(0.0,(FLOAT(IVIEW)-0.5)*YSPACE,-3)
        CALL NEWPEN(1)
        CALL GETCOLOR(ICOL0)
        CALL NEWCOLORNAME('orange')
C
C------ plot radial axis
        CALL PLOT(0.0,0.0,3)
        CALL PLOT(1.0,0.0,2)
C
C------ plot rotation axis
        CALL PLOT(0.0, 0.1,3)
        CALL PLOT(0.0,-0.2,2)
C
C------ plot radial tick marks
        DO NT=0, 5
          XT = 0.2*FLOAT(NT)
          CALL PLOT(XT,-.005,3)
          CALL PLOT(XT,0.005,2)
        ENDDO
C
C------ plot blade shape
        CALL NEWPEN(4)
        CALL XYLINE(NRL,XI,W1,0.0,1.0,0.0,1.0,1)
        CALL XYLINE(NRL,XI,W2,0.0,1.0,0.0,1.0,1)
C
C------ plot hub surface
        CALL NEWPEN(2)
        XLIN(1) = XIHUB
        XLIN(2) = XIHUB
        YLIN(1) = CHLMAX* XPAXIS     
        YLIN(2) = CHLMAX*(XPAXIS-1.0)
C
        CALL PLOT(XLIN(1),YLIN(1),3)
        CALL PLOT(XLIN(2),YLIN(2),2)
C
cc        CALL XYLINE(2,XLIN(2),YLIN(2),0.0,1.0,0.0,1.0,4)
C
C---- label for plot view
        XT = XI(1) + 1.0*CS
        YT = -0.35*YSPACE
        CALL PLCHAR(XT,YT,CS,'Lateral',0.0,-1)
C
        CALL PLOT(0.0,-(FLOAT(IVIEW)-0.5)*YSPACE,-3)
        CALL NEWCOLOR(ICOL0)
      ENDIF
C
C
C---- plot axial view if selected
      IVIEW = INDEX(VIEW,'A') + INDEX(VIEW,'a')
      IF(IVIEW.GT.0) THEN
        CALL PLOT(0.0,(FLOAT(IVIEW)-0.5)*YSPACE,-3)
        CALL NEWPEN(1)
        CALL GETCOLOR(ICOL0)
        CALL NEWCOLORNAME('green')
C
C------ plot radial axis
        CALL PLOT(0.0,0.0,3)
        CALL PLOT(1.0,0.0,2)
C
C------ for axial view, plot blade "spokes" (A)
        RSPOKE = MAX( 0.1 , 1.25*XIHUB )
        DO IB=2, NRBLD(N)
          ANG = 2.0*PI * FLOAT(IB-1)/FLOAT(NRBLD(N))
          XT = RSPOKE*COS(ANG)
          YT = RSPOKE*SIN(ANG)
          CALL PLOT(0.0,0.0,3)
          CALL PLOT(XT,YT,2)
        ENDDO
C
C------ plot radial tick marks
        DO NT=0, 5
          XT = 0.2*FLOAT(NT)
          CALL PLOT(XT,-.005,3)
          CALL PLOT(XT,0.005,2)
        ENDDO
C
C------ plot hub circle
        CALL NEWPEN(2)
        CALL PLCIRC(0.0,0.0,XIHUB,0)
C
C------ plot blade shape
        CALL NEWPEN(4)
        CALL XYLINE(NRL,XI,W3,0.0,1.0,0.0,1.0,1)
        CALL XYLINE(NRL,XI,W4,0.0,1.0,0.0,1.0,1)
C
C---- label for plot view
        XT = XI(1) + 1.0*CS
        YT = -0.35*YSPACE
        CALL PLCHAR(XT,YT,CS,'Axial',0.0,-1)
C
        CALL PLOT(0.0,-(FLOAT(IVIEW)-0.5)*YSPACE,-3)
        CALL NEWCOLOR(ICOL0)
      ENDIF
C
C
C---- plot end view if selected
      KVIEW = INDEX(VIEW,'E') + INDEX(VIEW,'e')
      IVIEW = INDEX(VIEW,'L') + INDEX(VIEW,'l')
      IF(KVIEW.GT.0 .AND. IVIEW.GT.0) THEN
        CALL GETCOLOR(ICOL0)
        CALL NEWCOLORNAME('red')
C
        DXEND = 1.25
        CALL PLOT(DXEND,(FLOAT(IVIEW)-0.5)*YSPACE,-3)
C
C------ plot blade shape
        CALL NEWPEN(4)
        CALL XYLINE(NRL,W3,W1,0.0,1.0,0.0,1.0,1)
        CALL XYLINE(NRL,W4,W2,0.0,1.0,0.0,1.0,1)
C
        I = 1
        CALL PLOT(W3(I),W1(I),3)
        CALL PLOT(W4(I),W2(I),2)
C
        I = NRL
        CALL PLOT(W3(I),W1(I),3)
        CALL PLOT(W4(I),W2(I),2)
C
C------ plot chord lines
        CALL NEWPEN(1)
        DO NT=1, 5
          XT = 0.2*FLOAT(NT)
          DO I=1, NRL-1
            IF(XT.GT.XI(I) .AND. XT.LE.XI(I+1)) THEN
             FRAC = (XT-XI(I)) / (XI(I+1)-XI(I))
             W1T = W1(I) + FRAC*(W1(I+1)-W1(I))
             W2T = W2(I) + FRAC*(W2(I+1)-W2(I))
             W3T = W3(I) + FRAC*(W3(I+1)-W3(I))
             W4T = W4(I) + FRAC*(W4(I+1)-W4(I))
             CALL PLOT(W3T,W1T,3)
             CALL PLOT(W4T,W2T,2)
            ENDIF
          ENDDO
        ENDDO
C
C------ plot rotation axis
        CALL PLOT(0.0, 0.1,3)
        CALL PLOT(0.0,-0.2,2)
C
C------ plot hub surface
        CALL NEWPEN(2)
        XLIN(1) = XIHUB
        XLIN(2) = XIHUB
        YLIN(1) = CHLMAX* XPAXIS     
        YLIN(2) = CHLMAX*(XPAXIS-1.0)
C
        CALL PLOT( XLIN(1),YLIN(1),3)
        CALL PLOT( XLIN(2),YLIN(2),2)
        CALL PLOT(-XLIN(1),YLIN(1),3)
        CALL PLOT(-XLIN(2),YLIN(2),2)
C
cc        CALL XYLINE(2,XLIN(2),YLIN(2),0.0, 1.0,0.0,1.0,4)
cc        CALL XYLINE(2,XLIN(2),YLIN(2),0.0,-1.0,0.0,1.0,4)
C
C------ plot prop-plane line
        DS = 0.2
        CALL NEWPEN(1)
        CALL NEWPAT(LMASK2)
        CALL PLOT(-DS,0.0,3)
        CALL PLOT( DS,0.0,2)
C
        CALL NEWPAT(LMASK1)
C
C------ plot hub lines
        DX = DS*COS(BETARP(1,N))
        DY = DS*SIN(BETARP(1,N))
        CALL PLOT(-DX,-DY,3)
        CALL PLOT( DX, DY,2)
C
C------ and tip-angle lines
        DX = DS*COS(BETARP(NRL,N))
        DY = DS*SIN(BETARP(NRL,N))
        CALL PLOT(-DX,-DY,3)
        CALL PLOT( DX, DY,2)
C
C---- label for plot view
        XT = -1.5*CS
        YT = -0.2 - 2.0*CS
        CALL PLCHAR(XT,YT,CS,'End',0.0,-1)
C
        CALL PLOT(-DXEND,-(FLOAT(IVIEW)-0.5)*YSPACE,-3)
        CALL NEWCOLOR(ICOL0)
      ENDIF
C
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! BLDPLT



      SUBROUTINE ROTRPLT(NR)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C-------------------------------------------
C     Plots rotor disk outlines (all blades)
C-------------------------------------------
      CHARACTER*80 PNAME
      EXTERNAL PLCHAR,PLMATH
      DIMENSION XLIN(4), YLIN(4)
      DIMENSION W1(IRX), W2(IRX), W3(IRX), 
     &          W4(IRX), W5(IRX), W6(IRX) 
      DIMENSION XI(IRX)
C
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
C---- set blade disk to plot
      N = NR
C
C----------------------------------------------------------------
C---- get physical blade geometry
      CALL GETROTG(N)
C
C---- mean radius
      YMNR = SQRT(0.5*ADISK(N)/PI + RHUB(N)**2)
C
C---- blade solidity
      CALL SPLINE(CHRP(1,N),W1,YRP,NRP)
      CMNR = SEVAL(YMNR,CHRP(1,N),W1,YRP,NRP)
      SIGMR = FLOAT(NRBLD(N))*CMNR/(2.0*PI*YMNR)
C
C---- blade angles and twist
      BHUB = BETARP(1,N)
      BTIP = BETARP(NRP,N)
C      B2   = BETAR(2,N)
C      BEXT = BHUB - YRC(1,N)*(BHUB-B2)/(YRC(1,N)-YRC(2,N))
      BTWIST = BHUB - BTIP
      CALL SPLINE(BETARP(1,N),W2,YRP,NRP)
      BETAMR = SEVAL(YMNR,BETARP(1,N),W2,YRP,NRP)
C----------------------------------------------------------------
C
C---- character size for axis numbers, labels
      CS  = 0.7*CSIZE
      CSL = CS*1.4
C
cc      CALL GETWINSIZE(XWIND,YWIND)
cc      write(*,*) 'xwind,ywind ',xwind,ywind
cc      CALL GETFACTORS(XSCALE,YSCALE)
cc      write(*,*) 'xscale,yscale ',xscale,yscale
C
      IF(NROTOR.LE.1) THEN
         PNAME=NAME
      ELSE
         CALL GETPNAME(PNAME,NAME,NR)
      ENDIF
C
C---- case title
      CALL PLOTABS(1.0,1.2,-3)
      CALL NEWPEN(4)
      XL1 = 0.0
      YL = 5.0*CS
      CALL PLCHAR(XL1,YL,CSL,PNAME,0.0,-1)
C
      CALL NEWPEN(3)
      XL1 = 0.0
      XL2 = XL1 + 17.0*CS
      XL3 = XL2 + 17.0*CS
      XL4 = XL3 + 17.0*CS
      XL5 = XL4 + 17.0*CS
C
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'blds = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, FLOAT(NRBLD(N)),0.0,-1)
      CALL PLSUBS(XL2 ,YL,0.8*CS,'mean' ,0.0,4,PLCHAR)
      CALL PLMATH(XL2 ,YL,CS,'s    = '  ,0.0,7)
      CALL PLNUMB(999.,YL,CS, SIGMR     ,0.0,4)
C
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'Rtip = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, RTIP(N)    ,0.0,4)
      CALL PLCHAR(XL2 ,YL,CS,'Ctip = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CHRP(NRP,N) ,0.0,4)
      CALL PLSUBS(XL3 ,YL,CS,'tip'       ,0.0,3,PLCHAR)
      CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BETARP(NRP,N)/DTR,0.0,2)
C
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'Rhub = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, RHUB(N)    ,0.0,4)
      CALL PLCHAR(XL2 ,YL,CS,'Chub = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CHRP(1,N)   ,0.0,4)
      CALL PLSUBS(XL3 ,YL,CS,'hub'       ,0.0,3,PLCHAR)
      CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BETARP(1,N)/DTR ,0.0,2)
C
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'Rmean= '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, YMNR       ,0.0,4)
      CALL PLCHAR(XL2 ,YL,CS,'Cmean= '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CMNR       ,0.0,4)
      CALL PLSUBS(XL3 ,YL,CS,'mean'      ,0.0,4,PLCHAR)
      CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BETAMR/DTR ,0.0,2)
C
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'Adisk= '  ,0.0,7)
      CALL PLNUMB(999.,YL,CS, ADISK(N)  ,0.0,4)
      CALL PLCHAR(XL2 ,YL,CS,'Atip = '  ,0.0,7)
      CALL PLNUMB(999.,YL,CS, ATIP(N)   ,0.0,4)
      CALL PLSUBS(XL3 ,YL,CS,'twist'     ,0.0,5,PLCHAR)
      CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BTWIST/DTR ,0.0,2)
C
C
cc      CALL PLSUBS(XL2 ,YL,0.8*CS,'wak' ,0.0,3,PLCHAR)
cc      CALL PLCHAR(XL2 ,YL,CS,'R    = ' ,0.0,7)
cc      CALL PLNUMB(999.,YL,CS, XW0*RAD  ,0.0,4)
C
C---- vertical space/R between projections
C
      XIHUB = YRP(1,N)/RTIP(N)
      CHUMAX = 0.
      CHLMAX = 0.
C
      IF(TGAP.GT.0.0 .AND. OMEGA(N).NE.0.0) THEN
         NRL = NRP-1
      ELSE
         NRL = NRP
      ENDIF
C
      DO I=1, NRL
        SINB = SIN(BETARP(I,N))
        COSB = COS(BETARP(I,N))
C
        XI(I) =  YRP(I,N)/RTIP(N)
C------ projection chords in front and behind radial axis
        W1(I) =  XPAXIS     *SINB*CHRP(I,N)/RTIP(N)
        W2(I) = (XPAXIS-1.0)*SINB*CHRP(I,N)/RTIP(N)
C
        W3(I) =  XPAXIS     *COSB*CHRP(I,N)/RTIP(N)
        W4(I) = (XPAXIS-1.0)*COSB*CHRP(I,N)/RTIP(N)
C
        W5(I) =  XPAXIS     *CHRP(I,N)/RTIP(N)
        W6(I) = (XPAXIS-1.0)*CHRP(I,N)/RTIP(N)
C
        CHUMAX = MAX(CHUMAX,     CHRP(I,N)/RTIP(N))
        CHLMAX = MAX(CHLMAX,SINB*CHRP(I,N)/RTIP(N))
      END DO
C
C
C---- plot axial view 
cc        CALL PLOT(0.5,0.75,-3)
        CALL PLOTABS(0.5*XWIND,0.5*XWIND+1.2,-3)
        CALL NEWPEN(1)
        CALL GETFACTORS(XSCALE,YSCALE)
        SCL = 0.8*(0.5*XWIND)
        CALL NEWFACTOR(SCL)
C
C------ plot radial axis
        CALL PLOT(0.0,0.0,3)
        CALL PLOT(1.0,0.0,2)
C
C------ for axial view, plot blade "spokes" (A)
        RSPOKE = MAX( 0.1 , 1.25*XIHUB )
        DO IB=2, NRBLD(N)
          ANG = 2.0*PI * FLOAT(IB-1)/FLOAT(NRBLD(N))
          XT = RSPOKE*COS(ANG)
          YT = RSPOKE*SIN(ANG)
          CALL PLOT(0.0,0.0,3)
          CALL PLOT(XT,YT,2)
        ENDDO
C
C------ plot radial tick marks on #1 (horizontal) blade
        DO NT=0, 5
          XT = 0.2*FLOAT(NT)
          CALL PLOT(XT,-.005,3)
          CALL PLOT(XT,0.005,2)
        ENDDO
C
C------ plot tip circle for duct
          CALL NEWPEN(2)
          CALL PLCIRC(0.0,0.0,1.0,72)
C
C------ plot hub circle
        CALL NEWPEN(2)
        CALL PLCIRC(0.0,0.0,XIHUB,0)
C
C------ plot blade shape
        CALL GETCOLOR(ICOL0)
        CALL NEWCOLORNAME('red')
        CALL NEWPEN(3)
        DO IB=1, NRBLD(N)
          ANG = 2.0*PI * FLOAT(IB-1)/FLOAT(NRBLD(N))
          SINA = SIN(ANG)
          COSA = COS(ANG)
          DO I = 1, NRL
c--- chordline defined on circular cylinder at radius XI(I)
             THET = W3(I)/XI(I)
             XX = XI(I)*COS(THET)
             YY = XI(I)*SIN(THET)
             XT = XX*COSA - YY*SINA
             YT = XX*SINA + YY*COSA
c--- chordline defined on plane normal to stacking line
c            XT = XI(I)*COSA - W3(I)*SINA
c            YT = XI(I)*SINA + W3(I)*COSA
            IF(I.EQ.1) THEN
              CALL PLOT(XT,YT,3)
             ELSE
              CALL PLOT(XT,YT,2)
            ENDIF
          ENDDO
          DO I = 1, NRL
c--- chordline defined on circular cylinder at radius XI
             THET = W4(I)/XI(I)
             XX = XI(I)*COS(THET)
             YY = XI(I)*SIN(THET)
             XT = XX*COSA - YY*SINA
             YT = XX*SINA + YY*COSA
c--- chordline defined on plane normal to stacking line
c            XT = XI(I)*COSA - W4(I)*SINA
c            YT = XI(I)*SINA + W4(I)*COSA
            IF(I.EQ.1) THEN
              CALL PLOT(XT,YT,3)
             ELSE
              CALL PLOT(XT,YT,2)
            ENDIF
          ENDDO
        ENDDO
cc        CALL XYLINE(NRC,XI,W3,0.0,1.0,0.0,1.0,1)
cc        CALL XYLINE(NRC,XI,W4,0.0,1.0,0.0,1.0,1)
        CALL NEWCOLOR(ICOL0)
C
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END ! PRPPLT



      SUBROUTINE BLDCLPLT(NR)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      DIMENSION XI(IRX), GAMB(IRX)
      DIMENSION W1(IRX), W2(IRX)
C---------------------------------------------------
C     Plots CL, eta, normalized Gamma vs r/R
C---------------------------------------------------
      LOGICAL LMAPLT
      CHARACTER*80 PNAME
C
      EXTERNAL PLCHAR,PLMATH
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
      IF(NROTOR.LE.1) THEN
         PNAME=NAME
      ELSE
         CALL GETPNAME(PNAME,NAME,NR)
      ENDIF
C
C---- set blade disk to plot
      N = NR
C
C----------------------------------------------------------------
C---- get physical blade geometry
      CALL GETROTG(N)
C
C---- mean radius
      YMNR = SQRT(0.5*ADISK(N)/PI + RHUB(N)**2)
C
C---- blade solidity
      CALL SPLINE(CHRP(1,N),W1,YRP,NRP)
      CMNR = SEVAL(YMNR,CHRP(1,N),W1,YRP,NRP)
      SIGMR = FLOAT(NRBLD(N))*CMNR/(2.0*PI*YMNR)
C
C---- blade angles and twist
      BHUB = BETARP(1,N)
      BTIP = BETARP(NRP,N)
C      B2   = BETAR(2,N)
C      BEXT = BHUB - YRC(1,N)*(BHUB-B2)/(YRC(1,N)-YRC(2,N))
      BTWIST = BHUB - BTIP
      CALL SPLINE(BETARP(1,N),W2,YRP,NRP)
      BETAMR = SEVAL(YMNR,BETARP(1,N),W2,YRP,NRP)
C----------------------------------------------------------------
C
C---- plot aspect ratio
      PLPAR = 0.6
C
C---- character size for axis numbers, labels
      CS  = 0.9*CSIZE
      CSL = CS*1.4
C
C---- system calcs
C
      RPM = 30.0*OMEGA(N)/PI
      IF(OMEGA(N).NE.0.0) THEN
        ADV = QINF/(OMEGA(N)*RTIP(N))
      ELSE
        ADV = 0.0
      ENDIF
C
C---- dimensional thrust, power, torque, rpm
      TDIM = TTOT + TDUCT
      QDIM = QTOT
      PDIM = PTOT
      TVDIM = TVIS
      PVDIM = PVIS
      TINV  = TTOT - TVIS 
      QINV  = QTOT - QVIS 
      PINV  = PTOT - PVIS 
C
C---- Thrust and power coefficients based on rotational speed
      DIA = 2.0*RTIP(N)
      EN = ABS(RPM/60.0)
      IF(EN.EQ.0.0) THEN
        CT = 0.0
        CP = 0.0
      ELSE
        CT = TDIM/(RHO*EN**2*DIA**4)
        CP = PDIM/(RHO*EN**3*DIA**5)
      ENDIF
C
C---- Thrust and power coefficients based on forward speed
      IF(QINF.GT.0.0) THEN
        TC = TDIM/(0.5*RHO*QINF**2 * PI*RTIP(N)**2)
        PC = PDIM/(0.5*RHO*QINF**3 * PI*RTIP(N)**2)
      ELSE
        TC = 0.0
        PC = 0.0
      ENDIF
C---- Thrust and power coefficients based on tip speed, disk area
      VTIP = ABS(OMEGA(N)*RTIP(N))
      IF(VTIP.NE.0.0) THEN
       CT0  = TDIM/(RHO*ADISK(N)*VTIP**2)
       CP0  = PDIM/(RHO*ADISK(N)*VTIP**3)
       CTOS = CT0 / SIGMR
       IF(CP0.NE.0.0 .AND. CT0.GE.0.0) THEN
         FOM = ABS(CT0)**1.5 / CP0 / 2.0
       ELSE
         FOM = 0.0
       ENDIF
      ELSE
       CT0  = 0.0
       CP0  = 0.0
       CTOS = 0.0
       FOM  = 0.0
      ENDIF
C
C---- overall efficiency
      IF(PDIM.NE.0.0) THEN
        EFFTOT = QINF*TDIM/PDIM
       ELSE
        EFFTOT = 0.0
      ENDIF
C---- induced efficiency (including nacelle thrust effect)
      IF(PINV.NE.0.0) THEN
        EFFIND = QINF*(TINV+TDUCT)/PINV
      ELSE
        EFFIND = 0.0
      ENDIF
C---- ideal (actuator disk) efficiency
      IF(TC.EQ.0.0) THEN
        EIDEAL = 0.0
      ELSE
        TCLIM = MAX( -1.0 , TC )
        EIDEAL = 2.0 / (1.0 + SQRT(TCLIM + 1.0))
      ENDIF
C
C
C---- rotor calcs
C---- dimensional thrust, power, torque, rpm
      TDIM2  = TTOTR(N)
      QDIM2  = QTOTR(N)
      PDIM2  = QDIM2    * OMEGA(N)
      TVDIM2 = TVISR(N)
      PVDIM2 = QVISR(N) * OMEGA(N)
      TINV2  = TINVR(N)
      QINV2  = QINVR(N)
      PINV2  = QINV2    * OMEGA(N)
C
C
C---- Define reference quantities for coefficients
      RREF = RTIP(N)
      AREF = ADISK(N)
      OMEGREF = OMEGA(N)
C
C---- standard thrust/power coefficients based on rotational speed
C
      IF(EN.EQ.0.0) THEN
        CT2 = 0.0
        CP2 = 0.0
      ELSE
        CT2 = TDIM2/(RHO*EN**2*DIA**4)
        CP2 = PDIM2/(RHO*EN**3*DIA**5)
      ENDIF
C
C---- standard thrust/power coefficients based on forward speed
      IF(QINF.GT.0.0) THEN
        TC2 = TDIM2/(0.5*RHO*QINF**2 * PI*RREF**2)
        PC2 = PDIM2/(0.5*RHO*QINF**3 * PI*RREF**2)
      ELSE
        TC2 = 0.0
        PC2 = 0.0
      ENDIF
C
C---- efficiency for rotor 
      IF(PDIM2.NE.0.0) THEN
         EFFTOT2 = QINF*TDIM2/PDIM2
      ELSE
         EFFTOT2 = 0.0
      ENDIF
C
C---------------------------------------------------------------
C
      CALL GETCOLOR(ICOL0)
C
      CALL NEWPEN(2)
      CALL XAXIS(0.0,0.0, 1.0,0.2    , 0.0,0.2, CS,1)
      CALL NEWPEN(3)
      CALL PLCHAR(0.7-1.5*CSL,-3.0*CSL,CSL,'r/R',0.0,3)
C
      CALL NEWPEN(1)
      NXG = 5
      NYG = 5
      CALL PLGRID(0.0,0.0, NXG,0.2, NYG,0.2*PLPAR, LMASK2 )
C
C---- set label locations
      YL = PLPAR + 15.5*CSL
      XL1 = 0.0
      XL2 = XL1 + 18.0*CS
      XL3 = XL2 + 18.0*CS
      XL4 = XL3 + 18.0*CS
      XL5 = XL4 + 18.0*CS
C
      CALL NEWPEN(4)
      CALL PLCHAR(XL1,YL,CSL,PNAME,0.0,-1)
C
      YL = YL - 1.0*CS
      CALL NEWPEN(1)
      CALL PLOT(0.0,YL,3)
      CALL PLOT(1.0,YL,2)
C
      CALL NEWPEN(3)
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'blds = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, FLOAT(NRBLD(N)),0.0,-1)
      CALL PLCHAR(XL2 ,YL,CS,'R(m) = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, RTIP(N)    ,0.0,3)
      CALL PLCHAR(XL3 ,YL,CS,'Adisk= '  ,0.0,7)
      CALL PLNUMB(999.,YL,CS, ADISK(N)  ,0.0,4)
      CALL PLSUBS(XL4 ,YL,0.8*CS,'mean',0.0,4,PLCHAR)
      CALL PLMATH(XL4 ,YL,CS,'s    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, SIGMR   ,0.0,4)
      CALL PLSUBS(XL5 ,YL,CS,'twist'  ,0.0,5,PLCHAR)
      CALL PLMATH(XL5 ,YL,CS,'b    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BTWIST/DTR  ,0.0,3)
C
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'Vm/s = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, QINF    ,0.0,3)
      CALL PLCHAR(XL2 ,YL,CS,'RPM  = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, RPM     ,0.0,1)
      CALL PLCHAR(XL3, YL,CS,'   R   ' ,0.0,7)
      CALL PLMATH(XL3, YL,CS,'  W    ' ,0.0,7)
      CALL PLCHAR(XL3, YL,CS,' /     ' ,0.0,7)
      CALL PLCHAR(XL3, YL,CS,'V    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, ADV     ,0.0,4)
      CALL PLCHAR(XL4 ,YL,CS,'h(km)= ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS,ALTH,0.0,3)
c      CALL PLCHAR(XL4 ,YL,CS,'hKft= ' ,0.0,6)
c      CALL PLNUMB(999.,YL,CS,3.281*ALT,0.0,0)
      CALL PLSUBS(XL5 ,YL,CS,'tip'    ,0.0,3,PLCHAR)
      CALL PLMATH(XL5 ,YL,CS,'b    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BTIP/DTR,0.0,3)
C
C
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'Vavg = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, VAAVG(N)   ,0.0,3)
      CALL PLCHAR(XL2 ,YL,CS,'T(N) = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, TDIM2      ,0.0,2)
      CALL PLCHAR(XL3 ,YL,CS,'P(kW)= '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, PDIM2/1000.,0.0,3)
      CALL PLCHAR(XL4 ,YL,CS,'Q(Nm)= '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, QDIM2      ,0.0,2)
      CALL PLSUBS(XL5 ,YL,CS,'hub'       ,0.0,3,PLCHAR)
      CALL PLMATH(XL5 ,YL,CS,'b    = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BHUB/DTR   ,0.0,3)
C
      YL = YL - 2.5*CS
      CALL PLSUBS(XL1 ,YL,CS,'C'       ,0.0,1,PLCHAR)
      CALL PLCHAR(XL1 ,YL,CS,'T    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, TC2      ,0.0,4)
      CALL PLSUBS(XL2 ,YL,CS,'C'       ,0.0,1,PLCHAR)
      CALL PLCHAR(XL2 ,YL,CS,'P    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, PC2      ,0.0,4)
      CALL PLSUBS(XL3 ,YL,CS,'T'       ,0.0,1,PLCHAR)
      CALL PLCHAR(XL3 ,YL,CS,'C    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CT2      ,0.0,4)
      CALL PLSUBS(XL4 ,YL,CS,'P'       ,0.0,1,PLCHAR)
      CALL PLCHAR(XL4 ,YL,CS,'C    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CP2      ,0.0,4)
      CALL PLMATH(XL5 ,YL,CS,'h    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, EFFTOT2  ,0.0,4)
C
      YL = YL - 1.0*CS
      CALL NEWPEN(1)
      CALL PLOT(0.0,YL,3)
      CALL PLOT(1.0,YL,2)
      CALL NEWPEN(3)
C
      YL = YL - 2.5*CS
c      CALL PLCHAR(XL1 ,YL,CS,'T lb= ' ,0.0,6)
c      CALL PLNUMB(999.,YL,CS, TDIM/4.45,0.0,1)
      CALL PLCHAR(XL1 ,YL,CS,'Ttot = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, TDIM,0.0,2)
      CALL PLCHAR(XL2 ,YL,CS,'Tduct= ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, TDUCT,0.0,2)
c      CALL PLCHAR(XL3 ,YL,CS,'P hp= ' ,0.0,6)
c      CALL PLNUMB(999.,YL,CS, TDIM/746.,0.0,1)
      CALL PLCHAR(XL3 ,YL,CS,'Ptot = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, PDIM/1000.,0.0,3)
      CALL PLCHAR(XL4 ,YL,CS,'Qtot = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, QDIM,0.0,3)
      CALL PLCHAR(XL5 ,YL,CS,'FOM  = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, FOM,0.0,4)
C
      YL = YL - 2.5*CS
      CALL PLSUBS(XL1 ,YL,CS,'C'      ,0.0,1,PLCHAR)
      CALL PLCHAR(XL1 ,YL,CS,'T    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS,TC,0.0,4)
      CALL PLSUBS(XL2 ,YL,CS,'C'      ,0.0,1,PLCHAR)
      CALL PLCHAR(XL2 ,YL,CS,'P    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS,PC,0.0,4)
      CALL PLSUBS(XL3 ,YL,CS,'T'      ,0.0,1,PLCHAR)
      CALL PLCHAR(XL3 ,YL,CS,'C    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CT      ,0.0,4)
      CALL PLSUBS(XL4 ,YL,CS,'P'      ,0.0,1,PLCHAR)
      CALL PLCHAR(XL4 ,YL,CS,'C    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CP      ,0.0,4)
      CALL PLMATH(XL5 ,YL,CS,'h    = ',0.0,7)
      CALL PLNUMB(999.,YL,CS,EFFTOT,0.0,4)
C
C
c      IF(ADV.LT.0.1) THEN
c       YL = YL - 2.5*CS
c       CALL PLCHAR(XL1 ,YL,CS,'Static:' ,0.0,-1)
c       CALL PLSUBS(XL2 ,YL,CS,'TH'     ,0.0,2,PLCHAR)
c       CALL PLCHAR(XL2 ,YL,CS,'C   = ' ,0.0,6)
c       CALL PLNUMB(999.,YL,CS, CTH     ,0.0,6)
c       CALL PLSUBS(XL3 ,YL,CS,'PH'     ,0.0,2,PLCHAR)
c       CALL PLCHAR(XL3 ,YL,CS,'C   = ' ,0.0,6)
c       CALL PLNUMB(999.,YL,CS, CPH     ,0.0,6)
c       CALL PLSUBS(XL4 ,YL,CS,'TH'     ,0.0,2,PLCHAR)
c       CALL PLCHAR(XL4 ,YL,CS,'C /'    ,0.0,3)
c       CALL PLMATH(999.,YL,CS,'s= '    ,0.0,3)
c       CALL PLNUMB(999.,YL,CS, CTOS    ,0.0,4)
C
c       YL = YL - 1.0*CS
c       CALL NEWPEN(1)
c       CALL PLOT(0.0,YL,3)
c       CALL PLOT(1.0,YL,2)
c       CALL NEWPEN(3)
c      ENDIF
C=====================================================================
C
C---- Blade radial station and circulation for plotting
      BLDS = FLOAT(NRBLD(N))
C
      IF(TGAP.GT.0.0 .AND. OMEGA(N).NE.0.0) THEN
         NRL = NRC-1
      ELSE
         NRL = NRC
      ENDIF
C
      DO I = 1, NRL
        XI(I) = YRC(I,N)/RTIP(N)
        GAMB(I) = BGAM(I,N)/BLDS
      END DO
      CALL MINMAX(NRL,GAMB,GMIN,GMAX)
      CALL MINMAX(NRL,CLR(1,N),CLMIN,CLMAX)
C
C---- CL-axis limit and increment
      CALL AXISADJ(CLMIN,CLMAX,CLSPAN,DCL,NCL)
      IF(CLMIN.LT.0.0 .AND. CLMAX.LE.0.0) THEN
        CLMIN = -2.0
        CLMAX = 0.0
      ELSEIF(CLMIN.GE.0.0 .AND. CLMAX.GT.0.0) THEN
        CLMAX = 2.0
        CLMIN = 0.0
      ENDIF
      CALL AXISADJ(CLMIN,CLMAX,CLSPAN,DCL,NCL)
C
      CALL NEWPEN(2)
c      CALL YAXIS(0.0,0.0, PLPAR,0.2*PLPAR, 0.0,0.2, CS,1)
      CLOFF = CLMIN
      CLFAC = PLPAR/CLSPAN
      YAX = PLPAR
      DYAX = YAX/FLOAT(NCL-1)
C
C---- CL vs r/R
      CALL NEWCOLORNAME('red')
      CALL NEWPEN(2)
      CALL YAXIS(1.1,0.0,YAX,DYAX,CLMIN,DCL,-CS,1)
      CALL NEWPEN(3)
      XPLT = 1.1 + 1.0*CSL
      YPLT = CLFAC*(CLMAX-1.5*DCL-CLOFF) - 0.5*CSL
      CALL PLCHAR(XPLT,YPLT,1.2*CSL,'c',0.0,1)
      CALL PLSUBS(XPLT,YPLT,1.2*CSL,'V',0.0,1,PLMATH)
C
      IF(IRTYPE(N).EQ.2) THEN
       CALL NEWPEN(4)
       CALL XYLINE(NRL,XI,CLR(1,N),0.0,1.0,CLOFF,CLFAC,1)     
       CALL NEWPEN(3)
       IL = NRC/2
       XPLT =                 XI(IL)
       YPLT = 1.0*CSL + CLFAC*(CLR(IL,N)-CLOFF)
       CALL PLCHAR(XPLT,YPLT,1.2*CSL,'c',0.0,1)
       CALL PLSUBS(XPLT,YPLT,1.2*CSL,'V',0.0,1,PLMATH)
      ENDIF
C
C---- normalized Gamma vs r/R
      CALL NEWCOLORNAME('cyan')
ccc      CALL SCALIT(NRC,GAMB,0.0,GFAC)
      CALL AXISADJ(GMIN,GMAX,GSPAN,GDEL,NGDEL)
      IF(GMIN.LT.0.0 .AND. GMAX.LE.0.0) THEN
        GMIN = MIN(-2.0,GMIN)
        GMAX = 0.0
      ELSEIF(GMIN.GE.0.0 .AND. GMAX.GT.0.0) THEN
        GMAX = MAX(2.0,GMAX)
        GMIN = 0.0
      ENDIF
      CALL AXISADJ(GMIN,GMAX,GSPAN,GDEL,NGDEL)
      GOFF = GMIN
      GFAC = PLPAR/GSPAN
      YAX = PLPAR
      DYAX = YAX/FLOAT(NGDEL-1)
C
      CALL NEWPEN(4)
      CALL XYLINE(NRL,XI,GAMB,0.0,1.0,GOFF,GFAC,2)
      IL = NRC/5
      CALL NEWPEN(2)
      CALL YAXIS(0.0,0.0,YAX,DYAX,GMIN,GDEL, CS,1)
      CALL NEWPEN(3)
      XPLT = 0.0 - 5.5*CSL
      YPLT = GFAC*(GMAX-1.5*GDEL-GOFF) - 0.5*CSL
      CALL PLMATH(XPLT,YPLT,CSL,'G',0.0,1)
      CALL PLCHAR(XPLT,YPLT,CS, '  blade',0.0,7)
C
      XPLT =  -2.0*CSL +             XI(IL)
      YPLT =   0.5*CSL + GFAC*(GAMB(IL)-GOFF)
      CALL PLMATH(XPLT,YPLT,CSL,'G',0.0,1)
      CALL PLCHAR(XPLT,YPLT,CS, '  blade',0.0,7)
C
C---- Mach vs r/R
      EOFF = 0.0
      EFAC = PLPAR/1.0
      EDEL = 0.2
      YAX = PLPAR
      DYAX = 0.2*YAX
C
      CALL NEWCOLORNAME('yellow')
      CALL NEWPEN(2)
      CALL YAXIS(1.0,0.0,-YAX,DYAX,0.0,EDEL,CS,1)
      CALL NEWPEN(3)
      XPLT = 1.0 - 2.0*CSL
      YPLT = EFAC*(1.0-1.5*EDEL-EOFF) - 0.5*CSL
      CALL PLCHAR(XPLT,YPLT,CSL,'M',0.0,1)
C
      CALL NEWPEN(4)
      CALL XYLINE(NRL,XI,MACHR(1,N),0.0,1.0,EOFF,EFAC,4)
      IL = 4*NRC/5
      CALL NEWPEN(3)
      XPLT =                 XI(IL)
      YPLT = 0.5*CSL + EFAC*(MACHR(IL,N)-EOFF)
      CALL PLCHAR(XPLT,YPLT,CSL,'M',0.0,1)
C
C---- local efficiency vs r/R
      CALL NEWPEN(2)
c      CALL YAXIS(0.0,0.0, PLPAR,0.2*PLPAR, 0.0,0.2, CS,1)
C
      CALL NEWCOLORNAME('green')
      CALL NEWPEN(2)
      CALL YAXIS(1.0,0.0,YAX,DYAX,0.0,EDEL,-CS,1)
      CALL NEWPEN(3)
      XPLT = 1.0 + 1.0*CSL
      YPLT = EFAC*(1.0-1.5*EDEL-EOFF) - 0.5*CSL
      CALL PLMATH(XPLT,YPLT,1.1*CSL,'h'    ,0.0,1)
      CALL PLSUBS(XPLT,YPLT,    CSL,'local',0.0,5,PLCHAR)
C
      DO I=1, NRL
C------ HHY define efficiency
        IF(DQII(I,N)*OMEGA(N).NE.0.) THEN
          EFFI = QINF*DTII(I,N)/(DQII(I,N)*OMEGA(N))
        ELSE
          EFFI = 0.0
        ENDIF
cc        W1(I) = EFFI*EFFP(I)
        W1(I) = EFFI
      ENDDO
C
      CALL NEWPEN(4)
      CALL XYLINE(NRL,XI,W1,0.0,1.0,0.0,EFAC,3)
      IL = 2*NRC/5
      CALL NEWPEN(3)
      XPLT =                 XI(IL)
      YPLT = 1.0*CSL + EFAC*W1(IL)
      CALL PLMATH(XPLT,YPLT,1.1*CSL,'h'    ,0.0,1)
      CALL PLSUBS(XPLT,YPLT,    CSL,'local',0.0,5,PLCHAR)
C
C
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      XYOFF(1) = 0.
      XYOFF(2) = 0.
      XYFAC(1) = 1.0
      XYFAC(2) = PLPAR
C
      RETURN
      END ! BLDCLPLT



      SUBROUTINE ROTVPLT(NR,LABS,LREL)
C-------------------------------------------
C     Plots velocities just downstream of rotor disk vs r/R
C-------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      CHARACTER*80 PNAME
      LOGICAL LABS,LREL
      EXTERNAL PLCHAR,PLMATH
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
      DIMENSION W1(IRX), W2(IRX), W3(IRX), 
     &          W4(IRX), W5(IRX), W6(IRX), 
     &          W7(IRX), W8(IRX), W9(IRX) 
      DIMENSION XI(IRX)
C
      IF(NROTOR.LE.1) THEN
         PNAME=NAME
      ELSE
         CALL GETPNAME(PNAME,NAME,NR)
      ENDIF
C
C---- set blade disk to plot
      N = NR
      RPM = 30.0*OMEGA(N)/PI
C
C---- plot aspect ratio
      PLPAR = 0.6
      PLFAC = 0.8
C
C---- character size for axis numbers, labels
      CS  = 0.9*CSIZE
      CSL = CS*1.4
C
C
      CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFAC*SIZE,LPLOT,LLAND)
      CALL PLOTABS(1.25,0.8,-3)
C
      CALL GETCOLOR(ICOL0)
C
C---- case title
c      CALL NEWPEN(4)
      XT = 0.0
c      YT = WFAC*WMAX + 1.5*CSL + 2.5*CS
C
C===================================================================
C
C---- get physical blade geometry
      CALL GETROTG(N)
C
C---- mean radius
      YMNR = SQRT(0.5*ADISK(N)/PI + RHUB(N)**2)
C
C---- blade solidity
      CALL SPLINE(CHRP(1,N),W1,YRP,NRP)
      CMNR = SEVAL(YMNR,CHRP(1,N),W1,YRP,NRP)
      SIGMR = FLOAT(NRBLD(N))*CMNR/(2.0*PI*YMNR)
C
C---- blade angles and twist
      BHUB = BETARP(1,N)
      BTIP = BETARP(NRP,N)
C      B2   = BETAR(2,N)
C      BEXT = BHUB - YRC(1,N)*(BHUB-B2)/(YRC(1,N)-YRC(2,N))
      BTWIST = BHUB - BTIP
      CALL SPLINE(BETARP(1,N),W2,YRP,NRP)
      BETAMR = SEVAL(YMNR,BETARP(1,N),W2,YRP,NRP)
C----------------------------------------------------------------
C
C---- system calcs
C
      RPM = 30.0*OMEGA(N)/PI
      IF(OMEGA(N).NE.0.0) THEN
        ADV = QINF/(OMEGA(N)*RTIP(N))
      ELSE
        ADV = 0.0
      ENDIF
C
C---- dimensional thrust, power, torque, rpm
      TDIM = TTOT + TDUCT
      QDIM = QTOT
      PDIM = PTOT
      TVDIM = TVIS
      PVDIM = PVIS
      TINV  = TTOT - TVIS 
      QINV  = QTOT - QVIS 
      PINV  = PTOT - PVIS 
C
C---- Thrust and power coefficients based on rotational speed
      DIA = 2.0*RTIP(N)
      EN = ABS(RPM/60.0)
      IF(EN.EQ.0.0) THEN
        CT = 0.0
        CP = 0.0
      ELSE
        CT = TDIM/(RHO*EN**2*DIA**4)
        CP = PDIM/(RHO*EN**3*DIA**5)
      ENDIF
C
C---- Thrust and power coefficients based on forward speed
      IF(QINF.GT.0.0) THEN
        TC = TDIM/(0.5*RHO*QINF**2 * PI*RTIP(N)**2)
        PC = PDIM/(0.5*RHO*QINF**3 * PI*RTIP(N)**2)
      ELSE
        TC = 0.0
        PC = 0.0
      ENDIF
C---- Thrust and power coefficients based on tip speed, disk area
      VTIP = ABS(OMEGA(N)*RTIP(N))
      IF(VTIP.NE.0.0) THEN
       CT0  = TDIM/(RHO*ADISK(N)*VTIP**2)
       CP0  = PDIM/(RHO*ADISK(N)*VTIP**3)
       CTOS = CT0 / SIGMR
       IF(CP0.NE.0.0 .AND. CT0.GE.0.0) THEN
         FOM = ABS(CT0)**1.5 / CP0 / 2.0
       ELSE
         FOM = 0.0
       ENDIF
      ELSE
       CT0  = 0.0
       CP0  = 0.0
       CTOS = 0.0
       FOM  = 0.0
      ENDIF
C
C---- overall efficiency
      IF(PDIM.NE.0.0) THEN
        EFFTOT = QINF*TDIM/PDIM
       ELSE
        EFFTOT = 0.0
      ENDIF
C---- induced efficiency (including nacelle thrust effect)
      IF(PINV.NE.0.0) THEN
        EFFIND = QINF*(TINV+TDUCT)/PINV
      ELSE
        EFFIND = 0.0
      ENDIF
C---- ideal (actuator disk) efficiency
      IF(TC.EQ.0.0) THEN
        EIDEAL = 0.0
      ELSE
        TCLIM = MAX( -1.0 , TC )
        EIDEAL = 2.0 / (1.0 + SQRT(TCLIM + 1.0))
      ENDIF
C
C
C---- rotor calcs
C---- dimensional thrust, power, torque, rpm
      TDIM2  = TTOTR(N)
      QDIM2  = QTOTR(N)
      PDIM2  = QDIM2    * OMEGA(N)
      TVDIM2 = TVISR(N)
      PVDIM2 = QVISR(N) * OMEGA(N)
      TINV2  = TINVR(N)
      QINV2  = QINVR(N)
      PINV2  = QINV2    * OMEGA(N)
C
C
C---- Define reference quantities for coefficients
      RREF = RTIP(N)
      AREF = ADISK(N)
      OMEGREF = OMEGA(N)
C
C---- standard thrust/power coefficients based on rotational speed
C
      IF(EN.EQ.0.0) THEN
        CT2 = 0.0
        CP2 = 0.0
      ELSE
        CT2 = TDIM2/(RHO*EN**2*DIA**4)
        CP2 = PDIM2/(RHO*EN**3*DIA**5)
      ENDIF
C
C---- standard thrust/power coefficients based on forward speed
      IF(QINF.GT.0.0) THEN
        TC2 = TDIM2/(0.5*RHO*QINF**2 * PI*RREF**2)
        PC2 = PDIM2/(0.5*RHO*QINF**3 * PI*RREF**2)
      ELSE
        TC2 = 0.0
        PC2 = 0.0
      ENDIF
C
C---- efficiency for rotor 
      IF(PDIM2.NE.0.0) THEN
         EFFTOT2 = QINF*TDIM2/PDIM2
      ELSE
         EFFTOT2 = 0.0
      ENDIF
C
C===================================================================
C---- set label locations
      YL = PLPAR + 15.5*CSL
      XL1 = 0.0
      XL2 = XL1 + 18.0*CS
      XL3 = XL2 + 18.0*CS
      XL4 = XL3 + 18.0*CS
      XL5 = XL4 + 18.0*CS
C
      CALL NEWPEN(4)
      CALL PLCHAR(XL1,YL,CSL,PNAME,0.0,-1)
C
      YL = YL - 1.0*CS
      CALL NEWPEN(1)
      CALL PLOT(0.0,YL,3)
      CALL PLOT(1.0,YL,2)
C
      CALL NEWPEN(3)
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'blds = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, FLOAT(NRBLD(N)),0.0,-1)
      CALL PLCHAR(XL2 ,YL,CS,'R(m) = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, RTIP(N)    ,0.0,3)
      CALL PLCHAR(XL3 ,YL,CS,'Adisk= '  ,0.0,7)
      CALL PLNUMB(999.,YL,CS, ADISK(N)  ,0.0,4)
      CALL PLSUBS(XL4 ,YL,0.8*CS,'mean',0.0,4,PLCHAR)
      CALL PLMATH(XL4 ,YL,CS,'s    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, SIGMR   ,0.0,4)
      CALL PLSUBS(XL5 ,YL,CS,'twist'  ,0.0,5,PLCHAR)
      CALL PLMATH(XL5 ,YL,CS,'b    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BTWIST/DTR  ,0.0,3)
C
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'Vm/s = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, QINF    ,0.0,3)
      CALL PLCHAR(XL2 ,YL,CS,'RPM  = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, RPM     ,0.0,1)
      CALL PLCHAR(XL3, YL,CS,'   R   ' ,0.0,7)
      CALL PLMATH(XL3, YL,CS,'  W    ' ,0.0,7)
      CALL PLCHAR(XL3, YL,CS,' /     ' ,0.0,7)
      CALL PLCHAR(XL3, YL,CS,'V    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, ADV     ,0.0,4)
      CALL PLCHAR(XL4 ,YL,CS,'h(km)= ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS,ALTH,0.0,3)
c      CALL PLCHAR(XL4 ,YL,CS,'hKft= ' ,0.0,6)
c      CALL PLNUMB(999.,YL,CS,3.281*ALT,0.0,0)
      CALL PLSUBS(XL5 ,YL,CS,'tip'    ,0.0,3,PLCHAR)
      CALL PLMATH(XL5 ,YL,CS,'b    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BTIP/DTR,0.0,3)
C
C
      YL = YL - 2.5*CS
      CALL PLCHAR(XL1 ,YL,CS,'Vavg = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, VAAVG(N)   ,0.0,3)
      CALL PLCHAR(XL2 ,YL,CS,'T(N) = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, TDIM2      ,0.0,2)
      CALL PLCHAR(XL3 ,YL,CS,'P(kW)= '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, PDIM2/1000.,0.0,3)
      CALL PLCHAR(XL4 ,YL,CS,'Q(Nm)= '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, QDIM2      ,0.0,2)
      CALL PLSUBS(XL5 ,YL,CS,'hub'       ,0.0,3,PLCHAR)
      CALL PLMATH(XL5 ,YL,CS,'b    = '   ,0.0,7)
      CALL PLNUMB(999.,YL,CS, BHUB/DTR   ,0.0,3)
C
      YL = YL - 2.5*CS
      CALL PLSUBS(XL1 ,YL,CS,'C'       ,0.0,1,PLCHAR)
      CALL PLCHAR(XL1 ,YL,CS,'T    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, TC2      ,0.0,4)
      CALL PLSUBS(XL2 ,YL,CS,'C'       ,0.0,1,PLCHAR)
      CALL PLCHAR(XL2 ,YL,CS,'P    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, PC2      ,0.0,4)
      CALL PLSUBS(XL3 ,YL,CS,'T'       ,0.0,1,PLCHAR)
      CALL PLCHAR(XL3 ,YL,CS,'C    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CT2      ,0.0,4)
      CALL PLSUBS(XL4 ,YL,CS,'P'       ,0.0,1,PLCHAR)
      CALL PLCHAR(XL4 ,YL,CS,'C    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CP2      ,0.0,4)
      CALL PLMATH(XL5 ,YL,CS,'h    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, EFFTOT2  ,0.0,4)
C
      YL = YL - 1.0*CS
      CALL NEWPEN(1)
      CALL PLOT(0.0,YL,3)
      CALL PLOT(1.0,YL,2)
      CALL NEWPEN(3)
C
      YL = YL - 2.5*CS
c      CALL PLCHAR(XL1 ,YL,CS,'T lb= ' ,0.0,6)
c      CALL PLNUMB(999.,YL,CS, TDIM/4.45,0.0,1)
      CALL PLCHAR(XL1 ,YL,CS,'Ttot = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, TDIM,0.0,2)
      CALL PLCHAR(XL2 ,YL,CS,'Tduct= ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, TDUCT,0.0,2)
c      CALL PLCHAR(XL3 ,YL,CS,'P hp= ' ,0.0,6)
c      CALL PLNUMB(999.,YL,CS, TDIM/746.,0.0,1)
      CALL PLCHAR(XL3 ,YL,CS,'Ptot = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, PDIM/1000.,0.0,3)
      CALL PLCHAR(XL4 ,YL,CS,'Qtot = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, QDIM,0.0,3)
      CALL PLCHAR(XL5 ,YL,CS,'FOM  = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, FOM,0.0,4)
C
      YL = YL - 2.5*CS
      CALL PLSUBS(XL1 ,YL,CS,'C'      ,0.0,1,PLCHAR)
      CALL PLCHAR(XL1 ,YL,CS,'T    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS,TC,0.0,4)
      CALL PLSUBS(XL2 ,YL,CS,'C'      ,0.0,1,PLCHAR)
      CALL PLCHAR(XL2 ,YL,CS,'P    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS,PC,0.0,4)
      CALL PLSUBS(XL3 ,YL,CS,'T'      ,0.0,1,PLCHAR)
      CALL PLCHAR(XL3 ,YL,CS,'C    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CT      ,0.0,4)
      CALL PLSUBS(XL4 ,YL,CS,'P'      ,0.0,1,PLCHAR)
      CALL PLCHAR(XL4 ,YL,CS,'C    = ' ,0.0,7)
      CALL PLNUMB(999.,YL,CS, CP      ,0.0,4)
      CALL PLMATH(XL5 ,YL,CS,'h    = ',0.0,7)
      CALL PLNUMB(999.,YL,CS,EFFTOT,0.0,4)
C
C
C================================================================
C
      BLDS = FLOAT(NRBLD(N))
      XIB = 0.67
      XIA = 0.33
      XIBX = 1.0
      XIAX = 1.0
C
C  Note: 
C     The blade relative velocity is defined just downstream of the blade
C     lifting line in a frame moving with the rotor.  
C     The absolute frame velocity is defined just downstream of the blade
C     lifting line.  
C     At the blade lifting line the tangential induced velocity is 1/2 of
C     the tangential velocity in the slipstream immediately downstream of
C     the blade.
C
      IF(TGAP.GT.0.0 .AND. OMEGA(N).NE.0.0) THEN
         NRL = NRC-1
      ELSE
         NRL = NRC
      ENDIF
C
      DO I = 1, NRL
        XI(I) =  YRC(I,N)/RTIP(N)
C---- blade relative velocity
        W1(I) = VREL(1,I,N)
        W2(I) = VREL(3,I,N)
        W3(I) = VREL(2,I,N)
C------ absolute frame velocity
        W4(I) = VABS(1,I,N)
        W5(I) = VABS(3,I,N)
        W6(I) = VABS(2,I,N)
C
ccc     write(*,*) i,xi(i),w1(i),w2(i),w3(i),w4(i)
C---- find radial stations closest to 1/2 and 3/4 radius for labels
        IF(ABS(XI(I)-XIA).LT.XIAX) THEN
          XIAX = ABS(XI(I)-XIA)
          IA = I
        ENDIF
        IF(ABS(XI(I)-XIB).LT.XIBX) THEN
          XIBX = ABS(XI(I)-XIB)
          IB = I
        ENDIF
      ENDDO
C
      WMIN = 0.
      WMAX = 0.
      DO I = 1, NRL
        IF(LREL) THEN
          WMIN = MIN( WMIN, W1(I), W2(I), W3(I) )
          WMAX = MAX( WMAX, W1(I), W2(I), W3(I) )
        ENDIF
        IF(LABS) THEN
          WMIN = MIN( WMIN, W4(I), W5(I), W6(I) )
          WMAX = MAX( WMAX, W4(I), W5(I), W6(I) )
        ENDIF
      ENDDO
C
      CALL SCALIT(1,WMIN,0.0,WMNFAC)
      CALL SCALIT(1,WMAX,0.0,WMXFAC)
      WFAC = MIN( WMNFAC , WMXFAC )
      WDEL = 1.0 / (5.0*WFAC)
      IF(WMIN .LT. 0.0) WMIN = -WDEL * AINT( -WMIN/WDEL + 0.99 )
      IF(WMAX .GT. 0.0) WMAX =  WDEL * AINT(  WMAX/WDEL + 0.99 )
C
      CALL AXISADJ(WMIN,WMAX,WSPAN,WDEL,ntics)
      WFAC = PLPAR / (WMAX - WMIN)
C
C====================================================================
      CALL PLOT(0.0,-WFAC*WMIN,-3)
      CALL NEWPEN(1)
      CALL PLOT(0.0,0.0,3)
      CALL PLOT(1.0,0.0,2)
C
      CALL NEWPEN(2)
      CALL XAXIS(0.0,WFAC*WMIN, 1.0,0.2    , 0.0,0.2, CS,1)
      CALL NEWPEN(3)
      CALL PLCHAR(0.7-1.5*CSL,WFAC*WMIN-3.0*CSL,CSL,'r/R',0.0,3)
C
cc      IF(LGRID) THEN
       CALL NEWPEN(1)
       NXG = 5
       NYG = INT( (WMAX-WMIN)/WDEL + 0.0001 )
       CALL PLGRID(0.0,WFAC*WMIN, NXG,0.2, NYG,WFAC*WDEL, LMASK2 )
cc      ENDIF
C
      CALL NEWPEN(2)
      CALL YAXIS(0.0,WFAC*WMIN,WFAC*(WMAX-WMIN),WFAC*WDEL,
     &            WMIN,WDEL, CS,-2)
      CALL NEWPEN(3)

      XL = -5.0*CSL
      YL = WFAC*(WMAX-1.5*WDEL) - 0.5*CSL
      CALL PLCHAR(XL,YL,CSL,'V'  ,0.0,1)
      CALL PLCHAR(XL+1.2*CSL,YL-0.0*CSL,0.8*CSL,'m/s'  ,0.0,3)
ccc      XL = -3.5*CSL
ccc      YL = WFAC*(WMAX-1.5*WDEL) - 0.5*CSL
ccc      CALL PLCHAR(XL,YL,CSL,'v/V'  ,0.0,3)

c      CALL PLSUBS(XL,YL+0.6*CSL,CSL,'axi',0.0,3,PLCHAR)
c      CALL PLOT(XL        ,YL,3)
c      CALL PLOT(XL+3.0*CSL,YL,2)
c      CALL PLCHAR(XL+1.0*CSL,YL-1.4*CSL,CSL,'V'    ,0.0,1)
C     
C---- Vaxial vs r/R
      CALL NEWCOLORNAME('cyan')
      CALL NEWPEN(4)
      XL = XI(IA) + 1.5*CS
      IF(LREL) THEN
       CALL XYLINE(NRL,XI,W1,0.0,1.0,0.0, WFAC,1)
       YL =  WFAC*W1(IA) + 1.5*CS
       CALL PLCHAR(XL,YL+0.6*CS,CS,'V'  ,0.0,1)
       CALL PLSUBS(XL,YL+0.6*CS,CS,'m',0.0,1,PLCHAR) 
      ENDIF
      IF(LABS) THEN
        ILIN = 1
        IF(LREL) ILIN = 2
        CALL XYLINE(NRL,XI,W4,0.0,1.0,0.0, WFAC,ILIN)
        YL =  WFAC*W4(IA) + 1.5*CS
        CALL PLCHAR(XL,YL+0.6*CS,CS,'V'  ,0.0,1)
        CALL PLSUBS(XL,YL+0.6*CS,CS,'m',0.0,1,PLCHAR)
      ENDIF
      CALL NEWPEN(3)
C
C---- Vduct (axial) vs r/R
      CALL NEWCOLORNAME('green')
      CALL NEWPEN(4)
      XL = XI(IB) + 1.5*CS
      IF(LREL) THEN
        CALL XYLINE(NRL,XI,W3,0.0,1.0,0.0, WFAC,1)
        YL =  WFAC*W3(IB) + 1.5*CS
        CALL PLCHAR(XL,YL+0.6*CS,CS,'V'  ,0.0,1)
        CALL PLSUBS(XL,YL+0.6*CS,CS,'r',0.0,1,PLCHAR)
      ENDIF
      IF(LABS) THEN
        ILIN = 1
        IF(LREL) ILIN = 2
        CALL XYLINE(NRL,XI,W6,0.0,1.0,0.0, WFAC,ILIN)
        YL =  WFAC*W6(IB) + 1.5*CS
        CALL PLCHAR(XL,YL+0.6*CS,CS,'V'  ,0.0,1)
        CALL PLSUBS(XL,YL+0.6*CS,CS,'r',0.0,1,PLCHAR)
      ENDIF
      CALL NEWPEN(3)
cc      CALL PLOT(XL       ,YL,3)
cc      CALL PLOT(XL+3.0*CS,YL,2)
cc      CALL PLCHAR(XL+1.0*CS,YL-1.4*CS,CS,'V'    ,0.0,1)
C
C---- Vtangential vs r/R
      CALL NEWCOLORNAME('red')
c      CALL NEWPEN(3)
c      XL = -3.5*CSL
c      YL = WFAC*(WMAX-0.5*WDEL) + 0.1*CSL
c      CALL PLCHAR(XL,YL+0.6*CSL,CSL,'v'  ,0.0,1)
c      CALL PLSUBS(XL,YL+0.6*CSL,CSL,'tan',0.0,3,PLCHAR)
c      CALL PLOT(XL        ,YL,3)
c      CALL PLOT(XL+3.0*CSL,YL,2)
c      CALL PLCHAR(XL+1.0*CSL,YL-1.4*CSL,CSL,'V'    ,0.0,1)
C     
      CALL NEWPEN(4)
      XL = XI(IB) + 1.5*CS
      IF(LREL) THEN
        CALL XYLINE(NRL,XI,W2,0.0,1.0,0.0, WFAC,1)
        YL =  WFAC*W2(IB) + 1.5*CS
        CALL PLCHAR(XL,YL+0.6*CS,CS,'V'  ,0.0,1)
        CALL PLSUBS(XL,YL+0.6*CS,CS,'tan',0.0,3,PLCHAR)
      ENDIF
      IF(LABS) THEN
       ILIN = 1
       IF(LREL) ILIN = 2
       CALL XYLINE(NRL,XI,W5,0.0,1.0,0.0, WFAC,ILIN)
       YL =  WFAC*W5(IB) + 1.5*CS
       CALL PLCHAR(XL,YL+0.6*CS,CS,'V'  ,0.0,1)
       CALL PLSUBS(XL,YL+0.6*CS,CS,'tan',0.0,3,PLCHAR)
      ENDIF
      CALL NEWPEN(3)
cc      CALL PLOT(XL       ,YL,3)
cc      CALL PLOT(XL+3.0*CS,YL,2)
cc      CALL PLCHAR(XL+1.0*CS,YL-1.4*CS,CS,'V'    ,0.0,1)
C
      CALL NEWCOLOR(ICOL0)
C
      W8(1) = -0.1    + 1.0*CSL
      W8(2) =         - 1.0*CSL
      W9(1) =  0.5*CSL
      W9(2) =  0.5*CSL
C
      XL = 0.05 
      YL = WFAC*WMAX - 2.0*CSL
      CALL NEWPEN(3)
      CALL PLCHAR(XL,YL,0.8*CSL,'Slipstream velocities',0.0,-1)
C
      XL = 0.1  + 0.5*CSL 
      IF(LREL) THEN
        YL = YL - 2.0*CSL
        CALL NEWPEN(4)
        CALL XYLINE(2,W8,W9,-XL,1.0,-YL,1.0,1)
        CALL NEWPEN(3)
        CALL PLCHAR(XL,YL,0.8*CSL,'Blade relative',0.0,-1)
      ENDIF
C
      IF(LABS) THEN
        YL = YL - 2.0*CSL
        CALL NEWPEN(4)
        CALL XYLINE(2,W8,W9,-XL,1.0,-YL,1.0,2)
        CALL NEWPEN(3)
        CALL PLCHAR(XL,YL,0.8*CSL,'Absolute frame',0.0,-1)
      ENDIF
C
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      XYOFF(1) = 0.
      XYOFF(2) = 0.
      XYFAC(1) = 1.0
      XYFAC(2) =  WFAC
C
      RETURN
      END ! RVIPLT


      SUBROUTINE PLOTBLDDATA(NR,PLTTYPE)
C----Plots data on blade disk # NR
C
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      LOGICAL LERROR, LXRELIMIT, LYRELIMIT
      LOGICAL LPLT_DATA, LEGEND2D
      CHARACTER*4 OPT
      CHARACTER*6 XVAR, YVAR
      CHARACTER PLTTYPE*(*)
      CHARACTER DATALBL*20, LEGENDDATA*80
      CHARACTER*80 LINE, FNAME*128
      CHARACTER*80 TITLE1, TITLE2, PNAME
C
      DIMENSION RINPUT(10), XLIM(2), YLIM(2)
C
C--- Temporary data for plotting 
      PARAMETER (IXP=500)
      DIMENSION  NPTDATA(IXP)
      DIMENSION  INFODATA(2,IXP), IFLW2DATA(0:IXP)
      DIMENSION  XDATA(IXP,IXP), YDATA(IXP,IXP)
      DIMENSION  LEGENDDATA(IXP)
C
      DATA XVAR, YVAR / 'XI', 'GAM' /
C
      PAR = 0.8
      LUTEMP = 3
      LXRELIMIT = .TRUE.
      LYRELIMIT = .TRUE.
      XMINSPEC = 999.
      XMAXSPEC = 999.
      YMINSPEC = 999.
      YMAXSPEC = 999.
      DATALBL = PLTTYPE
C
 1    FORMAT(A)
      NDATA = 0
      IF(IRTYPE(NR).EQ.1) THEN
        TITLE1 = 'DFDC Actuator Disk Data'
      ELSEIF(IRTYPE(NR).EQ.2) THEN
        TITLE1 = 'DFDC Blade Data'
      ENDIF
C
      IF(NROTOR.LE.1) THEN
         PNAME=NAME
      ELSE
         CALL GETPNAME(PNAME,NAME,NR)
      ENDIF
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
      WRITE(*,20) XVAR,YVAR,LEGEND
 20   FORMAT(/' ================================================'
     &       /'   A  abscissa select   currently: ',A,
     &       /'   O  ordinate select   currently: ',A,
     &       /'   L  limits for plot',
     &       /'   Z  zoom plot with cursor',
     &       /'   R  reset plot limits',
     &       /'   AN annotate plot',
     &       /'   LG legend for plot   currently: ',L2,
     &       /'   H  hardcopy plot',
     &       /'   W  write plot data to file',
     &      //' Select option (or <return>):  ', $)
C
      READ(*,1,ERR=10) OPT
      CALL LC2UC(OPT)
C
      IF(OPT.EQ.' ') THEN
        RETURN
C
C***************************************************
C---- Make hardcopy
       ELSE IF(OPT.EQ.'H') THEN
        CALL REPLOT(IDEVPS)
        GO TO 10
C
C***************************************************
C---- Write plot data to file
       ELSE IF(OPT.EQ.'W') THEN
        IF(NDATA.GT.0) THEN
         CALL ASKS('Enter plot data save filename^',FNAME)
         OPEN(LUTEMP,FILE=FNAME,STATUS='UNKNOWN',ERR=10)
         NDT = 1
         NPT = NPTDATA(NDT)
C--- Write out X data, Y data 
         DO I = 1, NPT
           WRITE(LUTEMP,51) XDATA(I,NDT),YDATA(I,NDT) 
         END DO
 51      FORMAT(2(G13.6,2X))
         CLOSE(LUTEMP)
        ENDIF
        GO TO 10
C
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
C---- Toggle legend flag
       ELSE IF(OPT.EQ.'LG') THEN
        LEGEND2D = .NOT.LEGEND2D
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
C***************************************************
C---- Change Abscissa 
       ELSE IF(OPT.EQ.'A') THEN
 40     WRITE(*,*) ' '
        WRITE(*,45) XVAR
 45     FORMAT(/' ================================================'
     &       /' Abscissa is: ',A,
     &       /'   R, XI '
     &      //' Select option (or <return>):  ', $)
        READ(*,1,ERR=40) OPT
        CALL LC2UC(OPT)
        IF(OPT.EQ.' ') THEN
          GO TO 10
         ELSE IF(OPT.EQ.'XI') THEN
          XVAR = 'XI'
         ELSE IF(OPT.EQ.'R') THEN
          XVAR = 'R'
        ENDIF
        XMINSPEC = 999.
        XMAXSPEC = 999.
        LXRELIMIT = .TRUE.
C
        IF(IRTYPE(NR).EQ.1) THEN
          TITLE1 = 'DFDC Actuator Disk Data'
        ELSEIF(IRTYPE(NR).EQ.2) THEN
          TITLE1 = 'DFDC Blade Data'
        ENDIF
C
        GO TO 100
C
C***************************************************
C---- Change Ordinate
       ELSE IF(OPT.EQ.'O') THEN
 50       WRITE(*,*) ' '
        WRITE(*,55) YVAR
 55     FORMAT(/' ================================================'
     &       /' Ordinate is: ',A,
     &       /'   CH, BE, ALF, SIG, STG',
     &       /'   GAM, CL, CD, RE, M, CLAlf',
     &       /'   TI, QI, TV, QV, T, Q',
     &       /'   CLSig, PSTAt',
     &       /'   VXA, VRA, VMA, VTA, VVA  ANGA (abs frame slipstrm)',
     &       /'   VXR, VRR, VMR, VTR, VVR  ANGR (rel frame slipstrm)',
     &       /'   VXB, VRB, VMB, VTB, VVB  ANGB (on lifting line)',
     &      //' Select option (or <return>):  ', $)
        READ(*,1,ERR=50) OPT
        CALL LC2UC(OPT)
        YMINSPEC = 999.
        YMAXSPEC = 999.
C
        IF(IRTYPE(NR).EQ.1) THEN
          TITLE1 = 'DFDC Actuator Disk Data'
        ELSEIF(IRTYPE(NR).EQ.2) THEN
          TITLE1 = 'DFDC Blade Data'
        ENDIF
C
        IF(OPT.EQ.' ') THEN
          GO TO 10
         ELSE IF(OPT.EQ.'CH') THEN
          YVAR = 'CHORD'
         ELSE IF(OPT.EQ.'BE') THEN
          YVAR = 'BETA'
         ELSE IF(OPT.EQ.'ALF') THEN
          YVAR = 'ALFA'
         ELSE IF(OPT.EQ.'SIG') THEN
          YVAR = 'SIGMA'
         ELSE IF(OPT.EQ.'STG') THEN
          YVAR = 'STAGR'
C
         ELSE IF(OPT.EQ.'GAM') THEN
          YVAR = 'GAM'
         ELSE IF(OPT.EQ.'CL') THEN
          YVAR = 'CL'
         ELSE IF(OPT.EQ.'CD') THEN
          YVAR = 'CD'
         ELSE IF(OPT.EQ.'RE') THEN
          YVAR = 'RE'
         ELSE IF(OPT.EQ.'M') THEN
          YVAR = 'MACH'
         ELSE IF(OPT.EQ.'CLA' .OR. OPT.EQ.'CLALF') THEN
          YVAR = 'CLalf'
C
         ELSE IF(OPT.EQ.'TI') THEN
          YVAR = 'Tinv'
         ELSE IF(OPT.EQ.'QI') THEN
          YVAR = 'Qinv'
         ELSE IF(OPT.EQ.'TV') THEN
          YVAR = 'Tvis'
         ELSE IF(OPT.EQ.'QV') THEN
          YVAR = 'Qvis'
         ELSE IF(OPT.EQ.'T') THEN
          YVAR = 'Ti+v'
         ELSE IF(OPT.EQ.'Q') THEN
          YVAR = 'Qi+v'
C
         ELSE IF(OPT.EQ.'CLSI' .OR. OPT.EQ.'CLSIG') THEN
          YVAR = 'CLsig'
         ELSE IF(OPT.EQ.'PSTA') THEN
          YVAR = 'Pstat'
C
         ELSE IF(OPT.EQ.'VXA') THEN
          YVAR = 'VXa'
          TITLE1 = 'Blade Absolute Frame Slipstream'
         ELSE IF(OPT.EQ.'VRA') THEN
          YVAR = 'VRa'
          TITLE1 = 'Blade Absolute Frame Slipstream'
         ELSE IF(OPT.EQ.'VMA') THEN
          YVAR = 'VMa'
          TITLE1 = 'Blade Absolute Frame Slipstream'
         ELSE IF(OPT.EQ.'VTA') THEN
          YVAR = 'VTa'
          TITLE1 = 'Blade Absolute Frame Slipstream'
         ELSE IF(OPT.EQ.'VVA') THEN
          YVAR = 'VVa'
          TITLE1 = 'Blade Absolute Frame Slipstream'
         ELSE IF(OPT.EQ.'ANGA') THEN
          YVAR = 'ANGa'
          TITLE1 = 'Blade Absolute Frame Slipstream'
C
         ELSE IF(OPT.EQ.'VXR') THEN
          YVAR = 'VXr'
          TITLE1 = 'Blade Relative Frame Slipstream'
         ELSE IF(OPT.EQ.'VRR') THEN
          YVAR = 'VRr'
          TITLE1 = 'Blade Relative Frame Slipstream'
         ELSE IF(OPT.EQ.'VMR') THEN
          YVAR = 'VMr'
          TITLE1 = 'Blade Relative Frame Slipstream'
         ELSE IF(OPT.EQ.'VTR') THEN
          YVAR = 'VTr'
          TITLE1 = 'Blade Relative Frame Slipstream'
         ELSE IF(OPT.EQ.'VVR') THEN
          YVAR = 'VVr'
          TITLE1 = 'Blade Relative Frame Slipstream'
         ELSE IF(OPT.EQ.'ANGR') THEN
          YVAR = 'ANGr'
          TITLE1 = 'Blade Relative Frame Slipstream'
C
         ELSE IF(OPT.EQ.'VXB') THEN
          YVAR = 'VXb'
          TITLE1 = 'Flow on Blade Lifting Line'
         ELSE IF(OPT.EQ.'VRB') THEN
          YVAR = 'VRb'
          TITLE1 = 'Flow on Blade Lifting Line'
         ELSE IF(OPT.EQ.'VMB') THEN
          YVAR = 'VMb'
          TITLE1 = 'Flow on Blade Lifting Line'
         ELSE IF(OPT.EQ.'VTB') THEN
          YVAR = 'VTb'
          TITLE1 = 'Flow on Blade Lifting Line'
         ELSE IF(OPT.EQ.'VVB') THEN
          YVAR = 'VVb'
          TITLE1 = 'Flow on Blade Lifting Line'
         ELSE IF(OPT.EQ.'ANGB') THEN
          YVAR = 'ANGb'
          TITLE1 = 'Flow on Blade Lifting Line'
         ELSE
            WRITE(*,*)'Command not recognized'
            GO TO 50
        ENDIF
        LYRELIMIT = .TRUE.
        GO TO 100
C
      ENDIF
      GO TO 10
C
C--- Put selected data into arrays
 100  CONTINUE
      NDT = 0
      NCASESX = 1
      DO NS= 1, NCASESX
        NDT = NDT + 1
C
        IF(TGAP.GT.0.0 .AND. OMEGA(NR).NE.0.0) THEN
          NPT = NRC-1
        ELSE
          NPT = NRC
        ENDIF
C
        IFLW = 1
ccc        IFLW = KCASE
        INFODATA(1,NDT) = IFLW
C---- set blade disk to plot
        N = NR
        BLDS = FLOAT(NRBLD(N))
        RPM = 30.0*OMEGA(NR)/PI
C
C--- Legend for each curve 
        CALL STRIP(DATALBL,NLBL)
        WRITE(LEGENDDATA(NDT),101) 'Case', NS,RPM,QINF,ALTH
 101    FORMAT(' ',A,I3,' @ ',F8.1,' rpm  ',F6.2,' m/s  ',
     &         F6.3,' km')
C
        DO I = 1, NPT
C--- Install X data
          IF(XVAR.EQ.'XI') THEN
            XDATA(I,NDT) = YRC(I,N)/RTIP(N)
           ELSEIF(XVAR.EQ.'R') THEN
            XDATA(I,NDT) = YRC(I,N)
          ENDIF
C--- Install Y data
          IF(YVAR.EQ.'CH') THEN
            YDATA(I,NDT) = CHR(I,N)
C
C--- Beta and Alpha in global or local coords
C
           ELSEIF(YVAR.EQ.'BETA') THEN
             IF(RPM.LE.0.0) THEN
                IF(LCOORD) THEN
                  YDATA(I,NDT) = (PI-BETAR(I,N))/DTR
                  TITLE1 = 'DFDC Blade Data (local coords)'
                ELSE
                  YDATA(I,NDT) = BETAR(I,N)/DTR
                  TITLE1 = 'DFDC Blade Data (global coords)'
                ENDIF
             ELSE
                  YDATA(I,NDT) = BETAR(I,N)/DTR                
             ENDIF
C
           ELSEIF(YVAR.EQ.'ALFA') THEN
             IF(RPM.LE.0.0) THEN
                IF(LCOORD) THEN
                  YDATA(I,NDT) = (2.0*AZERO(I,N)-ALFAR(I,N))/DTR
                  TITLE1 = 'DFDC Blade Data (local coords)'
                ELSE
                  YDATA(I,NDT) = ALFAR(I,N)/DTR
                  TITLE1 = 'DFDC Blade Data (global coords)'
                ENDIF
             ELSE
                  YDATA(I,NDT) = ALFAR(I,N)/DTR                
             ENDIF
C
C
           ELSEIF(YVAR.EQ.'SIGMA') THEN
            YDATA(I,NDT) = BLDS*CHR(I,N)/(2.0*PI*YRC(I,N))
           ELSEIF(YVAR.EQ.'STAGR') THEN
            YDATA(I,NDT) = (0.5*PI-BETAR(I,N))/DTR
C
           ELSEIF(YVAR.EQ.'GAM') THEN
            YDATA(I,NDT) = BGAM(I,N)/FLOAT(NRBLD(N))
           ELSEIF(YVAR.EQ.'CL') THEN
            YDATA(I,NDT) = CLR(I,N)
           ELSEIF(YVAR.EQ.'CD') THEN
            YDATA(I,NDT) = CDR(I,N)
           ELSEIF(YVAR.EQ.'RE') THEN
            YDATA(I,NDT) = RER(I,N)
           ELSEIF(YVAR.EQ.'MACH') THEN
            YDATA(I,NDT) = MACHR(I,N)
           ELSEIF(YVAR.EQ.'CLalf') THEN
            YDATA(I,NDT) = CLALF(I,N)
C
           ELSEIF(YVAR.EQ.'Tinv') THEN
            YDATA(I,NDT) = DTII(I,N)
           ELSEIF(YVAR.EQ.'Qinv') THEN
            YDATA(I,NDT) = DQII(I,N)
           ELSEIF(YVAR.EQ.'Tvis') THEN
            YDATA(I,NDT) = DTVI(I,N)
           ELSEIF(YVAR.EQ.'Qvis') THEN
            YDATA(I,NDT) = DQVI(I,N)
           ELSEIF(YVAR.EQ.'Ti+v') THEN
            YDATA(I,NDT) = DTII(I,N)+DTVI(I,N)
           ELSEIF(YVAR.EQ.'Qi+v') THEN
            YDATA(I,NDT) = DQII(I,N)+DQVI(I,N)
C
           ELSEIF(YVAR.EQ.'CLsig') THEN
            YDATA(I,NDT) = CLR(I,N)*BLDS*CHR(I,N)/(2.0*PI*YRC(I,N))
           ELSEIF(YVAR.EQ.'Pstat') THEN
            YDATA(I,NDT) = DPSI(I,N)
C
           ELSEIF(YVAR.EQ.'VXa') THEN
            YDATA(I,NDT) = VABS(1,I,N)
           ELSEIF(YVAR.EQ.'VRa') THEN
            YDATA(I,NDT) = VABS(2,I,N)
           ELSEIF(YVAR.EQ.'VMa') THEN
            YDATA(I,NDT) = SQRT(VABS(1,I,N)**2 + VABS(2,I,N)**2)
           ELSEIF(YVAR.EQ.'VTa') THEN
            YDATA(I,NDT) = VABS(3,I,N)
           ELSEIF(YVAR.EQ.'VVa') THEN
            VMA = SQRT(VABS(1,I,N)**2 + VABS(2,I,N)**2)
            YDATA(I,NDT) = SQRT(VMA**2 + VABS(3,I,N)**2)
           ELSEIF(YVAR.EQ.'ANGa') THEN
            VMA = SQRT(VABS(1,I,N)**2 + VABS(2,I,N)**2)
            YDATA(I,NDT) = ATAN2(VABS(3,I,N),VMA)/DTR
C
           ELSEIF(YVAR.EQ.'VXr') THEN
            YDATA(I,NDT) = VREL(1,I,N)
           ELSEIF(YVAR.EQ.'VRr') THEN
            YDATA(I,NDT) = VREL(2,I,N)
           ELSEIF(YVAR.EQ.'VMr') THEN
            YDATA(I,NDT) = SQRT(VREL(1,I,N)**2 + VREL(2,I,N)**2)
           ELSEIF(YVAR.EQ.'VTr') THEN
            YDATA(I,NDT) = VREL(3,I,N)
           ELSEIF(YVAR.EQ.'VVr') THEN
            VMR = SQRT(VREL(1,I,N)**2 + VREL(2,I,N)**2)
            YDATA(I,NDT) = SQRT(VMR**2 + VREL(3,I,N)**2)
           ELSEIF(YVAR.EQ.'ANGr') THEN
            VMR = SQRT(VREL(1,I,N)**2 + VREL(2,I,N)**2)
            YDATA(I,NDT) = ATAN2(VMR,-VREL(3,I,N))/DTR
C
           ELSEIF(YVAR.EQ.'VRr') THEN
            YDATA(I,NDT) = VREL(2,I,N)
           ELSEIF(YVAR.EQ.'VMr') THEN
            YDATA(I,NDT) = SQRT(VREL(1,I,N)**2 + VREL(2,I,N)**2)

           ELSEIF(YVAR.EQ.'VXb') THEN
            YDATA(I,NDT) = VREL(1,I,N)
           ELSEIF(YVAR.EQ.'VRb') THEN
            YDATA(I,NDT) = VREL(2,I,N)
           ELSEIF(YVAR.EQ.'VMb') THEN
            YDATA(I,NDT) = SQRT(VREL(1,I,N)**2 + VREL(2,I,N)**2)
           ELSEIF(YVAR.EQ.'VTb') THEN
            VTB = VREL(3,I,N) - 0.5*VIND(3,I,N)
            YDATA(I,NDT) = VTB
           ELSEIF(YVAR.EQ.'VVb') THEN
            VMB = SQRT(VREL(1,I,N)**2 + VREL(2,I,N)**2)
            VTB = VREL(3,I,N) - 0.5*VIND(3,I,N)
            YDATA(I,NDT) = SQRT(VMB**2 + VTB**2)
           ELSEIF(YVAR.EQ.'ANGb') THEN
            VTB = VREL(3,I,N) - 0.5*VIND(3,I,N)
            YDATA(I,NDT) = ATAN2(VREL(1,I,N),-VTB)/DTR
          ENDIF
        END DO
        NPTDATA(NDT) = NPT
C
      END DO
      NDATA = NDT
      GO TO 5
C    
      END


      SUBROUTINE PLOTXY2(NDIM,
     &                   NDATA,NPTDATA,
     &                   INFODATA,XDATA,YDATA,
     &                   XMINSPEC,XMAXSPEC,YMINSPEC,YMAXSPEC,
     &                   XVAR,YVAR,XLIM,YLIM,
     &                   LPLTLEGEND,LEGENDLABEL,VERSION,
     &                   TITLE1, TITLE2, LPLOT,
     &                   AR,CH,XMARG,YMARG,NCOLOR,SCRNFR,IDEV)
C
C--- Plot XY data in multiple data arrays
C
      PARAMETER ( ITMPX=500 )
C
      CHARACTER XVAR*(*), YVAR*(*)
      CHARACTER LINELABEL*80, LEGENDLABEL*(*),VERSION*(*)
      CHARACTER TITLE1*(*), TITLE2*(*)
      LOGICAL   LGRID, LPLTLEGEND, LPLOT
C
      DIMENSION  NPTDATA(*)
      DIMENSION  INFODATA(2,*)
      DIMENSION  XDATA(NDIM,*), YDATA(NDIM,*)
      DIMENSION  LEGENDLABEL(*)

      DIMENSION  XLIM(2), YLIM(2)
C
      DIMENSION  XX(ITMPX), YY(ITMPX), P1(3)
      DIMENSION  XCMX(ITMPX), XCMN(ITMPX), XCAV(ITMPX)
      DIMENSION  YCMX(ITMPX), YCMN(ITMPX), YCAV(ITMPX)
      DIMENSION  IFLW(ITMPX)
C
      COMMON / PLT2DATA / XOF2D, YOF2D, XSF2D, YSF2D
C
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
C--- Scaling functions to plot coordinates 
      XXMOD(XXX) = (XXX - XOF2D)*XSF2D
      YYMOD(YYY) = (YYY - YOF2D)*YSF2D
C
C
      IF(NDATA.LE.0) RETURN
C
C---- Character and symbol sizing
      CH2 = 0.95*CH
      CH3 = 0.85*CH2
      SH  = 0.7 *CH2
      CHV = 0.75*CH2
C---- local plot aspect ratio
      PAR = 0.90*AR
C---- number of grid intervals per axis annotation interval
      NGR = 2
C---- origin location / size
      XORG = 0.16
      YORG = 0.10
C---- Defaults
      LGRID = .TRUE.
      LPLTLEGEND = .TRUE.
      SIZE = 9.0
C
C
C***************************************************
C---- Plot the input data arrays
C
C--- Check for plot limits on selected data
        XPAV = -1.0E10
        XPMX = -1.0E10
        XPMN =  1.0E10
        YPAV = -1.0E10
        YPMX = -1.0E10
        YPMN =  1.0E10
        NSEL = 0
        NCHLEGEND = 0
        DO NC = 1, NDATA
           NSEL = NSEL + 1
C
           LINELABEL = LEGENDLABEL(NC)
           CALL STRIP(LINELABEL,NLBL)
           NCHLEGEND = MAX(NCHLEGEND,NLBL)
C--- Get case #s
           IFLW(NSEL) = INFODATA(1,NC)
           NP = NPTDATA(NC)
           XAV =  0.0
           XMX = -1.0E10
           XMN =  1.0E10
           YAV =  0.0
           YMX = -1.0E10
           YMN =  1.0E10
           DO I = 1, NP
             XMX = MAX(XMX,XDATA(I,NC))
             XMN = MIN(XMN,XDATA(I,NC))
             XAV = XAV + XDATA(I,NC)
             NAV = NAV + 1
             YMX = MAX(YMX,YDATA(I,NC))
             YMN = MIN(YMN,YDATA(I,NC))
             YAV = YAV + YDATA(I,NC)
           END DO
           XAV = XAV / FLOAT(NP)
           YAV = YAV / FLOAT(NP)
           XCMX(NC) = XMX   
           XCMN(NC) = XMN   
           YCMX(NC) = YMX   
           YCMN(NC) = YMN   
           XCAV(NC) = XAV   
           YCAV(NC) = YAV   
C
           XPMX = MAX(XPMX,XMX)
           XPMN = MIN(XPMN,XMN)
           YPMX = MAX(YPMX,YMX)
           YPMN = MIN(YPMN,YMN)
           XPAV = MAX(XPAV,XAV)
           YPAV = MAX(YPAV,YAV)
        END DO
C
        XMX = XPMX
        XMN = XPMN
        YMX = YPMX
        YMN = YPMN
C
        IF(NSEL.LE.0) THEN
         WRITE(*,*) 'No plot data selected!'
         RETURN
        ENDIF
C
C--- Adjustments to min/max values for specific quantities
        IF(YVAR.EQ.'CF') THEN  
          YMN = 0.0
          YMX = YPAV
        ENDIF
        IF(YVAR.EQ.'H') THEN  
          YMN = 0.0
        ENDIF
C
C--- Adjustments to min/max values for specific quantities
        IF(XMINSPEC.NE.999.)  XMN = XMINSPEC
        IF(XMAXSPEC.NE.999.)  XMX = XMAXSPEC
        IF(YMINSPEC.NE.999.)  YMN = YMINSPEC
        IF(YMAXSPEC.NE.999.)  YMX = YMAXSPEC
C
        IF(YMX.EQ.YMN) THEN
          YMX = YMN + 1.0
        ENDIF
C
C--- Check limits imposed (-999. is autorange)
        IF(XLIM(1).NE.-999.) XMN = XLIM(1)
        IF(XLIM(2).NE.-999.) XMX = XLIM(2)
        IF(YLIM(1).NE.-999.) YMN = YLIM(1)
        IF(YLIM(2).NE.-999.) YMX = YLIM(2)
C
        IF(ABS(XMX-XMN).LT.1.0E-10) THEN
         XMX = XMN + 1.0
        ENDIF
        IF(ABS(YMX-YMN).LT.1.0E-10) THEN
         YMX = YMN + 1.0
        ENDIF
C
        IF(XMX.LT.XMN) THEN
         TMP = XMN
         XMN = XMX
         XMX = TMP
         IF(XMX.EQ.XMN) XMX = XMN + 1.0
        ENDIF
        IF(YMX.LT.YMN) THEN
         TMP = YMN
         YMN = YMX
         YMX = TMP
         IF(YMX.EQ.YMN) YMX = YMN + 1.0
        ENDIF
C
C--- Scale the axes to plot
        CALL AXISADJ(XMN,XMX,XSPAN,DELTAX,NXTICS)
        CALL AXISADJ(YMN,YMX,YSPAN,DELTAY,NYTICS)
C
        XLIM(1) = XMN
        XLIM(2) = XMX
        YLIM(1) = YMN
        YLIM(2) = YMX
C
C---- set plot offsets and scaling factors
        XAXISLEN = 0.9
        YAXISLEN = 0.9*PAR
C
        XOF2D  = XMN
        XSF2D  = XAXISLEN / XSPAN
        YOF2D  = YMN
        YSF2D  = YAXISLEN / YSPAN
C
        XAXLEN = XSPAN*XSF2D
        YAXLEN = YSPAN*YSF2D
C
        XLEGND  = 0.98*XAXLEN - FLOAT(NCHLEGEND)*CH3
        YLEGND  = YAXLEN - 0.6*CH3 + 2.2*CH3*(NDATA-1)
C
        IF(YVAR.EQ.'D') THEN
          XLEGND = 0.1*XAXLEN
        ENDIF
C
        Y0     = YMN
        YSTRT  = YMN
        YLAB   = YMX
C
        DXANN  = DELTAX*XSF2D
        DYANN  = DELTAY*YSF2D
C
C
C------ open window and plot axes for current scale/offset
        CALL GETCOLOR(ICOL0)
        CALL NEWCOLORNAME('BLACK')
C
        IF (LPLOT) CALL PLEND
        CALL PLOPEN(SCRNFR,0,IDEV)
        LPLOT = .TRUE. 
        CALL PLOTABS(XORG*SIZE+XMARG,YORG*SIZE+YMARG,-3)
        CALL NEWPEN(2)
C
        CALL NEWFACTOR(SIZE)
C
C--- Title
        CALL NEWCOLORNAME('green')
        XPLT = 0.0
        YPLT = YAXLEN + 5.0*CH2
        CALL NEWPEN(2)
        CALL PLCHAR(XPLT,YPLT,CH2,TITLE1,0.0,-1)
        YPLT = YAXLEN + 2.1*CH2
        CALL PLCHAR(XPLT,YPLT,CH2,TITLE2,0.0,-1)
C
C--- X axis with first annotation suppressed
        CALL NEWCOLORNAME('black')
        CALL XAXIS2(XXMOD(XMN),YYMOD(Y0),
     &             XAXLEN,DXANN,XMN,DELTAX,1,CH2,-2)
        CALL STRIP(XVAR,NXCH)
        XPLT = XXMOD(XMX)-2.5*DXANN-0.5*CH*NXCH
        YPLT = YYMOD(Y0)-0.5*FLOAT(LEN(XVAR))*CH
        CALL PLCHAR(XPLT,YPLT,CH2,XVAR,0.0,NXCH)
C
C--- Y axis
        CALL YAXIS(XXMOD(XMN),YYMOD(YSTRT),
     &             YAXLEN,DYANN,YSTRT,DELTAY,CH2,-2)
        CALL STRIP(YVAR,NYCH)
        XPLT = XXMOD(XMN)  - CH*FLOAT(NYCH+1)
        YPLT = YYMOD(YLAB) - 0.5*CH
        IF(NYTICS.LE.3) THEN
           YPLT = YPLT - 0.5*DYANN
          ELSE
           YPLT = YPLT - 1.5*DYANN
        ENDIF
        CALL PLCHAR(XPLT,YPLT,CH2,YVAR,0.0,NYCH)
C
C--- Background grid
        CALL NEWCOLORNAME('cyan')
        IF(LGRID) THEN
         NXGR = NGR * INT(XAXISLEN/DXANN + 0.001)
         NYGR = NGR * INT(YAXISLEN/DYANN + 0.001)
         DXG = DXANN / FLOAT(NGR)
         DYG = DYANN / FLOAT(NGR)
         CALL NEWPEN(1)
         CALL PLGRID(0.0,0.0, NXGR,DXG, NYGR,DYG, LMASK2 )
        ENDIF
C
C--- Plot the data as line segments
        CALL NEWCOLORNAME('black')
        NCC = 0
        DO NC = 1, NDATA
          IF(NPTDATA(NC).GT.0) THEN
C
            NCC = NCC + 1
            NCMOD = NCOLOR/10
            IC = 10*MOD(NCC-1,NCMOD) + 1
            CALL NEWCOLOR(1-IC)
C
            NP = NPTDATA(NC)
            NPC = NP/2
            ID = INFODATA(1,NC)

              CHORD = 1.0
              XLE   = 0.0
C
C--- Plot the point pairs on the stored line segments
            DO I = 1, NP-1
               X1 = (XDATA(I,  NC)-XLE)/CHORD
               X2 = (XDATA(I+1,NC)-XLE)/CHORD
               XX1 = XXMOD(X1)
               XX2 = XXMOD(X2)
               YY1 = YYMOD(YDATA(I  ,NC))
               YY2 = YYMOD(YDATA(I+1,NC))
C 
               CALL PLOT(XX1,YY1,3)
               CALL PLOT(XX2,YY2,2)
C 
               CALL PLSYMB(XX1,YY1,SH,NC,0.0,0)
               CALL PLSYMB(XX2,YY2,SH,NC,0.0,0)
C               IF(I.EQ.NPC) THEN
C                CALL PLNUMB(XX2+CH2,YY2+CH2,CH2,FLOAT(NC),0.0,-1)
C               ENDIF
            END DO   
C
C--- Legend for curves
           IF(LPLTLEGEND) THEN
            LINELABEL = LEGENDLABEL(NC)
            CALL STRIP(LINELABEL,NLBL)
            XPLT = XLEGND
            YPLT = YLEGND - FLOAT(NCC-1)*2.2*CH3
            CALL PLSYMB(XPLT,YPLT+2.5*SH,SH,NC,0.0,0)
            CALL NEWCOLORNAME('BLACK')
            CALL PLCHAR(999.,YPLT+2.0*SH,CH3,'  ',0.0,2)
            CALL PLCHAR(999.,YPLT+2.0*SH,CH3,LINELABEL,0.0,NLBL)
           ENDIF
C
          ENDIF
C
        END DO
C
C--- Version info
C
      CALL STRIP(VERSION,NVERSION)
      XVERS = XAXLEN + 2.5*CHV
      YVERS = YAXLEN - FLOAT(NVERSION + 7)*CHV
C
      CALL NEWCOLORNAME('cyan')
      CALL PLCHAR(XVERS,YVERS,CHV,'DFDC v',90.,6)
      CALL PLCHAR(XVERS,999.,CHV,VERSION,90.,NVERSION)
C
      CALL PLFLUSH
C
      CALL NEWCOLOR(ICOL0)
      RETURN
      END  ! PLOTXY2
   


      SUBROUTINE OFFSET2D(XLIM,YLIM)
C---- Get zoom box from user mouse selection on 2D plot 
C
      DIMENSION XLIM(2), YLIM(2)
      CHARACTER*1 CHKEY
      COMMON / PLT2DATA / XOF2D, YOF2D, XSF2D, YSF2D
C
      SH = 2.0
C
      WRITE(*,*)
      WRITE(*,*) 'Mark off corners of blowup area'
      WRITE(*,*) '(2 spaces default to current area)'       
      CALL GETCURSORXY(XX1,YY1,CHKEY)
      CALL PLSYMB(XX1,YY1,SH,3,0.0,0)
      CALL GETCURSORXY(XX2,YY2,CHKEY)
      CALL PLSYMB(XX2,YY2,SH,3,0.0,0)
      IF(ABS(XX1-XX2).LT.0.01 .AND. ABS(YY1-YY2).LT.0.01) RETURN      
C
      XOF = MIN(XX1,XX2)/XSF2D + XOF2D
      YOF = MIN(YY1,YY2)/YSF2D + YOF2D
      XDIF = ABS(XX2 - XX1)/XSF2D
      YDIF = ABS(YY2 - YY1)/YSF2D   
      IF(XDIF.NE.0.0) THEN
        XLIM(1) = XOF
        XLIM(2) = XOF + XDIF
      ENDIF
      IF(YDIF.NE.0.0) THEN
        YLIM(1) = YOF
        YLIM(2) = YOF + YDIF
      ENDIF
C
      RETURN        
      END 




      SUBROUTINE XAXIS2(X1,Y1,XAXT,DXANN,FANN,DANN,IFLAG,CHT,NDIG)
C
C---XAXIS2 differs from libPlt XAXIS by having a flag to suppress both
C   end annotations rather than the zero annotation.
C.......................................................
C     X1,Y1  starting point of x axis
C     XAXT   length of x axis 
C     DXANN  distance between annotations
C     FANN   first annotation value
C     DANN   delta annotation value
C     IFLAG  flag to suppress end annotations 
C            = 0    all annotations
C            = 1    suppress first annotation
C            = 2    suppress last  annotation
C            = 3    suppress first and last annotations
C     CHT    character height   ( - = annotation above axis)
C     NDIG   number of digits to right of decimal point
C            = -1   no decimal point
C            = -2   number of digits determined internally
C.......................................................
C
      XAX = ABS(XAXT)
      IF(XAX.LE.0) RETURN
      CH = ABS(CHT)
C
C--- determine # of digits to use for annotations
      IF(NDIG.LE.-2) THEN
        ND = 1 - MAX( 0 , INT(LOG10(DANN)) )
        IF(DANN*10**ND - AINT(DANN*10**ND+0.01) .GT. 0.01) ND = ND + 1
        IF(DANN*10**ND - AINT(DANN*10**ND+0.01) .GT. 0.01) ND = ND + 1
      ELSE
        ND = NDIG
      ENDIF
      NDG = MAX(-1,ND)
ccc      write(*,*) 'xaxis2 dann,nd,ndig,ndg ',dann,nd,ndig,ndg
C
C---- x-axis
      CALL PLOT(X1,Y1,3)
      CALL PLOT(X1+XAX,Y1,2)
      NANN = 1 + IFIX(XAX/DXANN + 0.1)
ccc      write(*,*) 'nann ',nann
C
C---- annotate x-axis
      DO 10 NT=1, NANN
        XT = X1 + DXANN*FLOAT(NT-1)
C---- skip annotation for first or last annotation position as given by IFLG
        IF(MOD(IFLAG,2).EQ.1 .AND. NT.EQ.1)    GO TO 10
        IF(IFLAG.GT.1        .AND. NT.EQ.NANN) GO TO 10
C
        CALL PLOT(XT,Y1-0.2*CH,3)
        CALL PLOT(XT,Y1+0.2*CH,2)
        RN = FANN + DANN*FLOAT(NT-1)
        GRN = 0.
        IF(RN.NE.0.0) GRN = ALOG10(ABS(RN)+0.5/10.0**ND)
        GRN = MAX(GRN,0.0)
        NABC = INT(GRN) + 2 + ND
        WIDTH = 0.95*CH*FLOAT(NABC)
        IF(RN.LT.0.0) WIDTH = WIDTH + CH
        XNUM = XT - 0.5*WIDTH
        YNUM = Y1 - 2.1*CH
        IF(CHT.LT.0.0) YNUM = Y1 + 0.9*CH
C
        CALL PLNUMB(XNUM,YNUM,CH,RN,0.0,NDG)
   10 CONTINUE
C
      RETURN
      END ! XAXIS2
 
C
C---Esotec stuff ------------------------------------------------
C
      SUBROUTINE GETPNAME(PNAME,NAME,NR)
      CHARACTER*80 NAME,PNAME
      CHARACTER TNUM(9)*1
      DATA TNUM / '1','2','3','4','5','6','7','8','9' /
C
      L=LEN(NAME)
      DO J=L,1,-1
         IF(NAME(J:J).NE.' ')THEN
            L=J
            GO TO 10
         END IF
      ENDDO
 10   PNAME = NAME(1:L) // '-Dsk' // TNUM(NR)
C
      RETURN
      END  ! GETPNAME
C
C--------------------------------------------------------------------
C

      SUBROUTINE GETROTG(N)
      INCLUDE 'DFDC.INC'
      DIMENSION T1(IRX), T1S(IRX),
     &          T2(IRX), T2S(IRX)
C
      DO I=1,NRC
         T1(I) = CHR(I,N)
         T2(I) = BETAR(I,N)
      ENDDO
C
      CALL SEGSPL(T1,T1S,YRC(1,N),NRC)
      CALL SEGSPL(T2,T2S,YRC(1,N),NRC)
C
      DO IP=1,NRP
         CHRP  (IP,N) = SEVAL(YRP(IP,N),T1,T1S,YRC(1,N),NRC)
         BETARP(IP,N) = SEVAL(YRP(IP,N),T2,T2S,YRC(1,N),NRC)
      ENDDO
C
      RETURN
      END   ! GETROTG


























!

