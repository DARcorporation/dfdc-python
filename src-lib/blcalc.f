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

      SUBROUTINE BLOPER
C---- BL analysis
C---- Calculate strip BL's on surfaces
C
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      CHARACTER*4 COMAND, COMOLD
      CHARACTER*1 CIND
C
      CHARACTER*16 PROMPT
      CHARACTER*80 LINE
      CHARACTER*128 COMARG, ARGOLD
C
      LOGICAL LERROR, LXRELIMIT, LYRELIMIT
      LOGICAL LEGEND2D, LPLTMODE
C
      CHARACTER*6 XVAR, YVAR
      CHARACTER DATATYPE*40
      CHARACTER DATALBL*20, LEGENDDATA*80
      CHARACTER FNAME*128
      CHARACTER*80 TITLE1, TITLE2
C
      DIMENSION IINPUT(10), RINPUT(10)
      DIMENSION XLIM(2),    YLIM(2)
      DIMENSION XTRAN(2)
C
C--- Temporary data for plotting 
      PARAMETER (IXP=500,NXP=10)
      DIMENSION  NPTDATA(NXP)
      DIMENSION  INFODATA(2,NXP)
      DIMENSION  XDATA(IXP,NXP), YDATA(IXP,NXP)
      DIMENSION  LEGENDDATA(NXP)
C
      DATA XVAR, YVAR / 'X', 'Q' /
C
      PAR = 0.6
      LUTEMP = 3
      LXRELIMIT = .TRUE.
      LYRELIMIT = .TRUE.
      XMINSPEC = 999.
      XMAXSPEC = 999.
      YMINSPEC = 999.
      YMAXSPEC = 999.
      LEGEND2D = .TRUE.
      LPLTMODE = .FALSE.
C
      XTRAN(1) = 1.0
      XTRAN(2) = 1.0
      IF(AMPMAX.LE.0.0) THEN
        AMPMAX = 9.0
      ENDIF
C
 1    FORMAT(A)
C
      COMAND = '****'
      COMARG = ' '
C
      IEL = 0
C
      TITLE1 = 'DFDC Wall BL Data'
      UREN = RHO/RMU
C
C---- Display options 
c 10   WRITE(*,*) ' '
 10   URENQREF = QREF*RHO/RMU
C
      IF(LPLTMODE) THEN
        WRITE(*,22) XVAR,YVAR
        PROMPT = '.BLPL^'
      ELSE
        WRITE(*,20) URENQREF,XTRAN(1),XTRAN(2),AMPMAX 
        PROMPT = '.BL  ^'
      ENDIF
C
 20   FORMAT(/' ================================================'
     &       /'   BL  BL calculation',
     &       /'   RE  RN/m for V=Qref)        currently:  ',G12.6,
     &       /'   TR  transition points       currently: ',2(F6.2),
     &       /'   N   N-crit for transition   currently: ',F6.2,
     &       /'   P   plot  BL data',
     &       /'   W   write BL data to file')
cc     &       /' Select option (or <return>):  ', $)
C
 22   FORMAT(/' ================================================'
     &       /'   A   abscissa select   currently: ',A,
     &       /'   O   ordinate select   currently: ',A,
     &       /'   L   limits for plot',
     &       /'   R   reset plot limits',
     &       /'   AN  annotate plot'
     &       /'   H   hardcopy plot')
cc     &       /' Select option (or <return>):  ', $)
C
C====================================================================
C
 501  IF(IEL.GT.0 .AND. NBEL.GT.1) THEN
C------ add current element # to prompt
        IF    (IEL.LT. 10) THEN
         PROMPT(6:7) = '(' // PNUM(IEL)(2:2)
         K = 8
        ENDIF
        KEL = IEL
        PROMPT(K:K+1) = ')^'
      ENDIF
C
C---- Get command line from user
      CALL ASKC(PROMPT,COMAND,COMARG)
      CIND = COMAND
C---- extract numbers (if any) from command argument string
      DO I=1, 10
        IINPUT(I) = 0
        RINPUT(I) = 0.0
      ENDDO
      NINPUT = 0
      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
      NINPUT = 0
      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C
C---- Check command for element index
      IF(INDEX('0123456789',CIND) .NE. 0) THEN
        READ(COMAND,*,ERR=10) IE
        IF(IE.GE.0 .OR. IE.LE.NBEL) THEN
          IEL = IE
C---- If we are plotting already do a new plot...
          IF(LPLTMODE) GO TO 100
          GO TO 501
        ENDIF
      ENDIF
C
 505  IF(.NOT.LPLTMODE) THEN
C
C---- Basic BL analysis commands
       IF(COMAND.EQ.' ') THEN
        RETURN
C
       ELSEIF(COMAND.EQ.'?') THEN
         GO TO 10
C
C---- Calculate BL's on each element
       ELSEIF(COMAND.EQ.'BL') THEN
        IF(UREN.LE.0.0) THEN
          WRITE(*,*) '*** Set unit Reynolds number...'
          GO TO 10
        ENDIF
        CALL NEW_BL(XTRAN)
C
C---- Specify unit Reynold's number
       ELSEIF(COMAND.EQ.'RE') THEN
        IF(NINPUT.GE.1) THEN
          UREN = RINPUT(1)
        ELSE
         RINPUT(1) = UREN
         CALL ASKR('Enter unit RN^',RINPUT)
         UREN = RINPUT(1)
        ENDIF
C
C---- Specify N-crit for e^n transition 
       ELSEIF(COMAND.EQ.'N') THEN
        IF(NINPUT.GE.1) THEN
          AMPMAX = RINPUT(1)
        ELSE
         RINPUT(1) = AMPMAX
         CALL ASKR('Enter Ncrit for transition^',RINPUT)
         AMPMAX = RINPUT(1)
        ENDIF
C
C---- Specify Transition points
       ELSEIF(COMAND.EQ.'TR  ') THEN
        IF(NINPUT.GE.2) THEN
          XTRAN(1) = RINPUT(1)
          XTRAN(2) = RINPUT(2)
        ELSE
         RINPUT(1) = XTRAN(1)
         RINPUT(2) = XTRAN(2)
         CALL ASKRN('Enter S/Stot for fixed transition^',RINPUT,2)
         XTRAN(1) = RINPUT(1)
         XTRAN(2) = RINPUT(2)
        ENDIF
C
C---- Write out a file with the BL data
       ELSEIF(COMAND.EQ.'W') THEN
        NNEL = 0
        CALL WRITEBLDATA(NNEL)
C
C---- Plot the BL data
       ELSEIF(COMAND.EQ.'P' .OR. COMAND.EQ.'PLOT') THEN
        LPLTMODE = .TRUE.
        GO TO 100
C
       ENDIF
       GO TO 501
C
      ELSE
C
C***************************************************
C---- Plotting commands
       IF(COMAND.EQ.' ') THEN
         LPLTMODE = .FALSE. 
         GO TO 10
C
       ELSEIF(COMAND.EQ.'?') THEN
         GO TO 10
C
C---- Make hardcopy
       ELSEIF(COMAND.EQ.'H' .OR. COMAND.EQ.'HARD') THEN
        CALL REPLOT(IDEVH)
        GO TO 501
C
C---- Annotate plot
       ELSEIF(COMAND.EQ.'AN') THEN
        CALL ANNOT(CH)
        GO TO 501
C
C---- Toggle legend flag
       ELSEIF(COMAND.EQ.'LG') THEN
        LEGEND2D = .NOT.LEGEND2D
C
C---- Reset plot limits
       ELSEIF(COMAND.EQ.'R') THEN
        LXRELIMIT = .TRUE.
        LYRELIMIT = .TRUE.
C
C---- Set plot limits
       ELSEIF(COMAND.EQ.'L') THEN
        WRITE(*,*) 'Limits for plot:'
C
        WRITE(*,32) XLIM
 32     FORMAT('Xlimits  xmin: ',F12.6,' xmax: ',F12.6)
        IF(NINPUT.GE.2) THEN
          XLIM(1) = RINPUT(1)
          XLIM(2) = RINPUT(2)
        ELSE
         RINPUT(1) = XLIM(1)
         RINPUT(2) = XLIM(2)
         CALL ASKRN('Enter xmin, xmax for plot^',RINPUT,2)
         XLIM(1) = RINPUT(1)
         XLIM(2) = RINPUT(2)
        ENDIF
        LXRELIMIT = .FALSE.
C
        WRITE(*,34) YLIM
 34     FORMAT('Ylimits  ymin: ',F12.6,' ymax: ',F12.6)
        IF(NINPUT.GE.2) THEN
          YLIM(1) = RINPUT(1)
          YLIM(2) = RINPUT(2)
        ELSE
         RINPUT(1) = YLIM(1)
         RINPUT(2) = YLIM(2)
         CALL ASKRN('Enter ymin, ymax for plot^',RINPUT,2)
         YLIM(1) = RINPUT(1)
         YLIM(2) = RINPUT(2)
        ENDIF
        LYRELIMIT = .FALSE.
C
C---- Change Abscissa 
       ELSEIF(COMAND.EQ.'A') THEN
 40     WRITE(*,*) ' '
        WRITE(*,45) XVAR
 45     FORMAT(/' ================================================'
     &       /' Abscissa is: ',A,
     &       /'   X, Y, S'
     &      //' Select option (or <return>):  ', $)
        READ(*,1,ERR=40) COMAND
        CALL LC2UC(COMAND)
        IF(COMAND.EQ.' ') THEN
          GO TO 10
         ELSE IF(COMAND.EQ.'X') THEN
          XVAR = 'X'
         ELSE IF(COMAND.EQ.'Y') THEN
          XVAR = 'Y'
         ELSE IF(COMAND.EQ.'S') THEN
          XVAR = 'S'
        ENDIF
        LXRELIMIT = .TRUE.
C
C---- Change Ordinate
       ELSE IF(COMAND.EQ.'O') THEN
 50     WRITE(*,*) ' '
        WRITE(*,55) YVAR
 55     FORMAT(/' ================================================'
     &       /' Ordinate is: ',A,
     &       /'   X, Y, CP, Q, DIV, DS, TH, H, CF, N'
     &      //' Select option (or <return>):  ', $)
        READ(*,1,ERR=50) COMAND
        CALL LC2UC(COMAND)
        IF(COMAND.EQ.' ') THEN
          GO TO 10
         ELSE IF(COMAND.EQ.'X') THEN
          YVAR = 'X'
         ELSE IF(COMAND.EQ.'Y') THEN
          YVAR = 'Y'
         ELSE IF(COMAND.EQ.'S') THEN
          YVAR = 'S'
         ELSE IF(COMAND.EQ.'Q') THEN
          YVAR = 'Q'
         ELSE IF(COMAND.EQ.'CP') THEN
          YVAR = 'CP'
         ELSE IF(COMAND.EQ.'DIV') THEN
          YVAR = 'DIV'
         ELSE IF(COMAND.EQ.'DS') THEN
          YVAR = 'DS'
         ELSE IF(COMAND.EQ.'TH') THEN
          YVAR = 'TH'
         ELSE IF(COMAND.EQ.'H') THEN
          YVAR = 'H'
         ELSE IF(COMAND.EQ.'CF') THEN
          YVAR = 'CF'
         ELSE IF(COMAND.EQ.'N') THEN
          YVAR = 'N'
        ENDIF
        LYRELIMIT = .TRUE.
C
       ENDIF
C
      ENDIF
C
C***************************************************
C---- Plot the selected BL data
 100  IF(LXRELIMIT) THEN
       XLIM(1) = -999.
       XLIM(2) = -999.
      ENDIF
      IF(LYRELIMIT) THEN
       YLIM(1) = -999.
       YLIM(2) = -999.
      ENDIF
C--- Put element data into plotting arrays
 110  IF(IEL.EQ.0) THEN
        IE1 = 1
        IE2 = NBEL
      ELSE
        IE1 = IEL
        IE2 = IEL
      ENDIF
C
      NDT = 0
      DO IE = IE1, IE2
C
       DO IS = 1, 2
C--- Check for valid BL for element side BL's
        IF(ICBLBEG(IS,IE).GT.0) THEN
         NDT = NDT + 1
         ICB = ICBLBEG(IS,IE)
         ICE = ICBLEND(IS,IE)
         IDEL = 1
         IF(ICE.LT.ICB) IDEL = -1
         NPT = ABS(ICE - ICB) + 1
         INFODATA(1,NDT) = 2*(IE-1) + IS
C--- Legend for each curve 
         WRITE(LEGENDDATA(NDT),111) 'BL elem', IE,IS,QINF,QREF,ALTH
 111     FORMAT(' ',A,' ',I2,' side ',I2,' @ ',3(F8.2,','))
C
         NAMP = 0
         I    = 0
         DO IC = ICB, ICE, IDEL
          I = I + 1
C--- Install X data
          IF(XVAR.EQ.'X') THEN
           XDATA(I,NDT) = BLDATA(1,IC)
          ELSEIF(XVAR.EQ.'Y') THEN
           XDATA(I,NDT) = BLDATA(2,IC)
          ELSEIF(XVAR.EQ.'S') THEN
           XDATA(I,NDT) = BLDATA(3,IC)
          ENDIF
C--- Install Y data
          IF(YVAR.EQ.'X') THEN
           YDATA(I,NDT) = BLDATA(1,IC)
          ELSEIF(YVAR.EQ.'Y') THEN
           YDATA(I,NDT) = BLDATA(2,IC)
          ELSEIF(YVAR.EQ.'S') THEN
           YDATA(I,NDT) = BLDATA(3,IC)
          ELSEIF(YVAR.EQ.'Q') THEN
           YDATA(I,NDT) = BLDATA(4,IC)
          ELSEIF(YVAR.EQ.'CP') THEN
           YDATA(I,NDT) = BLDATA(5,IC)
          ELSEIF(YVAR.EQ.'Q') THEN
           YDATA(I,NDT) = BLDATA(6,IC)
          ELSEIF(YVAR.EQ.'DIV') THEN
           YDATA(I,NDT) = BLDATA(6,IC)
          ELSEIF(YVAR.EQ.'DS') THEN
           YDATA(I,NDT) = BLDATA(7,IC)
          ELSEIF(YVAR.EQ.'TH') THEN
           YDATA(I,NDT) = BLDATA(8,IC)
          ELSEIF(YVAR.EQ.'H') THEN
           YDATA(I,NDT) = BLDATA(9,IC)
          ELSEIF(YVAR.EQ.'CF') THEN
           YDATA(I,NDT) = BLDATA(10,IC)
          ELSEIF(YVAR.EQ.'N') THEN
           YDATA(I,NDT) = BLDATA(11,IC)
           IF(YDATA(I,NDT).GT.0.001) NAMP = I
          ENDIF
         END DO
         NPTDATA(NDT) = NPT
         IF(NAMP.NE.0) NPTDATA(NDT) = NAMP
        ENDIF
       END DO
C
      END DO
      NDATA = NDT
C
C---- Plot data
 120  IF(NDATA.GT.0) THEN 
       CALL PLOTXY2(IXP,
     &              NDATA,NPTDATA,
     &              INFODATA,XDATA,YDATA,
     &              XMINSPEC,XMAXSPEC,YMINSPEC,YMAXSPEC,
     &              XVAR,YVAR,XLIM,YLIM,
     &              LEGEND2D,LEGENDDATA,SVERSION,
     &              TITLE1, NAME, LPLOT,
     &              PAR,CSIZE,XMARG,YMARG,NCOLOR,SCRNFR,IDEV)
       LPLOT = .TRUE.
       CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
      ENDIF
      GO TO 10
C
      END

      
 
      SUBROUTINE WRITEBLDATA(IIEL)
      INCLUDE 'DFDC.INC'
      CHARACTER*80 FNAME
C
 1000 FORMAT(A)
C
      IF(.NOT.LVISC) THEN
        WRITE(*,*) 'BL data not defined...'
        RETURN
      ENDIF
C
      CALL ASKS('Enter BL data output filename^',FNAME)
      IF(FNAME.EQ.' ') THEN
        LU = 6
      ELSE
        LU = 3
        OPEN(UNIT=LU,FILE=FNAME,STATUS='UNKNOWN',ERR=5)
        GO TO 10
C
    5   WRITE(*,*) '***   File open error   ***'
        CLOSE(UNIT=LU)
        RETURN
      ENDIF
C
C--- Write header and titles for surface flow and BL data
 10   WRITE(LU,1000) NAME
      WRITE(LU,200)
C
      IE = IIEL
      IF(IE.LE.0) THEN
       IE1 = 1
       IE2 = NBEL
      ELSE
       IE1 = IE
       IE2 = IE
      ENDIF
C
C---- Dump BL's for elements
      DO IEL = IE1, IE2
C
       WRITE(LU,205) IEL
C
C--- First BL side
       IS = 1
       IC1 = ICBLBEG(IS,IEL)
       IC2 = ICBLEND(IS,IEL)
       ICS = ICSTG(IEL)
       IDEL = 1
       IF(IC2.LT.IC1) IDEL = -1
C
       N = ABS(IC2-IC1) + 1
       IF(N.GE.2) THEN
        WRITE(LU,210) IS
        DO IC = IC1, IC2, IDEL
         WRITE(LU,230) BLDATA(1,IC),  BLDATA(2,IC),
     &                 BLDATA(3,IC),  BLDATA(4,IC),
     &                 BLDATA(5,IC),  BLDATA(7,IC),
     &                 BLDATA(8,IC),  BLDATA(9,IC),
     &                 BLDATA(10,IC), BLDATA(11,IC)
        END DO
       ENDIF
C
C--- Second BL side
       IS = 2
       IC1 = ICBLBEG(IS,IEL)
       IC2 = ICBLEND(IS,IEL)
       ICS = ICSTG(IEL)
       IDEL = 1
       IF(IC2.LT.IC1) IDEL = -1
C
       N = ABS(IC2-IC1) + 1
       IF(N.GE.2) THEN
        WRITE(LU,220) IS
        DO IC = IC1, IC2, IDEL
         WRITE(LU,230) BLDATA(1,IC),  BLDATA(2,IC),
     &                 BLDATA(3,IC),  BLDATA(4,IC),
     &                 BLDATA(5,IC),  BLDATA(7,IC),
     &                 BLDATA(8,IC),  BLDATA(9,IC),
     &                 BLDATA(10,IC), BLDATA(11,IC)
        END DO
       ENDIF
C
      END DO
C
C
 200  FORMAT('#Surface Flow and BL data')
 205  FORMAT(/'# Element :',I3)
C
 210  FORMAT(/'#Side :',I3/
     &        '#X',T13,' Y',T26,' s/S',T39,
     &        ' Uedg',T52,' Cp',T65,
     &        ' Dstr',T78,' Theta',T91,' H12',T104,
     &        ' Cf',T117,' Namp',T140)
C
 220  FORMAT(/'#Side :',I3/
     &        '#X',T13,' Y',T26,' s/S',T39,
     &        ' Uedg',T52,' Cp',T65,
     &        ' Dstr',T78,' Theta',T91,' H12',T104,
     &        ' Cf',T117,' Namp',T140)
C
 230  FORMAT(10(G12.6,1X))
C
      CLOSE(UNIT=3)
      RETURN
      END



      SUBROUTINE NEW_BL(STRAN)
C
      INCLUDE 'DFDC.INC'
C
      DIMENSION NPBL(2),      ICBL(IPX,2)
      DIMENSION XBL(IPX,2),   YBL(IPX,2), 
     &          QBL(IPX,2),   SBL(IPX,2),
     &          DSTAR(IPX,2), THETA(IPX,2),
     &          HK(IPX,2),    CF(IPX,2),   
     &          BDIV(IPX,2),  AMP(IPX,2)
      DIMENSION STRAN(2)
C
      DO IEL = 1, NEL
        CXVIS(IEL) = 0.0
C
        IF(NETYPE(IEL).EQ.0) THEN
C
          IC1 = ICFRST(IEL)
          IC2 = ICLAST(IEL)
          ICS = ICSTG(IEL)
C---- Check locaton of stag. point vs control point on stagnation panel 
          XCS = XC(ICS)
          YCS = YC(ICS)
          DXS = XSTG(IEL)-XCS
          DYS = YSTG(IEL)-YCS
          DOTP = -ANC(2,ICS)*DXS + ANC(1,ICS)*DYS
          IF(DOTP.GE.0.0) THEN
            S1 =  DOTP
            S2 = 0.5*DSC(ICS) - DOTP
            IC0 = ICS
          ELSE
            S1 = 0.5*DSC(ICS) + DOTP
            S2 = -DOTP 
            IC0 = ICS-1
          ENDIF
C
C---- Assemble BL from stagnation point backwards
          IS = 0
          IF(IC0-IC1.GT.2) THEN
            IS = IS + 1
            I  = 1
            XBL(I,IS) = XSTG(IEL)  
            YBL(I,IS) = YSTG(IEL)  
            SBL(I,IS) = 0.0
            QBL(I,IS) = 0.0
            ICBL(I,IS) = 0
            DO IC = IC0, IC1, -1
              I = I + 1
              XBL(I,IS) = XC(IC) 
              YBL(I,IS) = YC(IC) 
              DY =      YBL(I,IS) - YBL(I-1,IS)
              YA = 0.5*(YBL(I,IS) + YBL(I-1,IS))
              IF(IC.EQ.2 .AND. DOTP.GT.0.0) THEN
                DS = S1
              ELSE
                DS = S1 + 0.5*DSC(IC)
              ENDIF
              SBL(I,IS) = SBL(I-1,IS) + DS
              QBL(I,IS) = SQRT(QCR(1,IC)**2 + QCR(2,IC)**2)   
              DS = SBL(I,IS) - SBL(I-1,IS)
              BDIV(I,IS) = DY/DS / YA
              ICBL(I,IS) = IC
              S1 = 0.5*DSC(IC)
            END DO
            BDIV(1,IS) = BDIV(2,IS)
c            BDIV(1,IS) = 0.0
            NPBL(IS) = I
          ENDIF
C
C---- Assemble BL from stagnation point forwards
          IF(IC2-(IC0+1).GT.2) THEN 
            IS = IS + 1
            I  = 1
            XBL(I,IS) = XSTG(IEL)  
            YBL(I,IS) = YSTG(IEL)  
            SBL(I,IS) = 0.0
            QBL(I,IS) = 0.0
            ICBL(I,IS) = 0
            DO IC = IC0+1, IC2
              I = I + 1
              XBL(I,IS) = XC(IC) 
              YBL(I,IS) = YC(IC) 
              DY =      YBL(I,IS) - YBL(I-1,IS)
              YA = 0.5*(YBL(I,IS) + YBL(I-1,IS))
              IF(IC.EQ.2 .AND. DOTP.LT.0.0) THEN
                DS = S2
              ELSE
                DS = S2 + 0.5*DSC(IC)
              ENDIF
              SBL(I,IS) = SBL(I-1,IS) + DS
              QBL(I,IS) = SQRT(QCR(1,IC)**2 + QCR(2,IC)**2)   
              DS = SBL(I,IS) - SBL(I-1,IS)
              BDIV(I,IS) = DY/DS/ YA
              ICBL(I,IS) = IC
              S2 = 0.5*DSC(IC)
            END DO
            BDIV(1,IS) = BDIV(2,IS)
c            BDIV(1,IS) = 0.0
            NPBL(IS) = I
          ENDIF
          NBL = IS
C
C--- Calculate BL for each side, could be 1 BL or 2 BL's
          FX = 0.0
          DO IS = 1, NBL
            SMAX = SBL(NPBL(IS),IS)
            STR  = STRAN(IS)*SMAX
            CMACH = MACH
C--- Calculate the boundary layer
            CALL BL (CMACH, UREN, STR, AMPMAX, 
     &               NPBL(IS), SBL(1,IS), QBL(1,IS), BDIV(1,IS),
     &               DSTAR(1,IS), THETA(1,IS), HK(1,IS), 
     &               CF(1,IS), AMP(1,IS), 
     &               CDP, CFALL, STRN(IS,IEL), SSEP(IS,IEL))
C
            LVISC  = .TRUE.
C
C--- Put results back into BLDATA arrays using ICBL index mapping array
            DO I = 1, NPBL(IS)
              IC = ICBL(I,IS)
              IF(IC.GT.0) THEN
                BLDATA( 1,IC) = XBL(I,IS)
                BLDATA( 2,IC) = YBL(I,IS)
                BLDATA( 3,IC) = SBL(I,IS)
                BLDATA( 4,IC) = QBL(I,IS)
                BLDATA( 5,IC) = CPR(IC)
                BLDATA( 6,IC) = BDIV(I,IS)
                BLDATA( 7,IC) = DSTAR(I,IS)
                BLDATA( 8,IC) = THETA(I,IS)
                BLDATA( 9,IC) = HK(I,IS)
                BLDATA(10,IC) = CF(I,IS)
                BLDATA(11,IC) = AMP(I,IS)
              ENDIF
            END DO
C--- Save start/end IC for element/side index for plotting
            ICBLBEG(IS,IEL) = ICBL(2,IS)
            ICBLEND(IS,IEL) = ICBL(NPBL(IS),IS)
C
C--- Integrate skin friction to generate viscous forces
            DO I = 1, NPBL(IS)
              TAU  = 0.5*RHO*QBL(I,IS)**2 * CF(I,IS)
              IF(I.GT.1) THEN
                TAUA = 0.5*(TAU + TAUO)
                DX  =      XBL(I,IS) - XBL(I-1,IS)
                YA  = 0.5*(YBL(I,IS) + YBL(I-1,IS))
                DFX = 2.0*PI*YA * DX * TAUA
                FX  = FX + DFX
              ENDIF
              TAUO = TAU
            END DO
            QUREF = 0.5*RHO*QREF**2
            CXVIS(IEL) = CXVIS(IEL) + FX / QUREF
C
C--- Print the results
            WRITE(*,230) IEL, IS, 
     &                   CMACH, UREN, 
     &                   AMPMAX, STRN(IS,IEL), SSEP(IS,IEL),
     &                   CXVIS(IEL)
            DO I = 1, NPBL(IS)
              WRITE(*,240)  SBL(I,IS)/SMAX, 
cc     &                      SBL(I,IS), XBL(I,IS), QBL(I,IS), 
     &                      SBL(I,IS), XBL(I,IS), QBL(I,IS), 
     &                      DSTAR(I,IS), THETA(I,IS), HK(I,IS), 
     &                      CF(I,IS), AMP(I,IS)
            END DO
C
          END DO
C
        ENDIF
C--- Update element and total forces with viscous force
        CX(IEL)  = CX(IEL)  + CXVIS(IEL)
        CX(0)    = CX(0)    + CXVIS(IEL)
        CXVIS(0) = CXVIS(0) + CXVIS(IEL)
      END DO
C
  230 FORMAT (/' ** B.L. Analysis on element ',I3,' side ',I3/
     *        '  Mach  =',F12.5/
     *        '  Uren  =',G12.3,'   Ampmax =',F10.5/  
     *        '  Strn  =',F12.5,'   Ssep   =',F12.5/
     *        '  CXvis =',F12.5//
cc     *        '  Cdp  =',F12.6,'   Cftot  =',F12.6//
     *        3X,'S/Smax',3X,'S',7X,'X',7X,'Ue',6X,'Dstar',5X,'Theta',
cc     *        3X,'S/Smax',4X,'S',8X,'Ue',7X,'Dstar',6X,'Theta',
     *        5X,'H12  ',5X,'Cf   ',3X,'Amp')
  240 FORMAT (1X,F6.3,F8.3,F8.3,F8.3,2F10.5,F8.3,F10.5,F8.3)
cc  240 FORMAT (1X,F7.3,F9.3,F9.3,2F11.5,F9.3,F11.5,F9.3)
c  123456712345678912345678912345678901123456789011234567123456789011234567
c   S/Smax    S       Ue        Dstar      Theta    H12      Cf       Amp
C
      RETURN
      END



      SUBROUTINE BL (MACH, UREN, SFIX, AMPMAX, NP, S, UE, UDIV,
     &               DSTAR, THETA, H12, CFE, AMPF, CDP, CFTOT, 
     &               STRN, SSEP)
C***********************************************************************
C   Copyright  May 1993    Harold Youngren - all rights reserved
C
C   Compressible 2D Integral BL Computation
C     This code computes a compressible boundary layer given an input
C     velocity distribution (U/Uinf) and a free stream Mach number. 
C     The method is valid for incompressible flow (Mach=0.) to roughly 
C     Mach=3.
C
C     A hybrid method is used:
C       Laminar flow   - uses two equation momentum/energy method 
C        correlations for H23-H12 lifted from Eppler method
C       Turbulent flow - Green's lag entrainment formulation (three 
C        equations)
C
C       Transition is computed using an (e**n) method, integrating the
C        disturbance amplification rate. This adds another ODE to the 
C        two laminar equations, a nice symmetry with the three equation
C        turbulent integration.  Alternatively, Eppler's much simpler  
C        H23/RT correlation can be used to calculate transition (if 
C        AMPMAX is input as 0). 
C
C     The BL ODE's are integrated using a variable step size second-order
C     Runge-Kutta scheme that uses the change in H23 (here HSBXXX) in 
C     laminar flow and H12 (here HBXXX) in turbulent flow as a measure
C     of error.  This scheme seems much more robust than higher order 
C     R-K with error estimation and allows adequate checking and sub-
C     stepping for transition and separation.
C
C   Assumptions for method:
C     Adiabatic wall conditions 
C     Prandtl number of unity (no weird wall temperatures)
C     Isentropic flow
C     Viscosity is evaluated by Sutherland's law.
C
C  HHY 1/16/98
C  Modified to include 3D effects from velocity divergence.
C  [This is the input data UDIV in the call statement above]
C  If the divergence term is not required set the elements of this array
C  to 0.0 and the BL will be 2D (no divergence).
C  NOTE: These divergence terms only act as a source term on the RHS
C        of the momentum equation and play no other role in the shape
C        factor or skin friction correlations.
C
C  HHY 9/6/05 
C  Note: for axisymmetric cases the term added to the VonKarman mom. eqn
C        is -(theta/r)*(dr/dx)
C        This implies that the UDIV term for axisymmetric cases is defined
C        by UDIV = 1/r * dr/dx
C
C   Note: Currently dimensioned for a max of 500 input flow stations
C***********************************************************************
C
      REAL  MACH, MMIN, MINF, MINFSQ, M, M1, M2, MSQ, ME, 
     &      MUEMUI, NUENUI, NU, NU1, NU2, NE
C
      PARAMETER (NBLX = 500)
      DIMENSION  S(NBLX),      UE(NBLX),     UDIV(NBLX)
      DIMENSION  ME(NBLX),     QE(NBLX),     NE(NBLX)
      DIMENSION  DSTAR(NBLX),  THETA(NBLX),  RTHETA(NBLX), 
     &           H12(NBLX),    CFE(NBLX),    AMPF(NBLX)
C
C...Gamma, Turbulent recovery factor, min. Mach for compressible calcs
      DATA  GAM   / 1.4   /
      DATA  TRF   / .88   /
      DATA  MMIN  / .01   /
C
C...Step size in H23, Laminar sep in H23   <- used for laminar integration
      DATA  HSBSTP, HSBSEP / .005, 1.51509 /
C
C   Step size in H12, Turbulent sep in H12 <- used for turbulent integration 
      DATA  HBSTP,  HBSEP  / .005, 2.8     /
C
C...Error fraction, min Rtheta for turbulent flow
      DATA  ERRFRC  / 0.2  /
      DATA  RTTRB  / 40. /
C
C...Check for max number of BL points
      IF(NP.GT.NBLX) THEN
        WRITE(*,*) '*** Too many input points for BL ',NP
        STOP
      ENDIF
C
      GF    = 0.5 * (GAM-1.)
      GEXP1  = (GAM-1.) / GAM
      GEXP2  = 1. / (GAM-1.)
C
      MINF = MACH
      IF(MACH.LT.MMIN) MINF = MMIN
      MINFSQ = MINF*MINF
C
C
C...Compute stagnation temp. ratio - T0TINF
C
      T0TINF = 1.0 + GF*MINFSQ
C
C...Compute flow parameters at edge of boundary layer.
C    Parameters are referenced to free stream values.
C    UE is limited to prevent close-to-zero divide
C
      UMAX = 0.
      DO 50  I   = 1, NP
        UMAX = MAX(UMAX,UE(I))
        USQ  = UE(I)**2
	MSQ  = USQ / ( GF*(1.0-USQ) + 1.0/MINFSQ )
C
        TETINF = T0TINF / (1.+GF*MSQ)
        RHOERI = TETINF**GEXP2
        MUEMUI = SQRT(TETINF)*(1.505/(1.+.505/TETINF))
        NUENUI = MUEMUI / RHOERI 
C
        ME(I) = SQRT(MSQ)
        QE(I) = RHOERI*UE(I)**2
        NE(I) = NUENUI / UREN
   50 CONTINUE
C
      ITURB = 0
      ISEP  = 0
      SSEP  = S(NP)
      CFTOT = 0.
      AMP   = 0.
      STRAN = SFIX
      STRN  = SFIX
C
      TH = 0.
      RT = 0.
      H  = 2.236
      CF = 0.
      EPTRN = 1.
C
C
C
C***********************************************************************
C...Calculate the boundary layer at each station
C***********************************************************************
C 
      DO 700   I = 1, NP
C
        S2  = S(I)
        M2  = ME(I)
        U2  = UE(I)
        NU2 = NE(I)
        Q2  = QE(I)
        B2  = UDIV(I)
C
        IF (I .EQ. 1)  GO TO 600
C
        IF(U2.LT..001*UMAX) U2 = .001*UMAX
        DELS = S2 - S1
        IF (DELS.LE.0.)  GO TO 600
C
        SBEG = S1
        SEND = S2
        IF (ITURB.NE.0)  GO TO 300
C
        IF (SBEG.GE.STRAN .AND. RT.GT.RTTRB)  GO TO 300
C
C***********************************************************************
C...Laminar BL (Two equation - momentum and energy method)
C***********************************************************************
C
C...Initialize BL. First point specially handled, calculate initial 
C    momentum thickness from stagnation value, else use flat plate BL.
C
      IF (I.EQ.2) THEN
        DS = .5*DELS
        SBEG = SBEG + DS
        CALL FLWTRP (SBEG,S1,M1,U1,Q1,NU1,B1,S2,M2,U2,Q2,NU2,B2,
     &               M,U,DUDS,Q,NU,B)
C
        IF (U1.LE.0. .OR. (U2-U1)/U1.GT.1.)  THEN
C...Stagnation point
          TH = SQRT(.084*NU/DUDS)
          HSB = 1.61998
          THETA(1) = TH
          DSTAR(1) = TH*H
        ELSE
C...Flat plate 
          TH = .66411*SQRT(DS*NU/U)
          HSB = 1.57258
          H12(1) = 2.597*(1.+.113*M*M) + .290*M*M
        ENDIF
      ENDIF
C
C
C...Integrate the laminar BL from SBEG to SEND
C
      SINT = SBEG
C
  110 DS = SEND - SINT
      IF (DS.LE.0.)    GO TO 600
C
C...Set initial step size using change in shape factor
      CALL FLWTRP (SINT,S1,M1,U1,Q1,NU1,B1,S2,M2,U2,Q2,NU2,B2,
     &             M,U,DUDS,Q0,NU,B)
      CALL LAMBL (M,U,DUDS,NU,B,TH,HSB,HS,HB,H,RT,CF0,
     &            DTHDX0,DHSDX0,DADX0)
      DHSB0 = DHSDX0*DS
C
      DELHSB = HSBSTP
      IF (HSB.LT.HSBSEP+2.*HSBSTP)  DELHSB = 0.1*DELHSB
      IF (ABS(DHSB0).GT.DELHSB)  DS = DS*ABS(DELHSB/DHSB0)
C
  120 DTH0  = DTHDX0*DS
      DHSB0 = DHSDX0*DS
C
      SSM  = SINT + 0.5*DS
      THM  = TH   + 0.5*DTH0
      HSBM = HSB  + 0.5*DHSB0
      CALL FLWTRP (SSM,S1,M1,U1,Q1,NU1,B1,S2,M2,U2,Q2,NU2,B2,
     &             M,U,DUDS,Q,NU,B)
      CALL LAMBL (M,U,DUDS,NU,B,THM,HSBM,HS,HB,H,RT,CF,
     &            DTHDXM,DHSDXM,DADXM)
      DHSB1 = DHSDXM*DS
      DTH1  = DTHDXM*DS
C
      ERRHSB = DHSB1 - DHSB0
      ERRTOL = ERRFRC*DELHSB
      IF (ABS(ERRHSB).GT.ERRTOL)  THEN
        DS = 0.5*DS
        GO TO 120
      ENDIF
C
      SS1  = SINT + DS
      HSB1 = HSB + DHSB1 
      TH1  = TH  + DTH1
      CALL FLWTRP (SS1,S1,M1,U1,Q1,NU1,B1,S2,M2,U2,Q2,NU2,B2,
     &             M,U,DUDS,Q,NU,B)
      CALL LAMBL (M,U,DUDS,NU,B,TH1,HSB1,HS,HB,H,RT,CF,
     &            DTHDX1,DHSDX1,DADX1)
C
C...Check for transition
C     triggered by laminar separation, 
C     Amplification ratio, or specified transition
C
      AMP1 = AMP + 0.5*DS*(DADX0 + DADX1)
      EPTRN1 = (18.4*HSB1 - 21.74) - ALOG(RT)
C
      FRAC = 1.
      IF (RT.GT.RTTRB) THEN
       FRAC1 = 1.
       FRAC2 = 1.
       FRAC3 = 1.
       FRAC4 = 1.
       IF (HSB1.LE.HSBSEP)           FRAC1 = (HSB1-HSBSEP)/DHSB1
       IF (AMPMAX.LE.0.)  THEN
cc         IF (EPTRN1.LT.0.)         FRAC2 = EPTRN/(EPTRN-EPTRN1)
       ELSE
         IF (AMP1.GE.AMPMAX)         FRAC3 = (AMPMAX-AMP)/(AMP1-AMP)
       ENDIF
       IF (STRAN-SS1.LT.0. .AND. 
     &     RT.GT.RTTRB)                FRAC4 = (STRAN-SS1)/DS
C--- Limiter on overshoot of FRAC4 fixed trans point  HHY 1/21/98
       FRAC4 = MAX(FRAC4,0.0)  
       FRAC = MIN(FRAC,FRAC1,FRAC2,FRAC3,FRAC4)
      ENDIF
C
      FRAC = MAX(0.0,FRAC) 
      DS  = FRAC*DS
      TH  = TH  + FRAC*DTH1
      HSB = HSB + FRAC*DHSB1 
      AMP = AMP + FRAC*(AMP1-AMP)
      EPTRN = EPTRN1
      CFTOT = CFTOT + 0.5*DS*(CF*Q + CF0*Q0)
      SINT = SINT + DS
      STRN = SINT
C
      IF (FRAC.GE.1.)  GO TO 110
C
C...Did the laminar BL complete the step from S1 to S2 before transition?
      ITURB = -1
      IF (SINT.GE.S2)  GO TO 600
        SBEG = SINT
        SEND = S2
C***********************************************************************
C
C
C
C
C***********************************************************************
C...Turbulent B.L. - Green's lag entrainment formulation
C***********************************************************************
C
C...If this is the first turbulent BL step set initial values for 
C    turbulent BL before proceeding.
C
  300 IF (ITURB.LE.0)  THEN
        CALL FLWTRP (SBEG,S1,M1,U1,Q1,NU1,B1,S2,M2,U2,Q2,NU2,B2,
     &               M,U,DUDS,Q,NU,B)
        CALL TRBINI (M,U,DUDS,NU,B,TH,HB,CE,H,RT,CF)
        HBLIM = AMAX1(HB*1.1,HBSEP)
        ITURB = 1
        STRAN = SBEG
        IWAK = 0
        AMP  = 0.
      ENDIF
C
C
C...Integrate the turbulent BL from SBEG to SEND
C
      SINT = SBEG
C
C...Integrate the separated BL using an empirical relation
  410 IF (ISEP.NE.0)  THEN
        CALL FLWTRP (SINT,S1,M1,U1,Q1,NU1,B1,S2,M2,U2,Q2,NU2,B2,
     &               M,U,DUDS,Q,NU,B)
        HB = HBSEP
        HBLIM = HBSEP
        H  = (HB+1.)*(1.+0.2*TRF*M*M) - 1.
        UCHG = U/U2
        IF(UCHG.LT.1.) UCHG = 1.
        TH = TH * (UCHG)**(0.5*(5.+H-(M*M+MINFSQ)))
        CF = 0.
      ELSE
C
C...Integrate the attached BL
      DS = SEND - SINT
      IF (DS.LE.0.) GO TO 600
C
C...Set initial step size using change in shape factor or entrainment
      CALL FLWTRP (SINT,S1,M1,U1,Q1,NU1,B1,S2,M2,U2,Q2,NU2,B2,
     &             M,U,DUDS,Q0,NU,B)
      CALL TRBBL (IWAK,M,U,DUDS,NU,B,TH,HB,CE,
     &            H,RT,CF0,DTHDX0,DHBDX0,DCEDX0)
      DHB0 = DHBDX0*DS
      STP1 = (HB-1.3)/.5
      DELHB = STP1
      if (DELHB .LT. 1.) DELHB = 1.
      DELHB = DELHB*HBSTP
      IF (ABS(DHB0).GT.DELHB)  DS = DS*ABS(DELHB/DHB0)
C
  420 DTH0 = DTHDX0*DS
      DHB0 = DHBDX0*DS
      DCE0 = DCEDX0*DS
C
      SSM = SINT + 0.5*DS
      THM = TH   + 0.5*DTH0
      HBM = HB   + 0.5*DHB0
      CEM = CE   + 0.5*DCE0
      CALL FLWTRP (SSM,S1,M1,U1,Q1,NU1,B1,S2,M2,U2,Q2,NU2,B2,
     &             M,U,DUDS,Q,NU,B)
      CALL TRBBL (IWAK,M,U,DUDS,NU,B,THM,HBM,CEM,
     &            H,RT,CF,DTHDXM,DHBDXM,DCEDXM)
      DTH1 = DTHDXM*DS
      DHB1 = DHBDXM*DS
      DCE1 = DCEDXM*DS
C
      ERRHB = DHB1 - DHB0
      ERRTOL = ERRFRC*DELHB
      IF (ABS(ERRHB).GT.ERRTOL)  THEN
        DS = 0.5*DS
        GO TO 420
      ENDIF
C
      SS1 = SINT + DS
      TH1 = TH + DTH1
      HB1 = HB + DHB1
      CE1 = CE + DCE1
      CALL FLWTRP (SS1,S1,M1,U1,Q1,NU1,B1,S2,M2,U2,Q2,NU2,B2,
     &             M,U,DUDS,Q,NU,B)
      CALL TRBBL (IWAK,M,U,DUDS,NU,B,TH1,HB1,CE1,
     &            H,RT,CF,DTHDX1,DHBDX1,DCEDX1)
C
      FRAC = 1.
      IF (HB1.GE.HBLIM) THEN
        ISEP = 1
        FRAC = (HBLIM-HB)/DHB1
        DS = FRAC*DS
        SSEP = SINT + DS
      ENDIF
      TH = TH + FRAC*DTH1
      HB = HB + FRAC*DHB1
      CE = CE + FRAC*DCE1
      CFTOT = CFTOT + 0.5*DS*(CF*Q + CF0*Q)
      SINT = SINT + DS
      GO TO 410
C
      ENDIF
C
C
C*********************************************************************
C...Find displacement thickness and skin friction
C*********************************************************************
C
  600 THETA (I) = TH
      RTHETA(I) = U2*TH/NU2
      H12   (I) = H
      DSTAR (I) = TH*H
      CFE   (I) = CF
      AMPF  (I) = AMP
cc      write(*,*) i,s2,amp
C
      S1  = S2
      M1  = M2
      U1  = U2
      NU1 = NU2
      Q1  = Q2
      B1  = B2
C
  700 CONTINUE
C
C
C*********************************************************************
C...Profile drag is computed by applying a compressible version
C    of the Squire/Young formula using values of BL parameters
C    from the trailing edge extrapolated to unity H12 in wake (from Eppler).
C
      H12TE = AMIN1(2.5,H)
      XPOTE = 0.5 *(H12TE + 5.0 - (M2**2+MINFSQ))
      CDP   = 2.*TH * U2**XPOTE
C
C...Normalize the profile drag (CDP) and total skin friction (CFTOT) 
C    based on total S (chord length?) - note that this is only 
C    approximately the chord for a curved surface (exact for flat plate)
C
      SCHORD = S(NP) - S(1)
      CDP   = CDP   / SCHORD
      CFTOT = CFTOT / SCHORD
C
      RETURN
      END
C
      SUBROUTINE FLWTRP (X,
     &                   X1,M1,UE1,QUE1,NUE1,B1,
     &                   X2,M2,UE2,QUE2,NUE2,B2,
     &                   M, UE, DUEDX, QUE, NUE, B) 
C
C...Interpolates flow conditions between input stations
C
C...Input    X                             !independent variable
C            X1, M1, UE1, QUE1, NUE1, B1   !external flow data at 1 (beg)
C            X2, M2, UE2, QUE2, NUE2, B2   !external flow data at 2 (end)
C
C...Output   M, UE, DUEDX, QUE, NUE, B     !external flow data at X
C
      REAL  M1, M2, M, NUE, NUE1, NUE2
C
C...Linear interpolation
      DX   = X -X1
      DX21 = X2-X1
      IF(DX21.LE.0.0) THEN
       WRITE(*,*) 'FLWTRP DX21=0, X,X1,X2 ',X,X1,X2
      ENDIF
      M   = M1   + (M2-M1)    *DX/DX21
      UE  = UE1  + (UE2-UE1)  *DX/DX21
      QUE = QUE1 + (QUE2-QUE1)*DX/DX21
      NUE = NUE1 + (NUE2-NUE1)*DX/DX21
      B   = B1   + (B2-B1)    *DX/DX21
      DUEDX = (UE2-UE1)/DX21
C
      RETURN
      END
C
      SUBROUTINE LAMBL (M,UE,DUEDX,NUE,UDIV,
     &                  THETA,HSK,
     &                  HS,HK,H,RT,CF,DTHDX,DHSKDX,DADX)
C
C   Copyright  May 1993    Harold Youngren - all rights reserved
C
C...Input    M, UE, NUE, DUEDX         !external flow data
C            UDIV                      !velocity divergence
C            THETA, HSK                !BL state variable
C...Output   HS,HK,H,RT,CF             !BL parameters
C            DTHDX,DHSKDX,DADX         !derivatives of variables
C
C...Laminar boundary layer method
C
C...Calculate quantities for derivative evaluation
C
C This system has two first order equations...
C     independent variable is X
C       dependent variables are   THETA    HSK
C              with derivatives  (DTHDX)  (DHSKDX)  
C
C The amplification rate derivative (DADX) is calculated for use
C     as a transition estimate
C
C--- HHY 1/20/98 added velocity divergence source term GRADB to mom. eqn
C
      REAL  M, MSQ, NUE
C
      MSQ = M*M
C
      GRAD  = THETA*DUEDX/UE
      GRADB = THETA*UDIV
      RT    = UE*THETA / NUE
C
C---- Compressible correction for energy thickness ( Whitfield )
      HS = (HSK+0.028*MSQ) / (1.+0.014*MSQ) 
      DHSKHS = (1.+0.014*MSQ)
C
C---- Laminar HS/HK correlations  ( from Eppler )
      IF (HSK.LT.1.51509) THEN
       HK = 4.02922
       SF = 0.
      ELSE
       IF (HSK.LT.1.57258) THEN
        HK = 4.02922 - SQRT(HSK-1.51509)*
     &       (583.60182 - 724.55916*HSK + 227.1822*HSK**2)
        SF = 2.512589 - 1.686095*HK + .391541*HK**2 - .031729*HK**3
       ELSE
        HK = 79.870845 - 89.58214*HSK + 25.715786*HSK**2
        SF = 1.372391 -   4.22625*HSK +  2.221687*HSK**2
       ENDIF
      ENDIF
C
C---- Shape factor (compressible - Whitfield) 
      H  = HK*(1.+.113*MSQ) + .290*MSQ
C
C---- Laminar skin friction function  ( Cf )
      CF = 2.*SF / RT
C
C---- density thickness shape parameter  ( H** )
      HSS = M * (0.064/(HK-0.8) + 0.251)
C
C---- Laminar dissipation function (CD)
      DS = 7.853976 - 10.26055*HSK + 3.418898*HSK**2
      CD = DS / RT
C
C...Amplification rate equation (from Drela)
C
      HMI = 1. / (HK-1.)
      RTLOG = 0.
      IF (RT.GT.0.) RTLOG = ALOG10(RT)
      GRCRIT = (1.415*HMI-0.489)*TANH(20.*HMI-12.9) + 3.295*HMI + 0.44
      DADX = 0.
      IF (RTLOG.GT.GRCRIT)  THEN
        T1 = 2.4*HK - 3.7 + 2.5*TANH(1.5*(HK-3.1))
        DADRT = .01*SQRT(T1**2 + .25)
        TFSQ  = (6.54*HK-14.07)/(HK**2)
        BUH   = (0.058*(4.-HK)**2/(HK-1.) - 0.068) / TFSQ
        DADX  = DADRT/THETA * TFSQ*0.5*(BUH+1.)
      ENDIF
C
C
C...Momentum equation
      DTHDX = 0.5*CF - (H+2.-MSQ)*GRAD - GRADB
C
C...Shape factor equation
      DHSDX = (2.*CD - 0.5*CF*HS - (2.*HSS/HS + 1.-H)*HS*GRAD)/THETA
      DHSKDX = DHSKHS*DHSDX
C
      RETURN
      END
C
      SUBROUTINE TRBINI (M,UE,DUEDX,NUE,UDIV,
     &                   THETA,HB,CE,H,RT,CF)
C
C   Copyright  May 1993    Harold Youngren - all rights reserved
C
C...Input    M, UE, DUEDX, NUE     !external flow data 
C            UDIV                  !velocity divergence
C            THETA                 !initial momentum thickness
C
C...Output   HB, CE                !initial state variables
C            H, RT, CF             !boundary layer parameters
C
C...Calculate initial quantities using equilibrium assumptions
C     a Newton-Raphson iteration is used to determine equilibrium
C     shape factor (HB) and entrainment (CE) corresponding to the 
C     initial conditions (THETA and DUEDX)
C
C     Uses Green's lag entrainment turbulent boundary layer method
C     follows R&M 3791
C
C
      REAL  M, MSQ, NUE
C
      DATA  TRF   / .88  /
      DATA  TOL   / .001 /
      DATA  ITER  / 30   /
      DATA  CEMIN / -.009/
      DATA  RTMIN / 40.  /
C
C
      MSQ = M*M
C
      GRAD =  THETA*DUEDX/UE
      GRADB = THETA*UDIV
      RT = UE*THETA / NUE
ccc   WRITE(*,*) 'RT,UE,THETA @ transition ',RT,UE,THETA
      IF(RT.LT.RTMIN) RT = RTMIN
C
      FB = 1 + 0.2*MSQ
      FC = SQRT(FB)
      FD = 1 + 0.04*MSQ
      FR = 1 + 0.056*MSQ
C
      CF0 = (0.01013/(ALOG10(FR*RT) - 1.02) - 0.00075) / FC
      HB0 = 1. / (1. - 6.55*SQRT(0.5*CF0*FD))
C
C...Initially assume HB(incompressible equiv. of H)=1.4
C
      HB = 1.4
C
C...Start iterative process to find HB and CEEQ
C
      DO 100 I = 1, ITER
C
      CF = CF0 * (0.9/(HB/HB0-0.4) -0.5)
      DCFDHB = -0.9*CF0 / (HB0*(HB/HB0-0.4)**2)
C
      H  = (HB+1.)*(1.+0.2*TRF*MSQ) - 1.
      DHDHB = (1.+0.2*TRF*MSQ)
      H1 = 3.15 + 1.72/(HB-1.) - 0.01*(HB-1.)**2
C
      GRDEQ0 = 1.25/H * (0.5*CF - ((HB-1.)/(6.432*HB))**2/FD)
      DGRDHB = .625/H * DCFDHB 
     &         - 1.25/H**2 * (0.5*CF-((HB-1.)/(6.432*HB))**2/FD)*DHDHB
     &         + 2.5/(H*FD) * ( ((HB-1.)/(6.432*HB))**2/HB 
     &                         - (HB-1.)/(6.432*HB)**2 )
      CEEQ0  = H1*(0.5*CF - (H+1)*GRDEQ0)
      IF (CEEQ0 .LT. CEMIN)  CEEQ0 = CEMIN
C
      REZ = GRAD - GRDEQ0
      DQDHB = - DGRDHB
C
      DHB = -REZ/DQDHB
C
      DRATIO = DHB/HB
      RLX = 1.
      IF (DRATIO .LT. -.5)  RLX = -.5/DRATIO
      IF (DRATIO .GT. 1.0)  RLX = 1.0/DRATIO
      HB = HB + RLX*DHB
C
      IF (ABS(DHB) .LT. TOL)  GO TO 200
C
  100 CONTINUE
C
C--Iteration did not converge
      WRITE (6,110)
  110 FORMAT (//' Iteration for HB, CE at transition failed')
C      STOP
C--Iteration did converge
  200 CE = CEEQ0
      RETURN
      END
C
      SUBROUTINE TRBBL (IWAKE,
     &                  M,UE,DUEDX,NUE,UDIV,
     &                  THETA,HB,CE,H,RT,CF,
     &                  DTHDX,DHBDX,DCEDX)
C
C   Copyright  May 1993    Harold Youngren - all rights reserved
C
C...Input    M, UE, DUEDX, NUE, UDIV  !external flow data at 1 (beg)
C            UDIV                     !velocity divergence at 1
C            THETA, HB, CE            !BL state variable
C            IWAKE                    !wake flag (set to 0 for normal BL,
C                                     !           set to 1 for wake BL)
C...Output   H, RT, CF                !BL parameters
C            DTHDX,DHBDX,DCEDX     !derivatives of variables
C
C...Green's lag entrainment turbulent boundary layer method
C     follows R&M 3791
C
C--- HHY 1/20/98 added velocity divergence source term GRADB to mom. eqn
C
C...Calculate quantities for derivative evaluation
C
C This system has three first order equations...
C     independent variable is X
C       dependent variables are    THETA    HB       CE
C              with derivatives  (DTHDX)  (DHBDX)  (DCEDX)
C
      REAL  M, MSQ, NUE
C
      DATA  TRF   / .88  /
      DATA  CEMIN / -.009/
      DATA  RTMIN / 40.  /
C
C...Special treatment if on a wake 
C
      if(UE.LT.0.0) then
        write(*,*) 'UE neg TRBBL ',UE
ccc        pause
      endif
      COEF = 1.
      IF (IWAKE .NE. 0)  COEF = 0.5
C
C...Limit CE to CEMIN to avoid excessive negative entrainment
C
      IF (CE .LT. CEMIN)  CE = CEMIN
C
      MSQ = M*M
C
      GRAD  = THETA*DUEDX/UE
      GRADB = THETA*UDIV
      RT = UE*THETA / NUE
      RT = AMAX1(RT,RTMIN)
C
      FA = 1 + 0.1*MSQ
      FB = 1 + 0.2*MSQ
      FC = SQRT(FB)
      FD = 1 + 0.04*MSQ
      FR = 1 + 0.056*MSQ
C
C
      CF0 = (0.01013/(ALOG10(FR*RT) - 1.02) - 0.00075) / FC
C
      IF (IWAKE .NE. 0)  CF0 = 0.
C
      HB0 = 1. / (1. - 6.55*SQRT(0.5*CF0*FD))
      CF = CF0 * (0.9/(HB/HB0-0.4) -0.5)
C
      H  = (HB+1.)*(1.+0.2*TRF*MSQ) - 1.
      H1 = 3.15 + 1.72/(HB-1.) - 0.01*(HB-1.)**2
      DHBDH1 = - (HB-1.)**2 / (1.72 + 0.02*(HB-1.)**3)
C
      CT = FA*(0.024*CE + 1.2*CE**2 + 0.32*CF0)
      F  = (0.02*CE + CE**2 + 0.8*CF0/3.) / (0.01 + CE)
C
      GRDEQ0 = 1.25/H * (0.5*CF - ((HB-1.)/(6.432*HB))**2/FD)
C
      CEEQ0  = H1*(0.5*CF - (H+1)*GRDEQ0)
      IF (CEEQ0 .LT. CEMIN)  CEEQ0 = CEMIN
C
      CTEQ0  = (0.024*CEEQ0 + 1.2*CEEQ0**2 + 0.32*CF0) * FA
      C = CTEQ0 / (FA*COEF**2) - 0.32*CF0
C
      CEEQ = SQRT(C/1.2 + .0001) - 0.01
      IF (CEEQ .LT. CEMIN)  CEEQ = CEMIN
      GRDEQ = (0.5*CF - CEEQ/H1) / (H + 1.)
C
C
C...Momentum equation
      DTHDX = 0.5*CF - (H+2.-MSQ)*GRAD - GRADB
C
C...Shape factor equation
      DHBDX = DHBDH1*(CE - H1*(0.5*CF - (H+1)*GRAD)) / THETA
C
C...Entrainment equation
      DCEDX = F/THETA * ( 2.8/(H+H1) * (SQRT(CTEQ0)-COEF*SQRT(CT)) 
     &                   + GRDEQ 
     &                   - GRAD*(1.+.075*MSQ*FB/FA) )
C
      RETURN
      END
