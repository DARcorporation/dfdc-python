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
C--- Aero data stored for one or more radial aerodynamic sections
C   
C-- aero data quantities for each defined radial aerodynamic section
C  NAERO     Number of aerodynamic datasets defined (NAERO>=1)
C  XIAERO    Radial station r/R where aero dataset is defined
C  AERODATA  Aerodynamic definition of the blade section at XIAERO
C            AERODATA( 1,x) = A0 (angle of zero lift)
C            AERODATA( 2,x) = CLMAX (Max CL)
C            AERODATA( 3,x) = CLMIN (Min CL)
C            AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
C            AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
C            AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
C            AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
C            AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
C            AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
C            AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
C            AERODATA(11,x) = REREF (reference Reynold's number)
C            AERODATA(12,x) = REXP (Reynold's number exponent Cd~Re^REXP)
C            AERODATA(13,x) = MCRIT (critical Mach #)
C            AERODATA(14,x) = TOC (thickness/chord)
C            AERODATA(15,x) = DCDCL2S (Scndary, annulus drag param d(Cd)/dCL^2)
C=========================================================================
C
C     Version 070-ES1
C     Philip Carter, Esotec Developments, February 2009
C     philip (at) esotec (dot) org
C
C     Changes from 0.70:
C
C     CL/CD plotting fixed (missing argument).
C     READ and WRIT commands fixed (LU).
C     FLIP command to toggle to/from neg-BGam parameters.
C     DISP formats repaired. Message when set for neg-BGam.
C     Interpolated A0 (AZERO) stored by GETCLCDCM.
C     Multi-Re plotting with constant Mach or vice-versa (AEROPLT2).
C     Modified plotting interface.
C     All plotting functionality duplicated in EDIT.
C     PLOTMACH and PLOTREYN subroutines to avoid code duplication.
C     HARD and ANNO fixed.
C     Disk index displayed with data for multi-disk cases.
C     Disk index added to plot titles for multi-disk cases.
C     Mach constant or Re constant written to plot legends.
C     Various fixes to control structure and cosmetics.
C
C=========================================================================
C
C
      SUBROUTINE DFAERO
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      CHARACTER*4 COMAND,COMOLD
      CHARACTER*132 COMARG,LINE, FNAME, SAVFIL,ARGOLD
      CHARACTER*1 CHKEY, ANS, CIND
      CHARACTER*16 PROMPT
      CHARACTER*8 PROMPT2
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR, LANS, RFIX
      PARAMETER (NMACH = 10)
      PARAMETER (NREYN = 10)
      DIMENSION AMACH(NMACH), AREYN(NREYN)
      REAL MCONST
C
      COMMON /AFPLOTS/ AMACH,AREYN,NM,NRE,MCONST,RCONST,RFIX
C
C--------------------------------------------------------------
C     Display, change, or write out airfoil section properties
C--------------------------------------------------------------
C---- check for loaded file
C
      IF(.NOT.LLOAD) THEN
         WRITE(*,*)
         WRITE(*,*)'No case file loaded'
         RETURN
      ENDIF
C
C---- start with disk 1
      NDSK = 1
      NR = NDSK
C
      IF(NPTOT.EQ.0) CALL GENGEOM
C
      LUSAVE=2
      LCHANGE = .FALSE.
      RFIX = .FALSE.
      NRE=0
      NM =0
      MCONST = 0.0
      RCONST = AERODATA(11,1,NR)
C
C---- "last" command (nothing yet)
      COMAND = '****'
      COMARG = ' '
      FNAME  = ' '
      COMOLD = COMAND
      ARGOLD = COMARG
C
C-------------------------------------------------------------
C---- begin user-interaction loop
C............................................................
C
 1    FORMAT(A)
C
 900  IF(NROTOR.GT.1)THEN
         WRITE(*,925) NR
      ELSE
         WRITE(*,920)
      ENDIF
C
      WRITE(*,927)
      DO N = 1, NAERO(NR)
         CALL GETAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &                DCLDA,DCLDA_STALL,DCL_STALL,
     &                CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &                CMCON,MCRIT,TOC,REREF,REXP)
        WRITE(*,930) N,XIAERO(N,NR),CLMAX,CLMIN,CDMIN,MCRIT,REXP,REREF
      END DO
C
 920  FORMAT(/' Defined aerodynamic sections: ')
 925  FORMAT(/' Defined aerodynamic sections: Disk',I2)
 927  FORMAT( '  N    r/R   CLmax   CLmin    CDmin',
     &        '   Mcrit   REexp       REref')
C
 930  FORMAT(I3,3(1X,F7.4),1X,F8.5,1X,F7.4,1X,F7.4,1X,G12.4)
C
C   N      r/R   CLmax   CLmin    CDmin   Mcrit   REexp        REref
Cxiiixffffffffxfffffffxfffffffxffffffffxfffffffxfffffffxffffffffffff
C
C
 950  CONTINUE
      PROMPT = '.AERO^'
C
      IF(NDSK.GT.0 .AND. NROTOR.GT.1) THEN
C------ add current target disk to prompt
        IF    (NDSK.LT. 10) THEN
         PROMPT(6:7) = '(' // PNUM(NDSK)(2:2)
         K = 8
        ELSEIF(NDSK.LT.100) THEN
         PROMPT(6:8) = '(' // PNUM(NDSK   )
         K = 9
        ENDIF
        PROMPT(K:K+1) = ')^'
      ENDIF
      CALL ASKC(PROMPT,COMAND,COMARG)
C
C---- process previous command ?
C
      IF(COMAND(1:1).EQ.'!') THEN
        IF(COMOLD.EQ.'****') THEN
          WRITE(*,*) 'Previous .AERO command not valid'
          GO TO 900
        ELSE
          COMAND = COMOLD
          COMARG = ARGOLD
        ENDIF
      ENDIF
C
      IF(COMAND.EQ.'    ') THEN
        IF(LPLOT) THEN
          CALL PLEND
          CALL CLRZOOM
        ENDIF
        IF(LCHANGE) CALL SETIAERO
        RETURN
      ENDIF
C
C---- can jump in here if COMARG string needs processing
 960  CONTINUE
C
C---- was only a new target-element index typed?
      CIND = COMAND(1:1)
      IF(INDEX('0123456789',CIND) .NE. 0) THEN
        READ(COMAND,*,ERR=970) ID
        IF(ID.LT.1 .OR. ID.GT.NROTOR) GO TO 970
C
C        IF(IRTYPE(ID).EQ.1) THEN
C           WRITE(*,1020) ID
C           GO TO 900
C        ELSEIF(IRTYPE(ID).EQ.2.AND.OMEGA(ID).NE.0.) THEN
C           WRITE(*,1022) ID
C        ELSEIF(IRTYPE(ID).EQ.2.AND.OMEGA(ID).EQ.0.) THEN
C           WRITE(*,1024) ID
C        ELSE
C           WRITE(*,1026) ID
C           GO TO 900
C        ENDIF
C
        NDSK = ID
        NR = NDSK
        GO TO 900
      ENDIF
C
 970  CONTINUE
C
C---- extract numbers (if any) from command argument string
      DO I=1, 20
        IINPUT(I) = 0
        RINPUT(I) = 0.0
      ENDDO
      NINPUT = 20
      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
      NINPUT = 20
      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C
 1000 FORMAT(A)
 1010 FORMAT(1X,A4,' command not recognized.  Type "?" for list')
 1020 FORMAT(/,1X, 'Disk',I2,' is an actuator disk')
 1022 FORMAT(/,1X, 'Disk',I2,' is a bladed rotor')
 1024 FORMAT(/,1X, 'Disk',I2,' is a bladed stator')
 1026 FORMAT(/,1X, 'Disk',I2,' is undefined')
C
C---- can jump in here if COMAND string has been set elsewhere
 980  CONTINUE
C
C------------------------------------------------------------------
C
      IF(COMAND.EQ.'HELP'.OR. COMAND.EQ.'?   ') THEN
        IF(NROTOR.GT.1)THEN
           WRITE(*,1100)
           WRITE(*,1110)
        ELSE
           WRITE(*,*)
           WRITE(*,1110)
        ENDIF
      GOTO 950
      ENDIF
C
      IF(COMAND.EQ.'DISP'.OR.COMAND.EQ.'D   ') GO TO 10
      IF(COMAND.EQ.'NEW ') GO TO 12
      IF(COMAND.EQ.'DEL ') GO TO 13
      IF(COMAND.EQ.'EDIT') GO TO 15
      IF(COMAND.EQ.'ROLD') GO TO 20
      IF(COMAND.EQ.'WRIT') GO TO 240
      IF(COMAND.EQ.'READ') GO TO 220
C
      IF(COMAND.EQ.'PLOM') GO TO 60
      IF(COMAND.EQ.'PLOR') GO TO 60
      IF(COMAND.EQ.'SETM') GO TO 80
      IF(COMAND.EQ.'SETR') GO TO 90
      IF(COMAND.EQ.'ANNO') GO TO 100
      IF(COMAND.EQ.'HARD') GO TO 110
      WRITE(*,1050) COMAND
      GO TO 900
C
 1050 FORMAT(1X,A4,' command not recognized.' //
     &             '  Type "?" for list, <Return> to exit menu.')
 1100 FORMAT(
     &       /'   <i>     Select disk index')
 1110 FORMAT(/'   DISP<d> Display airfoil parameters'
     &       /'   NEW     Create a new aero section'
     &       /'   DEL  i  Delete an aero section'
     &       /'   EDIT i  Edit section aero data'
     &       /'   READ    Read airfoil parameters from disk file'
     &       /'   ROLD    Read airfoil file (old dfdc format)'
     &       /'   WRIT    Write airfoil parameters to disk file'
C
     &      //'   PLOM i  Plot airfoil polars at multiple Mach'
     &       /'   PLOR i  Plot airfoil polars at multiple Re'
     &       /'   SETM r  Set constant Mach for Re plots'
     &       /'   SETR r  Set constant Re for Mach plots'
     &       /'   ANNO    Annotate current plot'
     &       /'   HARD    Hardcopy current plot')
C
C
C-------------------------------------------------------------
C--- Display aero data
 10   DO N = 1, NAERO(NR)
       CALL GETAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &              CMCON,MCRIT,TOC,REREF,REXP)
C
       IF(ABS(CLMAX).LT.ABS(CLMIN).AND.CLDMIN.LE.0.0) THEN
         WRITE(*,1425) N,XISECT
       ELSE
         WRITE(*,1420) N,XISECT
       ENDIF
C
       A0DEG = A0/DTR
       WRITE(*,1500) A0DEG,CDMIN
       WRITE(*,1510) DCLDA,CLDMIN
       WRITE(*,1515) DCLDA_STALL,DCDCL2
       WRITE(*,1520) CLMAX,REREF
       WRITE(*,1530) CLMIN,REXP
       WRITE(*,1540) DCL_STALL,CMCON
       WRITE(*,1550) MCRIT
C
C---- old format with DCDCL2S and TOC
c       A0DEG = A0/DTR
c       WRITE(*,1500) A0DEG,DCDCL2
c       WRITE(*,1510) DCLDA,DCDCL2S
c       WRITE(*,1515) DCLDA_STALL,REREF
c       WRITE(*,1520) CLMAX,REXP
c       WRITE(*,1530) CLMIN,CMCON
c       WRITE(*,1540) DCL_STALL,MCRIT
c       WRITE(*,1550) CDMIN,TOC
c       WRITE(*,1560) CLDMIN
C
      END DO
      GO TO 900
C
C-------------------------------------------------------------
C--- Define new aero section
 12   IF(NAERO(NR).GT.1) THEN
        WRITE(*,1200)
        READ(*,1) LINE
        IF(LINE.EQ.' ') GO TO 900
        READ(LINE,*,ERR=12) N
        IF(N.LT.0 .OR. N.GT.NAERO(NR)) GO TO 900
       ELSE
        N = 1
      ENDIF
      CALL GETAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &             DCLDA,DCLDA_STALL,DCL_STALL,
     &             CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &             CMCON,MCRIT,TOC,REREF,REXP)
 121  WRITE(*,1220) N,XISECT
      READ(*,1) LINE
      READ(LINE,*,ERR=121) XISECT
      IF(XISECT.LT.0.0 .OR. XISECT.GT.1.0) THEN
        WRITE(*,*) '*** Section r/R must be in range 0->1'
        GO TO 121
      ENDIF
      CALL PUTAERO(NR,NAERO(NR)+1,XISECT,A0,CLMAX,CLMIN,
     &             DCLDA,DCLDA_STALL,DCL_STALL,
     &             CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &             CMCON,MCRIT,TOC,REREF,REXP)
C
      NAERO(NR) = NAERO(NR)+1
C--- Sort aero data sections by XI location
      CALL SORTAR(NAERO(NR),XIAERO(1,NR),AERODATA(1,1,NR),NDX)
      LCHANGE = .TRUE.
      GO TO 900
C
 1200 FORMAT(/' Enter index of aero section to copy: ')
 1220 FORMAT(/' Section # ',I3,' @ r/R = ',F10.4,
     &       /' Enter r/R for new section: ')
C
C-------------------------------------------------------------
C--- Delete aero section
 13   IF(NAERO(NR).GT.1) THEN
        IF(NINPUT.GE.1) THEN
          NDEL = IINPUT(1)
        ELSE
          WRITE(*,1300)
          READ(*,1) LINE
          IF(LINE.EQ.' ') GO TO 900
          READ(LINE,*,ERR=13) NDEL
        ENDIF
        IF(NDEL.LT.0 .OR. NDEL.GT.NAERO(NR)) GO TO 900
       ELSE
        WRITE(*,*) '*** Cannot delete all aero sections'
        GO TO 900
      ENDIF
      WRITE(*,1320) NDEL,XIAERO(NDEL,NR)
      CALL ASKL('Confirm delete?^',LANS)
      IF(.NOT.LANS) GO TO 900
      DO N = NDEL+1, NAERO(NR)
        CALL GETAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &               DCLDA,DCLDA_STALL,DCL_STALL,
     &               CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &               CMCON,MCRIT,TOC,REREF,REXP)
        CALL PUTAERO(NR,N-1,XISECT,A0,CLMAX,CLMIN,
     &               DCLDA,DCLDA_STALL,DCL_STALL,
     &               CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &               CMCON,MCRIT,TOC,REREF,REXP)
      END DO
      NAERO(NR) = NAERO(NR)-1
      LCHANGE = .TRUE.
      GO TO 900
C
 1300 FORMAT(/' Enter index of aero section to delete: ')
 1320 FORMAT(/' Section # ',I3,' @ r/R = ',F10.4)
C
C-------------------------------------------------------------
C--- Edit aero data
C-------------------------------------------------------------
 15   IF(NAERO(NR).GT.1) THEN
        IF(NINPUT.GE.1) THEN
          N = IINPUT(1)
        ELSE
          WRITE(*,1400)
          READ(*,1) LINE
          IF(LINE.EQ.' ') GO TO 900
          READ(LINE,*,ERR=15) N
        ENDIF
C
        IF(N.LT.0 .OR. N.GT.NAERO(NR)) GO TO 900
      ELSE
        N = 1
      ENDIF
C
      CALL AEROEDIT(NR,N)
      GO TO 900
C
C-------------------------------------------------------------
C--- Read aero data from file - old format with DCDDCL2 and TOC
C
 20   CALL ASKS(' Enter aero data filename^',FNAME)
      IF(FNAME.EQ.' ') GO TO 900
      OPEN(LUSAVE,FILE=FNAME,STATUS='OLD',ERR=900)
      DO N = 1, NAX
C
       READ(LUSAVE,1,END=26,ERR=26) LINE
       READ(LUSAVE,1,END=26,ERR=26) LINE
       READ(LINE(7:10),*) NN
       READ(LINE(18:26),*) XISECT
       READ(LUSAVE,1) LINE
       READ(LUSAVE,1) LINE
       READ(LINE(27:35),*) A0DEG
       READ(LINE(66:75),*) DCDCL2
       READ(LUSAVE,1) LINE
       READ(LINE(27:35),*) DCLDA
       READ(LINE(66:75),*) DCDCL2S
       READ(LUSAVE,1) LINE
       READ(LINE(27:35),*) DCLDA_STALL
       READ(LINE(66:75),*) REREF
       READ(LUSAVE,1) LINE
       READ(LINE(27:35),*) CLMAX
       READ(LINE(66:75),*) REXP
       READ(LUSAVE,1) LINE
       READ(LINE(27:35),*) CLMIN
       READ(LINE(66:75),*) CMCON
       READ(LUSAVE,1) LINE
       READ(LINE(27:35),*) DCL_STALL
       READ(LINE(66:75),*) MCRIT
       READ(LUSAVE,1) LINE
       READ(LINE(27:35),*) CDMIN
       READ(LINE(66:75),*) TOC
C       READ(LINE(66:75),*) 
       READ(LUSAVE,1) LINE
       READ(LINE(27:35),*) CLDMIN
       READ(LUSAVE,1,END=26,ERR=26) LINE
C
       A0 = A0DEG*DTR
       CALL PUTAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &              CMCON,MCRIT,TOC,REREF,REXP)
      END DO
 26   CLOSE(LUSAVE)
      NAERO(NR) = N-1
      CALL SORTAR(NAERO(NR),XIAERO(1,NR),AERODATA(1,1,NR),NDX)
      LCHANGE = .TRUE.
      GO TO 900
C
C-------------------------------------------------------------
C--- Write aero data to file - old format
c 24   SAVFIL = ' '
c      CALL OPFILE(LUSAVE,SAVFIL)
c      DO N = 1, NAERO(NR)
c       CALL GETAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
c     &              DCLDA,DCLDA_STALL,DCL_STALL,
c     &              CDMIN,CLDMIN,DCDCL2,DCDCL2S,
c     &              CMCON,MCRIT,TOC,REREF,REXP)
c       A0DEG = A0/DTR
c       WRITE(LUSAVE,1420) N,XISECT
c       WRITE(LUSAVE,1500) A0DEG,DCDCL2
c       WRITE(LUSAVE,1510) DCLDA,DCDCL2S
c       WRITE(LUSAVE,1515) DCLDA_STALL,REREF
c       WRITE(LUSAVE,1520) CLMAX,REXP
c       WRITE(LUSAVE,1530) CLMIN,CMCON
c       WRITE(LUSAVE,1540) DCL_STALL,MCRIT
c       WRITE(LUSAVE,1550) CDMIN,TOC 
c       WRITE(LUSAVE,1560) CLDMIN
c      END DO
c      CLOSE(LUSAVE)
c      GO TO 900
C
C------------------------------------------------------------
C----Read aero data from file - xrotor format
C
 220  CALL ASKS(' Enter aero data filename^',FNAME)
      IF(FNAME.EQ.' ') GO TO 900
      OPEN(LUSAVE,FILE=FNAME,STATUS='OLD',ERR=900)
      DO N = 1, NAX
C
       READ(LUSAVE,1,END=226,ERR=226) LINE
       READ(LUSAVE,1,END=226,ERR=226) LINE
       READ(LINE(9:12),*) NN
       READ(LINE(20:30),*) XISECT
       READ(LUSAVE,1) LINE
       READ(LUSAVE,1) LINE
       READ(LINE(24:39),*) A0DEG
       READ(LINE(62:77),*) CDMIN
       READ(LUSAVE,1) LINE
       READ(LINE(24:39),*) DCLDA
       READ(LINE(62:77),*) CLDMIN
       READ(LUSAVE,1) LINE
       READ(LINE(24:39),*) DCLDA_STALL
       READ(LINE(62:77),*) DCDCL2
       READ(LUSAVE,1) LINE
       READ(LINE(24:39),*) CLMAX
       READ(LINE(62:77),*) REREF
       READ(LUSAVE,1) LINE
       READ(LINE(24:39),*) CLMIN
       READ(LINE(62:77),*) REXP
       READ(LUSAVE,1) LINE
       READ(LINE(24:39),*) DCL_STALL
       READ(LINE(62:77),*) CMCON
       READ(LUSAVE,1) LINE
       READ(LINE(62:77),*) MCRIT
       READ(LUSAVE,1,END=226,ERR=226) LINE
C
       A0 = A0DEG*DTR
       CALL PUTAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &              CMCON,MCRIT,TOC,REREF,REXP)
      END DO
 226  CLOSE(LUSAVE)
      NAERO(NR) = N-1
      CALL SORTAR(NAERO(NR),XIAERO(1,NR),AERODATA(1,1,NR),NDX)
      LCHANGE = .TRUE.
      GO TO 900
C
C-------------------------------------------------------------
C--- Write aero data to file - xrotor format
C
 240  SAVFIL = ' '
      CALL OPFILE(LUSAVE,SAVFIL)
      DO N = 1, NAERO(NR)
       CALL GETAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &              CMCON,MCRIT,TOC,REREF,REXP)
       A0DEG = A0/DTR
       WRITE(LUSAVE,1420) N,XISECT
       WRITE(LUSAVE,1500) A0DEG,CDMIN
       WRITE(LUSAVE,1510) DCLDA,CLDMIN
       WRITE(LUSAVE,1515) DCLDA_STALL,DCDCL2
       WRITE(LUSAVE,1520) CLMAX,REREF
       WRITE(LUSAVE,1530) CLMIN,REXP
       WRITE(LUSAVE,1540) DCL_STALL,CMCON
       WRITE(LUSAVE,1550) MCRIT
      END DO
      CLOSE(LUSAVE)
      GO TO 900
C
C-------------------------------------------------------------
C---- plotting: get section and load data
C
 60   IF(NAERO(NR).GT.1) THEN
        IF(NINPUT.GE.1) THEN
          N = IINPUT(1)
        ELSE
          WRITE(*,1400)
          READ(*,1) LINE
          IF(LINE.EQ.' ') GO TO 900
          READ(LINE,*,ERR=60) N
        ENDIF
        IF(N.LT.0 .OR. N.GT.NAERO(NR)) GO TO 900
       ELSE
        N = 1
      ENDIF
      CALL GETAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &             DCLDA,DCLDA_STALL,DCL_STALL,
     &             CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &             CMCON,MCRIT,TOC,REREF,REXP)
C
      IF(.NOT.RFIX) RCONST=REREF
C
C---- prompt for mach or re numbers and plot
C
      IF(COMAND.EQ.'PLOM') THEN
        CALL PLOTMACH (NR,N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,CDMIN,CLDMIN,
     &              DCDCL2,CMCON,MCRIT,REREF,REXP,DCDCL2S)
      ELSE
        CALL PLOTREYN (NR,N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,CDMIN,CLDMIN,
     &              DCDCL2,CMCON,MCRIT,REREF,REXP,DCDCL2S)
      ENDIF
C
      GO TO 900
C
C-------------------------------------------------------------
C---- set constant Mach for Re plots
C
 80   IF(NINPUT.GE.1) THEN
        MCONST = RINPUT(1)
      ELSE
        WRITE(*,2100) MCONST
        CALL READR(1,MCONST,ERROR)
        IF(ERROR) GO TO 80  
      ENDIF
C
      IF(MCONST.LT.0.0.OR.MCONST.GT.0.99) THEN
         WRITE(*,*)'Must be in range 0 --> 0.99'
         GO TO 80
      ENDIF
C
      GO TO 950
C
C--------------------------------------------------------------
C---- set constant Re for Mach plots
C
 90   IF(NINPUT.GE.1) THEN
        RCONST = RINPUT(1)
        RFIX = .TRUE.
      ELSE
        WRITE(*,2200) RCONST
        CALL READR(1,RCONST,ERROR)
        IF(ERROR) GO TO 90
        RFIX = .TRUE.
      ENDIF
C
      IF(RCONST.EQ.0.0) THEN
        RCONST = REREF
        RFIX = .FALSE.
      ENDIF
C
      GO TO 950
C
C-------------------------------------------------------------
 100  IF(LPLOT) THEN
        CALL ANNOT(CHGT)
      ELSE
        WRITE(*,*) 'No active plot to annotate'
      ENDIF
      GO TO 950
C
C-------------------------------------------------------------
 110  IF(LPLOT) CALL PLEND
      LPLOT = .FALSE.
      CALL REPLOT(IDEVPS)
      GO TO 950
C
C-------------------------------------------------------------
C
 1400 FORMAT(/' Enter index of aero section to process: ')
 1420 FORMAT(/' Section',I2,'   r/R = ',F5.3)
 1425 FORMAT(/' Section',I2,'   r/R = ',F5.3,
     &     24X,'Set for negative BGam')
C
C---- old formats
c 1500 FORMAT( 1X,72('=')/
c     &       ' 1) Zero-lift alpha (deg):',F6.2,7X,
c     &       ' 9) dCd/dCl**2           :',F8.4)
c 1510 FORMAT(' 2) dCl/dalpha           :',F7.3,6X,
c     &       '10) Secondary dCd/dCl**2 :',F8.4)
c 1515 FORMAT(' 3) dCl/dalpha@stall     :',F7.3,6X,
c     &       '11) Reference Re number  :',F9.0)
c 1520 FORMAT(' 4) Maximum Cl           :',F6.2,7X,
c     &       '12) Re scaling exponent  :',F8.4)
c 1530 FORMAT(' 5) Minimum Cl           :',F6.2,7X,
c     &       '13) Cm                   :',F7.3)
c 1540 FORMAT(' 6) Cl increment to stall:',F7.3,6X,
c     &       '14) Mcrit                :',F7.3)
c 1550 FORMAT(' 7) Minimum Cd           :',F8.4,5X,
c     &       '15) Section t/c          :',F7.3)
c 1560 FORMAT(' 8) Cl at minimum Cd     :',F7.3/
c     &         1X,72('='))
C
C---- the following formats identical to aero/xrotor
C
 1500 FORMAT( 1X,68('=')/
     &       ' Zero-lift alpha (deg):',F7.2,8X,
     &        'Minimum Cd           :',F7.4)
 1510 FORMAT(' d(Cl)/d(alpha)       :',F7.3,8X,
     &        'Cl at minimum Cd     :',F6.3)
 1515 FORMAT(' d(Cl)/d(alpha)@stall :',F7.3,8X,
     &        'd(Cd)/d(Cl**2)       :',F7.4)
 1520 FORMAT(' Maximum Cl           :',F6.2,9X,
     &        'Reference Re number  :',F9.0)
 1530 FORMAT(' Minimum Cl           :',F6.2,9X,
     &        'Re scaling exponent  :',F8.4)
 1540 FORMAT(' Cl increment to stall:',F7.3,8X,
     &        'Cm                   :',F7.3)
 1550 FORMAT('                       ',15X,
     &        'Mcrit                :',F7.3/
     &         1X,68('='))
C
 2100 FORMAT(/' Current value: ',F5.3, 
     &       /' Enter Mach constant for Re plotting')
 2200 FORMAT(/' Current value: ',E9.3,
     &       /' Enter Re constant for Mach plotting (0 for REref)')
C
      END ! AERO



      SUBROUTINE AEROEDIT(NR,NSEC)
C--------------------------------------------------------------
C     Edit airfoil section properties
C--------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      CHARACTER*4 COMAND
      CHARACTER*132 COMARG,LINE
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR, LANS, LNUMCMD, RFIX
      PARAMETER (NMACH = 10)
      PARAMETER (NREYN = 10)
      DIMENSION AMACH(NMACH),AREYN(NREYN)
      REAL MCONST
C
      COMMON /AFPLOTS/ AMACH,AREYN,NM,NRE,MCONST,RCONST,RFIX
C
 1    FORMAT(A)
C
      IF(NSEC.LT.0 .OR. NSEC.GT.NAERO(NR)) THEN
        WRITE(*,*) 'AEROEDIT: section index out of bounds: ',NSEC
        RETURN
      ENDIF
C
      CALL GETAERO(NR,NSEC,XISECT,A0,CLMAX,CLMIN,
     &             DCLDA,DCLDA_STALL,DCL_STALL,
     &             CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &             CMCON,MCRIT,TOC,REREF,REXP)
C
      IF(.NOT.RFIX) RCONST=REREF
      GO TO 10
C
C
 900  CALL ASKC('.EDIT^',COMAND,COMARG)
C
C---- extract numbers (if any) from command argument string
      DO I=1, 20
        IINPUT(I) = 0
        RINPUT(I) = 0.0
      ENDDO
      NINPUT = 20
      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
      NINPUT = 20
      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C
C
      IF(COMAND.EQ.'    ') THEN
       CALL PUTAERO(NR,NSEC,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &              CMCON,MCRIT,TOC,REREF,REXP)
C
C--- Sort aero data sections by XI location
C
       CALL SORTAR(NAERO(NR),XIAERO(1,NR),AERODATA(1,1,NR),NDX)
C
c       CALL SORTAR(NAERO(NR),XIAERO(NSEC,NR),AERODATA(1,1,NR),NDX)
C
       IF(LCHANGE) CALL SETIAERO
C
       RETURN
      ENDIF
C
      IF(COMAND.EQ.'?   ') WRITE(*,1100)
      IF(COMAND.EQ.'?   ') GO TO 900
      IF(COMAND.EQ.'DISP'.OR.COMAND.EQ.'D   ') GO TO 10
      IF(COMAND.EQ.'MOVE') GO TO 20
      IF(COMAND.EQ.'FLIP') GO TO 25
      IF(COMAND.EQ.'LIFT') GO TO 30
      IF(COMAND.EQ.'DRAG') GO TO 40
      IF(COMAND.EQ.'PLOM') GO TO 60
      IF(COMAND.EQ.'PLOR') GO TO 70
      IF(COMAND.EQ.'SETM') GO TO 80
      IF(COMAND.EQ.'SETR') GO TO 90
      IF(COMAND.EQ.'ANNO') GO TO 100
      IF(COMAND.EQ.'HARD') GO TO 110
C
      LNUMCMD = .TRUE.
      IF(COMAND.EQ.'1')    GO TO 31
      IF(COMAND.EQ.'2')    GO TO 32
      IF(COMAND.EQ.'3')    GO TO 33
      IF(COMAND.EQ.'4')    GO TO 34
      IF(COMAND.EQ.'5')    GO TO 35
      IF(COMAND.EQ.'6')    GO TO 36
      IF(COMAND.EQ.'7')    GO TO 41
      IF(COMAND.EQ.'8')    GO TO 42
      IF(COMAND.EQ.'9')    GO TO 43
      IF(COMAND.EQ.'10')   GO TO 44
      IF(COMAND.EQ.'11')   GO TO 45
      IF(COMAND.EQ.'12')   GO TO 37
      IF(COMAND.EQ.'13')   GO TO 46
C
      WRITE(*,1000) COMAND
      GO TO 900
C
 1000 FORMAT(1X,A4,' command not recognized.' //
     &             '  Type "?" for list, <Return> to exit menu.')
C
 1100 FORMAT(/'   <i>     Change individual parameter'
     &       /'   DISP<d> Display section aero characteristics'
     &       /'   LIFT    Change section lift characteristics'
     &       /'   DRAG    Change section drag characteristics'
     &       /'   MOVE    Change section r/R location'
     &       /'   FLIP    Toggle parameters for negative BGam'
C
     &      //'   PLOM    Plot airfoil polars at multiple Mach'
     &       /'   PLOR    Plot airfoil polars at multiple Re'
     &       /'   SETM r  Set constant Mach for Re plots'
     &       /'   SETR r  Set constant Re for Mach plots'
     &       /'   ANNO    Annotate current plot'
     &       /'   HARD    Hardcopy current plot')
C
C
C-------------------------------------------------------------
C--- Display aero data
C
 10   IF((CLMAX+CLMIN).LE.0.0 .AND. A0.GE.0.0) THEN
        WRITE(*,1425) NSEC,XISECT
      ELSE
        WRITE(*,1420) NSEC,XISECT
      ENDIF
C
      A0DEG = A0/DTR
c      WRITE(*,1500) A0DEG,DCDCL2
c      WRITE(*,1510) DCLDA,DCDCL2S
c      WRITE(*,1515) DCLDA_STALL,REREF
c      WRITE(*,1520) CLMAX,REXP
c      WRITE(*,1530) CLMIN,CMCON
c      WRITE(*,1540) DCL_STALL,MCRIT
c      WRITE(*,1550) CDMIN,TOC 
c      WRITE(*,1560) CLDMIN
C
      WRITE(*,1500) A0DEG,CDMIN
      WRITE(*,1510) DCLDA,CLDMIN
      WRITE(*,1515) DCLDA_STALL,DCDCL2
      WRITE(*,1520) CLMAX,REREF
      WRITE(*,1530) CLMIN,REXP
      WRITE(*,1540) DCL_STALL,CMCON
      WRITE(*,1550) MCRIT
      GO TO 900
C
C-------------------------------------------------------------
C--- Change section r/R location
 20   IF(NSEC.EQ.1) THEN
        WRITE(*,*) '*** Cannot move section #1'
        GO TO 900
      ENDIF
C
 21   WRITE(*,1220) NSEC,XISECT
      READ(*,1) LINE
      IF(LINE.EQ.' ') GO TO 900
      READ(LINE,*,ERR=21) XISECT
      IF(XISECT.LT.0.0 .OR. XISECT.GT.1.0) THEN
        WRITE(*,*) '*** Section r/R must be in range 0->1'
        GO TO 21
      ENDIF
      LCHANGE = .TRUE.
      GO TO 10
C
 1220 FORMAT(/'Section # ',I3,' @ r/R = ',F6.3,
     &       /'Enter new r/R for section: ')
C
C-------------------------------------------------------------
C--- Invert parameters for negative BGam
C
 25   A0    = -A0
      CLMAXS=  CLMAX
      CLMAX = -CLMIN
      CLMIN = -CLMAXS
      CLDMIN= -CLDMIN
      CMCON = -CMCON
C
      GO TO 10
C
C-------------------------------------------------------------
C--- Edit lift aero data
C-------------------------------------------------------------
 30   LNUMCMD = .FALSE.
C
 31   A0DEG = A0/DTR
      WRITE(*,2031) A0DEG
 2031 FORMAT(/' Zero-lift alpha (deg):',F9.4)
      CALL READR(1,A0DEG,ERROR)
      IF(ERROR) GO TO 31
      A0 = A0DEG * DTR
      IF(LNUMCMD) GO TO 39
C
 32   WRITE(*,2032) DCLDA
 2032 FORMAT(/' d(Cl)/d(alpha) (/rad):',F9.5)
      CALL READR(1,DCLDA,ERROR)
      IF(ERROR) GO TO 32
      IF(LNUMCMD) GO TO 39

 33   WRITE(*,2033) DCLDA_STALL
 2033 FORMAT(/' d(Cl)/d(alpha) stall (/rad):',F9.5)
      CALL READR(1,DCLDA_STALL,ERROR)
      IF(ERROR) GO TO 33
      IF(LNUMCMD) GO TO 39
C
 34   WRITE(*,2034) CLMAX
 2034 FORMAT(/' Maximum Cl:',F9.4)
      CALL READR(1,CLMAX,ERROR)
      IF(ERROR) GO TO 34
      IF(LNUMCMD) GO TO 39
C
 35   WRITE(*,2035) CLMIN
 2035 FORMAT(/' Minimum Cl:',F9.4)
      CALL READR(1,CLMIN,ERROR)
      IF(ERROR) GO TO 34
      IF(LNUMCMD) GO TO 39
C
 36   WRITE(*,2036) DCL_STALL
 2036 FORMAT(/' Cl increment to stall:',F9.5)
      CALL READR(1,DCL_STALL,ERROR)
      IF(ERROR) GO TO 36
      IF(LNUMCMD) GO TO 39
C
 37   WRITE(*,2037) CMCON
 2037 FORMAT(/' Cm:',F9.5)
      CALL READR(1,CMCON,ERROR)
      IF(ERROR) GO TO 37
C
 39   LCHANGE = .TRUE.
      GO TO 10
C
C-------------------------------------------------------------
C--- Edit drag aero data
C-------------------------------------------------------------
 40   LNUMCMD = .FALSE.
C
 41   WRITE(*,2041) CDMIN
 2041 FORMAT(/' Minimum Cd:',F10.6)
      CALL READR(1,CDMIN,ERROR)
      IF(ERROR) GO TO 41
      IF(LNUMCMD) GO TO 50
C
   42 WRITE(*,2042) CLDMIN
 2042 FORMAT(/' Cl at minimum Cd:',F10.5)
      CALL READR(1,CLDMIN,ERROR)
      IF(ERROR) GO TO 42
      IF(LNUMCMD) GO TO 50
C
   43 WRITE(*,2043) DCDCL2
 2043 FORMAT(/' d(Cd)/d(Cl**2):',F10.6)
      CALL READR(1,DCDCL2,ERROR)
      IF(ERROR) GO TO 43
      IF(LNUMCMD) GO TO 50
C
C----CODE CHANGE - DCDCL2S WAS CDSEC in 070
C
c 48   WRITE(*,2048) DCDCL2S
c 2048 FORMAT(/' Secondary drag dCd/dCL2:',F9.5)
c      CALL READR(1,DCDCL2S,ERROR)
c      IF(ERROR) GO TO 48
c      IF(LNUMCMD) GO TO 50
C
   44 WRITE(*,2044) REREF
 2044 FORMAT(/' Reference Re number:',F12.0)
      CALL READR(1,REREF,ERROR)
      IF(ERROR) GO TO 44
      IF(LNUMCMD) GO TO 50
C
   45 WRITE(*,2045) REXP
 2045 FORMAT(/' Re scaling exponent:',F10.4)
      CALL READR(1,REXP,ERROR)
      IF(ERROR) GO TO 45
      IF(LNUMCMD) GO TO 50
C
 46   WRITE(*,2046) MCRIT
 2046 FORMAT(/' Mcrit:',F9.5)
      CALL READR(1,MCRIT,ERROR)
      IF(ERROR) GO TO 46
c      IF(LNUMCMD) GO TO 50
C
c 47   WRITE(*,2047) TOC
c 2047 FORMAT(/' Section t/c:',F9.5)
c      CALL READR(1,TOC,ERROR)
c      IF(ERROR) GO TO 47
C
 50   LCHANGE = .TRUE.
      GO TO 10
C
C-------------------------------------------------------------
C---- plot multiple machs
C
 60   CALL PLOTMACH (NR,NSEC,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,CDMIN,CLDMIN,
     &              DCDCL2,CMCON,MCRIT,REREF,REXP,DCDCL2S)
      GO TO 900
C
C-------------------------------------------------------------
C---- plot multiple rns 

 70   CALL PLOTREYN (NR,NSEC,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,CDMIN,CLDMIN,
     &              DCDCL2,CMCON,MCRIT,REREF,REXP,DCDCL2S)
      GO TO 900
C
C-------------------------------------------------------------
C---- set constant Mach for Re plots
C
 80   IF(NINPUT.GE.1) THEN
        MCONST = RINPUT(1)
      ELSE
        WRITE(*,2100) MCONST
        CALL READR(1,MCONST,ERROR)
        IF(ERROR) GO TO 80  
      ENDIF
C
      IF(MCONST.LT.0.0.OR.MCONST.GT.0.99) THEN
         WRITE(*,*)'Must be in range 0 --> 0.99'
         GO TO 80
      ENDIF
C
      GO TO 900
C
C--------------------------------------------------------------
C---- set constant Re for Mach plots
C
 90   IF(NINPUT.GE.1) THEN
        RCONST = RINPUT(1)
        RFIX = .TRUE.
      ELSE
        WRITE(*,2200) RCONST
        CALL READR(1,RCONST,ERROR)
        IF(ERROR) GO TO 90
        RFIX = .TRUE.
      ENDIF
C
      IF(RCONST.EQ.0.0) THEN
        RCONST = REREF
        RFIX = .FALSE.
      ENDIF
C
      GO TO 900
C
C-------------------------------------------------------------
 100  IF(LPLOT) THEN
        CALL ANNOT(CHGT)
      ELSE
        WRITE(*,*) 'No active plot to annotate'
      ENDIF
      GO TO 900
C
C-------------------------------------------------------------
 110  IF(LPLOT) CALL PLEND
      LPLOT = .FALSE.
      CALL REPLOT(IDEVPS)
      GO TO 900
C
C--------------------------------------------------------------
C
 1400 FORMAT(/'Enter index of aero section to process: ')
 1420 FORMAT(/' Section',I2,'   r/R = ',F5.3)
 1425 FORMAT(/' Section',I2,'   r/R = ',F5.3,28X,
     &         'Set for negative BGam')
C
C---- old formats
c 1500 FORMAT( 1X,72('-')/
c     &       ' 1) Zero-lift alpha (deg):',F6.2,7X,
c     &       ' 9) dCd/dCl**2           :',F8.4)
c 1510 FORMAT(' 2) dCl/dalpha           :',F7.3,6X,
c     &       '10) Secondary dCd/dCl**2 :',F8.4)
c 1515 FORMAT(' 3) dCl/dalpha@stall     :',F7.3,6X,
c     &       '11) Reference Re number  :',F9.0)
c 1520 FORMAT(' 4) Maximum Cl           :',F6.2,7X,
c     &       '12) Re scaling exponent  :',F8.4)
c 1530 FORMAT(' 5) Minimum Cl           :',F6.2,7X,
c     &       '13) Cm                   :',F7.3)
c 1540 FORMAT(' 6) Cl increment to stall:',F7.3,6X,
c     &       '14) Mcrit                :',F7.3)
c 1550 FORMAT(' 7) Minimum Cd           :',F8.4,5X,
c     &       '15) Section t/c          :',F7.3)
c 1560 FORMAT(' 8) Cl at minimum Cd     :',F7.3/
c     &         1X,72('-'))
C
 1500 FORMAT( 1X,72('=')/
     &       ' 1) Zero-lift alpha (deg):',F7.2,6X,
     &       ' 7) Minimum Cd           :',F7.4)
 1510 FORMAT(' 2) d(Cl)/d(alpha)       :',F7.3,6X,
     &       ' 8) Cl at minimum Cd     :',F6.3)
 1515 FORMAT(' 3) d(Cl)/d(alpha)@stall :',F7.3,6X,
     &       ' 9) d(Cd)/d(Cl**2)       :',F7.4)
 1520 FORMAT(' 4) Maximum Cl           :',F6.2,7X,
     &       '10) Reference Re number  :',F9.0)
 1530 FORMAT(' 5) Minimum Cl           :',F6.2,7X,
     &       '11) Re scaling exponent  :',F8.4)
 1540 FORMAT(' 6) Cl increment to stall:',F7.3,6X,
     &       '12) Cm                   :',F7.3)
 1550 FORMAT('                          ',13X,
     &       '13) Mcrit                :',F7.3/
     &         1X,72('='))
C
 2100 FORMAT(/' Current value: ',F5.3, 
     &       /' Enter Mach constant for Re plotting')
 2200 FORMAT(/' Current value: ',E9.3,
     &       /' Enter Re constant for Mach plotting (0 for REref)')
C
      RETURN
      END ! AEROEDIT



      SUBROUTINE SETIAERO
C--------------------------------------------------
C     Sets up indices referring to aero section for 
C     each radial station
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C--- Find lower index of aero data sections XIAERO(N) bounding XI=YRC/RTIP
      DO NR = 1, NROTOR
       DO I=1, NRC
        IAERO(I,NR) = 1
        DO N = 1, NAERO(NR)
         XI = YRC(I,NR)/RTIP(NR)
         IF(XIAERO(N,NR).LE.XI) THEN
           IAERO(I,NR) = N
         ENDIF
        END DO
       END DO
      END DO
      RETURN
      END



      SUBROUTINE GETAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &                   DCLDA,DCLDA_STALL,DCL_STALL,
     &                   CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &                   CMCON,MCRIT,TOC,REREF,REXP)
C---------------------------------------------
C     Gets aero data from stored section array
C
C   AERODATA    Aerodynamic definition of the blade section at XIAERO
C               AERODATA( 1,x) = A0 (angle of zero lift)
C               AERODATA( 2,x) = CLMAX (Max CL)
C               AERODATA( 3,x) = CLMIN (Min CL)
C               AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
C               AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
C               AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
C               AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
C               AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
C               AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
C               AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
C               AERODATA(11,x) = REREF (reference Reynold's number)
C               AERODATA(12,x) = REXP (Reynold's number exponent Cd~Re^REXP)
C               AERODATA(13,x) = MCRIT (critical Mach #)
C               AERODATA(14,x) = TOC (thickness/chord)
C               AERODATA(15,x) = DCDCL2S (Secondary, annulus drag param d(Cd)/dCL^2)
C---------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(NR.LT.1 .OR. NR.GT.NROTOR) THEN
        WRITE(*,*) 'Error: blade index of aero section out of bounds'
        RETURN
      ENDIF
      IF(N.LT.1 .OR. N.GT.NAERO(NR)) THEN
        WRITE(*,*) 'Error: index of aero section out of bounds'
        RETURN
      ENDIF
C
      A0          = AERODATA( 1,N,NR)
      CLMAX       = AERODATA( 2,N,NR)
      CLMIN       = AERODATA( 3,N,NR)
      DCLDA       = AERODATA( 4,N,NR)
      DCLDA_STALL = AERODATA( 5,N,NR)
      DCL_STALL   = AERODATA( 6,N,NR)
      CDMIN       = AERODATA( 7,N,NR)
      CLDMIN      = AERODATA( 8,N,NR)
      DCDCL2      = AERODATA( 9,N,NR)
      CMCON       = AERODATA(10,N,NR)
      REREF       = AERODATA(11,N,NR)
      REXP        = AERODATA(12,N,NR)
      MCRIT       = AERODATA(13,N,NR)
      TOC         = AERODATA(14,N,NR)
      DCDCL2S     = AERODATA(15,N,NR)
      XISECT      = XIAERO(N,NR)
C
      RETURN
      END


      SUBROUTINE PUTAERO(NR,N,XISECT,A0,CLMAX,CLMIN,
     &                   DCLDA,DCLDA_STALL,DCL_STALL,
     &                   CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &                   CMCON,MCRIT,TOC,REREF,REXP)
C--------------------------------------------------------
C     Puts aero data into stored section array at index N
C
C   AERODATA    Aerodynamic definition of the blade section at XIAERO
C               AERODATA( 1,x) = A0 (angle of zero lift)
C               AERODATA( 2,x) = CLMAX (Max CL)
C               AERODATA( 3,x) = CLMIN (Min CL)
C               AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
C               AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
C               AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
C               AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
C               AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
C               AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
C               AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
C               AERODATA(11,x) = REREF (reference Reynold's number)
C               AERODATA(12,x) = REXP (Reynold's number exponent Cd~Re^REXP)
C               AERODATA(13,x) = MCRIT (critical Mach #)
C               AERODATA(14,x) = TOC (thickness/chord)
C               AERODATA(15,x) = DCDCL2S (Secondary, annulus drag param d(Cd)/dCL^2)
C--------------------------------------------------------
      INCLUDE 'DFDC.INC'
C 
      IF(NR.LT.1 .OR. NR.GT.NRX) THEN
        WRITE(*,*) 'Error: blade index of aero section out of bounds'
        RETURN
      ENDIF
      IF(N.LT.1) THEN
        WRITE(*,*) 'Error: index of aero section out of bounds'
        RETURN
      ENDIF
      IF(N.GT.NAX) THEN
        WRITE(*,*) 'Too many aero sections defined...'
        RETURN
      ENDIF
C
      AERODATA( 1,N,NR) = A0
      AERODATA( 2,N,NR) = CLMAX
      AERODATA( 3,N,NR) = CLMIN
      AERODATA( 4,N,NR) = DCLDA
      AERODATA( 5,N,NR) = DCLDA_STALL
      AERODATA( 6,N,NR) = DCL_STALL
      AERODATA( 7,N,NR) = CDMIN  
      AERODATA( 8,N,NR) = CLDMIN
      AERODATA( 9,N,NR) = DCDCL2
      AERODATA(10,N,NR) = CMCON 
      AERODATA(11,N,NR) = REREF
      AERODATA(12,N,NR) = REXP
      AERODATA(13,N,NR) = MCRIT
      AERODATA(14,N,NR) = TOC
      AERODATA(15,N,NR) = DCDCL2S
      XIAERO(N,NR)      = XISECT
C
      RETURN
      END



      SUBROUTINE SORTAR(NS,S,W,NDIM)
C----------------------------------------------------
C---- sort arrays by S values
C     Orders data monotonically increasing in S(i)
C----------------------------------------------------
      DIMENSION S(NS), W(NDIM,NS)
      LOGICAL DONE
C
      DO IPASS=1, 500
        DONE = .TRUE.
        DO N=1, NS-1
          NP = N+1
          IF(S(NP).LT.S(N)) THEN
           TEMP  = S(NP)
           S(NP) = S(N)
           S(N)  = TEMP
           DO L = 1, NDIM
             TEMP    = W(L,NP)
             W(L,NP) = W(L,N)
             W(L,N)  = TEMP
           END DO
           DONE = .FALSE.
          ENDIF
        END DO
        IF(DONE) GO TO 10
      END DO
      STOP 'SORTAR failed'
C
 10   RETURN
      END ! SORTAR


C*************************************************************************
C  Interpolated aero section properties functions
C  These routines implement a functional representation of the 
C  blade aero properties (CL,CD,CM) vs ALFA
C*************************************************************************


      SUBROUTINE GETCLCDCM(NR,IS,XI,ALF,W,REY,SECSIG,SECSTAGR,
     &                     CLIFT,CL_ALF,CL_W,
     &                     CLMAX,CLMIN,DCL_STALL,STALLF,
     &                     CDRAG,CD_ALF,CD_W,CD_REY,
     &                     CMOM,CM_AL,CM_W)
C-------------------------------------------------------------
C     CL(alpha),
C      CD(alpha), 
C       CM(alpha) interpolation function for blade at station IS at XI=r/R
C-------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL STALLF,STALLF2
C
      IF(XI.LT.0.0 .OR. XI.GT.1.0) THEN
        WRITE(*,*) 'Undefined section XI in GETCLCDCM ',XI
      ENDIF
C
C--- Check for installed aero data section index
      N = IAERO(IS,NR)
      IF(N.LT.1 .OR. N.GT.NAERO(NR)) THEN
C
       IF(NAERO(NR).GT.1) THEN
C--- Find lower index of aero data sections XIAERO(N) bounding XI
        DO N = 1, NAERO(NR)
         IF(XIAERO(N,NR).LE.XI) THEN
cc          write(*,*) 'getcl iaero= ',N,' is= ',is,xiaero(N),xi
           IAERO(IS,NR) = N
          ELSE
           GO TO 10
         ENDIF
        END DO
        WRITE(*,*) 'Aero section not found for station ',XI
       ENDIF
C
       N = 1
       IAERO(IS,NR) = N
      ENDIF
C
C--- Get section aero data from stored section array
 10   A0          = AERODATA( 1,N,NR)
      CLMAX       = AERODATA( 2,N,NR)
      CLMIN       = AERODATA( 3,N,NR)
      DCLDA       = AERODATA( 4,N,NR)
      DCLDA_STALL = AERODATA( 5,N,NR)
      DCL_STALL   = AERODATA( 6,N,NR)
      CDMIN       = AERODATA( 7,N,NR)
      CLDMIN      = AERODATA( 8,N,NR)
      DCDCL2      = AERODATA( 9,N,NR)
      CMCON       = AERODATA(10,N,NR)
      REREF       = AERODATA(11,N,NR)
      REXP        = AERODATA(12,N,NR)
      MCRIT       = AERODATA(13,N,NR)
      TOC         = AERODATA(14,N,NR)
      DCDCL2S     = AERODATA(15,N,NR)
      XISECT1     = XIAERO(N,NR)
C
C--- Get data for inner bounding aero section
      CALL CLCDCM(ALF,W,REY,VSO,SECSIG,SECSTAGR,
     &            CLIFT,CL_ALF,CL_W,STALLF,
     &            CDRAG,CD_ALF,CD_W,CD_REY,
     &            CMOM,CM_AL,CM_W,
     &            A0,CLMAX,CLMIN,DCLDA,DCLDA_STALL,DCL_STALL,
     &            CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP,
     &            TOC,DCDCL2S)
C
C--- Check for another bounding section, if not we are done, 
C    if we have another section linearly interpolate data to station IS
      IF(N.LT.NAERO(NR)) THEN
        XISECT2 = XIAERO(N+1,NR)
        FRAC = (XI-XISECT1)/(XISECT2-XISECT1)
        IF(FRAC.LE.0.0 .OR. FRAC.GT.1.0) THEN
cc         write(*,*) 'CL n,is,xi,frac = ',n,is,xi(is),frac
        ENDIF
C
C--- A02 sustituted for A0 in the following (2 places),
C    to permit A0 interpolation for storage in AZERO
C
        A02         = AERODATA( 1,N+1,NR)
        CLMAX2      = AERODATA( 2,N+1,NR)
        CLMIN2      = AERODATA( 3,N+1,NR)
        DCLDA       = AERODATA( 4,N+1,NR)
        DCLDA_STALL = AERODATA( 5,N+1,NR)
        DCL_STALL2  = AERODATA( 6,N+1,NR)
        CDMIN       = AERODATA( 7,N+1,NR)
        CLDMIN      = AERODATA( 8,N+1,NR)
        DCDCL2      = AERODATA( 9,N+1,NR)
        CMCON       = AERODATA(10,N+1,NR)
        REREF       = AERODATA(11,N+1,NR)
        REXP        = AERODATA(12,N+1,NR)
        MCRIT       = AERODATA(13,N+1,NR)
        TOC         = AERODATA(14,N+1,NR)
        DCDCL2S     = AERODATA(15,N+1,NR)
C
C--- Get data for outer bounding aero section
        CALL CLCDCM(ALF,W,REY,VSO,SECSIG,SECSTAGR,
     &              CLIFT2,CL_ALF2,CL_W2,STALLF2,
     &              CDRAG2,CD_ALF2,CD_W2,CD_REY2,
     &              CMOM2,CM_AL2,CM_W2,
     &              A02,CLMAX2,CLMIN2,DCLDA,DCLDA_STALL,DCL_STALL2,
     &              CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP,
     &              TOC,DCDCL2S)
C--- Interpolate aero data to blade station
        STALLF = STALLF .OR. STALLF2
        CLIFT  = (1.0-FRAC)*CLIFT  + FRAC*CLIFT2
        CL_ALF = (1.0-FRAC)*CL_ALF + FRAC*CL_ALF2
        CL_W   = (1.0-FRAC)*CL_W   + FRAC*CL_W2
        CLMAX  = (1.0-FRAC)*CLMAX  + FRAC*CLMAX2
        CLMIN  = (1.0-FRAC)*CLMIN  + FRAC*CLMIN2
        DCL_STALL = (1.0-FRAC)*DCL_STALL + FRAC*DCL_STALL2
C
        CMOM   = (1.0-FRAC)*CMOM   + FRAC*CMOM2
        CM_AL  = (1.0-FRAC)*CM_AL  + FRAC*CM_AL2
        CM_W   = (1.0-FRAC)*CM_W   + FRAC*CM_W2
C
        CDRAG  = (1.0-FRAC)*CDRAG  + FRAC*CDRAG2
        CD_ALF = (1.0-FRAC)*CD_ALF + FRAC*CD_ALF2
        CD_W   = (1.0-FRAC)*CD_W   + FRAC*CD_W2
        CD_REY = (1.0-FRAC)*CD_REY + FRAC*CD_REY2
        A0     = (1.0-FRAC)*A0     + FRAC*A02
      ENDIF
C
      AZERO(IS,NR) = A0
C
      RETURN
      END  !GETCLCDCM



      SUBROUTINE GETALF(NR,IS,XI,SECSIG,SECSTAGR,
     &                  CLIFT,W,ALF,ALF_CL,ALF_W,STALLF)
C------------------------------------------------------------
C     Inverse alpha(CL) function 
C     Uses Newton-Raphson iteration to get ALF from CL function
C------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL STALLF
      DATA NITER / 10 /
      DATA EPS   / 1.0E-5 /
C
      STALLF = .FALSE.
C
C---HHY A0 is now an aero section property
      A0  = AERODATA(1,IS,NR)
      REY = 0.0
C
      ALF = A0 
      DO ITER=1, NITER
        CALL GETCLCDCM(NR,IS,XI,ALF,W,REY,SECSIG,SECSTAGR,
     &                 CLTEMP,CL_ALF,CL_W,
     &                 CLMAX,CLMIN,DCL_STALL,STALLF,
     &                 CDRAG,CD_ALF,CD_W,CD_REY,
     &                 CMOM,CM_AL,CM_W)
cc      IF(STALLF) GO TO 20
        DALF = -(CLTEMP-CLIFT)/CL_ALF
        ALF = ALF + DALF
        ALF_CL =   1.0/CL_ALF
        ALF_W  = -CL_W/CL_ALF
        IF(ABS(DALF).LT.EPS) RETURN
      END DO
C
   20 WRITE(*,*) 'GETALF: alpha(CL) function inversion failed'
c      write(*,*) 'is,clift  ',is,clift
c      write(*,*) 'abs(dalf) ',abs(dalf)
c      write(*,*) 'cl_alf    ',cl_alf
C
      RETURN
      END ! GETALF



C*************************************************************************
C  Basic aero section properties functions
C  These routines implement a functional representation of the 
C  blade section aero properties (CL,CD,CM) vs ALF
C*************************************************************************

      SUBROUTINE CLCDCM(ALF,W,REY,VSO,SECSIG,SECSTAGR,
     &                  CLIFT,CL_ALF,CL_W,STALLF,
     &                  CDRAG,CD_ALF,CD_W,CD_REY,
     &                  CMOM,CM_AL,CM_W,
     &                  A0,CLMAX,CLMIN,DCLDA,DCLDA_STALL,DCL_STALL,
     &                  CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP,
     &                  TOC,DCDCL2S)
C------------------------------------------------------------
C     CL(alpha) function
C     Note that in addition to setting CLIFT and its derivatives
C     CLMAX and CLMIN (+ and - stall CL's) are set in this routine
C     In the compressible range the stall CL is reduced by a factor
C     proportional to Mcrit-Mach.  Stall limiting for compressible 
C     cases begins when the compressible drag added CDC > CDMstall
C------------------------------------------------------------
C     CD(alpha) function - presently CD is assumed to be a sum
C     of profile drag + stall drag + compressibility drag
C     In the linear lift range drag is CD0 + quadratic function of CL-CLDMIN
C     In + or - stall an additional drag is added that is proportional
C     to the extent of lift reduction from the linear lift value.
C     Compressible drag is based on adding drag proportional to 
C     (Mach-Mcrit_eff)^MEXP
C------------------------------------------------------------
C     CM(alpha) function - presently CM is assumed constant,
C     varying only with Mach by Prandtl-Glauert scaling
C------------------------------------------------------------
C
C      INCLUDE 'DFDC.INC'
      IMPLICIT REAL (M)
      LOGICAL STALLF
      DOUBLE PRECISION ECMIN, ECMAX
C
C---- Factors for compressibility drag model, HHY 10/23/00
C     Mcrit is set by user
C     Effective Mcrit is Mcrit_eff = Mcrit - CLMFACTOR*(CL-CLDmin) - DMDD
C     DMDD is the delta Mach to get CD=CDMDD (usually 0.0020)
C     Compressible drag is CDC = CDMFACTOR*(Mach-Mcrit_eff)^MEXP
C     CDMstall is the drag at which compressible stall begins
C
      CDMFACTOR = 10.0
      CLMFACTOR =  0.25
      MEXP      =  3.0
      CDMDD     =  0.0020
      CDMSTALL  =  0.1000
C
C---- Prandtl-Glauert compressibility factor
      MSQ   =   W*W/VSO**2
      MSQ_W = 2.0*W/VSO**2
      IF(MSQ.GE.1.0) THEN
       WRITE(*,*)
     &  'CLFUNC: Local Mach number limited to 0.99, was ', MSQ
       MSQ = 0.99
       MSQ_W = 0.
      ENDIF
      PGRT = 1.0 / SQRT(1.0 - MSQ)
      PGRT_W = 0.5*MSQ_W * PGRT**3
C
C---- Mach number and dependence on velocity
      MACH = SQRT(MSQ)
      MACH_W = 0.0
      IF(MACH.NE.0.0) MACH_W = 0.5*MSQ_W/MACH 
C
C------------------------------------------------------------
C--- Generate CLFACTOR for cascade effects from section solidity
      CLFACTOR = 1.0
      IF(SECSIG.GT.0.0) THEN
        CALL GETCLFACTOR(SECSIG,SECSTAGR,CLFACTOR)
      ENDIF
C
C------------------------------------------------------------
C--- Generate CL from dCL/dAlpha and Prandtl-Glauert scaling
      CLA     = DCLDA*PGRT  *(ALF-A0) * CLFACTOR
      CLA_ALF = DCLDA*PGRT            * CLFACTOR
      CLA_W   = DCLDA*PGRT_W*(ALF-A0) * CLFACTOR
C
Cccccccccc
C      WRITE(*,*)'CL Factor   ALF   A0   DCLDA  CLA'
C      WRITE(*,*) CLFACTOR,ALF,A0,DCLDA,CLA
c
CC--- Effective CLmax is limited by Mach effects
C    reduces CLmax to match the CL of onset of serious compressible drag
      CLMX = CLMAX
      CLMN = CLMIN
      DMSTALL  = (CDMSTALL/CDMFACTOR)**(1.0/MEXP)
      CLMAXM = MAX(0.0, (MCRIT+DMSTALL-MACH)/CLMFACTOR) + CLDMIN
      CLMAX  = MIN(CLMAX,CLMAXM)
      CLMINM = MIN(0.0,-(MCRIT+DMSTALL-MACH)/CLMFACTOR) + CLDMIN
      CLMIN  = MAX(CLMIN,CLMINM)
C
C--- CL limiter function (turns on after +-stall 
      ECMAX = DEXP( MIN(200.0D0,DBLE((CLA-CLMAX)/DCL_STALL)) )
      ECMIN = DEXP( MIN(200.0D0,DBLE((CLMIN-CLA)/DCL_STALL)) )
      CLLIM = DCL_STALL * DLOG( (1.0D0+ECMAX)/(1.0D0+ECMIN) )
      CLLIM_CLA = ECMAX/(1.0+ECMAX) + ECMIN/(1.0+ECMIN)
c
c      if(CLLIM.GT.0.001) then
c      write(*,999) 'cla,cllim,ecmax,ecmin ',cla,cllim,ecmax,ecmin
c      endif
c 999  format(a,2(1x,f10.6),3(1x,d12.6))
C
C--- Subtract off a (nearly unity) fraction of the limited CL function
C    This sets the dCL/dAlpha in the stalled regions to 1-FSTALL of that
C    in the linear lift range
      FSTALL = DCLDA_STALL/DCLDA
      CLIFT  = CLA     - (1.0-FSTALL)*CLLIM
      CL_ALF = CLA_ALF - (1.0-FSTALL)*CLLIM_CLA*CLA_ALF
      CL_W   = CLA_W   - (1.0-FSTALL)*CLLIM_CLA*CLA_W
C
      STALLF = .FALSE.
      IF(CLIFT.GT.CLMAX) STALLF = .TRUE.
      IF(CLIFT.LT.CLMIN) STALLF = .TRUE.
C
C
C------------------------------------------------------------
C--- CM from CMCON and Prandtl-Glauert scaling
      CMOM  = PGRT*CMCON
      CM_AL = 0.0
      CM_W  = PGRT_W*CMCON
C
C
C------------------------------------------------------------
C--- CD from profile drag, stall drag and compressibility drag 
C
C---- Reynolds number scaling factor
      IF(REY.LE.0) THEN
       RCORR = 1.0
       RCORR_REY = 0.0
      ELSE
       RCORR     = (REY/REREF)**REXP
       RCORR_REY =  REXP/REY
      ENDIF
C
C--- Include quadratic lift drag terms from airfoil and annulus
C
c      CDCL2 = DCDCL2 + DCDCL2S
      CDCL2 = DCDCL2  ! no chance of getting messed up...
C
C--- In the basic linear lift range drag is a function of lift
C    CD = CD0 (constant) + quadratic with CL)
      CDRAG  = (CDMIN + CDCL2*(CLIFT-CLDMIN)**2    ) * RCORR
      CD_ALF = (    2.0*CDCL2*(CLIFT-CLDMIN)*CL_ALF) * RCORR
      CD_W   = (    2.0*CDCL2*(CLIFT-CLDMIN)*CL_W  ) * RCORR
      CD_REY = CDRAG*RCORR_REY
C
C--- Post-stall drag added
      FSTALL = DCLDA_STALL/DCLDA
      DCDX    = (1.0-FSTALL)*CLLIM/(PGRT*DCLDA)
c      write(*,*) 'cla,cllim,fstall,pg,dclda ',cla,cllim,fstall,pg,dclda
      DCD     = 2.0* DCDX**2
      DCD_ALF = 4.0* DCDX * 
     &         (1.0-FSTALL)*CLLIM_CLA*CLA_ALF/(PGRT*DCLDA)
      DCD_W = 4.0* DCDX * 
     &       ( (1.0-FSTALL)*CLLIM_CLA*CLA_W/(PGRT*DCLDA)
     &          - DCD/PGRT*PGRT_W )
c      write(*,*) 'alf,cl,dcd,dcd_alf,dcd_w ',alf,clift,dcd,dcd_alf,dcd_w
C
C--- Compressibility drag (accounts for drag rise above Mcrit with CL effects
C    CDC is a function of a scaling factor*(M-Mcrit(CL))**MEXP
C    DMDD is the Mach difference corresponding to CD rise of CDMDD at MCRIT
      DMDD = (CDMDD/CDMFACTOR)**(1.0/MEXP)
      CRITMACH = MCRIT-CLMFACTOR*ABS(CLIFT-CLDMIN) - DMDD
      CRITMACH_ALF  = -CLMFACTOR*ABS(CL_ALF)
      CRITMACH_W    = -CLMFACTOR*ABS(CL_W)
      IF(MACH.LT.CRITMACH) THEN
       CDC     = 0.0
       CDC_ALF = 0.0
       CDC_W   = 0.0
      ELSE
       CDC = CDMFACTOR*(MACH-CRITMACH)**MEXP
       CDC_W   = MEXP*MACH_W*CDC/MACH - MEXP*CRITMACH_W  *CDC/CRITMACH
       CDC_ALF =                      - MEXP*CRITMACH_ALF*CDC/CRITMACH
      ENDIF
c      write(*,*) 'critmach,mach ',critmach,mach
c      write(*,*) 'cdc,cdc_w,cdc_alf ',cdc,cdc_w,cdc_alf
C
      FAC   = 1.0
      FAC_W = 0.0
C--- Although test data does not show profile drag increases due to Mach # 
C    you could use something like this to add increase drag by Prandtl-Glauert
C    (or any function you choose) 
cc      FAC   = PG
cc      FAC_W = PG_W
C--- Total drag terms
      CDRAG  = FAC*CDRAG              + DCD     + CDC
      CD_ALF = FAC*CD_ALF             + DCD_ALF + CDC_ALF
      CD_W   = FAC*CD_W + FAC_W*CDRAG + DCD_W   + CDC_ALF
      CD_REY = FAC*CD_REY
C
      RETURN
      END ! CLCDCM



      SUBROUTINE AEROPLT2(NR,NSEC,XISECT,NMACH,AMACH,NREYN,AREYN,
     &                  A0,CLMAX,CLMIN,DCLDA,DCLDA_STALL,DCL_STALL,
     &                  CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP,
     &                  DCDCL2S, MCONST,RCONST,IPL)
C------------------------------------------------------------
C     Plots section characteristics based on parametric model
C------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      PARAMETER ( NTMP=301 )
      LOGICAL LDUMMY, STALLF
      EXTERNAL PLCHAR,PLMATH
      REAL MCONST
      CHARACTER*80 AFNAME
C
      DIMENSION ATMP(NTMP), XTMP(NTMP), YTMP(NTMP)
      DIMENSION XLIN(2), YLIN(2), AMACH(NMACH),AREYN(NREYN)
C
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
C---- plot scale factor and aspect ratio
      PLFAC = 0.8
      PLPAR = 0.6
C
C---- character size for axis numbers, labels
      CS  = CSIZE*0.8
      CSL = CSIZE*1.0
C
C---- set solidity and stagger to 0.0 to give isolated foil CL vs alpha
      SECSIG   = 0.0
      SECSTAGR = 0.0
C
C---- CD, CL, alpha axis annotation increments
      DCD = 0.002
      DCL = 0.2
      DAL = 2.0
C
      DALFLIN = (CLMAX/DCLDA - CLMIN/DCLDA 
     &          + 8.0*DCL_STALL/DCLDA) * 180.0/PI
      IF(CLMAX-CLMIN .GT.  2.01)  DCL =  0.5
      IF(DALFLIN     .GT. 10.01)  DAL =  5.0
      IF(DALFLIN     .GT. 20.01)  DAL = 10.0
C
      IF(5.0*CDMIN   .GT. 0.016)  DCD = 0.005
      IF(5.0*CDMIN   .GT. 0.040)  DCD = 0.01
      IF(5.0*CDMIN   .GT. 0.080)  DCD = 0.02
      IF(5.0*CDMIN   .GT. 0.160)  DCD = 0.05
      IF(5.0*CDMIN   .GT. 0.400)  DCD = 0.10
C
C---- set plot limits
      CLMIN0 = DCL * AINT( MIN(CLMIN,0.0)/DCL - 0.5 )
      CLMAX0 = DCL * AINT( MAX(CLMAX,0.0)/DCL + 0.5 )
C
      CDMIN0 = 0.0
      CDMAX0 = DCD * AINT( 5.0*CDMIN/DCD + 0.5 )
      CDMAX0 = MAX( CDMAX0 , 0.001 )
C
      ALFMI = (CLMIN-8.0*DCL_STALL)/DCLDA
      ALFMA = (CLMAX+8.0*DCL_STALL)/DCLDA
      ALMIN0 = DAL * AINT( (ALFMI-0.02)*180.0/PI / DAL - 0.5 )
      ALMAX0 = DAL * AINT( (ALFMA+0.02)*180.0/PI / DAL + 0.5 )
C
C---- set CL, CD, alpha scaling factors
      CLWT = PLPAR/(CLMAX0-CLMIN0)
      CDWT = 0.60 /(CDMAX0-CDMIN0)
      ALWT = 0.30 /(ALMAX0-ALMIN0)
C
C---- get plot name
      IF(NROTOR.GT.1) THEN
         CALL GETPNAME(AFNAME,NAME,NR)
      ELSE
         AFNAME=NAME
      ENDIF
C
C
      CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFAC*SIZE,LPLOT,LLAND)
      CALL PLOTABS(1.0,0.75,-3)
C
      CALL GETCOLOR(ICOL0)
C
C---- re-origin for CD-CL plot
      CALL PLOT(0.0,-CLWT*CLMIN0,-3)
      CALL PLOT(5.0*CS,0.0,-3)
C
C---- plot case name, section # and r/R location
      CALL NEWPEN(2)
      XPLT = 0.0
      YPLT = CLWT*CLMAX0 + 1.0*CSL
      CALL PLCHAR(XPLT,YPLT,CSL,AFNAME,0.0,-1)
      XPLT = 0.0
      YPLT = CLWT*CLMIN0 - 2.5*CS
      CALL PLCHAR(XPLT,YPLT,CS,'Section ',0.0,8)
      CALL PLNUMB(999.,999.,CS, FLOAT(NSEC),0.0,-1)
      CALL PLCHAR(999.,999.,CS,'  r/R = ',0.0,8)
      CALL PLNUMB(999.,999.,CS, XISECT,0.0,3)
C
      IF(IPL.EQ.1) THEN
        CALL PLCHAR(999.,999.,CS,'  Re = ',0.0,7)
        CALL PLNUMB(999.,999.,CS, RCONST,0.0,-1)
      ELSE
        CALL PLCHAR(999.,999.,CS,'  Mach = ',0.0,9)
        CALL PLNUMB(999.,999.,CS, MCONST,0.0,3)
        CALL PLCHAR(999.,999.,CS,'  REexp = ',0.0,10)
        CALL PLNUMB(999.,999.,CS, REXP,0.0,3)
      ENDIF
C
C---- optional - to remind users this is not taken into account
C
c      CALL PLCHAR(999.,999.,CS,'  Solidity = 0',0.0,14)
C
C---- plot CD-CL axes
      CALL NEWPEN(2)
      CALL YAXIS(0.0,CLWT*CLMIN0, CLWT*(CLMAX0-CLMIN0),
     &               CLWT*DCL,CLMIN0,DCL,CS,1)
      CALL XAXIS(CDWT*CDMIN0,0.0,-CDWT*(CDMAX0-CDMIN0),
     &               CDWT*DCD,CDMIN0,DCD,CS,3)
      IF(LAGRID) THEN
       CALL NEWPEN(1)
       NXG = INT( (CDMAX0-CDMIN0)/DCD + 0.01 )
       NYG = INT( (CLMAX0-CLMIN0)/DCL + 0.01 )
       CALL PLGRID(CDWT*CDMIN0,CLWT*CLMIN0, 
     &             NXG,CDWT*DCD, NYG,CLWT*DCL, LMASK2 )
      ENDIF
C
C---- legend location
      XLEG    = CDWT*CDMIN0 + 10.0*CSL
      YLEG    = CLWT*CLMAX0 +  4.0*CSL
      XLIN(1) = CDWT*CDMIN0
      XLIN(2) = CDWT*CDMIN0 +  8.0*CSL
      YLIN(1) = CLWT*CLMAX0 +  4.4*CSL
      YLIN(2) = CLWT*CLMAX0 +  4.4*CSL
C
C---- CL label
      CALL NEWPEN(3)
      XPLT = -3.0*CSL
      YPLT = CLWT*(CLMAX0-1.5*DCL) - 0.3*CSL
      CALL PLCHAR(XPLT,YPLT,CSL,'c',0.0,1)
      CALL PLSUBS(XPLT,YPLT,CSL,'V',0.0,1,PLMATH)
C
C---- CD label
      CALL NEWPEN(2)
      XPLT = CDWT*(CDMAX0-1.5*DCD) - 0.7*CSL
      YPLT = -2.8*CSL
      CALL PLCHAR(XPLT,YPLT,CSL,'c',0.0,1)
      CALL PLSUBS(XPLT,YPLT,CSL,'d',0.0,1,PLCHAR)
C
C---- plot CL-CD polar curve 
      CALL NEWPEN(3)
      DALF = (ALMAX0-ALMIN0)/FLOAT(NTMP-1)
C
C-------------------------------------------------------------------
C
      IF(IPL.EQ.1) THEN
         NPLOT = NMACH
      ELSE
         NPLOT = NREYN
      ENDIF
C
      DO IPLOT=1, NPLOT
C
        IF(IPL.EQ.1) THEN
           TEMPM = AMACH(IPLOT)
           TEMPR = RCONST
        ELSE
           TEMPM = MCONST
           TEMPR = AREYN(IPLOT)
        ENDIF
C
C--- Limit MACH to 0.95 for sanity
        MACH = MIN(0.95,TEMPM)
        MSQ = MACH**2
        W = VSO*SQRT(MSQ)
        PGRT = 1.0 / SQRT(1.0 - MSQ)
C
        CALL NEWCOLOR(2+IPLOT)
        ILIN = IPLOT
C
        DO N=1, NTMP
C-------- set alpha
          ALFA = (ALMIN0 + DALF*FLOAT(N-1)) * PI/180.0
C-------- set corresponding CL and CD
          CLMX = CLMAX
          CLMN = CLMIN
          CALL CLCDCM(ALFA,W,TEMPR,VSO,SECSIG,SECSTAGR,
     &                CLIFT,CL_ALF,CL_W,STALLF,
     &                CDRAG,CD_ALF,CD_W,CD_REY,
     &                CMOM,CM_AL,CM_W,
     &                A0,CLMX,CLMN,DCLDA,DCLDA_STALL,DCL_STALL,
     &                CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP,
     &                TOC,DCDCL2S)
          YTMP(N) = CLIFT
          XTMP(N) = CDRAG
        END DO
        N = NTMP
        CALL CHKLIM(N,N1,N2,XTMP,1.1*CDMAX0)
        CALL XYLINE(N2-N1+1,XTMP(N1),YTMP(N1),0.0,CDWT,0.0,CLWT,ILIN)
C------ plot legend
        DELY = 2.5*CSL*FLOAT(IPLOT-1)
        CALL XYLINE(2,XLIN,YLIN,0.0,1.0,-DELY,1.0,ILIN)
C
        IF(IPL.EQ.1) THEN
          CALL PLCHAR(XLEG,YLEG+DELY,CSL,'M = ',0.0,4)
          CALL PLNUMB(999.,999.,CSL, MACH,0.0,2)
        ELSE
          CALL PLCHAR(XLEG,YLEG+DELY,CSL,'Re = ',0.0,5)
          CALL PLNUMB(999.,999.,CSL, TEMPR,0.0,-1)
        ENDIF
C
      END DO
C
      CALL NEWCOLOR(ICOL0)
C
C---- reset origin for alpha-CL plot
      CALL PLOT( CDWT*CDMAX0,0.0,-3)
      CALL PLOT( 0.10 , 0.0, -3)
      CALL PLOT(-ALWT*ALMIN0,0.0,-3)
C
C---- plot alpha-CL axes
      CALL NEWPEN(2)
      CALL YAXIS(0.0,CLWT*CLMIN0,-CLWT*(CLMAX0-CLMIN0),
     &               CLWT*DCL,CLMIN0,DCL,CS, 1)
      CALL XAXIS(ALWT*ALMIN0,0.0,-ALWT*(ALMAX0-ALMIN0),
     &               ALWT*DAL,ALMIN0,DAL,CS,-1)
      IF(LAGRID) THEN
       CALL NEWPEN(1)
       NXG = INT( (ALMAX0-ALMIN0)/DAL + 0.01 )
       NYG = INT( (CLMAX0-CLMIN0)/DCL + 0.01 )
       CALL PLGRID(ALWT*ALMIN0,CLWT*CLMIN0, 
     &              NXG,ALWT*DAL, NYG,CLWT*DCL, LMASK2 )
      ENDIF
C
C---- CL label
      CALL NEWPEN(3)
      XPLT = -3.0*CSL
      YPLT = CLWT*(CLMAX0-1.5*DCL) - 0.3*CSL
      CALL PLCHAR(XPLT,YPLT,CSL,'c',0.0,1)
      CALL PLSUBS(XPLT,YPLT,CSL,'V',0.0,1,PLMATH)
C
C---- alpha label
      CALL NEWPEN(3)
      XPLT = ALWT*(ALMAX0-0.5*DAL) - 0.8*CSL
      YPLT = -2.8*CSL
      CALL PLMATH(XPLT,YPLT,1.4*CS,'a"',0.0,2)
C
C---- plot alpha-CL curves
      CALL NEWPEN(3)
      DALF = (ALMAX0-ALMIN0)/FLOAT(NTMP-1)
C
C--------------------------------------------------------------------
C
      IF(IPL.EQ.1) THEN
         NPLOT = NMACH
      ELSE
C         NPLOT = NREYN  ! Turn on when Re influences CL
         NPLOT = 1       ! Only one CL/Alfa plot needed
      ENDIF
C
      DO IPLOT=1, NPLOT
C
        IF(IPL.EQ.1) THEN
           TEMPM = AMACH(IPLOT)
           TEMPR = RCONST
        ELSE
           TEMPM = MCONST
           TEMPR = AREYN(IPLOT)
        ENDIF
C
        MACH = MIN(0.95,TEMPM)
        MSQ = MACH**2
        W = VSO*SQRT(MSQ)
        PGRT = 1.0 / SQRT(1.0 - MSQ)
C
        CALL NEWCOLOR(2+IPLOT)
        ILIN = IPLOT
C
        DO N=1, NTMP
C-------- set alpha
          ALFA = (ALMIN0 + DALF*FLOAT(N-1)) * PI/180.0
C-------- set corresponding CL and CD
          CLMX = CLMAX
          CLMN = CLMIN
          CALL CLCDCM(ALFA,W,TEMPR,VSO,SECSIG,SECSTAGR,
     &                CLIFT,CL_ALF,CL_W,STALLF,
     &                CDRAG,CD_ALF,CD_W,CD_REY,
     &                CMOM,CM_AL,CM_W,
     &                A0,CLMX,CLMN,DCLDA,DCLDA_STALL,DCL_STALL,
     &                CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP,
     &                TOC,DCDCL2S)
          XTMP(N) = ALFA
          YTMP(N) = CLIFT
       END DO
        N = NTMP
        CALL CHKLIM(N,N1,N2,XTMP,1.1*ALMAX0*PI/180.)
        CALL XYLINE(N2-N1+1,XTMP(N1),YTMP(N1),0.0,ALWT*180.0/PI,0.0,
     &              CLWT,ILIN)
      END DO
      CALL NEWCOLOR(ICOL0)
C
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
C---- set factors and offsets for CL(CD) plot
      XYOFF(1) = 0.
      XYOFF(2) = 0.
      XYFAC(1) = CDWT
      XYFAC(2) = CLWT
C
      CALL PLOTABS(1.0,0.5,-3)
      CALL PLOT(0.0,-CLWT*CLMIN0,-3)
      CALL PLOT(5.0*CS,0.0,-3)
C
      RETURN
      END   ! AEROPLT2
C
C---------------------------------------------------------------------



      SUBROUTINE CHKLIM(N,NSTRT,NEND,F,FMAX)
C--- Get starting and end index for array values F(i) < FMAX
      DIMENSION F(N)
      NSTRT = 1
      NEND  = N
C--- Look for first point where F(i)<FMAX
      DO I=1,N
          IF(F(I).LT.FMAX) GO TO 10
      END DO
 10   NSTRT = MAX(I-1,1)
C--- Look for last point where F(i)<FMAX
      DO I=N,1,-1
          IF(F(I).LT.FMAX) GO TO 20
      END DO
 20   NEND = MIN(I+1,N)
C
      RETURN
      END



      SUBROUTINE OPFILE(LU,FNAME)
      CHARACTER*(*) FNAME
C
      CHARACTER*4 COMAND
      CHARACTER*128 COMARG,TMP
      CHARACTER*1 ANS, DUMMY
C
C---- get filename if it hasn't been already specified
      IF(FNAME.EQ.' ') CALL ASKS('Enter output filename^',FNAME)
C
C---- try to open file
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=50)
C
C---- file exists... ask how to proceed
      NF = INDEX(FNAME,' ') - 1
      TMP = 'File  '// FNAME(1:NF)//
     &      '  exists.  Overwrite / Append / New file ?^'
      CALL ASKC(TMP,COMAND,COMARG)
      ANS = COMAND(1:1)
C
C---- ask again if reply is invalid
      IF(INDEX('OoAaNn',ANS).EQ.0) THEN
        CALL ASKC(' O / A / N  ?^',COMAND,COMARG)
        ANS = COMAND(1:1)
C
        IF(INDEX('OoAaNn',ANS).EQ.0) THEN
C------- Still bad reply. Give up asking and just return
         WRITE(*,*) 'No action taken'
         RETURN
        ENDIF
      ENDIF
C
C---- at this point, file is open and reply is valid
      IF    (INDEX('Oo',ANS) .NE. 0) THEN
C------ go to beginning of file to overwrite
        REWIND(LU)
        GO TO 60
      ELSEIF(INDEX('Aa',ANS) .NE. 0) THEN
C------ go to end of file to append
        DO K=1, 12345678
          READ(LU,1000,END=40) DUMMY
 1000     FORMAT(A)
        ENDDO
 40     BACKSPACE(LU)
        GOTO 60
      ELSE
C------ new file... get filename from command argument, or ask if not supplied
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') CALL ASKS('Enter output filename^',FNAME)
      ENDIF
C
C---- at this point, file FNAME is new or is to be overwritten
 50   OPEN(LU,FILE=FNAME,STATUS='UNKNOWN',ERR=90)
      REWIND(LU)
C
 60   RETURN
C
 90   WRITE(*,*) 'Bad filename.'
      RETURN
      END ! OPFILE


      SUBROUTINE GETCLFACTOR(SIGMA,STAGGER,CLFACTOR)
C------------------------------------------------------------
C     Calculates multi-plane cascade effect on lift slope as a
C     function of solidity and stagger angle
C
C     Input:  SIGMA      solidity = Bc/(2*pi*r)
C             STAGGER    stagger angle (from axis to chordline, in rads)
C 
C     Output:  CLFACTOR  CLmultiplane/CL2D factor
C
C     Implements table-driven quadratic fit to Figure 6-29 in 
C     Wallis, Axial Flow Fans and Ducts.
C------------------------------------------------------------
      PARAMETER (PI=3.1415926535897932384)
      PARAMETER (DTR=PI/180.)
C
      DIMENSION X(11), A0(11), A1(11), A2(11)
C---- Table of quadratic fit coefficients
      DATA X      / 0.5, 
     &              0.6, 
     &              0.7, 
     &              0.8, 
     &              0.9, 
     &              1.0,
     &              1.1, 
     &              1.2, 
     &              1.3, 
     &              1.4, 
     &              1.5 /
      DATA A0     / 0.4755, 
     &              0.5255, 
     &              0.5722, 
     &              0.6142, 
     &              0.6647, 
     &              0.7016, 
     &              0.7643, 
     &              0.8302, 
     &              0.8932, 
     &              0.9366, 
     &              0.9814 /
      DATA A1     / -0.367495, 
     &              -0.341941, 
     &              -0.300058, 
     &              -0.255883, 
     &              -0.200593, 
     &              -0.114993, 
     &              -0.118602, 
     &              -0.130921, 
     &              -0.133442, 
     &              -0.077980, 
     &              -0.123071 /
      DATA A2     /  0.489466, 
     &               0.477648, 
     &               0.453027, 
     &               0.430048, 
     &               0.381462, 
     &               0.310028, 
     &               0.298309, 
     &               0.285309, 
     &               0.263084, 
     &               0.184165, 
     &               0.251594 /
C
      CLFACTOR = 1.0
      IF(SIGMA.LE.0.6) RETURN
C
C---- Interpolate quadratic fit coefficients by 1/solidity
      SIGI = 1.0/SIGMA
      CALL SEVLIN(SIGI,A0,X,11,AA0,DAA0)
      CALL SEVLIN(SIGI,A1,X,11,AA1,DAA1)
      CALL SEVLIN(SIGI,A2,X,11,AA2,DAA2)
C
C---- Only valid for stagger 20deg to 90deg, 
C     Limit low stagger to 20deg value to give constant lift ratio below that
      STAGR = STAGGER
      IF(STAGR.LT.20.0*DTR) STAGR = 20.0*DTR
      IF(STAGR.GT.90.0*DTR) STAGR = 90.0*DTR
C
C---- Quadratic fit for CLFACTOR at this SIGMA as function of STAGGER
      CLFACTOR = AA0 + AA1*STAGR + AA2*STAGR*STAGR
C---- maximum value of lift ratio should be limited to 1.0
      CLFACTOR = MIN(1.0,CLFACTOR)
C
      RETURN
      END


C
C--- New subroutines since 0.70 -----------------------------------
C
      SUBROUTINE GETAFNAME(AFNAME,NAME,NR,NSEC)
      CHARACTER*80 NAME,AFNAME
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
 10   AFNAME=NAME(1:L) //': Disk '// TNUM(NR)//': Sec '//TNUM(NSEC)
C
      RETURN
      END  ! GETAFNAME
C
C--------------------------------------------------------------------
C

      SUBROUTINE PLOTMACH(NR,NSEC,XISECT,A0,CLMAX,CLMIN,
     &                   DCLDA,DCLDA_STALL,DCL_STALL,CDMIN,CLDMIN,
     &                   DCDCL2,CMCON,MCRIT,REREF,REXP,DCDCL2S)
      INCLUDE 'DFDC.INC'
      CHARACTER*132 LINE
      LOGICAL RFIX, ERROR
      PARAMETER (NMACH = 10)
      PARAMETER (NREYN = 10)
      DIMENSION AMACH(NMACH),AREYN(NREYN)
      REAL MCONST
C
      COMMON /AFPLOTS/ AMACH,AREYN,NM,NRE,MCONST,RCONST,RFIX
C
 60   IPL=1
C
      WRITE(*,*)
C      WRITE(*,2005) RCONST
      IF(NM.GT.0) THEN
        WRITE(*,*) 'Current Mach set:'
        DO I = 1,NM
          WRITE(*,2010) AMACH(I)
        ENDDO
        WRITE(*,2020)
      ELSE
        WRITE(*,2030)
      ENDIF
C
      READ(*,1) LINE
C
      IF(LINE.EQ.'S'.OR. LINE.EQ.'s') GO TO 64
C
      IF(LINE.EQ.' ') THEN
        IF(NM.LE.0) THEN
           WRITE(*,*)'No Mach numbers specified'
           GO TO 60
        ELSE
           GO TO 68
        ENDIF
      ENDIF
C
      NM = NMACH
      CALL GETFLT(LINE,AMACH,NM,ERROR)          
      IF(.NOT.ERROR) GO TO 68
C
 64   CONTINUE
      AMACH(1) = 0.0
      AMACH(2) = 0.4
      AMACH(3) = 0.5
      AMACH(4) = 0.6
      AMACH(5) = 0.7
      AMACH(6) = 0.8
      AMACH(7) = 0.9
      NM = 7
C
 68   CALL AEROPLT2(NR,NSEC,XISECT,NM,AMACH,NRE,AREYN,
     &              A0,CLMAX,CLMIN,DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP,
     &              DCDCL2S, MCONST,RCONST,IPL)
C
 1    FORMAT(A)
C
 2005 FORMAT( ' Re constant: ',E9.3)
 2010 FORMAT(4X,F5.3)
 2020 FORMAT( ' <cr> to plot this set', 
     &       /' S for standard set',
     &      //' Enter Mach numbers')
 2030 FORMAT( ' Enter Mach numbers (S for standard set)')
C
      RETURN
      END    ! PLOTMACH
C
C-----------------------------------------------------------------------
C

      SUBROUTINE PLOTREYN (NR,NSEC,XISECT,A0,CLMAX,CLMIN,
     &                    DCLDA,DCLDA_STALL,DCL_STALL,CDMIN,CLDMIN,
     &                    DCDCL2,CMCON,MCRIT,REREF,REXP,DCDCL2S)
      INCLUDE 'DFDC.INC'
      CHARACTER*132 LINE
      LOGICAL RFIX, ERROR
      PARAMETER (NMACH = 10)
      PARAMETER (NREYN = 10)
      DIMENSION AMACH(NMACH),AREYN(NREYN)
      REAL MCONST
C
      COMMON /AFPLOTS/ AMACH,AREYN,NM,NRE,MCONST,RCONST,RFIX
C
 70   IPL=2
C
      WRITE(*,*)
C      WRITE(*,2005) MCONST
      IF(NRE.GT.0) THEN
        WRITE(*,*) 'Current Re set:'
        DO I = 1,NRE
          WRITE(*,2015) AREYN(I)
        ENDDO
        WRITE(*,2025)
      ELSE
        WRITE(*,2035)
      ENDIF
C
      READ(*,1) LINE
C
      IF(LINE.EQ.'S'.OR. LINE.EQ.'s') GO TO 74
C
      IF(LINE.EQ.' ') THEN
        IF(NRE.LE.0) THEN
           WRITE(*,*)'No Re numbers specified'
           GO TO 70
        ELSE
           GO TO 78
        ENDIF
      ENDIF
C
      NRE = NREYN
      CALL GETFLT(LINE,AREYN,NRE,ERROR)          
      IF(.NOT.ERROR) GO TO 78
C
 74   CONTINUE
      AREYN(1) = REREF * 0.5
      AREYN(2) = REREF * 0.75
      AREYN(3) = REREF
      AREYN(4) = REREF * 1.5
      AREYN(5) = REREF * 2.0
      NRE = 5
C
 78   CALL AEROPLT2(NR,NSEC,XISECT,NM,AMACH,NRE,AREYN,
     &              A0,CLMAX,CLMIN,DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP,
     &              DCDCL2S, MCONST,RCONST,IPL)
C
 1    FORMAT(A)
C
 2005 FORMAT( ' Mach constant: ',F5.3)
 2015 FORMAT(4X,E9.3)
 2025 FORMAT( ' <cr> to plot this set', 
     &       /' S for standard set',
     &      //' Enter Re numbers')
 2035 FORMAT( ' Enter Re numbers (S for standard set)')     
C
      RETURN
      END    ! PLOTREYN








 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
