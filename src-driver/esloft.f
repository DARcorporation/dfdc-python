C**********************************************************************
C    ESLOFT
C    Subroutine for rotor and stator lofting within DFDC.
C    Copyright (C) April 2009, 2013, Philip Carter, Esotec Developments.
C    philip (at) esotec (dot) org
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
C    Version 0.9, January 2014
C    - port back from ESLOFTX for CROTOR
C    - splined t and t/c
C    - persistent settings and airfoils across disks
C    - blade blockage geometry (velocity factor) calculation
C***********************************************************************
C
      SUBROUTINE ESLOFT
C--------------------------------------------------
C     Rotor and stator lofting in dfdc
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      CHARACTER*1 CHKEY, ANS, CIND
      CHARACTER*16 PROMPT
      CHARACTER*8 PROMPT2
      CHARACTER*4 COMAND,COMOLD
      CHARACTER*128 COMARG, ARGOLD, ARGPLT, FNAME, SAVFIL
      LOGICAL LOPMOD
C
C---- local arrays for calling GETINT,GETFLT...
      DIMENSION IINPUT(20),IPTEMP(NLSX)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
      CHARACTER*80 LINE,PLIN
C
C---- local arrays for repaneling and saving airfoils
      DIMENSION XXP(NPX,NAFX),YYP(NPX,NAFX),NNP(NAFX)
      DIMENSION IPFSAVE(NLSX)
C
C---- to prevent replots with changed airfoil numbers
      SAVE NAFOLD,NLOFTOLD
C
C---- "last" command (nothing yet)
      COMAND = '****'
      COMARG = ' '
      FNAME  = ' '
      ISPARS = ' '
C
      NDSK = 1
      N = NDSK
C
      COMOLD = COMAND
      ARGOLD = COMARG
C
C---- check if geometry loaded, process if available
C
      LLOFT = .TRUE.  !  GENGEOM will panel blades to geometric tip gap
      IF(LLOAD) THEN
         CALL GENGEOM
         IF(NAF.LT.2) THEN
           WRITE(*,*)'Load airfoils (LAF) or esloft set (LSET) to begin'
         ENDIF
      ELSE
         WRITE(*,*)
         WRITE(*,*)'No geometry loaded'
         WRITE(*,*)'You may work with airfoils only'
      ENDIF
C
C---- assume input data has changed
      LCALC  = .FALSE. !  Loft data is not calculated
      LBLEN  = .FALSE. !  Blended sections not calculated
      LTRAN  = .FALSE. !  Transformed sections not calculated
      LGEOPLX= .FALSE. !  No current xfoil plot
      ERROR  = .FALSE. !  Let things work smoothly...
      LTDEF  = .FALSE. !  No thickness data for spline routines
      LONAME =  NAME   !  Loft name begins as rotor name
C
      CALL GETLOFT(NDSK)
      IF(LLOAD.AND.NAF.GE.2) THEN
        CALL SHOWLOFT(NDSK)
      ENDIF
C
C-------------------------------------------------------------
C---- begin user-interaction loop
C............................................................
C
 500  CONTINUE
      LSMOD = .FALSE.  ! don't user-modify splined t or t/c
      PROMPT = '.LOFT^'
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
      IF(COMAND(1:1).EQ.'!') THEN
        IF(COMOLD.EQ.'****') THEN
          WRITE(*,*) 'Previous .LOFT command not valid'
          GO TO 500
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
        LPLOT = .FALSE.
        LLOFT = .FALSE.
        CALL PUTLOFT(NDSK)
        RETURN
      ENDIF
C
C---- can jump in here if COMARG string needs processing
 505  CONTINUE
C
C---- was only a new target-element index typed?
      CIND = COMAND(1:1)
      IF(INDEX('0123456789',CIND) .NE. 0) THEN
        READ(COMAND,*,ERR=508) ID
        IF(ID.LT.1 .OR. ID.GT.NROTOR) GO TO 508
C
        IF(IRTYPE(ID).EQ.1) THEN
           WRITE(*,1020) ID
           GO TO 500
        ELSEIF(IRTYPE(ID).EQ.2.AND.OMEGA(ID).NE.0.) THEN
           CALL PUTLOFT(NDSK)
           WRITE(*,1022) ID
        ELSEIF(IRTYPE(ID).EQ.2.AND.OMEGA(ID).EQ.0.) THEN
           CALL PUTLOFT(NDSK)
           WRITE(*,1024) ID
        ELSE
           WRITE(*,1026) ID
           GO TO 500
        ENDIF
C
        NDSK = ID
        N = NDSK
        LCALC = .FALSE.
        LBLEN = .FALSE.
        LTRAN = .FALSE.
C
        CALL GETLOFT(NDSK)
        IF(LLOAD.AND.NAF.GE.2) CALL SHOWLOFT(NDSK)
        GO TO 500
      ENDIF
 508  CONTINUE
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
 1010 FORMAT(1X,A4,' command not recognized.  Type a "?" for list')
 1020 FORMAT(/,1X, 'Disk',I2,' is an actuator disk')
 1022 FORMAT(/,1X, 'Disk',I2,' is a bladed rotor')
 1024 FORMAT(/,1X, 'Disk',I2,' is a bladed stator')
 1026 FORMAT(/,1X, 'Disk',I2,' is undefined')
C
C---- can jump in here if COMAND string has been set elsewhere
 510  CONTINUE
C
C------------------------------------------------------------------
C
      IF(    COMAND.EQ.'HELP'
     &  .OR. COMAND.EQ.'?   ') THEN
C
C------ Menu for commands
C
        IF(NROTOR.GT.1)THEN
           WRITE(*,1100)
           WRITE(*,1110)
        ELSE
           WRITE(*,*)
           WRITE(*,1110)
        ENDIF
C
 1100   FORMAT(
     &  /'    <i>     Select disk index')
 1110   FORMAT(
     &   '   SHOW <s> Show parent airfoil set',
     &  /'   LSET f   Load parent airfoil set from file',
     &  /'   SSET f   Save parent airfoil set to file',
     &  /'   NSET s   Name parent airfoil set',
     &  /'   DAF  i   Delete  airfoil from current set',
     &  /'   LAF  f   Load an airfoil from file',
     &  /'   SAF  i   Save an airfoil to file',
     &  /'   ALFZ i   Specify parent airfoil zero-lift alpha',
     &  /'   PANE ir  Specify parameters and repanel parents',
C
     & //'   DISP <d> Display current loft configuration',
     &  /'   THIC <t> Specify thickness distribution mode',
     &  /'   PARA r   Specify parabolic axis location',
     &  /'   THUB r   Specify blade thickness at hub',
     &  /'   TTIP r   Specify blade thickness at tip',
     &  /'   STN  ir  Specify station number and density (hub/tip)',
     &  /'   OSHO rr  Specify overshoot at hub and tip',     
     &  /'   PAX  rr  Specify pitch axes at hub and tip (x/c)',
     &  /'   TGAP r   Specify rotor geometric tip gap',
     &  /'   DIH  r   Specify dihedral at tip (circular arc)',
C
     & //'   DEF      Reset loft parameters to default settings',
     &  /'   BANG     Toggle local/global beta display (BGAM<0)',
     &  /'   NLOF s   Specify base name for loft output',
     &  /'   UNIT i   Specify units for loft output',
     &  /'   DIM      Toggle 2D/3D coordinate files',
     &  /'   ROTA     Toggle left/right hand rotation',
     &  /'   WRIT f   Save loft configuration data to file',
C
c     & //'   STOR     Store blade blockage data',
c     &  /'   PFAC     Specify blade blockage profile factor'
C
     & //'   SAVN ii  Save Normalized section(s) to file',
     &  /'   SAVT ii  Save Transformed section(s)   "  |',
     &  /'   SAVC     Save Centerbody coordinates   "  |',
     &  /'   SAVD     Save Duct coordinates         "  | output',
     &  /'   SAVR     Save loft station Radii       "  | units',
     &  /'   DISR     Display loft station Radii       |',
     &  /'   SAVB     Save ESBLADE loft data file       ',
C
     & //'   PLOP ii  Plot Parent section(s)',
     &  /'   PLON ii  Plot Normalized blended section(s)',
     &  /'   PLOT ii  Plot Transformed section(s)',
     &  /'  .DATA     Plot transformed section data vs. radius',
C
     & //'   BLOW <b> Blowup airfoil plot region',
     &  /'   RESE <r> Reset current airfoil plot scale and origin',
     &  /'   REPL     Replot current airfoil plot',
     &  /'   HARD     Hardcopy current plot',
     &  /'  .ANNO     Annotate plot',
     &  /'   SIZE r   Change absolute plot size',
     &  /'   Z        Zoom  ',
     &  /'   U        Unzoom',/)
C
C
C---------------------------------------------------------
C---- load airfoils from esloft file
C
      ELSEIF(COMAND.EQ.'LSET') THEN
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') THEN
          CALL ASKS(' Enter esloft airfoil filename^',FNAME)
        ENDIF
C
        IF(LLOAD) THEN
          DO NLD = 1,NROTOR
            IF(NLD.EQ.NDSK.OR.NAFQ(NLD).LT.2) THEN
              CALL AFLOAD(FNAME,ERROR,NLD)
            ENDIF
          ENDDO
        ELSE
          CALL AFLOAD(FNAME,ERROR,NDSK)
        ENDIF
C
        IF(.NOT.ERROR) THEN
          CALL SHOWAF(NDSK)
          LCALC = .FALSE.
          LBLEN = .FALSE.
          LTRAN = .FALSE.
        ENDIF
C
C---------------------------------------------------------
C---- save esloft airfoil file
C
      ELSEIF(COMAND.EQ.'SSET') THEN
        IF(NAF.LT.2) THEN
          WRITE(*,*)
          WRITE(*,*)'Minimum 2 airfoils required for esloft file'
          GO TO 500
        ENDIF
C
        CALL AFSAVE(COMARG,NDSK)
C
C--------------------------------------------------------------------
C---- load airfoil from generic file and process
C
      ELSEIF(COMAND.EQ.'LAF ') THEN
        IF(NAF.EQ.NAFX) THEN
          WRITE(*,1115) NAFX
          GO TO 500
        ENDIF
C
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') THEN
          CALL ASKS('Enter airfoil filename^',FNAME)
        ENDIF
C
        CALL DATLOAD(FNAME,ERROR,NDSK)
        IF(.NOT.ERROR) THEN
          CALL SHOWAF(NDSK)
          LCALC = .FALSE.
          LBLEN = .FALSE.
          LTRAN = .FALSE.
        ENDIF
C
 1115 FORMAT(/,' Maximum of ',I2,' parent airfoils',
     &       /,' Increase NAFX')
C
C--------------------------------------------------------------------
C---- save paneled airfoil from set
C
      ELSEIF(COMAND.EQ.'SAF ') THEN
        IF(NAF.LT.1) THEN
          WRITE(*,*)
          WRITE(*,*)'No airfoils loaded'
          GO TO 500
        ENDIF
C
        IF(NINPUT.EQ.1) THEN
          IDAT = IINPUT(1)
        ELSE
          WRITE(*,*)
          CALL ASKI('Enter airfoil index^',IDAT)
        ENDIF
C
        IF(IDAT.LT.1 .OR. IDAT.GT.NAF) THEN
          WRITE(*,*)
          WRITE(*,*)'Index out of range'
          GO TO 500
        ENDIF
C
        CALL DATSAVE(IDAT,NDSK)
C
C--------------------------------------------------------------------
C---- delete airfoil and process arrays
C
      ELSEIF(COMAND.EQ.'DAF ') THEN
        IF(NAF.LT.1) THEN
          WRITE(*,*)
          WRITE(*,*)'No airfoils loaded'
          GO TO 500
        ENDIF
C
        IF(NINPUT.EQ.1) THEN
          NDEL = IINPUT(1)
        ELSE
          WRITE(*,*)
          CALL ASKI('Enter airfoil index to delete^',NDEL)
        ENDIF
C
        IF(NDEL.LT.1 .OR. NDEL.GT.NAF) THEN
          WRITE(*,*)
          WRITE(*,*)'Index out of range'
          GO TO 500
        ENDIF
C
        DO I=NDEL,NAF-1
           DO K=1,NPP
             XXE(K,I)= XXE(K,I+1)
             YYE(K,I)= YYE(K,I+1)
           ENDDO
           PAREA (I,N)= PAREA (I+1,N)
           PLERAD(I,N)= PLERAD(I+1,N)
           PANGTE(I,N)= PANGTE(I+1,N)
           TOCE  (I,N)= TOCE  (I+1,N)
           TLOC  (I,N)= TLOC  (I+1,N)
           CAMLO (I,N)= CAMLO (I+1,N)
           CAMLOC(I,N)= CAMLOC(I+1,N)
           TETH  (I,N)= TETH  (I+1,N)
           ENAME (I,N)= ENAME (I+1,N)
           ALFZ  (I,N)= ALFZ  (I+1,N)
        ENDDO
C
        NAF=NAF-1
        LCALC = .FALSE.
        LBLEN = .FALSE.
        LTRAN = .FALSE.
C
        IF(NAF.GE.1) THEN
          CALL SHOWAF(NDSK)
        ELSE
          WRITE(*,*)
          WRITE(*,*)'No airfoils loaded'
        ENDIF
C
C--------------------------------------------------------------------
C---- specify A0
C
      ELSEIF(COMAND.EQ.'ALFZ') THEN
        IF(NAF.LT.1) THEN
          WRITE(*,*)
          WRITE(*,*)'No airfoils loaded'
          GO TO 500
        ENDIF
C
        IF(NINPUT.EQ.1) THEN
          NZER = IINPUT(1)
        ELSE
          WRITE(*,*)
          CALL ASKI('Enter airfoil index^',NZER)
        ENDIF
C
        IF(NZER.LT.1 .OR. NZER.GT.NAF) THEN
          WRITE(*,*)
          WRITE(*,*)'Index out of range'
          GO TO 500
        ENDIF
C
        WRITE(*,1120) NZER
 1120   FORMAT(/,' Airfoil ',I3)
        CALL ASKR('Enter zero-lift alpha (deg)^',ALFZ(NZER,N))
C
        LCALC = .FALSE.
        LBLEN = .FALSE.
        LTRAN = .FALSE.
        CALL SHOWAF(NDSK)     
C
C---------------------------------------------------------------
C---- display airfoil set
C
      ELSEIF(COMAND.EQ.'SHOW' .OR. COMAND.EQ.'S   ') THEN
        IF(NAF.LT.1) THEN
          WRITE(*,*)
          WRITE(*,*)'No airfoils loaded'
          GO TO 500
        ENDIF
C
        CALL SHOWAF(NDSK)
C
C---------------------------------------------------------------
C---- specify airfoil set name
C
      ELSEIF(COMAND.EQ.'NSET') THEN
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') THEN
          WRITE(*,1125) SNAME
          CALL ASKS(' Enter parent airfoil set name^',SNAME)
        ELSE
          SNAME = FNAME
        ENDIF
        CALL STRIP(SNAME,NSNAME)
C
 1125   FORMAT(/' Current airfoil set name is: ', A30)
C
C---------------------------------------------------------------
C---- change airfoil paneling
C
      ELSEIF(COMAND.EQ.'PANE') THEN
        IF(NAF.LT.1) THEN
          WRITE(*,*)
          WRITE(*,*)'No airfoils loaded'
          GO TO 500
        ENDIF
C
        NPP2T  = NPP2
        PAFT   = PAF
        NPMAX  = NPX/2
C
        IF(NINPUT.EQ.2)THEN
          NPP2T= IINPUT(1)
          PAFT = RINPUT(2)
        ELSEIF(NINPUT.EQ.1)THEN
          NPP2T= IINPUT(1)
        ELSE
          WRITE(*,*)
          CALL ASKI('Enter no. of points per side^',NPP2T)
          CALL ASKR('Enter panel density factor (le/te)^',PAFT)
        ENDIF
C
        IF(NPP2T.GT.NPMAX) THEN
          WRITE(*,1130) NPMAX
          GO TO 500
        ELSEIF(NPP2T.LT.20 .OR. PAFT.LT.1.0) THEN
          WRITE(*,*)'Parameters out of range'
          GO TO 500
        ELSEIF(NPP2T.EQ.NPP2.AND.PAFT.EQ.PAF) THEN
          WRITE(*,*)'Parameters unchanged'
          GO TO 500
        ENDIF
C
        DO I=1,NAF
          NNP(I)=NPP
        ENDDO
        NPP2= NPP2T
        NPP = NPP2*2
        PAF = PAFT
C
        DO I=1,NAF
          CALL AFPANEL(NPX,NAFX,I,NNP,XXE,YYE,NPP2,XXP,YYP,PAF)
          DO K=1,NPP
            XXE(K,I)= XXP(K,I)
            YYE(K,I)= YYP(K,I)
          ENDDO
        ENDDO
C
        LCALC = .FALSE.
        LBLEN = .FALSE.
        LTRAN = .FALSE.
        CALL AFDISP(6,NDSK)
C
 1130   FORMAT(/' Maximum ',I3,' points per side',
     &         /' Increase NPX')
C
C----------------------------------------------------------------
C---- change loft station parameters
C
      ELSEIF(COMAND.EQ.'STN ') THEN
        NLOFTT  = NLOFT
        PLOFTT  = PLOFT
        NLMAX   = NLSX-2
C
        IF(NINPUT.EQ.2) THEN
          NLOFTT = IINPUT(1)
          PLOFTT = RINPUT(2)
        ELSEIF(NINPUT.EQ.1) THEN
          NLOFTT = IINPUT(1)
        ELSE 
          WRITE(*,*)
          CALL ASKI('Enter no. of loft stations - hub->tip^',NLOFTT)
          CALL ASKR('Enter station density factor - hub/tip^',PLOFTT)
        ENDIF
C
        IF(NLOFTT.GT.NLMAX) THEN
           WRITE(*,1135) NLMAX
           GO TO 500
        ELSEIF(NLOFTT.LT.3 .OR. PLOFTT.LE.0.0) THEN
           WRITE(*,*)
           WRITE(*,*)'Parameters out of range'
           GO TO 500
        ELSEIF(NLOFTT.EQ.NLOFT.AND.PLOFTT.EQ.PLOFT) THEN
           WRITE(*,*)
           WRITE(*,*)'Parameters unchanged'
           GO TO 500
        ENDIF
C
        NLOFT = NLOFTT
        PLOFT = PLOFTT
C
        LCALC = .FALSE.
        LBLEN = .FALSE.
        LTRAN = .FALSE.
        CALL SHOWLOFT(NDSK)
C
 1135 FORMAT(/,' Maximum ',I2,' loft stations',
     &       /,' Increase NLSX')
C
C---------------------------------------------------------
C---- change tip gap - TGAPZL ignored for lofting
C
      ELSEIF(COMAND.EQ.'TGAP') THEN
        TGAPOLD = TGAP
        TGAPMM  = TGAP*1000.0
C        
        IF(NINPUT.EQ.1) THEN
          TGAPMM = RINPUT(1)
        ELSE
          WRITE(*,*)
          CALL ASKR('Enter rotor tip gap (mm)^',TGAPMM)
        ENDIF
        TGAP = TGAPMM/1000.0
C
        IF(TGAP.NE.TGAPOLD) THEN
           LCALC = .FALSE.
           LBLEN = .FALSE.
           LTRAN = .FALSE.
           CALL GENGEOM
           CALL SHOWLOFT(NDSK)
        ENDIF
C
C---------------------------------------------------------
C---- specify root thickness
C
      ELSEIF(COMAND.EQ.'THUB') THEN
        TKHUBT  = TKHUB
        TKHUBOLD= TKHUB
        IF(NINPUT.EQ.1) THEN
          TKHUBT = RINPUT(1)
        ELSE
          WRITE(*,*)
          WRITE(*,*)'0 to set thickest airfoil at hub'
          CALL ASKR('Enter blade thickness at hub (mm)^',TKHUBT)
        ENDIF
C
        IF(TKHUBT.EQ.TKHUBOLD) GO TO 500
C
        TKHUB = TKHUBT
        LCALC = .FALSE.
        LBLEN = .FALSE.
        LTRAN = .FALSE.
        CALL SHOWLOFT(NDSK)
C
C---------------------------------------------------------
C---- specify tip thickness
C
      ELSEIF(COMAND.EQ.'TTIP') THEN
        TKTIPT  = TKTIP
        TKTIPOLD= TKTIP
        IF(NINPUT.EQ.1) THEN
          TKTIPT = RINPUT(1)
        ELSE
          WRITE(*,*)
          WRITE(*,*)'0 to set thinnest airfoil at tip'
          CALL ASKR('Enter blade thickness at tip (mm)^',TKTIPT)
        ENDIF
C
        IF(TKTIPT.EQ.TKTIPOLD) GO TO 500
C
        TKTIP = TKTIPT
        LCALC = .FALSE.
        LBLEN = .FALSE.
        LTRAN = .FALSE.
        CALL SHOWLOFT(NDSK)
C
C---------------------------------------------------------
C---- specify thickness distribution
C
      ELSEIF(COMAND.EQ.'THIC' .OR. COMAND.EQ.'T   ') THEN
        ITTT    = ITTYPE
        ITTYPOLD= ITTYPE
        IF(NINPUT.EQ.1) THEN
         ITTT = IINPUT(1)
        ELSE
          WRITE(*,1140)
          CALL ASKI('Enter thickness distribution index^',ITTT)
        ENDIF
C
        IF(ITTT.LT.1 .OR. ITTT.GT.6) THEN
           WRITE(*,*)'Index out of range'
           GOTO 500
        ELSE
           ITTYPE=ITTT
           IF    (ITTYPE.EQ.1) THEN
              WRITE(*,*)
              WRITE(*,*)'Linear t/c selected'
           ELSEIF(ITTYPE.EQ.2) THEN
              WRITE(*,*)
              WRITE(*,*)'Parabolic t/c selected'
           ELSEIF(ITTYPE.EQ.3) THEN
              WRITE(*,*)
              WRITE(*,*)'Splined t/c selected'
              LSMOD = .TRUE.
           ELSEIF(ITTYPE.EQ.4) THEN
              WRITE(*,*)
              WRITE(*,*)'Linear thickness selected'
           ELSEIF(ITTYPE.EQ.5) THEN
              WRITE(*,*)
              WRITE(*,*)'Parabolic thickness selected'
           ELSEIF(ITTYPE.EQ.6) THEN
              WRITE(*,*)
              WRITE(*,*)'Splined thickness selected'
              LSMOD = .TRUE.
           ENDIF  
        ENDIF
C
        LCALC = .FALSE.
        LBLEN = .FALSE.
        LTRAN = .FALSE.
C
        IF(ITTYPE.EQ.3 .OR. ITTYPE.EQ.6) THEN
          CALL PRINTLOFT(6,NDSK)
        ELSE
          CALL SHOWLOFT(NDSK)
        ENDIF
C
 1140   FORMAT(
     &  /, ' 1  Linear    t/c',
     &  /, ' 2  Parabolic t/c',
     &  /, ' 3  Splined   t/c',
     &  /, ' 4  Linear    t ',
     &  /, ' 5  Parabolic t ',
     &  /, ' 6  Splined   t ',/)
C
C---------------------------------------------------------
C---- specify parabolic axis
C
      ELSEIF(COMAND.EQ.'PARA') THEN
        PAXIST =  PAXIS
        PAXISOLD= PAXIS
        IF(NINPUT.EQ.1) THEN
          PAXIST = RINPUT(1)
        ELSE
          WRITE(*,*)
          WRITE(*,*)'Units of bladelengths inboard of hub'
          CALL ASKR('Locate parabolic axis^',PAXIST)
        ENDIF
C
        IF(.NOT.LCALC) THEN
           IF(PAXIST.LE.0.0) THEN
             WRITE(*,*)
             WRITE(*,*)'PARA must be set greater than zero'
           ELSE
             PAXIS = PAXIST
           ENDIF
        ELSE
          YPARAT = PAXIST*(YLOFT(ITIP)-YLOFT(IHUB))
          PLOC   = YLOFT(IHUB)-YPARAT
          IF(PLOC.GE.YLOFT(1)) THEN
            WRITE(*,*)
            WRITE(*,*)'Axis must be inboard of Station 1'
            YPARAT= PAXIS*(YLOFT(ITIP)-YLOFT(IHUB))
          ELSE
            PAXIS = PAXIST
          ENDIF
        ENDIF
C
        WRITE(*,1145) PAXIS
        IF(LCALC) THEN
          YPARACM = YPARAT*100.0
          WRITE(*,1147) YPARACM
        ENDIF
C
        IF(PAXIS.NE.PAXISOLD) THEN
          LCALC = .FALSE.
          LBLEN = .FALSE.
          LTRAN = .FALSE.
          CALL SHOWLOFT(NDSK)
        ENDIF
C
 1145   FORMAT(/,1X,'Parabolic axis at',F6.3,' bladelengths')
 1147   FORMAT(1X,F7.3,' cm inboard of hub')
C
C---------------------------------------------------------
C---- specify overshoot 
C
      ELSEIF(COMAND.EQ.'OSHO') THEN
        OSHUBT = OSHUB
        OSTIPT = OSTIP
        OSHUBOLD= OSHUB
        OSTIPOLD= OSTIP
        IF(NINPUT.EQ.2) THEN
          OSHUBT = RINPUT(1)
          OSTIPT = RINPUT(2)
        ELSEIF(NINPUT.EQ.1) THEN
          OSHUBT = RINPUT(1)
          WRITE(*,*)
          CALL ASKR('Enter overshoot at tip (mm)^',OSTIPT)
        ELSE
          WRITE(*,*)
          CALL ASKR('Enter overshoot at hub (mm)^',OSHUBT)
          CALL ASKR('Enter overshoot at tip (mm)^',OSTIPT)
        ENDIF
C
        IF(OSHUBT.LT.0.0 .OR. OSTIPT.LT.0.0) THEN
          WRITE(*,*)
          WRITE(*,*)'Negative overshoot not allowed'
          GO TO 500
        ELSEIF(OSHUBT.EQ.OSHUBOLD.AND.OSTIPT.EQ.OSTIPOLD) THEN
          GO TO 500
        ENDIF
C
        OSHUB = OSHUBT
        OSTIP = OSTIPT
        LCALC = .FALSE.
        LBLEN = .FALSE.
        LTRAN = .FALSE.
        CALL SHOWLOFT(NDSK)
C
C---------------------------------------------------------
C---- change airfoil pitch axes x/c
C
      ELSEIF(COMAND.EQ.'PAX') THEN
        AXHUBOLD= AXHUB
        AXTIPOLD= AXTIP
        IF(NINPUT.EQ.2) THEN
          AXHUB = RINPUT(1)
          AXTIP = RINPUT(2)
        ELSEIF(NINPUT.EQ.1) THEN
          AXHUB = RINPUT(1)
          WRITE(*,*)
          CALL ASKR('Enter blade pitch axis at tip (x/c)^',AXTIP)
        ELSE
          WRITE(*,*)
          CALL ASKR('Enter blade pitch axis at hub (x/c)^',AXHUB)
          CALL ASKR('Enter blade pitch axis at tip (x/c)^',AXTIP)
        ENDIF
C
        IF(AXTIP.NE.AXTIPOLD.OR.AXHUB.NE.AXHUBOLD) THEN
          LCALC = .FALSE.
          LBLEN = .FALSE.
          LTRAN = .FALSE.
          CALL SHOWLOFT(NDSK)
        ENDIF
C
C---------------------------------------------------------
C---- change airfoil axes dihedral
C
      ELSEIF(COMAND.EQ.'DIH') THEN
        AZTIPOLD= AZTIP
        IF(NINPUT.EQ.1) THEN
          AZTIP = RINPUT(1)
        ELSE
          WRITE(*,*)
          CALL ASKR('Enter dihedral at tip (mm)^',AZTIP)
        ENDIF
C
        IF(AZTIP.NE.AZTIPOLD) THEN
          LTRAN = .FALSE.
          CALL SHOWLOFT(NDSK)
        ENDIF
C
C---------------------------------------------------------
C---- display current loft data
C
      ELSEIF(COMAND.EQ.'DISP'.OR.COMAND.EQ.'D   ') THEN
        CALL SHOWLOFT(NDSK)
C
C---------------------------------------------------------
C---- write loft data to disk
C
      ELSEIF(COMAND.EQ.'WRIT') THEN
        IF(.NOT.LCALC) THEN
           CALL LOFTGEOM(NDSK,ERROR)
           IF(ERROR) GO TO 500
        ENDIF
C
c        IF(COMARG(1:1).NE.' ') SAVFIL = COMARG
c        CALL OPFILE(LUSAVE,SAVFIL)
c        CALL PRINTLOFT(LUSAVE,NDSK)
c        CALL AFDISP(LUSAVE,NDSK)
c        CLOSE(LUSAVE)
C
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') THEN
          CALL ASKS(' Enter loft data filename^',FNAME)
        ENDIF
        CALL STRIP(FNAME,NF)
C
        LU = 19
        OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=605)
C
        WRITE(*,*)
        WRITE(*,*) 'Output file exists. Overwrite? (Y/n)'
        READ (*,1000) ANS
        IF(INDEX('Nn',ANS).EQ.0) GO TO 610
C
        CLOSE(LU)
        WRITE(*,*)
        WRITE(*,*) 'Output file not saved'
        GO TO 615
C
 605    OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=615)
 610    REWIND(LU)
C
        CALL PRINTLOFT(LU,NDSK)
        CALL AFDISP(LU,NDSK)
        CLOSE(LU)
C
 615    CONTINUE
C
C---------------------------------------------------------
C---- reset loft parameters to default settings
C
      ELSEIF(COMAND.EQ.'DEF ') THEN
         CALL LOFTINIT2(NDSK,NDSK)
         LONAME = NAME
         LCALC = .FALSE.
         LBLEN = .FALSE.
         LTRAN = .FALSE.
         CALL SHOWLOFT(NDSK)
C
         WRITE(*,*)
         WRITE(*,*) 'Loft parameters reset to defaults'
         WRITE(*,*) 'Airfoil set and output settings unchanged'
C
C---------------------------------------------------------
C---- toggle local/global beta coordinates
C
      ELSEIF(COMAND.EQ.'BANG') THEN
         LCOORD = .NOT.LCOORD
         WRITE(*,*)
         IF(LCOORD) THEN
            WRITE(*,*)'Local beta coordinates selected'
         ELSE
            WRITE(*,*)'Global beta coordinates selected'
         ENDIF
C
         CALL PRINTLOFT(6,NDSK)
C
C---------------------------------------------------------
C---- toggle between 2D and 3D points file output
C
      ELSEIF(COMAND.EQ.'DIM ') THEN
        LL2D = .NOT.LL2D
        WRITE(*,*)
        IF(LL2D) THEN
           WRITE(*,*)'2D points files will be saved'
        ELSE
           WRITE(*,*)'3D points files will be saved'
        ENDIF
C
C---------------------------------------------------------
C---- toggle between left and right hand rotation output
C
      ELSEIF(COMAND.EQ.'ROTA') THEN
        LROTATE = .NOT.LROTATE
        WRITE(*,*)
        IF(LROTATE) THEN
          WRITE(*,*)'Left hand rotation selected (when BGAM>0)'
        ELSE
         WRITE(*,*)'Right hand rotation selected (when BGAM>0)'
        ENDIF
        LTRAN= .FALSE.
C
C---------------------------------------------------------
C---- change loft output units
C
      ELSEIF(COMAND.EQ.'UNIT') THEN
        IF(NINPUT.EQ.1) THEN
           IUNTYPE = IINPUT(1)
        ELSE              
           WRITE(*,1150)
           IF(OUTFAC.EQ.1.) THEN
              IUNTYPE= 1
              WRITE(*,*)'Currently meters'
           ELSEIF(OUTFAC.EQ.100.) THEN
              IUNTYPE= 2
              WRITE(*,*)'Currently centimeters'
           ELSEIF(OUTFAC.EQ.1000.) THEN
              IUNTYPE= 3        
              WRITE(*,*)'Currently millimeters'
           ELSEIF(OUTFAC.EQ.MTI) THEN
              IUNTYPE= 4     
              WRITE(*,*)'Currently inches'
           ENDIF
           WRITE(*,*)
           CALL ASKI('Specify output units^',IUNTYPE)
        ENDIF
C
        WRITE(*,*)
        IF(IUNTYPE.EQ.1) THEN
           OUTFAC = 1.0
           WRITE(*,*) 'Output units set to meters'
        ELSEIF(IUNTYPE.EQ.2) THEN
           OUTFAC = 100.0
           WRITE(*,*) 'Output units set to centimeters'
        ELSEIF(IUNTYPE.EQ.3) THEN
           OUTFAC = 1000.0
           WRITE(*,*) 'Output units set to millimeters'
        ELSEIF(IUNTYPE.EQ.4) THEN
           OUTFAC = MTI
           WRITE(*,*) 'Output units set to inches'
        ELSE
           WRITE(*,*) 'Wrong index'
           WRITE(*,*) 'Output units set to meters'
           IUNTYPE= 1
           OUTFAC = 1.0
        ENDIF
C
 1150   FORMAT(
     &  /, '  1  meters',
     &  /, '  2  centimeters',
     &  /, '  3  millimeters',
     &  /, '  4  inches')
C
C---------------------------------------------------------
C---- change loft output base name
C
      ELSEIF(COMAND.EQ.'NLOF') THEN
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') THEN
          WRITE(*,*)
          WRITE(*,1160) LONAME
          CALL ASKS(' Enter new base name^',LONAME)
        ELSE
          LONAME = FNAME
        ENDIF
        CALL STRIP(LONAME,NLONAME)
C
 1160   FORMAT(1X,'Current loft output base name: ', A30)
C
C
C---------------------------------------------------------
C---- PLOTTING
C---- plot individual parent airfoils
C
      ELSEIF(COMAND.EQ.'PLOP') THEN
        IF(NAF.LT.1) THEN
          WRITE(*,*)
          WRITE(*,*)'No airfoils loaded'
          GO TO 500
        ENDIF
C
        DO I=1,NLSX
           IPTEMP(I)= 0
        ENDDO
C
        IF(NINPUT.GE.1) THEN
           DO I=1,NINPUT
              IPTEMP(I)= IINPUT(I)
           ENDDO
           NXF= NINPUT
        ELSE
           NXF = NAF
           CALL ASKS(' Enter airfoil indices - 0 for all^',PLIN)
           CALL GETINT(PLIN,IPTEMP,NXF,ERROR)
           IF(NXF.LT.1) GO TO 500
        ENDIF
C
        IF(IPTEMP(1).EQ.0) THEN
           DO I=1,NAF
             IXPLOT(I)=I
           ENDDO
           NXPLOT = NAF
        ELSE
           DO I=1,NXF
              IF(IPTEMP(I).LT.1 .OR. IPTEMP(I).GT.NAF) THEN
                 WRITE(*,*)
                 WRITE(*,*)'Airfoil index out of range'
                 GO TO 500
              ENDIF
           ENDDO
C
           NXPLOT = NXF
           DO I=1,NXPLOT
              IXPLOT(I)=IPTEMP(I)
           ENDDO
        ENDIF
C
        NAFOLD = NAF
        IXTYPE = 1
        CALL XFOILPLOT(NDSK, 1)  ! no blowup
C
C----------------------------------------------------------------
C---- plot blended or transformed sections
C
      ELSEIF(COMAND.EQ.'PLON' .OR. COMAND.EQ.'PLOT') THEN
C
        IF(COMAND.EQ.'PLON') THEN
          IF(.NOT.LBLEN) CALL AFBLEND(NDSK,ERROR)
          IXTYPE = 2
        ELSE
          IF(.NOT.LTRAN) CALL AFTRANS(NDSK,ERROR)
          IXTYPE = 3
        ENDIF
C
        IF(ERROR) THEN
          LGEOPLX= .FALSE.
          WRITE(*,*)
          WRITE(*,*) 'LCALC error'
          GO TO 500
        ENDIF
C
        DO I=1,NLSX
           IPTEMP(I)= 0
        ENDDO
C
        IF(NINPUT.GE.1) THEN
           DO I=1,NINPUT
              IPTEMP(I)= IINPUT(I)
           ENDDO
           NXF= NINPUT
        ELSE
           NXF = NLOFT2
           CALL ASKS(' Enter airfoil indices - 0 for all^',PLIN)
           CALL GETINT(PLIN,IPTEMP,NXF,ERROR)
           IF(NXF.LT.1) GO TO 500
        ENDIF
C
        IF(IPTEMP(1).EQ.0) THEN
           DO I=1,NLOFT2
             IXPLOT(I)=I
           ENDDO
           NXPLOT = NLOFT2
        ELSE
           DO I=1,NXF
              IF(IPTEMP(I).LT.1 .OR. IPTEMP(I).GT.NLOFT2) THEN
                 WRITE(*,*)
                 WRITE(*,*)'Airfoil index out of range'
                 GO TO 500
              ENDIF
           ENDDO
C
           NXPLOT = NXF
           DO I=1,NXPLOT
              IXPLOT(I)=IPTEMP(I)
           ENDDO
        ENDIF
C
        NLOFTOLD= NLOFT2
        CALL XFOILPLOT(NDSK, 1)  ! no blowup
C
C--------------------------------------------------------
C---- blow up current plot
C
      ELSEIF(COMAND.EQ.'BLOW' .OR.
     &       COMAND.EQ.'B   '      ) THEN
C
        IF(.NOT.LGEOPLX) THEN
          WRITE(*,*)
          WRITE(*,*)'No current xfoil plot to blowup'
          GO TO 500
        ENDIF
C
        IF(IXTYPE.EQ.1 .AND. NAFOLD.NE.NAF) THEN
          WRITE(*,*)
          WRITE(*,*)'Parent airfoil numbers have changed'
          WRITE(*,*)'Use PLOP to set airfoil indices'
          GO TO 500
C
        ELSEIF(IXTYPE.EQ.2) THEN
          IF(.NOT.LBLEN) THEN
            CALL AFBLEND(NDSK,ERROR)
            IF(ERROR) GO TO 500
          ENDIF
C
          IF(NLOFTOLD.NE.NLOFT2) THEN
            WRITE(*,*)
            WRITE(*,*)'Station numbers have changed'
            WRITE(*,*)'Use PLON to set station indices'
            GO TO 500
          ENDIF
C
        ELSEIF(IXTYPE.EQ.3) THEN
          IF(.NOT.LTRAN) THEN
            CALL AFTRANS(NDSK,ERROR)
            IF(ERROR) GO TO 500
          ENDIF
C
          IF(NLOFTOLD.NE.NLOFT2) THEN
            WRITE(*,*)
            WRITE(*,*)'Station numbers have changed'
            WRITE(*,*)'Use PLOT to set station indices'
            GO TO 500
          ENDIF
        ENDIF
C
        CALL OFFGETX
        CALL XFOILPLOT(NDSK, 2)  ! do blowup
C
C--------------------------------------------------------
C---- reset current airfoil plot scale and offset
C
      ELSEIF(COMAND.EQ.'RESE' .OR.
     &       COMAND.EQ.'R   '    ) THEN
C
        IF(.NOT.LGEOPLX) THEN
          WRITE(*,*)
          WRITE(*,*)'No current xfoil plot to reset'
          GO TO 500
        ENDIF
C
        IF(IXTYPE.EQ.1 .AND. NAFOLD.NE.NAF) THEN
          WRITE(*,*)
          WRITE(*,*)'Parent airfoil numbers have changed'
          WRITE(*,*)'Use PLOP to set airfoil indices'
          GO TO 500
C
        ELSEIF(IXTYPE.EQ.2) THEN
          IF(.NOT.LBLEN) THEN
            CALL AFBLEND(NDSK,ERROR)
            IF(ERROR) GO TO 500
          ENDIF
C
          IF(NLOFTOLD.NE.NLOFT2) THEN
            WRITE(*,*)
            WRITE(*,*)'Station numbers have changed'
            WRITE(*,*)'Use PLON to set station indices'
            GO TO 500
          ENDIF
C
        ELSEIF(IXTYPE.EQ.3) THEN
          IF(.NOT.LTRAN) THEN
            CALL AFTRANS(NDSK,ERROR)
            IF(ERROR) GO TO 500
          ENDIF
C
          IF(NLOFTOLD.NE.NLOFT2) THEN
            WRITE(*,*)
            WRITE(*,*)'Station numbers have changed'
            WRITE(*,*)'Use PLOT to set station indices'
            GO TO 500
          ENDIF
        ENDIF
C
        CALL XFOILPLOT(NDSK, 1)  ! no blowup
C
C--------------------------------------------------------
C---- replot current airfoil plot
C
      ELSEIF(COMAND.EQ.'REPL' ) THEN
C
        IF(.NOT.LGEOPLX) THEN
          WRITE(*,*)
          WRITE(*,*)'No current xfoil plot to replot'
          GO TO 500
        ENDIF
C
        IF(IXTYPE.EQ.1 .AND. NAFOLD.NE.NAF) THEN
          WRITE(*,*)
          WRITE(*,*)'Parent airfoil numbers have changed'
          WRITE(*,*)'Use PLOP to set airfoil indices'
          GO TO 500
C
        ELSEIF(IXTYPE.EQ.2) THEN
          IF(.NOT.LBLEN) THEN
            CALL AFBLEND(NDSK,ERROR)
            IF(ERROR) GO TO 500
          ENDIF
C
          IF(NLOFTOLD.NE.NLOFT2) THEN
            WRITE(*,*)
            WRITE(*,*)'Station numbers have changed'
            WRITE(*,*)'Use PLON to set station indices'
            GO TO 500
          ENDIF
C
        ELSEIF(IXTYPE.EQ.3) THEN
          IF(.NOT.LTRAN) THEN
            CALL AFTRANS(NDSK,ERROR)
            IF(ERROR) GO TO 500
          ENDIF
C
          IF(NLOFTOLD.NE.NLOFT2) THEN
            WRITE(*,*)
            WRITE(*,*)'Station numbers have changed'
            WRITE(*,*)'Use PLOT to set station indices'
            GO TO 500
          ENDIF
        ENDIF
C
        CALL XFOILPLOT(NDSK, 3)  ! use current scale and origin
C
C---------------------------------------------------------------
C---- plot section data
C
      ELSEIF(COMAND.EQ.'DATA') THEN
        IF(.NOT.LBLEN) THEN
           CALL AFBLEND(NDSK,ERROR)
           IF(ERROR) GO TO 500
        ENDIF
C
        CALL PLOTLOFDATA(NDSK)
        LGEOPLX = .FALSE.
C
C--------------------------------------------------------------
C---- save blended and transformed sections to disk
C
      ELSEIF(COMAND.EQ.'SAVN' .OR. COMAND.EQ.'SAVT') THEN
C
        IF(COMAND.EQ.'SAVN') THEN
          IF(.NOT.LBLEN) CALL AFBLEND(NDSK,ERROR)
          IF(ERROR) GO TO 500
          IPFTYPE = 1
        ELSE
          IF(.NOT.LTRAN) CALL AFTRANS(NDSK,ERROR)
          IF(ERROR) GO TO 500
          IPFTYPE = 2
        ENDIF
C
        DO I=1,NLSX
           IPTEMP(I)= 0
        ENDDO
C
        IF(NINPUT.GE.1) THEN
           DO I=1,NINPUT
              IPTEMP(I)= IINPUT(I)
           ENDDO
           NPF= NINPUT
        ELSE
           NPF = NLOFT2
           CALL ASKS(' Enter station indices - 0 for all^',PLIN)
           CALL GETINT(PLIN,IPTEMP,NPF,ERROR)
           IF(NPF.LT.1) GO TO 500
        ENDIF
C
        IF(IPTEMP(1).EQ.0) THEN
           DO I=1,NLOFT2
             IPFSAVE(I)=I
           ENDDO
           NPFSAVE = NLOFT2
        ELSE
           DO I=1,NPF
              IF(IPTEMP(I).LT.1 .OR. IPTEMP(I).GT.NLOFT2) THEN
                 WRITE(*,*)
                 WRITE(*,*)'Station index out of range'
                 GO TO 500
              ENDIF
           ENDDO
C
           NPFSAVE = NPF
           DO I=1,NPFSAVE
              IPFSAVE(I)= IPTEMP(I)
           ENDDO
        ENDIF
C
        CALL SAVLOFT(NDSK,NPFSAVE,IPFSAVE,IPFTYPE)
C
C--------------------------------------------------------------
C---- save centerbody coords to disk
C
      ELSEIF(COMAND.EQ.'SAVC') THEN
        IF(NPTOT.EQ.0) THEN
          WRITE(*,*)
          WRITE(*,*) 'No geometry loaded'
        ELSE
          CALL SAVGEOM(1)
        ENDIF
C
C--------------------------------------------------------------
C---- save duct coords to disk
C
      ELSEIF(COMAND.EQ.'SAVD') THEN
        IF(NPTOT.EQ.0) THEN
          WRITE(*,*)
          WRITE(*,*) 'No geometry loaded'
        ELSE
          CALL SAVGEOM(2)
        ENDIF
C
C--------------------------------------------------------------
C---- save station radii to disk
C
      ELSEIF(COMAND.EQ.'SAVR') THEN
        IF(.NOT.LCALC) THEN
          CALL LOFTGEOM(NDSK,ERROR)
          IF(ERROR) GO TO 500
        ENDIF
        LU = 19
        CALL SAVRAD(LU,NDSK,.FALSE.,.FALSE.)
C
C--------------------------------------------------------------
C---- display station radii
C
      ELSEIF(COMAND.EQ.'DISR') THEN
        IF(.NOT.LCALC) THEN
          CALL LOFTGEOM(NDSK,ERROR)
          IF(ERROR) GO TO 500
        ENDIF
        CALL SAVRAD(6,NDSK,.FALSE.,.FALSE.)
C
C--------------------------------------------------------------
C---- save ESBLADE file
C
      ELSEIF(COMAND.EQ.'SAVB') THEN
        IF(.NOT.LTRAN) THEN
          CALL AFTRANS(NDSK,ERROR)
          IF(ERROR) GO TO 500
        ENDIF
        CALL SAVEBLADE(NDSK)
C
C---------------------------------------------------------------
C---- store blade blockage data
C
c      ELSEIF(COMAND.EQ.'STOR') THEN
c        CALL GETBBFAC(NDSK)
c        WRITE(*,2100) N, BPLAST
c        WRITE(*,2150)
c        WRITE(*,2200)
c        DO I=1,NRC
c          WRITE(*,2300) I, YRC(I,N)*100.0, BBVFAC(I,N)
c        ENDDO
C
 2100 FORMAT(/,' Blade blockage data stored for Disk',I2,
     &       /,' Profile factor =',F5.2)
 2150 FORMAT(1X,37('-'))
 2200 FORMAT(12X,'r_cm',5X,'Vel_factor')
 2300 FORMAT(3X,I2,5X,F7.3,5X,F8.6)
C
C--------------------------------------------------------------------
C---- change blade blockage profile factor
C
      ELSEIF(COMAND.EQ.'PFAC') THEN
        IF(NINPUT.EQ.2) THEN
         BBPFAC = RINPUT(2)
        ELSEIF(NINPUT.EQ.1) THEN
         BBPFAC = RINPUT(1)
        ELSE
         WRITE(*,*)
      CALL ASKR('Enter blade blockage profile factor (1<->1.5)^',BBPFAC)
        ENDIF
C
C---------------------------------------------------------------
C---- standard plot controls
C---- annotate plot
C
      ELSEIF(COMAND.EQ.'ANNO') THEN
       IF(LPLOT) THEN
        CALL ANNOT(CHGT)
       ELSE
        WRITE(*,*) 'No active plot to annotate'
       ENDIF
C
C--------------------------------------------------------------------
C---- hard copy current plot
C
      ELSEIF(COMAND.EQ.'HARD') THEN
        IF(LPLOT) CALL PLEND
        LPLOT = .FALSE.
        CALL REPLOT(IDEVPS)
C
C--------------------------------------------------------------------
C---- change plot size
C
      ELSEIF(COMAND.EQ.'SIZE') THEN
        IF(NINPUT.GE.1) THEN
         SIZE = RINPUT(1)
        ELSE
         CALL ASKR('Enter plot size^',SIZE)
        ENDIF
C
C--------------------------------------------------------------------
C---- zoom plot
C
      ELSEIF(COMAND.EQ.'Z   ') THEN
        IF(LPLOT) THEN
         CALL USETZOOM(.TRUE.,.TRUE.)
         CALL REPLOT(IDEV)
        ENDIF
C
C--------------------------------------------------------------------
C---- unzoom plot
C
      ELSEIF(COMAND.EQ.'U   ') THEN
        IF(LPLOT) THEN
         CALL CLRZOOM
         CALL REPLOT(IDEV)
        ENDIF
C
C--------------------------------------------------------------------
C
C---- that's it!
C
      ELSE
        WRITE(*,1010) COMAND
      ENDIF
C
      GO TO 500
C
      END !  ESLOFT
C
C------------------------------------------------------------------




C------------------------------------------------------------------
C---- ESLOFT SUBROUTINES
C------------------------------------------------------------------


      SUBROUTINE SHOWLOFT(NDSK)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      LOGICAL ERROR
C
      PLFACD= 0.67  !  This controls size of LOFTPLT vs window
      XORG  = 0.12  !  These control the origin
      YORG  = 0.11
C
      IF(.NOT.LCALC) THEN
         CALL LOFTGEOM(NDSK,ERROR)
         IF(ERROR) RETURN
      ENDIF
C
      CALL PRINTLOFT(6,NDSK)
      CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFACD*SIZE,LPLOT,.FALSE.)
      CALL PLOT(XORG,YORG,-3)
      CALL LOFTPLT(NDSK)
C
      LGEOPLX = .FALSE.
C
      RETURN
      END  !  SHOWLOFT
C-------------------------------------------------------------------


      SUBROUTINE SHOWAF(NDSK)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      CALL AFDISP(6,NDSK)
      CALL PLTINI(SCRNFR,IPSLU,IDEV,SIZE,LPLOT,.FALSE.)
      CALL AFPLOT(NDSK)
C
      LGEOPLX = .FALSE.
C
      RETURN
      END  !  SHOWAF
C
C-------------------------------------------------------------------



      SUBROUTINE AFDISP(LU,NDSK)
C-------------------------------------------------------------------
C     Prints airfoil set data
C-------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      N = NDSK
C
      WRITE(LU,1050) N,NPP2,PAF
      WRITE(LU,1000)
      WRITE(LU,1100)SNAME
      DO I=1,NAF
        WRITE(LU,1120) I,ENAME(I,N),TOCE(I,N),TLOC(I,N),
     &  CAMLO(I,N),CAMLOC(I,N),PLERAD(I,N),TETH(I,N),ALFZ(I,N)
      ENDDO
      WRITE(LU,1000)
C
 1000 FORMAT(1X,78('-'))
 1050 FORMAT(/,1X,'Parent Airfoils: Disk',I2, 17X,
     &          'Points/side:',I3,3X,'Density-le/te:',F6.2)
 1100 FORMAT(  1X,A27,
     &   't/c  @  x/c   camber @  x/c    r_le    t_te   A0deg')
C
 1120 FORMAT(1X,I2,3X,A20,F6.4,2X,F6.4,2X,F6.4,2X,F6.4,2X,
     &       F6.4,2X,F6.4,1X,F6.2)
C
      RETURN
      END   ! AFDISP
C
C-------------------------------------------------------------------



      SUBROUTINE TOVC (NPX,NAFX,IAF,NP,XX,YY,THOVCD,THLOC)
C------------------------------------------------------------------
C    Generalized routine for calculating approx. max thickness/chord 
C    and location (x/c) of arbitrary airfoils. Assumes unit chord.
C------------------------------------------------------------------
      REAL XX(NPX,NAFX),YY(NPX,NAFX)
C
      NP2=NP/2
      THICKL=0.0
C
      DO 100 J=10,NP2
C            
         XPNT1=XX(J,IAF)
         YPNT1=YY(J,IAF)

         DO K=J+1,NP
            IF (XX(K,IAF) .GE. XPNT1)THEN
               GO TO 50
            END IF
         ENDDO
C
 50      XPNT2=XX(K,IAF)
         YPNT2=YY(K,IAF)
C
         THICK=YPNT1-YPNT2
C
         IF (THICK .LE. THICKL)THEN
            GO TO 200
         END IF
C           
         THICKL=THICK
C
 100  CONTINUE
C
 200  THOVCD = THICKL
      THLOC  = XX(J-1,IAF)
C
      RETURN
      END  ! TOVC
C
C-------------------------------------------------------------------




      SUBROUTINE AFPANEL(NPX,NAFX,IAF,NP,XX,YY,NPP2,XXP,YYP,PAF)
C-------------------------------------------------------------------
C    AFPANEL
C    Repaneling of arbitrary airfoils with equivalent no. of points
C    on upper and lower surfaces with spacing proportional to S.
C    Sets up airfoils for blending.
C
C    NPP2 -  no of points per side
C    PAF  -  coordinate spacing_TE/spacing_LE
C------------------------------------------------------------------- 
      REAL XX(NPX,NAFX),YY(NPX,NAFX),XXP(NPX,NAFX),YYP(NPX,NAFX),
     &XB(NPX),YB(NPX),XBP(NPX),YBP(NPX),SB(NPX),STINT(NPX),SPNT(NPX)
C
      INTEGER NP(NAFX)
C
      NB=NP(IAF)
      NPP=NPP2*2
C
      DO J=1,NB
         XB(J)=XX(J,IAF)
         YB(J)=YY(J,IAF)
      ENDDO
C
C---- initialize spline
C
      CALL SCALC (XB,YB, SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
C-----------------------------------------------------------------------
C---- old LE find routine
C---- find the Leading Edge vs. S
C
c      DO J=2,NB
c         IF (YB(J) .LT. 0.0) THEN
c            GO TO 100
c         END IF          
c      ENDDO  
C
c      WRITE(*,*)'Error: LE hunt did not converge'
c      GO TO 800
C
c 100  LE1=J-1
c      LE2=J
C
c      XLE1=XB(LE1)
c      YLE1=YB(LE1)
c      XLE2=XB(LE2)
c      YLE2=YB(LE2)
C
c      SLE1=SB(LE1)
c      SLE2=SB(LE2)
C
c      SLE=SLE1+(SLE2-SLE1)*YLE1/(YLE1-YLE2)
C
C---- confirm LE vs S
C
c      XLE=SEVAL(SLE,XB,XBP,SB,NB)
c      YLE=SEVAL(SLE,YB,YBP,SB,NB)
C
C      WRITE(UNIT=*,FMT=1000)IAF, XLE,YLE
C
 1000 FORMAT(/,' Airfoil ',I2,
     $           ' le at x = ',F6.4,', y = ', F6.4)
C
C-----------------------------------------------------------------------
C----xfoil LEFIND routine is better
C
      CALL LEFIND(SLE,XB,XBP,YB,YBP,SB,NB)
C
      XLE=SEVAL(SLE,XB,XBP,SB,NB)
      YLE=SEVAL(SLE,YB,YBP,SB,NB)
C      WRITE(*,1000)IAF,XLE,YLE
C
C----Work out raw panel intervals
C
      NINT2=2*NPP2-1
C
      DO I=1,NINT2
         STINT(I) = PAF - FLOAT(I-1)*(PAF-1.)/FLOAT(NINT2-1)
      ENDDO
C
C----Accumulate panel intervals to get unscaled S points
C
      SPNT(1)=0.0
C
      DO I=2,NPP2
         J=I-1
         K=2*I-3
         L=2*I-2
         SPNT(I)=SPNT(J)+STINT(K)+STINT(L)
      ENDDO
C
C----Calculate scale factors
C
      FACTUP = SLE/(SPNT(NPP2)+STINT(NINT2))
      SLOWER = SB(NB)-SLE
      FACTDN = SLOWER/(SPNT(NPP2)+STINT(NINT2))
C
C---- Scale lower S points
C
      DO I=2,NPP2
         K=NPP-I+1
         SPNT(K)=SB(NB)-FACTDN*SPNT(I)
      ENDDO
      SPNT(NPP)=SB(NB)
C
C---- Scale upper S points
C
      DO I=1,NPP2
         SPNT(I)= FACTUP*SPNT(I)        
      ENDDO
C
C----Calculate coordinates from S points
C
      DO I=2,NPP-1
         XXP(I,IAF)=SEVAL(SPNT(I),XB,XBP,SB,NB)
         YYP(I,IAF)=SEVAL(SPNT(I),YB,YBP,SB,NB)
      ENDDO
C
      XXP(1,IAF)=XB(1)
      YYP(1,IAF)=YB(1)
C
      XXP(NPP,IAF)=XB(NB)
      YYP(NPP,IAF)=YB(NB)
C
C----Airfoil is repaneled!
C
      WRITE(*,1010) IAF
C
 1010 FORMAT(' Airfoil ',I2,' repaneled')
C
 800  RETURN
      END  ! AFPANEL
C
C-------------------------------------------------------------------




      SUBROUTINE AFSAVE(FNAMEIN,NDSK)
C-------------------------------------------------------------------
C     Save airfoil set to esloft airfoil file    
C-------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*) FNAMEIN
      CHARACTER ANS*1, FNAME*128, FBLNK*80
      PARAMETER (NBUFX = 2000)
      INTEGER ILFRST(NAFX),ILLAST(NAFX)
      DIMENSION XT(NBUFX),YT(NBUFX)
C
      N = NDSK
      LU = 19
      FNAME = FNAMEIN
      CALL STRIP(FNAME,NF)
C
      IF(FNAME(1:1) .EQ. ' ') THEN
        CALL ASKS(' Enter airfoil set filename^',FNAME)
      ENDIF
C
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=5)
C
      WRITE(*,*)
      WRITE(*,*) 'Output file exists. Overwrite?  Y/n'
      READ (*,1000) ANS
      IF(INDEX('Nn',ANS).EQ.0) GO TO 6
C
      CLOSE(LU)
      WRITE(*,*) 'Current airfoil set not saved.'
      RETURN
C
 5    OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=90)
 6    REWIND(LU)
C
C      
C--- ESLOFT header and general data
C
      WRITE(LU,1100)
      WRITE(LU,1120) NAF
      WRITE(LU,1120) NPP2
      WRITE(LU,1130) PAF
C
      DO I=1,NAF
        WRITE(LU,1000) ENAME(I,N)
      ENDDO
C
      DO I=1,NAF
        WRITE(LU,1130) ALFZ(I,N)
      ENDDO
C
C--- Move coordinates to linear arrays and write
C
      DO I=1,NAF
        ILFRST(I)= 1 + NPP*(I-1)
        ILLAST(I)= NPP*I
        DO K=1,NPP
          L = ILFRST(I) + K-1
          XT(L) = XXE(K,I)
          YT(L) = YYE(K,I)
        ENDDO
      ENDDO
C
      IFTYPE = 2
      FBLNK = ' '
C
      CALL AWRITE(FBLNK,LU,
     &            NAF,ILFRST,ILLAST,XT,YT,
     &            SNAME, ISPARS, IFTYPE)
C
      CLOSE(LU)
      RETURN      
C
 90   WRITE(*,*) 'Bad filename.'
      WRITE(*,*) 'Current airfoil set not saved.'
      RETURN
C
 1000 FORMAT(A)
 1100 FORMAT('ESLOFT Airfoil File ')
 1120 FORMAT(I6)
 1130 FORMAT(F10.4)
C
      END  !  AFSAVE
C
C-----------------------------------------------------------



      SUBROUTINE AFLOAD(FNAMIN,ERROR,NDSK)
C-----------------------------------------------------------
C     Reads previously saved esloft airfoil file
C-----------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*)  FNAMIN
      CHARACTER*128  FNAME,LINE
      CHARACTER*80   NAMET
      LOGICAL LOPEN, ERROR
C
C---- local array for calling AREAD
      DIMENSION NT(NAFX)
C
C---- open file
C
      N = NDSK
      ERROR = .FALSE.
      LOPEN = .FALSE.
      LU = 19
      FNAME = FNAMIN
      CALL STRIP(FNAME,NF)
      IF(FNAME.NE.' ') THEN
        OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=98)
        LOPEN = .TRUE.
      ELSE
        RETURN
      ENDIF
C
      ICNT = 0
C
C---- Read header line from ESLOFT input file
C
      CALL RDLINE(LU,LINE,ICNT)
      IF(LINE.EQ.'END' .OR. LINE.EQ.'ERR') GO TO 210
C
      IF(LINE(1:6).NE.'ESLOFT') THEN
        WRITE(*,*) 'Not an ESLOFT airfoil file'
        CLOSE(LU)
        ERROR = .TRUE.
        RETURN
      ENDIF
C
C--- Read in general data
C
      CALL RDLINE(LU,LINE,ICNT)
      READ(LINE,*,ERR=210) NAF
C
      CALL RDLINE(LU,LINE,ICNT)
      READ(LINE,*,ERR=210) NPP2
      NPP = 2*NPP2
C
      CALL RDLINE(LU,LINE,ICNT)
      READ(LINE,*,ERR=210) PAF
C
      DO I=1,NAF
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) NAMET
        CALL STRIP(NAMET,NFT)
        ENAME(I,N)=NAMET
      ENDDO
C
      DO I=1,NAF
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) ALFZ(I,N)
      ENDDO
C
C---- read in airfoil coordinates
C
      CALL AREADNR(LU, NPX,NAFX, XXE,YYE,
     &             NT,NTEL,
     &             SNAME,ISPARS,IFTYPE)
C
      CALL STRIP(SNAME,NSNAME)
C
      IF(IFTYPE.EQ.0 .OR. NTEL.NE.NAF) THEN
         WRITE(*,1180) NTEL
         NAF = 0
         GOTO 300
      ENDIF  
C
      DO I=1,NAF
        IF(NT(I).NE.NPP) THEN
           WRITE(*,1190) I
           NAF = 0
           ERROR = .TRUE.
           GO TO 300
        ENDIF
      ENDDO
C
C---- save data
C
      SNAMEQ(N) = SNAME
      NAFQ(N)   = NAF   
      NPP2Q(N)  = NPP2
      PAFQ(N)   = PAF
C
      DO I=1,NAF
        DO J=1,NPP
          XXEQ(J,I,N) = XXE(J,I)
          YYEQ(J,I,N) = YYE(J,I)
        ENDDO
      ENDDO
C
C---- calculate geometry parameters
C
      NBB = NPP
      DO I=1,NAF
C
        DO J=1,NBB
          XBB(J)=XXE(J,I)
          YBB(J)=YYE(J,I)
        ENDDO
C
        CALL SCALC (XBB,YBB, SBB,NBB)
        CALL SEGSPL(XBB,XBBP,SBB,NBB)
        CALL SEGSPL(YBB,YBBP,SBB,NBB)
C
        CALL GEOPARX(XBB,XBBP,YBB,YBBP,SBB,NBB, WX1,
     &            SBLEX,CHORDBX,AREABX,RADBLEX,ANGBTEX,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKBX,CAMBRBX,XTHICKBX,XCAMBRBX)
C
        TOCE(I,N)  = THICKBX
        TLOC(I,N)  = XTHICKBX
        CAMLO(I,N) = CAMBRBX
        CAMLOC(I,N)= XCAMBRBX
        PAREA(I,N) = AREABX
        PLERAD(I,N)= RADBLEX
        PANGTE(I,N)= ANGBTEX
        TETH(I,N)  = ABS(YBB(1)-YBB(NPP))
      ENDDO
C
      DO I=1,NAF-1
        IF(TOCE(I,N).LT.TOCE(I+1,N)) THEN
          WRITE(*,1200)
          NAF = 0
          ERROR = .TRUE.
          GO TO 300
        ENDIF
      ENDDO
C
C---- end of input, close file
C
      GO TO 300
C...............................................................
   98 CONTINUE
      WRITE(*,1050) FNAME(1:NF)
      ERROR = .TRUE.
      GO TO 300
C
   99 CONTINUE
      WRITE(*,1100) FNAME(1:NF)
      NAF = 0
      ERROR = .TRUE.
      GO TO 300
C
 210  CONTINUE
      WRITE(*,1150) ICNT,FNAME(1:NF)
      NAF = 0
      ERROR = .TRUE.
C
 1050 FORMAT(/' File OPEN error:  ', A)
 1100 FORMAT(/' File READ error:  ', A)
 1150 FORMAT(/' File READ error on line ',I3,':  ', A)
 1180 FORMAT(/,
     $ ' Error reading airfoil data',/,
     $ ' Coordinates read for',I2,' foils')
 1190 FORMAT(/' Airfoil ',I2,' has differing no. of points',/,
     &        ' File read aborted')
 1200 FORMAT(/' File airfoils not of descending thickness',
     &       /' File read aborted')
C
 300  IF(LOPEN) CLOSE(LU)
      RETURN
      END  ! AFLOAD
C
C-----------------------------------------------------------




      SUBROUTINE DATLOAD(FNAMIN,ERROR,NDSK)
C-----------------------------------------------------------
C     Reads previously saved airfoil file, 
C     repanels airfoil and sorts airfoil array
C-----------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*)  FNAMIN
      CHARACTER*128  FNAME
      CHARACTER*80   NAMET
      CHARACTER*1    ANS
      DIMENSION XT(NPX,NAFX),YT(NPX,NAFX),XXP(NPX,NAFX),YYP(NPX,NAFX)
      INTEGER NT(NAFX)
      LOGICAL ERROR
C
      N = NDSK
      LU=19
      FNAME = FNAMIN
      ERROR = .FALSE.
C
      CALL AREAD(FNAME,LU,NPX,NAFX, XT,YT,
     &           NT,NRD,
     &           NAMET,ISPARS,IFTYPE)
C
      IF(IFTYPE.EQ.0) THEN
        WRITE(*,*) 'File read error. Aborted'
        ERROR = .TRUE.
        GO TO 500
      ENDIF
C
      IF(IFTYPE.EQ.1) THEN
        WRITE(*,*)'Plain airfoil file'
        CALL ASKS('Enter airfoil name^',NAMET)
      ENDIF
C
      IF(IFTYPE.EQ.4) THEN
        WRITE(*,*) 'Multi-element file, first element used'
      ENDIF
C
      CALL STRIP(NAMET,NNAMET)
C
C---- normalize and repanel airfoil
C
      IAF = 1
      CALL NORMX(NPX,NAFX,IAF,XT,YT,NT)
      CALL AFPANEL(NPX,NAFX,IAF,NT,XT,YT,NPP2,XXP,YYP,PAF)
C
C---- calculate geometry parameters
C
      NPP = 2*NPP2
      NBB = NPP
C
      DO J=1,NBB
        XBB(J)=XXP(J,IAF)
        YBB(J)=YYP(J,IAF)
      ENDDO
C
      TETHT = ABS(YBB(1)-YBB(NPP))
      IF(TETHT.LT.0.0001) THEN
        WRITE(*,*)'TE thickness/chord must exceed 0.0001'
        WRITE(*,*)'Airfoil not loaded'
        ERROR = .TRUE.
        GOTO 500
      ENDIF
C
      CALL SCALC (XBB,YBB, SBB,NBB)
      CALL SEGSPL(XBB,XBBP,SBB,NBB)
      CALL SEGSPL(YBB,YBBP,SBB,NBB)
C
      CALL GEOPARX(XBB,XBBP,YBB,YBBP,SBB,NBB, WX1,
     &            SBLEX,CHORDBX,AREABX,RADBLEX,ANGBTEX,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKBX,CAMBRBX,XTHICKBX,XCAMBRBX)
C
C---- deal with duplicate t/c
C---- adjust DELTC for min del t/c
C
      DELTC = 0.005
C
      DO I=1,NAF
         DTOC= ABS(THICKBX-TOCE(I,N))
         IF(DTOC.LT.DELTC) THEN
            WRITE(*,1230)I,I
            READ(*,1000) ANS
            IF(INDEX('Nn',ANS).EQ.0) THEN
              L=I
              NAF= NAF-1
              GO TO 50
            ELSE
              WRITE(*,*)'Airfoil not loaded'
              GO TO 500
            ENDIF
         ENDIF
      ENDDO
C
C---- shuffle airfoil array with decreasing t/c
C
      DO I=1,NAF
        IF(THICKBX.GT.TOCE(I,N)) THEN
          DO J=NAF,I,-1
            DO K=1,NPP
              XXE(K,J+1)= XXE(K,J)
              YYE(K,J+1)= YYE(K,J)
            ENDDO
C
            TOCE  (J+1,N)= TOCE(J,N)
            TLOC  (J+1,N)= TLOC(J,N)
            CAMLO (J+1,N)= CAMLO(J,N)
            CAMLOC(J+1,N)= CAMLOC(J,N)
            PAREA (J+1,N)= PAREA(J,N)
            PLERAD(J+1,N)= PLERAD(J,N)
            PANGTE(J+1,N)= PANGTE(J,N)
            TETH  (J+1,N)= TETH(J,N)
            ENAME (J+1,N)= ENAME(J,N)
            ALFZ  (J+1,N)= ALFZ(J,N)
          ENDDO
C
          L=I
          GO TO 50
        ENDIF
      ENDDO
C
      L=NAF+1
C
C---- load arrays
C
 50   DO K=1,NPP
        XXE(K,L)=XXP(K,IAF)
        YYE(K,L)=YYP(K,IAF)
      ENDDO
C
      TOCE(L,N)  = THICKBX
      TLOC(L,N)  = XTHICKBX
      CAMLO(L,N) = CAMBRBX
      CAMLOC(L,N)= XCAMBRBX
      PAREA(L,N) = AREABX
      PLERAD(L,N)= RADBLEX
      PANGTE(L,N)= ANGBTEX
      ENAME(L,N) = NAMET
      TETH(L,N)  = ABS(YYE(1,L)-YYE(NPP,L))
C
      ALFZ(L,N) = 0.0
      CALL ASKR('Enter zero-lift alpha (deg)^',ALFZ(L,N))
      NAF = NAF+1
C
C
 1000 FORMAT(A)          
 1020 FORMAT(/' Max thickness = ',F7.4,' @ ',F7.4 )
 1230 FORMAT(
     & /,' Imported airfoil has similar t/c as airfoil ',I2,
     & /,' Overwrite airfoil ',I2,' ?  Y/n')
C
 500  RETURN
      END  !  DATLOAD
C
C----------------------------------------------------------------
C


      SUBROUTINE DATSAVE(IDAT,NDSK)
C----------------------------------------------------------------
C     Saves section from parent airfoil set to .dat file
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER ANS*1, FNAME*128
      LU = 19
      N = NDSK
C
      CALL ASKS(' Enter airfoil filename^',FNAME)
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=5)
C
      WRITE(*,*) 'Output file exists. Overwrite?  Y/n'
      READ (*,1000) ANS
      IF(INDEX('Nn',ANS).EQ.0) GO TO 6
C
      CLOSE(LU)
      WRITE(*,*) 'Airfoil not saved.'
      RETURN
C
 5    OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=90)
 6    REWIND(LU)
C 
C--- Write data
C
      WRITE(LU,1040) ENAME(IDAT,N)
      DO I=1,NPP
         WRITE(LU,1045) XXE(I,IDAT), YYE(I,IDAT)
      ENDDO
C
      CLOSE(LU)
      RETURN
C
 90   WRITE(*,*) 'Bad filename. Airfoil not saved.'
      RETURN
C
 1000 FORMAT(A)
 1040 FORMAT(1X,A40)
 1045 FORMAT(2(F12.6))
C
      END   !  DATSAVE
C
C-------------------------------------------------------------------




      SUBROUTINE LOFTGEOM(NR,ERROR)
C-------------------------------------------------------------------
C    Panels blade to lofting station parameters
C    Splines chord and beta to lofting stations
C    Calculates blade thickness distribution
C      ITTYPE = 1   Linear thickness/chord
C      ITTYPE = 2   Parabolic    "
C      ITTYPE = 3   Splined      "
C      ITTYPE = 4   Linear   thickness
C      ITTYPE = 5   Parabolic    "
C      ITTYPE = 6   Splined      "
C    Interpolates aero and esloft A0 data to lofting stations
C    Corrects Beta using the A0 data
C    Sets flag: data is in place to display loft or generate sections
C--------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION DSTN(NLSX), TOCS(NLSX)
      DIMENSION T1(IRX), T1S(IRX),
     &          T2(IRX), T2S(IRX),
     &          T3(IRX), T3S(IRX)
      LOGICAL ERROR,SERROR
C
      N=NR
      ERROR=.FALSE.
C
C---- first make sure data is available
C
      IF(.NOT.LLOAD) THEN
        WRITE(*,*)
        WRITE(*,*) 'No case file loaded'
        ERROR=.TRUE.
        RETURN
      ENDIF
C
      IF(NPTOT.EQ.0) THEN
        WRITE(*,*)
        WRITE(*,*) '***  No paneling available  ***'
        ERROR=.TRUE.
        RETURN
      ENDIF
C
      IF(.NOT.LBLDEF) THEN
        WRITE(*,*)
        WRITE(*,*) 'No bladed disks defined in this case'
        ERROR=.TRUE.
        RETURN
      ENDIF
C
      IF(IRTYPE(N).NE.2) THEN
        WRITE(*,*)
        WRITE(*,*) 'Current disk is not a bladed rotor'
        ERROR=.TRUE.
        RETURN
      ENDIF
C
      IF(NAF.LT.2) THEN
         WRITE(*,*)
         WRITE(*,*) 'Minimum of 2 airfoils required'
         ERROR=.TRUE.
         RETURN
      ENDIF
C
C---- all is present and correct----------------------------------- 
C---- calculate lofting stations
C
      IF(TGAP.GT.0.0 .AND.OMEGA(N).NE.0.0)THEN
         ILTIP=NRP-1
      ELSE
         ILTIP=NRP
      ENDIF
C
      YTOT = YRP(ILTIP,N)- YRP(1,N)
      YTOT2= YTOT + (OSHUB+OSTIP)/1000.
      NLINT = NLOFT-1
      DTOT = 0.0
C
      DO I=1,NLINT
         DSTN(I)= 1.0 + FLOAT(I)*(PLOFT-1.0)/FLOAT(NLINT-1)
         DTOT = DTOT + DSTN(I)
      ENDDO
C
      FAC = YTOT/DTOT
      DO I=1,NLINT
         DSTN(I)=DSTN(I)*FAC
      ENDDO
C
      YLOFT(1)= YRP(1,N)
      DO I=2,NLOFT
         YLOFT(I)=YLOFT(I-1)+DSTN(I-1)
      ENDDO
C
C---- add overshoot stations (if any)
C---- store hub and root station indices
C
      NLOFT2 = NLOFT
      IF(OSHUB.NE.0.0) THEN
         IHUB = 2
         NLOFT2 = NLOFT2+1
         DO I=NLOFT2,2,-1
           YLOFT(I)=YLOFT(I-1)
         ENDDO
         YLOFT(1)=YLOFT(2)-OSHUB/1000.
      ELSE
         IHUB = 1
      ENDIF
C
      ITIP = IHUB+NLOFT-1
C
      IF(OSTIP.NE.0.0) THEN
         NLOFT2=NLOFT2+1
         YLOFT(NLOFT2)=YLOFT(NLOFT2-1)+OSTIP/1000.
      ENDIF       
C
C---- spline blade data to lofting stations
C
      DO I=1,NRC
         T1(I) = CHR  (I,N)
         T2(I) = BETAR(I,N)
      ENDDO
C
      CALL SEGSPL(T1,T1S,YRC(1,N),NRC)
      CALL SEGSPL(T2,T2S,YRC(1,N),NRC)
C
      DO I=1,NLOFT2
         CDLOFT(I) = SEVAL(YLOFT(I),T1,T1S,YRC(1,N),NRC)
         BELOFT(I) = SEVAL(YLOFT(I),T2,T2S,YRC(1,N),NRC)
      ENDDO
C
C--------------------------------------------------------------------
C---- thickness routines
C---- get thickness at root and tip
C
      IF(TKHUB.EQ.0.0) THEN
         THUB=TOCE(1,N)*CDLOFT(IHUB)
      ELSE
         THUB=TKHUB/1000.
      ENDIF
      TOCH=THUB/CDLOFT(IHUB)
C
      IF(TKTIP.EQ.0.0) THEN
         TTIP=TOCE(NAF,N)*CDLOFT(ITIP)
      ELSE
         TTIP=TKTIP/1000.
      ENDIF
      TOCT=TTIP/CDLOFT(ITIP)
C
C
C---- direct to thickness distribution types
C
C      IF(ITTYPE.EQ.1) GO TO 100
C      IF(ITTYPE.EQ.2) GO TO 200
C      IF(ITTYPE.EQ.3) GO TO 300
C      IF(ITTYPE.EQ.4) GO TO 400
C
      IF(ITTYPE.EQ.1) GO TO 100
      IF(ITTYPE.EQ.4) GO TO 200
      IF(ITTYPE.EQ.2) GO TO 300
      IF(ITTYPE.EQ.5) GO TO 400
      IF(ITTYPE.EQ.3) GO TO 500
      IF(ITTYPE.EQ.6) GO TO 550
C
      ITTYPE = 1
      WRITE(*,*) 'Thickness index out of bounds'
      WRITE(*,*) 'Set to linear t/c'
      GO TO 800
C
C---- linear thickness/chord---------------------------------
C
 100  CONTINUE
C
      DO I=1,NLOFT2
        TCLOFT(I)= TOCH - (TOCH-TOCT)*
     &    (YLOFT(I)-YLOFT(IHUB))/(YLOFT(ITIP)-YLOFT(IHUB))
C
        THLOFT(I) = TCLOFT(I)*CDLOFT(I)
      ENDDO
C
C---- locate airfoils
C
      DO I=1,NAF
        NTC(I) = 1
        TCLOC(I,1)= YLOFT(IHUB)+(TCLOFT(IHUB)-TOCE(I,N))*
     &  (YLOFT(ITIP)-YLOFT(IHUB))/(TCLOFT(IHUB)-TCLOFT(ITIP))
      ENDDO
C
      GO TO 600
C
C---- linear thickness---------------------------------------
C
 200  CONTINUE
C
      DO I=1,NLOFT2
        THLOFT(I)= THUB - (THUB-TTIP)*
     &  (YLOFT(I)-YLOFT(IHUB))/(YLOFT(ITIP)-YLOFT(IHUB))
C
        TCLOFT(I) = THLOFT(I)/CDLOFT(I)
      ENDDO
C
      CALL TCLOCATE(N)
C
      GO TO 600
C
C---- parabolic thickness/chord---------------------------------
C
 300  CONTINUE
C
C---- traps to prevent parabolic calcs blowing up
C
      BLENGTH= YLOFT(ITIP)-YLOFT(IHUB)
      YPARA  = PAXIS*BLENGTH
      PLOC   = YLOFT(IHUB)-YPARA   
C
      IF(PLOC.GE.YLOFT(1)) THEN
         WRITE(*,*)
         WRITE(*,*)'Parabolic axis must be inboard of Station 1'
         WRITE(*,*)'Increase the value of PARA'
         ERROR = .TRUE.
         GO TO 800
      ENDIF
C
      DTC = TOCH - TOCT
      IF(DTC.LT.0.0001) THEN
         WRITE(*,*)
         WRITE(*,*)'Parabolic thickness/chord:'
         WRITE(*,*)'Hub t/c must be greater than tip t/c'
         ERROR = .TRUE.
         GO TO 800
      ENDIF
C
C---- the big equation that took so long to figure out...
C      
      APARA = (BLENGTH + 2.0*YPARA 
     &         - 2.0*SQRT(YPARA*BLENGTH+YPARA**2))/DTC**2
C
      XPARA = SQRT(YPARA/APARA)
      DO I=1,NLOFT2
         YTC = YLOFT(I)+ YPARA - YLOFT(IHUB)
         XTC = SQRT(YTC/APARA)
         TCLOFT(I) = TOCH - XTC + XPARA
         THLOFT(I) = TCLOFT(I)*CDLOFT(I)
      ENDDO
C
      CALL TCLOCATE(N)
      GO TO 600
C
C---- parabolic thickness---------------------------------------
C
 400  CONTINUE
C
C---- traps
C
      BLENGTH= YLOFT(ITIP)-YLOFT(IHUB)
      YPARA  = PAXIS*BLENGTH
      PLOC   = YLOFT(IHUB)-YPARA
C
      IF(PLOC.GE.YLOFT(1)) THEN
         WRITE(*,*)
         WRITE(*,*)'Parabolic axis must be inboard of Station 1'
         WRITE(*,*)'Increase the value of PARA'
         ERROR = .TRUE.
         GO TO 800
      ENDIF
C
      DTH = THUB - TTIP
      IF(DTH.LT.0.0001) THEN
         WRITE(*,*)
         WRITE(*,*)'Parabolic thickness:'
         WRITE(*,*)'Blade must be thicker at hub than at tip'
         ERROR = .TRUE.
         GO TO 800
      ENDIF
C
C---- the equation again
C      
      APARA = (BLENGTH + 2.0*YPARA 
     &         - 2.0*SQRT(YPARA*BLENGTH+YPARA**2))/DTH**2
C
      XPARA = SQRT(YPARA/APARA)
      DO I=1,NLOFT2
         YTH = YLOFT(I)-YLOFT(IHUB)+YPARA
         XTH = SQRT(YTH/APARA)
         THLOFT(I) = THUB - XTH + XPARA
         TCLOFT(I) = THLOFT(I)/CDLOFT(I)
      ENDDO
C
      CALL TCLOCATE(N)
      GO TO 600
C
C
C---- splined t/c------------------------------------------------
C
 500  CONTINUE
C
      IF(.NOT.LTDEF) THEN
        DO I=1,NLOFT2
          TCLOFT(I) = TOCH - (TOCH-TOCT)*
     &    (YLOFT(I)-YLOFT(IHUB))/(YLOFT(ITIP)-YLOFT(IHUB))
          THLOFT(I) = TCLOFT(I)*CDLOFT(I)
        ENDDO
      ENDIF
C
      IF(LSMOD) THEN
        CALL MODTC
        DO I=1,NLOFT2
          THLOFT(I) = TCLOFT(I)*CDLOFT(I)
        ENDDO
        TKHUB = THLOFT(IHUB)*1000.
        TKTIP = THLOFT(ITIP)*1000.
      ELSE 
        DO I=1,NLOFT2
          TCLOFT(I) = SEVAL(YLOFT(I),TT1,TT1S,TT3,NLOLD) 
        ENDDO
C
        DHUB = TOCH-TCLOFT(IHUB)
        DTIP = TOCT-TCLOFT(ITIP)
        DO I=1,NLOFT2
          TCLOFT(I) = TCLOFT(I)
     &    +DHUB*(YLOFT(ITIP)-YLOFT(I))/(YLOFT(ITIP)-YLOFT(IHUB))
     &    +DTIP*(YLOFT(I)-YLOFT(IHUB))/(YLOFT(ITIP)-YLOFT(IHUB))
          THLOFT(I) = TCLOFT(I)*CDLOFT(I)
        ENDDO
      ENDIF
C
C---- store splines for future use
C
      DO I=1,NLOFT2
        TT1(I) = TCLOFT(I)
        TT2(I) = THLOFT(I)
        TT3(I) = YLOFT(I)
      ENDDO
C
      NLOLD = NLOFT2
      CALL SEGSPL(TT1,TT1S,TT3,NLOLD)
      CALL SEGSPL(TT2,TT2S,TT3,NLOLD)
C
      CALL TCLOCATE(N)
      GO TO 600
C
C
C---- splined t---------------------------------------------------
C
 550  CONTINUE
C
      IF(.NOT.LTDEF) THEN
        DO I=1,NLOFT2
          THLOFT(I)= THUB - (THUB-TTIP)*
     &    (YLOFT(I)-YLOFT(IHUB))/(YLOFT(ITIP)-YLOFT(IHUB))
          TCLOFT(I) = THLOFT(I)/CDLOFT(I)
        ENDDO
      ENDIF
C
      IF(LSMOD) THEN
        CALL MODTH
        DO I=1,NLOFT2
          TCLOFT(I) = THLOFT(I)/CDLOFT(I)
        ENDDO
        TKHUB = THLOFT(IHUB)*1000.
        TKTIP = THLOFT(ITIP)*1000.
      ELSE 
        DO I=1,NLOFT2
          THLOFT(I) = SEVAL(YLOFT(I),TT2,TT2S,TT3,NLOLD) 
        ENDDO
C
        DHUB = THUB-THLOFT(IHUB)
        DTIP = TTIP-THLOFT(ITIP)
        DO I=1,NLOFT2
          THLOFT(I) = THLOFT(I)
     &    +DHUB*(YLOFT(ITIP)-YLOFT(I))/(YLOFT(ITIP)-YLOFT(IHUB))
     &    +DTIP*(YLOFT(I)-YLOFT(IHUB))/(YLOFT(ITIP)-YLOFT(IHUB))
          TCLOFT(I) = THLOFT(I)/CDLOFT(I)
        ENDDO
      ENDIF
C
C---- store splines for future use
C
      DO I=1,NLOFT2
        TT1(I) = TCLOFT(I)
        TT2(I) = THLOFT(I)
        TT3(I) = YLOFT(I)
      ENDDO
C
      NLOLD = NLOFT2
      CALL SEGSPL(TT1,TT1S,TT3,NLOLD)
      CALL SEGSPL(TT2,TT2S,TT3,NLOLD)
C
      CALL TCLOCATE(N)
C
C
C---- Beta calcs--------------------------------------------------
C---- interpolate DFDC and ESLOFT A0 data to lofting stations 
C
 600  CALL GETA0LOFT(N)
      CALL GETA0DFDC(N)
C
C---- perform Beta correction
C
      DO I=1,NLOFT2
         IF(OMEGA(N).LE.0.0) THEN
            BDEL = (AZDFDC(I)-AZLOFT(I))*DTR
         ELSE
            BDEL = (AZLOFT(I)-AZDFDC(I))*DTR
         ENDIF
         BECORR(I) = BELOFT(I)+BDEL
      ENDDO
C
C---- all data required to generate sections is now in place
C
      LCALC=.TRUE.
      LBLEN=.FALSE.
      LTRAN=.FALSE.
      LTDEF=.TRUE.
C
 800  RETURN
      END   ! LOFTGEOM
C
C--------------------------------------------------------------------



      SUBROUTINE GETA0LOFT(NDSK)
C--------------------------------------------------------------------
C    Linear interpolation of ESLOFT airfoils A0 to lofting stations
C    Interpolation is according to t/c, not blade radial location
C    Allows extrapolation to meet lofting stations
C    Loft AZERO inverted for negative Bgam disks
C--------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      N = NDSK
C
      IF(OMEGA(N).LE.0.0) THEN
        AZFAC = -1.0
      ELSE
        AZFAC =  1.0
      ENDIF
C
      DO 200 L=1,NLOFT2
        DO 100 I=1,NAF
          IF(TOCE(I,N).EQ.TCLOFT(L)) THEN
             AZLOFT(L)= ALFZ(I,N) * AZFAC
             GO TO 200
          ELSEIF(TOCE(I,N).LT.TCLOFT(L)) THEN
             IF(I.EQ.1) THEN
                I1 = 1
                I2 = 2
             ELSE
                I1 = I-1
                I2 = I
             ENDIF
             GO TO 150
          ELSEIF(I.EQ.NAF) THEN
             I1 = NAF-1
             I2 = NAF
          ENDIF
C
 100  CONTINUE
 150  CONTINUE
C
      FAC = (TOCE(I1,N)-TCLOFT(L))/(TOCE(I1,N)-TOCE(I2,N))
      AZLOFT(L)=(ALFZ(I1,N)-(ALFZ(I1,N)-ALFZ(I2,N))*FAC)*AZFAC
C
 200  CONTINUE
C
      RETURN
      END   ! GETA0LOFT
C
C-------------------------------------------------------------------




      SUBROUTINE GETA0DFDC(NR)
C--------------------------------------------------------------------
C    Linear interpolation of AERO airfoils A0 to lofting stations
C    Does the same calculation as the (fixed) GETCLCDCM in AERO,but
C    allows extrapolation if required to meet the lofting stations
C--------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION YADIM(NAX)
C
      N = NAERO(NR)
C
      IF(N.LT.1) THEN
         WRITE(*,*) 'AERO sections not defined'
         GO TO 300
      ENDIF
C
C---- single section needs no interpolation
C
      IF(N.EQ.1) THEN
         DO I=1,NLOFT2
           AZDFDC(I)= AERODATA(1,1,NR)/DTR
         ENDDO
         GO TO 300
      ENDIF
C
C---- put XI into dimensioned Y
C
      DO I=1,N
         YADIM(I)= XIAERO(I,NR)*RTIP(NR)
      ENDDO
C
C---- interpolate A0 to lofting stations
C
      DO 200 I=1,NLOFT2
        DO J=1,N
          IF(YADIM(J).GT.YLOFT(I))GO TO 100
        ENDDO
        J=N+1
C
 100    CONTINUE
        IF(J.EQ.1) THEN
          I1 = 1
          I2 = 2
        ELSE IF(J.GT.N) THEN
          I1 = N-1
          I2 = N
        ELSE
          I1 = J-1
          I2 = J
        ENDIF
C
        A01 = AERODATA(1,I1,NR)
        A02 = AERODATA(1,I2,NR)
C
        FAC = (YLOFT(I)-YADIM(I1))/(YADIM(I2)-YADIM(I1))
        AZDFDC(I)= (A01-(A01-A02)*FAC)/DTR  
C
 200  CONTINUE
C
 300  CONTINUE
      RETURN
      END   ! GETA0DFDC
C
C-------------------------------------------------------------------




      SUBROUTINE GETLNAMES(LONAME,NR,NDSK,ISTN,NAMTYPE,NAMOUT)
C---------------------------------------------------------------------
C     Assembles titles and filenames for various types of loft output
C     NAMTYPE= 1  Returns root name (LONAME) + disk index if NR>1
C     NAMTYPE= 2  Returns root name (+ disk index) + loft station index
C     NAMTYPE= 3  Returns root name (+ disk index) + stn index + '.txt'
C     NAMTYPE= 4  Returns root name (+ disk index) + stn index + '.dat'
C     NAMTYPE= 5  Returns root name (+ disk index) + '-radii.txt'
C     NAMTYPE= 6  Returns 'Stn ' + stn index
C     NAMTYPE= 7  Returns root name + '-cbody.txt'
C     NAMTYPE= 8  Returns root name + '-duct.txt'
C     NAMTYPE= 9  Returns root name (+ disk index) + '-ESBLADE.txt'
C
C     NR     - number of disks
C     NDSK   - current disk index
C     ISTN   - current loft station index
C     NAMOUT - output string
C---------------------------------------------------------------------
      CHARACTER*(*)LONAME
      CHARACTER*80 NAMOUT
      CHARACTER  SINDEX(99)*3, DINDEX(4)*5, SIND*3, DIND*5
C
      DATA SINDEX/ 
     $      '-01','-02','-03','-04','-05','-06','-07','-08','-09',
     $'-10','-11','-12','-13','-14','-15','-16','-17','-18','-19',
     $'-20','-21','-22','-23','-24','-25','-26','-27','-28','-29',
     $'-30','-31','-32','-33','-34','-35','-36','-37','-38','-39',
     $'-40','-41','-42','-43','-44','-45','-46','-47','-48','-49',
     $'-50','-51','-52','-53','-54','-55','-56','-57','-58','-59',
     $'-60','-61','-62','-63','-64','-65','-66','-67','-68','-69',
     $'-70','-71','-72','-73','-74','-75','-76','-77','-78','-79',
     $'-80','-81','-82','-83','-84','-85','-86','-87','-88','-89',
     $'-90','-91','-92','-93','-94','-95','-96','-97','-98','-99'/
C
      DATA DINDEX/'-Dsk1','-Dsk2','-Dsk3','-Dsk4'/
C
      IF(ISTN.LT.1 .OR. ISTN.GT.99)THEN
         WRITE(*,*) 'Station index outside range (1-99)'
         SIND = '-**'
      ELSE
         SIND = SINDEX(ISTN)
      END IF
C
      IF(NDSK.LT.1 .OR. NDSK.GT.4)THEN
         WRITE(*,*) 'Disk index outside range (1-4)'
         DIND = '-Dsk*'
      ELSE
         DIND = DINDEX(NDSK)
      END IF
C
      L=LEN(LONAME)
      DO J=L,1,-1
         IF(LONAME(J:J).NE.' ')THEN
            L=J
            GO TO 20
         ENDIF
      ENDDO
C
 20   IF(NAMTYPE.EQ.1) THEN
        IF(NR.GT.1) THEN
          NAMOUT = LONAME(1:L)//DIND
        ELSE
          NAMOUT = LONAME(1:L)
        ENDIF
C
      ELSEIF(NAMTYPE.EQ.2) THEN
        IF(NR.GT.1) THEN
          NAMOUT = LONAME(1:L)//DIND//SIND
        ELSE
          NAMOUT = LONAME(1:L)//SIND
        ENDIF
C
      ELSEIF(NAMTYPE.EQ.3) THEN
        IF(NR.GT.1) THEN
          NAMOUT=LONAME(1:L)//DIND//SIND//'.txt'
        ELSE
          NAMOUT=LONAME(1:L)//SIND//'.txt'
        ENDIF
C
      ELSEIF(NAMTYPE.EQ.4) THEN
        IF(NR.GT.1) THEN
          NAMOUT=LONAME(1:L)//DIND//SIND//'.dat'
        ELSE
          NAMOUT=LONAME(1:L)//SIND//'.dat'
        ENDIF
C
      ELSEIF(NAMTYPE.EQ.5) THEN
        IF(NR.GT.1) THEN
          NAMOUT=LONAME(1:L)//DIND//'-radii.txt'
        ELSE
          NAMOUT=LONAME(1:L)//'-radii.txt'
        ENDIF
C
      ELSEIF(NAMTYPE.EQ.6) THEN
          NAMOUT= 'Stn'//SIND
C
      ELSEIF(NAMTYPE.EQ.7) THEN
          NAMOUT=LONAME(1:L)//'-cbody.txt'
C
      ELSEIF(NAMTYPE.EQ.8) THEN
          NAMOUT=LONAME(1:L)//'-duct.txt'
C
      ELSEIF(NAMTYPE.EQ.9) THEN
        IF(NR.GT.1) THEN
          NAMOUT=LONAME(1:L)//DIND//'-ESBLADE.txt'
        ELSE
          NAMOUT=LONAME(1:L)//'-ESBLADE.txt'
        ENDIF
C
      ELSE
         WRITE(*,*) 'Name type index is out of range (GETLNAMES)'
         NAMOUT = 'ERROR'
      ENDIF
C
      RETURN
      END  ! GETLNAMES
C
C----------------------------------------------------------------------




      SUBROUTINE PRINTLOFT(LU,NDSK)
C--------------------------------------------------------------------
C    Writes loft data to terminal or disk
C--------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*80 NAMOUT,OUTSET,TTYPE,ROTAT,COORDSET
      LOGICAL ERROR, LXHUB, LXTIP
      CHARACTER*2 HUBTIP
C
      IF(.NOT.LCALC) THEN
         CALL LOFTGEOM(NDSK,ERROR)
         IF(ERROR) RETURN
      ENDIF
C
C---- header
C
      N = NDSK
      ISTN=1
      CALL GETLNAMES(LONAME,NROTOR,NDSK,ISTN,1,NAMOUT)
C
      IF(OMEGA(NDSK).GT.0.0) THEN
         ROTAT = '    Rotor: Positive rotation'      
      ELSEIF(OMEGA(NDSK).LT.0.0) THEN
         ROTAT = '    Rotor: Negative rotation'
      ELSEIF(OMEGA(NDSK).EQ.0.0) THEN
         ROTAT = 'Stator: Negative circulation'
      ENDIF
C
      WRITE(LU,*)
C      WRITE(LU,1000)
      WRITE(LU,1010) NAMOUT,ROTAT
      WRITE(LU,1000) 
C
C      IF(LTERROR) WRITE(LU,1015)
C
C---- airfoil data
C
      WRITE(LU,1035) SNAME,NPP2,PAF
C
C---- paneling and thickness data
C
      IF(ITTYPE.EQ.1) THEN
         TTYPE = 'Linear_t/c'
      ELSEIF(ITTYPE.EQ.2) THEN
         TTYPE = 'Parabolic_t/c'
      ELSEIF(ITTYPE.EQ.3) THEN
         TTYPE = 'Splined_t/c'
      ELSEIF(ITTYPE.EQ.4) THEN
         TTYPE = 'Linear_t'
      ELSEIF(ITTYPE.EQ.5) THEN
         TTYPE = 'Parabolic_t'
      ELSEIF(ITTYPE.EQ.6) THEN
         TTYPE = 'Splined_t'
      ELSE
         TTYPE = 'Error'
      ENDIF
C
      TGAPMM = TGAP*1000.
      WRITE(LU,1000)
      WRITE(LU,1104) NLOFT,PLOFT,TGAPMM,
     &               TKHUB,OSHUB,AXHUB,
     &               TKTIP,OSTIP,AXTIP,
     &               AZTIP
C
C      IF(ITTYPE.LE.2) THEN
C         WRITE(LU,1110) TTYPE
C      ELSE
C         WRITE(LU,1120) TTYPE,PAXIS,APARA
C      ENDIF
C
      IF(ITTYPE.EQ.2 .OR. ITTYPE.EQ.5) THEN
         WRITE(LU,1120) TTYPE,PAXIS,APARA
      ELSE
         WRITE(LU,1110) TTYPE
      ENDIF
C
C---- station data
C
      WRITE(LU,1000)
      WRITE(LU,1200)
C
      DO I=1,NLOFT2
        YCM  = YLOFT(I) *100.
        YIN  = YLOFT(I) *MTI
        CDCM = CDLOFT(I)*100.
        CDIN = CDLOFT(I)*MTI
        THMM = THLOFT(I)*1000.
        THIN = THLOFT(I)*MTI
C
        IF(LCOORD.AND.OMEGA(NDSK).LE.0.0) THEN
          BLDEG= (PI-BELOFT(I))/DTR
          BCDEG= (PI-BECORR(I))/DTR
        ELSE
          BLDEG= BELOFT(I)/DTR
          BCDEG= BECORR(I)/DTR
        ENDIF
C
        IF(I.EQ.IHUB) THEN
           HUBTIP='H '
        ELSEIF(I.EQ.ITIP) THEN
           HUBTIP='T '
        ELSE
           HUBTIP='  '
        ENDIF
C
        WRITE(LU,1220)I,HUBTIP,YCM,YIN,CDCM,CDIN,BLDEG,BCDEG,
     &                AZDFDC(I),AZLOFT(I),THMM,THIN,TCLOFT(I)
      ENDDO
C
C---- output settings
C
      IF(LL2D) THEN
         IF(OUTFAC.EQ.1.) THEN
           OUTSET=' Output: 2D meters'
         ELSEIF(OUTFAC.EQ.100.) THEN
           OUTSET=' Output: 2D centimeters'
         ELSEIF(OUTFAC.EQ.1000.) THEN
           OUTSET=' Output: 2D millimeters'
         ELSEIF(OUTFAC.EQ.MTI) THEN
           OUTSET=' Output: 2D inches' 
         ELSE
           OUTSET=' Output units incorrectly configured'
         ENDIF
       ELSE
         IF(OUTFAC.EQ.1.) THEN
           OUTSET=' Output: 3D meters'
         ELSEIF(OUTFAC.EQ.100.) THEN
           OUTSET=' Output: 3D centimeters'
         ELSEIF(OUTFAC.EQ.1000.) THEN
           OUTSET=' Output: 3D millimeters'
         ELSEIF(OUTFAC.EQ.MTI) THEN
           OUTSET=' Output: 3D inches' 
         ELSE
           OUTSET=' Output units incorrectly configured'
         ENDIF
       ENDIF
C
       WRITE(LU,1000)
       IF(OMEGA(NDSK).LE.0.0) THEN
          IF(LCOORD) THEN
             COORDSET= '   Local Beta Coords'
          ELSE
             COORDSET= '  Global Beta Coords'
          ENDIF
          WRITE(LU,1230) OUTSET,COORDSET
       ELSE
          WRITE(LU,1240) OUTSET
       ENDIF
C
       WRITE(LU,1000)
C
C---- check for extrapolations and warn
C
      LXHUB = .FALSE.
      LXTIP = .FALSE.
      DO I=IHUB,ITIP
        DTCHUB = TCLOFT(I)-TOCE(1,N)
        DTCTIP = TOCE(NAF,N)-TCLOFT(I)
        IF(DTCHUB.GT.0.0001) LXHUB=.TRUE.
        IF(DTCTIP.GT.0.0001) LXTIP=.TRUE.
      ENDDO
C
      IF(LXHUB) THEN
        WRITE(LU,*)
     & 'Warning: sections extrapolated beyond thickest parent'
      ENDIF
C
      IF(LXTIP) THEN
        WRITE(LU,*)
     & 'Warning: sections extrapolated beyond thinnest parent'
      ENDIF
C
C
 1000 FORMAT(1X,78('-'))
 1010 FORMAT( ' Lofted Blade: ',A32,4X,A28)
 1035 FORMAT( 
     &    ' Parents: ',A18,          'Points per side:',I4,7X,
     &    'Density -le/te  :',F6.2)
C
 1104 FORMAT(
     &    ' StationSpec    :',I4,7X, 'Density-hub/tip:',F6.2,5X,
     &    'RotorTipGap-mm  : ',F6.3,
     &  /,' ThickSpecHub-mm:',F6.2,5X,'OverShootHub-mm:',F6.2,5X,
     &    'PitchAxisHub-x/c: ',F6.3,
     &  /,' ThickSpecTip   :',F6.2,5X,'OverShootTip   :',F6.2,5X,
     &    'PitchAxisTip    : ',F6.3,
     &  /,' TipDihedral-mm :',F6.2)
C
 1110 FORMAT(' Distribution   : ',A15)
C
 1120 FORMAT(
     &       ' Distribution: ',A14,2X, 'AxisLocatn-blds:',F7.3,
     &       '      Coefficient:',G13.4)
C
 1200 FORMAT(  ' Stn Rad_cm Rad_in Cd_cm Cd_in  Bdfdc',
C               XIIXXFFFFFFFXFFFFFFFXFFFFFFXFFFFFFFFFFFFF
     &    '  Bloft A0aero A0loft  t_mm   t_in    t/c')
C          FFFFFFFFFFFFFFFFFFFFFXFFFFFFXFFFFFFXXFFFFFF
C
 1220 FORMAT(1X,I2,A2,F6.2,1X,F6.2,1X,F5.2,1X,F5.2,1X,F6.2,
     &    1X,F6.2,1X,F6.2,F6.2,1X,F6.2,1X,F6.3,2X,F6.4)
 1230 FORMAT(A50,9X,A20)
 1240 FORMAT(A50)
C
       RETURN
       END  ! PRINTLOFT
C
C------------------------------------------------------------------------



      SUBROUTINE AFBLEND(NDSK,ERROR)
C--------------------------------------------------------------------
C    Interpolates airfoils to specified thickness/chord
C--------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION XTEMP(NPX),YTEMP(NPX),STEMP(NPX),XPTMP(NPX),
     &          YPTMP(NPX)
      LOGICAL ERROR
C
      ERROR=.FALSE.
      N = NDSK
C
      IF(.NOT.LCALC) THEN
         CALL LOFTGEOM(NDSK,ERROR)
         IF(ERROR) RETURN
      ENDIF
C
      NL = NLOFT2
      DO 200 L=1,NL
        DO 100 I=1,NAF
          IF(TOCE(I,N).EQ.TCLOFT(L)) THEN
             DO J=1,NPP
               BLXX(J,L)= XXE(J,I)
               BLYY(J,L)= YYE(J,I)
             ENDDO
             GO TO 200
          ELSEIF(TOCE(I,N).LT.TCLOFT(L)) THEN
             IF(I.EQ.1) THEN
                I1 = 1
                I2 = 2
             ELSE
                I1 = I-1
                I2 = I
             ENDIF
             GO TO 150
          ELSEIF(I.EQ.NAF) THEN
             I1 = NAF-1
             I2 = NAF
          ENDIF
C          
 100  CONTINUE
 150  CONTINUE
C
      FAC = (TOCE(I1,N)-TCLOFT(L))/(TOCE(I1,N)-TOCE(I2,N))
C
      DO J=1,NPP
        BLXX(J,L) = XXE(J,I1) - FAC*(XXE(J,I1)-XXE(J,I2))
        BLYY(J,L) = YYE(J,I1) - FAC*(YYE(J,I1)-YYE(J,I2))
      ENDDO
C
 200  CONTINUE
C
C---- load geometry data arrays
C
      DO 300 I=1,NL
         DO J=1,NPP
           XTEMP(J) = BLXX(J,I)
           YTEMP(J) = BLYY(J,I)
         ENDDO
C        
         CALL SCALC (XTEMP,YTEMP,STEMP,NPP)
         CALL SEGSPL(XTEMP,XPTMP,STEMP,NPP)
         CALL SEGSPL(YTEMP,YPTMP,STEMP,NPP)
C
         CALL GEOPARX(XTEMP,XPTMP,YTEMP,YPTMP,STEMP,NPP,WX1,
     &               SBLEX,CHORDBX,AREABX,RADBLEX,ANGBTEX,
     &               EI11BA,EI22BA,APX1BA,APX2BA,
     &               EI11BT,EI22BT,APX1BT,APX2BT,
     &               THICKBX,CAMBRBX,XTHICKBX,XCAMBRBX)
C
         BLENDATA(1,I)= AREABX
         BLENDATA(2,I)= RADBLEX
         BLENDATA(3,I)= ANGBTEX
         BLENDATA(4,I)= THICKBX
         BLENDATA(5,I)= XTHICKBX
         BLENDATA(6,I)= CAMBRBX
         BLENDATA(7,I)= XCAMBRBX
         BLENDATA(8,I)= YTEMP(1)-YTEMP(NPP)
C
 300  CONTINUE
C
C---- all done
C
      LBLEN= .TRUE.
      LTRAN= .FALSE.
C
      RETURN
      END   ! AFBLEND
C
C------------------------------------------------------------------------



      SUBROUTINE AFTRANS(NDSK,ERROR)
C---------------------------------------------------------------------
C     Transforms blended airfoils to specified chord, chord axis, beta
C---------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL ERROR
      DIMENSION XTEMP(NPX),YTEMP(NPX)
C
      ERROR=.FALSE.
C
      IF(.NOT.LCALC) THEN
         CALL LOFTGEOM(NDSK,ERROR)
         IF(ERROR) RETURN
      ENDIF
C
      IF(.NOT.LBLEN) THEN
         CALL AFBLEND(NDSK,ERROR)
         IF(ERROR) RETURN
      ENDIF
C
      IF(LROTATE) THEN
        ROTFAC1=  1.0
        ROTFAC2= -1.0
      ELSE
        ROTFAC1= -1.0
        ROTFAC2=  1.0
      ENDIF
C
      NL = NLOFT2
      DO 100 L=1,NL
C
        CDAXIS= AXHUB-(AXHUB-AXTIP)*
     &         (YLOFT(L)-YLOFT(IHUB))/(YLOFT(ITIP)-YLOFT(IHUB))
C
C---- additions start here
c---- find thickness at pitch axis
c
        DO J=2,NPP2
	  IF(BLXX(J,L).LT.CDAXIS) THEN
	     NT2 = J
	     NT1 = J-1
	     YAXU = BLYY(NT1,L) + (BLXX(NT1,L)-CDAXIS) * 
     &            (BLYY(NT2,L)-BLYY(NT1,L))/(BLXX(NT1,L)-BLXX(NT2,L))
	     GOTO 20
          ENDIF
        ENDDO
C
        WRITE(*,*)'Blade thickness calculation failed'
c
 20	DO J=NPP-1,NPP2+2, -1
	  IF(BLXX(J,L).LT.CDAXIS) THEN
	     NT2 = J
             NT1 = J+1
	     YAXL = BLYY(NT1,L) + (BLXX(NT1,L)-CDAXIS) * 
     &            (BLYY(NT2,L)-BLYY(NT1,L))/(BLXX(NT1,L)-BLXX(NT2,L))
	     GOTO 30
          ENDIF
        ENDDO
c
        WRITE(*,*)'Blade thickness calculation failed'
C
 30     CONTINUE
C
C---- calculate axis, translate and scale section
C
	AXY  = (YAXU+YAXL)/2.0 
        DO J=1,NPP
          XTEMP(J) = (BLXX(J,L)-CDAXIS)*CDLOFT(L)
          YTEMP(J) = (BLYY(J,L)-AXY)   *CDLOFT(L)
        ENDDO
C
C---- rotate
C
        IF(OMEGA(NDSK).LE.0.0) THEN
          SBETA = SIN(PI-BECORR(L))
          CBETA = COS(PI-BECORR(L))
          DO J=1,NPP
            TRXX(J,L)= (CBETA*XTEMP(J) + SBETA*YTEMP(J))*ROTFAC2
            TRYY(J,L)=  CBETA*YTEMP(J) - SBETA*XTEMP(J)
          ENDDO
C
        ELSE
          SBETA = SIN(BECORR(L))
          CBETA = COS(BECORR(L))
          DO J=1,NPP
            TRXX(J,L)= (CBETA*XTEMP(J) + SBETA*YTEMP(J))*ROTFAC1
            TRYY(J,L)=  CBETA*YTEMP(J) - SBETA*XTEMP(J)
          ENDDO
        ENDIF
C
C---- translate if dihedral is not zero
C
	IF(AZTIP.NE.0.0) THEN
          AZTIPM = AZTIP/1000.
          YRAD= ((YLOFT(ITIP)-YLOFT(IHUB))**2 + AZTIPM**2) /
     &           (AZTIPM * 2.0)
          ZTHETA = ASIN((YLOFT(L)-YLOFT(IHUB))/YRAD)
          CZAXIS = YRAD*(1.-COS(ZTHETA))   
C
	  DO J=1,NPP
	    TRYY(J,L) = TRYY(J,L) + CZAXIS
          ENDDO    
        ENDIF
C
 100  CONTINUE
C
      RETURN
      END   !  AFTRANS
C
C----------------------------------------------------------------------
        

      SUBROUTINE SAVRAD(LU,NDSK,LAUTO,LOVERW)
      INCLUDE 'DFDC.INC'
      LOGICAL LAUTO,LOPEN,LOVERW
      CHARACTER ANS*1,FNAME*80,LNAMT*80,STNT*1
      CHARACTER*30 OUTUN,OUTROT,OUTDIM
C
      LOPEN=.FALSE.
      IF(LU.EQ.6) GO TO 100
C
      IF(.NOT.LAUTO) THEN
        WRITE(*,1010) LONAME
        CALL ASKS(' Enter loft radii filename^',FNAME)
        IF(FNAME.EQ.'A' .OR. FNAME.EQ.'a') RETURN
      ENDIF
C
      IF(FNAME(1:1).EQ.' ' .OR.LAUTO) THEN
        CALL GETLNAMES(LONAME,NROTOR,NDSK,1,5,FNAME)
      ENDIF
C
C---- open file
C
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=5)
      IF(LOVERW) GO TO 6
C
      WRITE(*,*)
      WRITE(*,*) 'Output file exists. Overwrite?  Y/n'
      READ (*,1000) ANS
      IF(INDEX('Nn',ANS).EQ.0) GO TO 6
C
      CLOSE(LU)
      WRITE(*,*) 'Loft radii file not saved'
      RETURN
C
 5    OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=190)
 6    REWIND(LU)
      LOPEN = .TRUE.
C
C---- write to terminal or file
C
 100  CONTINUE
      CALL GETLNAMES(LONAME,NROTOR,NDSK,1,1,LNAMT)
C
      IF(OUTFAC.EQ.MTI) THEN
         OUTUN = 'Units: inches'
      ELSEIF(OUTFAC.EQ.1000.) THEN
         OUTUN = 'Units: millimeters'
      ELSEIF(OUTFAC.EQ.100.) THEN
         OUTUN = 'Units: centimeters'
      ELSE
         OUTUN = 'Units: meters'
      ENDIF
C
      IF(OMEGA(NDSK).GT.0.0) THEN
         IF(LROTATE) THEN
           OUTROT = 'Left handed rotor'
         ELSE
           OUTROT = 'Right handed rotor'
         ENDIF
      ELSEIF(OMEGA(NDSK).EQ.0.0) THEN
         IF(LROTATE) THEN
           OUTROT = 'Right handed stator'
         ELSE
           OUTROT = 'Left handed stator'
         ENDIF         
      ELSE
         IF(LROTATE) THEN
           OUTROT = 'Right handed rotor'
         ELSE
           OUTROT = 'Left handed rotor'
         ENDIF         
      ENDIF
C
      IF(LL2D) THEN
         OUTDIM= '2D points files'
      ELSE
         OUTDIM= '3D points files'
      ENDIF
C
      WRITE(LU,1030)
      WRITE(LU,1020)
      WRITE(LU,1035) LNAMT,OUTROT,OUTDIM,OUTUN,AZTIP
      WRITE(LU,1020)
      WRITE(LU,1040)
C
      DO I=1,NLOFT2
         IF(I.EQ.IHUB) THEN
            STNT = 'H'
         ELSEIF(I.EQ.ITIP) THEN
            STNT = 'T'
         ELSE
            STNT = ' '
         ENDIF
         RADUN = YLOFT(I)*OUTFAC
         WRITE(LU,1050) I,STNT,RADUN
      ENDDO
C
      WRITE(LU,1020)
      WRITE(LU,*)
C
      IF(LOPEN) THEN
        CLOSE(LU)
        WRITE(*,1060) FNAME
      ENDIF
      RETURN
C
 190  WRITE(*,*) 'Bad filename'
      WRITE(*,*) 'Loft radii file not saved'
      RETURN
C
 1000 FORMAT(A)
 1010 FORMAT(/,' Current loft output base name: ',A30,
     &       /,' <enter> to use base name with std suffix',
     &       /,' <a>  to abort')
C
 1020 FORMAT(1X,29('-'))
 1030 FORMAT(/,'  ESLOFT Lofted Station Radii')
 1035 FORMAT(2X,A40,/, 2X,A20, /, 2X,A20, /, 2X,A20,/,
     &       2X,'Dihedral-mm:',F5.2)
C
 1040 FORMAT(1X,' Station             Radius')
C
c 1050 FORMAT(3X,I3,A1,10X,G16.7)   ! exponential output
 1050 FORMAT(3X,I3,A1,7X,F15.5)  ! decimal output
C
 1060 FORMAT(/,' Station radii written to disk: ',A30)
C
      END   ! SAVRAD
C
C----------------------------------------------------------------------




      SUBROUTINE SAVLOFT(NDSK,NPFSAVE,IPFSAVE,IPFTYPE)
C------------------------------------------------------------
C     Saves arbitrary station points-files to disk
C     IPFTYPE = 1   Blended normalized sections saved
C     IPFTYPE = 2   Transformed sections saved
C------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION IPFSAVE(NLSX)
      CHARACTER ANS*1,FNAME*80,PFNAME*80
      LOGICAL LAUTO,LOVERW
C
      LU = 19
      LAUTO = .FALSE.
      LOVERW= .FALSE.
C
      IF(IPFTYPE.EQ.1) THEN
         NAMTYPE = 4
      ELSE
         NAMTYPE = 3
      ENDIF
C
      DO 100 I= 1,NPFSAVE
        IPF = IPFSAVE(I)
C
        IF(.NOT.LAUTO) THEN
          WRITE(*,1010) IPF, LONAME
          CALL ASKS(' Enter points-file filename^',FNAME)
          IF(FNAME.EQ.'A' .OR. FNAME.EQ.'a') RETURN
          IF(FNAME(1:1).EQ.' ') LAUTO = .TRUE.
        ENDIF
C
        IF(LAUTO) THEN
          CALL GETLNAMES(LONAME,NROTOR,NDSK,IPF,NAMTYPE,FNAME)
        ENDIF
C
C---- open file
C
        OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=5)
        IF(LOVERW) GO TO 6
C
        WRITE(*,1015)
        READ (*,1000) ANS
C
        IF(ANS.EQ.'N' .OR. ANS.EQ.'n') THEN
          CLOSE(LU)
          WRITE(*,*) 'Points file not saved'
          RETURN
        ELSEIF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
          GO TO 6
        ELSEIF(ANS.EQ.' ') THEN
          LOVERW = .TRUE.
          GO TO 6
        ENDIF
C
 5      OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=190)
 6      REWIND(LU)
C
C---- write data
C
        IF(IPFTYPE.EQ.1) THEN
          CALL GETLNAMES(LONAME,NROTOR,NDSK,IPF,2,PFNAME)
          WRITE(LU,1040) PFNAME
          DO J=1,NPP
            WRITE(LU,1045) BLXX(J,IPF),BLYY(J,IPF)
          ENDDO
        ELSE
          ZTEMP = YLOFT(IPF)* OUTFAC
          DO J=1,NPP
            XTEMP= TRXX(J,IPF)* OUTFAC
            YTEMP= TRYY(J,IPF)* OUTFAC
            IF(LL2D) THEN
              WRITE(LU,1050) XTEMP,YTEMP
            ELSE
              WRITE(LU,1060) XTEMP,YTEMP,ZTEMP
            ENDIF
          ENDDO
        ENDIF
C
        CLOSE(LU)
        WRITE(*,1080) IPF,FNAME
C
 100  CONTINUE
C
C---- save station radii file if all stations are saved
C
      IF(NPFSAVE.EQ.NLOFT2 .AND. IPFTYPE.EQ.2) THEN
        CALL SAVRAD(LU,NDSK,LAUTO,LOVERW)
      ENDIF
C
      RETURN
C
 190  WRITE(*,*) 'Bad filename'
      WRITE(*,*) 'Points file not saved'
      RETURN
C
C
C---- formats
C
 1000 FORMAT(A)
 1010 FORMAT(/,' Loft station ',I2,
     &       /,' Current loft output base name: ',A30,
     &       /,' <enter> to use base name with std suffixes',
     &       /,' <a>  to abort')
C
 1015 FORMAT(/,' Output file exists. Overwrite?  y/n',
     &       /,' <enter> to overwrite all')
C
 1040 FORMAT(1X,A40)
 1045 FORMAT(2(F12.6))
C
 1050 FORMAT(2(F13.7))
 1060 FORMAT(3(F13.7))
C 1050 FORMAT(2(G16.7))
C 1060 FORMAT(3(G16.7))
C
 1080 FORMAT(' Station ',I2,' written to disk: ',A30)
C
      END  !  SAVLOFT
C
c-----------------------------------------------------------------------


      SUBROUTINE MODTC
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C------------------------------------------------
C     Takes user cursor input to modify 
C     current thickness/chord array.
C     Modified from MODCH
C     PJC, Esotec Developements, Sept 2011, Aug 2013
C------------------------------------------------
      DIMENSION YLIMS(2)      
      EXTERNAL PLCHAR,PLMATH
C
      PLFAC = 0.95
      PLPAR = 1.20*PAR
      SH = 0.3*CSIZE
      NSPLT = 20
C
C---- work with temporary arrays
      TCMAX = TCLOFT(1)
      DO I=1, NLOFT2
        WE1(I) = YLOFT(I)/YLOFT(ITIP)
        WE2(I) = TCLOFT(I)
        TCMAX = MAX(TCLOFT(I),TCMAX)
      ENDDO
      CALL SPLINE(WE2,WE3,WE1,NLOFT2)
C
      DY = 0.01
      YLIMS(1) = 0.
      YLIMS(2) = 1.1*TCMAX
      CALL PLTMODL(NLOFT2,WE1,WE2,DY,YLIMS,
     &            PLFAC,PLPAR,NSPLT,XOFF,XSF,YOFF,YSF)
C
      XPLT = -5.0*CSIZE
      YPLT = PLPAR - 0.5*DY*YSF - 0.7*CSIZE
c      YPLT = 1.5*DY*YSF - 0.7*CSIZE
      CALL PLCHAR(XPLT,YPLT,1.4*CSIZE,'t/c',0.0,3)
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
C---- get new WE2 array
      CALL CRSMOD(NLOFT2,WE1,WE2,WE3,
     &            XOFF,XSF,YOFF,YSF, SH, NSPLT,
     &            LSLOPE, IMOD1,IMOD2 )
C
C---- store as thickness/chords
      DO I = IMOD1, IMOD2
        TCLOFT(I) = WE2(I)
      ENDDO
C
      RETURN
      END   ! MODTC


C------------------------------------------------------------

      SUBROUTINE MODTH
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C------------------------------------------------
C     Takes user cursor input to modify 
C     current thickness array.
C     Modified from MODCH
C     PJC, Esotec Developements, Sept 2011, Aug 2013
C------------------------------------------------
      DIMENSION YLIMS(2)
      EXTERNAL PLCHAR,PLMATH
C
      PLFAC = 0.95
      PLPAR = 1.2*PAR
      SH = 0.3*CSIZE
      NSPLT = 20
C
C---- work with temporary arrays
      THMAX = THLOFT(1)*1000.
      DO I=1, NLOFT2
        WE1(I) = YLOFT(I)/YLOFT(ITIP)
        WE2(I) = THLOFT(I)*1000.
        THMAX = MAX(THLOFT(I)*1000.,THMAX)
      ENDDO
      CALL SPLINE(WE2,WE3,WE1,NLOFT2)
C
      DY = 1.0
      YLIMS(1) = 0.
      YLIMS(2) = 1.1*THMAX
      CALL PLTMODL(NLOFT2,WE1,WE2,DY,YLIMS,
     &            PLFAC,PLPAR,NSPLT,XOFF,XSF,YOFF,YSF)
C
      XPLT = -6.0*CSIZE
      YPLT = PLPAR - 0.5*DY*YSF - 0.7*CSIZE
c      YPLT = 1.5*DY*YSF - 0.7*CSIZE
      CALL PLCHAR(XPLT,YPLT,1.4*CSIZE,'t mm',0.0,4)
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
C---- get new WE2 array
      CALL CRSMOD(NLOFT2,WE1,WE2,WE3,
     &            XOFF,XSF,YOFF,YSF, SH, NSPLT,
     &            LSLOPE, IMOD1,IMOD2 )
C
C---- store as thickness
      DO I = IMOD1, IMOD2
        THLOFT(I) = WE2(I)/1000.
      ENDDO
C
      RETURN
      END   ! MODTH


C-------------------------------------------------------------

      SUBROUTINE PLTMODL(N,X,Y,DYMIN,YLIMS,
     &                  PLFAC,PLPAR,NSPLT,XOFF,XSF,YOFF,YSF)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      DIMENSION X(N),Y(N)
      DIMENSION YLIMS(2)
C-------------------------------------------------------------
C     Plots Y(X) array with grid overlay, assumes 0<=X<=1
C
C     user can optionally specify Y limits(min,max) or autoscaling 
C     in arrays YLIMS(1) = ymin (999. for autoscale on ymin)
C               YLIMS(2) = ymax (999. for autoscale on ymax)
C
C     Returns the scaling factors and offsets.
C
C     Intended for use with user-driven cursor interaction.
C-------------------------------------------------------------
C     Modified for use by ESLOFT
C     PJC, Esotec Developments, Sept 2011, Aug 2013
C-------------------------------------------------------------
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
      YMIN = Y(1)
      YMAX = Y(1)
      DO I=2, NLOFT2
        YMIN = MIN(YMIN,Y(I))
        YMAX = MAX(YMAX,Y(I))
      ENDDO
      IF(YLIMS(1).NE.999.) YMIN = YLIMS(1) 
      IF(YLIMS(2).NE.999.) YMAX = YLIMS(2) 
C
      IF(ABS(YMAX-YMIN).EQ.0.) THEN
        IF(YMAX.NE.0.) THEN
          YMIN = 0.5*YMAX
          YMAX = 1.5*YMAX
         ELSE
          YMAX = YMIN + 1.0
        ENDIF
      ENDIF
C
      CALL SCALIT(1,YMAX,YMIN,YFAC)
C
      YDEL = MAX( 1.0/(5.0*YFAC) , DYMIN )
      YMIN = YDEL*(AINT(YMIN/YDEL + 1000.8) - 1001.0)
      YMAX = YDEL*(AINT(YMAX/YDEL + 1001.2) - 1000.0)
C
      IF(YMIN.LT.0.0) YMIN = 0.
C
C---- increase y range if it is too narrow for useful editing
c      IF((YMAX-YMIN) .LT. 0.5*(YMAX+YMIN)) THEN
c        YMAX = YMAX + YDEL
c        YMIN = YMIN - YDEL
c      ENDIF
c      IF((YMAX-YMIN) .LT. 0.5*(YMAX+YMIN)) THEN
c        YMAX = YMAX + YDEL
c        YMIN = YMIN - YDEL
c      ENDIF
C
      DYMIN = YDEL
C
      YSF = PLPAR / (YMAX-YMIN)
      YOFF = YMIN
C
      XSF = 1./1.1
      XOFF = 0.0
C
C
      CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFAC*SIZE,LPLOT,LLAND)
      CALL PLOTABS(1.0,1.0,-3)
C
      CALL GETCOLOR(ICOL0)
C
      CS = 1.2*CSIZE
      CALL NEWPEN(2)
C
      CALL XAXIS(0.0,0.0,1.0,0.2*XSF, 0.0,0.2, CSIZE,1)
      XPLT = 0.5*XSF - 1.2*CS
      YPLT =     - 2.5*CS
      CALL PLCHAR(XPLT,YPLT,CS,'r/R',0.0,3)
C
      CALL YAXIS(0.0,0.0,PLPAR,YDEL*YSF, YMIN,YDEL, CSIZE,-2)
      IF(LGRID) THEN
       CALL NEWPEN(1)
       CALL NEWCOLORNAME('cyan')
       NXG = 11
       NYG = INT( (YMAX-YMIN)/(0.5*YDEL) + 0.01 )
       CALL PLGRID(0.0,0.0, NXG,0.1*XSF, NYG,0.5*YDEL*YSF, LMASK2 )
      ENDIF
C
      SH = 0.4*CSIZE
      CALL NEWCOLORNAME('blue')
      CALL XYSYMB(NLOFT2,X,Y,XOFF,XSF,YOFF,YSF,SH,1)
C
      CALL NEWCOLOR(ICOL0)
      CALL NEWPEN(2)
      CALL XYLINE(NLOFT2,X,Y,XOFF,XSF,YOFF,YSF,1)
C
      CALL PLFLUSH
C
      RETURN
      END   ! PLOTMODL
C
C-----------------------------------------------------


      SUBROUTINE TCLOCATE(NDSK)
      INCLUDE 'DFDC.INC'
      DIMENSION TOCS(NLSX)
      LOGICAL LHI,SERROR
C
      N = NDSK
      NYLOC = 100  ! number of intervals calculated on y
      YRX = 0.15   ! searched bladelengths inbd of root
      YTX = 0.15   ! searched bladelengths outbd of tip
C
      BLENGTH = YLOFT(ITIP)-YLOFT(IHUB)
      YTC1 = YLOFT(IHUB) - YRX*BLENGTH
      YTC2 = YLOFT(ITIP) + YTX*BLENGTH
      YINT = (YTC2-YTC1)/FLOAT(NYLOC)
C
      DO I=1,NAF
        NTC(I) = 0
      ENDDO
C
      CALL SEGSPL(TCLOFT,TOCS,YLOFT,NLOFT2)
C
      DO 200 J=1,NAF
        TCAF = TOCE(J,N)
        YLOC = YTC1 - YINT
        TCT = 0.0
        DO 100 I=1,NYLOC+1
          YLOC = YLOC + YINT
          TCTL=TCT
          TCT = SEVAL(YLOC,TCLOFT,TOCS,YLOFT,NLOFT2)
C
          IF(I.GE.2 .AND. TCT.EQ.TCTL) GOTO 100
C 
          IF(I.EQ.1) THEN
            TCT2=SEVAL(2.*YLOC,TCLOFT,TOCS,YLOFT,NLOFT2) 
            IF(TCT.LT.TCAF) THEN
              LHI = .FALSE.
            ELSEIF(TCT.GT.TCAF) THEN
              LHI = .TRUE.
            ELSE
              IF(TCT2.LE.TCT) THEN
                LHI = .TRUE.
              ELSE
                LHI = .FALSE.
              ENDIF
            ENDIF
          ENDIF
C
          IF(TCT.EQ.TCAF) THEN
            NTC(J)=NTC(J)+1
            K = NTC(J)
            TCLOC(J,K) = YLOC
            LHI = .NOT.LHI
          ENDIF
C
          IF(TCT.GT.TCAF.AND..NOT.LHI .OR. TCT.LT.TCAF.AND.LHI)THEN
            YAF = YLOC
            CALL SINVRTER(YAF,TCAF,TCLOFT,TOCS,YLOFT,NLOFT2,SERROR)
            IF(.NOT.SERROR .AND. NTC(J)+1.LE.NTCX) THEN
              NTC(J) = NTC(J)+1
              K = NTC(J)
              TCLOC(J,K) = YAF
            ENDIF
            LHI = .NOT.LHI
          ENDIF
C
 100    CONTINUE
 200  CONTINUE
C
      RETURN
      END   ! TCLOCATE


C---------------------------------------------------------------------

      SUBROUTINE GETLOFT(NDSK)
      INCLUDE 'DFDC.INC'
      N = NDSK
C
      SNAME  = SNAMEQ(N)
      NAF    = NAFQ(N)
      NPP2   = NPP2Q(N)
      NPP    = 2*NPP2
      PAF    = PAFQ(N)
C
      DO I=1,NAF
        DO J=1,NPP
          XXE(J,I) = XXEQ(J,I,N)
          YYE(J,I) = YYEQ(J,I,N)
        ENDDO
      ENDDO
C
      TKHUB = TKHUBQ(N)
      TKTIP = TKTIPQ(N)
      AXHUB = AXHUBQ(N)
      AXTIP = AXTIPQ(N)
      AZTIP = AZTIPQ(N)
      NLOFT = NLOFTQ(N)
      NLOFT2= NLOFT2Q(N)
      PLOFT = PLOFTQ(N)
      OSHUB = OSHUBQ(N)
      OSTIP = OSTIPQ(N)
      ITTYPE= ITTYPEQ(N)
      PAXIS = PAXISQ(N)
      LTDEF = LTDEFQ(N)
C
      IF(LTDEF) THEN
        DO I=1,NLOFT2
          TCLOFT(I) = TCLOFTQ(I,N)
          THLOFT(I) = THLOFTQ(I,N)
        ENDDO
      ENDIF
C
      RETURN
      END  !  GETLOFT
C
C---------------------------------------------------------------------

      SUBROUTINE PUTLOFT(NDSK)
      INCLUDE 'DFDC.INC'
      N = NDSK
C
      SNAMEQ(N) = SNAME
      NAFQ(N)   = NAF   
      NPP2Q(N)  = NPP2
      PAFQ(N)   = PAF
C
      DO I=1,NAF
        DO J=1,NPP
          XXEQ(J,I,N) = XXE(J,I)
          YYEQ(J,I,N) = YYE(J,I)
        ENDDO
      ENDDO
C
      TKHUBQ(N) = TKHUB
      TKTIPQ(N) = TKTIP
      AXHUBQ(N) = AXHUB
      AXTIPQ(N) = AXTIP
      AZTIPQ(N) = AZTIP
      NLOFTQ(N) = NLOFT
      NLOFT2Q(N)= NLOFT2
      PLOFTQ(N) = PLOFT
      OSHUBQ(N) = OSHUB
      OSTIPQ(N) = OSTIP
      ITTYPEQ(N)= ITTYPE
      PAXISQ(N) = PAXIS
      LTDEFQ(N) = LTDEF
C
      IF(LTDEF) THEN
        DO I=1,NLOFT2
          TCLOFTQ(I,N) = TCLOFT(I)
          THLOFTQ(I,N) = THLOFT(I)
        ENDDO
      ENDIF
C
      CALL GETBBFAC(N)
      WRITE(*,1010) N
C
 1010 FORMAT(/,' Blade blockage data stored for Disk',I2)
C
      RETURN
      END  !  PUTLOFT
C
C---------------------------------------------------------------------


      SUBROUTINE GETBBFAC(NDSK)
      INCLUDE 'DFDC.INC'
      DIMENSION T1(NLSX),T1S(NLSX)
      LOGICAL ERROR
C
      N=NDSK
      BLDS = FLOAT(NRBLD(N))
C
      IF(.NOT.LBLEN) THEN
         CALL AFBLEND(NDSK,ERROR)
         IF(ERROR) THEN
           WRITE(*,*) 'AFBLEND: Error calculating section areas'
           RETURN
         ENDIF
      ENDIF
C
      DO I=1,NLOFT2
        SECAR= BLENDATA(1,I)*CDLOFT(I)*CDLOFT(I)*BLDS*BBPFAC
        CYLAR= YLOFT(I)*2.0*PI * CDLOFT(I)*SIN(BECORR(I))
        T1(I)= CYLAR/(CYLAR-SECAR)
      ENDDO
C
      CALL SEGSPL(T1,T1S,YLOFT,NLOFT2)
      DO I=1,NRC
        BBVFAC(I,N) = SEVAL(YRC(I,N),T1,T1S,YLOFT,NLOFT2)
      ENDDO
C
      BPLAST = BBPFAC
      LBBLOFT(N) = .TRUE.
      RETURN
      END  ! GETBBFAC
C
C---------------------------------------------------------------------


      SUBROUTINE BBUPDATE(NDSK)
      INCLUDE 'DFDC.INC'
      LOGICAL ERROR
C
      NR = NDSK
C
      CALL GETLOFT(NR)
      CALL LOFTGEOM(NR,ERROR)
      IF(ERROR) RETURN
C
      CALL AFBLEND(NR,ERROR)
      IF(ERROR) THEN
        WRITE(*,*) 'AFBLEND: Error calculating section areas'
        RETURN
      ENDIF
C     
      CALL GETBBFAC(NR)
      WRITE(*,1010) NR
C
 1010 FORMAT(' Blade blockage data updated for Disk',I2)
C
      RETURN
      END    ! BBUPDATE
C
C----------------------------------------------------------------------
        

      SUBROUTINE SAVGEOM(IGEOM)
      INCLUDE 'DFDC.INC'
      CHARACTER ANS*1,FNAME*80
C
      LU = 19
      IG = IGEOM
C
      WRITE(*,1010) LONAME
      IF(IG.EQ.1) THEN
        CALL ASKS(' Enter centerbody filename^',FNAME)
      ELSE
        CALL ASKS(' Enter duct filename^',FNAME)
      ENDIF
C
      IF(FNAME.EQ.'A' .OR. FNAME.EQ.'a') RETURN
C
      IF(FNAME(1:1).EQ.' ') THEN
        IF(IG.EQ.1) THEN
          CALL GETLNAMES(LONAME,NROTOR,1,1,7,FNAME)
        ELSEIF(IG.EQ.2) THEN
          CALL GETLNAMES(LONAME,NROTOR,1,1,8,FNAME)
        ELSE
          FNAME = 'Name generation error'
        ENDIF
      ENDIF
C
      WRITE(*,*)
      CALL ASKI('Save CB points: 1 side or 2 ?^',ICBS)      
C
C---- open file
C
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=5)
      WRITE(*,*)
      WRITE(*,*) 'Output file exists. Overwrite?  Y/n'
      READ (*,1000) ANS
      IF(INDEX('Nn',ANS).EQ.0) GO TO 6
C
      CLOSE(LU)
      WRITE(*,*) 'Geometry file not saved'
      RETURN
C
 5    OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=190)
 6    REWIND(LU)
C
C---- Write out coordinates
C
      ZTEMP = 0.
      DO J = IPFRST(IG),IPLAST(IG)
        XTEMP= XP(J) * OUTFAC
        YTEMP= YP(J) * OUTFAC
        IF(LL2D) THEN
          WRITE(LU,1070) XTEMP,YTEMP
        ELSE
          WRITE(LU,1080) XTEMP,YTEMP,ZTEMP
        ENDIF
      ENDDO
C
C---- append second side
C
      IF(IG.EQ.1 .AND. ICBS.EQ.2) THEN
        DO J = IPLAST(IG)-1,IPFRST(IG),-1
          XTEMP=  XP(J) * OUTFAC
          YTEMP= -YP(J) * OUTFAC
          IF(LL2D) THEN
            WRITE(LU,1070) XTEMP,YTEMP
          ELSE
            WRITE(LU,1080) XTEMP,YTEMP,ZTEMP
          ENDIF
        ENDDO
      ENDIF
C
      CLOSE(LU)
C
      IF(IG.EQ.1) THEN
        WRITE(*,1090) FNAME
      ELSE
        WRITE(*,1095) FNAME
      ENDIF
C
      RETURN
C
 190  WRITE(*,*) 'Bad filename'
      WRITE(*,*) 'Geometry file not saved'
      RETURN
C
 1000 FORMAT(A)
 1010 FORMAT(/,' Current loft output base name: ',A30,
     &       /,' <enter> to use base name with suffix',
     &       /,' <a> to abort')
C
 1070 FORMAT(2(F13.7))
 1080 FORMAT(3(F13.7))
 1090 FORMAT(' Centerbody coordinates written to disk: ',A30)
 1095 FORMAT(' Duct coordinates written to disk: ',A30)
C
      END   ! SAVGEOM

C
C---------------------------------------------------------------------


      SUBROUTINE SAVEBLADE(NDSK)
      INCLUDE 'DFDC.INC'
C------------------------------------------------------------
C     Saves ESBLADE loft data file to disk
C------------------------------------------------------------
      CHARACTER ANS*1,FNAM*80,BNAM*80
C
      LU = 19
      NAMTYPE = 9
C
      WRITE(*,1010) LONAME
      CALL ASKS(' Enter ESBLADE filename^',FNAM)
      IF(FNAM.EQ.'A' .OR. FNAM.EQ.'a') RETURN
      IF(FNAM(1:1).EQ.' ') THEN
        CALL GETLNAMES(LONAME,NROTOR,NDSK,1,NAMTYPE,FNAM)
      ENDIF
C
C---- open file
C
      OPEN(LU,FILE=FNAM,STATUS='OLD',ERR=5)
      WRITE(*,1015)
      READ (*,1000) ANS
C
      IF(ANS.EQ.'N' .OR. ANS.EQ.'n') THEN
        CLOSE(LU)
        WRITE(*,*) 'File not saved'
        RETURN
      ELSEIF(ANS.EQ.'y' .OR. ANS.EQ.' ') THEN
        GO TO 6
      ENDIF
C
 5    OPEN(LU,FILE=FNAM,STATUS='NEW',ERR=190)
 6    REWIND(LU)
C
C---- write header data
C
c      IF(LCIRC) THEN
c        ICIRC = 1
c      ELSE
c        ICIRC = 0
c      ENDIF
C
      CALL GETLNAMES(LONAME,NROTOR,NDSK,1,1,BNAM)
C
      ICIRC = 0
      WRITE(LU,1100) BNAM
      WRITE(LU,1120) ICIRC,NLOFT2,NPP
C
C---- write coordinates
C
      DO IPF = 1,NLOFT2
        WRITE(LU,1130) YLOFT(IPF)
        DO J=1,NPP
          WRITE(LU,1140) TRXX(J,IPF), TRYY(J,IPF)
        ENDDO
      ENDDO
C
      WRITE(LU,*)
      CLOSE(LU)
C
      WRITE(*,1150) FNAM
      RETURN
C
 190  WRITE(*,*) 'Bad filename'
      WRITE(*,*) 'File not saved'
      RETURN
C
C
C---- formats
C
 1000 FORMAT(A)
 1010 FORMAT(/,' Current loft base name: ',A30,
     &       /,' <enter> to use base name with std suffix',
     &       /,' <a>  to abort')
 1015 FORMAT(/,' Output file exists. Overwrite?  Y/n')
C
 1100 FORMAT('ESBLADE Loft Data File generated by ESLOFT',
     &     /,'Loft name: ',A32)
 1120 FORMAT(3(I6))
 1130 FORMAT(F14.7)
 1140 FORMAT(3(F14.7))
 1150 FORMAT(/,' ESBLADE file written to disk: ', A40)
C
      END  !SAVEBLADE



C---------------------------------------------------------------------


      SUBROUTINE LOFTINIT1
      INCLUDE 'DFDC.INC'
C
      SNAME  = 'Unnamed Airfoil Set'  
      NAF    = 0      ! no airfoils loaded
      NPP2   = 60     ! points per side for paneled airfoils
      PAF    = 8.0    ! airfoil point spacing parameter (TE/LE)
      OUTFAC = 1.0    ! lofted points file units are meters
      LL2D   =.FALSE. ! lofted points files are 3D
      LLOFT  =.FALSE. ! TGAPZL is subtracted when paneling tip gap
      LROTATE=.TRUE.  ! left hand rotation selected (when BGAM>0)
      ICBS   = 1      ! 1 side included in CB points output
C
      BLFAC  = 0.5    ! default interpolation factor for MODI
      LBLBL  =.FALSE. ! blade blockage corrections not applied
      BBPFAC = 1.3    ! blade blockage profile factor (1 <--> 1.5)
      BPLAST = BBPFAC ! profile factor of last stored BB data
C
C--- plotting stuff
C
      LGGRIDX =.TRUE.  ! xfoil grid plotting is on
      LGTICKX =.TRUE.  ! xfoil ticks are on
      LGPARMX =.TRUE.  ! airfoil geo parameters will be plotted
C
      DO N=1,NRX
        SNAMEQ(N) = SNAME
        NAFQ(N)   = NAF   
        NPP2Q(N)  = NPP2
        PAFQ(N)   = PAF
      ENDDO
C
      RETURN
      END  ! LOFTINIT1
C
C--------------------------------------------------------------------


      SUBROUTINE LOFTINIT2(NR1,NR2)
      INCLUDE 'DFDC.INC'
C
      TKHUB = 0.0  ! root thickness defined by thickest airfoil
      TKTIP = 0.0  ! tip thickness defined by thinnest airfoil
      AXHUB = 0.35 ! pitch axis at root (x/c)
      AXTIP = 0.35 ! pitch axis at tip 
      AZTIP = 0.0  ! dihedral at tip (mm)
      NLOFT = 16   ! number of loft stations (not including overshoot)
      PLOFT = 1.5  ! loft station spacing parameter (tip/root)
      OSHUB = 0.0  ! no overshoot at hub
      OSTIP = 0.0  ! no overshoot at tip
      ITTYPE= 1    ! thickness defined as linear thickness/chord
      PAXIS = 0.3  ! parabolic axis, units of bladelengths inbd of hub
      GTICKX= 0.0005 ! tick length on xfoil plots, fraction of arc
      LTDEF =.FALSE. ! t undefined -initial splined t & t/c are linear
C
      DO N=NR1,NR2
        LBBLOFT(N)= .FALSE.  ! blade blockage factors not stored
        LTDEFQ(N) = LTDEF
        TKHUBQ(N) = TKHUB
        TKTIPQ(N) = TKTIP
        AXHUBQ(N) = AXHUB
        AXTIPQ(N) = AXTIP
        AZTIPQ(N) = AZTIP
        NLOFTQ(N) = NLOFT
        PLOFTQ(N) = PLOFT
        OSHUBQ(N) = OSHUB
        OSTIPQ(N) = OSTIP
        ITTYPEQ(N)= ITTYPE
        PAXISQ(N) = PAXIS
        LTDEFQ(N) = LTDEF
      ENDDO
C
      RETURN
      END  !  LOFTINIT2
C
C
C----------------------------------------------------------------------


