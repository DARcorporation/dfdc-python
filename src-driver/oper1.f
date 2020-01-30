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

      SUBROUTINE DFOPER
C--------------------------------------------------
C     Direct analysis driver routine
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      CHARACTER*1 CHKEY, ANS, CIND
      CHARACTER*16 PROMPT
      CHARACTER*8 PROMPT2
      CHARACTER*4 COMAND, COMOLD, COMPLT, CIDOTS
      LOGICAL LRECALC, LOPMOD, LANS, LMODG, LPALLWK
C
      CHARACTER*128 COMARG, ARGOLD, ARGPLT, FNAME, PLTTYPE
C
      LOGICAL LPLTEL(NEX)
C
C---- local arrays for calling GETINT,GETFLT...
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
      CHARACTER*80 LINE
C
C---- local arrays for flow survey
      PARAMETER (NFX=1)
      DIMENSION QF(2,NFX), CPF(NFX), CPF_QF(2,NFX), CPF_QINF(NFX)
C
      LOGICAL LGUI
C
C---- retain last-command info if OPER is exited and then re-entered
      SAVE COMOLD, ARGOLD
C
C---- retain last u,v sampling location for convenience
      SAVE XFVELS, YFVELS
      DATA XFVELS, YFVELS / 0.0 , 0.0 /
C
C---- generate paneled geometry for case
      CALL GENGEOM
C
      IF(NPTOT.EQ.0) THEN
       WRITE(*,*)
       WRITE(*,*) '***  No paneling available  ***'
       RETURN
      ENDIF
C
C---- "last" command (nothing yet)
      COMAND = '****'
      COMARG = ' '
      FNAME = ' '
      LRECALC = .FALSE.
      LCONV   = .FALSE.
C
C---- startup default: plot all elements
      DO IEL=1, NEL
        LPLTEL(IEL) = .TRUE.
      ENDDO
      LPALLWK = .FALSE.
C
C---- clear imposed viscous strengths
      DO IP = 1, NPTOT
        GAMVSP(IP) = 0.
        SIGVSP(IP) = 0.
      ENDDO
C
C---- plot Cp(x) or Q(x)
      ICPQXY = 1
C---- drive total thrust (ISPEC=2), drive rotor thrust for ISPEC=1
      ISPEC = 2
C---- default working disk
      NDSK = 1
C
      COMOLD = COMAND
      ARGOLD = COMARG
C
C---- "last" plot command  (start off with Cp(x))
      COMPLT = 'CPX '
      ARGPLT = ' '
C
C---- begin user-interaction loop
C............................................................
C
 500  CONTINUE
      PROMPT = '.OPER'
      PROMPT(6:7) = 'i^'
      IF(LPACC) THEN
       PROMPT(7:8) = 'a^'
      ENDIF
C
      IF(NDSK.GT.0 .AND. NROTOR.GT.1) THEN
C------ add current target disk to prompt
        IF    (NDSK.LT. 10) THEN
         PROMPT(7:8) = '(' // PNUM(NDSK)(2:2)
         K = 9
        ELSEIF(NDSK.LT.100) THEN
         PROMPT(7:9) = '(' // PNUM(NDSK   )
         K = 10
        ENDIF
        PROMPT(K:K+1) = ')^'
      ENDIF
      CALL ASKC(PROMPT,COMAND,COMARG)
C
C---- process previous command ?
      IF(COMAND(1:1).EQ.'!') THEN
        IF(COMOLD.EQ.'****') THEN
          WRITE(*,*) 'Previous .OPER command not valid'
          GO TO 500
        ELSE
          COMAND = COMOLD
          COMARG = ARGOLD
          LRECALC = .TRUE.
        ENDIF
      ELSE
        LRECALC = .FALSE.
      ENDIF
C
      IF(COMAND.EQ.'    ') THEN
       IF(LPLOT) THEN
        CALL PLEND
        CALL CLRZOOM
       ENDIF
       LPLOT = .FALSE.
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
       IF(ID.LT.0 .OR. ID.GT.NROTOR) GO TO 508
       NDSK = ID
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
C---- can jump in here if COMAND string has been set elsewhere
 510  CONTINUE
C
C--------------------------------------------------------------------
C
      IF(    COMAND.EQ.'HELP'
     &  .OR. COMAND.EQ.'?   ') THEN
C------ display menu variants depending on flow type
C
        IF(NEL.GT.1) THEN
         CIDOTS = 'ii. '
        ELSE
         CIDOTS = '    '
        ENDIF
C
C------ Menu for disk commands
        WRITE(*,2010)
 2010   FORMAT(
     &  /'   EXEC     (or E) converge current operating point',
     &  /'   THRU r   Prescribe thrust and converge case',
     &  /'   RPM  r   Prescribe disk RPM',
     & //'  .AERO     Set or edit blade section aero data',
     &  /'   DESI     design blade chord to current BGAM,CL',
     &  /'   PITC r,r Set blade pitch angle (constant or linear)',
     &  /'   CL   r,r Set blade CL (constant or linear)',
     &  /'   NBLD i   Set # blades in disk',
     &  /'   RFLG     Toggle bladed disk/actuator disk flag'
     &  /'   TSPE     Toggle thrust spec (total/rotor)',
     &  /'   TGAP     Set tip gap for rotor at shroud'
     &  /'   XDSK r   Set disk axial location',
     &  /'   NRS  i   Set # radial stations')
C
C------ Menu for flow conditions
        WRITE(*,2030)
 2030   FORMAT(
     &  /'   VINF r   Prescribe freestream speed'
     &  /'   VREF r   Change reference velocity'
     &  /'   ATMO r   Set fluid properties'
     &  /'   VSOU r   Change speed of sound')
c     & /'   VSEQ rrr Prescribe freestream speed sequence'
c     & /'   SING irr Prescribe sheet,line,point strengths'
c     & /'  .SEQ      General sequence specification')
C
C------ Menu for debugging or exotic commands
        IF(LDBG) WRITE(*,2040)
 2040   FORMAT(
     &  /'   NDOF i   Normal-velocity DOF toggle',
     &  /'  .ETYP     Toggle element-type (solid/prescribed)',
     &  /'  .TPAN i   Toggle TE panel inclusion')
C
C------ all
        WRITE(*,2050) CIDOTS,CIDOTS,CIDOTS,CIDOTS,CIDOTS,CIDOTS,CIDOTS
 2050   FORMAT(
     & /'  .DRAG     Drag area menu',
     & /'  .BL       BL analysis',
     & /'   FORC     Force summary',
     & /'   SHOW     Show current geometry,dragareas,rotor status'
     & /'   OUT      Solution output listing',
C
     &//'   CPX  ',A,'Plot Cp vs x',
     & /'   CPY  ',A,'Plot Cp vs y',
     & /'   CPAX rrr Change min,max,delta  Cp annotation',
     & /'   QX   ',A,'Plot Q vs x',
     & /'   QY   ',A,'Plot Q vs y',
     & /'   QAX  rrr Change min,max,delta  Q annotation',
     & /'   CPV  ',A,'Plot foils with pressure vectors |',
     & /'   QV   ',A,'Plot foils with velocity vectors | gee wiz',
     & /'   QN   ',A,'Plot foils with normal velocity  |',
cc     &//'   WAKE i   Generate wake from trailing edge of element i',
     & /'   WRLX     Align free vortex slipstreams with flow',
     & /'   WITR     Realign wakes with solver iterations (slooow!)',
     & /'   WMOV i   Align wake element i with flow',
C
     &//'   SIDE i   Select element side for plotting',
     & /'   Opt      Toggle plotting options',
     & /'  .ANNO     Annotate plot',
     & /'   HARD     Hardcopy current plot',
     & /'   SIZE r   Change absolute plot size',
C
     &//'   UV   rr  Get flow data at a specified location',
     & /'   UVC      Get flow data at cursor-specified locations',
     & /'   SLC      Trace streamlines at cursor-specified locations',
     & /'   CPWR f   Output x vs Cp to file',
C
     &//'   Z        Zoom  ',
     & /'   U        Unzoom',
     & /'  .OOPT     Set OPER plot and display options',
     & /'   WBOX     Set X,Y box for streamlines',
     & /'   GBOX     Set X,Y box for geometry plotting',
     & /'   VGET     Get inflow axial and tangential velocities')
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'DEBU') THEN
        LDBG = .NOT.LDBG
        IF(LDBG) THEN
          WRITE(*,*) 'DEBUG output is enabled'
        ELSE
          WRITE(*,*) 'DEBUG output is disabled'
        ENDIF
C
C--------------------------------------------------------------------
C---- Rotor data related commands
C--------------------------------------------------------------------
C
      ELSEIF(COMAND.EQ.'NRS') THEN
        CALL ASKI('Enter number of rotor radial stations^',NRSTA)
        CALL GENGEOM
C---- invalidate any existing solution
        LNCVP = .FALSE.
        LQAIC = .FALSE.
        LQGIC = .FALSE.
        LQCNT = .FALSE.
        LSYSP = .FALSE.
        LGSYS = .FALSE.
        LGAMU = .FALSE.
        LGAMA = .FALSE.
        LSIGP = .FALSE.
        LSIGM = .FALSE.
C
        NQSP = 0
        NSEG = 0
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'RPM'
     &  .OR. COMAND.EQ.'R   ') THEN
 12     CONTINUE
        NR = NDSK
        IF(NINPUT.EQ.2) THEN
         RPM = RINPUT(2)
        ELSEIF(NINPUT.EQ.1) THEN
         RPM = RINPUT(1)
        ELSE
         RPM = 30.0*OMEGA(NR)/PI
         CALL ASKR('Enter rotor RPM^',RPM)
        ENDIF
        OMEGA(NR) = RPM*PI/30.
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'NBLD') THEN
        N = NDSK
        IF(NINPUT.EQ.1) THEN
         NRBLD(N) = IINPUT(1)
        ELSE
         CALL ASKI('Enter number of rotor blades^',NRBLD(N))
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'CL') THEN
       IF(NINPUT.EQ.1) THEN
         CLBLD = RINPUT(1)
         DO IR = 1, NRC
           CLDES(IR) = CLBLD
         END DO
       ELSEIF(NINPUT.EQ.2) THEN
         CL1 = RINPUT(1)
         CL2 = RINPUT(2)
         DO IR = 1, NRC
           CLDES(IR) = CL1 + FLOAT(IR-1)/FLOAT(NRC-1) * (CL2-CL1)
         END DO
       ELSE
         CL1 = CLDES(1)
         CALL ASKR('Enter rotor blade root CL^',CL1)
         CL2 = CLDES(NRC)
         CALL ASKR('Enter rotor blade tip  CL^',CL2)
         DO IR = 1, NRC
           CLDES(IR) = CL1 + FLOAT(IR-1)/FLOAT(NRC-1) * (CL2-CL1)
         END DO
       ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'PITC') THEN
       N = 1
       IF(NINPUT.EQ.1) THEN
         ANG = RINPUT(1)
         IF(ABS(ANG).GT.90.) THEN
           CALL ASKR('Enter change in blade pitch (deg)^',ANG)
         ENDIF
         DO IR = 1, NRC
           BETAR(IR,N) = BETAR(IR,N) + ANG*DTR
         END DO
       ELSEIF(NINPUT.EQ.2) THEN
         ANG1 = RINPUT(1)
         ANG2 = RINPUT(2)
         IF(ABS(ANG1).GT.90. .OR. ABS(ANG2).GT.90.) THEN
           CALL ASKR('Enter change in hub pitch (deg)^',ANG1)
           CALL ASKR('Enter change in tip pitch (deg)^',ANG2)
         ENDIF
         DO IR = 1, NRC
           BETAR(IR,N) = BETAR(IR,N) 
     &                 + ANG1*DTR 
     &                 + FLOAT(IR-1)/FLOAT(NRC-1) * (ANG2-ANG1)*DTR
         END DO
       ELSE
         ANG = 0.0
         CALL ASKR('Enter change in blade pitch^',ANG)
         DO IR = 1, NRC
           BETAR(IR,N) = BETAR(IR,N) + ANG*DTR
         END DO
       ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'BGAM') THEN
        N = NDSK
C---- constant BGAM if one value specified
        IF(NINPUT.EQ.1) THEN
         BGAMI = RINPUT(1)
         DO IR = 1, NRC
           BGAM(IR,N) = BGAMI
         END DO
C---- linear variation of BGAM if two values specified
        ELSEIF(NINPUT.EQ.2) THEN
         BGAM1 = RINPUT(1)
         BGAM2 = RINPUT(2)
         DO IR = 1, NRC
           BGAM(IR,N) = BGAM1 + FLOAT(IR-1)/FLOAT(NRC-1) * (BGAM2-BGAM1)
         END DO
C---- three values specify disk # and linear variation of BGAM
        ELSEIF(NINPUT.EQ.3) THEN
         N = IINPUT(1)
         BGAM1 = RINPUT(2)
         BGAM2 = RINPUT(3)
         DO IR = 1, NRC
           BGAM(IR,N) = BGAM1 + FLOAT(IR-1)/FLOAT(NRC-1) * (BGAM2-BGAM1)
         END DO
C---- query user for BGAM value (constant)
        ELSE
         LINE = ' '
         CALL ASKS('Enter B*GAMMA (or <ret> to use current)^',LINE)
         IF(LINE.NE.' ') THEN
           READ(LINE,*,ERR=505) BGAMI
           DO IR = 1, NRC
             BGAM(IR,N) = BGAMI
           END DO
         ENDIF
        ENDIF
C---- set circulation to 0.0 in tip gap         
        IF(TGAP.GT.0.0) THEN
          IR = NRC
          BGAM(IR,N) = 0.0
        ENDIF
C----- initialize rotor operating point from B*GAMMA
        CALL ROTINITBGAM
C---- Reset rotor type back to actuator disk
        IRTYPE(N) = 1
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'TGAP') THEN
C---- get tip gap if one value specified
        WRITE(*,*) 'Set rotor tip gap (in input length units)'
        TGAPOLD   = TGAP
        TGAPZLOLD = TGAPZL
        IF(NINPUT.EQ.2) THEN
         TGAP   = RINPUT(1)
         TGAPZL = RINPUT(2)
        ELSEIF(NINPUT.EQ.1) THEN
         TGAP = RINPUT(1)
         CALL ASKR('Enter tip gap for zero loss^',TGAPZL)
        ELSE
         CALL ASKR('Enter rotor blade tip gap^',TGAP)
         CALL ASKR('Enter tip gap for zero loss^',TGAPZL)
        ENDIF
        TG    = TGAP   -TGAPZL
        TGOLD = TGAPOLD-TGAPZLOLD
        IF(TG.NE.TGOLD) CALL GENGEOM
C---- set circulation to 0.0 in tip gap         
        N = 1
        IF(TGAP.GT.0.0) THEN
          IR = NRC
          BGAM(IR,N) = 0.0
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'AERO') THEN
C---- modify section aerodynamic properties
        N = NDSK
        CALL AERO(N)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'RFLG') THEN
C---- switch between rotor and actuator disk for analysis
        N = NDSK
        IF(IRTYPE(N).EQ.1) THEN
          IF(LBLDEF) THEN
            IRTYPE(N) = 2
            WRITE(*,*) 'Rotor blade geometry used for analysis'
          ELSE
            WRITE(*,*) 'Define blade geometry first...'
          ENDIF
        ELSEIF(IRTYPE(N).EQ.2) THEN
          IRTYPE(N) = 1
          WRITE(*,*) 'Actuator disk used for analysis'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'TSPE') THEN
C---- toggle total/rotor thrust spec flag
        IF(ISPEC.EQ.1) THEN
          ISPEC = 2
          WRITE(*,*) 'Total thrust used for target thrust'
        ELSE
          ISPEC = 1
          WRITE(*,*) 'Rotor thrust used for target thrust'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'XDSK') THEN
        NR = NDSK
        IF(NINPUT.EQ.1) THEN
         XDISK(NR) = RINPUT(1)
        ELSEIF(NINPUT.EQ.2) THEN
         NR = IINPUT(1)
         XDISK(NR) = RINPUT(2)
        ELSE
         CALL ASKR('Enter rotor disk X location^',XDISK(NR))
        ENDIF
        CALL GENGEOM
C---- invalidate any existing solution
        LNCVP = .FALSE.
        LQAIC = .FALSE.
        LQGIC = .FALSE.
        LQCNT = .FALSE.
        LSYSP = .FALSE.
        LGSYS = .FALSE.
        LGAMU = .FALSE.
        LGAMA = .FALSE.
        LSIGP = .FALSE.
        LSIGM = .FALSE.
C
        NQSP = 0
        NSEG = 0
C
        COMOLD = COMAND
        ARGOLD = COMARG
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
       ELSEIF(COMAND.EQ.'BL') THEN
C---- BL analysis of CB and duct walls, generate viscous forces
        CALL BLOPER
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'DRAG') THEN
C---- BL analysis of CB and duct walls, generate viscous forces
       CALL DRGOPER(LMODG)
       IF(LMODG) CALL GENGEOM
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ATMO') THEN
       IF(NINPUT.EQ.1) THEN
         ALTH = RINPUT(1)
       ELSE
         CALL ASKR('Enter altitude in km^',ALTH)
       ENDIF
       CALL ATMO(ALTH,DELTAT,VSO,RHO,RMU)
       CALL FLOSHO(6, VSO, RHO, RMU)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'RHO') THEN
        IF(NINPUT.GE.1) THEN
         RHO = RINPUT(1)
        ELSE
         CALL ASKR('Enter density (rho)^',RHO)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'VSOU') THEN
        IF(NINPUT.GE.1) THEN
         VSO = RINPUT(1)
        ELSE
         CALL ASKR('Enter speed of sound^',VSO)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ITER') THEN
        IF(NINPUT.GE.1) THEN
         ITRMAXSOLV = IINPUT(1)
        ELSE
         CALL ASKI('Enter iteration limit^',ITRMAXSOLV)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'RLX') THEN
        IF(NINPUT.GE.1) THEN
         RLXSOLV = RINPUT(1)
        ELSE
         CALL ASKR('Enter velocity relaxation factor^',RLXSOLV)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'EPS') THEN
        IF(NINPUT.GE.1) THEN
         EPSSOLV = RINPUT(1)
        ELSE
         CALL ASKR('Enter convergence fraction of WX^',EPSSOLV)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'VINF'
     &  .OR. COMAND.EQ.'V   ') THEN
         QINFOLD = QINF
         IF(NINPUT.GE.1) THEN
          QINF = RINPUT(1)
         ELSE
          CALL ASKR('Enter freestream velocity^',QINF)
         ENDIF
         IF(QINF.NE.QINFOLD) LGAMA = .FALSE.
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'VREF') THEN
        IF(NINPUT.GE.1) THEN
         QREF = RINPUT(1)
        ELSE
         CALL ASKR('Enter reference velocity Vref^', QREF)
        ENDIF
C
        IF(QREF .LE. 0.0) THEN
         WRITE(*,*) 'Setting Vref = 1'
         QREF = 1.0
        ENDIF
C
        CALL QCPFOR
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'VGET') THEN
        CALL GETVEL(COMARG)
C
C--------------------------------------------------------------------
C----- design a blade from current circulation
      ELSEIF(COMAND.EQ.'DESI') THEN
        NR = NDSK
        IF(NINPUT.GE.1) NR = IINPUT(1)
        CALL DESBLADE(NR)
        LBLDEF = .TRUE.
C
C--------------------------------------------------------------------
C---- Rotor data and initialization
C
C---- initialize slipstream and circulation using blade CH,BETA geometry
      ELSEIF(COMAND.EQ.'RINB') THEN
       CALL ROTINITBLD
C
C---- initialize slipstream using B*GAM specified or current circulation
      ELSEIF(COMAND.EQ.'RING') THEN
       CALL ROTINITBGAM
C
C---- initialize slipstream using specified thrust
      ELSEIF(COMAND.EQ.'RINT') THEN
        IF(NINPUT.GE.1) THEN
         THR = RINPUT(1)
        ELSE
         THR = TTOT
         CALL ASKR('Enter rotor thrust^',THR)
        ENDIF
       CALL ROTINITTHR(THR)
C
C--------------------------------------------------------------------
C----- display slipstream velocities just downstream of rotor
      ELSEIF(COMAND.EQ.'PVEL') THEN
       NR = NDSK
       IF(NINPUT.GE.1) NR = IINPUT(1)
       CALL PRTVEL(6,'Rotor velocities',.TRUE.,.TRUE.,.TRUE.,NR)
C
C--------------------------------------------------------------------
C----- print slipstream velocities just downstream of rotor
      ELSEIF(COMAND.EQ.'WVEL') THEN
       NR = NDSK
       IF(NINPUT.GE.1) NR = IINPUT(1)
       LU = 19
       CALL OPFILE(LU,FNAME)
       CALL PRTVEL(LU,'Rotor velocities',.TRUE.,.TRUE.,.TRUE.,NR)
       CLOSE(LU)
C--------------------------------------------------------------------
C----- display blade geometry
      ELSEIF(COMAND.EQ.'DISP') THEN
       CALL ROTRPRT(6)
C
C--------------------------------------------------------------------
C----- display duct, actuator disk or blade geometry, and drag objects
      ELSEIF(COMAND.EQ.'SHOW') THEN
        CALL FLOSHO(6, VSO, RHO, RMU)
        CALL SHOWDUCT
        CALL SHOWACTDSK
        CALL SHOWBLADE
        CALL SHOWDRAGOBJ(0)
C
C--------------------------------------------------------------------
C----- recalculate solution
      ELSEIF(COMAND.EQ.'RECO') THEN
        LGAMA = .FALSE.
        CALL GAMSOLV
        COMOLD = COMAND
        ARGOLD = COMARG
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
C----- call TQCALC to get blade forces and display
      ELSEIF(COMAND.EQ.'TQCA'
     &  .OR. COMAND.EQ.'TQ  ') THEN
       IF(NINPUT.EQ.1) THEN
         ITYP = IFIX(RINPUT(1))
        ELSE
         ITYP = 1
       ENDIF
       CALL TQCALC(ITYP)
       CALL ROTRPRT(6)
C
C--------------------------------------------------------------------
C----- print out blade flow and forces
      ELSEIF(COMAND.EQ.'TQW') THEN
       ITYP = 1
       CALL TQCALC(ITYP)
       LU = 19
       CALL OPFILE(LU,FNAME)
       CALL ROTRPRT(LU)
C
C--------------------------------------------------------------------
C----- converge solution for current blade circulation BGAM (sets GTH)
      ELSEIF(COMAND.EQ.'CONV') THEN
        RLX    = RLXSOLV
        WXEPS  = EPSSOLV
        ITRMAX = ITRMAXSOLV
        RLXF = RLX
        CALL CONVGTH(ITRMAX,RLXF,WXEPS)
        ITYP = 1
        CALL TQCALC(ITYP)
        CALL ROTRPRT(6)
C
        COMOLD = COMAND
        ARGOLD = COMARG
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
C----- Calculate solution for current actuator or blade using new solver
      ELSEIF(COMAND.EQ.'EXEC'
     &  .OR. COMAND.EQ.'E   ') THEN
C
        RLX    = RLXSOLV
        WXEPS  = EPSSOLV
        ITRMAX = ITRMAXSOLV
C
        RLXF = RLX
        IF(NROTOR.GT.0) THEN
         N = 1
         IF(IRTYPE(N).EQ.1) THEN
           WRITE(*,*) 'Using actuator disk...'
           IF(.NOT.LVMAV) THEN
             CALL ROTINITBGAM
C---- Put flow data into wake grid
             CALL SETGRDFLW
           ENDIF
           CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
cc           CALL CONVGTH(ITRMAX,RLXF,WXEPS)
         ELSEIF(IRTYPE(N).EQ.2) THEN
           WRITE(*,*) 'Using current blade...'
C---- check for uninitialized BGAM for blade
           DO IR = 1, NRC
             IF(BGAM(IR,N).NE.0.0) GO TO 32
           END DO
           CALL ROTINITBLD
C---- Put flow data into wake grid
           CALL SETGRDFLW
C
 32        CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
         ENDIF
C
         RLXF = RLX
         IF(LWRLX) THEN 
           CALL WAKERESET
           IF(IRTYPE(N).EQ.1) THEN
             CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
cc             CALL CONVGTH(ITRMAX,RLXF,WXEPS)
           ELSE
             CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
           ENDIF
         ENDIF
C
         ITYP = 1
         CALL TQCALC(ITYP)
         CALL ROTRPRT(6)
        ENDIF
C
        COMOLD = COMAND
        ARGOLD = COMARG
        COMAND = COMPLT
        COMARG = ARGPLT
        IF(LCONV) GO TO 505
        GO TO 500
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'THRU') THEN
        IF(NINPUT.GE.1) THEN
         THRUST = RINPUT(1)
        ELSE
         IF(ISPEC.EQ.1) THEN
           THRUST = TTOT
           CALL ASKR('Enter rotor thrust^',THRUST)
         ELSE
           ITYP = 1
           CALL TQCALC(ITYP)
           THRUST = TTOT + TDUCT
           CALL ASKR('Enter total thrust^',THRUST)
         ENDIF
        ENDIF
C
        RLX    = RLXSOLV
        WXEPS  = EPSSOLV
        ITRMAX = ITRMAXSOLV
C
        N = 1
        IF(IRTYPE(N).EQ.2) THEN
          WRITE(*,*) 'Driving PITCH to get specified thrust'
C----- initialize rotor from blade definition
cc          IF(.NOT.LCONV) CALL ROTINITBLD
C----- if rotor is defined converge operating point from rotor geometry
          RLXF = RLX
          CALL CONVGTHBGT(ITRMAX,RLXF,WXEPS,THRUST,ISPEC)
        ELSE
          WRITE(*,*) 'Driving BGAM to get specified thrust'
C----- initialize rotor BGAM from thrust spec (constant BGAM)
          IF(.NOT.LCONV) CALL ROTINITTHR(THRUST)
C---- Put flow data into wake grid
          CALL SETGRDFLW
C----- converge operating point from specified BGAM
          RLXF = RLX
          CALL CONVGTHBGT(ITRMAX,RLXF,WXEPS,THRUST,ISPEC)
cc          CALL CONVGTHT(ITRMAX,RLXF,WXEPS,THRUST,ISPEC)
        ENDIF
        ITYP = 1
        CALL TQCALC(ITYP)
        CALL ROTRPRT(6)
C
        COMOLD = COMAND
        ARGOLD = COMARG
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SEQ ') THEN
        CALL GETSEQ
        IF(NPSEQ.LE.0) GO TO 500
C
        IF(COMPLT(1:2).EQ.'CP') THEN
          CALL CPPINI(LPLTEL,ICPQXY)
        ELSEIF(COMPLT(1:1).EQ.'Q') THEN
          CALL QPINI(LPLTEL,ICPQXY)
        ENDIF
C
C---- step through sequence of flow conditions
        DO IPSEQ = 1, NPSEQ
C
          QINF = QINFP1 + QINFPD*FLOAT(IPSEQ-1)
          DO IEL = 1, NEL
            GAMSET(IEL) = GAMP1(IEL) + GAMPD(IEL)*FLOAT(IPSEQ-1)
            SIGSET(IEL) = SIGP1(IEL) + SIGPD(IEL)*FLOAT(IPSEQ-1)
          ENDDO
C
          LGAMA = .FALSE.
          RLX    = RLXSOLV
          WXEPS  = EPSSOLV
          ITRMAX = ITRMAXSOLV
C
          RLXF = RLX
          IF(NROTOR.GT.0) THEN
           N = 1
           IF(IRTYPE(N).EQ.1) THEN
             WRITE(*,*) 'Using actuator disk...'
             IF(.NOT.LVMAV) THEN
               CALL ROTINITBGAM
             ENDIF
             CALL CONVGTH(ITRMAX,RLXF,WXEPS)
           ELSEIF(IRTYPE(N).EQ.2) THEN
             WRITE(*,*) 'Using current blade...'
C---- check for uninitialized BGAM for blade
             DO IR = 1, NRC
               IF(BGAM(IR,N).NE.0.0) GO TO 332
             END DO
             CALL ROTINITBLD
 332         CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
           ENDIF
C
           RLXF = RLX
           IF(LWRLX) THEN 
             CALL WAKERESET
             IF(IRTYPE(N).EQ.1) THEN
               CALL CONVGTH(ITRMAX,RLXF,WXEPS)
             ELSE
               CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
             ENDIF
           ENDIF
C
           ITYP = 1
           CALL TQCALC(ITYP)
           CALL ROTRPRT(6)
          ENDIF
C
c          LGAMA = .FALSE.
c          CALL GAMSOLVE
C
C---- plot case Cp or Q to display
          IF(COMPLT(1:2).EQ.'CP') THEN
            CALL CPLINE(LPLTEL,ICPQXY)
          ELSEIF(COMPLT(1:1).EQ.'Q') THEN
            CALL QPLINE(LPLTEL,ICPQXY)
          ENDIF
        ENDDO
        CALL PLFLUSH
        ARGPLT = '    '
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'FORC') THEN
        IF(NINPUT.EQ.0) THEN
C------- no element index specified... display all elements
         IEL = 0
        ELSE
C------- display only specified element
         IEL = IINPUT(1)
        ENDIF
C
        CALL PRFORC(IEL)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'OUT ') THEN
        CALL PRDUMP
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'XY  ') THEN
        CALL XYPLOT(LPLTEL)
C
        COMPLT = COMAND
        ARGPLT = COMARG
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'CPX '
     &  .OR. COMAND.EQ.'CPY ') THEN
C
        IF(NINPUT.EQ.0) THEN
C------- no element index specified... plot all elements
         DO IEL=1, NEL
           LPLTEL(IEL) = .TRUE.
           IF(.NOT.LPALLWK) THEN
            IF(NETYPE(IEL).EQ.7 .AND.
     &        (IEL2IR(IEL).NE.1 .AND. IEL2IR(IEL).NE.NRP)) THEN
             LPLTEL(IEL) = .FALSE.
            ENDIF
           ENDIF
         ENDDO
        ELSE
C------- mark only specified elements for plotting
         DO IEL=1, NEL
           LPLTEL(IEL) = .FALSE.
         ENDDO
         DO K=1, NINPUT
           IEL = IINPUT(K)
           IF(IEL.GT.0 .AND. IEL.LE.NEL) LPLTEL(IEL) = .TRUE.
         ENDDO
        ENDIF
C
        IF(COMAND.EQ.'CPX ') THEN
         ICPQXY = 1
        ELSE
         ICPQXY = 2
        ENDIF
C
        IPTYPE = JCPLT
        CALL CPPLOT(LPLTEL,ICPQXY)
C
        COMPLT = COMAND
        ARGPLT = COMARG
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'CPAX') THEN
        IF(NINPUT.LT.3) THEN
 71      RINPUT(1) = CPLMIN
         RINPUT(2) = CPLMAX
         RINPUT(3) = CPLDEL
         WRITE(*,1710) RINPUT(1), RINPUT(2), RINPUT(3)
 1710    FORMAT(/'  Enter   CPmin, CPmax, CPdel:', 3F10.4)
         CALL READR(3,RINPUT,ERROR)
         IF(ERROR) GO TO 71
        ENDIF
C
        CPLMIN = RINPUT(1)
        CPLMAX = RINPUT(2)
        CPLDEL = RINPUT(3)
C
        CPLDEL = -MAX( ABS(CPLDEL) , 0.001 )
C
        IF    (ICPQXY.EQ.1) THEN
         COMAND = 'CPX '
        ELSEIF(ICPQXY.EQ.2) THEN
         COMAND = 'CPY '
        ELSE
         COMAND = 'CPX '
        ENDIF
        COMARG = ' '
        GO TO 505
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'QX '
     &  .OR. COMAND.EQ.'QY ') THEN
C
        IF(NINPUT.EQ.0) THEN
C------- no element index specified... plot all elements
         DO IEL=1, NEL
           LPLTEL(IEL) = .TRUE.
           IF(.NOT.LPALLWK) THEN
            IF(NETYPE(IEL).EQ.7 .AND.
     &        (IEL2IR(IEL).NE.1 .AND. IEL2IR(IEL).NE.NRP)) THEN
             LPLTEL(IEL) = .FALSE.
            ENDIF
           ENDIF
         ENDDO
        ELSE
C------- mark only specified elements for plotting
         DO IEL=1, NEL
           LPLTEL(IEL) = .FALSE.
         ENDDO
         DO K=1, NINPUT
           IEL = IINPUT(K)
           IF(IEL.GT.0 .AND. IEL.LE.NEL) LPLTEL(IEL) = .TRUE.
         ENDDO
        ENDIF
C
        IF(COMAND.EQ.'QX ') THEN
         ICPQXY = 1
        ELSE
         ICPQXY = 2
        ENDIF
C
        IPTYPE = JQPLT
        CALL QPLOT(LPLTEL,ICPQXY)
C
        COMPLT = COMAND
        ARGPLT = COMARG
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'QAX') THEN
        IF(NINPUT.LT.3) THEN
 72       RINPUT(1) = QPLMIN
          RINPUT(2) = QPLMAX
          RINPUT(3) = QPLDEL
          WRITE(*,1720) RINPUT(1), RINPUT(2), RINPUT(3)
          CALL READR(3,RINPUT,ERROR)
          IF(ERROR) GO TO 72
        ENDIF
 1720   FORMAT(/'  Enter   Qmin, Qmax, Qdel:', 3F10.4)
C
        QPLMIN = RINPUT(1)
        QPLMAX = RINPUT(2)
        QPLDEL = RINPUT(3)
C
        QPLDEL = -MAX( ABS(QPLDEL) , 0.001 )
C
        IF    (ICPQXY.EQ.1) THEN
         COMAND = 'QX '
        ELSEIF(ICPQXY.EQ.2) THEN
         COMAND = 'QY '
        ELSE
         COMAND = 'QX '
        ENDIF
        COMARG = ' '
        GO TO 505
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'CPV '
     &  .OR. COMAND.EQ.'QV  '
     &  .OR. COMAND.EQ.'QN  ') THEN
C
        IF(NINPUT.EQ.0) THEN
         DO IEL=1, NEL
           LPLTEL(IEL) = .TRUE.
           IF(.NOT.LPALLWK) THEN
            IF(NETYPE(IEL).EQ.7 .AND.
     &        (IEL2IR(IEL).NE.1 .AND. IEL2IR(IEL).NE.NRP)) THEN
             LPLTEL(IEL) = .FALSE.
            ENDIF
           ENDIF
         ENDDO
        ELSE
         DO IEL=1, NEL
           LPLTEL(IEL) = .FALSE.
         ENDDO
         DO K=1, NINPUT
           IEL = IINPUT(K)
           IF(IEL.GT.0 .AND. IEL.LE.NEL) LPLTEL(IEL) = .TRUE.
         ENDDO
        ENDIF
C
        IF    (COMAND.EQ.'CPV ') THEN
         IVEC = 1
        ELSEIF(COMAND.EQ.'QV  ') THEN
         IVEC = 2
        ELSE
         IVEC = 3
        ENDIF
        IPTYPE = JAPLT
        CALL PQVEC(LPLTEL,IVEC)
C
        COMPLT = COMAND
        ARGPLT = COMARG
C
C---------------------------------------------------------------------
C--- Set options for plots and displays
      ELSEIF(COMAND.EQ.'OOPT') THEN
       WRITE(*,1998) PVFAC,QVFAC,LPALLWK
       COMOLD = COMAND
       ARGOLD = COMARG
       PROMPT2 = '.OOPT^'
       CALL ASKC(PROMPT2,COMAND,COMARG)
       IF(    COMAND.EQ.'PV') THEN
         CALL ASKR('Enter pressure vector scale factor^',PVFAC)
       ELSEIF(COMAND.EQ.'QV') THEN
        CALL ASKR('Enter velocity  vector scale factor^',QVFAC)
       ELSEIF(COMAND.EQ.'PLW') THEN
        CALL ASKL('Enter flag for all wakes plotted ^',LPALLWK)
       ENDIF
C
 1998  FORMAT(
     &  /'  Set options for vector plot displays',
     &  /'  PVfac =',G12.5,
     &  /'  QVfac =',G12.5,
     &  /'  LPWAK =',L4,
     & //'  PV        set pressure vector scaling',
     &  /'  QV        set velocity vector scaling',
     &  /'  PLW       plot all wakes flag')
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'WMOV') THEN
 77     CONTINUE
        IF(NINPUT.GE.2) THEN
         IELA = IINPUT(1)
         ISID = IINPUT(2)
        ELSEIF(NINPUT.GE.1) THEN
         IELA = IINPUT(1)
        ELSE
         IELA = 0
         CALL ASKI('Enter index of element to get wake moved^',IELA)
         ISID = ISPLOT(IELA)
         CALL ASKI('Enter side index of wake to get moved^',ISID)
        ENDIF
C
C------ align wake attached to element IELA with velocity vectors
C       on current specified side (0/1/2) from ISPLOT(IEL)
        CALL WAKMOV(IELA,ISID)
C
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
C----- toggle automatic wake iteration flag 
      ELSEIF(COMAND.EQ.'WITR') THEN
        LWRLX = .NOT.LWRLX
        IF(LWRLX) THEN
          WRITE(*,*) 'Duct wakes will be re-aligned to flow'
        ELSE
          WRITE(*,*) 'Duct wakes will be fixed'
        ENDIF
C
C--------------------------------------------------------------------
C------ align vortex wake elements bounding grid to current flow
      ELSEIF(COMAND.EQ.'WRLX') THEN
        CALL WAKERESET
C
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SIDE') THEN
       CALL ELPLSD
C
       COMAND = COMPLT
       COMARG = ARGPLT
       GO TO 505
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'OPT '
     &  .OR. COMAND.EQ.'O   ') THEN
        CALL PLTOGG(LOPMOD)
C
        IF(LOPMOD) THEN
C------- go replot current plot
         COMAND = COMPLT
         COMARG = ARGPLT
         GO TO 505
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ANNO') THEN
       IF(LPLOT) THEN
        CALL ANNOT(CHGT)
       ELSE
        WRITE(*,*) 'No active plot to annotate'
       ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'HARD') THEN
        IF(LPLOT) CALL PLEND
        LPLOT = .FALSE.
        CALL REPLOT(IDEVPS)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SIZE') THEN
        IF(NINPUT.GE.1) THEN
         SIZE = RINPUT(1)
        ELSE
         CALL ASKR('Enter plot size^',SIZE)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'CPWR') THEN
       CALL CPDUMP(COMARG)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SLC') THEN
       DX = 0.01*(XWBOX(2)-XWBOX(1))
       DY = 0.01*(YWBOX(2)-YWBOX(1))
       DS = MIN(DX,DY)
       CALL GETSLC(XS,YS,DS)
C
C--------------------------------------------------------------------
C---- Survey flow field for flow properties 
      ELSEIF(COMAND.EQ.'UV  '
     &  .OR. COMAND.EQ.'UVC ') THEN
        IF(.NOT.LPLOT) THEN
         WRITE(*,*) 'No active plot'
         GO TO 500
        ENDIF
        COMOLD = COMAND
        ARGOLD = COMARG
C
 90     CONTINUE
        IF(COMAND.EQ.'UV  ') THEN
C---- keyboard input of x,y coordinates
         IF(NINPUT.LT.2) THEN
C---- no command arguments... get x,y coordinates from prompted input
 91       RINPUT(1) = XFVELS
          RINPUT(2) = YFVELS
          WRITE(*,1791) RINPUT(1), RINPUT(2)
 1791     FORMAT(/'  Enter  x, y :  ', 2G12.5)
          CALL READR(2,RINPUT,ERROR)
          IF(ERROR) GO TO 91
         ENDIF
         XFVELS = RINPUT(1)
         YFVELS = RINPUT(2)
C
        ELSE
C---- Interactive input, set up offset and scale for type of plot 
         IF(IPTYPE.EQ.JCPLT .OR. IPTYPE.EQ.JQPLT) THEN
           XOFF = XOFA
           YOFF = YOFA
           FAC  = FACA
         ELSE
           XOFF = -XOFG
           YOFF = -YOFG
           FAC  = GSF
         ENDIF
cc         KQUIT = 1
cc         X1 = XMARG + 0.93*(XPAGE-2.0*XMARG)
cc         X2 = XMARG + 0.99*(XPAGE-2.0*XMARG)
cc         Y1 = YMARG + 0.01*(YPAGE-2.0*YMARG)
cc         Y2 = YMARG + 0.04*(YPAGE-2.0*YMARG)
cc         CALL GUIBOX(KQUIT, X1,X2,Y1,Y2, 'GREEN', ' Quit ')
C
         WRITE(*,*) 'Click on points (type "Q" to finish)...'
C---- interactive input of points with cursor (on plot)
         CALL GETCURSORXY(XPLT,YPLT,CHKEY)
         IF(INDEX('Qq',CHKEY).NE.0) GO TO 500
cc       IF(LGUI(KQUIT,XPLT,YPLT))  GO TO 500
C------- undo plot offset/scaling
         XFVELS = XPLT/FAC - XOFF
         YFVELS = YPLT/FAC - YOFF
        ENDIF
C
        CALL GETUV(XFVELS,YFVELS,UFS,VFS)
        QF(1,1) = UFS
        QF(2,1) = VFS
        QMAG = SQRT(UFS**2 + VFS**2)
C------ evaluate Cp from EIF
        NF = 1
        CALL CPCALC(NF, QF, QINF,QREF, CPF, CPF_QF,CPF_QINF)
C----- add rotor enthalpy, entropy, pressure
        CALL CPADDHS(XFVELS,YFVELS,CPF,DHF,DSF,VTF)
C----------------------------------------------------------------
C------ add point to current geometry and Cp(x) plot
        CALL GETCOLOR(ICOL0)
        CALL NEWCOLORNAME('MAGENTA')
C------ put a "+" on selected spot
        XPLT = (XFVELS + XOFF)*FAC
        YPLT = (YFVELS + YOFF)*FAC
        SSIZ = 0.8*SHGT
        CALL PLSYMB(XPLT,YPLT,SSIZ,3,0.0,0)
C
C------ put on a vector arrow for flow direction and magnitude 
        VAFAC = 0.3*QVFAC
        IP = 1
        DX = QF(1,IP)/QREF*VAFAC*FAC
        DY = QF(2,IP)/QREF*VAFAC*FAC
        X0 = (XFVELS + XOFF)*FAC
        Y0 = (YFVELS + YOFF)*FAC
        CALL ARROW(X0-0.5*DX,Y0-0.5*DY,DX,DY)
C
C------ put Cp point symbol on Cp(x) or Q(x) plot
        IF(IPTYPE.EQ.JCPLT) THEN
          CALL XYSYMB(1,XFVELS,CPF,-XOFF,FAC,0.0,-CPLFAC,SHGT,0)
        ELSEIF(IPTYPE.EQ.JQPLT) THEN
          CALL XYSYMB(1,XFVELS,QMAG/QREF,-XOFF,FAC,0.0,QPLFAC,SHGT,0)
        ENDIF
        CALL NEWCOLOR(ICOL0)
        CALL PLFLUSH
C----------------------------------------------------------------
        IP = 1
        VXF = QF(1,IP)
        VRF = QF(1,IP)
        ALF  = ATAN2(VRF,VXF)
        SWRL = ATAN2(VTF,VXF)
        WRITE(*,1800) XFVELS, YFVELS,
     &              VXF, VRF, VTF, QMAG,
     &              VXF/QREF, VRF/QREF, VTF/QREF, QMAG/QREF,
     &              CPF(IP), ALF/DTR, SWRL/DTR,
     &              DHF, DSF 
C
 1800   FORMAT(/' x     =', G11.5,'  y     =', G11.5,
     &         /' u     =', G11.5,'  v     =', G11.5,
     &         '  w     =', G11.5,'  q     =', G11.5,
     &         /' u/Qref=', G11.5,'  v/Qref=', G11.5,
     &         '  w/Qref=', G11.5,'  q/Qref=', G11.5,
     &         /' Cp    =', G11.5,'  alfa  =', G11.5,' deg',
     &                            '  swirl =', G11.5,' deg',
     &         /' delH  =', G11.5,'  delS  =', G11.5/ )
C
C------ request another point
        IF(COMAND.EQ.'UVC ') GO TO 90
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'Z   ') THEN
        IF(LPLOT) THEN
         CALL USETZOOM(.TRUE.,.TRUE.)
         CALL REPLOT(IDEV)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'U   ') THEN
        IF(LPLOT) THEN
         CALL CLRZOOM
         CALL REPLOT(IDEV)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'GBOX') THEN
        WRITE(*,1990) XGBOX(1),YGBOX(1),XGBOX(2),YGBOX(2)
C
        IF(NINPUT.GE.4) THEN
          XGBOX(1) = RINPUT(1)
          YGBOX(1) = RINPUT(2)
          XGBOX(2) = RINPUT(3)
          YGBOX(2) = RINPUT(4)
        ELSE
          WRITE(*,1991)
          PROMPT2 = '.GBOX^'
          CALL ASKC(PROMPT2,COMAND,COMARG)
C
          IF(COMAND.EQ.'D') THEN
            XGBOX(1) = 0.0
            YGBOX(1) = 0.0
            XGBOX(2) = 0.0
            YGBOX(2) = 0.0
          ELSEIF(COMAND.EQ.'K') THEN
 190        RINPUT(1) = XGBOX(1)
            RINPUT(2) = YGBOX(1)
            WRITE(*,1992) RINPUT(1), RINPUT(2)
            CALL READR(2,RINPUT,ERROR)
            IF(ERROR) GO TO 190
            XGBOX(1) = RINPUT(1)
            YGBOX(1) = RINPUT(2)
C
 192        RINPUT(1) = XGBOX(2)
            RINPUT(2) = YGBOX(2)
            WRITE(*,1993) RINPUT(1), RINPUT(2)
            CALL READR(2,RINPUT,ERROR)
            IF(ERROR) GO TO 192
            XGBOX(2) = RINPUT(1)
            YGBOX(2) = RINPUT(2)
          ELSEIF(COMAND.EQ.'C') THEN
C---- Interactive box input with cursor on plot 
            IF(.NOT.LPLOT) THEN
              WRITE(*,*) 'No active plot'
              GO TO 500
            ENDIF
            WRITE(*,*) 'Click on corner point Xmin,Ymin'
            CALL GETCURSORXY(XPLT,YPLT,CHKEY)
            IF(INDEX('Qq',CHKEY).NE.0) GO TO 500
            XGBOX(1) = XPLT/FACA - XOFA
            YGBOX(1) = YPLT/FACA - YOFA
C
            WRITE(*,*) 'Click on corner point Xmax,Ymax'
            CALL GETCURSORXY(XPLT,YPLT,CHKEY)
            IF(INDEX('Qq',CHKEY).NE.0) GO TO 500
            XGBOX(2) = XPLT/FACA - XOFA
            YGBOX(2) = YPLT/FACA - YOFA
          ENDIF
C
        ENDIF
        COMAND = COMPLT
        COMARG = ARGPLT
        IF(LPLOT) THEN
          GO TO 505
        ELSE
          GO TO 500
        ENDIF
C
 1990   FORMAT(/' Geometry box Xmin, Ymin :  ', 2(G12.5,2X)
     &         /'              Xmax, Ymax :  ', 2(G12.5,2X))
 1991   FORMAT(
     &  /'  Enter option for geometry box input',
     &  /'   C  cursor   input of corner points',
     &  /'   K  keyboard input of corner points',
     &  /'   D  default box')
 1992   FORMAT(/'  Enter box Xmin, Ymin :  ', 2(G12.5,2X))
 1993   FORMAT(/'  Enter box Xmax, Ymax :  ', 2(G12.5,2X))
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'WBOX') THEN
        WRITE(*,1994) XWBOX(1),YWBOX(1),XWBOX(2),YWBOX(2)
C
        IF(NINPUT.GE.4) THEN
          XWBOX(1) = RINPUT(1)
          YWBOX(1) = RINPUT(2)
          XWBOX(2) = RINPUT(3)
          YWBOX(2) = RINPUT(4)
        ELSE
          WRITE(*,1995)
          PROMPT2 = '.WBOX^'
          CALL ASKC(PROMPT2,COMAND,COMARG)
C
          IF(COMAND.EQ.'D') THEN
            CALL WAKEBOX
          ELSEIF(COMAND.EQ.'K') THEN
 194         RINPUT(1) = XWBOX(1)
            RINPUT(2) = YWBOX(1)
            WRITE(*,1992) RINPUT(1), RINPUT(2)
            CALL READR(2,RINPUT,ERROR)
            IF(ERROR) GO TO 194
            XWBOX(1) = RINPUT(1)
            YWBOX(1) = RINPUT(2)
C
 196        RINPUT(1) = XWBOX(2)
            RINPUT(2) = YWBOX(2)
            WRITE(*,1993) RINPUT(1), RINPUT(2)
            CALL READR(2,RINPUT,ERROR)
            IF(ERROR) GO TO 196
            XWBOX(2) = RINPUT(1)
            YWBOX(2) = RINPUT(2)
          ELSEIF(COMAND.EQ.'C') THEN
C---- Interactive box input with cursor on plot 
            IF(.NOT.LPLOT) THEN
              WRITE(*,*) 'No active plot'
              GO TO 500
            ENDIF
            WRITE(*,*) 'Click on corner point Xmin,Ymin'
            CALL GETCURSORXY(XPLT,YPLT,CHKEY)
            IF(INDEX('Qq',CHKEY).NE.0) GO TO 500
            XWBOX(1) = XPLT/FACA - XOFA
            YWBOX(1) = YPLT/FACA - YOFA
C
            WRITE(*,*) 'Click on corner point Xmax,Ymax'
            CALL GETCURSORXY(XPLT,YPLT,CHKEY)
            IF(INDEX('Qq',CHKEY).NE.0) GO TO 500
            XWBOX(2) = XPLT/FACA - XOFA
            YWBOX(2) = YPLT/FACA - YOFA
          ENDIF
C
        ENDIF
        COMAND = COMPLT
        COMARG = ARGPLT
        IF(LPLOT) THEN
          GO TO 505
        ELSE
          GO TO 500
        ENDIF
C
 1994   FORMAT(/' Wake box Xmin, Ymin :  ', 2(G12.5,2X)
     &         /'          Xmax, Ymax :  ', 2(G12.5,2X))
 1995   FORMAT(
     &  /'  Enter option for wake box input',
     &  /'   C  cursor   input of corner points',
     &  /'   K  keyboard input of corner points',
     &  /'   D  default box')
C
C---------------------------------------------------------------------
C--- Plots for blade geometry and operating point (from XROTOR)
      ELSEIF(COMAND.EQ.'PLOT') THEN
        NR = NDSK
        IF(NINPUT.GE.1) THEN
         NPLOT = IINPUT(1)
        ELSE
         WRITE(*,2000)
         NPLOT = 3
         CALL ASKI('select plot number^',NPLOT)
        ENDIF
C
        PLFAC1 = 0.6
        PLFAC2 = 0.8
        PLFACD = 0.6
        XORG = 0.15
        YORG = 0.10
C
      IF(NPLOT.EQ.0) THEN
       GO TO 505
C--- 3 view geometry plot of single blade
      ELSEIF(NPLOT.EQ.1) THEN
       IF(IRTYPE(NR).EQ.2) THEN
        CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFAC1*SIZE,LPLOT,LLAND)
        CALL PLOT(XORG,YORG,-3)
        CALL BLDPLT(NR,'ALUE')
       ELSE
        WRITE(*,*) 'Blade must be defined before plotting'
       ENDIF
C--- Geometry of all blades, axial view
      ELSEIF(NPLOT.EQ.2) THEN
       IF(IRTYPE(NR).EQ.2) THEN
        CALL PLTINI(SCRNFR,IPSLU,IDEV,SIZE,LPLOT,.NOT.LLAND)
        CALL ROTRPLT(NR)
       ELSE
        WRITE(*,*) 'Blade must be defined before plotting'
       ENDIF
C--- Plot of operating point (Gam, CL, M, eff) + data
      ELSEIF(NPLOT.EQ.3) THEN
       CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFAC2*SIZE,LPLOT,LLAND)
       CALL PLOT(XORG,YORG,-3)
       CALL BLDCLPLT(NR)
C--- Combined geometry and operating point
      ELSEIF(NPLOT.EQ.4) THEN
       IF(IRTYPE(NR).EQ.2) THEN
        CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFACD*SIZE,LPLOT,.NOT.LLAND)
        CALL PLOT(0.175,0.075,-3)
        CALL BLDPLT(NR,'AL')
        CALL PLOTABS(0.0,0.0,-3)
        CALL PLOT(0.175,0.875,-3)
        CALL BLDCLPLT(NR)
       ELSE
        WRITE(*,*) 'Blade must be defined before plotting'
       ENDIF
C--- Relative frame velocities on blade
      ELSEIF(NPLOT.EQ.7) THEN
       CALL ROTVPLT(NR,.FALSE.,.TRUE.)
C--- Absolute frame velocities on blade
      ELSEIF(NPLOT.EQ.8) THEN
       CALL ROTVPLT(NR,.TRUE.,.FALSE.)
C--- Velocity triangles
c      ELSE IF(NPLOT.EQ.9) THEN
c       CALL TRIPLT
C--- Imposed external slipstream velocities
c      ELSE IF(NPLOT.EQ.10) THEN
c       IF(NADD.LT.2) THEN
c        WRITE(*,*) 'No slipstream profiles present'
c        GO TO 900
c       ENDIF
c       CALL VELPLT
C--- Plot reference x,y data
c      ELSE IF(NPLOT.EQ.11) THEN
c       FNAME = ' '
c       CALL REFPLT(FNAME, XYOFF(1),XYOFF(2),XYFAC(1),XYFAC(2),
c     &             0.5*CSIZE, 1)
C--- Plot blade parameters vs r/R
      ELSEIF(NPLOT.EQ.12) THEN
        CALL PLOTBLDDATA(NR,PLTTYPE)
        GO TO 500
C
      ENDIF
C
 2000 FORMAT(/'  0   CANCEL'
     &       /'  1   Rotor blade geometry'
     &       /'  2   Rotor disk axial view (all blades)'
     &       /'  3   Radial distributions on blade (GAM,CL,M)'
     &       /'  4   Radial distributions plus geometry'
     &       /'  7   Blade-relative slipstream velocities vs r/R'
     &       /'  8   Absolute frame slipstream velocities vs r/R'
c     &       /'  9   Velocity triangles'
c     &       /' 10   External slipstream velocity profiles'
c     &       /' 11   Reference x,y data'
     &       /' 12   Plot blade data (Gam,CL,CD,etc) vs r/R')
C
C
C
C--------------------------------------------------------------------
C  Relatively comprehensible commands, used for cool stuff
C--------------------------------------------------------------------
C---------------------------------------------------------
C---- Command to generate panel geometry from buffer (REMOVE LATER)
C
C
C--------------------------------------------------------------------
C---- Regenerate paneled geometry and grid
      ELSEIF(COMAND.EQ.'GENG') THEN
        CALL GENGEOM
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'NDOF') THEN
        IF(NINPUT.GE.1) THEN
         IEL = IINPUT(1)
        ELSE
         CALL ASKI('Enter element to toggle^',IEL)
        ENDIF
C
        IF(IEL.GE.1 .AND. IEL.LE.NEL) THEN
         LNDOF(IEL) = .NOT.LNDOF(IEL)
         IF(LNDOF(IEL)) THEN
          WRITE(*,*) 'Zero normal velocity enforced'
         ELSE
          WRITE(*,*) 'TE gamma regularity enforced'
         ENDIF
        ENDIF
C
        LSYSP = .FALSE.
        LGSYS = .FALSE.
        LGAMU = .FALSE.
        LGAMA = .FALSE.
C
        CALL GAMSOLV
C
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ETYP') THEN
        IEL = 0
        CALL ELTYPE(IEL)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'TPAN') THEN
        IF(NINPUT.EQ.0) THEN
         IEL = 0
        ELSE
         IEL = IINPUT(1)
        ENDIF
        CALL ELTPAN(IEL)
      ELSEIF(COMAND.EQ.'SING') THEN
        ERROR = .TRUE.
        DO IEL = 1, NEL
          IF(NETYPE(IEL).NE.0) ERROR = .FALSE.
        ENDDO
        IF(ERROR) THEN
         WRITE(*,*) 'No prescribable elements present'
         GO TO 500
        ENDIF
C
        IF(NINPUT.GE.3) THEN
         IEL = INT(RINPUT(1)+0.01)
         IF(IEL.LT.1 .OR. IEL.GT.NEL) THEN
          WRITE(*,*) 'Element index out of range'
         ELSEIF(NETYPE(IEL).EQ.7) THEN
cc         ELSEIF(NETYPE(IEL).EQ.0) THEN
          WRITE(*,*) 'Not a prescribable element'
          GO TO 500
         ELSE
          GAMSET(IEL) = RINPUT(2)
          SIGSET(IEL) = RINPUT(3)
         ENDIF
C
        ELSE
         CALL GETGSP
C
        ENDIF
C
        LGAMA = .FALSE.
cc        CALL GAMSOLV
C
        COMOLD = COMAND
        ARGOLD = COMARG
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'FF') THEN
         CALL NFCALC
         CALL FFCALC
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'BLOW') THEN
C---- Specify surface blowing in region on element
        IF(NINPUT.GE.1) THEN
         IEL = INT(RINPUT(1)+0.01)
         IF(IEL.LT.1 .OR. IEL.GT.NEL) THEN
          WRITE(*,*) 'Element index out of range'
          GO TO 500
         ELSEIF(NETYPE(IEL).NE.0) THEN
          WRITE(*,*) 'Not a prescribable element'
          GO TO 500
         ENDIF
        ENDIF
C
        IF(NINPUT.GE.4) THEN
         X1 = RINPUT(2)
         X2 = RINPUT(3)
         SRC = RINPUT(4)
        ELSE
          WRITE(*,*) 'No blowing X1,X2,VN prescribed'
          GO TO 500
        ENDIF
C
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        DO IP = IP1, IP2
          SIGVSP(IP) = 0.0
          IF(XP(IP).GE.X1 .AND. XP(IP).LE.X2) THEN
           SIGVSP(IP) = SRC
          ENDIF
        END DO
C
        LGAMA = .FALSE.
C
        COMOLD = COMAND
        ARGOLD = COMARG
        COMAND = COMPLT
        COMARG = ARGPLT
        GO TO 505
C
C--------------------------------------------------------------------
C  Relatively incomprehensible commands, used for debugging
C--------------------------------------------------------------------
C
C---- Solver hacks
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'GSOL') THEN
        IF(NINPUT.GE.1) THEN
         ITR = IINPUT(1)
        ELSE
         ITR = 1
        ENDIF
        DO I = 1, ITR
C---- Set rotor wake gamma
          CALL GTHCALC(GTH)
C---- Solution for current flow condition
          CALL GAMSOLV
        END DO
c        ITYP = 1
c        CALL TQCALC(ITYP)
c        CALL ROTRPRT(6)
C
      ELSEIF(COMAND.EQ.'GTHI') THEN
        IF(NINPUT.GE.1) THEN
         VAIN = RINPUT(1)
        ELSE
         CALL ASKR('Enter average Vaxial^',VAIN)
        ENDIF
        CALL VMAVGINIT(VAIN)
C---- Set rotor wake gamma
        CALL GTHCALC(GTH)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'RES ') THEN
c        CALL QAIC(.FALSE.)
c        CALL SYSP
        CALL GSYS(.TRUE.,.TRUE.,.TRUE., 1,NP, 1,0)
          
        do k = 1, nsys
          write(*,*) k, res(k)
        enddo
C
C---- Velocity updates
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'VMAI') THEN
        CALL VMAVGINIT(VAAVG)
C
C---- Grid updates
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ING') THEN
       CALL INIGRD
       CALL UPDROTWAK
C
      ELSEIF(COMAND.EQ.'UPG') THEN
       CALL UPDGRD
cc       CALL RLXGRD

      ELSEIF(COMAND.EQ.'RLXG') THEN
       CALL RLXGRD2
C
      ELSEIF(COMAND.EQ.'GETG') THEN
       CALL GETGRDFLW
C
C---- Solution updates and recalculation
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'QUP') THEN
       CALL ASKL('Call QAIC^',LANS)
       IF(LANS) CALL QAIC(.FALSE.)
C
       CALL ASKL('Call GAMSOLV^',LANS)
       IF(LANS) THEN
         LGAMU = .FALSE.
         CALL GAMSOLV
       ENDIF
C
       CALL ASKL('Call GUSUM^',LANS)
       IF(LANS) CALL GUSUM
C
       CALL ASKL('Call GCSUM^',LANS)
       IF(LANS) CALL QCSUM
C
       CALL ASKL('Call GCPFOR^',LANS)
       IF(LANS) CALL QCPFOR
C
       COMAND = COMPLT
       COMARG = ARGPLT
       GO TO 505
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'PSIN') THEN
        CALL GAMSOL
c        IF(.NOT.LSIGP) CALL SYSPS
c        IF(.NOT.LSIGM) CALL GSCALC
C
        WRITE(*,*)
        WRITE(*,*) '  n    ip1   ip2'
        DO IEL = 1, NEL
          WRITE(*,8100) IEL, IPFRST(IEL), IPLAST(IEL)
 8100     FORMAT(1X, I3, 2I6)
        ENDDO
C
 81     WRITE(*,*)
        WRITE(*,*) 'Enter ip, sigma'
        READ (*,*,ERR=81) IP, SIGIP
C
        IF(IP.EQ.0) THEN
         CALL QCSUM
C
         CALL QQCALC(NCTOT,ANC, IPCO,IPCP, GAM,SIG,GTH, QC, QCL,QCR )
C         CALL CPCALC(NCTOT, QCL, QINF,QREF, CPL, CPL_QCL,CPL_QINF)
         CALL CPCALC(NCTOT, QCR, QINF,QREF, CPR, CPR_QCR,CPR_QINF)
C
         COSA = 1.0
         SINA = 0.0
C
         CX(0) = 0.
         CY(0) = 0.
         CM(0) = 0.
         CD(0) = 0.
         DO IEL = 1, NEL
           I = ICFRST(IEL)
           N = ICLAST(IEL) - ICFRST(IEL) + 1
           IF(N.GE.1) THEN
            CALL FCOEFF(N,CPL(I),CPR(I),XC(I),YC(I),ANC(1,I),DSC(I), 
     &                  CX(IEL),CY(IEL),CM(IEL))
            CD(IEL) = CX(IEL)*COSA + CY(IEL)*SINA
C
            CX(0) = CX(0) + CX(IEL)
            CY(0) = CY(0) + CY(IEL)
            CM(0) = CM(0) + CM(IEL)
            CD(0) = CD(0) + CD(IEL)
           ENDIF
         ENDDO
C
         CALL CPPLOT(LPLTEL,ICPQXY)
        ELSE
C
c         SIG(IP) = SIGIP
c         SIGSP(IP) = SIGIP
c         LINF = LINFSIG(IP)
c         DO JP = 1, NPTOT
c           JSYS = JSYSGAM(JP)
c           GAM(JP) = GAM(JP) + GAM_SIG(JSYS,LINF)*SIGIP
c         ENDDO
         GO TO 81
        ENDIF
C
C--------------------------------------------------------------------
      ELSE
        WRITE(*,1010) COMAND

      ENDIF
C
      GO TO 500
C...................................................................
C
 1000 FORMAT(A)
 1010 FORMAT(1X,A4,' command not recognized.  Type a "?" for list')
C
 1200 FORMAT(/' Sonic Cp =', F10.2, '      Sonic Q/Qinf =', F10.3/)
 1300 FORMAT(/' Enter ', A,'  min, max, increment: ', $)
      END ! OPER



      SUBROUTINE GETGSP
C-------------------------------------------------------
C     Obtains from the user the prescribed singularity 
C     strengths for each element.
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*80 LINE
      DIMENSION RINPUT(20)
      LOGICAL ERROR
C
 1000 FORMAT(A)
C
 10   WRITE(*,1050)
 1050 FORMAT(
     &  ' _________________________________________________'
     &//'                       n          strengths       '
     & /'                      --  ------------------------')
C
      DO IEL = 1, NEL
        IP = IPFRST(IEL)
C
        IF    (NETYPE(IEL).EQ.1) THEN
         WRITE(*,1200) 'Sheet  gamma sigma:',
     &     IEL,GAMSET(IEL),SIGSET(IEL), XP(IP),YP(IP)
C
        ELSEIF(NETYPE(IEL).EQ.2) THEN
         WRITE(*,1200) 'Axis   doubl sourc:',
     &     IEL,GAMSET(IEL),SIGSET(IEL), XP(IP),YP(IP)
C
        ELSEIF(NETYPE(IEL).EQ.3) THEN
          WRITE(*,1200) 'Ring   Gamma Sigma:',
     &     IEL,GAMSET(IEL),SIGSET(IEL), XP(IP),YP(IP)
C
        ELSEIF(NETYPE(IEL).EQ.4) THEN
         WRITE(*,1200) 'Point  Doubl Sourc:',
     &     IEL,GAMSET(IEL),SIGSET(IEL), XP(IP),YP(IP)
C
        ELSEIF(NETYPE(IEL).EQ.5) THEN
         WRITE(*,1200) 'Source line  Gamma Sourc:',
     &     IEL,GAMSET(IEL),SIGSET(IEL), XP(IP),YP(IP)
C
        ELSEIF(NETYPE(IEL).EQ.6) THEN
         WRITE(*,1200) 'Rotor sources  Gamma Sourc:',
     &     IEL,GAMSET(IEL),SIGSET(IEL), XP(IP),YP(IP)
C
        ENDIF
      ENDDO
 1200 FORMAT( 1X, A, I4, 2G13.5,'     at  x,y = (', 2G12.4,' )' )
C
      WRITE(*,*)
 20   WRITE(*,2000)
 2000 FORMAT( ' Enter  n, strengths:  ', $)
      READ (*,1000) LINE
      IF(INDEX(LINE,'?').NE.0) GO TO 10
C
      NINPUT = 20
      CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
      IF(ERROR) GO TO 20
C
      IF(NINPUT.LT.3) THEN
       RETURN
C
      ELSE
       IEL = INT( RINPUT(1) + 0.001 )
       IF(IEL.LT.1 .OR. IEL.GT.NEL) THEN
        WRITE(*,*) 'Index out of range'
        GO TO 20
       ELSEIF(NETYPE(IEL).EQ.0) THEN
        WRITE(*,*) 'Not a prescribable element'
        GO TO 20
       ENDIF
C
      ENDIF
C
C---- OK, we have a valid element index...
      GAMSET(IEL) = RINPUT(2)
      SIGSET(IEL) = RINPUT(3)
C
      GO TO 20
      END ! GETGSP



      SUBROUTINE GETSEQ
C-------------------------------------------------------
C     Obtains from the user a sequence of prescribed 
C     singularity strengths for each element.
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*80 LINE
      DIMENSION RINPUT(20)
      DIMENSION IINPUT(20)
      LOGICAL ERROR
C
 1000 FORMAT(A)
C
c      IF(NPSEQ.LE.0) THEN
       CALL ASKI('Enter number of operating points^',NPSEQ)
       IF(NPSEQ.LT.2) THEN
        WRITE(*,*) 'Need at least two points for a sequence'
        RETURN
       ENDIF
c      ENDIF
C
      IF(QINFP1.EQ.0.0 .AND. QINFP2.EQ.0.0) THEN
       QINFP1 = QINF
       QINFP2 = QINF
      ENDIF
C
 5    CONTINUE
      QINFPD = (QINFP2-QINFP1)/FLOAT(NPSEQ-1)
      DO IEL = 1, NEL
        GAMPD(IEL) = (GAMP2(IEL)-GAMP1(IEL))/FLOAT(NPSEQ-1)
        SIGPD(IEL) = (SIGP2(IEL)-SIGP1(IEL))/FLOAT(NPSEQ-1)
      ENDDO
C
C
 10   WRITE(*,1010) NPSEQ
 1010 FORMAT(
     &  ' ______________________',
     &  '_____________________________________________'
     &//'         number of points:   N', I5
     &//'                       ',
     &  '      #     first       last          delta'
     & /'                             ',
     &  '   --  ---------- -----------   -----------')
C
CCC      'Sheet  gam1 gam2 dgam:
CCC        II123456789012123456789012  123456789012
C
      WRITE(*,1100) 'Qinf       Qinf1 Qinf2 dQinf:   Q',
     &     QINFP1,QINFP2,QINFPD
C
      DO IEL = 1, NEL
        IP = IPFRST(IEL)
C
        IF    (NETYPE(IEL).EQ.1) THEN
          WRITE(*,1200) 'Sheet      gam1  gam2  dgam :',
     &     IEL,GAMP1(IEL),GAMP2(IEL),GAMPD(IEL)
          WRITE(*,1200) '           sig1  sig2  dsig :',
     &    -IEL,SIGP1(IEL),SIGP2(IEL),SIGPD(IEL)
C
        ELSEIF(NETYPE(IEL).EQ.2) THEN
          WRITE(*,1200) 'Axis       dbl1  dbl2  ddbl :',
     &     IEL,GAMP1(IEL),GAMP2(IEL),GAMPD(IEL)
          WRITE(*,1200) '           src1  src2  dsrc :',
     &    -IEL,SIGP1(IEL),SIGP2(IEL),SIGPD(IEL)
C
        ELSEIF(NETYPE(IEL).EQ.3) THEN
          WRITE(*,1200) 'Ring       Gam1  Gam2  dGam :',
     &     IEL,GAMP1(IEL),GAMP2(IEL),GAMPD(IEL)
          WRITE(*,1200) '           Sig1  Sig2  dSig :',
     &    -IEL,SIGP1(IEL),SIGP2(IEL),SIGPD(IEL)
C
        ELSEIF(NETYPE(IEL).EQ.4) THEN
          WRITE(*,1200) 'Point      Dbl1  Dbl2  dDbl :',
     &     IEL,GAMP1(IEL),GAMP2(IEL),GAMPD(IEL)
          WRITE(*,1200) '           Src1  Src2  dSrc :',
     &    -IEL,SIGP1(IEL),SIGP2(IEL),SIGPD(IEL)
C
        ELSEIF(NETYPE(IEL).EQ.5) THEN
          WRITE(*,1200) 'Src line   Sig1  Sig2  dSig :',
     &    -IEL,SIGP1(IEL),SIGP2(IEL),SIGPD(IEL)
C
        ELSEIF(NETYPE(IEL).EQ.6) THEN
          WRITE(*,1200) 'Rotor src  Sig1  Sig2  dSig :',
     &    -IEL,SIGP1(IEL),SIGP2(IEL),SIGPD(IEL)
        ENDIF
      ENDDO
 1100 FORMAT( 1X, A,     2G12.4, 2X, G12.4)
 1200 FORMAT( 1X, A, I4, 2G12.4, 2X, G12.4)
C
      WRITE(*,*)
 20   WRITE(*,2000)
 2000 FORMAT( ' Enter  # (or Q), first, last:  ', $)
      READ (*,1000) LINE
C
      DO K = 1, 80
        IF(LINE(K:K).NE.' ') THEN
         LINE = LINE(K:80)
         GO TO 25
        ENDIF
      ENDDO
      RETURN
C
 25   CONTINUE
C
C---- if there's a question mark in the line, clear it and reprint menu later
      IQUEST = INDEX(LINE,'?')
      IF(IQUEST.GT.0) LINE(IQUEST:IQUEST) = ' '
C
C---- check for special cases...
      IF    (INDEX('Nn',LINE(1:1)).NE.0) THEN
        LINE(1:1) = ' '
        NINPUT = 20
        CALL GETINT(LINE,IINPUT,NINPUT,ERROR)
        IF(ERROR .OR. NINPUT.LT.1) GO TO 20
C
        NPSEQ = IINPUT(1)
        IF(NPSEQ.LT.1) THEN
         RETURN
        ELSE
C------- go recalculate increments
         GO TO 5
        ENDIF
C
      ELSEIF(INDEX('Qq',LINE(1:1)).NE.0) THEN
        LINE(1:1) = ' '
        NINPUT = 20
        CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
        IF(ERROR .OR. NINPUT.LT.2) GO TO 20
C
        QINFP1 = RINPUT(1)
        QINFP2 = RINPUT(2)
        QINFPD = (QINFP2-QINFP1)/FLOAT(NPSEQ-1)
C
      ELSEIF(INDEX('Ss',LINE(1:1)).NE.0) THEN
        LINE(1:1) = ' '
        NINPUT = 20
        CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
        IF(ERROR .OR. NINPUT.LT.2) GO TO 20
C
        SIGP1(IEL) = RINPUT(1)
        SIGP2(IEL) = RINPUT(2)
        SIGPD(IEL) = (SIGP2(IEL)-SIGP1(IEL))/FLOAT(NPSEQ-1)
C
      ELSE
        NINPUT = 20
        CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
        IF(ERROR) GO TO 20
        IF(NINPUT.LT.3) GO TO 90
C
        IELS = INT( RINPUT(1) + SIGN(0.001,RINPUT(1)) )
        IEL = IABS(IELS)
        IF(IEL.EQ.0 .OR. IEL.GT.NEL) THEN
         WRITE(*,*) 'Index out of range'
         GO TO 20
        ELSEIF(NETYPE(IEL).EQ.0) THEN
         WRITE(*,*) 'Not a prescribable element'
         GO TO 20
        ENDIF
C
C------ OK, we have a valid element index IEL...
        IF(IELS.GT.0) THEN
         GAMP1(IEL) = RINPUT(2)
         GAMP2(IEL) = RINPUT(3)
         GAMPD(IEL) = (GAMP2(IEL)-GAMP1(IEL))/FLOAT(NPSEQ-1)
C
        ELSE
         SIGP1(IEL) = RINPUT(2)
         SIGP2(IEL) = RINPUT(3)
         SIGPD(IEL) = (SIGP2(IEL)-SIGP1(IEL))/FLOAT(NPSEQ-1)
        ENDIF
      ENDIF
C
 90   CONTINUE
      IF(IQUEST.GT.0) THEN
       GO TO 10
      ELSE
       GO TO 20
      ENDIF
      END ! GETSEQ



      SUBROUTINE PLTOGG(LTMOD)
      INCLUDE 'DFDC.INC'
      LOGICAL LTMOD
C
      CHARACTER*1 ANS
C
      LTMOD = .FALSE.
C
 1000 FORMAT(A)
C
 10   WRITE(*,1010) LPGRID, LPSTAG, LSPLT, LVPLT, LEPLT, 
     &              LPREF, LFREF
 1010 FORMAT(/'     ----------------------------'
     &       /'      G rid on Cp plots   ',L4
     &       /'      S tagnation pt plot ',L4
     &       /'      P oints on line plot',L4
     &       /'      V ertices on panels ',L4
     &       /'      N umbers of elements',L4
     &       /'      C p data overlay    ',L4
     &       /'      F orce data overlay ',L4 )
      WRITE(*,2100)
 2100 FORMAT(/'     Toggle: ',$)
      READ(*,1000) ANS
C
      IF(ANS.EQ.' ') THEN
       RETURN
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Gg',ANS).NE.0) THEN
       LPGRID = .NOT.LPGRID
       LTMOD = .TRUE.
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Ss',ANS).NE.0) THEN
       LPSTAG = .NOT.LPSTAG
       LTMOD = .TRUE.
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Pp',ANS).NE.0) THEN
       LSPLT = .NOT.LSPLT
       LTMOD = .TRUE.
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Vv',ANS).NE.0) THEN
       LVPLT = .NOT.LVPLT
       LTMOD = .TRUE.
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Nn',ANS).NE.0) THEN
       LEPLT = .NOT.LEPLT
       LTMOD = .TRUE.
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Cc',ANS).NE.0) THEN
       LPREF = .NOT.LPREF
       LTMOD = .TRUE.
C
C--------------------------------------------------------------------
      ELSEIF(INDEX('Ff',ANS).NE.0) THEN
       LFREF = .NOT.LFREF
       LTMOD = .TRUE.
C
      ELSE
       WRITE(*,*) '* Unrecognized command'
C
      ENDIF
C
      GO TO 10
      END



      SUBROUTINE CPDUMP(FNAME1)
      INCLUDE 'DFDC.INC'
      CHARACTER*(*) FNAME1
C
      CHARACTER*80 FNAME, FILDEF
C
 1000 FORMAT(A)
C
      IF(FNAME1(1:1).NE.' ') THEN
       FNAME = FNAME1
      ELSE
C----- no argument... get it somehow
       IF(NPREFIX.GT.0) THEN
C------ offer default using existing prefix
        FILDEF = PREFIX(1:NPREFIX) // '.cp'
        WRITE(*,1100) FILDEF
 1100   FORMAT(/' Enter filename:  ', A)
        READ(*,1000) FNAME
        CALL STRIP(FNAME,NFN)
        IF(NFN.EQ.0) FNAME = FILDEF
       ELSE
C------ nothing available... just ask for filename
        CALL ASKS('Enter filename^',FNAME)
       ENDIF
      ENDIF
C
C
      LU = 19
      OPEN(LU,FILE=FNAME,STATUS='UNKNOWN')
      REWIND(LU)
C
      WRITE(LU,1000)
     & '#     x          y         Cp         Q  '
C         0.23451  0.23451  0.23451  0.23451
c             123456789012345678901234567890
C
ccc      CALL COMSET
C
ccc      BETA = SQRT(1.0 - MINF**2)
ccc      BFAC = 0.5*MINF**2 / (1.0 + BETA)
cc      BETA = 1.0
cc      BFAC = 0.0
C
      DO 10 IEL=1, NEL
        IF(ICFRST(IEL).LE.0) GO TO 10
C
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        N = IC2 - IC1 + 1
C
        IF(ISPLOT(IEL).EQ. 0 .OR.
     &     ISPLOT(IEL).EQ. 1      ) THEN
C-------- dump Cp on right side
          WRITE(LU,1000) ' '
          WRITE(LU,1000) ' '
          DO IC = IC1, IC2
            QMAG = SQRT(QCR(1,IC)**2 + QCR(2,IC)**2) 
            WRITE(LU,8500) XC(IC), YC(IC), CPR(IC), QMAG
          ENDDO
        ENDIF

        IF(ISPLOT(IEL).EQ. 0 .OR.
     &     ISPLOT(IEL).EQ. 2      ) THEN
C-------- dump Cp on left side
          WRITE(LU,1000) ' '
          WRITE(LU,1000) ' '
          DO IC = IC1, IC2
            QMAG = SQRT(QCL(1,IC)**2 + QCL(2,IC)**2) 
            WRITE(LU,8500) XC(IC), YC(IC), CPL(IC), QMAG
          ENDDO
        ENDIF
 10   CONTINUE
C
 8500 FORMAT(5(1X,F10.5))
C
      CLOSE(LU)
      RETURN
      END ! CPDUMP



      SUBROUTINE OPFIL(LU,FNAME)
      CHARACTER*(*) FNAME
C
      CHARACTER*4 COMAND
      CHARACTER*128 COMARG,TMP
      CHARACTER*1 ANS, DUMMY
C
C---- get filename if it hasn't been already specified
      IF(FNAME(1:1).EQ.' ') CALL ASKS('Enter output filename^',FNAME)
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
          READ(LU,1000,END=60) DUMMY
 1000     FORMAT(A)
        ENDDO
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
      END ! OPFIL


      SUBROUTINE DRGOPER(LMODG)
C------------------------------------------------------
C     Drag area options
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
cc      INCLUDE 'PLOT.INC'
      CHARACTER*4 COMAND, COMOLD
      CHARACTER*1 CIND
C
      CHARACTER*16 PROMPT
      CHARACTER*80 LINE
      CHARACTER*128 COMARG, ARGOLD
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
      LOGICAL LMODG
C
      SAVE COMOLD, ARGOLD
C
 1000 FORMAT(A)
C
      COMAND = '****'
      COMARG = ' '
C
c      IF(NDOBJ.EQ.0) THEN
c       WRITE(*,*)
c       WRITE(*,*) '***  No drag areas defined  ***'
c       RETURN
c      ENDIF
C
      WRITE(*,*)
      CALL SHOWDRAGOBJ(0)
C
      IDRG = 0
      DCDA = 0.
      DELX = 0.
      DELY = 0.
      XFAC = 1.0
      YFAC = 1.0
      FAC  = 1.0
      LMODG = .FALSE.
C
C---- make sure we have a valid target element
      IF(IDRG.LT.1 .OR. IDRG.GT.NDOBJ) THEN
       IF(NDOBJ.LE.1) THEN
        IDRG = 0
       ELSE
        WRITE(*,*) ' '
 2      CALL ASKI('Specify target drag area index^',IDRG)
        IF(IDRG.LT.0 .OR. IDRG.GT.NDOBJ) THEN
         WRITE(*,*) 'Number of drag areas present:', NDOBJ
         IDRG = 1
         GO TO 2
        ENDIF
       ENDIF
      ENDIF

C---- go plot geometry with scale,offset initialization
c      LOFINI = .TRUE.
c      COMAND = 'PLOT'
c      GO TO 510
CC
C---- begin user-interaction loop
C................................................................
C
 500  CONTINUE
      COMOLD = COMAND
      ARGOLD = COMARG
C
 501  CONTINUE
      PROMPT = '.DRAG ^'
C
      IF(IDRG.GT.0 .AND. NDOBJ.GT.1) THEN
C------ add current target drag area to prompt
        IF    (IDRG.LT. 10) THEN
         PROMPT(7:8) = '(' // PNUM(IDRG)(2:2)
         K = 9
        ELSEIF(IDRG.LT.100) THEN
         PROMPT(7:9) = '(' // PNUM(IDRG   )
         K = 10
        ENDIF
        KEL = IDRG
        PROMPT(K:K+1) = ')^'
      ENDIF
C
      CALL ASKC(PROMPT,COMAND,COMARG)
C
C---- Quit with blank line at prompt
      IF(COMAND.EQ.'    ') THEN
c       IF(LPLOT) THEN
c        CALL PLEND
c        CALL CLRZOOM
c       ENDIF
c       LPLOT = .FALSE.
       RETURN
      ENDIF
C
 505  CONTINUE
C
C---- was only a new target-element index typed?
      CIND = COMAND(1:1)
      IF(INDEX('0123456789',CIND) .NE. 0) THEN
       READ(COMAND,*,ERR=508) ID
       IF(ID.LT.0 .OR. ID.GT.NDOBJ) GO TO 508
       IDRG = ID
       GO TO 501
      ENDIF
 508  CONTINUE
C
C---- extract numbers (if any) from command argument string
      DO I=1, 20
        IINPUT(I) = 0
        RINPUT(I) = 0.0
      ENDDO
      NINPUT = 0
      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
      NINPUT = 0
      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C
C---- can enter here if COMAND string is set to desired command
 510  CONTINUE
C--------------------------------------------------------------------
      IF(    COMAND.EQ.'HELP'
     &  .OR. COMAND.EQ.'?   ') THEN
        WRITE(*,*) ' '
        IF(NDOBJ.GT.1) WRITE(*,1104) IDRG
 1104   FORMAT(
     &  '        i   Change target drag area (currently =', I3, ')')
C
        WRITE(*,1106) 
 1106  FORMAT(' ...............................................'
     & /'   List     List drag areas'
     & /'   Togg     Toggle drag areas on/off'
     & /'   New      New drag area'
     & /'   Edit     Edit x,y,CDA values'
     & /'   DXY  rr  Add delta dx,dy to drag area x,y'
     & /'   SXY  rr  Scale drag area x,y'
     & /'   DCD   r  Add delta CDA to drag area CDA'
     & /'   SCD   r  Scale drag area CDA values')
C
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'LIST' 
     &  .OR. COMAND.EQ.'L') THEN
       CALL SHOWDRAGOBJ(0)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'TOGG' 
     &  .OR. COMAND.EQ.'T') THEN
C---- toggle drag area objects flag
       LDRGOBJ = .NOT.LDRGOBJ
       IF(LDRGOBJ) THEN
         WRITE(*,*) 'Drag areas are active'
       ELSE
         WRITE(*,*) 'Drag areas de-activated'
       ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'NEW' 
     &  .OR. COMAND.EQ.'N') THEN
C
       IF(NDOBJ.GE.NDRGX) THEN
         WRITE(*,*) '# drag objects exceeds dimension ',NDRGX
         GO TO 500
       ENDIF
C
       ND = NDOBJ + 1
       WRITE(*,1108)
 1108  FORMAT(/' Input x, y, CDA points (terminate with blank line)')
       ID = 0
       DO I = 1, IDX
 11       READ(*,1000,ERR=500,END=500) LINE
         IF(LINE.EQ.' ') GO TO 20
         READ(LINE,*,ERR=14,END=14) XX,YY,CDAIN
         ID = ID + 1
         XDDEF(ID,ND) = XX
         YDDEF(ID,ND) = YY
         CDADEF(ID,ND) = CDAIN
         GO TO 15
C
   14    WRITE(*,*) 'Try again'
         GO TO 11
   15    CONTINUE
       END DO
C
 20    IF(ID.GE.2) THEN
         NDDEF(ND) = ID
         NDOBJ = ND
         LMODG = .TRUE.
         WRITE(*,*) 'New drag area defined ',NDOBJ
         IDRG = ND
       ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'EDIT' 
     &  .OR. COMAND.EQ.'E') THEN
        IF(IDRG.LE.0) THEN
         WRITE(*,*) 'Specify drag object # first...'
         GO TO 500
        ENDIF
C
       CALL SHOWDRAGOBJ(IDRG)
 22    WRITE(*,1112)
 1112  FORMAT(/' Input drag area pt index to edit (blank to end)')
       READ(*,1000,ERR=500,END=500) LINE
       IF(LINE.EQ.' ') GO TO 30
       READ(LINE,*,ERR=22) ID
       IF(ID.LT.0 .OR. ID.GT.NDDEF(IDRG)) GO TO 22
C
 24    WRITE(*,1114)       
 1114  FORMAT(/' Input new x, y, CDA (skip with blank line)')
       READ(*,1000,ERR=24,END=22) LINE
       IF(LINE.NE.' ') THEN
        READ(LINE,*,ERR=24,END=24) XX,YY,CDAIN
        XDDEF(ID,IDRG)  = XX
        YDDEF(ID,IDRG)  = YY
        CDADEF(ID,IDRG) = CDAIN
        LMODG = .TRUE.
        GO TO 22
       ENDIF
C
 30    GO TO 500
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'DXY') THEN
        IF(IDRG.LE.0) THEN
         WRITE(*,*) 'Specify drag object # first...'
         GO TO 500
        ENDIF
        IF    (NINPUT.LE.0) THEN
         RINPUT(1) = DELX
         RINPUT(2) = DELY
         CALL ASKRN('Enter delta(x), delta(y)^',RINPUT,2)
        ELSEIF(NINPUT.LE.1) THEN
         RINPUT(2) = DELY
         CALL ASKR('Enter delta(y)^',RINPUT(2))
        ENDIF
        DELX = RINPUT(1)
        DELY = RINPUT(2)
C
        ND = IDRG
        DO ID = 1, NDDEF(ND)
          XDDEF(ID,ND) = XDDEF(ID,ND) + DELX
          YDDEF(ID,ND) = YDDEF(ID,ND) + DELY
        END DO
        LMODG = .TRUE.
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SXY') THEN
        IF(IDRG.LE.0) THEN
         WRITE(*,*) 'Specify drag object # first...'
         GO TO 500
        ENDIF
        IF(NINPUT.GE.1) THEN
         FAC = RINPUT(1)
         XFAC = FAC
         YFAC = FAC
        ELSE
        CALL ASKR('Enter scale factor (0 for separate x,y scales)^',FAC)
         XFAC = FAC
         YFAC = FAC
        ENDIF
C
        IF(FAC .EQ. 0.0) THEN
         IF(NINPUT.GE.3) THEN
          XFAC = RINPUT(2)
          YFAC = RINPUT(3)
         ELSE
          CALL ASKR('Enter x scale factor^',XFAC)
          CALL ASKR('Enter y scale factor^',YFAC)
         ENDIF
        ENDIF
C
        ND = IDRG
        DO ID = 1, NDDEF(ND)
          XDDEF(ID,ND) = XDDEF(ID,ND) * XFAC
          YDDEF(ID,ND) = YDDEF(ID,ND) * YFAC
        END DO
        LMODG = .TRUE.
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'DCD') THEN
        IF(IDRG.LE.0) THEN
         WRITE(*,*) 'Specify drag object # first...'
         GO TO 500
        ENDIF
        IF(NINPUT.GE.1) THEN
         DCDA = RINPUT(1)
        ELSE
         CALL ASKR('Enter CDA increment^',DCDA)
        ENDIF
C
        ND = IDRG
        DO ID = 1, NDDEF(ND)
          CDADEF(ID,ND) = CDADEF(ID,ND) + DCDA
        END DO
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SCD') THEN
        IF(IDRG.LE.0) THEN
         WRITE(*,*) 'Specify drag object # first...'
         GO TO 500
        ENDIF
        IF(NINPUT.GE.1) THEN
         FAC = RINPUT(1)
        ELSE
         CALL ASKR('Enter CDA scale factor^',FAC)
        ENDIF
C
        ND = IDRG
        DO ID = 1, NDDEF(ND)
          CDADEF(ID,ND) = CDADEF(ID,ND) * FAC
        END DO
C
      ELSE
       WRITE(*,1050) COMAND
C
      ENDIF
C
      GO TO 500
C
C....................................................
 1050 FORMAT(' Command ',A4,' not recognized.  Type a " ? " for list.')
      END


