C***********************************************************************
C    MODIDF
C    Blade geometry modification in DFDC
C    August 2013, Esotec Developments
C
C    Some code from module xmodi.f - Xrotor 7.55 
C    Copyright (C) 2011 Mark Drela 
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


      SUBROUTINE DFMODI
C---------------------------------------------------
C     Rotor and stator geometry modification in dfdc
C---------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      CHARACTER*1 CHKEY, ANS, CIND
      CHARACTER*16 PROMPT
      CHARACTER*8 PROMPT2
      CHARACTER*4 COMAND,COMOLD
      CHARACTER*128 COMARG, ARGOLD, ARGPLT, FNAME
C
C---- local arrays for calling GETINT,GETFLT...
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
C
C--- temporary arrays
      DIMENSION T1(IRX), T1S(IRX),
     &          T2(IRX), T2S(IRX),
     &          T3(IRX), T3S(IRX)
C
c      PLFACD = 0.6
c      PLFAC1 = 0.7
c      PLFAC2 = 0.85
C
      PLFAC1 = 0.6
      PLFAC2 = 0.8
      PLFACD = 0.6
      XORG = 0.15
      YORG = 0.10
C
      CFACA = 1.0
      CFACB = 0.0
      TDEG  = 0.0
C
C---- "last" command (nothing yet)
      COMAND = '****'
      COMARG = ' '
      FNAME  = ' '
      ISPARS = ' '
C
      COMOLD = COMAND
      ARGOLD = COMARG
C
C---- start with disk 1
      NDSK = 1
      N = NDSK
C
C---- check if geometry loaded, process if required
C
      IF(LLOAD) THEN
        IF(NPTOT.EQ.0) CALL GENGEOM
        IF(.NOT.LBLDEF) THEN
          WRITE(*,*)
          WRITE(*,*) 'No bladed disks defined in this case'
          RETURN
        ENDIF
      ELSE
        WRITE(*,*)
        WRITE(*,*)'No geometry loaded'
        RETURN
      ENDIF
C
C---- load arrays for bladed disks
C
      DO NR=1,NROTOR
        IF(IRTYPE(NR).NE.2) THEN
          WRITE(*,1050) NR
        ELSE
          CALL GETBLADE(NR)
        ENDIF
      ENDDO
C
      WRITE(*,*)
      WRITE(*,*) 'You are working with the buffer blade'
C
C-------------------------------------------------------------
C---- begin user-interaction loop
C............................................................
C
 500  CONTINUE
      PROMPT = '.MODI^'
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
          WRITE(*,*) 'Previous .MODI command not valid'
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
           WRITE(*,1022) ID
        ELSEIF(IRTYPE(ID).EQ.2.AND.OMEGA(ID).EQ.0.) THEN
           WRITE(*,1024) ID
        ELSE
           WRITE(*,1026) ID
           GO TO 500
        ENDIF
C
        NDSK = ID
        N = NDSK
        GO TO 500
      ENDIF
C
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
 1010 FORMAT(1X,A4,' command not recognized.  Type "?" for list')
 1020 FORMAT(/,1X, 'Disk',I2,' is an actuator disk')
 1022 FORMAT(/,1X, 'Disk',I2,' is a bladed rotor')
 1024 FORMAT(/,1X, 'Disk',I2,' is a bladed stator')
 1026 FORMAT(/,1X, 'Disk',I2,' is undefined')
 1050 FORMAT(/,1x, 'Disk',I2,' is not a bladed rotor')
C
C---- can jump in here if COMAND string has been set elsewhere
 510  CONTINUE
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
      GOTO 500
      ENDIF
C
 1100 FORMAT(
     &  /'   <i>     Select disk index')
C
 1110 FORMAT(
     &  /'   GSET    Set buffer  geometry <== current geometry'
     &  /'   EXEC    Set current geometry <== buffer  geometry'
C
     & //'   MODC    Modify blade chord distribution'
     &  /'   MODB    Modify blade twist distribution'
     &  /'   SCAL rr Scale blade chords'
     &  /'   TLIN r  Add linear blade twist (proportional to r/R)'
     &  /'   INTE    Interpolate between current rotors'
C
     & //'   CURR    Plot geometry - current blade'
     &  /'   BUFF    Plot geometry - buffer blade'
     &  /'   CHOR    Plot chord - current and buffer blades'
     &  /'   BETA    Plot beta  - current and buffer blades'
     &  /'   ANNO    Annotate current plot'
     &  /'   HARD    Hardcopy current plot'
     &  /'   SIZE r  Change plot-object size')
C
C
      IF(COMAND.EQ.'EXEC'.OR.COMAND.EQ.'E   ') GO TO 26
      IF(COMAND.EQ.'GSET') GO TO 28
C
      IF(COMAND.EQ.'MODB') GO TO 30
      IF(COMAND.EQ.'MODC') GO TO 32
      IF(COMAND.EQ.'SCAL') GO TO 34
      IF(COMAND.EQ.'TLIN') GO TO 36
      IF(COMAND.EQ.'INTE') GO TO 38
C
      IF(COMAND.EQ.'CHOR') GO TO 40
      IF(COMAND.EQ.'BETA') GO TO 50
      IF(COMAND.EQ.'BUFF') GO TO 60
      IF(COMAND.EQ.'CURR') GO TO 70
C
      IF(COMAND.EQ.'ANNO') GO TO 92
      IF(COMAND.EQ.'HARD') GO TO 94
      IF(COMAND.EQ.'SIZE') GO TO 96
C
      WRITE(*,1010) COMAND
      GO TO 500
C
C---------------------------------------------------------------------
 26   CONTINUE
      CALL PUTBLADE(NDSK)
      GO TO 500
C---------------------------------------------------------------------
 28   CONTINUE
      CALL GETBLADE(NDSK)
      WRITE(*,*)
      WRITE(*,*) 'Current blade written to buffer blade'
      GO TO 500
C---------------------------------------------------------------------
 30   CONTINUE
      CALL MODBE(NDSK)
      GO TO 500
C---------------------------------------------------------------------
 32   CONTINUE
      CALL MODCH(NDSK)
      GO TO 500
C---------------------------------------------------------------------
 34   IF(NINPUT.GE.2) THEN
       CFACA = RINPUT(1)
       CFACB = RINPUT(2)
      ELSE
       WRITE(*,*) ' '
       WRITE(*,*) '(Chord)new = (Chord)old * [A + B(r/R)]'
       CALL ASKR('Enter constant chord scaling factor A^',CFACA)
       CALL ASKR('Enter linear chord scaling factor B^',CFACB)
      ENDIF
C
      DO I=1, NRP
       BUFCHD(I,N)=BUFCHD(I,N)*(CFACA+CFACB*YRP(I,N)/YRP(NRP,N))
      ENDDO
C
      GO TO 500
C
C---------------------------------------------------------------------
 36   IF(NINPUT.GE.1) THEN
       TDEG = RINPUT(1)
      ELSE
       CALL ASKR('Tip angle change (deg)^',TDEG)
      ENDIF
C
      TRAD = TDEG * DTR
      DO I=1, NRP
        IF(OMEGA(N).GT.0.0)THEN
          BUFBET(I,N) = BUFBET(I,N) + TRAD*YRP(I,N)/YRP(NRP,N)
        ELSE
          BUFBET(I,N) = BUFBET(I,N) - TRAD*YRP(I,N)/YRP(NRP,N)
        ENDIF
      ENDDO
C
      GO TO 500
C
C---------------------------------------------------------------------
 38   CONTINUE
      IF(NROTOR.LT.2) THEN
        WRITE(*,*) 'Interpolation requires two rotors'
        GO TO 500
      ENDIF
C
      DO NR = 1,2
        IF(IRTYPE(NR).NE.2 .OR. OMEGA(NR).EQ.0.0) THEN
          WRITE(*,*) 'Interpolation requires two rotors'
          GO TO 500
        ENDIF
      ENDDO
C
      IF(NINPUT.GE.1) THEN
        BLFAC = RINPUT(1)
      ELSE
        WRITE(*,*)
        WRITE(*,*)'   Disk1    0 <------> 1    Disk2'
        CALL ASKR('Enter rotor interpolation factor^',BLFAC)
      ENDIF
C
      LBA = .FALSE.
      LBC = .FALSE.
C
 39   CALL ASKS(' Blend angles, chords, or both? (a,c,b)^',CHKEY)
      WRITE(*,*)
      IF (CHKEY.EQ.' ') THEN
         WRITE(*,*) 'Interpolation aborted - geometry unchanged'
         GO TO 500
      ENDIF
C         
      IF    (CHKEY.EQ.'A' .OR. CHKEY.EQ.'a') THEN
         LBA = .TRUE.
      ELSEIF(CHKEY.EQ.'C' .OR. CHKEY.EQ.'c') THEN
         LBC = .TRUE.
      ELSEIF(CHKEY.EQ.'B' .OR. CHKEY.EQ.'b') THEN
         LBA = .TRUE.
         LBC = .TRUE.
      ELSE
         GO TO 39
      ENDIF
C                  
      CALL DFBLEND(NDSK)
      GO TO 500
C
C---- Plotting -----------------------------------------------------
C---- Chord plot
 40   CONTINUE
      CALL PLOTCH(NDSK)
      GO TO 500
C
C---- Beta plot
 50   CONTINUE
      CALL PLOTBE(NDSK)
      GO TO 500
C
C--- 3 view geometry plot of single buffer blade
 60   CONTINUE
      CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFAC1*SIZE,LPLOT,LLAND)
      CALL PLOT(XORG,YORG,-3)
      CALL BLDPLTM(N,'ALUE')
      GO TO 500
C
C--- 3 view geometry plot of single current blade
 70   CONTINUE
      CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFAC1*SIZE,LPLOT,LLAND)
      CALL PLOT(XORG,YORG,-3)
      CALL BLDPLT(N,'ALUE')
      GO TO 500
C
C---------------------------------------------------------------------
 92   IF(LPLOT) THEN
       CALL ANNOT(1.2*CSIZE)
      ELSE
       WRITE(*,*) 'No current plot'
      ENDIF
      GO TO 500
C
C---------------------------------------------------------------------
 94   IF(LPLOT) THEN
       CALL PLEND
       CALL REPLOT(IDEVPS)
      ELSE
       WRITE(*,*) 'No current plot'
      ENDIF
      GO TO 500
C
C---------------------------------------------------------------------
 96   IF(NINPUT.GE.1) THEN
       SIZE = RINPUT(1)
      ELSE
       WRITE(*,*) 'Current plot-object size =', SIZE
       CALL ASKR('Enter new plot-object size^',SIZE)
      ENDIF
      GO TO 500
C.....................................................................
C
      END ! MODIDF



      SUBROUTINE DFBLEND(NDSK)
      INCLUDE 'DFDC.INC'
      DIMENSION TBET(NRX)
C----------------------------------------------------------------
C     interpolates geometry between bladed rotors Disk1 and Disk2
C----------------------------------------------------------------
      N = NDSK
      IF(LBC) THEN
        DO I = 1,NRP
          BUFCHD(I,N)=CURCHD(I,1)+BLFAC*(CURCHD(I,2)-CURCHD(I,1))
        ENDDO
        WRITE(*,*) 'Interpolated chords written to buffer blade'
      ENDIF
C
      IF(LBA) THEN
        DO I=1,NRP
          DO NR=1,2
            IF(OMEGA(NR).GT.0.0) THEN
              TBET(NR) = CURBET(I,NR)
            ELSE
              TBET(NR) = PI-CURBET(I,NR)
            ENDIF
          ENDDO
C
          TMPBET=TBET(1)+BLFAC*(TBET(2)-TBET(1))
C
          IF(OMEGA(N).GT.0.0) THEN
            BUFBET(I,N) = TMPBET
          ELSE
            BUFBET(I,N) = PI-TMPBET
          ENDIF
        ENDDO
C
C---- rotate blades to match current beta
C
        IRPNT = 1  ! beta match point (points inboard of tip)
        NRPNT = NRP - IRPNT
C
        DBET = CURBET(NRPNT,N) - BUFBET(NRPNT,N)
        DO I=1,NRP
          BUFBET(I,N)=BUFBET(I,N) + DBET
        ENDDO
        WRITE(*,*) 'Interpolated angles written to buffer blade'
C
      ENDIF
C
      RETURN
      END  ! DFBLEND

C------------------------------------------------------------


      SUBROUTINE GETBLADE(NDSK)
      INCLUDE 'DFDC.INC'
C---- imports current blade, splines to rotor points,
C     stores chord and beta to buffer and current arrays
C
      DIMENSION T1(IRX), T1S(IRX),
     &          T2(IRX), T2S(IRX)
C
      N = NDSK
C
C---- spline blade data to rotor points
C
      DO I=1,NRC
         T1(I) = CHR  (I,N)
         T2(I) = BETAR(I,N)
      ENDDO
C
      CALL SEGSPL(T1,T1S,YRC(1,N),NRC)
      CALL SEGSPL(T2,T2S,YRC(1,N),NRC)
C
      DO I=1,NRP
         BUFCHD(I,N) = SEVAL(YRP(I,N),T1,T1S,YRC(1,N),NRC)
         BUFBET(I,N) = SEVAL(YRP(I,N),T2,T2S,YRC(1,N),NRC)
         CURCHD(I,N) = BUFCHD(I,N)
         CURBET(I,N) = BUFBET(I,N)
      ENDDO
C
      RETURN
      END  ! GETBLADE
C-------------------------------------------------------------

      SUBROUTINE PUTBLADE(NDSK)
      INCLUDE 'DFDC.INC'
      LOGICAL ERROR
C---- exports buffer rotor to dfdc centerpoints
C
      DIMENSION T1(IRX), T1S(IRX),
     &          T2(IRX), T2S(IRX)
C
      N = NDSK
C
      DO I=1,NRP
         T1(I) = BUFCHD(I,N)
         T2(I) = BUFBET(I,N)
      ENDDO

      CALL SEGSPL(T1,T1S,YRP(1,N),NRP)
      CALL SEGSPL(T2,T2S,YRP(1,N),NRP)
C
      DO I=1,NRC
         CHR  (I,N) = SEVAL(YRC(I,N),T1,T1S,YRP(1,N),NRP)
         BETAR(I,N) = SEVAL(YRC(I,N),T2,T2S,YRP(1,N),NRP)
      ENDDO
C
      DO I=1,NRP
         CURCHD(I,N) = BUFCHD(I,N)
         CURBET(I,N) = BUFBET(I,N)
      ENDDO
C
      WRITE(*,*)
      WRITE(*,*) 'Buffer blade written to current blade'
      IF(LBBLOFT(N)) CALL BBUPDATE(N)
C
      RETURN
      END  ! PUTBLADE
C ---------------------------------------------------------------


      SUBROUTINE MODBE(NDSK)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C------------------------------------------------
C     Takes user cursor input to modify 
C     buffer blade angle array.
C------------------------------------------------
      DIMENSION YLIMS(2)
      EXTERNAL PLCHAR,PLMATH
C
      N = NDSK
      RTD = 180./PI
      PLFAC = 0.95
      PLPAR = 1.20*PAR
      SH = 0.3*CSIZE
      CSL = 1.0*CSIZE
      NSPLT = 20
C
C---- work with temporary arrays
C
      IF(OMEGA(N).GT.0.0) THEN
        BETMIN = BUFBET(1,N)*RTD
        BETMAX = BUFBET(1,N)*RTD
        DO I=1, NRP
          WB1(I) = YRP(I,N)/YRP(NRP,N)
          WB2(I) = BUFBET(I,N) * RTD
          BETMIN = MIN(WB2(I),BETMIN)
          BETMAX = MAX(WB2(I),BETMAX)
        ENDDO
      ELSE
        BETMIN = (PI-BUFBET(1,N))*RTD
        BETMAX = (PI-BUFBET(1,N))*RTD
        DO I=1, NRP
          WB1(I) = YRP(I,N)/YRP(NRP,N)
          WB2(I) = (PI-BUFBET(I,N)) * RTD
          BETMIN = MIN(WB2(I),BETMIN)
          BETMAX = MAX(WB2(I),BETMAX)
        ENDDO
      ENDIF
C
      CALL SPLINE(WB2,WB3,WB1,NRP)
C
      DY = 1.0
      DBET = BETMAX - BETMIN
      YLIMS(1) = BETMIN-0.15*DBET
      YLIMS(2) = BETMAX+0.15*DBET
      CALL PLTMOD(NRP,WB1,WB2,DY,YLIMS,
     &            PLFAC,PLPAR,NSPLT,XOFF,XSF,YOFF,YSF)
C
      XPLT = -4.0*CSIZE
      YPLT = PLPAR - 0.5*DY*YSF - 0.7*CSIZE
      CALL PLMATH(XPLT,YPLT,1.4*CSIZE,'b"',0.0,2)
C
      XPLTL = 2.0*CSL
      YPLTL = PLPAR - 2.5*CSL
      IF(OMEGA(N).LE.0.0) THEN
        CALL PLCHAR(XPLTL,YPLTL,CSL,
     &      'Local Coordinates',0.0,17)
      ENDIF
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
C---- get new W2 array
      CALL CRSMOD(NRP,WB1,WB2,WB3,
     &            XOFF,XSF,YOFF,YSF, SH, NSPLT,
     &            LSLOPE, IMOD1,IMOD2 )
C
C---- store buffer blade angles
C
      IF(OMEGA(N).GT.0.0) THEN
        DO I = IMOD1, IMOD2
          DBET = WB2(I) * DTR - BUFBET(I,N)
          BUFBET(I,N)  = BUFBET(I,N) + DBET
        ENDDO
      ELSE
        DO I = IMOD1, IMOD2
          DBET = WB2(I) * DTR - (PI-BUFBET(I,N))
          BUFBET(I,N)  = BUFBET(I,N) - DBET
        ENDDO
      ENDIF
C
      RETURN
      END  !MODBE



      SUBROUTINE MODCH(NDSK)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C------------------------------------------------
C     Takes user cursor input to modify 
C     buffer blade chord array.
C------------------------------------------------
      DIMENSION YLIMS(2)
      EXTERNAL PLCHAR,PLMATH
C
      N = NDSK
      PLFAC = 0.95
      PLPAR = 1.20*PAR
      SH = 0.3*CSIZE
      NSPLT = 20
C
C---- work with temporary arrays
C
      CHMAX = BUFCHD(1,N) * 100.
      DO I=1, NRP
        WB1(I) = YRP(I,N)/YRP(NRP,N)
        WB2(I) = BUFCHD(I,N) * 100.
        CHMAX = MAX(WB2(I),CHMAX)
      ENDDO
      CALL SPLINE(WB2,WB3,WB1,NRP)
C
      DY = 0.1
      YLIMS(1) = 0.
      YLIMS(2) = 1.15*CHMAX
      CALL PLTMOD(NRP,WB1,WB2,DY,YLIMS,
     &            PLFAC,PLPAR,NSPLT,XOFF,XSF,YOFF,YSF)
C
      XPLT = -6.0*CSIZE
      YPLT = PLPAR - 0.5*DY*YSF - 0.7*CSIZE
      CALL PLCHAR(XPLT,YPLT,1.4*CSIZE,'c cm',0.0,4)
C
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
C---- get new W2 array
      CALL CRSMOD(NRP,WB1,WB2,WB3,
     &            XOFF,XSF,YOFF,YSF, SH, NSPLT,
     &            LSLOPE, IMOD1,IMOD2 )
C
C---- store as buffer blade chords
      DO I = IMOD1, IMOD2
        BUFCHD(I,N) = WB2(I)/100.
      ENDDO
C
      RETURN
      END  !MODCH



      SUBROUTINE PLTMOD(N,X,Y,DYMIN,YLIMS,
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
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
      YMIN = Y(1)
      YMAX = Y(1)
      DO I=2, NRP
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
C        YMIN = YMIN - YDEL
c      ENDIF
c      IF((YMAX-YMIN) .LT. 0.5*(YMAX+YMIN)) THEN
c        YMAX = YMAX + YDEL
C        YMIN = YMIN - YDEL
c      ENDIF
C
      DYMIN = YDEL
C
      YSF = PLPAR / (YMAX-YMIN)
      YOFF = YMIN
C
      XSF = 1.0
      XOFF = 0.0
C
C
      CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFAC*SIZE,LPLOT,LLAND)
      CALL PLOTABS(1.0,1.0,-3)
C
      CALL GETCOLOR(ICOL0)
C
      CS = 1.4*CSIZE
      CALL NEWPEN(2)
C
      CALL XAXIS(0.0,0.0,1.0,0.2    , 0.0,0.2, CSIZE,1)
      XPLT = 0.5 - 1.2*CS
      YPLT =     - 2.5*CS
      CALL PLCHAR(XPLT,YPLT,CS,'r/R',0.0,3)
C
      CALL YAXIS(0.0,0.0,PLPAR,YDEL*YSF, YMIN,YDEL, CSIZE,-2)
      IF(LGRID) THEN
       CALL NEWPEN(1)
       NXG = 10
       NYG = INT( (YMAX-YMIN)/(0.5*YDEL) + 0.01 )
       CALL PLGRID(0.0,0.0, NXG,0.1, NYG,0.5*YDEL*YSF, LMASK2 )
      ENDIF
C
      SH = 0.4*CSIZE
      CALL NEWCOLORNAME('blue')
      CALL XYSYMB(NRP,X,Y,XOFF,XSF,YOFF,YSF,SH,1)
C
      CALL NEWCOLOR(ICOL0)
      CALL NEWPEN(2)
      CALL XYLINE(NRP,X,Y,XOFF,XSF,YOFF,YSF,1)
C
      CALL PLFLUSH
C
      RETURN
      END ! PLTMOD


C -----------------------------------------------------------


      SUBROUTINE PLOTBE(NDSK)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C------------------------------------------------
C     Plots current and buffer blade angle arrays
C------------------------------------------------
      DIMENSION YLIMS(2)
      EXTERNAL PLCHAR,PLMATH
C
      N = NDSK
      RTD = 180./PI
      PLFAC = 0.95
      PLPAR = 1.20*PAR
      SH = 0.3*CSIZE
      CSL= 1.0*CSIZE
      NSPLT = 20
C
C---- work with temporary arrays
C
      IF(OMEGA(N).GT.0.0) THEN
        BETMIN = BUFBET(1,N)*RTD
        BETMAX = BUFBET(1,N)*RTD
        DO I=1, NRP
          WB1(I) = YRP(I,N)/YRP(NRP,N)
          WB2(I) = BUFBET(I,N) * RTD
          WC2(I) = CURBET(I,N) * RTD
          BETMIN = MIN(WB2(I),BETMIN)
          BETMAX = MAX(WB2(I),BETMAX)
          BETMIN = MIN(WC2(I),BETMIN)
          BETMAX = MAX(WC2(I),BETMAX)
        ENDDO
      ELSE
        BETMIN = (PI-BUFBET(1,N))*RTD
        BETMAX = (PI-BUFBET(1,N))*RTD
        DO I=1, NRP
          WB1(I) = YRP(I,N)/YRP(NRP,N)
          WB2(I) = (PI-BUFBET(I,N)) * RTD
          WC2(I) = (PI-CURBET(I,N)) * RTD
          BETMIN = MIN(WB2(I),BETMIN)
          BETMAX = MAX(WB2(I),BETMAX)
          BETMIN = MIN(WC2(I),BETMIN)
          BETMAX = MAX(WC2(I),BETMAX)
        ENDDO
      ENDIF
C
      DY = 1.0
      DBET = BETMAX - BETMIN
      YLIMS(1) = BETMIN-0.15*DBET
      YLIMS(2) = BETMAX+0.15*DBET
      CALL PLTMODG(NRP,WB1,WB2,WC2,DY,YLIMS,
     &            PLFAC,PLPAR,NSPLT,XOFF,XSF,YOFF,YSF)
C
      XPLT = -4.0*CSIZE
      YPLT = PLPAR - 0.5*DY*YSF - 0.7*CSIZE
      CALL PLMATH(XPLT,YPLT,1.4*CSIZE,'b"',0.0,2)
C
      XPLTL = 2.0*CSL
      YPLTL = PLPAR - 2.5 * CSL
      IF(OMEGA(N).LE.0.0) THEN
        CALL PLCHAR(XPLTL,YPLTL,CSL,
     &      'Local Coordinates',0.0,17)
      ENDIF
C
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END  !PLOTBE



      SUBROUTINE PLOTCH(NDSK)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C------------------------------------------------
C     Plots current and buffer blade chord arrays
C------------------------------------------------
      DIMENSION YLIMS(2)
      EXTERNAL PLCHAR,PLMATH
C
      N = NDSK
      PLFAC = 0.95
      PLPAR = 1.20*PAR
      SH = 0.3*CSIZE
      NSPLT = 20
C
C---- work with temporary arrays
C
      CHMAX = BUFCHD(1,N) * 100.
      DO I=1, NRP
        WB1(I) = YRP(I,N)/YRP(NRP,N)
        WB2(I) = BUFCHD(I,N) * 100.
        WC2(I) = CURCHD(I,N) * 100.
        CHMAX = MAX(WB2(I),CHMAX)
        CHMAX = MAX(WC2(I),CHMAX)
      ENDDO
C
      DY = 0.1
      YLIMS(1) = 0.
      YLIMS(2) = 1.15*CHMAX
      CALL PLTMODG(NRP,WB1,WB2,WC2,DY,YLIMS,
     &            PLFAC,PLPAR,NSPLT,XOFF,XSF,YOFF,YSF)
C
      CALL NEWPEN(2)
      CALL NEWCOLORNAME('black')
      XPLT = -6.0*CSIZE
      YPLT = PLPAR - 0.5*DY*YSF - 0.7*CSIZE
      CALL PLCHAR(XPLT,YPLT,1.4*CSIZE,'c cm',0.0,4)
C
      CALL PLFLUSH
      CALL REPLOT(IDEV) ! Crude hack to fix win graphics 
C
      RETURN
      END  !PLOTCH



      SUBROUTINE PLTMODG(N,X,Y,Z,DYMIN,YLIMS,
     &                  PLFAC,PLPAR,NSPLT,XOFF,XSF,YOFF,YSF)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      DIMENSION X(N),Y(N),Z(N)
      DIMENSION YLIMS(2)
      DIMENSION XLEG(2),YLEG(2)
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
C     Modified for plotting of buffer and current chord & beta
C     Esotec Developments 2013
C ------------------------------------------------------------
      DATA LMASK1, LMASK2, LMASK3 / -32640, -30584, -21846 /
C
      YMIN = Y(1)
      YMAX = Y(1)
      DO I=2, NRP
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
C      IF((YMAX-YMIN) .LT. 0.5*(YMAX+YMIN)) THEN
C        YMAX = YMAX + YDEL
C        YMIN = YMIN - YDEL
C      ENDIF
C      IF((YMAX-YMIN) .LT. 0.5*(YMAX+YMIN)) THEN
C        YMAX = YMAX + YDEL
C        YMIN = YMIN - YDEL
C      ENDIF
C
      DYMIN = YDEL
C
      YSF = PLPAR / (YMAX-YMIN)
      YOFF = YMIN
C
      XSF = 1.0
      XOFF = 0.0
C
C
      CALL PLTINI(SCRNFR,IPSLU,IDEV,PLFAC*SIZE,LPLOT,LLAND)
      CALL PLOTABS(1.0,1.0,-3)
C
      CALL GETCOLOR(ICOL0)
C
      CS = 1.4*CSIZE
      CALL NEWPEN(2)
C
      CALL XAXIS(0.0,0.0,1.0,0.2    , 0.0,0.2, CSIZE,1)
      XPLT = 0.5 - 1.2*CS
      YPLT =     - 2.5*CS
      CALL PLCHAR(XPLT,YPLT,CS,'r/R',0.0,3)
C
      CALL YAXIS(0.0,0.0,PLPAR,YDEL*YSF, YMIN,YDEL, CSIZE,-2)
      IF(LGRID) THEN
       CALL NEWPEN(1)
       NXG = 10
       NYG = INT( (YMAX-YMIN)/(0.5*YDEL) + 0.01 )
       CALL PLGRID(0.0,0.0, NXG,0.1, NYG,0.5*YDEL*YSF, LMASK2 )
      ENDIF
C
      SH = 0.4*CSIZE
      CSL = 1.0*CSIZE
      XPLT = 2.0*CSL
      YPLT = 0.06
      YPLT2= YPLT - 2.5*CSL
C
      XLEG(1)= XPLT + 14.*CSL
      XLEG(2)= XLEG(1) + 0.08
      YLEG(1)= YPLT + 0.4*CSL
      YLEG(2)= YLEG(1)
C
C----buffer blade
C
      CALL NEWCOLORNAME('green')
      CALL XYSYMB(NRP,X,Y,XOFF,XSF,YOFF,YSF,SH,1)
C
      CALL NEWPEN(2)
      CALL XYLINE(NRP,X,Y,XOFF,XSF,YOFF,YSF,1)
C
      CALL PLCHAR(XPLT,YPLT,CSL,'Buffer  blade',0.0,13)
      CALL XYLINE(2,XLEG,YLEG,0.,1.,0.,1.,1)
C
C----current blade
C
      YLEG(1) = YLEG(1) - 2.5*CSL
      YLEG(2) = YLEG(1)
C
      CALL NEWCOLORNAME('blue')
      CALL NEWPEN(1)
      CALL XYSYMB(NRP,X,Z,XOFF,XSF,YOFF,YSF,SH,1)
C
      CALL NEWPEN(2)
      CALL XYLINE(NRP,X,Z,XOFF,XSF,YOFF,YSF,2)
C
      CALL PLCHAR(XPLT,YPLT2,CSL,'Current blade',0.0,13)
      CALL XYLINE(2,XLEG,YLEG,0.,1.,0.,1.,2)
C
      CALL NEWCOLORNAME('black')
      CALL PLFLUSH
C
      RETURN
      END  ! PLTMODG
C
C-------------------------------------------------------------------


      SUBROUTINE BLDPLTM(NR,VIEW)
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      EXTERNAL PLCHAR,PLMATH
C
      CHARACTER*(*) VIEW
      CHARACTER*80  PNAME
C------------------------------------------------
C     Plots buffer blade planform and projections
C     Modified from BLDPLT
C------------------------------------------------
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
      CALL STRIP(PNAME,NPN)
C
C---- case title
      CALL NEWPEN(2)
      XT = 0.5
      YT = -0.030
      CALL PLCHAR(XT,YT,CSL,PNAME,0.0,NPN)
      CALL PLCHAR(999.,YT,CSL,': Buffer Blade',0.0,14)
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
C      CALL GETROTG(N)
C
C---- mean radius
      YMNR = SQRT(0.5*ADISK(N)/PI + RHUB(N)**2)
C
C---- blade solidity
      CALL SPLINE(BUFCHD(1,N),W1,YRP,NRP)
      CMNR = SEVAL(YMNR,BUFCHD(1,N),W1,YRP,NRP)
      SIGMR = FLOAT(NRBLD(N))*CMNR/(2.0*PI*YMNR)
C
C---- blade angles and twist
      BHUB = BUFBET(1,N)
      BTIP = BUFBET(NRP,N)
C      B2   = BETAR(2,N)
C      BEXT = BHUB - YRC(1,N)*(BHUB-B2)/(YRC(1,N)-YRC(2,N))
      BTWIST = BHUB - BTIP
      CALL SPLINE(BUFBET(1,N),W2,YRP,NRP)
      BETAMR = SEVAL(YMNR,BUFBET(1,N),W2,YRP,NRP)
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
       CALL PLNUMB(999.,YL,CS, BUFCHD(NRP,N) ,0.0,4)
       CALL PLSUBS(XL3 ,YL,CS,'tip'       ,0.0,3,PLCHAR)
       CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, BUFBET(NRP,N)/DTR,0.0,2)
C
       YL = YL - 2.5*CS
       CALL PLCHAR(XL1 ,YL,CS,'Rhub = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, RHUB(N)    ,0.0,4)
       CALL PLCHAR(XL2 ,YL,CS,'Chub = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, BUFCHD(1,N)   ,0.0,4)
       CALL PLSUBS(XL3 ,YL,CS,'hub'       ,0.0,3,PLCHAR)
       CALL PLMATH(XL3 ,YL,CS,'b    = '   ,0.0,7)
       CALL PLNUMB(999.,YL,CS, BUFBET(1,N)/DTR ,0.0,2)
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
        SINB = SIN(BUFBET(I,N))
        COSB = COS(BUFBET(I,N))
C
        XI(I) =  YRP(I,N)/RTIP(N)
C------ projection chords in front and behind radial axis
        W1(I) =  XPAXIS     *SINB*BUFCHD(I,N)/RTIP(N)
        W2(I) = (XPAXIS-1.0)*SINB*BUFCHD(I,N)/RTIP(N)
C
        W3(I) =  XPAXIS     *COSB*BUFCHD(I,N)/RTIP(N)
        W4(I) = (XPAXIS-1.0)*COSB*BUFCHD(I,N)/RTIP(N)
C
        W5(I) =  XPAXIS     *BUFCHD(I,N)/RTIP(N)
        W6(I) = (XPAXIS-1.0)*BUFCHD(I,N)/RTIP(N)
C
        CHUMAX = MAX(CHUMAX,     BUFCHD(I,N)/RTIP(N))
        CHLMAX = MAX(CHLMAX,SINB*BUFCHD(I,N)/RTIP(N))
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
        DX = DS*COS(BUFBET(1,N))
        DY = DS*SIN(BUFBET(1,N))
        CALL PLOT(-DX,-DY,3)
        CALL PLOT( DX, DY,2)
C
C------ and tip-angle lines
        DX = DS*COS(BUFBET(NRL,N))
        DY = DS*SIN(BUFBET(NRL,N))
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
      END ! BLDPLTM




