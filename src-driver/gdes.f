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
C     GNAM command to edit geometry name
C     
C========================================================================
C
C
      SUBROUTINE DFGDES(LMODG)
C------------------------------------------------------
C     Geometry design routine. 
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      CHARACTER*4 COMAND, COMOLD
      CHARACTER*1 CIND
      LOGICAL LRECALC, LOFINI, LMODPL, LMODG
C
      CHARACTER*16 PROMPT, PROMPT2
      CHARACTER*80 LINE
      CHARACTER*128 COMARG, ARGOLD
C
      LOGICAL LPLTEL(NEX)
      DIMENSION XBNEW(IBX), YBNEW(IBX)
C
      DIMENSION XBOX(2), YBOX(2)
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
      LOGICAL LSAME, OK
C
      LOGICAL INSID
      LOGICAL INSIDE
C
      SAVE COMOLD, ARGOLD
C
      COMAND = '****'
      COMARG = ' '
      LRECALC = .FALSE.
      LMODG   = .FALSE.
C
c      IF(NBTOT.EQ.0) THEN
c       WRITE(*,*)
c       WRITE(*,*) '***  No airfoil available  ***'
c       RETURN
c      ENDIF
C
C---- startup default: plot all elements
      DO IEL=1, NBEL
        LPLTEL(IEL) = .TRUE.
      ENDDO
C
ccc      LSYM = .TRUE.
C
      WRITE(*,*)
      WRITE(*,*) 'You are working with the buffer airfoil'
C
      NTGSPL = 10
C
      ARAD = 0.
      DELX = 0.
      DELY = 0.
      XFAC = 1.0
      YFAC = 1.0
      FAC = 1.0
C
C---- make sure we have a valid target element
      IF(IELGDES.LT.1 .OR. IELGDES.GT.NBEL) THEN
       IF(NBEL.LE.1) THEN
        IELGDES = 1
       ELSE
 5      CALL ASKI('Specify target element index^',IELGDES)
        IF(IELGDES.LT.1 .OR. IELGDES.GT.NBEL) THEN
         WRITE(*,*) 'Number of elements present:', NBEL
         IELGDES = 1
         GO TO 5
        ENDIF
       ENDIF
      ENDIF
C
C---- set plot limits
      IELPLT = 0
      DO IEL = 1, NBEL
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL MINMAX(N,XB(I),XBMINE(IEL),XBMAXE(IEL))
        CALL MINMAX(N,YB(I),YBMINE(IEL),YBMAXE(IEL))
      ENDDO
C
C---- go plot geometry with scale,offset initialization
      LOFINI = .TRUE.
      COMAND = 'PLOT'
      GO TO 510
CC
C---- begin user-interaction loop
C................................................................
C
 500  CONTINUE
      COMOLD = COMAND
      ARGOLD = COMARG
C
 501  CONTINUE
      IF(LGSYMM) THEN
       PROMPT = '.GDESs^'
      ELSE
       PROMPT = '.GDES ^'
      ENDIF
C
      IF(NBEL.GT.1) THEN
C------ add current target element to prompt
        IF    (IELGDES.LT. 10) THEN
         PROMPT(7:8) = '(' // PNUM(IELGDES)(2:2)
         K = 9
        ELSEIF(IELGDES.LT.100) THEN
         PROMPT(7:9) = '(' // PNUM(IELGDES)
         K = 10
        ENDIF
        KEL = IELGDES
        PROMPT(K:K+1) = ')^'
      ENDIF
C
      CALL ASKC(PROMPT,COMAND,COMARG)
C
      IF(COMAND.EQ.'    ') THEN
       IF(LPLOT) THEN
        CALL PLEND
        CALL CLRZOOM
       ENDIF
       LPLOT = .FALSE.
C
       IF(.NOT.LGSAME) THEN
        WRITE(*,*)
        WRITE(*,*) 'Buffer geometry not identical to current geometry'
       ENDIF
       RETURN
      ENDIF
C
C---- process previous command ?
      IF(COMAND(1:1).EQ.'!') THEN
        IF(COMOLD.EQ.'****') THEN
          WRITE(*,*) 'Previous .GDES command not valid'
          GO TO 501
        ELSE
          COMAND = COMOLD
          COMARG  = ARGOLD
          LRECALC = .TRUE.
        ENDIF
      ELSE
        LRECALC = .FALSE.
      ENDIF
C
 505  CONTINUE
C
C---- was only a new target-element index typed?
      CIND = COMAND(1:1)
      IF(INDEX('0123456789',CIND) .NE. 0) THEN
       READ(COMAND,*,ERR=508) IEL
       IF(IEL.LT.1 .OR. IEL.GT.NBEL) GO TO 508
       IELGDES = IEL
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
        IEL = IELGDES
        WRITE(*,1102)
 1102   FORMAT(
     &//'   GSET     Set buffer  geometry <== current geometry'
     & /'   EXEC     Set current geometry <== buffer  geometry')
C
        IF(NBEL.GT.1) WRITE(*,1104) IEL
 1104   FORMAT(
     & /'        i   Change target element (currently =', I3, ')')
C
        WRITE(*,1106) XBREFE(IEL),YBREFE(IEL)
 1106   FORMAT(
     & /'   REF  rr  Set element reference pt (currently =', 2G12.5,')'
C
     &//'   TRAN rr  Translate'
     & /'   SCAL r   Scale about reference pt'
     & /'   ADEG r   Rotate (degrees) about reference pt'
     & /'   ARAD r   Rotate (radians) about reference pt'
     & /'   HSET     Set element home to current position'
     & /'   HOME     Put element back to home'
C
     &//'   TGAP rr  Change trailing edge gap'
     & /'   LERA rr  Change leading edge radius'
     & /'   TCAM rr  Change thickness, camber'
     & /'   CAMP     Add to camber line from added delta(Cp)'
     & /'   CAMS     Add to camber line from splined input points'
     & /'   FLAP rrr Deflect trailing edge flap'
C
     &//'   MODG     Modify geometry via cursor'
     & /'   SLOP     Toggle modified-contour slope matching flag'
     & /'   SYMM     Toggle y-symmetry flag'
C
     &//'  .POIN     Pointwise modification options'
     & /'   DIST     Determine distance between 2 cursor points'
C
     &//'   CLIS     List curvatures'
     & /'   CANG     List panel corner angles'
C
     &//'   TICK     Toggle node tick-mark plotting'
     & /'   GRID     Toggle grid plotting'
     & /'   GPAR     Toggle geometric parameter plotting'
     & /'   HPLO     Toggle home-position plotting'
     & /'   OVER f   Overlay disk file airfoil'
     & /'   GNAM     Change duct/centerbody geometry name'
C
     &//'   PLOT     Plot buffer airfoil'
     & /'   BLOW     Blowup plot'
     & /'   RESE     Reset to original plot scale'
     & /'  .ANNO     Annotate plot'
     & /'   HARD     Hardcopy current plot'
     & /'   SIZE r   Change absolute plot size'
     & /'   Z        Zoom  '
     & /'   U        Unzoom')
C
C
C--------------------------------------------------------------------

      ELSEIF(COMAND.EQ.'GNAM') THEN
        WRITE(*,2020) ANAME
 2020   FORMAT(/'  Current geometry name is :  ',A)
C
        CALL ASKS(' Enter new geometry name^',ANAME)
C
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'GSET') THEN
C---- Only reset geometry for first NBEL elements in buffer geometry
        DO IEL = 1, NBEL
          IP1 = IPFRST(IEL)
          IP2 = IPLAST(IEL)
          DO IP=IP1, IP2
            XB(IP) = XP(IP)
            YB(IP) = YP(IP)
          ENDDO
          NBTOT = IP2
          IBFRST(IEL) = IP1
          IBLAST(IEL) = IP2
C------ spline new coordinates
          N = IP2 - IP1 + 1
          I = IP1
          CALL SCALC(XB(I),YB(I),SB(I),N)
          CALL SEGSPL(XB(I),XBS(I),SB(I),N)
          CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
          CALL CLRHOM(IEL)
        ENDDO
C
        LGSAME = .TRUE.
C
        LOFINI = .FALSE.
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'EXEC') THEN
cc        CALL PANCOP
C
        LOFINI = .FALSE.
        LMODG  = .TRUE.
        LGSAME = .FALSE.
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ELEM'
     &  .OR. COMAND.EQ.'E   ') THEN
    8   IF(NINPUT.GE.1) THEN
         IELGDES = IINPUT(1)
        ELSE
         CALL ASKI('Enter target element index^',IELGDES)
        ENDIF
C
        IF(IELGDES.LT.1 .OR. IELGDES.GT.NBEL) THEN
         WRITE(*,*) 'Number of elements present:', NBEL
         NINPUT = 0
         GO TO 8
        ENDIF
C
        LOFINI = .TRUE.
        CALL GPLINI(LOFINI,NBEL,LPLTEL,
     &              XBMINE,XBMAXE,YBMINE,YBMAXE,
     &              XBOFF,XBSF,YBOFF,YBSF)
        CALL GBPLOT(LPLTEL)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'REF ') THEN
        IEL = IELGDES
        IF(NINPUT.GE.2) THEN
          XBREFE(IEL) = RINPUT(1)
          YBREFE(IEL) = RINPUT(2)
        ELSEIF(NINPUT.GE.1) THEN
         RINPUT(2) = YBREFE(IEL)
         CALL ASKR('Enter Yref^',RINPUT(2))
         XBREFE(IEL) = RINPUT(1)
         YBREFE(IEL) = RINPUT(2) 
        ELSE
         WRITE(*,1108) XBREFE(IEL),YBREFE(IEL)
         COMOLD = COMAND
         ARGOLD = COMARG
         PROMPT2 = '.REF^'
         CALL ASKC(PROMPT2,COMAND,COMARG)
         IF(    COMAND(1:1).EQ.'L') THEN
          XBREFE(IEL) = XBLE(IEL)
          YBREFE(IEL) = YBLE(IEL) 
         ELSEIF(COMAND(1:1).EQ.'T') THEN
          XBREFE(IEL) = XBTE(IEL)
          YBREFE(IEL) = YBTE(IEL) 
         ELSEIF(COMAND(1:1).EQ.'Q') THEN
          XBREFE(IEL) = XBLE(IEL) + 0.25*(XBTE(IEL)-XBLE(IEL))
          YBREFE(IEL) = YBLE(IEL) + 0.25*(YBTE(IEL)-YBLE(IEL))
         ELSEIF(COMAND(1:1).EQ.'C') THEN
          WRITE(*,*) 'Click on new reference point location...'
          CALL GETCURSORXY(XCRS,YCRS,CHKEY)
          XBREFE(IEL) = XCRS/XBSF + XBOFF
          YBREFE(IEL) = YCRS/YBSF + YBOFF
         ELSEIF(COMAND(1:1).EQ.'K') THEN
          RINPUT(1) = XBREFE(IEL)
          RINPUT(2) = YBREFE(IEL)
          CALL ASKRN('Enter Xref,Yref^',RINPUT,2)
          XBREFE(IEL) = RINPUT(1)
          YBREFE(IEL) = RINPUT(2) 
         ENDIF
       ENDIF
C
        COMAND = 'PLOT'
        GO TO 510
C
 1108   FORMAT(
     & /'  Set element reference pt (currently =', 2G12.5,')'
     & /'   L E       set reference pt to LE'
     & /'   T E       set reference pt to LE'
     & /'   Q uarter  set reference pt to 1/4 pt on chordline'
     & /'   C urs     set reference pt from cursor'
     & /'   K ey      set reference pt from keyboard input')
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'TRAN'
     &  .OR. COMAND.EQ.'T   ') THEN
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
        IEL = IELGDES
        IF(IEL.EQ.1) THEN
          WRITE(*,*) 'DeltaY set to 0.0 for centerbody'
          DELY = 0.0
        ENDIF
        CALL TRAN(IEL,DELX,DELY)
C
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL MINMAX(N,XB(I),XBMINE(IEL),XBMAXE(IEL))
        CALL MINMAX(N,YB(I),YBMINE(IEL),YBMAXE(IEL))
C
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SCAL'
     &  .OR. COMAND.EQ.'S   ') THEN
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
        IEL = IELGDES
        CALL SCAL(IEL,XFAC,YFAC,XBREFE(IEL),YBREFE(IEL))
C
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL MINMAX(N,XB(I),XBMINE(IEL),XBMAXE(IEL))
        CALL MINMAX(N,YB(I),YBMINE(IEL),YBMAXE(IEL))
C
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ADEG'
     &  .OR. COMAND.EQ.'A   ') THEN
        IF(NINPUT.GE.1) THEN
         ADEG = RINPUT(1)
        ELSE
         ADEG = 0.
         CALL ASKR('Enter angle change (+clockwise, deg)^',ADEG)
        ENDIF
C
        IEL = IELGDES
        IF(IEL.EQ.1) THEN
          WRITE(*,*) 'Angle set to 0.0 for centerbody'
          ADEG = 0.0
        ENDIF
        CALL ANGL(IEL,ADEG,XBREFE(IEL),YBREFE(IEL))
C
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL MINMAX(N,XB(I),XBMINE(IEL),XBMAXE(IEL))
        CALL MINMAX(N,YB(I),YBMINE(IEL),YBMAXE(IEL))
C
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ARAD') THEN
        IF(NINPUT.GE.1) THEN
         ARAD = RINPUT(1)
        ELSE
         ARAD = 0.
         CALL ASKR('Enter angle change (+clockwise, radians)^',ADEG)
        ENDIF
        ADEG = ARAD * PI/180.0
C
        IEL = IELGDES
        IF(IEL.EQ.1) THEN
          WRITE(*,*) 'Angle set to 0.0 for centerbody'
          ADEG = 0.0
        ENDIF
        CALL ANGL(IEL,ADEG,XBREFE(IEL),YBREFE(IEL))
C
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL MINMAX(N,XB(I),XBMINE(IEL),XBMAXE(IEL))
        CALL MINMAX(N,YB(I),YBMINE(IEL),YBMAXE(IEL))
C
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'HSET') THEN
        IEL = IELGDES
        CALL CLRHOM(IEL)
C
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'HOME') THEN
        IEL = IELGDES
        DELX = -DXBSUM(IEL)
        DELY = -DYBSUM(IEL)
        ADEG = -AGBSUM(IEL)
        XFAC = 1.0/XFBSUM(IEL)
        YFAC = 1.0/YFBSUM(IEL)
C
        CALL TRAN(IEL,DELX,DELY)
        CALL ANGL(IEL,ADEG,     XBREFE(IEL),YBREFE(IEL))
        CALL SCAL(IEL,XFAC,YFAC,XBREFE(IEL),YBREFE(IEL))
C
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL MINMAX(N,XB(I),XBMINE(IEL),XBMAXE(IEL))
        CALL MINMAX(N,YB(I),YBMINE(IEL),YBMAXE(IEL))
C
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'MODG'
     &  .OR. COMAND.EQ.'M   ') THEN
        CALL GETCOLOR(ICOL0)
C
C------ put up Qspec(s) plot if there isn't one on the screen
        IF(.NOT.LGPLOT) THEN
         LOFINI = .TRUE.
         CALL GPLINI(LOFINI,NBEL,LPLTEL,
     &               XBMINE,XBMAXE,YBMINE,YBMAXE,
     &               XBOFF,XBSF,YBOFF,YBSF)
         CALL GBPLOT(LPLTEL)
        ENDIF
C
C------ box in which modification cursor inputs are recognized
        CALL GETZOOMABS(XZOFF,YZOFF,XZFAC,YZFAC)
        XBOX(1) = MIN(0.0  ,       1.0001*XMARG)/XZFAC - XZOFF
        XBOX(2) = MIN(XWIND, XPAGE-1.0001*XMARG)/XZFAC - XZOFF
        YBOX(1) = MAX(0.0  ,       1.0001*YMARG)/YZFAC - YZOFF
        YBOX(2) = MAX(YWIND, YPAGE-1.0001*YMARG)/YZFAC - YZOFF
C
C------ modify geometry of target element
        LMODPL = .FALSE.
        CALL MODIXY(IBX,IBFRST,IBLAST,NBEL,
     &              XB,YB,XBS,YBS,SB, LGSLOP,
     &              IBMOD1,IBMOD2,IELMOD,
     &              XBOX,YBOX, XBOX,YBOX,
     &              XBOFF,YBOFF,XBSF,YBSF, LMODPL)
C
        IF(IELMOD.GE.1     .AND.
     &     IELMOD.LE.NBEL   .AND.
     &     IBMOD1.LT.IBMOD2     ) THEN
C------- process changed element
         CALL ELPROC(IELMOD)
C
C------- replot just the modified piece
         CALL NEWCOLORNAME('MAGENTA')
         I = IBMOD1
         N = IBMOD2 - IBMOD1 + 1
         CALL PLTGSP(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, 
     &               XBOFF,XBSF,YBOFF,YBSF, 20)
C
C------- add symbols at blend points
         CALL PLSYMB((XB(IBMOD1)-XBOFF)*XBSF,
     &               (YB(IBMOD1)-YBOFF)*YBSF,SHGT,1,0.0,0)
         CALL PLSYMB((XB(IBMOD2)-XBOFF)*XBSF,
     &               (YB(IBMOD2)-YBOFF)*YBSF,SHGT,1,0.0,0)
C
         CALL NEWCOLOR(ICOL0)
         CALL PLFLUSH
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SLOP') THEN
C------ input segment endpoint tangency toggle
        LGSLOP = .NOT.LGSLOP
        IF(     LGSLOP) THEN
         WRITE(*,*)
     &     'Modified segment will be made tangent at endpoints'
        ELSE
         WRITE(*,*)
     &     'Modified segment will not be made tangent at endpoints'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SYMM') THEN
        LGSYMM = .NOT.LGSYMM
        IF(LGSYMM) THEN
         WRITE(*,*) 'y-symmetry forcing enabled.'
         CALL ZERCAM(IELGDES)
        ELSE
         WRITE(*,*) 'y-symmetry forcing disabled.'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'TGAP') THEN
C------ set new trailing edge gap
        IEL = IELGDES
        IF(.NOT.LBODY(IEL)) THEN
         WRITE(*,*) 'TE gap is defined only for closed elements'
         GO TO 500
        ENDIF
C
        IB1 = IBFRST(IEL)
        IB2 = IBLAST(IEL)
        DXB = XB(IB1) - XB(IB2)
        DYB = YB(IB1) - YB(IB2)
        GAP = SQRT(DXB**2 + DYB**2)
C
        IF    (NINPUT .GE. 2) THEN
         GAPNEW = RINPUT(1)
         DOC    = RINPUT(2)
        ELSEIF(NINPUT .GE. 1) THEN
         GAPNEW = RINPUT(1)
         CALL ASKR('Enter blending distance/c (0..1)^',DOC)
        ELSE
         GAPNEW = GAP
         CALL ASKR('Enter new TE gap^',GAPNEW)
         CALL ASKR('Enter blending distance/c (0..1)^',DOC)
        ENDIF
        DOC = MIN( MAX( DOC , 0.0 ) , 1.0 )
C
        I = IB1
        N = IB2 - IB1 + 1
        CALL TEGAP(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, DOC,GAPNEW)
        LGSAME = .FALSE.
C
        CALL NEWCOLORNAME('MAGENTA')
        CALL NEWPEN(2)
        CALL PLTGSP(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, 
     &              XBOFF,XBSF,YBOFF,YBSF, NTGSPL)
        CALL NEWCOLOR(ICOL0)
        CALL PLFLUSH
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'LERA') THEN
C------ set new LE radius
        IEL = IELGDES
C
        IF(.NOT.LBODY(IEL)) THEN
         WRITE(*,*) 'LE radius is defined only for closed elements'
         GO TO 500
        ENDIF
C
        IB1 = IBFRST(IEL)
        IB2 = IBLAST(IEL)
C
        IF    (NINPUT .GE. 2) THEN
         RFAC = RINPUT(1)
         DOC  = RINPUT(2)
        ELSEIF(NINPUT .GE. 1) THEN
         RFAC = RINPUT(1)
         CALL ASKR('Enter blending distance/c from LE^',DOC)
        ELSE
         CALL ASKR('Enter new/old LE radius scaling ratio^',RFAC)
         CALL ASKR('Enter blending distance/c from LE^',DOC)
        ENDIF
        DOC = MAX( DOC , 0.001 )
C
        I = IB1
        N = IB2 - IB1 + 1
        CALL LERSCL(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, 
     &              DOC,RFAC, XBNEW(I),YBNEW(I))
C
        DO IB = IB1, IB2
          XB(IB) = XBNEW(IB)
          YB(IB) = YBNEW(IB)
        ENDDO
C
C------ spline new coordinates
        CALL SCALC(XB(I),YB(I),SB(I),N)
        CALL SEGSPL(XB(I),XBS(I),SB(I),N)
        CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
C------ find max curvature
        CVMAX = 0.
        DO IB = IB1+1, IB2-1
          CV = CURV(SB(IB),XB(I),XBS(I),YB(I),YBS(I),SB(I),N)
          CVMAX = MAX( ABS(CV) , CVMAX )
        ENDDO
C
        RADIUS = 1.0/CVMAX
C
        WRITE(*,4200) 100.0*RADIUS
 4200   FORMAT(/' New LE radius = ',F7.3,' %')
C
        LGSAME = .FALSE.
C
        CALL NEWCOLORNAME('MAGENTA')
        CALL NEWPEN(2)
        CALL PLTGSP(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, 
     &              XBOFF,XBSF,YBOFF,YBSF, NTGSPL)
        CALL NEWCOLOR(ICOL0)
        CALL PLFLUSH
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'CLIS') THEN
C------ display curvatures
        IEL = IELGDES
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL CLIS(XB(I),XBS(I),YB(I),YBS(I),SB(I),N)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'CANG') THEN
C------ display corner angles
        IEL = IELGDES
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL CANG(XB(I),YB(I),N,.TRUE.)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'TCAM') THEN
C------ show/change thickness and camber
        IEL = IELGDES
C
        IF(.NOT.LBODY(IEL)) THEN
         WRITE(*,*) 'Thickness is defined only for closed elements'
         GO TO 500
        ENDIF
C
        IB1 = IBFRST(IEL)
        IB2 = IBLAST(IEL)
C
        I = IB1
        N = IB2 - IB1 + 1
C
C------ first see what the current thickness,camber are
        TFAC = 1.0
        CFAC = 1.0
        CALL TCSET(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, 
     &             TMAX,CMAX, TFAC,CFAC, XBNEW(I),YBNEW(I), .FALSE.)
        WRITE(*,4300) TMAX, CMAX
 4300     FORMAT(/' Current thickness:', F13.6,
     &           /' Current    camber:', F13.6)
C
        IF    (NINPUT .GE. 2) THEN
         TNEW = RINPUT(1)
         CNEW = RINPUT(2)
        ELSEIF(NINPUT .GE. 1) THEN
         TNEW = RINPUT(1)
         CNEW = 999.0
         IF(LGSYMM) THEN
          WRITE(*,*) 'Symmetry enforced:  Maintaining zero camber.'
         ELSE
          CALL ASKR('Enter new max  camber  (ret = no change)^',CNEW)
         ENDIF
        ELSE
         TNEW = TMAX
         CNEW = CMAX
         CALL ASKR('Enter new max thickness (ret = no change)^',TNEW)
         IF(LGSYMM) THEN
          WRITE(*,*) 'Symmetry enforced:  Maintaining zero camber.'
         ELSE
          CALL ASKR('Enter new max  camber   (ret = no change)^',CNEW)
         ENDIF
        ENDIF
C
        TFAC = 1.0
        CFAC = 1.0
        IF(TMAX.NE.0.0 .AND. TNEW.NE.999.0) TFAC = TNEW/TMAX
        IF(CMAX.NE.0.0 .AND. CNEW.NE.999.0) CFAC = CNEW/CMAX
        IF(LGSYMM) CFAC = 0.0
C
C------ sanity checks on scaling factors
        IF(ABS(TFAC) .GT. 1.0 .OR. ABS(CFAC) .GT. 1.0) THEN
          WRITE(*,4310) TFAC, CFAC
 4310     FORMAT(/' Implied scaling factors are:', F10.3,' x thickness'
     &           /'                             ', F10.3,' x camber   ')
          CALL ASKL('Apply scaling factors?^',OK)
          IF(.NOT.OK) THEN
            WRITE(*,*) 'No action taken'
            GO TO 500
          ENDIF
        ENDIF
C
        CALL TCSET(XB(I),XBS(I),YB(I),YBS(I),SB(I),N,
     &             TMAX,CMAX, TFAC,CFAC, XBNEW(I),YBNEW(I), .TRUE.)
C
        DO IB = IB1, IB2
          XB(IB) = XBNEW(IB)
          YB(IB) = YBNEW(IB)
        ENDDO
C
C------ spline new coordinates
        CALL SCALC(XB(I),YB(I),SB(I),N)
        CALL SEGSPL(XB(I),XBS(I),SB(I),N)
        CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
        LGSAME = .FALSE.
C
        LOFINI = .FALSE.
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'CAMP') THEN
C------ modify camber via change to loading delta(Cp)
        IF(LGSYMM) THEN
         WRITE(*,*) 'Disabling symmetry enforcement.'
         LGSYMM = .FALSE.
        ENDIF
C
        IEL = IELGDES
        CALL CAMP(IEL)
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'CAMS') THEN
C------ modify camber via splined input 
        IF(LGSYMM) THEN
         WRITE(*,*) 'Disabling symmetry enforcement.'
         LGSYMM = .FALSE.
        ENDIF
C
        IEL = IELGDES
        CALL CAMS(IEL)
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'FLAP') THEN
C------ deflect flap
        IEL = IELGDES
        IF(NINPUT.GE.2) THEN
         XBF = RINPUT(1)
         YBF = RINPUT(2)
        ELSE
         XBF = -999.0
         YBF = -999.0
        ENDIF
C
        IB1 = IBFRST(IEL)
        IB2 = IBLAST(IEL)
C
        I = IB1
        N = IB2 - IB1 + 1
        CALL GETXYF(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, 
     &              TOPS,BOTS,XBF,YBF)
        INSID = INSIDE(XB(I),YB(I),N,XBF,YBF)
C
        WRITE(*,4500) XBF, YBF
 4500   FORMAT(/' Flap hinge: x,y =', 2F9.5 )
C
        IF(NINPUT.GE.3) THEN
         DDEF = RINPUT(3)
        ELSE
         DDEF = 0.0
         CALL ASKR('Enter flap deflection in degrees (+ down)^',DDEF)
        ENDIF
        RDEF = DDEF*PI/180.0
        IF(RDEF .EQ. 0.0) RETURN
C
        I = IB1
        N = IB2 - IB1 + 1
        CALL FLAP(XB(I),XBS(I),YB(I),YBS(I),SB(I),N,
     &            XBF,YBF,RDEF,INSID, 
     &            XBNEW,YBNEW,NNEW,
     &            TOPS,ATOP,XTOP,YTOP,
     &            BOTS,ABOT,XBOT,YBOT )
C
C------ move array slots for all higher elements to fit any new points
        NADD = NNEW - N
        CALL BSHIFT(IEL+1,IBFRST(IEL+1),NADD)
C------ remember to change pointer to first point
        IBFRST(IEL+1) = IBFRST(IEL+1) + NADD
C
C------ insert updated element with flap
        IBLAST(IEL) = IBFRST(IEL) + NNEW - 1
        INEW  = 0
        DO IB = IBFRST(IEL), IBLAST(IEL)
          INEW = INEW + 1
          XB(IB) = XBNEW(INEW)
          YB(IB) = YBNEW(INEW)
        ENDDO
C
        I = IBFRST(IEL)
        N = IBLAST(IEL) - IBFRST(IEL) + 1
        CALL SCALC(XB(I),YB(I),SB(I),N)
        CALL SEGSPL(XB(I),XBS(I),SB(I),N)
        CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
C---- replot 
        CALL GPLINI(LOFINI,NBEL,LPLTEL,
     &              XBMINE,XBMAXE,YBMINE,YBMAXE,
     &              XBOFF,XBSF,YBOFF,YBSF)
        LOFINI = .FALSE.
        CALL GBPLOT(LPLTEL)
C------ save current color and set new color
        CALL GETCOLOR(ICOL0)
C
        CALL NEWCOLORNAME('green')
        CALL PLOT((XBF -XBOFF)*XBSF,(YBF -YBOFF)*YBSF,3)
        CALL PLOT((XTOP-XBOFF)*XBSF,(YTOP-YBOFF)*YBSF,2)
        CALL PLOT((XBF -XBOFF)*XBSF,(YBF -YBOFF)*YBSF,3)
        CALL PLOT((XBOT-XBOFF)*XBSF,(YBOT-YBOFF)*YBSF,2)
C
        COSD = COS(RDEF)
        SIND = SIN(RDEF)
        IF(ATOP .EQ. 0.0) THEN
         XBAR = XTOP - XBF
         YBAR = YTOP - YBF
         XT1C = XBF  +  XBAR*COSD + YBAR*SIND
         YT1C = YBF  -  XBAR*SIND + YBAR*COSD
         CALL PLOT((XBF -XBOFF)*XBSF,(YBF -YBOFF)*YBSF,3)
         CALL PLOT((XT1C-XBOFF)*XBSF,(YT1C-YBOFF)*YBSF,2)
        ENDIF
C
        IF(ABOT .EQ. 0.0) THEN
         XBAR = XBOT - XBF
         YBAR = YBOT - YBF
         XB1C = XBF  +  XBAR*COSD + YBAR*SIND
         YB1C = YBF  -  XBAR*SIND + YBAR*COSD
         CALL PLOT((XBF -XBOFF)*XBSF,(YBF -YBOFF)*YBSF,3)
         CALL PLOT((XB1C-XBOFF)*XBSF,(YB1C-YBOFF)*YBSF,2)
        ENDIF
C
        SHB = SHGT*XBSF
        CALL NEWCOLORNAME('red')
        CALL PLSYMB((XBF-XBOFF)*XBSF,(YBF-YBOFF)*YBSF,SHB,1,0.0,0)
C
        CALL NEWCOLOR(ICOL0)
        LGSAME = .FALSE.
C
        CALL PLFLUSH
c        COMAND = 'PLOT'
c        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'POIN'
     &  .OR. COMAND.EQ.'P   ') THEN
C
        CALL PNTMOD(IELGDES,LPLTEL)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'DIST') THEN
C------ distance between cursor clicks
        CALL DIST(XBOFF,YBOFF,XBSF,YBSF,SHB)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'PLOT') THEN
        CALL GPLINI(LOFINI,NBEL,LPLTEL,
     &              XBMINE,XBMAXE,YBMINE,YBMAXE,
     &              XBOFF,XBSF,YBOFF,YBSF)
        LOFINI = .FALSE.
        CALL GBPLOT(LPLTEL)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'BLOW'
     &  .OR. COMAND.EQ.'B   ') THEN
C------ blowup
        XWS = XWIND/SIZE
        YWS = YWIND/SIZE
        LSAME = .TRUE.
        LCURS = .TRUE.
        SHO = 2.0*XWS
        CALL OFFGET(XBOFF,YBOFF,XBSF,YBSF,XWS,YWS, LSAME,LCURS,SHO)
C
        LOFINI = .FALSE.
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'RESE'
     &  .OR. COMAND.EQ.'R   ') THEN
C------ reset default blowup parameters
        IF(IELPLT.EQ.0) THEN
         IEL1 = 1
         IEL2 = NBEL
        ELSE
         IEL1 = IELPLT
         IEL2 = IELPLT
        ENDIF
        DO IEL = IEL1, IEL2
          I = IBFRST(IEL)
          N = IBLAST(IEL) - IBFRST(IEL) + 1
          CALL MINMAX(N,XB(I),XBMINE(IEL),XBMAXE(IEL))
          CALL MINMAX(N,YB(I),YBMINE(IEL),YBMAXE(IEL))
        ENDDO
C
C------ go replot airfoil with initialization
        LOFINI = .TRUE.
        COMAND = 'PLOT'
        GO TO 510
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
        COMAND = 'PLOT'
        GO TO 510
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
      ELSEIF(COMAND.EQ.'TICK') THEN
C------ toggle tick mark plotting
        LGTICK = .NOT.LGTICK
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'GRID') THEN
C------ toggle grid plotting
        LGGRID = .NOT.LGGRID
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'GPAR') THEN
C------ toggle geometric parameter plotting
        LGPARM = .NOT.LGPARM
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'HPLO') THEN
C------ toggle geometric parameter plotting
        LHOMPL = .NOT.LHOMPL
        COMAND = 'PLOT'
        GO TO 510
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'OVER') THEN
C------ overlay disk airfoil file
        CALL OVER(COMARG)
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





      SUBROUTINE ZERCAM(IEL)
C-----------------------------------------
C     Zeros out camber of buffer airfoil
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION XBNEW(IBX), YBNEW(IBX)
C
      IF(LXBOD(IEL)) THEN
       WRITE(*,*) 'Axis body has zero camber by definition'
       RETURN
      ENDIF
C
      IF(.NOT.LBODY(IEL)) THEN
       WRITE(*,*) 'Auto y-symmetry is only for closed-body elements'
       RETURN
      ENDIF
C
      WRITE(*,*) 'Setting current camber to zero.'
C
      I = IBFRST(IEL)
      N = IBLAST(IEL) - IBFRST(IEL) + 1
C
C---- zero out camber of closed body
      TFAC = 1.0
      CFAC = 0.0
      CALL TCSET(XB(I),XBS(I),YB(I),YBS(I),SB(I),N, 
     &            TMAX,CMAX, TFAC,CFAC,XBNEW(I),YBNEW(I),.TRUE.)
      DO IB = IBFRST(IEL), IBLAST(IEL)
        XB(IB) = XBNEW(IB)
        YB(IB) = YBNEW(IB)
      ENDDO
      CALL SCALC(XB(I),YB(I),SB(I),N)
      CALL SEGSPL(XB(I),XBS(I),SB(I),N)
      CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
C---- make points exact mirror images
      CALL YSYM(XB(I),XBS(I),YB(I),YBS(I),SB(I),IBX,N,1,
     &          XBNEW(I),YBNEW(I),NNEW)
C
C---- move array slots for all higher elements to fit any new points
      NADD = NNEW - N
      CALL BSHIFT(IEL+1,IBFRST(IEL+1),NADD)
C
C---- insert new element possibly modified by YSYM
      IBLAST(IEL) = IBFRST(IEL) + NNEW - 1
      DO IB = IBFRST(IEL), IBLAST(IEL)
        XB(IB) = XBNEW(IB)
        YB(IB) = YBNEW(IB)
      ENDDO
      N = NNEW
C
      CALL SCALC(XB(I),YB(I),SB(I),N)
      CALL SEGSPL(XB(I),XBS(I),SB(I),N)
      CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
      LGSAME = .FALSE.
C
      RETURN
      END


      SUBROUTINE CAMS(IEL)
C-----------------------------------------
C     Adds to buffer airfoil camber line
C     another shape splined in x/c.
C-----------------------------------------
      INCLUDE 'DFDC.INC'
C
      PARAMETER (NIX = 300)
      DIMENSION XI(NIX), YI(NIX), YIX(NIX)
C
      PARAMETER (KK = 101)
      DIMENSION XPL(KK), YPL(KK)
      DIMENSION XPLAX(2), YPLAX(2)
      LOGICAL OK
C
      N = IBLAST(IEL) - IBFRST(IEL) + 1
      IF(N.GT.NIX) THEN
        WRITE(*,*) 'CAMS # input points exceeds NIX: ',N
        STOP
      ENDIF
C
    2 WRITE(*,1000)
 1000 FORMAT(/' Input x/c, y/c pairs for x/c = 0..1'
     &       /' Identical successive points enable a slope break')
C
      DO 5 I = 1, NIX
    3   READ(*,*,ERR=4) XI(I), YI(I)
        IM1 = MAX(1,I-1)
        IF(XI(I) .LT. XI(IM1)) GO TO 4
C
        IF(XI(I) .EQ. 1.0) GO TO 6
        GO TO 5
C
    4   WRITE(*,*) 'try again'
        GO TO 3
    5 CONTINUE
C
C---- spline camber line y(x)
    6 CONTINUE
      NI = I
      CALL SEGSPL(YI,YIX,XI,NI)
C
C
      CHXB = XBTE(IEL) - XBLE(IEL)
      CHYB = YBTE(IEL) - YBLE(IEL)
      CHBSQ = CHXB**2 + CHYB**2
C
C---- set new camber line for plotting and finding max
      XMAX = 0.0
      YMAX = 0.0
      DO K = 1, KK
        XOC = FLOAT(K-1)/FLOAT(KK-1)
        YOC = SEVAL(XOC,YI,YIX,XI,NI)
        IF(ABS(YOC) .GT. ABS(YMAX)) THEN
         XMAX = XOC
         YMAX = YOC
        ENDIF
        XPL(K) = XBLE(IEL) + CHXB*XOC - CHYB*YOC
        YPL(K) = YBLE(IEL) + CHXB*YOC + CHYB*XOC
      ENDDO
C
C---- define chord line
      XPLAX(1) = XPL(1)
      XPLAX(2) = XPL(KK)
      YPLAX(1) = YPL(1)
      YPLAX(2) = YPL(KK)
C
C
      CALL GETCOLOR(ICOL0)
C
C---- plot chord line, change in camber line
      CALL NEWPEN(1)
      CALL NEWCOLORNAME('magenta')
      CALL XYLINE(KK,XPL,YPL,XBOFF,XBSF,YBOFF,YBSF,2)
C
      CALL NEWCOLORNAME('blue')
      CALL XYLINE( 2,XPLAX,YPLAX,XBOFF,XBSF,YBOFF,YBSF,4)
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
C
C
      ALE = ATAN( DEVAL(0.0,YI,YIX,XI,NI) ) * 180.0/PI
      ATE = ATAN( DEVAL(1.0,YI,YIX,XI,NI) ) * 180.0/PI
C
      WRITE(*,1100) ALE, ATE, YMAX, XMAX
 1100 FORMAT(/' Added camber line incidence at LE =', F8.3, '  deg.',
     &       /' Added camber line incidence at TE =', F8.3, '  deg.',
     &       /' Max added camber/c =', F8.4, '  at x/c =', F7.3  )
C
      CALL ASKL('Is this acceptable?^',OK)
      IF(.NOT.OK) GO TO 2
C
C---- go over each point, changing the camber line appropriately
      IB1 = IBFRST(IEL)
      IB2 = IBLAST(IEL)
      DO IB = IB1, IB2
        XOC = (  (XB(IB)-XBLE(IEL))*CHXB
     &         + (YB(IB)-YBLE(IEL))*CHYB ) / CHBSQ
        TOC = SEVAL(XOC,YI,YIX,XI,NI)
        XB(IB) = XB(IB) - CHYB*TOC
        YB(IB) = YB(IB) + CHXB*TOC
      ENDDO
C
      I = IB1
      N = IB2 - IB1 + 1
      CALL SCALC(XB(I),YB(I),SB(I),N)
      CALL SEGSPL(XB(I),XBS(I),SB(I),N)
      CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
      LGSAME = .FALSE.
C
      RETURN
      END ! CAMS


      SUBROUTINE CAMP(IEL)
C-------------------------------------------------
C     Adds to buffer airfoil camber line another
C     shape derived from a specified delta(Cp).
C-------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      PARAMETER (NIX = 200)
      DIMENSION XI(NIX), DCPI(NIX), DCPIX(NIX)
C
      PARAMETER (KK = 101)
      DIMENSION XPL(KK), YPL(KK)
      DIMENSION XPLAX(2), YPLAX(2)
C
      LOGICAL OK
C
      DIMENSION X(KK), Y(KK), DYDX(KK), DCP(KK)
C
C---- singular part of camber line from finite delta(Cp) at LE and TE
      YSING(XX) = QOPI*DCP1*((XX-1.0)*LOG(MAX(1.0-XX,1.0E-6)) - XX)
     &          - QOPI*DCP0*( XX     *LOG(MAX(    XX,1.0E-6)) - XX)
C
      QOPI = 1.0/(4.0*PI)
C
C
    2 WRITE(*,1000)
 1000 FORMAT(/' Input x/c, delta(Cp) pairs for x/c = 0..1'
     &       /' Identical successive points enable a slope break')
C
      DO 5 I = 1, NIX
    3   READ(*,*,ERR=4) XI(I), DCPI(I)
C
        IM1 = MAX(1,I-1)
        IF(XI(I) .LT. XI(IM1)) GO TO 4
C
        IF(XI(I) .EQ. 1.0) GO TO 6
        GO TO 5
C
    4   WRITE(*,*) 'try again'
        GO TO 3
    5 CONTINUE
C
C---- spline dCp(x)
    6 CONTINUE
      NI = I
      CALL SEGSPL(DCPI,DCPIX,XI,NI)
C
C---- fill cosine-spaced delta(Cp) array from splined input points
      DO K=1, KK
        FRAC = FLOAT(K-1)/FLOAT(KK-1)
        X(K) = 0.5 - 0.5*COS(PI*FRAC)
        DCP(K) = SEVAL(X(K),DCPI,DCPIX,XI,NI)
      ENDDO
C
C---- ensure exact endpoint values
      X(1)  = 0.0
      X(KK) = 1.0
C
C
      DCP0 = DCP(1)
      DCP1 = DCP(KK)
C
      X0 = X(1)
      X1 = X(KK)
C         
C---- calculate Cauchy integral for y'(x) with removed singularity
      DO 10 K = 1, KK
        DYDX(K) = 0.0
C
        J = 1
        IF(K.EQ.J) THEN
         YP1 = DEVAL(X(J),DCPI,DCPIX,XI,NI)
        ELSE
         YP1 = (DCP(J) - DCP(K)) / (X(J) - X(K))
        ENDIF
C
        DO J = 2, KK
          IF(K.EQ.J) THEN
           YP2 = DEVAL(X(J),DCPI,DCPIX,XI,NI)
          ELSE
           YP2 = (DCP(J) - DCP(K)) / (X(J) - X(K))
          ENDIF
          DYDX(K) = DYDX(K) + 0.5*(YP1+YP2)*(X(J)-X(J-1))
C
          YP1 = YP2
        ENDDO
C
        DYDX(K) = QOPI*DYDX(K)
C
C------ add on removed part of Cauchy integral, further leaving out the
C-      possible infinities at LE,TE, so that y(x) can be safely splined. 
C-      The infinities are analytically integrated, and added on to y(x)
C-      with the statement function YSING.
        IF(K.NE.1) THEN
         DYDX(K) = DYDX(K)
     &           - QOPI*(DCP(K) - DCP0)*LOG(X(K) - X0)
        ENDIF
C
        IF(K.NE.KK) THEN
         DYDX(K) = DYDX(K)
     &           + QOPI*(DCP(K) - DCP1)*LOG(X1 - X(K))
        ENDIF
C
   10 CONTINUE
C
C---- integrate regular part of y'(x) from LE
      Y(1) = 0.
      DO K = 2, KK
        Y(K) = Y(K-1)
     &       + 0.5*(DYDX(K) + DYDX(K-1))*(X(K) - X(K-1))
      ENDDO
C
C---- add offset and angle of attack to get y(0) = y(1) = 0
      Y0 = Y(1)  + YSING(X(1) )
      Y1 = Y(KK) + YSING(X(KK))
      DO K = 1, KK
        Y(K) = Y(K)
     &       - Y0*(X1  -X(K))/(X1-X0)
     &       - Y1*(X(K)-X0  )/(X1-X0)
      ENDDO
C
      XMAX = 0.0
      YMAX = 0.0
      DO K = 1, KK
        IF(ABS(Y(K)) .GT. ABS(YMAX)) THEN
         XMAX = X(K)
         YMAX = Y(K)
        ENDIF
      ENDDO
C
      CHXB = XBTE(IEL) - XBLE(IEL)
      CHYB = YBTE(IEL) - YBLE(IEL)
      CHBSQ = CHXB**2 + CHYB**2
C
C---- set camber line in cartesian space
      DO K = 1, KK
        XT = X(K)
        YT = Y(K) + YSING(X(K))
        XPL(K) = XBLE(IEL) + CHXB*XT - CHYB*YT
        YPL(K) = YBLE(IEL) + CHYB*XT + CHXB*YT
      ENDDO
C
C---- define chord line
      XPLAX(1) = X(1)
      XPLAX(2) = X(KK)
      YPLAX(1) = Y(1)
      YPLAX(2) = Y(KK)
C
      CALL GETCOLOR(ICOL0)
C
C---- plot chord line, change in camber line
      CALL NEWPEN(1)
      CALL NEWCOLORNAME('magenta')
      CALL XYLINE(KK,XPL,YPL,XBOFF,XBSF,YBOFF,YBSF,2)
      CALL NEWCOLORNAME('blue')
      CALL XYLINE( 2,XPLAX,YPLAX,XBOFF,XBSF,YBOFF,YBSF,4)
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
C
C
      WRITE(*,1100) YMAX, XMAX
 1100 FORMAT(/' Max added camber/c =', F8.4, '  at x/c =', F7.3  )
C
      CALL ASKL('Is this acceptable?^',OK)
      IF(.NOT.OK) GO TO 2
C
C---- spline regular part of y(x)
      CALL SEGSPL(Y,DYDX,X,KK)
C
C---- go over each point, changing the camber line appropriately
      IB1 = IBFRST(IEL)
      IB2 = IBLAST(IEL)
      DO IB = IB1, IB2
        XOC = (  (XB(IB)-XBLE(IEL))*CHXB
     &         + (YB(IB)-YBLE(IEL))*CHYB ) / CHBSQ
        YOC = SEVAL(XOC,Y,DYDX,X,KK) + YSING(XOC)
        XB(IB) = XB(IB) - CHYB*YOC
        YB(IB) = YB(IB) + CHXB*YOC
      ENDDO
C
      I = IB1
      N = IB2 - IB1 + 1
      CALL SCALC(XB(I),YB(I),SB(I),N)
      CALL SEGSPL(XB(I),XBS(I),SB(I),N)
      CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
      LGSAME = .FALSE.
C
      RETURN
      END ! CAMP



      SUBROUTINE OVER(FNAME1)
C----------------------------------------------------
C     Overlays plot of airfoil from coordinate file.
C----------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*) FNAME1
C
      CHARACTER*80 FNAME
C
C---- local arrays for calling AREAD
      PARAMETER (ITX=300)
      DIMENSION XT(ITX,NEX), YT(ITX,NEX)
      DIMENSION ST(ITX), XTS(ITX), YTS(ITX)
      DIMENSION NT(NEX)
C
      CHARACTER*32 ONAME
      CHARACTER*80 OPARS
C
      DATA KK / 20 /
C
      FNAME = FNAME1
      IF(FNAME(1:1).EQ.' ') CALL ASKS('Enter filename^',FNAME)
C
      LU = 7
      CALL AREAD(FNAME,LU, ITX,NEX, XT,YT,
     &           NT,NTEL,
     &           ONAME,OPARS,IOTYPE)
C
      IF(IOTYPE.EQ.0) THEN
C----- read error occurred
       RETURN
      ENDIF
      CALL STRIP(ONAME,NONAME)
C
      CALL GETCOLOR(ICOL0)
C
      CALL NEWPEN(2)
      CALL NEWCOLORNAME('cyan')
      DO IEL = 1, NTEL
        N = NT(IEL)
        CALL SCALC(XT(1,IEL),YT(1,IEL),ST,N)
        CALL SEGSPL(XT(1,IEL),XTS,ST,N)
        CALL SEGSPL(YT(1,IEL),YTS,ST,N)
C
        CALL PLTGSP(XT(1,IEL),XTS,YT(1,IEL),YTS,ST,N, 
     &              XBOFF,XBSF,YBOFF,YBSF, KK)
      ENDDO
      CALL NEWCOLOR(ICOL0)
      CALL PLFLUSH
C
      RETURN
      END




      SUBROUTINE CLRHOM(IEL)
      INCLUDE 'DFDC.INC'
C
      DXBSUM(IEL) = 0.
      DYBSUM(IEL) = 0.
      AGBSUM(IEL) = 0.
      XFBSUM(IEL) = 1.0
      YFBSUM(IEL) = 1.0
C
cc      XBREFE(IEL) = 0.
cc      YBREFE(IEL) = 0.
C
      RETURN
      END ! CLRHOM


      SUBROUTINE TRAN(IELOPT,DX,DY)
C-------------------------------------------
C     Translates element IEL 
C     elements by distance DX,DY.
C-------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LDONE(NEX)
C
      DO IEL = 1, NBEL
        LDONE(IEL) = .FALSE.
      ENDDO
C
      IEL = IELOPT
C
 10   IF(LDONE(IEL)) RETURN
C
cC---- send changes to panel routines
c      CALL ADDXYA(IEL,DX,DY,0.0,0.0)
C
      IB1 = IBFRST(IEL)
      IB2 = IBLAST(IEL)
C
C---- shift element nodes
      DO IB = IB1, IB2
        XB(IB) = XB(IB) + DX
        YB(IB) = YB(IB) + DY
      ENDDO
C
C---- shift other associated points
      XBLE(IEL) = XBLE(IEL) + DX
      YBLE(IEL) = YBLE(IEL) + DY
      XBTE(IEL) = XBTE(IEL) + DX
      YBTE(IEL) = YBTE(IEL) + DY
      XBREFE(IEL) = XBREFE(IEL) + DX
      YBREFE(IEL) = YBREFE(IEL) + DY
      XBMINE(IEL) = XBMINE(IEL) + DX
      YBMINE(IEL) = YBMINE(IEL) + DY
      XBMAXE(IEL) = XBMAXE(IEL) + DX
      YBMAXE(IEL) = YBMAXE(IEL) + DY
      XELNUM(IEL) = XELNUM(IEL) + DX
      YELNUM(IEL) = YELNUM(IEL) + DY
C
      XBCEN2DA(IEL) = XBCEN2DA(IEL) + DX
      YBCEN2DA(IEL) = YBCEN2DA(IEL) + DY
      XBCEN2DT(IEL) = XBCEN2DT(IEL) + DX
      YBCEN2DT(IEL) = YBCEN2DT(IEL) + DY
      XBCENV(IEL)   = XBCENV(IEL) + DX
      YBCENV(IEL)   = YBCENV(IEL) + DY
      XBCENVT(IEL)  = XBCENVT(IEL) + DX
      YBCENVT(IEL)  = YBCENVT(IEL) + DY
C
C---- accumulate shifts
      DXBSUM(IEL) = DXBSUM(IEL) + DX
      DYBSUM(IEL) = DYBSUM(IEL) + DY
C
cc    LHOME(IEL) = .FALSE.
      LDONE(IEL) = .TRUE.
C
      RETURN
      END ! TRAN


      SUBROUTINE SCAL(IELOPT,XSCL,YSCL,XCT,YCT)
C------------------------------------------------
C     Scales element NN 
C     by XSCL,YSCL, about point XCT,YCT.
C------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LDONE(NEX)
C
      DO IEL = 1, NEX
        LDONE(IEL) = .FALSE.
      ENDDO
C
      XCEN = XCT
      YCEN = YCT
C
      IEL = IELOPT
C
 10   IF(LDONE(IEL)) RETURN
C
C---- set changes about origin from changes about reference point
      DF = 0.5*LOG(ABS(XSCL*YSCL))
      DA = 0.0
      DX = (1.0-XSCL)*XCEN
      DY = (1.0-YSCL)*YCEN
C
cC---- send changes to panel routines
c      CALL ADDXYA(IEL,DX,DY,DF,DA)
C
      IB1 = IBFRST(IEL)
      IB2 = IBLAST(IEL)
C
C---- scale element nodes
      DO IB = IB1, IB2
        XB(IB) = XSCL*XB(IB) + DX
        YB(IB) = YSCL*YB(IB) + DY
      ENDDO
C
C---- scale other associated points
      XBLE(IEL) = XSCL*XBLE(IEL) + DX
      YBLE(IEL) = YSCL*YBLE(IEL) + DY
      XBTE(IEL) = XSCL*XBTE(IEL) + DX
      YBTE(IEL) = YSCL*YBTE(IEL) + DY
      XBMINE(IEL) = XSCL*XBMINE(IEL) + DX
      YBMINE(IEL) = YSCL*YBMINE(IEL) + DY
      XBMAXE(IEL) = XSCL*XBMAXE(IEL) + DX
      YBMAXE(IEL) = YSCL*YBMAXE(IEL) + DY
      XBREFE(IEL) = XSCL*XBREFE(IEL) + DX
      YBREFE(IEL) = YSCL*YBREFE(IEL) + DY
      XELNUM(IEL) = XSCL*XELNUM(IEL) + DX
      YELNUM(IEL) = YSCL*YELNUM(IEL) + DY
C
      XBCEN2DA(IEL) = XSCL*XBCEN2DA(IEL) + DX
      YBCEN2DA(IEL) = YSCL*YBCEN2DA(IEL) + DY
      XBCEN2DT(IEL) = XSCL*XBCEN2DT(IEL) + DX
      YBCEN2DT(IEL) = YSCL*YBCEN2DT(IEL) + DY
      XBCENV(IEL)   = XSCL*XBCENV(IEL)   + DX
      YBCENV(IEL)   = YSCL*YBCENV(IEL)   + DY
      XBCENVT(IEL)  = XSCL*XBCENVT(IEL)  + DX
      YBCENVT(IEL)  = YSCL*YBCENVT(IEL)  + DY
C
C---- add to total shift assumulators
      DXBSUM(IEL) = XSCL*DXBSUM(IEL) + DX
      DYBSUM(IEL) = YSCL*DYBSUM(IEL) + DY
C
      IF(LBODY(IEL) .AND. XSCL*YSCL .LT. 0.0) THEN
C----- closed element was flipped, so reverse node ordering
       DO IB = IB1, (IB2+IB1)/2
         XTMP = XB(IB2-IB+IB1)
         YTMP = YB(IB2-IB+IB1)
         XB(IB2-IB+IB1) = XB(IB)
         YB(IB2-IB+IB1) = YB(IB)
         XB(IB) = XTMP
         YB(IB) = YTMP
       ENDDO
      ENDIF
C
      I = IB1
      N = IB2 - IB1 + 1
      CALL SCALC(XB(I),YB(I),SB(I),N)
      CALL SEGSPL(XB(I),XBS(I),SB(I),N)
      CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
C---- accumulate scales
      XFBSUM(IEL) = XFBSUM(IEL)*XSCL
      YFBSUM(IEL) = YFBSUM(IEL)*YSCL
C
cc    LHOME(IEL) = .FALSE.
      LDONE(IEL) = .TRUE.
C
      RETURN
      END ! SCAL

 
      SUBROUTINE ANGL(IELOPT,ADEG,XCT,YCT)
C------------------------------------------------
C     Rotates element IEL 
C     by ADEG (degrees) about the point XCT,YCT.
C------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LDONE(NEX)
C
      DO IEL = 1, NBEL
        LDONE(IEL) = .FALSE.
      ENDDO
C
      XCEN = XCT
      YCEN = YCT
C
      ARAD = ADEG * PI/180.0
      SA = SIN(ARAD)
      CA = COS(ARAD)
C
      IEL = IELOPT
C
 10   IF(LDONE(IEL)) RETURN
C
C---- set changes about origin from changes about reference point
      DF = 0.0
      DA = ARAD
      DX = XCEN - (CA*XCEN + SA*YCEN)
      DY = YCEN - (CA*YCEN - SA*XCEN)
C
cC---- send changes to panel routines
c      CALL ADDXYA(N,DX,DY,DF,DA)
C
      IB1 = IBFRST(IEL)
      IB2 = IBLAST(IEL)
C
C---- rotate element nodes
      DO IB = IB1, IB2
        XX = XB(IB)
        YY = YB(IB)
        XB(IB) = CA*XX + SA*YY + DX
        YB(IB) = CA*YY - SA*XX + DY
      ENDDO
C
C---- rotate other associated points
      XX = XBLE(IEL)
      YY = YBLE(IEL)
      XBLE(IEL) = CA*XX + SA*YY + DX
      YBLE(IEL) = CA*YY - SA*XX + DY
C
      XX = XBTE(IEL)
      YY = YBTE(IEL)
      XBTE(IEL) = CA*XX + SA*YY + DX
      YBTE(IEL) = CA*YY - SA*XX + DY
C
      XX = XELNUM(IEL)
      YY = YELNUM(IEL)
      XELNUM(IEL) = CA*XX + SA*YY + DX
      YELNUM(IEL) = CA*YY - SA*XX + DY
C
      XX = XBREFE(IEL)
      YY = YBREFE(IEL)
      XBREFE(IEL) = CA*XX + SA*YY + DX
      YBREFE(IEL) = CA*YY - SA*XX + DY
C
      XX = XBCEN2DA(IEL)
      YY = YBCEN2DA(IEL)
      XBCEN2DA(IEL) = CA*XX + SA*YY + DX
      YBCEN2DA(IEL) = CA*YY - SA*XX + DY
C
      XX = XBCEN2DT(IEL)
      YY = YBCEN2DT(IEL)
      XBCEN2DT(IEL) = CA*XX + SA*YY + DX
      YBCEN2DT(IEL) = CA*YY - SA*XX + DY
C
      XX = DXBSUM(IEL)
      YY = DYBSUM(IEL)
      DXBSUM(IEL) = CA*XX + SA*YY + DX
      DYBSUM(IEL) = CA*YY - SA*XX + DY
C
      I = IB1
      N = IB2 - IB1 + 1
      CALL SEGSPL(XB(I),XBS(I),SB(I),N)
      CALL SEGSPL(YB(I),YBS(I),SB(I),N)
C
C---- accumulate rotation
      AGBSUM(IEL) = AGBSUM(IEL) + ADEG
C
cc    LHOME(IEL) = .FALSE.
      LDONE(IEL) = .TRUE.
C
      RETURN
      END ! ANGL
 
