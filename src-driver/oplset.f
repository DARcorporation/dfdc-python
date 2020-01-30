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

      SUBROUTINE OPLSET
      INCLUDE 'PLOT.INC'
C--------------------------------------------------------------
C     Allows user modification of plot parameters in PLOT.INC
C--------------------------------------------------------------
      CHARACTER*1 VAR
      CHARACTER*4 COMAND
      CHARACTER*128 COMARG
      CHARACTER*10 CHCURS, CHLAND
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR, LCOLOR
C
 1000 FORMAT(A)
C
 1    CONTINUE
      IF(LCURS) THEN
       CHCURS = 'Cursor    '
      ELSE
       CHCURS = 'Keyboard  '
      ENDIF
C
      IF(LLAND) THEN
       CHLAND = 'Landscape '
      ELSE
       CHLAND = 'Portrait  '
      ENDIF
C
      LCOLOR = IDEVPS.EQ.4
C
      WRITE(*,2000) SIZE, XPAGE,YPAGE, XMARG,YMARG, 
     &              CHGT, SHGT, SCRNFR,
     &              CHCURS, CHLAND, LCOLOR
 2000 FORMAT(' ...............................................'
     &     //'  S ize of plot object         ', F6.2,'"'
     &      /'  D imensions of page          ', F6.2,' x',F6.2,'"'
     &      /'  M argins from page edges     ', F6.2,'",',F6.2,'"'
     &      /'  F ont size (relative)        ', F8.4
     &      /'  P oint-symbol size (relative)', F8.4
     &      /'  W indow/screen size fraction ', F8.4
     &      /'  B lowup input method:        ', A 
     &      /'  O rientation of plot:        ', A 
     &      /'  C olor PostScript output?    ', L2 )
C
C   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
C   x x x x   x             x   x x     x       x
C
 5    CALL ASKC('      Option, Value   (or <Return>) ^',COMAND,COMARG)
C
      DO I=1, 20
        IINPUT(I) = 0.0
        RINPUT(I) = 0.0
      ENDDO
      NINPUT = 0
      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
      NINPUT = 0
      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C
      VAR = COMAND(1:1)
      IF (VAR.EQ.'0' .OR. VAR.EQ.' ') THEN
       RETURN
C
      ELSEIF (INDEX('Ss',VAR).NE.0) THEN
        IF(NINPUT.GE.1) THEN
          SIZE = RINPUT(1)
        ELSE
          CALL ASKR('Enter size (in)^',SIZE)
        ENDIF
C
      ELSEIF (INDEX('Dd',VAR).NE.0) THEN
        IF(NINPUT.GE.2) THEN
          XPAGE = RINPUT(1)
          YPAGE = RINPUT(2)
        ELSEIF(NINPUT.GE.1) THEN
          XPAGE = RINPUT(1)
          CALL ASKR('Enter page Y dimension (in)^',YPAGE)
        ELSE
          CALL ASKR('Enter page X dimension (in)^',XPAGE)
          CALL ASKR('Enter page Y dimension (in)^',YPAGE)
        ENDIF
C
      ELSEIF (INDEX('Mm',VAR).NE.0) THEN
        IF(NINPUT.GE.2) THEN
          XMARG = RINPUT(1)
          YMARG = RINPUT(2)
        ELSEIF(NINPUT.GE.1) THEN
          XMARG = RINPUT(1)
          CALL ASKR('Enter page Y margin (in)^',YMARG)
        ELSE
          CALL ASKR('Enter page X margin (in)^',XMARG)
          CALL ASKR('Enter page Y margin (in)^',YMARG)
        ENDIF
C
      ELSEIF (INDEX('Ff',VAR).NE.0) THEN
        IF(NINPUT.GE.1) THEN
          CHGT = RINPUT(1)
        ELSE
          CALL ASKR('Enter character font size^',CHGT)
        ENDIF
C
      ELSEIF (INDEX('Pp',VAR).NE.0) THEN
        IF(NINPUT.GE.1) THEN
          SHGT = RINPUT(1)
        ELSE
          CALL ASKR('Enter point-symbol size^',SHGT)
        ENDIF
C
      ELSEIF (INDEX('Ww',VAR).NE.0) THEN
        IF(NINPUT.GE.1) THEN
          SCRNFR = RINPUT(1)
        ELSE
          CALL ASKR('Enter window/screen size fraction^',SCRNFR)
        ENDIF
C
      ELSEIF (INDEX('Bb',VAR).NE.0) THEN
        LCURS = .NOT. LCURS
C
      ELSEIF (INDEX('Oo',VAR).NE.0) THEN
        LLAND = .NOT. LLAND
        WRITE(*,*)
        WRITE(*,*) 'Swapping X,Y page dimensions'
        XTMP = XPAGE
        YTMP = YPAGE
        XPAGE = YTMP
        YPAGE = XTMP
C
      ELSEIF (INDEX('Cc',VAR).NE.0) THEN
        LCOLOR = .NOT. LCOLOR
        IF(     LCOLOR) IDEVPS = 4
        IF(.NOT.LCOLOR) IDEVPS = 2
C
      ELSE
        WRITE(*,*) '*** Item not recognized ***'
      ENDIF
      GO TO 1
C
      END ! OPLSET
