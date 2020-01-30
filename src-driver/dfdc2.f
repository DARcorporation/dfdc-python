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


      SUBROUTINE WRGNULIC
C---- Print out a copy of the GNU license for DFDC
      WRITE(*,10)
 10   FORMAT(
     &/'==============================================================='
     &,'===============',
     &/' DFDC (Ducted Fan Design Code) is an aerodynamic and aeroacoust'
     &,'ic design',
     &/' and analysis tool for aircraft with propulsors in ducted fan'
     &,' configurations.',
     &/' ',
     &/' This software was developed under the auspices and sponsorship'
     &,' of the',
     &/' Tactical Technology Office (TTO) of the Defense Advanced Resea'
     &,'rch',
     &/' Projects Agency (DARPA).',
     &/' ',
     &/' Copyright (c) 2004, 2005, Booz Allen Hamilton Inc., All Rights'
     &,' Reserved'
     &/' ',
     &/' This program is free software; you can redistribute it and/or '
     &,'modify it',
     &/' under the terms of the GNU General Public License as published'
     &,' by the',
     &/' Free Software Foundation; either version 2 of the License, or '
     &,'(at your',
     &/' option) any later version.',
     &/' ',
     &/' This program is distributed in the hope that it will be useful'
     &,', but',
     &/' WITHOUT ANY WARRANTY; without even the implied warranty of'
     &/' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the '
     &,'GNU',
     &/' General Public License for more details.',
     &/' ',
     &/' You should have received a copy of the GNU General Public Lice'
     &,'nse along',
     &/' with this program; if not, write to the Free Software Foundati'
     &,'on, Inc.,',
     &/' 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA',
     &/' ',
     &/' Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (dre'
     &,'la@mit.edu)',
     &/' Program Management: Brad Tousley, Paul Eremenko (eremenko@alum'
     &,'.mit.edu)',
     &/'==============================================================='
     &,'===============')
      RETURN
      END



      SUBROUTINE ELTYPE(IELDEF)
C--------------------------------------
C     Toggles solid/wake element type
C--------------------------------------
      INCLUDE 'DFDC.INC'
C
      CHARACTER*80 LINE
      LOGICAL ERROR, CHANGE
      DIMENSION IINPUT(10)
C
 1000 FORMAT(A)
C
      CHANGE = .FALSE.
C
      IF(IELDEF.GT.0) THEN
C----- use nonzero argument as the first input and jump to process it
       NINPUT = 1
       IINPUT(1) = IELDEF
       GO TO 20
      ENDIF
C
 10   CONTINUE
      WRITE(*,*)
      WRITE(*,*) '   Surface types...'
      WRITE(*,*) '   ------------------------------------------'
      DO IEL = 1, NEL
        IF    (NETYPE(IEL).EQ.0) THEN
         WRITE(*,2100) IEL
 2100    FORMAT('     element',I3, ':  strength unknown (solid)')
        ELSEIF(NETYPE(IEL).EQ.1) THEN
         WRITE(*,2110) IEL
 2110    FORMAT('     element',I3, ':  strength prescribed')
        ENDIF
      ENDDO
C
      WRITE(*,*)
 12   WRITE(*,3000) 
 3000 FORMAT('    Enter elements to toggle:  ',$)
      READ(*,1000) LINE
      NINPUT = 10
      CALL GETINT(LINE,IINPUT,NINPUT,ERROR)
      IF(ERROR) GO TO 12
      IF(NINPUT.EQ.0) THEN
       IF(CHANGE) THEN
        LSYSP = .FALSE.
        LGSYS = .FALSE.
        LGAMU = .FALSE.
        LGAMA = .FALSE.
        LQCNT = .FALSE.
C
C------ see if TE panels need to be removed and AIC matrices recomputed
        DO IEL = 1, NEL
          IF(NETYPE(IEL).EQ.1 .AND. LTPAN(IEL)) THEN
           LTPAN(IEL) = .FALSE.
           LQAIC = .FALSE.
           LQGIC = .FALSE.
          ENDIF
        ENDDO
       ENDIF
       RETURN
      ENDIF
C
 20   CONTINUE
      DO I = 1, NINPUT
        IEL = IINPUT(I)
        IF(IEL.LT.0 .OR. IEL.GT.NEL) THEN
         WRITE(*,*) 'Element index out of range'
        ELSEIF(NETYPE(IEL).EQ.0) THEN
         NETYPE(IEL) = 1
         CHANGE = .TRUE.
        ELSEIF(NETYPE(IEL).EQ.1) THEN
         NETYPE(IEL) = 0
         CHANGE = .TRUE.
        ELSE
         WRITE(*,*) 'Cannot change line or point element-type'
        ENDIF
      ENDDO
C
C---- go back to menu
      GO TO 10
      END ! ELTYPE



      SUBROUTINE ELTPAN(IELDEF)
C--------------------------------------
C     Toggles TE panel inclusion
C--------------------------------------
      INCLUDE 'DFDC.INC'
C
      CHARACTER*80 LINE
      LOGICAL ERROR, CHANGE
      DIMENSION IINPUT(10)
C
 1000 FORMAT(A)
C
      CHANGE = .FALSE.
C
      IF(IELDEF.GT.0) THEN
C----- use nonzero argument as the first input and jump to process it
       NINPUT = 1
       IINPUT(1) = IELDEF
       GO TO 20
      ENDIF
C
 10   CONTINUE
      WRITE(*,*)
      WRITE(*,*) '   TE panel presence...'
      WRITE(*,*) '   ------------------------------'
      DO IEL = 1, NEL
        IF(NETYPE(IEL).EQ.0 .OR. NETYPE(IEL).EQ.7) THEN
          WRITE(*,1105) IEL, LTPAN(IEL)
        ENDIF
      ENDDO
 1105 FORMAT( '     element',I3,':    ', L1)
C
      WRITE(*,*)
 12   WRITE(*,1200) 
 1200 FORMAT('    Enter elements to toggle:  ',$)
      READ(*,1000) LINE
      NINPUT = 10
      CALL GETINT(LINE,IINPUT,NINPUT,ERROR)
      IF(ERROR) GO TO 12
      IF(NINPUT.EQ.0) THEN
       IF(CHANGE) THEN
cccHHY
        LQAIC = .FALSE.
cccHHY
        LGSYS = .FALSE.
        LGAMU = .FALSE.
        LGAMA = .FALSE.
       ENDIF
       RETURN
      ENDIF
C
 20   CONTINUE
      DO I = 1, NINPUT
        IEL = IINPUT(I)
        IF(IEL.LE.0 .OR. IEL.GT.NEL) THEN
         WRITE(*,*) 'Element index out of range'
        ELSEIF(NETYPE(IEL).NE.0 .AND. NETYPE(IEL).NE.7) THEN
         WRITE(*,*) 'Must be solid body or vortex wake panel'
        ELSEIF(LXBOD(IEL)) THEN
         WRITE(*,*) 'Closed axis-body elements cannot have TE panels'
        ELSE
         LTPAN(IEL) = .NOT.LTPAN(IEL)
         CHANGE = .TRUE.
        ENDIF
      ENDDO
C
C---- go back to menu
      GO TO 10
      END ! ELTPAN


      SUBROUTINE ELPLSD
C--------------------------------------
C     Toggles element-side plot flags
C--------------------------------------
      INCLUDE 'DFDC.INC'
C
      CHARACTER*80 LINE
      LOGICAL ERROR
C
      CHARACTER*3 CSPLOT(-1:2)
      DATA CSPLOT /  '---', 'L R', '  R', 'L  ' /
C
 1000 FORMAT(A)
C
 10   CONTINUE
      WRITE(*,*)
      WRITE(*,*) '   Element-side plot toggles...'
      WRITE(*,*) '   ------------------------------'
      DO IEL = 1, NEL
        IF( NETYPE(IEL).EQ.0 
     &  .OR.NETYPE(IEL).EQ.1
     &  .OR.NETYPE(IEL).EQ.2
     &  .OR.NETYPE(IEL).EQ.5
     &  .OR.NETYPE(IEL).EQ.6
     &  .OR.NETYPE(IEL).EQ.7 ) THEN
         WRITE(*,1105) IEL, CSPLOT(ISPLOT(IEL))
        ENDIF
      ENDDO
 1105 FORMAT( '     element',I3,':    ', A)
C
      WRITE(*,*)
 12   WRITE(*,1200)
 1200 FORMAT('    Enter  element, plotted-sides (L,R): ', $)
      READ(*,1000) LINE
      IF(LINE.EQ.' ') THEN
       RETURN
      ENDIF
C
      READ(LINE,*,ERR=12) IEL
      IF(IEL.LT.1 .OR. IEL.GT.NEL) THEN
       WRITE(*,*) 'Element index out of range'
      ELSE
       KR = INDEX(LINE,'R') + INDEX(LINE,'r')
       KL = INDEX(LINE,'L') + INDEX(LINE,'l')
C
       IF    (KR.GT.0 .AND. KL.GT.0) THEN
        ISPLOT(IEL) = 0
       ELSEIF(KR.GT.0) THEN
        ISPLOT(IEL) = 1
       ELSEIF(KL.GT.0) THEN
        ISPLOT(IEL) = 2
       ELSE
        ISPLOT(IEL) = -1
       ENDIF
      ENDIF
C
      GO TO 10
      END ! ELPLSD


