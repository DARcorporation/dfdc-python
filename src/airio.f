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
C     Trivial reporting format adjustment (line 275)
C
C=========================================================================
C
C
      SUBROUTINE AREAD(FNAME,LU,IBX,NBX,XB,YB,
     &                 NB,NBL,
     &                 NAME,ISPARS,IFTYPE)
      CHARACTER*(*) FNAME
      CHARACTER*(*) NAME
      CHARACTER*(*) ISPARS
      DIMENSION XB(IBX,NBX) ,YB(IBX,NBX)
      DIMENSION NB(NBX)
C--------------------------------------------------------
C     Reads in several types of airfoil coordinate file.
C
C  Input:
C       FNAME     name of coordinate file to be read
C                   if blank, unit LU is assumed to be open
C       LU        logical unit for file read
C
C  Output:
C       XB,YB(..)  coordinates
C       NB(.)      number of XB,YB coordinates
C
C       NAME    character name string        (if IFTYPE > 1)
C       ISPARS  ISES/MSES domain-size string (if IFTYPE > 2)
C       IFTYPE  returns type of file:
C            0    None.  Read error occurred.
C            1    Generic.
C            2    Labeled generic.
C            3    MSES single element.
C            4    MSES multi-element.
C--------------------------------------------------------
      CHARACTER*128 LINE
      LOGICAL LOPEN, ERROR
C
      DIMENSION AINPUT(10)
C
C---- permitted characters in string containing Fortran-readable numbers
      CHARACTER*20 NCHARS
      DATA NCHARS / '0123456789-+.,EDed  ' /
C
C---- tab character is also permitted
      NCHARS(20:20) = CHAR(9)
C
C
C---- first assume that there will be a read error
      IFTYPE = 0
C
      LOPEN = FNAME(1:1) .NE. ' '
      IF(LOPEN) THEN
       NF = INDEX(FNAME,' ') + 1
       OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=98)
       REWIND(LU)
      ENDIF
C
      READ(LU,1000) LINE
C
C---- if first line has any non-numeric character, go treat it as the name
      DO K=1, 80
        IF(INDEX(NCHARS,LINE(K:K)) .EQ. 0) GO TO 20
      ENDDO
C
C---- plain unlabeled file: rewind, and just read in X,Y coordinates
      NAME = ' '
C
      REWIND(LU)
      WRITE(*,*)
      WRITE(*,*) 'Reading plain coordinate file'
      IFTYPE = 1
      GO TO 40
C
C---- first line interpreted as label string
   20 CONTINUE
      NAME = LINE
C
C---- read second line
      READ(LU,1000) LINE
      NINPUT = 10
      CALL GETFLT(LINE,AINPUT,NINPUT,ERROR)
      IF(ERROR) GO TO 99
C
      IF(NINPUT.LT.4) THEN
C------ no domain parameters: re-read name string and then read X,Y coordinates
        REWIND(LU)
        READ(LU,1000) NAME
        IFTYPE = 2
      ENDIF
C
      WRITE(*,1010) NAME
 1010 FORMAT(/1X,'Reading airfoil: ', A32)
C
 40   CONTINUE
C
C---- read in airfoil coordinates
      DO 55 N=1, NBX+1
        DO IB=1, IBX+1
          READ(LU,*,END=56,ERR=99) XBT, YBT
C
          IF(XBT.EQ.999.0) THEN
C--------- save number of points and type for element which was just read in
           NB(N) = IB-1
C--------- go try to read next element
           GO TO 55
          ENDIF
C
          IF(N.GT.NBX) THEN
           WRITE(*,*) 'AREAD: Too many airfoil elements. Increase NBX.'
           STOP
          ENDIF
C
C-------- uneventful XBT,YBT coordinates... just store them
          XB(IB,N) = XBT
          YB(IB,N) = YBT
C
        END DO
        WRITE(*,*) 'AREAD: Too many airfoil points.  Increase IBX.'
        STOP
C
   55 CONTINUE
      N = NBX
C
   56 IF(IB.EQ.1) THEN
C----- coordinate file has "999.0 999.0" at the end ...
       NBL = N-1
      ELSE
C----- coordinate file has no special ending line
       NBL = N
       NB(N) = IB-1
      ENDIF
C
C---- set type if coordinate file this was
      IF(IFTYPE.EQ.0) THEN
       IF(NBL.EQ.1) THEN
        IFTYPE = 3
       ELSE
        IFTYPE = 4
       ENDIF
      ENDIF
C
      IF(LOPEN) CLOSE(LU)
      RETURN
C
   98 CONTINUE
      WRITE(*,1050) FNAME(1:NF)
      IFTYPE = 0
      RETURN
C
   99 CONTINUE
      WRITE(*,1100) FNAME(1:NF)
      IFTYPE = 0
      IF(LOPEN) CLOSE(LU)
      RETURN
C...............................................................
 1000 FORMAT(A)
 1050 FORMAT(/' File OPEN error:  ', A)
 1100 FORMAT(/' File READ error:  ', A)
      END ! AREAD



      SUBROUTINE AREADNR(LU,IBX,NBX,XB,YB,
     &                   NB,NBL,
     &                   NAME,ISPARS,IFTYPE)
      CHARACTER*80 NAME
      CHARACTER*(*) ISPARS
      DIMENSION XB(IBX,NBX) ,YB(IBX,NBX)
      DIMENSION NB(NBX)
C--------------------------------------------------------
C     Reads in several types of airfoil coordinate file
C     without performing file open or rewinds.  
C     Assumes input airfoil file is open on entry
C
C  Input:
C       LU        logical unit for file read
C
C  Output:
C       XB,YB(..)  coordinates
C       NB(.)      number of XB,YB coordinates
C
C       NAME    character name string        (if IFTYPE > 1)
C       ISPARS  ISES/MSES domain-size string (if IFTYPE > 2)
C       IFTYPE  returns type of file:
C            0    None.  Read error occurred.
C            1    Generic.
C            2    Labeled generic.
C            3    MSES single element.
C            4    MSES multi-element.
C--------------------------------------------------------
      CHARACTER*128 LINE
      LOGICAL ERROR
C
      DIMENSION AINPUT(10)
C
C---- permitted characters in string containing Fortran-readable numbers
      CHARACTER*20 NCHARS
      DATA NCHARS / '0123456789-+.,EDed  ' /
C
C---- tab character is also permitted
      NCHARS(20:20) = CHAR(9)
C
C
C---- first assume that there will be a read error
      IFTYPE = 0
      READ(LU,1000) LINE
C---- if first line has any non-numeric character, go treat it as the name
      DO K=1, 80
        IF(INDEX(NCHARS,LINE(K:K)) .EQ. 0) GO TO 20
      ENDDO
C
C---- plain unlabeled file: rewind, and just read in X,Y coordinates
      NAME = ' '
      WRITE(*,*)
      WRITE(*,*) 'Reading plain coordinate file'
      IFTYPE = 1
      GO TO 40
C
C---- first line interpreted as label string
   20 CONTINUE
      NAME = LINE(1:60)
C      CALL STRIP(NAME,NNM)
C
C---- read second line
      READ(LU,1000) LINE
      NINPUT = 10
      CALL GETFLT(LINE,AINPUT,NINPUT,ERROR)
      IF(ERROR) GO TO 99
C
      IF(NINPUT.LT.4) THEN
C------ no domain parameters: read X,Y coordinates directly with this line
        IFTYPE = 2
      ELSE
C---- second line had parameters, read third line (coordinates line)
        READ(LU,1000) LINE
      ENDIF
C
      WRITE(*,1010) NAME
 1010 FORMAT(/1X,'Reading airfoil: ', A32)
C
 40   CONTINUE
C
C---- read in airfoil coordinates assuming first line has been read into
C     character variable LINE
      DO 55 N=1, NBX+1
        DO IB=1, IBX+1
          IF(N.EQ.1 .AND. IB.EQ.1) THEN
            READ(LINE,*,END=56,ERR=99) XBT, YBT
          ELSE
            READ(LU,1000,END=56,ERR=99) LINE
            IF(INDEX(LINE,'END').NE.0) GO TO 56
            READ(LINE,*,END=56,ERR=99) XBT, YBT
          ENDIF
C
          IF(XBT.EQ.999.0) THEN
C--------- save number of points and type for element which was just read in
           NB(N) = IB-1
C--------- go try to read next element
           GO TO 55
          ENDIF
C
          IF(N.GT.NBX) THEN
           WRITE(*,*) 'AREAD: Too many airfoil elements. Increase NBX.'
           STOP
          ENDIF
C
C-------- uneventful XBT,YBT coordinates... just store them
          XB(IB,N) = XBT
          YB(IB,N) = YBT
C
        END DO
        WRITE(*,*) 'AREAD: Too many airfoil points.  Increase IBX.'
        STOP
C
   55 CONTINUE
      N = NBX
C
   56 IF(IB.EQ.1) THEN
C----- coordinate file has "999.0 999.0" at the end ...
       NBL = N-1
      ELSE
C----- coordinate file has no special ending line
       NBL = N
       NB(N) = IB-1
      ENDIF
C
C---- set type if coordinate file this was
      IF(IFTYPE.EQ.0) THEN
       IF(NBL.EQ.1) THEN
        IFTYPE = 3
       ELSE
        IFTYPE = 4
       ENDIF
      ENDIF
C
      RETURN
C
   99 CONTINUE
      WRITE(*,1100)
      IFTYPE = 0
      RETURN
C...............................................................
 1000 FORMAT(A)
 1100 FORMAT(/' File READ error:  ')
      END ! AREADNR



      SUBROUTINE AWRITE(FNAME,LU,
     &                  NEL,IFRST,ILAST,X,Y,
     &                  NAME, ISPARS, IFTYPE)
      CHARACTER*(*) FNAME
      CHARACTER*(*) NAME
      CHARACTER*(*) ISPARS
      DIMENSION IFRST(*), ILAST(*)
      DIMENSION X(*), Y(*)
C--------------------------------------------------------
C     Writes one of several types of airfoil coordinate files
C
C  Input:
C       FNAME     name of coordinate file to be read
C                   if blank, unit LU is assumed to be open
C       LU        logical unit for file write
C
C       NEL       number of elements
C       IFRST(.)
C       ILAST(.)
C       X(.)       x coordinates
C       Y(.)       r coordinates
C
C       NAME    character name string        (if IFTYPE > 1)
C       ISPARS  ISES/MSES domain-size string (if IFTYPE > 2)
C       IFTYPE  type of file to be written
C            1    Generic.
C            2    Labeled generic.
C            3    MSES single-element
C            4    MSES multi-element.
C--------------------------------------------------------
      CHARACTER*1 ANS
      LOGICAL LOPEN
C
      LOPEN = FNAME(1:1) .NE. ' '
      IF(LOPEN) THEN
       NF = INDEX(FNAME,' ') + 1
       OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=98)
       REWIND(LU)
      ENDIF
C
      IF(IFTYPE.NE.1) WRITE(LU,1000) NAME
      IF(IFTYPE.EQ.3 .OR.
     &   IFTYPE.EQ.4) WRITE(LU,1000) ISPARS
C
      DO IEL = 1, NEL
        DO I = IFRST(IEL), ILAST(IEL)
          WRITE(LU,1500) X(I), Y(I)
        ENDDO
C
        IF(IEL.LT.NEL) WRITE(LU,1550) 999.0, 999.0
      ENDDO
C
      IF(LOPEN) CLOSE(LU)
      RETURN
C
   98 CONTINUE
      WRITE(*,1050) FNAME(1:NF)
      RETURN
C...............................................................
 1000 FORMAT(A)
 1050 FORMAT(/' File OPEN error:  ', A)
C
 1500 FORMAT(1X,2F12.6)
 1550 FORMAT(1X,2F6.1)
      END ! AWRITE






