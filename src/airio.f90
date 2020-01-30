!=========================================================================
! DFDC (Ducted Fan Design Code) is an aerodynamic and aeroacoustic design
! and analysis tool for aircraft with propulsors in ducted fan
! configurations.
! 
! This software was developed under the auspices and sponsorship of the
! Tactical Technology Office (TTO) of the Defense Advanced Research
! Projects Agency (DARPA).
! 
! Copyright (c) 2004, 2005, Booz Allen Hamilton Inc., All Rights Reserved
!
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation; either version 2 of the License, or (at your
! option) any later version.
! 
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (drela@mit.edu)
! Program Management: Brad Tousley, Paul Eremenko (eremenko@alum.mit.edu)
!
!=========================================================================
!
!     Version 070-ES1
!     Philip Carter, Esotec Developments, February 2009
!     philip (at) esotec (dot) org
!
!     Changes from 0.70:
!
!     Trivial reporting format adjustment (line 275)
!
!=========================================================================
!
!
SUBROUTINE AREAD(FNAME, LU, IBX, NBX, XB, YB, &
        NB, NBL, &
        NAME, ISPARS, IFTYPE)
    CHARACTER*(*) FNAME
    CHARACTER*(*) NAME
    CHARACTER*(*) ISPARS
    DIMENSION XB(IBX, NBX), YB(IBX, NBX)
    DIMENSION NB(NBX)
    !--------------------------------------------------------
    !     Reads in several types of airfoil coordinate file.
    !
    !  Input:
    !       FNAME     name of coordinate file to be read
    !                   if blank, unit LU is assumed to be open
    !       LU        logical unit for file read
    !
    !  Output:
    !       XB,YB(..)  coordinates
    !       NB(.)      number of XB,YB coordinates
    !
    !       NAME    character name string        (if IFTYPE > 1)
    !       ISPARS  ISES/MSES domain-size string (if IFTYPE > 2)
    !       IFTYPE  returns type of file:
    !            0    None.  Read error occurred.
    !            1    Generic.
    !            2    Labeled generic.
    !            3    MSES single element.
    !            4    MSES multi-element.
    !--------------------------------------------------------
    CHARACTER*128 LINE
    LOGICAL LOPEN, ERROR
    !
    DIMENSION AINPUT(10)
    !
    !---- permitted characters in string containing Fortran-readable numbers
    CHARACTER*20 NCHARS
    DATA NCHARS / '0123456789-+.,EDed  ' /
    !
    !---- tab character is also permitted
    NCHARS(20:20) = CHAR(9)
    !
    !
    !---- first assume that there will be a read error
    IFTYPE = 0
    !
    LOPEN = FNAME(1:1) .NE. ' '
    IF(LOPEN) THEN
        NF = INDEX(FNAME, ' ') + 1
        OPEN(LU, FILE = FNAME, STATUS = 'OLD', ERR = 98)
        REWIND(LU)
    ENDIF
    !
    READ(LU, 1000) LINE
    !
    !---- if first line has any non-numeric character, go treat it as the name
    DO K = 1, 80
        IF(INDEX(NCHARS, LINE(K:K)) .EQ. 0) GO TO 20
    ENDDO
    !
    !---- plain unlabeled file: rewind, and just read in X,Y coordinates
    NAME = ' '
    !
    REWIND(LU)
    WRITE(*, *)
    WRITE(*, *) 'Reading plain coordinate file'
    IFTYPE = 1
    GO TO 40
    !
    !---- first line interpreted as label string
    20 CONTINUE
    NAME = LINE
    !
    !---- read second line
    READ(LU, 1000) LINE
    NINPUT = 10
    CALL GETFLT(LINE, AINPUT, NINPUT, ERROR)
    IF(ERROR) GO TO 99
    !
    IF(NINPUT.LT.4) THEN
        !------ no domain parameters: re-read name string and then read X,Y coordinates
        REWIND(LU)
        READ(LU, 1000) NAME
        IFTYPE = 2
    ENDIF
    !
    WRITE(*, 1010) NAME
    1010 FORMAT(/1X, 'Reading airfoil: ', A32)
    !
    40   CONTINUE
    !
    !---- read in airfoil coordinates
    DO 55 N = 1, NBX + 1
        DO IB = 1, IBX + 1
            READ(LU, *, END = 56, ERR = 99) XBT, YBT
            !
            IF(XBT.EQ.999.0) THEN
                !--------- save number of points and type for element which was just read in
                NB(N) = IB - 1
                !--------- go try to read next element
                GO TO 55
            ENDIF
            !
            IF(N.GT.NBX) THEN
                WRITE(*, *) 'AREAD: Too many airfoil elements. Increase NBX.'
                STOP
            ENDIF
            !
            !-------- uneventful XBT,YBT coordinates... just store them
            XB(IB, N) = XBT
            YB(IB, N) = YBT
            !
        END DO
        WRITE(*, *) 'AREAD: Too many airfoil points.  Increase IBX.'
        STOP
        !
    55 CONTINUE
    N = NBX
    !
    56 IF(IB.EQ.1) THEN
        !----- coordinate file has "999.0 999.0" at the end ...
        NBL = N - 1
    ELSE
        !----- coordinate file has no special ending line
        NBL = N
        NB(N) = IB - 1
    ENDIF
    !
    !---- set type if coordinate file this was
    IF(IFTYPE.EQ.0) THEN
        IF(NBL.EQ.1) THEN
            IFTYPE = 3
        ELSE
            IFTYPE = 4
        ENDIF
    ENDIF
    !
    IF(LOPEN) CLOSE(LU)
    RETURN
    !
    98 CONTINUE
    WRITE(*, 1050) FNAME(1:NF)
    IFTYPE = 0
    RETURN
    !
    99 CONTINUE
    WRITE(*, 1100) FNAME(1:NF)
    IFTYPE = 0
    IF(LOPEN) CLOSE(LU)
    RETURN
    !...............................................................
    1000 FORMAT(A)
    1050 FORMAT(/' File OPEN error:  ', A)
    1100 FORMAT(/' File READ error:  ', A)
END
! AREAD



SUBROUTINE AREADNR(LU, IBX, NBX, XB, YB, &
        NB, NBL, &
        NAME, ISPARS, IFTYPE)
    CHARACTER*80 NAME
    CHARACTER*(*) ISPARS
    DIMENSION XB(IBX, NBX), YB(IBX, NBX)
    DIMENSION NB(NBX)
    !--------------------------------------------------------
    !     Reads in several types of airfoil coordinate file
    !     without performing file open or rewinds.
    !     Assumes input airfoil file is open on entry
    !
    !  Input:
    !       LU        logical unit for file read
    !
    !  Output:
    !       XB,YB(..)  coordinates
    !       NB(.)      number of XB,YB coordinates
    !
    !       NAME    character name string        (if IFTYPE > 1)
    !       ISPARS  ISES/MSES domain-size string (if IFTYPE > 2)
    !       IFTYPE  returns type of file:
    !            0    None.  Read error occurred.
    !            1    Generic.
    !            2    Labeled generic.
    !            3    MSES single element.
    !            4    MSES multi-element.
    !--------------------------------------------------------
    CHARACTER*128 LINE
    LOGICAL ERROR
    !
    DIMENSION AINPUT(10)
    !
    !---- permitted characters in string containing Fortran-readable numbers
    CHARACTER*20 NCHARS
    DATA NCHARS / '0123456789-+.,EDed  ' /
    !
    !---- tab character is also permitted
    NCHARS(20:20) = CHAR(9)
    !
    !
    !---- first assume that there will be a read error
    IFTYPE = 0
    READ(LU, 1000) LINE
    !---- if first line has any non-numeric character, go treat it as the name
    DO K = 1, 80
        IF(INDEX(NCHARS, LINE(K:K)) .EQ. 0) GO TO 20
    ENDDO
    !
    !---- plain unlabeled file: rewind, and just read in X,Y coordinates
    NAME = ' '
    WRITE(*, *)
    WRITE(*, *) 'Reading plain coordinate file'
    IFTYPE = 1
    GO TO 40
    !
    !---- first line interpreted as label string
    20 CONTINUE
    NAME = LINE(1:60)
    !      CALL STRIP(NAME,NNM)
    !
    !---- read second line
    READ(LU, 1000) LINE
    NINPUT = 10
    CALL GETFLT(LINE, AINPUT, NINPUT, ERROR)
    IF(ERROR) GO TO 99
    !
    IF(NINPUT.LT.4) THEN
        !------ no domain parameters: read X,Y coordinates directly with this line
        IFTYPE = 2
    ELSE
        !---- second line had parameters, read third line (coordinates line)
        READ(LU, 1000) LINE
    ENDIF
    !
    WRITE(*, 1010) NAME
    1010 FORMAT(/1X, 'Reading airfoil: ', A32)
    !
    40   CONTINUE
    !
    !---- read in airfoil coordinates assuming first line has been read into
    !     character variable LINE
    DO 55 N = 1, NBX + 1
        DO IB = 1, IBX + 1
            IF(N.EQ.1 .AND. IB.EQ.1) THEN
                READ(LINE, *, END = 56, ERR = 99) XBT, YBT
            ELSE
                READ(LU, 1000, END = 56, ERR = 99) LINE
                IF(INDEX(LINE, 'END').NE.0) GO TO 56
                READ(LINE, *, END = 56, ERR = 99) XBT, YBT
            ENDIF
            !
            IF(XBT.EQ.999.0) THEN
                !--------- save number of points and type for element which was just read in
                NB(N) = IB - 1
                !--------- go try to read next element
                GO TO 55
            ENDIF
            !
            IF(N.GT.NBX) THEN
                WRITE(*, *) 'AREAD: Too many airfoil elements. Increase NBX.'
                STOP
            ENDIF
            !
            !-------- uneventful XBT,YBT coordinates... just store them
            XB(IB, N) = XBT
            YB(IB, N) = YBT
            !
        END DO
        WRITE(*, *) 'AREAD: Too many airfoil points.  Increase IBX.'
        STOP
        !
    55 CONTINUE
    N = NBX
    !
    56 IF(IB.EQ.1) THEN
        !----- coordinate file has "999.0 999.0" at the end ...
        NBL = N - 1
    ELSE
        !----- coordinate file has no special ending line
        NBL = N
        NB(N) = IB - 1
    ENDIF
    !
    !---- set type if coordinate file this was
    IF(IFTYPE.EQ.0) THEN
        IF(NBL.EQ.1) THEN
            IFTYPE = 3
        ELSE
            IFTYPE = 4
        ENDIF
    ENDIF
    !
    RETURN
    !
    99 CONTINUE
    WRITE(*, 1100)
    IFTYPE = 0
    RETURN
    !...............................................................
    1000 FORMAT(A)
    1100 FORMAT(/' File READ error:  ')
END
! AREADNR



SUBROUTINE AWRITE(FNAME, LU, &
        NEL, IFRST, ILAST, X, Y, &
        NAME, ISPARS, IFTYPE)
    CHARACTER*(*) FNAME
    CHARACTER*(*) NAME
    CHARACTER*(*) ISPARS
    DIMENSION IFRST(*), ILAST(*)
    DIMENSION X(*), Y(*)
    !--------------------------------------------------------
    !     Writes one of several types of airfoil coordinate files
    !
    !  Input:
    !       FNAME     name of coordinate file to be read
    !                   if blank, unit LU is assumed to be open
    !       LU        logical unit for file write
    !
    !       NEL       number of elements
    !       IFRST(.)
    !       ILAST(.)
    !       X(.)       x coordinates
    !       Y(.)       r coordinates
    !
    !       NAME    character name string        (if IFTYPE > 1)
    !       ISPARS  ISES/MSES domain-size string (if IFTYPE > 2)
    !       IFTYPE  type of file to be written
    !            1    Generic.
    !            2    Labeled generic.
    !            3    MSES single-element
    !            4    MSES multi-element.
    !--------------------------------------------------------
    CHARACTER*1 ANS
    LOGICAL LOPEN
    !
    LOPEN = FNAME(1:1) .NE. ' '
    IF(LOPEN) THEN
        NF = INDEX(FNAME, ' ') + 1
        OPEN(LU, FILE = FNAME, STATUS = 'OLD', ERR = 98)
        REWIND(LU)
    ENDIF
    !
    IF(IFTYPE.NE.1) WRITE(LU, 1000) NAME
    IF(IFTYPE.EQ.3 .OR.&
            IFTYPE.EQ.4) WRITE(LU, 1000) ISPARS
    !
    DO IEL = 1, NEL
        DO I = IFRST(IEL), ILAST(IEL)
            WRITE(LU, 1500) X(I), Y(I)
        ENDDO
        !
        IF(IEL.LT.NEL) WRITE(LU, 1550) 999.0, 999.0
    ENDDO
    !
    IF(LOPEN) CLOSE(LU)
    RETURN
    !
    98 CONTINUE
    WRITE(*, 1050) FNAME(1:NF)
    RETURN
    !...............................................................
    1000 FORMAT(A)
    1050 FORMAT(/' File OPEN error:  ', A)
    !
    1500 FORMAT(1X, 2F12.6)
    1550 FORMAT(1X, 2F6.1)
END
! AWRITE





