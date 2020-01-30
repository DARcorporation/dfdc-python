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
!
!==== user input routines with prompting and error trapping
!
SUBROUTINE ASKI1(PROMPT, IINPUT)
    !
    !---- integer input
    !
    CHARACTER*(*) PROMPT
    INTEGER IINPUT
    !
    NP = INDEX(PROMPT, '^') - 1
    IF(NP.EQ.0) NP = LEN(PROMPT)
    !
    10   WRITE(*, 1000) PROMPT(1:NP)
    READ (*, *, ERR = 10) IINPUT
    RETURN
    !
    1000 FORMAT(/A, '   i>  ', $)
END
! ASKI

SUBROUTINE ASKR1(PROMPT, RINPUT)
    !
    !---- real input
    !
    CHARACTER*(*) PROMPT
    REAL RINPUT
    !
    NP = INDEX(PROMPT, '^') - 1
    IF(NP.EQ.0) NP = LEN(PROMPT)
    !
    10   WRITE(*, 1000) PROMPT(1:NP)
    READ (*, *, ERR = 10) RINPUT
    RETURN
    !
    1000 FORMAT(/A, '   r>  ', $)
END
! ASKR


SUBROUTINE ASKI(PROMPT, IINPUT)
    !---- integer input
    CHARACTER*(*) PROMPT
    INTEGER IINPUT
    !
    LOGICAL ERROR
    !
    NP = INDEX(PROMPT, '^') - 1
    IF(NP.EQ.0) NP = LEN(PROMPT)
    !
    10   WRITE(*, 1010) PROMPT(1:NP), IINPUT
    1010 FORMAT(1X, A, ':  ', I6)
    !
    CALL READI(1, IINPUT, ERROR)
    IF(ERROR) GO TO 10
    !
    RETURN
END
! ASKI


SUBROUTINE ASKR(PROMPT, RINPUT)
    !---- real input
    CHARACTER*(*) PROMPT
    REAL RINPUT
    !
    LOGICAL ERROR
    !
    NP = INDEX(PROMPT, '^') - 1
    IF(NP.EQ.0) NP = LEN(PROMPT)
    !
    10   WRITE(*, 1010) PROMPT(1:NP), RINPUT
    1010 FORMAT(1X, A, ':  ', G12.6)
    !
    CALL READR(1, RINPUT, ERROR)
    IF(ERROR) GO TO 10
    !
    RETURN
END
! ASKR


SUBROUTINE ASKIN(PROMPT, IINPUT, NINP1)
    !---- integer input
    CHARACTER*(*) PROMPT
    INTEGER IINPUT(*)
    !
    LOGICAL ERROR
    !
    NP = INDEX(PROMPT, '^') - 1
    IF(NP.EQ.0) NP = LEN(PROMPT)
    !
    NINP = MIN(NINP1, 20)
    !
    10   WRITE(*, 1010) PROMPT(1:NP), (IINPUT(I), I = 1, NINP)
    1010 FORMAT(1X, A, ':  ', 20I6)
    !
    CALL READI(NINP, IINPUT, ERROR)
    IF(ERROR) GO TO 10
    !
    RETURN
END
! ASKIN


SUBROUTINE ASKRN(PROMPT, RINPUT, NINP1)
    !---- real input
    CHARACTER*(*) PROMPT
    REAL RINPUT(*)
    !
    LOGICAL ERROR
    !
    NP = INDEX(PROMPT, '^') - 1
    IF(NP.EQ.0) NP = LEN(PROMPT)
    !
    NINP = MIN(NINP1, 20)
    !
    10   WRITE(*, 1010) PROMPT(1:NP), (RINPUT(I), I = 1, NINP)
    1010 FORMAT(1X, A, ':  ', 20G12.6)
    !
    CALL READR(NINP, RINPUT, ERROR)
    IF(ERROR) GO TO 10
    !
    RETURN
END
! ASKRN



SUBROUTINE ASKL(PROMPT, LINPUT)
    !
    !---- logical input
    !
    CHARACTER*(*) PROMPT
    LOGICAL LINPUT
    CHARACTER*1 CHAR
    !
    NP = INDEX(PROMPT, '^') - 1
    IF(NP.EQ.0) NP = LEN(PROMPT)
    !
    10   WRITE(*, 1000) PROMPT(1:NP)
    READ (*, 1010) CHAR
    IF(CHAR.EQ.'y') CHAR = 'Y'
    IF(CHAR.EQ.'n') CHAR = 'N'
    IF(CHAR.NE.'Y' .AND. CHAR.NE.'N') GO TO 10
    !
    LINPUT = CHAR .EQ. 'Y'
    RETURN
    !
    1000 FORMAT(/A, ' y/n>  ', $)
    1010 FORMAT(A)
END
! ASKL


SUBROUTINE ASKS(PROMPT, INPUT)
    !
    !---- string of arbitrary length input
    !
    CHARACTER*(*) PROMPT
    CHARACTER*(*) INPUT
    !
    NP = INDEX(PROMPT, '^') - 1
    IF(NP.EQ.0) NP = LEN(PROMPT)
    !
    WRITE(*, 1000) PROMPT(1:NP)
    READ (*, 1010) INPUT
    !
    RETURN
    !
    1000 FORMAT(/A, '   s>  ', $)
    1010 FORMAT(A)
END
! ASKS


SUBROUTINE ASKC(PROMPT, COMAND, CARGS)
    !
    !---- returns 4-byte character string input converted to uppercase
    !---- also returns rest of input characters in CARGS string
    !
    CHARACTER*(*) PROMPT
    CHARACTER*(*) COMAND, CARGS
    !
    CHARACTER*128 LINE
    LOGICAL ERROR
    !
    IZERO = ICHAR('0')
    !
    NP = INDEX(PROMPT, '^') - 1
    IF(NP.EQ.0) NP = LEN(PROMPT)
    !
    WRITE(*, 1000) PROMPT(1:NP)
    READ (*, 1020) LINE
    !
    !---- strip off leading blanks
    DO K = 1, 128
        IF(LINE(1:1) .EQ. ' ') THEN
            LINE = LINE(2:128)
        ELSE
            GO TO 5
        ENDIF
    ENDDO
    5    CONTINUE
    !
    !---- find position of first blank, "+", "-", ".", ",", or numeral
    K = INDEX(LINE, ' ')
    KI = INDEX(LINE, '-')
    IF(KI.NE.0) K = MIN(K, KI)
    KI = INDEX(LINE, '+')
    IF(KI.NE.0) K = MIN(K, KI)
    KI = INDEX(LINE, '.')
    IF(KI.NE.0) K = MIN(K, KI)
    KI = INDEX(LINE, ',')
    IF(KI.NE.0) K = MIN(K, KI)
    DO I = 0, 9
        KI = INDEX(LINE, CHAR(IZERO + I))
        IF(KI.NE.0) K = MIN(K, KI)
    ENDDO
    !
    !---- there is no blank between command and argument... use first 4 characters
    IF(K.LE.0) K = 5
    !
    IF(K.EQ.1) THEN
        !------ the "command" is a number... set entire COMAND string with it
        COMAND = LINE
    ELSE
        !------ the "command" is some string... just use the part up to the argument
        COMAND = LINE(1:K - 1)
    ENDIF
    !
    !---- convert it to uppercase
    CALL LC2UC(COMAND)
    !
    CARGS = LINE(K:128)
    CALL STRIP(CARGS, NCARGS)
    RETURN
    !
    1000 FORMAT(/A, '   c>  ', $)
    1020 FORMAT(A)
END
! ASKC


SUBROUTINE LC2UC(INPUT)
    CHARACTER*(*) INPUT
    !
    CHARACTER*26 LCASE, UCASE
    DATA LCASE / 'abcdefghijklmnopqrstuvwxyz' /
    DATA UCASE / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
    !
    N = LEN(INPUT)
    !
    DO 10 I = 1, N
        K = INDEX(LCASE, INPUT(I:I))
        IF(K.GT.0) INPUT(I:I) = UCASE(K:K)
    10   CONTINUE
    !
    RETURN
END
! LC2UC



SUBROUTINE READI(N, IVAR, ERROR)
    DIMENSION IVAR(N)
    LOGICAL ERROR
    !--------------------------------------------------
    !     Reads N integer variables, leaving unchanged
    !     if only <return> is entered.
    !--------------------------------------------------
    DIMENSION IVTMP(40)
    CHARACTER*80 LINE
    !
    READ(*, 1000) LINE
    1000 FORMAT(A80)
    !
    DO 10 I = 1, N
        IVTMP(I) = IVAR(I)
    10   CONTINUE
    !
    NTMP = 40
    CALL GETINT(LINE, IVTMP, NTMP, ERROR)
    !
    IF(ERROR) RETURN
    !
    DO 20 I = 1, N
        IVAR(I) = IVTMP(I)
    20   CONTINUE
    !
    RETURN
END
! READI



SUBROUTINE READR(N, VAR, ERROR)
    DIMENSION VAR(N)
    LOGICAL ERROR
    !-------------------------------------------------
    !     Reads N real variables, leaving unchanged
    !     if only <return> is entered.
    !-------------------------------------------------
    DIMENSION VTMP(40)
    CHARACTER*80 LINE
    !
    READ(*, 1000) LINE
    1000 FORMAT(A80)
    !
    DO 10 I = 1, N
        VTMP(I) = VAR(I)
    10   CONTINUE
    !
    NTMP = 40
    CALL GETFLT(LINE, VTMP, NTMP, ERROR)
    !
    IF(ERROR) RETURN
    !
    DO 20 I = 1, N
        VAR(I) = VTMP(I)
    20   CONTINUE
    !
    RETURN
END
! READR




SUBROUTINE GETINT(INPUT, A, N, ERROR)
    CHARACTER*(*) INPUT
    INTEGER A(*)
    LOGICAL ERROR
    !----------------------------------------------------------
    !     Parses character string INPUT into an array
    !     of integer numbers returned in A(1...N)
    !
    !     Will attempt to extract no more than N numbers,
    !     unless N = 0, in which case all numbers present
    !     in INPUT will be extracted.
    !
    !     N returns how many numbers were actually extracted.
    !----------------------------------------------------------
    CHARACTER*130 REC
    !
    !---- only first 128 characters in INPUT will be parsed
    ILEN = MIN(LEN(INPUT), 128)
    ILENP = ILEN + 2
    !
    !---- put input into local work string (which will be munched)
    REC(1:ILENP) = INPUT(1:ILEN) // ' ,'
    !
    !---- ignore everything after a "!" character
    K = INDEX(REC, '!')
    IF(K.GT.0) REC(1:ILEN) = REC(1:K - 1)
    !
    NINP = N
    !
    !---- count up how many numbers are to be extracted
    N = 0
    K = 1
    DO 10 IPASS = 1, ILEN
        !------ search for next space or comma starting with current index K
        KSPACE = INDEX(REC(K:ILENP), ' ') + K - 1
        KCOMMA = INDEX(REC(K:ILENP), ',') + K - 1
        !
        IF(K.EQ.KSPACE) THEN
            !------- just skip this space
            K = K + 1
            GO TO 9
        ENDIF
        !
        IF(K.EQ.KCOMMA) THEN
            !------- comma found.. increment number count and keep looking
            N = N + 1
            K = K + 1
            GO TO 9
        ENDIF
        !
        !------ neither space nor comma found, so we ran into a number...
        !-    ...increment number counter and keep looking after next space or comma
        N = N + 1
        K = MIN(KSPACE, KCOMMA) + 1
        !
        9     IF(K.GE.ILEN) GO TO 11
    10   CONTINUE
    !
    !---- decide on how many numbers to read, and go ahead and read them
    11   IF(NINP.GT.0) N = MIN(N, NINP)
    READ(REC(1:ILEN), *, ERR = 20) (A(I), I = 1, N)
    ERROR = .FALSE.
    RETURN
    !
    !---- bzzzt !!!
    20   CONTINUE
    !cc   WRITE(*,*) 'GETINT: String-to-integer conversion error.'
    N = 0
    ERROR = .TRUE.
    RETURN
END


SUBROUTINE GETFLT(INPUT, A, N, ERROR)
    CHARACTER*(*) INPUT
    REAL A(*)
    LOGICAL ERROR
    !----------------------------------------------------------
    !     Parses character string INPUT into an array
    !     of real numbers returned in A(1...N)
    !
    !     Will attempt to extract no more than N numbers,
    !     unless N = 0, in which case all numbers present
    !     in INPUT will be extracted.
    !
    !     N returns how many numbers were actually extracted.
    !----------------------------------------------------------
    CHARACTER*130 REC
    !
    !---- only first 128 characters in INPUT will be parsed
    ILEN = MIN(LEN(INPUT), 128)
    ILENP = ILEN + 2
    !
    !---- put input into local work string (which will be munched)
    REC(1:ILENP) = INPUT(1:ILEN) // ' ,'
    !
    !---- ignore everything after a "!" character
    K = INDEX(REC, '!')
    IF(K.GT.0) REC(1:ILEN) = REC(1:K - 1)
    !
    NINP = N
    !
    !---- count up how many numbers are to be extracted
    N = 0
    K = 1
    DO 10 IPASS = 1, ILEN
        !------ search for next space or comma starting with current index K
        KSPACE = INDEX(REC(K:ILENP), ' ') + K - 1
        KCOMMA = INDEX(REC(K:ILENP), ',') + K - 1
        !
        IF(K.EQ.KSPACE) THEN
            !------- just skip this space
            K = K + 1
            GO TO 9
        ENDIF
        !
        IF(K.EQ.KCOMMA) THEN
            !------- comma found.. increment number count and keep looking
            N = N + 1
            K = K + 1
            GO TO 9
        ENDIF
        !
        !------ neither space nor comma found, so we ran into a number...
        !-    ...increment number counter and keep looking after next space or comma
        N = N + 1
        K = MIN(KSPACE, KCOMMA) + 1
        !
        9     IF(K.GE.ILEN) GO TO 11
    10   CONTINUE
    !
    !---- decide on how many numbers to read, and go ahead and read them
    11   IF(NINP.GT.0) N = MIN(N, NINP)
    READ(REC(1:ILEN), *, ERR = 20) (A(I), I = 1, N)
    ERROR = .FALSE.
    RETURN
    !
    !---- bzzzt !!!
    20   CONTINUE
    !cc   WRITE(*,*) 'GETFLT: String-to-integer conversion error.'
    N = 0
    ERROR = .TRUE.
    RETURN
END


SUBROUTINE STRIP(STRING, NS)
    CHARACTER*(*) STRING
    !-------------------------------------------
    !     Strips leading blanks off string
    !     and returns length of non-blank part.
    !-------------------------------------------
    N = LEN(STRING)
    !
    !---- find last non-blank character
    DO 10 K2 = N, 1, -1
        IF(STRING(K2:K2).NE.' ') GO TO 11
    10 CONTINUE
    K2 = 0
    11 CONTINUE
    !
    !---- find first non-blank character
    DO 20 K1 = 1, K2
        IF(STRING(K1:K1).NE.' ') GO TO 21
    20 CONTINUE
    21 CONTINUE
    !
    !---- number of non-blank characters
    NS = K2 - K1 + 1
    IF(NS.EQ.0) RETURN
    !
    !---- shift STRING so first character is non-blank
    STRING(1:NS) = STRING(K1:K2)
    !
    !---- pad tail of STRING with blanks
    DO 30 K = NS + 1, N
        STRING(K:K) = ' '
    30 CONTINUE
    !
    RETURN
END

SUBROUTINE GETARG0(IARG, ARG)
    !------------------------------------------------
    !     Same as GETARG, but...
    !
    !     ...in the case of Intel Fortran, this one
    !     doesn't barf if there's no Unix argument
    !      (just returns blank string instead)
    !------------------------------------------------
    CHARACTER*(*) ARG
    !
    NARG = IARGC()
    IF(NARG.GE.IARG) THEN
        CALL GETARG(IARG, ARG)
    ELSE
        ARG = ' '
    ENDIF
    !
    RETURN
END
! GETARG0


SUBROUTINE RDLINE(LUN, LINE, ICNT)
    !...Purpose  Read a non-comment line from the input file
    !
    !...Input    Data read from unit LUN
    !...Output   LINE  Character string with input line
    !                  LINE is set to 'END' for end or errors
    !
    !   Comment lines are assumed to start with ! or #
    !
    CHARACTER*(*) LINE
    LOGICAL LECHO
    DATA LECHO / .FALSE. /
    !
    1000 FORMAT(A)
    1010 FORMAT(I4, 1X, A)
    !
    20 ICNT = ICNT + 1
    READ (LUN, 1000, END = 80, ERR = 90) LINE
    IF(LECHO) WRITE(*, 1010) ICNT, LINE(1:60)
    !
    !---- skip comment line
    IF(INDEX('!#', LINE(1:1)) .NE. 0) GO TO 20
    !
    !---- skip blank line
    IF(LINE.EQ.' ') GO TO 20
    !
    !---- normal return after significant line
    RETURN
    !
    80 LINE = 'END '
    RETURN
    !
    90 LINE = 'ERR '
    RETURN
END




