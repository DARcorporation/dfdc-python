!*==ASKI1.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
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
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER :: IINPUT
    CHARACTER(*) :: PROMPT
    !
    ! Local variables
    !
    INTEGER :: NP
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !---- integer input
    !
    !
    NP = INDEX(PROMPT, '^') - 1
    IF (NP==0) NP = LEN(PROMPT)
    !
    10   WRITE (*, 1000) PROMPT(1:NP)
    READ (*, *, ERR = 10) IINPUT
    RETURN
    !
    1000 FORMAT (/A, '   i>  ', $)
END SUBROUTINE ASKI1
!*==ASKR1.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! ASKI

SUBROUTINE ASKR1(PROMPT, RINPUT)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    CHARACTER(*) :: PROMPT
    REAL :: RINPUT
    !
    ! Local variables
    !
    INTEGER :: NP
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !---- real input
    !
    !
    NP = INDEX(PROMPT, '^') - 1
    IF (NP==0) NP = LEN(PROMPT)
    !
    10   WRITE (*, 1000) PROMPT(1:NP)
    READ (*, *, ERR = 10) RINPUT
    RETURN
    !
    1000 FORMAT (/A, '   r>  ', $)
END SUBROUTINE ASKR1
!*==ASKI.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! ASKR


SUBROUTINE ASKI(PROMPT, IINPUT)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER :: IINPUT
    CHARACTER(*) :: PROMPT
    !
    ! Local variables
    !
    LOGICAL :: ERROR
    INTEGER :: NP
    !
    !*** End of declarations rewritten by SPAG
    !
    !---- integer input
    !
    !
    NP = INDEX(PROMPT, '^') - 1
    IF (NP==0) NP = LEN(PROMPT)
    DO
        !
        WRITE (*, 1010) PROMPT(1:NP), IINPUT
        1010    FORMAT (1X, A, ':  ', I6)
        !
        CALL READI(1, IINPUT, ERROR)
        IF (.NOT.(ERROR)) EXIT
    ENDDO
    !
END SUBROUTINE ASKI
!*==ASKR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! ASKI


SUBROUTINE ASKR(PROMPT, RINPUT)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    CHARACTER(*) :: PROMPT
    REAL :: RINPUT
    !
    ! Local variables
    !
    LOGICAL :: ERROR
    INTEGER :: NP
    !
    !*** End of declarations rewritten by SPAG
    !
    !---- real input
    !
    !
    NP = INDEX(PROMPT, '^') - 1
    IF (NP==0) NP = LEN(PROMPT)
    DO
        !
        WRITE (*, 1010) PROMPT(1:NP), RINPUT
        1010    FORMAT (1X, A, ':  ', G12.6)
        !
        CALL READR(1, RINPUT, ERROR)
        IF (.NOT.(ERROR)) EXIT
    ENDDO
    !
END SUBROUTINE ASKR
!*==ASKIN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! ASKR


SUBROUTINE ASKIN(PROMPT, IINPUT, NINP1)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER :: NINP1
    CHARACTER(*) :: PROMPT
    INTEGER, DIMENSION(*) :: IINPUT
    !
    ! Local variables
    !
    LOGICAL :: ERROR
    INTEGER :: I, NINP, NP
    !
    !*** End of declarations rewritten by SPAG
    !
    !---- integer input
    !
    !
    NP = INDEX(PROMPT, '^') - 1
    IF (NP==0) NP = LEN(PROMPT)
    !
    NINP = MIN(NINP1, 20)
    DO
        !
        WRITE (*, 1010) PROMPT(1:NP), (IINPUT(I), I = 1, NINP)
        1010    FORMAT (1X, A, ':  ', 20I6)
        !
        CALL READI(NINP, IINPUT, ERROR)
        IF (.NOT.(ERROR)) EXIT
    ENDDO
    !
END SUBROUTINE ASKIN
!*==ASKRN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! ASKIN


SUBROUTINE ASKRN(PROMPT, RINPUT, NINP1)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER :: NINP1
    CHARACTER(*) :: PROMPT
    REAL, DIMENSION(*) :: RINPUT
    !
    ! Local variables
    !
    LOGICAL :: ERROR
    INTEGER :: I, NINP, NP
    !
    !*** End of declarations rewritten by SPAG
    !
    !---- real input
    !
    !
    NP = INDEX(PROMPT, '^') - 1
    IF (NP==0) NP = LEN(PROMPT)
    !
    NINP = MIN(NINP1, 20)
    DO
        !
        WRITE (*, 1010) PROMPT(1:NP), (RINPUT(I), I = 1, NINP)
        1010    FORMAT (1X, A, ':  ', 20G12.6)
        !
        CALL READR(NINP, RINPUT, ERROR)
        IF (.NOT.(ERROR)) EXIT
    ENDDO
    !
END SUBROUTINE ASKRN
!*==ASKL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! ASKRN



SUBROUTINE ASKL(PROMPT, LINPUT)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    LOGICAL :: LINPUT
    CHARACTER(*) :: PROMPT
    !
    ! Local variables
    !
    CHARACTER(1) :: CHAR
    INTEGER :: NP
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !---- logical input
    !
    !
    NP = INDEX(PROMPT, '^') - 1
    IF (NP==0) NP = LEN(PROMPT)
    DO
        !
        WRITE (*, 1000) PROMPT(1:NP)
        READ (*, 1010) CHAR
        IF (CHAR=='y') CHAR = 'Y'
        IF (CHAR=='n') CHAR = 'N'
        IF (CHAR=='Y' .OR. CHAR=='N') THEN
            !
            LINPUT = CHAR=='Y'
            RETURN
            !
            1000       FORMAT (/A, ' y/n>  ', $)
            1010       FORMAT (A)
        ENDIF
    ENDDO
END SUBROUTINE ASKL
!*==ASKS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! ASKL


SUBROUTINE ASKS(PROMPT, INPUT)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    CHARACTER(*) :: INPUT, PROMPT
    !
    ! Local variables
    !
    INTEGER :: NP
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !---- string of arbitrary length input
    !
    !
    NP = INDEX(PROMPT, '^') - 1
    IF (NP==0) NP = LEN(PROMPT)
    !
    WRITE (*, 1000) PROMPT(1:NP)
    READ (*, 1010) INPUT
    !
    RETURN
    !
    1000 FORMAT (/A, '   s>  ', $)
    1010 FORMAT (A)
END SUBROUTINE ASKS
!*==ASKC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! ASKS


SUBROUTINE ASKC(PROMPT, COMAND, CARGS)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    CHARACTER(*) :: CARGS, COMAND, PROMPT
    !
    ! Local variables
    !
    INTEGER :: I, IZERO, K, KI, NCARGS, NP
    CHARACTER(128) :: LINE
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !---- returns 4-byte character string input converted to uppercase
    !---- also returns rest of input characters in CARGS string
    !
    !
    !
    IZERO = ICHAR('0')
    !
    NP = INDEX(PROMPT, '^') - 1
    IF (NP==0) NP = LEN(PROMPT)
    !
    WRITE (*, 1000) PROMPT(1:NP)
    READ (*, 1020) LINE
    !
    !---- strip off leading blanks
    DO K = 1, 128
        IF (LINE(1:1)/=' ') EXIT
        LINE = LINE(2:128)
    ENDDO
    !
    !---- find position of first blank, "+", "-", ".", ",", or numeral
    K = INDEX(LINE, ' ')
    KI = INDEX(LINE, '-')
    IF (KI/=0) K = MIN(K, KI)
    KI = INDEX(LINE, '+')
    IF (KI/=0) K = MIN(K, KI)
    KI = INDEX(LINE, '.')
    IF (KI/=0) K = MIN(K, KI)
    KI = INDEX(LINE, ',')
    IF (KI/=0) K = MIN(K, KI)
    DO I = 0, 9
        KI = INDEX(LINE, CHAR(IZERO + I))
        IF (KI/=0) K = MIN(K, KI)
    ENDDO
    !
    !---- there is no blank between command and argument... use first 4 characters
    IF (K<=0) K = 5
    !
    IF (K==1) THEN
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
    1000 FORMAT (/A, '   c>  ', $)
    1020 FORMAT (A)
END SUBROUTINE ASKC
!*==LC2UC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! ASKC


SUBROUTINE LC2UC(INPUT)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    CHARACTER(*) :: INPUT
    !
    ! Local variables
    !
    INTEGER :: I, K, N
    CHARACTER(26), SAVE :: LCASE, UCASE
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    DATA LCASE/'abcdefghijklmnopqrstuvwxyz'/
    DATA UCASE/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
    !
    N = LEN(INPUT)
    !
    DO I = 1, N
        K = INDEX(LCASE, INPUT(I:I))
        IF (K>0) INPUT(I:I) = UCASE(K:K)
    ENDDO
    !
END SUBROUTINE LC2UC
!*==READI.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! LC2UC



SUBROUTINE READI(N, IVAR, ERROR)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    LOGICAL :: ERROR
    INTEGER :: N
    INTEGER, DIMENSION(N) :: IVAR
    !
    ! Local variables
    !
    INTEGER :: I, NTMP
    INTEGER, DIMENSION(40) :: IVTMP
    CHARACTER(80) :: LINE
    !
    !*** End of declarations rewritten by SPAG
    !
    !--------------------------------------------------
    !     Reads N integer variables, leaving unchanged
    !     if only <return> is entered.
    !--------------------------------------------------
    !
    READ (*, 1000) LINE
    1000 FORMAT (A80)
    !
    DO I = 1, N
        IVTMP(I) = IVAR(I)
    ENDDO
    !
    NTMP = 40
    CALL GETINT(LINE, IVTMP, NTMP, ERROR)
    !
    IF (ERROR) RETURN
    !
    DO I = 1, N
        IVAR(I) = IVTMP(I)
    ENDDO
    !
END SUBROUTINE READI
!*==READR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! READI



SUBROUTINE READR(N, VAR, ERROR)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    LOGICAL :: ERROR
    INTEGER :: N
    REAL, DIMENSION(N) :: VAR
    !
    ! Local variables
    !
    INTEGER :: I, NTMP
    CHARACTER(80) :: LINE
    REAL, DIMENSION(40) :: VTMP
    !
    !*** End of declarations rewritten by SPAG
    !
    !-------------------------------------------------
    !     Reads N real variables, leaving unchanged
    !     if only <return> is entered.
    !-------------------------------------------------
    !
    READ (*, 1000) LINE
    1000 FORMAT (A80)
    !
    DO I = 1, N
        VTMP(I) = VAR(I)
    ENDDO
    !
    NTMP = 40
    CALL GETFLT(LINE, VTMP, NTMP, ERROR)
    !
    IF (ERROR) RETURN
    !
    DO I = 1, N
        VAR(I) = VTMP(I)
    ENDDO
    !
END SUBROUTINE READR
!*==GETINT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! READR




SUBROUTINE GETINT(INPUT, A, N, ERROR)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    LOGICAL :: ERROR
    CHARACTER(*) :: INPUT
    INTEGER :: N
    INTEGER, DIMENSION(*) :: A
    !
    ! Local variables
    !
    INTEGER :: I, ILEN, ILENP, IPASS, K, KCOMMA, KSPACE, NINP
    CHARACTER(130) :: REC
    !
    !*** End of declarations rewritten by SPAG
    !
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
    IF (K>0) REC(1:ILEN) = REC(1:K - 1)
    !
    NINP = N
    !
    !---- count up how many numbers are to be extracted
    N = 0
    K = 1
    DO IPASS = 1, ILEN
        !------ search for next space or comma starting with current index K
        KSPACE = INDEX(REC(K:ILENP), ' ') + K - 1
        KCOMMA = INDEX(REC(K:ILENP), ',') + K - 1
        !
        IF (K==KSPACE) THEN
            !------- just skip this space
            K = K + 1
            GOTO 9
        ENDIF
        !
        IF (K==KCOMMA) THEN
            !------- comma found.. increment number count and keep looking
            N = N + 1
            K = K + 1
            GOTO 9
        ENDIF
        !
        !------ neither space nor comma found, so we ran into a number...
        !-    ...increment number counter and keep looking after next space or comma
        N = N + 1
        K = MIN(KSPACE, KCOMMA) + 1
        !
        9       IF (K>=ILEN) EXIT
    ENDDO
    !
    !---- decide on how many numbers to read, and go ahead and read them
    IF (NINP>0) N = MIN(N, NINP)
    READ (REC(1:ILEN), *, ERR = 20) (A(I), I = 1, N)
    ERROR = .FALSE.
    RETURN
    !
    !---- bzzzt !!!
    !cc   WRITE(*,*) 'GETINT: String-to-integer conversion error.'
    20   N = 0
    ERROR = .TRUE.
END SUBROUTINE GETINT
!*==GETFLT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


SUBROUTINE GETFLT(INPUT, A, N, ERROR)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    LOGICAL :: ERROR
    CHARACTER(*) :: INPUT
    INTEGER :: N
    REAL, DIMENSION(*) :: A
    !
    ! Local variables
    !
    INTEGER :: I, ILEN, ILENP, IPASS, K, KCOMMA, KSPACE, NINP
    CHARACTER(130) :: REC
    !
    !*** End of declarations rewritten by SPAG
    !
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
    IF (K>0) REC(1:ILEN) = REC(1:K - 1)
    !
    NINP = N
    !
    !---- count up how many numbers are to be extracted
    N = 0
    K = 1
    DO IPASS = 1, ILEN
        !------ search for next space or comma starting with current index K
        KSPACE = INDEX(REC(K:ILENP), ' ') + K - 1
        KCOMMA = INDEX(REC(K:ILENP), ',') + K - 1
        !
        IF (K==KSPACE) THEN
            !------- just skip this space
            K = K + 1
            GOTO 9
        ENDIF
        !
        IF (K==KCOMMA) THEN
            !------- comma found.. increment number count and keep looking
            N = N + 1
            K = K + 1
            GOTO 9
        ENDIF
        !
        !------ neither space nor comma found, so we ran into a number...
        !-    ...increment number counter and keep looking after next space or comma
        N = N + 1
        K = MIN(KSPACE, KCOMMA) + 1
        !
        9       IF (K>=ILEN) EXIT
    ENDDO
    !
    !---- decide on how many numbers to read, and go ahead and read them
    IF (NINP>0) N = MIN(N, NINP)
    READ (REC(1:ILEN), *, ERR = 20) (A(I), I = 1, N)
    ERROR = .FALSE.
    RETURN
    !
    !---- bzzzt !!!
    !cc   WRITE(*,*) 'GETFLT: String-to-integer conversion error.'
    20   N = 0
    ERROR = .TRUE.
END SUBROUTINE GETFLT
!*==STRIP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


SUBROUTINE STRIP(STRING, NS)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER :: NS
    CHARACTER(*) :: STRING
    !
    ! Local variables
    !
    INTEGER :: K, K1, K2, N
    !
    !*** End of declarations rewritten by SPAG
    !
    !-------------------------------------------
    !     Strips leading blanks off string
    !     and returns length of non-blank part.
    !-------------------------------------------
    N = LEN(STRING)
    !
    !---- find last non-blank character
    DO K2 = N, 1, -1
        IF (STRING(K2:K2)/=' ') GOTO 11
    ENDDO
    K2 = 0
    !
    !---- find first non-blank character
    11   DO K1 = 1, K2
        IF (STRING(K1:K1)/=' ') EXIT
    ENDDO
    !
    !---- number of non-blank characters
    NS = K2 - K1 + 1
    IF (NS==0) RETURN
    !
    !---- shift STRING so first character is non-blank
    STRING(1:NS) = STRING(K1:K2)
    !
    !---- pad tail of STRING with blanks
    DO K = NS + 1, N
        STRING(K:K) = ' '
    ENDDO
    !
END SUBROUTINE STRIP
!*==GETARG0.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020

SUBROUTINE GETARG0(IARG, ARG)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    CHARACTER(*) :: ARG
    INTEGER :: IARG
    !
    ! Local variables
    !
    INTEGER :: NARG
    !
    !*** End of declarations rewritten by SPAG
    !
    !------------------------------------------------
    !     Same as GETARG, but...
    !
    !     ...in the case of Intel Fortran, this one
    !     doesn't barf if there's no Unix argument
    !      (just returns blank string instead)
    !------------------------------------------------
    !
    NARG = IARGC()
    IF (NARG>=IARG) THEN
        CALL GETARG(IARG, ARG)
    ELSE
        ARG = ' '
    ENDIF
    !
END SUBROUTINE GETARG0
!*==RDLINE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! GETARG0


SUBROUTINE RDLINE(LUN, LINE, ICNT)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER :: ICNT, LUN
    CHARACTER(*) :: LINE
    !
    ! Local variables
    !
    LOGICAL, SAVE :: LECHO
    !
    !*** End of declarations rewritten by SPAG
    !
    !...Purpose  Read a non-comment line from the input file
    !
    !...Input    Data read from unit LUN
    !...Output   LINE  Character string with input line
    !                  LINE is set to 'END' for end or errors
    !
    !   Comment lines are assumed to start with ! or #
    !
    DATA LECHO/.FALSE./
    !
    1000 FORMAT (A)
    1010 FORMAT (I4, 1X, A)
    DO
        !
        ICNT = ICNT + 1
        READ (LUN, 1000, END = 80, ERR = 90) LINE
        IF (LECHO) WRITE (*, 1010) ICNT, LINE(1:60)
        !
        !---- skip comment line
        IF (INDEX('!#', LINE(1:1))==0) THEN
            !
            !---- skip blank line
            !
            !---- normal return after significant line
            IF (LINE/=' ') RETURN
        ENDIF
    ENDDO
    !
    80   LINE = 'END '
    RETURN
    !
    90   LINE = 'ERR '
END SUBROUTINE RDLINE
