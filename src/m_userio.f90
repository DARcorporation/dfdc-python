module m_userio
    implicit none
contains
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
    subroutine aski1(prompt, iinput)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iinput
        character(*) :: prompt
        !
        ! Local variables
        !
        integer :: np
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- integer input
        !
        !
        np = index(prompt, '^') - 1
        if (np==0) np = len(prompt)
        !
        10   write (*, 1000) prompt(1:np)
        read (*, *, err = 10) iinput
        return
        !
        1000 format (/a, '   i>  ', $)
    end subroutine aski1
    !*==ASKR1.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ASKI

    subroutine askr1(prompt, rinput)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: prompt
        real :: rinput
        !
        ! Local variables
        !
        integer :: np
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- real input
        !
        !
        np = index(prompt, '^') - 1
        if (np==0) np = len(prompt)
        !
        10   write (*, 1000) prompt(1:np)
        read (*, *, err = 10) rinput
        return
        !
        1000 format (/a, '   r>  ', $)
    end subroutine askr1
    !*==ASKI.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ASKR


    subroutine aski(prompt, iinput)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iinput
        character(*) :: prompt
        !
        ! Local variables
        !
        logical :: error
        integer :: np
        integer, dimension(1) :: iinput_temp
        !
        !*** End of declarations rewritten by SPAG
        !
        !---- integer input
        !
        !
        np = index(prompt, '^') - 1
        if (np==0) np = len(prompt)
        do
            !
            write (*, 1010) prompt(1:np), iinput
            1010    format (1x, a, ':  ', i6)
            !
            call readi(1, iinput_temp, error)
            if (.not.(error)) then
                iinput = iinput_temp(1)
                exit
            end if
        enddo
        !
    end subroutine aski
    !*==ASKR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ASKI


    subroutine askr(prompt, rinput)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: prompt
        real :: rinput
        !
        ! Local variables
        !
        logical :: error
        integer :: np
        real, dimension(1) :: rinput_temp
        !
        !*** End of declarations rewritten by SPAG
        !
        !---- real input
        !
        !
        np = index(prompt, '^') - 1
        if (np==0) np = len(prompt)
        do
            !
            write (*, 1010) prompt(1:np), rinput
            1010    format (1x, a, ':  ', g12.6)
            !
            call readr(1, rinput_temp, error)
            if (.not.(error)) then
                rinput = rinput_temp(1)
                exit
            end if
        enddo
        !
    end subroutine askr
    !*==ASKIN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ASKR


    subroutine askin(prompt, iinput, ninp1)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ninp1
        character(*) :: prompt
        integer, dimension(*) :: iinput
        !
        ! Local variables
        !
        logical :: error
        integer :: i, ninp, np
        !
        !*** End of declarations rewritten by SPAG
        !
        !---- integer input
        !
        !
        np = index(prompt, '^') - 1
        if (np==0) np = len(prompt)
        !
        ninp = min(ninp1, 20)
        do
            !
            write (*, 1010) prompt(1:np), (iinput(i), i = 1, ninp)
            1010    format (1x, a, ':  ', 20i6)
            !
            call readi(ninp, iinput, error)
            if (.not.(error)) exit
        enddo
        !
    end subroutine askin
    !*==ASKRN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ASKIN


    subroutine askrn(prompt, rinput, ninp1)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ninp1
        character(*) :: prompt
        real, dimension(*) :: rinput
        !
        ! Local variables
        !
        logical :: error
        integer :: i, ninp, np
        !
        !*** End of declarations rewritten by SPAG
        !
        !---- real input
        !
        !
        np = index(prompt, '^') - 1
        if (np==0) np = len(prompt)
        !
        ninp = min(ninp1, 20)
        do
            !
            write (*, 1010) prompt(1:np), (rinput(i), i = 1, ninp)
            1010    format (1x, a, ':  ', 20g12.6)
            !
            call readr(ninp, rinput, error)
            if (.not.(error)) exit
        enddo
        !
    end subroutine askrn
    !*==ASKL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ASKRN



    subroutine askl(prompt, linput)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: linput
        character(*) :: prompt
        !
        ! Local variables
        !
        character(1) :: char
        integer :: np
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- logical input
        !
        !
        np = index(prompt, '^') - 1
        if (np==0) np = len(prompt)
        do
            !
            write (*, 1000) prompt(1:np)
            read (*, 1010) char
            if (char=='y') char = 'Y'
            if (char=='n') char = 'N'
            if (char=='Y' .or. char=='N') then
                !
                linput = char=='Y'
                return
                !
                1000       format (/a, ' y/n>  ', $)
                1010       format (a)
            endif
        enddo
    end subroutine askl
    !*==ASKS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ASKL


    subroutine asks(prompt, input)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: input, prompt
        !
        ! Local variables
        !
        integer :: np
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- string of arbitrary length input
        !
        !
        np = index(prompt, '^') - 1
        if (np==0) np = len(prompt)
        !
        write (*, 1000) prompt(1:np)
        read (*, 1010) input
        !
        return
        !
        1000 format (/a, '   s>  ', $)
        1010 format (a)
    end subroutine asks
    !*==ASKC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ASKS


    subroutine askc(prompt, comand, cargs)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: cargs, comand, prompt
        !
        ! Local variables
        !
        integer :: i, izero, k, ki, ncargs, np
        character(128) :: line
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- returns 4-byte character string input converted to uppercase
        !---- also returns rest of input characters in CARGS string
        !
        !
        !
        izero = ichar('0')
        !
        np = index(prompt, '^') - 1
        if (np==0) np = len(prompt)
        !
        write (*, 1000) prompt(1:np)
        read (*, 1020) line
        !
        !---- strip off leading blanks
        do k = 1, 128
            if (line(1:1)/=' ') exit
            line = line(2:128)
        enddo
        !
        !---- find position of first blank, "+", "-", ".", ",", or numeral
        k = index(line, ' ')
        ki = index(line, '-')
        if (ki/=0) k = min(k, ki)
        ki = index(line, '+')
        if (ki/=0) k = min(k, ki)
        ki = index(line, '.')
        if (ki/=0) k = min(k, ki)
        ki = index(line, ',')
        if (ki/=0) k = min(k, ki)
        do i = 0, 9
            ki = index(line, char(izero + i))
            if (ki/=0) k = min(k, ki)
        enddo
        !
        !---- there is no blank between command and argument... use first 4 characters
        if (k<=0) k = 5
        !
        if (k==1) then
            !------ the "command" is a number... set entire COMAND string with it
            comand = line
        else
            !------ the "command" is some string... just use the part up to the argument
            comand = line(1:k - 1)
        endif
        !
        !---- convert it to uppercase
        call lc2uc(comand)
        !
        cargs = line(k:128)
        call strip(cargs, ncargs)
        return
        !
        1000 format (/a, '   c>  ', $)
        1020 format (a)
    end subroutine askc
    !*==LC2UC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ASKC


    subroutine lc2uc(input)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: input
        !
        ! Local variables
        !
        integer :: i, k, n
        character(26), save :: lcase, ucase
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        data lcase/'abcdefghijklmnopqrstuvwxyz'/
        data ucase/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
        !
        n = len(input)
        !
        do i = 1, n
            k = index(lcase, input(i:i))
            if (k>0) input(i:i) = ucase(k:k)
        enddo
        !
    end subroutine lc2uc
    !*==READI.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! LC2UC



    subroutine readi(n, ivar, error)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: error
        integer :: n
        integer, dimension(n) :: ivar
        !
        ! Local variables
        !
        integer :: i, ntmp
        integer, dimension(40) :: ivtmp
        character(80) :: line
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Reads N integer variables, leaving unchanged
        !     if only <return> is entered.
        !--------------------------------------------------
        !
        read (*, 1000) line
        1000 format (a80)
        !
        do i = 1, n
            ivtmp(i) = ivar(i)
        enddo
        !
        ntmp = 40
        call getint(line, ivtmp, ntmp, error)
        !
        if (error) return
        !
        do i = 1, n
            ivar(i) = ivtmp(i)
        enddo
        !
    end subroutine readi
    !*==READR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! READI



    subroutine readr(n, var, error)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: error
        integer :: n
        real, dimension(n) :: var
        !
        ! Local variables
        !
        integer :: i, ntmp
        character(80) :: line
        real, dimension(40) :: vtmp
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------
        !     Reads N real variables, leaving unchanged
        !     if only <return> is entered.
        !-------------------------------------------------
        !
        read (*, 1000) line
        1000 format (a80)
        !
        do i = 1, n
            vtmp(i) = var(i)
        enddo
        !
        ntmp = 40
        call getflt(line, vtmp, ntmp, error)
        !
        if (error) return
        !
        do i = 1, n
            var(i) = vtmp(i)
        enddo
        !
    end subroutine readr
    !*==GETINT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! READR




    subroutine getint(input, a, n, error)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: error
        character(*) :: input
        integer :: n
        integer, dimension(*) :: a
        !
        ! Local variables
        !
        integer :: i, ilen, ilenp, ipass, k, kcomma, kspace, ninp
        character(130) :: rec
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
        ilen = min(len(input), 128)
        ilenp = ilen + 2
        !
        !---- put input into local work string (which will be munched)
        rec(1:ilenp) = input(1:ilen) // ' ,'
        !
        !---- ignore everything after a "!" character
        k = index(rec, '!')
        if (k>0) rec(1:ilen) = rec(1:k - 1)
        !
        ninp = n
        !
        !---- count up how many numbers are to be extracted
        n = 0
        k = 1
        do ipass = 1, ilen
            !------ search for next space or comma starting with current index K
            kspace = index(rec(k:ilenp), ' ') + k - 1
            kcomma = index(rec(k:ilenp), ',') + k - 1
            !
            if (k==kspace) then
                !------- just skip this space
                k = k + 1
                goto 9
            endif
            !
            if (k==kcomma) then
                !------- comma found.. increment number count and keep looking
                n = n + 1
                k = k + 1
                goto 9
            endif
            !
            !------ neither space nor comma found, so we ran into a number...
            !-    ...increment number counter and keep looking after next space or comma
            n = n + 1
            k = min(kspace, kcomma) + 1
            !
            9       if (k>=ilen) exit
        enddo
        !
        !---- decide on how many numbers to read, and go ahead and read them
        if (ninp>0) n = min(n, ninp)
        read (rec(1:ilen), *, err = 20) (a(i), i = 1, n)
        error = .false.
        return
        !
        !---- bzzzt !!!
        !cc   WRITE(*,*) 'GETINT: String-to-integer conversion error.'
        20   n = 0
        error = .true.
    end subroutine getint
    !*==GETFLT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine getflt(input, a, n, error)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: error
        character(*) :: input
        integer :: n
        real, dimension(*) :: a
        !
        ! Local variables
        !
        integer :: i, ilen, ilenp, ipass, k, kcomma, kspace, ninp
        character(130) :: rec
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
        ilen = min(len(input), 128)
        ilenp = ilen + 2
        !
        !---- put input into local work string (which will be munched)
        rec(1:ilenp) = input(1:ilen) // ' ,'
        !
        !---- ignore everything after a "!" character
        k = index(rec, '!')
        if (k>0) rec(1:ilen) = rec(1:k - 1)
        !
        ninp = n
        !
        !---- count up how many numbers are to be extracted
        n = 0
        k = 1
        do ipass = 1, ilen
            !------ search for next space or comma starting with current index K
            kspace = index(rec(k:ilenp), ' ') + k - 1
            kcomma = index(rec(k:ilenp), ',') + k - 1
            !
            if (k==kspace) then
                !------- just skip this space
                k = k + 1
                goto 9
            endif
            !
            if (k==kcomma) then
                !------- comma found.. increment number count and keep looking
                n = n + 1
                k = k + 1
                goto 9
            endif
            !
            !------ neither space nor comma found, so we ran into a number...
            !-    ...increment number counter and keep looking after next space or comma
            n = n + 1
            k = min(kspace, kcomma) + 1
            !
            9       if (k>=ilen) exit
        enddo
        !
        !---- decide on how many numbers to read, and go ahead and read them
        if (ninp>0) n = min(n, ninp)
        read (rec(1:ilen), *, err = 20) (a(i), i = 1, n)
        error = .false.
        return
        !
        !---- bzzzt !!!
        !cc   WRITE(*,*) 'GETFLT: String-to-integer conversion error.'
        20   n = 0
        error = .true.
    end subroutine getflt
    !*==STRIP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine strip(string, ns)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ns
        character(*) :: string
        !
        ! Local variables
        !
        integer :: k, k1, k2, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------
        !     Strips leading blanks off string
        !     and returns length of non-blank part.
        !-------------------------------------------
        n = len(string)
        !
        !---- find last non-blank character
        do k2 = n, 1, -1
            if (string(k2:k2)/=' ') goto 11
        enddo
        k2 = 0
        !
        !---- find first non-blank character
        11   do k1 = 1, k2
            if (string(k1:k1)/=' ') exit
        enddo
        !
        !---- number of non-blank characters
        ns = k2 - k1 + 1
        if (ns==0) return
        !
        !---- shift STRING so first character is non-blank
        string(1:ns) = string(k1:k2)
        !
        !---- pad tail of STRING with blanks
        do k = ns + 1, n
            string(k:k) = ' '
        enddo
        !
    end subroutine strip
    !*==GETARG0.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020

    subroutine getarg0(iarg, arg)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: arg
        integer :: iarg
        !
        ! Local variables
        !
        integer :: narg
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
        narg = iargc()
        if (narg>=iarg) then
            call getarg(iarg, arg)
        else
            arg = ' '
        endif
        !
    end subroutine getarg0
    !*==RDLINE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GETARG0


    subroutine rdline(lun, line, icnt)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: icnt, lun
        character(*) :: line
        !
        ! Local variables
        !
        logical, save :: lecho
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
        data lecho/.false./
        !
        1000 format (a)
        1010 format (i4, 1x, a)
        do
            !
            icnt = icnt + 1
            read (lun, 1000, end = 80, err = 90) line
            if (lecho) write (*, 1010) icnt, line(1:60)
            !
            !---- skip comment line
            if (index('!#', line(1:1))==0) then
                !
                !---- skip blank line
                !
                !---- normal return after significant line
                if (line/=' ') return
            endif
        enddo
        !
        80   line = 'END '
        return
        !
        90   line = 'ERR '
    end subroutine rdline
end module m_userio
