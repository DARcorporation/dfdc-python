module m_airio
    implicit none
contains
    !*==AREAD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PNTAICD
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
    subroutine aread(fname, lu, ibx, nbx, xb, yb, nb, nbl, name, ispars, iftype)
        use m_userio, only : getflt
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: fname, ispars, name
        integer :: ibx, iftype, lu, nbl, nbx
        integer, dimension(nbx) :: nb
        real, dimension(ibx, nbx) :: xb, yb
        !
        ! Local variables
        !
        real, dimension(10) :: ainput
        logical :: error, lopen
        integer :: ib, k, n, nf, ninput
        character(128) :: line
        character(20), save :: nchars
        real :: xbt, ybt
        !
        !*** End of declarations rewritten by SPAG
        !
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
        !
        !
        !---- permitted characters in string containing Fortran-readable numbers
        data nchars/'0123456789-+.,EDed  '/
        !
        !---- tab character is also permitted
        nchars(20:20) = char(9)
        !
        !
        !---- first assume that there will be a read error
        iftype = 0
        !
        lopen = fname(1:1)/=' '
        if (lopen) then
            nf = index(fname, ' ') + 1
            open (lu, file = fname, status = 'OLD', err = 98)
            rewind (lu)
        endif
        !
        read (lu, 1000) line
        !
        !---- if first line has any non-numeric character, go treat it as the name
        do k = 1, 80
            if (index(nchars, line(k:k))==0) goto 20
        enddo
        !
        !---- plain unlabeled file: rewind, and just read in X,Y coordinates
        name = ' '
        !
        rewind (lu)
        write (*, *)
        write (*, *) 'Reading plain coordinate file'
        iftype = 1
        goto 40
        !
        !---- first line interpreted as label string
        20   name = line
        !
        !---- read second line
        read (lu, 1000) line
        ninput = 10
        call getflt(line, ainput, ninput, error)
        if (error) goto 99
        !
        if (ninput<4) then
            !------ no domain parameters: re-read name string and then read X,Y coordinates
            rewind (lu)
            read (lu, 1000) name
            iftype = 2
        endif
        !
        write (*, 1010) name
        1010 format (/1x, 'Reading airfoil: ', a32)
        !
        !
        !---- read in airfoil coordinates
        40   do n = 1, nbx + 1
            do ib = 1, ibx + 1
                read (lu, *, end = 56, err = 99) xbt, ybt
                !
                if (xbt==999.0) then
                    !--------- save number of points and type for element which was just read in
                    nb(n) = ib - 1
                    !--------- go try to read next element
                    goto 55
                endif
                !
                if (n>nbx) then
                    write (*, *)                                              &
                            &'AREAD: Too many airfoil elements. Increase NBX.'
                    stop
                endif
                !
                !-------- uneventful XBT,YBT coordinates... just store them
                xb(ib, n) = xbt
                yb(ib, n) = ybt
                !
            enddo
            write (*, *) 'AREAD: Too many airfoil points.  Increase IBX.'
            stop
            !
        55   enddo
        n = nbx
        !
        56   if (ib==1) then
            !----- coordinate file has "999.0 999.0" at the end ...
            nbl = n - 1
        else
            !----- coordinate file has no special ending line
            nbl = n
            nb(n) = ib - 1
        endif
        !
        !---- set type if coordinate file this was
        if (iftype==0) then
            if (nbl==1) then
                iftype = 3
            else
                iftype = 4
            endif
        endif
        !
        if (lopen) close (lu)
        return
        !
        98   write (*, 1050) fname(1:nf)
        iftype = 0
        return
        !
        99   write (*, 1100) fname(1:nf)
        iftype = 0
        if (lopen) close (lu)
        return
        !...............................................................
        1000 format (a)
        1050 format (/' File OPEN error:  ', a)
        1100 format (/' File READ error:  ', a)
    end subroutine aread
    !*==AREADNR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! AREAD



    subroutine areadnr(lu, ibx, nbx, xb, yb, nb, nbl, name, ispars, iftype)
        use m_userio, only : getflt
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ibx, iftype, lu, nbl, nbx
        character(*) :: ispars
        character(80) :: name
        integer, dimension(nbx) :: nb
        real, dimension(ibx, nbx) :: xb, yb
        !
        ! Local variables
        !
        real, dimension(10) :: ainput
        logical :: error
        integer :: ib, k, n, ninput
        character(128) :: line
        character(20), save :: nchars
        real :: xbt, ybt
        !
        !*** End of declarations rewritten by SPAG
        !
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
        !
        !
        !---- permitted characters in string containing Fortran-readable numbers
        data nchars/'0123456789-+.,EDed  '/
        !
        !---- tab character is also permitted
        nchars(20:20) = char(9)
        !
        !
        !---- first assume that there will be a read error
        iftype = 0
        read (lu, 1000) line
        !---- if first line has any non-numeric character, go treat it as the name
        do k = 1, 80
            if (index(nchars, line(k:k))==0) goto 20
        enddo
        !
        !---- plain unlabeled file: rewind, and just read in X,Y coordinates
        name = ' '
        write (*, *)
        write (*, *) 'Reading plain coordinate file'
        iftype = 1
        goto 40
        !
        !---- first line interpreted as label string
        20   name = line(1:60)
        !      CALL STRIP(NAME,NNM)
        !
        !---- read second line
        read (lu, 1000) line
        ninput = 10
        call getflt(line, ainput, ninput, error)
        if (error) goto 99
        !
        if (ninput<4) then
            !------ no domain parameters: read X,Y coordinates directly with this line
            iftype = 2
        else
            !---- second line had parameters, read third line (coordinates line)
            read (lu, 1000) line
        endif
        !
        write (*, 1010) name
        1010 format (/1x, 'Reading airfoil: ', a32)
        !
        !
        !---- read in airfoil coordinates assuming first line has been read into
        !     character variable LINE
        40   do n = 1, nbx + 1
            do ib = 1, ibx + 1
                if (n==1 .and. ib==1) then
                    read (line, *, end = 56, err = 99) xbt, ybt
                else
                    read (lu, 1000, end = 56, err = 99) line
                    if (index(line, 'END')/=0) goto 56
                    read (line, *, end = 56, err = 99) xbt, ybt
                endif
                !
                if (xbt==999.0) then
                    !--------- save number of points and type for element which was just read in
                    nb(n) = ib - 1
                    !--------- go try to read next element
                    goto 55
                endif
                !
                if (n>nbx) then
                    write (*, *)                                              &
                            &'AREAD: Too many airfoil elements. Increase NBX.'
                    stop
                endif
                !
                !-------- uneventful XBT,YBT coordinates... just store them
                xb(ib, n) = xbt
                yb(ib, n) = ybt
                !
            enddo
            write (*, *) 'AREAD: Too many airfoil points.  Increase IBX.'
            stop
            !
        55   enddo
        n = nbx
        !
        56   if (ib==1) then
            !----- coordinate file has "999.0 999.0" at the end ...
            nbl = n - 1
        else
            !----- coordinate file has no special ending line
            nbl = n
            nb(n) = ib - 1
        endif
        !
        !---- set type if coordinate file this was
        if (iftype==0) then
            if (nbl==1) then
                iftype = 3
            else
                iftype = 4
            endif
        endif
        !
        return
        !
        99   write (*, 1100)
        iftype = 0
        return
        !...............................................................
        1000 format (a)
        1100 format (/' File READ error:  ')
    end subroutine areadnr
    !*==AWRITE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! AREADNR



    subroutine awrite(fname, lu, nel, ifrst, ilast, x, y, name, ispars, iftype)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: fname, ispars, name
        integer :: iftype, lu, nel
        integer, dimension(*) :: ifrst, ilast
        real, dimension(*) :: x, y
        !
        ! Local variables
        !
        integer :: i, iel, nf
        logical :: lopen
        !
        !*** End of declarations rewritten by SPAG
        !
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
        !
        lopen = fname(1:1)/=' '
        if (lopen) then
            nf = index(fname, ' ') + 1
            open (lu, file = fname, status = 'OLD', err = 98)
            rewind (lu)
        endif
        !
        if (iftype/=1) write (lu, 1000) name
        if (iftype==3 .or. iftype==4) write (lu, 1000) ispars
        !
        do iel = 1, nel
            do i = ifrst(iel), ilast(iel)
                write (lu, 1500) x(i), y(i)
            enddo
            !
            if (iel<nel) write (lu, 1550) 999.0, 999.0
        enddo
        !
        if (lopen) close (lu)
        return
        !
        98   write (*, 1050) fname(1:nf)
        return
        !...............................................................
        1000 format (a)
        1050 format (/' File OPEN error:  ', a)
        !
        1500 format (1x, 2f12.6)
        1550 format (1x, 2f6.1)
    end subroutine awrite
end module m_airio
