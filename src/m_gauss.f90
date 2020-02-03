module m_gauss
    implicit none
contains
    !*==LUDCMP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! FCOEFF
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

    subroutine ludcmp(nsiz, n, a)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n, nsiz
        real, dimension(nsiz, nsiz) :: a
        !
        ! Local variables
        !
        real :: dum, sum
        integer :: i, j, k
        !
        !*** End of declarations rewritten by SPAG
        !
        !     *******************************************************
        !     *                                                     *
        !     *   Factors a full NxN matrix A into an LU form.      *
        !     *   Subr. BAKSUB can back-substitute it with some RHS.*
        !     *   Assumes matrix is non-singular...                 *
        !     *    ...if it isn't, a divide by zero will result.    *
        !     *                                                     *
        !     *   A is the matrix...                                *
        !     *     ...replaced with its LU factors.                *
        !     *                                                     *
        !     *   Stolen from Numerical Recipes, removed pivoting.  *
        !     *                                                     *
        !     *                              Mark Drela  1988       *
        !     *******************************************************
        !
        !
        do j = 1, n
            do i = 1, j - 1
                sum = a(i, j)
                do k = 1, i - 1
                    sum = sum - a(i, k) * a(k, j)
                enddo
                a(i, j) = sum
            enddo
            !
            do i = j, n
                sum = a(i, j)
                do k = 1, j - 1
                    sum = sum - a(i, k) * a(k, j)
                enddo
                a(i, j) = sum
            enddo
            !
            dum = 1.0 / a(j, j)
            do i = j + 1, n
                a(i, j) = a(i, j) * dum
            enddo
            !
        enddo
        !
    end subroutine ludcmp
    !*==BAKSUB.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! LUDCMP



    subroutine baksub(nsiz, n, a, b)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n, nsiz
        real, dimension(nsiz, nsiz) :: a
        real, dimension(nsiz) :: b
        !
        ! Local variables
        !
        integer :: i, ii, j
        real :: sum
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        do ii = 1, n
            if (b(ii)/=0.0) exit
        enddo
        !
        do i = ii + 1, n
            sum = b(i)
            do j = ii, i - 1
                sum = sum - a(i, j) * b(j)
            enddo
            b(i) = sum
        enddo
        !
        b(n) = b(n) / a(n, n)
        !
        do i = n - 1, 1, -1
            sum = b(i)
            do j = i + 1, n
                sum = sum - a(i, j) * b(j)
            enddo
            b(i) = sum / a(i, i)
        enddo
        !
    end subroutine baksub
end module m_gauss
