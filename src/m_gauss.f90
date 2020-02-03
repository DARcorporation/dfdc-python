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

    SUBROUTINE LUDCMP(NSIZ, N, A)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N, NSIZ
        REAL, DIMENSION(NSIZ, NSIZ) :: A
        !
        ! Local variables
        !
        REAL :: DUM, SUM
        INTEGER :: I, J, K
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
        DO J = 1, N
            DO I = 1, J - 1
                SUM = A(I, J)
                DO K = 1, I - 1
                    SUM = SUM - A(I, K) * A(K, J)
                ENDDO
                A(I, J) = SUM
            ENDDO
            !
            DO I = J, N
                SUM = A(I, J)
                DO K = 1, J - 1
                    SUM = SUM - A(I, K) * A(K, J)
                ENDDO
                A(I, J) = SUM
            ENDDO
            !
            DUM = 1.0 / A(J, J)
            DO I = J + 1, N
                A(I, J) = A(I, J) * DUM
            ENDDO
            !
        ENDDO
        !
    END SUBROUTINE LUDCMP
    !*==BAKSUB.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! LUDCMP



    SUBROUTINE BAKSUB(NSIZ, N, A, B)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N, NSIZ
        REAL, DIMENSION(NSIZ, NSIZ) :: A
        REAL, DIMENSION(NSIZ) :: B
        !
        ! Local variables
        !
        INTEGER :: I, II, J
        REAL :: SUM
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        DO II = 1, N
            IF (B(II)/=0.0) EXIT
        ENDDO
        !
        DO I = II + 1, N
            SUM = B(I)
            DO J = II, I - 1
                SUM = SUM - A(I, J) * B(J)
            ENDDO
            B(I) = SUM
        ENDDO
        !
        B(N) = B(N) / A(N, N)
        !
        DO I = N - 1, 1, -1
            SUM = B(I)
            DO J = I + 1, N
                SUM = SUM - A(I, J) * B(J)
            ENDDO
            B(I) = SUM / A(I, I)
        ENDDO
        !
    END SUBROUTINE BAKSUB
end module m_gauss
