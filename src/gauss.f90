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
    DIMENSION A(NSIZ, NSIZ)
    !
    DO 19 J = 1, N
        DO 14 I = 1, J - 1
            SUM = A(I, J)
            DO 13 K = 1, I - 1
                SUM = SUM - A(I, K) * A(K, J)
            13     CONTINUE
            A(I, J) = SUM
        14   CONTINUE
        !
        DO 16 I = J, N
            SUM = A(I, J)
            DO 15 K = 1, J - 1
                SUM = SUM - A(I, K) * A(K, J)
            15     CONTINUE
            A(I, J) = SUM
        16   CONTINUE
        !
        DUM = 1.0 / A(J, J)
        DO 18 I = J + 1, N
            A(I, J) = A(I, J) * DUM
        18   CONTINUE
        !
    19 CONTINUE
    !
    RETURN
END
! LUDCMP



SUBROUTINE BAKSUB(NSIZ, N, A, B)
    DIMENSION A(NSIZ, NSIZ), B(NSIZ)
    !
    DO II = 1, N
        IF(B(II).NE.0.0) GO TO 5
    ENDDO
    5    CONTINUE
    !
    DO 12 I = II + 1, N
        SUM = B(I)
        DO 11 J = II, I - 1
            SUM = SUM - A(I, J) * B(J)
        11   CONTINUE
        B(I) = SUM
    12 CONTINUE
    !
    B(N) = B(N) / A(N, N)
    !
    DO 14 I = N - 1, 1, -1
        SUM = B(I)
        DO 13 J = I + 1, N
            SUM = SUM - A(I, J) * B(J)
        13   CONTINUE
        B(I) = SUM / A(I, I)
    14 CONTINUE
    !
    RETURN
END
! BAKSUB