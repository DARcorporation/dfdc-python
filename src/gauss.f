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

      SUBROUTINE LUDCMP(NSIZ,N,A)
C     *******************************************************
C     *                                                     *
C     *   Factors a full NxN matrix A into an LU form.      *
C     *   Subr. BAKSUB can back-substitute it with some RHS.*
C     *   Assumes matrix is non-singular...                 *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   A is the matrix...                                *
C     *     ...replaced with its LU factors.                *
C     *                                                     *
C     *   Stolen from Numerical Recipes, removed pivoting.  *
C     *                                                     *
C     *                              Mark Drela  1988       *
C     *******************************************************
C
      DIMENSION A(NSIZ,NSIZ)
C
      DO 19 J=1, N
        DO 14 I=1, J-1
          SUM = A(I,J)
          DO 13 K=1, I-1
            SUM = SUM - A(I,K)*A(K,J)
   13     CONTINUE
          A(I,J) = SUM
   14   CONTINUE
C
        DO 16 I=J, N
          SUM = A(I,J)
          DO 15 K=1, J-1
            SUM = SUM - A(I,K)*A(K,J)
   15     CONTINUE
          A(I,J) = SUM
   16   CONTINUE
C
        DUM = 1.0/A(J,J)
        DO 18 I=J+1, N
          A(I,J) = A(I,J)*DUM
   18   CONTINUE
C
   19 CONTINUE
C
      RETURN
      END ! LUDCMP



      SUBROUTINE BAKSUB(NSIZ,N,A,B)
      DIMENSION A(NSIZ,NSIZ), B(NSIZ)
C
      DO II=1, N
        IF(B(II).NE.0.0) GO TO 5
      ENDDO
 5    CONTINUE
C
      DO 12 I=II+1, N
        SUM = B(I)
        DO 11 J=II, I-1
          SUM = SUM - A(I,J)*B(J)
   11   CONTINUE
        B(I) = SUM
   12 CONTINUE
C
      B(N) = B(N)/A(N,N)
C
      DO 14 I=N-1, 1, -1
        SUM = B(I)
        DO 13 J=I+1, N
          SUM = SUM - A(I,J)*B(J)
   13   CONTINUE
        B(I) = SUM/A(I,I)
   14 CONTINUE
C
      RETURN
      END ! BAKSUB
