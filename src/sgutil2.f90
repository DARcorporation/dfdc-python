module m_sgutil2
    implicit none
contains
    !*==ISGFIND.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETCRV
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
    
    FUNCTION ISGFIND(SS, S, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SS
        INTEGER :: ISGFIND
        REAL, DIMENSION(N) :: S
        !
        ! Local variables
        !
        INTEGER :: I, ILOW, IMID
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- Find closest value in S1 array to SS and return index
        !     Assumes array S in ascending value order and uses binary search to
        !     locate interval with S(I-1)<=SS<S(I)
        ILOW = 1
        I = N
        !
        DO WHILE (I - ILOW>1)
            !
            IMID = (I + ILOW) / 2
            IF (SS<S(IMID)) THEN
                I = IMID
            ELSE
                ILOW = IMID
            ENDIF
        ENDDO
        !
        !---- Check S(I-1) and S(I) for closest value to SS
        ISGFIND = I
        IF (I>1) THEN
            IF (ABS(SS - S(I))>ABS(SS - S(I - 1))) ISGFIND = I - 1
        ENDIF
        !
    END FUNCTION ISGFIND
    !*==SGCOPY.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE SGCOPY(S1, S2, N)
        !---- Copy N elements from S1 to S2
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(N) :: S1, S2
        !
        ! Local variables
        !
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        DO I = 1, N
            S2(I) = S1(I)
        ENDDO
        !
    END SUBROUTINE SGCOPY
    !*==SGCOPF.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE SGCOPF(S1, S2, N, SOFF1, SWT1, SOFF2, SWT2, FOFF, F1, FX1, X1, N1, &
            & F2, FX2, X2, N2)
        use m_spline, only: seval, deval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: FOFF, SOFF1, SOFF2, SWT1, SWT2
        INTEGER :: N, N1, N2
        REAL, DIMENSION(N1) :: F1, FX1, X1
        REAL, DIMENSION(N2) :: F2, FX2, X2
        REAL, DIMENSION(N) :: S1, S2
        !
        ! Local variables
        !
        REAL :: DELS1, DEV1, DS1, FEV1, FEV2, RES, RES_S1
        INTEGER :: I, ITER
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- reset array S1 so that  f1(s1) = f2(s2) + foff ,
        !-    with  s1 = S1*SWT1 + SOFF1  and  s2 = S2*SWT2 + SOFF2
        !
        DELS1 = S1(N) - S1(1)
        DO I = 2, N - 1
            FEV2 = SEVAL((S2(I) * SWT2 + SOFF2), F2, FX2, X2, N2) + FOFF
            !
            DO ITER = 1, 10
                FEV1 = SEVAL((S1(I) * SWT1 + SOFF1), F1, FX1, X1, N1)
                DEV1 = DEVAL((S1(I) * SWT1 + SOFF1), F1, FX1, X1, N1)
                RES = FEV1 - FEV2
                RES_S1 = DEV1 * SWT1
                DS1 = -RES / RES_S1
                S1(I) = S1(I) + DS1
                IF (ABS(DS1 / DELS1)<1.0E-5) GOTO 20
            ENDDO
            WRITE (*, *) 'SGCOPF: Convergence failed.  dS/Smax =', &
                    & DS1 / DELS1
            !
        20   ENDDO
        !
    END SUBROUTINE SGCOPF
    !*==SGAVG.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE SGAVG(S1, S2, N, C1)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: C1
        INTEGER :: N
        REAL, DIMENSION(N) :: S1, S2
        !
        ! Local variables
        !
        REAL :: F1, F2
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- average arrays S1 S2,
        !-    preserving S1 spacing at right point and S2 spacing at left endpoint
        DO I = 1, N
            F1 = FLOAT(I - 1) * C1
            F2 = FLOAT(N - I)
            S1(I) = (S1(I) * F1 + S2(I) * F2) / (F1 + F2)
            S2(I) = S1(I)
        ENDDO
        !
    END SUBROUTINE SGAVG
    !*==SGAVG1.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE SGAVG1(S1, S2, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(N) :: S1, S2
        !
        ! Local variables
        !
        REAL :: F1, F2
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- impose spacing of array S1 on S2 at right endpoint
        DO I = 1, N
            F1 = FLOAT(I - 1)
            F2 = FLOAT(N - I)
            S2(I) = (S1(I) * F1 + S2(I) * F2) / (F1 + F2)
        ENDDO
        !
    END SUBROUTINE SGAVG1
    !*==SGAVG2.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE SGAVG2(S1, S2, N)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(N) :: S1, S2
        !
        ! Local variables
        !
        REAL :: F1, F2
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- impose spacing of array S2 on S1 at left endpoint
        DO I = 1, N
            F1 = FLOAT(I - 1)
            F2 = FLOAT(N - I)
            S1(I) = (S1(I) * F1 + S2(I) * F2) / (F1 + F2)
        ENDDO
        !
    END SUBROUTINE SGAVG2
    !*==SGSHFT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE SGSHFT(SBEG, SEND, SOLD, SNEW, S1, S2, N)
        !
        !---- Shift S1 values between limits SBEG and SEND such that
        !     value at SOLD is shifted to SNEW in S2.  S1 values
        !     between SBEG and SOLD are scaled to correspond to
        !     the range SBEG to SNEW. S1 values between SOLD and
        !     SEND are scaled to the range SNEW to SEND.
        !
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SBEG, SEND, SNEW, SOLD
        REAL, DIMENSION(N) :: S1, S2
        !
        ! Local variables
        !
        REAL :: FAC1, FAC2
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        IF (SBEG>=SEND) THEN
            WRITE (*, *) 'SGSHFT: Bad input range SBEG,SEND ', SBEG, SEND
            STOP
        ENDIF
        !
        IF (SBEG>SOLD .OR. SEND<SOLD) THEN
            WRITE (*, *) 'SGSHFT: Bad input range SBEG/SOLD/SEND ', SBEG, &
                    & SOLD, SEND
            STOP
        ENDIF
        !
        IF (SBEG>SNEW .OR. SEND<SNEW) THEN
            WRITE (*, *) 'SGSHFT: Bad input range SBEG/SNEW/SEND ', SBEG, &
                    & SNEW, SEND
            STOP
        ENDIF
        !
        !---- Shift S1 values in SBEG-SEND range to new values such
        !     that SOLD becomes SNEW (linearly interpolated from both
        !     ends of range).
        FAC1 = (SNEW - SBEG) / (SOLD - SBEG)
        FAC2 = (SEND - SNEW) / (SEND - SOLD)
        DO I = 1, N
            IF (S1(I)<=SBEG) THEN
                S2(I) = S1(I)
            ELSEIF (S1(I)>SBEG .AND. S1(I)<=SOLD) THEN
                S2(I) = SBEG + (S1(I) - SBEG) * FAC1
            ELSEIF (S1(I)>SOLD .AND. S1(I)<SEND) THEN
                S2(I) = SEND - (SEND - S1(I)) * FAC2
            ELSE
                S2(I) = S1(I)
            ENDIF
        ENDDO
        !
    END SUBROUTINE SGSHFT
    !*==SGRENUM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE SGRENUM(S, N1, N2)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        INTEGER, PARAMETER :: NMAX = 2000
        !
        ! Dummy arguments
        !
        INTEGER :: N1, N2
        REAL, DIMENSION(2000) :: S
        !
        ! Local variables
        !
        REAL, DIMENSION(NMAX) :: C, X, XI
        REAL :: CX1, CX2, FRAC, RI, T
        INTEGER :: I, J
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Interpolates input array S(1:N1) into new number of points 1:N2.
        !     The interpolation is a spline in the array index.
        !     The first and last elements in S will remain the same.
        !----------------------------------------------------------------------
        !
        IF (N1>NMAX) STOP 'NEWNUM:  Array overflow'
        !
        !---- save old input array
        DO I = 1, N1
            X(I) = S(I)
        ENDDO
        !
        !---- spline X(i)  (set up and solve tridiagonal system on the fly)
        C(1) = 0.5
        XI(1) = 1.5 * (X(2) - X(1))
        DO I = 2, N1 - 1
            C(I) = 1.0 / (4.0 - C(I - 1))
            XI(I) = (3.0 * (X(I + 1) - X(I - 1)) - XI(I - 1)) * C(I)
        ENDDO
        I = N1
        XI(I) = (3.0 * (X(I) - X(I - 1)) - XI(I - 1)) / (2.0 - C(I - 1))
        !
        DO I = N1 - 1, 1, -1
            XI(I) = XI(I) - C(I) * XI(I + 1)
        ENDDO
        !
        !---- evaluate s(i) spline at new points
        DO J = 1, N2
            FRAC = FLOAT(J - 1) / FLOAT(N2 - 1)
            RI = 1.0 + FRAC * FLOAT(N1 - 1)
            I = MIN(INT(RI), N1 - 1)
            !
            T = RI - FLOAT(I)
            CX1 = XI(I) - X(I + 1) + X(I)
            CX2 = XI(I + 1) - X(I + 1) + X(I)
            S(J) = T * X(I + 1) + (1.0 - T) * X(I) + (T - T * T) * ((1.0 - T) * CX1 - T * CX2)
        ENDDO
        !
        !---- make sure new endpoints are exactly the same as old ones
        S(1) = X(1)
        S(N2) = X(N1)
        !
    END SUBROUTINE SGRENUM
end module m_sgutil2
