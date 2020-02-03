module m_aero
    implicit none
contains
    !*==SETIAERO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
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
    !--- Aero data stored for one or more radial aerodynamic sections
    !
    !-- aero data quantities for each defined radial aerodynamic section
    !  NAERO     Number of aerodynamic datasets defined (NAERO>=1)
    !  XIAERO    Radial station r/R where aero dataset is defined
    !  AERODATA  Aerodynamic definition of the blade section at XIAERO
    !            AERODATA( 1,x) = A0 (angle of zero lift)
    !            AERODATA( 2,x) = CLMAX (Max CL)
    !            AERODATA( 3,x) = CLMIN (Min CL)
    !            AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
    !            AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
    !            AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
    !            AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
    !            AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
    !            AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
    !            AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
    !            AERODATA(11,x) = REREF (reference Reynold's number)
    !            AERODATA(12,x) = REXP (Reynold's number exponent Cd~Re^REXP)
    !            AERODATA(13,x) = MCRIT (critical Mach #)
    !            AERODATA(14,x) = TOC (thickness/chord)
    !            AERODATA(15,x) = DCDCL2S (Scndary, annulus drag param d(Cd)/dCL^2)
    !=========================================================================
    !
    !     Version 070-ES1
    !     Philip Carter, Esotec Developments, February 2009
    !     philip (at) esotec (dot) org
    !
    !     Changes from 0.70:
    !
    !     CL/CD plotting fixed (missing argument).
    !     READ and WRIT commands fixed (LU).
    !     FLIP command to toggle to/from neg-BGam parameters.
    !     DISP formats repaired. Message when set for neg-BGam.
    !     Interpolated A0 (AZERO) stored by GETCLCDCM.
    !     Multi-Re plotting with constant Mach or vice-versa (AEROPLT2).
    !     Modified plotting interface.
    !     All plotting functionality duplicated in EDIT.
    !     PLOTMACH and PLOTREYN subroutines to avoid code duplication.
    !     HARD and ANNO fixed.
    !     Disk index displayed with data for multi-disk cases.
    !     Disk index added to plot titles for multi-disk cases.
    !     Mach constant or Re constant written to plot legends.
    !     Various fixes to control structure and cosmetics.
    !
    !=========================================================================
    !
    !
    
    SUBROUTINE SETIAERO
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: I, N, NR
        REAL :: XI
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Sets up indices referring to aero section for
        !     each radial station
        !--------------------------------------------------
        !
        !--- Find lower index of aero data sections XIAERO(N) bounding XI=YRC/RTIP
        DO NR = 1, NROTOR
            DO I = 1, NRC
                IAERO(I, NR) = 1
                DO N = 1, NAERO(NR)
                    XI = YRC(I, NR) / RTIP(NR)
                    IF (XIAERO(N, NR)<=XI) IAERO(I, NR) = N
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE SETIAERO
    !*==GETAERO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE GETAERO(NR, N, XISECT, A0, CLMAX, CLMIN, DCLDA, DCLDA_STALL, &
            & DCL_STALL, CDMIN, CLDMIN, DCDCL2, DCDCL2S, CMCON, &
            & MCRIT, TOC, REREF, REXP)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: A0, CDMIN, CLDMIN, CLMAX, CLMIN, CMCON, DCDCL2, &
                & DCDCL2S, DCLDA, DCLDA_STALL, DCL_STALL, MCRIT, &
                & REREF, REXP, TOC, XISECT
        INTEGER :: N, NR
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------
        !     Gets aero data from stored section array
        !
        !   AERODATA    Aerodynamic definition of the blade section at XIAERO
        !               AERODATA( 1,x) = A0 (angle of zero lift)
        !               AERODATA( 2,x) = CLMAX (Max CL)
        !               AERODATA( 3,x) = CLMIN (Min CL)
        !               AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
        !               AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
        !               AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
        !               AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
        !               AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
        !               AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
        !               AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
        !               AERODATA(11,x) = REREF (reference Reynold's number)
        !               AERODATA(12,x) = REXP (Reynold's number exponent Cd~Re^REXP)
        !               AERODATA(13,x) = MCRIT (critical Mach #)
        !               AERODATA(14,x) = TOC (thickness/chord)
        !               AERODATA(15,x) = DCDCL2S (Secondary, annulus drag param d(Cd)/dCL^2)
        !---------------------------------------------
        !
        IF (NR<1 .OR. NR>NROTOR) THEN
            WRITE (*, *) 'Error: blade index of aero section out of bounds'
            RETURN
        ENDIF
        IF (N<1 .OR. N>NAERO(NR)) THEN
            WRITE (*, *) 'Error: index of aero section out of bounds'
            RETURN
        ENDIF
        !
        A0 = AERODATA(1, N, NR)
        CLMAX = AERODATA(2, N, NR)
        CLMIN = AERODATA(3, N, NR)
        DCLDA = AERODATA(4, N, NR)
        DCLDA_STALL = AERODATA(5, N, NR)
        DCL_STALL = AERODATA(6, N, NR)
        CDMIN = AERODATA(7, N, NR)
        CLDMIN = AERODATA(8, N, NR)
        DCDCL2 = AERODATA(9, N, NR)
        CMCON = AERODATA(10, N, NR)
        REREF = AERODATA(11, N, NR)
        REXP = AERODATA(12, N, NR)
        MCRIT = AERODATA(13, N, NR)
        TOC = AERODATA(14, N, NR)
        DCDCL2S = AERODATA(15, N, NR)
        XISECT = XIAERO(N, NR)
        !
    END SUBROUTINE GETAERO
    !*==PUTAERO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE PUTAERO(NR, N, XISECT, A0, CLMAX, CLMIN, DCLDA, DCLDA_STALL, &
            & DCL_STALL, CDMIN, CLDMIN, DCDCL2, DCDCL2S, CMCON, &
            & MCRIT, TOC, REREF, REXP)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: A0, CDMIN, CLDMIN, CLMAX, CLMIN, CMCON, DCDCL2, &
                & DCDCL2S, DCLDA, DCLDA_STALL, DCL_STALL, MCRIT, &
                & REREF, REXP, TOC, XISECT
        INTEGER :: N, NR
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------
        !     Puts aero data into stored section array at index N
        !
        !   AERODATA    Aerodynamic definition of the blade section at XIAERO
        !               AERODATA( 1,x) = A0 (angle of zero lift)
        !               AERODATA( 2,x) = CLMAX (Max CL)
        !               AERODATA( 3,x) = CLMIN (Min CL)
        !               AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
        !               AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
        !               AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
        !               AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
        !               AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
        !               AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
        !               AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
        !               AERODATA(11,x) = REREF (reference Reynold's number)
        !               AERODATA(12,x) = REXP (Reynold's number exponent Cd~Re^REXP)
        !               AERODATA(13,x) = MCRIT (critical Mach #)
        !               AERODATA(14,x) = TOC (thickness/chord)
        !               AERODATA(15,x) = DCDCL2S (Secondary, annulus drag param d(Cd)/dCL^2)
        !--------------------------------------------------------
        !
        IF (NR<1 .OR. NR>NRX) THEN
            WRITE (*, *) 'Error: blade index of aero section out of bounds'
            RETURN
        ENDIF
        IF (N<1) THEN
            WRITE (*, *) 'Error: index of aero section out of bounds'
            RETURN
        ENDIF
        IF (N>NAX) THEN
            WRITE (*, *) 'Too many aero sections defined...'
            RETURN
        ENDIF
        !
        AERODATA(1, N, NR) = A0
        AERODATA(2, N, NR) = CLMAX
        AERODATA(3, N, NR) = CLMIN
        AERODATA(4, N, NR) = DCLDA
        AERODATA(5, N, NR) = DCLDA_STALL
        AERODATA(6, N, NR) = DCL_STALL
        AERODATA(7, N, NR) = CDMIN
        AERODATA(8, N, NR) = CLDMIN
        AERODATA(9, N, NR) = DCDCL2
        AERODATA(10, N, NR) = CMCON
        AERODATA(11, N, NR) = REREF
        AERODATA(12, N, NR) = REXP
        AERODATA(13, N, NR) = MCRIT
        AERODATA(14, N, NR) = TOC
        AERODATA(15, N, NR) = DCDCL2S
        XIAERO(N, NR) = XISECT
        !
    END SUBROUTINE PUTAERO
    !*==SORTAR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE SORTAR(NS, S, W, NDIM)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: NDIM, NS
        REAL, DIMENSION(NS) :: S
        REAL, DIMENSION(NDIM, NS) :: W
        !
        ! Local variables
        !
        LOGICAL :: DONE
        INTEGER :: IPASS, L, N, NP
        REAL :: TEMP
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !---- sort arrays by S values
        !     Orders data monotonically increasing in S(i)
        !----------------------------------------------------
        !
        DO IPASS = 1, 500
            DONE = .TRUE.
            DO N = 1, NS - 1
                NP = N + 1
                IF (S(NP)<S(N)) THEN
                    TEMP = S(NP)
                    S(NP) = S(N)
                    S(N) = TEMP
                    DO L = 1, NDIM
                        TEMP = W(L, NP)
                        W(L, NP) = W(L, N)
                        W(L, N) = TEMP
                    ENDDO
                    DONE = .FALSE.
                ENDIF
            ENDDO
            IF (DONE) GOTO 10
        ENDDO
        STOP 'SORTAR failed'
        !
    10   END SUBROUTINE SORTAR
    !*==GETCLCDCM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SORTAR
    
    
    !*************************************************************************
    !  Interpolated aero section properties functions
    !  These routines implement a functional representation of the
    !  blade aero properties (CL,CD,CM) vs ALFA
    !*************************************************************************
    
    
    SUBROUTINE GETCLCDCM(NR, IS, XI, ALF, W, REY, SECSIG, SECSTAGR, CLIFT, &
            & CL_ALF, CL_W, CLMAX, CLMIN, DCL_STALL, STALLF, &
            & CDRAG, CD_ALF, CD_W, CD_REY, CMOM, CM_AL, CM_W)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: ALF, CDRAG, CD_ALF, CD_REY, CD_W, CLIFT, CLMAX, &
                & CLMIN, CL_ALF, CL_W, CMOM, CM_AL, CM_W, &
                & DCL_STALL, REY, SECSIG, SECSTAGR, W, XI
        INTEGER :: IS, NR
        LOGICAL :: STALLF
        !
        ! Local variables
        !
        REAL :: A0, A02, CDMIN, CDRAG2, CD_ALF2, CD_REY2, CD_W2, &
                & CLDMIN, CLIFT2, CLMAX2, CLMIN2, CL_ALF2, CL_W2, &
                & CMCON, CMOM2, CM_AL2, CM_W2, DCDCL2, DCDCL2S, &
                & DCLDA, DCLDA_STALL, DCL_STALL2, FRAC, MCRIT, &
                & REREF, REXP, TOC, XISECT1, XISECT2
        INTEGER :: N
        LOGICAL :: STALLF2
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     CL(alpha),
        !      CD(alpha),
        !       CM(alpha) interpolation function for blade at station IS at XI=r/R
        !-------------------------------------------------------------
        !
        IF (XI<0.0 .OR. XI>1.0) WRITE (*, *)                             &
                &'Undefined section XI in GETCLCDCM '&
                &, XI
        !
        !--- Check for installed aero data section index
        N = IAERO(IS, NR)
        IF (N<1 .OR. N>NAERO(NR)) THEN
            !
            IF (NAERO(NR)>1) THEN
                !--- Find lower index of aero data sections XIAERO(N) bounding XI
                DO N = 1, NAERO(NR)
                    IF (XIAERO(N, NR)>XI) GOTO 10
                    !c          write(*,*) 'getcl iaero= ',N,' is= ',is,xiaero(N),xi
                    IAERO(IS, NR) = N
                ENDDO
                WRITE (*, *) 'Aero section not found for station ', XI
            ENDIF
            !
            N = 1
            IAERO(IS, NR) = N
        ENDIF
        !
        !--- Get section aero data from stored section array
        10   A0 = AERODATA(1, N, NR)
        CLMAX = AERODATA(2, N, NR)
        CLMIN = AERODATA(3, N, NR)
        DCLDA = AERODATA(4, N, NR)
        DCLDA_STALL = AERODATA(5, N, NR)
        DCL_STALL = AERODATA(6, N, NR)
        CDMIN = AERODATA(7, N, NR)
        CLDMIN = AERODATA(8, N, NR)
        DCDCL2 = AERODATA(9, N, NR)
        CMCON = AERODATA(10, N, NR)
        REREF = AERODATA(11, N, NR)
        REXP = AERODATA(12, N, NR)
        MCRIT = AERODATA(13, N, NR)
        TOC = AERODATA(14, N, NR)
        DCDCL2S = AERODATA(15, N, NR)
        XISECT1 = XIAERO(N, NR)
        !
        !--- Get data for inner bounding aero section
        CALL CLCDCM(ALF, W, REY, VSO, SECSIG, SECSTAGR, CLIFT, CL_ALF, CL_W, &
                & STALLF, CDRAG, CD_ALF, CD_W, CD_REY, CMOM, CM_AL, CM_W, A0, &
                & CLMAX, CLMIN, DCLDA, DCLDA_STALL, DCL_STALL, CDMIN, CLDMIN, &
                & DCDCL2, CMCON, MCRIT, REREF, REXP, TOC, DCDCL2S)
        !
        !--- Check for another bounding section, if not we are done,
        !    if we have another section linearly interpolate data to station IS
        IF (N<NAERO(NR)) THEN
            XISECT2 = XIAERO(N + 1, NR)
            FRAC = (XI - XISECT1) / (XISECT2 - XISECT1)
            IF (FRAC<=0.0 .OR. FRAC>1.0) THEN
                !c         write(*,*) 'CL n,is,xi,frac = ',n,is,xi(is),frac
            ENDIF
            !
            !--- A02 sustituted for A0 in the following (2 places),
            !    to permit A0 interpolation for storage in AZERO
            !
            A02 = AERODATA(1, N + 1, NR)
            CLMAX2 = AERODATA(2, N + 1, NR)
            CLMIN2 = AERODATA(3, N + 1, NR)
            DCLDA = AERODATA(4, N + 1, NR)
            DCLDA_STALL = AERODATA(5, N + 1, NR)
            DCL_STALL2 = AERODATA(6, N + 1, NR)
            CDMIN = AERODATA(7, N + 1, NR)
            CLDMIN = AERODATA(8, N + 1, NR)
            DCDCL2 = AERODATA(9, N + 1, NR)
            CMCON = AERODATA(10, N + 1, NR)
            REREF = AERODATA(11, N + 1, NR)
            REXP = AERODATA(12, N + 1, NR)
            MCRIT = AERODATA(13, N + 1, NR)
            TOC = AERODATA(14, N + 1, NR)
            DCDCL2S = AERODATA(15, N + 1, NR)
            !
            !--- Get data for outer bounding aero section
            CALL CLCDCM(ALF, W, REY, VSO, SECSIG, SECSTAGR, CLIFT2, CL_ALF2, CL_W2, &
                    & STALLF2, CDRAG2, CD_ALF2, CD_W2, CD_REY2, CMOM2, CM_AL2, &
                    & CM_W2, A02, CLMAX2, CLMIN2, DCLDA, DCLDA_STALL, &
                    & DCL_STALL2, CDMIN, CLDMIN, DCDCL2, CMCON, MCRIT, REREF, &
                    & REXP, TOC, DCDCL2S)
            !--- Interpolate aero data to blade station
            STALLF = STALLF .OR. STALLF2
            CLIFT = (1.0 - FRAC) * CLIFT + FRAC * CLIFT2
            CL_ALF = (1.0 - FRAC) * CL_ALF + FRAC * CL_ALF2
            CL_W = (1.0 - FRAC) * CL_W + FRAC * CL_W2
            CLMAX = (1.0 - FRAC) * CLMAX + FRAC * CLMAX2
            CLMIN = (1.0 - FRAC) * CLMIN + FRAC * CLMIN2
            DCL_STALL = (1.0 - FRAC) * DCL_STALL + FRAC * DCL_STALL2
            !
            CMOM = (1.0 - FRAC) * CMOM + FRAC * CMOM2
            CM_AL = (1.0 - FRAC) * CM_AL + FRAC * CM_AL2
            CM_W = (1.0 - FRAC) * CM_W + FRAC * CM_W2
            !
            CDRAG = (1.0 - FRAC) * CDRAG + FRAC * CDRAG2
            CD_ALF = (1.0 - FRAC) * CD_ALF + FRAC * CD_ALF2
            CD_W = (1.0 - FRAC) * CD_W + FRAC * CD_W2
            CD_REY = (1.0 - FRAC) * CD_REY + FRAC * CD_REY2
            A0 = (1.0 - FRAC) * A0 + FRAC * A02
        ENDIF
        !
        AZERO(IS, NR) = A0
        !
    END SUBROUTINE GETCLCDCM
    !*==GETALF.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    !GETCLCDCM
    
    
    
    SUBROUTINE GETALF(NR, IS, XI, SECSIG, SECSTAGR, CLIFT, W, ALF, ALF_CL, &
            & ALF_W, STALLF)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: ALF, ALF_CL, ALF_W, CLIFT, SECSIG, SECSTAGR, W, &
                & XI
        INTEGER :: IS, NR
        LOGICAL :: STALLF
        !
        ! Local variables
        !
        REAL :: A0, CDRAG, CD_ALF, CD_REY, CD_W, CLMAX, CLMIN, &
                & CLTEMP, CL_ALF, CL_W, CMOM, CM_AL, CM_W, DALF, &
                & DCL_STALL, REY
        REAL, SAVE :: EPS
        INTEGER :: ITER
        INTEGER, SAVE :: NITER
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------
        !     Inverse alpha(CL) function
        !     Uses Newton-Raphson iteration to get ALF from CL function
        !------------------------------------------------------------
        DATA NITER/10/
        DATA EPS/1.0E-5/
        !
        STALLF = .FALSE.
        !
        !---HHY A0 is now an aero section property
        A0 = AERODATA(1, IS, NR)
        REY = 0.0
        !
        ALF = A0
        DO ITER = 1, NITER
            CALL GETCLCDCM(NR, IS, XI, ALF, W, REY, SECSIG, SECSTAGR, CLTEMP, &
                    & CL_ALF, CL_W, CLMAX, CLMIN, DCL_STALL, STALLF, CDRAG, &
                    & CD_ALF, CD_W, CD_REY, CMOM, CM_AL, CM_W)
            !c      IF(STALLF) GO TO 20
            DALF = -(CLTEMP - CLIFT) / CL_ALF
            ALF = ALF + DALF
            ALF_CL = 1.0 / CL_ALF
            ALF_W = -CL_W / CL_ALF
            IF (ABS(DALF)<EPS) RETURN
        ENDDO
        !
        WRITE (*, *) 'GETALF: alpha(CL) function inversion failed'
        !      write(*,*) 'is,clift  ',is,clift
        !      write(*,*) 'abs(dalf) ',abs(dalf)
        !      write(*,*) 'cl_alf    ',cl_alf
        !
    END SUBROUTINE GETALF
    !*==CLCDCM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GETALF
    
    
    
    !*************************************************************************
    !  Basic aero section properties functions
    !  These routines implement a functional representation of the
    !  blade section aero properties (CL,CD,CM) vs ALF
    !*************************************************************************
    
    SUBROUTINE CLCDCM(ALF, W, REY, VSO, SECSIG, SECSTAGR, CLIFT, CL_ALF, CL_W, &
            & STALLF, CDRAG, CD_ALF, CD_W, CD_REY, CMOM, CM_AL, CM_W, &
            & A0, CLMAX, CLMIN, DCLDA, DCLDA_STALL, DCL_STALL, &
            & CDMIN, CLDMIN, DCDCL2, CMCON, MCRIT, REREF, REXP, TOC, &
            & DCDCL2S)
        !------------------------------------------------------------
        !     CL(alpha) function
        !     Note that in addition to setting CLIFT and its derivatives
        !     CLMAX and CLMIN (+ and - stall CL's) are set in this routine
        !     In the compressible range the stall CL is reduced by a factor
        !     proportional to Mcrit-Mach.  Stall limiting for compressible
        !     cases begins when the compressible drag added CDC > CDMstall
        !------------------------------------------------------------
        !     CD(alpha) function - presently CD is assumed to be a sum
        !     of profile drag + stall drag + compressibility drag
        !     In the linear lift range drag is CD0 + quadratic function of CL-CLDMIN
        !     In + or - stall an additional drag is added that is proportional
        !     to the extent of lift reduction from the linear lift value.
        !     Compressible drag is based on adding drag proportional to
        !     (Mach-Mcrit_eff)^MEXP
        !------------------------------------------------------------
        !     CM(alpha) function - presently CM is assumed constant,
        !     varying only with Mach by Prandtl-Glauert scaling
        !------------------------------------------------------------
        !
        !C    INCLUDE 'DFDC.inc'
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: A0, ALF, CDMIN, CDRAG, CD_ALF, CD_REY, CD_W, &
                & CLDMIN, CLIFT, CLMAX, CLMIN, CL_ALF, CL_W, CMCON, &
                & CMOM, CM_AL, CM_W, DCDCL2, DCDCL2S, DCLDA, &
                & DCLDA_STALL, DCL_STALL, MCRIT, REREF, REXP, REY, &
                & SECSIG, SECSTAGR, TOC, VSO, W
        LOGICAL :: STALLF
        !
        ! Local variables
        !
        REAL :: CDC, CDCL2, CDC_ALF, CDC_W, CDMDD, CDMFACTOR, &
                & CDMSTALL, CLA, CLA_ALF, CLA_W, CLFACTOR, CLLIM, &
                & CLLIM_CLA, CLMAXM, CLMFACTOR, CLMINM, CLMN, CLMX, &
                & CRITMACH, CRITMACH_ALF, CRITMACH_W, DCD, DCDX, &
                & DCD_ALF, DCD_W, DMDD, DMSTALL, FAC, FAC_W, &
                & FSTALL, MACH, MACH_W, MEXP, MSQ, MSQ_W, PGRT, &
                & PGRT_W, RCORR, RCORR_REY
        REAL(R8KIND) :: ECMAX, ECMIN
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- Factors for compressibility drag model, HHY 10/23/00
        !     Mcrit is set by user
        !     Effective Mcrit is Mcrit_eff = Mcrit - CLMFACTOR*(CL-CLDmin) - DMDD
        !     DMDD is the delta Mach to get CD=CDMDD (usually 0.0020)
        !     Compressible drag is CDC = CDMFACTOR*(Mach-Mcrit_eff)^MEXP
        !     CDMstall is the drag at which compressible stall begins
        !
        CDMFACTOR = 10.0
        CLMFACTOR = 0.25
        MEXP = 3.0
        CDMDD = 0.0020
        CDMSTALL = 0.1000
        !
        !---- Prandtl-Glauert compressibility factor
        MSQ = W * W / VSO**2
        MSQ_W = 2.0 * W / VSO**2
        IF (MSQ>=1.0) THEN
            WRITE (*, *) 'CLFUNC: Local Mach number limited to 0.99, was ', &
                    & MSQ
            MSQ = 0.99
            MSQ_W = 0.
        ENDIF
        PGRT = 1.0 / SQRT(1.0 - MSQ)
        PGRT_W = 0.5 * MSQ_W * PGRT**3
        !
        !---- Mach number and dependence on velocity
        MACH = SQRT(MSQ)
        MACH_W = 0.0
        IF (MACH/=0.0) MACH_W = 0.5 * MSQ_W / MACH
        !
        !------------------------------------------------------------
        !--- Generate CLFACTOR for cascade effects from section solidity
        CLFACTOR = 1.0
        IF (SECSIG>0.0) CALL GETCLFACTOR(SECSIG, SECSTAGR, CLFACTOR)
        !
        !------------------------------------------------------------
        !--- Generate CL from dCL/dAlpha and Prandtl-Glauert scaling
        CLA = DCLDA * PGRT * (ALF - A0) * CLFACTOR
        CLA_ALF = DCLDA * PGRT * CLFACTOR
        CLA_W = DCLDA * PGRT_W * (ALF - A0) * CLFACTOR
        !
        !ccccccccc
        !      WRITE(*,*)'CL Factor   ALF   A0   DCLDA  CLA'
        !      WRITE(*,*) CLFACTOR,ALF,A0,DCLDA,CLA
        !
        !C--- Effective CLmax is limited by Mach effects
        !    reduces CLmax to match the CL of onset of serious compressible drag
        CLMX = CLMAX
        CLMN = CLMIN
        DMSTALL = (CDMSTALL / CDMFACTOR)**(1.0 / MEXP)
        CLMAXM = MAX(0.0, (MCRIT + DMSTALL - MACH) / CLMFACTOR) + CLDMIN
        CLMAX = MIN(CLMAX, CLMAXM)
        CLMINM = MIN(0.0, -(MCRIT + DMSTALL - MACH) / CLMFACTOR) + CLDMIN
        CLMIN = MAX(CLMIN, CLMINM)
        !
        !--- CL limiter function (turns on after +-stall
        ECMAX = DEXP(MIN(200.0D0, DBLE((CLA - CLMAX) / DCL_STALL)))
        ECMIN = DEXP(MIN(200.0D0, DBLE((CLMIN - CLA) / DCL_STALL)))
        CLLIM = DCL_STALL * DLOG((1.0D0 + ECMAX) / (1.0D0 + ECMIN))
        CLLIM_CLA = ECMAX / (1.0 + ECMAX) + ECMIN / (1.0 + ECMIN)
        !
        !      if(CLLIM.GT.0.001) then
        !      write(*,999) 'cla,cllim,ecmax,ecmin ',cla,cllim,ecmax,ecmin
        !      endif
        ! 999  format(a,2(1x,f10.6),3(1x,d12.6))
        !
        !--- Subtract off a (nearly unity) fraction of the limited CL function
        !    This sets the dCL/dAlpha in the stalled regions to 1-FSTALL of that
        !    in the linear lift range
        FSTALL = DCLDA_STALL / DCLDA
        CLIFT = CLA - (1.0 - FSTALL) * CLLIM
        CL_ALF = CLA_ALF - (1.0 - FSTALL) * CLLIM_CLA * CLA_ALF
        CL_W = CLA_W - (1.0 - FSTALL) * CLLIM_CLA * CLA_W
        !
        STALLF = .FALSE.
        IF (CLIFT>CLMAX) STALLF = .TRUE.
        IF (CLIFT<CLMIN) STALLF = .TRUE.
        !
        !
        !------------------------------------------------------------
        !--- CM from CMCON and Prandtl-Glauert scaling
        CMOM = PGRT * CMCON
        CM_AL = 0.0
        CM_W = PGRT_W * CMCON
        !
        !
        !------------------------------------------------------------
        !--- CD from profile drag, stall drag and compressibility drag
        !
        !---- Reynolds number scaling factor
        IF (REY<=0) THEN
            RCORR = 1.0
            RCORR_REY = 0.0
        ELSE
            RCORR = (REY / REREF)**REXP
            RCORR_REY = REXP / REY
        ENDIF
        !
        !--- Include quadratic lift drag terms from airfoil and annulus
        !
        !      CDCL2 = DCDCL2 + DCDCL2S
        CDCL2 = DCDCL2
        ! no chance of getting messed up...
        !
        !--- In the basic linear lift range drag is a function of lift
        !    CD = CD0 (constant) + quadratic with CL)
        CDRAG = (CDMIN + CDCL2 * (CLIFT - CLDMIN)**2) * RCORR
        CD_ALF = (2.0 * CDCL2 * (CLIFT - CLDMIN) * CL_ALF) * RCORR
        CD_W = (2.0 * CDCL2 * (CLIFT - CLDMIN) * CL_W) * RCORR
        CD_REY = CDRAG * RCORR_REY
        !
        !--- Post-stall drag added
        FSTALL = DCLDA_STALL / DCLDA
        DCDX = (1.0 - FSTALL) * CLLIM / (PGRT * DCLDA)
        !      write(*,*) 'cla,cllim,fstall,pg,dclda ',cla,cllim,fstall,pg,dclda
        DCD = 2.0 * DCDX**2
        DCD_ALF = 4.0 * DCDX * (1.0 - FSTALL) * CLLIM_CLA * CLA_ALF / (PGRT * DCLDA)
        DCD_W = 4.0 * DCDX * ((1.0 - FSTALL) * CLLIM_CLA * CLA_W / (PGRT * DCLDA)       &
                & - DCD / PGRT * PGRT_W)
        !      write(*,*) 'alf,cl,dcd,dcd_alf,dcd_w ',alf,clift,dcd,dcd_alf,dcd_w
        !
        !--- Compressibility drag (accounts for drag rise above Mcrit with CL effects
        !    CDC is a function of a scaling factor*(M-Mcrit(CL))**MEXP
        !    DMDD is the Mach difference corresponding to CD rise of CDMDD at MCRIT
        DMDD = (CDMDD / CDMFACTOR)**(1.0 / MEXP)
        CRITMACH = MCRIT - CLMFACTOR * ABS(CLIFT - CLDMIN) - DMDD
        CRITMACH_ALF = -CLMFACTOR * ABS(CL_ALF)
        CRITMACH_W = -CLMFACTOR * ABS(CL_W)
        IF (MACH<CRITMACH) THEN
            CDC = 0.0
            CDC_ALF = 0.0
            CDC_W = 0.0
        ELSE
            CDC = CDMFACTOR * (MACH - CRITMACH)**MEXP
            CDC_W = MEXP * MACH_W * CDC / MACH - MEXP * CRITMACH_W * CDC / CRITMACH
            CDC_ALF = -MEXP * CRITMACH_ALF * CDC / CRITMACH
        ENDIF
        !      write(*,*) 'critmach,mach ',critmach,mach
        !      write(*,*) 'cdc,cdc_w,cdc_alf ',cdc,cdc_w,cdc_alf
        !
        FAC = 1.0
        FAC_W = 0.0
        !--- Although test data does not show profile drag increases due to Mach #
        !    you could use something like this to add increase drag by Prandtl-Glauert
        !    (or any function you choose)
        !c      FAC   = PG
        !c      FAC_W = PG_W
        !--- Total drag terms
        CDRAG = FAC * CDRAG + DCD + CDC
        CD_ALF = FAC * CD_ALF + DCD_ALF + CDC_ALF
        CD_W = FAC * CD_W + FAC_W * CDRAG + DCD_W + CDC_ALF
        CD_REY = FAC * CD_REY
        !
    END SUBROUTINE CLCDCM
    !*==CHKLIM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CLCDCM
    
    
    
    SUBROUTINE CHKLIM(N, NSTRT, NEND, F, FMAX)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: FMAX
        INTEGER :: N, NEND, NSTRT
        REAL, DIMENSION(N) :: F
        !
        ! Local variables
        !
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !--- Get starting and end index for array values F(i) < FMAX
        NSTRT = 1
        NEND = N
        !--- Look for first point where F(i)<FMAX
        DO I = 1, N
            IF (F(I)<FMAX) EXIT
        ENDDO
        NSTRT = MAX(I - 1, 1)
        !--- Look for last point where F(i)<FMAX
        DO I = N, 1, -1
            IF (F(I)<FMAX) EXIT
        ENDDO
        NEND = MIN(I + 1, N)
        !
    END SUBROUTINE CHKLIM
    !*==OPFILE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE OPFILE(LU, FNAME)
        use m_userio, only: askc, asks
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        CHARACTER(*) :: FNAME
        INTEGER :: LU
        !
        ! Local variables
        !
        CHARACTER(1) :: ANS, DUMMY
        CHARACTER(4) :: COMAND
        CHARACTER(128) :: COMARG, TMP
        INTEGER :: K, NF
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !
        !---- get filename if it hasn't been already specified
        IF (FNAME==' ') CALL ASKS('Enter output filename^', FNAME)
        !
        !---- try to open file
        OPEN (LU, FILE = FNAME, STATUS = 'OLD', ERR = 50)
        !
        !---- file exists... ask how to proceed
        NF = INDEX(FNAME, ' ') - 1
        TMP = 'File  ' // FNAME(1:NF)                                       &
                & // '  exists.  Overwrite / Append / New file ?^'
        CALL ASKC(TMP, COMAND, COMARG)
        ANS = COMAND(1:1)
        !
        !---- ask again if reply is invalid
        IF (INDEX('OoAaNn', ANS)==0) THEN
            CALL ASKC(' O / A / N  ?^', COMAND, COMARG)
            ANS = COMAND(1:1)
            !
            IF (INDEX('OoAaNn', ANS)==0) THEN
                !------- Still bad reply. Give up asking and just return
                WRITE (*, *) 'No action taken'
                RETURN
            ENDIF
        ENDIF
        !
        !---- at this point, file is open and reply is valid
        IF (INDEX('Oo', ANS)/=0) THEN
            !------ go to beginning of file to overwrite
            REWIND (LU)
            GOTO 60
        ELSEIF (INDEX('Aa', ANS)/=0) THEN
            !------ go to end of file to append
            DO K = 1, 12345678
                READ (LU, 1000, END = 40) DUMMY
                1000       FORMAT (A)
            ENDDO
            40      BACKSPACE (LU)
            GOTO 60
        ELSE
            !------ new file... get filename from command argument, or ask if not supplied
            FNAME = COMARG
            IF (FNAME(1:1)==' ') CALL ASKS('Enter output filename^', &
                    & FNAME)
        ENDIF
        !
        !---- at this point, file FNAME is new or is to be overwritten
        50   OPEN (LU, FILE = FNAME, STATUS = 'UNKNOWN', ERR = 90)
        REWIND (LU)
        !
        60   RETURN
        !
        90   WRITE (*, *) 'Bad filename.'
    END SUBROUTINE OPFILE
    !*==GETCLFACTOR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! OPFILE
    
    
    SUBROUTINE GETCLFACTOR(SIGMA, STAGGER, CLFACTOR)
        use m_spline, only: sevlin
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        REAL, PARAMETER :: PI = 3.1415926535897932384, DTR = PI / 180.
        !
        ! Dummy arguments
        !
        REAL :: CLFACTOR, SIGMA, STAGGER
        !
        ! Local variables
        !
        REAL, DIMENSION(11), SAVE :: A0, A1, A2, X
        REAL :: AA0, AA1, AA2, DAA0, DAA1, DAA2, SIGI, STAGR
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------
        !     Calculates multi-plane cascade effect on lift slope as a
        !     function of solidity and stagger angle
        !
        !     Input:  SIGMA      solidity = Bc/(2*pi*r)
        !             STAGGER    stagger angle (from axis to chordline, in rads)
        !
        !     Output:  CLFACTOR  CLmultiplane/CL2D factor
        !
        !     Implements table-driven quadratic fit to Figure 6-29 in
        !     Wallis, Axial Flow Fans and Ducts.
        !------------------------------------------------------------
        !
        !---- Table of quadratic fit coefficients
        DATA X/0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, &
                & 1.5/
        DATA A0/0.4755, 0.5255, 0.5722, 0.6142, 0.6647, 0.7016, &
                & 0.7643, 0.8302, 0.8932, 0.9366, 0.9814/
        DATA A1/ - 0.367495, -0.341941, -0.300058, -0.255883, &
                & -0.200593, -0.114993, -0.118602, -0.130921, -0.133442, &
                & -0.077980, -0.123071/
        DATA A2/0.489466, 0.477648, 0.453027, 0.430048, 0.381462, &
                & 0.310028, 0.298309, 0.285309, 0.263084, 0.184165, &
                & 0.251594/
        !
        CLFACTOR = 1.0
        IF (SIGMA<=0.6) RETURN
        !
        !---- Interpolate quadratic fit coefficients by 1/solidity
        SIGI = 1.0 / SIGMA
        CALL SEVLIN(SIGI, A0, X, 11, AA0, DAA0)
        CALL SEVLIN(SIGI, A1, X, 11, AA1, DAA1)
        CALL SEVLIN(SIGI, A2, X, 11, AA2, DAA2)
        !
        !---- Only valid for stagger 20deg to 90deg,
        !     Limit low stagger to 20deg value to give constant lift ratio below that
        STAGR = STAGGER
        IF (STAGR<20.0 * DTR) STAGR = 20.0 * DTR
        IF (STAGR>90.0 * DTR) STAGR = 90.0 * DTR
        !
        !---- Quadratic fit for CLFACTOR at this SIGMA as function of STAGGER
        CLFACTOR = AA0 + AA1 * STAGR + AA2 * STAGR * STAGR
        !---- maximum value of lift ratio should be limited to 1.0
        CLFACTOR = MIN(1.0, CLFACTOR)
        !
    END SUBROUTINE GETCLFACTOR
end module m_aero
