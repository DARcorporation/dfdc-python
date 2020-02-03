module m_atmo
    implicit none
contains
    !*==ATMO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
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

    SUBROUTINE ATMO(ALSPEC, DELTAT, VSOALT, RHOALT, RMUALT)
        use m_spline, only : seval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        INTEGER, PARAMETER :: N = 44
        !
        ! Dummy arguments
        !
        REAL :: ALSPEC, DELTAT, RHOALT, RMUALT, VSOALT
        !
        ! Local variables
        !
        REAL :: ALFRAC, DALT, DRHO, DRMU, DVSO
        REAL, DIMENSION(N), SAVE :: ALT, RHO, RMU, VSO
        LOGICAL, SAVE :: FIRST
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Returns speed of sound (VSO) in m/s, density (RHO)
        !     in kg/m^3, and dynamic viscosity (RMU) in kg/m-s
        !     of standard atmosphere at specified altitude ALSPEC
        !     (in kilometers).  If ALSPEC=-1, water properties
        !     at 15 Celsius are returned.
        !
        !     Reference:  "U.S. Standard Atmosphere", NOAA.
        !---------------------------------------------------------
        !
        !
        DATA FIRST/.TRUE./
        DATA ALT/0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, &
                & 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, &
                & 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, &
                & 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, &
                & 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 45.0, &
                & 60.0, 75.0/
        DATA VSO/340.0, 336.0, 332.0, 329.0, 325.0, 320.0, 316.0, &
                & 312.0, 308.0, 304.0, 299.0, 295.0, 295.0, 295.0, &
                & 295.0, 295.0, 295.0, 295.0, 295.0, 295.0, 295.0, &
                & 295.8, 296.4, 297.1, 297.8, 298.5, 299.1, 299.8, &
                & 300.5, 301.1, 301.8, 302.5, 303.1, 305.0, 306.8, &
                & 308.7, 310.5, 312.3, 314.0, 316.0, 318.0, 355.0, &
                & 372.0, 325.0/
        DATA RHO/1.226, 1.112, 1.007, 0.909, 0.820, 0.737, 0.660, &
                & 0.589, 0.526, 0.467, 0.413, 0.364, 0.311, 0.265, &
                & 0.227, 0.194, 0.163, 0.141, 0.121, 0.103, .0880, &
                & .0749, .0637, .0543, .0463, .0395, .0338, .0288, &
                & .0246, .0210, .0180, .0154, .0132, .0113, .0096, &
                & .0082, .0070, .0060, .0052, .0044, 0.004, 0.002, &
                & 3.9E-4, 8.0E-5/
        DATA RMU/1.780, 1.749, 1.717, 1.684, 1.652, 1.619, 1.586, &
                & 1.552, 1.517, 1.482, 1.447, 1.418, 1.418, 1.418, &
                & 1.418, 1.418, 1.418, 1.418, 1.418, 1.418, 1.418, &
                & 1.427, 1.433, 1.438, 1.444, 1.449, 1.454, 1.460, &
                & 1.465, 1.471, 1.476, 1.481, 1.487, 1.502, 1.512, &
                & 1.532, 1.546, 1.561, 1.580, 1.600, 1.700, 1.912, &
                & 2.047, 1.667/
        !
        !---- special case: Water at STP
        IF (ALSPEC==-1.0) THEN
            VSOALT = 1500.
            RHOALT = 1000.
            RMUALT = 1.15E-3
            WRITE (*, *) '                              o        '
            WRITE (*, *) 'ATMO: You are underwater at 15  Celsius'
            RETURN
        ENDIF
        !
        !---- linearly interpolate quantities from tabulated values
        DO I = 2, N
            IF (ALSPEC<=ALT(I)) THEN
                !
                DALT = ALT(I) - ALT(I - 1)
                DVSO = VSO(I) - VSO(I - 1)
                DRHO = RHO(I) - RHO(I - 1)
                DRMU = RMU(I) - RMU(I - 1)
                !
                ALFRAC = (ALSPEC - ALT(I - 1)) / DALT
                !
                VSOALT = VSO(I - 1) + DVSO * ALFRAC
                RHOALT = RHO(I - 1) + DRHO * ALFRAC
                RMUALT = RMU(I - 1) + DRMU * ALFRAC
                RMUALT = RMUALT * 1.0E-5
                !
                RETURN
            ENDIF
        ENDDO
        !
        !
        IF (ALSPEC>ALT(N)) THEN
            WRITE (*, *) ' '
            WRITE (*, *) 'ATMO: You''re in low earth orbit.  Good luck.'
            VSOALT = VSO(N)
            RHOALT = RHO(N)
            RMUALT = RMU(N) * 1.0E-5
            RETURN
        ENDIF
        !
        !      IF(FIRST) THEN
        !       DO 20 I=1, N
        !         RHO(I) = ALOG(RHO(I))
        ! 20    CONTINUE
        !       CALL SPLINE(VSO,VSOH,ALT,N)
        !       CALL SPLIND(RHO,RHOH,ALT,N,999.0,0.0)
        !       CALL SPLINE(RMU,RMUH,ALT,N)
        !       FIRST = .FALSE.
        !      ENDIF
        !C
        !C---- interpolate quantities from splines
        !      VSOALT = SEVAL(ALSPEC,VSO,VSOH,ALT,N)
        !      RHOALT = SEVAL(ALSPEC,RHO,RHOH,ALT,N)
        !      RMUALT = SEVAL(ALSPEC,RMU,RMUH,ALT,N) * 1.0E-5
        !      RHOALT = EXP(RHOALT)
        !C
    END SUBROUTINE ATMO
    !*==FLOSHO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ATMO


    SUBROUTINE FLOSHO(LU, VSO, RHO, RMU)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: LU
        REAL :: RHO, RMU, VSO
        !
        !*** End of declarations rewritten by SPAG
        !
        WRITE (LU, 10) VSO, RHO, RMU
        10   FORMAT (/' Speed of sound (m/s)  :', &
                &F10.3/' Air density   (kg/m^3):', &
                &F10.5/' Air viscosity (kg/m-s):', E11.4)
    END SUBROUTINE FLOSHO
end module m_atmo
