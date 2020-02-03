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

    subroutine atmo(alspec, deltat, vsoalt, rhoalt, rmualt)
        use m_spline, only : seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: n = 44
        !
        ! Dummy arguments
        !
        real :: alspec, deltat, rhoalt, rmualt, vsoalt
        !
        ! Local variables
        !
        real :: alfrac, dalt, drho, drmu, dvso
        real, dimension(n), save :: alt, rho, rmu, vso
        logical, save :: first
        integer :: i
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
        data first/.true./
        data alt/0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, &
                & 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, &
                & 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, &
                & 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, &
                & 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 45.0, &
                & 60.0, 75.0/
        data vso/340.0, 336.0, 332.0, 329.0, 325.0, 320.0, 316.0, &
                & 312.0, 308.0, 304.0, 299.0, 295.0, 295.0, 295.0, &
                & 295.0, 295.0, 295.0, 295.0, 295.0, 295.0, 295.0, &
                & 295.8, 296.4, 297.1, 297.8, 298.5, 299.1, 299.8, &
                & 300.5, 301.1, 301.8, 302.5, 303.1, 305.0, 306.8, &
                & 308.7, 310.5, 312.3, 314.0, 316.0, 318.0, 355.0, &
                & 372.0, 325.0/
        data rho/1.226, 1.112, 1.007, 0.909, 0.820, 0.737, 0.660, &
                & 0.589, 0.526, 0.467, 0.413, 0.364, 0.311, 0.265, &
                & 0.227, 0.194, 0.163, 0.141, 0.121, 0.103, .0880, &
                & .0749, .0637, .0543, .0463, .0395, .0338, .0288, &
                & .0246, .0210, .0180, .0154, .0132, .0113, .0096, &
                & .0082, .0070, .0060, .0052, .0044, 0.004, 0.002, &
                & 3.9e-4, 8.0e-5/
        data rmu/1.780, 1.749, 1.717, 1.684, 1.652, 1.619, 1.586, &
                & 1.552, 1.517, 1.482, 1.447, 1.418, 1.418, 1.418, &
                & 1.418, 1.418, 1.418, 1.418, 1.418, 1.418, 1.418, &
                & 1.427, 1.433, 1.438, 1.444, 1.449, 1.454, 1.460, &
                & 1.465, 1.471, 1.476, 1.481, 1.487, 1.502, 1.512, &
                & 1.532, 1.546, 1.561, 1.580, 1.600, 1.700, 1.912, &
                & 2.047, 1.667/
        !
        !---- special case: Water at STP
        if (alspec==-1.0) then
            vsoalt = 1500.
            rhoalt = 1000.
            rmualt = 1.15e-3
            write (*, *) '                              o        '
            write (*, *) 'ATMO: You are underwater at 15  Celsius'
            return
        endif
        !
        !---- linearly interpolate quantities from tabulated values
        do i = 2, n
            if (alspec<=alt(i)) then
                !
                dalt = alt(i) - alt(i - 1)
                dvso = vso(i) - vso(i - 1)
                drho = rho(i) - rho(i - 1)
                drmu = rmu(i) - rmu(i - 1)
                !
                alfrac = (alspec - alt(i - 1)) / dalt
                !
                vsoalt = vso(i - 1) + dvso * alfrac
                rhoalt = rho(i - 1) + drho * alfrac
                rmualt = rmu(i - 1) + drmu * alfrac
                rmualt = rmualt * 1.0e-5
                !
                return
            endif
        enddo
        !
        !
        if (alspec>alt(n)) then
            write (*, *) ' '
            write (*, *) 'ATMO: You''re in low earth orbit.  Good luck.'
            vsoalt = vso(n)
            rhoalt = rho(n)
            rmualt = rmu(n) * 1.0e-5
            return
        endif
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
    end subroutine atmo
    !*==FLOSHO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ATMO


    subroutine flosho(lu, vso, rho, rmu)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: lu
        real :: rho, rmu, vso
        !
        !*** End of declarations rewritten by SPAG
        !
        write (lu, 10) vso, rho, rmu
        10   format (/' Speed of sound (m/s)  :', &
                &f10.3/' Air density   (kg/m^3):', &
                &f10.5/' Air viscosity (kg/m-s):', e11.4)
    end subroutine flosho
end module m_atmo
