!*==GAMSOLV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! SGRENUM

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


SUBROUTINE GAMSOLV
    USE I_DFDC
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    !*** End of declarations rewritten by SPAG
    !
    !------------------------------------------------------------------
    !     Generates inviscid panel solution with current flow condition
    !------------------------------------------------------------------
    !
    IF (LDBG) THEN
        WRITE (*, *) 'Entering GAMSOLV'
        WRITE (*, *) ' LNCVP ', LNCVP
        WRITE (*, *) ' LQAIC ', LQAIC
        WRITE (*, *) ' LSYSP ', LSYSP
        WRITE (*, *) ' LGSYS ', LGSYS
        WRITE (*, *) ' LGAMA ', LGAMA
        WRITE (*, *) ' LVMAV ', LVMAV
    ENDIF
    !
    !---- set control points and assign their pointers
    IF (.NOT.LNCVP) CALL CVPGEN
    !
    !---- set AIC matrix for velocities at control points
    IF (.NOT.LQAIC) CALL QAIC(.FALSE.)
    !
    !---- set up pointers for direct problem
    IF (.NOT.LSYSP) CALL SYSP
    !
    !---- set up and factor system for direct problem
    IF (.NOT.LGSYS) THEN
        CALL GSYS(.TRUE., .TRUE., .TRUE., 1, NPTOT, 1, 0)
        IF (LDBG) WRITE (*, *) 'Factoring system matrix...'
        CALL LUDCMP(NSYSX, NSYS, SYS(1, 1))
    ENDIF
    !
    !---- Set rotor drag source strengths from current velocities
    CALL SETROTORSRC
    !---- Set (known) drag object source strengths from current velocities
    CALL SETDRGOBJSRC
    !
    !---- Initialize VMAVG for wake gamma relaxatio
    IF (.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
    !
    !---- Solve for current RHS knowns (Qinf,sources,GTH)
    IF (.NOT.LGAMA) CALL GSOLVE
    !
    !---- Calculate velocities at control points
    CALL QCSUM
    !---- Calculate Cp's on surfaces, forces
    CALL QCPFOR
    !---- Put flow data into wake grid
    CALL SETGRDFLW
    !---- Find stagnation points
    CALL STGFIND
    !
    !c      write(*,*) 'Qn =', (QNDOF(IEL), IEL=1, 2)

END SUBROUTINE GAMSOLV
!*==GAMSOL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! GAMSOLV



SUBROUTINE GAMSOL
    USE I_DFDC
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Local variables
    !
    INTEGER :: IEL, IP, IP1, IP2, IR, IU
    !
    !*** End of declarations rewritten by SPAG
    !
    !------------------------------------------------------------------
    !     Generates inviscid panel solution by superposing unit-strength
    !     solutions
    !--------------------------------------------
    !
    !---- set control points and assign their pointers
    IF (.NOT.LNCVP) CALL CVPGEN
    !
    !---- set AIC matrix for velocities at control points
    IF (.NOT.LQAIC) CALL QAIC(.FALSE.)
    !
    !---- set up pointers for direct problem
    IF (.NOT.LSYSP) CALL SYSP
    !
    !---- set up and factor system for direct problem
    IF (.NOT.LGSYS) THEN
        CALL GSYS(.TRUE., .TRUE., .TRUE., 1, NPTOT, 1, 0)
        IF (LDBG) WRITE (*, *) 'Factoring system matrix...'
        CALL LUDCMP(NSYSX, NSYS, SYS(1, 1))
    ENDIF
    !
    IF (.NOT.LGAMU) THEN
        !
        !----- all r.h.s. vectors are to be considered invalid...
        IF (LDBG) WRITE (*, *) 'Clearing unit-strength solutions...'
        NU = 0
        IUQINF = 0
        DO IR = 1, NRP
            IUVWK(IR) = 0
        ENDDO
        DO IP = 1, NPTOT
            IUGAM(IP) = 0
            IUSIG(IP) = 0
        ENDDO
        DO IU = 1, NUX
            LUSET(IU) = .FALSE.
        ENDDO
        LGAMU = .TRUE.
        !
    ENDIF
    !
    !---- New unit solution for Qinf
    IF (QINF/=0.0) THEN
        IF (LDBG) WRITE (*, *) 'Qinf unit calc'
        CALL GUCALC(.TRUE., .FALSE., 0, -1)
    ENDIF
    !
    !---- New unit solutions for B*Circ rotor circulation distribution
    !      LBGAM = .FALSE.
    !      DO IR = 1, NRP
    !        LBGAM = (LBGAM .OR. WBGAM(IR).NE.0.0)
    !      END DO
    !      IF(LBGAM) THEN
    !        IF(LDBG) WRITE(*,*) 'Rotor circulation unit calc'
    !        CALL GUCALC(.FALSE.,.TRUE., 0,-1)
    !      ENDIF
    !
    !---- New unit solutions for RHS singularities
    DO IEL = 1, NEL
        IF (GAMSET(IEL)/=0.0 .OR. SIGSET(IEL)/=0.0 .OR. NETYPE(IEL)   &
                & ==5 .OR. NETYPE(IEL)==6) THEN
            IP1 = IPFRST(IEL)
            IP2 = IPLAST(IEL)
            IF (LDBG) WRITE (*, *) 'Element set unit calc IEL ', IEL
            CALL GUCALC(.FALSE., .FALSE., IP1, IP2)
        ENDIF
    ENDDO
    !
    !---- no imposed viscous strengths
    DO IP = 1, NPTOT
        GAMVSP(IP) = 0.
        SIGVSP(IP) = 0.
    ENDDO
    !
    !---- Set rotor drag source strengths from current velocities
    CALL SETROTORSRC
    !---- Set (known) drag object source strengths from current velocities
    CALL SETDRGOBJSRC
    !
    !---- Initialize VMAVG for wake gamma relaxatio
    IF (.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
    !
    !---- combine unit vorticities weighted by current singularities and freestream
    IF (.NOT.LGAMA) CALL GUSUM
    !ccHHY
    !---- Calculate velocities at control points
    CALL QCSUM
    !ccHHY
    !
    !---- Calculate Cp's on surfaces, forces
    CALL QCPFOR
    !
    !---- Put flow data into wake grid
    CALL SETGRDFLW
    !---- Find stagnation points
    CALL STGFIND
    !
    !c      write(*,*) 'Qn =', (QNDOF(IEL), IEL=1, 2)

END SUBROUTINE GAMSOL
