module m_solve
    implicit none
contains
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


    subroutine gamsolv
        use i_dfdc
        use m_inigrd, only : setgrdflw
        use m_gauss, only : ludcmp
        use m_geom, only : cvpgen
        use m_rotoper, only : vmavginit, setrotorsrc, stgfind, setdrgobjsrc
        use m_system, only : qcsum, gsolve, gsys, sysp
        use m_qaic, only : qaic
        use m_vels, only : qcpfor
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------------
        !     Generates inviscid panel solution with current flow condition
        !------------------------------------------------------------------
        !
        if (ldbg) then
            write (*, *) 'Entering GAMSOLV'
            write (*, *) ' LNCVP ', lncvp
            write (*, *) ' LQAIC ', lqaic
            write (*, *) ' LSYSP ', lsysp
            write (*, *) ' LGSYS ', lgsys
            write (*, *) ' LGAMA ', lgama
            write (*, *) ' LVMAV ', lvmav
        endif
        !
        !---- set control points and assign their pointers
        if (.not.lncvp) call cvpgen
        !
        !---- set AIC matrix for velocities at control points
        if (.not.lqaic) call qaic(.false.)
        !
        !---- set up pointers for direct problem
        if (.not.lsysp) call sysp
        !
        !---- set up and factor system for direct problem
        if (.not.lgsys) then
            call gsys(.true., .true., .true., 1, nptot, 1, 0)
            if (ldbg) write (*, *) 'Factoring system matrix...'
            call ludcmp(nsysx, nsys, sys(1, 1))
        endif
        !
        !---- Set rotor drag source strengths from current velocities
        call setrotorsrc
        !---- Set (known) drag object source strengths from current velocities
        call setdrgobjsrc
        !
        !---- Initialize VMAVG for wake gamma relaxatio
        if (.not.lvmav) call vmavginit(vavginit)
        !
        !---- Solve for current RHS knowns (Qinf,sources,GTH)
        if (.not.lgama) call gsolve
        !
        !---- Calculate velocities at control points
        call qcsum
        !---- Calculate Cp's on surfaces, forces
        call qcpfor
        !---- Put flow data into wake grid
        call setgrdflw
        !---- Find stagnation points
        call stgfind
        !
        !c      write(*,*) 'Qn =', (QNDOF(IEL), IEL=1, 2)

    end subroutine gamsolv
    !*==GAMSOL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GAMSOLV



    subroutine gamsol
        use i_dfdc
        use m_inigrd, only : setgrdflw
        use m_gauss, only : ludcmp
        use m_geom, only : cvpgen
        use m_rotoper, only : vmavginit, setrotorsrc, stgfind, setdrgobjsrc
        use m_system, only : gusum, qcsum, gucalc, gsys, sysp
        use m_qaic, only : qaic
        use m_vels, only : qcpfor
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: iel, ip, ip1, ip2, ir, iu
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------------
        !     Generates inviscid panel solution by superposing unit-strength
        !     solutions
        !--------------------------------------------
        !
        !---- set control points and assign their pointers
        if (.not.lncvp) call cvpgen
        !
        !---- set AIC matrix for velocities at control points
        if (.not.lqaic) call qaic(.false.)
        !
        !---- set up pointers for direct problem
        if (.not.lsysp) call sysp
        !
        !---- set up and factor system for direct problem
        if (.not.lgsys) then
            call gsys(.true., .true., .true., 1, nptot, 1, 0)
            if (ldbg) write (*, *) 'Factoring system matrix...'
            call ludcmp(nsysx, nsys, sys(1, 1))
        endif
        !
        if (.not.lgamu) then
            !
            !----- all r.h.s. vectors are to be considered invalid...
            if (ldbg) write (*, *) 'Clearing unit-strength solutions...'
            nu = 0
            iuqinf = 0
            do ir = 1, nrp
                iuvwk(ir) = 0
            enddo
            do ip = 1, nptot
                iugam(ip) = 0
                iusig(ip) = 0
            enddo
            do iu = 1, nux
                luset(iu) = .false.
            enddo
            lgamu = .true.
            !
        endif
        !
        !---- New unit solution for Qinf
        if (qinf/=0.0) then
            if (ldbg) write (*, *) 'Qinf unit calc'
            call gucalc(.true., .false., 0, -1)
        endif
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
        do iel = 1, nel
            if (gamset(iel)/=0.0 .or. sigset(iel)/=0.0 .or. netype(iel)   &
                    & ==5 .or. netype(iel)==6) then
                ip1 = ipfrst(iel)
                ip2 = iplast(iel)
                if (ldbg) write (*, *) 'Element set unit calc IEL ', iel
                call gucalc(.false., .false., ip1, ip2)
            endif
        enddo
        !
        !---- no imposed viscous strengths
        do ip = 1, nptot
            gamvsp(ip) = 0.
            sigvsp(ip) = 0.
        enddo
        !
        !---- Set rotor drag source strengths from current velocities
        call setrotorsrc
        !---- Set (known) drag object source strengths from current velocities
        call setdrgobjsrc
        !
        !---- Initialize VMAVG for wake gamma relaxatio
        if (.not.lvmav) call vmavginit(vavginit)
        !
        !---- combine unit vorticities weighted by current singularities and freestream
        if (.not.lgama) call gusum
        !ccHHY
        !---- Calculate velocities at control points
        call qcsum
        !ccHHY
        !
        !---- Calculate Cp's on surfaces, forces
        call qcpfor
        !
        !---- Put flow data into wake grid
        call setgrdflw
        !---- Find stagnation points
        call stgfind
        !
        !c      write(*,*) 'Qn =', (QNDOF(IEL), IEL=1, 2)

    end subroutine gamsol
end module m_solve
