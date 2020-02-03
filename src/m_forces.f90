module m_forces
    implicit none
contains
    !*==CPCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! DFSAVE


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


    subroutine cpcalc(n, q, qinf, qref, cp, cp_q, cp_qinf)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: qinf, qref
        real, dimension(*) :: cp, cp_qinf
        real, dimension(2, *) :: cp_q, q
        !
        ! Local variables
        !
        integer :: i
        real :: qsq, qsqref
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Sets Cp from speed.
        !-------------------------------------------------------------
        !
        if (qref/=0.0) then
            qsqref = qref**2
        else
            qsqref = 1.0
        endif
        !
        do i = 1, n
            qsq = q(1, i)**2 + q(2, i)**2
            cp(i) = (qinf**2 - qsq) / qsqref
            cp_q(1, i) = -2.0 * q(1, i) / qsqref
            cp_q(2, i) = -2.0 * q(2, i) / qsqref
            cp_qinf(i) = 2.0 * qinf / qsqref
        enddo
        !
    end subroutine cpcalc
    !*==CPADDHS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CPCALC



    subroutine cpaddhs(xx, yy, cp, dhh, dss, vtt)
        use i_dfdc
        use m_grdutils, only : xygrdfind
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: cp, dhh, dss, vtt, xx, yy
        !
        ! Local variables
        !
        real :: circ, qrsq, vtsq
        real, save :: eps
        integer :: ic, jc
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Add enthalpy from rotor to Cp at XX,YY if point
        !     falls within the slipstream grid
        !     CP is modified and DH,DS and Vtheta are returned
        !-------------------------------------------------------
        !
        data eps/1.0e-3/
        !
        vtt = 0.0
        dhh = 0.0
        dss = 0.0
        !
        !---- Check for point upstream of rotor
        if (xx<xgmin) return
        !
        !---- check for point in rotor wake grid
        call xygrdfind(xx, yy, ix, xg, yg, ii, jj, ic, jc)
        !
        if (ic/=0 .and. jc/=0) then
            !c        write(*,99) 'cpaddhs grid point i,j for x,y ',ic,jc,xx,yy
            qrsq = qref**2
            circ = bgamg(ic, jc)
            vtt = bgamg(ic, jc) * pi2i / yy
            vtsq = vtt**2
            dhh = dhg(ic, jc)
            dss = dsg(ic, jc)
            cp = cp + 2.0 * (dhh - dss) / qrsq - vtsq / qrsq
            !c      ELSE
            !c        write(*,99) 'cpaddhs no grid point i,j for x,y ',ic,jc,xx,yy
        endif
        99   format (a, 2i5, 5(1x, f12.6))
        !
    end subroutine cpaddhs
    !*==FCOEFF.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine fcoeff(nc, cpl, cpr, xc, yc, anc, dsc, cx, cy, cm)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        real, parameter :: pi = 3.1415926535897932384
        !
        ! Dummy arguments
        !
        real :: cm, cx, cy
        integer :: nc
        real, dimension(2, *) :: anc
        real, dimension(*) :: cpl, cpr, dsc, xc, yc
        !
        ! Local variables
        !
        integer :: ic
        real :: twopir
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- Integrate inviscid panel forces from Cp's at control points
        !     using the difference of pressures on both sides of sheet
        !
        !
        !
        cx = 0.
        cy = 0.
        cm = 0.
        do ic = 1, nc
            twopir = 2.0 * pi * yc(ic)
            cx = cx + (cpl(ic) - cpr(ic)) * anc(1, ic) * dsc(ic) * twopir
            !---- only axisymmetric forces on axisymetric geometry !!
            !cc     CY = CY + (CPL(IC) - CPR(IC))*ANC(2,IC)*DSC(IC)*TWOPIR
            !cc     CM = CM + (CPL(IC) - CPR(IC))
            !cc  &           *(  XC(IC)*ANC(2,IC)
            !cc  &             - YC(IC)*ANC(1,IC) )*DSC(IC)*TWOPIR
        enddo
        !
    end subroutine fcoeff
end module m_forces
