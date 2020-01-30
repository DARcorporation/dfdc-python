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


SUBROUTINE CPCALC(N, Q, QINF, QREF, CP, CP_Q, CP_QINF)
    !-------------------------------------------------------------
    !     Sets Cp from speed.
    !-------------------------------------------------------------
    DIMENSION Q(2, *), CP(*), CP_Q(2, *), CP_QINF(*)
    !
    IF(QREF.NE.0.0) THEN
        QSQREF = QREF**2
    ELSE
        QSQREF = 1.0
    ENDIF
    !
    DO I = 1, N
        QSQ = Q(1, I)**2 + Q(2, I)**2
        CP(I) = (QINF**2 - QSQ) / QSQREF
        CP_Q(1, I) = -2.0 * Q(1, I) / QSQREF
        CP_Q(2, I) = -2.0 * Q(2, I) / QSQREF
        CP_QINF(I) = 2.0 * QINF / QSQREF
    ENDDO
    !
    RETURN
END
! CPCALC



SUBROUTINE CPADDHS(XX, YY, CP, DHH, DSS, VTT)
    !-------------------------------------------------------
    !     Add enthalpy from rotor to Cp at XX,YY if point
    !     falls within the slipstream grid
    !     CP is modified and DH,DS and Vtheta are returned
    !-------------------------------------------------------
    INCLUDE 'DFDC.INC'
    !
    DATA EPS /1.0E-3/
    !
    VTT = 0.0
    DHH = 0.0
    DSS = 0.0
    !
    !---- Check for point upstream of rotor
    IF(XX.LT.XGMIN) RETURN
    !
    !---- check for point in rotor wake grid
    CALL XYGRDFIND(XX, YY, IX, XG, YG, II, JJ, IC, JC)
    !
    IF(IC.NE.0 .AND. JC.NE.0) THEN
        !c        write(*,99) 'cpaddhs grid point i,j for x,y ',ic,jc,xx,yy
        QRSQ = QREF**2
        CIRC = BGAMG(IC, JC)
        VTT = BGAMG(IC, JC) * PI2I / YY
        VTSQ = VTT**2
        DHH = DHG(IC, JC)
        DSS = DSG(IC, JC)
        CP = CP + 2.0 * (DHH - DSS) / QRSQ - VTSQ / QRSQ
        !c      ELSE
        !c        write(*,99) 'cpaddhs no grid point i,j for x,y ',ic,jc,xx,yy
    ENDIF
    99   format(A, 2I5, 5(1X, F12.6))
    !
    RETURN
END


SUBROUTINE FCOEFF(NC, CPL, CPR, XC, YC, ANC, DSC, &
        CX, CY, CM)
    !
    !---- Integrate inviscid panel forces from Cp's at control points
    !     using the difference of pressures on both sides of sheet
    !
    DIMENSION CPL(*), CPR(*), XC(*), YC(*), ANC(2, *), DSC(*)
    !
    PARAMETER (PI = 3.1415926535897932384)
    !
    CX = 0.
    CY = 0.
    CM = 0.
    DO IC = 1, NC
        TWOPIR = 2.0 * PI * YC(IC)
        CX = CX + (CPL(IC) - CPR(IC)) * ANC(1, IC) * DSC(IC) * TWOPIR
        !---- only axisymmetric forces on axisymetric geometry !!
        !cc     CY = CY + (CPL(IC) - CPR(IC))*ANC(2,IC)*DSC(IC)*TWOPIR
        !cc     CM = CM + (CPL(IC) - CPR(IC))
        !cc  &           *(  XC(IC)*ANC(2,IC)
        !cc  &             - YC(IC)*ANC(1,IC) )*DSC(IC)*TWOPIR
    ENDDO
    !
END
! FCOEFF
