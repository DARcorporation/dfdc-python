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


      SUBROUTINE CPCALC(N, Q,QINF,QREF, CP,CP_Q,CP_QINF )
C-------------------------------------------------------------
C     Sets Cp from speed.
C-------------------------------------------------------------
      DIMENSION Q(2,*), CP(*), CP_Q(2,*), CP_QINF(*)
C
      IF(QREF.NE.0.0) THEN
       QSQREF = QREF**2
      ELSE
       QSQREF = 1.0
      ENDIF
C
      DO I = 1, N
        QSQ = Q(1,I)**2 + Q(2,I)**2
        CP(I) = (QINF**2 - QSQ)/QSQREF
        CP_Q(1,I) = -2.0*Q(1,I)/QSQREF
        CP_Q(2,I) = -2.0*Q(2,I)/QSQREF
        CP_QINF(I) = 2.0*QINF  /QSQREF
      ENDDO
C
      RETURN
      END ! CPCALC



      SUBROUTINE CPADDHS(XX,YY,CP,DHH,DSS,VTT)
C-------------------------------------------------------
C     Add enthalpy from rotor to Cp at XX,YY if point 
C     falls within the slipstream grid
C     CP is modified and DH,DS and Vtheta are returned
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DATA EPS /1.0E-3/
C
      VTT = 0.0
      DHH = 0.0
      DSS = 0.0
C
C---- Check for point upstream of rotor
      IF(XX.LT.XGMIN) RETURN
C
C---- check for point in rotor wake grid
      CALL XYGRDFIND(XX,YY,IX,XG,YG,II,JJ,IC,JC)
C
      IF(IC.NE.0 .AND. JC.NE.0) THEN
cc        write(*,99) 'cpaddhs grid point i,j for x,y ',ic,jc,xx,yy
        QRSQ = QREF**2
        CIRC = BGAMG(IC,JC)
        VTT  = BGAMG(IC,JC)*PI2I/YY
        VTSQ = VTT**2
        DHH = DHG(IC,JC)
        DSS = DSG(IC,JC)
        CP = CP + 2.0*(DHH-DSS)/QRSQ - VTSQ/QRSQ
cc      ELSE
cc        write(*,99) 'cpaddhs no grid point i,j for x,y ',ic,jc,xx,yy
      ENDIF
 99   format(A,2I5,5(1X,F12.6))
C
      RETURN
      END


      SUBROUTINE FCOEFF(NC, CPL,CPR, XC,YC,ANC,DSC,
     &                  CX,CY,CM)
C
C---- Integrate inviscid panel forces from Cp's at control points 
C     using the difference of pressures on both sides of sheet
C
      DIMENSION CPL(*), CPR(*), XC(*), YC(*), ANC(2,*), DSC(*)
C
      PARAMETER (PI=3.1415926535897932384)
C
      CX = 0.
      CY = 0.
      CM = 0.
      DO IC = 1, NC
        TWOPIR = 2.0*PI*YC(IC)
        CX = CX + (CPL(IC) - CPR(IC))*ANC(1,IC)*DSC(IC)*TWOPIR
C---- only axisymmetric forces on axisymetric geometry !!
ccc     CY = CY + (CPL(IC) - CPR(IC))*ANC(2,IC)*DSC(IC)*TWOPIR
ccc     CM = CM + (CPL(IC) - CPR(IC))
ccc  &           *(  XC(IC)*ANC(2,IC)
ccc  &             - YC(IC)*ANC(1,IC) )*DSC(IC)*TWOPIR
      ENDDO
C
      END ! FCOEFF

