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

SUBROUTINE XSIETA(XM, X1, X2, X3, X4, YM, Y1, Y2, Y3, Y4, XSI, ETA, &
        XSI_X, XSI_Y, ETA_X, ETA_Y)
    !............................................................
    !     Subroutine for calculating the local bi-linear
    !     coordinates of a marker inside a quadrilateral cell.
    !     Also calculates the local coordinate derivatives.
    !
    !   Input:
    !   ------
    !     XM,YM          marker coordinates
    !     X1,X2,X3,X4    cell vertex x-coordinates
    !     Y1,Y2,Y3,Y4    cell vertex y-coordinates
    !
    !   Output:
    !   -------
    !     XSI,ETA        local cell coordinates at marker
    !     XSI_X          dXSI/dX  at marker
    !     XSI_Y          dXSI/dY  at marker
    !     ETA_X          dETA/dX  at marker
    !     ETA_Y          dETA/dY  at marker
    !
    !............................................................
    !
    !---- initial guess for xsi, eta  (center of cell)
    XS = 0.
    ET = 0.
    !
    !---- perform 12 or less Newton iterations for improving the guess
    DO 10 IXY = 1, 12
        Z1 = (1.0 - XS) * (1.0 - ET)
        Z2 = (1.0 + XS) * (1.0 - ET)
        Z3 = (1.0 + XS) * (1.0 + ET)
        Z4 = (1.0 - XS) * (1.0 + ET)
        DR1 = -(Z1 * X1 + Z2 * X2 + Z3 * X3 + Z4 * X4 - 4.0 * XM)
        DR2 = -(Z1 * Y1 + Z2 * Y2 + Z3 * Y3 + Z4 * Y4 - 4.0 * YM)
        A11 = -(1.0 - ET) * X1 + (1.0 - ET) * X2 + (1.0 + ET) * X3 - (1.0 + ET) * X4
        A12 = -(1.0 - XS) * X1 - (1.0 + XS) * X2 + (1.0 + XS) * X3 + (1.0 - XS) * X4
        A21 = -(1.0 - ET) * Y1 + (1.0 - ET) * Y2 + (1.0 + ET) * Y3 - (1.0 + ET) * Y4
        A22 = -(1.0 - XS) * Y1 - (1.0 + XS) * Y2 + (1.0 + XS) * Y3 + (1.0 - XS) * Y4
        DETINV = 1.0 / (A11 * A22 - A12 * A21)
        DXS = DETINV * (DR1 * A22 - A12 * DR2)
        DET = DETINV * (A11 * DR2 - DR1 * A21)
        XS = XS + DXS
        ET = ET + DET
        IF(AMAX1(ABS(DXS), ABS(DET)) .LT. 1.0E-4) GO TO 11
    10 CONTINUE
    WRITE(6, *) 'Xsi-Eta conv failed at mrkr', XM, YM, DXS, DET
    !
    11 XSI = XS
    ETA = ET
    !
    !---- derivatives are a free bonus from Newton procedure
    XSI_X = 4.0 * A22 * DETINV
    ETA_X = -4.0 * A21 * DETINV
    XSI_Y = -4.0 * A12 * DETINV
    ETA_Y = 4.0 * A11 * DETINV
    !
    RETURN
END
! XSIETA



FUNCTION FINT(XSI, ETA, Q1, Q2, Q3, Q4)
    !...................................................
    !     Interpolates cell corner values Q1,Q2,Q3,Q4
    !     to marker at local cell coodinates XSI,ETA
    !...................................................
    !
    Z1 = (1.0 - XSI) * (1.0 - ETA)
    Z2 = (1.0 + XSI) * (1.0 - ETA)
    Z3 = (1.0 + XSI) * (1.0 + ETA)
    Z4 = (1.0 - XSI) * (1.0 + ETA)
    !
    FINT = (Q1 * Z1 + Q2 * Z2 + Q3 * Z3 + Q4 * Z4) * 0.25
    !
    RETURN
END
! FINT



SUBROUTINE XYGRDFIND(XF, YF, IX, X, Y, II, JJ, IC, JC)
    DIMENSION X(IX, *), Y(IX, *)
    DIMENSION XQ(5), YQ(5)
    !
    !..............................................
    !     Last-resort search through all
    !     cells to locate point XF,YF in the grid
    !
    !     Returns indices IC,JC (to IC+1,JC+1) of cell
    !     containing point or 0,0 for point outside
    !     of grid.
    !..............................................
    !
    IC = 0
    JC = 0
    !
    !---- sweep over all cells in grid
    DO IO = 1, II - 1
        IP = IO + 1
        DO JO = 1, JJ - 1
            JP = JO + 1
            !
            XQ(1) = X(IO, JO)
            XQ(2) = X(IP, JO)
            XQ(3) = X(IP, JP)
            XQ(4) = X(IO, JP)
            YQ(1) = Y(IO, JO)
            YQ(2) = Y(IP, JO)
            YQ(3) = Y(IP, JP)
            YQ(4) = Y(IO, JP)
            !
            !------ initial check to see if point is outside cell limits
            !       (this disqualifies nearly all cells)
            XMAX = AMAX1(XQ(1), XQ(2), XQ(3), XQ(4))
            XMIN = AMIN1(XQ(1), XQ(2), XQ(3), XQ(4))
            YMAX = AMAX1(YQ(1), YQ(2), YQ(3), YQ(4))
            YMIN = AMIN1(YQ(1), YQ(2), YQ(3), YQ(4))
            !
            IF(XF.GT.XMAX .OR. XF.LT.XMIN .OR.&
                    YF.GT.YMAX .OR. YF.LT.YMIN) GO TO 10
            !
            !------ Inside test for polygon
            NQ = 4
            CALL PNPOLY(XF, YF, XQ, YQ, NQ, INOUT)
            IF(INOUT.GE.0) THEN
                !------ Found point in cell
                IC = IO
                JC = JO
                RETURN
            ENDIF
            !
            10     CONTINUE
        END DO
    END DO
    !
    !c      WRITE(6,*) 'Point grid location failed. x y =',XF,YF
    RETURN
    !
END
! XYGRDFIND


!                                                                       
!     ..................................................................
!                                                                       
!        SUBROUTINE PNPOLY                                              
!                                                                       
!        PURPOSE                                                        
!           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
!                                                                       
!        USAGE                                                          
!           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
!                                                                       
!        DESCRIPTION OF THE PARAMETERS                                  
!           PX      - X-COORDINATE OF POINT IN QUESTION.                
!           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!                     VERTICES OF POLYGON.                              
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!                     VERTICES OF POLYGON.                              
!           N       - NUMBER OF VERTICES IN THE POLYGON.                
!           INOUT   - THE SIGNAL RETURNED:                              
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
!                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
!                                                                       
!        REMARKS                                                        
!           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
!           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
!           OPTIONALLY BE INCREASED BY 1.                               
!           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
!           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
!           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
!           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
!           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
!           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
!
!           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
!           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
!           POINT IS INSIDE OF THE POLYGON.                             
!                                                                       
!     ..................................................................
!                                                                       
SUBROUTINE PNPOLY(PX, PY, XX, YY, N, INOUT)
    PARAMETER (MAXDIM = 200)
    REAL X(MAXDIM), Y(MAXDIM), XX(N), YY(N)
    LOGICAL MX, MY, NX, NY

    !--- Output unit for printed messages
    IF(N.GT.MAXDIM) THEN
        WRITE(*, 1)
        1      FORMAT('WARNING:', I5, ' too many polygon sides in PNPOLY')
        RETURN
    ENDIF
    !
    DO I = 1, N
        X(I) = XX(I) - PX
        Y(I) = YY(I) - PY
    END DO
    !
    INOUT = -1
    DO 2 I = 1, N
        J = 1 + MOD(I, N)
        MX = X(I).GE.0.0
        NX = X(J).GE.0.0
        MY = Y(I).GE.0.0
        NY = Y(J).GE.0.0
        IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX))      GO TO 2
        IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3
        INOUT = -INOUT
        GO TO 2
        !
        3        TEST = (Y(I) * X(J) - X(I) * Y(J)) / (X(J) - X(I))
        IF(TEST.LT.0.0) THEN
            GO TO 2
        ELSEIF(TEST.EQ.0.0) THEN
            INOUT = 0
            RETURN
        ELSEIF(TEST.GT.0.0) THEN
            INOUT = -INOUT
        ENDIF
        !
    2    CONTINUE
    !
    RETURN
END


SUBROUTINE UVGRDC(IX, II, JJ, X, Y, XPOS, YPOS, UC, VC)
    IMPLICIT REAL (A-H, M, O-Z)
    DIMENSION X(IX, *), Y(IX, *)
    DIMENSION XPOS(*), YPOS(*)
    DIMENSION UC(IX, *), VC(IX, *)
    !-------------------------------------------------------------
    !     Calculates velocity on cell centers for axisymmetric
    !     streamfunction grid by finte difference of streamfunction
    !     values.
    !
    !   Input:
    !     IX        first grid-array dimension
    !     II,JJ     grid i,j size
    !     X(i,j)    grid z coordinates
    !     Y(i,j)    grid r coordinates
    !     XPOS(i)   grid s values
    !     YPOS(j)   grid e values  (axisymmetric streamfunction)
    !
    !   Output:
    !     UC(i,j)   u velocity at grid cell centers (II-1 x JJ-1)
    !     VC(i,j)   v velocity at grid cell centers (II-1 x JJ-1)
    !
    !.............................................................
    !
    !   Uses solution to Thompson's grid-generation equations with
    !   e(x,r) satisfying the axisymmetric streamfunction equation.
    !   Hence,  e = constant  lines are streamlines.
    !
    !       s_zz + s_rr = 0
    !       e_zz + e_rr = e_r / r
    !
    !   The equations for u(s,e), v(s,e) are
    !
    !       u = z_e / (r J)
    !
    !       v = r_e / (r J)
    !
    !   where
    !
    !         J  =  z_s r_e  -  z_e r_s
    !
    !-------------------------------------------------------------
    !
    !
    !------ go over all streamline cell centers
    DO 100 JO = 1, JJ - 1
        JP = JO + 1
        !
        !-------- for all cell centers on this streamline
        DO 10 IO = 1, II - 1
            IP = IO + 1
            !
            XOO = X(IO, JO)
            XPO = X(IP, JO)
            XOP = X(IO, JP)
            XPP = X(IP, JP)
            YOO = Y(IO, JO)
            YPO = Y(IP, JO)
            YOP = Y(IO, JP)
            YPP = Y(IP, JP)
            !
            XAV = 0.25 * (XOO + XOP + XPO + XPP)
            YAV = 0.25 * (YOO + YOP + YPO + YPP)
            !
            DXAXI = 0.5 * (-XOO - XOP + XPO + XPP)
            DYAXI = 0.5 * (-YOO - YOP + YPO + YPP)
            DXAET = 0.5 * (-XOO + XOP - XPO + XPP)
            DYAET = 0.5 * (-YOO + YOP - YPO + YPP)
            !
            DXI = XPOS(IP) - XPOS(IO)
            DET = YPOS(JP) - YPOS(JO)
            !
            DXDXI = DXAXI / DXI
            DYDXI = DYAXI / DXI
            DXDET = DXAET / DET
            DYDET = DYAET / DET
            !
            AJA = DYDET * DXDXI - DXDET * DYDXI
            !
            UC(IO, JO) = DXDXI / YAV / AJA
            VC(IO, JO) = DYDXI / YAV / AJA
            !
        10       CONTINUE
        !
    100    CONTINUE
    !
    RETURN
END
! UVGRDC


SUBROUTINE UVGRD(IX, II, JJ, X, Y, XPOS, YPOS, UG, VG)
    IMPLICIT REAL (A-H, M, O-Z)
    DIMENSION X(IX, *), Y(IX, *)
    DIMENSION XPOS(*), YPOS(*)
    DIMENSION UG(IX, *), VG(IX, *)
    !-------------------------------------------------------------
    !     Calculates velocity on axisymmetric streamfunction grid
    !     by finte difference of streamfunction values.
    !
    !   Input:
    !     IX        first grid-array dimension
    !     II,JJ     grid i,j size
    !     X(i,j)    grid z coordinates
    !     Y(i,j)    grid r coordinates
    !     XPOS(i)   grid s values
    !     YPOS(j)   grid e values  (axisymmetric streamfunction)
    !
    !   Output:
    !     UG(i,j)   u velocity at grid points
    !     VG(i,j)   v velocity at grid points
    !
    !.............................................................
    !
    !   Uses solution to Thompson's grid-generation equations with
    !   e(x,r) satisfying the axisymmetric streamfunction equation.
    !   Hence,  e = constant  lines are streamlines.
    !
    !       s_zz + s_rr = 0
    !       e_zz + e_rr = e_r / r
    !
    !   The equations for u(s,e), v(s,e) are
    !
    !       u = z_e / (rho r J)
    !
    !       v = r_e / (rho r J)
    !
    !   where
    !
    !         J  =  z_s r_e  -  z_e r_s
    !
    !-------------------------------------------------------------
    !
    !------ go over all interior streamlines
    DO 100 JO = 2, JJ - 1
        JM = JO - 1
        JP = JO + 1
        !
        !-------- for all interior points on this streamline
        DO 10 IO = 2, II - 1
            IM = IO - 1
            IP = IO + 1
            !
            XMM = X(IM, JM)
            XOM = X(IO, JM)
            XPM = X(IP, JM)
            XMO = X(IM, JO)
            XOO = X(IO, JO)
            XPO = X(IP, JO)
            XMP = X(IM, JP)
            XOP = X(IO, JP)
            XPP = X(IP, JP)
            YMM = Y(IM, JM)
            YOM = Y(IO, JM)
            YPM = Y(IP, JM)
            YMO = Y(IM, JO)
            YOO = Y(IO, JO)
            YPO = Y(IP, JO)
            YMP = Y(IM, JP)
            YOP = Y(IO, JP)
            YPP = Y(IP, JP)
            !
            DXIM = XPOS(IO) - XPOS(IM)
            DXIP = XPOS(IP) - XPOS(IO)
            DXIAV = 0.5 * (DXIM + DXIP)
            !
            DETM = YPOS(JO) - YPOS(JM)
            DETP = YPOS(JP) - YPOS(JO)
            DETAV = 0.5 * (DETM + DETP)
            !
            DXDET = 0.5 * (XOP - XOM) / DETAV
            DYDET = 0.5 * (YOP - YOM) / DETAV
            DXDXI = 0.5 * (XPO - XMO) / DXIAV
            DYDXI = 0.5 * (YPO - YMO) / DXIAV
            !
            AJA = DYDET * DXDXI - DXDET * DYDXI
            YAV = YOO
            !
            UG(IO, JO) = DXDXI / YAV / AJA
            VG(IO, JO) = DYDXI / YAV / AJA
            !
        10       CONTINUE
        !
    100    CONTINUE
    !
    !
    !-------- Calculate velocities on grid boundaries (uses backward differences)
    !
    !-------- Velocities on lower streamline
    JO = 1
    DO IO = 2, II - 1
        IP = IO + 1
        IM = IO - 1
        !
        JP = JO + 1
        JQ = JO + 2
        !
        XOO = X(IO, JO)
        XOP = X(IO, JP)
        XOQ = X(IO, JQ)
        YOO = Y(IO, JO)
        YOP = Y(IO, JP)
        YOQ = Y(IO, JQ)
        !
        DETP = YPOS(JP) - YPOS(JO)
        DETQ = YPOS(JQ) - YPOS(JP)
        !
        DXDE1 = (XOP - XOO) / DETP
        DXDE2 = (XOQ - XOP) / DETQ
        DYDE1 = (YOP - YOO) / DETP
        DYDE2 = (YOQ - YOP) / DETQ
        !--- backwards difference at boundary
        DXDET = DXDE1 - DETP * (DXDE2 - DXDE1) / (DETP + DETQ)
        DYDET = DYDE1 - DETP * (DYDE2 - DYDE1) / (DETP + DETQ)
        !
        XMO = X(IM, JO)
        XOO = X(IO, JO)
        XPO = X(IP, JO)
        YMO = Y(IM, JO)
        YOO = Y(IO, JO)
        YPO = Y(IP, JO)
        !
        DXIM = XPOS(IO) - XPOS(IM)
        DXIP = XPOS(IP) - XPOS(IO)
        DXIAV = 0.5 * (DXIM + DXIP)
        DXDXI = 0.5 * (XPO - XMO) / DXIAV
        DYDXI = 0.5 * (YPO - YMO) / DXIAV
        !
        AJA = DYDET * DXDXI - DXDET * DYDXI
        YAV = YOO
        !
        UG(IO, JO) = DXDXI / YAV / AJA
        VG(IO, JO) = DYDXI / YAV / AJA
    ENDDO
    !
    !-------- Velocities on upper streamline
    JO = JJ
    DO IO = 2, II - 1
        IP = IO + 1
        IM = IO - 1
        !
        JM = JO - 1
        JL = JO - 2
        !
        XOO = X(IO, JO)
        XOM = X(IO, JM)
        XOL = X(IO, JL)
        YOO = Y(IO, JO)
        YOM = Y(IO, JM)
        YOL = Y(IO, JL)
        !
        DETM = YPOS(JO) - YPOS(JM)
        DETL = YPOS(JM) - YPOS(JL)
        DXDE1 = (XOO - XOM) / DETM
        DYDE1 = (YOO - YOM) / DETM
        DXDE2 = (XOM - XOL) / DETL
        DYDE2 = (YOM - YOL) / DETL
        !--- backwards difference at boundary
        DXDET = DXDE1 + DETM * (DXDE1 - DXDE2) / (DETM + DETL)
        DYDET = DYDE1 + DETM * (DYDE1 - DYDE2) / (DETM + DETL)
        !
        XMO = X(IM, JO)
        XOO = X(IO, JO)
        XPO = X(IP, JO)
        YMO = Y(IM, JO)
        YOO = Y(IO, JO)
        YPO = Y(IP, JO)
        !
        DXIM = XPOS(IO) - XPOS(IM)
        DXIP = XPOS(IP) - XPOS(IO)
        DXIAV = 0.5 * (DXIM + DXIP)
        DXDXI = 0.5 * (XPO - XMO) / DXIAV
        DYDXI = 0.5 * (YPO - YMO) / DXIAV
        !
        AJA = DYDET * DXDXI - DXDET * DYDXI
        YAV = YOO
        !
        UG(IO, JO) = DXDXI / YAV / AJA
        VG(IO, JO) = DYDXI / YAV / AJA
    ENDDO
    !
    !-------- Velocities on(i=1) inlet plane (includes corner points J=1,J=JJ)
    IO = 1
    IP = IO + 1
    IQ = IO + 2
    !
    DO JO = 1, JJ
        JL = JO - 2
        JM = JO - 1
        JP = JO + 1
        JQ = JO + 2
        !
        IF(JO.EQ.1) THEN
            !---- lower streamline inlet corner
            XOO = X(IO, JO)
            XOP = X(IO, JP)
            XOQ = X(IO, JQ)
            YOO = Y(IO, JO)
            YOP = Y(IO, JP)
            YOQ = Y(IO, JQ)
            !
            DETP = YPOS(JP) - YPOS(JO)
            DETQ = YPOS(JQ) - YPOS(JP)
            DXDE1 = (XOP - XOO) / DETP
            DXDE2 = (XOQ - XOP) / DETQ
            DYDE1 = (YOP - YOO) / DETP
            DYDE2 = (YOQ - YOP) / DETQ
            !--------- backwards difference at boundary
            DXDET = DXDE1 - DETP * (DXDE2 - DXDE1) / (DETP + DETQ)
            DYDET = DYDE1 - DETP * (DYDE2 - DYDE1) / (DETP + DETQ)
        ELSEIF(JO.EQ.JJ) THEN
            !---- upper streamline inlet corner
            XOO = X(IO, JO)
            XOM = X(IO, JM)
            XOL = X(IO, JL)
            YOO = Y(IO, JO)
            YOM = Y(IO, JM)
            YOL = Y(IO, JL)
            !
            DETM = YPOS(JO) - YPOS(JM)
            DETL = YPOS(JM) - YPOS(JL)
            DXDE1 = (XOO - XOM) / DETM
            DYDE1 = (YOO - YOM) / DETM
            DXDE2 = (XOM - XOL) / DETL
            DYDE2 = (YOM - YOL) / DETL
            !--------- backwards difference at boundary
            DXDET = DXDE1 + DETM * (DXDE1 - DXDE2) / (DETM + DETL)
            DYDET = DYDE1 + DETM * (DYDE1 - DYDE2) / (DETM + DETL)

        ELSE
            XOO = X(IO, JO)
            XPO = X(IP, JO)
            XQO = X(IQ, JO)
            YOO = Y(IO, JO)
            YPO = Y(IP, JO)
            YQO = Y(IQ, JO)
            !
            DXIP = XPOS(IP) - XPOS(IO)
            DXIQ = XPOS(IQ) - XPOS(IP)
            DXDX1 = (XPO - XOO) / DXIP
            DXDX2 = (XQO - XPO) / DXIQ
            DYDX1 = (YPO - YOO) / DXIP
            DYDX2 = (YQO - YPO) / DXIQ
            !---------- 2nd-order backward 3-point difference at boundary
            DXDXI = DXDX1 - DXIP * (DXDX2 - DXDX1) / (DXIP + DXIQ)
            DYDXI = DYDX1 - DXIP * (DYDX2 - DYDX1) / (DXIP + DXIQ)
        ENDIF
        !
        XOM = X(IO, JM)
        XOO = X(IO, JO)
        XOP = X(IO, JP)
        YOM = Y(IO, JM)
        YOO = Y(IO, JO)
        YOP = Y(IO, JP)
        !
        DETM = YPOS(JO) - YPOS(JM)
        DETP = YPOS(JP) - YPOS(JO)
        DETAV = 0.5 * (DETM + DETP)
        DXDET = 0.5 * (XOP - XOM) / DETAV
        DYDET = 0.5 * (YOP - YOM) / DETAV
        !
        AJA = DYDET * DXDXI - DXDET * DYDXI
        YAV = YOO
        !
        UG(IO, JO) = DXDXI / YAV / AJA
        VG(IO, JO) = DYDXI / YAV / AJA
    ENDDO
    !
    !-------- Velocities on (i=II) outlet plane (includes corner points J=1,J=JJ)
    IO = II
    IM = IO - 1
    IL = IO - 2
    DO JO = 1, JJ
        JL = JO - 2
        JM = JO - 1
        JP = JO + 1
        JQ = JO + 2
        !
        IF(JO.EQ.1) THEN
            !---- lower streamline outlet corner
            XOO = X(IO, JO)
            XOP = X(IO, JP)
            XOQ = X(IO, JQ)
            YOO = Y(IO, JO)
            YOP = Y(IO, JP)
            YOQ = Y(IO, JQ)
            !
            DETP = YPOS(JP) - YPOS(JO)
            DETQ = YPOS(JQ) - YPOS(JP)
            DXDE1 = (XOP - XOO) / DETP
            DXDE2 = (XOQ - XOP) / DETQ
            DYDE1 = (YOP - YOO) / DETP
            DYDE2 = (YOQ - YOP) / DETQ
            !--------- backwards difference at boundary
            DXDET = DXDE1 - DETP * (DXDE2 - DXDE1) / (DETP + DETQ)
            DYDET = DYDE1 - DETP * (DYDE2 - DYDE1) / (DETP + DETQ)
        ELSEIF(JO.EQ.JJ) THEN
            !---- upper streamline outlet corner
            XOO = X(IO, JO)
            XOM = X(IO, JM)
            XOL = X(IO, JL)
            YOO = Y(IO, JO)
            YOM = Y(IO, JM)
            YOL = Y(IO, JL)
            !
            DETM = YPOS(JO) - YPOS(JM)
            DETL = YPOS(JM) - YPOS(JL)
            DXDE1 = (XOO - XOM) / DETM
            DYDE1 = (YOO - YOM) / DETM
            DXDE2 = (XOM - XOL) / DETL
            DYDE2 = (YOM - YOL) / DETL
            !--------- backwards difference at boundary
            DXDET = DXDE1 + DETM * (DXDE1 - DXDE2) / (DETM + DETL)
            DYDET = DYDE1 + DETM * (DYDE1 - DYDE2) / (DETM + DETL)

        ELSE
            !
            XOO = X(IO, JO)
            XMO = X(IM, JO)
            XLO = X(IL, JO)
            YOO = Y(IO, JO)
            YMO = Y(IM, JO)
            YLO = Y(IL, JO)
            !
            DXIM = XPOS(IO) - XPOS(IM)
            DXIL = XPOS(IM) - XPOS(IL)
            DXDX1 = (XOO - XMO) / DXIM
            DXDX2 = (XMO - XLO) / DXIL
            DYDX1 = (YOO - YMO) / DXIM
            DYDX2 = (YMO - YLO) / DXIL
            !---------- 2nd-order 3-point difference for tangential velocity
            DXDXI = DXDX1 + DXIM * (DXDX1 - DXDX2) / (DXIM + DXIL)
            DYDXI = DYDX1 + DXIM * (DYDX1 - DYDX2) / (DXIM + DXIL)
        ENDIF
        !
        XOM = X(IO, JM)
        XOO = X(IO, JO)
        XOP = X(IO, JP)
        YOM = Y(IO, JM)
        YOO = Y(IO, JO)
        YOP = Y(IO, JP)
        !
        DETM = YPOS(JO) - YPOS(JM)
        DETP = YPOS(JP) - YPOS(JO)
        DETAV = 0.5 * (DETM + DETP)
        DXDET = 0.5 * (XOP - XOM) / DETAV
        DYDET = 0.5 * (YOP - YOM) / DETAV
        !
        AJA = DYDET * DXDXI - DXDET * DYDXI
        YAV = YOO
        !
        UG(IO, JO) = DXDXI / YAV / AJA
        VG(IO, JO) = DYDXI / YAV / AJA
    ENDDO
    !
    RETURN
END
! UVGRD
