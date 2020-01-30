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

      SUBROUTINE XSIETA(XM,X1,X2,X3,X4,YM,Y1,Y2,Y3,Y4,XSI,ETA,
     &                  XSI_X, XSI_Y, ETA_X, ETA_Y)
C............................................................
C     Subroutine for calculating the local bi-linear
C     coordinates of a marker inside a quadrilateral cell.
C     Also calculates the local coordinate derivatives.
C
C   Input:
C   ------
C     XM,YM          marker coordinates
C     X1,X2,X3,X4    cell vertex x-coordinates
C     Y1,Y2,Y3,Y4    cell vertex y-coordinates
C
C   Output:
C   -------
C     XSI,ETA        local cell coordinates at marker
C     XSI_X          dXSI/dX  at marker
C     XSI_Y          dXSI/dY  at marker
C     ETA_X          dETA/dX  at marker
C     ETA_Y          dETA/dY  at marker
C
C............................................................
C
C---- initial guess for xsi, eta  (center of cell)
      XS = 0.
      ET = 0.
C
C---- perform 12 or less Newton iterations for improving the guess
      DO 10 IXY=1, 12
        Z1 = (1.0 - XS)*(1.0 - ET)
        Z2 = (1.0 + XS)*(1.0 - ET)
        Z3 = (1.0 + XS)*(1.0 + ET)
        Z4 = (1.0 - XS)*(1.0 + ET)
        DR1 = -(Z1*X1 + Z2*X2 + Z3*X3 + Z4*X4  -  4.0*XM)
        DR2 = -(Z1*Y1 + Z2*Y2 + Z3*Y3 + Z4*Y4  -  4.0*YM)
        A11 = -(1.0-ET)*X1 + (1.0-ET)*X2 + (1.0+ET)*X3 - (1.0+ET)*X4
        A12 = -(1.0-XS)*X1 - (1.0+XS)*X2 + (1.0+XS)*X3 + (1.0-XS)*X4
        A21 = -(1.0-ET)*Y1 + (1.0-ET)*Y2 + (1.0+ET)*Y3 - (1.0+ET)*Y4
        A22 = -(1.0-XS)*Y1 - (1.0+XS)*Y2 + (1.0+XS)*Y3 + (1.0-XS)*Y4
        DETINV = 1.0 / (A11*A22 - A12*A21)
        DXS = DETINV * (DR1*A22 - A12*DR2)
        DET = DETINV * (A11*DR2 - DR1*A21)
        XS = XS + DXS
        ET = ET + DET
        IF(AMAX1(ABS(DXS),ABS(DET)) .LT. 1.0E-4) GO TO 11
   10 CONTINUE
      WRITE(6,*) 'Xsi-Eta conv failed at mrkr',XM,YM,DXS,DET
C
   11 XSI = XS
      ETA = ET
C
C---- derivatives are a free bonus from Newton procedure
      XSI_X =  4.0*A22*DETINV
      ETA_X = -4.0*A21*DETINV
      XSI_Y = -4.0*A12*DETINV
      ETA_Y =  4.0*A11*DETINV
C
      RETURN
      END ! XSIETA



      FUNCTION FINT(XSI,ETA,Q1,Q2,Q3,Q4)
C...................................................
C     Interpolates cell corner values Q1,Q2,Q3,Q4
C     to marker at local cell coodinates XSI,ETA
C...................................................
C
      Z1 = (1.0 - XSI)*(1.0 - ETA)
      Z2 = (1.0 + XSI)*(1.0 - ETA)
      Z3 = (1.0 + XSI)*(1.0 + ETA)
      Z4 = (1.0 - XSI)*(1.0 + ETA)
C
      FINT = (Q1*Z1 + Q2*Z2 + Q3*Z3 + Q4*Z4) * 0.25
C
      RETURN
      END ! FINT



      SUBROUTINE XYGRDFIND(XF,YF,IX,X,Y,II,JJ,IC,JC)
      DIMENSION X(IX,*), Y(IX,*)
      DIMENSION XQ(5), YQ(5)
C
C..............................................
C     Last-resort search through all
C     cells to locate point XF,YF in the grid
C
C     Returns indices IC,JC (to IC+1,JC+1) of cell
C     containing point or 0,0 for point outside
C     of grid.
C..............................................
C
      IC = 0
      JC = 0
C
C---- sweep over all cells in grid
      DO IO=1, II-1
        IP = IO+1
        DO JO=1, JJ-1
          JP = JO+1
C
          XQ(1) = X(IO,JO)
          XQ(2) = X(IP,JO)
          XQ(3) = X(IP,JP)
          XQ(4) = X(IO,JP)
          YQ(1) = Y(IO,JO)
          YQ(2) = Y(IP,JO)
          YQ(3) = Y(IP,JP)
          YQ(4) = Y(IO,JP)
C
C------ initial check to see if point is outside cell limits
C       (this disqualifies nearly all cells)
          XMAX = AMAX1( XQ(1),XQ(2),XQ(3),XQ(4) )
          XMIN = AMIN1( XQ(1),XQ(2),XQ(3),XQ(4) )
          YMAX = AMAX1( YQ(1),YQ(2),YQ(3),YQ(4) )
          YMIN = AMIN1( YQ(1),YQ(2),YQ(3),YQ(4) )
C
          IF( XF.GT.XMAX .OR. XF.LT.XMIN .OR. 
     &        YF.GT.YMAX .OR. YF.LT.YMIN      ) GO TO 10
C
C------ Inside test for polygon 
          NQ = 4
          CALL PNPOLY(XF,YF,XQ,YQ,NQ,INOUT)                            
          IF( INOUT.GE.0) THEN
C------ Found point in cell
            IC = IO
            JC = JO
            RETURN
          ENDIF
C
   10     CONTINUE
        END DO
      END DO
C
cc      WRITE(6,*) 'Point grid location failed. x y =',XF,YF
      RETURN
C
      END ! XYGRDFIND


C                                                                       
C     ..................................................................
C                                                                       
C        SUBROUTINE PNPOLY                                              
C                                                                       
C        PURPOSE                                                        
C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
C                                                                       
C        USAGE                                                          
C           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
C                                                                       
C        DESCRIPTION OF THE PARAMETERS                                  
C           PX      - X-COORDINATE OF POINT IN QUESTION.                
C           PY      - Y-COORDINATE OF POINT IN QUESTION.                
C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
C                     VERTICES OF POLYGON.                              
C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
C                     VERTICES OF POLYGON.                              
C           N       - NUMBER OF VERTICES IN THE POLYGON.                
C           INOUT   - THE SIGNAL RETURNED:                              
C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
C                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
C                                                                       
C        REMARKS                                                        
C           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
C           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
C           OPTIONALLY BE INCREASED BY 1.                               
C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
C           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
C
C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
C                                                                       
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
C           NONE                                                        
C                                                                       
C        METHOD                                                         
C           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
C           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
C           POINT IS INSIDE OF THE POLYGON.                             
C                                                                       
C     ..................................................................
C                                                                       
      SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)                            
      PARAMETER (MAXDIM=200)                                                      C
      REAL X(MAXDIM),Y(MAXDIM),XX(N),YY(N)                                    
      LOGICAL MX,MY,NX,NY                                               

C--- Output unit for printed messages                                 
      IF(N.GT.MAXDIM) THEN
        WRITE(*,1)                                                        
 1      FORMAT('WARNING:',I5,' too many polygon sides in PNPOLY')
        RETURN
      ENDIF
C                                                            
      DO I=1,N                                                        
        X(I) = XX(I)-PX                                                     
        Y(I) = YY(I)-PY
      END DO                                                     
C
      INOUT = -1                                                          
      DO 2 I = 1,N                                                        
        J = 1+MOD(I,N)                                                      
        MX = X(I).GE.0.0                                                    
        NX = X(J).GE.0.0                                                    
        MY = Y(I).GE.0.0                                                    
        NY = Y(J).GE.0.0                                                    
        IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX))      GO TO 2       
        IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
          INOUT = -INOUT                                                      
          GO TO 2                                                           
C
 3        TEST = (Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))
        IF(TEST.LT.0.0) THEN
          GO TO 2
        ELSEIF(TEST.EQ.0.0) THEN
          INOUT = 0                                                           
          RETURN
        ELSEIF(TEST.GT.0.0) THEN
          INOUT = -INOUT
        ENDIF                                                      
C
 2    CONTINUE                                                          
C
      RETURN                                                            
      END                                                               




      SUBROUTINE UVGRDC(IX,II,JJ,X,Y,XPOS,YPOS,UC,VC)
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION X(IX,*), Y(IX,*)
      DIMENSION XPOS(*), YPOS(*)
      DIMENSION UC(IX,*), VC(IX,*)
C-------------------------------------------------------------
C     Calculates velocity on cell centers for axisymmetric 
C     streamfunction grid by finte difference of streamfunction 
C     values.
C
C   Input:
C     IX        first grid-array dimension
C     II,JJ     grid i,j size
C     X(i,j)    grid z coordinates
C     Y(i,j)    grid r coordinates
C     XPOS(i)   grid s values
C     YPOS(j)   grid e values  (axisymmetric streamfunction)
C
C   Output:
C     UC(i,j)   u velocity at grid cell centers (II-1 x JJ-1)
C     VC(i,j)   v velocity at grid cell centers (II-1 x JJ-1)
C
C.............................................................
C
C   Uses solution to Thompson's grid-generation equations with
C   e(x,r) satisfying the axisymmetric streamfunction equation.
C   Hence,  e = constant  lines are streamlines.
C
C       s_zz + s_rr = 0
C       e_zz + e_rr = e_r / r
C
C   The equations for u(s,e), v(s,e) are 
C
C       u = z_e / (r J)
C
C       v = r_e / (r J)
C
C   where
C
C         J  =  z_s r_e  -  z_e r_s
C
C-------------------------------------------------------------
C
C
C------ go over all streamline cell centers
        DO 100 JO = 1, JJ-1
          JP = JO+1
C
C-------- for all cell centers on this streamline 
          DO 10 IO=1, II-1
            IP = IO+1
C
            XOO = X(IO,JO)
            XPO = X(IP,JO)
            XOP = X(IO,JP)
            XPP = X(IP,JP)
            YOO = Y(IO,JO)
            YPO = Y(IP,JO)
            YOP = Y(IO,JP)
            YPP = Y(IP,JP)
C
            XAV  = 0.25*( XOO + XOP + XPO + XPP)
            YAV  = 0.25*( YOO + YOP + YPO + YPP)
C
            DXAXI = 0.5*(-XOO - XOP + XPO + XPP)
            DYAXI = 0.5*(-YOO - YOP + YPO + YPP)
            DXAET = 0.5*(-XOO + XOP - XPO + XPP)
            DYAET = 0.5*(-YOO + YOP - YPO + YPP)
C
            DXI = XPOS(IP) - XPOS(IO)
            DET = YPOS(JP) - YPOS(JO)
C
            DXDXI = DXAXI / DXI
            DYDXI = DYAXI / DXI
            DXDET = DXAET / DET
            DYDET = DYAET / DET
C
            AJA = DYDET*DXDXI - DXDET*DYDXI
C
            UC(IO,JO) = DXDXI / YAV / AJA
            VC(IO,JO) = DYDXI / YAV / AJA
C
 10       CONTINUE
C
 100    CONTINUE
C
      RETURN
      END ! UVGRDC


      SUBROUTINE UVGRD(IX,II,JJ,X,Y,XPOS,YPOS,UG,VG)
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION X(IX,*), Y(IX,*)
      DIMENSION XPOS(*), YPOS(*)
      DIMENSION UG(IX,*), VG(IX,*)
C-------------------------------------------------------------
C     Calculates velocity on axisymmetric streamfunction grid
C     by finte difference of streamfunction values.
C
C   Input:
C     IX        first grid-array dimension
C     II,JJ     grid i,j size
C     X(i,j)    grid z coordinates
C     Y(i,j)    grid r coordinates
C     XPOS(i)   grid s values
C     YPOS(j)   grid e values  (axisymmetric streamfunction)
C
C   Output:
C     UG(i,j)   u velocity at grid points
C     VG(i,j)   v velocity at grid points
C
C.............................................................
C
C   Uses solution to Thompson's grid-generation equations with
C   e(x,r) satisfying the axisymmetric streamfunction equation.
C   Hence,  e = constant  lines are streamlines.
C
C       s_zz + s_rr = 0
C       e_zz + e_rr = e_r / r
C
C   The equations for u(s,e), v(s,e) are 
C
C       u = z_e / (rho r J)
C
C       v = r_e / (rho r J)
C
C   where
C
C         J  =  z_s r_e  -  z_e r_s
C
C-------------------------------------------------------------
C
C------ go over all interior streamlines
        DO 100 JO = 2, JJ-1
          JM = JO-1
          JP = JO+1
C
C-------- for all interior points on this streamline 
          DO 10 IO=2, II-1
            IM = IO-1
            IP = IO+1
C
            XMM = X(IM,JM)
            XOM = X(IO,JM)
            XPM = X(IP,JM)
            XMO = X(IM,JO)
            XOO = X(IO,JO)
            XPO = X(IP,JO)
            XMP = X(IM,JP)
            XOP = X(IO,JP)
            XPP = X(IP,JP)
            YMM = Y(IM,JM)
            YOM = Y(IO,JM)
            YPM = Y(IP,JM)
            YMO = Y(IM,JO)
            YOO = Y(IO,JO)
            YPO = Y(IP,JO)
            YMP = Y(IM,JP)
            YOP = Y(IO,JP)
            YPP = Y(IP,JP)
C
            DXIM = XPOS(IO) - XPOS(IM)
            DXIP = XPOS(IP) - XPOS(IO)
            DXIAV = 0.5*(DXIM+DXIP)
C
            DETM = YPOS(JO) - YPOS(JM)
            DETP = YPOS(JP) - YPOS(JO)
            DETAV = 0.5*(DETM+DETP)
C
            DXDET = 0.5*(XOP - XOM) / DETAV
            DYDET = 0.5*(YOP - YOM) / DETAV
            DXDXI = 0.5*(XPO - XMO) / DXIAV
            DYDXI = 0.5*(YPO - YMO) / DXIAV
C
            AJA = DYDET*DXDXI - DXDET*DYDXI
            YAV = YOO
C
            UG(IO,JO) = DXDXI / YAV / AJA
            VG(IO,JO) = DYDXI / YAV / AJA
C
 10       CONTINUE
C
 100    CONTINUE
C
C
C-------- Calculate velocities on grid boundaries (uses backward differences)
C
C-------- Velocities on lower streamline
        JO = 1
        DO IO = 2, II-1
          IP = IO+1
          IM = IO-1
C
          JP = JO+1
          JQ = JO+2
C
          XOO = X(IO,JO)
          XOP = X(IO,JP)
          XOQ = X(IO,JQ)
          YOO = Y(IO,JO)
          YOP = Y(IO,JP)
          YOQ = Y(IO,JQ)
C
          DETP = YPOS(JP)-YPOS(JO)
          DETQ = YPOS(JQ)-YPOS(JP)
C
          DXDE1 = ( XOP - XOO ) / DETP
          DXDE2 = ( XOQ - XOP ) / DETQ
          DYDE1 = ( YOP - YOO ) / DETP
          DYDE2 = ( YOQ - YOP ) / DETQ
C--- backwards difference at boundary
          DXDET = DXDE1 - DETP*(DXDE2 - DXDE1)/(DETP+DETQ)
          DYDET = DYDE1 - DETP*(DYDE2 - DYDE1)/(DETP+DETQ)
C
          XMO = X(IM,JO)
          XOO = X(IO,JO)
          XPO = X(IP,JO)
          YMO = Y(IM,JO)
          YOO = Y(IO,JO)
          YPO = Y(IP,JO)
C
          DXIM = XPOS(IO) - XPOS(IM)
          DXIP = XPOS(IP) - XPOS(IO)
          DXIAV = 0.5*(DXIM+DXIP)
          DXDXI = 0.5*(XPO - XMO) / DXIAV
          DYDXI = 0.5*(YPO - YMO) / DXIAV
C
          AJA = DYDET*DXDXI - DXDET*DYDXI
          YAV = YOO
C
          UG(IO,JO) = DXDXI / YAV / AJA
          VG(IO,JO) = DYDXI / YAV / AJA
        ENDDO
C
C-------- Velocities on upper streamline
        JO = JJ
        DO IO = 2, II-1
          IP = IO+1
          IM = IO-1
C
          JM = JO-1
          JL = JO-2
C
          XOO = X(IO,JO)
          XOM = X(IO,JM)
          XOL = X(IO,JL)
          YOO = Y(IO,JO)
          YOM = Y(IO,JM)
          YOL = Y(IO,JL)
C
          DETM = YPOS(JO)-YPOS(JM)
          DETL = YPOS(JM)-YPOS(JL)
          DXDE1 = ( XOO - XOM ) / DETM
          DYDE1 = ( YOO - YOM ) / DETM
          DXDE2 = ( XOM - XOL ) / DETL
          DYDE2 = ( YOM - YOL ) / DETL
C--- backwards difference at boundary
          DXDET = DXDE1 + DETM*(DXDE1 - DXDE2)/(DETM+DETL)
          DYDET = DYDE1 + DETM*(DYDE1 - DYDE2)/(DETM+DETL)
C
          XMO = X(IM,JO)
          XOO = X(IO,JO)
          XPO = X(IP,JO)
          YMO = Y(IM,JO)
          YOO = Y(IO,JO)
          YPO = Y(IP,JO)
C
          DXIM = XPOS(IO) - XPOS(IM)
          DXIP = XPOS(IP) - XPOS(IO)
          DXIAV = 0.5*(DXIM+DXIP)
          DXDXI = 0.5*(XPO - XMO) / DXIAV
          DYDXI = 0.5*(YPO - YMO) / DXIAV
C
          AJA = DYDET*DXDXI - DXDET*DYDXI
          YAV = YOO
C
          UG(IO,JO) = DXDXI / YAV / AJA
          VG(IO,JO) = DYDXI / YAV / AJA
        ENDDO
C
C-------- Velocities on(i=1) inlet plane (includes corner points J=1,J=JJ)
        IO = 1
        IP = IO+1
        IQ = IO+2
C
        DO JO = 1, JJ
          JL = JO-2
          JM = JO-1
          JP = JO+1
          JQ = JO+2
C
          IF(JO.EQ.1) THEN
C---- lower streamline inlet corner
           XOO = X(IO,JO)
           XOP = X(IO,JP)
           XOQ = X(IO,JQ)
           YOO = Y(IO,JO)
           YOP = Y(IO,JP)
           YOQ = Y(IO,JQ)
C
           DETP = YPOS(JP)-YPOS(JO)
           DETQ = YPOS(JQ)-YPOS(JP)
           DXDE1 = ( XOP - XOO ) / DETP
           DXDE2 = ( XOQ - XOP ) / DETQ
           DYDE1 = ( YOP - YOO ) / DETP
           DYDE2 = ( YOQ - YOP ) / DETQ
C--------- backwards difference at boundary
           DXDET = DXDE1 - DETP*(DXDE2 - DXDE1)/(DETP+DETQ)
           DYDET = DYDE1 - DETP*(DYDE2 - DYDE1)/(DETP+DETQ)
          ELSEIF(JO.EQ.JJ) THEN
C---- upper streamline inlet corner
           XOO = X(IO,JO)
           XOM = X(IO,JM)
           XOL = X(IO,JL)
           YOO = Y(IO,JO)
           YOM = Y(IO,JM)
           YOL = Y(IO,JL)
C
           DETM = YPOS(JO)-YPOS(JM)
           DETL = YPOS(JM)-YPOS(JL)
           DXDE1 = ( XOO - XOM ) / DETM
           DYDE1 = ( YOO - YOM ) / DETM
           DXDE2 = ( XOM - XOL ) / DETL
           DYDE2 = ( YOM - YOL ) / DETL
C--------- backwards difference at boundary
           DXDET = DXDE1 + DETM*(DXDE1 - DXDE2)/(DETM+DETL)
           DYDET = DYDE1 + DETM*(DYDE1 - DYDE2)/(DETM+DETL)

          ELSE
           XOO = X(IO,JO)
           XPO = X(IP,JO)
           XQO = X(IQ,JO)
           YOO = Y(IO,JO)
           YPO = Y(IP,JO)
           YQO = Y(IQ,JO)
C
           DXIP = XPOS(IP) - XPOS(IO)
           DXIQ = XPOS(IQ) - XPOS(IP)
           DXDX1 = ( XPO - XOO ) / DXIP
           DXDX2 = ( XQO - XPO ) / DXIQ
           DYDX1 = ( YPO - YOO ) / DXIP
           DYDX2 = ( YQO - YPO ) / DXIQ
C---------- 2nd-order backward 3-point difference at boundary
           DXDXI = DXDX1 - DXIP*(DXDX2 - DXDX1)/(DXIP+DXIQ)
           DYDXI = DYDX1 - DXIP*(DYDX2 - DYDX1)/(DXIP+DXIQ)
          ENDIF
C
          XOM = X(IO,JM)
          XOO = X(IO,JO)
          XOP = X(IO,JP)
          YOM = Y(IO,JM)
          YOO = Y(IO,JO)
          YOP = Y(IO,JP)
C
          DETM = YPOS(JO) - YPOS(JM)
          DETP = YPOS(JP) - YPOS(JO)
          DETAV = 0.5*(DETM+DETP)
          DXDET = 0.5*(XOP - XOM) / DETAV
          DYDET = 0.5*(YOP - YOM) / DETAV
C
          AJA = DYDET*DXDXI - DXDET*DYDXI
          YAV = YOO
C
          UG(IO,JO) = DXDXI / YAV / AJA
          VG(IO,JO) = DYDXI / YAV / AJA
        ENDDO
C
C-------- Velocities on (i=II) outlet plane (includes corner points J=1,J=JJ)
        IO = II
        IM = IO-1
        IL = IO-2
        DO JO = 1, JJ
          JL = JO-2
          JM = JO-1
          JP = JO+1
          JQ = JO+2
C
          IF(JO.EQ.1) THEN
C---- lower streamline outlet corner
           XOO = X(IO,JO)
           XOP = X(IO,JP)
           XOQ = X(IO,JQ)
           YOO = Y(IO,JO)
           YOP = Y(IO,JP)
           YOQ = Y(IO,JQ)
C
           DETP = YPOS(JP)-YPOS(JO)
           DETQ = YPOS(JQ)-YPOS(JP)
           DXDE1 = ( XOP - XOO ) / DETP
           DXDE2 = ( XOQ - XOP ) / DETQ
           DYDE1 = ( YOP - YOO ) / DETP
           DYDE2 = ( YOQ - YOP ) / DETQ
C--------- backwards difference at boundary
           DXDET = DXDE1 - DETP*(DXDE2 - DXDE1)/(DETP+DETQ)
           DYDET = DYDE1 - DETP*(DYDE2 - DYDE1)/(DETP+DETQ)
          ELSEIF(JO.EQ.JJ) THEN
C---- upper streamline outlet corner
           XOO = X(IO,JO)
           XOM = X(IO,JM)
           XOL = X(IO,JL)
           YOO = Y(IO,JO)
           YOM = Y(IO,JM)
           YOL = Y(IO,JL)
C
           DETM = YPOS(JO)-YPOS(JM)
           DETL = YPOS(JM)-YPOS(JL)
           DXDE1 = ( XOO - XOM ) / DETM
           DYDE1 = ( YOO - YOM ) / DETM
           DXDE2 = ( XOM - XOL ) / DETL
           DYDE2 = ( YOM - YOL ) / DETL
C--------- backwards difference at boundary
           DXDET = DXDE1 + DETM*(DXDE1 - DXDE2)/(DETM+DETL)
           DYDET = DYDE1 + DETM*(DYDE1 - DYDE2)/(DETM+DETL)

          ELSE
C
           XOO = X(IO,JO)
           XMO = X(IM,JO)
           XLO = X(IL,JO)
           YOO = Y(IO,JO)
           YMO = Y(IM,JO)
           YLO = Y(IL,JO)
C
           DXIM = XPOS(IO) - XPOS(IM)
           DXIL = XPOS(IM) - XPOS(IL)
           DXDX1 = ( XOO - XMO ) / DXIM
           DXDX2 = ( XMO - XLO ) / DXIL
           DYDX1 = ( YOO - YMO ) / DXIM
           DYDX2 = ( YMO - YLO ) / DXIL
C---------- 2nd-order 3-point difference for tangential velocity
           DXDXI = DXDX1 + DXIM*(DXDX1 - DXDX2)/(DXIM+DXIL)
           DYDXI = DYDX1 + DXIM*(DYDX1 - DYDX2)/(DXIM+DXIL)
          ENDIF
C
          XOM = X(IO,JM)
          XOO = X(IO,JO)
          XOP = X(IO,JP)
          YOM = Y(IO,JM)
          YOO = Y(IO,JO)
          YOP = Y(IO,JP)
C
          DETM = YPOS(JO) - YPOS(JM)
          DETP = YPOS(JP) - YPOS(JO)
          DETAV = 0.5*(DETM+DETP)
          DXDET = 0.5*(XOP - XOM) / DETAV
          DYDET = 0.5*(YOP - YOM) / DETAV
C
          AJA = DYDET*DXDXI - DXDET*DYDXI
          YAV = YOO
C
          UG(IO,JO) = DXDXI / YAV / AJA
          VG(IO,JO) = DYDXI / YAV / AJA
        ENDDO
C
      RETURN
      END ! UVGRD

