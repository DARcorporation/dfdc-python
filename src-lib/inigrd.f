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


      SUBROUTINE INIGRD
C------------------------------------------------------
C     Sets up vortex wake grid system from rotor and wake
C     geometry 
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION S(IPX)
C
      IF(LDBG) THEN
        WRITE(*,*) ' '
        WRITE(*,*) 'INIGRD setup grid for slipstream'
C
        WRITE(LUNDBG,*) ' '
        WRITE(LUNDBG,*) 'INIGRD setup grid for slipstream'
      ENDIF
C
C---- Find most-upstream rotor point on duct wall and set as start of grid
      NRUPSTRM = 1
      IPUP = IPROTDW(NRUPSTRM)
      DO NR = 2, NROTOR
        IF(IPROTDW(NR).LT.IPUP) THEN
          NRUPSTRM = NR
          IPUP = IPROTDW(NR)
        ENDIF
      END DO
C---- Setup rotor at upstream end of grid
      DO NR = 1, NROTOR
        IGROTOR(NR) = 1 + IPROTDW(NR)-IPUP
      END DO      
C
      NR = NRUPSTRM
C---- Indices of rotor line on CB and duct define upstream grid boundary
C---- find axial location of rotor and TE's on CB and duct walls
      IEL = 1       
      IP1CB = IPFRST(IEL)
      IP2CB = IPLAST(IEL)
      IPCB = IPROTCB(NR)
      IF(IPCB.EQ.0) THEN
        WRITE(*,*) 'Error locating rotor on CB wall'
        STOP
      ENDIF
      IF(LDBG) WRITE(LUNDBG,*) 'Rotor on CB @ ',IPCB
      XCB = XP(IPCB)
      YCB = YP(IPCB)
      XTECB = XP(IP1CB)
cc    IF(LDBG) WRITE(*,*) 'XR, XTECB ',XCB,XTECB
C
      IEL = 2       
      IP1DW = IPFRST(IEL)
      IP2DW = IPLAST(IEL)
      IPDW = IPROTDW(NR)
      IF(IPDW.EQ.0) THEN
        WRITE(*,*) 'Error locating rotor on Duct wall'
        STOP
      ENDIF
      IF(LDBG) WRITE(LUNDBG,*) 'Rotor on Duct wall @ ',IPDW
      XDW = XP(IPDW)
      YDW = YP(IPDW)
      XTEDW = XP(IP2DW)
cc    IF(LDBG) WRITE(*,*) 'XR XTEDW ',XDW,XTEDW
C
C
C---- Define lower streamline (J=1)
      I = IGROTOR(NR) - 1
      J = 1
      DO IP = IPCB, IP1CB, -1
        I = I + 1
        XG(I,J) = XP(IP)
        YG(I,J) = YP(IP)
        IP2IG(IP) = I
      END DO
      XLAST = XG(I,J)
      YLAST = YG(I,J)
      IGTECB = I
cc      write(*,*) 'IGTECB,II defined as ',IGTECB,II
cc    IF(LDBG) WRITE(*,*) 'L IGTECB,XYLAST ',IGTECB,XLAST,YLAST
C
C---- Add points at X values of duct wall if it extends downstream of CB
      IPD = IPDW + IGTECB - 1
      IF(IPD.LT.IP2DW) THEN
       DO IP = IPD+1, IP2DW
         I = I + 1
         XG(I,J) = XP(IP)
         YG(I,J) = YLAST
       END DO
       XLAST = XG(I,J)
       YLAST = YG(I,J)
       ITE = I
cc     IF(LDBG) WRITE(*,*) 'L ITE,XYLAST ',ITE,XLAST,YLAST
      ENDIF
C
C---- Check downstream wake end location
      IF(LDBG) WRITE(LUNDBG,*) 'INIGRD XWAKE ',XWAKE
      IF(XDWKLEN.LT.0.0 .OR. XDWKLEN.GT.5.0) THEN
        WRITE(*,*) 'Rotor wake length out of bounds ',XDWKLEN
        WRITE(*,*) '  CB   TE at X,Y =', XBTE(1),YBTE(1)
        WRITE(*,*) '  Duct TE at X,Y =', XBTE(2),YBTE(2)
        WRITE(*,*) '  Rotor Diameter =',2.0*RTIP(1)
        CALL ASKR('Enter length of wake from TE dXlen/D^',XDWKLEN)
      ENDIF
      XWAKE = MAX(XBTE(1),XBTE(2)) + XDWKLEN*2.0*RTIP(1)
C
C---- Add points downstream to end of wake
      DX = XG(I,J)-XG(I-1,J)
      DY = YG(I,J)-YG(I-1,J)
      DS = SQRT(DX*DX + DY*DY)
C
      SMAX = XWAKE-XLAST
      DXW  = SMAX / FLOAT(NWAKE-1) 
      DS1  = 0.2*DXW
      CALL SETEXP(S,DS1,SMAX,NWAKE)
      DO K = 2, NWAKE
        I = I + 1
        XG(I,J) = XLAST + S(K)
        YG(I,J) = YLAST
      END DO
      ILWR = I
C
C
C---- Define upper streamline (J=JJ=NRP from rotor line points)
      JJ = NRP
C
      I = IGROTOR(NR) - 1
      J = JJ
      DO IP = IPDW, IP2DW
        I = I + 1
        XG(I,J) = XP(IP)
        YG(I,J) = YP(IP)
        IP2IG(IP) = I
      END DO
      XLAST = XG(I,J)
      YLAST = YG(I,J)
      IGTEDW = I
cc      write(*,*) 'IGTEDW,II defined as ',IGTEDW,II
cc    IF(LDBG) WRITE(*,*) 'U IGTEDW,XYLAST ',IGTEDW,XLAST,YLAST
C
C---- Add points at X values of CB wall if it extends downstream of duct
      IPC = IPCB - IGTEDW + 1
      IF(IPC.GT.IP1CB) THEN
       DO IP = IPC-1, IP1CB, -1
         I = I + 1
         XG(I,J) = XP(IP)
         YG(I,J) = YLAST
       END DO
       XLAST = XG(I,J)
       YLAST = YG(I,J)
       ITE = I
cc     IF(LDBG) WRITE(*,*) 'U ITE,XYLAST ',ITE,XLAST,YLAST
      ENDIF
C
C---- Add points downstream to end of wake
      DX = XG(I,J)-XG(I-1,J)
      DY = YG(I,J)-YG(I-1,J)
      DS = SQRT(DX*DX + DY*DY)
C
      SMAX = XWAKE-XLAST
      DXW  = SMAX / FLOAT(NWAKE-1) 
      DS1  = 0.2*DXW
      CALL SETEXP(S,DS1,SMAX,NWAKE)
      DO K = 2, NWAKE
        I = I + 1
        XG(I,J) = XLAST + S(K)
        YG(I,J) = YLAST
      END DO
      IUPR = I
C
      IF(ILWR.NE.IUPR) THEN
        WRITE(*,*) 'INIGRD: incompatible # points on J=1 and J=JJ',
     &             ILWR, IUPR
        WRITE(*,*) 'IP1CB,IP2CB,IPCB ',IP1CB,IP2CB,IPCB
        WRITE(*,*) 'IP1DW,IP2DW,IPDW ',IP1DW,IP2DW,IPDW
        STOP
      ENDIF
C---- Set total # of streamwise points
      II = IUPR
C
C
C---- Add points at inlet and outlet boundaries
C---- Inlet points are simply the rotor line points XR,YR
      I = 1
      XGMIN = XRP(1,NR)
      DO J = 1, JJ
        XG(I,J) = XRP(J,NR)
        YG(I,J) = YRP(J,NR)
        XGMIN = MIN(XGMIN,XG(I,J))
      END DO
C
C---- Sweep downstream setting intermediate J lines preserving relative 
C     massflow between upper and lower walls 
      AMASS1 = PI*(YG(1,JJ)**2-YG(1,1)**2)
      DO I = 2, II
        XL = XG(I,1)
        YL = YG(I,1)
        XU = XG(I,JJ)
        YU = YG(I,JJ)
        AMASS = PI*(YU**2-YL**2)
C---- Create radial spacing by preserving relative area 
        DO J = 2, JJ-1
          DMASS1 = PI*(YG(1,J)**2-YG(1,J-1)**2)
          FRC = DMASS1/AMASS1
          YY = SQRT(FRC*AMASS/PI + YG(I,J-1)**2)
          FRAC = (YY-YL)/(YU-YL)
          XG(I,J) = XL + FRAC*(XU-XL)
          YG(I,J) = YY
        END DO
      END DO
C
C---- Done, initial grid defined!
      IF(LDBG) THEN
        WRITE(*,*) 'INIGRD grid defined with II,JJ ',II,JJ
        WRITE(LUNDBG,*) 'INIGRD grid defined with II,JJ ',II,JJ
      ENDIF
C
C     Add initial grid smoothing relaxation?
C---- Grid smoothing by elliptic solver
      CALL RLXGRD
C
C---- invalidate existing solution to incorporate new vortex wakes
      LNCVP = .FALSE.
      LQAIC = .FALSE.
      LQGIC = .FALSE.
      LQCNT = .FALSE.
      LGSYS = .FALSE.
      LGAMU = .FALSE.
      LGAMA = .FALSE.
C
      RETURN
      END ! INIGRD


      SUBROUTINE UPDGRD
C------------------------------------------------------
C     Updates grid upper and lower boundaries with 
C     points from current upper and lower vortex wakes
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) THEN
        WRITE(*,*) ' '
        WRITE(*,*) 'Updating grid system boundaries'
      ENDIF
C
C---- Redefine lower streamline (J=1) with IR=1 wake points
      J = 1
      IEL = IR2IEL(1)
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      DO IP = IP1, IP2
        I = IGTECB + (IP-IP1)
        XG(I,J) = XP(IP)
        YG(I,J) = YP(IP)
      END DO
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'Updated lower boundary from element ',IEL
      ENDIF
C
C---- Redefine upper streamline (J=JJ) with IR=NRP wake points
      J = JJ
      IEL = IR2IEL(NRP)
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      DO IP = IP1, IP2
        I = IGTEDW + (IP-IP1)
        XG(I,J) = XP(IP)
        YG(I,J) = YP(IP)
      END DO
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'Updated upper boundary from element ',IEL
      ENDIF
C
C---- Grid boundaries redefined!
      IF(LDBG) WRITE(*,*) 'UPDGRD grid boundaries updated'
C
C---- Grid smoothing by elliptic solver
      CALL RLXGRD
C
      RETURN
      END ! UPDGRD




      SUBROUTINE RLXGRD
C------------------------------------------------------
C     Runs elliptic grid solution to smooth grid
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(II.LE.0 .OR. JJ.LE.0) THEN
        WRITE(*,*) 'RLXGRD:  no grid defined II,JJ ',II,JJ
        RETURN
      ENDIF
C
C---- Set up streamfunction array
      YPOS(1) = 0.0
      DO J = 2, JJ
        YPOS(J) = YG(1,J)**2 - YG(1,J-1)**2 + YPOS(J-1)
      END DO
C---- Set up streamwise spacing array
      DO I = 1, II
        XPOS(I) = XG(I,1)
      END DO
C
C---- Run grid solver
      ITMAXS = -400
      TOLER = 1.0E-9
      KBCINL = 0   ! inlet plane Dirichlet, fix nodes
      KBCOUT = 1   ! outlet plane Neumann, move nodes on boundary
      KBCBOT = 0   ! bottom streamline Dirichlet, fix nodes
      KBCTOP = 0   ! top streamline Dirichlet, fix nodes
c      WRITE(*,*) 
c      WRITE(*,*) 'Grid relaxation iteration'
      CALL AXELL(IX,II,JJ,XG,YG,XPOS,YPOS,ITMAXS,TOLER,
     &           KBCINL,KBCOUT,KBCBOT,KBCTOP)
C
      RETURN
      END ! RLXGRD

c     &  XG(IX,JX),    YG(IX,JX),
c     &  QG(IX,JX),    QXG(IX,JX),   QYG(IX,JX),  QTG(IX,JX),  
c     &  RG(IX,JX),    PG(IX,JX),    POG(IX,JX),
c     &  BGAMG(IX,JX), DSG(IX,JX),   DHG(IX,JX),
c     &  XPOS(IX),     YPOS(JX)


      SUBROUTINE SETGRDFLW
C---------------------------------------------------------
C     Sets grid flow data from circulation and entropy on 
C     rotor lines
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Clear grid data for accumulation over sources
      CALL CLRGRDFLW
C
C---- Set B*GAM circulation, delta entropy, and delta enthalpy on gridC
C     from rotor and drag sources

      DO N = 1, NROTOR
C
       IEL = IELROTOR(N)
       IC1 = ICFRST(IEL)
       IC2 = ICLAST(IEL)
C---- Set values at rotor I station due to blade row 
       IG = IGROTOR(N)
C
       NNC = IC2-IC1+1
       IF(NNC.NE.JJ-1) THEN
         WRITE(*,*) 'Rotor source element-grid mismatch ',NNC,JJ-1
         STOP
       ENDIF
C
       DO J = 1, JJ-1
        IR = J
C---- Circulation and enthalpy
        DBGAMG = BGAM(IR,N)
        DDHG   = OMEGA(N)*BGAM(IR,N)*PI2I
C---- Entropy from drag sources on rotor
        IC = IC1 + IR-1
        IP1 = IPCO(IC)
        IP2 = IPCP(IC)
        VMC  = SQRT(QC(1,IC)**2 + QC(2,IC)**2)
        SIGC = 0.5*(SIG(IP1)+SIG(IP2))
        DDSG = VMC*SIGC
C
C---- Convect values downstream from rotor 
        DO I = IG, II-1
          BGAMG(I,J) = BGAMG(I,J) + DBGAMG
          DHG(I,J)   = DHG(I,J)   + DDHG
          DSG(I,J)   = DSG(I,J)   + DDSG
cc        RG(I,J)    = RG(1,J)
        END DO
C
       END DO
      END DO
C
      RETURN
      END ! SETGRDFLW


      SUBROUTINE CLRGRDFLW
C---------------------------------------------------------
C     Clears (initializes to zero) grid flow BGAMG,DHG,DSG
C     This is normally done before accumulating disk 
C     contributions
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Clear grid data for accumulation of circulation, enthalpy, entropy
      DO J = 1, JJ
        DO I = 1, II
          BGAMG(I,J) = 0.0
          DHG(I,J)   = 0.0
          DSG(I,J)   = 0.0
          RG(I,J)    = RHO
        END DO
      END DO
C
      RETURN
      END ! CLRGRDFLW


      SUBROUTINE ROTBG2GRD(N)
C---------------------------------------------------------
C     Updates grid flow (circulation) from BGAM on rotor N
C     rotor lines
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Set B*GAM circulation, delta enthalpy on grid from blade/rotors
C---- Set values at grid IG station due to blade row 
      IG = IGROTOR(N)
C
      DO J = 1, JJ-1
       IR = J
C---- Circulation and enthalpy
       DBGAMG = BGAM(IR,N)
       DDHG   = OMEGA(N)*BGAM(IR,N)*PI2I
C
C---- Convect values downstream from rotor 
       DO I = IG, II-1
         BGAMG(I,J) = BGAMG(I,J) + DBGAMG
         DHG(I,J)   = DHG(I,J)   + DDHG
       END DO
C
      END DO
C
      RETURN
      END ! ROTBG2GRD
  


      SUBROUTINE GETGRDFLW
C---------------------------------------------------------
C     Sets grid initial flow data from rotor line conditions
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      CALL UVGRDC(IX,II,JJ,XG,YG,XPOS,YPOS,QXG,QYG)
C---- Print grid velocities
      WRITE(44,*) 'I,J, XC,YC, U, V'
      DO I = 1, II-1
        IO = I
        IP = I + 1
        DO J = 1, JJ-1
          JO = J
          JP = J + 1
          XOO = XG(IO,JO)
          XPO = XG(IP,JO)
          XOP = XG(IO,JP)
          XPP = XG(IP,JP)
          YOO = YG(IO,JO)
          YPO = YG(IP,JO)
          YOP = YG(IO,JP)
          YPP = YG(IP,JP)
          XCG = 0.25*(XOO+XPO+XOP+XPP)
          YCG = 0.25*(YOO+YPO+YOP+YPP)
          WRITE(44,100) I,J,XCG,YCG,QXG(I,J),QYG(I,J)
        END DO
      END DO
 100  FORMAT(2I5,5(1X,F11.6))
C
      RETURN
      END ! GETGRDFLW



      SUBROUTINE INLSFCN
C---------------------------------------------------------
C     Calculates streamfunction (YPOS) at grid inlet using
C     velocities at rotor and current grid points
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      NR = NRUPSTRM
C---- Set up streamfunction array using grid and rotor line velocity
      YPOS(1) = 0.0
      DO J = 1, JJ
        VX = VABS(1,J-1,NR)
        YPOS(J) = 0.5*VX*(YG(1,J)**2 - YG(1,J-1)**2) + YPOS(J-1) 
cc        YPOS(J) = 0.5*VX*YG(1,J)**2
      END DO
      YSCL = 1.0/YPOS(JJ)
cc      DO J = 1, JJ
cc        write(*,*) 'j ypos ',j,ypos(j)
cc        YPOS(J) = YSCL*YPOS(J)
cc      END DO
C
      RETURN
      END ! INLSFCN


      SUBROUTINE RLXGRD2
C------------------------------------------------------
C     Runs elliptic grid solution to smooth grid
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(II.LE.0 .OR. JJ.LE.0) THEN
        WRITE(*,*) 'RLXGRD:  no grid defined II,JJ ',II,JJ
        RETURN
      ENDIF
C
C---- Set up streamfunction array
      CALL INLSFCN
C
C---- Set up streamwise spacing array
      DO I = 1, II
        XPOS(I) = XG(I,1)
      END DO
C
C---- Run grid solver
      ITMAXS = -400
      TOLER = 1.0E-9
      KBCINL = 0   ! inlet plane Dirichlet, fix nodes
      KBCOUT = 1   ! outlet plane Neumann, move nodes on boundary
      KBCBOT = 0   ! bottom streamline Dirichlet, fix nodes
      KBCTOP = 0   ! top streamline Dirichlet, fix nodes
      CALL AXELL(IX,II,JJ,XG,YG,XPOS,YPOS,ITMAXS,TOLER,
     &           KBCINL,KBCOUT,KBCBOT,KBCTOP)
C
      RETURN
      END ! RLXGRD2


      SUBROUTINE PLTGRID(SSIZEL,
     &                   XG,YG,IX,II,JJ,
     &                   XOFA,YOFA,FACA)
C----------------------------------------------
C     Plots grid in Cp vs x plot
C----------------------------------------------
      DIMENSION XG(IX,*),YG(IX,*)
      DIMENSION SSIZEL(0:7)
C
      INCLUDE 'MASKS.INC'
C 
      IF(II.LE.1 .OR. JJ.LE.1) THEN
        WRITE(*,*) 'No grid defined in PLTGRID II,JJ ',II,JJ
        RETURN
      ENDIF
C
      CALL GETCOLOR(ICOL0)
      CALL GETPAT(IPAT0)
      CALL NEWPEN(2)
      CALL NEWCOLORNAME('GRAY80')
      CALL NEWPAT(LMASK1)
C
C------ plot I lines of grid
      DO J = 1, JJ
        CALL PLOT((XG(1,J)+XOFA)*FACA,(YG(1,J)+YOFA)*FACA,3)
        DO I = 2, II
          CALL PLOT((XG(I,J)+XOFA)*FACA,(YG(I,J)+YOFA)*FACA,2)
        END DO
      END DO
C
C------ plot J lines of grid
      DO I = 1, II
        CALL PLOT((XG(I,1)+XOFA)*FACA,(YG(I,1)+YOFA)*FACA,3)
        DO J = 2, JJ
          CALL PLOT((XG(I,J)+XOFA)*FACA,(YG(I,J)+YOFA)*FACA,2)
        END DO
      END DO
      CALL NEWCOLOR(ICOL0)
      CALL NEWPAT(IPAT0)
C
      RETURN
      END ! PLTGRID


      SUBROUTINE AXELL(IX,II,JJ,X,Y,XPOS,YPOS,ITMAXS,TOLER,
     &                 KBCINL,KBCOUT,KBCBOT,KBCTOP)
      IMPLICIT REAL (A-H,M,O-Z)
      DIMENSION X(IX,*), Y(IX,*)
      DIMENSION XPOS(*), YPOS(*)
C-------------------------------------------------------------
C     Axisymmetric elliptic grid smoother.
C     Uses SLOR implicit along i direction.
C
C   Input:
C     IX        first grid-array dimension
C     II,JJ     grid i,j size
C     X(i,j)    grid z coordinates
C     Y(i,j)    grid r coordinates
C     XPOS(i)   grid s values
C     YPOS(j)   grid e values  (axisymmetric streamfunction)
C                ( YPOS(j) = (j-1)**2 gives uniform spacing in r )
C     ITMAXS    max number of SLOR smoothing passes
C               if ITMAXS < 0, overrelaxation is NOT used
C     TOLER     convergence tolerance (movement distance)
C     KBCINL    1 = Neumann, 0 = Dirichlet   BC applied at i=1 
C     KBCOUT    1 = Neumann, 0 = Dirichlet   BC applied at i=II
C     KBCBOT    1 = Neumann, 0 = Dirichlet   BC applied at j=1 
C     KBCTOP    1 = Neumann, 0 = Dirichlet   BC applied at j=JJ
C
C   Output:
C     X(i,j)    smoothed-grid z coordinates
C     Y(i,j)    smoothed-grid r coordinates
C
C.............................................................
C
C   Solves Thompson's grid-generation equations with source term 
C   which makes e(x,r) satisfy the axisymmetric streamfunction 
C   equation.  Hence,  e = constant  lines are streamlines.
C
C       s_zz + s_rr = 0
C       e_zz + e_rr = e_r / r
C
C   These equations for s(z,r), e(z,r) are inverted to
C
C       a z_ss  -  2b z_se  +  c z_ee  =  (J/r) z_e z_s
C
C       a r_ss  -  2b r_se  +  c r_ee  =  (J/r) r_e z_s
C
C   where
C                    2        2                    2        2
C         a  =  (z_e)  + (r_e)           c  = (z_s)  + (r_s)
C
C         b  =  z_e z_s  +  r_e r_s
C
C         J  =  z_s r_e  -  z_e r_s
C
C
C   The equivalent form which is actually differenced is
C
C       a z_ss  -  2b z_se  +  c/r (r z_e)_e  -  b/r r_s z_e  =  0
C
C       a r_ss  -  2b r_se  +  c/r (r r_e)_e  -  b/r r_s r_e  =  0
C
C-------------------------------------------------------------
      PARAMETER (IDIM=500)
      DIMENSION C(IDIM),D(2,IDIM)
C
      DIMENSION XT(IDIM),YT(IDIM),ST(IDIM),XTS(IDIM),YTS(IDIM),STI(IDIM)
      DIMENSION XB(IDIM),YB(IDIM),SB(IDIM),XBS(IDIM),YBS(IDIM),SBI(IDIM)
C
      IF(II.GT.IDIM) STOP 'ELLIP: Array overflow.  Increase IDIM.'
C
      ITMAX = IABS(ITMAXS)
C
C---- convergence tolerance for each phase
      DSET1 = 1.0E-1
      DSET2 = 5.0E-3
      DSET3 = 5.0E-7
C
      IF(ITMAXS.GT.0) THEN
C----- over-relaxation parameter for each phase
       RLX1 = 1.00       !          DMAX > DSET1
       RLX2 = 1.10       !  DSET1 > DMAX > DSET2
       RLX3 = 1.40       !  DSET2 > DMAX > DSET3
CCC    STOP              !  DSET3 > DMAX
      ELSE
       RLX1 = 1.0
       RLX2 = 1.0
       RLX3 = 1.0
      ENDIF
C
C---- spline bottom & top wall shapes
      JO = 1
      DO IO = 1, II
        XB(IO) = X(IO,JO)
        YB(IO) = Y(IO,JO)
      ENDDO
C
      DX = ABS(XB(II)-XB(1))
      DY = ABS(YB(II)-YB(1))
      DXY = MAX(XB(1),YB(1),DX,DY)
C
      CALL SCALC(XB,YB,SB,II)
      CALL SPLINE(XB,XBS,SB,II)
      CALL SPLINE(YB,YBS,SB,II)
      DO IO = 1, II
        SBI(IO) = SB(IO)
      ENDDO
C
      JO = JJ
      DO IO = 1, II
        XT(IO) = X(IO,JO)
        YT(IO) = Y(IO,JO)
      ENDDO
      CALL SCALC(XT,YT,ST,II)
      CALL SPLINE(XT,XTS,ST,II)
      CALL SPLINE(YT,YTS,ST,II)
      DO IO = 1, II
        STI(IO) = ST(IO)
      ENDDO
C
      RLX = RLX1
C
      DO 1000 IPASS = 1, ITMAX
C
        DMAX = 0.
C
        JO = 1
        IF(KBCBOT.EQ.1 .AND. JJ.GT.2) THEN
C-------- relax to get dX/dj = 0 on lower streamline
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
C
            YOO = Y(IO,JO)
            YOP = Y(IO,JP)
            YOQ = Y(IO,JQ)
C
            XS = DEVAL(SBI(IO),XB,XBS,SB,II)
            YS = DEVAL(SBI(IO),YB,YBS,SB,II)
C
            DETP = YPOS(JP)-YPOS(JO)
            DETQ = YPOS(JQ)-YPOS(JP)
C
            DXDE1 = ( XOP - XOO ) / DETP
            DXDE2 = ( XOQ - XOP ) / DETQ
C
            DYDE1 = ( YOP - YOO ) / DETP
            DYDE2 = ( YOQ - YOP ) / DETQ
C
            REZ   = XS*(DXDE1 - DETP*(DXDE2 - DXDE1)/(DETP+DETQ))
     &            + YS*(DYDE1 - DETP*(DYDE2 - DYDE1)/(DETP+DETQ))
C
            Z_XOO = XS*(-1.0/DETP - 1.0/(DETP+DETQ))
            Z_YOO = YS*(-1.0/DETP - 1.0/(DETP+DETQ))
C
            Z_S = Z_XOO*XS + Z_YOO*YS
C
            DS = -REZ/Z_S
C
            SBI(IO) = SBI(IO) + DS
C
            X(IO,JO) = SEVAL(SBI(IO),XB,XBS,SB,II)
            Y(IO,JO) = SEVAL(SBI(IO),YB,YBS,SB,II)
          ENDDO
        ENDIF
C
        JO = JJ
        IF(KBCTOP.EQ.1 .AND. JJ.GT.2) THEN
C-------- relax to get dX/dj = 0 on upper streamline
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
C
            YOO = Y(IO,JO)
            YOM = Y(IO,JM)
            YOL = Y(IO,JL)
C
            XS = DEVAL(STI(IO),XT,XTS,ST,II)
            YS = DEVAL(STI(IO),YT,YTS,ST,II)
C
            DETM = YPOS(JO)-YPOS(JM)
            DETL = YPOS(JM)-YPOS(JL)
C
            DXDE1 = ( XOO - XOM ) / DETM
            DXDE2 = ( XOM - XOL ) / DETL
C
            DYDE1 = ( YOO - YOM ) / DETM
            DYDE2 = ( YOM - YOL ) / DETL
C
            REZ   = XS*(DXDE1 + DETM*(DXDE1 - DXDE2)/(DETM+DETL))
     &            + YS*(DYDE1 + DETM*(DYDE1 - DYDE2)/(DETM+DETL))
C
            Z_XOO = XS*(1.0/DETM + 1.0/(DETM+DETL))
            Z_YOO = YS*(1.0/DETM + 1.0/(DETM+DETL))
C
            Z_S = Z_XOO*XS + Z_YOO*YS
C
            DS = -REZ/Z_S
C
            STI(IO) = STI(IO) + DS
C
            X(IO,JO) = SEVAL(STI(IO),XT,XTS,ST,II)
            Y(IO,JO) = SEVAL(STI(IO),YT,YTS,ST,II)
          ENDDO
        ENDIF
C
C
C------ go over all interior streamlines
        DO 100 JO = 2, JJ-1
          JM = JO-1
          JP = JO+1
C
          IF(KBCINL.EQ.1) THEN
C---------- Neumann BC is specified on i=1 plane:  relax node positions
C
            IO = 1
            IP = IO+1
            IQ = IO+2
C
            XOO = X(IO,JO)
            XPO = X(IP,JO)
            XQO = X(IQ,JO)
C
            YOO = Y(IO,JO)
            YPO = Y(IP,JO)
            YQO = Y(IQ,JO)
C
            XS = X(IO,JP) - X(IO,JM)
            YS = Y(IO,JP) - Y(IO,JM)
C
            XSN = XS / SQRT(XS**2 + YS**2)
            YSN = YS / SQRT(XS**2 + YS**2)
C
            DXIP = XPOS(IP) - XPOS(IO)
            DXIQ = XPOS(IQ) - XPOS(IP)
C
            DXDX1 = ( XPO - XOO ) / DXIP
            DXDX2 = ( XQO - XPO ) / DXIQ
C
            DYDX1 = ( YPO - YOO ) / DXIP
            DYDX2 = ( YQO - YPO ) / DXIQ
C
C---------- 2nd-order 3-point difference for tangential velocity
            REZ   = XSN*(DXDX1 - DXIP*(DXDX2 - DXDX1)/(DXIP+DXIQ))
     &            + YSN*(DYDX1 - DXIP*(DYDX2 - DYDX1)/(DXIP+DXIQ))
C
            Z_XOO = XSN*(-1.0/DXIP - 1.0/(DXIP+DXIQ))
            Z_YOO = YSN*(-1.0/DXIP - 1.0/(DXIP+DXIQ))
C
            Z_S = Z_XOO*XSN + Z_YOO*YSN
C
            SRLX = 1.00*RLX
            DS = -SRLX*REZ/Z_S
C
            X(IO,JO) = X(IO,JO) + DS*XSN
            Y(IO,JO) = Y(IO,JO) + DS*YSN
          ENDIF
C
          IF(KBCOUT.EQ.1) THEN
C---------- Neumann BC is specified on i=II plane:  relax node positions
C
            IO = II
            IM = IO-1
            IL = IO-2
C
            XOO = X(IO,JO)
            XMO = X(IM,JO)
            XLO = X(IL,JO)
C
            YOO = Y(IO,JO)
            YMO = Y(IM,JO)
            YLO = Y(IL,JO)
C
            XS = X(IO,JP) - X(IO,JM)
            YS = Y(IO,JP) - Y(IO,JM)
C
            XSN = XS / SQRT(XS**2 + YS**2)
            YSN = YS / SQRT(XS**2 + YS**2)
C
            DXIM = XPOS(IO) - XPOS(IM)
            DXIL = XPOS(IM) - XPOS(IL)
C
            DXDX1 = ( XOO - XMO ) / DXIM
            DXDX2 = ( XMO - XLO ) / DXIL
C
            DYDX1 = ( YOO - YMO ) / DXIM
            DYDX2 = ( YMO - YLO ) / DXIL
C
C---------- 2nd-order 3-point difference for tangential velocity
            REZ   = XSN*(DXDX1 + DXIM*(DXDX1 - DXDX2)/(DXIM+DXIL))
     &            + YSN*(DYDX1 + DXIM*(DYDX1 - DYDX2)/(DXIM+DXIL))
C
            Z_XOO = XSN*(1.0/DXIM + 1.0/(DXIM+DXIL))
            Z_YOO = YSN*(1.0/DXIM + 1.0/(DXIM+DXIL))
C
            Z_S = Z_XOO*XSN + Z_YOO*YSN
C
            SRLX = 1.00*RLX
            DS = -SRLX*REZ/Z_S
C
            X(IO,JO) = X(IO,JO) + DS*XSN
            Y(IO,JO) = Y(IO,JO) + DS*YSN
          ENDIF

C-------- relax all points on this streamline by SLOR
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
            YAP = 0.5*(YOO + YOP)
            YAM = 0.5*(YOO + YOM)
            YAA = 0.5*(YAP + YAM)
C
            DXIM = XPOS(IO) - XPOS(IM)
            DXIP = XPOS(IP) - XPOS(IO)
            DXIAV = 0.5*(DXIM+DXIP)
C
            DETM = YPOS(JO) - YPOS(JM)
            DETP = YPOS(JP) - YPOS(JO)
            DETAV = 0.5*(DETM+DETP)
C
c----------------
c            DXDETP = (XOP - XOO)/DETP * YAP/YAA
c            DYDETP = (YOP - YOO)/DETP * YAP/YAA
c            ALFP = DXDETP**2 + DYDETP**2
cC
c            DXDETM = (XOO - XOM)/DETM * YAM/YAA
c            DYDETM = (YOO - YOM)/DETM * YAM/YAA
c            ALFM = DXDETM**2 + DYDETM**2
cC
c            DXDET = 0.5*(DXDETP + DXDETM)
c            DYDET = 0.5*(DYDETP + DYDETM)
c            ALF = 0.5*(ALFP + ALFM)
cC
cC
c            DXDXIP = (XPO - XOO)/DXIP
c            DYDXIP = (YPO - YOO)/DXIP
c            GAMP = DXDXIP**2 + DYDXIP**2
cC
c            DXDXIM = (XOO - XMO)/DXIM
c            DYDXIM = (YOO - YMO)/DXIM
c            GAMM = DXDXIM**2 + DYDXIM**2
cC
c            DXDXI = 0.5*(DXDXIP + DXDXIM)
c            DYDXI = 0.5*(DYDXIP + DYDXIM)
c            GAM = 0.5*(GAMP + GAMM) 
cC
cC
c            BET = DXDET*DXDXI + DYDET*DYDXI
C
C----------------
c            DXDET = 0.5*( (XOP - XOO)/DETP * YAP/YAA
c     &                  + (XOO - XOM)/DETM * YAM/YAA )
c            DYDET = 0.5*( (YOP - YOO)/DETP * YAP/YAA
c     &                  + (YOO - YOM)/DETM * YAM/YAA )
c            DXDXI = 0.5*( (XPO - XOO)/DXIP
c     &                  + (XOO - XMO)/DXIM )
c            DYDXI = 0.5*( (YPO - YOO)/DXIP
c     &                  + (YOO - YMO)/DXIM )

C
            DXDET = 0.5*(XOP - XOM) / DETAV
            DYDET = 0.5*(YOP - YOM) / DETAV
            DXDXI = 0.5*(XPO - XMO) / DXIAV
            DYDXI = 0.5*(YPO - YMO) / DXIAV
C
c            DXDET2 = 0.5*( (XOP - XOO) * YAP/YAA
c     &                   + (XOO - XOM) * YAM/YAA ) / DETAV
c            DYDET2 = 0.5*( (YOP - YOO) * YAP/YAA
c     &                   + (YOO - YOM) * YAM/YAA ) / DETAV
c            ALF = DXDET2**2 + DYDET2**2 
c            GAM = DXDXI**2 + DYDXI**2
c            BET = DXDET2*DXDXI + DYDET2*DYDXI
c            AJA = DYDET2*DXDXI - DXDET2*DYDXI
C
            ALF = DXDET**2 + DYDET**2
            GAM = DXDXI**2 + DYDXI**2
            BET = DXDET*DXDXI + DYDET*DYDXI
            AJA = DYDET*DXDXI - DXDET*DYDXI


C----------------
C
c            CXIM = 1.0/(DXIM*DXIAV)
c            CXIP = 1.0/(DXIP*DXIAV)
c            CETM = 1.0/(DETM*DETAV) * YAM/YAA
c            CETP = 1.0/(DETP*DETAV) * YAP/YAA
C
            CXIM = DETM*DETP/(DXIM*DXIAV)
            CXIP = DETM*DETP/(DXIP*DXIAV)
            CETM =      DETP/      DETAV * YAM/YAA
            CETP = DETM     /      DETAV * YAP/YAA
C
            B =      -ALF*CXIM
            A = ALF*(CXIM+CXIP) + GAM*(CETM + CETP)
            C(IO) =  -ALF*CXIP
            IF(IO.EQ.2) B = 0.0
C
            XDSE = DETM*DETP*(XPP - XMP - XPM + XMM) / (4.0*DXIAV*DETAV)
            YDSE = DETM*DETP*(YPP - YMP - YPM + YMM) / (4.0*DXIAV*DETAV)
C
            D(1,IO) = ALF*((XPO-XOO)*CXIP - (XOO-XMO)*CXIM)
     &              - 2.0*BET*XDSE  !!! (XPP-XMP-XPM+XMM) / (4.0*DXIAV*DETAV)
     &              + GAM*((XOP-XOO)*CETP - (XOO-XOM)*CETM)
     &              - BET*DYDXI*DXDET*DETM*DETP/YAA
C
            D(2,IO) = ALF*((YPO-YOO)*CXIP - (YOO-YMO)*CXIM)
     &              - 2.0*BET*YDSE  !!! (YPP-YMP-YPM+YMM) / (4.0*DXIAV*DETAV)
     &              + GAM*((YOP-YOO)*CETP - (YOO-YOM)*CETM)
     &              - BET*DYDXI*DYDET*DETM*DETP/YAA
C
            AINV = 1.0/(A - B*C(IM))
            C(IO) = C(IO) * AINV
            D(1,IO) = ( D(1,IO) - B*D(1,IM) ) * AINV
            D(2,IO) = ( D(2,IO) - B*D(2,IM) ) * AINV
C
 10       CONTINUE
C
          D(1,II) = 0.
          D(2,II) = 0.
C
          IFIN = II-1
          DO 20 IBACK=2, IFIN
            IO = II-IBACK+1
            IP = IO+1
            D(1,IO) = D(1,IO) - C(IO)*D(1,IP)
            D(2,IO) = D(2,IO) - C(IO)*D(2,IP)
C
            X(IO,JO) = X(IO,JO) + RLX*D(1,IO)
            Y(IO,JO) = Y(IO,JO) + RLX*D(2,IO)
            AD1 = ABS(D(1,IO))
            AD2 = ABS(D(2,IO))
            DMAX = MAX(DMAX,AD1,AD2)
 20       CONTINUE
C
 100    CONTINUE
C
c        IF(MOD(IPASS,10).EQ.0) THEN
c          WRITE(*,*) IPASS, '  Dmax = ', DMAX, RLX
c        ENDIF
C
        IF(DMAX.LT.TOLER*DXY) RETURN
C
        RLX = RLX1
        IF(DMAX.LT.DSET1*DXY) RLX = RLX2
        IF(DMAX.LT.DSET2*DXY) RLX = RLX3
        IF(DMAX.LT.DSET3*DXY) RETURN
C
 1000 CONTINUE
C
      RETURN
      END ! AXELL

