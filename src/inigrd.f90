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


SUBROUTINE INIGRD
    !------------------------------------------------------
    !     Sets up vortex wake grid system from rotor and wake
    !     geometry
    !------------------------------------------------------
    INCLUDE 'DFDC.INC'
    DIMENSION S(IPX)
    !
    IF(LDBG) THEN
        WRITE(*, *) ' '
        WRITE(*, *) 'INIGRD setup grid for slipstream'
        !
        WRITE(LUNDBG, *) ' '
        WRITE(LUNDBG, *) 'INIGRD setup grid for slipstream'
    ENDIF
    !
    !---- Find most-upstream rotor point on duct wall and set as start of grid
    NRUPSTRM = 1
    IPUP = IPROTDW(NRUPSTRM)
    DO NR = 2, NROTOR
        IF(IPROTDW(NR).LT.IPUP) THEN
            NRUPSTRM = NR
            IPUP = IPROTDW(NR)
        ENDIF
    END DO
    !---- Setup rotor at upstream end of grid
    DO NR = 1, NROTOR
        IGROTOR(NR) = 1 + IPROTDW(NR) - IPUP
    END DO
    !
    NR = NRUPSTRM
    !---- Indices of rotor line on CB and duct define upstream grid boundary
    !---- find axial location of rotor and TE's on CB and duct walls
    IEL = 1
    IP1CB = IPFRST(IEL)
    IP2CB = IPLAST(IEL)
    IPCB = IPROTCB(NR)
    IF(IPCB.EQ.0) THEN
        WRITE(*, *) 'Error locating rotor on CB wall'
        STOP
    ENDIF
    IF(LDBG) WRITE(LUNDBG, *) 'Rotor on CB @ ', IPCB
    XCB = XP(IPCB)
    YCB = YP(IPCB)
    XTECB = XP(IP1CB)
    !c    IF(LDBG) WRITE(*,*) 'XR, XTECB ',XCB,XTECB
    !
    IEL = 2
    IP1DW = IPFRST(IEL)
    IP2DW = IPLAST(IEL)
    IPDW = IPROTDW(NR)
    IF(IPDW.EQ.0) THEN
        WRITE(*, *) 'Error locating rotor on Duct wall'
        STOP
    ENDIF
    IF(LDBG) WRITE(LUNDBG, *) 'Rotor on Duct wall @ ', IPDW
    XDW = XP(IPDW)
    YDW = YP(IPDW)
    XTEDW = XP(IP2DW)
    !c    IF(LDBG) WRITE(*,*) 'XR XTEDW ',XDW,XTEDW
    !
    !
    !---- Define lower streamline (J=1)
    I = IGROTOR(NR) - 1
    J = 1
    DO IP = IPCB, IP1CB, -1
        I = I + 1
        XG(I, J) = XP(IP)
        YG(I, J) = YP(IP)
        IP2IG(IP) = I
    END DO
    XLAST = XG(I, J)
    YLAST = YG(I, J)
    IGTECB = I
    !c      write(*,*) 'IGTECB,II defined as ',IGTECB,II
    !c    IF(LDBG) WRITE(*,*) 'L IGTECB,XYLAST ',IGTECB,XLAST,YLAST
    !
    !---- Add points at X values of duct wall if it extends downstream of CB
    IPD = IPDW + IGTECB - 1
    IF(IPD.LT.IP2DW) THEN
        DO IP = IPD + 1, IP2DW
            I = I + 1
            XG(I, J) = XP(IP)
            YG(I, J) = YLAST
        END DO
        XLAST = XG(I, J)
        YLAST = YG(I, J)
        ITE = I
        !c     IF(LDBG) WRITE(*,*) 'L ITE,XYLAST ',ITE,XLAST,YLAST
    ENDIF
    !
    !---- Check downstream wake end location
    IF(LDBG) WRITE(LUNDBG, *) 'INIGRD XWAKE ', XWAKE
    IF(XDWKLEN.LT.0.0 .OR. XDWKLEN.GT.5.0) THEN
        WRITE(*, *) 'Rotor wake length out of bounds ', XDWKLEN
        WRITE(*, *) '  CB   TE at X,Y =', XBTE(1), YBTE(1)
        WRITE(*, *) '  Duct TE at X,Y =', XBTE(2), YBTE(2)
        WRITE(*, *) '  Rotor Diameter =', 2.0 * RTIP(1)
        CALL ASKR('Enter length of wake from TE dXlen/D^', XDWKLEN)
    ENDIF
    XWAKE = MAX(XBTE(1), XBTE(2)) + XDWKLEN * 2.0 * RTIP(1)
    !
    !---- Add points downstream to end of wake
    DX = XG(I, J) - XG(I - 1, J)
    DY = YG(I, J) - YG(I - 1, J)
    DS = SQRT(DX * DX + DY * DY)
    !
    SMAX = XWAKE - XLAST
    DXW = SMAX / FLOAT(NWAKE - 1)
    DS1 = 0.2 * DXW
    CALL SETEXP(S, DS1, SMAX, NWAKE)
    DO K = 2, NWAKE
        I = I + 1
        XG(I, J) = XLAST + S(K)
        YG(I, J) = YLAST
    END DO
    ILWR = I
    !
    !
    !---- Define upper streamline (J=JJ=NRP from rotor line points)
    JJ = NRP
    !
    I = IGROTOR(NR) - 1
    J = JJ
    DO IP = IPDW, IP2DW
        I = I + 1
        XG(I, J) = XP(IP)
        YG(I, J) = YP(IP)
        IP2IG(IP) = I
    END DO
    XLAST = XG(I, J)
    YLAST = YG(I, J)
    IGTEDW = I
    !c      write(*,*) 'IGTEDW,II defined as ',IGTEDW,II
    !c    IF(LDBG) WRITE(*,*) 'U IGTEDW,XYLAST ',IGTEDW,XLAST,YLAST
    !
    !---- Add points at X values of CB wall if it extends downstream of duct
    IPC = IPCB - IGTEDW + 1
    IF(IPC.GT.IP1CB) THEN
        DO IP = IPC - 1, IP1CB, -1
            I = I + 1
            XG(I, J) = XP(IP)
            YG(I, J) = YLAST
        END DO
        XLAST = XG(I, J)
        YLAST = YG(I, J)
        ITE = I
        !c     IF(LDBG) WRITE(*,*) 'U ITE,XYLAST ',ITE,XLAST,YLAST
    ENDIF
    !
    !---- Add points downstream to end of wake
    DX = XG(I, J) - XG(I - 1, J)
    DY = YG(I, J) - YG(I - 1, J)
    DS = SQRT(DX * DX + DY * DY)
    !
    SMAX = XWAKE - XLAST
    DXW = SMAX / FLOAT(NWAKE - 1)
    DS1 = 0.2 * DXW
    CALL SETEXP(S, DS1, SMAX, NWAKE)
    DO K = 2, NWAKE
        I = I + 1
        XG(I, J) = XLAST + S(K)
        YG(I, J) = YLAST
    END DO
    IUPR = I
    !
    IF(ILWR.NE.IUPR) THEN
        WRITE(*, *) 'INIGRD: incompatible # points on J=1 and J=JJ', &
                ILWR, IUPR
        WRITE(*, *) 'IP1CB,IP2CB,IPCB ', IP1CB, IP2CB, IPCB
        WRITE(*, *) 'IP1DW,IP2DW,IPDW ', IP1DW, IP2DW, IPDW
        STOP
    ENDIF
    !---- Set total # of streamwise points
    II = IUPR
    !
    !
    !---- Add points at inlet and outlet boundaries
    !---- Inlet points are simply the rotor line points XR,YR
    I = 1
    XGMIN = XRP(1, NR)
    DO J = 1, JJ
        XG(I, J) = XRP(J, NR)
        YG(I, J) = YRP(J, NR)
        XGMIN = MIN(XGMIN, XG(I, J))
    END DO
    !
    !---- Sweep downstream setting intermediate J lines preserving relative
    !     massflow between upper and lower walls
    AMASS1 = PI * (YG(1, JJ)**2 - YG(1, 1)**2)
    DO I = 2, II
        XL = XG(I, 1)
        YL = YG(I, 1)
        XU = XG(I, JJ)
        YU = YG(I, JJ)
        AMASS = PI * (YU**2 - YL**2)
        !---- Create radial spacing by preserving relative area
        DO J = 2, JJ - 1
            DMASS1 = PI * (YG(1, J)**2 - YG(1, J - 1)**2)
            FRC = DMASS1 / AMASS1
            YY = SQRT(FRC * AMASS / PI + YG(I, J - 1)**2)
            FRAC = (YY - YL) / (YU - YL)
            XG(I, J) = XL + FRAC * (XU - XL)
            YG(I, J) = YY
        END DO
    END DO
    !
    !---- Done, initial grid defined!
    IF(LDBG) THEN
        WRITE(*, *) 'INIGRD grid defined with II,JJ ', II, JJ
        WRITE(LUNDBG, *) 'INIGRD grid defined with II,JJ ', II, JJ
    ENDIF
    !
    !     Add initial grid smoothing relaxation?
    !---- Grid smoothing by elliptic solver
    CALL RLXGRD
    !
    !---- invalidate existing solution to incorporate new vortex wakes
    LNCVP = .FALSE.
    LQAIC = .FALSE.
    LQGIC = .FALSE.
    LQCNT = .FALSE.
    LGSYS = .FALSE.
    LGAMU = .FALSE.
    LGAMA = .FALSE.
    !
    RETURN
END
! INIGRD


SUBROUTINE UPDGRD
    !------------------------------------------------------
    !     Updates grid upper and lower boundaries with
    !     points from current upper and lower vortex wakes
    !------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    IF(LDBG) THEN
        WRITE(*, *) ' '
        WRITE(*, *) 'Updating grid system boundaries'
    ENDIF
    !
    !---- Redefine lower streamline (J=1) with IR=1 wake points
    J = 1
    IEL = IR2IEL(1)
    IP1 = IPFRST(IEL)
    IP2 = IPLAST(IEL)
    DO IP = IP1, IP2
        I = IGTECB + (IP - IP1)
        XG(I, J) = XP(IP)
        YG(I, J) = YP(IP)
    END DO
    IF(LDBG) THEN
        WRITE(LUNDBG, *) 'Updated lower boundary from element ', IEL
    ENDIF
    !
    !---- Redefine upper streamline (J=JJ) with IR=NRP wake points
    J = JJ
    IEL = IR2IEL(NRP)
    IP1 = IPFRST(IEL)
    IP2 = IPLAST(IEL)
    DO IP = IP1, IP2
        I = IGTEDW + (IP - IP1)
        XG(I, J) = XP(IP)
        YG(I, J) = YP(IP)
    END DO
    IF(LDBG) THEN
        WRITE(LUNDBG, *) 'Updated upper boundary from element ', IEL
    ENDIF
    !
    !---- Grid boundaries redefined!
    IF(LDBG) WRITE(*, *) 'UPDGRD grid boundaries updated'
    !
    !---- Grid smoothing by elliptic solver
    CALL RLXGRD
    !
    RETURN
END
! UPDGRD




SUBROUTINE RLXGRD
    !------------------------------------------------------
    !     Runs elliptic grid solution to smooth grid
    !------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    IF(II.LE.0 .OR. JJ.LE.0) THEN
        WRITE(*, *) 'RLXGRD:  no grid defined II,JJ ', II, JJ
        RETURN
    ENDIF
    !
    !---- Set up streamfunction array
    YPOS(1) = 0.0
    DO J = 2, JJ
        YPOS(J) = YG(1, J)**2 - YG(1, J - 1)**2 + YPOS(J - 1)
    END DO
    !---- Set up streamwise spacing array
    DO I = 1, II
        XPOS(I) = XG(I, 1)
    END DO
    !
    !---- Run grid solver
    ITMAXS = -400
    TOLER = 1.0E-9
    KBCINL = 0   ! inlet plane Dirichlet, fix nodes
    KBCOUT = 1   ! outlet plane Neumann, move nodes on boundary
    KBCBOT = 0   ! bottom streamline Dirichlet, fix nodes
    KBCTOP = 0   ! top streamline Dirichlet, fix nodes
    !      WRITE(*,*)
    !      WRITE(*,*) 'Grid relaxation iteration'
    CALL AXELL(IX, II, JJ, XG, YG, XPOS, YPOS, ITMAXS, TOLER, &
            KBCINL, KBCOUT, KBCBOT, KBCTOP)
    !
    RETURN
END
! RLXGRD

!     &  XG(IX,JX),    YG(IX,JX),
!     &  QG(IX,JX),    QXG(IX,JX),   QYG(IX,JX),  QTG(IX,JX),  
!     &  RG(IX,JX),    PG(IX,JX),    POG(IX,JX),
!     &  BGAMG(IX,JX), DSG(IX,JX),   DHG(IX,JX),
!     &  XPOS(IX),     YPOS(JX)


SUBROUTINE SETGRDFLW
    !---------------------------------------------------------
    !     Sets grid flow data from circulation and entropy on
    !     rotor lines
    !---------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    !---- Clear grid data for accumulation over sources
    CALL CLRGRDFLW
    !
    !---- Set B*GAM circulation, delta entropy, and delta enthalpy on gridC
    !     from rotor and drag sources

    DO N = 1, NROTOR
        !
        IEL = IELROTOR(N)
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        !---- Set values at rotor I station due to blade row
        IG = IGROTOR(N)
        !
        NNC = IC2 - IC1 + 1
        IF(NNC.NE.JJ - 1) THEN
            WRITE(*, *) 'Rotor source element-grid mismatch ', NNC, JJ - 1
            STOP
        ENDIF
        !
        DO J = 1, JJ - 1
            IR = J
            !---- Circulation and enthalpy
            DBGAMG = BGAM(IR, N)
            DDHG = OMEGA(N) * BGAM(IR, N) * PI2I
            !---- Entropy from drag sources on rotor
            IC = IC1 + IR - 1
            IP1 = IPCO(IC)
            IP2 = IPCP(IC)
            VMC = SQRT(QC(1, IC)**2 + QC(2, IC)**2)
            SIGC = 0.5 * (SIG(IP1) + SIG(IP2))
            DDSG = VMC * SIGC
            !
            !---- Convect values downstream from rotor
            DO I = IG, II - 1
                BGAMG(I, J) = BGAMG(I, J) + DBGAMG
                DHG(I, J) = DHG(I, J) + DDHG
                DSG(I, J) = DSG(I, J) + DDSG
                !c        RG(I,J)    = RG(1,J)
            END DO
            !
        END DO
    END DO
    !
    RETURN
END
! SETGRDFLW


SUBROUTINE CLRGRDFLW
    !---------------------------------------------------------
    !     Clears (initializes to zero) grid flow BGAMG,DHG,DSG
    !     This is normally done before accumulating disk
    !     contributions
    !---------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    !---- Clear grid data for accumulation of circulation, enthalpy, entropy
    DO J = 1, JJ
        DO I = 1, II
            BGAMG(I, J) = 0.0
            DHG(I, J) = 0.0
            DSG(I, J) = 0.0
            RG(I, J) = RHO
        END DO
    END DO
    !
    RETURN
END
! CLRGRDFLW


SUBROUTINE ROTBG2GRD(N)
    !---------------------------------------------------------
    !     Updates grid flow (circulation) from BGAM on rotor N
    !     rotor lines
    !---------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    !---- Set B*GAM circulation, delta enthalpy on grid from blade/rotors
    !---- Set values at grid IG station due to blade row
    IG = IGROTOR(N)
    !
    DO J = 1, JJ - 1
        IR = J
        !---- Circulation and enthalpy
        DBGAMG = BGAM(IR, N)
        DDHG = OMEGA(N) * BGAM(IR, N) * PI2I
        !
        !---- Convect values downstream from rotor
        DO I = IG, II - 1
            BGAMG(I, J) = BGAMG(I, J) + DBGAMG
            DHG(I, J) = DHG(I, J) + DDHG
        END DO
        !
    END DO
    !
    RETURN
END
! ROTBG2GRD



SUBROUTINE GETGRDFLW
    !---------------------------------------------------------
    !     Sets grid initial flow data from rotor line conditions
    !---------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    CALL UVGRDC(IX, II, JJ, XG, YG, XPOS, YPOS, QXG, QYG)
    !---- Print grid velocities
    WRITE(44, *) 'I,J, XC,YC, U, V'
    DO I = 1, II - 1
        IO = I
        IP = I + 1
        DO J = 1, JJ - 1
            JO = J
            JP = J + 1
            XOO = XG(IO, JO)
            XPO = XG(IP, JO)
            XOP = XG(IO, JP)
            XPP = XG(IP, JP)
            YOO = YG(IO, JO)
            YPO = YG(IP, JO)
            YOP = YG(IO, JP)
            YPP = YG(IP, JP)
            XCG = 0.25 * (XOO + XPO + XOP + XPP)
            YCG = 0.25 * (YOO + YPO + YOP + YPP)
            WRITE(44, 100) I, J, XCG, YCG, QXG(I, J), QYG(I, J)
        END DO
    END DO
    100  FORMAT(2I5, 5(1X, F11.6))
    !
    RETURN
END
! GETGRDFLW



SUBROUTINE INLSFCN
    !---------------------------------------------------------
    !     Calculates streamfunction (YPOS) at grid inlet using
    !     velocities at rotor and current grid points
    !---------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    NR = NRUPSTRM
    !---- Set up streamfunction array using grid and rotor line velocity
    YPOS(1) = 0.0
    DO J = 1, JJ
        VX = VABS(1, J - 1, NR)
        YPOS(J) = 0.5 * VX * (YG(1, J)**2 - YG(1, J - 1)**2) + YPOS(J - 1)
        !c        YPOS(J) = 0.5*VX*YG(1,J)**2
    END DO
    YSCL = 1.0 / YPOS(JJ)
    !c      DO J = 1, JJ
    !c        write(*,*) 'j ypos ',j,ypos(j)
    !c        YPOS(J) = YSCL*YPOS(J)
    !c      END DO
    !
    RETURN
END
! INLSFCN


SUBROUTINE RLXGRD2
    !------------------------------------------------------
    !     Runs elliptic grid solution to smooth grid
    !------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    IF(II.LE.0 .OR. JJ.LE.0) THEN
        WRITE(*, *) 'RLXGRD:  no grid defined II,JJ ', II, JJ
        RETURN
    ENDIF
    !
    !---- Set up streamfunction array
    CALL INLSFCN
    !
    !---- Set up streamwise spacing array
    DO I = 1, II
        XPOS(I) = XG(I, 1)
    END DO
    !
    !---- Run grid solver
    ITMAXS = -400
    TOLER = 1.0E-9
    KBCINL = 0   ! inlet plane Dirichlet, fix nodes
    KBCOUT = 1   ! outlet plane Neumann, move nodes on boundary
    KBCBOT = 0   ! bottom streamline Dirichlet, fix nodes
    KBCTOP = 0   ! top streamline Dirichlet, fix nodes
    CALL AXELL(IX, II, JJ, XG, YG, XPOS, YPOS, ITMAXS, TOLER, &
            KBCINL, KBCOUT, KBCBOT, KBCTOP)
    !
    RETURN
END
! RLXGRD2


SUBROUTINE AXELL(IX, II, JJ, X, Y, XPOS, YPOS, ITMAXS, TOLER, &
        KBCINL, KBCOUT, KBCBOT, KBCTOP)
    IMPLICIT REAL (A-H, M, O-Z)
    DIMENSION X(IX, *), Y(IX, *)
    DIMENSION XPOS(*), YPOS(*)
    !-------------------------------------------------------------
    !     Axisymmetric elliptic grid smoother.
    !     Uses SLOR implicit along i direction.
    !
    !   Input:
    !     IX        first grid-array dimension
    !     II,JJ     grid i,j size
    !     X(i,j)    grid z coordinates
    !     Y(i,j)    grid r coordinates
    !     XPOS(i)   grid s values
    !     YPOS(j)   grid e values  (axisymmetric streamfunction)
    !                ( YPOS(j) = (j-1)**2 gives uniform spacing in r )
    !     ITMAXS    max number of SLOR smoothing passes
    !               if ITMAXS < 0, overrelaxation is NOT used
    !     TOLER     convergence tolerance (movement distance)
    !     KBCINL    1 = Neumann, 0 = Dirichlet   BC applied at i=1
    !     KBCOUT    1 = Neumann, 0 = Dirichlet   BC applied at i=II
    !     KBCBOT    1 = Neumann, 0 = Dirichlet   BC applied at j=1
    !     KBCTOP    1 = Neumann, 0 = Dirichlet   BC applied at j=JJ
    !
    !   Output:
    !     X(i,j)    smoothed-grid z coordinates
    !     Y(i,j)    smoothed-grid r coordinates
    !
    !.............................................................
    !
    !   Solves Thompson's grid-generation equations with source term
    !   which makes e(x,r) satisfy the axisymmetric streamfunction
    !   equation.  Hence,  e = constant  lines are streamlines.
    !
    !       s_zz + s_rr = 0
    !       e_zz + e_rr = e_r / r
    !
    !   These equations for s(z,r), e(z,r) are inverted to
    !
    !       a z_ss  -  2b z_se  +  c z_ee  =  (J/r) z_e z_s
    !
    !       a r_ss  -  2b r_se  +  c r_ee  =  (J/r) r_e z_s
    !
    !   where
    !                    2        2                    2        2
    !         a  =  (z_e)  + (r_e)           c  = (z_s)  + (r_s)
    !
    !         b  =  z_e z_s  +  r_e r_s
    !
    !         J  =  z_s r_e  -  z_e r_s
    !
    !
    !   The equivalent form which is actually differenced is
    !
    !       a z_ss  -  2b z_se  +  c/r (r z_e)_e  -  b/r r_s z_e  =  0
    !
    !       a r_ss  -  2b r_se  +  c/r (r r_e)_e  -  b/r r_s r_e  =  0
    !
    !-------------------------------------------------------------
    PARAMETER (IDIM = 500)
    DIMENSION C(IDIM), D(2, IDIM)
    !
    DIMENSION XT(IDIM), YT(IDIM), ST(IDIM), XTS(IDIM), YTS(IDIM), STI(IDIM)
    DIMENSION XB(IDIM), YB(IDIM), SB(IDIM), XBS(IDIM), YBS(IDIM), SBI(IDIM)
    !
    IF(II.GT.IDIM) STOP 'ELLIP: Array overflow.  Increase IDIM.'
    !
    ITMAX = IABS(ITMAXS)
    !
    !---- convergence tolerance for each phase
    DSET1 = 1.0E-1
    DSET2 = 5.0E-3
    DSET3 = 5.0E-7
    !
    IF(ITMAXS.GT.0) THEN
        !----- over-relaxation parameter for each phase
        RLX1 = 1.00       !          DMAX > DSET1
        RLX2 = 1.10       !  DSET1 > DMAX > DSET2
        RLX3 = 1.40       !  DSET2 > DMAX > DSET3
        !CC    STOP              !  DSET3 > DMAX
    ELSE
        RLX1 = 1.0
        RLX2 = 1.0
        RLX3 = 1.0
    ENDIF
    !
    !---- spline bottom & top wall shapes
    JO = 1
    DO IO = 1, II
        XB(IO) = X(IO, JO)
        YB(IO) = Y(IO, JO)
    ENDDO
    !
    DX = ABS(XB(II) - XB(1))
    DY = ABS(YB(II) - YB(1))
    DXY = MAX(XB(1), YB(1), DX, DY)
    !
    CALL SCALC(XB, YB, SB, II)
    CALL SPLINE(XB, XBS, SB, II)
    CALL SPLINE(YB, YBS, SB, II)
    DO IO = 1, II
        SBI(IO) = SB(IO)
    ENDDO
    !
    JO = JJ
    DO IO = 1, II
        XT(IO) = X(IO, JO)
        YT(IO) = Y(IO, JO)
    ENDDO
    CALL SCALC(XT, YT, ST, II)
    CALL SPLINE(XT, XTS, ST, II)
    CALL SPLINE(YT, YTS, ST, II)
    DO IO = 1, II
        STI(IO) = ST(IO)
    ENDDO
    !
    RLX = RLX1
    !
    DO 1000 IPASS = 1, ITMAX
        !
        DMAX = 0.
        !
        JO = 1
        IF(KBCBOT.EQ.1 .AND. JJ.GT.2) THEN
            !-------- relax to get dX/dj = 0 on lower streamline
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
                !
                YOO = Y(IO, JO)
                YOP = Y(IO, JP)
                YOQ = Y(IO, JQ)
                !
                XS = DEVAL(SBI(IO), XB, XBS, SB, II)
                YS = DEVAL(SBI(IO), YB, YBS, SB, II)
                !
                DETP = YPOS(JP) - YPOS(JO)
                DETQ = YPOS(JQ) - YPOS(JP)
                !
                DXDE1 = (XOP - XOO) / DETP
                DXDE2 = (XOQ - XOP) / DETQ
                !
                DYDE1 = (YOP - YOO) / DETP
                DYDE2 = (YOQ - YOP) / DETQ
                !
                REZ = XS * (DXDE1 - DETP * (DXDE2 - DXDE1) / (DETP + DETQ))&
                        + YS * (DYDE1 - DETP * (DYDE2 - DYDE1) / (DETP + DETQ))
                !
                Z_XOO = XS * (-1.0 / DETP - 1.0 / (DETP + DETQ))
                Z_YOO = YS * (-1.0 / DETP - 1.0 / (DETP + DETQ))
                !
                Z_S = Z_XOO * XS + Z_YOO * YS
                !
                DS = -REZ / Z_S
                !
                SBI(IO) = SBI(IO) + DS
                !
                X(IO, JO) = SEVAL(SBI(IO), XB, XBS, SB, II)
                Y(IO, JO) = SEVAL(SBI(IO), YB, YBS, SB, II)
            ENDDO
        ENDIF
        !
        JO = JJ
        IF(KBCTOP.EQ.1 .AND. JJ.GT.2) THEN
            !-------- relax to get dX/dj = 0 on upper streamline
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
                !
                YOO = Y(IO, JO)
                YOM = Y(IO, JM)
                YOL = Y(IO, JL)
                !
                XS = DEVAL(STI(IO), XT, XTS, ST, II)
                YS = DEVAL(STI(IO), YT, YTS, ST, II)
                !
                DETM = YPOS(JO) - YPOS(JM)
                DETL = YPOS(JM) - YPOS(JL)
                !
                DXDE1 = (XOO - XOM) / DETM
                DXDE2 = (XOM - XOL) / DETL
                !
                DYDE1 = (YOO - YOM) / DETM
                DYDE2 = (YOM - YOL) / DETL
                !
                REZ = XS * (DXDE1 + DETM * (DXDE1 - DXDE2) / (DETM + DETL))&
                        + YS * (DYDE1 + DETM * (DYDE1 - DYDE2) / (DETM + DETL))
                !
                Z_XOO = XS * (1.0 / DETM + 1.0 / (DETM + DETL))
                Z_YOO = YS * (1.0 / DETM + 1.0 / (DETM + DETL))
                !
                Z_S = Z_XOO * XS + Z_YOO * YS
                !
                DS = -REZ / Z_S
                !
                STI(IO) = STI(IO) + DS
                !
                X(IO, JO) = SEVAL(STI(IO), XT, XTS, ST, II)
                Y(IO, JO) = SEVAL(STI(IO), YT, YTS, ST, II)
            ENDDO
        ENDIF
        !
        !
        !------ go over all interior streamlines
        DO 100 JO = 2, JJ - 1
            JM = JO - 1
            JP = JO + 1
            !
            IF(KBCINL.EQ.1) THEN
                !---------- Neumann BC is specified on i=1 plane:  relax node positions
                !
                IO = 1
                IP = IO + 1
                IQ = IO + 2
                !
                XOO = X(IO, JO)
                XPO = X(IP, JO)
                XQO = X(IQ, JO)
                !
                YOO = Y(IO, JO)
                YPO = Y(IP, JO)
                YQO = Y(IQ, JO)
                !
                XS = X(IO, JP) - X(IO, JM)
                YS = Y(IO, JP) - Y(IO, JM)
                !
                XSN = XS / SQRT(XS**2 + YS**2)
                YSN = YS / SQRT(XS**2 + YS**2)
                !
                DXIP = XPOS(IP) - XPOS(IO)
                DXIQ = XPOS(IQ) - XPOS(IP)
                !
                DXDX1 = (XPO - XOO) / DXIP
                DXDX2 = (XQO - XPO) / DXIQ
                !
                DYDX1 = (YPO - YOO) / DXIP
                DYDX2 = (YQO - YPO) / DXIQ
                !
                !---------- 2nd-order 3-point difference for tangential velocity
                REZ = XSN * (DXDX1 - DXIP * (DXDX2 - DXDX1) / (DXIP + DXIQ))&
                        + YSN * (DYDX1 - DXIP * (DYDX2 - DYDX1) / (DXIP + DXIQ))
                !
                Z_XOO = XSN * (-1.0 / DXIP - 1.0 / (DXIP + DXIQ))
                Z_YOO = YSN * (-1.0 / DXIP - 1.0 / (DXIP + DXIQ))
                !
                Z_S = Z_XOO * XSN + Z_YOO * YSN
                !
                SRLX = 1.00 * RLX
                DS = -SRLX * REZ / Z_S
                !
                X(IO, JO) = X(IO, JO) + DS * XSN
                Y(IO, JO) = Y(IO, JO) + DS * YSN
            ENDIF
            !
            IF(KBCOUT.EQ.1) THEN
                !---------- Neumann BC is specified on i=II plane:  relax node positions
                !
                IO = II
                IM = IO - 1
                IL = IO - 2
                !
                XOO = X(IO, JO)
                XMO = X(IM, JO)
                XLO = X(IL, JO)
                !
                YOO = Y(IO, JO)
                YMO = Y(IM, JO)
                YLO = Y(IL, JO)
                !
                XS = X(IO, JP) - X(IO, JM)
                YS = Y(IO, JP) - Y(IO, JM)
                !
                XSN = XS / SQRT(XS**2 + YS**2)
                YSN = YS / SQRT(XS**2 + YS**2)
                !
                DXIM = XPOS(IO) - XPOS(IM)
                DXIL = XPOS(IM) - XPOS(IL)
                !
                DXDX1 = (XOO - XMO) / DXIM
                DXDX2 = (XMO - XLO) / DXIL
                !
                DYDX1 = (YOO - YMO) / DXIM
                DYDX2 = (YMO - YLO) / DXIL
                !
                !---------- 2nd-order 3-point difference for tangential velocity
                REZ = XSN * (DXDX1 + DXIM * (DXDX1 - DXDX2) / (DXIM + DXIL))&
                        + YSN * (DYDX1 + DXIM * (DYDX1 - DYDX2) / (DXIM + DXIL))
                !
                Z_XOO = XSN * (1.0 / DXIM + 1.0 / (DXIM + DXIL))
                Z_YOO = YSN * (1.0 / DXIM + 1.0 / (DXIM + DXIL))
                !
                Z_S = Z_XOO * XSN + Z_YOO * YSN
                !
                SRLX = 1.00 * RLX
                DS = -SRLX * REZ / Z_S
                !
                X(IO, JO) = X(IO, JO) + DS * XSN
                Y(IO, JO) = Y(IO, JO) + DS * YSN
            ENDIF

            !-------- relax all points on this streamline by SLOR
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
                YAP = 0.5 * (YOO + YOP)
                YAM = 0.5 * (YOO + YOM)
                YAA = 0.5 * (YAP + YAM)
                !
                DXIM = XPOS(IO) - XPOS(IM)
                DXIP = XPOS(IP) - XPOS(IO)
                DXIAV = 0.5 * (DXIM + DXIP)
                !
                DETM = YPOS(JO) - YPOS(JM)
                DETP = YPOS(JP) - YPOS(JO)
                DETAV = 0.5 * (DETM + DETP)
                !
                !----------------
                !            DXDETP = (XOP - XOO)/DETP * YAP/YAA
                !            DYDETP = (YOP - YOO)/DETP * YAP/YAA
                !            ALFP = DXDETP**2 + DYDETP**2
                !C
                !            DXDETM = (XOO - XOM)/DETM * YAM/YAA
                !            DYDETM = (YOO - YOM)/DETM * YAM/YAA
                !            ALFM = DXDETM**2 + DYDETM**2
                !C
                !            DXDET = 0.5*(DXDETP + DXDETM)
                !            DYDET = 0.5*(DYDETP + DYDETM)
                !            ALF = 0.5*(ALFP + ALFM)
                !C
                !C
                !            DXDXIP = (XPO - XOO)/DXIP
                !            DYDXIP = (YPO - YOO)/DXIP
                !            GAMP = DXDXIP**2 + DYDXIP**2
                !C
                !            DXDXIM = (XOO - XMO)/DXIM
                !            DYDXIM = (YOO - YMO)/DXIM
                !            GAMM = DXDXIM**2 + DYDXIM**2
                !C
                !            DXDXI = 0.5*(DXDXIP + DXDXIM)
                !            DYDXI = 0.5*(DYDXIP + DYDXIM)
                !            GAM = 0.5*(GAMP + GAMM)
                !C
                !C
                !            BET = DXDET*DXDXI + DYDET*DYDXI
                !
                !----------------
                !            DXDET = 0.5*( (XOP - XOO)/DETP * YAP/YAA
                !     &                  + (XOO - XOM)/DETM * YAM/YAA )
                !            DYDET = 0.5*( (YOP - YOO)/DETP * YAP/YAA
                !     &                  + (YOO - YOM)/DETM * YAM/YAA )
                !            DXDXI = 0.5*( (XPO - XOO)/DXIP
                !     &                  + (XOO - XMO)/DXIM )
                !            DYDXI = 0.5*( (YPO - YOO)/DXIP
                !     &                  + (YOO - YMO)/DXIM )

                !
                DXDET = 0.5 * (XOP - XOM) / DETAV
                DYDET = 0.5 * (YOP - YOM) / DETAV
                DXDXI = 0.5 * (XPO - XMO) / DXIAV
                DYDXI = 0.5 * (YPO - YMO) / DXIAV
                !
                !            DXDET2 = 0.5*( (XOP - XOO) * YAP/YAA
                !     &                   + (XOO - XOM) * YAM/YAA ) / DETAV
                !            DYDET2 = 0.5*( (YOP - YOO) * YAP/YAA
                !     &                   + (YOO - YOM) * YAM/YAA ) / DETAV
                !            ALF = DXDET2**2 + DYDET2**2
                !            GAM = DXDXI**2 + DYDXI**2
                !            BET = DXDET2*DXDXI + DYDET2*DYDXI
                !            AJA = DYDET2*DXDXI - DXDET2*DYDXI
                !
                ALF = DXDET**2 + DYDET**2
                GAM = DXDXI**2 + DYDXI**2
                BET = DXDET * DXDXI + DYDET * DYDXI
                AJA = DYDET * DXDXI - DXDET * DYDXI


                !----------------
                !
                !            CXIM = 1.0/(DXIM*DXIAV)
                !            CXIP = 1.0/(DXIP*DXIAV)
                !            CETM = 1.0/(DETM*DETAV) * YAM/YAA
                !            CETP = 1.0/(DETP*DETAV) * YAP/YAA
                !
                CXIM = DETM * DETP / (DXIM * DXIAV)
                CXIP = DETM * DETP / (DXIP * DXIAV)
                CETM = DETP / DETAV * YAM / YAA
                CETP = DETM / DETAV * YAP / YAA
                !
                B = -ALF * CXIM
                A = ALF * (CXIM + CXIP) + GAM * (CETM + CETP)
                C(IO) = -ALF * CXIP
                IF(IO.EQ.2) B = 0.0
                !
                XDSE = DETM * DETP * (XPP - XMP - XPM + XMM) / (4.0 * DXIAV * DETAV)
                YDSE = DETM * DETP * (YPP - YMP - YPM + YMM) / (4.0 * DXIAV * DETAV)
                !
                D(1, IO) = ALF * ((XPO - XOO) * CXIP - (XOO - XMO) * CXIM)&
                        - 2.0 * BET * XDSE  & !!! (XPP-XMP-XPM+XMM) / (4.0*DXIAV*DETAV)&
                + GAM*((XOP-XOO)*CETP - (XOO-XOM)*CETM)&
                        - BET*DYDXI*DXDET*DETM*DETP/YAA
                !
                D(2, IO) = ALF * ((YPO - YOO) * CXIP - (YOO - YMO) * CXIM)&
                        - 2.0 * BET * YDSE & !!! (YPP-YMP-YPM+YMM) / (4.0*DXIAV*DETAV)&
                + GAM*((YOP-YOO)*CETP - (YOO-YOM)*CETM)&
                        - BET*DYDXI*DYDET*DETM*DETP/YAA
                !
                AINV = 1.0 / (A - B * C(IM))
                C(IO) = C(IO) * AINV
                D(1, IO) = (D(1, IO) - B * D(1, IM)) * AINV
                D(2, IO) = (D(2, IO) - B * D(2, IM)) * AINV
                !
            10       CONTINUE
            !
            D(1, II) = 0.
            D(2, II) = 0.
            !
            IFIN = II - 1
            DO 20 IBACK = 2, IFIN
                IO = II - IBACK + 1
                IP = IO + 1
                D(1, IO) = D(1, IO) - C(IO) * D(1, IP)
                D(2, IO) = D(2, IO) - C(IO) * D(2, IP)
                !
                X(IO, JO) = X(IO, JO) + RLX * D(1, IO)
                Y(IO, JO) = Y(IO, JO) + RLX * D(2, IO)
                AD1 = ABS(D(1, IO))
                AD2 = ABS(D(2, IO))
                DMAX = MAX(DMAX, AD1, AD2)
            20       CONTINUE
            !
        100    CONTINUE
        !
        !        IF(MOD(IPASS,10).EQ.0) THEN
        !          WRITE(*,*) IPASS, '  Dmax = ', DMAX, RLX
        !        ENDIF
        !
        IF(DMAX.LT.TOLER * DXY) RETURN
        !
        RLX = RLX1
        IF(DMAX.LT.DSET1 * DXY) RLX = RLX2
        IF(DMAX.LT.DSET2 * DXY) RLX = RLX3
        IF(DMAX.LT.DSET3 * DXY) RETURN
        !
    1000 CONTINUE
    !
    RETURN
END
! AXELL
