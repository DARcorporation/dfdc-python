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

SUBROUTINE WAKERESET
    !-------------------------------------------------------------
    !     Realigns the upper and lower wakes bounding the duct
    !     vortex wake grid with the current flow.
    !
    !     A new wake grid system is defined by updating and re-relaxing
    !     the grid.  All wake elements are re-initialized with the new
    !     geometry.
    !-------------------------------------------------------------
    INCLUDE 'DFDC.INC'
    !
    !---- CB wake
    IEL = IR2IEL(1)
    CALL WAKMOV(IEL, ISPLOT(IEL))
    !
    !---- duct wake (check for tip gap...)
    IF(TGAP.LE.0.0) THEN
        IEL = IR2IEL(NRP)
        CALL WAKMOV(IEL, ISPLOT(IEL))
    ELSE
        !---- if tip gap move next inner wake (the one with circulation)
        IEL = IR2IEL(NRP - 1)
        CALL WAKMOV(IEL, ISPLOT(IEL))
        !---- move outer wake streamline to match moved streamline
        IP1 = IPFRST(IEL)
        IELO = IR2IEL(NRP)
        IP1O = IPFRST(IELO)
        IP2O = IPLAST(IELO)
        IPTE = IP1 + IGTEDW - 1
        DX0 = XP(IP1O) - XP(IPTE)
        DY0 = YP(IP1O) - YP(IPTE)
        !c        write(*,*) 'WAKERESET first dx0,dy0 ',dx0,dy0
        DO IP = IP1O, IP2O
            IPM = IPTE + (IP - IP1O)
            DX = XP(IP) - XP(IPM)
            DY = YP(IP) - YP(IPM)
            !c          write(*,*) 'WAKERESET ip,dx,dy ',ip,dx,dy
            YP(IP) = YP(IP) - (DY - DY0)
        END DO
    ENDIF
    !
    !---- change grid boundaries, relax grid and update vortex wakes
    CALL UPDGRD
    CALL UPDROTWAK
    !
    !---- invalidate existing solution
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
! WAKERESET




SUBROUTINE WAKMOV(IELA, ISQ)
    !-------------------------------------------------------------
    !     Moves existing wake (element IELA) by
    !     aligning panels with current velocities.
    !     Panel side flag ISQ indicates which velocity to use.
    !       ISQ=0 uses mean surface velocities QC
    !       ISQ=1 uses right surface velocities QCR
    !       ISQ=2 uses left  surface velocities QCL
    !
    !     Panel points are moved preserving existing panel lengths
    !     while aligning the panels with the current velocity
    !     vectors on the panel.
    !-------------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    IF(NETYPE(IELA).EQ.0) THEN
        !----- this element is a body panel
        WRITE(*, *) 'WAKMOV: Element', IELA, ' is body, cannot be moved!'
        RETURN
    ENDIF
    !
    IEL = IELA
    !
    !---- Save X location for end of wake
    IP1 = IPFRST(IEL)
    IP2 = IPLAST(IEL)
    XP2 = XP(IP2)
    !
    CALL WAKMOVR(IEL, XP, YP, ISQ)
    !
    !---- Reset last point to initial (unchanged) XWAKE location
    XP(IP2) = XP2
    !
    !---- respline node powgrints
    CALL XYPSPL(IEL)
    !---- reset control points and normal vectors
    CALL XYCSET(IEL)
    CALL ANPSET(IEL)
    !
    !      WRITE(*,*) 'Generating dq/d(gam,sig) influence of moved wake...'
    !      CALL QAIC1(.FALSE., 1,NCTOT, IEL,IEL)
    !      DO KEL = 1, NEL
    !        ICE = NCX + KEL
    !        IF(LBODY(KEL)) CALL QAIC1(.FALSE.,ICE,ICE,IEL,IEL)
    !      ENDDO
    !
    !      IP1 = IPFRST(IEL)
    !      IP2 = IPLAST(IEL)
    !      CALL GSYS(.FALSE., .FALSE., .TRUE., IP1,IP2, 1,0)
    !
    !---- free existing unit wake strength vectors (if any), mark vectors invalid
    DO IP = IP1, IP2
        LUSET(IUGAM(IP)) = .FALSE.
        LUSET(IUSIG(IP)) = .FALSE.
        IUGAM(IP) = 0
        IUSIG(IP) = 0
    ENDDO
    DO IR = 1, NRP
        LUSET(IUVWK(IR)) = .FALSE.
        IUVWK(IR) = 0
    ENDDO
    !
    RETURN
END
! WAKMOV



SUBROUTINE WAKMOVR(IEL, XPNEW, YPNEW, ISQ)
    !----------------------------------------------------------
    !---- Routine to move points on element to new positions
    !     to align with velocity vector on side ISQ (0/1/2)
    !----------------------------------------------------------
    INCLUDE 'DFDC.inc'
    DIMENSION XPNEW(*), YPNEW(*)
    !
    DIMENSION DSP(IPX)
    !
    !---- speed below this is treated as stagnation
    QSTAG = 0.0001 * QREF
    !
    IP1 = IPFRST(IEL)
    IP2 = IPLAST(IEL)
    !
    IC1 = ICFRST(IEL)
    IC2 = ICLAST(IEL)
    !
    !---- calculate velocity and its Jacobians at all wake points
    !cc      CALL QAIC1(.FALSE., IC1,IC2, 1,NEL)
    !---- set velocities on both sides of panels at all control points
    !cc      CALL QQCALC(NCTOT,ANC, IPCO,IPCP, GAM,SIG, QC, QCL,QCR )
    !
    DO IP = IP1, IP2 - 1
        DXP = XP(IP + 1) - XP(IP)
        DYP = YP(IP + 1) - YP(IP)
        DSP(IP) = SQRT(DXP**2 + DYP**2)
    ENDDO
    !
    !---- march down wake setting wake points
    IP = IPCO(IC1)
    XPNEW(IP) = XP(IP)
    YPNEW(IP) = YP(IP)
    !
    DO IC = IC1, IC2
        IPO = IPCO(IC)
        IPP = IPCP(IC)
        !
        IF(ISQ .EQ. 0) THEN
            QX = QC(1, IC)
            QY = QC(2, IC)
        ELSEIF(ISQ .EQ. 1) THEN
            QX = QCR(1, IC)
            QY = QCR(2, IC)
        ELSEIF(ISQ .EQ. 2) THEN
            QX = QCL(1, IC)
            QY = QCL(2, IC)
        ELSE
            QX = QC(1, IC)
            QY = QC(2, IC)
        ENDIF
        !
        QMAG = SQRT(QX**2 + QY**2)
        IF(QMAG.LT.QSTAG) THEN
            !------- stagnation... don't try to move points
            RETURN
        ENDIF
        !
        UN = QX / QMAG
        VN = QY / QMAG
        XPNEW(IPP) = XPNEW(IPO) + UN * DSP(IPO)
        YPNEW(IPP) = YPNEW(IPO) + VN * DSP(IPO)
        !
    END DO
    !
    RETURN
END
! WAKMOVR


