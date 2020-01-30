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

      SUBROUTINE WAKERESET
C-------------------------------------------------------------
C     Realigns the upper and lower wakes bounding the duct
C     vortex wake grid with the current flow.  
C
C     A new wake grid system is defined by updating and re-relaxing
C     the grid.  All wake elements are re-initialized with the new
C     geometry.
C-------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- CB wake
      IEL = IR2IEL(1)
      CALL WAKMOV(IEL,ISPLOT(IEL))
C
C---- duct wake (check for tip gap...)
      IF(TGAP.LE.0.0) THEN
        IEL = IR2IEL(NRP)
        CALL WAKMOV(IEL,ISPLOT(IEL))
      ELSE
C---- if tip gap move next inner wake (the one with circulation)
        IEL = IR2IEL(NRP-1)
        CALL WAKMOV(IEL,ISPLOT(IEL))
C---- move outer wake streamline to match moved streamline
        IP1 = IPFRST(IEL)
        IELO = IR2IEL(NRP)
        IP1O = IPFRST(IELO)
        IP2O = IPLAST(IELO)
        IPTE = IP1+IGTEDW-1
        DX0 = XP(IP1O)-XP(IPTE)
        DY0 = YP(IP1O)-YP(IPTE)
cc        write(*,*) 'WAKERESET first dx0,dy0 ',dx0,dy0
        DO IP = IP1O, IP2O
          IPM = IPTE + (IP-IP1O)
          DX = XP(IP)-XP(IPM)
          DY = YP(IP)-YP(IPM)
cc          write(*,*) 'WAKERESET ip,dx,dy ',ip,dx,dy
          YP(IP) = YP(IP) - (DY-DY0)
        END DO
      ENDIF
C
C---- change grid boundaries, relax grid and update vortex wakes
      CALL UPDGRD
      CALL UPDROTWAK
C
C---- invalidate existing solution
      LNCVP = .FALSE.
      LQAIC = .FALSE.
      LQGIC = .FALSE.
      LQCNT = .FALSE.
      LGSYS = .FALSE.
      LGAMU = .FALSE.
      LGAMA = .FALSE.
C
      RETURN
      END ! WAKERESET




      SUBROUTINE WAKMOV(IELA,ISQ)
C-------------------------------------------------------------
C     Moves existing wake (element IELA) by
C     aligning panels with current velocities.
C     Panel side flag ISQ indicates which velocity to use.
C       ISQ=0 uses mean surface velocities QC
C       ISQ=1 uses right surface velocities QCR
C       ISQ=2 uses left  surface velocities QCL
C
C     Panel points are moved preserving existing panel lengths 
C     while aligning the panels with the current velocity 
C     vectors on the panel.
C-------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(NETYPE(IELA).EQ.0) THEN
C----- this element is a body panel
       WRITE(*,*) 'WAKMOV: Element',IELA,' is body, cannot be moved!'
       RETURN
      ENDIF
C
      IEL = IELA
C
C---- Save X location for end of wake
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      XP2 = XP(IP2)
C
      CALL WAKMOVR(IEL,XP,YP,ISQ)
C
C---- Reset last point to initial (unchanged) XWAKE location
      XP(IP2) = XP2
C
C---- respline node powgrints
      CALL XYPSPL(IEL)
C---- reset control points and normal vectors
      CALL XYCSET(IEL)
      CALL ANPSET(IEL)
C
c      WRITE(*,*) 'Generating dq/d(gam,sig) influence of moved wake...'
c      CALL QAIC1(.FALSE., 1,NCTOT, IEL,IEL)
c      DO KEL = 1, NEL
c        ICE = NCX + KEL
c        IF(LBODY(KEL)) CALL QAIC1(.FALSE.,ICE,ICE,IEL,IEL)
c      ENDDO
C
c      IP1 = IPFRST(IEL)
c      IP2 = IPLAST(IEL)
c      CALL GSYS(.FALSE., .FALSE., .TRUE., IP1,IP2, 1,0)
C
C---- free existing unit wake strength vectors (if any), mark vectors invalid
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
C
      RETURN
      END ! WAKMOV



      SUBROUTINE WAKMOVR(IEL,XPNEW,YPNEW,ISQ)
C----------------------------------------------------------
C---- Routine to move points on element to new positions
C     to align with velocity vector on side ISQ (0/1/2)
C----------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION XPNEW(*), YPNEW(*)
C
      DIMENSION DSP(IPX)
C
C---- speed below this is treated as stagnation
      QSTAG = 0.0001*QREF
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
C
C---- calculate velocity and its Jacobians at all wake points
ccc      CALL QAIC1(.FALSE., IC1,IC2, 1,NEL)
C---- set velocities on both sides of panels at all control points
ccc      CALL QQCALC(NCTOT,ANC, IPCO,IPCP, GAM,SIG, QC, QCL,QCR )
C
      DO IP = IP1, IP2-1
        DXP = XP(IP+1) - XP(IP)
        DYP = YP(IP+1) - YP(IP)
        DSP(IP) = SQRT(DXP**2 + DYP**2)
      ENDDO
C
C---- march down wake setting wake points
      IP = IPCO(IC1)
      XPNEW(IP) = XP(IP)
      YPNEW(IP) = YP(IP)
C
      DO IC = IC1, IC2
        IPO = IPCO(IC)
        IPP = IPCP(IC)
C
        IF(ISQ .EQ. 0) THEN
          QX = QC(1,IC)
          QY = QC(2,IC)
        ELSEIF(ISQ .EQ. 1) THEN
          QX = QCR(1,IC)
          QY = QCR(2,IC)
        ELSEIF(ISQ .EQ. 2) THEN
          QX = QCL(1,IC)
          QY = QCL(2,IC)
        ELSE
          QX = QC(1,IC)
          QY = QC(2,IC)
        ENDIF
C
        QMAG = SQRT(QX**2 + QY**2)
        IF(QMAG.LT.QSTAG) THEN
C------- stagnation... don't try to move points
         RETURN
        ENDIF
C
        UN = QX/QMAG
        VN = QY/QMAG
        XPNEW(IPP) = XPNEW(IPO) + UN*DSP(IPO)
        YPNEW(IPP) = YPNEW(IPO) + VN*DSP(IPO)
C
      END DO
C
      RETURN
      END ! WAKMOVR



