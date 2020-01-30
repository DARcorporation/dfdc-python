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
!

SUBROUTINE ADJPANL
    !-----------------------------------------------------------------
    !     Readjusts panel nodes for CB and duct to give
    !     spacing compatible with a vortex wake grid between
    !     the elements.
    !-----------------------------------------------------------------
    INCLUDE 'DFDC.INC'
    LOGICAL LOVRMAX
    !
    PARAMETER (ITX = 600)
    DIMENSION XT1(ITX), YT1(ITX), ST1(ITX), XTS1(ITX), YTS1(ITX)
    DIMENSION XT2(ITX), YT2(ITX), ST2(ITX), XTS2(ITX), YTS2(ITX)
    DIMENSION XX(ITX), YY(ITX), SS(ITX)
    !
    !---- Flag for spacing in overlap area to TE (T uses max # of points from CB
    !     or duct for overlap area, F simply copies spacing from opposite body)
    LOVRMAX = .TRUE.
    !
    !---- Copy CB panel geometry to temporary array
    IEL = 1
    IP1 = IPFRST(IEL)
    IP2 = IPLAST(IEL)
    N = IP2 - IP1 + 1
    IF(N.GT.ITX) THEN
        WRITE(*, *) 'ADJPANL: ITX too small for N ', ITX, N
        STOP
    ENDIF
    DO I = 1, N
        IP = IP1 + I - 1
        XT1(I) = XP(IP)
        YT1(I) = YP(IP)
        ST1(I) = SP(IP)
        XTS1(I) = XPS(IP)
        YTS1(I) = YPS(IP)
    END DO
    NT1 = N
    !
    XLE1 = XPLE(IEL)
    YLE1 = YPLE(IEL)
    SLE1 = SPLE(IEL)
    !---- Use TE point for CB foil
    XTE1 = XP(IP1)
    YTE1 = YP(IP1)
    !
    !---- Copy duct panel geometry to temporary array
    IEL = 2
    IP1 = IPFRST(IEL)
    IP2 = IPLAST(IEL)
    N = IP2 - IP1 + 1
    IF(N.GT.ITX) THEN
        WRITE(*, *) 'ADJPANL: ITX too small for N ', ITX, N
        STOP
    ENDIF
    DO I = 1, N
        IP = IP1 + I - 1
        XT2(I) = XP(IP)
        YT2(I) = YP(IP)
        ST2(I) = SP(IP)
        XTS2(I) = XPS(IP)
        YTS2(I) = YPS(IP)
    END DO
    NT2 = N
    !
    XLE2 = XPLE(IEL)
    YLE2 = YPLE(IEL)
    SLE2 = SPLE(IEL)
    !---- Use inner TE point for duct foil
    XTE2 = XP(IP2)
    YTE2 = YP(IP2)
    !
    !---- Axial region for rotor(s) defined by overlap of LE/TE's
    XLEMAX = MAX(XLE1, XLE2)
    XTEMIN = MIN(XTE1, XTE2)
    !
    !cc      GO TO 100
    !
    !
    !
    !
    !-----------------------------------------------------------------
    !-----------------------------------------------------------------
    !---- start process with upstream rotor
    !
    NR = 1
    !---- Check rotor station for sanity
    IF(XDISK(NR).LT.XLEMAX .OR.&
            XDISK(NR).GT.XTEMIN) THEN
        WRITE(*, *) '*Rotor station out of bounds ', XDISK(NR)
        WRITE(*, *) '  XLE1 =', XLE1, '  XTE1 =', XTE1
        WRITE(*, *) '  XLE2 =', XLE2, '  XTE2 =', XTE2
        !c        CALL ASKR('Enter rotor disk X location^',XDISK(NR))
        !---- default rotor is halfway along passage
        XDISK(NR) = 0.5 * (XLEMAX + XTEMIN)
        WRITE(*, *) '*Rotor station set to ', XDISK(NR)
    ENDIF
    IF(LDBG) WRITE(*, *) ' '
    !
    !-----------------------------------------------------------------
    !---- find location of rotor on CB wall
    STCB = ST1(1)
    CALL SINVRT(STCB, XDISK(NR), XT1, XTS1, ST1, NT1)
    YCB = SEVAL(STCB, YT1, YTS1, ST1, NT1)
    RHUB(NR) = YCB
    IF(LDBG) WRITE(*, *) 'Rhub,S on CB @ ', RHUB(NR), STCB
    !---- Find nearest point on CB corresponding to rotor
    IPRCB = ISGFIND(STCB, ST1, NT1)
    IPROTCB(NR) = IPRCB
    IF(LDBG) WRITE(*, *) 'Nearest CB point to rotor @ ', &
            IPRCB, ST1(IPRCB)
    !---- Readjust point spacing to put grid point at rotor intersection
    SBEG = ST1(1)
    SEND = ST1(NT1)
    SOLD = ST1(IPRCB)
    CALL SGSHFT(SBEG, SEND, SOLD, STCB, ST1, SS, NT1)
    !---- Generate new X,Y coordinates
    CALL SGSPLUPD(XT1, XTS1, YT1, YTS1, ST1, NT1, &
            XT1, XTS1, YT1, YTS1, SS, NT1)
    CALL SGCOPY(SS, ST1, NT1)
    !
    !
    !-----------------------------------------------------------------
    !---- find location of rotor on duct wall
    STDW = ST2(NT2)
    CALL SINVRT(STDW, XDISK(NR), XT2, XTS2, ST2, NT2)
    YDW = SEVAL(STDW, YT2, YTS2, ST2, NT2)
    RTIP(NR) = YDW
    IF(LDBG) WRITE(*, *) 'Rtip,S on duct @ ', RTIP(NR), STDW
    !---- Find nearest point on duct corresponding to rotor
    IPRDW = ISGFIND(STDW, ST2, NT2)
    IPROTDW(NR) = IPRDW
    IF(LDBG) WRITE(*, *) 'Nearest Duct point to rotor @ ', &
            IPRDW, ST2(IPRDW)
    !---- Readjust point spacing to put grid point at rotor intersection
    SBEG = SLE2
    SEND = ST2(NT2)
    SOLD = ST2(IPRDW)
    CALL SGSHFT(SBEG, SEND, SOLD, STDW, ST2, SS, NT2)
    !---- Generate new X,Y coordinates
    CALL SGSPLUPD(XT2, XTS2, YT2, YTS2, ST2, NT2, &
            XT2, XTS2, YT2, YTS2, SS, NT2)
    CALL SGCOPY(SS, ST2, NT2)
    !
    !c      GO TO 100
    !
    !-----------------------------------------------------------------
    !---- Match spacing in overlap area downstream of rotor (to TE's)
    !     to define I direction slipstream grid
    !
    !---- Duct TE upstream of CB TE
    IF(XTE2.LE.XTE1) THEN
        !
        IF(LDBG) WRITE(*, *) 'XTEduct <= XTECB ', XTE2, XTE1
        !---- Find point on CB corresponding to duct TE
        STE = ST1(1)
        CALL SINVRT(STE, XTE2, XT1, XTS1, ST1, NT1)
        IF(LDBG) WRITE(*, *) 'CB point @ Xduct TE @ ', STE
        !---- Find nearest point on CB corresponding to duct TE
        ITE = 0
        ITE = ISGFIND(STE, ST1, NT1)
        IF(LDBG) WRITE(*, *) 'Nearest CB point to duct TE @ ', &
                ITE, STE, ST1(ITE)
        !---- If this point can be moved (not TE point!!), align to duct TE
        IF(ITE.GT.1) THEN
            !---- Readjust CB point spacing to put grid point at duct TE X location
            SBEG = ST1(1)
            SEND = ST1(IPRCB)
            SOLD = ST1(ITE)
            CALL SGSHFT(SBEG, SEND, SOLD, STE, ST1, SS, NT1)
            !---- Update CB spline definition
            CALL SGSPLUPD(XT1, XTS1, YT1, YTS1, ST1, NT1, &
                    XT1, XTS1, YT1, YTS1, SS, NT1)
            CALL SGCOPY(SS, ST1, NT1)
        ENDIF
        !
        !---- Respace overlap zone from rotor to CB TE point
        NDUCT = NT2 - IPRDW + 1
        NCB = IPRCB - ITE + 1
        !---- Check # points in overlap zone, use densest paneling from CB or duct
        IF(NCB.LT.NDUCT .OR. (.NOT.LOVRMAX)) THEN
            !---- Respace CB points from rotor to duct TE point with # from duct
            NNEW = NDUCT
            NOLD = NCB
            CALL SGCOPY(ST1, SS, NT1)
            CALL SGRENUM(SS(ITE), NOLD, NNEW)
            CALL SGCOPY(ST1(IPRCB), SS(ITE + NNEW - 1), NT1 - IPRCB + 1)
            NT1NEW = NT1 + NNEW - NOLD
            !---- Generate new X,Y coordinates
            CALL SGSPLUPD(XT1, XTS1, YT1, YTS1, ST1, NT1, &
                    XT1, XTS1, YT1, YTS1, SS, NT1NEW)
            CALL SGCOPY(SS, ST1, NT1NEW)
            NT1 = NT1NEW
            IPRCB = IPROTCB(NR) + NNEW - NOLD
            IPROTCB(NR) = IPRCB
            !
        ELSEIF(NCB.GT.NDUCT) THEN
            !---- Respace duct points from rotor to duct TE point with # from CB
            NNEW = NCB
            NOLD = NDUCT
            CALL SGCOPY(ST2, SS, NT2)
            CALL SGRENUM(SS(IPRDW), NOLD, NNEW)
            NT2NEW = NT2 + NNEW - NOLD
            !---- Generate new X,Y coordinates
            CALL SGSPLUPD(XT2, XTS2, YT2, YTS2, ST2, NT2, &
                    XT2, XTS2, YT2, YTS2, SS, NT2NEW)
            CALL SGCOPY(SS, ST2, NT2NEW)
            NT2 = NT2NEW
        ENDIF
        !
        !
        !---- CB TE upstream of duct TE
    ELSEIF(XTE2.GT.XTE1) THEN
        !
        IF(LDBG) WRITE(*, *) 'XTEduct > XTECB ', XTE2, XTE1
        !---- Find point on duct corresponding to CB TE
        STE = ST2(NT2)
        CALL SINVRT(STE, XTE1, XT2, XTS2, ST2, NT2)
        !---- Find nearest point on duct corresponding to CB TE
        ITE = ISGFIND(STE, ST2, NT2)
        IF(LDBG) WRITE(*, *) 'Nearest Duct point to CB TE @ ', &
                ITE, STE, ST2(ITE)
        !---- If this point can be moved (not TE point!!), align to CB TE
        IF(ITE.LT.NT2) THEN
            !---- Readjust duct point spacing to put grid point at CB TE X location
            SBEG = ST2(IPRDW)
            SEND = ST2(NT2)
            SOLD = ST2(ITE)
            CALL SGSHFT(SBEG, SEND, SOLD, STE, ST2, SS, NT2)
            CALL SGSPLUPD(XT2, XTS2, YT2, YTS2, ST2, NT2, &
                    XT2, XTS2, YT2, YTS2, SS, NT2)
            CALL SGCOPY(SS, ST2, NT2)
        ENDIF
        !
        !---- Respace overlap zone from rotor to CB TE point
        NDUCT = ITE - IPRDW + 1
        NCB = IPRCB
        !---- Check # points in overlap zone, use densest paneling from CB or duct
        IF(NCB.GT.NDUCT .OR. (.NOT.LOVRMAX)) THEN
            !---- Respace duct points from rotor to CB TE point with # from CB
            NNEW = NCB
            NOLD = NDUCT
            CALL SGCOPY(ST2, SS, NT2)
            CALL SGRENUM(SS(IPRDW), NOLD, NNEW)
            CALL SGCOPY(ST2(ITE), SS(IPRDW + NNEW - 1), NT2 - ITE + 1)
            NT2NEW = NT2 + NNEW - NOLD
            !---- Generate new X,Y coordinates
            CALL SGSPLUPD(XT2, XTS2, YT2, YTS2, ST2, NT2, &
                    XT2, XTS2, YT2, YTS2, SS, NT2NEW)
            CALL SGCOPY(SS, ST2, NT2NEW)
            NT2 = NT2NEW
        ELSEIF(NCB.LT.NDUCT) THEN
            !---- Respace CB points from rotor to CB TE point with # from duct
            NNEW = NDUCT
            NOLD = NCB
            CALL SGCOPY(ST1, SS, NT1)
            CALL SGRENUM(SS(1), NOLD, NNEW)
            CALL SGCOPY(ST1(NOLD + 1), SS(NNEW + 1), NT1 - NOLD)
            NT1NEW = NT1 + NNEW - NOLD
            !---- Generate new X,Y coordinates
            CALL SGSPLUPD(XT1, XTS1, YT1, YTS1, ST1, NT1, &
                    XT1, XTS1, YT1, YTS1, SS, NT1NEW)
            CALL SGCOPY(SS, ST1, NT1NEW)
            NT1 = NT1NEW
            IPRCB = IPROTCB(NR) + NNEW - NOLD
            IPROTCB(NR) = IPRCB
        ENDIF
        !
    ENDIF
    !
    !      GO TO 100
    !
    !-----------------------------------------------------------------
    !-----------------------------------------------------------------
    !---- Process downstream rotor(s), fiddling grid to put move points
    !     to match downstream rotor line
    DO NR = 2, NROTOR
        !
        !---- Check rotor station for sanity
        IF(LDBG) WRITE(*, *) ' '
        IF(XDISK(NR).LT.XLEMAX .OR.&
                XDISK(NR).GT.XTEMIN) THEN
            WRITE(*, *) '*Rotor#N station out of bounds ', XDISK(NR)
            WRITE(*, *) '  XLE1 =', XLE1, '  XTE1 =', XTE1
            WRITE(*, *) '  XLE2 =', XLE2, '  XTE2 =', XTE2
            !c        CALL ASKR('Enter rotor disk X location^',XDISK)
            !---- default rotor #2 is halfway along passage
            XDISK(NR) = 0.5 * (XDISK(1) + XTEMIN)
            WRITE(*, *) '*Rotor#N station set to ', XDISK(NR)
        ENDIF
        !
        !-----------------------------------------------------------------
        !---- find location of rotor on CB wall
        STCB = ST1(1)
        CALL SINVRT(STCB, XDISK(NR), XT1, XTS1, ST1, NT1)
        YCB = SEVAL(STCB, YT1, YTS1, ST1, NT1)
        RHUB(NR) = YCB
        IF(LDBG) WRITE(*, *) 'Rhub,S #N on CB @ ', RHUB(NR), STCB
        !---- Find nearest point on CB corresponding to rotor
        IPRCB2 = ISGFIND(STCB, ST1, NT1)
        IPROTCB(NR) = IPRCB2
        IF(LDBG) WRITE(*, *) 'Nearest CB point to rotor @ ', &
                IPRCB2, ST1(IPRCB2)
        !---- Readjust point spacing to put grid point at rotor intersection
        SBEG = ST1(1)
        SEND = ST1(IPRCB)
        SOLD = ST1(IPRCB2)
        CALL SGSHFT(SBEG, SEND, SOLD, STCB, ST1, SS, NT1)
        !---- Generate new X,Y coordinates
        CALL SGSPLUPD(XT1, XTS1, YT1, YTS1, ST1, NT1, &
                XT1, XTS1, YT1, YTS1, SS, NT1)
        CALL SGCOPY(SS, ST1, NT1)
        !
        !-----------------------------------------------------------------
        !---- find location of rotor on duct wall
        STDW = ST2(IPRDW)
        CALL SINVRT(STDW, XDISK(NR), XT2, XTS2, ST2, NT2)
        YDW = SEVAL(STDW, YT2, YTS2, ST2, NT2)
        RTIP(NR) = YDW
        IF(LDBG) WRITE(*, *) 'Rtip,S #N on duct @ ', RTIP(NR), STDW
        !---- index of rotor#2 from rotor#1 must match offset from rotor#1 on CB
        IPRDW2 = IPRDW + (IPRCB - IPRCB2)
        IPROTDW(NR) = IPRDW2
        IF(LDBG) WRITE(*, *) 'Nearest Duct point to rotor @ ', &
                IPRDW2, ST2(IPRDW2)
        !---- Readjust point spacing to put grid point at rotor intersection
        SBEG = ST2(IPRDW)
        SEND = ST2(NT2)
        SOLD = ST2(IPRDW2)
        CALL SGSHFT(SBEG, SEND, SOLD, STDW, ST2, SS, NT2)
        !---- Generate new X,Y coordinates
        CALL SGSPLUPD(XT2, XTS2, YT2, YTS2, ST2, NT2, &
                XT2, XTS2, YT2, YTS2, SS, NT2)
        CALL SGCOPY(SS, ST2, NT2)
        !
    END DO  ! loop over NROTOR=>2
    !
    !---- Paneling (points) modified to match on CB and duct and
    !     to match rotor lines
    !
    !-----------------------------------------------------------------
    !---- Before installing new panel geometry check IEL>2 to save old
    !     geometry
    100  IF(NEL.GE.3) THEN
        NPNEW = NT1 + NT2
        NMOV = NPNEW - IPFRST(3) + 1
        IF(NMOV.GT.0) THEN
            !---- Move existing elements up to give room for modified elements 1 and 2
            DO IEL = NEL, 3, -1
                IP1 = IPFRST(IEL)
                IP2 = IPLAST(IEL)
                DO IP = IP2, IP1, -1
                    IPMOV = IP + NMOV
                    XP(IPMOV) = XP(IP)
                    YP(IPMOV) = YP(IP)
                END DO
                IPFRST(IEL) = IP1 + NMOV
                IPLAST(IEL) = IP2 + NMOV
            END DO
        ENDIF
    ENDIF
    !
    !-----------------------------------------------------------------
    !---- Put adjusted body wall points back into paneled point arrays
    IP = 0
    !---- CB element
    IEL = 1
    IP1 = IP + 1
    DO I = 1, NT1
        IP = IP + 1
        XP(IP) = XT1(I)
        YP(IP) = YT1(I)
    END DO
    IP2 = IP
    IPFRST(IEL) = IP1
    IPLAST(IEL) = IP2
    I = IP1
    N = IP2 - IP1 + 1
    !---- spline and set other stuff for element IEL
    CALL XYPSPL(IEL)
    !---- set indices of TE panel, if any
    CALL ITPSET(IEL)
    !---- update pointers to CB wall point at rotor line(s)
    DO NR = 1, NROTOR
        IPROTCB(NR) = IP1 - 1 + IPROTCB(NR)
    END DO
    !
    !---- Put new duct wall points back into paneled point arrays
    IEL = 2
    IP1 = IP + 1
    DO I = 1, NT2
        IP = IP + 1
        XP(IP) = XT2(I)
        YP(IP) = YT2(I)
    END DO
    IP2 = IP
    IPFRST(IEL) = IP1
    IPLAST(IEL) = IP2
    I = IP1
    N = IP2 - IP1 + 1
    !---- spline and set other stuff for element IEL
    CALL XYPSPL(IEL)
    !---- set indices of TE panel, if any
    CALL ITPSET(IEL)
    !---- update pointer to duct wall point at rotor line
    DO NR = 1, NROTOR
        IPROTDW(NR) = IP1 - 1 + IPROTDW(NR)
        IF(LDBG) WRITE(*, *) 'NR IPROTCB IPROTDW  ', &
                NR, IPROTCB(NR), IPROTDW(NR)
    END DO
    !
    !---- Shift old points in any remaining defined elements to match new ordering
    IF(NEL.GE.3) THEN
        NPNEW = NT1 + NT2
        NMOV = NPNEW + 1 - IPFRST(3)
        IF(NMOV.LT.0) THEN
            !---- Move existing elements down to give room for modified elements 1 and 2
            DO IEL = 3, NEL
                IP1 = IPFRST(IEL)
                IP2 = IPLAST(IEL)
                DO IP = IP1, IP2
                    IPMOV = IP + NMOV
                    XP(IPMOV) = XP(IP)
                    YP(IPMOV) = YP(IP)
                END DO
                IPFRST(IEL) = IP1 + NMOV
                IPLAST(IEL) = IP2 + NMOV
            END DO
        ENDIF
    ENDIF
    !---- Reset total # of paneled points
    NPTOT = IPLAST(NEL)
    !
    !cc        write(14,99) i,xt1(i),xx(i),xdisk(nr)
    99     format(i5, 5(1x, f12.6))
    !
    !---- Set flag for rotor defined (at least geometrically)
    !c      LROTOR = .TRUE.
    !
    !---- invalidate any existing solution
    LNCVP = .FALSE.
    LQAIC = .FALSE.
    LQGIC = .FALSE.
    LQCNT = .FALSE.
    LSYSP = .FALSE.
    LGSYS = .FALSE.
    LGAMU = .FALSE.
    LGAMA = .FALSE.
    LSIGP = .FALSE.
    LSIGM = .FALSE.
    !
    RETURN
END
! ADJPANL


SUBROUTINE SGSPLUPD(X1, XS1, Y1, YS1, S1, N1, &
        X2, XS2, Y2, YS2, S2, N2)
    !-----------------------------------------------------------------
    !     Generates points and spline definitions for spline curve
    !     X1,XS1,Y1,YS1,S1,N1 (with N1 points at S1) to new spline
    !     with spacing S2 with N2 points.
    !     New curve with spline and points returned in X2,XS2,Y2,YS2,S2,N2
    !-----------------------------------------------------------------
    !
    DIMENSION X1(*), Y1(*), S1(*), XS1(*), YS1(*)
    DIMENSION X2(*), Y2(*), S2(*), XS2(*), YS2(*)
    !
    PARAMETER (IX = 600)
    DIMENSION XX(IX), YY(IX)
    !
    IF(N1.GT.IX .OR. N2.GT.IX) THEN
        WRITE(*, *) 'SGSPLUPD number of input points exceeds IX ', N1, N2
        STOP
    ENDIF
    !
    !---- Generate new X,Y coordinates at S2 locations
    DO I = 1, N2
        XX(I) = SEVAL(S2(I), X1, XS1, S1, N1)
        YY(I) = SEVAL(S2(I), Y1, YS1, S1, N1)
    END DO
    !
    !---- Replace points in output arrays for redefined curve
    DO I = 1, N2
        X2(I) = XX(I)
        Y2(I) = YY(I)
    END DO
    !---- Respline new points
    CALL SEGSPL(X2, XS2, S2, N2)
    CALL SEGSPL(Y2, YS2, S2, N2)
    !
    RETURN
END