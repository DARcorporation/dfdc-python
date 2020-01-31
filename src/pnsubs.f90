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

SUBROUTINE PANDEF(IELDEF)
    !-----------------------------------------------
    !     Sets default paneling parameters
    !-----------------------------------------------
    INCLUDE 'DFDC.INC'
    DIMENSION SBTOT(NEX)
    !
    IF(IELDEF.EQ.0) THEN
        IEL1 = 1
        IEL2 = NBEL
    ELSE
        IEL1 = IELDEF
        IEL2 = IELDEF
    ENDIF
    !
    SBMIN = 0.0
    SBMAX = 0.0
    DO IEL = IEL1, IEL2
        SBTOT(IEL) = SB(IBLAST(IEL)) - SB(IBFRST(IEL))
        IF(NETYPE(IEL).EQ.0) THEN
            SBMIN = MIN(SBMIN, SBTOT(IEL))
            SBMAX = MAX(SBMAX, SBTOT(IEL))
        ENDIF
    ENDDO
    !---- use mean element arclength for scaling
    SBMAG = 0.5 * (SBMIN + SBMAX)
    !
    ANPAN = FLOAT(NPANDEF)
    !---- set min, max # panels per element
    NPANMIN = NPANDEF / 2
    NPANMAX = 2 * NPANDEF
    !
    DO IEL = IEL1, IEL2
        CVEX(IEL) = 1.0
        SMOF(IEL) = 1.0
        FSLE(IEL) = 0.6
        FSTE(IEL) = 0.6
        !
        IF(SBTOT(IEL).EQ.0.0) THEN
            NPAN(IEL) = 1
        ELSE
            IF (NETYPE(IEL).EQ.0) THEN
                !-------- solid airfoil element
                NPAN(IEL) = INT(ANPAN * (SBTOT(IEL) / SBMAG)**FPANDEF)
                NPAN(IEL) = MAX(NPANMIN, NPAN(IEL))
                NPAN(IEL) = MIN(NPANMAX, NPAN(IEL))
                CVEX(IEL) = 1.0
                FSLE(IEL) = 0.6
                FSTE(IEL) = 0.6
                !
            ELSEIF(NETYPE(IEL).EQ.1&
                    .OR. NETYPE(IEL).EQ.2&
                    .OR. NETYPE(IEL).EQ.5) THEN
                !-------- wake sheet, wake line or source line element
                NPAN(IEL) = INT(ANPAN / 4.0 * (SBTOT(IEL) / SBMAG)**FPANDEF)
                !---- not less than 2 panels for these...
                NPAN(IEL) = MAX(2, NPAN(IEL))
                CVEX(IEL) = 1.0
                FSLE(IEL) = 0.6
                FSTE(IEL) = 5.0
            ELSE
                !-------- point or ring element (these won't really be used)
                NPAN(IEL) = 1
                CVEX(IEL) = 1.0
                FSLE(IEL) = 1.0
                FSTE(IEL) = 1.0
            ENDIF
        ENDIF
        !
        !------ refinement-station parameters
        DO IRPN = 1, NRPNX
            CRRAT(IRPN, IEL) = 1.0
            SRPN1(IRPN, IEL) = 0.0
            SRPN2(IRPN, IEL) = 0.0
        ENDDO
    ENDDO
    !
    !---- Set flag for spacing defined, unset user-defined flag
    LRSPCDEF = .TRUE.
    LRSPCUSR = .FALSE.
    !
    RETURN
END
! PANDEF


SUBROUTINE PANCOP
    !------------------------------------------------------------
    !     Copies XB,YB arrays directly into panel arrays.
    !------------------------------------------------------------
    INCLUDE 'DFDC.inc'
    !
    IF(NBTOT.GT.IPX) STOP 'PANCOP: Array overflow on IPX.'
    !
    !---- total number of panel nodes
    NPTOT = NBTOT
    NEL = NBEL
    !c      write(60,*) 'PANCOP '
    !c      write(60,*) 'NBEL,NBTOT ',NBEL,NBTOT
    !
    !---- panel node locations
    DO IB = 1, NBTOT
        IP = IB
        !c        write(60,99) ib,xb(ip),yb(ip),xp(ip),yp(ip)
        XP(IP) = XB(IB)
        YP(IP) = YB(IB)
    ENDDO
    99   format(i5, 6F12.6)
    !
    !---- reference point location for each element
    DO IEL = 1, NEL
        XPREFE(IEL) = XBREFE(IEL)
        YPREFE(IEL) = YBREFE(IEL)
        !
        XPCENT(IEL) = XBCEN2DT(IEL)
        XPCENT(IEL) = XBCEN2DT(IEL)
    ENDDO
    !
    !---- set new stuff for each element
    DO IEL = 1, NEL
        !c       write(60,*) IEL,IPFRST(IEL),IPLAST(IEL),IBFRST(IEL),IBLAST(IEL)
        IPFRST(IEL) = IBFRST(IEL)
        IPLAST(IEL) = IBLAST(IEL)
        IF(IPLAST(IEL).GT.IPX) STOP 'PANCOP: Array overflow on IPX'
        !
        !------ set prescribed panel number to actual new panel number
        NPAN(IEL) = IPLAST(IEL) - IPFRST(IEL) + 1
        !
        !------ spline and set other stuff for element IEL
        CALL XYPSPL(IEL)
        !
        !------ set indices of TE panel, if any
        CALL ITPSET(IEL)
        !
        I = IPFRST(IEL)
        N = IPLAST(IEL) - IPFRST(IEL) + 1
        CALL MINMAX(N, XP(I), XPMINE(IEL), XPMAXE(IEL))
        CALL MINMAX(N, YP(I), YPMINE(IEL), YPMAXE(IEL))
    ENDDO
    !
    !---- set overall min,max over all elements
    CALL MINMAX(NEL, XPMINE, XPMIN, DUMMY)
    CALL MINMAX(NEL, YPMINE, YPMIN, DUMMY)
    CALL MINMAX(NEL, XPMAXE, DUMMY, XPMAX)
    CALL MINMAX(NEL, YPMAXE, DUMMY, YPMAX)
    !
    !---- we now have a new geometry... invalidate any existing solution
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
    LGSAME = .TRUE.
    !
    RETURN
END
! PANCOP

SUBROUTINE PANGEN(LQUERY)
    !------------------------------------------------------------------
    !     Generates panel nodes X,Y from buffer-airfoil XB,YB arrays.
    !     If LQUERY=T, then user's interaction is requested.
    !------------------------------------------------------------------
    INCLUDE 'DFDC.INC'
    LOGICAL LQUERY
    !
    LOGICAL LCPLOT
    DIMENSION SG(IPX)
    !
    DATA LCPLOT / .FALSE. /
    !
    IF(NBEL.EQ.0) RETURN
    NEL = NBEL
    !
    !---- process each element
    DO 100 IEL = 1, NBEL
        !
        IB1 = IBFRST(IEL)
        IB2 = IBLAST(IEL)
        NBPTS = IB2 - IB1 + 1
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        NPPTS = IP2 - IP1 + 1
        IF(LDBG) THEN
            WRITE(*, *) 'PANGEN IEL ', IEL
            WRITE(*, *) 'IB1,IB2 ', IB1, IB2
            WRITE(*, *) 'IP1,IP2 ', IP1, IP2
        ENDIF
        !
        IF(NPPTS.LE.1 .AND. NBPTS.LE.1) THEN
            !------- don't try to panel point or ring element
            NG = 1
            SG(1) = 0.
            !
        ELSE
            !
            IF(LQUERY) THEN
                !-------- set current spacing so SGCURV has something to plot
                IF(NPPTS.GT.1) THEN
                    !--------- set SG from existing paneled geometry, if defined
                    NG = NPPTS
                    DO IP = IP1, IP2
                        IG = IP - IP1 + 1
                        SG(IG) = (SP(IP) - SP(IP1)) / (SP(IP2) - SP(IP1))
                    ENDDO
                ELSEIF(NBPTS.GT.1) THEN
                    !--------- or, set SG from buffer geometry
                    NG = NBPTS
                    DO IB = IB1, IB2
                        IG = IB - IB1 + 1
                        SG(IG) = (SB(IB) - SB(IB1)) / (SB(IB2) - SB(IB1))
                    ENDDO
                    IF(NPAN(IEL).LE.0) NPAN(IEL) = NBPTS
                ELSE
                    WRITE(*, *)
                    WRITE(*, *) '  No current paneling to modify'
                    WRITE(*, *) '  Input case and/or Execute PANE command'
                    RETURN
                ENDIF
                !
            ELSE
                !--------- no initial paneling for SGCURV
                SG(1) = 0.
            ENDIF
            !
            !         WRITE(*,1150) IEL
            ! 1150    FORMAT(
            !     & /' ============================================================='
            !     &//'  Element', I4)
            !
            !------- surface element... generate spacing array SG
            !
            CALL SGCURV(LQUERY, LCPLOT, NBPTS, &
                    XB(IB1), XBS(IB1), &
                    YB(IB1), YBS(IB1), &
                    SB(IB1), &
                    IPX, NG, SG, &
                    NPAN(IEL), CVEX(IEL), SMOF(IEL), FSLE(IEL), FSTE(IEL), &
                    NRPNX, SRPN1(1, IEL), SRPN2(1, IEL), CRRAT(1, IEL))
        ENDIF
        !
        IF(LQUERY) THEN
            !---- Set flag for spacing defined, unset user-defined flag
            LRSPCDEF = .TRUE.
            LRSPCUSR = .TRUE.
        ENDIF
        !
        !------ starting and ending indices for this element IEL
        IF(IEL.EQ.1) THEN
            IP1 = 1
            IP2 = NG
        ELSE
            IP1 = IPLAST(IEL - 1) + 1
            IP2 = IPLAST(IEL - 1) + NG
        ENDIF
        IF(IP2.GT.IPX) STOP 'PANGEN: Array overflow on IPX'
        !
        IF(LQUERY) THEN
            IPDEL = IP2 - IPLAST(IEL)
            IF(NPTOT + IPDEL.GT.IPX) STOP 'PANGEN: Array overflow on IPX'
            !
            !-------- move current points in subsequent elements to make/take-up room
            IF(IPDEL.GT.0) THEN
                KP1 = NPTOT
                KP2 = IPLAST(IEL) + 1
                KPD = -1
            ELSE
                KP1 = IPLAST(IEL) + 1
                KP2 = NPTOT
                KPD = +1
            ENDIF
            DO KP = KP1, KP2, KPD
                SP(KP + IPDEL) = SP(KP)
                XP(KP + IPDEL) = XP(KP)
                YP(KP + IPDEL) = YP(KP)
                XPS(KP + IPDEL) = XPS(KP)
                YPS(KP + IPDEL) = YPS(KP)
            ENDDO
            DO KEL = IEL + 1, NEL
                IPFRST(KEL) = IPFRST(KEL) + IPDEL
                IPLAST(KEL) = IPLAST(KEL) + IPDEL
            ENDDO
        ENDIF
        !
        !------ set panel nodes at fractional arc length locations SG
        DO IP = IP1, IP2
            IG = IP - IP1 + 1
            SBI = SB(IB1) + (SB(IB2) - SB(IB1)) * SG(IG)
            !
            XP(IP) = SEVAL(SBI, XB(IB1), XBS(IB1), SB(IB1), NBPTS)
            YP(IP) = SEVAL(SBI, YB(IB1), YBS(IB1), SB(IB1), NBPTS)
        ENDDO
        IPFRST(IEL) = IP1
        IPLAST(IEL) = IP2
        !
        !------ centroid location
        XPCENT(IEL) = XBCEN2DT(IEL)
        YPCENT(IEL) = YBCEN2DT(IEL)

        !------ spline and set other stuff for element IEL
        CALL XYPSPL(IEL)
        !
        !------ set indices of TE panel, if any
        CALL ITPSET(IEL)
        !
        I = IPFRST(IEL)
        N = IPLAST(IEL) - IPFRST(IEL) + 1
        CALL MINMAX(N, XP(I), XPMINE(IEL), XPMAXE(IEL))
        CALL MINMAX(N, YP(I), YPMINE(IEL), YPMAXE(IEL))
    100  CONTINUE
    !
    !---- set total number of panel vertices
    NPTOT = IPLAST(NEL)
    !
    !---- set overall min,max over all elements
    CALL MINMAX(NEL, XPMINE, XPMIN, DUMMY)
    CALL MINMAX(NEL, YPMINE, YPMIN, DUMMY)
    CALL MINMAX(NEL, XPMAXE, DUMMY, XPMAX)
    CALL MINMAX(NEL, YPMAXE, DUMMY, YPMAX)
    !
    !---- default reference point location for each element
    DO IEL = 1, NEL
        XPREFE(IEL) = 0.0
        YPREFE(IEL) = 0.0
    ENDDO
    !
    !---- we now have a new geometry... invalidate any existing solution
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
    LQSPEC = .FALSE.
    LGSAME = .FALSE.
    !
    RETURN
END
! PANGEN



SUBROUTINE PANWRT
    !-----------------------------------------
    !     Writes paneling-parameter file
    !-----------------------------------------
    INCLUDE 'DFDC.inc'
    CHARACTER*80 PFNDEF
    !
    1000 FORMAT(A)
    !
    LU = 8
    !
    !---- use existing paneling-parameter filename (if any) as the default
    20   WRITE(*, 1200) PFILE(1:60)
    1200 FORMAT(/' Enter output filename:  ', A)
    !
    READ(*, 1000) PFNDEF
    IF(PFNDEF(1:1).EQ.' ') THEN
        IF(PFILE(1:1).EQ.' ') GO TO 20
    ELSE
        PFILE = PFNDEF
    ENDIF
    !
    OPEN(LU, FILE = PFILE, STATUS = 'UNKNOWN', ERR = 20)
    REWIND(LU)
    !
    WRITE(LU, *) NEL, NRPNX
    DO IEL = 1, NEL
        WRITE(LU, *) NPAN(IEL)
        WRITE(LU, *) CVEX(IEL), SMOF(IEL), FSLE(IEL), FSTE(IEL)
        DO IRPN = 1, NRPNX
            WRITE(LU, *) SRPN1(IRPN, IEL), &
                    SRPN2(IRPN, IEL), &
                    CRRAT(IRPN, IEL)
        ENDDO
    ENDDO
    CLOSE(LU)
    !
    RETURN
END
! PANWRT



SUBROUTINE PANGET(FNAME, ERROR)
    !-----------------------------------------
    !     Reads paneling-parameter file FNAME
    !-----------------------------------------
    INCLUDE 'DFDC.inc'
    CHARACTER*(*) FNAME
    LOGICAL ERROR
    !
    LU = 8
    !
    IF(FNAME(1:1).EQ.' ') THEN
        CALL ASKS('Enter paneling-parameter filename^', FNAME)
    ENDIF
    !
    !---- strip off leading blanks and get number of characters
    CALL STRIP(FNAME, NFNAME)
    !
    WRITE(*, *)
    !
    OPEN(LU, FILE = FNAME, STATUS = 'OLD', ERR = 90)
    REWIND(LU)
    !
    READ(LU, *, ERR = 80) NEL0, NRPN0
    NEL0 = MIN(NEL0, NEX)
    NRPN0 = MIN(NRPN0, NRPNX)
    !
    DO IEL = 1, NEL0
        READ(LU, *, ERR = 80) NPAN(IEL)
        READ(LU, *, ERR = 80) CVEX(IEL), SMOF(IEL), FSLE(IEL), FSTE(IEL)
        DO IRPN = 1, NRPN0
            READ(LU, *, ERR = 80) SRPN1(IRPN, IEL), &
                    SRPN2(IRPN, IEL), &
                    CRRAT(IRPN, IEL)
        ENDDO
    ENDDO
    CLOSE(LU)
    WRITE(*, *) 'Paneling parameters read from file ', FNAME(1:NFNAME)
    LRSPCDEF = .TRUE.
    LRSPCUSR = .TRUE.
    ERROR = .FALSE.
    RETURN
    !
    80   CONTINUE
    CLOSE(LU)
    WRITE(*, *) 'READ error on panel-parameter file ', FNAME(1:NFNAME)
    ERROR = .TRUE.
    RETURN
    !
    90   CONTINUE
    !      WRITE(*,*) 'OPEN error on panel-parameter file ', FNAME(1:NFNAME)
    ERROR = .TRUE.
    RETURN
END
! PANGET