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

      SUBROUTINE PANDEF(IELDEF)
C-----------------------------------------------
C     Sets default paneling parameters
C-----------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION SBTOT(NEX)
C
      IF(IELDEF.EQ.0) THEN
       IEL1 = 1
       IEL2 = NBEL
      ELSE
       IEL1 = IELDEF
       IEL2 = IELDEF
      ENDIF
C
      SBMIN = 0.0
      SBMAX = 0.0
      DO IEL = IEL1, IEL2
        SBTOT(IEL) = SB(IBLAST(IEL)) - SB(IBFRST(IEL))
        IF(NETYPE(IEL).EQ.0) THEN
         SBMIN = MIN( SBMIN , SBTOT(IEL) )
         SBMAX = MAX( SBMAX , SBTOT(IEL) )
        ENDIF
      ENDDO
C---- use mean element arclength for scaling
      SBMAG = 0.5*(SBMIN+SBMAX)
C
      ANPAN = FLOAT(NPANDEF)
C---- set min, max # panels per element 
      NPANMIN =   NPANDEF/2
      NPANMAX = 2*NPANDEF
C
      DO IEL = IEL1, IEL2
        CVEX(IEL) = 1.0
        SMOF(IEL) = 1.0
        FSLE(IEL) = 0.6
        FSTE(IEL) = 0.6
C
        IF(SBTOT(IEL).EQ.0.0) THEN
         NPAN(IEL) = 1
        ELSE
         IF (NETYPE(IEL).EQ.0) THEN
C-------- solid airfoil element
          NPAN(IEL) = INT( ANPAN * (SBTOT(IEL)/SBMAG)**FPANDEF )
          NPAN(IEL) = MAX(NPANMIN,NPAN(IEL))
          NPAN(IEL) = MIN(NPANMAX,NPAN(IEL))
          CVEX(IEL) = 1.0
          FSLE(IEL) = 0.6
          FSTE(IEL) = 0.6
C
         ELSEIF(NETYPE(IEL).EQ.1
     &     .OR. NETYPE(IEL).EQ.2
     &     .OR. NETYPE(IEL).EQ.5) THEN
C-------- wake sheet, wake line or source line element
          NPAN(IEL) = INT( ANPAN/4.0 * (SBTOT(IEL)/SBMAG)**FPANDEF )
C---- not less than 2 panels for these...
          NPAN(IEL) = MAX(2,NPAN(IEL))
          CVEX(IEL) = 1.0
          FSLE(IEL) = 0.6
          FSTE(IEL) = 5.0
         ELSE
C-------- point or ring element (these won't really be used)
          NPAN(IEL) = 1
          CVEX(IEL) = 1.0
          FSLE(IEL) = 1.0
          FSTE(IEL) = 1.0
         ENDIF
        ENDIF
C
C------ refinement-station parameters
        DO IRPN = 1, NRPNX
          CRRAT(IRPN,IEL) = 1.0
          SRPN1(IRPN,IEL) = 0.0
          SRPN2(IRPN,IEL) = 0.0
        ENDDO
      ENDDO
C
C---- Set flag for spacing defined, unset user-defined flag
      LRSPCDEF = .TRUE.
      LRSPCUSR = .FALSE.
C
      RETURN
      END ! PANDEF


      SUBROUTINE PANCOP
C------------------------------------------------------------
C     Copies XB,YB arrays directly into panel arrays.
C------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(NBTOT.GT.IPX) STOP 'PANCOP: Array overflow on IPX.'
C
C---- total number of panel nodes
      NPTOT = NBTOT
      NEL   = NBEL
cc      write(60,*) 'PANCOP '
cc      write(60,*) 'NBEL,NBTOT ',NBEL,NBTOT
C
C---- panel node locations
      DO IB = 1, NBTOT
        IP = IB
cc        write(60,99) ib,xb(ip),yb(ip),xp(ip),yp(ip)
        XP(IP) = XB(IB)
        YP(IP) = YB(IB)
      ENDDO
 99   format(i5,6F12.6)
C
C---- reference point location for each element
      DO IEL = 1, NEL
        XPREFE(IEL) = XBREFE(IEL)
        YPREFE(IEL) = YBREFE(IEL)
C
        XPCENT(IEL) = XBCEN2DT(IEL)
        XPCENT(IEL) = XBCEN2DT(IEL)
      ENDDO
C
C---- set new stuff for each element
      DO IEL=1, NEL
cc       write(60,*) IEL,IPFRST(IEL),IPLAST(IEL),IBFRST(IEL),IBLAST(IEL)
        IPFRST(IEL) = IBFRST(IEL)
        IPLAST(IEL) = IBLAST(IEL)
        IF(IPLAST(IEL).GT.IPX) STOP 'PANCOP: Array overflow on IPX'
C
C------ set prescribed panel number to actual new panel number
        NPAN(IEL) = IPLAST(IEL) - IPFRST(IEL) + 1
C
C------ spline and set other stuff for element IEL
        CALL XYPSPL(IEL)
C
C------ set indices of TE panel, if any
        CALL ITPSET(IEL)
C
        I = IPFRST(IEL)
        N = IPLAST(IEL) - IPFRST(IEL) + 1
        CALL MINMAX(N,XP(I),XPMINE(IEL),XPMAXE(IEL))
        CALL MINMAX(N,YP(I),YPMINE(IEL),YPMAXE(IEL))
      ENDDO
C
C---- set overall min,max over all elements
      CALL MINMAX(NEL,XPMINE,XPMIN,DUMMY)
      CALL MINMAX(NEL,YPMINE,YPMIN,DUMMY)
      CALL MINMAX(NEL,XPMAXE,DUMMY,XPMAX)
      CALL MINMAX(NEL,YPMAXE,DUMMY,YPMAX)
C
C---- we now have a new geometry... invalidate any existing solution
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
C
      LGSAME = .TRUE.
C
      RETURN
      END ! PANCOP



      SUBROUTINE PANGEN(LQUERY)
C------------------------------------------------------------------
C     Generates panel nodes X,Y from buffer-airfoil XB,YB arrays.
C     If LQUERY=T, then user's interaction is requested.
C------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LQUERY
C
      LOGICAL LCPLOT
      DIMENSION SG(IPX)
C
      DATA LCPLOT / .FALSE. /
C
      IF(NBEL.EQ.0) RETURN
      NEL = NBEL
C
C---- process each element
      DO 100 IEL = 1, NBEL
C
        IB1 = IBFRST(IEL)
        IB2 = IBLAST(IEL)
        NBPTS = IB2 - IB1 + 1
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        NPPTS = IP2 - IP1 + 1
        IF(LDBG) THEN
          WRITE(*,*) 'PANGEN IEL ',IEL
          WRITE(*,*) 'IB1,IB2 ',IB1,IB2
          WRITE(*,*) 'IP1,IP2 ',IP1,IP2
        ENDIF
C
        IF(NPPTS.LE.1 .AND. NBPTS.LE.1) THEN
C------- don't try to panel point or ring element
         NG = 1
         SG(1) = 0.
C
        ELSE
C
         IF(LQUERY) THEN
C-------- set current spacing so SGCURV has something to plot
           IF(NPPTS.GT.1) THEN
C--------- set SG from existing paneled geometry, if defined
            NG = NPPTS
            DO IP = IP1, IP2
              IG = IP - IP1 + 1
              SG(IG) = (SP(IP)-SP(IP1))/(SP(IP2)-SP(IP1))
            ENDDO
           ELSEIF(NBPTS.GT.1) THEN
C--------- or, set SG from buffer geometry
            NG = NBPTS
            DO IB = IB1, IB2
              IG = IB - IB1 + 1
              SG(IG) = (SB(IB)-SB(IB1))/(SB(IB2)-SB(IB1))
            ENDDO
            IF(NPAN(IEL).LE.0) NPAN(IEL) = NBPTS
           ELSE
            WRITE(*,*)
            WRITE(*,*) '  No current paneling to modify'
            WRITE(*,*) '  Input case and/or Execute PANE command'
            RETURN
           ENDIF
C
         ELSE
C--------- no initial paneling for SGCURV
           SG(1) = 0.
         ENDIF
C
c         WRITE(*,1150) IEL
c 1150    FORMAT(
c     & /' ============================================================='
c     &//'  Element', I4)
C
C------- surface element... generate spacing array SG 
C
         CALL SGCURV(LQUERY,LCPLOT,NBPTS,
     &               XB(IB1),XBS(IB1),
     &               YB(IB1),YBS(IB1),
     &               SB(IB1),
     &               IPX, NG,SG,
     &               NPAN(IEL),CVEX(IEL),SMOF(IEL),FSLE(IEL),FSTE(IEL),
     &               NRPNX,SRPN1(1,IEL),SRPN2(1,IEL),CRRAT(1,IEL))
        ENDIF
C
        IF(LQUERY) THEN
C---- Set flag for spacing defined, unset user-defined flag
         LRSPCDEF = .TRUE.
         LRSPCUSR = .TRUE.
        ENDIF       
C
C------ starting and ending indices for this element IEL
        IF(IEL.EQ.1) THEN
         IP1 = 1
         IP2 = NG
        ELSE
         IP1 = IPLAST(IEL-1) + 1
         IP2 = IPLAST(IEL-1) + NG
        ENDIF
        IF(IP2.GT.IPX) STOP 'PANGEN: Array overflow on IPX'
C
        IF(LQUERY) THEN
          IPDEL = IP2 - IPLAST(IEL)
          IF(NPTOT+IPDEL.GT.IPX) STOP 'PANGEN: Array overflow on IPX'
C
C-------- move current points in subsequent elements to make/take-up room
          IF(IPDEL.GT.0) THEN
           KP1 = NPTOT
           KP2 = IPLAST(IEL)+1
           KPD = -1
          ELSE
           KP1 = IPLAST(IEL)+1
           KP2 = NPTOT
           KPD = +1
          ENDIF
          DO KP = KP1, KP2, KPD
            SP(KP+IPDEL) = SP(KP)
            XP(KP+IPDEL) = XP(KP)
            YP(KP+IPDEL) = YP(KP)
            XPS(KP+IPDEL) = XPS(KP)
            YPS(KP+IPDEL) = YPS(KP)
          ENDDO
          DO KEL = IEL+1, NEL
            IPFRST(KEL) = IPFRST(KEL) + IPDEL
            IPLAST(KEL) = IPLAST(KEL) + IPDEL
          ENDDO
        ENDIF
C
C------ set panel nodes at fractional arc length locations SG
        DO IP = IP1, IP2
          IG = IP - IP1 + 1
          SBI = SB(IB1) + (SB(IB2)-SB(IB1))*SG(IG)
C
          XP(IP) = SEVAL(SBI,XB(IB1),XBS(IB1),SB(IB1),NBPTS)
          YP(IP) = SEVAL(SBI,YB(IB1),YBS(IB1),SB(IB1),NBPTS)
        ENDDO
        IPFRST(IEL) = IP1
        IPLAST(IEL) = IP2
C
C------ centroid location
        XPCENT(IEL) = XBCEN2DT(IEL)
        YPCENT(IEL) = YBCEN2DT(IEL)

C------ spline and set other stuff for element IEL
        CALL XYPSPL(IEL)
C
C------ set indices of TE panel, if any
        CALL ITPSET(IEL)
C
        I = IPFRST(IEL)
        N = IPLAST(IEL) - IPFRST(IEL) + 1
        CALL MINMAX(N,XP(I),XPMINE(IEL),XPMAXE(IEL))
        CALL MINMAX(N,YP(I),YPMINE(IEL),YPMAXE(IEL))
 100  CONTINUE
C
C---- set total number of panel vertices
      NPTOT = IPLAST(NEL)
C
C---- set overall min,max over all elements
      CALL MINMAX(NEL,XPMINE,XPMIN,DUMMY)
      CALL MINMAX(NEL,YPMINE,YPMIN,DUMMY)
      CALL MINMAX(NEL,XPMAXE,DUMMY,XPMAX)
      CALL MINMAX(NEL,YPMAXE,DUMMY,YPMAX)
C
C---- default reference point location for each element
      DO IEL = 1, NEL
        XPREFE(IEL) = 0.0
        YPREFE(IEL) = 0.0
      ENDDO
C
C---- we now have a new geometry... invalidate any existing solution
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
C
      LQSPEC = .FALSE.
      LGSAME = .FALSE.
C
      RETURN
      END ! PANGEN



      SUBROUTINE PANWRT
C-----------------------------------------
C     Writes paneling-parameter file
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*80 PFNDEF
C
 1000 FORMAT(A)
C
      LU = 8
C
C---- use existing paneling-parameter filename (if any) as the default
 20   WRITE(*,1200) PFILE(1:60)
 1200 FORMAT(/' Enter output filename:  ', A)
C
      READ(*,1000) PFNDEF
      IF(PFNDEF(1:1).EQ.' ') THEN
        IF(PFILE(1:1).EQ.' ') GO TO 20
      ELSE
        PFILE = PFNDEF
      ENDIF
C
      OPEN(LU,FILE=PFILE,STATUS='UNKNOWN',ERR=20)
      REWIND(LU)
C
      WRITE(LU,*) NEL, NRPNX
      DO IEL = 1, NEL
        WRITE(LU,*) NPAN(IEL)
        WRITE(LU,*) CVEX(IEL), SMOF(IEL), FSLE(IEL), FSTE(IEL)
        DO IRPN = 1, NRPNX
          WRITE(LU,*) SRPN1(IRPN,IEL),
     &                SRPN2(IRPN,IEL),
     &                CRRAT(IRPN,IEL)
        ENDDO
      ENDDO
      CLOSE(LU)
C
      RETURN
      END ! PANWRT



      SUBROUTINE PANGET(FNAME,ERROR)
C-----------------------------------------
C     Reads paneling-parameter file FNAME
C-----------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*) FNAME
      LOGICAL ERROR
C
      LU = 8
C
      IF(FNAME(1:1).EQ.' ') THEN
       CALL ASKS('Enter paneling-parameter filename^',FNAME)
      ENDIF
C
C---- strip off leading blanks and get number of characters
      CALL STRIP(FNAME,NFNAME)
C
      WRITE(*,*)
C
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=90)
      REWIND(LU)
C
      READ(LU,*,ERR=80) NEL0, NRPN0
      NEL0  = MIN( NEL0  , NEX   )
      NRPN0 = MIN( NRPN0 , NRPNX )
C
      DO IEL = 1, NEL0
        READ(LU,*,ERR=80) NPAN(IEL)
        READ(LU,*,ERR=80) CVEX(IEL), SMOF(IEL), FSLE(IEL), FSTE(IEL)
        DO IRPN = 1, NRPN0
          READ(LU,*,ERR=80) SRPN1(IRPN,IEL),
     &                      SRPN2(IRPN,IEL),
     &                      CRRAT(IRPN,IEL)
        ENDDO
      ENDDO
      CLOSE(LU)
      WRITE(*,*) 'Paneling parameters read from file ', FNAME(1:NFNAME)
      LRSPCDEF = .TRUE.
      LRSPCUSR = .TRUE.
      ERROR = .FALSE.
      RETURN
C
 80   CONTINUE
      CLOSE(LU)
      WRITE(*,*) 'READ error on panel-parameter file ', FNAME(1:NFNAME)
      ERROR = .TRUE.
      RETURN
C
 90   CONTINUE
c      WRITE(*,*) 'OPEN error on panel-parameter file ', FNAME(1:NFNAME)
      ERROR = .TRUE.
      RETURN
      END ! PANGET
