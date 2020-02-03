module m_geom
    implicit none
contains
    !*==XYPSPL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! BAKSUB
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
    
    SUBROUTINE XYPSPL(IEL)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IEL
        !
        ! Local variables
        !
        INTEGER :: I, IP1, IP2, N
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------
        !     Splines panel-node coordinates
        !     and sets other element-related stuff.
        !---------------------------------------------
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        !
        I = IP1
        N = IP2 - IP1 + 1
        IF (N>1) THEN
            CALL SCALC(XP(I), YP(I), SP(I), N)
            CALL SEGSPL(XP(I), XPS(I), SP(I), N)
            CALL SEGSPL(YP(I), YPS(I), SP(I), N)
        ENDIF
        !
        IF (LXBOD(IEL)) THEN
            !----- axisymmetric body on axis
            IF (XP(IP1)<XP(IP2)) THEN
                SPLE(IEL) = SP(IP1)
                XPLE(IEL) = XP(IP1)
                YPLE(IEL) = YP(IP1)
                XPTE(IEL) = XP(IP2)
                YPTE(IEL) = YP(IP2)
            ELSE
                SPLE(IEL) = SP(IP2)
                XPLE(IEL) = XP(IP2)
                YPLE(IEL) = YP(IP2)
                XPTE(IEL) = XP(IP1)
                YPTE(IEL) = YP(IP1)
            ENDIF
            !
        ELSEIF (LBODY(IEL) .AND. (LV1ZR(IEL) .OR. LV2ZR(IEL))) THEN
            !----- axisymmetric body not closed on axis
            IF (XP(IP1)<XP(IP2)) THEN
                SPLE(IEL) = SP(IP1)
                XPLE(IEL) = XP(IP1)
                YPLE(IEL) = YP(IP1)
                XPTE(IEL) = XP(IP2)
                YPTE(IEL) = YP(IP2)
            ELSE
                SPLE(IEL) = SP(IP2)
                XPLE(IEL) = XP(IP2)
                YPLE(IEL) = YP(IP2)
                XPTE(IEL) = XP(IP1)
                YPTE(IEL) = YP(IP1)
            ENDIF
            !
        ELSEIF (LBODY(IEL) .AND. .NOT.LXBOD(IEL)) THEN
            !----- body off axis
            CALL LEFIND(SPLE(IEL), XP(I), XPS(I), YP(I), YPS(I), SP(I), N)
            XPLE(IEL) = SEVAL(SPLE(IEL), XP(I), XPS(I), SP(I), N)
            YPLE(IEL) = SEVAL(SPLE(IEL), YP(I), YPS(I), SP(I), N)
            XPTE(IEL) = 0.5 * (XP(IP1) + XP(IP2))
            YPTE(IEL) = 0.5 * (YP(IP1) + YP(IP2))
            !
        ELSE
            !----- zero-thickness body
            SPLE(IEL) = SP(IP1)
            XPLE(IEL) = XP(IP1)
            YPLE(IEL) = YP(IP1)
            XPTE(IEL) = XP(IP2)
            YPTE(IEL) = YP(IP2)
        ENDIF
        !
        !---- set location for plotted element index
        IF (NETYPE(IEL)==0 .OR. NETYPE(IEL)==1 .OR. NETYPE(IEL)==2) THEN
            !----- surface or axis line... set plot location at skin centroid
            XELNUM(IEL) = XPCENT(IEL)
            YELNUM(IEL) = YPCENT(IEL)
        ELSE
            !----- point,ring,source line or vortex wake element...
            !      just set plot location at the point itself
            XELNUM(IEL) = XP(IP1)
            YELNUM(IEL) = YP(IP1)
        ENDIF
        !
    END SUBROUTINE XYPSPL
    !*==CVPGEN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! XYPSPL
    
    
    
    
    SUBROUTINE CVPGEN
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        REAL :: DSP, DXP, DYP
        INTEGER :: IC, IEL, IP
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Sets control-point locations on paneled airfoil.
        !     Sets locations of point sources.
        !---------------------------------------------------------
        !
        IF (LDBG) WRITE (*, *) 'Setting control points for elements'
        !
        DO IEL = 1, NEX
            ICFRST(IEL) = 0
            ICLAST(IEL) = 0
        ENDDO
        !
        DO IC = 1, ICX
            IPCO(IC) = 0
            IPCP(IC) = 0
        ENDDO
        !
        !
        !---- initialize counters for pointer accumulation
        IC = 0
        !
        DO IEL = 1, NEL
            IF (NETYPE(IEL)==3 .OR. NETYPE(IEL)==4) THEN
                !------- ring or point singularity... no control points
                IP = IPFRST(IEL)
                !
                !------- set limits so as to skip over IC do-loops
                ICFRST(IEL) = 0
                ICLAST(IEL) = -1
                !
                CYCLE
            ENDIF
            !
            !------ first control point in element IEL
            ICFRST(IEL) = IC + 1
            !
            !------ go over all panels on element IEL
            DO IP = IPFRST(IEL), IPLAST(IEL) - 1
                IC = IC + 1
                IF (IC>NCX) STOP 'CVPGEN: Array overflow on NCX'
                !
                DXP = XP(IP + 1) - XP(IP)
                DYP = YP(IP + 1) - YP(IP)
                DSP = SQRT(DXP**2 + DYP**2)
                !
                IF (DSP==0.0) THEN
                    !--------- zero-length panel gets dummy control point
                    ICTYPE(IC) = -1
                    !
                ELSE
                    !--------- ordinary panel
                    ICTYPE(IC) = 0
                    IF (IEL2IR(IEL)==1) THEN
                        ICTYPE(IC) = 0
                    ELSEIF (IEL2IR(IEL)>1 .AND. IEL2IR(IEL)<NRP) THEN
                        !c             ICTYPE(IC) = -2
                    ELSEIF (IEL2IR(IEL)==NRP) THEN
                        ICTYPE(IC) = 0
                    ENDIF
                    !
                ENDIF
                !
                !-------- panel nodes defining sheet strength at IC
                IPCO(IC) = IP
                IPCP(IC) = IP + 1
            ENDDO
            !
            !------ last control point in element IEL
            ICLAST(IEL) = IC
            !
        ENDDO
        !
        !---- total number of control points
        NCTOT = IC
        !
        !---- compute normal vectors at vortex nodes, compute geometric quantities
        DO IEL = 1, NEL
            !------ don't do point singularities
            IF (NETYPE(IEL)==0 .OR. NETYPE(IEL)==1 .OR. NETYPE(IEL)       &
                    & ==2 .OR. NETYPE(IEL)==5 .OR. NETYPE(IEL)==6 .OR.          &
                    & NETYPE(IEL)==7) THEN
                CALL XYCSET(IEL)
                CALL ANPSET(IEL)
            ENDIF
        ENDDO
        !
        !---- invalidate any existing solution
        LQAIC = .FALSE.
        LSYSP = .FALSE.
        LGSYS = .FALSE.
        LGAMU = .FALSE.
        LGAMA = .FALSE.
        !
        LNCVP = .TRUE.
        !
    END SUBROUTINE CVPGEN
    !*==XYCSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CVPGEN
    
    
    
    SUBROUTINE XYCSET(IEL)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IEL
        !
        ! Local variables
        !
        REAL, SAVE :: DCFRAC
        REAL :: DSCA, DSP, DSPINV, DXP, DYP, STAN, XTAN, YTAN
        INTEGER :: IC, IC1, IC2, ICE, IP, IP1, IP2, IPO, IPP
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Sets control points and their unit normal vectors
        !     for element IEL.
        !----------------------------------------------------------------
        !
        !c    DATA DCFRAC / 0.1 /
        DATA DCFRAC/0.05/
        !c    DATA DCFRAC / 0.02 /
        !
        DO IC = ICFRST(IEL), ICLAST(IEL)
            IPO = IPCO(IC)
            IPP = IPCP(IC)
            !
            DXP = XP(IPP) - XP(IPO)
            DYP = YP(IPP) - YP(IPO)
            DSP = SQRT(DXP**2 + DYP**2)
            !
            IF (DSP==0.0) THEN
                !------- zero-length panel
                DSPINV = 0.0
            ELSE
                !------- ordinary panel
                DSPINV = 1.0 / DSP
            ENDIF
            !
            !------ assign control point at center of panel
            XC(IC) = XP(IPO) + 0.5 * DXP
            YC(IC) = YP(IPO) + 0.5 * DYP
            DSC(IC) = DSP
            DSC_DXY(1, IC) = DXP * DSPINV
            DSC_DXY(2, IC) = DYP * DSPINV
            !
            !------ normal vector at control point
            ANC(1, IC) = DYP * DSPINV
            ANC(2, IC) = -DXP * DSPINV
            !
            !------ geometry sensitivities of normal vector for design/optimization stuff
            ANC_DXY(1, 1, IC) = ANC(1, IC) * ANC(2, IC) * DSPINV
            ANC_DXY(1, 2, IC) = ANC(2, IC)**2 * DSPINV
            ANC_DXY(2, 1, IC) = -ANC(1, IC)**2 * DSPINV
            ANC_DXY(2, 2, IC) = -ANC(2, IC) * ANC(1, IC) * DSPINV
        ENDDO
        !
        !
        !------ set up extra QNDOF control points for each body element
        IF (.NOT.LBODY(IEL)) THEN
            ICE = NCX + IEL
            XC(ICE) = 0.
            YC(ICE) = 0.
            !
        ELSEIF (LXBOD(IEL)) THEN
            IC = ICLAST(IEL)
            IP = IPLAST(IEL)
            !
            ICE = NCX + IEL
            XC(ICE) = XP(IP) - DCFRAC * DSC(IC)
            YC(ICE) = 0.
            ANC(1, ICE) = 1.0
            ANC(2, ICE) = 0.
            !
        ELSE
            IC1 = ICFRST(IEL)
            IC2 = ICLAST(IEL)
            !
            IP1 = IPFRST(IEL)
            IP2 = IPLAST(IEL)
            !
            XTAN = XP(IP1) - XP(IP1 + 1) + XP(IP2) - XP(IP2 - 1)
            YTAN = YP(IP1) - YP(IP1 + 1) + YP(IP2) - YP(IP2 - 1)
            STAN = SQRT(XTAN**2 + YTAN**2)
            !
            DSCA = 0.5 * (DSC(IC1) + DSC(IC2))
            !
            ICE = NCX + IEL
            XC(ICE) = 0.5 * (XP(IP1) + XP(IP2)) - DCFRAC * DSCA * XTAN / STAN
            YC(ICE) = 0.5 * (YP(IP1) + YP(IP2)) - DCFRAC * DSCA * YTAN / STAN
            ANC(1, ICE) = XTAN / STAN
            ANC(2, ICE) = YTAN / STAN
            !
        ENDIF
        !
    END SUBROUTINE XYCSET
    !*==ANPSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! XYCSET
    
    
    
    SUBROUTINE ANPSET(IEL)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IEL
        !
        ! Local variables
        !
        REAL :: DS, DSQ, DSQM, DSQP, DX, DXM, DXP, DX_DXM, &
                & DX_DXP, DX_DYM, DX_DYP, DY, DYM, DYP, DY_DXM, &
                & DY_DXP, DY_DYM, DY_DYP, TMP1, TMP2
        INTEGER :: IP, IPM, IPP
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Sets unit normal vectors on vortex node locations
        !     for element IEL.
        !----------------------------------------------------------------
        !
        !---- don't process single-point elements
        IF (IPFRST(IEL)==IPLAST(IEL)) RETURN
        !
        DO IP = IPFRST(IEL), IPLAST(IEL)
            !
            !------ examine the two panels IPM..IP, IP..IPP adjoining vortex node IP
            IPM = MAX(IP - 1, IPFRST(IEL))
            IPP = MIN(IP + 1, IPLAST(IEL))
            !
            DXM = XP(IP) - XP(IPM)
            DYM = YP(IP) - YP(IPM)
            !
            DXP = XP(IPP) - XP(IP)
            DYP = YP(IPP) - YP(IP)
            !
            DSQM = DXM * DXM + DYM * DYM
            DSQP = DXP * DXP + DYP * DYP
            !
            IF (DSQM==0.0) THEN
                !------- first point on element, or panel IPM..IP has zero length
                DX = DXP
                DY = DYP
                DX_DXM = 0.
                DX_DYM = 0.
                DX_DXP = 1.0
                DX_DYP = 0.
                DY_DXM = 0.
                DY_DYM = 0.
                DY_DXP = 0.
                DY_DYP = 1.0
            ELSEIF (DSQP==0.0) THEN
                !------- last point on element, or panel IP..IPP has zero length
                DX = DXM
                DY = DYM
                DX_DXM = 1.0
                DX_DYM = 0.
                DX_DXP = 0.
                DX_DYP = 0.
                DY_DXM = 0.
                DY_DYM = 1.0
                DY_DXP = 0.
                DY_DYP = 0.
            ELSE
                !------- usual interior point... normal vector is weighted average
                DX = DXM * DSQP + DXP * DSQM
                DY = DYM * DSQP + DYP * DSQM
                DX_DXM = DSQP + DXP * 2.0 * DXM
                DX_DYM = +DXP * 2.0 * DYM
                DX_DXP = DSQM + DXM * 2.0 * DXP
                DX_DYP = +DXM * 2.0 * DYP
                DY_DXM = +DYP * 2.0 * DXM
                DY_DYM = DSQP + DYP * 2.0 * DYM
                DY_DXP = +DYM * 2.0 * DXP
                DY_DYP = DSQM + DYM * 2.0 * DYP
            ENDIF
            !
            DSQ = DX * DX + DY * DY
            IF (DSQ==0.0) THEN
                !------- internal error
                WRITE (*, *) '? ANPSET: zero avg panel length. IEL IP =', &
                        & IEL, IP
                RETURN
            ENDIF
            !
            !------ set unit normal ANP
            DS = SQRT(DSQ)
            ANP(1, IP) = DY / DS
            ANP(2, IP) = -DX / DS
            !
            TMP1 = -ANP(1, IP) / DSQ
            ANP_XYM(1, 1, IP) = -DY_DXM / DS - TMP1 * (DX * DX_DXM + DY * DY_DXM)
            ANP_XYM(1, 2, IP) = -DY_DYM / DS - TMP1 * (DX * DX_DYM + DY * DY_DYM)
            !
            ANP_XYO(1, 1, IP) = DY_DXM / DS + TMP1 * (DX * DX_DXM + DY * DY_DXM)       &
                    & - DY_DXP / DS - TMP1 * (DX * DX_DXP + DY * DY_DXP)
            ANP_XYO(1, 2, IP) = DY_DYM / DS + TMP1 * (DX * DX_DYM + DY * DY_DYM)       &
                    & - DY_DYP / DS - TMP1 * (DX * DX_DYP + DY * DY_DYP)
            !
            ANP_XYP(1, 1, IP) = DY_DXP / DS + TMP1 * (DX * DX_DXP + DY * DY_DXP)
            ANP_XYP(1, 2, IP) = DY_DYP / DS + TMP1 * (DX * DX_DYP + DY * DY_DYP)
            !
            TMP2 = -ANP(2, IP) / DSQ
            ANP_XYM(2, 1, IP) = DX_DXM / DS - TMP2 * (DX * DX_DXM + DY * DY_DXM)
            ANP_XYM(2, 2, IP) = DX_DYM / DS - TMP2 * (DX * DX_DYM + DY * DY_DYM)
            !
            ANP_XYO(2, 1, IP) = -DX_DXM / DS + TMP2 * (DX * DX_DXM + DY * DY_DXM)      &
                    & + DX_DXP / DS - TMP2 * (DX * DX_DXP + DY * DY_DXP)
            ANP_XYO(2, 2, IP) = -DX_DYM / DS + TMP2 * (DX * DX_DYM + DY * DY_DYM)      &
                    & + DX_DYP / DS - TMP2 * (DX * DX_DYP + DY * DY_DYP)
            !
            ANP_XYP(2, 1, IP) = -DX_DXP / DS + TMP2 * (DX * DX_DXP + DY * DY_DXP)
            ANP_XYP(2, 2, IP) = -DX_DYP / DS + TMP2 * (DX * DX_DYP + DY * DY_DYP)
        ENDDO
        !
    END SUBROUTINE ANPSET
    !*==ITPSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ANPSET
    
    
    SUBROUTINE ITPSET(IEL)
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IEL
        !
        ! Local variables
        !
        INTEGER :: IP1, IP2
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------
        !     Sets TE panel endpoint indices
        !--------------------------------------
        !
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        !
        !---- default case... no TE panel
        IPTE1(IEL) = -999
        IPTE2(IEL) = -999
        !
        !---- no TE panel on closed axis body or non-solid body
        IF (LXBOD(IEL) .OR. NETYPE(IEL)/=0) RETURN
        !
        !---- TE panel must be chosen...
        IF (LTPAN(IEL)) THEN
            !----- usual TE panel closing off body
            IPTE1(IEL) = IP1
            IPTE2(IEL) = IP2
            IF (YP(IP1)==0.0 .AND. YP(IP2)>0.0) THEN
                IPTE1(IEL) = 0
                IPTE2(IEL) = IP2
            ELSEIF (YP(IP2)==0.0 .AND. YP(IP1)>0.0) THEN
                IPTE1(IEL) = IP1
                IPTE2(IEL) = 0
            ENDIF
        ENDIF
        !
    END SUBROUTINE ITPSET
end module m_geom
