module m_vels
    implicit none
contains
    !*==QCPFOR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    
    
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
    
    
    SUBROUTINE QCPFOR
        use m_forces, only: cpcalc, fcoeff
        use m_system, only: qcsum
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        REAL, DIMENSION(2) :: ANT
        REAL :: CIRC, CMTE, COSA, CPTL, CPTR, CXTE, CYTE, &
                & DELTACP, DSCT, DXT, DYT, QRSQ, QUE, SINA, VTSQ, &
                & VTT, XCT, YCT
        INTEGER :: I, IC, ICT1, ICT2, IEL, IG, IP1, IP2, IPT1, &
                & IPT2, IR, IR1, IR2, IU, JG, N
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------
        !     Computes velocities, Cp, forces, etc.
        !     for current singularity distribitions
        !---------------------------------------------
        !
        !---- compute velocities at control points (average between two panel sides)
        IF (.NOT.LQCNT) CALL QCSUM
        !
        !---- set velocities on both sides of panels at all control points
        CALL QQCALC(NCTOT, ANC, IPCO, IPCP, GAM, SIG, GTH, QC, QCL, QCR)
        !
        !---- also set QCL,QCR sensitivities to imposed freestream, singularities, etc.
        DO IU = 1, NU
            CALL QQCALC(NCTOT, ANC, IPCO, IPCP, GAMU(1, IU), SIGU(1, IU), &
                    & GTHU(1, IU), QCU(1, 1, IU), QCLU(1, 1, IU), QCRU(1, 1, IU))
        ENDDO
        !
        !---- set velocities and Cp on both sides of panels at all control points
        CALL CPCALC(NCTOT, QCL, QINF, QREF, CPL, CPL_QCL, CPL_QINF)
        CALL CPCALC(NCTOT, QCR, QINF, QREF, CPR, CPR_QCR, CPR_QINF)
        !
        !---- add in contribution from rotor enthalpy and entropy losses to Cp's
        !     CB and duct walls and on upper and lower wakes
        QRSQ = QREF**2
        DO IC = 1, NCTOT
            IP1 = IPCO(IC)
            IP2 = IPCP(IC)
            IR1 = IP2IR(IP1)
            IR2 = IP2IR(IP2)
            IG = IC2IG(IC)
            !---- check for control point in rotor slipstream (in grid)
            IF (IG>0) THEN
                IR = IR1
                IF (IR==1) THEN
                    JG = IR
                    IF (IP2<=IPLAST(1)) THEN
                        !---- on CB downstream of rotor
                        CIRC = BGAMG(IG, JG)
                        VTT = CIRC * PI2I / YC(IC)
                        VTSQ = VTT**2
                        DELTACP = 2.0 * (DHG(IG, JG) - DSG(IG, JG)) / QRSQ - VTSQ / QRSQ
                        CPR(IC) = CPR(IC) + DELTACP
                    ELSE
                        !---- on CB wake element
                        CIRC = BGAMG(IG, JG)
                        VTT = CIRC * PI2I / YC(IC)
                        VTSQ = VTT**2
                        DELTACP = 2.0 * (DHG(IG, JG) - DSG(IG, JG)) / QRSQ - VTSQ / QRSQ
                        CPL(IC) = CPL(IC) + DELTACP
                    ENDIF
                    !
                ELSEIF (IR==NRP) THEN
                    JG = IR - 1
                    IF (IP1<IPLAST(2)) THEN
                        !---- on DUCT downstream of rotor
                        CIRC = BGAMG(IG, JG)
                        VTT = CIRC * PI2I / YC(IC)
                        VTSQ = VTT**2
                        DELTACP = 2.0 * (DHG(IG, JG) - DSG(IG, JG)) / QRSQ - VTSQ / QRSQ
                        CPR(IC) = CPR(IC) + DELTACP
                    ELSE
                        !---- on DUCT wake element
                        CIRC = BGAMG(IG, JG)
                        VTT = CIRC * PI2I / YC(IC)
                        VTSQ = VTT**2
                        DELTACP = 2.0 * (DHG(IG, JG) - DSG(IG, JG)) / QRSQ - VTSQ / QRSQ
                        CPR(IC) = CPR(IC) + DELTACP
                    ENDIF
                ELSE
                    JG = IR - 1
                    CIRC = BGAMG(IG, JG)
                    VTT = CIRC * PI2I / YC(IC)
                    VTSQ = VTT**2
                    DELTACP = 2.0 * (DHG(IG, JG) - DSG(IG, JG)) / QRSQ - VTSQ / QRSQ
                    CPR(IC) = CPR(IC) + DELTACP
                    !
                    JG = IR
                    CIRC = BGAMG(IG, JG)
                    VTT = CIRC * PI2I / YC(IC)
                    VTSQ = VTT**2
                    DELTACP = 2.0 * (DHG(IG, JG) - DSG(IG, JG)) / QRSQ - VTSQ / QRSQ
                    CPL(IC) = CPL(IC) + DELTACP
                ENDIF
            ENDIF
        ENDDO
        !
        !---- also set CPL,CPR sensitivities
        DO IU = 1, NU
            DO IC = 1, NCTOT
                CPLU(IC, IU) = CPL_QCL(1, IC) * QCLU(1, IC, IU) + CPL_QCL(2, IC)   &
                        & * QCLU(2, IC, IU)
                CPRU(IC, IU) = CPR_QCR(1, IC) * QCRU(1, IC, IU) + CPR_QCR(2, IC)   &
                        & * QCRU(2, IC, IU)
            ENDDO
        ENDDO
        !
        !---- integrate Cp to get forces and moment, and all sensitivities
        COSA = 1.0
        SINA = 0.0
        CX(0) = 0.
        CY(0) = 0.
        CM(0) = 0.
        CD(0) = 0.
        CXVIS(0) = 0.
        DO IU = 1, NU
            CXU(0, IU) = 0.
            CYU(0, IU) = 0.
            CMU(0, IU) = 0.
            CDU(0, IU) = 0.
        ENDDO
        DO IEL = 1, NEL
            I = ICFRST(IEL)
            N = ICLAST(IEL) - ICFRST(IEL) + 1
            IF (N>=1) THEN
                CALL FCOEFF(N, CPL(I), CPR(I), XC(I), YC(I), ANC(1, I), DSC(I), &
                        & CX(IEL), CY(IEL), CM(IEL))
                !
                CD(IEL) = CX(IEL)
                !
                !---- Totals go in index 0
                CX(0) = CX(0) + CX(IEL)
                !         CY(0) = CY(0) + CY(IEL)
                !         CM(0) = CM(0) + CM(IEL)
                !         CD(0) = CD(0) + CD(IEL)
                !
                !        cosa = u / sqrt(u**2 + v**2)
                !        sina = v / sqrt(u**2 + v**2)
                !        cosa_u =  v**2 / q**3 = sina**2/q
                !        cosa_v =  -u*v / q**3 = -sina cosa/q
                !        sina_u =  -v*u / q**3 = -sina cosa/q
                !        sina_v =  u**2 / q**3 = cosa**2/q
                !
                !------- also set sensitivities due to freestream
                IU = IUQINF
                IF (LUSET(IU)) THEN
                    CALL FCOEFF(N, CPLU(I, IU), CPRU(I, IU), XC(I), YC(I), ANC(1, I), &
                            & DSC(I), CXU(IEL, IU), CYU(IEL, IU), CMU(IEL, IU))
                    CDU(IEL, IU) = CXU(IEL, IU)
                    !---- Totals go in index 0
                    CXU(0, IU) = CXU(0, IU) + CXU(IEL, IU)
                    !           CYU(0,IU) = CYU(0,IU) + CYU(IEL,IU)
                    !           CMU(0,IU) = CMU(0,IU) + CMU(IEL,IU)
                    !           CDU(0,IU) = CDU(0,IU) + CDU(IEL,IU)
                ENDIF
                !
                !---- Add in viscous forces if BL has been calculated
                IF (LVISC) THEN
                    CX(IEL) = CX(IEL) + CXVIS(IEL)
                    CD(IEL) = CX(IEL)
                    !---- Add viscous forces to total forces in index 0
                    CX(0) = CX(0) + CXVIS(IEL)
                    CXVIS(0) = CXVIS(0) + CXVIS(IEL)
                ENDIF
                !
                !
                !---- Elements with TE panels
                !     Add forces due to pressures on TE panel
                IF (LTPAN(IEL) .AND. NETYPE(IEL)==0) THEN
                    IPT1 = IPTE1(IEL)
                    IPT2 = IPTE2(IEL)
                    ICT1 = 0
                    ICT2 = 0
                    IF (IPT1==IPFRST(IEL)) ICT1 = ICFRST(IEL)
                    IF (IPT1==IPLAST(IEL)) ICT1 = ICLAST(IEL)
                    IF (IPT2==IPFRST(IEL)) ICT2 = ICFRST(IEL)
                    IF (IPT2==IPLAST(IEL)) ICT2 = ICLAST(IEL)
                    !          write(12,*) 'iel ',iel
                    !          write(12,*) 'ipt1,ipt2 ',ipt1,ipt2
                    !          write(12,*) 'ipfrst,iplast ',ipfrst(iel),iplast(iel)
                    !          write(12,*) 'ict1,ict2 ',ict1,ict2
                    !          write(12,*) 'icfrst,iclast ',icfrst(iel),iclast(iel)
                    !
                    IF (ICT1/=0 .AND. ICT2/=0) THEN
                        CPTL = 0.5 * (CPL(ICT1) + CPL(ICT2))
                        CPTR = 0.5 * (CPR(ICT1) + CPR(ICT2))
                    ELSEIF (ICT1==0) THEN
                        CPTL = CPL(ICT2)
                        CPTR = CPR(ICT2)
                    ELSEIF (ICT2==0) THEN
                        CPTL = CPL(ICT1)
                        CPTR = CPR(ICT1)
                    ENDIF
                    XCT = 0.5 * (XPT(1, IEL) + XPT(2, IEL))
                    YCT = 0.5 * (YPT(1, IEL) + YPT(2, IEL))
                    DXT = XPT(1, IEL) - XPT(2, IEL)
                    DYT = YPT(1, IEL) - YPT(2, IEL)
                    DSCT = SQRT(DXT * DXT + DYT * DYT)
                    ANT(1) = DYT / DSCT
                    ANT(2) = -DXT / DSCT
                    !          write(12,*) 'cptl,cptr ',cptl,cptr
                    !          write(12,*) 'xct,yct ',xct,yct
                    !          write(12,*) 'dsct ',dsct
                    !          write(12,*) 'ant ',ant(1),ant(2)
                    !
                    CALL FCOEFF(1, CPTL, CPTR, XCT, YCT, ANT, DSCT, CXTE, CYTE, CMTE)
                    !          write(12,*) 'CXTE ',cxte
                    !          write(12,*) 'CYTE ',cyte
                    !          write(12,*) 'CMTE ',cmte
                    !
                    !---- Add in TE forces
                    CX(IEL) = CX(IEL) + CXTE
                    !          CY(IEL) = CY(IEL) + CYTE
                    !          CD(IEL) = CX(IEL)
                    !---- Add TE forces to total forces in index 0
                    CX(0) = CX(0) + CXTE
                    !          CY(0) = CY(0) + CYTE
                    !          CM(0) = CM(0) + CMTE
                    !          CD(0) = CD(0) + CXTE
    
                    !------- also set sensitivities due to freestream
                    IU = IUQINF
                    IF (LUSET(IU)) THEN
                        IF (ICT1/=0 .AND. ICT2/=0) THEN
                            CPTL = 0.5 * (CPLU(ICT1, IU) + CPLU(ICT2, IU))
                            CPTR = 0.5 * (CPRU(ICT1, IU) + CPRU(ICT2, IU))
                        ELSEIF (ICT1==0) THEN
                            CPTL = CPLU(ICT2, IU)
                            CPTR = CPRU(ICT2, IU)
                        ELSEIF (ICT2==0) THEN
                            CPTL = CPLU(ICT1, IU)
                            CPTR = CPRU(ICT1, IU)
                        ENDIF
                        CALL FCOEFF(1, CPTL, CPTR, XCT, YCT, ANT, DSCT, CXTE, CYTE, &
                                & CMTE)
                        CXU(IEL, IU) = CXU(IEL, IU) + CXTE
                        !            CYU(IEL,IU) = CYU(IEL,IU) + CYTE
                        !            CDU(IEL,IU) = CXU(IEL,IU)
                        !---- Totals go in index 0
                        CXU(0, IU) = CXU(0, IU) + CXTE
                        !            CYU(0,IU) = CYU(0,IU) + CYTE
                        !            CMU(0,IU) = CMU(0,IU) + CMTE
                        !            CDU(0,IU) = CDU(0,IU) + CXTE
                    ENDIF
                ENDIF
                !
            ENDIF
        ENDDO
        !
        !---- Save duct + centerbody (IEL=1,2) X force (thrust)
        QUE = 0.5 * RHO * QREF**2
        TDUCT = -(CX(1) + CX(2)) * QUE
        !
    END SUBROUTINE QCPFOR
    !*==QQCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCPFOR
    
    
    
    SUBROUTINE QQCALC(N, AN, IPO, IPP, GAM, SIG, GTH, Q, QL, QR)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(2, *) :: AN, Q, QL, QR
        REAL, DIMENSION(*) :: GAM, GTH, SIG
        INTEGER, DIMENSION(*) :: IPO, IPP
        !
        ! Local variables
        !
        REAL :: DQN, DQT, GAMC, SIGC
        INTEGER :: I, IO, IP
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Sets velocity on both sides of vortex,source sheet.
        !-------------------------------------------------------------
        !
        DO I = 1, N
            IO = IPO(I)
            IP = IPP(I)
            !
            GAMC = 0.5 * (GAM(IO) + GAM(IP))
            SIGC = 0.5 * (SIG(IO) + SIG(IP))
            !---- Add vorticity contribution from additional vortex wake gamma
            GAMC = GAMC + 0.5 * (GTH(IO) + GTH(IP))
            !
            !------ 1/2 of sheet strength adds,subtracts to left,right speed
            DQT = 0.5 * GAMC
            DQN = 0.5 * SIGC
            !
            QR(1, I) = Q(1, I) + AN(1, I) * DQN + AN(2, I) * DQT
            QR(2, I) = Q(2, I) + AN(2, I) * DQN - AN(1, I) * DQT
            !
            QL(1, I) = Q(1, I) - AN(1, I) * DQN - AN(2, I) * DQT
            QL(2, I) = Q(2, I) - AN(2, I) * DQN + AN(1, I) * DQT
        ENDDO
        !
    END SUBROUTINE QQCALC
    !*==QTCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QQCALC
    
    
    
    SUBROUTINE QTCALC(N, AN, IPO, IPP, GAM, GTH, Q, QT)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(2, *) :: AN, Q
        REAL, DIMENSION(*) :: GAM, GTH, QT
        INTEGER, DIMENSION(*) :: IPO, IPP
        !
        ! Local variables
        !
        REAL, DIMENSION(2) :: AT
        REAL :: DQT, GAMC
        INTEGER :: I, IO, IP
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Sets tangential velocity component
        !-------------------------------------------------------------
        !
        !
        DO I = 1, N
            IO = IPO(I)
            IP = IPP(I)
            !
            GAMC = 0.5 * (GAM(IO) + GAM(IP))
            !---- Add vorticity from vortex wake gamma
            GAMC = GAMC + 0.5 * (GTH(IO) + GTH(IP))
            !
            !------ set tangential vector along sheet (AN  points to right side)
            AT(1) = -AN(2, I)
            AT(2) = AN(1, I)
            !
            !------ 1/2 of sheet strength subtracts to right-side speed (KVSIDE=+1)
            DQT = -0.5 * GAMC
            !
            QT(I) = Q(1, I) * AT(1) + Q(2, I) * AT(2) + DQT
        ENDDO
        !
    END SUBROUTINE QTCALC
    !*==QJAC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QTCALC
    
    
    
    SUBROUTINE QJAC(DX, DY, DU, DV, QJ, QJ_DX, QJ_DY, QJ_DU, QJ_DV)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: DU, DV, DX, DY
        REAL, DIMENSION(2, 2) :: QJ, QJ_DU, QJ_DV, QJ_DX, QJ_DY
        !
        ! Local variables
        !
        REAL :: DSQ, DSQI
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------------
        !     Computes velocity jacobian from position,velocity deltas,
        !     and zero velocity div and curl conditions.
        !------------------------------------------------------------------
        !
        DSQ = DX**2 + DY**2
        IF (DSQ==0.0) THEN
            DSQI = 1.0
        ELSE
            DSQI = 1.0 / DSQ
        ENDIF
        !
        QJ(1, 1) = (DU * DX - DV * DY) * DSQI
        QJ(1, 2) = (DU * DY + DV * DX) * DSQI
        !
        QJ_DX(1, 1) = (DU - 2.0 * QJ(1, 1) * DX) * DSQI
        QJ_DY(1, 1) = (-DV - 2.0 * QJ(1, 1) * DY) * DSQI
        !
        QJ_DX(1, 2) = (DV - 2.0 * QJ(1, 2) * DX) * DSQI
        QJ_DY(1, 2) = (DU - 2.0 * QJ(1, 2) * DY) * DSQI
        !
        QJ_DU(1, 1) = DX * DSQI
        QJ_DV(1, 1) = -DY * DSQI
        !
        QJ_DU(1, 2) = DY * DSQI
        QJ_DV(1, 2) = DX * DSQI
        !
        !
        QJ(2, 1) = QJ(1, 2)
        QJ_DX(2, 1) = QJ_DX(1, 2)
        QJ_DY(2, 1) = QJ_DY(1, 2)
        QJ_DU(2, 1) = QJ_DU(1, 2)
        QJ_DV(2, 1) = QJ_DV(1, 2)
        !
        QJ(2, 2) = -QJ(1, 1)
        QJ_DX(2, 2) = -QJ_DX(1, 1)
        QJ_DY(2, 2) = -QJ_DY(1, 1)
        QJ_DU(2, 2) = -QJ_DU(1, 1)
        QJ_DV(2, 2) = -QJ_DV(1, 1)
        !
    END SUBROUTINE QJAC
    !*==QCURV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QJAC
    
    
    SUBROUTINE QCURV(U, V, QJ, CRV, CRV_U, CRV_V, CRV_QJ)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: CRV, CRV_U, CRV_V, U, V
        REAL, DIMENSION(2, 2) :: CRV_QJ, QJ
        !
        ! Local variables
        !
        REAL :: CRV_QI, CRV_UN, CRV_VN, Q, QI, QI_U, QI_V, QSQ, &
                & UN, UN_U, UN_V, VN, VN_U, VN_V
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Computes streamline curvature from velocity
        !     and velocity Jacobian.
        !-----------------------------------------------------
        !
        QSQ = U**2 + V**2
        Q = SQRT(QSQ)
        !
        IF (Q==0.0) THEN
            QI = 1.0
        ELSE
            QI = 1.0 / Q
        ENDIF
        !
        UN = U * QI
        VN = V * QI
        UN_U = VN * VN * QI
        UN_V = -UN * VN * QI
        VN_U = -UN * VN * QI
        VN_V = UN * UN * QI
        !
        QI_U = -UN * QI**2
        QI_V = -VN * QI**2
        !
        !
        CRV = (UN * UN * QJ(2, 1) - UN * VN * QJ(1, 1) + UN * VN * QJ(2, 2) - VN * VN * QJ(1, 2)) * QI
        !
        !---- CRV( QJ UN(U V) VN(U V) QI(U V) )
        CRV_QJ(1, 1) = -UN * VN * QI
        CRV_QJ(1, 2) = -VN * VN * QI
        CRV_QJ(2, 1) = UN * UN * QI
        CRV_QJ(2, 2) = UN * VN * QI
        !
        CRV_UN = (2.0 * UN * QJ(2, 1) - VN * QJ(1, 1) + VN * QJ(2, 2)) * QI
        CRV_VN = (-UN * QJ(1, 1) + UN * QJ(2, 2) - 2.0 * VN * QJ(1, 2)) * QI
        CRV_QI = UN * UN * QJ(2, 1) - UN * VN * QJ(1, 1) + UN * VN * QJ(2, 2)            &
                & - VN * VN * QJ(1, 2)
        !
        !---- CRV( QJ U V )
        CRV_U = CRV_UN * UN_U + CRV_VN * VN_U + CRV_QI * QI_U
        CRV_V = CRV_UN * UN_V + CRV_VN * VN_V + CRV_QI * QI_V
        !
    END SUBROUTINE QCURV
    !*==CV2SET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCURV
    
    
    
    SUBROUTINE CV2SET(ICM, ICO, QC, XC, YC, ANC, GAMC, SIGC, CV, CV_QC, CV_RC, &
            & CV_ANC, CV_GAMC, CV_SIGC)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: CV
        INTEGER :: ICM, ICO
        REAL, DIMENSION(2, *) :: ANC, QC
        REAL, DIMENSION(2, -1:0) :: CV_ANC, CV_QC, CV_RC
        REAL, DIMENSION(-1:0) :: CV_GAMC, CV_SIGC
        REAL, DIMENSION(*) :: GAMC, SIGC, XC, YC
        !
        ! Local variables
        !
        REAL :: CV_DQ, CV_DR, GAMCM, GAMCO, QSGN, SIGCM, SIGCO
        REAL, DIMENSION(2) :: CV_Q, CV_QSM, CV_QSO, DQ, DR, QSA, &
                & QSM, QSO
        REAL, DIMENSION(2, 2) :: CV_QJ, QJ
        INTEGER :: K
        REAL, DIMENSION(2, 2, 2) :: QJ_DQ, QJ_DR
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !
        GAMCM = GAMC(ICM)
        SIGCM = SIGC(ICM)
        !
        GAMCO = GAMC(ICO)
        SIGCO = SIGC(ICO)
        !
        !---- 1/2 of sheet strength adds,subtracts to left,right speed
        QSGN = +0.5
        !
        QSM(1) = QC(1, ICM) + QSGN * (ANC(1, ICM) * SIGCM + ANC(2, ICM) * GAMCM)
        QSM(2) = QC(2, ICM) + QSGN * (ANC(2, ICM) * SIGCM - ANC(1, ICM) * GAMCM)
        !
        QSO(1) = QC(1, ICO) + QSGN * (ANC(1, ICO) * SIGCO + ANC(2, ICO) * GAMCO)
        QSO(2) = QC(2, ICO) + QSGN * (ANC(2, ICO) * SIGCO - ANC(1, ICO) * GAMCO)
        !
        DQ(1) = QSO(1) - QSM(1)
        DQ(2) = QSO(2) - QSM(2)
        !
        QSA(1) = (QSO(1) + QSM(1)) * 0.5
        QSA(2) = (QSO(2) + QSM(2)) * 0.5
        !
        DR(1) = XC(ICO) - XC(ICM)
        DR(2) = YC(ICO) - YC(ICM)
        !
        !
        CALL QJAC(DR(1), DR(2), DQ(1), DQ(2), QJ, QJ_DR(1, 1, 1), QJ_DR(1, 1, 2), &
                & QJ_DQ(1, 1, 1), QJ_DQ(1, 1, 2))
        !
        CALL QCURV(QSA(1), QSA(2), QJ, CV, CV_Q(1), CV_Q(2), CV_QJ)
        !
        DO K = 1, 2
            CV_DQ = CV_QJ(1, 1) * QJ_DQ(1, 1, K) + CV_QJ(1, 2) * QJ_DQ(1, 2, K)      &
                    & + CV_QJ(2, 1) * QJ_DQ(2, 1, K) + CV_QJ(2, 2) * QJ_DQ(2, 2, K)
            CV_DR = CV_QJ(1, 1) * QJ_DR(1, 1, K) + CV_QJ(1, 2) * QJ_DR(1, 2, K)      &
                    & + CV_QJ(2, 1) * QJ_DR(2, 1, K) + CV_QJ(2, 2) * QJ_DR(2, 2, K)
            !
            CV_QSM(K) = -CV_DQ + 0.5 * CV_Q(K)
            CV_QSO(K) = CV_DQ + 0.5 * CV_Q(K)
            !
            CV_QC(K, -1) = CV_QSM(K)
            CV_QC(K, 0) = CV_QSO(K)
            !
            CV_RC(K, -1) = -CV_DR
            CV_RC(K, 0) = CV_DR
            !
            CV_ANC(K, -1) = CV_QSM(K) * QSGN * SIGCM
            CV_ANC(K, 0) = CV_QSO(K) * QSGN * SIGCO
        ENDDO
        !
        CV_ANC(1, -1) = CV_ANC(1, -1) - CV_QSM(2) * QSGN * GAMCM
        CV_ANC(1, 0) = CV_ANC(1, 0) - CV_QSO(2) * QSGN * GAMCO
        CV_ANC(2, -1) = CV_ANC(2, -1) + CV_QSM(1) * QSGN * GAMCM
        CV_ANC(2, 0) = CV_ANC(2, 0) + CV_QSO(1) * QSGN * GAMCO
        !
        CV_SIGC(-1) = QSGN * (CV_QSM(1) * ANC(1, ICM) + CV_QSM(2) * ANC(2, ICM))
        CV_SIGC(0) = QSGN * (CV_QSO(1) * ANC(1, ICO) + CV_QSO(2) * ANC(2, ICO))
        !
        CV_GAMC(-1) = QSGN * (CV_QSM(1) * ANC(2, ICM) - CV_QSM(2) * ANC(1, ICM))
        CV_GAMC(0) = QSGN * (CV_QSO(1) * ANC(2, ICO) - CV_QSO(2) * ANC(1, ICO))
        !
    END SUBROUTINE CV2SET
    !*==CV3SET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CV2SET
    
    
    SUBROUTINE CV3SET(ICM, ICO, ICP, QC, XC, YC, ANC, ANSGN, GAMC, SIGC, CV, &
            & CV_QC, CV_RC, CV_ANC, CV_GAMC, CV_SIGC)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: ANSGN, CV
        INTEGER :: ICM, ICO, ICP
        REAL, DIMENSION(2, *) :: ANC, QC
        REAL, DIMENSION(2, -1:1) :: CV_ANC, CV_QC, CV_RC
        REAL, DIMENSION(-1:1) :: CV_GAMC, CV_SIGC
        REAL, DIMENSION(*) :: GAMC, SIGC, XC, YC
        !
        ! Local variables
        !
        REAL :: CSGN, CV_DQ, CV_DR, GAMCM, GAMCO, GAMCP, QSGN, &
                & QXN, SIGCM, SIGCO, SIGCP
        REAL, DIMENSION(2) :: CV_Q, CV_QSM, CV_QSO, CV_QSP, DQ, &
                & DR, QSM, QSO, QSP
        REAL, DIMENSION(2, 2) :: CV_QJ, QJ
        INTEGER :: K
        REAL, DIMENSION(2, 2, 2) :: QJ_DQ, QJ_DR
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !
        GAMCM = GAMC(ICM)
        SIGCM = SIGC(ICM)
        !
        GAMCO = GAMC(ICO)
        SIGCO = SIGC(ICO)
        !
        GAMCP = GAMC(ICP)
        SIGCP = SIGC(ICP)
        !
        !---- 1/2 of sheet strength adds,subtracts to left,right speed
        QSGN = 0.5 * ANSGN
        !
        QSM(1) = QC(1, ICM) + QSGN * (ANC(1, ICM) * SIGCM + ANC(2, ICM) * GAMCM)
        QSM(2) = QC(2, ICM) + QSGN * (ANC(2, ICM) * SIGCM - ANC(1, ICM) * GAMCM)
        !
        QSO(1) = QC(1, ICO) + QSGN * (ANC(1, ICO) * SIGCO + ANC(2, ICO) * GAMCO)
        QSO(2) = QC(2, ICO) + QSGN * (ANC(2, ICO) * SIGCO - ANC(1, ICO) * GAMCO)
        !
        QSP(1) = QC(1, ICP) + QSGN * (ANC(1, ICP) * SIGCP + ANC(2, ICP) * GAMCP)
        QSP(2) = QC(2, ICP) + QSGN * (ANC(2, ICP) * SIGCP - ANC(1, ICP) * GAMCP)
        !
        DQ(1) = QSP(1) - QSM(1)
        DQ(2) = QSP(2) - QSM(2)
        !
        DR(1) = XC(ICP) - XC(ICM)
        DR(2) = YC(ICP) - YC(ICM)
        !
        QXN = ANSGN * (QSO(1) * ANC(2, ICO) - QSO(2) * ANC(1, ICO))
        CSGN = SIGN(1.0, QXN)
        !
        CALL QJAC(DR(1), DR(2), DQ(1), DQ(2), QJ, QJ_DR(1, 1, 1), QJ_DR(1, 1, 2), &
                & QJ_DQ(1, 1, 1), QJ_DQ(1, 1, 2))
        !
        CALL QCURV(QSO(1), QSO(2), QJ, CV, CV_Q(1), CV_Q(2), CV_QJ)
        !
        CV = CSGN * CV
        DO K = 1, 2
            CV_Q(K) = CSGN * CV_Q(K)
            CV_QJ(K, 1) = CSGN * CV_QJ(K, 1)
            CV_QJ(K, 2) = CSGN * CV_QJ(K, 2)
        ENDDO
        !
        DO K = 1, 2
            CV_DQ = CV_QJ(1, 1) * QJ_DQ(1, 1, K) + CV_QJ(1, 2) * QJ_DQ(1, 2, K)      &
                    & + CV_QJ(2, 1) * QJ_DQ(2, 1, K) + CV_QJ(2, 2) * QJ_DQ(2, 2, K)
            CV_DR = CV_QJ(1, 1) * QJ_DR(1, 1, K) + CV_QJ(1, 2) * QJ_DR(1, 2, K)      &
                    & + CV_QJ(2, 1) * QJ_DR(2, 1, K) + CV_QJ(2, 2) * QJ_DR(2, 2, K)
            !
            CV_QSM(K) = -CV_DQ
            CV_QSO(K) = CV_Q(K)
            CV_QSP(K) = CV_DQ
            !
            CV_QC(K, -1) = CV_QSM(K)
            CV_QC(K, 0) = CV_QSO(K)
            CV_QC(K, +1) = CV_QSP(K)
            !
            CV_RC(K, -1) = -CV_DR
            CV_RC(K, 0) = 0.
            CV_RC(K, +1) = CV_DR
            !
            CV_ANC(K, -1) = CV_QSM(K) * QSGN * SIGCM
            CV_ANC(K, 0) = CV_QSO(K) * QSGN * SIGCO
            CV_ANC(K, +1) = CV_QSP(K) * QSGN * SIGCP
        ENDDO
        !
        CV_ANC(1, -1) = CV_ANC(1, -1) - CV_QSM(2) * QSGN * GAMCM
        CV_ANC(1, 0) = CV_ANC(1, 0) - CV_QSO(2) * QSGN * GAMCO
        CV_ANC(1, +1) = CV_ANC(1, +1) - CV_QSP(2) * QSGN * GAMCP
        CV_ANC(2, -1) = CV_ANC(2, -1) + CV_QSM(1) * QSGN * GAMCM
        CV_ANC(2, 0) = CV_ANC(2, 0) + CV_QSO(1) * QSGN * GAMCO
        CV_ANC(2, +1) = CV_ANC(2, +1) + CV_QSP(1) * QSGN * GAMCP
        !
        CV_SIGC(-1) = QSGN * (CV_QSM(1) * ANC(1, ICM) + CV_QSM(2) * ANC(2, ICM))
        CV_SIGC(0) = QSGN * (CV_QSO(1) * ANC(1, ICO) + CV_QSO(2) * ANC(2, ICO))
        CV_SIGC(+1) = QSGN * (CV_QSP(1) * ANC(1, ICP) + CV_QSP(2) * ANC(2, ICP))
        !
        CV_GAMC(-1) = QSGN * (CV_QSM(1) * ANC(2, ICM) - CV_QSM(2) * ANC(1, ICM))
        CV_GAMC(0) = QSGN * (CV_QSO(1) * ANC(2, ICO) - CV_QSO(2) * ANC(1, ICO))
        CV_GAMC(+1) = QSGN * (CV_QSP(1) * ANC(2, ICP) - CV_QSP(2) * ANC(1, ICP))
        !
    END SUBROUTINE CV3SET
    !*==GETUV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CV3SET
    
    
    SUBROUTINE GETUV(XX, YY, US, VS)
        use m_qaic, only: qfcalc
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: US, VS, XX, YY
        !
        ! Local variables
        !
        INTEGER :: IFTYPE, IP1, IP2, IPFO, IPFP, NF
        REAL, DIMENSION(2) :: QF
        REAL, DIMENSION(2, IPX) :: QF_GAM, QF_GTH, QF_SIG
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Get EIF velocity US,VS at point XX,YY
        !-----------------------------------------------
        !
        !---- local arrays for calling QFCALC
        !
        !------ we will evaluate only this point
        IP1 = 1
        IP2 = 1
        NF = 1
        !------ assume the point is not on a panel
        IPFO = 0
        IPFP = 0
        IFTYPE = 0
        !------ evaluate velocity components
        CALL QFCALC(IP1, IP2, XX, YY, IPFO, IPFP, IFTYPE, QF, QF_GAM, QF_SIG, &
                & QF_GTH)
        !------ add freestream
        US = QF(1) + QINF
        VS = QF(2)
        !
    END SUBROUTINE GETUV
end module m_vels
