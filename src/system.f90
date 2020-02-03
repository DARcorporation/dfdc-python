module m_system
    implicit none
contains
    !*==SYSP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! INTERS
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
    ! This version assumes delta(B*GAM)*WT/WM defined at each rotor center
    
    SUBROUTINE SYSP
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: IC, IEL, IP, IR, JAICG, JAICS, JAICV, JAICW, &
                & JAICZ, JP, JSYS, KSYS, NDOF
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------
        !     Sets up pointers associating....
        !        matrix rows with control points,
        !        matrix columns with unknowns
        !-------------------------------------------
        !
        IF (LDBG) WRITE (*, *) 'Generating flow-tangency system pointers'
        !
        !---- set what type of equation gets assigned to each control point,
        !-    and which unknown vortex gets associated with each control point
        !
        !---- -999 indicates un-enforced control point
        !----    0 indicates undefined variable
        DO IC = 1, NCTOT
            KSYSVNC(IC) = -999
            KSYSGMG(IC) = -999
            IELDOFC(IC) = 0
        ENDDO
        DO IP = 1, NPTOT
            KSYSGAM(IP) = -999
            JSYSGAM(IP) = 0
            JSYSDXY(IP) = 0
            JAICGAM(IP) = 0
            JAICSIG(IP) = 0
            JAICXYP(IP) = 0
            JAICGTH(IP) = 0
        ENDDO
        DO IEL = 1, NEL
            KSYSKUT(IEL) = -999
            KSYSQNC(IEL) = -999
            KSYSGSS(IEL) = -999
            KSYSQTT(IEL) = -999
            JSYSQNC(IEL) = 0
            JAICVWK(IEL) = 0
        ENDDO
        !
        !
        !---- initialize system-row and system-column counters
        KSYS = 0
        JSYS = 0
        !
        !---- initialize gamma,sigma,xy column counters
        JAICG = 0
        JAICS = 0
        JAICZ = 0
        !---- initialize wake gamma (point) column counters
        JAICW = 0
        !---- initialize vortex wake row column counters
        JAICV = 0
        !
        !---- go over all elements
        DO IEL = 1, NEL
            !
            !cC------ LVeZR = f if first,last vortex node is not an unknown
            !c        LV1ZR(IEL) = YP(IPFRST(IEL)) .EQ. 0.0
            !c        LV2ZR(IEL) = YP(IPLAST(IEL)) .EQ. 0.0
            !
            IF (NETYPE(IEL)==0) THEN
                !-------- normal solid-wall airfoil...
                !
                !--------------------------------------------------
                !-------- set up variable indices
                !
                DO JP = IPFRST(IEL), IPLAST(IEL)
                    !---------- wall gamma will be unknown
                    JSYS = JSYS + 1
                    JSYSGAM(JP) = JSYS
                    JPSYS(JSYS) = JP
                    !---------- wall sigma will be imposed externally (e.g. BL blowing model)
                    JAICS = JAICS + 1
                    JAICSIG(JP) = JAICS
                ENDDO
                !
                IF (LBODY(IEL)) THEN
                    !--------- declare normal-velocity DOF, and set element which it influences
                    JSYS = JSYS + 1
                    JSYSQNC(IEL) = JSYS
                    JPSYS(JSYS) = 0
                    NDOF = IEL
                ELSE
                    NDOF = 0
                ENDIF
                !
                IF (LXYMOV) THEN
                    DO JP = IPFRST(IEL), IPLAST(IEL)
                        !----------- wall x,y movement will be imposed externally (e.g. optimization)
                        JAICZ = JAICZ + 1
                        JAICXYP(JP) = JAICZ
                    ENDDO
                ENDIF
                !
                !--------------------------------------------------
                !-------- set up equation indices
                !
                IF (LV1ZR(IEL)) THEN
                    !---------- enforce zero-gamma condition at first point
                    KSYS = KSYS + 1
                    KSYSGAM(IPFRST(IEL)) = KSYS
                ENDIF
                !
                !-------- enforce flow tangency at all control points, note Qn DOF if any
                DO IC = ICFRST(IEL), ICLAST(IEL)
                    KSYS = KSYS + 1
                    !
                    IF (ICTYPE(IC)>=0) THEN
                        !----------- usual control point
                        KSYSVNC(IC) = KSYS
                        IELDOFC(IC) = NDOF
                    ELSE
                        !----------- zero-length panel... set gamma-continuity equation
                        KSYSGMG(IC) = KSYS
                    ENDIF
                ENDDO
                !
                IF (LV2ZR(IEL)) THEN
                    !---------- enforce zero-gamma condition at last point
                    KSYS = KSYS + 1
                    KSYSGAM(IPLAST(IEL)) = KSYS
                ENDIF
                !
                IF (.NOT.(LV1ZR(IEL) .OR. LV2ZR(IEL))) THEN
                    !--------- enforce Kutta condition
                    KSYS = KSYS + 1
                    KSYSKUT(IEL) = KSYS
                ENDIF
                !
                IF (LBODY(IEL) .AND. .NOT.(LV1ZR(IEL) .AND. LV2ZR(IEL)))  &
                        & THEN
                    !--------- will have normal-velocity DOF... enforce some constraint on it
                    KSYS = KSYS + 1
                    IF (LQNZR(IEL)) THEN
                        !---------- direct constraint
                        KSYSQNC(IEL) = KSYS
                    ELSE
                        !ccC---------- indirect constraint: set zero loading curvature at TE
                        !cc            KSYSGSS(IEL) = KSYS
                        !
                        !---------- indirect constraint: set zero velocity at interior point near TE
                        KSYSQTT(IEL) = KSYS
                    ENDIF
                ENDIF
                !
                !-------- check for vortex wake gamma on body element
                DO JP = IPFRST(IEL), IPLAST(IEL)
                    IF (IP2IR(JP)>0) THEN
                        !---------- GTH will be imposed or solved for at this point
                        JAICW = JAICW + 1
                        JAICGTH(JP) = JAICW
                    ENDIF
                ENDDO
                !
            ELSEIF (NETYPE(IEL)==1) THEN
                !-------- wake element...
                !---------- wake gamma will be imposed externally (e.g. BL curvature relation)
                DO JP = IPFRST(IEL), IPLAST(IEL)
                    JAICG = JAICG + 1
                    JAICGAM(JP) = JAICG
                    !---------- wake sigma will be imposed externally (e.g. BL blowing model)
                    !c            JAICS = JAICS + 1
                    !c            JAICSIG(JP) = JAICS
                ENDDO
                !
                IF (LXYMOV) THEN
                    DO JP = IPFRST(IEL), IPLAST(IEL)
                        !----------- wake x,y movement will be imposed externally (e.g. optimization)
                        JAICZ = JAICZ + 1
                        JAICXYP(JP) = JAICZ
                    ENDDO
                ENDIF
                !
            ELSEIF (NETYPE(IEL)>=2 .AND. NETYPE(IEL)<=4) THEN
                !-------- axis line, ring or point singularities...(NETYPE=2,3,4)
                !-------- sigma, doublet will be imposed externally
                DO JP = IPFRST(IEL), IPLAST(IEL)
                    JAICG = JAICG + 1
                    JAICGAM(JP) = JAICG
                    JAICS = JAICS + 1
                    JAICSIG(JP) = JAICS
                ENDDO
                !
            ELSEIF (NETYPE(IEL)==5) THEN
                !-------- drag area source element NETYPE=5 (source RHS only)
                !-------- sigma will be imposed externally (e.g. drag source)
                DO JP = IPFRST(IEL), IPLAST(IEL)
                    JAICS = JAICS + 1
                    JAICSIG(JP) = JAICS
                ENDDO
                !
            ELSEIF (NETYPE(IEL)==6) THEN
                !-------- rotor source element NETYPE=6 (source RHS only)
                !-------- sigma will be imposed externally (e.g. drag source)
                DO JP = IPFRST(IEL), IPLAST(IEL)
                    JAICS = JAICS + 1
                    JAICSIG(JP) = JAICS
                ENDDO
                !
            ELSEIF (NETYPE(IEL)==7) THEN
                !-------- vortex wake element (RHS)
                !-------- wake gamma will be imposed externally
                IF (JAICV<NRP) THEN
                    JAICV = JAICV + 1
                    IR = IEL2IR(IEL)
                    JAICVWK(IR) = JAICV
                    !c            write(19,*) 'iel,ir,jaicv ',iel,ir,jaicv
                ENDIF
                !-------- set up points with vortex wake gamma
                DO JP = IPFRST(IEL), IPLAST(IEL)
                    !---------- GTH will be imposed or solved for at this point
                    JAICW = JAICW + 1
                    JAICGTH(JP) = JAICW
                ENDDO
                !
            ELSE
                WRITE (*, *) 'Unknown NETYPE in SYSP ', NETYPE(IEL)
                !
            ENDIF
        ENDDO
        !
        !---- set system size
        NSYS = KSYS
        !
        IF (KSYS/=JSYS) WRITE (*, *) 'Error in SYSP, KSYS JSYS =', &
                & KSYS, JSYS
        !
        NAICGAM = JAICG
        NAICSIG = JAICS
        NAICXYP = JAICZ
        NAICVWK = JAICV
        NAICGTH = JAICW
        !
        IF (LDBG) THEN
            WRITE (*, *) ' '
            WRITE (*, *) 'SYSP system setup and pointers'
            WRITE (*, *) '  NP      = ', NPTOT
            WRITE (*, *) '  NC      = ', NCTOT
            WRITE (*, *) '  NSYS    = ', NSYS
            WRITE (*, *) ' '
            WRITE (*, *) '  NAICGAM = ', NAICGAM
            WRITE (*, *) '  NAICSIG = ', NAICSIG
            WRITE (*, *) '  NAICXYP = ', NAICXYP
            WRITE (*, *) '  NAICVWK = ', NAICVWK
            WRITE (*, *) '  NAICGTH = ', NAICGTH
            WRITE (*, *) ' '
            WRITE (*, *) ' '
            !
            WRITE (LUNDBG, *) 'SYSP system setup and pointers'
            WRITE (LUNDBG, *) '  NP      = ', NPTOT
            WRITE (LUNDBG, *) '  NC      = ', NCTOT
            WRITE (LUNDBG, *) '  NSYS    = ', NSYS
            WRITE (LUNDBG, *) ' '
            WRITE (LUNDBG, *) '  NAICGAM = ', NAICGAM
            WRITE (LUNDBG, *) '  NAICSIG = ', NAICSIG
            WRITE (LUNDBG, *) '  NAICXYP = ', NAICXYP
            WRITE (LUNDBG, *) '  NAICVWK = ', NAICVWK
            WRITE (LUNDBG, *) '  NAICGTH = ', NAICGTH
            WRITE (LUNDBG, *) ' '
        ENDIF
        !
        LSYSP = .TRUE.
        !
    END SUBROUTINE SYSP
    !*==CLRSYS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SYSP
    
    
    SUBROUTINE CLRSYS
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: J, K
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- clear all system arrays
        DO K = 1, NSYS
            RES(K) = 0.
            DO J = 0, NSYS
                SYS(K, J) = 0.
            ENDDO
            AICQFF(K, 1) = 0.
            DO J = 0, NAICGAM
                AICGAM(K, J) = 0.
            ENDDO
            DO J = 0, NAICSIG
                AICSIG(K, J) = 0.
            ENDDO
            DO J = 0, NAICXYP
                AICXYP(K, 1, J) = 0.
                AICXYP(K, 2, J) = 0.
            ENDDO
            DO J = 0, NAICVWK
                AICVWK(K, J) = 0.
            ENDDO
        ENDDO
        !
    END SUBROUTINE CLRSYS
    !*==GSYS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CLRSYS
    
    
    SUBROUTINE GSYS(LSYS, LQFF, LVWK, IP1, IP2, IZ1, IZ2)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IP1, IP2, IZ1, IZ2
        LOGICAL :: LQFF, LSYS, LVWK
        !
        ! Local variables
        !
        REAL :: DN1SW, DN2SW, DS1M, DS1O, DS1SW, DS2M, DS2O, &
                & DS2SW, F1L, F1M, F1O, F2L, F2M, F2O, QTTF, &
                & RES_DX, RES_DY
        INTEGER :: IC, IE, IEL, IELW, IPO, IPP, IR, J, J1, &
                & J1L, J1M, J1O, J2, J2L, J2M, J2O, JP, JP1, &
                & JP1L, JP1M, JP1O, JP2, JP2L, JP2M, JP2O, JPO, &
                & JPP, JPW, JW, K
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------
        !     Generates influence matrices of unknown gamma
        !     and freestream on flow tangency and other residuals.
        !
        !     Input:  LSYS     if T, set up surface-gamma Jacobian
        !             LQFF     if T, set up freestream Jacobian
        !             IP1,IP2  set up Jacobian for gamma,sigma of nodes IP1..IP2
        !             IZ1,IZ2  set up Jacobian for x,y of nodes IZ1..IZ2
        !
        !     Input:  NSYS         system size
        !             GAM(j)       panel-node vortex sheet strengths
        !             QNDOF(e)     element normalwash DOF variables
        !             ANC(.i)      control-point normal vectors
        !             QC(.i)       control-point velocities
        !             QC_GAM(.ij)  AIC matrix  dQC(.i)/dGAM(j) on surfaces
        !             QC_SIG(.ij)  AIC matrix  dQC(.i)/dSIG(j) on surfaces
        !             QC_XYP(.i.j) AIC matrix  dQC(.i)/dXP(j),dYP(j)
        !             QC_GTH(.ij)  AIC matrix  dQC(.i)/dGTH(j) on surfaces
        !             KSYS...      active-equation pointers
        !             JSYS...      active-variable pointers
        !
        !     Output: linear system components
        !
        !        LHS system
        !             RES(k)       Qn(k)  for current QC(j), QNDOF(e)
        !                      or  GAM-continuity residual
        !                      or  Kutta condition residual
        !                      or  TE loading regularity residual
        !             SYS(kj)      dRES(k)/dGAM(j) , dRES(k)/dQNDOF(e)   if LSYS=t
        !
        !        RHS columns
        !             AICQFF(k)    dRES(k)/dQinf                         if LQFF=t
        !             AICVWK(e)    dRES(k)/dWBGAM on rotor               if LVWK=t
        !
        !             AICGAM(kj)   dRES(k)/dGAM(j)    j = IP1..IP2
        !             AICSIG(kj)   dRES(k)/dSIG(j)    j = IP1..IP2
        !
        !             AICXYP(k.j)  dRES(k)/dXP(j)     j = IZ1..IZ2
        !                          dRES(k)/dYP(j)     j = IZ1..IZ2
        !
        !-----------------------------------------------------------
        !
        IF (NSYS==0) RETURN
        !
        IF (LDBG) THEN
            WRITE (*, *) 'Setting up direct-problem Newton system...'
            WRITE (LUNDBG, *) 'GSYS direct problem system setup'
            WRITE (LUNDBG, *) '   LSYS,LQFF,LVWK ', LSYS, LQFF, LVWK
            WRITE (LUNDBG, *) '   IP1,IP2        ', IP1, IP2
            WRITE (LUNDBG, *) '   IZ1,IZ2        ', IZ1, IZ2
        ENDIF
        !
        !----------------------------------------------------------
        !---- clear system arrays
        IF (LSYS) THEN
            DO K = 1, NSYS
                RES(K) = 0.
                DO J = 1, NSYS
                    SYS(K, J) = 0.
                ENDDO
            ENDDO
        ENDIF
        !
        !---- clear RHS arrays
        IF (LQFF) THEN
            DO K = 1, NSYS
                AICQFF(K, 1) = 0.
            ENDDO
        ENDIF
        !
        IF (LVWK) THEN
            DO IR = 1, NRP
                J = JAICVWK(IR)
                DO K = 1, NSYS
                    AICVWK(K, J) = 0.
                ENDDO
            ENDDO
        ENDIF
        !
        DO K = 1, NSYS
            DO JP = IP1, IP2
                J = JAICGAM(JP)
                AICGAM(K, J) = 0.
                !
                J = JAICSIG(JP)
                AICSIG(K, J) = 0.
                !
                J = JAICGTH(JP)
                AICGTH(K, J) = 0.
            ENDDO
            !
            DO JP = IZ1, IZ2
                J = JAICXYP(JP)
                AICXYP(K, 1, J) = 0.
                AICXYP(K, 2, J) = 0.
            ENDDO
        ENDDO
        !
        !----------------------------------------------------------
        !----------------------------------------------------------
        !---- set flow tangency at all control points
        DO IC = 1, NCTOT
            !
            !        write(77,*) 'IC ',IC
            !        do jp = 1, np
            !          write(77,*) 'JP, QC1-2 ',JP,QC_GTH(1,JP,IC),QC_GTH(2,JP,IC)
            !        end do
            !------ K = system row index for control point IC
            K = KSYSVNC(IC)
            !
            !------ skip un-enforced control points
            IF (K>=0) THEN
                !
                !----------------------------------------------------------
                IF (LSYS) THEN
                    !------- flow-tangency residual
                    RES(K) = ANC(1, IC) * QC(1, IC) + ANC(2, IC) * QC(2, IC)
                    !
                    !
                    DO JP = 1, NPTOT
                        !--------- dRES/dgamma_surface
                        !-         J = system column index of GAM(JP)  (0 if GAM is on r.h.s.)
                        J = JSYSGAM(JP)
                        SYS(K, J) = ANC(1, IC) * QC_GAM(1, JP, IC) + ANC(2, IC)      &
                                & * QC_GAM(2, JP, IC)
                    ENDDO
                ENDIF
                !
                !----------------------------------------------------------
                !------- dRES/dQinf
                IF (LQFF) AICQFF(K, 1) = ANC(1, IC)
                !
                !----------------------------------------------------------
                IF (LVWK) THEN
                    !c          write(19,*) 'ic ',ic
                    !------- dRES/dBCirc rotor vortex wake influence
                    DO JP = 1, NPTOT
                        IR = IP2IR(JP)
                        IF (IR/=0) THEN
                            !---- we are on a vortex wake element
                            !
                            !---- accumulate wake influences by wake rows using wake index
                            J = JAICVWK(IR)
                            AICVWK(K, J) = AICVWK(K, J)                          &
                                    & + (ANC(1, IC) * QC_GTH(1, JP, IC)         &
                                            & + ANC(2, IC) * QC_GTH(2, JP, IC)) / YP(JP)
                            !c              write(19,*) 'jp,ir,j ',jp,ir,j
                            !
                            !---- set wake point influence individually
                            J = JAICGTH(JP)
                            AICGTH(K, J) = ANC(1, IC) * QC_GTH(1, JP, IC) + ANC(2, IC)&
                                    & * QC_GTH(2, JP, IC)
                        ENDIF
                    ENDDO
                ENDIF
                !
                !----------------------------------------------------------
                DO JP = IP1, IP2
                    !-------- dRES/dgamma, dRES/dsigma
                    !-        J = r.h.s. column index  (0 if GAM,SIG is a l.h.s. unknown)
                    J = JAICGAM(JP)
                    AICGAM(K, J) = ANC(1, IC) * QC_GAM(1, JP, IC) + ANC(2, IC)      &
                            & * QC_GAM(2, JP, IC)
                    !
                    J = JAICSIG(JP)
                    AICSIG(K, J) = ANC(1, IC) * QC_SIG(1, JP, IC) + ANC(2, IC)      &
                            & * QC_SIG(2, JP, IC)
                ENDDO
                !
                !----------------------------------------------------------
                DO JP = IZ1, IZ2
                    !-------- dRES/dxy
                    !-        J = r.h.s. column index
                    J = JAICXYP(JP)
                    AICXYP(K, 1, J) = ANC(1, IC) * QC_XP(1, JP, IC) + ANC(2, IC)     &
                            & * QC_XP(2, JP, IC)
                    AICXYP(K, 2, J) = ANC(1, IC) * QC_YP(1, JP, IC) + ANC(2, IC)     &
                            & * QC_YP(2, JP, IC)
                ENDDO
                !
                !----------------------------------------------------------
                !------ local source contribution to flow-tangency residual
                IPO = IPCO(IC)
                IPP = IPCP(IC)
                !
                IF (LSYS) RES(K) = RES(K) - 0.25 * (SIG(IPO) + SIG(IPP))
                !
                IF (IPO>=IP1 .AND. IPO<=IP2) THEN
                    J = JAICSIG(IPO)
                    AICSIG(K, J) = AICSIG(K, J) - 0.25
                ENDIF
                IF (IPP>=IP1 .AND. IPP<=IP2) THEN
                    J = JAICSIG(IPP)
                    AICSIG(K, J) = AICSIG(K, J) - 0.25
                ENDIF
                !
                !------ derivatives w.r.t. dn movement of the two nodes of panel IC
                RES_DX = ANC_DXY(1, 1, IC) * QC(1, IC) + ANC_DXY(2, 1, IC) * QC(2, IC)
                RES_DY = ANC_DXY(1, 2, IC) * QC(1, IC) + ANC_DXY(2, 2, IC) * QC(2, IC)
                !
                IF (IPO>=IZ1 .AND. IPO<=IZ2) THEN
                    J = JAICXYP(IPO)
                    AICXYP(K, 1, J) = AICXYP(K, 1, J) - RES_DX
                    AICXYP(K, 2, J) = AICXYP(K, 2, J) - RES_DY
                ENDIF
                IF (IPP>=IZ1 .AND. IPP<=IZ2) THEN
                    J = JAICXYP(IPP)
                    AICXYP(K, 1, J) = AICXYP(K, 1, J) + RES_DX
                    AICXYP(K, 2, J) = AICXYP(K, 2, J) + RES_DY
                ENDIF
                !
                !----------------------------------------------------------
                IF (LSYS .AND. IELDOFC(IC)>0) THEN
                    !------- this control point has a normal-velocity QNDOF contributing to Q.n
                    IE = IELDOFC(IC)
                    RES(K) = RES(K) - QNDOF(IE)
                    !
                    !------- J = system column index of QNDOF(IE)
                    J = JSYSQNC(IE)
                    IF (J<=0) THEN
                        WRITE (*, *) '? Unassigned normal-velocity DOF', IE
                        WRITE (LUNDBG, *) '? Unassigned normal-velocity DOF', &
                                & IE
                    ENDIF
                    !
                    SYS(K, J) = -1.0
                ENDIF
                IF (LDBG) WRITE (LUNDBG, *) K, 'V.n=Qn ', RES(K)
                !@@@
            ENDIF
        ENDDO
        !
        !---- set gamma-continuity at all zero-length panels
        DO IC = 1, NCTOT
            !------ K = system row index of gamma-continuity equation (if any)
            K = KSYSGMG(IC)
            !
            IF (LSYS .AND. K>0) THEN
                JPO = IPCO(IC)
                JPP = IPCP(IC)
                RES(K) = GAM(JPP) - GAM(JPO)
                !
                SYS(K, JPO) = -1.0
                SYS(K, JPP) = 1.0
                IF (LDBG) WRITE (LUNDBG, *) K, 'dGam=0 ', RES(K)
                !@@@
            ENDIF
        ENDDO
        !
        !
        DO IEL = 1, NEL
            !------ zero net curvature of gamma at TE (triangular loading towards TE)
            K = KSYSGSS(IEL)
            IF (LSYS .AND. K>0) THEN
    
                JP1O = IPFRST(IEL) + 2
                JP1M = IPFRST(IEL) + 1
                JP1L = IPFRST(IEL)
                JP2O = IPLAST(IEL)
                JP2M = IPLAST(IEL) - 1
                JP2L = IPLAST(IEL) - 2
                !
                DS1O = SP(JP1O) - SP(JP1M)
                DS1M = SP(JP1M) - SP(JP1L)
                DS2O = SP(JP2O) - SP(JP2M)
                DS2M = SP(JP2M) - SP(JP2L)
                !
                !-------- 1st derivative
                !          F1O = 1.0/DS1O + 1.0/(DS1O+DS1M)
                !          F1L = DS1O/((DS1O+DS1M)*DS1M)
                !          F1M = -(F1O+F1L)
                !          F2O = 1.0/DS2O + 1.0/(DS2O+DS2M)
                !          F2L = DS2O/((DS2O+DS2M)*DS2M)
                !          F2M = -(F2O+F2L)
                !
                !-------- 2nd derivative
                F1O = 1.0 / DS1O
                F1M = -(1.0 / DS1O + 1.0 / DS1M)
                F1L = 1.0 / DS1M
                F2O = 1.0 / DS2O
                F2M = -(1.0 / DS2O + 1.0 / DS2M)
                F2L = 1.0 / DS2M
                !
                J1O = JSYSGAM(JP1O)
                J1M = JSYSGAM(JP1M)
                J1L = JSYSGAM(JP1L)
                J2O = JSYSGAM(JP2O)
                J2M = JSYSGAM(JP2M)
                J2L = JSYSGAM(JP2L)
                !
                IF (LV1ZR(IEL)) THEN
                    RES(K) = -F2O * GAM(JP2O) - F2M * GAM(JP2M) - F2L * GAM(JP2L)
                    SYS(K, J2O) = -F2O
                    SYS(K, J2M) = -F2M
                    SYS(K, J2L) = -F2L
                    !
                ELSEIF (LV2ZR(IEL)) THEN
                    RES(K) = F1O * GAM(JP1O) + F1M * GAM(JP1M) + F1L * GAM(JP1L)
                    SYS(K, J1O) = F1O
                    SYS(K, J1M) = F1M
                    SYS(K, J1L) = F1L
                    !
                ELSE
                    RES(K) = F1O * GAM(JP1O) + F1M * GAM(JP1M) + F1L * GAM(JP1L)   &
                            & - F2O * GAM(JP2O) - F2M * GAM(JP2M) - F2L * GAM(JP2L)
                    SYS(K, J1O) = F1O
                    SYS(K, J1M) = F1M
                    SYS(K, J1L) = F1L
                    SYS(K, J2O) = -F2O
                    SYS(K, J2M) = -F2M
                    SYS(K, J2L) = -F2L
                    !
                ENDIF
                IF (LDBG) WRITE (LUNDBG, *) K, 'Qss=0 ', RES(K)
                !@@@
            ENDIF
            !
            !------ Kutta condition
            K = KSYSKUT(IEL)
            IF (K>0) THEN
                IF (LSYS) THEN
                    !
                    IF (LBODY(IEL)) THEN
                        !------------ closed body
                        CALL SWDOTS(IEL, DS1SW, DS2SW, DN1SW, DN2SW)
                        JP1 = IPFRST(IEL)
                        JP2 = IPLAST(IEL)
                        !
                        !@@@
                        ds1sw = -1.0
                        ds2sw = 1.0
                        dn1sw = 0.
                        dn2sw = 0.
                        !------------ basic Kutta condition
                        RES(K) = DS2SW * GAM(JP2) - DS1SW * GAM(JP1)              &
                                & + DN2SW * SIG(JP2) - DN1SW * SIG(JP1)
                        !cc     &               + CIRDOT*2.0/(GAM(JP1)-GAM(JP2))
                        !
                        J1 = JSYSGAM(JP1)
                        J2 = JSYSGAM(JP2)
                        SYS(K, J1) = -DS1SW
                        SYS(K, J2) = DS2SW
                        !
                        IF (JP1>=IP1 .AND. JP1<=IP2) THEN
                            J1 = JAICSIG(JP1)
                            AICSIG(K, J1) = -DN1SW
                        ENDIF
                        IF (JP2>=IP1 .AND. JP2<=IP2) THEN
                            !cc               J2 = JAICGAM(JP2)   !!! bug
                            J2 = JAICSIG(JP2)
                            AICSIG(K, J2) = DN2SW
                        ENDIF
                        !
                        IF (IEWAKE(IEL)/=0) THEN
                            !------------- add wake panel contribution
                            IELW = IEWAKE(IEL)
                            JPW = IPFRST(IELW)
                            RES(K) = RES(K) - GAM(JPW)
                            RES(K) = RES(K) - GTH(JPW)      !!HHY
                            IF (JPW>=IP1 .AND. JP2<=IP2) THEN
                                JW = JAICGAM(JPW)
                                AICGAM(K, JW) = -1.0
                                !                JW = JAICGTH(JPW)
                                !                AICGTH(K,JW) = -1.0
                            ENDIF
                        ENDIF
                        !
                    ELSE
                        !------------ zero-thickness plate
                        JP = IPLAST(IEL)
                        !
                        RES(K) = GAM(JP)
                        !
                        J = JSYSGAM(JP)
                        SYS(K, J) = 1.0
                        !
                        IF (IEWAKE(IEL)/=0) THEN
                            !------------- add wake panel contribution
                            IELW = IEWAKE(IEL)
                            JPW = IPFRST(IELW)
                            RES(K) = RES(K) - GAM(JPW)
                            RES(K) = RES(K) - GTH(JPW)      !!HHY
                            !
                            IF (JPW>=IP1 .AND. JP2<=IP2) THEN
                                JW = JAICGAM(JPW)
                                AICGAM(K, JW) = -1.0
                                !                JW = JAICGTH(JPW)
                                !                AICGTH(K,JW) = -1.0
                            ENDIF
                        ENDIF
                        !
                    ENDIF
                ENDIF
                IF (LDBG) WRITE (LUNDBG, *) K, 'Kutta ', RES(K)
                !@@@
            ENDIF
            !
            !------ normal-velocity QNDOF constraint
            K = KSYSQNC(IEL)
            IF (LSYS .AND. K>0) THEN
    
                J = JSYSQNC(IEL)
                RES(K) = QNDOF(IEL)
                !CC  &          - QNDOFSP(IEL)
                SYS(K, J) = 1.0
                IF (LDBG) WRITE (LUNDBG, *) K, 'Qn=0 ', RES(K)
                !@@@
            ENDIF
            !
            !------ zero internal tangential velocity constraint
            K = KSYSQTT(IEL)
            IF (K>0) THEN
                !
                IC = NCX + IEL
                !
                IF (LSYS) THEN
    
                    !cc debug stuff to fiddle with net massflow from QNdof
                    !cc            write(*,*) 'enter Qtt'   !@@@
                    !cc            read(*,*) qttf    !@@@
                    qttf = 0.0
    
                    !-------- flow-tangency residual
                    RES(K) = ANC(1, IC) * QC(1, IC) + ANC(2, IC) * QC(2, IC)         &
                            & - qttf * QINF                              !@@@
                    !
                    IF (LDBG) WRITE (LUNDBG, *) K, 'Qtt=0 ', RES(K)
                    !@@@
    
                    DO JP = 1, NPTOT
                        !---------- dRES/dgamma_surface
                        !-          J = system column index of GAM(JP)  (0 if GAM is on r.h.s.)
                        J = JSYSGAM(JP)
                        SYS(K, J) = ANC(1, IC) * QC_GAM(1, JP, IC) + ANC(2, IC)      &
                                & * QC_GAM(2, JP, IC)
                    ENDDO
                ENDIF
                !
                !-------- dRES/dQinf
                IF (LQFF) AICQFF(K, 1) = ANC(1, IC) - qttf
                !@@@
                !
                !----------------------------------------------------------
                IF (LVWK) THEN
                    !c          write(19,*) 'icx ',ic
                    !------- dRES/dBCirc rotor vortex wake influence
                    DO JP = 1, NPTOT
                        IR = IP2IR(JP)
                        IF (IR/=0) THEN
                            !---- we are on a vortex wake element
                            !
                            !---- accumulate wake influences by wake rows using wake index
                            J = JAICVWK(IR)
                            AICVWK(K, J) = AICVWK(K, J)                          &
                                    & + (ANC(1, IC) * QC_GTH(1, JP, IC)         &
                                            & + ANC(2, IC) * QC_GTH(2, JP, IC)) / YP(JP)
                            !c              write(19,*) 'jp,ir,j ',jp,ir,j
                            !
                            !---- set wake point influence individually
                            J = JAICGTH(JP)
                            AICGTH(K, J) = ANC(1, IC) * QC_GTH(1, JP, IC) + ANC(2, IC)&
                                    & * QC_GTH(2, JP, IC)
                        ENDIF
                    ENDDO
                ENDIF
                !
                !----------------------------------------------------------
                DO JP = IP1, IP2
                    !--------- dRES/dgamma, dRES/dsigma
                    !-         J = r.h.s. column index  (0 if GAM,SIG is a l.h.s. unknown)
                    J = JAICGAM(JP)
                    AICGAM(K, J) = ANC(1, IC) * QC_GAM(1, JP, IC) + ANC(2, IC)      &
                            & * QC_GAM(2, JP, IC)
                    !
                    J = JAICSIG(JP)
                    AICSIG(K, J) = ANC(1, IC) * QC_SIG(1, JP, IC) + ANC(2, IC)      &
                            & * QC_SIG(2, JP, IC)
                ENDDO
            ENDIF
            !
        ENDDO
        !
        !---- set zero vortex strengths explicitly (if selected)
        IF (LSYS) THEN
            DO JP = IP1, IP2
                K = KSYSGAM(JP)
                IF (K>0) THEN
                    J = JSYSGAM(JP)
                    RES(K) = GAM(JP)
                    SYS(K, J) = 1.0
                    IF (LDBG) WRITE (LUNDBG, *) K, 'gam=0 ', RES(K)
                    !@@@
                ENDIF
            ENDDO
        ENDIF
        !
        IF (LSYS) THEN
            LGSYS = .TRUE.
            LGAMU = .FALSE.
            LQGIC = .FALSE.
        ENDIF
        !
    END SUBROUTINE GSYS
    !*==SYSPQ.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GSYS
    
    
    
    SUBROUTINE SYSPQ
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: IEL, IP, IPS1, IPS2, ISEG, JSYS, K, KSYS
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------
        !     Modifies pointers from SYSP to solve
        !     mixed-inverse problem
        !-------------------------------------------
        !
        IF (LDBG) THEN
            WRITE (*, *) 'Modifying system pointers for mixed-inverse...'
            WRITE (LUNDBG, *)                                               &
                    &'SYSPQ modify system pointers for mixed inverse'
            WRITE (LUNDBG, *) '   NSYS =', NSYS, ' (before augmentation)'
        ENDIF
        !
        !---- -999 indicates un-enforced control point
        !----    0 indicates undefined variable
        DO ISEG = 0, NSEG
            DO K = 1, 4
                KSYSDNS(K, ISEG) = -999
                JSYSQSP(K, ISEG) = 0
            ENDDO
        ENDDO
        !
        !---- assume that no panel node lies inside an inverse segment
        DO IP = 1, NPTOT
            KSEGP(IP) = 0
        ENDDO
        !
        !---- initialize system-row and system-column counters to current count
        KSYS = NSYS
        JSYS = NSYS
        !
        DO ISEG = 1, NSEG
            IPS1 = IPSEG1(ISEG)
            IPS2 = IPSEG2(ISEG)
            !
            !------ for nodes inside target segment, change variable from dgam to dn
            DO IP = IPS1, IPS2
                JSYSDXY(IP) = JSYSGAM(IP)
                JSYSGAM(IP) = 0
                !
                KSEGP(IP) = ISEG
            ENDDO
            !
            IEL = IELSEG(ISEG)
            !
            IF (LNFIX1(ISEG) .OR. LNFIX2(ISEG)) THEN
                !------- use constant term if at least one endpoint is fixed
                JSYS = JSYS + 1
                JSYSQSP(1, ISEG) = JSYS
                !
                KSYS = KSYS + 1
                IF (LNFIX1(ISEG)) THEN
                    KSYSDNS(1, ISEG) = KSYS
                ELSE
                    KSYSDNS(2, ISEG) = KSYS
                ENDIF
            ENDIF
            !
            IF (LNFIX1(ISEG) .AND. LNFIX2(ISEG)) THEN
                !------- also use linear term if the other endpoint is also fixed
                JSYS = JSYS + 1
                JSYSQSP(2, ISEG) = JSYS
                !
                KSYS = KSYS + 1
                IF (LNFIX2(ISEG)) THEN
                    KSYSDNS(2, ISEG) = KSYS
                ELSE
                    KSYSDNS(1, ISEG) = KSYS
                ENDIF
            ENDIF
            !
            IF (LDFIX1(IEL)) THEN
                !------- match slope on the left
                KSYS = KSYS + 1
                JSYS = JSYS + 1
                KSYSDNS(3, ISEG) = KSYS
                JSYSQSP(3, ISEG) = JSYS
            ENDIF
            !
            IF (LDFIX2(IEL)) THEN
                !------- match slope on the right
                KSYS = KSYS + 1
                JSYS = JSYS + 1
                KSYSDNS(4, ISEG) = KSYS
                JSYSQSP(4, ISEG) = JSYS
            ENDIF
            !
        ENDDO
        !
        !---- new augmented-system size
        NSYS = KSYS
        IF (LDBG) WRITE (LUNDBG, *) '   NSYS =', NSYS, &
                &' augmented for inverse'
        !
        IF (KSYS/=JSYS) WRITE (*, *) 'Error in SYSPQ, KSYS JSYS =', &
                & KSYS, JSYS
        !
        !---- standard direct-problem pointers no longer valid
        LSYSP = .FALSE.
        !
    END SUBROUTINE SYSPQ
    !*==QSYS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SYSPQ
    
    
    SUBROUTINE QSYS
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        REAL :: DS1M, DS1O, DS2M, DS2O, DSM, DSP, F1L, F1M, &
                & F1O, F2L, F2M, F2O, FRM, FRO, FRP, QTTF, &
                & RES_ANC1, RES_ANC2, RES_DX, RES_DY, RES_GAM, &
                & RES_GAMO, RES_GAMP, RES_QC1, RES_QC2, RES_SIG, &
                & RES_SIGO, RES_SIGP, RES_UINF, RES_XP, RES_XPO, &
                & RES_XPP, RES_YP, RES_YPO, RES_YPP
        INTEGER :: IC, IEL, IELW, ISEG, J, JO, JP, JP1, JP1L, &
                & JP1M, JP1O, JP2, JP2L, JP2M, JP2O, JPM, JPO, &
                & JPP, K, KQSP
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------
        !     Sets up Newton system for mixed-inverse problem.
        !     The equation residuals are the same as the
        !     direct problem set up in GSYS.  The GAM variables
        !     in the inverse segments ISEG = 1..NSEG are replaced
        !     by DNP normal-movement variables.
        !
        !     In addition, the system is augmented by 1,2,3, or 4
        !     geometry-regularity conditions per inverse segment.
        !
        !     Unlike the direct problem, this system is nonlinear,
        !     so the Jacobian matrix depends on the current solution.
        !-----------------------------------------------------------
        !
        IF (NSYS==0) RETURN
        !
        IF (LDBG) THEN
            WRITE (*, *) 'Setting up mixed-inverse Newton system...'
            WRITE (LUNDBG, *) 'QSYS mixed-inverse problem system setup'
            WRITE (LUNDBG, *) '   NSYS =', NSYS
        ENDIF
        !
        !---- clear system arrays
        DO K = 1, NSYS
            RES(K) = 0.
            DO J = 0, NSYS
                SYS(K, J) = 0.
            ENDDO
            DO J = 0, NAICGAM
                AICGAM(K, J) = 0.
            ENDDO
            DO J = 0, NAICSIG
                AICSIG(K, J) = 0.
            ENDDO
            AICQFF(K, 1) = 0.
        ENDDO
        !
        !---- set flow tangency at all control points
        DO IC = 1, NCTOT
            K = KSYSVNC(IC)
            !
            !------ skip un-enforced control points
            IF (K>=0) THEN
                !
                !------- flow tangency at control point IC
                RES(K) = ANC(1, IC) * QC(1, IC) + ANC(2, IC) * QC(2, IC)
                !
                RES_QC1 = ANC(1, IC)
                RES_QC2 = ANC(2, IC)
                RES_ANC1 = QC(1, IC)
                RES_ANC2 = QC(2, IC)
                RES_UINF = ANC(1, IC)
                !
                !------- derivatives w.r.t. dn movement of the two nodes of panel IC
                JPO = IPCO(IC)
                JPP = IPCP(IC)
                RES_DX = RES_ANC1 * ANC_DXY(1, 1, IC) + RES_ANC2 * ANC_DXY(2, 1, IC)
                RES_DY = RES_ANC1 * ANC_DXY(1, 2, IC) + RES_ANC2 * ANC_DXY(2, 2, IC)
                RES_XPO = -RES_DX
                RES_YPO = -RES_DY
                RES_XPP = RES_DX
                RES_YPP = RES_DY
                !
                JO = JSYSDXY(JPO)
                JP = JSYSDXY(JPP)
                SYS(K, JO) = SYS(K, JO) + RES_XPO * ANP(1, JPO)                  &
                        & + RES_YPO * ANP(2, JPO)
                SYS(K, JP) = SYS(K, JP) + RES_XPP * ANP(1, JPP)                  &
                        & + RES_YPP * ANP(2, JPP)
                !
                !------- derivatives w.r.t. vorticity and position of all nodes
                DO JP = 1, NPTOT
                    RES_GAM = RES_QC1 * QC_GAM(1, JP, IC)                        &
                            & + RES_QC2 * QC_GAM(2, JP, IC)
                    RES_SIG = RES_QC1 * QC_SIG(1, JP, IC)                        &
                            & + RES_QC2 * QC_SIG(2, JP, IC)
                    !
                    RES_XP = RES_QC1 * QC_XP(1, JP, IC) + RES_QC2 * QC_XP(2, JP, IC)
                    RES_YP = RES_QC1 * QC_YP(1, JP, IC) + RES_QC2 * QC_YP(2, JP, IC)
                    !
                    CALL GXYLIN(K, JP, RES_GAM, RES_SIG, RES_XP, RES_YP)
                ENDDO
                !
                IF (IELDOFC(IC)>0) THEN
                    !-------- flow tangency contribution of normal-velocity DOF
                    IEL = IELDOFC(IC)
                    RES(K) = RES(K) - QNDOF(IEL)
                    !
                    J = JSYSQNC(IEL)
                    SYS(K, J) = -1.0
                ENDIF
                IF (LDBG) WRITE (LUNDBG, *) K, 'V.n=Qn ', RES(K)
                !@@@
                !
                AICQFF(K, 1) = RES_UINF
            ENDIF
        ENDDO
        !
        !---- set gamma or dn-continuity at all zero-length panels
        DO IC = 1, NCTOT
            K = KSYSGMG(IC)
            !
            !------ skip if there's no constraint
            IF (K>=0) THEN
                !
                JPO = IPCO(IC)
                JPP = IPCP(IC)
                !
                IF (JSYSGAM(JPO)>0 .AND. JSYSGAM(JPP)>0) THEN
                    RES(K) = GAM(JPP) - GAM(JPO)
                    RES_GAMP = 1.0
                    RES_GAMO = -1.0
                    RES_SIGP = 0.
                    RES_SIGO = 0.
                    RES_XPP = 0.
                    RES_XPO = 0.
                    RES_YPP = 0.
                    RES_YPO = 0.
                ELSE
                    !CC      RES(K) = DN(JPP) - DN(JPO)
                    RES(K) = 0.
                    RES_GAMP = 0.
                    RES_GAMO = 0.
                    RES_XPP = ANP(1, JPP)
                    RES_XPO = ANP(1, JPO)
                    RES_YPP = ANP(2, JPP)
                    RES_YPO = ANP(2, JPO)
                ENDIF
                IF (LDBG) WRITE (LUNDBG, *) K, 'dGam=0 ', RES(K)
                !@@@
                !
                CALL GXYLIN(K, JPO, RES_GAMO, RES_SIGO, RES_XPO, RES_YPO)
                CALL GXYLIN(K, JPP, RES_GAMP, RES_SIGP, RES_XPP, RES_YPP)
            ENDIF
        ENDDO
        !
        DO IEL = 1, NEL
            !------ zero net curvature of gamma at TE (triangular loading towards TE)
            K = KSYSGSS(IEL)
            IF (K>0) THEN
    
                JP1O = IPFRST(IEL) + 2
                JP1M = IPFRST(IEL) + 1
                JP1L = IPFRST(IEL)
                JP2O = IPLAST(IEL)
                JP2M = IPLAST(IEL) - 1
                JP2L = IPLAST(IEL) - 2
                !
                DS1O = SP(JP1O) - SP(JP1M)
                DS1M = SP(JP1M) - SP(JP1L)
                DS2O = SP(JP2O) - SP(JP2M)
                DS2M = SP(JP2M) - SP(JP2L)
                !
                !-------- 1st derivative
                !          F1O = 1.0/DS1O + 1.0/(DS1O+DS1M)
                !          F1L = DS1O/((DS1O+DS1M)*DS1M)
                !          F1M = -(F1O+F1L)
                !          F2O = 1.0/DS2O + 1.0/(DS2O+DS2M)
                !          F2L = DS2O/((DS2O+DS2M)*DS2M)
                !          F2M = -(F2O+F2L)
                !
                !-------- 2nd derivative
                F1O = 1.0 / DS1O
                F1M = -(1.0 / DS1O + 1.0 / DS1M)
                F1L = 1.0 / DS1M
                F2O = 1.0 / DS2O
                F2M = -(1.0 / DS2O + 1.0 / DS2M)
                F2L = 1.0 / DS2M
                !
                IF (LV1ZR(IEL)) THEN
                    RES(K) = -F2O * GAM(JP2O) - F2M * GAM(JP2M) - F2L * GAM(JP2L)
                    CALL GXYLIN(K, JP2O, -F2O, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP2M, -F2M, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP2L, -F2L, 0.0, 0.0, 0.0)
                    !
                ELSEIF (LV2ZR(IEL)) THEN
                    RES(K) = F1O * GAM(JP1O) + F1M * GAM(JP1M) + F1L * GAM(JP1L)
                    CALL GXYLIN(K, JP1O, F1O, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP1M, F1M, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP1L, F1L, 0.0, 0.0, 0.0)
                    !
                ELSE
                    RES(K) = F1O * GAM(JP1O) + F1M * GAM(JP1M) + F1L * GAM(JP1L)   &
                            & - F2O * GAM(JP2O) - F2M * GAM(JP2M) - F2L * GAM(JP2L)
                    CALL GXYLIN(K, JP1O, F1O, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP1M, F1M, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP1L, F1L, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP2O, -F2O, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP2M, -F2M, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP2L, -F2L, 0.0, 0.0, 0.0)
                    !
                ENDIF
                IF (LDBG) WRITE (LUNDBG, *) K, 'Qss=0 ', RES(K)
                !@@@
            ENDIF
            !
            !------ Kutta condition
            K = KSYSKUT(IEL)
            IF (K>0) THEN
                !
                IF (LBODY(IEL)) THEN
                    JP1 = IPFRST(IEL)
                    JP2 = IPLAST(IEL)
                    RES(K) = GAM(JP1) + GAM(JP2)
                    CALL GXYLIN(K, JP1, 1.0, 0.0, 0.0, 0.0)
                    CALL GXYLIN(K, JP2, 1.0, 0.0, 0.0, 0.0)
                    !
                ELSE
                    JP = IPLAST(IEL)
                    RES(K) = GAM(JP)
                    CALL GXYLIN(K, JP, 1.0, 0.0, 0.0, 0.0)
                    !
                ENDIF
                !
                IF (IEWAKE(IEL)/=0) THEN
                    IELW = IEWAKE(IEL)
                    JP = IPFRST(IELW)
                    RES(K) = RES(K) - GAM(JP)
                    !           RES(K) = RES(K) - GTH(JP)  !!! HHY
                    CALL GXYLIN(K, JP, -1.0, 0.0, 0.0, 0.0)
                ENDIF
                IF (LDBG) WRITE (LUNDBG, *) K, 'Kutta ', RES(K)
                !@@@
            ENDIF
            !
            !------ normal-velocity DOF
            K = KSYSQNC(IEL)
            IF (K>0) THEN
                !
                J = JSYSQNC(IEL)
                RES(K) = QNDOF(IEL)
                !CC  &          - QNDOFSP(IEL)
                SYS(K, J) = 1.0
                IF (LDBG) WRITE (LUNDBG, *) K, 'Qn=0 ', RES(K)
                !@@@
            ENDIF
            !
            !------ zero internal tangential velocity constraint
            K = KSYSQTT(IEL)
            IF (K>0) THEN
                !
                IC = NCX + IEL
                !
                !            write(*,*) 'enter Qtan/Qinf'   !@@@
                !            read(*,*) qttf    !@@@
                qttf = 0.0
    
                !-------- flow-tangency residual
                RES(K) = ANC(1, IC) * QC(1, IC) + ANC(2, IC) * QC(2, IC) - qttf * QINF
                !@@@
                !
                IF (LDBG) WRITE (LUNDBG, *) K, 'Qtan=0 ', RES(K)
                !@@@
                RES_QC1 = ANC(1, IC)
                RES_QC2 = ANC(2, IC)
                RES_ANC1 = QC(1, IC)
                RES_ANC2 = QC(2, IC)
                RES_UINF = ANC(1, IC)
                !
                !C-------- derivatives w.r.t. dn movement of the two nodes of panel IC
                !          JPO = IPCO(IC)
                !          JPP = IPCP(IC)
                !          RES_DX = RES_ANC1*ANC_DXY(1,1,IC) + RES_ANC2*ANC_DXY(2,1,IC)
                !          RES_DY = RES_ANC1*ANC_DXY(1,2,IC) + RES_ANC2*ANC_DXY(2,2,IC)
                !          RES_XPO = -RES_DX
                !          RES_YPO = -RES_DY
                !          RES_XPP =  RES_DX
                !          RES_YPP =  RES_DY
                !
                !          JO = JSYSDXY(JPO)
                !          JP = JSYSDXY(JPP)
                !          SYS(K,JO) = SYS(K,JO) + RES_XPO*ANP(1,JPO) + RES_YPO*ANP(2,JPO)
                !          SYS(K,JP) = SYS(K,JP) + RES_XPP*ANP(1,JPP) + RES_YPP*ANP(2,JPP)
                !
                !-------- derivatives w.r.t. vorticity and position of all nodes
                DO JP = 1, NPTOT
                    RES_GAM = RES_QC1 * QC_GAM(1, JP, IC)                        &
                            & + RES_QC2 * QC_GAM(2, JP, IC)
                    RES_SIG = RES_QC1 * QC_SIG(1, JP, IC)                        &
                            & + RES_QC2 * QC_SIG(2, JP, IC)
                    !
                    RES_XP = RES_QC1 * QC_XP(1, JP, IC) + RES_QC2 * QC_XP(2, JP, IC)
                    RES_YP = RES_QC1 * QC_YP(1, JP, IC) + RES_QC2 * QC_YP(2, JP, IC)
                    !
                    CALL GXYLIN(K, JP, RES_GAM, RES_SIG, RES_XP, RES_YP)
                ENDDO
                !
            ENDIF
        ENDDO
        !
        !
        !---- set zero vortex strengths explicitly
        DO JP = 1, NPTOT
            K = KSYSGAM(JP)
            !
            !------ skip non-specified vortices
            IF (K>=0) THEN
                !
                RES(K) = GAM(JP)
                CALL GXYLIN(K, JP, 1.0, 0.0, 0.0, 0.0)
                IF (LDBG) WRITE (LUNDBG, *) K, 'gam=0 ', RES(K)
                !@@@
            ENDIF
        ENDDO
        !
        !
        !---- segment-fixing conditions
        DO ISEG = 1, NSEG
            K = KSYSDNS(1, ISEG)
            IF (K>0) THEN
    
                JP = IPSEG1(ISEG)
                !CC      RES(K) = DN(JP)
                RES(K) = 0.
                J = JSYSDXY(JP)
                SYS(K, J) = 1.0
                IF (LDBG) WRITE (LUNDBG, *) K, 'dn1=0 ', RES(K)
                !@@@
            ENDIF
            !
            K = KSYSDNS(2, ISEG)
            IF (K>0) THEN
    
                JP = IPSEG2(ISEG)
                !CC      RES(K) = DN(JP)
                RES(K) = 0.
                J = JSYSDXY(JP)
                SYS(K, J) = 1.0
                IF (LDBG) WRITE (LUNDBG, *) K, 'dn2=0 ', RES(K)
                !@@@
            ENDIF
            !
            K = KSYSDNS(3, ISEG)
            IF (K>0) THEN
    
                !
                !         JP = IPSEG1(ISEG) + 1
                !CCC      RES(K) = DN(JP)
                !         RES(K) = 0.
                !         J = JSYSDXY(JP)
                !         SYS(K,J) = 1.0
                !
                KQSP = KQTARG
                JPM = IPSEG1(ISEG) - 1
                JPO = IPSEG1(ISEG)
                JPP = IPSEG1(ISEG) + 1
                DSM = SSPEC(JPO) - SSPEC(JPM)
                DSP = SSPEC(JPP) - SSPEC(JPO)
                FRM = 1.0 / DSM
                FRO = -(1.0 / DSM + 1.0 / DSP)
                FRP = 1.0 / DSP
                RES(K) = FRM * GAM(JPM) + FRO * GAM(JPO) + FRP * GAM(JPP)         &
                        & - FRM * QSPEC(JPM, KQSP) - FRO * QSPEC(JPO, KQSP)        &
                        & - FRP * QSPEC(JPP, KQSP)
                CALL GXYLIN(K, JPM, FRM, 0.0, 0.0, 0.0)
                CALL GXYLIN(K, JPO, FRO, 0.0, 0.0, 0.0)
                CALL GXYLIN(K, JPP, FRP, 0.0, 0.0, 0.0)
                IF (LDBG) WRITE (LUNDBG, *) K, 'gam_ss1=0 ', RES(K)
                !@@@
            ENDIF
            !
            K = KSYSDNS(4, ISEG)
            IF (K>0) THEN
                !
                !
                !         JP = IPSEG2(ISEG) - 1
                !CCC      RES(K) = DN(JP)
                !         RES(K) = 0.
                !         J = JSYSDXY(JP)
                !         SYS(K,J) = 1.0
                !
                KQSP = KQTARG
                JPM = IPSEG2(ISEG) - 1
                JPO = IPSEG2(ISEG)
                JPP = IPSEG2(ISEG) + 1
                DSM = SSPEC(JPO) - SSPEC(JPM)
                DSP = SSPEC(JPP) - SSPEC(JPO)
                FRM = 1.0 / DSM
                FRO = -(1.0 / DSM + 1.0 / DSP)
                FRP = 1.0 / DSP
                RES(K) = FRM * GAM(JPM) + FRO * GAM(JPO) + FRP * GAM(JPP)         &
                        & - FRM * QSPEC(JPM, KQSP) - FRO * QSPEC(JPO, KQSP)        &
                        & - FRP * QSPEC(JPP, KQSP)
                CALL GXYLIN(K, JPM, FRM, 0.0, 0.0, 0.0)
                CALL GXYLIN(K, JPO, FRO, 0.0, 0.0, 0.0)
                CALL GXYLIN(K, JPP, FRP, 0.0, 0.0, 0.0)
                IF (LDBG) WRITE (LUNDBG, *) K, 'gam_ss2=0 ', RES(K)
                !@@@
            ENDIF
        ENDDO
        !
        LQSYS = .TRUE.
        !
    END SUBROUTINE QSYS
    !*==GXYLIN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QSYS
    
    
    
    
    SUBROUTINE GXYLIN(K, JP, RES_GAM, RES_SIG, RES_XP, RES_YP)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: JP, K
        REAL :: RES_GAM, RES_SIG, RES_XP, RES_YP
        !
        ! Local variables
        !
        INTEGER :: ISEG, J, J1, J2, J3, J4
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        J = JSYSGAM(JP)
        SYS(K, J) = SYS(K, J) + RES_GAM
        !
        J = JAICGAM(JP)
        AICGAM(K, J) = AICGAM(K, J) + RES_GAM
        !
        J = JAICSIG(JP)
        AICSIG(K, J) = AICSIG(K, J) + RES_SIG
        !
        ISEG = KSEGP(JP)
        J1 = JSYSQSP(1, ISEG)
        J2 = JSYSQSP(2, ISEG)
        J3 = JSYSQSP(3, ISEG)
        J4 = JSYSQSP(4, ISEG)
        SYS(K, J1) = SYS(K, J1) + RES_GAM * FSPEC(1, JP)
        SYS(K, J2) = SYS(K, J2) + RES_GAM * FSPEC(2, JP)
        SYS(K, J3) = SYS(K, J3) + RES_GAM * FSPEC(3, JP)
        SYS(K, J4) = SYS(K, J4) + RES_GAM * FSPEC(4, JP)
        !
        J = JSYSDXY(JP)
        SYS(K, J) = SYS(K, J) + RES_XP * ANP(1, JP) + RES_YP * ANP(2, JP)
        !
    END SUBROUTINE GXYLIN
    !*==GUCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GXYLIN
    
    
    
    SUBROUTINE GUCALC(LGUQI, LVWK, IP1, IP2)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IP1, IP2
        LOGICAL :: LGUQI, LVWK
        !
        ! Local variables
        !
        INTEGER :: IP, IR, J
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Sets up V.n residuals at all control points.
        !     Calculates unit-freestream vortex sheet strengths
        !     GAMU by back-substituting into already-factored
        !     Jacobian matrix.
        !-----------------------------------------------------
        !
        IF (LDBG) WRITE (*, *)                                           &
                &'Generating unit vorticity distributions...'
        !
        !---------------------------------------------------------------------
        !---- assign unit-solution indices
        IF (LGUQI) THEN
            IF (IUQINF==0) THEN
                !------ freestream x-velocity
                IF (LDBG) WRITE (*, *) 'Generating solution for unit Uinf'
                CALL GUMAKE(IUQINF, AICQFF(1, 1))
            ENDIF
        ENDIF
        !
        IF (LVWK) THEN
            DO IR = 1, NRP
                IF (IUVWK(IR)==0) THEN
                    !------ rotor circulation, radial station IR
                    IF (LDBG) WRITE (*, *)                                  &
                            &'Generating solution for unit GAM IR '&
                            &, IR
                    J = JAICVWK(IR)
                    CALL GUMAKE(IUVWK(IR), AICVWK(1, J))
                    CALL GUVWK(IR)
                ENDIF
            ENDDO
        ENDIF
        !
        DO IP = IP1, IP2
            J = JAICGAM(IP)
            IF (IUGAM(IP)==0 .AND. J>0) THEN
                IF (LDBG) WRITE (*, *) 'Generating solution for unit gamma'&
                        &, IP
                CALL GUMAKE(IUGAM(IP), AICGAM(1, J))
                GAMU(IP, IUGAM(IP)) = 1.0
            ENDIF
            !
            J = JAICSIG(IP)
            IF (IUSIG(IP)==0 .AND. J>0) THEN
                IF (LDBG) WRITE (*, *) 'Generating solution for unit sigma'&
                        &, IP
                CALL GUMAKE(IUSIG(IP), AICSIG(1, J))
                SIGU(IP, IUSIG(IP)) = 1.0
            ENDIF
        ENDDO
        !
    END SUBROUTINE GUCALC
    !*==GUMAKE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GUCALC
    
    
    
    SUBROUTINE GUMAKE(IU, AIC)
        use m_gauss, only: baksub
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IU
        REAL, DIMENSION(*) :: AIC
        !
        ! Local variables
        !
        INTEGER :: IEL, IP, J, KU
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     Finds the index IU of the first available unit r.h.s. vector
        !     Generates GAMU(.IU),SIGU(.IU),QNDOFU(.IU) unit solutions.
        !-------------------------------------------------------------------
        !
        !---- find first available vector index KU
        DO KU = 1, NU
            IF (.NOT.LUSET(KU)) GOTO 10
        ENDDO
        !
        !---- no free vectors found... add new one to list
        NU = NU + 1
        IF (NU>NUX) WRITE (*, *)                                         &
                &'GUMAKE: Array limit exceeded. Increase NUX to'&
                &, NU
        !
        KU = MIN(NU, NUX)
        !
        !-----------------------------------
        !
        !---- set Jacobian vector
        10   DO J = 1, NSYS
            RES(J) = AIC(J)
        ENDDO
        !
        !---- calculate response
        IF (NSYS>0) CALL BAKSUB(NSYSX, NSYS, SYS(1, 1), RES(1))
        RES(0) = 0.
        !
        !---- first clear entire response vectors
        DO IP = 1, NPTOT
            GAMU(IP, KU) = 0.
            SIGU(IP, KU) = 0.
        ENDDO
        DO IEL = 1, NEL
            QNDOFU(IEL, KU) = 0.
        ENDDO
        !
        !---- store response of surface gamma
        DO IP = 1, NPTOT
            J = JSYSGAM(IP)
            GAMU(IP, KU) = -RES(J)
        ENDDO
        !
        !---- store response of normal-velocity DOFs
        DO IEL = 1, NEL
            J = JSYSQNC(IEL)
            QNDOFU(IEL, KU) = -RES(J)
        ENDDO
        !
        LUSET(KU) = .TRUE.
        !
        IU = KU
        !
    END SUBROUTINE GUMAKE
    !*==GSOLVE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GUMAKE
    
    
    
    SUBROUTINE GSOLVE
        use m_gauss, only: baksub
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        REAL, DIMENSION(IPX) :: GAMK, SIGK
        INTEGER :: IEL, IP, J, K
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     Does a direct solve from current RHS (knowns)
        !     Generates GAM(.),SIG(.),QNDOF(.)
        !-------------------------------------------------------------------
        !
        !---- set imposed (known) gamma and sigma strengths
        DO IP = 1, NPTOT
            GAMK(IP) = GAMVSP(IP)
            SIGK(IP) = SIGVSP(IP)
        ENDDO
        !
        DO IEL = 1, NEL
            DO IP = IPFRST(IEL), IPLAST(IEL)
                GAMK(IP) = GAMK(IP) + GAMSET(IEL)
                SIGK(IP) = SIGK(IP) + SIGSET(IEL)
            ENDDO
        ENDDO
        !
        !---- set Jacobian vector (RHS)
        DO K = 1, NSYS
            RES(K) = AICQFF(K, 1) * QINF
            DO IP = 1, NPTOT
                J = JAICGTH(IP)
                RES(K) = RES(K) + AICGTH(K, J) * GTH(IP)
            ENDDO
            DO IP = 1, NPTOT
                J = JAICGAM(IP)
                RES(K) = RES(K) + AICGAM(K, J) * GAMK(IP)
            ENDDO
            DO IP = 1, NPTOT
                J = JAICSIG(IP)
                RES(K) = RES(K) + AICSIG(K, J) * SIGK(IP)
            ENDDO
        ENDDO
        !
        !---- calculate solution
        IF (NSYS>0) THEN
            CALL BAKSUB(NSYSX, NSYS, SYS(1, 1), RES(1))
        ELSE
            WRITE (*, *) 'GSOLVE error NSYS=', NSYS
        ENDIF
        RES(0) = 0.
        !
        !---- first clear unknown solution vectors
        DO IP = 1, NPTOT
            GAM(IP) = 0.
            SIG(IP) = 0.
        ENDDO
        DO IEL = 1, NEL
            QNDOF(IEL) = 0.
        ENDDO
        !
        !---- store solution for surface gamma + known gamma values
        !---- store solution for surface sigma + known sigma values
        DO IP = 1, NPTOT
            J = JSYSGAM(IP)
            GAM(IP) = -RES(J) + GAMK(IP)
            SIG(IP) = SIG(IP) + SIGK(IP)
        ENDDO
        !
        !---- store solution for normal-velocity DOFs
        DO IEL = 1, NEL
            J = JSYSQNC(IEL)
            QNDOF(IEL) = -RES(J)
        ENDDO
        !
        !---- Set strength of TE panels (for debugging use)
        CALL SETGSTE
        !
        LGAMA = .TRUE.
        LQCNT = .FALSE.
        !
    END SUBROUTINE GSOLVE
    !*==GSOLVE0.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GSOLVE
    
    
    SUBROUTINE GSOLVE0
        use m_gauss, only: baksub
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        REAL, DIMENSION(IPX) :: GAMK, SIGK
        INTEGER :: IEL, IP, J, K
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     Does a direct solve from current RHS (knowns)
        !     Generates GAM(.),SIG(.),QNDOF(.)
        !-------------------------------------------------------------------
        !
        !---- set known gamma and sigma strengths
        DO IP = 1, NPTOT
            GAMK(IP) = GAMVSP(IP)
            SIGK(IP) = SIGVSP(IP)
        ENDDO
        !
        DO IEL = 1, NEL
            DO IP = IPFRST(IEL), IPLAST(IEL)
                GAMK(IP) = GAMK(IP) + GAMSET(IEL)
                SIGK(IP) = SIGK(IP) + SIGSET(IEL)
            ENDDO
        ENDDO
        !
        !---- set Jacobian vector (RHS)
        DO K = 1, NSYS
            RES(K) = AICQFF(K, 1) * QINF
            !        DO IR = 1, NRP
            !          RES(K) = RES(K) + AICVWK(K,IR)*WBGAM(IR)
            !        ENDDO
            DO IP = 1, NPTOT
                J = JAICGAM(IP)
                RES(K) = RES(K) + AICGAM(K, J) * GAMK(IP)
            ENDDO
            DO IP = 1, NPTOT
                J = JAICSIG(IP)
                RES(K) = RES(K) + AICSIG(K, J) * SIGK(IP)
            ENDDO
        ENDDO
        !
        !---- calculate response
        IF (NSYS>0) CALL BAKSUB(NSYSX, NSYS, SYS(1, 1), RES(1))
        RES(0) = 0.
        !
        !---- first clear unknown solution vectors
        DO IP = 1, NPTOT
            GAM(IP) = 0.
            GTH(IP) = 0.
            SIG(IP) = 0.
        ENDDO
        DO IEL = 1, NEL
            QNDOF(IEL) = 0.
        ENDDO
        !
        !---- store solution for surface gamma + known gamma values
        !---- store solution for surface sigma + known sigma values
        DO IP = 1, NPTOT
            J = JSYSGAM(IP)
            GAM(IP) = -RES(J) + GAMK(IP)
            SIG(IP) = SIG(IP) + SIGK(IP)
            !---- set GTH for wake point from WBGAM for wake row
            !        IF(IP2IR(IP).NE.0) THEN
            !          IR = IP2IR(IP)
            !          GTH(IP) = WBGAM(IR)/YP(IP)
            !        ENDIF
        ENDDO
        !
        !---- store solution for normal-velocity DOFs
        DO IEL = 1, NEL
            J = JSYSQNC(IEL)
            QNDOF(IEL) = -RES(J)
        ENDDO
        !
        !---- Set strength of TE panels (for debugging use)
        CALL SETGSTE
        LGAMA = .TRUE.
        LQCNT = .FALSE.
        !
    END SUBROUTINE GSOLVE0
    !*==GUVWK.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GSOLVE0
    
    
    SUBROUTINE GUVWK(IRADD)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: IRADD
        !
        ! Local variables
        !
        INTEGER :: IR, JP, KU
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     Creates wake gamma "mode" shape corresponding to unit
        !     rotor row IRADD circulation (actually Wt/Wm *delta BGAM /2*pi)
        !
        !     Creates GTHU(.IU) unit vorticity distribution
        !-------------------------------------------------------------------
        !
        !c      WRITE(*,*) 'GAMADDVWK: add rotor wake gammas to unit soln ',IRADD
        !
        !---- store imposed unit rotor WBGAM circulation gamma distribution (on wake
        !     and body nodes) as GTHU
        DO JP = 1, NPTOT
            !---- clear unit response
            KU = IUVWK(IRADD)
            GTHU(JP, KU) = 0.
            !
            IR = IP2IR(JP)
            IF (IR/=0) THEN
                !---- we are on a rotor wake element, use rotor index to locate influence
                IF (IR==IRADD) THEN
                    IF (IUVWK(IR)>0) THEN
                        KU = IUVWK(IR)
                        GTHU(JP, KU) = 1.0 / YP(JP)
                    ELSE
                        WRITE (*, *) 'GUVWK: Unassigned rotor RHS IR =', IR
                    ENDIF
                ENDIF
                !
            ENDIF
        ENDDO
        !
    END SUBROUTINE GUVWK
    !*==GUSUM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GUVWK
    
    
    
    SUBROUTINE GUSUM
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: IEL, IP, IU
        REAL, DIMENSION(0:NUX) :: UVAL
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Superimposes unit solutions for current...
        !       * freestream Qinf
        !       * rotor circulations
        !       * specified wake  gamma,sigma
        !       * specified ring  Gamma,Sigma
        !       * specified point Doublet,Source
        !       * specified source-only elements
        !-----------------------------------------------------
        !
        IF (LDBG) WRITE (*, *)                                           &
                &'Superimposing unit-strength distributions'
        !
        DO IU = 0, NU
            UVAL(IU) = 0.
        ENDDO
        !
        !---- put prescribed strengths into indexed array (inactive quantity has IU=0)
        !
        !----- axisymmetric case
        UVAL(IUQINF) = QINF
        !
        !---- Rotor wake circulation (actually Wt/Wm *delta BGAM /2*pi)
        !      DO IR = 1, NRP
        !        UVAL(IUVWK(IR)) = WBGAM(IR)
        !c       write(18,*) ir,iuvwk(ir),wbgam(ir)
        !      ENDDO
        !
        !---- Additional gamma specified from viscous response (from BL effects)
        DO IP = 1, NPTOT
            !c        write(24,*) IP,IUGAM(IP),GAMVSP(IP)
            UVAL(IUGAM(IP)) = GAMVSP(IP)
            UVAL(IUSIG(IP)) = SIGVSP(IP)
        ENDDO
        !
        DO IEL = 1, NEL
            DO IP = IPFRST(IEL), IPLAST(IEL)
                UVAL(IUGAM(IP)) = UVAL(IUGAM(IP)) + GAMSET(IEL)
                UVAL(IUSIG(IP)) = UVAL(IUSIG(IP)) + SIGSET(IEL)
            ENDDO
        ENDDO
        !
        !---- clear all singularities for accumulation
        DO IP = 1, NPTOT
            GAM(IP) = 0.
            SIG(IP) = 0.
            GTH(IP) = 0.
        ENDDO
        DO IEL = 1, NEL
            QNDOF(IEL) = 0.
        ENDDO
        !
        !---- accumulate unit contributions weighted by prescribed strengths
        DO IU = 1, NU
            !------ no sense accumulating zeros
            IF (UVAL(IU)/=0.0) THEN
                !c          write(20,*) iu,uval(iu)
                !
                DO IP = 1, NPTOT
                    GAM(IP) = GAM(IP) + GAMU(IP, IU) * UVAL(IU)
                    SIG(IP) = SIG(IP) + SIGU(IP, IU) * UVAL(IU)
                    GTH(IP) = GTH(IP) + GTHU(IP, IU) * UVAL(IU)
                    !c          write(20,99) ip,iu,gamu(ip,iu),GTHU(IP,IU)
                ENDDO
                DO IEL = 1, NEL
                    QNDOF(IEL) = QNDOF(IEL) + QNDOFU(IEL, IU) * UVAL(IU)
                ENDDO
            ENDIF
            !
        ENDDO
        !
        !      DO IP = 1, NPTOT
        !        write(24,*) IP,GAM(IP),GAMVSP(IP)
        !      ENDDO
        !
        !---- Set strength of TE panels (for debugging use)
        CALL SETGSTE
        !
        99   FORMAT (2I5, 5(1x, F12.6))
        !
        LGAMA = .TRUE.
        LQCNT = .FALSE.
        !
    END SUBROUTINE GUSUM
    !*==QCSUM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GUSUM
    
    
    
    SUBROUTINE QCSUM
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: IC, JP
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Computes velocities at control points for
        !     current panel and point-singularity strengths.
        !-----------------------------------------------------
        !
        IF (LDBG) WRITE (*, *) 'Computing control-point velocities'
        !
        DO IC = 1, NCTOT
            IF (ICTYPE(IC)==-999) THEN
                WRITE (*, *) '? unassigned control point ', IC
                !
            ELSEIF (ICTYPE(IC)==-1) THEN
                !------- zero-length panel (skip)
                !
            ELSEIF (ICTYPE(IC)/=-2) THEN
                !------- vortex wake panels w/o control points (skip)
                !
                !------- usual velocity contributions for all vortex and source strengths
                QC(1, IC) = 0.
                QC(2, IC) = 0.
                !
                DO JP = 1, NPTOT
                    QC(1, IC) = QC(1, IC) + QC_GAM(1, JP, IC) * GAM(JP)            &
                            & + QC_SIG(1, JP, IC) * SIG(JP)
                    QC(2, IC) = QC(2, IC) + QC_GAM(2, JP, IC) * GAM(JP)            &
                            & + QC_SIG(2, JP, IC) * SIG(JP)
                    !---- Add contribution from rotor wake vorticity
                    QC(1, IC) = QC(1, IC) + QC_GTH(1, JP, IC) * GTH(JP)
                    QC(2, IC) = QC(2, IC) + QC_GTH(2, JP, IC) * GTH(JP)
                ENDDO
                !
                !---- Add freestream flow contribution
                QC(1, IC) = QC(1, IC) + QINF
            ENDIF
            !
        ENDDO
        !
        LQCNT = .TRUE.
        !
    END SUBROUTINE QCSUM
    !*==QCUSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCSUM
    
    
    
    
    SUBROUTINE QCUSET
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: IC, IP, IU, JP
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Computes velocities at control points for
        !     current unit-forcing solutions
        !-----------------------------------------------------
        !
        IF (LDBG) WRITE (*, *)                                           &
                &'Computing control-pt unit-forcing velocities'
        !
        DO IC = 1, NCTOT
            IF (ICTYPE(IC)==-999) THEN
                WRITE (*, *) '? unassigned control point ', IC
                !
            ELSEIF (ICTYPE(IC)==-1) THEN
                !------- zero-length panel (skip)
                !
            ELSEIF (ICTYPE(IC)/=-2) THEN
                !------- vortex wake panel w/o control point (skip)
                !
                !------- usual velocity contributions for all vortex strengths
                DO IU = 1, NU
                    QCU(1, IC, IU) = 0.
                    QCU(2, IC, IU) = 0.
                    !
                    DO JP = 1, NPTOT
                        !c             IF(JSYSGAM(JP).GT.0) THEN
                        !---- Contribution from solution unit responses
                        QCU(1, IC, IU) = QCU(1, IC, IU) + QC_GAM(1, JP, IC)         &
                                & * GAMU(JP, IU)
                        QCU(2, IC, IU) = QCU(2, IC, IU) + QC_GAM(2, JP, IC)         &
                                & * GAMU(JP, IU)
                        !---- Contribution from rotor wake vorticity unit forcing
                        QCU(1, IC, IU) = QCU(1, IC, IU) + QC_GTH(1, JP, IC)         &
                                & * GTHU(JP, IU)
                        QCU(2, IC, IU) = QCU(2, IC, IU) + QC_GTH(2, JP, IC)         &
                                & * GTHU(JP, IU)
                        !c             ENDIF
                    ENDDO
                ENDDO
                !
                !---- Contribution from freestream
                IU = IUQINF
                QCU(1, IC, IU) = QCU(1, IC, IU) + 1.0
                !
                !---- Contribution from specified vortex and source panels
                DO IP = 1, NPTOT
                    IU = IUGAM(IP)
                    QCU(1, IC, IU) = QCU(1, IC, IU) + QC_GAM(1, IP, IC)
                    QCU(2, IC, IU) = QCU(2, IC, IU) + QC_GAM(2, IP, IC)
                    !
                    IU = IUSIG(IP)
                    QCU(1, IC, IU) = QCU(1, IC, IU) + QC_SIG(1, IP, IC)
                    QCU(2, IC, IU) = QCU(2, IC, IU) + QC_SIG(2, IP, IC)
                ENDDO
                !
            ENDIF
        ENDDO
        !
    END SUBROUTINE QCUSET
    !*==QCUSUM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCUSET
    
    
    
    SUBROUTINE QCUSUM
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: IC, IEL, IP, IU
        REAL, DIMENSION(0:NUX) :: UVAL
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Computes velocities at control points for
        !     current panel and point-singularity strengths.
        !-----------------------------------------------------
        !
        IF (LDBG) WRITE (*, *) 'Superimposing unit-strength velocities'
        !
        DO IU = 0, NU
            UVAL(IU) = 0.
        ENDDO
        !
        !---- put prescribed strengths into indexed array (inactive quantity has IU=0)
        !
        !----- freestream
        UVAL(IUQINF) = QINF
        !
        !---- Rotor wake circulation (actually Wt/Wm *delta BGAM /2*pi)
        !      DO IR = 1, NRP
        !        UVAL(IUVWK(IR)) = WBGAM(IR)
        !c       write(18,*) ir,iuvwk(ir),wbgam(ir)
        !      ENDDO
        !
        DO IP = 1, NPTOT
            UVAL(IUGAM(IP)) = GAMVSP(IP)
            UVAL(IUSIG(IP)) = SIGVSP(IP)
        ENDDO
        !---- Specified singularities
        DO IEL = 1, NEL
            DO IP = IPFRST(IEL), IPLAST(IEL)
                UVAL(IUGAM(IP)) = UVAL(IUGAM(IP)) + GAMSET(IEL)
                UVAL(IUSIG(IP)) = UVAL(IUSIG(IP)) + SIGSET(IEL)
            ENDDO
        ENDDO
        !
        !---- accumulate unit contributions weighted by prescribed strengths
        DO IC = 1, NCTOT
            QC(1, IC) = 0.
            QC(2, IC) = 0.
        ENDDO
        !
        DO IU = 1, NU
            !------ no sense accumulating zeros
            IF (UVAL(IU)/=0.0) THEN
                !
                DO IC = 1, NCTOT
                    QC(1, IC) = QC(1, IC) + QCU(1, IC, IU) * UVAL(IU)
                    QC(2, IC) = QC(2, IC) + QCU(2, IC, IU) * UVAL(IU)
                ENDDO
            ENDIF
            !
        ENDDO
        !
        LQCNT = .TRUE.
        !
    END SUBROUTINE QCUSUM
    !*==SETGSTE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCUSUM
    
    
    
    
    SUBROUTINE SETGSTE
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        INTEGER :: IEL, IP1, IP2, K
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Sets TE vortex and source strength
        !     from vortex and source strength on body panels
        !-----------------------------------------------------
        !
        IF (LDBG) WRITE (*, *) 'Setting TE GAMT,SIGT'
        !
        !---- Set TE panel strengths
        !c      write(30,*) 'SETGSTE IEL,GAMT1-2,SIGT1-2'
        !
        DO IEL = 1, NEL
            DO K = 1, 2
                GAMT(K, IEL) = 0.
                SIGT(K, IEL) = 0.
            ENDDO
            !
            !---- Body elements with GAM
            IF (LTPAN(IEL) .AND. NETYPE(IEL)==0) THEN
                IP1 = IPTE1(IEL)
                IP2 = IPTE2(IEL)
                !
                IF (IP1/=0) THEN
                    DO K = 1, 2
                        GAMT(K, IEL) = GAMT(K, IEL) + GAMT_GAM(K, 1, IEL) * GAM(IP1)&
                                & + GAMT_SIG(K, 1, IEL) * SIG(IP1)
                        SIGT(K, IEL) = SIGT(K, IEL) + SIGT_GAM(K, 1, IEL) * GAM(IP1)&
                                & + SIGT_SIG(K, 1, IEL) * SIG(IP1)
                    ENDDO
                ENDIF
                IF (IP2/=0) THEN
                    DO K = 1, 2
                        GAMT(K, IEL) = GAMT(K, IEL) + GAMT_GAM(K, 2, IEL) * GAM(IP2)&
                                & + GAMT_SIG(K, 2, IEL) * SIG(IP2)
                        SIGT(K, IEL) = SIGT(K, IEL) + SIGT_GAM(K, 2, IEL) * GAM(IP2)&
                                & + SIGT_SIG(K, 2, IEL) * SIG(IP2)
                    ENDDO
                ENDIF
            ENDIF
            !
            !---- vortex wake elements with GTH specified
            IF (LTPAN(IEL) .AND. NETYPE(IEL)==7) THEN
                IP1 = IPTE1(IEL)
                IP2 = IPTE2(IEL)
                !
                IF (IP1/=0) THEN
                    DO K = 1, 2
                        GAMT(K, IEL) = GAMT(K, IEL) + GAMT_GAM(K, 1, IEL) * GTH(IP1)
                        SIGT(K, IEL) = SIGT(K, IEL) + SIGT_GAM(K, 1, IEL) * GTH(IP1)
                    ENDDO
                ENDIF
                IF (IP2/=0) THEN
                    DO K = 1, 2
                        GAMT(K, IEL) = GAMT(K, IEL) + GAMT_GAM(K, 2, IEL) * GTH(IP2)
                        SIGT(K, IEL) = SIGT(K, IEL) + SIGT_GAM(K, 2, IEL) * GTH(IP2)
                    ENDDO
                ENDIF
            ENDIF
            !
            !c        write(30,99) iel,GAMT(1,IEL),GAMT(2,IEL),SIGT(1,IEL),SIGT(2,IEL)
            99      FORMAT (i5, 4(1x, f12.6))
        ENDDO
        !
    END SUBROUTINE SETGSTE
    !*==NWDOTS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETGSTE
    
    
    
    SUBROUTINE NWDOTS(IEL, DS1NW, DS2NW, DN1NW, DN2NW)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: DN1NW, DN2NW, DS1NW, DS2NW
        INTEGER :: IEL
        !
        ! Local variables
        !
        REAL, DIMENSION(2) :: ANW
        REAL :: ANWABS
        INTEGER :: IC1, IC2
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        !
        !      IELW = IEWAKE(IEL)
        !      IF(IELW .NE. 0) THEN
        !       ICW = ICFRST(IELW)
        !       ANW(1) = ANC(1,ICW)
        !       ANW(2) = ANC(2,ICW)
        !      ELSE
        ANW(1) = ANC(1, IC2) - ANC(1, IC1)
        ANW(2) = ANC(2, IC2) - ANC(2, IC1)
        ANWABS = SQRT(ANW(1)**2 + ANW(2)**2)
        ANW(1) = ANW(1) / ANWABS
        ANW(2) = ANW(2) / ANWABS
        !      ENDIF
        !
        DS1NW = ANC(2, IC1) * ANW(1) - ANC(1, IC1) * ANW(2)
        DS2NW = ANC(2, IC2) * ANW(1) - ANC(1, IC2) * ANW(2)
        !
        DN1NW = ANC(1, IC1) * ANW(1) + ANC(2, IC1) * ANW(2)
        DN2NW = ANC(1, IC2) * ANW(1) + ANC(2, IC2) * ANW(2)
        !
    END SUBROUTINE NWDOTS
    !*==SWDOTS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! NWDOTS
    
    
    SUBROUTINE SWDOTS(IEL, DS1SW, DS2SW, DN1SW, DN2SW)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: DN1SW, DN2SW, DS1SW, DS2SW
        INTEGER :: IEL
        !
        ! Local variables
        !
        REAL, DIMENSION(2) :: ANW
        REAL :: ANWABS
        INTEGER :: IC1, IC2
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        !
        !      IELW = IEWAKE(IEL)
        !      IF(IELW .NE. 0) THEN
        !       ICW = ICFRST(IELW)
        !       ANW(1) = ANC(1,ICW)
        !       ANW(2) = ANC(2,ICW)
        !      ELSE
        ANW(1) = ANC(1, IC2) - ANC(1, IC1)
        ANW(2) = ANC(2, IC2) - ANC(2, IC1)
        ANWABS = SQRT(ANW(1)**2 + ANW(2)**2)
        ANW(1) = ANW(1) / ANWABS
        ANW(2) = ANW(2) / ANWABS
        !      ENDIF
        !
        DS1SW = ANC(2, IC1) * ANW(2) + ANC(1, IC1) * ANW(1)
        DS2SW = ANC(2, IC2) * ANW(2) + ANC(1, IC2) * ANW(1)
        !
        DN1SW = ANC(1, IC1) * ANW(2) - ANC(2, IC1) * ANW(1)
        DN2SW = ANC(1, IC2) * ANW(2) - ANC(2, IC2) * ANW(1)
        !
    END SUBROUTINE SWDOTS
end module m_system
