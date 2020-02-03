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

    subroutine sysp
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: ic, iel, ip, ir, jaicg, jaics, jaicv, jaicw, &
                & jaicz, jp, jsys, ksys, ndof
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------
        !     Sets up pointers associating....
        !        matrix rows with control points,
        !        matrix columns with unknowns
        !-------------------------------------------
        !
        if (ldbg) write (*, *) 'Generating flow-tangency system pointers'
        !
        !---- set what type of equation gets assigned to each control point,
        !-    and which unknown vortex gets associated with each control point
        !
        !---- -999 indicates un-enforced control point
        !----    0 indicates undefined variable
        do ic = 1, nctot
            ksysvnc(ic) = -999
            ksysgmg(ic) = -999
            ieldofc(ic) = 0
        enddo
        do ip = 1, nptot
            ksysgam(ip) = -999
            jsysgam(ip) = 0
            jsysdxy(ip) = 0
            jaicgam(ip) = 0
            jaicsig(ip) = 0
            jaicxyp(ip) = 0
            jaicgth(ip) = 0
        enddo
        do iel = 1, nel
            ksyskut(iel) = -999
            ksysqnc(iel) = -999
            ksysgss(iel) = -999
            ksysqtt(iel) = -999
            jsysqnc(iel) = 0
            jaicvwk(iel) = 0
        enddo
        !
        !
        !---- initialize system-row and system-column counters
        ksys = 0
        jsys = 0
        !
        !---- initialize gamma,sigma,xy column counters
        jaicg = 0
        jaics = 0
        jaicz = 0
        !---- initialize wake gamma (point) column counters
        jaicw = 0
        !---- initialize vortex wake row column counters
        jaicv = 0
        !
        !---- go over all elements
        do iel = 1, nel
            !
            !cC------ LVeZR = f if first,last vortex node is not an unknown
            !c        LV1ZR(IEL) = YP(IPFRST(IEL)) .EQ. 0.0
            !c        LV2ZR(IEL) = YP(IPLAST(IEL)) .EQ. 0.0
            !
            if (netype(iel)==0) then
                !-------- normal solid-wall airfoil...
                !
                !--------------------------------------------------
                !-------- set up variable indices
                !
                do jp = ipfrst(iel), iplast(iel)
                    !---------- wall gamma will be unknown
                    jsys = jsys + 1
                    jsysgam(jp) = jsys
                    jpsys(jsys) = jp
                    !---------- wall sigma will be imposed externally (e.g. BL blowing model)
                    jaics = jaics + 1
                    jaicsig(jp) = jaics
                enddo
                !
                if (lbody(iel)) then
                    !--------- declare normal-velocity DOF, and set element which it influences
                    jsys = jsys + 1
                    jsysqnc(iel) = jsys
                    jpsys(jsys) = 0
                    ndof = iel
                else
                    ndof = 0
                endif
                !
                if (lxymov) then
                    do jp = ipfrst(iel), iplast(iel)
                        !----------- wall x,y movement will be imposed externally (e.g. optimization)
                        jaicz = jaicz + 1
                        jaicxyp(jp) = jaicz
                    enddo
                endif
                !
                !--------------------------------------------------
                !-------- set up equation indices
                !
                if (lv1zr(iel)) then
                    !---------- enforce zero-gamma condition at first point
                    ksys = ksys + 1
                    ksysgam(ipfrst(iel)) = ksys
                endif
                !
                !-------- enforce flow tangency at all control points, note Qn DOF if any
                do ic = icfrst(iel), iclast(iel)
                    ksys = ksys + 1
                    !
                    if (ictype(ic)>=0) then
                        !----------- usual control point
                        ksysvnc(ic) = ksys
                        ieldofc(ic) = ndof
                    else
                        !----------- zero-length panel... set gamma-continuity equation
                        ksysgmg(ic) = ksys
                    endif
                enddo
                !
                if (lv2zr(iel)) then
                    !---------- enforce zero-gamma condition at last point
                    ksys = ksys + 1
                    ksysgam(iplast(iel)) = ksys
                endif
                !
                if (.not.(lv1zr(iel) .or. lv2zr(iel))) then
                    !--------- enforce Kutta condition
                    ksys = ksys + 1
                    ksyskut(iel) = ksys
                endif
                !
                if (lbody(iel) .and. .not.(lv1zr(iel) .and. lv2zr(iel)))  &
                        & then
                    !--------- will have normal-velocity DOF... enforce some constraint on it
                    ksys = ksys + 1
                    if (lqnzr(iel)) then
                        !---------- direct constraint
                        ksysqnc(iel) = ksys
                    else
                        !ccC---------- indirect constraint: set zero loading curvature at TE
                        !cc            KSYSGSS(IEL) = KSYS
                        !
                        !---------- indirect constraint: set zero velocity at interior point near TE
                        ksysqtt(iel) = ksys
                    endif
                endif
                !
                !-------- check for vortex wake gamma on body element
                do jp = ipfrst(iel), iplast(iel)
                    if (ip2ir(jp)>0) then
                        !---------- GTH will be imposed or solved for at this point
                        jaicw = jaicw + 1
                        jaicgth(jp) = jaicw
                    endif
                enddo
                !
            elseif (netype(iel)==1) then
                !-------- wake element...
                !---------- wake gamma will be imposed externally (e.g. BL curvature relation)
                do jp = ipfrst(iel), iplast(iel)
                    jaicg = jaicg + 1
                    jaicgam(jp) = jaicg
                    !---------- wake sigma will be imposed externally (e.g. BL blowing model)
                    !c            JAICS = JAICS + 1
                    !c            JAICSIG(JP) = JAICS
                enddo
                !
                if (lxymov) then
                    do jp = ipfrst(iel), iplast(iel)
                        !----------- wake x,y movement will be imposed externally (e.g. optimization)
                        jaicz = jaicz + 1
                        jaicxyp(jp) = jaicz
                    enddo
                endif
                !
            elseif (netype(iel)>=2 .and. netype(iel)<=4) then
                !-------- axis line, ring or point singularities...(NETYPE=2,3,4)
                !-------- sigma, doublet will be imposed externally
                do jp = ipfrst(iel), iplast(iel)
                    jaicg = jaicg + 1
                    jaicgam(jp) = jaicg
                    jaics = jaics + 1
                    jaicsig(jp) = jaics
                enddo
                !
            elseif (netype(iel)==5) then
                !-------- drag area source element NETYPE=5 (source RHS only)
                !-------- sigma will be imposed externally (e.g. drag source)
                do jp = ipfrst(iel), iplast(iel)
                    jaics = jaics + 1
                    jaicsig(jp) = jaics
                enddo
                !
            elseif (netype(iel)==6) then
                !-------- rotor source element NETYPE=6 (source RHS only)
                !-------- sigma will be imposed externally (e.g. drag source)
                do jp = ipfrst(iel), iplast(iel)
                    jaics = jaics + 1
                    jaicsig(jp) = jaics
                enddo
                !
            elseif (netype(iel)==7) then
                !-------- vortex wake element (RHS)
                !-------- wake gamma will be imposed externally
                if (jaicv<nrp) then
                    jaicv = jaicv + 1
                    ir = iel2ir(iel)
                    jaicvwk(ir) = jaicv
                    !c            write(19,*) 'iel,ir,jaicv ',iel,ir,jaicv
                endif
                !-------- set up points with vortex wake gamma
                do jp = ipfrst(iel), iplast(iel)
                    !---------- GTH will be imposed or solved for at this point
                    jaicw = jaicw + 1
                    jaicgth(jp) = jaicw
                enddo
                !
            else
                write (*, *) 'Unknown NETYPE in SYSP ', netype(iel)
                !
            endif
        enddo
        !
        !---- set system size
        nsys = ksys
        !
        if (ksys/=jsys) write (*, *) 'Error in SYSP, KSYS JSYS =', &
                & ksys, jsys
        !
        naicgam = jaicg
        naicsig = jaics
        naicxyp = jaicz
        naicvwk = jaicv
        naicgth = jaicw
        !
        if (ldbg) then
            write (*, *) ' '
            write (*, *) 'SYSP system setup and pointers'
            write (*, *) '  NP      = ', nptot
            write (*, *) '  NC      = ', nctot
            write (*, *) '  NSYS    = ', nsys
            write (*, *) ' '
            write (*, *) '  NAICGAM = ', naicgam
            write (*, *) '  NAICSIG = ', naicsig
            write (*, *) '  NAICXYP = ', naicxyp
            write (*, *) '  NAICVWK = ', naicvwk
            write (*, *) '  NAICGTH = ', naicgth
            write (*, *) ' '
            write (*, *) ' '
            !
            write (lundbg, *) 'SYSP system setup and pointers'
            write (lundbg, *) '  NP      = ', nptot
            write (lundbg, *) '  NC      = ', nctot
            write (lundbg, *) '  NSYS    = ', nsys
            write (lundbg, *) ' '
            write (lundbg, *) '  NAICGAM = ', naicgam
            write (lundbg, *) '  NAICSIG = ', naicsig
            write (lundbg, *) '  NAICXYP = ', naicxyp
            write (lundbg, *) '  NAICVWK = ', naicvwk
            write (lundbg, *) '  NAICGTH = ', naicgth
            write (lundbg, *) ' '
        endif
        !
        lsysp = .true.
        !
    end subroutine sysp
    !*==CLRSYS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SYSP


    subroutine clrsys
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: j, k
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- clear all system arrays
        do k = 1, nsys
            res(k) = 0.
            do j = 0, nsys
                sys(k, j) = 0.
            enddo
            aicqff(k, 1) = 0.
            do j = 0, naicgam
                aicgam(k, j) = 0.
            enddo
            do j = 0, naicsig
                aicsig(k, j) = 0.
            enddo
            do j = 0, naicxyp
                aicxyp(k, 1, j) = 0.
                aicxyp(k, 2, j) = 0.
            enddo
            do j = 0, naicvwk
                aicvwk(k, j) = 0.
            enddo
        enddo
        !
    end subroutine clrsys
    !*==GSYS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CLRSYS


    subroutine gsys(lsys, lqff, lvwk, ip1, ip2, iz1, iz2)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ip1, ip2, iz1, iz2
        logical :: lqff, lsys, lvwk
        !
        ! Local variables
        !
        real :: dn1sw, dn2sw, ds1m, ds1o, ds1sw, ds2m, ds2o, &
                & ds2sw, f1l, f1m, f1o, f2l, f2m, f2o, qttf, &
                & res_dx, res_dy
        integer :: ic, ie, iel, ielw, ipo, ipp, ir, j, j1, &
                & j1l, j1m, j1o, j2, j2l, j2m, j2o, jp, jp1, &
                & jp1l, jp1m, jp1o, jp2, jp2l, jp2m, jp2o, jpo, &
                & jpp, jpw, jw, k
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
        if (nsys==0) return
        !
        if (ldbg) then
            write (*, *) 'Setting up direct-problem Newton system...'
            write (lundbg, *) 'GSYS direct problem system setup'
            write (lundbg, *) '   LSYS,LQFF,LVWK ', lsys, lqff, lvwk
            write (lundbg, *) '   IP1,IP2        ', ip1, ip2
            write (lundbg, *) '   IZ1,IZ2        ', iz1, iz2
        endif
        !
        !----------------------------------------------------------
        !---- clear system arrays
        if (lsys) then
            do k = 1, nsys
                res(k) = 0.
                do j = 1, nsys
                    sys(k, j) = 0.
                enddo
            enddo
        endif
        !
        !---- clear RHS arrays
        if (lqff) then
            do k = 1, nsys
                aicqff(k, 1) = 0.
            enddo
        endif
        !
        if (lvwk) then
            do ir = 1, nrp
                j = jaicvwk(ir)
                do k = 1, nsys
                    aicvwk(k, j) = 0.
                enddo
            enddo
        endif
        !
        do k = 1, nsys
            do jp = ip1, ip2
                j = jaicgam(jp)
                aicgam(k, j) = 0.
                !
                j = jaicsig(jp)
                aicsig(k, j) = 0.
                !
                j = jaicgth(jp)
                aicgth(k, j) = 0.
            enddo
            !
            do jp = iz1, iz2
                j = jaicxyp(jp)
                aicxyp(k, 1, j) = 0.
                aicxyp(k, 2, j) = 0.
            enddo
        enddo
        !
        !----------------------------------------------------------
        !----------------------------------------------------------
        !---- set flow tangency at all control points
        do ic = 1, nctot
            !
            !        write(77,*) 'IC ',IC
            !        do jp = 1, np
            !          write(77,*) 'JP, QC1-2 ',JP,QC_GTH(1,JP,IC),QC_GTH(2,JP,IC)
            !        end do
            !------ K = system row index for control point IC
            k = ksysvnc(ic)
            !
            !------ skip un-enforced control points
            if (k>=0) then
                !
                !----------------------------------------------------------
                if (lsys) then
                    !------- flow-tangency residual
                    res(k) = anc(1, ic) * qc(1, ic) + anc(2, ic) * qc(2, ic)
                    !
                    !
                    do jp = 1, nptot
                        !--------- dRES/dgamma_surface
                        !-         J = system column index of GAM(JP)  (0 if GAM is on r.h.s.)
                        j = jsysgam(jp)
                        sys(k, j) = anc(1, ic) * qc_gam(1, jp, ic) + anc(2, ic)      &
                                & * qc_gam(2, jp, ic)
                    enddo
                endif
                !
                !----------------------------------------------------------
                !------- dRES/dQinf
                if (lqff) aicqff(k, 1) = anc(1, ic)
                !
                !----------------------------------------------------------
                if (lvwk) then
                    !c          write(19,*) 'ic ',ic
                    !------- dRES/dBCirc rotor vortex wake influence
                    do jp = 1, nptot
                        ir = ip2ir(jp)
                        if (ir/=0) then
                            !---- we are on a vortex wake element
                            !
                            !---- accumulate wake influences by wake rows using wake index
                            j = jaicvwk(ir)
                            aicvwk(k, j) = aicvwk(k, j)                          &
                                    & + (anc(1, ic) * qc_gth(1, jp, ic)         &
                                            & + anc(2, ic) * qc_gth(2, jp, ic)) / yp(jp)
                            !c              write(19,*) 'jp,ir,j ',jp,ir,j
                            !
                            !---- set wake point influence individually
                            j = jaicgth(jp)
                            aicgth(k, j) = anc(1, ic) * qc_gth(1, jp, ic) + anc(2, ic)&
                                    & * qc_gth(2, jp, ic)
                        endif
                    enddo
                endif
                !
                !----------------------------------------------------------
                do jp = ip1, ip2
                    !-------- dRES/dgamma, dRES/dsigma
                    !-        J = r.h.s. column index  (0 if GAM,SIG is a l.h.s. unknown)
                    j = jaicgam(jp)
                    aicgam(k, j) = anc(1, ic) * qc_gam(1, jp, ic) + anc(2, ic)      &
                            & * qc_gam(2, jp, ic)
                    !
                    j = jaicsig(jp)
                    aicsig(k, j) = anc(1, ic) * qc_sig(1, jp, ic) + anc(2, ic)      &
                            & * qc_sig(2, jp, ic)
                enddo
                !
                !----------------------------------------------------------
                do jp = iz1, iz2
                    !-------- dRES/dxy
                    !-        J = r.h.s. column index
                    j = jaicxyp(jp)
                    aicxyp(k, 1, j) = anc(1, ic) * qc_xp(1, jp, ic) + anc(2, ic)     &
                            & * qc_xp(2, jp, ic)
                    aicxyp(k, 2, j) = anc(1, ic) * qc_yp(1, jp, ic) + anc(2, ic)     &
                            & * qc_yp(2, jp, ic)
                enddo
                !
                !----------------------------------------------------------
                !------ local source contribution to flow-tangency residual
                ipo = ipco(ic)
                ipp = ipcp(ic)
                !
                if (lsys) res(k) = res(k) - 0.25 * (sig(ipo) + sig(ipp))
                !
                if (ipo>=ip1 .and. ipo<=ip2) then
                    j = jaicsig(ipo)
                    aicsig(k, j) = aicsig(k, j) - 0.25
                endif
                if (ipp>=ip1 .and. ipp<=ip2) then
                    j = jaicsig(ipp)
                    aicsig(k, j) = aicsig(k, j) - 0.25
                endif
                !
                !------ derivatives w.r.t. dn movement of the two nodes of panel IC
                res_dx = anc_dxy(1, 1, ic) * qc(1, ic) + anc_dxy(2, 1, ic) * qc(2, ic)
                res_dy = anc_dxy(1, 2, ic) * qc(1, ic) + anc_dxy(2, 2, ic) * qc(2, ic)
                !
                if (ipo>=iz1 .and. ipo<=iz2) then
                    j = jaicxyp(ipo)
                    aicxyp(k, 1, j) = aicxyp(k, 1, j) - res_dx
                    aicxyp(k, 2, j) = aicxyp(k, 2, j) - res_dy
                endif
                if (ipp>=iz1 .and. ipp<=iz2) then
                    j = jaicxyp(ipp)
                    aicxyp(k, 1, j) = aicxyp(k, 1, j) + res_dx
                    aicxyp(k, 2, j) = aicxyp(k, 2, j) + res_dy
                endif
                !
                !----------------------------------------------------------
                if (lsys .and. ieldofc(ic)>0) then
                    !------- this control point has a normal-velocity QNDOF contributing to Q.n
                    ie = ieldofc(ic)
                    res(k) = res(k) - qndof(ie)
                    !
                    !------- J = system column index of QNDOF(IE)
                    j = jsysqnc(ie)
                    if (j<=0) then
                        write (*, *) '? Unassigned normal-velocity DOF', ie
                        write (lundbg, *) '? Unassigned normal-velocity DOF', &
                                & ie
                    endif
                    !
                    sys(k, j) = -1.0
                endif
                if (ldbg) write (lundbg, *) k, 'V.n=Qn ', res(k)
                !@@@
            endif
        enddo
        !
        !---- set gamma-continuity at all zero-length panels
        do ic = 1, nctot
            !------ K = system row index of gamma-continuity equation (if any)
            k = ksysgmg(ic)
            !
            if (lsys .and. k>0) then
                jpo = ipco(ic)
                jpp = ipcp(ic)
                res(k) = gam(jpp) - gam(jpo)
                !
                sys(k, jpo) = -1.0
                sys(k, jpp) = 1.0
                if (ldbg) write (lundbg, *) k, 'dGam=0 ', res(k)
                !@@@
            endif
        enddo
        !
        !
        do iel = 1, nel
            !------ zero net curvature of gamma at TE (triangular loading towards TE)
            k = ksysgss(iel)
            if (lsys .and. k>0) then

                jp1o = ipfrst(iel) + 2
                jp1m = ipfrst(iel) + 1
                jp1l = ipfrst(iel)
                jp2o = iplast(iel)
                jp2m = iplast(iel) - 1
                jp2l = iplast(iel) - 2
                !
                ds1o = sp(jp1o) - sp(jp1m)
                ds1m = sp(jp1m) - sp(jp1l)
                ds2o = sp(jp2o) - sp(jp2m)
                ds2m = sp(jp2m) - sp(jp2l)
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
                f1o = 1.0 / ds1o
                f1m = -(1.0 / ds1o + 1.0 / ds1m)
                f1l = 1.0 / ds1m
                f2o = 1.0 / ds2o
                f2m = -(1.0 / ds2o + 1.0 / ds2m)
                f2l = 1.0 / ds2m
                !
                j1o = jsysgam(jp1o)
                j1m = jsysgam(jp1m)
                j1l = jsysgam(jp1l)
                j2o = jsysgam(jp2o)
                j2m = jsysgam(jp2m)
                j2l = jsysgam(jp2l)
                !
                if (lv1zr(iel)) then
                    res(k) = -f2o * gam(jp2o) - f2m * gam(jp2m) - f2l * gam(jp2l)
                    sys(k, j2o) = -f2o
                    sys(k, j2m) = -f2m
                    sys(k, j2l) = -f2l
                    !
                elseif (lv2zr(iel)) then
                    res(k) = f1o * gam(jp1o) + f1m * gam(jp1m) + f1l * gam(jp1l)
                    sys(k, j1o) = f1o
                    sys(k, j1m) = f1m
                    sys(k, j1l) = f1l
                    !
                else
                    res(k) = f1o * gam(jp1o) + f1m * gam(jp1m) + f1l * gam(jp1l)   &
                            & - f2o * gam(jp2o) - f2m * gam(jp2m) - f2l * gam(jp2l)
                    sys(k, j1o) = f1o
                    sys(k, j1m) = f1m
                    sys(k, j1l) = f1l
                    sys(k, j2o) = -f2o
                    sys(k, j2m) = -f2m
                    sys(k, j2l) = -f2l
                    !
                endif
                if (ldbg) write (lundbg, *) k, 'Qss=0 ', res(k)
                !@@@
            endif
            !
            !------ Kutta condition
            k = ksyskut(iel)
            if (k>0) then
                if (lsys) then
                    !
                    if (lbody(iel)) then
                        !------------ closed body
                        call swdots(iel, ds1sw, ds2sw, dn1sw, dn2sw)
                        jp1 = ipfrst(iel)
                        jp2 = iplast(iel)
                        !
                        !@@@
                        ds1sw = -1.0
                        ds2sw = 1.0
                        dn1sw = 0.
                        dn2sw = 0.
                        !------------ basic Kutta condition
                        res(k) = ds2sw * gam(jp2) - ds1sw * gam(jp1)              &
                                & + dn2sw * sig(jp2) - dn1sw * sig(jp1)
                        !cc     &               + CIRDOT*2.0/(GAM(JP1)-GAM(JP2))
                        !
                        j1 = jsysgam(jp1)
                        j2 = jsysgam(jp2)
                        sys(k, j1) = -ds1sw
                        sys(k, j2) = ds2sw
                        !
                        if (jp1>=ip1 .and. jp1<=ip2) then
                            j1 = jaicsig(jp1)
                            aicsig(k, j1) = -dn1sw
                        endif
                        if (jp2>=ip1 .and. jp2<=ip2) then
                            !cc               J2 = JAICGAM(JP2)   !!! bug
                            j2 = jaicsig(jp2)
                            aicsig(k, j2) = dn2sw
                        endif
                        !
                        if (iewake(iel)/=0) then
                            !------------- add wake panel contribution
                            ielw = iewake(iel)
                            jpw = ipfrst(ielw)
                            res(k) = res(k) - gam(jpw)
                            res(k) = res(k) - gth(jpw)      !!HHY
                            if (jpw>=ip1 .and. jp2<=ip2) then
                                jw = jaicgam(jpw)
                                aicgam(k, jw) = -1.0
                                !                JW = JAICGTH(JPW)
                                !                AICGTH(K,JW) = -1.0
                            endif
                        endif
                        !
                    else
                        !------------ zero-thickness plate
                        jp = iplast(iel)
                        !
                        res(k) = gam(jp)
                        !
                        j = jsysgam(jp)
                        sys(k, j) = 1.0
                        !
                        if (iewake(iel)/=0) then
                            !------------- add wake panel contribution
                            ielw = iewake(iel)
                            jpw = ipfrst(ielw)
                            res(k) = res(k) - gam(jpw)
                            res(k) = res(k) - gth(jpw)      !!HHY
                            !
                            if (jpw>=ip1 .and. jp2<=ip2) then
                                jw = jaicgam(jpw)
                                aicgam(k, jw) = -1.0
                                !                JW = JAICGTH(JPW)
                                !                AICGTH(K,JW) = -1.0
                            endif
                        endif
                        !
                    endif
                endif
                if (ldbg) write (lundbg, *) k, 'Kutta ', res(k)
                !@@@
            endif
            !
            !------ normal-velocity QNDOF constraint
            k = ksysqnc(iel)
            if (lsys .and. k>0) then

                j = jsysqnc(iel)
                res(k) = qndof(iel)
                !CC  &          - QNDOFSP(IEL)
                sys(k, j) = 1.0
                if (ldbg) write (lundbg, *) k, 'Qn=0 ', res(k)
                !@@@
            endif
            !
            !------ zero internal tangential velocity constraint
            k = ksysqtt(iel)
            if (k>0) then
                !
                ic = ncx + iel
                !
                if (lsys) then

                    !cc debug stuff to fiddle with net massflow from QNdof
                    !cc            write(*,*) 'enter Qtt'   !@@@
                    !cc            read(*,*) qttf    !@@@
                    qttf = 0.0

                    !-------- flow-tangency residual
                    res(k) = anc(1, ic) * qc(1, ic) + anc(2, ic) * qc(2, ic)         &
                            & - qttf * qinf                              !@@@
                    !
                    if (ldbg) write (lundbg, *) k, 'Qtt=0 ', res(k)
                    !@@@

                    do jp = 1, nptot
                        !---------- dRES/dgamma_surface
                        !-          J = system column index of GAM(JP)  (0 if GAM is on r.h.s.)
                        j = jsysgam(jp)
                        sys(k, j) = anc(1, ic) * qc_gam(1, jp, ic) + anc(2, ic)      &
                                & * qc_gam(2, jp, ic)
                    enddo
                endif
                !
                !-------- dRES/dQinf
                if (lqff) aicqff(k, 1) = anc(1, ic) - qttf
                !@@@
                !
                !----------------------------------------------------------
                if (lvwk) then
                    !c          write(19,*) 'icx ',ic
                    !------- dRES/dBCirc rotor vortex wake influence
                    do jp = 1, nptot
                        ir = ip2ir(jp)
                        if (ir/=0) then
                            !---- we are on a vortex wake element
                            !
                            !---- accumulate wake influences by wake rows using wake index
                            j = jaicvwk(ir)
                            aicvwk(k, j) = aicvwk(k, j)                          &
                                    & + (anc(1, ic) * qc_gth(1, jp, ic)         &
                                            & + anc(2, ic) * qc_gth(2, jp, ic)) / yp(jp)
                            !c              write(19,*) 'jp,ir,j ',jp,ir,j
                            !
                            !---- set wake point influence individually
                            j = jaicgth(jp)
                            aicgth(k, j) = anc(1, ic) * qc_gth(1, jp, ic) + anc(2, ic)&
                                    & * qc_gth(2, jp, ic)
                        endif
                    enddo
                endif
                !
                !----------------------------------------------------------
                do jp = ip1, ip2
                    !--------- dRES/dgamma, dRES/dsigma
                    !-         J = r.h.s. column index  (0 if GAM,SIG is a l.h.s. unknown)
                    j = jaicgam(jp)
                    aicgam(k, j) = anc(1, ic) * qc_gam(1, jp, ic) + anc(2, ic)      &
                            & * qc_gam(2, jp, ic)
                    !
                    j = jaicsig(jp)
                    aicsig(k, j) = anc(1, ic) * qc_sig(1, jp, ic) + anc(2, ic)      &
                            & * qc_sig(2, jp, ic)
                enddo
            endif
            !
        enddo
        !
        !---- set zero vortex strengths explicitly (if selected)
        if (lsys) then
            do jp = ip1, ip2
                k = ksysgam(jp)
                if (k>0) then
                    j = jsysgam(jp)
                    res(k) = gam(jp)
                    sys(k, j) = 1.0
                    if (ldbg) write (lundbg, *) k, 'gam=0 ', res(k)
                    !@@@
                endif
            enddo
        endif
        !
        if (lsys) then
            lgsys = .true.
            lgamu = .false.
            lqgic = .false.
        endif
        !
    end subroutine gsys
    !*==SYSPQ.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GSYS



    subroutine syspq
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: iel, ip, ips1, ips2, iseg, jsys, k, ksys
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------
        !     Modifies pointers from SYSP to solve
        !     mixed-inverse problem
        !-------------------------------------------
        !
        if (ldbg) then
            write (*, *) 'Modifying system pointers for mixed-inverse...'
            write (lundbg, *)                                               &
                    &'SYSPQ modify system pointers for mixed inverse'
            write (lundbg, *) '   NSYS =', nsys, ' (before augmentation)'
        endif
        !
        !---- -999 indicates un-enforced control point
        !----    0 indicates undefined variable
        do iseg = 0, nseg
            do k = 1, 4
                ksysdns(k, iseg) = -999
                jsysqsp(k, iseg) = 0
            enddo
        enddo
        !
        !---- assume that no panel node lies inside an inverse segment
        do ip = 1, nptot
            ksegp(ip) = 0
        enddo
        !
        !---- initialize system-row and system-column counters to current count
        ksys = nsys
        jsys = nsys
        !
        do iseg = 1, nseg
            ips1 = ipseg1(iseg)
            ips2 = ipseg2(iseg)
            !
            !------ for nodes inside target segment, change variable from dgam to dn
            do ip = ips1, ips2
                jsysdxy(ip) = jsysgam(ip)
                jsysgam(ip) = 0
                !
                ksegp(ip) = iseg
            enddo
            !
            iel = ielseg(iseg)
            !
            if (lnfix1(iseg) .or. lnfix2(iseg)) then
                !------- use constant term if at least one endpoint is fixed
                jsys = jsys + 1
                jsysqsp(1, iseg) = jsys
                !
                ksys = ksys + 1
                if (lnfix1(iseg)) then
                    ksysdns(1, iseg) = ksys
                else
                    ksysdns(2, iseg) = ksys
                endif
            endif
            !
            if (lnfix1(iseg) .and. lnfix2(iseg)) then
                !------- also use linear term if the other endpoint is also fixed
                jsys = jsys + 1
                jsysqsp(2, iseg) = jsys
                !
                ksys = ksys + 1
                if (lnfix2(iseg)) then
                    ksysdns(2, iseg) = ksys
                else
                    ksysdns(1, iseg) = ksys
                endif
            endif
            !
            if (ldfix1(iel)) then
                !------- match slope on the left
                ksys = ksys + 1
                jsys = jsys + 1
                ksysdns(3, iseg) = ksys
                jsysqsp(3, iseg) = jsys
            endif
            !
            if (ldfix2(iel)) then
                !------- match slope on the right
                ksys = ksys + 1
                jsys = jsys + 1
                ksysdns(4, iseg) = ksys
                jsysqsp(4, iseg) = jsys
            endif
            !
        enddo
        !
        !---- new augmented-system size
        nsys = ksys
        if (ldbg) write (lundbg, *) '   NSYS =', nsys, &
                &' augmented for inverse'
        !
        if (ksys/=jsys) write (*, *) 'Error in SYSPQ, KSYS JSYS =', &
                & ksys, jsys
        !
        !---- standard direct-problem pointers no longer valid
        lsysp = .false.
        !
    end subroutine syspq
    !*==QSYS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SYSPQ


    subroutine qsys
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: ds1m, ds1o, ds2m, ds2o, dsm, dsp, f1l, f1m, &
                & f1o, f2l, f2m, f2o, frm, fro, frp, qttf, &
                & res_anc1, res_anc2, res_dx, res_dy, res_gam, &
                & res_gamo, res_gamp, res_qc1, res_qc2, res_sig, &
                & res_sigo, res_sigp, res_uinf, res_xp, res_xpo, &
                & res_xpp, res_yp, res_ypo, res_ypp
        integer :: ic, iel, ielw, iseg, j, jo, jp, jp1, jp1l, &
                & jp1m, jp1o, jp2, jp2l, jp2m, jp2o, jpm, jpo, &
                & jpp, k, kqsp
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
        if (nsys==0) return
        !
        if (ldbg) then
            write (*, *) 'Setting up mixed-inverse Newton system...'
            write (lundbg, *) 'QSYS mixed-inverse problem system setup'
            write (lundbg, *) '   NSYS =', nsys
        endif
        !
        !---- clear system arrays
        do k = 1, nsys
            res(k) = 0.
            do j = 0, nsys
                sys(k, j) = 0.
            enddo
            do j = 0, naicgam
                aicgam(k, j) = 0.
            enddo
            do j = 0, naicsig
                aicsig(k, j) = 0.
            enddo
            aicqff(k, 1) = 0.
        enddo
        !
        !---- set flow tangency at all control points
        do ic = 1, nctot
            k = ksysvnc(ic)
            !
            !------ skip un-enforced control points
            if (k>=0) then
                !
                !------- flow tangency at control point IC
                res(k) = anc(1, ic) * qc(1, ic) + anc(2, ic) * qc(2, ic)
                !
                res_qc1 = anc(1, ic)
                res_qc2 = anc(2, ic)
                res_anc1 = qc(1, ic)
                res_anc2 = qc(2, ic)
                res_uinf = anc(1, ic)
                !
                !------- derivatives w.r.t. dn movement of the two nodes of panel IC
                jpo = ipco(ic)
                jpp = ipcp(ic)
                res_dx = res_anc1 * anc_dxy(1, 1, ic) + res_anc2 * anc_dxy(2, 1, ic)
                res_dy = res_anc1 * anc_dxy(1, 2, ic) + res_anc2 * anc_dxy(2, 2, ic)
                res_xpo = -res_dx
                res_ypo = -res_dy
                res_xpp = res_dx
                res_ypp = res_dy
                !
                jo = jsysdxy(jpo)
                jp = jsysdxy(jpp)
                sys(k, jo) = sys(k, jo) + res_xpo * anp(1, jpo)                  &
                        & + res_ypo * anp(2, jpo)
                sys(k, jp) = sys(k, jp) + res_xpp * anp(1, jpp)                  &
                        & + res_ypp * anp(2, jpp)
                !
                !------- derivatives w.r.t. vorticity and position of all nodes
                do jp = 1, nptot
                    res_gam = res_qc1 * qc_gam(1, jp, ic)                        &
                            & + res_qc2 * qc_gam(2, jp, ic)
                    res_sig = res_qc1 * qc_sig(1, jp, ic)                        &
                            & + res_qc2 * qc_sig(2, jp, ic)
                    !
                    res_xp = res_qc1 * qc_xp(1, jp, ic) + res_qc2 * qc_xp(2, jp, ic)
                    res_yp = res_qc1 * qc_yp(1, jp, ic) + res_qc2 * qc_yp(2, jp, ic)
                    !
                    call gxylin(k, jp, res_gam, res_sig, res_xp, res_yp)
                enddo
                !
                if (ieldofc(ic)>0) then
                    !-------- flow tangency contribution of normal-velocity DOF
                    iel = ieldofc(ic)
                    res(k) = res(k) - qndof(iel)
                    !
                    j = jsysqnc(iel)
                    sys(k, j) = -1.0
                endif
                if (ldbg) write (lundbg, *) k, 'V.n=Qn ', res(k)
                !@@@
                !
                aicqff(k, 1) = res_uinf
            endif
        enddo
        !
        !---- set gamma or dn-continuity at all zero-length panels
        do ic = 1, nctot
            k = ksysgmg(ic)
            !
            !------ skip if there's no constraint
            if (k>=0) then
                !
                jpo = ipco(ic)
                jpp = ipcp(ic)
                !
                if (jsysgam(jpo)>0 .and. jsysgam(jpp)>0) then
                    res(k) = gam(jpp) - gam(jpo)
                    res_gamp = 1.0
                    res_gamo = -1.0
                    res_sigp = 0.
                    res_sigo = 0.
                    res_xpp = 0.
                    res_xpo = 0.
                    res_ypp = 0.
                    res_ypo = 0.
                else
                    !CC      RES(K) = DN(JPP) - DN(JPO)
                    res(k) = 0.
                    res_gamp = 0.
                    res_gamo = 0.
                    res_xpp = anp(1, jpp)
                    res_xpo = anp(1, jpo)
                    res_ypp = anp(2, jpp)
                    res_ypo = anp(2, jpo)
                endif
                if (ldbg) write (lundbg, *) k, 'dGam=0 ', res(k)
                !@@@
                !
                call gxylin(k, jpo, res_gamo, res_sigo, res_xpo, res_ypo)
                call gxylin(k, jpp, res_gamp, res_sigp, res_xpp, res_ypp)
            endif
        enddo
        !
        do iel = 1, nel
            !------ zero net curvature of gamma at TE (triangular loading towards TE)
            k = ksysgss(iel)
            if (k>0) then

                jp1o = ipfrst(iel) + 2
                jp1m = ipfrst(iel) + 1
                jp1l = ipfrst(iel)
                jp2o = iplast(iel)
                jp2m = iplast(iel) - 1
                jp2l = iplast(iel) - 2
                !
                ds1o = sp(jp1o) - sp(jp1m)
                ds1m = sp(jp1m) - sp(jp1l)
                ds2o = sp(jp2o) - sp(jp2m)
                ds2m = sp(jp2m) - sp(jp2l)
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
                f1o = 1.0 / ds1o
                f1m = -(1.0 / ds1o + 1.0 / ds1m)
                f1l = 1.0 / ds1m
                f2o = 1.0 / ds2o
                f2m = -(1.0 / ds2o + 1.0 / ds2m)
                f2l = 1.0 / ds2m
                !
                if (lv1zr(iel)) then
                    res(k) = -f2o * gam(jp2o) - f2m * gam(jp2m) - f2l * gam(jp2l)
                    call gxylin(k, jp2o, -f2o, 0.0, 0.0, 0.0)
                    call gxylin(k, jp2m, -f2m, 0.0, 0.0, 0.0)
                    call gxylin(k, jp2l, -f2l, 0.0, 0.0, 0.0)
                    !
                elseif (lv2zr(iel)) then
                    res(k) = f1o * gam(jp1o) + f1m * gam(jp1m) + f1l * gam(jp1l)
                    call gxylin(k, jp1o, f1o, 0.0, 0.0, 0.0)
                    call gxylin(k, jp1m, f1m, 0.0, 0.0, 0.0)
                    call gxylin(k, jp1l, f1l, 0.0, 0.0, 0.0)
                    !
                else
                    res(k) = f1o * gam(jp1o) + f1m * gam(jp1m) + f1l * gam(jp1l)   &
                            & - f2o * gam(jp2o) - f2m * gam(jp2m) - f2l * gam(jp2l)
                    call gxylin(k, jp1o, f1o, 0.0, 0.0, 0.0)
                    call gxylin(k, jp1m, f1m, 0.0, 0.0, 0.0)
                    call gxylin(k, jp1l, f1l, 0.0, 0.0, 0.0)
                    call gxylin(k, jp2o, -f2o, 0.0, 0.0, 0.0)
                    call gxylin(k, jp2m, -f2m, 0.0, 0.0, 0.0)
                    call gxylin(k, jp2l, -f2l, 0.0, 0.0, 0.0)
                    !
                endif
                if (ldbg) write (lundbg, *) k, 'Qss=0 ', res(k)
                !@@@
            endif
            !
            !------ Kutta condition
            k = ksyskut(iel)
            if (k>0) then
                !
                if (lbody(iel)) then
                    jp1 = ipfrst(iel)
                    jp2 = iplast(iel)
                    res(k) = gam(jp1) + gam(jp2)
                    call gxylin(k, jp1, 1.0, 0.0, 0.0, 0.0)
                    call gxylin(k, jp2, 1.0, 0.0, 0.0, 0.0)
                    !
                else
                    jp = iplast(iel)
                    res(k) = gam(jp)
                    call gxylin(k, jp, 1.0, 0.0, 0.0, 0.0)
                    !
                endif
                !
                if (iewake(iel)/=0) then
                    ielw = iewake(iel)
                    jp = ipfrst(ielw)
                    res(k) = res(k) - gam(jp)
                    !           RES(K) = RES(K) - GTH(JP)  !!! HHY
                    call gxylin(k, jp, -1.0, 0.0, 0.0, 0.0)
                endif
                if (ldbg) write (lundbg, *) k, 'Kutta ', res(k)
                !@@@
            endif
            !
            !------ normal-velocity DOF
            k = ksysqnc(iel)
            if (k>0) then
                !
                j = jsysqnc(iel)
                res(k) = qndof(iel)
                !CC  &          - QNDOFSP(IEL)
                sys(k, j) = 1.0
                if (ldbg) write (lundbg, *) k, 'Qn=0 ', res(k)
                !@@@
            endif
            !
            !------ zero internal tangential velocity constraint
            k = ksysqtt(iel)
            if (k>0) then
                !
                ic = ncx + iel
                !
                !            write(*,*) 'enter Qtan/Qinf'   !@@@
                !            read(*,*) qttf    !@@@
                qttf = 0.0

                !-------- flow-tangency residual
                res(k) = anc(1, ic) * qc(1, ic) + anc(2, ic) * qc(2, ic) - qttf * qinf
                !@@@
                !
                if (ldbg) write (lundbg, *) k, 'Qtan=0 ', res(k)
                !@@@
                res_qc1 = anc(1, ic)
                res_qc2 = anc(2, ic)
                res_anc1 = qc(1, ic)
                res_anc2 = qc(2, ic)
                res_uinf = anc(1, ic)
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
                do jp = 1, nptot
                    res_gam = res_qc1 * qc_gam(1, jp, ic)                        &
                            & + res_qc2 * qc_gam(2, jp, ic)
                    res_sig = res_qc1 * qc_sig(1, jp, ic)                        &
                            & + res_qc2 * qc_sig(2, jp, ic)
                    !
                    res_xp = res_qc1 * qc_xp(1, jp, ic) + res_qc2 * qc_xp(2, jp, ic)
                    res_yp = res_qc1 * qc_yp(1, jp, ic) + res_qc2 * qc_yp(2, jp, ic)
                    !
                    call gxylin(k, jp, res_gam, res_sig, res_xp, res_yp)
                enddo
                !
            endif
        enddo
        !
        !
        !---- set zero vortex strengths explicitly
        do jp = 1, nptot
            k = ksysgam(jp)
            !
            !------ skip non-specified vortices
            if (k>=0) then
                !
                res(k) = gam(jp)
                call gxylin(k, jp, 1.0, 0.0, 0.0, 0.0)
                if (ldbg) write (lundbg, *) k, 'gam=0 ', res(k)
                !@@@
            endif
        enddo
        !
        !
        !---- segment-fixing conditions
        do iseg = 1, nseg
            k = ksysdns(1, iseg)
            if (k>0) then

                jp = ipseg1(iseg)
                !CC      RES(K) = DN(JP)
                res(k) = 0.
                j = jsysdxy(jp)
                sys(k, j) = 1.0
                if (ldbg) write (lundbg, *) k, 'dn1=0 ', res(k)
                !@@@
            endif
            !
            k = ksysdns(2, iseg)
            if (k>0) then

                jp = ipseg2(iseg)
                !CC      RES(K) = DN(JP)
                res(k) = 0.
                j = jsysdxy(jp)
                sys(k, j) = 1.0
                if (ldbg) write (lundbg, *) k, 'dn2=0 ', res(k)
                !@@@
            endif
            !
            k = ksysdns(3, iseg)
            if (k>0) then

                !
                !         JP = IPSEG1(ISEG) + 1
                !CCC      RES(K) = DN(JP)
                !         RES(K) = 0.
                !         J = JSYSDXY(JP)
                !         SYS(K,J) = 1.0
                !
                kqsp = kqtarg
                jpm = ipseg1(iseg) - 1
                jpo = ipseg1(iseg)
                jpp = ipseg1(iseg) + 1
                dsm = sspec(jpo) - sspec(jpm)
                dsp = sspec(jpp) - sspec(jpo)
                frm = 1.0 / dsm
                fro = -(1.0 / dsm + 1.0 / dsp)
                frp = 1.0 / dsp
                res(k) = frm * gam(jpm) + fro * gam(jpo) + frp * gam(jpp)         &
                        & - frm * qspec(jpm, kqsp) - fro * qspec(jpo, kqsp)        &
                        & - frp * qspec(jpp, kqsp)
                call gxylin(k, jpm, frm, 0.0, 0.0, 0.0)
                call gxylin(k, jpo, fro, 0.0, 0.0, 0.0)
                call gxylin(k, jpp, frp, 0.0, 0.0, 0.0)
                if (ldbg) write (lundbg, *) k, 'gam_ss1=0 ', res(k)
                !@@@
            endif
            !
            k = ksysdns(4, iseg)
            if (k>0) then
                !
                !
                !         JP = IPSEG2(ISEG) - 1
                !CCC      RES(K) = DN(JP)
                !         RES(K) = 0.
                !         J = JSYSDXY(JP)
                !         SYS(K,J) = 1.0
                !
                kqsp = kqtarg
                jpm = ipseg2(iseg) - 1
                jpo = ipseg2(iseg)
                jpp = ipseg2(iseg) + 1
                dsm = sspec(jpo) - sspec(jpm)
                dsp = sspec(jpp) - sspec(jpo)
                frm = 1.0 / dsm
                fro = -(1.0 / dsm + 1.0 / dsp)
                frp = 1.0 / dsp
                res(k) = frm * gam(jpm) + fro * gam(jpo) + frp * gam(jpp)         &
                        & - frm * qspec(jpm, kqsp) - fro * qspec(jpo, kqsp)        &
                        & - frp * qspec(jpp, kqsp)
                call gxylin(k, jpm, frm, 0.0, 0.0, 0.0)
                call gxylin(k, jpo, fro, 0.0, 0.0, 0.0)
                call gxylin(k, jpp, frp, 0.0, 0.0, 0.0)
                if (ldbg) write (lundbg, *) k, 'gam_ss2=0 ', res(k)
                !@@@
            endif
        enddo
        !
        lqsys = .true.
        !
    end subroutine qsys
    !*==GXYLIN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QSYS




    subroutine gxylin(k, jp, res_gam, res_sig, res_xp, res_yp)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: jp, k
        real :: res_gam, res_sig, res_xp, res_yp
        !
        ! Local variables
        !
        integer :: iseg, j, j1, j2, j3, j4
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        j = jsysgam(jp)
        sys(k, j) = sys(k, j) + res_gam
        !
        j = jaicgam(jp)
        aicgam(k, j) = aicgam(k, j) + res_gam
        !
        j = jaicsig(jp)
        aicsig(k, j) = aicsig(k, j) + res_sig
        !
        iseg = ksegp(jp)
        j1 = jsysqsp(1, iseg)
        j2 = jsysqsp(2, iseg)
        j3 = jsysqsp(3, iseg)
        j4 = jsysqsp(4, iseg)
        sys(k, j1) = sys(k, j1) + res_gam * fspec(1, jp)
        sys(k, j2) = sys(k, j2) + res_gam * fspec(2, jp)
        sys(k, j3) = sys(k, j3) + res_gam * fspec(3, jp)
        sys(k, j4) = sys(k, j4) + res_gam * fspec(4, jp)
        !
        j = jsysdxy(jp)
        sys(k, j) = sys(k, j) + res_xp * anp(1, jp) + res_yp * anp(2, jp)
        !
    end subroutine gxylin
    !*==GUCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GXYLIN



    subroutine gucalc(lguqi, lvwk, ip1, ip2)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ip1, ip2
        logical :: lguqi, lvwk
        !
        ! Local variables
        !
        integer :: ip, ir, j
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
        if (ldbg) write (*, *)                                           &
                &'Generating unit vorticity distributions...'
        !
        !---------------------------------------------------------------------
        !---- assign unit-solution indices
        if (lguqi) then
            if (iuqinf==0) then
                !------ freestream x-velocity
                if (ldbg) write (*, *) 'Generating solution for unit Uinf'
                call gumake(iuqinf, aicqff(1, 1))
            endif
        endif
        !
        if (lvwk) then
            do ir = 1, nrp
                if (iuvwk(ir)==0) then
                    !------ rotor circulation, radial station IR
                    if (ldbg) write (*, *)                                  &
                            &'Generating solution for unit GAM IR '&
                            &, ir
                    j = jaicvwk(ir)
                    call gumake(iuvwk(ir), aicvwk(1, j))
                    call guvwk(ir)
                endif
            enddo
        endif
        !
        do ip = ip1, ip2
            j = jaicgam(ip)
            if (iugam(ip)==0 .and. j>0) then
                if (ldbg) write (*, *) 'Generating solution for unit gamma'&
                        &, ip
                call gumake(iugam(ip), aicgam(1, j))
                gamu(ip, iugam(ip)) = 1.0
            endif
            !
            j = jaicsig(ip)
            if (iusig(ip)==0 .and. j>0) then
                if (ldbg) write (*, *) 'Generating solution for unit sigma'&
                        &, ip
                call gumake(iusig(ip), aicsig(1, j))
                sigu(ip, iusig(ip)) = 1.0
            endif
        enddo
        !
    end subroutine gucalc
    !*==GUMAKE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GUCALC



    subroutine gumake(iu, aic)
        use i_dfdc
        use m_gauss, only : baksub
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iu
        real, dimension(*) :: aic
        !
        ! Local variables
        !
        integer :: iel, ip, j, ku
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     Finds the index IU of the first available unit r.h.s. vector
        !     Generates GAMU(.IU),SIGU(.IU),QNDOFU(.IU) unit solutions.
        !-------------------------------------------------------------------
        !
        !---- find first available vector index KU
        do ku = 1, nu
            if (.not.luset(ku)) goto 10
        enddo
        !
        !---- no free vectors found... add new one to list
        nu = nu + 1
        if (nu>nux) write (*, *)                                         &
                &'GUMAKE: Array limit exceeded. Increase NUX to'&
                &, nu
        !
        ku = min(nu, nux)
        !
        !-----------------------------------
        !
        !---- set Jacobian vector
        10   do j = 1, nsys
            res(j) = aic(j)
        enddo
        !
        !---- calculate response
        if (nsys>0) call baksub(nsysx, nsys, sys(1, 1), res(1))
        res(0) = 0.
        !
        !---- first clear entire response vectors
        do ip = 1, nptot
            gamu(ip, ku) = 0.
            sigu(ip, ku) = 0.
        enddo
        do iel = 1, nel
            qndofu(iel, ku) = 0.
        enddo
        !
        !---- store response of surface gamma
        do ip = 1, nptot
            j = jsysgam(ip)
            gamu(ip, ku) = -res(j)
        enddo
        !
        !---- store response of normal-velocity DOFs
        do iel = 1, nel
            j = jsysqnc(iel)
            qndofu(iel, ku) = -res(j)
        enddo
        !
        luset(ku) = .true.
        !
        iu = ku
        !
    end subroutine gumake
    !*==GSOLVE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GUMAKE



    subroutine gsolve
        use i_dfdc
        use m_gauss, only : baksub
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real, dimension(ipx) :: gamk, sigk
        integer :: iel, ip, j, k
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     Does a direct solve from current RHS (knowns)
        !     Generates GAM(.),SIG(.),QNDOF(.)
        !-------------------------------------------------------------------
        !
        !---- set imposed (known) gamma and sigma strengths
        do ip = 1, nptot
            gamk(ip) = gamvsp(ip)
            sigk(ip) = sigvsp(ip)
        enddo
        !
        do iel = 1, nel
            do ip = ipfrst(iel), iplast(iel)
                gamk(ip) = gamk(ip) + gamset(iel)
                sigk(ip) = sigk(ip) + sigset(iel)
            enddo
        enddo
        !
        !---- set Jacobian vector (RHS)
        do k = 1, nsys
            res(k) = aicqff(k, 1) * qinf
            do ip = 1, nptot
                j = jaicgth(ip)
                res(k) = res(k) + aicgth(k, j) * gth(ip)
            enddo
            do ip = 1, nptot
                j = jaicgam(ip)
                res(k) = res(k) + aicgam(k, j) * gamk(ip)
            enddo
            do ip = 1, nptot
                j = jaicsig(ip)
                res(k) = res(k) + aicsig(k, j) * sigk(ip)
            enddo
        enddo
        !
        !---- calculate solution
        if (nsys>0) then
            call baksub(nsysx, nsys, sys(1, 1), res(1))
        else
            write (*, *) 'GSOLVE error NSYS=', nsys
        endif
        res(0) = 0.
        !
        !---- first clear unknown solution vectors
        do ip = 1, nptot
            gam(ip) = 0.
            sig(ip) = 0.
        enddo
        do iel = 1, nel
            qndof(iel) = 0.
        enddo
        !
        !---- store solution for surface gamma + known gamma values
        !---- store solution for surface sigma + known sigma values
        do ip = 1, nptot
            j = jsysgam(ip)
            gam(ip) = -res(j) + gamk(ip)
            sig(ip) = sig(ip) + sigk(ip)
        enddo
        !
        !---- store solution for normal-velocity DOFs
        do iel = 1, nel
            j = jsysqnc(iel)
            qndof(iel) = -res(j)
        enddo
        !
        !---- Set strength of TE panels (for debugging use)
        call setgste
        !
        lgama = .true.
        lqcnt = .false.
        !
    end subroutine gsolve
    !*==GSOLVE0.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GSOLVE


    subroutine gsolve0
        use i_dfdc
        use m_gauss, only : baksub
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real, dimension(ipx) :: gamk, sigk
        integer :: iel, ip, j, k
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     Does a direct solve from current RHS (knowns)
        !     Generates GAM(.),SIG(.),QNDOF(.)
        !-------------------------------------------------------------------
        !
        !---- set known gamma and sigma strengths
        do ip = 1, nptot
            gamk(ip) = gamvsp(ip)
            sigk(ip) = sigvsp(ip)
        enddo
        !
        do iel = 1, nel
            do ip = ipfrst(iel), iplast(iel)
                gamk(ip) = gamk(ip) + gamset(iel)
                sigk(ip) = sigk(ip) + sigset(iel)
            enddo
        enddo
        !
        !---- set Jacobian vector (RHS)
        do k = 1, nsys
            res(k) = aicqff(k, 1) * qinf
            !        DO IR = 1, NRP
            !          RES(K) = RES(K) + AICVWK(K,IR)*WBGAM(IR)
            !        ENDDO
            do ip = 1, nptot
                j = jaicgam(ip)
                res(k) = res(k) + aicgam(k, j) * gamk(ip)
            enddo
            do ip = 1, nptot
                j = jaicsig(ip)
                res(k) = res(k) + aicsig(k, j) * sigk(ip)
            enddo
        enddo
        !
        !---- calculate response
        if (nsys>0) call baksub(nsysx, nsys, sys(1, 1), res(1))
        res(0) = 0.
        !
        !---- first clear unknown solution vectors
        do ip = 1, nptot
            gam(ip) = 0.
            gth(ip) = 0.
            sig(ip) = 0.
        enddo
        do iel = 1, nel
            qndof(iel) = 0.
        enddo
        !
        !---- store solution for surface gamma + known gamma values
        !---- store solution for surface sigma + known sigma values
        do ip = 1, nptot
            j = jsysgam(ip)
            gam(ip) = -res(j) + gamk(ip)
            sig(ip) = sig(ip) + sigk(ip)
            !---- set GTH for wake point from WBGAM for wake row
            !        IF(IP2IR(IP).NE.0) THEN
            !          IR = IP2IR(IP)
            !          GTH(IP) = WBGAM(IR)/YP(IP)
            !        ENDIF
        enddo
        !
        !---- store solution for normal-velocity DOFs
        do iel = 1, nel
            j = jsysqnc(iel)
            qndof(iel) = -res(j)
        enddo
        !
        !---- Set strength of TE panels (for debugging use)
        call setgste
        lgama = .true.
        lqcnt = .false.
        !
    end subroutine gsolve0
    !*==GUVWK.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GSOLVE0


    subroutine guvwk(iradd)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iradd
        !
        ! Local variables
        !
        integer :: ir, jp, ku
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
        do jp = 1, nptot
            !---- clear unit response
            ku = iuvwk(iradd)
            gthu(jp, ku) = 0.
            !
            ir = ip2ir(jp)
            if (ir/=0) then
                !---- we are on a rotor wake element, use rotor index to locate influence
                if (ir==iradd) then
                    if (iuvwk(ir)>0) then
                        ku = iuvwk(ir)
                        gthu(jp, ku) = 1.0 / yp(jp)
                    else
                        write (*, *) 'GUVWK: Unassigned rotor RHS IR =', ir
                    endif
                endif
                !
            endif
        enddo
        !
    end subroutine guvwk
    !*==GUSUM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GUVWK



    subroutine gusum
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: iel, ip, iu
        real, dimension(0:nux) :: uval
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
        if (ldbg) write (*, *)                                           &
                &'Superimposing unit-strength distributions'
        !
        do iu = 0, nu
            uval(iu) = 0.
        enddo
        !
        !---- put prescribed strengths into indexed array (inactive quantity has IU=0)
        !
        !----- axisymmetric case
        uval(iuqinf) = qinf
        !
        !---- Rotor wake circulation (actually Wt/Wm *delta BGAM /2*pi)
        !      DO IR = 1, NRP
        !        UVAL(IUVWK(IR)) = WBGAM(IR)
        !c       write(18,*) ir,iuvwk(ir),wbgam(ir)
        !      ENDDO
        !
        !---- Additional gamma specified from viscous response (from BL effects)
        do ip = 1, nptot
            !c        write(24,*) IP,IUGAM(IP),GAMVSP(IP)
            uval(iugam(ip)) = gamvsp(ip)
            uval(iusig(ip)) = sigvsp(ip)
        enddo
        !
        do iel = 1, nel
            do ip = ipfrst(iel), iplast(iel)
                uval(iugam(ip)) = uval(iugam(ip)) + gamset(iel)
                uval(iusig(ip)) = uval(iusig(ip)) + sigset(iel)
            enddo
        enddo
        !
        !---- clear all singularities for accumulation
        do ip = 1, nptot
            gam(ip) = 0.
            sig(ip) = 0.
            gth(ip) = 0.
        enddo
        do iel = 1, nel
            qndof(iel) = 0.
        enddo
        !
        !---- accumulate unit contributions weighted by prescribed strengths
        do iu = 1, nu
            !------ no sense accumulating zeros
            if (uval(iu)/=0.0) then
                !c          write(20,*) iu,uval(iu)
                !
                do ip = 1, nptot
                    gam(ip) = gam(ip) + gamu(ip, iu) * uval(iu)
                    sig(ip) = sig(ip) + sigu(ip, iu) * uval(iu)
                    gth(ip) = gth(ip) + gthu(ip, iu) * uval(iu)
                    !c          write(20,99) ip,iu,gamu(ip,iu),GTHU(IP,IU)
                enddo
                do iel = 1, nel
                    qndof(iel) = qndof(iel) + qndofu(iel, iu) * uval(iu)
                enddo
            endif
            !
        enddo
        !
        !      DO IP = 1, NPTOT
        !        write(24,*) IP,GAM(IP),GAMVSP(IP)
        !      ENDDO
        !
        !---- Set strength of TE panels (for debugging use)
        call setgste
        !
        99   format (2i5, 5(1x, f12.6))
        !
        lgama = .true.
        lqcnt = .false.
        !
    end subroutine gusum
    !*==QCSUM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GUSUM



    subroutine qcsum
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: ic, jp
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Computes velocities at control points for
        !     current panel and point-singularity strengths.
        !-----------------------------------------------------
        !
        if (ldbg) write (*, *) 'Computing control-point velocities'
        !
        do ic = 1, nctot
            if (ictype(ic)==-999) then
                write (*, *) '? unassigned control point ', ic
                !
            elseif (ictype(ic)==-1) then
                !------- zero-length panel (skip)
                !
            elseif (ictype(ic)/=-2) then
                !------- vortex wake panels w/o control points (skip)
                !
                !------- usual velocity contributions for all vortex and source strengths
                qc(1, ic) = 0.
                qc(2, ic) = 0.
                !
                do jp = 1, nptot
                    qc(1, ic) = qc(1, ic) + qc_gam(1, jp, ic) * gam(jp)            &
                            & + qc_sig(1, jp, ic) * sig(jp)
                    qc(2, ic) = qc(2, ic) + qc_gam(2, jp, ic) * gam(jp)            &
                            & + qc_sig(2, jp, ic) * sig(jp)
                    !---- Add contribution from rotor wake vorticity
                    qc(1, ic) = qc(1, ic) + qc_gth(1, jp, ic) * gth(jp)
                    qc(2, ic) = qc(2, ic) + qc_gth(2, jp, ic) * gth(jp)
                enddo
                !
                !---- Add freestream flow contribution
                qc(1, ic) = qc(1, ic) + qinf
            endif
            !
        enddo
        !
        lqcnt = .true.
        !
    end subroutine qcsum
    !*==QCUSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCSUM




    subroutine qcuset
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: ic, ip, iu, jp
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Computes velocities at control points for
        !     current unit-forcing solutions
        !-----------------------------------------------------
        !
        if (ldbg) write (*, *)                                           &
                &'Computing control-pt unit-forcing velocities'
        !
        do ic = 1, nctot
            if (ictype(ic)==-999) then
                write (*, *) '? unassigned control point ', ic
                !
            elseif (ictype(ic)==-1) then
                !------- zero-length panel (skip)
                !
            elseif (ictype(ic)/=-2) then
                !------- vortex wake panel w/o control point (skip)
                !
                !------- usual velocity contributions for all vortex strengths
                do iu = 1, nu
                    qcu(1, ic, iu) = 0.
                    qcu(2, ic, iu) = 0.
                    !
                    do jp = 1, nptot
                        !c             IF(JSYSGAM(JP).GT.0) THEN
                        !---- Contribution from solution unit responses
                        qcu(1, ic, iu) = qcu(1, ic, iu) + qc_gam(1, jp, ic)         &
                                & * gamu(jp, iu)
                        qcu(2, ic, iu) = qcu(2, ic, iu) + qc_gam(2, jp, ic)         &
                                & * gamu(jp, iu)
                        !---- Contribution from rotor wake vorticity unit forcing
                        qcu(1, ic, iu) = qcu(1, ic, iu) + qc_gth(1, jp, ic)         &
                                & * gthu(jp, iu)
                        qcu(2, ic, iu) = qcu(2, ic, iu) + qc_gth(2, jp, ic)         &
                                & * gthu(jp, iu)
                        !c             ENDIF
                    enddo
                enddo
                !
                !---- Contribution from freestream
                iu = iuqinf
                qcu(1, ic, iu) = qcu(1, ic, iu) + 1.0
                !
                !---- Contribution from specified vortex and source panels
                do ip = 1, nptot
                    iu = iugam(ip)
                    qcu(1, ic, iu) = qcu(1, ic, iu) + qc_gam(1, ip, ic)
                    qcu(2, ic, iu) = qcu(2, ic, iu) + qc_gam(2, ip, ic)
                    !
                    iu = iusig(ip)
                    qcu(1, ic, iu) = qcu(1, ic, iu) + qc_sig(1, ip, ic)
                    qcu(2, ic, iu) = qcu(2, ic, iu) + qc_sig(2, ip, ic)
                enddo
                !
            endif
        enddo
        !
    end subroutine qcuset
    !*==QCUSUM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCUSET



    subroutine qcusum
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: ic, iel, ip, iu
        real, dimension(0:nux) :: uval
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Computes velocities at control points for
        !     current panel and point-singularity strengths.
        !-----------------------------------------------------
        !
        if (ldbg) write (*, *) 'Superimposing unit-strength velocities'
        !
        do iu = 0, nu
            uval(iu) = 0.
        enddo
        !
        !---- put prescribed strengths into indexed array (inactive quantity has IU=0)
        !
        !----- freestream
        uval(iuqinf) = qinf
        !
        !---- Rotor wake circulation (actually Wt/Wm *delta BGAM /2*pi)
        !      DO IR = 1, NRP
        !        UVAL(IUVWK(IR)) = WBGAM(IR)
        !c       write(18,*) ir,iuvwk(ir),wbgam(ir)
        !      ENDDO
        !
        do ip = 1, nptot
            uval(iugam(ip)) = gamvsp(ip)
            uval(iusig(ip)) = sigvsp(ip)
        enddo
        !---- Specified singularities
        do iel = 1, nel
            do ip = ipfrst(iel), iplast(iel)
                uval(iugam(ip)) = uval(iugam(ip)) + gamset(iel)
                uval(iusig(ip)) = uval(iusig(ip)) + sigset(iel)
            enddo
        enddo
        !
        !---- accumulate unit contributions weighted by prescribed strengths
        do ic = 1, nctot
            qc(1, ic) = 0.
            qc(2, ic) = 0.
        enddo
        !
        do iu = 1, nu
            !------ no sense accumulating zeros
            if (uval(iu)/=0.0) then
                !
                do ic = 1, nctot
                    qc(1, ic) = qc(1, ic) + qcu(1, ic, iu) * uval(iu)
                    qc(2, ic) = qc(2, ic) + qcu(2, ic, iu) * uval(iu)
                enddo
            endif
            !
        enddo
        !
        lqcnt = .true.
        !
    end subroutine qcusum
    !*==SETGSTE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCUSUM




    subroutine setgste
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: iel, ip1, ip2, k
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Sets TE vortex and source strength
        !     from vortex and source strength on body panels
        !-----------------------------------------------------
        !
        if (ldbg) write (*, *) 'Setting TE GAMT,SIGT'
        !
        !---- Set TE panel strengths
        !c      write(30,*) 'SETGSTE IEL,GAMT1-2,SIGT1-2'
        !
        do iel = 1, nel
            do k = 1, 2
                gamt(k, iel) = 0.
                sigt(k, iel) = 0.
            enddo
            !
            !---- Body elements with GAM
            if (ltpan(iel) .and. netype(iel)==0) then
                ip1 = ipte1(iel)
                ip2 = ipte2(iel)
                !
                if (ip1/=0) then
                    do k = 1, 2
                        gamt(k, iel) = gamt(k, iel) + gamt_gam(k, 1, iel) * gam(ip1)&
                                & + gamt_sig(k, 1, iel) * sig(ip1)
                        sigt(k, iel) = sigt(k, iel) + sigt_gam(k, 1, iel) * gam(ip1)&
                                & + sigt_sig(k, 1, iel) * sig(ip1)
                    enddo
                endif
                if (ip2/=0) then
                    do k = 1, 2
                        gamt(k, iel) = gamt(k, iel) + gamt_gam(k, 2, iel) * gam(ip2)&
                                & + gamt_sig(k, 2, iel) * sig(ip2)
                        sigt(k, iel) = sigt(k, iel) + sigt_gam(k, 2, iel) * gam(ip2)&
                                & + sigt_sig(k, 2, iel) * sig(ip2)
                    enddo
                endif
            endif
            !
            !---- vortex wake elements with GTH specified
            if (ltpan(iel) .and. netype(iel)==7) then
                ip1 = ipte1(iel)
                ip2 = ipte2(iel)
                !
                if (ip1/=0) then
                    do k = 1, 2
                        gamt(k, iel) = gamt(k, iel) + gamt_gam(k, 1, iel) * gth(ip1)
                        sigt(k, iel) = sigt(k, iel) + sigt_gam(k, 1, iel) * gth(ip1)
                    enddo
                endif
                if (ip2/=0) then
                    do k = 1, 2
                        gamt(k, iel) = gamt(k, iel) + gamt_gam(k, 2, iel) * gth(ip2)
                        sigt(k, iel) = sigt(k, iel) + sigt_gam(k, 2, iel) * gth(ip2)
                    enddo
                endif
            endif
            !
            !c        write(30,99) iel,GAMT(1,IEL),GAMT(2,IEL),SIGT(1,IEL),SIGT(2,IEL)
            99      format (i5, 4(1x, f12.6))
        enddo
        !
    end subroutine setgste
    !*==NWDOTS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETGSTE



    subroutine nwdots(iel, ds1nw, ds2nw, dn1nw, dn2nw)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: dn1nw, dn2nw, ds1nw, ds2nw
        integer :: iel
        !
        ! Local variables
        !
        real, dimension(2) :: anw
        real :: anwabs
        integer :: ic1, ic2
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        !
        !      IELW = IEWAKE(IEL)
        !      IF(IELW .NE. 0) THEN
        !       ICW = ICFRST(IELW)
        !       ANW(1) = ANC(1,ICW)
        !       ANW(2) = ANC(2,ICW)
        !      ELSE
        anw(1) = anc(1, ic2) - anc(1, ic1)
        anw(2) = anc(2, ic2) - anc(2, ic1)
        anwabs = sqrt(anw(1)**2 + anw(2)**2)
        anw(1) = anw(1) / anwabs
        anw(2) = anw(2) / anwabs
        !      ENDIF
        !
        ds1nw = anc(2, ic1) * anw(1) - anc(1, ic1) * anw(2)
        ds2nw = anc(2, ic2) * anw(1) - anc(1, ic2) * anw(2)
        !
        dn1nw = anc(1, ic1) * anw(1) + anc(2, ic1) * anw(2)
        dn2nw = anc(1, ic2) * anw(1) + anc(2, ic2) * anw(2)
        !
    end subroutine nwdots
    !*==SWDOTS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! NWDOTS


    subroutine swdots(iel, ds1sw, ds2sw, dn1sw, dn2sw)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: dn1sw, dn2sw, ds1sw, ds2sw
        integer :: iel
        !
        ! Local variables
        !
        real, dimension(2) :: anw
        real :: anwabs
        integer :: ic1, ic2
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        !
        !      IELW = IEWAKE(IEL)
        !      IF(IELW .NE. 0) THEN
        !       ICW = ICFRST(IELW)
        !       ANW(1) = ANC(1,ICW)
        !       ANW(2) = ANC(2,ICW)
        !      ELSE
        anw(1) = anc(1, ic2) - anc(1, ic1)
        anw(2) = anc(2, ic2) - anc(2, ic1)
        anwabs = sqrt(anw(1)**2 + anw(2)**2)
        anw(1) = anw(1) / anwabs
        anw(2) = anw(2) / anwabs
        !      ENDIF
        !
        ds1sw = anc(2, ic1) * anw(2) + anc(1, ic1) * anw(1)
        ds2sw = anc(2, ic2) * anw(2) + anc(1, ic2) * anw(1)
        !
        dn1sw = anc(1, ic1) * anw(2) - anc(2, ic1) * anw(1)
        dn2sw = anc(1, ic2) * anw(2) - anc(2, ic2) * anw(1)
        !
    end subroutine swdots
end module m_system
