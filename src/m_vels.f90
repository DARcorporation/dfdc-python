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


    subroutine qcpfor
        use i_dfdc
        use m_forces, only : cpcalc, fcoeff
        use m_system, only : qcsum
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real, dimension(2) :: ant
        real :: circ, cmte, cosa, cxte, cyte, &
                & deltacp, dxt, dyt, qrsq, que, sina, vtsq, &
                & vtt
        integer :: i, ic, ict1, ict2, iel, ig, ip1, ip2, ipt1, &
                & ipt2, ir, ir1, ir2, iu, jg, n
        real, dimension(1) :: cptl, cptr, dsct, xct, yct
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------
        !     Computes velocities, Cp, forces, etc.
        !     for current singularity distribitions
        !---------------------------------------------
        !
        !---- compute velocities at control points (average between two panel sides)
        if (.not.lqcnt) call qcsum
        !
        !---- set velocities on both sides of panels at all control points
        call qqcalc(nctot, anc, ipco, ipcp, gam, sig, gth, qc, qcl, qcr)
        !
        !---- also set QCL,QCR sensitivities to imposed freestream, singularities, etc.
        do iu = 1, nu
            call qqcalc(nctot, anc, ipco, ipcp, gamu(1, iu), sigu(1, iu), &
                    & gthu(1, iu), qcu(1, 1, iu), qclu(1, 1, iu), qcru(1, 1, iu))
        enddo
        !
        !---- set velocities and Cp on both sides of panels at all control points
        call cpcalc(nctot, qcl, qinf, qref, cpl, cpl_qcl, cpl_qinf)
        call cpcalc(nctot, qcr, qinf, qref, cpr, cpr_qcr, cpr_qinf)
        !
        !---- add in contribution from rotor enthalpy and entropy losses to Cp's
        !     CB and duct walls and on upper and lower wakes
        qrsq = qref**2
        do ic = 1, nctot
            ip1 = ipco(ic)
            ip2 = ipcp(ic)
            ir1 = ip2ir(ip1)
            ir2 = ip2ir(ip2)
            ig = ic2ig(ic)
            !---- check for control point in rotor slipstream (in grid)
            if (ig>0) then
                ir = ir1
                if (ir==1) then
                    jg = ir
                    if (ip2<=iplast(1)) then
                        !---- on CB downstream of rotor
                        circ = bgamg(ig, jg)
                        vtt = circ * pi2i / yc(ic)
                        vtsq = vtt**2
                        deltacp = 2.0 * (dhg(ig, jg) - dsg(ig, jg)) / qrsq - vtsq / qrsq
                        cpr(ic) = cpr(ic) + deltacp
                    else
                        !---- on CB wake element
                        circ = bgamg(ig, jg)
                        vtt = circ * pi2i / yc(ic)
                        vtsq = vtt**2
                        deltacp = 2.0 * (dhg(ig, jg) - dsg(ig, jg)) / qrsq - vtsq / qrsq
                        cpl(ic) = cpl(ic) + deltacp
                    endif
                    !
                elseif (ir==nrp) then
                    jg = ir - 1
                    if (ip1<iplast(2)) then
                        !---- on DUCT downstream of rotor
                        circ = bgamg(ig, jg)
                        vtt = circ * pi2i / yc(ic)
                        vtsq = vtt**2
                        deltacp = 2.0 * (dhg(ig, jg) - dsg(ig, jg)) / qrsq - vtsq / qrsq
                        cpr(ic) = cpr(ic) + deltacp
                    else
                        !---- on DUCT wake element
                        circ = bgamg(ig, jg)
                        vtt = circ * pi2i / yc(ic)
                        vtsq = vtt**2
                        deltacp = 2.0 * (dhg(ig, jg) - dsg(ig, jg)) / qrsq - vtsq / qrsq
                        cpr(ic) = cpr(ic) + deltacp
                    endif
                else
                    jg = ir - 1
                    circ = bgamg(ig, jg)
                    vtt = circ * pi2i / yc(ic)
                    vtsq = vtt**2
                    deltacp = 2.0 * (dhg(ig, jg) - dsg(ig, jg)) / qrsq - vtsq / qrsq
                    cpr(ic) = cpr(ic) + deltacp
                    !
                    jg = ir
                    circ = bgamg(ig, jg)
                    vtt = circ * pi2i / yc(ic)
                    vtsq = vtt**2
                    deltacp = 2.0 * (dhg(ig, jg) - dsg(ig, jg)) / qrsq - vtsq / qrsq
                    cpl(ic) = cpl(ic) + deltacp
                endif
            endif
        enddo
        !
        !---- also set CPL,CPR sensitivities
        do iu = 1, nu
            do ic = 1, nctot
                cplu(ic, iu) = cpl_qcl(1, ic) * qclu(1, ic, iu) + cpl_qcl(2, ic)   &
                        & * qclu(2, ic, iu)
                cpru(ic, iu) = cpr_qcr(1, ic) * qcru(1, ic, iu) + cpr_qcr(2, ic)   &
                        & * qcru(2, ic, iu)
            enddo
        enddo
        !
        !---- integrate Cp to get forces and moment, and all sensitivities
        cosa = 1.0
        sina = 0.0
        cx(0) = 0.
        cy(0) = 0.
        cm(0) = 0.
        cd(0) = 0.
        cxvis(0) = 0.
        do iu = 1, nu
            cxu(0, iu) = 0.
            cyu(0, iu) = 0.
            cmu(0, iu) = 0.
            cdu(0, iu) = 0.
        enddo
        do iel = 1, nel
            i = icfrst(iel)
            n = iclast(iel) - icfrst(iel) + 1
            if (n>=1) then
                call fcoeff(n, cpl(i), cpr(i), xc(i), yc(i), anc(1, i), dsc(i), &
                        & cx(iel), cy(iel), cm(iel))
                !
                cd(iel) = cx(iel)
                !
                !---- Totals go in index 0
                cx(0) = cx(0) + cx(iel)
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
                iu = iuqinf
                if (luset(iu)) then
                    call fcoeff(n, cplu(i, iu), cpru(i, iu), xc(i), yc(i), anc(1, i), &
                            & dsc(i), cxu(iel, iu), cyu(iel, iu), cmu(iel, iu))
                    cdu(iel, iu) = cxu(iel, iu)
                    !---- Totals go in index 0
                    cxu(0, iu) = cxu(0, iu) + cxu(iel, iu)
                    !           CYU(0,IU) = CYU(0,IU) + CYU(IEL,IU)
                    !           CMU(0,IU) = CMU(0,IU) + CMU(IEL,IU)
                    !           CDU(0,IU) = CDU(0,IU) + CDU(IEL,IU)
                endif
                !
                !---- Add in viscous forces if BL has been calculated
                if (lvisc) then
                    cx(iel) = cx(iel) + cxvis(iel)
                    cd(iel) = cx(iel)
                    !---- Add viscous forces to total forces in index 0
                    cx(0) = cx(0) + cxvis(iel)
                    cxvis(0) = cxvis(0) + cxvis(iel)
                endif
                !
                !
                !---- Elements with TE panels
                !     Add forces due to pressures on TE panel
                if (ltpan(iel) .and. netype(iel)==0) then
                    ipt1 = ipte1(iel)
                    ipt2 = ipte2(iel)
                    ict1 = 0
                    ict2 = 0
                    if (ipt1==ipfrst(iel)) ict1 = icfrst(iel)
                    if (ipt1==iplast(iel)) ict1 = iclast(iel)
                    if (ipt2==ipfrst(iel)) ict2 = icfrst(iel)
                    if (ipt2==iplast(iel)) ict2 = iclast(iel)
                    !          write(12,*) 'iel ',iel
                    !          write(12,*) 'ipt1,ipt2 ',ipt1,ipt2
                    !          write(12,*) 'ipfrst,iplast ',ipfrst(iel),iplast(iel)
                    !          write(12,*) 'ict1,ict2 ',ict1,ict2
                    !          write(12,*) 'icfrst,iclast ',icfrst(iel),iclast(iel)
                    !
                    if (ict1/=0 .and. ict2/=0) then
                        cptl = 0.5 * (cpl(ict1) + cpl(ict2))
                        cptr = 0.5 * (cpr(ict1) + cpr(ict2))
                    elseif (ict1==0) then
                        cptl = cpl(ict2)
                        cptr = cpr(ict2)
                    elseif (ict2==0) then
                        cptl = cpl(ict1)
                        cptr = cpr(ict1)
                    endif
                    xct = 0.5 * (xpt(1, iel) + xpt(2, iel))
                    yct = 0.5 * (ypt(1, iel) + ypt(2, iel))
                    dxt = xpt(1, iel) - xpt(2, iel)
                    dyt = ypt(1, iel) - ypt(2, iel)
                    dsct = sqrt(dxt * dxt + dyt * dyt)
                    ant(1) = dyt / dsct(1)
                    ant(2) = -dxt / dsct(1)
                    !          write(12,*) 'cptl,cptr ',cptl,cptr
                    !          write(12,*) 'xct,yct ',xct,yct
                    !          write(12,*) 'dsct ',dsct
                    !          write(12,*) 'ant ',ant(1),ant(2)
                    !
                    call fcoeff(1, cptl, cptr, xct, yct, ant, dsct, cxte, cyte, cmte)
                    !          write(12,*) 'CXTE ',cxte
                    !          write(12,*) 'CYTE ',cyte
                    !          write(12,*) 'CMTE ',cmte
                    !
                    !---- Add in TE forces
                    cx(iel) = cx(iel) + cxte
                    !          CY(IEL) = CY(IEL) + CYTE
                    !          CD(IEL) = CX(IEL)
                    !---- Add TE forces to total forces in index 0
                    cx(0) = cx(0) + cxte
                    !          CY(0) = CY(0) + CYTE
                    !          CM(0) = CM(0) + CMTE
                    !          CD(0) = CD(0) + CXTE

                    !------- also set sensitivities due to freestream
                    iu = iuqinf
                    if (luset(iu)) then
                        if (ict1/=0 .and. ict2/=0) then
                            cptl = 0.5 * (cplu(ict1, iu) + cplu(ict2, iu))
                            cptr = 0.5 * (cpru(ict1, iu) + cpru(ict2, iu))
                        elseif (ict1==0) then
                            cptl = cplu(ict2, iu)
                            cptr = cpru(ict2, iu)
                        elseif (ict2==0) then
                            cptl = cplu(ict1, iu)
                            cptr = cpru(ict1, iu)
                        endif
                        call fcoeff(1, cptl, cptr, xct, yct, ant, dsct, cxte, cyte, &
                                & cmte)
                        cxu(iel, iu) = cxu(iel, iu) + cxte
                        !            CYU(IEL,IU) = CYU(IEL,IU) + CYTE
                        !            CDU(IEL,IU) = CXU(IEL,IU)
                        !---- Totals go in index 0
                        cxu(0, iu) = cxu(0, iu) + cxte
                        !            CYU(0,IU) = CYU(0,IU) + CYTE
                        !            CMU(0,IU) = CMU(0,IU) + CMTE
                        !            CDU(0,IU) = CDU(0,IU) + CXTE
                    endif
                endif
                !
            endif
        enddo
        !
        !---- Save duct + centerbody (IEL=1,2) X force (thrust)
        que = 0.5 * rho * qref**2
        tduct = -(cx(1) + cx(2)) * que
        !
    end subroutine qcpfor
    !*==QQCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCPFOR



    subroutine qqcalc(n, an, ipo, ipp, gam, sig, gth, q, ql, qr)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(2, *) :: an, q, ql, qr
        real, dimension(*) :: gam, gth, sig
        integer, dimension(*) :: ipo, ipp
        !
        ! Local variables
        !
        real :: dqn, dqt, gamc, sigc
        integer :: i, io, ip
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Sets velocity on both sides of vortex,source sheet.
        !-------------------------------------------------------------
        !
        do i = 1, n
            io = ipo(i)
            ip = ipp(i)
            !
            gamc = 0.5 * (gam(io) + gam(ip))
            sigc = 0.5 * (sig(io) + sig(ip))
            !---- Add vorticity contribution from additional vortex wake gamma
            gamc = gamc + 0.5 * (gth(io) + gth(ip))
            !
            !------ 1/2 of sheet strength adds,subtracts to left,right speed
            dqt = 0.5 * gamc
            dqn = 0.5 * sigc
            !
            qr(1, i) = q(1, i) + an(1, i) * dqn + an(2, i) * dqt
            qr(2, i) = q(2, i) + an(2, i) * dqn - an(1, i) * dqt
            !
            ql(1, i) = q(1, i) - an(1, i) * dqn - an(2, i) * dqt
            ql(2, i) = q(2, i) - an(2, i) * dqn + an(1, i) * dqt
        enddo
        !
    end subroutine qqcalc
    !*==QTCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QQCALC



    subroutine qtcalc(n, an, ipo, ipp, gam, gth, q, qt)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(2, *) :: an, q
        real, dimension(*) :: gam, gth, qt
        integer, dimension(*) :: ipo, ipp
        !
        ! Local variables
        !
        real, dimension(2) :: at
        real :: dqt, gamc
        integer :: i, io, ip
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Sets tangential velocity component
        !-------------------------------------------------------------
        !
        !
        do i = 1, n
            io = ipo(i)
            ip = ipp(i)
            !
            gamc = 0.5 * (gam(io) + gam(ip))
            !---- Add vorticity from vortex wake gamma
            gamc = gamc + 0.5 * (gth(io) + gth(ip))
            !
            !------ set tangential vector along sheet (AN  points to right side)
            at(1) = -an(2, i)
            at(2) = an(1, i)
            !
            !------ 1/2 of sheet strength subtracts to right-side speed (KVSIDE=+1)
            dqt = -0.5 * gamc
            !
            qt(i) = q(1, i) * at(1) + q(2, i) * at(2) + dqt
        enddo
        !
    end subroutine qtcalc
    !*==QJAC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QTCALC



    subroutine qjac(dx, dy, du, dv, qj, qj_dx, qj_dy, qj_du, qj_dv)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: du, dv, dx, dy
        real, dimension(2, 2) :: qj, qj_du, qj_dv, qj_dx, qj_dy
        !
        ! Local variables
        !
        real :: dsq, dsqi
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------------
        !     Computes velocity jacobian from position,velocity deltas,
        !     and zero velocity div and curl conditions.
        !------------------------------------------------------------------
        !
        dsq = dx**2 + dy**2
        if (dsq==0.0) then
            dsqi = 1.0
        else
            dsqi = 1.0 / dsq
        endif
        !
        qj(1, 1) = (du * dx - dv * dy) * dsqi
        qj(1, 2) = (du * dy + dv * dx) * dsqi
        !
        qj_dx(1, 1) = (du - 2.0 * qj(1, 1) * dx) * dsqi
        qj_dy(1, 1) = (-dv - 2.0 * qj(1, 1) * dy) * dsqi
        !
        qj_dx(1, 2) = (dv - 2.0 * qj(1, 2) * dx) * dsqi
        qj_dy(1, 2) = (du - 2.0 * qj(1, 2) * dy) * dsqi
        !
        qj_du(1, 1) = dx * dsqi
        qj_dv(1, 1) = -dy * dsqi
        !
        qj_du(1, 2) = dy * dsqi
        qj_dv(1, 2) = dx * dsqi
        !
        !
        qj(2, 1) = qj(1, 2)
        qj_dx(2, 1) = qj_dx(1, 2)
        qj_dy(2, 1) = qj_dy(1, 2)
        qj_du(2, 1) = qj_du(1, 2)
        qj_dv(2, 1) = qj_dv(1, 2)
        !
        qj(2, 2) = -qj(1, 1)
        qj_dx(2, 2) = -qj_dx(1, 1)
        qj_dy(2, 2) = -qj_dy(1, 1)
        qj_du(2, 2) = -qj_du(1, 1)
        qj_dv(2, 2) = -qj_dv(1, 1)
        !
    end subroutine qjac
    !*==QCURV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QJAC


    subroutine qcurv(u, v, qj, crv, crv_u, crv_v, crv_qj)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: crv, crv_u, crv_v, u, v
        real, dimension(2, 2) :: crv_qj, qj
        !
        ! Local variables
        !
        real :: crv_qi, crv_un, crv_vn, q, qi, qi_u, qi_v, qsq, &
                & un, un_u, un_v, vn, vn_u, vn_v
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Computes streamline curvature from velocity
        !     and velocity Jacobian.
        !-----------------------------------------------------
        !
        qsq = u**2 + v**2
        q = sqrt(qsq)
        !
        if (q==0.0) then
            qi = 1.0
        else
            qi = 1.0 / q
        endif
        !
        un = u * qi
        vn = v * qi
        un_u = vn * vn * qi
        un_v = -un * vn * qi
        vn_u = -un * vn * qi
        vn_v = un * un * qi
        !
        qi_u = -un * qi**2
        qi_v = -vn * qi**2
        !
        !
        crv = (un * un * qj(2, 1) - un * vn * qj(1, 1) + un * vn * qj(2, 2) - vn * vn * qj(1, 2)) * qi
        !
        !---- CRV( QJ UN(U V) VN(U V) QI(U V) )
        crv_qj(1, 1) = -un * vn * qi
        crv_qj(1, 2) = -vn * vn * qi
        crv_qj(2, 1) = un * un * qi
        crv_qj(2, 2) = un * vn * qi
        !
        crv_un = (2.0 * un * qj(2, 1) - vn * qj(1, 1) + vn * qj(2, 2)) * qi
        crv_vn = (-un * qj(1, 1) + un * qj(2, 2) - 2.0 * vn * qj(1, 2)) * qi
        crv_qi = un * un * qj(2, 1) - un * vn * qj(1, 1) + un * vn * qj(2, 2)            &
                & - vn * vn * qj(1, 2)
        !
        !---- CRV( QJ U V )
        crv_u = crv_un * un_u + crv_vn * vn_u + crv_qi * qi_u
        crv_v = crv_un * un_v + crv_vn * vn_v + crv_qi * qi_v
        !
    end subroutine qcurv
    !*==CV2SET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QCURV



    subroutine cv2set(icm, ico, qc, xc, yc, anc, gamc, sigc, cv, cv_qc, cv_rc, &
            & cv_anc, cv_gamc, cv_sigc)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: cv
        integer :: icm, ico
        real, dimension(2, *) :: anc, qc
        real, dimension(2, -1:0) :: cv_anc, cv_qc, cv_rc
        real, dimension(-1:0) :: cv_gamc, cv_sigc
        real, dimension(*) :: gamc, sigc, xc, yc
        !
        ! Local variables
        !
        real :: cv_dq, cv_dr, gamcm, gamco, qsgn, sigcm, sigco
        real, dimension(2) :: cv_q, cv_qsm, cv_qso, dq, dr, qsa, &
                & qsm, qso
        real, dimension(2, 2) :: cv_qj, qj
        integer :: k
        real, dimension(2, 2, 2) :: qj_dq, qj_dr
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !
        gamcm = gamc(icm)
        sigcm = sigc(icm)
        !
        gamco = gamc(ico)
        sigco = sigc(ico)
        !
        !---- 1/2 of sheet strength adds,subtracts to left,right speed
        qsgn = +0.5
        !
        qsm(1) = qc(1, icm) + qsgn * (anc(1, icm) * sigcm + anc(2, icm) * gamcm)
        qsm(2) = qc(2, icm) + qsgn * (anc(2, icm) * sigcm - anc(1, icm) * gamcm)
        !
        qso(1) = qc(1, ico) + qsgn * (anc(1, ico) * sigco + anc(2, ico) * gamco)
        qso(2) = qc(2, ico) + qsgn * (anc(2, ico) * sigco - anc(1, ico) * gamco)
        !
        dq(1) = qso(1) - qsm(1)
        dq(2) = qso(2) - qsm(2)
        !
        qsa(1) = (qso(1) + qsm(1)) * 0.5
        qsa(2) = (qso(2) + qsm(2)) * 0.5
        !
        dr(1) = xc(ico) - xc(icm)
        dr(2) = yc(ico) - yc(icm)
        !
        !
        call qjac(dr(1), dr(2), dq(1), dq(2), qj, qj_dr(1, 1, 1), qj_dr(1, 1, 2), &
                & qj_dq(1, 1, 1), qj_dq(1, 1, 2))
        !
        call qcurv(qsa(1), qsa(2), qj, cv, cv_q(1), cv_q(2), cv_qj)
        !
        do k = 1, 2
            cv_dq = cv_qj(1, 1) * qj_dq(1, 1, k) + cv_qj(1, 2) * qj_dq(1, 2, k)      &
                    & + cv_qj(2, 1) * qj_dq(2, 1, k) + cv_qj(2, 2) * qj_dq(2, 2, k)
            cv_dr = cv_qj(1, 1) * qj_dr(1, 1, k) + cv_qj(1, 2) * qj_dr(1, 2, k)      &
                    & + cv_qj(2, 1) * qj_dr(2, 1, k) + cv_qj(2, 2) * qj_dr(2, 2, k)
            !
            cv_qsm(k) = -cv_dq + 0.5 * cv_q(k)
            cv_qso(k) = cv_dq + 0.5 * cv_q(k)
            !
            cv_qc(k, -1) = cv_qsm(k)
            cv_qc(k, 0) = cv_qso(k)
            !
            cv_rc(k, -1) = -cv_dr
            cv_rc(k, 0) = cv_dr
            !
            cv_anc(k, -1) = cv_qsm(k) * qsgn * sigcm
            cv_anc(k, 0) = cv_qso(k) * qsgn * sigco
        enddo
        !
        cv_anc(1, -1) = cv_anc(1, -1) - cv_qsm(2) * qsgn * gamcm
        cv_anc(1, 0) = cv_anc(1, 0) - cv_qso(2) * qsgn * gamco
        cv_anc(2, -1) = cv_anc(2, -1) + cv_qsm(1) * qsgn * gamcm
        cv_anc(2, 0) = cv_anc(2, 0) + cv_qso(1) * qsgn * gamco
        !
        cv_sigc(-1) = qsgn * (cv_qsm(1) * anc(1, icm) + cv_qsm(2) * anc(2, icm))
        cv_sigc(0) = qsgn * (cv_qso(1) * anc(1, ico) + cv_qso(2) * anc(2, ico))
        !
        cv_gamc(-1) = qsgn * (cv_qsm(1) * anc(2, icm) - cv_qsm(2) * anc(1, icm))
        cv_gamc(0) = qsgn * (cv_qso(1) * anc(2, ico) - cv_qso(2) * anc(1, ico))
        !
    end subroutine cv2set
    !*==CV3SET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CV2SET


    subroutine cv3set(icm, ico, icp, qc, xc, yc, anc, ansgn, gamc, sigc, cv, &
            & cv_qc, cv_rc, cv_anc, cv_gamc, cv_sigc)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: ansgn, cv
        integer :: icm, ico, icp
        real, dimension(2, *) :: anc, qc
        real, dimension(2, -1:1) :: cv_anc, cv_qc, cv_rc
        real, dimension(-1:1) :: cv_gamc, cv_sigc
        real, dimension(*) :: gamc, sigc, xc, yc
        !
        ! Local variables
        !
        real :: csgn, cv_dq, cv_dr, gamcm, gamco, gamcp, qsgn, &
                & qxn, sigcm, sigco, sigcp
        real, dimension(2) :: cv_q, cv_qsm, cv_qso, cv_qsp, dq, &
                & dr, qsm, qso, qsp
        real, dimension(2, 2) :: cv_qj, qj
        integer :: k
        real, dimension(2, 2, 2) :: qj_dq, qj_dr
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !
        gamcm = gamc(icm)
        sigcm = sigc(icm)
        !
        gamco = gamc(ico)
        sigco = sigc(ico)
        !
        gamcp = gamc(icp)
        sigcp = sigc(icp)
        !
        !---- 1/2 of sheet strength adds,subtracts to left,right speed
        qsgn = 0.5 * ansgn
        !
        qsm(1) = qc(1, icm) + qsgn * (anc(1, icm) * sigcm + anc(2, icm) * gamcm)
        qsm(2) = qc(2, icm) + qsgn * (anc(2, icm) * sigcm - anc(1, icm) * gamcm)
        !
        qso(1) = qc(1, ico) + qsgn * (anc(1, ico) * sigco + anc(2, ico) * gamco)
        qso(2) = qc(2, ico) + qsgn * (anc(2, ico) * sigco - anc(1, ico) * gamco)
        !
        qsp(1) = qc(1, icp) + qsgn * (anc(1, icp) * sigcp + anc(2, icp) * gamcp)
        qsp(2) = qc(2, icp) + qsgn * (anc(2, icp) * sigcp - anc(1, icp) * gamcp)
        !
        dq(1) = qsp(1) - qsm(1)
        dq(2) = qsp(2) - qsm(2)
        !
        dr(1) = xc(icp) - xc(icm)
        dr(2) = yc(icp) - yc(icm)
        !
        qxn = ansgn * (qso(1) * anc(2, ico) - qso(2) * anc(1, ico))
        csgn = sign(1.0, qxn)
        !
        call qjac(dr(1), dr(2), dq(1), dq(2), qj, qj_dr(1, 1, 1), qj_dr(1, 1, 2), &
                & qj_dq(1, 1, 1), qj_dq(1, 1, 2))
        !
        call qcurv(qso(1), qso(2), qj, cv, cv_q(1), cv_q(2), cv_qj)
        !
        cv = csgn * cv
        do k = 1, 2
            cv_q(k) = csgn * cv_q(k)
            cv_qj(k, 1) = csgn * cv_qj(k, 1)
            cv_qj(k, 2) = csgn * cv_qj(k, 2)
        enddo
        !
        do k = 1, 2
            cv_dq = cv_qj(1, 1) * qj_dq(1, 1, k) + cv_qj(1, 2) * qj_dq(1, 2, k)      &
                    & + cv_qj(2, 1) * qj_dq(2, 1, k) + cv_qj(2, 2) * qj_dq(2, 2, k)
            cv_dr = cv_qj(1, 1) * qj_dr(1, 1, k) + cv_qj(1, 2) * qj_dr(1, 2, k)      &
                    & + cv_qj(2, 1) * qj_dr(2, 1, k) + cv_qj(2, 2) * qj_dr(2, 2, k)
            !
            cv_qsm(k) = -cv_dq
            cv_qso(k) = cv_q(k)
            cv_qsp(k) = cv_dq
            !
            cv_qc(k, -1) = cv_qsm(k)
            cv_qc(k, 0) = cv_qso(k)
            cv_qc(k, +1) = cv_qsp(k)
            !
            cv_rc(k, -1) = -cv_dr
            cv_rc(k, 0) = 0.
            cv_rc(k, +1) = cv_dr
            !
            cv_anc(k, -1) = cv_qsm(k) * qsgn * sigcm
            cv_anc(k, 0) = cv_qso(k) * qsgn * sigco
            cv_anc(k, +1) = cv_qsp(k) * qsgn * sigcp
        enddo
        !
        cv_anc(1, -1) = cv_anc(1, -1) - cv_qsm(2) * qsgn * gamcm
        cv_anc(1, 0) = cv_anc(1, 0) - cv_qso(2) * qsgn * gamco
        cv_anc(1, +1) = cv_anc(1, +1) - cv_qsp(2) * qsgn * gamcp
        cv_anc(2, -1) = cv_anc(2, -1) + cv_qsm(1) * qsgn * gamcm
        cv_anc(2, 0) = cv_anc(2, 0) + cv_qso(1) * qsgn * gamco
        cv_anc(2, +1) = cv_anc(2, +1) + cv_qsp(1) * qsgn * gamcp
        !
        cv_sigc(-1) = qsgn * (cv_qsm(1) * anc(1, icm) + cv_qsm(2) * anc(2, icm))
        cv_sigc(0) = qsgn * (cv_qso(1) * anc(1, ico) + cv_qso(2) * anc(2, ico))
        cv_sigc(+1) = qsgn * (cv_qsp(1) * anc(1, icp) + cv_qsp(2) * anc(2, icp))
        !
        cv_gamc(-1) = qsgn * (cv_qsm(1) * anc(2, icm) - cv_qsm(2) * anc(1, icm))
        cv_gamc(0) = qsgn * (cv_qso(1) * anc(2, ico) - cv_qso(2) * anc(1, ico))
        cv_gamc(+1) = qsgn * (cv_qsp(1) * anc(2, icp) - cv_qsp(2) * anc(1, icp))
        !
    end subroutine cv3set
    !*==GETUV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CV3SET


    subroutine getuv(xx, yy, us, vs)
        use i_dfdc
        use m_qaic, only : qfcalc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: us, vs, xx, yy
        !
        ! Local variables
        !
        integer :: ip1, ip2, nf
        integer, dimension(1) :: iftype, ipfo, ipfp
        real, dimension(2) :: qf
        real, dimension(2, ipx) :: qf_gam, qf_gth, qf_sig
        real, dimension(1) :: xx_temp, yy_temp
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
        ip1 = 1
        ip2 = 1
        nf = 1
        !------ assume the point is not on a panel
        ipfo(1) = 0
        ipfp(1) = 0
        iftype(1) = 0
        !------ evaluate velocity components
        call qfcalc(ip1, ip2, xx_temp, yy_temp, ipfo, ipfp, iftype, qf, qf_gam, qf_sig, &
                & qf_gth)
        xx = xx_temp(1)
        yy = yy_temp(1)
        !------ add freestream
        us = qf(1) + qinf
        vs = qf(2)
        !
    end subroutine getuv
end module m_vels
