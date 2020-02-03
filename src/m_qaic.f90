module m_qaic
    implicit none
contains
    !*==QAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PANGET
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

    subroutine qaic(lxyjac)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: lxyjac
        !
        ! Local variables
        !
        integer :: ice, kel
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        if (ldbg) then
            if (lxyjac) then
                write (*, *)                                                 &
                        &'Generating dq/d(gam,sig,x,y) Jacobian matrices...'
            else
                write (*, *) 'Generating dq/d(gam,sig) Jacobian matrices...'
            endif
        endif
        !
        call qaic1(lxyjac, 1, nctot, 1, nel)
        do kel = 1, nel
            ice = ncx + kel
            if (lbody(kel)) call qaic1(lxyjac, ice, ice, 1, nel)
        enddo
        !
        lqcnt = .true.
        lqaic = .true.
        lqgic = lxyjac
        !
    end subroutine qaic
    !*==QAIC1.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QAIC



    subroutine qaic1(lxyjac, ic1, ic2, iel1, iel2)
        use i_dfdc
        use m_aic, only : pntaic, axlaic, panaic, panaicd, linaic
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ic1, ic2, iel1, iel2
        logical :: lxyjac
        !
        ! Local variables
        !
        real :: dspt, dxpt, dypt, splen
        integer :: ic, iel, ip, ip1, ip2, ipo, ipp, j, jpm, &
                & jpo, k, kc, kp, numic
        real, dimension(2, 2, icx) :: qc_gamte, qc_sigte, qc_xpt, &
                & qc_ypt
        real, dimension(ipx) :: rcore
        real, save :: tfrac
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------
        !     Computes velocities and Jacobians at control points IC1...IC2
        !     w.r.t aero quantities at nodes on elements IEL1..IEL2
        !
        !     If LXYJAC=T, also computes Jacobians w.r.t. panel geometry
        !
        !     Input:  XC,YC(i)  control points
        !             XP,YP(j)  panel geometry
        !             GAM(j)    panel vorticity  (or line circ. , or point doublet)
        !             SIG(j)    panel source     (or line source, or point source )
        !             GTH(j)    vortex wake panel vorticity
        !             QINF      freestream speed
        !             ALFA      freestream angle
        !
        !     Output: QC(.i)         velocity at XC(i),YC(i)
        !             QC_GAM(.ji)    dQC(.i)/dGAM(j)
        !             QC_SIG(.ji)
        !             QC_GTH(.ji)    dQC(.i)/dGTH(j)
        !
        !     Additonal output if LXYJAC=T:
        !             QC_XP(.ji)
        !             QC_YP(.ji)
        !             QC_XC(.i)
        !             QC_YC(.i)
        !
        !     Note:   TE panel information and strengths are set by
        !             this routine
        !---------------------------------------------------------------
        !
        !
        !
        !---- maximum TE gap (relative to perimeter) to give closed airfoil a TE panel
        data tfrac/0.00001/
        !
        !---- clear control-point velocities for accumulation in PANAIC and PNTAIC
        do ic = ic1, ic2
            qc(1, ic) = 0.
            qc(2, ic) = 0.
            do iel = iel1, iel2
                do ip = ipfrst(iel), iplast(iel)
                    qc_gam(1, ip, ic) = 0.
                    qc_gam(2, ip, ic) = 0.
                    qc_sig(1, ip, ic) = 0.
                    qc_sig(2, ip, ic) = 0.
                    qc_xp(1, ip, ic) = 0.
                    qc_xp(2, ip, ic) = 0.
                    qc_yp(1, ip, ic) = 0.
                    qc_yp(2, ip, ic) = 0.
                    qc_gth(1, ip, ic) = 0.
                    qc_gth(2, ip, ic) = 0.
                enddo
                !HHY- Uncomment to save TE panel influences separately
                !c          DO KP = 1, 2
                !c            QC_GAMT(1,KP,IEL,IC) = 0.
                !c            QC_GAMT(2,KP,IEL,IC) = 0.
                !c            QC_SIGT(1,KP,IEL,IC) = 0.
                !c            QC_SIGT(2,KP,IEL,IC) = 0.
                !c          ENDDO
            enddo
        enddo
        !
        kc = ic1
        numic = ic2 - ic1 + 1
        !
        !---- Set AIC matrix for all singularity types
        do iel = iel1, iel2
            if (netype(iel)==0 .or. netype(iel)==1 .or. netype(iel)       &
                    & ==5 .or. netype(iel)==6 .or. netype(iel)==7) then
                !-------- wall,wake,source-only line and vortex wake sheets
                if (lxyjac) then
                    !--------- linearize QC w.r.t GAM,SIG,XP,YP,XC,YC
                    call panaicd(ipfrst(iel), iplast(iel), xp, yp, numic, xc(kc), &
                            & yc(kc), ipco(kc), ipcp(kc), ictype(kc), gam, sig, &
                            & ipx, qc(1, kc), qc_gam(1, 1, kc), qc_sig(1, 1, kc), &
                            & qc_xp(1, 1, kc), qc_yp(1, 1, kc), qc_xc(1, kc), &
                            & qc_yc(1, kc))
                    do ic = ic1, ic2
                        ipo = ipco(ic)
                        ipp = ipcp(ic)
                        if (ipo/=0 .and. ipp/=0) then
                            qc_xp(1, ipo, ic) = qc_xp(1, ipo, ic) + qc_xc(1, ic) * 0.5
                            qc_xp(2, ipo, ic) = qc_xp(2, ipo, ic) + qc_xc(2, ic) * 0.5
                            qc_xp(1, ipp, ic) = qc_xp(1, ipp, ic) + qc_xc(1, ic) * 0.5
                            qc_xp(2, ipp, ic) = qc_xp(2, ipp, ic) + qc_xc(2, ic) * 0.5
                            qc_yp(1, ipo, ic) = qc_yp(1, ipo, ic) + qc_yc(1, ic) * 0.5
                            qc_yp(2, ipo, ic) = qc_yp(2, ipo, ic) + qc_yc(2, ic) * 0.5
                            qc_yp(1, ipp, ic) = qc_yp(1, ipp, ic) + qc_yc(1, ic) * 0.5
                            qc_yp(2, ipp, ic) = qc_yp(2, ipp, ic) + qc_yc(2, ic) * 0.5
                        endif
                    enddo
                else
                    !--------- linearize QC w.r.t GAM,SIG only
                    call panaic(ipfrst(iel), iplast(iel), xp, yp, numic, xc(kc), &
                            & yc(kc), ipco(kc), ipcp(kc), ictype(kc), gam, sig, &
                            & ipx, qc(1, kc), qc_gam(1, 1, kc), qc_sig(1, 1, kc))
                endif
                !
            elseif (netype(iel)==2) then
                !------- line on axis
                do ip = ipfrst(iel), iplast(iel)
                    rcore(ip) = 0.01
                enddo
                call axlaic(ipfrst(iel), iplast(iel), xp, yp, numic, xc(kc), &
                        & yc(kc), ipco(kc), ipcp(kc), ictype(kc), gam, sig, &
                        & rcore, ipx, qc(1, kc), qc_gam(1, 1, kc), qc_sig(1, 1, kc)&
                        &)
                !
            elseif (netype(iel)==3) then
                !------- ring
                call linaic(ipfrst(iel), iplast(iel), xp, yp, numic, xc(kc), &
                        & yc(kc), gam, sig, ipx, qc(1, kc), qc_gam(1, 1, kc), &
                        & qc_sig(1, 1, kc))
                !
            elseif (netype(iel)==4) then
                !------- point
                call pntaic(ipfrst(iel), iplast(iel), xp, yp, numic, xc(kc), &
                        & yc(kc), gam, sig, ipx, qc(1, kc), qc_gam(1, 1, kc), &
                        & qc_sig(1, 1, kc))
            endif
            !
        enddo
        !
        !---- Accumulate velocity from wake vorticity, GTH
        !---- Save velocity AIC matrix (w/o TE panel influences) for wake vorticity
        do iel = iel1, iel2
            if (netype(iel)==0 .or. netype(iel)==7) then
                !
                ip1 = ipfrst(iel)
                ip2 = iplast(iel)
                do ic = ic1, ic2
                    do ip = ip1, ip2
                        qc_gth(1, ip, ic) = qc_gam(1, ip, ic)
                        qc_gth(2, ip, ic) = qc_gam(2, ip, ic)
                        qc(1, ic) = qc(1, ic) + qc_gth(1, ip, ic) * gth(ip)
                        qc(2, ic) = qc(2, ic) + qc_gth(2, ip, ic) * gth(ip)
                    enddo
                enddo
                !
            endif
        enddo
        !
        !---- Add contributions of TE panels...
        do iel = iel1, iel2
            if (.not.(ltpan(iel) .and. (netype(iel)==0 .or. netype(iel)==7&
                    &))) then
                !------- no TE panel... zero out panel strengths
                do k = 1, 2
                    xpt(k, iel) = 0.
                    ypt(k, iel) = 0.
                    gamt(k, iel) = 0.
                    sigt(k, iel) = 0.
                    do j = 1, 2
                        gamt_gam(k, j, iel) = 0.
                        gamt_sig(k, j, iel) = 0.
                        sigt_gam(k, j, iel) = 0.
                        sigt_sig(k, j, iel) = 0.
                    enddo
                enddo
                cycle
            endif
            !---- clear accumulating control-point velocities for each TE panel
            do ic = ic1, ic2
                do kp = 1, 2
                    qc_gamte(1, kp, ic) = 0.
                    qc_gamte(2, kp, ic) = 0.
                    qc_sigte(1, kp, ic) = 0.
                    qc_sigte(2, kp, ic) = 0.
                    qc_xpt(1, kp, ic) = 0.
                    qc_xpt(2, kp, ic) = 0.
                    qc_ypt(1, kp, ic) = 0.
                    qc_ypt(2, kp, ic) = 0.
                enddo
                qc_xc(1, ic) = 0.
                qc_xc(2, ic) = 0.
                qc_yc(1, ic) = 0.
                qc_yc(2, ic) = 0.
            enddo
            !
            !------ nodes on element spanned by TE panel
            ip1 = ipte1(iel)
            ip2 = ipte2(iel)
            !
            !------ set panel coordinates, strengths, and all Jacobians
            call tpanset(ip1, ip2, gam, sig, xp, yp, xpt(1, iel), xpt_xp(1, 1, iel), &
                    & ypt(1, iel), ypt_yp(1, 1, iel), gamt(1, iel), &
                    & gamt_gam(1, 1, iel), gamt_sig(1, 1, iel), &
                    & gamt_xp(1, 1, iel), gamt_yp(1, 1, iel), gamt_dx(1, iel), &
                    & gamt_dy(1, iel), sigt(1, iel), sigt_gam(1, 1, iel), &
                    & sigt_sig(1, 1, iel), sigt_xp(1, 1, iel), &
                    & sigt_yp(1, 1, iel), sigt_dx(1, iel), sigt_dy(1, iel))
            !
            !------ skip influence computation if TE panel has negligible length
            dxpt = xpt(2, iel) - xpt(1, iel)
            dypt = ypt(2, iel) - ypt(1, iel)
            dspt = sqrt(dxpt**2 + dypt**2)
            !
            splen = abs(sp(iplast(iel)) - sp(ipfrst(iel)))
            if (dspt<tfrac * splen) then
                write (*, *) 'Negligible TE panel ignored for element', iel
                cycle
            endif
            !
            !
            if (lxyjac) then
                !-------- linearize QC w.r.t GAM,SIG,XP,YP,XC,YC
                call panaicd(1, 2, xpt(1, iel), ypt(1, iel), numic, xc(kc), yc(kc), &
                        & izero, izero, ictype(kc), gamt(1, iel), sigt(1, iel), &
                        & 2, qc(1, kc), qc_gamte(1, 1, kc), qc_sigte(1, 1, kc), &
                        & qc_xpt(1, 1, kc), qc_ypt(1, 1, kc), qc_xc(1, kc), &
                        & qc_yc(1, kc))
                do ic = ic1, ic2
                    ipo = ipco(ic)
                    ipp = ipcp(ic)
                    if (ipo/=0 .and. ipp/=0) then
                        qc_xp(1, ipo, ic) = qc_xp(1, ipo, ic) + qc_xc(1, ic) * 0.5
                        qc_xp(2, ipo, ic) = qc_xp(2, ipo, ic) + qc_xc(2, ic) * 0.5
                        qc_xp(1, ipp, ic) = qc_xp(1, ipp, ic) + qc_xc(1, ic) * 0.5
                        qc_xp(2, ipp, ic) = qc_xp(2, ipp, ic) + qc_xc(2, ic) * 0.5
                        qc_yp(1, ipo, ic) = qc_yp(1, ipo, ic) + qc_yc(1, ic) * 0.5
                        qc_yp(2, ipo, ic) = qc_yp(2, ipo, ic) + qc_yc(2, ic) * 0.5
                        qc_yp(1, ipp, ic) = qc_yp(1, ipp, ic) + qc_yc(1, ic) * 0.5
                        qc_yp(2, ipp, ic) = qc_yp(2, ipp, ic) + qc_yc(2, ic) * 0.5
                    endif
                enddo
                !
            else
                !-------- linearize QC w.r.t GAM,SIG only
                call panaic(1, 2, xpt(1, iel), ypt(1, iel), numic, xc(kc), yc(kc), &
                        & izero, izero, ictype(kc), gamt(1, iel), sigt(1, iel), &
                        & 2, qc(1, kc), qc_gamte(1, 1, kc), qc_sigte(1, 1, kc))
            endif
            !
            !HHY- Uncomment to save TE panel influences separately
            !---- Save TE panel velocity influences
            !c        DO K = 1, 2
            !c          DO IC = IC1, IC2
            !c            QC_GAMT(K,1,IEL,IC) = QC_GAMTE(K,1,IC)
            !c            QC_GAMT(K,2,IEL,IC) = QC_GAMTE(K,2,IC)
            !c            QC_SIGT(K,1,IEL,IC) = QC_SIGTE(K,1,IC)
            !c            QC_SIGT(K,2,IEL,IC) = QC_SIGTE(K,2,IC)
            !c          END DO
            !c        END DO
            !
            do ic = ic1, ic2
                do kp = 1, 2
                    if (kp==1) then
                        ip = ip1
                        jpo = ip1 + 1
                        jpm = ip1
                    else
                        ip = ip2
                        jpo = ip2
                        jpm = ip2 - 1
                    endif
                    if (ip>0) then
                        do k = 1, 2
                            qc_gam(k, ip, ic) = qc_gam(k, ip, ic)                  &
                                    & + qc_gamte(k, 1, ic)               &
                                            & * gamt_gam(1, kp, iel)              &
                                    & + qc_gamte(k, 2, ic)               &
                                            & * gamt_gam(2, kp, iel)              &
                                    & + qc_sigte(k, 1, ic)               &
                                            & * sigt_gam(1, kp, iel)              &
                                    & + qc_sigte(k, 2, ic)               &
                                            & * sigt_gam(2, kp, iel)
                            qc_sig(k, ip, ic) = qc_sig(k, ip, ic)                  &
                                    & + qc_gamte(k, 1, ic)               &
                                            & * gamt_sig(1, kp, iel)              &
                                    & + qc_gamte(k, 2, ic)               &
                                            & * gamt_sig(2, kp, iel)              &
                                    & + qc_sigte(k, 1, ic)               &
                                            & * sigt_sig(1, kp, iel)              &
                                    & + qc_sigte(k, 2, ic)               &
                                            & * sigt_sig(2, kp, iel)
                            !
                            if (lxyjac) then
                                qc_xp(k, ip, ic) = qc_xp(k, ip, ic) + qc_xpt(k, 1, ic)&
                                        & * xpt_xp(1, kp, iel) + qc_xpt(k, 2, ic)           &
                                        & * xpt_xp(2, kp, iel) + qc_gamte(k, 1, ic)         &
                                        & * gamt_xp(1, kp, iel) + qc_gamte(k, 2, ic)        &
                                        & * gamt_xp(2, kp, iel) + qc_sigte(k, 1, ic)        &
                                        & * sigt_xp(1, kp, iel) + qc_sigte(k, 2, ic)        &
                                        & * sigt_xp(2, kp, iel)
                                qc_yp(k, ip, ic) = qc_yp(k, ip, ic) + qc_ypt(k, 1, ic)&
                                        & * ypt_yp(1, kp, iel) + qc_ypt(k, 2, ic)           &
                                        & * ypt_yp(2, kp, iel) + qc_gamte(k, 1, ic)         &
                                        & * gamt_yp(1, kp, iel) + qc_gamte(k, 2, ic)        &
                                        & * gamt_yp(2, kp, iel) + qc_sigte(k, 1, ic)        &
                                        & * sigt_yp(1, kp, iel) + qc_sigte(k, 2, ic)        &
                                        & * sigt_yp(2, kp, iel)
                                qc_xp(k, jpo, ic) = qc_xp(k, jpo, ic)               &
                                        & + qc_gamte(k, kp, ic) * gamt_dx(kp, iel)          &
                                        & + qc_sigte(k, kp, ic) * sigt_dx(kp, iel)
                                qc_yp(k, jpo, ic) = qc_yp(k, jpo, ic)               &
                                        & + qc_gamte(k, kp, ic) * gamt_dy(kp, iel)          &
                                        & + qc_sigte(k, kp, ic) * sigt_dy(kp, iel)
                                qc_xp(k, jpm, ic) = qc_xp(k, jpm, ic)               &
                                        & - qc_gamte(k, kp, ic) * gamt_dx(kp, iel)          &
                                        & - qc_sigte(k, kp, ic) * sigt_dx(kp, iel)
                                qc_yp(k, jpm, ic) = qc_yp(k, jpm, ic)               &
                                        & - qc_gamte(k, kp, ic) * gamt_dy(kp, iel)          &
                                        & - qc_sigte(k, kp, ic) * sigt_dy(kp, iel)
                            endif
                        enddo   ! K loop
                    endif
                enddo  ! KP loop
            enddo ! IC loop
            !
            !---- Special treatment for vortex wake last GTH point to include TE panel
            if (ltpan(iel) .and. netype(iel)==7) then
                ip = iplast(iel)
                do ic = ic1, ic2
                    !---- Remove velocity w/o TE panel
                    qc(1, ic) = qc(1, ic) - qc_gth(1, ip, ic) * gth(ip)
                    qc(2, ic) = qc(2, ic) - qc_gth(2, ip, ic) * gth(ip)
                    !---- Save influence with TE panel
                    qc_gth(1, ip, ic) = qc_gam(1, ip, ic)
                    qc_gth(2, ip, ic) = qc_gam(2, ip, ic)
                    !---- Add back in velocity with TE panel
                    qc(1, ic) = qc(1, ic) + qc_gth(1, ip, ic) * gth(ip)
                    qc(2, ic) = qc(2, ic) + qc_gth(2, ip, ic) * gth(ip)
                enddo
            endif
            !
        enddo
        !
        !---- add on freestream
        do ic = ic1, ic2
            qc(1, ic) = qc(1, ic) + qinf
        enddo
        !
    end subroutine qaic1
    !*==TPANSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QAIC1



    subroutine tpanset(ipte1, ipte2, gam, sig, xp, yp, xpt, xpt_xp, ypt, &
            & ypt_yp, gamt, gamt_gam, gamt_sig, gamt_xp, gamt_yp, &
            & gamt_dx, gamt_dy, sigt, sigt_gam, sigt_sig, sigt_xp, &
            & sigt_yp, sigt_dx, sigt_dy)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ipte1, ipte2
        real, dimension(*) :: gam, sig, xp, yp
        real, dimension(2) :: gamt, gamt_dx, gamt_dy, sigt, &
                & sigt_dx, sigt_dy, xpt, ypt
        real, dimension(2, 2) :: gamt_gam, gamt_sig, gamt_xp, &
                & gamt_yp, sigt_gam, sigt_sig, &
                & sigt_xp, sigt_yp, xpt_xp, ypt_yp
        !
        ! Local variables
        !
        real, dimension(2) :: anx, anxt_xp, anxt_yp, anx_dx, &
                & anx_dy, any, anyt_xp, anyt_yp, &
                & any_dx, any_dy, crs, crs_dx, crs_dy, &
                & dot, dot_dx, dot_dy
        real :: anxt, anxt_dxt, anxt_dyt, anyaxi, anyt, anyt_dxt, &
                & anyt_dyt, ds, dsti, dstsq, dx, dxt, dy, dyt
        real, dimension(2, 2) :: crs_xp, crs_yp, dot_xp, dot_yp
        integer :: ip, ip1, ip2, ipopp, j, jpm, jpo, k, kopp
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------
        !     Sets geometry, strengths, and associated Jacobians
        !     of TE panel on element IEL
        !----------------------------------------------------------
        !
        !
        !---- local setup arrays
        !
        !
        do k = 1, 2
            do j = 1, 2
                xpt_xp(k, j) = 0.
                ypt_yp(k, j) = 0.
                gamt_gam(k, j) = 0.
                gamt_sig(k, j) = 0.
                sigt_gam(k, j) = 0.
                sigt_sig(k, j) = 0.
                gamt_xp(k, j) = 0.
                gamt_yp(k, j) = 0.
                sigt_xp(k, j) = 0.
                sigt_yp(k, j) = 0.
            enddo
            anx_dx(k) = 0.
            any_dx(k) = 0.
            anx_dy(k) = 0.
            any_dy(k) = 0.
        enddo
        !
        !---- geometry nodes spanned by TE panel
        !-    (if node index = 0, then that TE panel node lies on axis)
        ip1 = ipte1
        ip2 = ipte2
        !
        !---- set node x,y of TE panel, and dx,dy of adjacent panels
        !
        do k = 1, 2
            if (k==1) then
                ip = ip1
                ipopp = ip2
                kopp = 2
                jpo = ip + 1
                jpm = ip
                anyaxi = 1.0
            else
                ip = ip2
                ipopp = ip1
                kopp = 1
                jpo = ip
                jpm = ip - 1
                anyaxi = -1.0
            endif
            !
            if (ip>0) then
                !------- usual TE panel
                xpt(k) = xp(ip)
                ypt(k) = yp(ip)
                xpt_xp(k, k) = 1.0
                ypt_yp(k, k) = 1.0
                dx = xp(jpo) - xp(jpm)
                dy = yp(jpo) - yp(jpm)
                ds = sqrt(dx**2 + dy**2)
                anx(k) = dy / ds
                any(k) = -dx / ds
                anx_dx(k) = -dx * dy / ds**3
                anx_dy(k) = dx * dx / ds**3
                any_dx(k) = -dy * dy / ds**3
                any_dy(k) = dy * dx / ds**3
            else
                !------- node IP is on axis, at same x location as node IPOPP
                xpt(k) = xp(ipopp)
                ypt(k) = 0.
                xpt_xp(k, kopp) = 1.0
                anx(k) = 0.
                any(k) = anyaxi
            endif
        enddo
        !
        !---- now set TE panel normal vector ANXT,ANYT
        dxt = xpt(1) - xpt(2)
        dyt = ypt(1) - ypt(2)
        dstsq = dxt**2 + dyt**2
        if (dstsq==0.0) then
            dsti = 1.0
        else
            dsti = 1.0 / sqrt(dstsq)
        endif
        anxt = dyt * dsti
        anyt = -dxt * dsti
        anxt_dxt = -dxt * dyt * dsti**3
        anxt_dyt = dxt * dxt * dsti**3
        anyt_dxt = -dyt * dyt * dsti**3
        anyt_dyt = dyt * dxt * dsti**3
        !
        do j = 1, 2
            anxt_xp(j) = anxt_dxt * (xpt_xp(1, j) - xpt_xp(2, j))
            anxt_yp(j) = anxt_dyt * (ypt_yp(1, j) - ypt_yp(2, j))
            !
            anyt_xp(j) = anyt_dxt * (xpt_xp(1, j) - xpt_xp(2, j))
            anyt_yp(j) = anyt_dyt * (ypt_yp(1, j) - ypt_yp(2, j))
        enddo
        !
        !
        !---- form dot,cross products with adjacent-panel normal vectors
        do k = 1, 2
            dot(k) = anxt * anx(k) + anyt * any(k)
            crs(k) = anxt * any(k) - anyt * anx(k)
            dot_dx(k) = anxt * anx_dx(k) + anyt * any_dx(k)
            dot_dy(k) = anxt * anx_dy(k) + anyt * any_dy(k)
            crs_dx(k) = anxt * any_dx(k) - anyt * anx_dx(k)
            crs_dy(k) = anxt * any_dy(k) - anyt * anx_dy(k)
            do j = 1, 2
                dot_xp(k, j) = anxt_xp(j) * anx(k) + anyt_xp(j) * any(k)
                dot_yp(k, j) = anxt_yp(j) * anx(k) + anyt_yp(j) * any(k)
                crs_xp(k, j) = anxt_xp(j) * any(k) - anyt_xp(j) * anx(k)
                crs_yp(k, j) = anxt_yp(j) * any(k) - anyt_yp(j) * anx(k)
            enddo
        enddo
        !
        !
        do k = 1, 2
            if (k==1) then
                ip = ip1
                ipopp = ip2
                kopp = 2
            else
                ip = ip2
                ipopp = ip1
                kopp = 1
            endif
            !
            if (ip>0) then
                gamt(k) = dot(k) * gam(ip) - crs(k) * sig(ip)
                sigt(k) = crs(k) * gam(ip) + dot(k) * sig(ip)
                gamt_sig(k, k) = -crs(k)
                gamt_gam(k, k) = dot(k)
                sigt_sig(k, k) = dot(k)
                sigt_gam(k, k) = crs(k)
                gamt_dx(k) = dot_dx(k) * gam(ip) - crs_dx(k) * sig(ip)
                gamt_dy(k) = dot_dy(k) * gam(ip) - crs_dy(k) * sig(ip)
                sigt_dx(k) = crs_dx(k) * gam(ip) + dot_dx(k) * sig(ip)
                sigt_dy(k) = crs_dy(k) * gam(ip) + dot_dy(k) * sig(ip)
                do j = 1, 2
                    gamt_xp(k, j) = dot_xp(k, j) * gam(ip) - crs_xp(k, j) * sig(ip)
                    gamt_yp(k, j) = dot_yp(k, j) * gam(ip) - crs_yp(k, j) * sig(ip)
                    sigt_xp(k, j) = crs_xp(k, j) * gam(ip) + dot_xp(k, j) * sig(ip)
                    sigt_yp(k, j) = crs_yp(k, j) * gam(ip) + dot_yp(k, j) * sig(ip)
                enddo
            else
                gamt(k) = 0.
                sigt(k) = crs(kopp) * gam(ipopp) + dot(kopp) * sig(ipopp)
                sigt_sig(k, kopp) = dot(kopp)
                sigt_gam(k, kopp) = crs(kopp)
                sigt_dx(k) = crs_dx(k) * gam(ipopp) + dot_dx(k) * sig(ipopp)
                sigt_dy(k) = crs_dy(k) * gam(ipopp) + dot_dy(k) * sig(ipopp)
                do j = 1, 2
                    sigt_xp(k, j) = crs_xp(k, j) * gam(ipopp) + dot_xp(k, j)      &
                            & * sig(ipopp)
                    sigt_yp(k, j) = crs_yp(k, j) * gam(ipopp) + dot_yp(k, j)      &
                            & * sig(ipopp)
                enddo
            endif
        enddo
        !
    end subroutine tpanset
    !*==QFCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! TPANSET


    subroutine qfcalc(if1, if2, xf, yf, ipfo, ipfp, iftype, qf, qf_gam, qf_sig, &
            & qf_gth)
        use i_dfdc
        use m_aic, only : axlaic, linaic, panaic, pntaic
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: if1, if2
        integer, dimension(*) :: iftype, ipfo, ipfp
        real, dimension(2, *) :: qf
        real, dimension(2, ipx, *) :: qf_gam, qf_gth, qf_sig
        real, dimension(*) :: xf, yf
        !
        ! Local variables
        !
        integer :: i, iel, if, ip, ip1, ip2, k, kf, kp, n, nf
        real, dimension(2, 2, ipx) :: qf_gamte, qf_sigte
        real, dimension(ipx) :: rcore
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------
        !     Calculates velocity components QF(.) at field points
        !     IF1...IF2 in coordinate arrays XF,YF w.r.t aero quantities
        !     at all nodes on all elements (1..NEL)
        !---------------------------------------------------------------
        !     Input:  XF,YF(i)  field points
        !             IPFO(.)   first corner point index (0 if not on panel)
        !             IPFP(.)   second corner point index (0 if not on panel)
        !             IFTYPE(.) type for field point (0 if not on panel)
        !
        !     Assumes:
        !             GAM(j)    panel vorticity  (or line circ. , or point doublet)
        !             SIG(j)    panel source     (or line source, or point source )
        !             GTH(j)    vortex wake panel vorticity
        !             QINF      freestream speed
        !
        !             GAMT(j)   TE panel vorticity
        !             SIGT(j)   TE panel source
        !
        !     Output: QF(.i)         velocity at XC(i),YC(i)
        !             QF_GAM(.ji)    dQC(.i)/dGAM(j)
        !             QF_SIG(.ji)
        !             QF_GTH(.ji)    dQC(.i)/dGTH(j)
        !---------------------------------------------------------------
        !    Note:    This routine assumes that TE panel strengths are set
        !             before entry.
        !---------------------------------------------------------------
        !HHY- Use these to save TE panel influences separately
        !c      SUBROUTINE QFCALC(IF1,IF2, XF,YF, IPFO,IPFP,IFTYPE,
        !c     &                  QF,QF_GAM,QF_SIG, QF_GAMT,QF_SIGT)
        !c      DIMENSION          QF_GAMT(2,2,NEX,*), QF_SIGT(2,2,NEX,*)
        !--------------------------------------------------------
        !
        !
        do k = 1, 2
            do if = if1, if2
                qf(k, if) = 0.
                do ip = 1, nptot
                    qf_gam(k, ip, if) = 0.
                    qf_sig(k, ip, if) = 0.
                    qf_gth(k, ip, if) = 0.
                enddo
                !HHY- Uncomment to save TE panel influences separately
                !c          DO IEL = 1, NEL
                !c            QF_GAMT(K,1,IEL,IF) = 0.
                !c            QF_SIGT(K,2,IEL,IF) = 0.
                !c          ENDDO
            enddo
        enddo
        !
        nf = if2 - if1 + 1
        kf = if1
        !
        !---- Set contribution of all vortex/source lines, rings, and points
        do iel = 1, nel
            if (netype(iel)==0 .or. netype(iel)==1 .or. netype(iel)       &
                    & ==5 .or. netype(iel)==6 .or. netype(iel)==7) then
                !------- sheets
                call panaic(ipfrst(iel), iplast(iel), xp, yp, nf, xf(kf), yf(kf), &
                        & ipfo(kf), ipfp(kf), iftype(kf), gam, sig, ipx, &
                        & qf(1, kf), qf_gam(1, 1, kf), qf_sig(1, 1, kf))
                !
            elseif (netype(iel)==2) then
                !------- line on axis
                do ip = ipfrst(iel), iplast(iel)
                    rcore(ip) = 0.01
                enddo
                call axlaic(ipfrst(iel), iplast(iel), xp, yp, nf, xf(kf), yf(kf), &
                        & ipfo(kf), ipfp(kf), iftype(kf), gam, sig, rcore, ipx, &
                        & qf(1, kf), qf_gam(1, 1, kf), qf_sig(1, 1, kf))
                !
            elseif (netype(iel)==3) then
                !------- lines, rings
                call linaic(ipfrst(iel), iplast(iel), xp, yp, nf, xf(kf), yf(kf), &
                        & gam, sig, ipx, qf(1, kf), qf_gam(1, 1, kf), &
                        & qf_sig(1, 1, kf))
                !
            elseif (netype(iel)==4) then
                !------- point
                i = ipfrst(iel)
                n = 1
                call pntaic(ipfrst(iel), iplast(iel), xp, yp, nf, xf(kf), yf(kf), &
                        & gam, sig, ipx, qf(1, kf), qf_gam(1, 1, kf), &
                        & qf_sig(1, 1, kf))
            endif
        enddo
        !
        !
        !
        !---- Accumulate velocity from wake vorticity, GTH
        !---- Save velocity AIC matrix (w/o TE panel influences) for wake vorticity
        do iel = 1, nel
            if (netype(iel)==0 .or. netype(iel)==7) then
                !
                ip1 = ipfrst(iel)
                ip2 = iplast(iel)
                do if = if1, if2
                    do ip = ip1, ip2
                        qf_gth(1, ip, if) = qf_gam(1, ip, if)
                        qf_gth(2, ip, if) = qf_gam(2, ip, if)
                        qf(1, if) = qf(1, if) + qf_gth(1, ip, if) * gth(ip)
                        qf(2, if) = qf(2, if) + qf_gth(2, ip, if) * gth(ip)
                    enddo
                enddo
                !
            endif
        enddo

        !
        !---- Add contribution of all TE panels
        do iel = 1, nel
            if (ltpan(iel) .and. (netype(iel)==0 .or. netype(iel)==7))   &
                    & then
                !---- clear field-point velocities for each TE panel
                do if = if1, if2
                    do kp = 1, 2
                        qf_gamte(1, kp, if) = 0.
                        qf_gamte(2, kp, if) = 0.
                        qf_sigte(1, kp, if) = 0.
                        qf_sigte(2, kp, if) = 0.
                    enddo
                enddo
                !
                call panaic(1, 2, xpt(1, iel), ypt(1, iel), nf, xf(kf), yf(kf), &
                        & izero, izero, iftype(kf), gamt(1, iel), sigt(1, iel), &
                        & 2, qf(1, kf), qf_gamte(1, 1, kf), qf_sigte(1, 1, kf))
                !
                !HHY- Uncomment to save TE panel influences separately
                !---- Save TE velocity influences
                !c          DO IF = IF1, IF2
                !c            DO K = 1, 2
                !c              QF_GAMT(K,1,IEL,IF) = QF_GAMTE(K,1,IF)
                !c              QF_GAMT(K,2,IEL,IF) = QF_GAMTE(K,2,IF)
                !c              QF_SIGT(K,1,IEL,IF) = QF_SIGTE(K,1,IF)
                !c              QF_SIGT(K,2,IEL,IF) = QF_SIGTE(K,2,IF)
                !c            END DO
                !c          END DO
                !
                !---- Add TE influences to panel nodes (assumes TPANSET has been called)
                do if = if1, if2
                    do kp = 1, 2
                        if (kp==1) then
                            ip = ipte1(iel)
                        else
                            ip = ipte2(iel)
                        endif
                        if (ip>0) then
                            do k = 1, 2
                                qf_gam(k, ip, if) = qf_gam(k, ip, if)               &
                                        & + qf_gamte(k, 1, if) * gamt_gam(1, kp, iel)        &
                                        & + qf_gamte(k, 2, if) * gamt_gam(2, kp, iel)        &
                                        & + qf_sigte(k, 1, if) * sigt_gam(1, kp, iel)        &
                                        & + qf_sigte(k, 2, if) * sigt_gam(2, kp, iel)
                                qf_sig(k, ip, if) = qf_sig(k, ip, if)               &
                                        & + qf_gamte(k, 1, if) * gamt_sig(1, kp, iel)        &
                                        & + qf_gamte(k, 2, if) * gamt_sig(2, kp, iel)        &
                                        & + qf_sigte(k, 1, if) * sigt_sig(1, kp, iel)        &
                                        & + qf_sigte(k, 2, if) * sigt_sig(2, kp, iel)
                            enddo    ! K loop
                        endif
                    enddo   ! KP loop
                enddo  ! IF loop
                !
            endif
            !
            !---- Special treatment for vortex wake last GTH point to include TE panel
            if (ltpan(iel) .and. netype(iel)==7) then
                ip = iplast(iel)
                do if = if1, if2
                    !---- Remove velocity w/o TE panel
                    !c            QF(1,IF) = QF(1,IF) - QF_GTH(1,IP,IF)*GTH(IP)
                    !c            QF(2,IF) = QF(2,IF) - QF_GTH(2,IP,IF)*GTH(IP)
                    !---- Save influence with TE panel
                    qf_gth(1, ip, if) = qf_gam(1, ip, if)
                    qf_gth(2, ip, if) = qf_gam(2, ip, if)
                    !---- Add back in velocity with TE panel
                    !c            QF(1,IF) = QF(1,IF) + QF_GTH(1,IP,IF)*GTH(IP)
                    !c            QF(2,IF) = QF(2,IF) + QF_GTH(2,IP,IF)*GTH(IP)
                enddo
            endif

        enddo
        !
        !---- Add freestream
        !      DO IF = 1, NF
        !        QF(1,IF) = QF(1,IF) + QINF
        !      ENDDO
        !
    end subroutine qfcalc
end module m_qaic
