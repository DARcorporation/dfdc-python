module m_sgutil
    implicit none
contains
    !*==SGCURV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020

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

    subroutine sgcurv(lquery, lcplot, nb, xb, xpb, yb, ypb, sb, ngx, ng, sg, &
            & npts, cvex, smof, fsl, fsr, nref, sref1, sref2, crrat)
        use m_spline, only : trisol, segspl, seval, curv
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: ix = 1601, nrefx = 20
        !
        ! Dummy arguments
        !
        real :: cvex, fsl, fsr, smof
        logical :: lcplot, lquery
        integer :: nb, ng, ngx, npts, nref
        real, dimension(nref) :: crrat, sref1, sref2
        real, dimension(nb) :: sb, xb, xpb, yb, ypb
        real, dimension(ngx) :: sg
        !
        ! Local variables
        !
        real, dimension(ix) :: aa, bb, cc, cv, cv0, cvint, cvs, &
                & s, snew
        real :: curvk, curvl, cvasq, cvavg, cvle, cvm, cvmax, &
                & cvo, cvte, dle, dsa, dsavg, dsb, dsle, dsm, &
                & dsmax, dsmin, dsp, dste, dte, frac, gcv, sbtot, &
                & sgi, sm, so, stot
        real, dimension(nrefx) :: cvr
        integer :: i, ib, ig, ipass, iref, k, kk, kmult, n
        integer, save :: npass
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------------
        !     Sets spacing array SG based on curvature and a few parameters
        !
        !  Input:  LQUERY   if T, requests user input for parameters before paneling
        !          LCPLOT   if T, spacing parameters are plotted
        !          NB       number of coordinates
        !          XB(i)    x coordinates
        !          XPB(i)   dXB/dSB spline derivative
        !          YB(i)    y coordinates
        !          YPB(i)   dYB/dSB spline derivative
        !          SB(i)    spline parameter (typically arc length)
        !          NGX      max number of points which can be generated
        !
        !  Input, also Output if interactively modified:
        !          NG       number of points actually generated
        !          SG(k)    fractional spacing array  0.0 .. 1.0
        !          NPTS     number of points to be generated
        !          CVEX     curvature-exponent parameter
        !          SMOF     smoothing-length factor
        !          FSL      fractional spacing at left  endpoint
        !          FSR      fractional spacing at right endpoint
        !          NREF     number of local-refinement regions
        !          SREF1(r) first point of refinement region
        !          SREF2(r) last  point of refinement region
        !          CRRAT(r) fractional spacing over refinement region
        !
        !
        !---------------------------------------------------------------------
        !
        !
        !
        !
        !
        !---- max number of point-setting passes
        data npass/1/
        !
        if (nb>ix) stop 'SGCURV: Local array overflow on IX'
        if (nref>nrefx) stop 'SGCURV: Local array overflow on NREFX'
        !
        !
        sbtot = sb(nb) - sb(1)
        if (sbtot<=0.0) then
            write (*, *) '? SGCURV: Zero element arc length'
            return
        endif
        !
        !---- if no interaction or no spacing is present, go compute it first
        if (.not.lquery .or. ng==0) then
        endif
        !================================================================
        !==== Set surface panel node distribution based on curvature
        !
        ng = npts
        !
        if (ng<2) then
            write (*, *) '? SGCURV: Must specify at least two panel nodes'
            return
        endif
        !
        kmult = min(8, (ix - 1) / nb)
        !
        !---- set up working coordinate array
        i = 1
        ib = 1
        s(i) = sb(ib)
        do ib = 1, nb - 1
            dsb = sb(ib + 1) - sb(ib)
            !
            if (dsb>0.0) then
                kk = kmult
            else
                kk = 1
            endif
            !
            do k = 1, kk
                frac = float(k) / float(kk)
                i = i + 1
                if (i>ix) stop 'SGCURV:  Array overflow on IX'
                s(i) = sb(ib + 1) * frac + sb(ib) * (1.0 - frac)
            enddo
        enddo
        !
        !---- number of points in working array
        n = i
        !
        stot = s(n) - s(1)
        do ipass = 1, npass
            !
            dsavg = stot / float(ng - 1)
            !
            dsle = fsl * dsavg
            dste = fsr * dsavg
            !
            !---- set up curvature array (nondimensionalized with perimeter)
            do i = 1, n
                cv(i) = curv(s(i), xb, xpb, yb, ypb, sb, nb) * stot
                cv(i) = abs(cv(i)) + 1.0
            enddo
            !
            !---- reset curvature at corner nodes from adjacent nodes
            do i = 2, n - 2
                if (s(i)==s(i + 1)) then
                    cv(i) = 0.5 * (cv(i - 1) + cv(i + 2))
                    cv(i + 1) = cv(i)
                endif
            enddo
            !
            !---- raise curvature to power sqrt(CVEX), limiting to prevent overflows
            do i = 1, n
                gcv = min(0.5 * cvex * log(cv(i)), 12.0)
                cv(i) = exp(gcv)
            enddo
            !
            !---- set max and average curvature
            cvmax = 0.01
            cvavg = 0.
            do i = 1, n - 1
                cvasq = 0.5 * (cv(i)**2 + cv(i + 1)**2)
                cvmax = max(cvmax, sqrt(cvasq))
                cvavg = cvavg + cvasq * (s(i + 1) - s(i))
            enddo
            cvavg = sqrt(cvavg / (s(n) - s(1)))
            !
            cvavg = max(cvavg, 1.0)
            cvmax = max(cvmax, 1.0)
            !
            !---- set artificial curvature at ends to get approximate ds/c there
            cvle = cvavg * dsavg / dsle
            cvte = cvavg * dsavg / dste
            cv(1) = cvle
            cv(n) = cvte
            !
            do i = 1, n
                cv0(i) = cv(i)
            enddo
            !
            !---- set curvature smoothing length
            curvl = 0.008 * stot
            !
            !---- set up implicit system for smoothed curvature array CV
            curvk = curvl**2 * smof
            aa(1) = 1.0
            cc(1) = 0.
            do i = 2, n - 1
                dsm = s(i) - s(i - 1)
                dsp = s(i + 1) - s(i)
                dsa = 0.5 * (dsm + dsp)
                if (dsm==0.0 .or. dsp==0.0) then
                    bb(i) = 0.0
                    aa(i) = 1.0
                    cc(i) = 0.0
                else
                    bb(i) = -curvk / (dsm * dsa)
                    aa(i) = curvk / (dsm * dsa) + curvk / (dsp * dsa) + 1.0
                    cc(i) = -curvk / (dsp * dsa)
                endif
            enddo
            aa(n) = 1.0
            bb(n) = 0.
            !
            !---- set artificial curvature at the bunching points
            do iref = 1, nref
                if (crrat(iref)>0.0) then
                    cvr(iref) = cvavg / max(crrat(iref), 0.0001)
                else
                    cvr(iref) = 0.0
                endif
                !
                do i = 2, n - 1
                    if (crrat(iref)>0.0) then
                        !--------- set modified curvature if point is in refinement area
                        sgi = (s(i) - s(1)) / stot
                        !c            BB(I) = 0.
                        !c            AA(I) = 1.0
                        !c            CC(I) = 0.
                        if (sgi>sref1(iref) .and. sgi<sref2(iref)) cv(i)    &
                                & = cvr(iref)
                    endif
                enddo
            enddo
            !
            !---- calculate smoothed curvature array
            call trisol(aa, bb, cc, cv, n)
            !
            !---- spline curvature array
            call segspl(cv, cvs, s, n)
            !
            !---- Integrate exponentiated curvature with arc length
            cvint(1) = 0.
            do i = 2, n
                sm = s(i - 1)
                so = s(i)
                cvm = seval(sm, cv, cvs, s, n)
                cvo = seval(so, cv, cvs, s, n)
                cvm = max(cvm, 1.0)
                cvo = max(cvo, 1.0)
                !
                cvint(i) = cvint(i - 1) + 0.5 * (cvo + cvm) * (so - sm)
            enddo
            !
            do i = 1, n
                cvint(i) = cvint(i) + s(i)
            enddo
            !
            !---- Calculate normalized surface spacing distribution arrays
            call setcrv(snew, ng, s, cvint, n)
            !
            do ig = 1, ng
                sg(ig) = (snew(ig) - snew(1)) / (snew(ng) - snew(1))
            enddo
            !
            !---- calculate actual obtained grid spacings
            dle = abs(snew(2) - snew(1))
            dte = abs(snew(ng) - snew(ng - 1))
            dsmin = dle
            dsmax = 0.0
            do ig = 2, ng
                dsmin = min(dsmin, abs(snew(ig) - snew(ig - 1)))
                dsmax = max(dsmax, abs(snew(ig) - snew(ig - 1)))
            enddo
            !
        enddo
    end subroutine sgcurv
    !*==SETCRV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SGCURV



    subroutine setcrv(snew, npts, sb, hh, n)
        use m_spline, only : seval, splina
        !------------------------------------------------------------------
        !     Sets spacing array SNEW(i) based on surface curvature.
        !
        !     Input: NPTS    number of points defining surface
        !            SB(.)   surface arc length
        !            HH(.)   integrated spacing function
        !            N       number of spacing array points to be generated
        !
        !     Output: SNEW(.)  spacing array to be generated
        !             LFUDGL   T if curvature was augmented al left endpoint
        !             LFUDGT   T if curvature was augmented al left endpoint
        !------------------------------------------------------------------
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: nmax = 1601
        !
        ! Dummy arguments
        !
        integer :: n, npts
        real, dimension(n) :: hh, sb
        real, dimension(npts) :: snew
        !
        ! Local variables
        !
        real :: dh, hhi, rn
        integer :: i
        real, dimension(nmax) :: sbp
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !
        if (n>nmax) stop 'SETCRV: array overflow'
        !
        rn = float(npts - 1)
        !
        !
        !---- spline  HH(i)  =  h(s)  =  s + c(s)
        call splina(hh, sbp, sb, n)
        !
        !---- invert derivative from dh/ds to ds/dh
        do i = 1, n
            sbp(i) = 1.0 / sbp(i)
        enddo
        !
        !---- now,  SB(i),SBP(i)  =  s(h), s'(h)
        !
        !
        dh = hh(n) / rn
        !
        snew(1) = 0.0
        do i = 2, npts - 1
            hhi = float(i - 1) * dh
            snew(i) = seval(hhi, sb, sbp, hh, n)
        enddo
        snew(npts) = sb(n)
        !
    end subroutine setcrv
end module m_sgutil
