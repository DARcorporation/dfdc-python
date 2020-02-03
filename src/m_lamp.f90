module m_lamp
    implicit none
contains
    !*==LAMP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! AXELL
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

    subroutine lamp(x1, r1, x2, r2, xf, rf, ug1, vg1, ug2, vg2, us1, vs1, us2, vs2)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: nromx = 10
        !
        ! Dummy arguments
        !
        real :: r1, r2, rf, ug1, ug2, us1, us2, vg1, vg2, vs1, &
                & vs2, x1, x2, xf
        !
        ! Local variables
        !
        real :: dels, delsq, dt, err, errug1, errug2, errus1, &
                & errus2, errvg1, errvg2, errvs1, errvs2, refl, rt, &
                & t, tb, ugt, ust, vgt, vst, w, xt
        integer :: irom, it, krom, nt
        real, save :: romtol
        real, dimension(nromx) :: ug1i, ug2i, us1i, us2i, vg1i, &
                & vg2i, vs1i, vs2i
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------------
        !
        !     Computes the velocities at XF,RF induced by a vortex + source
        !     sheet "lampshade"  extending from X1,R1 to X2,R2.
        !
        !      Vortex sheet density = gamma  (positive counterclockwise)
        !      Source sheet density = sigma
        !
        !     Both densities are assumed to be linear in meridional
        !     arc length over the panel.
        !
        !
        !  Input:
        !  ------
        !     X1,R1   x,r of one   lampshade edge
        !     X2,R2   x,r of other lampshade edge
        !     XF,RF   x,r of the field point
        !
        !  Output:
        !  -------
        !     UG1,VG1   x,r velocities at XF,RF for unit gamma at X1,R1
        !     UG2,VG2   x,r velocities at XF,RF for unit gamma at X2,R2
        !     US1,VS1   x,r velocities at XF,RF for unit sigma at X1,R1
        !     US2,VS2   x,r velocities at XF,RF for unit sigma at X2,R2
        !
        !     Total x,r velocities for given endpoint sheet densities are:
        !
        !       U = UG1*gamma1 + UG2*gamma2 + US1*sigma1 + US2*sigma2
        !       V = VG1*gamma1 + VG2*gamma2 + VS1*sigma1 + VS2*sigma2
        !
        !-----------------------------------------------------------------------
        !
        !     Integrates point vortex/source ring velocities over
        !     the lampshade using a Romberg sequence applied to the
        !     simple midpoint rule:
        !
        !        |       *       |  ->  I1  -   I21  -   I321  -   I4321
        !                                   /        /         /     .
        !        |   *       *   |  ->  I2  -   I32  -   I432        .
        !                                   /        /     .       8th order
        !        | *   *   *   * |  ->  I3  -   I43        .
        !                                   /    .       6th order
        !        |* * * * * * * *|  ->  I4       .
        !                .               .      4th order
        !                .               .
        !                               2nd order
        !               etc
        !
        !  The first column I1,I2... are the 2nd-order integral approximations
        !  computed with the stock midpoint rule.  The subsequent columns
        !  are the Richardson extrapolations to higher-order accuracy.
        !
        !
        !  Algorithm for four Romberg stages:
        !
        !    I21 = (4*I2 - I1) / 3           | extrapolation of 2nd-order results
        !    I32 = (4*I3 - I2) / 3           |
        !    I43 = (4*I4 - I3) / 3           |
        !
        !    I321 = (16*I32 - I21) / 15      | extrapolation of 4th-order results
        !    I432 = (16*I43 - I32) / 15      |
        !
        !    I4321 = (64*I432 - I321) / 63   | extrapolation of 6th-order results
        !
        !
        !  Quantities stored in the NROMX arrays after each IROM stage:
        !
        !    IROM  =   1      2       3       4
        !            ----   ----    -----   ------
        !    U(1)  =  I1     I21     I321    I4321
        !    U(2)  =         I2      I32     I432
        !    U(3)  =                 I3      I43
        !    U(4)  =                         I4
        !
        !-----------------------------------------------------------------------
        !
        !---- NROMX = max number of Romberg stages
        !-      This limits the integration resolution and cost.
        !-      The influence of this panel on a field point
        !-      which is closer than  O(panel_length)/2**NROMX
        !-      will not be accurately represented.
        !
        !
        !
        !
        !---- Romberg convergence tolerance
        data romtol/1.0e-6/
        !cc      DATA ROMTOL / 1.0E-12 /
        !
        !---- reference length for convergence tolerance
        refl = 0.5 * (r1 + r2)
        !
        !---- evaluate integrals over  0..t..1  on increasingly fine grids
        do irom = 1, nromx
            nt = 2**irom / 2
            !
            ug1i(irom) = 0.
            vg1i(irom) = 0.
            ug2i(irom) = 0.
            vg2i(irom) = 0.
            us1i(irom) = 0.
            vs1i(irom) = 0.
            us2i(irom) = 0.
            vs2i(irom) = 0.
            !
            !------ visit the midpoints of each of the NT intervals
            do it = 1, nt
                t = (float(it) - 0.5) / float(nt)
                tb = 1.0 - t
                !
                dt = 1.0 / float(nt)
                !
                xt = x1 * tb + x2 * t
                rt = r1 * tb + r2 * t
                !
                !-------- get induced velocities for vortex,source ring at XT,RT
                call ring(xt, rt, xf, rf, ugt, vgt, ust, vst)
                !
                !-------- accumulate the separate unit-gamma, unit-sigma integrals
                ug1i(irom) = ug1i(irom) + dt * ugt * tb
                vg1i(irom) = vg1i(irom) + dt * vgt * tb
                ug2i(irom) = ug2i(irom) + dt * ugt * t
                vg2i(irom) = vg2i(irom) + dt * vgt * t
                us1i(irom) = us1i(irom) + dt * ust * tb
                vs1i(irom) = vs1i(irom) + dt * vst * tb
                us2i(irom) = us2i(irom) + dt * ust * t
                vs2i(irom) = vs2i(irom) + dt * vst * t
            enddo
            !
            !------ Romberg sequence using all previous grid results
            do krom = irom, 2, -1
                !-------- weight needed to cancel lowest-order error terms in KROM level
                w = 2.0**(2 * (irom - krom + 1))
                !
                !-------- put Richardson extrapolation for KROM level into KROM-1 level
                ug1i(krom - 1) = (w * ug1i(krom) - ug1i(krom - 1)) / (w - 1.0)
                vg1i(krom - 1) = (w * vg1i(krom) - vg1i(krom - 1)) / (w - 1.0)
                ug2i(krom - 1) = (w * ug2i(krom) - ug2i(krom - 1)) / (w - 1.0)
                vg2i(krom - 1) = (w * vg2i(krom) - vg2i(krom - 1)) / (w - 1.0)
                us1i(krom - 1) = (w * us1i(krom) - us1i(krom - 1)) / (w - 1.0)
                vs1i(krom - 1) = (w * vs1i(krom) - vs1i(krom - 1)) / (w - 1.0)
                us2i(krom - 1) = (w * us2i(krom) - us2i(krom - 1)) / (w - 1.0)
                vs2i(krom - 1) = (w * vs2i(krom) - vs2i(krom - 1)) / (w - 1.0)
            enddo
            !
            if (irom>1) then
                !------- compare the best-current and best-previous integrals
                errug1 = ug1i(1) - ug1i(2)
                errvg1 = vg1i(1) - vg1i(2)
                errug2 = ug2i(1) - ug2i(2)
                errvg2 = vg2i(1) - vg2i(2)
                errus1 = us1i(1) - us1i(2)
                errvs1 = vs1i(1) - vs1i(2)
                errus2 = us2i(1) - us2i(2)
                errvs2 = vs2i(1) - vs2i(2)
                !
                err = max(abs(errug1), abs(errvg1), abs(errug2), abs(errvg2), &
                        & abs(errus1), abs(errvs1), abs(errus2), abs(errvs2))

                !         write(13,1200)
                !     &     ABS(ERRUG1),
                !     &     ABS(ERRVG1),
                !     &     ABS(ERRUG2),
                !     &     ABS(ERRVG2),
                !     &     ABS(ERRUS1),
                !     &     ABS(ERRVS1),
                !     &     ABS(ERRUS2),
                !     &     ABS(ERRVS2)
                ! 1200    format(8e16.9)

                if (err * refl<romtol) goto 101
            endif
        enddo
        write (*, *) 'LAMP: Romberg convergence failed.  Error =', err
        !

        !cc      write(*,*) IROM, ERR
        !
        !---- return best final results
        101  delsq = (x1 - x2)**2 + (r1 - r2)**2
        dels = sqrt(delsq)
        !
        ug1 = ug1i(1) * dels
        vg1 = vg1i(1) * dels
        ug2 = ug2i(1) * dels
        vg2 = vg2i(1) * dels
        us1 = us1i(1) * dels
        vs1 = vs1i(1) * dels
        us2 = us2i(1) * dels
        vs2 = vs2i(1) * dels
        !
    end subroutine lamp
    !*==LAMPC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! LAMP



    subroutine lampc(x1, r1, x2, r2, ug1, vg1, ug2, vg2, us1, vs1, us2, vs2)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: nromx = 9
        real, parameter :: pi = 3.14159265358979
        !
        ! Dummy arguments
        !
        real :: r1, r2, ug1, ug2, us1, us2, vg1, vg2, vs1, &
                & vs2, x1, x2
        !
        ! Local variables
        !
        real :: dels, delsq, dsq, dt, err, errug1, errug2, &
                & errus1, errus2, errvg1, errvg2, errvs1, errvs2, &
                & refl, rf, rt, t, tb, uga, ugai, ugt, usa, &
                & usai, ust, vga, vgai, vgt, vsa, vsai, vst, w, &
                & wm1, xf, xt
        integer :: irom, it, krom, nt
        real, save :: romtol
        real, dimension(nromx) :: ug1i, ug2i, us1i, us2i, vg1i, &
                & vg2i, vs1i, vs2i
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Same as LAMP, but the field point is assumed
        !     to be at the lampshade-panel midpoint.
        !
        !     The 1/r and log(r) singularities in the integrands
        !     are removed from the numerical integration,
        !     and are computed analytically.
        !
        !     The induced velocites returned by this routine
        !     are the average of the two values on each side
        !     of the sheet.  The sheet jumps are not included.
        !---------------------------------------------------------
        !
        !---- max number of Romberg integration stages
        !
        !
        !
        !
        !---- Romberg convergence tolerance (actual error may be much less than this)
        data romtol/1.0e-6/
        !cc      DATA ROMTOL / 1.0E-12 /
        !
        !---- lampshade meridional length**2
        delsq = (x1 - x2)**2 + (r1 - r2)**2
        !
        !---- reference length for convergence tolerance
        refl = 0.5 * (r1 + r2)
        !
        !
        !---- field point is assumed to be at midpoint
        xf = 0.5 * (x1 + x2)
        rf = 0.5 * (r1 + r2)
        !
        !---- evaluate integrals on increasingly fine grids
        !-     (start with two intervals to avoid landing right on the midpoint)
        do irom = 1, nromx
            nt = 2**irom
            !
            ug1i(irom) = 0.
            vg1i(irom) = 0.
            ug2i(irom) = 0.
            vg2i(irom) = 0.
            us1i(irom) = 0.
            vs1i(irom) = 0.
            us2i(irom) = 0.
            vs2i(irom) = 0.
            !
            !------ visit the midpoints of each of the NT intervals
            do it = 1, nt
                t = (float(it) - 0.5) / float(nt)
                tb = 1.0 - t
                !
                dt = 1.0 / float(nt)
                !
                xt = x1 * tb + x2 * t
                rt = r1 * tb + r2 * t
                !
                call ring(xt, rt, xf, rf, ugt, vgt, ust, vst)
                !
                !-------- singular parts of velocities in the limit  XT,RT -> XF,RF
                dsq = (xt - xf)**2 + (rt - rf)**2
                uga = (rt - rf) / (4.0 * pi * dsq) - 0.5 * log(dsq / (64.0 * rf**2))      &
                        & / (8.0 * pi * rf)
                vga = -(xt - xf) / (4.0 * pi * dsq)
                usa = -(xt - xf) / (4.0 * pi * dsq)
                vsa = -(rt - rf) / (4.0 * pi * dsq) - 0.5 * log(dsq / rf**2) / (8.0 * pi * rf)
                !
                !-------- accumulate integrals, with singular parts (at t=0.5) removed
                ug1i(irom) = ug1i(irom) + dt * (ugt * tb - uga)
                vg1i(irom) = vg1i(irom) + dt * (vgt * tb - vga)
                ug2i(irom) = ug2i(irom) + dt * (ugt * t - uga)
                vg2i(irom) = vg2i(irom) + dt * (vgt * t - vga)
                us1i(irom) = us1i(irom) + dt * (ust * tb - usa)
                vs1i(irom) = vs1i(irom) + dt * (vst * tb - vsa)
                us2i(irom) = us2i(irom) + dt * (ust * t - usa)
                vs2i(irom) = vs2i(irom) + dt * (vst * t - vsa)
            enddo
            !
            !------ Romberg sequence using all previous grid results
            do krom = irom, 2, -1
                !-------- weight needed to cancel lowest-order error terms in KROM level
                w = 2.0**(2 * (irom - krom + 1))
                wm1 = w - 1.0
                !
                !-------- put Richardson extrapolation for KROM level into KROM-1 level
                ug1i(krom - 1) = (w * ug1i(krom) - ug1i(krom - 1)) / wm1
                vg1i(krom - 1) = (w * vg1i(krom) - vg1i(krom - 1)) / wm1
                ug2i(krom - 1) = (w * ug2i(krom) - ug2i(krom - 1)) / wm1
                vg2i(krom - 1) = (w * vg2i(krom) - vg2i(krom - 1)) / wm1
                us1i(krom - 1) = (w * us1i(krom) - us1i(krom - 1)) / wm1
                vs1i(krom - 1) = (w * vs1i(krom) - vs1i(krom - 1)) / wm1
                us2i(krom - 1) = (w * us2i(krom) - us2i(krom - 1)) / wm1
                vs2i(krom - 1) = (w * vs2i(krom) - vs2i(krom - 1)) / wm1
            enddo
            !
            if (irom>1) then
                !------- compare the best-current and best-previous integrals
                errug1 = ug1i(1) - ug1i(2)
                errvg1 = vg1i(1) - vg1i(2)
                errug2 = ug2i(1) - ug2i(2)
                errvg2 = vg2i(1) - vg2i(2)
                errus1 = us1i(1) - us1i(2)
                errvs1 = vs1i(1) - vs1i(2)
                errus2 = us2i(1) - us2i(2)
                errvs2 = vs2i(1) - vs2i(2)
                !
                err = max(abs(errug1), abs(errvg1), abs(errug2), abs(errvg2), &
                        & abs(errus1), abs(errvs1), abs(errus2), abs(errvs2))
                !
                if (err * refl<romtol) goto 101
            endif
        enddo
        write (*, *) 'LAMPC: Romberg convergence failed.  Error =', err
        !
        !
        101  delsq = (x1 - x2)**2 + (r1 - r2)**2
        dels = sqrt(delsq)
        !
        !---- analytically-integrated singular parts which were removed
        ugai = (1.0 + log(16.0 * rf / dels)) / (4.0 * pi * rf)
        vgai = 0.
        usai = 0.
        vsai = (1.0 + log(2.0 * rf / dels)) / (4.0 * pi * rf)
        !
        !---- return final results, with removed parts added back on
        ug1 = (ug1i(1) + ugai * 0.5) * dels
        vg1 = (vg1i(1) + vgai * 0.5) * dels
        ug2 = (ug2i(1) + ugai * 0.5) * dels
        vg2 = (vg2i(1) + vgai * 0.5) * dels
        us1 = (us1i(1) + usai * 0.5) * dels
        vs1 = (vs1i(1) + vsai * 0.5) * dels
        us2 = (us2i(1) + usai * 0.5) * dels
        vs2 = (vs2i(1) + vsai * 0.5) * dels
        !
    end subroutine lampc
    !*==GLAMP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! LAMPC



    subroutine glamp(r1, r2, rf, qg1, qg1_rf, qg2, qg2_rf, qs1, qs1_rf, qs2, &
            & qs2_rf)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: nromx = 10
        !
        ! Dummy arguments
        !
        real, dimension(2) :: qg1, qg2, qs1, qs2, r1, r2, rf
        real, dimension(2, 2) :: qg1_rf, qg2_rf, qs1_rf, qs2_rf
        !
        ! Local variables
        !
        real :: dels, delsq, dt, err, errug1, errug2, errus1, &
                & errus2, errvg1, errvg2, errvs1, errvs2, refl, t, &
                & tb, w, wm1
        integer :: irom, it, j, k, krom, nt
        real, dimension(2, nromx) :: qg1i, qg2i, qs1i, qs2i
        real, dimension(2, 2, nromx) :: qg1i_rf, qg2i_rf, qs1i_rf, &
                & qs2i_rf
        real, dimension(2) :: qgt, qst, rt
        real, dimension(2, 2) :: qgt_rf, qgt_rt, qst_rf, qst_rt
        real, save :: romtol
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------------
        !     Same as LAMP, but also returns velocity gradient
        !-----------------------------------------------------------------------
        !
        !
        !
        !
        !
        !---- Romberg convergence tolerance
        data romtol/1.0e-6/
        !cc   DATA ROMTOL / 1.0E-12 /
        !
        !---- reference length for convergence tolerance
        refl = 0.5 * (r1(2) + r2(2))
        !
        !---- evaluate integrals over  0..t..1  on increasingly fine grids
        do irom = 1, nromx
            nt = 2**irom / 2
            !
            do k = 1, 2
                qg1i(k, irom) = 0.
                qg2i(k, irom) = 0.
                qs1i(k, irom) = 0.
                qs2i(k, irom) = 0.
                do j = 1, 2
                    qg1i_rf(k, j, irom) = 0.
                    qg2i_rf(k, j, irom) = 0.
                    qs1i_rf(k, j, irom) = 0.
                    qs2i_rf(k, j, irom) = 0.
                enddo
            enddo
            !
            !------ visit the midpoints of each of the NT intervals
            do it = 1, nt
                t = (float(it) - 0.5) / float(nt)
                tb = 1.0 - t
                !
                dt = 1.0 / float(nt)
                !
                rt(1) = r1(1) * tb + r2(1) * t
                rt(2) = r1(2) * tb + r2(2) * t
                !
                !-------- get induced velocities for vortex,source ring at XT,RT
                call dring(rt(1), rt(2), rf(1), rf(2), qgt(1), qgt_rt(1, 1), &
                        & qgt_rt(1, 2), qgt_rf(1, 1), qgt_rf(1, 2), qgt(2), &
                        & qgt_rt(2, 1), qgt_rt(2, 2), qgt_rf(2, 1), qgt_rf(2, 2), &
                        & qst(1), qst_rt(1, 1), qst_rt(1, 2), qst_rf(1, 1), &
                        & qst_rf(1, 2), qst(2), qst_rt(2, 1), qst_rt(2, 2), &
                        & qst_rf(2, 1), qst_rf(2, 2))
                !
                !-------- accumulate the separate unit-gamma, unit-sigma integrals
                do k = 1, 2
                    qg1i(k, irom) = qg1i(k, irom) + dt * qgt(k) * tb
                    qg2i(k, irom) = qg2i(k, irom) + dt * qgt(k) * t
                    qs1i(k, irom) = qs1i(k, irom) + dt * qst(k) * tb
                    qs2i(k, irom) = qs2i(k, irom) + dt * qst(k) * t
                    do j = 1, 2
                        qg1i_rf(k, j, irom) = qg1i_rf(k, j, irom) + dt * qgt_rf(k, j)&
                                & * tb
                        qg2i_rf(k, j, irom) = qg2i_rf(k, j, irom) + dt * qgt_rf(k, j)&
                                & * t
                        qs1i_rf(k, j, irom) = qs1i_rf(k, j, irom) + dt * qst_rf(k, j)&
                                & * tb
                        qs2i_rf(k, j, irom) = qs2i_rf(k, j, irom) + dt * qst_rf(k, j)&
                                & * t
                    enddo
                enddo
            enddo
            !
            !------ Romberg sequence using all previous grid results
            do krom = irom, 2, -1
                !-------- weight needed to cancel lowest-order error terms in KROM level
                w = 2.0**(2 * (irom - krom + 1))
                wm1 = w - 1.0
                !
                !-------- put Richardson extrapolation for KROM level into KROM-1 level
                do k = 1, 2
                    qg1i(k, krom - 1) = (w * qg1i(k, krom) - qg1i(k, krom - 1)) / wm1
                    qg2i(k, krom - 1) = (w * qg2i(k, krom) - qg2i(k, krom - 1)) / wm1
                    qs1i(k, krom - 1) = (w * qs1i(k, krom) - qs1i(k, krom - 1)) / wm1
                    qs2i(k, krom - 1) = (w * qs2i(k, krom) - qs2i(k, krom - 1)) / wm1
                    do j = 1, 2
                        qg1i_rf(k, j, krom - 1) = (w * qg1i_rf(k, j, krom) - qg1i_rf(k, j&
                                &, krom - 1)) / wm1
                        qg2i_rf(k, j, krom - 1) = (w * qg2i_rf(k, j, krom) - qg2i_rf(k, j&
                                &, krom - 1)) / wm1
                        qs1i_rf(k, j, krom - 1) = (w * qs1i_rf(k, j, krom) - qs1i_rf(k, j&
                                &, krom - 1)) / wm1
                        qs2i_rf(k, j, krom - 1) = (w * qs2i_rf(k, j, krom) - qs2i_rf(k, j&
                                &, krom - 1)) / wm1
                    enddo
                enddo
            enddo
            !
            if (irom>1) then
                !------- compare the best-current and best-previous integrals
                errug1 = qg1i(1, 1) - qg1i(1, 2)
                errvg1 = qg1i(2, 1) - qg1i(2, 2)
                errug2 = qg2i(1, 1) - qg2i(1, 2)
                errvg2 = qg2i(2, 1) - qg2i(2, 2)
                errus1 = qs1i(1, 1) - qs1i(1, 2)
                errvs1 = qs1i(2, 1) - qs1i(2, 2)
                errus2 = qs2i(1, 1) - qs2i(1, 2)
                errvs2 = qs2i(2, 1) - qs2i(2, 2)
                !
                err = max(abs(errug1), abs(errvg1), abs(errug2), abs(errvg2), &
                        & abs(errus1), abs(errvs1), abs(errus2), abs(errvs2))

                !         write(13,1200)
                !     &     ABS(ERRUG1),
                !     &     ABS(ERRVG1),
                !     &     ABS(ERRUG2),
                !     &     ABS(ERRVG2),
                !     &     ABS(ERRUS1),
                !     &     ABS(ERRVS1),
                !     &     ABS(ERRUS2),
                !     &     ABS(ERRVS2)
                ! 1200    format(8e16.9)

                if (err * refl<romtol) goto 101
            endif
        enddo
        write (*, *) 'GLAMP: Romberg convergence failed.  Error =', err
        !

        !cc      write(*,*) IROM, ERR
        !
        !---- return best final results
        101  delsq = (r1(1) - r2(1))**2 + (r1(2) - r2(2))**2
        dels = sqrt(delsq)
        !
        do k = 1, 2
            qg1(k) = qg1i(k, 1) * dels
            qg2(k) = qg2i(k, 1) * dels
            qs1(k) = qs1i(k, 1) * dels
            qs2(k) = qs2i(k, 1) * dels
            do j = 1, 2
                qg1_rf(k, j) = qg1i_rf(k, j, 1) * dels
                qg2_rf(k, j) = qg2i_rf(k, j, 1) * dels
                qs1_rf(k, j) = qs1i_rf(k, j, 1) * dels
                qs2_rf(k, j) = qs2i_rf(k, j, 1) * dels
            enddo
        enddo
        !
    end subroutine glamp
    !*==GLAMPC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GLAMP


    subroutine glampc(r1, r2, qg1, qg1_rf, qg2, qg2_rf, qs1, qs1_rf, qs2, &
            & qs2_rf)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: nromx = 10
        real, parameter :: pi = 3.14159265358979
        !
        ! Dummy arguments
        !
        real, dimension(2) :: qg1, qg2, qs1, qs2, r1, r2
        real, dimension(2, 2) :: qg1_rf, qg2_rf, qs1_rf, qs2_rf
        !
        ! Local variables
        !
        real :: dels, delsq, dr1, dr2, drsq, dsq, dt, err, &
                & errug1, errug2, errus1, errus2, errvg1, errvg2, &
                & errvs1, errvs2, pid4, pir8, refl, t, tb, w, wm1
        real, dimension(2) :: dsq_rf, qga, qgai, qgt, qsa, qsai, &
                & qst, rf, rt
        integer :: irom, it, j, k, krom, nt
        real, dimension(2, nromx) :: qg1i, qg2i, qs1i, qs2i
        real, dimension(2, 2, nromx) :: qg1i_rf, qg2i_rf, qs1i_rf, &
                & qs2i_rf
        real, dimension(2, 2) :: qgai_rf, qga_rf, qgt_rf, qgt_rt, &
                & qsai_rf, qsa_rf, qst_rf, qst_rt
        real, save :: romtol
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------------
        !     Same as LAMPC, but also returns velocity gradient
        !-----------------------------------------------------------------------
        !
        !
        !
        !
        !
        !
        !---- Romberg convergence tolerance
        data romtol/1.0e-6/
        !cc   DATA ROMTOL / 1.0E-12 /
        !
        !---- reference length for convergence tolerance
        refl = 0.5 * (r1(2) + r2(2))
        !
        rf(1) = 0.5 * (r1(1) + r2(1))
        rf(2) = 0.5 * (r1(2) + r2(2))
        !
        !---- evaluate integrals over  0..t..1  on increasingly fine grids
        do irom = 1, nromx
            nt = 2**irom
            !
            do k = 1, 2
                qg1i(k, irom) = 0.
                qg2i(k, irom) = 0.
                qs1i(k, irom) = 0.
                qs2i(k, irom) = 0.
                do j = 1, 2
                    qg1i_rf(k, j, irom) = 0.
                    qg2i_rf(k, j, irom) = 0.
                    qs1i_rf(k, j, irom) = 0.
                    qs2i_rf(k, j, irom) = 0.
                enddo
            enddo
            !
            !------ visit the midpoints of each of the NT intervals
            do it = 1, nt
                t = (float(it) - 0.5) / float(nt)
                tb = 1.0 - t
                !
                dt = 1.0 / float(nt)
                !
                rt(1) = r1(1) * tb + r2(1) * t
                rt(2) = r1(2) * tb + r2(2) * t
                !
                !-------- get induced velocities for vortex,source ring at XT,RT
                call dring(rt(1), rt(2), rf(1), rf(2), qgt(1), qgt_rt(1, 1), &
                        & qgt_rt(1, 2), qgt_rf(1, 1), qgt_rf(1, 2), qgt(2), &
                        & qgt_rt(2, 1), qgt_rt(2, 2), qgt_rf(2, 1), qgt_rf(2, 2), &
                        & qst(1), qst_rt(1, 1), qst_rt(1, 2), qst_rf(1, 1), &
                        & qst_rf(1, 2), qst(2), qst_rt(2, 1), qst_rt(2, 2), &
                        & qst_rf(2, 1), qst_rf(2, 2))
                !
                !-------- singular parts of velocities in the limit  XT,RT -> XF,RF
                dsq = (rt(1) - rf(1))**2 + (rt(2) - rf(2))**2
                dsq_rf(1) = -2.0 * (rt(1) - rf(1))
                dsq_rf(2) = -2.0 * (rt(2) - rf(2))
                !
                !
                pid4 = 4.0 * pi * dsq
                pir8 = 8.0 * pi * rf(2)
                drsq = dsq / rf(2)**2
                !
                dr1 = (rt(1) - rf(1)) / pid4
                dr2 = (rt(2) - rf(2)) / pid4
                !
                qga(1) = dr2 - 0.5 * log(drsq / 64.0) / pir8
                qga(2) = -dr1
                qsa(1) = -dr1
                qsa(2) = -dr2 - 0.5 * log(drsq) / pir8
                !
                qga_rf(1, 1) = (-dr2 - 0.5 / pir8) / dsq * dsq_rf(1)
                qga_rf(1, 2) = (-dr2 - 0.5 / pir8) / dsq * dsq_rf(2) - 1.0 / pid4 + &
                        & (1.0 + 0.5 * log(drsq / 64.0)) / (rf(2) * pir8)
                !
                qga_rf(2, 1) = dr1 / dsq * dsq_rf(1) + 1.0 / pid4
                qga_rf(2, 2) = dr1 / dsq * dsq_rf(2)
                !
                qsa_rf(1, 1) = dr1 / dsq * dsq_rf(1) + 1.0 / pid4
                qsa_rf(1, 2) = dr1 / dsq * dsq_rf(2)
                !
                qsa_rf(2, 1) = (dr2 - 0.5 / pir8) / dsq * dsq_rf(1)
                qsa_rf(2, 2) = (dr2 - 0.5 / pir8) / dsq * dsq_rf(2) + 1.0 / pid4 + &
                        & (1.0 + 0.5 * log(drsq)) / (rf(2) * pir8)
                !
                !-------- accumulate integrals, with singular parts (at t=0.5) removed
                do k = 1, 2
                    qg1i(k, irom) = qg1i(k, irom) + dt * (qgt(k) * tb - qga(k))
                    qg2i(k, irom) = qg2i(k, irom) + dt * (qgt(k) * t - qga(k))
                    qs1i(k, irom) = qs1i(k, irom) + dt * (qst(k) * tb - qsa(k))
                    qs2i(k, irom) = qs2i(k, irom) + dt * (qst(k) * t - qsa(k))
                    do j = 1, 2
                        qg1i_rf(k, j, irom) = qg1i_rf(k, j, irom)                 &
                                & + dt * (qgt_rf(k, j) * tb - qga_rf(k, j))
                        qg2i_rf(k, j, irom) = qg2i_rf(k, j, irom)                 &
                                & + dt * (qgt_rf(k, j) * t - qga_rf(k, j))
                        qs1i_rf(k, j, irom) = qs1i_rf(k, j, irom)                 &
                                & + dt * (qst_rf(k, j) * tb - qsa_rf(k, j))
                        qs2i_rf(k, j, irom) = qs2i_rf(k, j, irom)                 &
                                & + dt * (qst_rf(k, j) * t - qsa_rf(k, j))
                    enddo
                enddo
            enddo
            !
            !------ Romberg sequence using all previous grid results
            do krom = irom, 2, -1
                !-------- weight needed to cancel lowest-order error terms in KROM level
                w = 2.0**(2 * (irom - krom + 1))
                wm1 = w - 1.0
                !
                !-------- put Richardson extrapolation for KROM level into KROM-1 level
                do k = 1, 2
                    qg1i(k, krom - 1) = (w * qg1i(k, krom) - qg1i(k, krom - 1)) / wm1
                    qg2i(k, krom - 1) = (w * qg2i(k, krom) - qg2i(k, krom - 1)) / wm1
                    qs1i(k, krom - 1) = (w * qs1i(k, krom) - qs1i(k, krom - 1)) / wm1
                    qs2i(k, krom - 1) = (w * qs2i(k, krom) - qs2i(k, krom - 1)) / wm1
                    do j = 1, 2
                        qg1i_rf(k, j, krom - 1) = (w * qg1i_rf(k, j, krom) - qg1i_rf(k, j&
                                &, krom - 1)) / wm1
                        qg2i_rf(k, j, krom - 1) = (w * qg2i_rf(k, j, krom) - qg2i_rf(k, j&
                                &, krom - 1)) / wm1
                        qs1i_rf(k, j, krom - 1) = (w * qs1i_rf(k, j, krom) - qs1i_rf(k, j&
                                &, krom - 1)) / wm1
                        qs2i_rf(k, j, krom - 1) = (w * qs2i_rf(k, j, krom) - qs2i_rf(k, j&
                                &, krom - 1)) / wm1
                    enddo
                enddo
            enddo
            !
            if (irom>1) then
                !------- compare the best-current and best-previous integrals
                errug1 = qg1i(1, 1) - qg1i(1, 2)
                errvg1 = qg1i(2, 1) - qg1i(2, 2)
                errug2 = qg2i(1, 1) - qg2i(1, 2)
                errvg2 = qg2i(2, 1) - qg2i(2, 2)
                errus1 = qs1i(1, 1) - qs1i(1, 2)
                errvs1 = qs1i(2, 1) - qs1i(2, 2)
                errus2 = qs2i(1, 1) - qs2i(1, 2)
                errvs2 = qs2i(2, 1) - qs2i(2, 2)
                !
                err = max(abs(errug1), abs(errvg1), abs(errug2), abs(errvg2), &
                        & abs(errus1), abs(errvs1), abs(errus2), abs(errvs2))

                !         write(13,1200)
                !     &     ABS(ERRUG1),
                !     &     ABS(ERRVG1),
                !     &     ABS(ERRUG2),
                !     &     ABS(ERRVG2),
                !     &     ABS(ERRUS1),
                !     &     ABS(ERRVS1),
                !     &     ABS(ERRUS2),
                !     &     ABS(ERRVS2)
                ! 1200    format(8e16.9)

                if (err * refl<romtol) goto 101
            endif
        enddo
        write (*, *) 'GLAMPC: Romberg convergence failed.  Error =', err
        !

        !cc      write(*,*) IROM, ERR
        !
        !---- return best final results
        101  delsq = (r1(1) - r2(1))**2 + (r1(2) - r2(2))**2
        dels = sqrt(delsq)
        !
        !---- analytically-integrated singular parts which were removed
        qgai(1) = (1.0 + log(16.0 * rf(2) / dels)) / (8.0 * pi * rf(2))
        qgai(2) = 0.
        qsai(1) = 0.
        qsai(2) = (1.0 + log(2.0 * rf(2) / dels)) / (8.0 * pi * rf(2))
        !
        qgai_rf(1, 1) = 0.
        qgai_rf(2, 1) = 0.
        qsai_rf(1, 1) = 0.
        qsai_rf(2, 1) = 0.
        !
        qgai_rf(1, 2) = (1.0 / rf(2)) / (8.0 * pi * rf(2)) - qgai(1) / rf(2)
        qgai_rf(2, 2) = 0.
        qsai_rf(1, 2) = 0.
        qsai_rf(2, 2) = (1.0 / rf(2)) / (8.0 * pi * rf(2)) - qsai(2) / rf(2)
        !
        !
        do k = 1, 2
            qg1(k) = (qg1i(k, 1) + qgai(k)) * dels
            qg2(k) = (qg2i(k, 1) + qgai(k)) * dels
            qs1(k) = (qs1i(k, 1) + qsai(k)) * dels
            qs2(k) = (qs2i(k, 1) + qsai(k)) * dels
            do j = 1, 2
                qg1_rf(k, j) = (qg1i_rf(k, j, 1) + qgai_rf(k, j)) * dels
                qg2_rf(k, j) = (qg2i_rf(k, j, 1) + qgai_rf(k, j)) * dels
                qs1_rf(k, j) = (qs1i_rf(k, j, 1) + qsai_rf(k, j)) * dels
                qs2_rf(k, j) = (qs2i_rf(k, j, 1) + qsai_rf(k, j)) * dels
            enddo
        enddo
        !
    end subroutine glampc
    !*==DLAMP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GLAMPC



    subroutine dlamp(r1, r2, rf, qg1, qg1_r1, qg1_r2, qg1_rf, qg2, qg2_r1, &
            & qg2_r2, qg2_rf, qs1, qs1_r1, qs1_r2, qs1_rf, qs2, &
            & qs2_r1, qs2_r2, qs2_rf)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: nromx = 10
        !
        ! Dummy arguments
        !
        real, dimension(2) :: qg1, qg2, qs1, qs2, r1, r2, rf
        real, dimension(2, 2) :: qg1_r1, qg1_r2, qg1_rf, qg2_r1, &
                & qg2_r2, qg2_rf, qs1_r1, qs1_r2, &
                & qs1_rf, qs2_r1, qs2_r2, qs2_rf
        !
        ! Local variables
        !
        real :: dels, delsq, dt, err, errug1, errug2, errus1, &
                & errus2, errvg1, errvg2, errvs1, errvs2, refl, t, &
                & tb, w, wm1
        real, dimension(2) :: dels_r1, dels_r2, qgt, qst, rt
        integer :: irom, it, j, k, krom, nt
        real, dimension(2, nromx) :: qg1i, qg2i, qs1i, qs2i
        real, dimension(2, 2, nromx) :: qg1i_r1, qg1i_r2, qg1i_rf, &
                & qg2i_r1, qg2i_r2, qg2i_rf, &
                & qs1i_r1, qs1i_r2, qs1i_rf, &
                & qs2i_r1, qs2i_r2, qs2i_rf
        real, dimension(2, 2) :: qgt_rf, qgt_rt, qst_rf, qst_rt
        real, save :: romtol
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------------
        !     Same as LAMP, but also returns derivatives.
        !-----------------------------------------------------------------------
        !
        !
        !
        !
        !
        !---- Romberg convergence tolerance
        data romtol/1.0e-6/
        !cc   DATA ROMTOL / 1.0E-12 /
        !
        !---- reference length for convergence tolerance
        refl = 0.5 * (r1(2) + r2(2))
        !
        !---- evaluate integrals over  0..t..1  on increasingly fine grids
        do irom = 1, nromx
            nt = 2**irom / 2
            !
            do k = 1, 2
                qg1i(k, irom) = 0.
                qg2i(k, irom) = 0.
                qs1i(k, irom) = 0.
                qs2i(k, irom) = 0.
                do j = 1, 2
                    qg1i_r1(k, j, irom) = 0.
                    qg2i_r1(k, j, irom) = 0.
                    qs1i_r1(k, j, irom) = 0.
                    qs2i_r1(k, j, irom) = 0.
                    qg1i_r2(k, j, irom) = 0.
                    qg2i_r2(k, j, irom) = 0.
                    qs1i_r2(k, j, irom) = 0.
                    qs2i_r2(k, j, irom) = 0.
                    qg1i_rf(k, j, irom) = 0.
                    qg2i_rf(k, j, irom) = 0.
                    qs1i_rf(k, j, irom) = 0.
                    qs2i_rf(k, j, irom) = 0.
                enddo
            enddo
            !
            !------ visit the midpoints of each of the NT intervals
            do it = 1, nt
                t = (float(it) - 0.5) / float(nt)
                tb = 1.0 - t
                !
                dt = 1.0 / float(nt)
                !
                rt(1) = r1(1) * tb + r2(1) * t
                rt(2) = r1(2) * tb + r2(2) * t
                !
                !-------- get induced velocities for vortex,source ring at XT,RT
                call dring(rt(1), rt(2), rf(1), rf(2), qgt(1), qgt_rt(1, 1), &
                        & qgt_rt(1, 2), qgt_rf(1, 1), qgt_rf(1, 2), qgt(2), &
                        & qgt_rt(2, 1), qgt_rt(2, 2), qgt_rf(2, 1), qgt_rf(2, 2), &
                        & qst(1), qst_rt(1, 1), qst_rt(1, 2), qst_rf(1, 1), &
                        & qst_rf(1, 2), qst(2), qst_rt(2, 1), qst_rt(2, 2), &
                        & qst_rf(2, 1), qst_rf(2, 2))
                !
                !-------- accumulate the separate unit-gamma, unit-sigma integrals
                do k = 1, 2
                    qg1i(k, irom) = qg1i(k, irom) + dt * qgt(k) * tb
                    qg2i(k, irom) = qg2i(k, irom) + dt * qgt(k) * t
                    qs1i(k, irom) = qs1i(k, irom) + dt * qst(k) * tb
                    qs2i(k, irom) = qs2i(k, irom) + dt * qst(k) * t
                    do j = 1, 2
                        qg1i_r1(k, j, irom) = qg1i_r1(k, j, irom) + dt * qgt_rt(k, j)&
                                & * tb * tb
                        qg2i_r1(k, j, irom) = qg2i_r1(k, j, irom) + dt * qgt_rt(k, j)&
                                & * t * tb
                        qs1i_r1(k, j, irom) = qs1i_r1(k, j, irom) + dt * qst_rt(k, j)&
                                & * tb * tb
                        qs2i_r1(k, j, irom) = qs2i_r1(k, j, irom) + dt * qst_rt(k, j)&
                                & * t * tb
                        !
                        qg1i_r2(k, j, irom) = qg1i_r2(k, j, irom) + dt * qgt_rt(k, j)&
                                & * tb * t
                        qg2i_r2(k, j, irom) = qg2i_r2(k, j, irom) + dt * qgt_rt(k, j)&
                                & * t * t
                        qs1i_r2(k, j, irom) = qs1i_r2(k, j, irom) + dt * qst_rt(k, j)&
                                & * tb * t
                        qs2i_r2(k, j, irom) = qs2i_r2(k, j, irom) + dt * qst_rt(k, j)&
                                & * t * t
                        !
                        qg1i_rf(k, j, irom) = qg1i_rf(k, j, irom) + dt * qgt_rf(k, j)&
                                & * tb
                        qg2i_rf(k, j, irom) = qg2i_rf(k, j, irom) + dt * qgt_rf(k, j)&
                                & * t
                        qs1i_rf(k, j, irom) = qs1i_rf(k, j, irom) + dt * qst_rf(k, j)&
                                & * tb
                        qs2i_rf(k, j, irom) = qs2i_rf(k, j, irom) + dt * qst_rf(k, j)&
                                & * t
                    enddo
                enddo
            enddo
            !
            !------ Romberg sequence using all previous grid results
            do krom = irom, 2, -1
                !-------- weight needed to cancel lowest-order error terms in KROM level
                w = 2.0**(2 * (irom - krom + 1))
                wm1 = w - 1.0
                !
                !-------- put Richardson extrapolation for KROM level into KROM-1 level
                do k = 1, 2
                    qg1i(k, krom - 1) = (w * qg1i(k, krom) - qg1i(k, krom - 1)) / wm1
                    qg2i(k, krom - 1) = (w * qg2i(k, krom) - qg2i(k, krom - 1)) / wm1
                    qs1i(k, krom - 1) = (w * qs1i(k, krom) - qs1i(k, krom - 1)) / wm1
                    qs2i(k, krom - 1) = (w * qs2i(k, krom) - qs2i(k, krom - 1)) / wm1
                    do j = 1, 2
                        qg1i_r1(k, j, krom - 1) = (w * qg1i_r1(k, j, krom) - qg1i_r1(k, j&
                                &, krom - 1)) / wm1
                        qg2i_r1(k, j, krom - 1) = (w * qg2i_r1(k, j, krom) - qg2i_r1(k, j&
                                &, krom - 1)) / wm1
                        qs1i_r1(k, j, krom - 1) = (w * qs1i_r1(k, j, krom) - qs1i_r1(k, j&
                                &, krom - 1)) / wm1
                        qs2i_r1(k, j, krom - 1) = (w * qs2i_r1(k, j, krom) - qs2i_r1(k, j&
                                &, krom - 1)) / wm1
                        !
                        qg1i_r2(k, j, krom - 1) = (w * qg1i_r2(k, j, krom) - qg1i_r2(k, j&
                                &, krom - 1)) / wm1
                        qg2i_r2(k, j, krom - 1) = (w * qg2i_r2(k, j, krom) - qg2i_r2(k, j&
                                &, krom - 1)) / wm1
                        qs1i_r2(k, j, krom - 1) = (w * qs1i_r2(k, j, krom) - qs1i_r2(k, j&
                                &, krom - 1)) / wm1
                        qs2i_r2(k, j, krom - 1) = (w * qs2i_r2(k, j, krom) - qs2i_r2(k, j&
                                &, krom - 1)) / wm1
                        !
                        qg1i_rf(k, j, krom - 1) = (w * qg1i_rf(k, j, krom) - qg1i_rf(k, j&
                                &, krom - 1)) / wm1
                        qg2i_rf(k, j, krom - 1) = (w * qg2i_rf(k, j, krom) - qg2i_rf(k, j&
                                &, krom - 1)) / wm1
                        qs1i_rf(k, j, krom - 1) = (w * qs1i_rf(k, j, krom) - qs1i_rf(k, j&
                                &, krom - 1)) / wm1
                        qs2i_rf(k, j, krom - 1) = (w * qs2i_rf(k, j, krom) - qs2i_rf(k, j&
                                &, krom - 1)) / wm1
                    enddo
                enddo
            enddo
            !
            if (irom>1) then
                !------- compare the best-current and best-previous integrals
                errug1 = qg1i(1, 1) - qg1i(1, 2)
                errvg1 = qg1i(2, 1) - qg1i(2, 2)
                errug2 = qg2i(1, 1) - qg2i(1, 2)
                errvg2 = qg2i(2, 1) - qg2i(2, 2)
                errus1 = qs1i(1, 1) - qs1i(1, 2)
                errvs1 = qs1i(2, 1) - qs1i(2, 2)
                errus2 = qs2i(1, 1) - qs2i(1, 2)
                errvs2 = qs2i(2, 1) - qs2i(2, 2)
                !
                err = max(abs(errug1), abs(errvg1), abs(errug2), abs(errvg2), &
                        & abs(errus1), abs(errvs1), abs(errus2), abs(errvs2))

                !         write(13,1200)
                !     &     ABS(ERRUG1),
                !     &     ABS(ERRVG1),
                !     &     ABS(ERRUG2),
                !     &     ABS(ERRVG2),
                !     &     ABS(ERRUS1),
                !     &     ABS(ERRVS1),
                !     &     ABS(ERRUS2),
                !     &     ABS(ERRVS2)
                ! 1200    format(8e16.9)

                if (err * refl<romtol) goto 101
            endif
        enddo
        write (*, *) 'DLAMP: Romberg convergence failed.  Error =', err
        !

        !cc      write(*,*) IROM, ERR
        !
        !---- return best final results
        101  delsq = (r1(1) - r2(1))**2 + (r1(2) - r2(2))**2
        dels = sqrt(delsq)
        dels_r1(1) = (r1(1) - r2(1)) / dels
        dels_r2(1) = -(r1(1) - r2(1)) / dels
        dels_r1(2) = (r1(2) - r2(2)) / dels
        dels_r2(2) = -(r1(2) - r2(2)) / dels
        !
        do k = 1, 2
            qg1(k) = qg1i(k, 1) * dels
            qg2(k) = qg2i(k, 1) * dels
            qs1(k) = qs1i(k, 1) * dels
            qs2(k) = qs2i(k, 1) * dels
            do j = 1, 2
                qg1_r1(k, j) = qg1i_r1(k, j, 1) * dels + qg1i(k, 1) * dels_r1(j)
                qg2_r1(k, j) = qg2i_r1(k, j, 1) * dels + qg2i(k, 1) * dels_r1(j)
                qs1_r1(k, j) = qs1i_r1(k, j, 1) * dels + qs1i(k, 1) * dels_r1(j)
                qs2_r1(k, j) = qs2i_r1(k, j, 1) * dels + qs2i(k, 1) * dels_r1(j)
                !
                qg1_r2(k, j) = qg1i_r2(k, j, 1) * dels + qg1i(k, 1) * dels_r2(j)
                qg2_r2(k, j) = qg2i_r2(k, j, 1) * dels + qg2i(k, 1) * dels_r2(j)
                qs1_r2(k, j) = qs1i_r2(k, j, 1) * dels + qs1i(k, 1) * dels_r2(j)
                qs2_r2(k, j) = qs2i_r2(k, j, 1) * dels + qs2i(k, 1) * dels_r2(j)
                !
                qg1_rf(k, j) = qg1i_rf(k, j, 1) * dels
                qg2_rf(k, j) = qg2i_rf(k, j, 1) * dels
                qs1_rf(k, j) = qs1i_rf(k, j, 1) * dels
                qs2_rf(k, j) = qs2i_rf(k, j, 1) * dels
            enddo
        enddo
        !
    end subroutine dlamp
    !*==DLAMPC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! DLAMP



    subroutine dlampc(r1, r2, qg1, qg1_r1, qg1_r2, qg2, qg2_r1, qg2_r2, qs1, &
            & qs1_r1, qs1_r2, qs2, qs2_r1, qs2_r2)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: nromx = 10
        real, parameter :: pi = 3.14159265358979
        !
        ! Dummy arguments
        !
        real, dimension(2) :: qg1, qg2, qs1, qs2, r1, r2
        real, dimension(2, 2) :: qg1_r1, qg1_r2, qg2_r1, qg2_r2, &
                & qs1_r1, qs1_r2, qs2_r1, qs2_r2
        !
        ! Local variables
        !
        real :: dels, delsq, dr1, dr2, drsq, dsq, dt, err, &
                & errug1, errug2, errus1, errus2, errvg1, errvg2, &
                & errvs1, errvs2, pid4, pir8, refl, t, tb, w, wm1
        real, dimension(2) :: dels_r1, dels_r2, dsq_rf, dsq_rt, &
                & qga, qgai, qgt, qsa, qsai, qst, &
                & rf, rt
        integer :: irom, it, j, k, krom, nt
        real, dimension(2, nromx) :: qg1i, qg2i, qs1i, qs2i
        real, dimension(2, 2, nromx) :: qg1i_r1, qg1i_r2, qg1i_rf, &
                & qg2i_r1, qg2i_r2, qg2i_rf, &
                & qs1i_r1, qs1i_r2, qs1i_rf, &
                & qs2i_r1, qs2i_r2, qs2i_rf
        real, dimension(2, 2) :: qgai_r1, qgai_r2, qgai_rf, qga_rf, &
                & qga_rt, qgt_rf, qgt_rt, qsai_r1, &
                & qsai_r2, qsai_rf, qsa_rf, qsa_rt, &
                & qst_rf, qst_rt
        real, save :: romtol
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------------
        !     Same as LAMPC, but also returns derivatives.
        !-----------------------------------------------------------------------
        !
        !
        !
        !
        !
        !
        !---- Romberg convergence tolerance
        data romtol/1.0e-6/
        !cc   DATA ROMTOL / 1.0E-12 /
        !
        !---- reference length for convergence tolerance
        refl = 0.5 * (r1(2) + r2(2))
        !
        rf(1) = 0.5 * (r1(1) + r2(1))
        rf(2) = 0.5 * (r1(2) + r2(2))
        !
        !---- evaluate integrals over  0..t..1  on increasingly fine grids
        do irom = 1, nromx
            nt = 2**irom
            !
            do k = 1, 2
                qg1i(k, irom) = 0.
                qg2i(k, irom) = 0.
                qs1i(k, irom) = 0.
                qs2i(k, irom) = 0.
                do j = 1, 2
                    qg1i_r1(k, j, irom) = 0.
                    qg2i_r1(k, j, irom) = 0.
                    qs1i_r1(k, j, irom) = 0.
                    qs2i_r1(k, j, irom) = 0.
                    qg1i_r2(k, j, irom) = 0.
                    qg2i_r2(k, j, irom) = 0.
                    qs1i_r2(k, j, irom) = 0.
                    qs2i_r2(k, j, irom) = 0.
                    qg1i_rf(k, j, irom) = 0.
                    qg2i_rf(k, j, irom) = 0.
                    qs1i_rf(k, j, irom) = 0.
                    qs2i_rf(k, j, irom) = 0.
                enddo
            enddo
            !
            !------ visit the midpoints of each of the NT intervals
            do it = 1, nt
                t = (float(it) - 0.5) / float(nt)
                tb = 1.0 - t
                !
                dt = 1.0 / float(nt)
                !
                rt(1) = r1(1) * tb + r2(1) * t
                rt(2) = r1(2) * tb + r2(2) * t
                !
                !-------- get induced velocities for vortex,source ring at XT,RT
                call dring(rt(1), rt(2), rf(1), rf(2), qgt(1), qgt_rt(1, 1), &
                        & qgt_rt(1, 2), qgt_rf(1, 1), qgt_rf(1, 2), qgt(2), &
                        & qgt_rt(2, 1), qgt_rt(2, 2), qgt_rf(2, 1), qgt_rf(2, 2), &
                        & qst(1), qst_rt(1, 1), qst_rt(1, 2), qst_rf(1, 1), &
                        & qst_rf(1, 2), qst(2), qst_rt(2, 1), qst_rt(2, 2), &
                        & qst_rf(2, 1), qst_rf(2, 2))
                !
                !-------- singular parts of velocities in the limit  XT,RT -> XF,RF
                dsq = (rt(1) - rf(1))**2 + (rt(2) - rf(2))**2
                dsq_rt(1) = 2.0 * (rt(1) - rf(1))
                dsq_rt(2) = 2.0 * (rt(2) - rf(2))
                dsq_rf(1) = -2.0 * (rt(1) - rf(1))
                dsq_rf(2) = -2.0 * (rt(2) - rf(2))
                !
                !
                pid4 = 4.0 * pi * dsq
                pir8 = 8.0 * pi * rf(2)
                drsq = dsq / rf(2)**2
                !
                dr1 = (rt(1) - rf(1)) / pid4
                dr2 = (rt(2) - rf(2)) / pid4
                !
                qga(1) = dr2 - 0.5 * log(drsq / 64.0) / pir8
                qga(2) = -dr1
                qsa(1) = -dr1
                qsa(2) = -dr2 - 0.5 * log(drsq) / pir8
                !
                qga_rt(1, 1) = (-dr2 - 0.5 / pir8) / dsq * dsq_rt(1)
                qga_rt(1, 2) = (-dr2 - 0.5 / pir8) / dsq * dsq_rt(2) + 1.0 / pid4
                qga_rf(1, 1) = (-dr2 - 0.5 / pir8) / dsq * dsq_rf(1)
                qga_rf(1, 2) = (-dr2 - 0.5 / pir8) / dsq * dsq_rf(2) - 1.0 / pid4 + &
                        & (1.0 + 0.5 * log(drsq / 64.0)) / (rf(2) * pir8)
                !
                qga_rt(2, 1) = dr1 / dsq * dsq_rt(1) - 1.0 / pid4
                qga_rt(2, 2) = dr1 / dsq * dsq_rt(2)
                qga_rf(2, 1) = dr1 / dsq * dsq_rf(1) + 1.0 / pid4
                qga_rf(2, 2) = dr1 / dsq * dsq_rf(2)
                !
                qsa_rt(1, 1) = dr1 / dsq * dsq_rt(1) - 1.0 / pid4
                qsa_rt(1, 2) = dr1 / dsq * dsq_rt(2)
                qsa_rf(1, 1) = dr1 / dsq * dsq_rf(1) + 1.0 / pid4
                qsa_rf(1, 2) = dr1 / dsq * dsq_rf(2)
                !
                qsa_rt(2, 1) = (dr2 - 0.5 / pir8) / dsq * dsq_rt(1)
                qsa_rt(2, 2) = (dr2 - 0.5 / pir8) / dsq * dsq_rt(2) - 1.0 / pid4
                qsa_rf(2, 1) = (dr2 - 0.5 / pir8) / dsq * dsq_rf(1)
                qsa_rf(2, 2) = (dr2 - 0.5 / pir8) / dsq * dsq_rf(2) + 1.0 / pid4 + &
                        & (1.0 + 0.5 * log(drsq)) / (rf(2) * pir8)
                !
                !-------- accumulate integrals, with singular parts (at t=0.5) removed
                do k = 1, 2
                    qg1i(k, irom) = qg1i(k, irom) + dt * (qgt(k) * tb - qga(k))
                    qg2i(k, irom) = qg2i(k, irom) + dt * (qgt(k) * t - qga(k))
                    qs1i(k, irom) = qs1i(k, irom) + dt * (qst(k) * tb - qsa(k))
                    qs2i(k, irom) = qs2i(k, irom) + dt * (qst(k) * t - qsa(k))
                    do j = 1, 2
                        qg1i_r1(k, j, irom) = qg1i_r1(k, j, irom)                 &
                                & + dt * (qgt_rt(k, j) * tb - qga_rt(k, j)) &
                                        & * tb
                        qg2i_r1(k, j, irom) = qg2i_r1(k, j, irom)                 &
                                & + dt * (qgt_rt(k, j) * t - qga_rt(k, j))  &
                                        & * tb
                        qs1i_r1(k, j, irom) = qs1i_r1(k, j, irom)                 &
                                & + dt * (qst_rt(k, j) * tb - qsa_rt(k, j)) &
                                        & * tb
                        qs2i_r1(k, j, irom) = qs2i_r1(k, j, irom)                 &
                                & + dt * (qst_rt(k, j) * t - qsa_rt(k, j))  &
                                        & * tb
                        !
                        qg1i_r2(k, j, irom) = qg1i_r2(k, j, irom)                 &
                                & + dt * (qgt_rt(k, j) * tb - qga_rt(k, j)) &
                                        & * t
                        qg2i_r2(k, j, irom) = qg2i_r2(k, j, irom)                 &
                                & + dt * (qgt_rt(k, j) * t - qga_rt(k, j)) * t
                        qs1i_r2(k, j, irom) = qs1i_r2(k, j, irom)                 &
                                & + dt * (qst_rt(k, j) * tb - qsa_rt(k, j)) &
                                        & * t
                        qs2i_r2(k, j, irom) = qs2i_r2(k, j, irom)                 &
                                & + dt * (qst_rt(k, j) * t - qsa_rt(k, j)) * t
                        !
                        qg1i_rf(k, j, irom) = qg1i_rf(k, j, irom)                 &
                                & + dt * (qgt_rf(k, j) * tb - qga_rf(k, j))
                        qg2i_rf(k, j, irom) = qg2i_rf(k, j, irom)                 &
                                & + dt * (qgt_rf(k, j) * t - qga_rf(k, j))
                        qs1i_rf(k, j, irom) = qs1i_rf(k, j, irom)                 &
                                & + dt * (qst_rf(k, j) * tb - qsa_rf(k, j))
                        qs2i_rf(k, j, irom) = qs2i_rf(k, j, irom)                 &
                                & + dt * (qst_rf(k, j) * t - qsa_rf(k, j))
                    enddo
                enddo
            enddo
            !
            !------ Romberg sequence using all previous grid results
            do krom = irom, 2, -1
                !-------- weight needed to cancel lowest-order error terms in KROM level
                w = 2.0**(2 * (irom - krom + 1))
                wm1 = w - 1.0
                !
                !-------- put Richardson extrapolation for KROM level into KROM-1 level
                do k = 1, 2
                    qg1i(k, krom - 1) = (w * qg1i(k, krom) - qg1i(k, krom - 1)) / wm1
                    qg2i(k, krom - 1) = (w * qg2i(k, krom) - qg2i(k, krom - 1)) / wm1
                    qs1i(k, krom - 1) = (w * qs1i(k, krom) - qs1i(k, krom - 1)) / wm1
                    qs2i(k, krom - 1) = (w * qs2i(k, krom) - qs2i(k, krom - 1)) / wm1
                    do j = 1, 2
                        qg1i_r1(k, j, krom - 1) = (w * qg1i_r1(k, j, krom) - qg1i_r1(k, j&
                                &, krom - 1)) / wm1
                        qg2i_r1(k, j, krom - 1) = (w * qg2i_r1(k, j, krom) - qg2i_r1(k, j&
                                &, krom - 1)) / wm1
                        qs1i_r1(k, j, krom - 1) = (w * qs1i_r1(k, j, krom) - qs1i_r1(k, j&
                                &, krom - 1)) / wm1
                        qs2i_r1(k, j, krom - 1) = (w * qs2i_r1(k, j, krom) - qs2i_r1(k, j&
                                &, krom - 1)) / wm1
                        !
                        qg1i_r2(k, j, krom - 1) = (w * qg1i_r2(k, j, krom) - qg1i_r2(k, j&
                                &, krom - 1)) / wm1
                        qg2i_r2(k, j, krom - 1) = (w * qg2i_r2(k, j, krom) - qg2i_r2(k, j&
                                &, krom - 1)) / wm1
                        qs1i_r2(k, j, krom - 1) = (w * qs1i_r2(k, j, krom) - qs1i_r2(k, j&
                                &, krom - 1)) / wm1
                        qs2i_r2(k, j, krom - 1) = (w * qs2i_r2(k, j, krom) - qs2i_r2(k, j&
                                &, krom - 1)) / wm1
                        !
                        qg1i_rf(k, j, krom - 1) = (w * qg1i_rf(k, j, krom) - qg1i_rf(k, j&
                                &, krom - 1)) / wm1
                        qg2i_rf(k, j, krom - 1) = (w * qg2i_rf(k, j, krom) - qg2i_rf(k, j&
                                &, krom - 1)) / wm1
                        qs1i_rf(k, j, krom - 1) = (w * qs1i_rf(k, j, krom) - qs1i_rf(k, j&
                                &, krom - 1)) / wm1
                        qs2i_rf(k, j, krom - 1) = (w * qs2i_rf(k, j, krom) - qs2i_rf(k, j&
                                &, krom - 1)) / wm1
                    enddo
                enddo
            enddo
            !
            if (irom>1) then
                !------- compare the best-current and best-previous integrals
                errug1 = qg1i(1, 1) - qg1i(1, 2)
                errvg1 = qg1i(2, 1) - qg1i(2, 2)
                errug2 = qg2i(1, 1) - qg2i(1, 2)
                errvg2 = qg2i(2, 1) - qg2i(2, 2)
                errus1 = qs1i(1, 1) - qs1i(1, 2)
                errvs1 = qs1i(2, 1) - qs1i(2, 2)
                errus2 = qs2i(1, 1) - qs2i(1, 2)
                errvs2 = qs2i(2, 1) - qs2i(2, 2)
                !
                err = max(abs(errug1), abs(errvg1), abs(errug2), abs(errvg2), &
                        & abs(errus1), abs(errvs1), abs(errus2), abs(errvs2))

                !         write(13,1200)
                !     &     ABS(ERRUG1),
                !     &     ABS(ERRVG1),
                !     &     ABS(ERRUG2),
                !     &     ABS(ERRVG2),
                !     &     ABS(ERRUS1),
                !     &     ABS(ERRVS1),
                !     &     ABS(ERRUS2),
                !     &     ABS(ERRVS2)
                ! 1200    format(8e16.9)

                if (err * refl<romtol) goto 101
            endif
        enddo
        write (*, *) 'DLAMPC: Romberg convergence failed.  Error =', err
        !

        !cc      write(*,*) IROM, ERR
        !
        !---- return best final results
        101  delsq = (r1(1) - r2(1))**2 + (r1(2) - r2(2))**2
        dels = sqrt(delsq)
        dels_r1(1) = (r1(1) - r2(1)) / dels
        dels_r2(1) = -(r1(1) - r2(1)) / dels
        dels_r1(2) = (r1(2) - r2(2)) / dels
        dels_r2(2) = -(r1(2) - r2(2)) / dels
        !
        !---- analytically-integrated singular parts which were removed
        qgai(1) = (1.0 + log(16.0 * rf(2) / dels)) / (8.0 * pi * rf(2))
        qgai(2) = 0.
        qsai(1) = 0.
        qsai(2) = (1.0 + log(2.0 * rf(2) / dels)) / (8.0 * pi * rf(2))
        do j = 1, 2
            qgai_r1(1, j) = (-dels_r1(j) / dels) / (8.0 * pi * rf(2))
            qgai_r1(2, j) = 0.
            qsai_r1(1, j) = 0.
            qsai_r1(2, j) = (-dels_r1(j) / dels) / (8.0 * pi * rf(2))
            !
            qgai_r2(1, j) = (-dels_r2(j) / dels) / (8.0 * pi * rf(2))
            qgai_r2(2, j) = 0.
            qsai_r2(1, j) = 0.
            qsai_r2(2, j) = (-dels_r2(j) / dels) / (8.0 * pi * rf(2))
        enddo
        !
        qgai_rf(1, 1) = 0.
        qgai_rf(2, 1) = 0.
        qsai_rf(1, 1) = 0.
        qsai_rf(2, 1) = 0.
        !
        qgai_rf(1, 2) = (1.0 / rf(2)) / (8.0 * pi * rf(2)) - qgai(1) / rf(2)
        qgai_rf(2, 2) = 0.
        qsai_rf(1, 2) = 0.
        qsai_rf(2, 2) = (1.0 / rf(2)) / (8.0 * pi * rf(2)) - qsai(2) / rf(2)
        !
        !
        do k = 1, 2
            qg1(k) = (qg1i(k, 1) + qgai(k)) * dels
            qg2(k) = (qg2i(k, 1) + qgai(k)) * dels
            qs1(k) = (qs1i(k, 1) + qsai(k)) * dels
            qs2(k) = (qs2i(k, 1) + qsai(k)) * dels
            do j = 1, 2
                qg1_r1(k, j) = (qg1i_r1(k, j, 1) + qgai_r1(k, j))                 &
                        & * dels + (qg1i(k, 1) + qgai(k)) * dels_r1(j)        &
                        & + 0.5 * (qg1i_rf(k, j, 1) + qgai_rf(k, j)) * dels
                qg2_r1(k, j) = (qg2i_r1(k, j, 1) + qgai_r1(k, j))                 &
                        & * dels + (qg2i(k, 1) + qgai(k)) * dels_r1(j)        &
                        & + 0.5 * (qg2i_rf(k, j, 1) + qgai_rf(k, j)) * dels
                qs1_r1(k, j) = (qs1i_r1(k, j, 1) + qsai_r1(k, j))                 &
                        & * dels + (qs1i(k, 1) + qsai(k)) * dels_r1(j)        &
                        & + 0.5 * (qs1i_rf(k, j, 1) + qsai_rf(k, j)) * dels
                qs2_r1(k, j) = (qs2i_r1(k, j, 1) + qsai_r1(k, j))                 &
                        & * dels + (qs2i(k, 1) + qsai(k)) * dels_r1(j)        &
                        & + 0.5 * (qs2i_rf(k, j, 1) + qsai_rf(k, j)) * dels
                !
                qg1_r2(k, j) = (qg1i_r2(k, j, 1) + qgai_r2(k, j))                 &
                        & * dels + (qg1i(k, 1) + qgai(k)) * dels_r2(j)        &
                        & + 0.5 * (qg1i_rf(k, j, 1) + qgai_rf(k, j)) * dels
                qg2_r2(k, j) = (qg2i_r2(k, j, 1) + qgai_r2(k, j))                 &
                        & * dels + (qg2i(k, 1) + qgai(k)) * dels_r2(j)        &
                        & + 0.5 * (qg2i_rf(k, j, 1) + qgai_rf(k, j)) * dels
                qs1_r2(k, j) = (qs1i_r2(k, j, 1) + qsai_r2(k, j))                 &
                        & * dels + (qs1i(k, 1) + qsai(k)) * dels_r2(j)        &
                        & + 0.5 * (qs1i_rf(k, j, 1) + qsai_rf(k, j)) * dels
                qs2_r2(k, j) = (qs2i_r2(k, j, 1) + qsai_r2(k, j))                 &
                        & * dels + (qs2i(k, 1) + qsai(k)) * dels_r2(j)        &
                        & + 0.5 * (qs2i_rf(k, j, 1) + qsai_rf(k, j)) * dels
            enddo
        enddo
        !
    end subroutine dlampc
    !*==RING.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! DLAMPC




    subroutine ring(xv, rv, xf, rf, ux, ur, sx, sr)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        real, parameter :: pi = 3.14159265358979
        !
        ! Dummy arguments
        !
        real :: rf, rv, sr, sx, ur, ux, xf, xv
        !
        ! Local variables
        !
        real :: ak, ele, elk, f, r, srp, x, xrm, xrp
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------------
        !     Computes the velocities induced by a vortex/source ring
        !     located at XV with radius RV with unit circulation
        !     and unit source/length density.
        !
        !     Adapted from routines provided by J. Kerwin.
        !-----------------------------------------------------------------------
        !  Input:
        !     XV,RV  x,r of the ring
        !     XF,RF  x,r of the field point
        !
        !  Output:
        !     UX,UR  velocity at XF,RF for unit circulation (+ counterclockwise)
        !     SX,SR  velocity at XF,RF for unit source/perimeter
        !-----------------------------------------------------------------------
        !
        if (rv<=0.0) then
            !------ zero-radius ring
            ux = 0.
            ur = 0.
            sx = 0.
            sr = 0.
            return
        endif
        !
        !---- this fails if R=1 and X=0  (on the ring itself)
        r = rf / rv
        x = (xv - xf) / rv
        !
        if (r==1.0 .and. x==0.0) then
            ux = 0.
            ur = 0.
            sx = 0.
            sr = 0.
            return
        endif
        !
        if (rf==0.0) then
            !----- Control Point on the axis
            ux = 1.0 / sqrt(1.0 + x**2)**3 / (2.0 * rv)
            ur = 0.0
            sx = -x / sqrt(1.0 + x**2)**3 / (2.0 * rv)
            sr = 0.0
            !
        else
            !----- Control Point not on X-axis
            xrp = x**2 + (1.0 + r)**2
            xrm = x**2 + (1.0 - r)**2
            !
            srp = sqrt(xrp)
            !
            ak = xrm / xrp
            call ellek(ak, ele, elk)
            !
            f = 2.0 / xrm
            !
            ux = (1.0 / srp) * (elk - ele * (1.0 + f * (r - 1.0))) / (2.0 * pi * rv)
            ur = (x / (srp * r)) * (elk - ele * (1.0 + f * r)) / (2.0 * pi * rv)
            !
            sx = (x / srp) * (-ele * f) / (2.0 * pi * rv)
            sr = (1.0 / (srp * r)) * (elk - ele * (1.0 + f * (r - r * r))) / (2.0 * pi * rv)
            !
            !cC----- streamfunction due to vortex
            !c       PSI = ((1.0 - 2.0*XRP*R)*ELK - ELE)*RV / (2.0*PI*SQRT(XRP))
        endif
        !
    end subroutine ring
    !*==DRING.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine dring(xv, rv, xf, rf, ux, ux_xv, ux_rv, ux_xf, ux_rf, ur, ur_xv, &
            & ur_rv, ur_xf, ur_rf, sx, sx_xv, sx_rv, sx_xf, sx_rf, sr, &
            & sr_xv, sr_rv, sr_xf, sr_rf)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        real, parameter :: pi = 3.14159265358979
        !
        ! Dummy arguments
        !
        real :: rf, rv, sr, sr_rf, sr_rv, sr_xf, sr_xv, sx, &
                & sx_rf, sx_rv, sx_xf, sx_xv, ur, ur_rf, ur_rv, &
                & ur_xf, ur_xv, ux, ux_rf, ux_rv, ux_xf, ux_xv, &
                & xf, xv
        !
        ! Local variables
        !
        real :: ak, ak_r, ak_x, dele, delk, ele, ele_r, ele_x, &
                & elk, elk_r, elk_x, f, f_r, f_x, hrip, r, rvi, &
                & r_rf, r_rv, srp, srp_r, srp_x, sr_r, sr_x, sx_r, &
                & sx_x, ur_r, ur_x, ux_r, ux_x, x, xrm, xrm_r, &
                & xrm_x, xrp, xrp_r, xrp_x, x_rv, x_xf, x_xv
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------------
        !     Same as RING, but also returns AIC derivatives w.r.t. geometry
        !-----------------------------------------------------------------------
        !
        ux = 0.
        ur = 0.
        sx = 0.
        sr = 0.
        !
        ux_x = 0.
        ur_x = 0.
        sx_x = 0.
        sr_x = 0.
        !
        ux_r = 0.
        ur_r = 0.
        sx_r = 0.
        sr_r = 0.
        !
        rvi = 0.
        !
        !----- zero-radius ring
        if (rv>0.0) then
            !
            rvi = 1.0 / rv
            !
            !---- this fails if R=1 and X=0  (on the ring itself)
            r = rf * rvi
            x = (xv - xf) * rvi
            !
            if (r/=1.0 .or. x/=0.0) then
                !
                if (rf==0.0) then
                    !----- Control Point on the axis
                    ux = 1.0 / sqrt(1.0 + x**2)**3 * pi ! / (2.0*PI*RV)
                    sx = -x / sqrt(1.0 + x**2)**3 * pi ! / (2.0*PI*RV)
                    !
                    ux_x = -3.0 * x * ux / (1.0 + x**2)
                    sx_x = -3.0 * x * ux / (1.0 + x**2) - ux
                    !
                else
                    !----- Control Point not on X-axis
                    xrp = x**2 + (1.0 + r)**2
                    xrm = x**2 + (1.0 - r)**2
                    !
                    xrp_x = 2.0 * x
                    xrm_x = 2.0 * x
                    xrp_r = 2.0 * (1.0 + r)
                    xrm_r = -2.0 * (1.0 - r)
                    !
                    srp = sqrt(xrp)
                    srp_x = (0.5 / srp) * xrp_x
                    srp_r = (0.5 / srp) * xrp_r
                    !
                    ak = xrm / xrp
                    ak_x = (xrm_x - ak * xrp_x) / xrp
                    ak_r = (xrm_r - ak * xrp_r) / xrp
                    !
                    call dellek(ak, ele, dele, elk, delk)
                    ele_x = dele * ak_x
                    ele_r = dele * ak_r
                    elk_x = delk * ak_x
                    elk_r = delk * ak_r
                    !
                    f = 2.0 / xrm
                    f_x = (-f / xrm) * xrm_x
                    f_r = (-f / xrm) * xrm_r
                    !
                    ux = (1.0 / srp) * (elk - ele * (1.0 + f * (r - 1.0)))
                    ux_x = (1.0 / srp)                                         &
                            & * (elk_x - ele_x * (1.0 + f * (r - 1.0)) - ele * (f_x * (r - 1.0)))  &
                            & - (ux / srp) * srp_x
                    ux_r = (1.0 / srp)                                         &
                            & * (elk_r - ele_r * (1.0 + f * (r - 1.0)) - ele * (f_r * (r - 1.0))   &
                                    & - ele * f) - (ux / srp) * srp_r
                    !
                    ur = (x / (srp * r)) * (elk - ele * (1.0 + f * r))
                    ur_x = (x / (srp * r)) * (elk_x - ele_x * (1.0 + f * r) - ele * (f_x * r))   &
                            & + (1.0 / (srp * r)) * (elk - ele * (1.0 + f * r)) - (ur / srp)    &
                            & * srp_x
                    ur_r = (x / (srp * r)) * (elk_r - ele_r * (1.0 + f * r) - ele * (f_r * r + f)) &
                            & - (ur / srp) * srp_r - ur / r
                    !
                    sx = (x / srp) * (-ele * f)
                    sx_x = (x / srp) * (-ele_x * f - ele * f_x) + (1.0 / srp) * (-ele * f)   &
                            & - (sx / srp) * srp_x
                    sx_r = (x / srp) * (-ele_r * f - ele * f_r) - (sx / srp) * srp_r
                    !
                    sr = (1.0 / (srp * r)) * (elk - ele * (1.0 + f * (r - r * r)))
                    sr_x = (1.0 / (srp * r))                                     &
                            & * (elk_x - ele_x * (1.0 + f * (r - r * r)) - ele * (f_x * (r - r * r)))  &
                            & - (sr / srp) * srp_x
                    sr_r = (1.0 / (srp * r))                                     &
                            & * (elk_r - ele_r * (1.0 + f * (r - r * r)) - ele * (f_r * (r - r * r)    &
                                    & + f * (1.0 - 2.0 * r))) - (sr / srp) * srp_r - sr / r
                    !
                endif
            endif
        endif
        !
        !c    X = (XV-XF)/RV
        x_xv = rvi
        x_rv = -x * rvi
        x_xf = -rvi
        !
        !c    R = RF/RV
        r_rv = -r * rvi
        r_rf = rvi
        !
        hrip = rvi / (2.0 * pi)
        !
        ux = hrip * ux
        ur = hrip * ur
        sx = hrip * sx
        sr = hrip * sr
        !
        ux_xv = hrip * ux_x * x_xv
        ur_xv = hrip * ur_x * x_xv
        sx_xv = hrip * sx_x * x_xv
        sr_xv = hrip * sr_x * x_xv
        ux_rv = hrip * (ux_x * x_rv + ux_r * r_rv) - ux * rvi
        ur_rv = hrip * (ur_x * x_rv + ur_r * r_rv) - ur * rvi
        sx_rv = hrip * (sx_x * x_rv + sx_r * r_rv) - sx * rvi
        sr_rv = hrip * (sr_x * x_rv + sr_r * r_rv) - sr * rvi
        ux_xf = hrip * ux_x * x_xf
        ur_xf = hrip * ur_x * x_xf
        sx_xf = hrip * sx_x * x_xf
        sr_xf = hrip * sr_x * x_xf
        ux_rf = hrip * ux_r * r_rf
        ur_rf = hrip * ur_r * r_rf
        sx_rf = hrip * sx_r * r_rf
        sr_rf = hrip * sr_r * r_rf
        !
    end subroutine dring
    !*==ELLEK.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine ellek(ak, ele, elk)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: ak, ele, elk
        !
        ! Local variables
        !
        real :: alk
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------------
        !              ELLIPTIC FUNCTIONS ROUTINE
        !
        !     Adapted from routines provided by J. Kerwin.
        !-----------------------------------------------------------------------
        !   Input
        !     AK     elliptic-integral argument
        !
        !   Output
        !     ELE    complete elliptic integral of the second kind
        !     ELK    complete elliptic integral of the first  kind
        !_______________________________________________________________________
        !
        alk = -log(ak)
        !
        ele = 1.00000000000 + &
                & (0.44325141463 + (0.06260601220 + (0.04757383546 + &
                        & 0.01736506451 * ak) * ak) * ak)                                   &
                        & * ak + ((0.24998368310 + (0.09200180037 + &
                & (0.04069697526 + 0.00526449639 * ak) * ak) * ak) * ak) * alk
        !
        elk = 1.38629436112 + &
                & (0.09666344259 + (0.03590092383 + (0.03742563713 + &
                        & 0.01451196212 * ak) * ak) * ak)                                   &
                        & * ak + (0.50000000000 + (0.12498593597 + &
                & (0.06880248576 + (0.03328355346 + 0.00441787012 * ak) * ak) * ak) * ak) &
                & * alk
        !
    end subroutine ellek
    !*==DELLEK.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine dellek(ak, ele, ele_ak, elk, elk_ak)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: ak, ele, ele_ak, elk, elk_ak
        !
        ! Local variables
        !
        real :: alk
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------------
        !              ELLIPTIC FUNCTIONS + DERIVATIVE ROUTINE
        !
        !     Adapted from routines provided by J. Kerwin.
        !-----------------------------------------------------------------------
        !   Input
        !     AK     elliptic-integral argument
        !
        !   Output
        !     ELE     complete elliptic integral of the second kind
        !     ELK     complete elliptic integral of the first  kind
        !     ELE_AK  d(ELE)/d(AK)
        !     ELK_AK  d(ELK)/d(AK)
        !_______________________________________________________________________
        !
        alk = -log(ak)
        !
        ele = 1.00000000000 + &
                & (0.44325141463 + (0.06260601220 + (0.04757383546 + &
                        & 0.01736506451 * ak) * ak) * ak)                                   &
                        & * ak + ((0.24998368310 + (0.09200180037 + &
                & (0.04069697526 + 0.00526449639 * ak) * ak) * ak) * ak) * alk
        !
        elk = 1.38629436112 + &
                & (0.09666344259 + (0.03590092383 + (0.03742563713 + &
                        & 0.01451196212 * ak) * ak) * ak)                                   &
                        & * ak + (0.50000000000 + (0.12498593597 + &
                & (0.06880248576 + (0.03328355346 + 0.00441787012 * ak) * ak) * ak) * ak) &
                & * alk
        !
        ele_ak = 0.19326773153 + &
                & (0.03321022403 + (0.10202453112 + 0.06419576165 * ak) * ak)      &
                        & * ak + (0.24998368310 + &
                & (0.18400360074 + (0.12209092578 + 0.02105798556 * ak) * ak) * ak)  &
                & * alk
        !
        elk_ak = -0.028322493380 + &
                & (0.002999361900 + (0.078993357930 + 0.053629978360 * ak) * ak)   &
                        & * ak - 0.50000000000 / ak + &
                & (0.12498593597 + (0.13760497152 + (0.09985066038 + &
                        & 0.01767148048 * ak) * ak) * ak) * alk
        !
    end subroutine dellek
end module m_lamp
