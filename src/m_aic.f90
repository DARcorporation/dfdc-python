module m_aic
    implicit none
contains
    !*==PANAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
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

    subroutine panaic(jfrst, jlast, xp, yp, nf, xf, yf, jpan1, jpan2, iftype, &
            & gam, sig, jdim, qf, qf_gam, qf_sig)
        use m_lamp, only : lamp, lampc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: jdim, jfrst, jlast, nf
        real, dimension(*) :: gam, sig, xf, xp, yf, yp
        integer, dimension(*) :: iftype, jpan1, jpan2
        real, dimension(2, *) :: qf
        real, dimension(2, jdim, *) :: qf_gam, qf_sig
        !
        ! Local variables
        !
        real :: anwx, anwy, aswx, aswy, ug1, ug2, us1, us2, &
                & vg1, vg2, vs1, vs2
        integer :: i, j, j1, j2
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Computes field-point velocities due to vortex/source panels.
        !     Computes the velocities' jacobians w.r.t. the panel strengths.
        !
        !  Input:
        !     JFRST     first panel node
        !     JLAST     last panel node
        !     XP,YP(j)  panel node x,y locations
        !
        !     NF        number of field points
        !     XF,YF(i)  field-point x locations  i = 1..NF
        !     JPAN1(i)  j indices of the two nodes defining the panel
        !     JPAN2(i)      which contains field point i at its midpoint
        !     IFTYPE(i) dummy field point if IFTYPE < 0  (skip)
        !
        !     GAM(j)   sheet vortex strength at X,Y(j)  (+ clockwise)
        !     SIG(j)   sheet source strength at X,Y(j)
        !
        !     JDIM     Second dimension of QF_XXX arrays
        !
        !  Output:
        !     QF(1:2,i)     x,y-velocities at XF,YF(i)
        !     QF_GAM(.j.)   dQF(..)/dGAM(j)
        !     QF_SIG(.j.)   dQF(..)/dSIG(j)
        !----------------------------------------------------------------------
        !
        !---- unit vectors tangent and normal to wall, if any
        aswx = 1.
        aswy = 0.
        !
        anwx = -aswy
        anwy = aswx
        !
        do i = 1, nf
            if (iftype(i)>=0) then
                !
                do j = jfrst, jlast - 1
                    j1 = j
                    j2 = j + 1
                    !
                    if (j1==jpan1(i) .and. j2==jpan2(i)) then
                        !--------- field point is at the center of the J..J+1 panel
                        call lampc(xp(j1), yp(j1), xp(j2), yp(j2), ug1, vg1, ug2, &
                                & vg2, us1, vs1, us2, vs2)
                    else
                        !--------- field point is somewhere else
                        call lamp(xp(j1), yp(j1), xp(j2), yp(j2), xf(i), yf(i), ug1, &
                                & vg1, ug2, vg2, us1, vs1, us2, vs2)
                    endif
                    !
                    !-------- accumulate velocity components and their sensitivities
                    !         Note:  GAM is defined positive clockwise
                    !                PANE,LAMP assume positive counterclockwise
                    qf(1, i) = qf(1, i) - ug1 * gam(j1) - ug2 * gam(j2)            &
                            & + us1 * sig(j1) + us2 * sig(j2)
                    !
                    qf(2, i) = qf(2, i) - vg1 * gam(j1) - vg2 * gam(j2)            &
                            & + vs1 * sig(j1) + vs2 * sig(j2)
                    !
                    qf_gam(1, j1, i) = qf_gam(1, j1, i) - ug1
                    qf_gam(1, j2, i) = qf_gam(1, j2, i) - ug2
                    qf_sig(1, j1, i) = qf_sig(1, j1, i) + us1
                    qf_sig(1, j2, i) = qf_sig(1, j2, i) + us2
                    !
                    qf_gam(2, j1, i) = qf_gam(2, j1, i) - vg1
                    qf_gam(2, j2, i) = qf_gam(2, j2, i) - vg2
                    qf_sig(2, j1, i) = qf_sig(2, j1, i) + vs1
                    qf_sig(2, j2, i) = qf_sig(2, j2, i) + vs2
                enddo
            endif
        enddo
        !
    end subroutine panaic
    !*==LINAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine linaic(jfrst, jlast, xp, yp, nf, xf, yf, vor, sou, jdim, qf, &
            & qf_vor, qf_sou)
        use m_lamp, only : ring
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: jdim, jfrst, jlast, nf
        real, dimension(2, *) :: qf
        real, dimension(2, jdim, *) :: qf_sou, qf_vor
        real, dimension(*) :: sou, vor, xf, xp, yf, yp
        !
        ! Local variables
        !
        real :: anwx, anwy, aswx, aswy, ug, us, vg, vs
        integer :: i, j
        real, save :: pi
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------
        !     Generates line-singularity AIC matrices.
        !
        !  Input:
        !     JFRST
        !     JLAST
        !     XP,YP(j)  Singularity x,y locations       j = JFRST..JLAST
        !
        !     NF        Number of field points
        !     XF,YF(i)  Field-point x,y locations  i = 1..NF
        !
        !     VOR      Vortex line (or ring) strength (+ clockwise)
        !     SOU      Source line (or ring) strength
        !
        !     JDIM     Second dimension of VIC array
        !
        !  Output:
        !     QF(1:2,i)     x,y-velocities at XF,YF(i)
        !     QF_VOR(.j.)   dQF(..)/dVOR(j)
        !     QF_SOU(.j.)   dQF(..)/dSOU(j)
        !--------------------------------------------------------
        data pi/3.1415926535897932384/
        !
        !---- unit vectors tangent and normal to wall, if any
        aswx = 0.
        aswy = 1.
        anwx = 1.
        anwy = 0.
        !
        do i = 1, nf
            do j = jfrst, jlast
                call ring(xp(j), yp(j), xf(i), yf(i), ug, vg, us, vs)
                !
                !-------- accumulate velocity components and their sensitivities
                !-        Note:  GAM is defined positive clockwise
                !-               PANE,LAMP assume positive counterclockwise
                qf(1, i) = qf(1, i) - ug * vor(j) + us * sou(j)
                qf(2, i) = qf(2, i) - vg * vor(j) + vs * sou(j)
                !
                qf_vor(1, j, i) = qf_vor(1, j, i) - ug
                qf_sou(1, j, i) = qf_sou(1, j, i) + us
                !
                qf_vor(2, j, i) = qf_vor(2, j, i) - vg
                qf_sou(2, j, i) = qf_sou(2, j, i) + vs
            enddo
        enddo
        !
    end subroutine linaic
    !*==AXLAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! LINAIC



    subroutine axlaic(jfrst, jlast, xp, yp, nf, xf, yf, jpan1, jpan2, iftype, &
            & dbl, src, rcore, jdim, qf, qf_dbl, qf_src)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: jdim, jfrst, jlast, nf
        real, dimension(*) :: dbl, rcore, src, xf, xp, yf, yp
        integer, dimension(*) :: iftype, jpan1, jpan2
        real, dimension(2, *) :: qf
        real, dimension(2, jdim, *) :: qf_dbl, qf_src
        !
        ! Local variables
        !
        real :: dx, eps, eps1, eps2, r1, r2, rcore1, rcore2, &
                & rmx1, rmx2, rsq1, rsq2, ug1, ug2, us1, us2, &
                & vg1, vg2, vs1, vs2, x1, x2, y
        integer :: i, j, j1, j2
        real, save :: pi
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------
        !     Generates axis-line singularity AIC matrices.
        !
        !  Input:
        !     JFRST
        !     JLAST
        !     XP,YP(j)  Singularity x,r locations       j = JFRST, JLAST
        !
        !     NF        Number of field points
        !     XF,YF(i)  Field-point x,r locations  i = 1..NF
        !
        !     DBL(j)    Line-doublet strength (positive when flow is to -x on axis)
        !     SRC(j)    Line-source  strength
        !
        !     JDIM      Second dimension of VIC array
        !
        !  Output:
        !     QF(1:2,i)     x,r-velocities at XF,YF(i)
        !     QF_DBL(.j.)   dQF(..)/dDBL(j)
        !     QF_SRC(.j.)   dQF(..)/dSRC(j)
        !
        !     Special case:  if YF(i)=0, then QF(i) is actually q*r
        !--------------------------------------------------------
        data pi/3.1415926535897932384/
        !
        do i = 1, nf
            if (iftype(i)>=0) then
                !
                do j = jfrst, jlast - 1
                    j1 = j
                    j2 = j + 1
                    !
                    dx = xp(j2) - xp(j1)
                    !
                    x1 = xf(i) - xp(j1)
                    x2 = xf(i) - xp(j2)
                    y = yf(i)
                    !
                    rcore1 = rcore(j1)
                    rcore2 = rcore(j2)
                    !
                    !cc       IF(J1.EQ.JPAN1(I) .AND. J2.EQ.JPAN2(I)) THEN
                    !--------- field point is at the center of the J..J+1 line segment
                    rsq1 = x1 * x1 + y * y + 0.25 * rcore1**2
                    rsq2 = x2 * x2 + y * y + 0.25 * rcore2**2
                    !          ELSE
                    !C--------- field point is somewhere else
                    !           RSQ1 = X1*X1 + Y*Y
                    !           RSQ2 = X2*X2 + Y*Y
                    !          ENDIF
                    !
                    r1 = sqrt(rsq1)
                    r2 = sqrt(rsq2)
                    !
                    eps1 = y * y / rsq1
                    eps2 = y * y / rsq2
                    !
                    eps = min(eps1, eps2)
                    !
                    !-------- compute r-x using exact expression,
                    !-          or 4th-order Taylor series about Y=0 to avoid roundoff with sqrt()
                    !          IF(EPS .GT. 0.02) THEN
                    rmx1 = r1 - x1
                    rmx2 = r2 - x2
                    !          ELSE
                    !           RSG1 = SIGN(R1,X1)
                    !           RSG2 = SIGN(R2,X2)
                    !           RMX1 = (R1 - RSG1) + RSG1*EPS1*(0.5 + 0.125*EPS1)
                    !           RMX2 = (R2 - RSG2) + RSG2*EPS2*(0.5 + 0.125*EPS2)
                    !          ENDIF
                    !
                    !          IF(MIN(EPS1,EPS2) .LT. 0.0001
                    !     &       .AND. RSG1.GT.0.0
                    !     &       .AND. RSG2.GT.0.0) THEN
                    !           RMX1 = R2
                    !           RMX2 = R1
                    !          ENDIF
                    !C
                    us1 = -(1.0 / r1 + log(rmx1 / rmx2) / dx) / (4.0 * pi)
                    us2 = (1.0 / r2 + log(rmx1 / rmx2) / dx) / (4.0 * pi)
                    vs1 = y * ((1.0 / r1 - 1.0 / dx) / rmx1 + (1.0 / dx) / rmx2) / (4.0 * pi)
                    vs2 = y * ((1.0 / dx) / rmx1 - (1.0 / r2 + 1.0 / dx) / rmx2) / (4.0 * pi)
                    !
                    !-------- no doublet line currently implemented
                    ug1 = 0.
                    ug2 = 0.
                    vg1 = 0.
                    vg2 = 0.
                    !
                    !-------- accumulate velocity components and their sensitivities
                    qf(1, i) = qf(1, i) - ug1 * dbl(j1) - ug2 * dbl(j2)            &
                            & + us1 * src(j1) + us2 * src(j2)
                    !
                    qf(2, i) = qf(2, i) - vg1 * dbl(j1) - vg2 * dbl(j2)            &
                            & + vs1 * src(j1) + vs2 * src(j2)
                    !
                    qf_dbl(1, j1, i) = qf_dbl(1, j1, i) - ug1
                    qf_dbl(1, j2, i) = qf_dbl(1, j2, i) - ug2
                    qf_src(1, j1, i) = qf_src(1, j1, i) + us1
                    qf_src(1, j2, i) = qf_src(1, j2, i) + us2
                    !
                    qf_dbl(2, j1, i) = qf_dbl(2, j1, i) - vg1
                    qf_dbl(2, j2, i) = qf_dbl(2, j2, i) - vg2
                    qf_src(2, j1, i) = qf_src(2, j1, i) + vs1
                    qf_src(2, j2, i) = qf_src(2, j2, i) + vs2
                enddo
            endif
        enddo
        !
    end subroutine axlaic
    !*==PNTAIC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! AXLAIC



    subroutine pntaic(jfrst, jlast, xp, yp, nf, xf, yf, dbl, src, jdim, qf, &
            & qf_dbl, qf_src)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: jdim, jfrst, jlast, nf
        real, dimension(*) :: dbl, src, xf, xp, yf, yp
        real, dimension(2, *) :: qf
        real, dimension(2, jdim, *) :: qf_dbl, qf_src
        !
        ! Local variables
        !
        real :: dx, dy, r32, rsq, ug, us, vg, vs
        integer :: i, j
        real, save :: pi
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------
        !     Generates point-source and point-doublet AIC matrices
        !     (axisymmetric case only)
        !
        !  Input:
        !     JFRST
        !     JLAST
        !     XP,YP(j)  Singularity x,r locations       j = JFRST..JLAST
        !
        !     NF        Number of field points
        !     XF,YF(i)  Field-point x,r locations  i = 1..NF
        !
        !     DBL       Point-doublet strength (positive when flow is to -x on axis)
        !     SRC       Point-source  strength
        !
        !     JDIM      Second dimension of VIC array
        !
        !  Output:
        !     QF(1:2,i)     x,r-velocities at XF,YF(i)
        !     QF_DBL(.j.)   dQF(..)/dDBL(j)
        !     QF_SRC(.j.)   dQF(..)/dSRC(j)
        !--------------------------------------------------------
        data pi/3.1415926535897932384/
        !
        do i = 1, nf
            do j = jfrst, jlast
                !-------- point singularity on axis... DBL is x-doublet, SRC is point source
                dx = xf(i) - xp(j)
                dy = yf(i)
                rsq = dx * dx + dy * dy
                r32 = sqrt(rsq)**3
                ug = (2.0 * dx * dx - dy * dy) / (4.0 * pi * r32 * rsq)
                vg = (3.0 * dx * dy) / (4.0 * pi * r32 * rsq)
                us = dx / (4.0 * pi * r32)
                vs = dy / (4.0 * pi * r32)
                !
                qf(1, i) = qf(1, i) - ug * dbl(j) + us * src(j)
                qf(2, i) = qf(2, i) - vg * dbl(j) + vs * src(j)
                !
                qf_dbl(1, j, i) = qf_dbl(1, j, i) - ug
                qf_dbl(2, j, i) = qf_dbl(2, j, i) - vg
                qf_src(1, j, i) = qf_src(1, j, i) + us
                qf_src(2, j, i) = qf_src(2, j, i) + vs
            enddo
        enddo
        !
    end subroutine pntaic
    !*==PANAICD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PNTAIC



    subroutine panaicd(jfrst, jlast, xp, yp, nf, xf, yf, jpan1, jpan2, iftype, &
            & gam, sig, jdim, qf, qf_gam, qf_sig, qf_xp, qf_yp, &
            & qf_xf, qf_yf)
        use m_lamp, only : dlampc, dlamp
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: jdim, jfrst, jlast, nf
        real, dimension(*) :: gam, sig, xf, xp, yf, yp
        integer, dimension(*) :: iftype, jpan1, jpan2
        real, dimension(2, *) :: qf, qf_xf, qf_yf
        real, dimension(2, jdim, *) :: qf_gam, qf_sig, qf_xp, qf_yp
        !
        ! Local variables
        !
        real :: anwx, anwy, aswx, aswy
        integer :: i, j, j1, j2, kq, kr
        real, dimension(2) :: qg1, qg2, qs1, qs2, r1, r2, rf
        real, dimension(2, 2) :: qg1_r1, qg1_r2, qg1_rf, qg2_r1, &
                & qg2_r2, qg2_rf, qs1_r1, qs1_r2, &
                & qs1_rf, qs2_r1, qs2_r2, qs2_rf
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Sames as PANAIC, but also returns Jacobians w.r.t. the geometry
        !----------------------------------------------------------------------
        !
        !---- local arrays for calling DLAMP
        !
        !---- unit vectors tangent and normal to wall, if any
        aswx = 1.
        aswy = 0.
        !
        anwx = -aswy
        anwy = aswx
        !
        do i = 1, nf
            if (iftype(i)>=0) then
                !
                rf(1) = xf(i)
                rf(2) = yf(i)
                !
                do j = jfrst, jlast - 1
                    j1 = j
                    j2 = j + 1
                    !
                    if (j1==jpan1(i) .and. j2==jpan2(i)) then
                        !--------- field point is at the center of the J..J+1 panel
                        r1(1) = xp(j1)
                        r1(2) = yp(j1)
                        r2(1) = xp(j2)
                        r2(2) = yp(j2)
                        call dlampc(r1, r2, qg1, qg1_r1, qg1_r2, qg2, qg2_r1, qg2_r2, &
                                & qs1, qs1_r1, qs1_r2, qs2, qs2_r1, qs2_r2)
                        do kq = 1, 2
                            do kr = 1, 2
                                qg1_rf(kq, kr) = 0.
                                qg2_rf(kq, kr) = 0.
                                qs1_rf(kq, kr) = 0.
                                qs2_rf(kq, kr) = 0.
                            enddo
                        enddo
                    else
                        !--------- field point is somewhere else
                        r1(1) = xp(j1)
                        r1(2) = yp(j1)
                        r2(1) = xp(j2)
                        r2(2) = yp(j2)
                        call dlamp(r1, r2, rf, qg1, qg1_r1, qg1_r2, qg1_rf, qg2, &
                                & qg2_r1, qg2_r2, qg2_rf, qs1, qs1_r1, qs1_r2, &
                                & qs1_rf, qs2, qs2_r1, qs2_r2, qs2_rf)
                    endif
                    !
                    !
                    !-------- accumulate velocity components and their sensitivities
                    !-        Note:  GAM is defined positive clockwise
                    !-               PANE,LAMP assume positive counterclockwise
                    do kq = 1, 2
                        qf(kq, i) = qf(kq, i) - qg1(kq) * gam(j1) - qg2(kq)       &
                                & * gam(j2) + qs1(kq) * sig(j1) + qs2(kq)       &
                                & * sig(j2)
                        !
                        qf_gam(kq, j1, i) = qf_gam(kq, j1, i) - qg1(kq)
                        qf_gam(kq, j2, i) = qf_gam(kq, j2, i) - qg2(kq)
                        qf_sig(kq, j1, i) = qf_sig(kq, j1, i) + qs1(kq)
                        qf_sig(kq, j2, i) = qf_sig(kq, j2, i) + qs2(kq)
                        !
                        qf_xp(kq, j1, i) = qf_xp(kq, j1, i) - qg1_r1(kq, 1) * gam(j1)&
                                & - qg2_r1(kq, 1) * gam(j2) + qs1_r1(kq, 1)&
                                & * sig(j1) + qs2_r1(kq, 1) * sig(j2)
                        qf_xp(kq, j2, i) = qf_xp(kq, j2, i) - qg1_r2(kq, 1) * gam(j1)&
                                & - qg2_r2(kq, 1) * gam(j2) + qs1_r2(kq, 1)&
                                & * sig(j1) + qs2_r2(kq, 1) * sig(j2)
                        qf_xf(kq, i) = qf_xf(kq, i) - qg1_rf(kq, 1) * gam(j1)      &
                                & - qg2_rf(kq, 1) * gam(j2) + qs1_rf(kq, 1)   &
                                & * sig(j1) + qs2_rf(kq, 1) * sig(j2)
                        !
                        qf_yp(kq, j1, i) = qf_yp(kq, j1, i) - qg1_r1(kq, 2) * gam(j1)&
                                & - qg2_r1(kq, 2) * gam(j2) + qs1_r1(kq, 2)&
                                & * sig(j1) + qs2_r1(kq, 2) * sig(j2)
                        qf_yp(kq, j2, i) = qf_yp(kq, j2, i) - qg1_r2(kq, 2) * gam(j1)&
                                & - qg2_r2(kq, 2) * gam(j2) + qs1_r2(kq, 2)&
                                & * sig(j1) + qs2_r2(kq, 2) * sig(j2)
                        qf_yf(kq, i) = qf_yf(kq, i) - qg1_rf(kq, 2) * gam(j1)      &
                                & - qg2_rf(kq, 2) * gam(j2) + qs1_rf(kq, 2)   &
                                & * sig(j1) + qs2_rf(kq, 2) * sig(j2)
                    enddo
                enddo
            endif
        enddo
        !
    end subroutine panaicd
    !*==LINAICD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! PANAICD



    subroutine linaicd(jfrst, jlast, xp, yp, nf, xf, yf, vor, sou, jdim, qf, &
            & qf_vor, qf_sou, qf_xp, qf_yp, qf_xf, qf_yf)
        use m_lamp, only : dring
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: jdim, jfrst, jlast, nf
        real, dimension(2, *) :: qf, qf_xf, qf_yf
        real, dimension(2, jdim, *) :: qf_sou, qf_vor, qf_xp, qf_yp
        real, dimension(*) :: sou, vor, xf, xp, yf, yp
        !
        ! Local variables
        !
        real :: anwx, anwy, aswx, aswy, ug, ug_xf, ug_xp, ug_yf, &
                & ug_yp, us, us_xf, us_xp, us_yf, us_yp, vg, &
                & vg_xf, vg_xp, vg_yf, vg_yp, vs, vs_xf, vs_xp, &
                & vs_yf, vs_yp
        integer :: i, j
        real, save :: pi
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Sames as LINAIC, but also returns Jacobians w.r.t. the geometry
        !----------------------------------------------------------------------
        data pi/3.1415926535897932384/
        !
        !---- unit vectors tangent and normal to wall, if any
        aswx = 0.
        aswy = 1.
        anwx = 1.
        anwy = 0.
        !
        do i = 1, nf
            do j = jfrst, jlast
                call dring(xp(j), yp(j), xf(i), yf(i), ug, ug_xp, ug_yp, ug_xf, &
                        & ug_yf, vg, vg_xp, vg_yp, vg_xf, vg_yf, us, us_xp, us_yp, &
                        & us_xf, us_yf, vs, vs_xp, vs_yp, vs_xf, vs_yf)
                !
                !-------- accumulate velocity components and their sensitivities
                !-        Note:  GAM is defined positive clockwise
                !-               PANE,LAMP assume positive counterclockwise
                qf(1, i) = qf(1, i) - ug * vor(j) + us * sou(j)
                qf(2, i) = qf(2, i) - vg * vor(j) + vs * sou(j)
                !
                qf_vor(1, i, j) = qf_vor(1, i, j) - ug
                qf_sou(1, i, j) = qf_sou(1, i, j) + us
                !
                qf_vor(2, i, j) = qf_vor(2, i, j) - vg
                qf_sou(2, i, j) = qf_sou(2, i, j) + vs
                !
                qf_xp(1, i, j) = qf_xp(1, i, j) - ug_xp * vor(j) + us_xp * sou(j)
                qf_yp(1, i, j) = qf_yp(1, i, j) - ug_yp * vor(j) + us_yp * sou(j)
                qf_xf(1, i) = qf_xf(1, i) - ug_xf * vor(j) + us_xf * sou(j)
                qf_yf(1, i) = qf_yf(1, i) - ug_yf * vor(j) + us_yf * sou(j)
                !
                qf_xp(2, i, j) = qf_xp(2, i, j) - vg_xp * vor(j) + vs_xp * sou(j)
                qf_yp(2, i, j) = qf_yp(2, i, j) - vg_yp * vor(j) + vs_yp * sou(j)
                qf_xf(2, i) = qf_xf(2, i) - vg_xf * vor(j) + vs_xf * sou(j)
                qf_yf(2, i) = qf_yf(2, i) - vg_yf * vor(j) + vs_yf * sou(j)
                !
            enddo
        enddo
        !
    end subroutine linaicd
    !*==PNTAICD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! LINAICD



    subroutine pntaicd(jfrst, jlast, xp, yp, nf, xf, yf, dbl, src, jdim, qf, &
            & qf_dbl, qf_src, qf_xp, qf_yp, qf_xf, qf_yf)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: jdim, jfrst, jlast, nf
        real, dimension(*) :: dbl, src, xf, xp, yf, yp
        real, dimension(2, *) :: qf, qf_xf, qf_yf
        real, dimension(2, jdim, *) :: qf_dbl, qf_src, qf_xp, qf_yp
        !
        ! Local variables
        !
        real :: dx, dy, r32, r52, rsq, ug, ug_xf, ug_xp, ug_yf, &
                & ug_yp, us, us_xf, us_xp, us_yf, us_yp, vg, &
                & vg_xf, vg_xp, vg_yf, vg_yp, vs, vs_xf, vs_xp, &
                & vs_yf, vs_yp
        integer :: i, j
        real, save :: pi
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Sames as PNTAIC, but also returns Jacobians w.r.t. the geometry
        !----------------------------------------------------------------------
        data pi/3.1415926535897932384/
        !
        do i = 1, nf
            do j = jfrst, jlast
                !-------- point singularity on axis... DBL is x-doublet, SRC is point source
                dx = xf(i) - xp(j)
                dy = yf(i)
                !
                rsq = dx * dx + dy * dy
                r32 = sqrt(rsq)**3
                r52 = r32 * rsq
                !
                ug = (2.0 * dx * dx - dy * dy) / (4.0 * pi * r52)
                vg = (3.0 * dx * dy) / (4.0 * pi * r52)
                us = dx / (4.0 * pi * r32)
                vs = dy / (4.0 * pi * r32)
                !
                ug_xf = 4.0 * dx / (4.0 * pi * r52) - 5.0 * dx * ug / rsq
                ug_yf = -2.0 * dy / (4.0 * pi * r52) - 5.0 * dy * ug / rsq
                vg_xf = 3.0 * dy / (4.0 * pi * r52) - 5.0 * dx * vg / rsq
                vg_yf = 3.0 * dx / (4.0 * pi * r52) - 5.0 * dy * vg / rsq
                !
                us_xf = 1.0 / (4.0 * pi * r32) - 3.0 * dx * us / rsq
                us_yf = -3.0 * dy * us / rsq
                vs_xf = -3.0 * dx * vs / rsq
                vs_yf = 1.0 / (4.0 * pi * r32) - 3.0 * dy * vs / rsq
                !
                ug_xp = -ug_xf
                ug_yp = 0.
                vg_xp = -vg_xf
                vg_yp = 0.
                !
                us_xp = -us_xf
                us_yp = 0.
                vs_xp = -vs_xf
                vs_yp = 0.
                !
                !
                qf(1, i) = qf(1, i) - ug * dbl(j) + us * src(j)
                qf(2, i) = qf(2, i) - vg * dbl(j) + vs * src(j)
                !
                qf_dbl(1, i, j) = qf_dbl(1, i, j) - ug
                qf_dbl(2, i, j) = qf_dbl(2, i, j) - vg
                qf_src(1, i, j) = qf_src(1, i, j) + us
                qf_src(2, i, j) = qf_src(2, i, j) + vs
                !
                qf_xp(1, i, j) = qf_xp(1, i, j) - ug_xp * dbl(j) + us_xp * src(j)
                qf_yp(1, i, j) = qf_yp(1, i, j) - ug_yp * dbl(j) + us_yp * src(j)
                qf_xf(1, i) = qf_xf(1, i) - ug_xf * dbl(j) + us_xf * src(j)
                qf_yf(1, i) = qf_yf(1, i) - ug_yf * dbl(j) + us_yf * src(j)
                !
                qf_xp(2, i, j) = qf_xp(2, i, j) - vg_xp * dbl(j) + vs_xp * src(j)
                qf_yp(2, i, j) = qf_yp(2, i, j) - vg_yp * dbl(j) + vs_yp * src(j)
                qf_xf(2, i) = qf_xf(2, i) - vg_xf * dbl(j) + vs_xf * src(j)
                qf_yf(2, i) = qf_yf(2, i) - vg_yf * dbl(j) + vs_yf * src(j)
                !
            enddo
        enddo
        !
    end subroutine pntaicd
end module m_aic
