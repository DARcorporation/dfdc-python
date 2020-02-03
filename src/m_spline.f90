module m_spline
    implicit none
contains
    !*==SPLINE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GAMSOL
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
    !    1-D Spline Package.
    !    Interpolates a function x(s) over an s interval
    !
    !                                              Mark Drela
    !                                              1985
    !
    !    Usage:
    !
    !C---- fill S(i), X(i) arrays
    !      S(i) = ...
    !      X(i) = ...
    !
    !C---- calculate spline coefficients for x(s) from discrete values
    !C-    (can use SPLIND or SPLINA instead)
    !      CALL SPLINE(X,XS,S,N)
    !
    !C---- evaluate splined x(s) and/or its derivatives
    !C-    at any number of s points SS
    !      XX   = SEVAL(SS,X,XS,S,N)
    !      XXS  = DEVAL(SS,X,XS,S,N)
    !      XXSS = D2VAL(SS,X,XS,S,N)
    !
    !C---- alternative to calling SEVAL,DEVAL,D2VAL separately
    !C-    (slightly more efficient if all three quantities are needed)
    !      CALL SEVALL(SS,X,XS,S,N, XX,XXS,XXSS)
    !
    !


    subroutine spline(x, xs, s, n)
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
        integer :: n
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        real, dimension(nmax) :: a, b, c
        real :: dsm, dsp
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Calculates spline coefficients for X(S).          |
        !     Natural end conditions are used (zero 3rd         |
        !      derivative over first, last intervals).          |
        !                                                       |
        !     To evaluate the spline at some value of S,        |
        !     use SEVAL and/or DEVAL.                           |
        !                                                       |
        !     S        independent variable array (input)       |
        !     X        dependent variable array   (input)       |
        !     XS       dX/dS array                (calculated)  |
        !     N        number of points           (input)       |
        !                                                       |
        !-------------------------------------------------------
        if (n>nmax) stop 'SPLINE: array overflow, increase NMAX'
        !
        if (n==1) then
            xs(1) = 0.
            return
        endif
        !
        do i = 2, n - 1
            dsm = s(i) - s(i - 1)
            dsp = s(i + 1) - s(i)
            b(i) = dsp
            a(i) = 2.0 * (dsm + dsp)
            c(i) = dsm
            xs(i) = 3.0 * ((x(i + 1) - x(i)) * dsm / dsp + (x(i) - x(i - 1)) * dsp / dsm)
        enddo
        !
        !---- set zero 3rd derivative end conditions
        a(1) = 1.0
        c(1) = 1.0
        xs(1) = 2.0 * (x(2) - x(1)) / (s(2) - s(1))
        !
        b(n) = 1.0
        a(n) = 1.0
        xs(n) = 2.0 * (x(n) - x(n - 1)) / (s(n) - s(n - 1))
        !
        if (n==2) then
            !----- if only two points are present, specify zero 2nd derivative instead
            !-     (straight line interpolation will result)
            b(n) = 1.0
            a(n) = 2.0
            xs(n) = 3.0 * (x(n) - x(n - 1)) / (s(n) - s(n - 1))
        endif
        !
        !---- solve for derivative array XS
        call trisol(a, b, c, xs, n)
        !
    end subroutine spline
    !*==SPLIND.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SPLINE



    subroutine splind(x, xs, s, n, xs1, xs2)
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
        integer :: n
        real :: xs1, xs2
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        real, dimension(nmax) :: a, b, c
        real :: dsm, dsp
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Calculates spline coefficients for X(S).          |
        !     Same as SPLINE, but also allows specified-slope   |
        !     or zero-curvature end conditions to be imposed.   |
        !                                                       |
        !     To evaluate the spline at some value of S,        |
        !     use SEVAL and/or DEVAL.                           |
        !                                                       |
        !     S        independent variable array (input)       |
        !     X        dependent variable array   (input)       |
        !     XS       dX/dS array                (calculated)  |
        !     N        number of points           (input)       |
        !     XS1,XS2  endpoint derivatives       (input)       |
        !              If =  999.0, then usual zero second      |
        !              derivative end condition(s) are used     |
        !              If = -999.0, then zero third             |
        !              derivative end condition(s) are used     |
        !                                                       |
        !     Note: specifying both XS1,XS2 = -999.0            |
        !           is equivalent to using SPLINE.              |
        !                                                       |
        !-------------------------------------------------------
        if (n>nmax) stop 'SPLIND: array overflow, increase NMAX'
        !
        if (n==1) then
            xs(1) = 0.
            return
        endif
        !
        do i = 2, n - 1
            dsm = s(i) - s(i - 1)
            dsp = s(i + 1) - s(i)
            b(i) = dsp
            a(i) = 2.0 * (dsm + dsp)
            c(i) = dsm
            xs(i) = 3.0 * ((x(i + 1) - x(i)) * dsm / dsp + (x(i) - x(i - 1)) * dsp / dsm)
        enddo
        !
        if (xs1==999.0) then
            !----- set zero second derivative end condition
            a(1) = 2.0
            c(1) = 1.0
            xs(1) = 3.0 * (x(2) - x(1)) / (s(2) - s(1))
        elseif (xs1==-999.0) then
            !----- set zero third derivative end condition
            a(1) = 1.0
            c(1) = 1.0
            xs(1) = 2.0 * (x(2) - x(1)) / (s(2) - s(1))
        else
            !----- set specified first derivative end condition
            a(1) = 1.0
            c(1) = 0.
            xs(1) = xs1
        endif
        !
        if (xs2==999.0) then
            b(n) = 1.0
            a(n) = 2.0
            xs(n) = 3.0 * (x(n) - x(n - 1)) / (s(n) - s(n - 1))
        elseif (xs2==-999.0) then
            b(n) = 1.0
            a(n) = 1.0
            xs(n) = 2.0 * (x(n) - x(n - 1)) / (s(n) - s(n - 1))
        else
            a(n) = 1.0
            b(n) = 0.
            xs(n) = xs2
        endif
        !
        if (n==2 .and. xs1==-999.0 .and. xs2==-999.0) then
            !----- cannot have zero 3rd derivative at both endpoints of a single interval
            b(n) = 1.0
            a(n) = 2.0
            xs(n) = 3.0 * (x(n) - x(n - 1)) / (s(n) - s(n - 1))
        endif
        !
        !---- solve for derivative array XS
        call trisol(a, b, c, xs, n)
        !
    end subroutine splind
    !*==SPLINA.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SPLIND


    subroutine splina(x, xs, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        real :: ds, dx, xs1, xs2
        integer :: i
        logical :: lend
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Calculates spline coefficients for X(S) by a      |
        !     simple averaging of adjacent segment slopes.      |
        !                                                       |
        !     Interpolated X(S) is less likely to oscillate     |
        !     than with SPLINE, but does not have continuity    |
        !     in curvature.                                     |
        !                                                       |
        !     To evaluate the spline at some value of S,        |
        !     use SEVAL and/or DEVAL.                           |
        !                                                       |
        !     S        independent variable array (input)       |
        !     X        dependent variable array   (input)       |
        !     XS       dX/dS array                (calculated)  |
        !     N        number of points           (input)       |
        !                                                       |
        !-------------------------------------------------------
        !
        if (n==1) then
            xs(1) = 0.
            return
        endif
        !
        lend = .true.
        do i = 1, n - 1
            ds = s(i + 1) - s(i)
            if (ds==0.) then
                xs(i) = xs1
                lend = .true.
            else
                dx = x(i + 1) - x(i)
                xs2 = dx / ds
                if (lend) then
                    xs(i) = xs2
                    lend = .false.
                else
                    xs(i) = 0.5 * (xs1 + xs2)
                endif
            endif
            xs1 = xs2
        enddo
        xs(n) = xs1
        !
    end subroutine splina
    !*==TRISOL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SPLINA


    subroutine trisol(a, b, c, d, kk)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: kk
        real, dimension(kk) :: a, b, c, d
        !
        ! Local variables
        !
        integer :: k, km
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------
        !     Solves KK long, tri-diagonal system |
        !                                         |
        !             A C          D              |
        !             B A C        D              |
        !               B A .      .              |
        !                 . . C    .              |
        !                   B A    D              |
        !                                         |
        !     The righthand side D is replaced by |
        !     the solution.  A, C are destroyed.  |
        !-----------------------------------------
        !
        do k = 2, kk
            km = k - 1
            c(km) = c(km) / a(km)
            d(km) = d(km) / a(km)
            a(k) = a(k) - b(k) * c(km)
            d(k) = d(k) - b(k) * d(km)
        enddo
        !
        d(kk) = d(kk) / a(kk)
        !
        do k = kk - 1, 1, -1
            d(k) = d(k) - c(k) * d(k + 1)
        enddo
        !
    end subroutine trisol
    !*==GEVAL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! TRISOL


    real function geval(ss, x, xs, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: ss
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        real :: cx1, cx2, dgev, ds, t
        integer :: i, ilow, imid, k
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates int( X(SS) ) dS                   |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        if (n==1) then
            geval = x(1) * (ss - s(1))
            return
        endif
        !
        ilow = 1
        i = n
        !
        do while (i - ilow>1)
            !
            imid = (i + ilow) / 2
            if (ss<s(imid)) then
                i = imid
            else
                ilow = imid
            endif
        enddo
        !
        !
        !---- first integrate up to I-1 point
        geval = 0.
        do k = 2, i - 1
            ds = s(k) - s(k - 1)
            !
            !------ Int X(t) dt  for t = 0..1
            dgev = 0.5 * (x(k) + x(k - 1)) + (xs(k - 1) - xs(k)) * ds / 12.0
            !
            geval = geval + dgev * ds
        enddo
        !
        !---- now integrate up to SS value in I-1..I interval
        ds = s(i) - s(i - 1)
        t = (ss - s(i - 1)) / ds
        cx1 = ds * xs(i - 1) - x(i) + x(i - 1)
        cx2 = ds * xs(i) - x(i) + x(i - 1)
        !
        dgev = 0.5 * t * t * x(i) + (t - 0.5 * t * t) * x(i - 1) + (6.0 - 8.0 * t + 3.0 * t * t)    &
                & * t * t * cx1 / 12.0 + (-4.0 * t + 3.0 * t * t) * t * t * cx2 / 12.0
        !
        geval = geval + dgev * ds
        !
    end function geval
    !*==SEVAL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GEVAL


    real function seval(ss, x, xs, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: ss
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        real :: cx1, cx2, ds, t
        integer :: i, ilow, imid
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates X(SS)                             |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        if (n==1) then
            seval = x(1) + xs(1) * (ss - s(1))
            return
        endif
        !
        ilow = 1
        i = n
        !
        do while (i - ilow>1)
            !
            imid = (i + ilow) / 2
            if (ss<s(imid)) then
                i = imid
            else
                ilow = imid
            endif
        enddo
        !
        ds = s(i) - s(i - 1)
        t = (ss - s(i - 1)) / ds
        cx1 = ds * xs(i - 1) - x(i) + x(i - 1)
        cx2 = ds * xs(i) - x(i) + x(i - 1)
        seval = t * x(i) + (1.0 - t) * x(i - 1) + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)
    end function seval
    !*==DEVAL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEVAL


    real function deval(ss, x, xs, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: ss
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        real :: cx1, cx2, ds, t
        integer :: i, ilow, imid
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates dX/dS(SS)                         |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        if (n==1) then
            deval = xs(1)
            return
        endif
        !
        ilow = 1
        i = n
        !
        do while (i - ilow>1)
            !
            imid = (i + ilow) / 2
            if (ss<s(imid)) then
                i = imid
            else
                ilow = imid
            endif
        enddo
        !
        ds = s(i) - s(i - 1)
        t = (ss - s(i - 1)) / ds
        cx1 = ds * xs(i - 1) - x(i) + x(i - 1)
        cx2 = ds * xs(i) - x(i) + x(i - 1)
        deval = x(i) - x(i - 1) + (1. - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.) * cx2
        deval = deval / ds
    end function deval
    !*==D2VAL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! DEVAL


    real function d2val(ss, x, xs, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: ss
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        real :: cx1, cx2, ds, t
        integer :: i, ilow, imid
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates d2X/dS2(SS)                       |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        if (n==1) then
            d2val = 0.0
            return
        endif
        !
        ilow = 1
        i = n
        !
        do while (i - ilow>1)
            !
            imid = (i + ilow) / 2
            if (ss<s(imid)) then
                i = imid
            else
                ilow = imid
            endif
        enddo
        !
        ds = s(i) - s(i - 1)
        t = (ss - s(i - 1)) / ds
        cx1 = ds * xs(i - 1) - x(i) + x(i - 1)
        cx2 = ds * xs(i) - x(i) + x(i - 1)
        d2val = (6. * t - 4.) * cx1 + (6. * t - 2.0) * cx2
        d2val = d2val / ds**2
    end function d2val
    !*==SEVALL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! D2VAL


    subroutine sevall(ss, x, xs, s, n, xx, xxs, xxss)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: ss, xx, xxs, xxss
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        real :: ds, f0, f1, f2, f3, t
        integer :: i, ilow, imid
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Calculates all spline derivatives.           |
        !     (Combines SEVAL, DEVAL, D2VAL)               |
        !     XS array must have been calculated by SPLINE |
        !--------------------------------------------------
        if (n==1) then
            xx = x(1) + xs(1) * (ss - s(1))
            xxs = xs(1)
            xxss = 0.
            return
        endif
        !
        ilow = 1
        i = n
        !
        do while (i - ilow>1)
            !
            imid = (i + ilow) / 2
            if (ss<s(imid)) then
                i = imid
            else
                ilow = imid
            endif
        enddo
        !
        ds = s(i) - s(i - 1)
        t = (ss - s(i - 1)) / ds
        !
        f0 = x(i - 1)
        f1 = ds * xs(i - 1)
        f2 = -ds * (2.0 * xs(i - 1) + xs(i)) + 3.0 * (x(i) - x(i - 1))
        f3 = ds * (xs(i - 1) + xs(i)) - 2.0 * (x(i) - x(i - 1))
        !
        xx = f0 + t * (f1 + t * (f2 + t * f3))
        xxs = f1 + t * (2.0 * f2 + t * 3.0 * f3)
        xxss = 2.0 * f2 + t * 6.0 * f3
        !
        xxs = xxs / ds
        xxss = xxss / ds**2
        !
    end subroutine sevall
    !*==SEVLIN.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEVALL



    subroutine sevlin(ss, x, s, n, xx, xxs)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: ss, xx, xxs
        real, dimension(n) :: s, x
        !
        ! Local variables
        !
        real :: ds, t
        integer :: i, ilow, imid
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------
        !     Calculates X(SS) and dX/ds(SS) using piecewise-linear  |
        !     interpolation. This is intended for intepolating very  |
        !     noisy data for which a cubic spline is inappropriate.  |
        !------------------------------------------------------------
        if (n==1) then
            xx = x(1)
            xxs = 0.
            return
        endif
        !
        ilow = 1
        i = n
        !
        do while (i - ilow>1)
            !
            imid = (i + ilow) / 2
            if (ss<s(imid)) then
                i = imid
            else
                ilow = imid
            endif
        enddo
        !
        ds = s(i) - s(i - 1)
        t = (ss - s(i - 1)) / ds
        xx = t * x(i) + (1.0 - t) * x(i - 1)
        xxs = (x(i) - x(i - 1)) / ds
        !
    end subroutine sevlin
    !*==CURV.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEVLIN



    real function curv(ss, x, xs, y, ys, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: ss
        real, dimension(n) :: s, x, xs, y, ys
        !
        ! Local variables
        !
        real :: ds, f1, f2, f3, g1, g2, g3, t, xd, xdd, yd, &
                & ydd
        integer :: i, ilow, imid
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Calculates curvature of splined 2-D curve |
        !     at S = SS                                 |
        !                                               |
        !     S        arc length array of curve        |
        !     X, Y     coordinate arrays of curve       |
        !     XS,YS    derivative arrays                |
        !              (calculated earlier by SPLINE)   |
        !-----------------------------------------------
        if (n==1) then
            curv = 0.
            return
        endif
        !
        ilow = 1
        i = n
        !
        do while (i - ilow>1)
            !
            imid = (i + ilow) / 2
            if (ss<s(imid)) then
                i = imid
            else
                ilow = imid
            endif
        enddo
        !
        ds = s(i) - s(i - 1)
        t = (ss - s(i - 1)) / ds
        !
        f1 = ds * xs(i - 1)
        f2 = -ds * (2.0 * xs(i - 1) + xs(i)) + 3.0 * (x(i) - x(i - 1))
        f3 = ds * (xs(i - 1) + xs(i)) - 2.0 * (x(i) - x(i - 1))
        !
        xd = f1 + t * (2.0 * f2 + t * 3.0 * f3)
        xdd = 2.0 * f2 + t * 6.0 * f3
        !
        !
        g1 = ds * ys(i - 1)
        g2 = -ds * (2.0 * ys(i - 1) + ys(i)) + 3.0 * (y(i) - y(i - 1))
        g3 = ds * (ys(i - 1) + ys(i)) - 2.0 * (y(i) - y(i - 1))
        !
        yd = g1 + t * (2.0 * g2 + t * 3.0 * g3)
        ydd = 2.0 * g2 + t * 6.0 * g3
        !
        !
        curv = (xd * ydd - yd * xdd) / sqrt((xd * xd + yd * yd)**3)
        !
    end function curv
    !*==CURVS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CURV


    real function curvs(ss, x, xs, y, ys, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: ss
        real, dimension(n) :: s, x, xs, y, ys
        !
        ! Local variables
        !
        real :: bot, cx1, cx2, cy1, cy2, dbotdt, ds, dtopdt, &
                & f1, f2, f3, g1, g2, g3, sqrtb, t, top, xd, &
                & xdd, xddd, yd, ydd, yddd
        integer :: i, ilow, imid
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Calculates curvature derivative of        |
        !     splined 2-D curve at S = SS               |
        !                                               |
        !     S        arc length array of curve        |
        !     X, Y     coordinate arrays of curve       |
        !     XS,YS    derivative arrays                |
        !              (calculated earlier by SPLINE)   |
        !-----------------------------------------------
        if (n==1) then
            curvs = 0.
            return
        endif
        !
        ilow = 1
        i = n
        !
        do while (i - ilow>1)
            !
            imid = (i + ilow) / 2
            if (ss<s(imid)) then
                i = imid
            else
                ilow = imid
            endif
        enddo
        !
        ds = s(i) - s(i - 1)
        t = (ss - s(i - 1)) / ds
        !
        cx1 = ds * xs(i - 1) - x(i) + x(i - 1)
        cx2 = ds * xs(i) - x(i) + x(i - 1)
        xd = x(i) - x(i - 1) + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2
        xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2
        xddd = 6.0 * cx1 + 6.0 * cx2
        !
        cy1 = ds * ys(i - 1) - y(i) + y(i - 1)
        cy2 = ds * ys(i) - y(i) + y(i - 1)
        yd = y(i) - y(i - 1) + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2
        ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2
        yddd = 6.0 * cy1 + 6.0 * cy2
        !

        f1 = ds * xs(i - 1)
        f2 = -ds * (2.0 * xs(i - 1) + xs(i)) + 3.0 * (x(i) - x(i - 1))
        f3 = ds * (xs(i - 1) + xs(i)) - 2.0 * (x(i) - x(i - 1))
        !
        xd = f1 + t * (2.0 * f2 + t * 3.0 * f3)
        xdd = 2.0 * f2 + t * 6.0 * f3
        xddd = 6.0 * f3
        !
        !
        g1 = ds * ys(i - 1)
        g2 = -ds * (2.0 * ys(i - 1) + ys(i)) + 3.0 * (y(i) - y(i - 1))
        g3 = ds * (ys(i - 1) + ys(i)) - 2.0 * (y(i) - y(i - 1))
        !
        yd = g1 + t * (2.0 * g2 + t * 3.0 * g3)
        ydd = 2.0 * g2 + t * 6.0 * g3
        yddd = 6.0 * g3
        !
        sqrtb = sqrt(xd * xd + yd * yd)
        bot = sqrtb**3
        dbotdt = 3.0 * sqrtb * (xd * xdd + yd * ydd)
        !
        top = xd * ydd - yd * xdd
        dtopdt = xd * yddd - yd * xddd
        !
        curvs = (dtopdt * bot - dbotdt * top) / bot**2 / ds
        !
    end function curvs
    !*==SINVRT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CURVS


    subroutine sinvrt(si, xi, x, xs, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: si, xi
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        real :: ds, xx, xxs, xxss
        integer :: iter
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !     Calculates the "inverse" spline function S(X). |
        !     Since S(X) can be multi-valued or not defined, |
        !      this is not a "black-box" routine.  The call- |
        !      ing program must pass via SI a sufficiently   |
        !      good initial guess for S(XI).                 |
        !                                                    |
        !     XI      specified X value       (input)        |
        !     SI      calculated S(XI) value  (input,output) |
        !     X,XS,S  usual spline arrays     (input)        |
        !                                                    |
        !----------------------------------------------------
        !
        do iter = 1, 10
            call sevall(si, x, xs, s, n, xx, xxs, xxss)
            if (xxs==0.0) exit
            !
            ds = (xi - xx) / xxs
            si = si + ds
            if (abs(ds / (s(n) - s(1)))<1.0e-5) return
        enddo
        write (*, *) 'SINVRT: spline inversion failed.  Continuing...'
        !
    end subroutine sinvrt
    !*==SCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SINVRT


    subroutine scalc(x, y, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(n) :: s, x, y
        !
        ! Local variables
        !
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------
        !     Calculates the arc length array S  |
        !     for a 2-D array of points (X,Y).   |
        !----------------------------------------
        !
        s(1) = 0.
        do i = 2, n
            s(i) = s(i - 1) + sqrt((x(i) - x(i - 1))**2 + (y(i) - y(i - 1))**2)
        enddo
        !
    end subroutine scalc
    !*==SEGSPL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SCALC


    subroutine segspl(x, xs, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        integer :: iseg, iseg0, nseg
        real :: xx
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Splines X(S) array just like SPLINE,      |
        !     but allows derivative discontinuities     |
        !     at segment joints.  Segment joints are    |
        !     defined by identical successive S values. |
        !-----------------------------------------------
        if (n==1) then
            xs(1) = 0.
            return
        endif
        !
        xx = 1.0 / (s(2) - s(1))
        !c        STOP 'SEGSPL:  First input point duplicated'
        if (s(1)==s(2)) xx = 1.0 / (s(2) - s(1))
        !c        STOP 'SEGSPL:  Last  input point duplicated'
        if (s(n)==s(n - 1)) xx = 1.0 / (s(n) - s(n - 1))
        !
        iseg0 = 1
        do iseg = 2, n - 2
            if (s(iseg)==s(iseg + 1)) then
                nseg = iseg - iseg0 + 1
                call spline(x(iseg0), xs(iseg0), s(iseg0), nseg)
                iseg0 = iseg + 1
            endif
        enddo
        !
        nseg = n - iseg0 + 1
        call spline(x(iseg0), xs(iseg0), s(iseg0), nseg)
        !
    end subroutine segspl
    !*==SEGSPD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEGSPL


    subroutine segspd(x, xs, s, n, xs1, xs2)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: xs1, xs2
        real, dimension(n) :: s, x, xs
        !
        ! Local variables
        !
        integer :: iseg, iseg0, nseg
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------
        !     Splines X(S) array just like SPLIND,      |
        !     but allows derivative discontinuities     |
        !     at segment joints.  Segment joints are    |
        !     defined by identical successive S values. |
        !-----------------------------------------------
        if (n==1) then
            xs(1) = 0.
            return
        endif
        !
        if (s(1)==s(2)) stop 'SEGSPD:  First input point duplicated'
        if (s(n)==s(n - 1)) stop 'SEGSPD:  Last  input point duplicated'
        !
        iseg0 = 1
        do iseg = 2, n - 2
            if (s(iseg)==s(iseg + 1)) then
                nseg = iseg - iseg0 + 1
                call splind(x(iseg0), xs(iseg0), s(iseg0), nseg, xs1, xs2)
                iseg0 = iseg + 1
            endif
        enddo
        !
        nseg = n - iseg0 + 1
        call splind(x(iseg0), xs(iseg0), s(iseg0), nseg, xs1, xs2)
        !
    end subroutine segspd
    !*==INTERS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SEGSPD



    subroutine inters(ok, ss1, ss2, x1, xs1, y1, ys1, s1, n1, x2, xs2, y2, ys2, s2, &
            & n2)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n1, n2
        logical :: ok
        real :: ss1, ss2
        real, dimension(n1) :: s1, x1, xs1, y1, ys1
        real, dimension(n2) :: s2, x2, xs2, y2, ys2
        !
        ! Local variables
        !
        real :: a11, a12, a21, a22, det, ds1, ds2, rlx, rs1, &
                & rs1old, rs2, rs2old, ss1old, ss2old, xx1, xx2, &
                & yy1, yy2
        logical :: clip1, clip2
        real, save :: eps
        integer :: irlx, iter
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Finds spline coordinate values SS1, SS2 at the
        !     intersection of two space curves (X1,Y1), (X2,Y2).
        !-------------------------------------------------------
        data eps/1.0e-5/
        !
        ok = .true.
        !cc      SS1 = S1(1)
        !cc      SS2 = S2(1)
        rs1 = 1.0e12
        rs2 = 1.0e12
        ds1 = 0.0
        ds2 = 0.0
        !
        do iter = 1, 12
            !
            rlx = 1.0
            ss1old = ss1
            ss2old = ss2
            rs1old = abs(rs1)
            rs2old = abs(rs2)
            !
            do irlx = 1, 16
                !
                clip1 = .false.
                clip2 = .false.
                ss1 = ss1old + rlx * ds1
                ss2 = ss2old + rlx * ds2
                !
                if (ss1<s1(1) .or. ss1>s1(n1)) then
                    clip1 = .true.
                    ss1 = max(ss1, s1(1))
                    ss1 = min(ss1, s1(n1))
                endif
                if (ss2<s2(1) .or. ss2>s2(n2)) then
                    clip2 = .true.
                    ss2 = max(ss2, s2(1))
                    ss2 = min(ss2, s2(n2))
                endif
                !
                xx1 = seval(ss1, x1, xs1, s1, n1)
                xx2 = seval(ss2, x2, xs2, s2, n2)
                yy1 = seval(ss1, y1, ys1, s1, n1)
                yy2 = seval(ss2, y2, ys2, s2, n2)
                !
                rs1 = xx1 - xx2
                rs2 = yy1 - yy2
                !
                if (abs(rs1)<rs1old .and. abs(rs2)<rs2old) goto 11
                !
                rlx = 0.5 * rlx
                !
            enddo
            write (*, *) 'INTERS: Under-relaxation loop failed.'
            !
            11      a11 = deval(ss1, x1, xs1, s1, n1)
            a12 = -deval(ss2, x2, xs2, s2, n2)
            a21 = deval(ss1, y1, ys1, s1, n1)
            a22 = -deval(ss2, y2, ys2, s2, n2)
            !
            det = a11 * a22 - a12 * a21
            ds1 = -(rs1 * a22 - a12 * rs2) / det
            ds2 = -(a11 * rs2 - rs1 * a21) / det
            !
            if (abs(ds1)<eps * (s1(n1) - s1(1)) .and. abs(ds2)                &
                    & <eps * (s2(n2) - s2(1))) return
            !
        enddo
        write (*, *) 'INTERS: Convergence failed. Res =', rs1, rs2
        if (clip1) write (*, *) '        S1 clip:', s1(1), s1(n1), &
                & ss1, ds1
        if (clip2) write (*, *) '        S2 clip:', s2(1), s2(n2), &
                & ss2, ds2
        ok = .false.
        !
    end subroutine inters
end module m_spline
