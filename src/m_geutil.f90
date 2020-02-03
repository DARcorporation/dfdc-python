module m_geutil
    implicit none
contains
    !*==MINMAX.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ITPSET
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

    subroutine minmax(n, x, xmin, xmax)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: xmax, xmin
        real, dimension(*) :: x
        !
        ! Local variables
        !
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------
        !     Calculates min,max of array X
        !--------------------------------------
        !
        i = 1
        xmin = x(i)
        xmax = x(i)
        do i = 2, n
            xmin = min(xmin, x(i))
            xmax = max(xmax, x(i))
        enddo
        !
    end subroutine minmax
    !*==NCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine ncalc(x, y, s, n, xn, yn)
        use m_spline, only : segspl
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(*) :: s, x, xn, y, yn
        !
        ! Local variables
        !
        integer :: i
        real :: smod, sx, sy
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------
        !     Calculates normal unit vector
        !     components at airfoil panel nodes
        !---------------------------------------
        !
        if (n<=1) return
        !
        call segspl(x, xn, s, n)
        call segspl(y, yn, s, n)
        do i = 1, n
            sx = yn(i)
            sy = -xn(i)
            smod = sqrt(sx * sx + sy * sy)
            xn(i) = sx / smod
            yn(i) = sy / smod
        enddo
        !
        !---- average normal vectors at corner points
        do i = 1, n - 1
            if (s(i)==s(i + 1)) then
                sx = 0.5 * (xn(i) + xn(i + 1))
                sy = 0.5 * (yn(i) + yn(i + 1))
                smod = sqrt(sx * sx + sy * sy)
                xn(i) = sx / smod
                yn(i) = sy / smod
                xn(i + 1) = sx / smod
                yn(i + 1) = sy / smod
            endif
        enddo
        !
    end subroutine ncalc
    !*==LEFIND.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine lefind(sle, x, xp, y, yp, s, n)
        use m_spline, only : d2val, seval, deval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: sle
        real, dimension(*) :: s, x, xp, y, yp
        !
        ! Local variables
        !
        real :: dotp, dseps, dsle, dx, dxdd, dxds, dxte, dy, &
                & dydd, dyds, dyte, res, ress, xchord, xle, xte, &
                & ychord, yle, yte
        integer :: i, iter
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Locates leading edge spline-parameter value SLE
        !
        !     The defining condition is
        !
        !      (X-XTE,Y-YTE) . (X',Y') = 0     at  S = SLE
        !
        !     i.e. the surface tangent is normal to the chord
        !     line connecting X(SLE),Y(SLE) and the TE point.
        !------------------------------------------------------
        if (n<=1) then
            sle = s(1)
            return
        endif
        !
        !---- convergence tolerance
        dseps = (s(n) - s(1)) * 1.0e-5
        !
        !---- set trailing edge point coordinates
        xte = 0.5 * (x(1) + x(n))
        yte = 0.5 * (y(1) + y(n))
        !
        !---- get first guess for SLE
        do i = 3, n - 2
            dxte = x(i) - xte
            dyte = y(i) - yte
            dx = x(i + 1) - x(i)
            dy = y(i + 1) - y(i)
            dotp = dxte * dx + dyte * dy
            if (dotp<0.0) exit
        enddo
        !
        sle = s(i)
        !
        !---- check for sharp LE case
        if (s(i)==s(i - 1)) return
        !
        !---- Newton iteration to get exact SLE value
        do iter = 1, 50
            xle = seval(sle, x, xp, s, n)
            yle = seval(sle, y, yp, s, n)
            dxds = deval(sle, x, xp, s, n)
            dyds = deval(sle, y, yp, s, n)
            dxdd = d2val(sle, x, xp, s, n)
            dydd = d2val(sle, y, yp, s, n)
            !
            xchord = xle - xte
            ychord = yle - yte
            !
            !------ drive dot product between chord line and LE tangent to zero
            res = xchord * dxds + ychord * dyds
            ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd
            !
            if (ress==0.0) exit
            !
            !------ Newton delta for SLE
            dsle = -res / ress
            !
            dsle = max(dsle, -0.02 * abs(xchord + ychord))
            dsle = min(dsle, 0.02 * abs(xchord + ychord))
            sle = sle + dsle
            if (abs(dsle)<dseps) return
        enddo
        write (*, *) 'LEFIND:  LE point not found.  Continuing...'
        !
        sle = s(i)
    end subroutine lefind
    !*==CNORM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine cnorm(x, xp, y, yp, s, n, xtran, ytran, scale, angle)
        use m_spline, only : seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: angle, scale, xtran, ytran
        integer :: n
        real, dimension(*) :: s, x, xp, y, yp
        !
        ! Local variables
        !
        real :: cinv, csqinv, cssq, cx, cy, sle, xbar, xle, &
                & xte, ybar, yle, yte
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------
        !     Scales and rotates coordinates to get unit chordline
        !----------------------------------------------------------
        call lefind(sle, x, xp, y, yp, s, n)
        !
        xte = 0.5 * (x(1) + x(n))
        yte = 0.5 * (y(1) + y(n))
        xle = seval(sle, x, xp, s, n)
        yle = seval(sle, y, yp, s, n)
        !
        cx = xte - xle
        cy = yte - yle
        cssq = cx * cx + cy * cy
        !
        if (cssq==0.0) then
            csqinv = 1.0
            cinv = 1.0
        else
            csqinv = 1.0 / cssq
            cinv = 1.0 / sqrt(cssq)
        endif
        !
        do i = 1, n
            xbar = x(i) - xle
            ybar = y(i) - yle
            x(i) = (xbar * cx + ybar * cy) * csqinv
            y(i) = (-xbar * cy + ybar * cx) * csqinv
            s(i) = s(i) * cinv
        enddo
        !
        xtran = -xle
        ytran = -yle
        scale = cinv
        angle = atan2(cy, cx)
        !
    end subroutine cnorm
    !*==GEOPAR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CNORM


    subroutine geopar(x, xp, y, yp, s, n, t, sle, chord, radle, angte, area, &
            & eixxa, eiyya, eixya, askn, eixxt, eiyyt, eixyt, thick, &
            & cambr, volm, vskn, asrf, rgxv, rgyv)
        use m_spline, only : seval, curv
        use m_xutils, only : atanc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: angte, area, askn, asrf, cambr, chord, eixxa, &
                & eixxt, eixya, eixyt, eiyya, eiyyt, radle, rgxv, &
                & rgyv, sle, thick, volm, vskn
        integer :: n
        real, dimension(*) :: s, t, x, xp, y, yp
        !
        ! Local variables
        !
        real :: ang1, ang2, cfac, chsq, curvle, perim, rgxt, &
                & rgyt, tfac, xcena, xcent, xcenv, xcenvt, xle, &
                & xnew, xte, ycena, ycent, ycenv, ycenvt, yle, &
                & ynew, yte
        real, dimension(1) :: xnew_temp, ynew_temp
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Sets geometric parameters for airfoil shape
        !------------------------------------------------------
        !
        call lefind(sle, x, xp, y, yp, s, n)
        !
        xle = seval(sle, x, xp, s, n)
        yle = seval(sle, y, yp, s, n)
        xte = 0.5 * (x(1) + x(n))
        yte = 0.5 * (y(1) + y(n))
        !
        chsq = (xte - xle)**2 + (yte - yle)**2
        chord = sqrt(chsq)
        !
        curvle = curv(sle, x, xp, y, yp, s, n)
        !
        if (abs(curvle)>0.001 * (s(n) - s(1))) then
            radle = 1.0 / curvle
        else
            radle = 0.
        endif
        !
        ang1 = atan2(-yp(1), -xp(1))
        ang2 = atanc(yp(n), xp(n), ang1)
        angte = ang2 - ang1
        !
        !---- 2D properties
        call aecalc(n, x, y, t, 0, perim, area, xcena, ycena, eixxa, eiyya, eixya)
        call aecalc(n, x, y, t, 2, perim, askn, xcent, ycent, eixxt, eiyyt, eixyt)
        tfac = 1.0
        cfac = 1.0
        call tcset(x, xp, y, yp, s, n, thick, cambr, tfac, cfac, xnew_temp, ynew_temp, .false.)
        xnew = xnew_temp(1)
        ynew = ynew_temp(1)
        !
        !---- Axisymmetric properties
        call axcalc(n, x, y, t, 0, volm, asrf, xcenv, ycenv, rgxv, rgyv)
        call axcalc(n, x, y, t, 2, vskn, asrf, xcenvt, ycenvt, rgxt, rgyt)
        !
    end subroutine geopar
    !*==AXCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine axcalc(n, x, y, t, itype, volm, area, xcen, ycen, rgx, rgy)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: area, rgx, rgy, volm, xcen, ycen
        integer :: itype, n
        real, dimension(*) :: t, x, y
        !
        ! Local variables
        !
        real :: aint, da, ds, dste, dv, dx, dy, rgxsq, sint, &
                & ta, vint, xa, xint, xxint, ya, yint, yyint
        integer :: io, ip
        logical :: laxibod, lopen
        real, save :: pi
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------
        !     Calculates geometric properties of axisymmetric shape X,Y
        !     where X is axial coordinate and Y is radial coordinate
        !
        !     Input:
        !       N      number of points
        !       X(.)   shape coordinate point arrays
        !       Y(.)
        !       T(.)   skin-thickness array, used only if ITYPE = 2
        !       ITYPE  = 0 ...   integration is over volume
        !              = 2 ...   integration is over skin with thickness specified
        !
        !     Output:
        !       VOLM  volume
        !       AREA  surface area
        !       XCEN  centroid location
        !       YCEN
        !       RGX   radius of gyration about X centroid
        !       RGY   radius of gyration about Y centroid (Y=0)
        !---------------------------------------------------------------
        data pi/3.141592653589793238/
        !
        !--- check for shape symmetric about X axis
        !    if first or last point is on centerline it is an axisymmetric body
        laxibod = (y(1)==0.0 .or. y(n)==0.0)
        !
        !--- check shape for "open" or closed body by comparing distance between
        !    first and last points to 0.1 of perimeter length
        !
        sint = 0.0
        do io = 2, n
            dx = x(io) - x(io - 1)
            dy = y(io) - y(io - 1)
            ds = sqrt(dx * dx + dy * dy)
            sint = sint + ds
        enddo
        dx = x(n) - x(1)
        dy = y(n) - y(1)
        dste = sqrt(dx * dx + dy * dy)
        lopen = (dste>0.1 * sint)
        !
        !--- integrate around contour to get geometric quantities
        vint = 0.0
        aint = 0.0
        xint = 0.0
        yint = 0.0
        xxint = 0.0
        yyint = 0.0
        !
        do io = 1, n
            if (io<n) then
                ip = io + 1
            else
                !--- if last point is on centerline or body is open end here...
                if (laxibod .or. lopen) cycle
                !--- otherwise close body with last segment between first and last points
                ip = 1
            endif
            !
            !        write(54,199) io,x(io),y(io),x(ip),y(ip)
            dx = x(io) - x(ip)
            dy = y(io) - y(ip)
            xa = (x(io) + x(ip)) * 0.5
            ya = (y(io) + y(ip)) * 0.5
            ds = sqrt(dx * dx + dy * dy)
            199     format (i4, 7f10.6)
            !---- surface area
            da = 2.0 * pi * ya * ds
            aint = aint + da

            if (itype==0) then
                !-------- treat as a axisymmetric solid body (YCEN=0,intXY=0)
                dv = pi * dx * ya**2
                !        write(55,199) io,dx,dy,xa,ya,ds,da,dv
                vint = vint + dv
                xint = xint + xa * dv
                yint = 0.0
                xxint = xxint + xa * xa * dv
                yyint = yyint + ya * ya * dv / 4.0
            else
                !--------- integrate over skin thickness
                ta = (t(io) + t(ip)) * 0.50
                dv = 2.0 * pi * ya * ta * ds
                vint = vint + dv
                xint = xint + xa * dv
                yint = 0.0
                xxint = xxint + xa * xa * dv
                yyint = yyint + ya * ya * dv
            endif
        enddo
        !
        !---- volume and surface area
        if (lopen .and. .not.laxibod) vint = 0.0
        volm = vint
        area = aint
        !
        !---- calculate centroid location
        if (vint<=0.0) then
            xcen = 0.
            ycen = 0.
            rgx = 0.
            rgy = 0.
        else
            xcen = xint / vint
            ycen = yint / vint
            !---- calculate radii of gyration
            rgxsq = xxint / vint - xcen**2
            if (rgxsq>=0.0) rgx = sqrt(rgxsq)
            rgy = sqrt(yyint / vint)
        endif
        !
    end subroutine axcalc
    !*==AECALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! AXCALC



    subroutine aecalc(n, x, y, t, itype, perim, area, xcen, ycen, eixx, eiyy, &
            & eixy)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: area, eixx, eixy, eiyy, perim, xcen, ycen
        integer :: itype, n
        real, dimension(*) :: t, x, y
        !
        ! Local variables
        !
        real :: aint, da, ds, dx, dy, sint, ta, vint, xa, &
                & xint, xxint, xyint, ya, yint, yyint
        integer :: io, ip
        logical :: laxis
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------
        !     Calculates geometric properties of 2D shape X,Y
        !
        !     Input:
        !       N      number of points
        !       X(.)   shape coordinate point arrays
        !       Y(.)
        !       T(.)   skin-thickness array, used only if ITYPE = 2
        !       ITYPE  = 0 ...   integration is over whole area  dA = dx dy
        !              = 1 ...   integration is over perimeter   dA =   ds
        !              = 2 ...   integration is over skin area   dA = t ds
        !
        !     Output:
        !       PERIM perimeter length
        !       AREA  area (assumes 2D)
        !       XCEN  centroid location
        !       YCEN
        !       EIXX  moments of inertia
        !       EIYY
        !       EIXY
        !---------------------------------------------------------------
        !
        !--- check for shape with symmetry about X axis
        !    if first or last point is on centerline it is an axisymmetric body
        laxis = (y(1)==0.0 .or. y(n)==0.0)
        !
        sint = 0.0
        vint = 0.0
        aint = 0.0
        xint = 0.0
        yint = 0.0
        xxint = 0.0
        xyint = 0.0
        yyint = 0.0
        !
        do io = 1, n
            !
            if (io<n) then
                ip = io + 1
            else
                !--- if last point is on axis (symmetrical body)
                if (laxis) cycle
                !--- otherwise close body with last segment between first and last points
                ip = 1
            endif
            !
            dx = x(io) - x(ip)
            dy = y(io) - y(ip)
            xa = (x(io) + x(ip)) * 0.5
            ya = (y(io) + y(ip)) * 0.5
            !
            ds = sqrt(dx * dx + dy * dy)
            sint = sint + ds

            if (itype==0) then
                !-------- integrate over airfoil cross-section
                da = ya * dx
                aint = aint + da
                xint = xint + xa * da
                yint = yint + ya * da / 2.0
                xxint = xxint + xa * xa * da
                xyint = xyint + xa * ya * da / 2.0
                yyint = yyint + ya * ya * da / 3.0
            else
                if (itype==1) then
                    !--------- integrate over perimeter
                    da = ds
                else
                    !--------- integrate over skin thickness
                    ta = (t(io) + t(ip)) * 0.50
                    da = ta * ds
                endif
                aint = aint + da
                xint = xint + xa * da
                yint = yint + ya * da
                xxint = xxint + xa * xa * da
                xyint = xyint + xa * ya * da
                yyint = yyint + ya * ya * da
            endif
        enddo
        !
        perim = sint
        area = aint
        !
        !---- calculate centroid location
        if (aint==0.0) then
            xcen = 0.
            ycen = 0.
        else
            xcen = xint / aint
            ycen = yint / aint
        endif
        !
        !---- calculate inertias
        eixx = yyint - ycen * ycen * aint
        eiyy = xxint - xcen * xcen * aint
        eixy = xyint - xcen * ycen * aint
        !
    end subroutine aecalc
    !*==PRAXIS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! AECALC



    subroutine praxis(fxx, fyy, fxy, f11, f22, ap1, ap2)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: ap1, ap2, f11, f22, fxx, fxy, fyy
        !
        ! Local variables
        !
        real :: c1, c2, fsq, s1, s2, sgn
        real, save :: pi
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------
        !     Sets prinicipal components F11,F22 of tensor FXX,FYY,FXY.
        !     Also sets corresponding principal direction angls AP1,AP2.
        !-----------------------------------------------------------------
        !
        data pi/3.141592653589793238/
        !
        fsq = 0.25 * (fxx - fyy)**2 + fxy**2
        sgn = sign(1.0, fyy - fxx)
        f11 = 0.5 * (fxx + fyy) - sgn * sqrt(fsq)
        f22 = 0.5 * (fxx + fyy) + sgn * sqrt(fsq)
        !
        if (f11==0.0 .or. f22==0.0) then
            !----- vanishing tensor limit
            ap1 = 0.0
            ap2 = atan2(1.0, 0.0)
            !
        elseif (fsq / (f11 * f22)<0.0001) then
            !----- rotationally-invariant tensor (circle, square, etc.)
            ap1 = 0.0
            ap2 = atan2(1.0, 0.0)
            !
        else
            !----- normal general tensor case
            c1 = fxy
            s1 = fxx - f11
            !
            c2 = fxy
            s2 = fxx - f22
            !
            if (abs(s1)>abs(s2)) then
                ap1 = atan2(s1, c1)
                ap2 = ap1 + 0.5 * pi
            else
                ap2 = atan2(s2, c2)
                ap1 = ap2 - 0.5 * pi
            endif

            if (ap1<-0.5 * pi) ap1 = ap1 + pi
            if (ap1>+0.5 * pi) ap1 = ap1 - pi
            if (ap2<-0.5 * pi) ap2 = ap2 + pi
            if (ap2>+0.5 * pi) ap2 = ap2 - pi
            !
        endif

    end subroutine praxis
    !*==TEGAP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine tegap(x, xs, y, ys, s, n, doc, gapnew)
        use m_spline, only : scalc, seval, segspl
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: doc, gapnew
        integer :: n
        real, dimension(*) :: s, x, xs, y, ys
        !
        ! Local variables
        !
        real :: arg, chsq, dgap, dxn, dxu, dyn, dyu, gap, sle, &
                & tfac, xle, xoc, xte, yle, yte
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------
        !     Used to set buffer airfoil
        !     trailing edge gap
        !----------------------------------
        !
        call lefind(sle, x, xs, y, ys, s, n)
        xle = seval(sle, x, xs, s, n)
        yle = seval(sle, y, ys, s, n)
        xte = 0.5 * (x(1) + x(n))
        yte = 0.5 * (y(1) + y(n))
        chsq = (xte - xle)**2 + (yte - yle)**2
        !
        dxn = x(1) - x(n)
        dyn = y(1) - y(n)
        gap = sqrt(dxn**2 + dyn**2)
        !
        !---- components of unit vector parallel to TE gap
        if (gap>0.0) then
            dxu = dxn / gap
            dyu = dyn / gap
        else
            dxu = -.5 * (ys(n) - ys(1))
            dyu = 0.5 * (xs(n) - xs(1))
        endif
        !
        dgap = gapnew - gap
        !
        !---- go over each point, changing the y-thickness appropriately
        do i = 1, n
            !------ chord-based x/c
            xoc = ((x(i) - xle) * (xte - xle) + (y(i) - yle) * (yte - yle)) / chsq
            !
            !------ thickness factor tails off exponentially away from trailing edge
            if (doc==0.0) then
                tfac = 0.0
                if (i==1 .or. i==n) tfac = 1.0
            else
                arg = min((1.0 - xoc) * (1.0 / doc - 1.0), 15.0)
                tfac = exp(-arg)
            endif
            !
            if (s(i)<=sle) then
                x(i) = x(i) + 0.5 * dgap * xoc * tfac * dxu
                y(i) = y(i) + 0.5 * dgap * xoc * tfac * dyu
            else
                x(i) = x(i) - 0.5 * dgap * xoc * tfac * dxu
                y(i) = y(i) - 0.5 * dgap * xoc * tfac * dyu
            endif
        enddo
        !
        call scalc(x, y, s, n)
        call segspl(x, xs, s, n)
        call segspl(y, ys, s, n)
        !
    end subroutine tegap
    !*==TCSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! TEGAP



    subroutine tcset(x, xp, y, yp, s, n, tmax, cmax, tfac, cfac, xnew, ynew, &
            & lnewset)
        use m_spline, only : deval, seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: cfac, cmax, tfac, tmax
        logical :: lnewset
        integer :: n
        real, dimension(*) :: s, x, xnew, xp, y, ynew, yp
        !
        ! Local variables
        !
        real :: chbsq, cmin, dsopp, dxc, dyc, res, resd, sle, &
                & sopp, stot, xbar, xle, xopp, xoppd, xte, ybar, &
                & ybarct, ybarop, yle, yopp, yoppd, yte
        real, save :: eps
        integer :: i, itopp
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Determines max thickness and camber.
        !     Scales thickness and camber by TFAC,CFAC
        !-----------------------------------------------------
        !
        data eps/1.0e-5/
        !
        call lefind(sle, x, xp, y, yp, s, n)
        xle = seval(sle, x, xp, s, n)
        yle = seval(sle, y, yp, s, n)
        xte = 0.5 * (x(1) + x(n))
        yte = 0.5 * (y(1) + y(n))
        chbsq = (xte - xle)**2 + (yte - yle)**2
        !
        if (chbsq==0.0) then
            tmax = 0.
            cmax = 0.
            return
        endif
        !
        !---- set unit chord-line vector
        dxc = (xte - xle) / sqrt(chbsq)
        dyc = (yte - yle) / sqrt(chbsq)
        !
        stot = s(n) - s(1)
        !
        tmax = 0.
        cmax = 0.
        cmin = 0.
        !
        !---- go over each point, finding and/or changing the y-thickness
        !-    (defined normal to chord line)
        do i = 1, n
            !
            !------ coordinates in chord-line axes
            xbar = (x(i) - xle) * dxc + (y(i) - yle) * dyc
            ybar = (y(i) - yle) * dxc - (x(i) - xle) * dyc
            !
            sopp = 2.0 * sle - s(i)
            sopp = max(sopp, s(1))
            sopp = min(sopp, s(n))
            !
            if (abs(sopp - sle)<=eps * stot) then
                sopp = sle
            else
                !------- converge on exact opposite point with same chord x value
                do itopp = 1, 12
                    xopp = seval(sopp, x, xp, s, n)
                    yopp = seval(sopp, y, yp, s, n)
                    xoppd = deval(sopp, x, xp, s, n)
                    yoppd = deval(sopp, y, yp, s, n)
                    !
                    res = (xopp - xle) * dxc + (yopp - yle) * dyc - xbar
                    resd = xoppd * dxc + yoppd * dyc
                    !
                    if (abs(res) / (s(n) - s(1))<eps) goto 305
                    if (resd==0.0) exit
                    !
                    dsopp = -res / resd
                    sopp = sopp + dsopp
                    !
                    if (abs(dsopp / (s(n) - s(1)))<eps) goto 305
                enddo
                !
                write (*, *)                                                 &
                        &'TCSET: Opposite-point location failed. Continuing...'
                sopp = 2.0 * sle - s(i)
                !
            endif
            !
            !------ set point on the opposite side with the same chord x value
            305     xopp = seval(sopp, x, xp, s, n)
            yopp = seval(sopp, y, yp, s, n)
            !
            ybarop = (yopp - yle) * dxc - (xopp - xle) * dyc
            !
            if (lnewset) then
                !------- set new chord x,y coordinates by changing camber & thickness
                ybarct = cfac * 0.5 * (ybar + ybarop) + tfac * 0.5 * (ybar - ybarop)
                xnew(i) = xle + xbar * dxc - ybarct * dyc
                ynew(i) = yle + ybarct * dxc + xbar * dyc
            endif
            !
            tmax = max(tmax, abs(ybar - ybarop))
            cmax = max(cmax, 0.5 * (ybar + ybarop))
            cmin = min(cmin, 0.5 * (ybar + ybarop))
        enddo
        !
        if (-cmin>cmax) cmax = cmin
        !
    end subroutine tcset
    !*==YSYM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine ysym(x, xp, y, yp, s, nx, n, iside, xnew, ynew, nnew)
        use m_spline, only : seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: iside, n, nnew, nx
        real, dimension(nx) :: s, x, xnew, xp, y, ynew, yp
        !
        ! Local variables
        !
        real :: chsq, ds, dxc, dyc, sle, xbar, xle, xte, ybar, &
                & yle, yte
        integer :: i, i1, i2, ig1, ig2, igdir, ile, ile1, ile2
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Makes passed-in airfoil symmetric about chord line.
        !---------------------------------------------------------
        !
        call lefind(sle, x, xp, y, yp, s, n)
        xle = seval(sle, x, xp, s, n)
        yle = seval(sle, y, yp, s, n)
        xte = 0.5 * (x(1) + x(n))
        yte = 0.5 * (y(1) + y(n))
        chsq = (xte - xle)**2 + (yte - yle)**2
        !
        if (chsq==0.0) return
        !
        !---- set unit chord-line vector
        dxc = (xte - xle) / sqrt(chsq)
        dyc = (yte - yle) / sqrt(chsq)
        !
        !---- find index of node ILE which is just before leading edge point
        do i = 1, n
            if (s(i)>=sle) exit
        enddo
        ile = i - 1
        !
        ds = s(ile + 1) - s(ile)
        if (sle - s(ile - 1)<0.1 * ds) then
            !------ point is just before LE, we will move it ahead to LE
            ile1 = ile - 1
            ile2 = ile + 1
        elseif (s(ile + 1) - sle<0.1 * ds) then
            !------ point is just after LE, we will move it back to LE
            ile1 = ile
            ile2 = ile + 2
        else
            !------ no point is near LE ... we will add new point
            ile1 = ile
            ile2 = ile + 1
        endif
        !
        !---- set index limits of side which will set symmetric geometry
        if (iside==1) then
            ig1 = 1
            ig2 = ile1
            igdir = +1
        else
            ig1 = n
            ig2 = ile2
            igdir = -1
        endif
        !
        !---- set new number of points, including LE point
        nnew = 2 * (iabs(ig2 - ig1) + 1) + 1
        if (nnew>nx) stop 'YSYM:  Array overflow on passed arrays.'
        !
        !---- set symmetric geometry
        do i = ig1, ig2, igdir
            !
            !------ coordinates in chord-line axes
            xbar = (x(i) - xle) * dxc + (y(i) - yle) * dyc
            ybar = (y(i) - yle) * dxc - (x(i) - xle) * dyc
            !
            i1 = 1 + (i - ig1) * igdir
            i2 = nnew - (i - ig1) * igdir
            !
            xnew(i1) = xle + xbar * dxc - ybar * dyc
            xnew(i2) = xle + xbar * dxc + ybar * dyc
            !
            ynew(i1) = yle + ybar * dxc + xbar * dyc
            ynew(i2) = yle - ybar * dxc + xbar * dyc
        enddo
        !
        !---- set new LE point
        xnew(nnew / 2 + 1) = xle
        ynew(nnew / 2 + 1) = yle
        !
    end subroutine ysym
    !*==LERSCL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! YSYM


    subroutine lerscl(x, xp, y, yp, s, n, doc, rfac, xnew, ynew)
        use m_spline, only : deval, seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: doc, rfac
        integer :: n
        real, dimension(*) :: s, x, xnew, xp, y, ynew, yp
        !
        ! Local variables
        !
        real :: arg, chord, dsopp, dxc, dyc, res, resd, sle, &
                & sopp, srfac, tfac, xbar, xle, xoc, xopp, xoppd, &
                & xte, ybar, ybarct, ybarop, yle, yopp, yoppd, yte
        integer :: i, itopp
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Adjusts airfoil to scale LE radius by factor RFAC.
        !     Blending of new shape is done with decay length DOC.
        !---------------------------------------------------------
        !
        call lefind(sle, x, xp, y, yp, s, n)
        xle = seval(sle, x, xp, s, n)
        yle = seval(sle, y, yp, s, n)
        xte = 0.5 * (x(1) + x(n))
        yte = 0.5 * (y(1) + y(n))
        chord = sqrt((xte - xle)**2 + (yte - yle)**2)
        !
        !---- set unit chord-line vector
        dxc = (xte - xle) / chord
        dyc = (yte - yle) / chord
        !
        srfac = sqrt(abs(rfac))
        !
        !---- go over each point, changing the y-thickness appropriately
        do i = 1, n
            !
            !------ coordinates in chord-line axes
            xbar = (x(i) - xle) * dxc + (y(i) - yle) * dyc
            ybar = (y(i) - yle) * dxc - (x(i) - xle) * dyc
            !
            sopp = 2.0 * sle - s(i)
            sopp = max(sopp, s(1))
            sopp = min(sopp, s(n))
            !
            if (abs(sopp / sle - 1.0)<=1.0e-5) then
                sopp = sle
            else
                !------- converge on exact opposite point with same chord x value
                do itopp = 1, 12
                    xopp = seval(sopp, x, xp, s, n)
                    yopp = seval(sopp, y, yp, s, n)
                    xoppd = deval(sopp, x, xp, s, n)
                    yoppd = deval(sopp, y, yp, s, n)
                    !
                    res = (xopp - xle) * dxc + (yopp - yle) * dyc - xbar
                    resd = xoppd * dxc + yoppd * dyc
                    !
                    if (abs(res) / (s(n) - s(1))<1.0e-5) goto 305
                    if (resd==0.0) exit
                    !
                    dsopp = -res / resd
                    sopp = sopp + dsopp
                    !
                    if (abs(dsopp / (s(n) - s(1)))<1.0e-5) goto 305
                enddo
                write (*, *)                                                 &
                        &'LERSCL: Opposite-point location failed. Continuing...'
                sopp = 2.0 * sle - s(i)
            endif
            !
            !------ set point on the opposite side with the same chord x value
            305     xopp = seval(sopp, x, xp, s, n)
            yopp = seval(sopp, y, yp, s, n)
            !
            ybarop = (yopp - yle) * dxc - (xopp - xle) * dyc
            !
            !------ thickness factor tails off exponentially towards trailing edge
            xoc = xbar / chord
            arg = min(xoc / doc, 15.0)
            tfac = 1.0 - (1.0 - srfac) * exp(-arg)
            !
            !------ set new chord x,y coordinates by changing thickness locally
            ybarct = 0.5 * (ybar + ybarop) + tfac * 0.5 * (ybar - ybarop)
            !
            xnew(i) = xle + xbar * dxc - ybarct * dyc
            ynew(i) = yle + ybarct * dxc + xbar * dyc
        enddo
        !
    end subroutine lerscl
    !*==FLAP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine flap(x, xs, y, ys, s, n, xf, yf, adef, insid, xnew, ynew, nnew, &
            & tops, atop, xtop, ytop, bots, abot, xbot, ybot)
        use m_spline, only : deval, seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: abot, adef, atop, bots, tops, xbot, xf, xtop, &
                & ybot, yf, ytop
        logical :: insid
        integer :: n, nnew
        real, dimension(*) :: s, x, xnew, xs, y, ynew, ys
        !
        ! Local variables
        !
        real :: ang, ca, chx, chy, cosd, crsp, dang, dsavg, &
                & dsnew, fvx, fvy, sa, sb1, sb1p, sb1q, sb2, &
                & sb2p, sb2q, sfrac, sind, st1, st1p, st1q, st2, &
                & st2p, st2q, xb1, xb1new, xb2, xb2new, xbar, xt1, &
                & xt1new, xt2, xt2new, yb1, yb1new, yb2, yb2new, &
                & ybar, yt1, yt1new, yt2, yt2new
        integer :: i, ib1, ib2, ib2q, inew, ip, it1, it2, it2q, &
                & npadd
        logical :: lb1new, lb2new, lt1new, lt2new
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !     Sets modified airfoil with a deflected flap.
        !     Points may be added/subtracted in the flap
        !     break vicinity to clean things up.
        !     Note: angle deflection ADEF is in radians
        !----------------------------------------------------
        !
        !
        !
        if (insid) then
            atop = max(0.0, -adef)
            abot = max(0.0, adef)
        else
            chx = deval(bots, x, xs, s, n) - deval(tops, x, xs, s, n)
            chy = deval(bots, y, ys, s, n) - deval(tops, y, ys, s, n)
            fvx = seval(bots, x, xs, s, n) + seval(tops, x, xs, s, n)
            fvy = seval(bots, y, ys, s, n) + seval(tops, y, ys, s, n)
            crsp = chx * (yf - 0.5 * fvy) - chy * (xf - 0.5 * fvx)
            if (crsp>0.0) then
                !-------- flap hinge is above airfoil
                atop = max(0.0, adef)
                abot = max(0.0, adef)
            else
                !-------- flap hinge is below airfoil
                atop = max(0.0, -adef)
                abot = max(0.0, -adef)
            endif
        endif
        !
        !---- find upper and lower surface break arc length values...
        call sss(tops, st1, st2, atop, xf, yf, x, xs, y, ys, s, n, 1)
        call sss(bots, sb1, sb2, abot, xf, yf, x, xs, y, ys, s, n, 2)
        !
        !---- ... and x,y coordinates
        xt1 = seval(st1, x, xs, s, n)
        yt1 = seval(st1, y, ys, s, n)
        xt2 = seval(st2, x, xs, s, n)
        yt2 = seval(st2, y, ys, s, n)
        xb1 = seval(sb1, x, xs, s, n)
        yb1 = seval(sb1, y, ys, s, n)
        xb2 = seval(sb2, x, xs, s, n)
        yb2 = seval(sb2, y, ys, s, n)
        !
        !
        write (*, 1100) xt1, yt1, xt2, yt2, xb1, yb1, xb2, yb2
        1100 format (/' Top breaks: x,y =  ', 2f9.5, 4x, &
                &2f9.5/' Bot breaks: x,y =  ', 2f9.5, 4x, 2f9.5)
        !
        !---- find points adjacent to breaks
        do i = 1, n - 1
            if (s(i)<=st1 .and. s(i + 1)>st1) it1 = i + 1
            if (s(i)<st2 .and. s(i + 1)>=st2) it2 = i
            if (s(i)<=sb1 .and. s(i + 1)>sb1) ib1 = i
            if (s(i)<sb2 .and. s(i + 1)>=sb2) ib2 = i + 1
        enddo
        !
        dsavg = (s(n) - s(1)) / float(n - 1)
        !
        !---- smallest fraction of s increments i+1 and i+2 away from break point
        sfrac = 0.33333
        !
        if (atop/=0.0) then
            st1p = st1 + sfrac * (s(it1) - st1)
            st1q = st1 + sfrac * (s(it1 + 1) - st1)
            if (s(it1)<st1q) then
                !-------- simply move adjacent point to ideal SFRAC location
                xt1new = seval(st1q, x, xs, s, n)
                yt1new = seval(st1q, y, ys, s, n)
                lt1new = .false.
            else
                !-------- make new point at SFRAC location
                xt1new = seval(st1p, x, xs, s, n)
                yt1new = seval(st1p, y, ys, s, n)
                lt1new = .true.
            endif
            !
            st2p = st2 + sfrac * (s(it2) - st2)
            it2q = max(it2 - 1, 1)
            st2q = st2 + sfrac * (s(it2q) - st2)
            if (s(it2)>st2q) then
                !-------- simply move adjacent point
                xt2new = seval(st2q, x, xs, s, n)
                yt2new = seval(st2q, y, ys, s, n)
                lt2new = .false.
            else
                !-------- make new point
                xt2new = seval(st2p, x, xs, s, n)
                yt2new = seval(st2p, y, ys, s, n)
                lt2new = .true.
            endif
        endif
        !
        if (abot/=0.0) then
            sb1p = sb1 + sfrac * (s(ib1) - sb1)
            sb1q = sb1 + sfrac * (s(ib1 - 1) - sb1)
            if (s(ib1)>sb1q) then
                !-------- simply move adjacent point
                xb1new = seval(sb1q, x, xs, s, n)
                yb1new = seval(sb1q, y, ys, s, n)
                lb1new = .false.
            else
                !-------- make new point
                xb1new = seval(sb1p, x, xs, s, n)
                yb1new = seval(sb1p, y, ys, s, n)
                lb1new = .true.
            endif
            !
            sb2p = sb2 + sfrac * (s(ib2) - sb2)
            ib2q = min(ib2 + 1, n)
            sb2q = sb2 + sfrac * (s(ib2q) - sb2)
            if (s(ib2)<sb2q) then
                !-------- simply move adjacent point
                xb2new = seval(sb2q, x, xs, s, n)
                yb2new = seval(sb2q, y, ys, s, n)
                lb2new = .false.
            else
                !-------- make new point
                xb2new = seval(sb2p, x, xs, s, n)
                yb2new = seval(sb2p, y, ys, s, n)
                lb2new = .true.
            endif
        endif
        !
        !c      DSTOP = ABS(S(IT2)-S(IT1))
        !c      DSBOT = ABS(S(IB2)-S(IB1))
        !
        sind = sin(adef)
        cosd = cos(adef)
        !
        !-------------------------------------------------------------------
        !---- initialize accumulator index for new airfoil
        inew = 0
        !
        !---- upper flap surface
        do i = 1, it2
            xbar = x(i) - xf
            ybar = y(i) - yf
            !
            inew = inew + 1
            xnew(inew) = xf + xbar * cosd + ybar * sind
            ynew(inew) = yf - xbar * sind + ybar * cosd
        enddo
        !
        !
        if (atop==0.0) then
            !------ arc length of newly created surface on top of airfoil
            dsnew = abs(adef) * sqrt((xt1 - xf)**2 + (yt1 - yf)**2)
            !
            !------ number of points to be added to define newly created surface
            npadd = int(1.5 * dsnew / dsavg + 1.0)
            !
            if (npadd>0) then
                !------- add new points along the new surface circular arc segment
                dang = adef / float(npadd)
                xbar = xt1 - xf
                ybar = yt1 - yf
                do ip = 1, npadd
                    ang = dang * (float(ip) - 0.5)
                    ca = cos(ang)
                    sa = sin(ang)
                    !
                    inew = inew + 1
                    xnew(inew) = xf + xbar * ca + ybar * sa
                    ynew(inew) = yf - xbar * sa + ybar * ca
                enddo
            endif
            !
        else
            !------ set point in the corner and possibly two adjacent points
            if (lt2new) then
                xbar = xt2new - xf
                ybar = yt2new - yf
                inew = inew + 1
                xnew(inew) = xf + xbar * cosd + ybar * sind
                ynew(inew) = yf - xbar * sind + ybar * cosd
            endif
            !
            inew = inew + 1
            xnew(inew) = xt1
            ynew(inew) = yt1
            !
            if (lt1new) then
                inew = inew + 1
                xnew(inew) = xt1new
                ynew(inew) = yt1new
            endif
            !
        endif
        !
        !
        !---- unchanged portion ahead of flap breaks
        do i = it1, ib1
            inew = inew + 1
            xnew(inew) = x(i)
            ynew(inew) = y(i)
        enddo
        !
        !
        if (abot==0.0) then
            !------ arc length of newly created surface on top of airfoil
            dsnew = abs(adef) * sqrt((xb1 - xf)**2 + (yb1 - yf)**2)
            !
            !------ number of points to be added to define newly created surface
            npadd = int(1.5 * dsnew / dsavg + 1.0)
            !
            if (npadd>0) then
                !------- add new points along the new surface circular arc segment
                dang = adef / float(npadd)
                xbar = xb1 - xf
                ybar = yb1 - yf
                do ip = 1, npadd
                    ang = dang * (float(ip) - 0.5)
                    ca = cos(ang)
                    sa = sin(ang)
                    !
                    inew = inew + 1
                    xnew(inew) = xf + xbar * ca + ybar * sa
                    ynew(inew) = yf - xbar * sa + ybar * ca
                enddo
            endif
            !
        else
            if (lb1new) then
                inew = inew + 1
                xnew(inew) = xb1new
                ynew(inew) = yb1new
            endif
            !
            inew = inew + 1
            xnew(inew) = xb1
            ynew(inew) = yb1
            !
            if (lb2new) then
                xbar = xb2new - xf
                ybar = yb2new - yf
                inew = inew + 1
                xnew(inew) = xf + xbar * cosd + ybar * sind
                ynew(inew) = yf - xbar * sind + ybar * cosd
            endif
            !
        endif
        !
        !
        !---- rotate flap points about the hinge point (XF,YF)
        do i = ib2, n
            !------ rotated part on flap
            xbar = x(i) - xf
            ybar = y(i) - yf
            !
            inew = inew + 1
            xnew(inew) = xf + xbar * cosd + ybar * sind
            ynew(inew) = yf - xbar * sind + ybar * cosd
        enddo
        !
        !
        !---- total number of new points
        nnew = inew
        !
        xtop = xt1
        ytop = yt1
        xbot = xb1
        ybot = yb1
        !
    end subroutine flap
    !*==INSIDE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! FLAP



    function inside(x, y, n, xf, yf)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: xf, yf
        logical :: inside
        real, dimension(*) :: x, y
        !
        ! Local variables
        !
        real :: angle, xb1, xb2, yb1, yb2
        integer :: i, ip
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------
        !     Returns .TRUE. if point XF,YF
        !     is inside contour X(i),Y(i).
        !-------------------------------------
        !
        !---- integrate subtended angle around airfoil perimeter
        angle = 0.0
        do i = 1, n
            ip = i + 1
            if (i==n) ip = 1
            xb1 = x(i) - xf
            yb1 = y(i) - yf
            xb2 = x(ip) - xf
            yb2 = y(ip) - yf
            angle = angle + (xb1 * yb2 - yb1 * xb2)                              &
                    & / sqrt((xb1**2 + yb1**2) * (xb2**2 + yb2**2))
        enddo
        !
        !---- angle = 0 if XF,YF is outside, angle = +/- 2 pi  if XF,YF is inside
        inside = abs(angle)>1.0
        !
    end function inside
    !*==GETXYF.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine getxyf(x, xp, y, yp, s, n, tops, bots, xf, yf)
        use m_spline, only : seval, sinvrt
        use m_userio, only : askr
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: bots, tops, xf, yf
        integer :: n
        real, dimension(*) :: s, x, xp, y, yp
        !
        ! Local variables
        !
        real :: boty, topy, yrel
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        if (xf==-999.0) then
            xf = 0.
            call askr('Enter flap hinge x location^', xf)
        endif
        !
        !---- find top and bottom y at hinge x location
        tops = s(1) + (x(1) - xf)
        bots = s(n) - (x(n) - xf)
        call sinvrt(tops, xf, x, xp, s, n)
        call sinvrt(bots, xf, x, xp, s, n)
        topy = seval(tops, y, yp, s, n)
        boty = seval(bots, y, yp, s, n)
        !
        write (*, 1000) topy, boty
        1000 format (/'  Top    surface:  y =', f8.4, &
                &'     y/t = 1.0'/'  Bottom surface:  y =', f8.4, &
                &'     y/t = 0.0')
        !
        if (yf==-999.0) then
            yf = 999.0
            call askr(&
                    &'Enter flap hinge y location (or 999 to specify y/t)^'&
                    &, yf)
        endif
        !
        if (yf==999.0) then
            yrel = 0.5
            call askr('Enter flap hinge relative y/t location^', yrel)
            yf = topy * yrel + boty * (1.0 - yrel)
        endif
        !
    end subroutine getxyf
    !*==SSS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    !

    !      IF(LGPARM) THEN
    !        CALL NEWPEN(3)
    !        CALL GPARPL(XL0,YL0,0.9*CH,.TRUE.,
    !     &              CHORDB,AREAB,RADBLE,ANGBTE,
    !     &              EI11BA,EI22BA,APX1BA,APX2BA,
    !     &              EI11BT,EI22BT,APX1BT,APX2BT,
    !     &              THICKB,CAMBRB)
    !      ENDIF
    !C
    !      CALL PLFLUSH
    !C
    !      LGEOPL = .TRUE.
    !      NOVER = 0
    !C
    !      RETURN
    !      END





    subroutine sss(ss, s1, s2, del, xbf, ybf, x, xp, y, yp, s, n, iside)
        use m_spline, only : d2val, seval, deval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: del, s1, s2, ss, xbf, ybf
        integer :: iside, n
        real, dimension(*) :: s, x, xp, y, yp
        !
        ! Local variables
        !
        real :: a11, a12, a21, a22, det, ds1, ds2, r1, r1sq, &
                & r1_s1, r2, r2sq, r2_s2, rr, rrsq, rr_s1, rr_s2, &
                & rs1, rs2, rsq, sind, ssgn, stot, x1, x1p, x1pp, &
                & x2, x2p, x2pp, xtot, y1, y1p, y1pp, y2, y2p, &
                & y2pp, ytot
        real, save :: eps
        integer :: iter
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Returns arc length points S1,S2 at flap surface break
        !     locations.  S1 is on fixed airfoil part, S2 is on flap.
        !     The points are defined according to two cases:
        !
        !
        !     If DEL > 0:  Surface will be eliminated in S1 < s < S2
        !
        !     Returns the arc length values S1,S2 of the endpoints
        !     of the airfoil surface segment which "disappears" as a
        !     result of the flap deflection.  The line segments between
        !     these enpoints and the flap hinge point (XBF,YBF) have
        !     an included angle of DEL.  DEL is therefore the flap
        !     deflection which will join up the points at S1,S2.
        !     SS is an approximate arc length value near S1 and S2.
        !     It is used as an initial guess for the Newton loop
        !     for S1 and S2.
        !
        !
        !     If DEL = 0:  Surface will be created at s = S1 = S2
        !
        !     If DEL=0, then S1,S2 will cooincide, and will be located
        !     on the airfoil surface where the segment joining the
        !     point at S1,S2 and the hinge point is perpendicular to
        !     the airfoil surface.  This will be the point where the
        !     airfoil surface must be broken to permit a gap to open
        !     as a result of the flap deflection.
        !----------------------------------------------------------------
        !
        !---- convergence epsilon
        data eps/1.0e-5/
        !
        stot = abs(s(n) - s(1))
        !
        write (*, *) 'sss del ', del
        sind = sin(0.5 * abs(del))
        !
        ssgn = 1.0
        if (iside==1) ssgn = -1.0
        !
        !---- initial guesses for S1, S2
        rsq = (seval(ss, x, xp, s, n) - xbf)**2 + (seval(ss, y, yp, s, n) - ybf)**2
        s1 = ss - (sind * sqrt(rsq) + eps * stot) * ssgn
        s2 = ss + (sind * sqrt(rsq) + eps * stot) * ssgn
        !
        !---- Newton iteration loop
        do iter = 1, 10
            x1 = seval(s1, x, xp, s, n)
            x1p = deval(s1, x, xp, s, n)
            y1 = seval(s1, y, yp, s, n)
            y1p = deval(s1, y, yp, s, n)
            !
            x2 = seval(s2, x, xp, s, n)
            x2p = deval(s2, x, xp, s, n)
            y2 = seval(s2, y, yp, s, n)
            y2p = deval(s2, y, yp, s, n)
            !
            r1sq = (x1 - xbf)**2 + (y1 - ybf)**2
            r2sq = (x2 - xbf)**2 + (y2 - ybf)**2
            r1 = sqrt(r1sq)
            r2 = sqrt(r2sq)
            !
            rrsq = (x1 - x2)**2 + (y1 - y2)**2
            rr = sqrt(rrsq)
            !
            if (r1<=eps * stot .or. r2<=eps * stot) then
                s1 = ss
                s2 = ss
                return
            endif
            !
            r1_s1 = (x1p * (x1 - xbf) + y1p * (y1 - ybf)) / r1
            r2_s2 = (x2p * (x2 - xbf) + y2p * (y2 - ybf)) / r2
            !
            if (sind>0.01) then
                !
                if (rr==0.0) return
                !
                rr_s1 = (x1p * (x1 - x2) + y1p * (y1 - y2)) / rr
                rr_s2 = -(x2p * (x1 - x2) + y2p * (y1 - y2)) / rr
                !
                !------- Residual 1: set included angle via dot product
                rs1 = ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) / rr - sind * r1
                a11 = ((xbf - x1) * (-x1p) + (ybf - y1) * (-y1p))                     &
                        & / rr + ((-x1p) * (x2 - x1) + (-y1p) * (y2 - y1))                 &
                        & / rr - ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) * rr_s1 / rrsq - &
                        & sind * r1_s1
                a12 = ((xbf - x1) * (x2p) + (ybf - y1) * (y2p))                       &
                        & / rr - ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) * rr_s2 / rrsq
                !
                !------- Residual 2: set equal length segments
                rs2 = r1 - r2
                a21 = r1_s1
                a22 = -r2_s2
            else
                !
                !------- Residual 1: set included angle via small angle approximation
                rs1 = (r1 + r2) * sind + (s1 - s2) * ssgn
                a11 = r1_s1 * sind + ssgn
                a12 = r2_s2 * sind - ssgn
                !
                !------- Residual 2: set vector sum of line segments beteen the
                !-       endpoints and flap hinge to be perpendicular to airfoil surface.
                x1pp = d2val(s1, x, xp, s, n)
                y1pp = d2val(s1, y, yp, s, n)
                x2pp = d2val(s2, x, xp, s, n)
                y2pp = d2val(s2, y, yp, s, n)
                !
                xtot = x1 + x2 - 2.0 * xbf
                ytot = y1 + y2 - 2.0 * ybf
                !
                rs2 = xtot * (x1p + x2p) + ytot * (y1p + y2p)
                a21 = x1p * (x1p + x2p) + y1p * (y1p + y2p) + xtot * x1pp + ytot * y1pp
                a22 = x2p * (x1p + x2p) + y2p * (y1p + y2p) + xtot * x2pp + ytot * y2pp
            endif
            !
            det = a11 * a22 - a12 * a21
            ds1 = -(rs1 * a22 - a12 * rs2) / det
            ds2 = -(a11 * rs2 - rs1 * a21) / det
            !
            ds1 = min(ds1, 0.01 * stot)
            ds1 = max(ds1, -.01 * stot)
            ds2 = min(ds2, 0.01 * stot)
            ds2 = max(ds2, -.01 * stot)
            !
            s1 = s1 + ds1
            s2 = s2 + ds2
            if (abs(ds1) + abs(ds2)<eps * stot) goto 11
        enddo
        write (*, *) 'SSS: failed to converge subtending angle points'
        s1 = ss
        s2 = ss
        !
        !
        !---- make sure points are identical if included angle is zero.
        11   if (del==0.0) then
            s1 = 0.5 * (s1 + s2)
            s2 = s1
        endif
        !
    end subroutine sss
    !*==CLIS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine clis(x, xp, y, yp, s, n)
        use m_spline, only : curv
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(*) :: s, x, xp, y, yp
        !
        ! Local variables
        !
        real :: cmax, cv
        integer :: i, imax
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     Displays curvatures at panel nodes.
        !-------------------------------------------------------------------
        !
        cmax = 0.0
        imax = 1
        !
        !---- go over each point, calculating curvature
        write (*, 1050)
        do i = 1, n
            cv = curv(s(i), x, xp, y, yp, s, n)
            write (*, 1100) i, x(i), y(i), cv
            if (abs(cv)>abs(cmax)) then
                cmax = cv
                imax = i
            endif
        enddo
        !
        write (*, 1200) cmax, imax, x(imax), y(imax)
        !
        return
        !
        1050 format (/'  i       x        y           curv')
        !CC             120   0.2134  -0.0234      2025.322
        1100 format (1x, i3, 2f9.4, f14.3)
        1200 format (/' Maximum curvature =', f14.3, '   at  i,x,y  = ', i3, 2f9.4)
    end subroutine clis
    !*==CANG.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CLIS


    subroutine cang(x, y, n, lprint)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: lprint
        integer :: n
        real, dimension(*) :: x, y
        !
        ! Local variables
        !
        real :: amax, angl, crossp, dx1, dx2, dy1, dy2
        integer :: i, imax
        real, save :: pi
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     LPRINT=t:   Displays all panel node corner angles
        !     LPRINT=f:   Displays max panel node corner angle
        !-------------------------------------------------------------------
        data pi/3.1415926535897932384/
        !
        amax = 0.0
        imax = 1
        !
        !---- go over each point, calculating corner angle
        if (lprint) write (*, 1050)
        do i = 2, n - 1
            dx1 = x(i) - x(i - 1)
            dy1 = y(i) - y(i - 1)
            dx2 = x(i) - x(i + 1)
            dy2 = y(i) - y(i + 1)
            !
            !------ allow for doubled points
            if (dx1==0.0 .and. dy1==0.0) then
                dx1 = x(i) - x(i - 2)
                dy1 = y(i) - y(i - 2)
            endif
            if (dx2==0.0 .and. dy2==0.0) then
                dx2 = x(i) - x(i + 2)
                dy2 = y(i) - y(i + 2)
            endif
            !
            crossp = (dx2 * dy1 - dy2 * dx1)                                     &
                    & / sqrt((dx1**2 + dy1**2) * (dx2**2 + dy2**2))
            angl = asin(crossp) * (180.0 / pi)
            if (lprint) write (*, 1100) i, x(i), y(i), angl
            if (abs(angl)>abs(amax)) then
                amax = angl
                imax = i
            endif
        enddo
        !
        write (*, 1200) amax, imax, x(imax), y(imax)
        !
        return
        !
        1050 format (/'  i       x        y      angle')
        !CC             120   0.2134  -0.0234   25.322
        1100 format (1x, i3, 2f9.4, f9.3)
        1200 format (/' Maximum panel corner angle =', f7.3, '   at  i,x,y  = ', &
                & i3, 2f9.4)
    end subroutine cang
end module m_geutil
