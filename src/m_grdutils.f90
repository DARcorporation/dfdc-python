module m_grdutils
    implicit none
contains
    !*==XSIETA.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CANG


    !      SUBROUTINE ADDXYA(IEL,DX,DY,DF,DA)
    !      INCLUDE 'AIRPAN.INC'
    !C
    !      EFE = EXP(DF)
    !      SAE = SIN(-DA)
    !      CAE = COS(-DA)
    !      DXE = DXEL(IEL)
    !      DYE = DYEL(IEL)
    !C
    !      DXEL(IEL) = DXE*CAE - DYE*SAE + DX
    !      DYEL(IEL) = DXE*SAE + DYE*CAE + DY
    !      DFEL(IEL) = DFEL(IEL) + DF
    !      DAEL(IEL) = DAEL(IEL) - DA
    !C
    !      RETURN
    !      END
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

    subroutine xsieta(xm, x1, x2, x3, x4, ym, y1, y2, y3, y4, xsi, eta, xsi_x, &
            & xsi_y, eta_x, eta_y)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: eta, eta_x, eta_y, x1, x2, x3, x4, xm, xsi, &
                & xsi_x, xsi_y, y1, y2, y3, y4, ym
        !
        ! Local variables
        !
        real :: a11, a12, a21, a22, det, detinv, dr1, dr2, dxs, &
                & et, xs, z1, z2, z3, z4
        integer :: ixy
        !
        !*** End of declarations rewritten by SPAG
        !
        !............................................................
        !     Subroutine for calculating the local bi-linear
        !     coordinates of a marker inside a quadrilateral cell.
        !     Also calculates the local coordinate derivatives.
        !
        !   Input:
        !   ------
        !     XM,YM          marker coordinates
        !     X1,X2,X3,X4    cell vertex x-coordinates
        !     Y1,Y2,Y3,Y4    cell vertex y-coordinates
        !
        !   Output:
        !   -------
        !     XSI,ETA        local cell coordinates at marker
        !     XSI_X          dXSI/dX  at marker
        !     XSI_Y          dXSI/dY  at marker
        !     ETA_X          dETA/dX  at marker
        !     ETA_Y          dETA/dY  at marker
        !
        !............................................................
        !
        !---- initial guess for xsi, eta  (center of cell)
        xs = 0.
        et = 0.
        !
        !---- perform 12 or less Newton iterations for improving the guess
        do ixy = 1, 12
            z1 = (1.0 - xs) * (1.0 - et)
            z2 = (1.0 + xs) * (1.0 - et)
            z3 = (1.0 + xs) * (1.0 + et)
            z4 = (1.0 - xs) * (1.0 + et)
            dr1 = -(z1 * x1 + z2 * x2 + z3 * x3 + z4 * x4 - 4.0 * xm)
            dr2 = -(z1 * y1 + z2 * y2 + z3 * y3 + z4 * y4 - 4.0 * ym)
            a11 = -(1.0 - et) * x1 + (1.0 - et) * x2 + (1.0 + et) * x3 - (1.0 + et) * x4
            a12 = -(1.0 - xs) * x1 - (1.0 + xs) * x2 + (1.0 + xs) * x3 + (1.0 - xs) * x4
            a21 = -(1.0 - et) * y1 + (1.0 - et) * y2 + (1.0 + et) * y3 - (1.0 + et) * y4
            a22 = -(1.0 - xs) * y1 - (1.0 + xs) * y2 + (1.0 + xs) * y3 + (1.0 - xs) * y4
            detinv = 1.0 / (a11 * a22 - a12 * a21)
            dxs = detinv * (dr1 * a22 - a12 * dr2)
            det = detinv * (a11 * dr2 - dr1 * a21)
            xs = xs + dxs
            et = et + det
            if (amax1(abs(dxs), abs(det))<1.0e-4) goto 11
        enddo
        write (6, *) 'Xsi-Eta conv failed at mrkr', xm, ym, dxs, det
        !
        11   xsi = xs
        eta = et
        !
        !---- derivatives are a free bonus from Newton procedure
        xsi_x = 4.0 * a22 * detinv
        eta_x = -4.0 * a21 * detinv
        xsi_y = -4.0 * a12 * detinv
        eta_y = 4.0 * a11 * detinv
        !
    end subroutine xsieta
    !*==FINT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! XSIETA



    function fint(xsi, eta, q1, q2, q3, q4)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: eta, q1, q2, q3, q4, xsi
        real :: fint
        !
        ! Local variables
        !
        real :: z1, z2, z3, z4
        !
        !*** End of declarations rewritten by SPAG
        !
        !...................................................
        !     Interpolates cell corner values Q1,Q2,Q3,Q4
        !     to marker at local cell coodinates XSI,ETA
        !...................................................
        !
        z1 = (1.0 - xsi) * (1.0 - eta)
        z2 = (1.0 + xsi) * (1.0 - eta)
        z3 = (1.0 + xsi) * (1.0 + eta)
        z4 = (1.0 - xsi) * (1.0 + eta)
        !
        fint = (q1 * z1 + q2 * z2 + q3 * z3 + q4 * z4) * 0.25
        !
    end function fint
    !*==XYGRDFIND.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! FINT



    subroutine xygrdfind(xf, yf, ix, x, y, ii, jj, ic, jc)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ic, ii, ix, jc, jj
        real :: xf, yf
        real, dimension(ix, *) :: x, y
        !
        ! Local variables
        !
        integer :: inout, io, ip, jo, jp, nq
        real :: xmax, xmin, ymax, ymin
        real, dimension(5) :: xq, yq
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !..............................................
        !     Last-resort search through all
        !     cells to locate point XF,YF in the grid
        !
        !     Returns indices IC,JC (to IC+1,JC+1) of cell
        !     containing point or 0,0 for point outside
        !     of grid.
        !..............................................
        !
        ic = 0
        jc = 0
        !
        !---- sweep over all cells in grid
        do io = 1, ii - 1
            ip = io + 1
            do jo = 1, jj - 1
                jp = jo + 1
                !
                xq(1) = x(io, jo)
                xq(2) = x(ip, jo)
                xq(3) = x(ip, jp)
                xq(4) = x(io, jp)
                yq(1) = y(io, jo)
                yq(2) = y(ip, jo)
                yq(3) = y(ip, jp)
                yq(4) = y(io, jp)
                !
                !------ initial check to see if point is outside cell limits
                !       (this disqualifies nearly all cells)
                xmax = amax1(xq(1), xq(2), xq(3), xq(4))
                xmin = amin1(xq(1), xq(2), xq(3), xq(4))
                ymax = amax1(yq(1), yq(2), yq(3), yq(4))
                ymin = amin1(yq(1), yq(2), yq(3), yq(4))
                !
                if (xf<=xmax .and. xf>=xmin .and. yf<=ymax .and. yf>=ymin)&
                        & then
                    !
                    !------ Inside test for polygon
                    nq = 4
                    call pnpoly(xf, yf, xq, yq, nq, inout)
                    if (inout>=0) then
                        !------ Found point in cell
                        ic = io
                        jc = jo
                        return
                    endif
                endif
                !
            enddo
        enddo
        !
        !c      WRITE(6,*) 'Point grid location failed. x y =',XF,YF
        !
    end subroutine xygrdfind
    !*==PNPOLY.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! XYGRDFIND


    !
    !     ..................................................................
    !
    !        SUBROUTINE PNPOLY
    !
    !        PURPOSE
    !           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON
    !
    !        USAGE
    !           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )
    !
    !        DESCRIPTION OF THE PARAMETERS
    !           PX      - X-COORDINATE OF POINT IN QUESTION.
    !           PY      - Y-COORDINATE OF POINT IN QUESTION.
    !           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF
    !                     VERTICES OF POLYGON.
    !           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF
    !                     VERTICES OF POLYGON.
    !           N       - NUMBER OF VERTICES IN THE POLYGON.
    !           INOUT   - THE SIGNAL RETURNED:
    !                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,
    !                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,
    !                      1 IF THE POINT IS INSIDE OF THE POLYGON.
    !
    !        REMARKS
    !           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.
    !           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY
    !           OPTIONALLY BE INCREASED BY 1.
    !           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING
    !           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX
    !           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING
    !           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.
    !           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.
    !           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM
    !
    !           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.
    !
    !        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
    !           NONE
    !
    !        METHOD
    !           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT
    !           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE
    !           POINT IS INSIDE OF THE POLYGON.
    !
    !     ..................................................................
    !
    subroutine pnpoly(px, py, xx, yy, n, inout)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: maxdim = 200
        !
        ! Dummy arguments
        !
        integer :: inout, n
        real :: px, py
        real, dimension(n) :: xx, yy
        !
        ! Local variables
        !
        integer :: i, j
        logical :: mx, my, nx, ny
        real :: test
        real, dimension(maxdim) :: x, y
        !
        !*** End of declarations rewritten by SPAG
        !

        !--- Output unit for printed messages
        if (n>maxdim) then
            write (*, 1)
            1       format ('WARNING:', i5, ' too many polygon sides in PNPOLY')
            return
        endif
        !
        do i = 1, n
            x(i) = xx(i) - px
            y(i) = yy(i) - py
        enddo
        !
        inout = -1
        do i = 1, n
            j = 1 + mod(i, n)
            mx = x(i)>=0.0
            nx = x(j)>=0.0
            my = y(i)>=0.0
            ny = y(j)>=0.0
            if (.not.(.not.((my.or.ny) .and. (mx.or.nx)) .or.             &
                    & (mx .and. nx))) then
                if (.not.(my .and. ny .and. (mx .or. nx) .and.             &
                        & .not.(mx .and. nx))) then
                    !
                    test = (y(i) * x(j) - x(i) * y(j)) / (x(j) - x(i))
                    if (test<0.0) then
                    elseif (test==0.0) then
                        inout = 0
                        return
                    elseif (test>0.0) then
                        inout = -inout
                    endif
                else
                    inout = -inout
                endif
            endif
            !
        enddo
        !
    end subroutine pnpoly
    !*==UVGRDC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine uvgrdc(ix, ii, jj, x, y, xpos, ypos, uc, vc)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ii, ix, jj
        real, dimension(ix, *) :: uc, vc, x, y
        real, dimension(*) :: xpos, ypos
        !
        ! Local variables
        !
        real :: aja, det, dxaet, dxaxi, dxdet, dxdxi, dxi, &
                & dyaet, dyaxi, dydet, dydxi, xav, xoo, xop, xpo, &
                & xpp, yav, yoo, yop, ypo, ypp
        integer :: io, ip, jo, jp
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Calculates velocity on cell centers for axisymmetric
        !     streamfunction grid by finte difference of streamfunction
        !     values.
        !
        !   Input:
        !     IX        first grid-array dimension
        !     II,JJ     grid i,j size
        !     X(i,j)    grid z coordinates
        !     Y(i,j)    grid r coordinates
        !     XPOS(i)   grid s values
        !     YPOS(j)   grid e values  (axisymmetric streamfunction)
        !
        !   Output:
        !     UC(i,j)   u velocity at grid cell centers (II-1 x JJ-1)
        !     VC(i,j)   v velocity at grid cell centers (II-1 x JJ-1)
        !
        !.............................................................
        !
        !   Uses solution to Thompson's grid-generation equations with
        !   e(x,r) satisfying the axisymmetric streamfunction equation.
        !   Hence,  e = constant  lines are streamlines.
        !
        !       s_zz + s_rr = 0
        !       e_zz + e_rr = e_r / r
        !
        !   The equations for u(s,e), v(s,e) are
        !
        !       u = z_e / (r J)
        !
        !       v = r_e / (r J)
        !
        !   where
        !
        !         J  =  z_s r_e  -  z_e r_s
        !
        !-------------------------------------------------------------
        !
        !
        !------ go over all streamline cell centers
        do jo = 1, jj - 1
            jp = jo + 1
            !
            !-------- for all cell centers on this streamline
            do io = 1, ii - 1
                ip = io + 1
                !
                xoo = x(io, jo)
                xpo = x(ip, jo)
                xop = x(io, jp)
                xpp = x(ip, jp)
                yoo = y(io, jo)
                ypo = y(ip, jo)
                yop = y(io, jp)
                ypp = y(ip, jp)
                !
                xav = 0.25 * (xoo + xop + xpo + xpp)
                yav = 0.25 * (yoo + yop + ypo + ypp)
                !
                dxaxi = 0.5 * (-xoo - xop + xpo + xpp)
                dyaxi = 0.5 * (-yoo - yop + ypo + ypp)
                dxaet = 0.5 * (-xoo + xop - xpo + xpp)
                dyaet = 0.5 * (-yoo + yop - ypo + ypp)
                !
                dxi = xpos(ip) - xpos(io)
                det = ypos(jp) - ypos(jo)
                !
                dxdxi = dxaxi / dxi
                dydxi = dyaxi / dxi
                dxdet = dxaet / det
                dydet = dyaet / det
                !
                aja = dydet * dxdxi - dxdet * dydxi
                !
                uc(io, jo) = dxdxi / yav / aja
                vc(io, jo) = dydxi / yav / aja
                !
            enddo
            !
        enddo
        !
    end subroutine uvgrdc
    !*==UVGRD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! UVGRDC


    subroutine uvgrd(ix, ii, jj, x, y, xpos, ypos, ug, vg)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ii, ix, jj
        real, dimension(ix, *) :: ug, vg, x, y
        real, dimension(*) :: xpos, ypos
        !
        ! Local variables
        !
        real :: aja, detav, detl, detm, detp, detq, dxde1, &
                & dxde2, dxdet, dxdx1, dxdx2, dxdxi, dxiav, dxil, &
                & dxim, dxip, dxiq, dyde1, dyde2, dydet, dydx1, &
                & dydx2, dydxi, xlo, xmm, xmo, xmp, xol, xom, &
                & xoo, xop, xoq, xpm, xpo, xpp, xqo, yav, ylo, &
                & ymm, ymo, ymp, yol, yom, yoo, yop, yoq, ypm, &
                & ypo, ypp, yqo
        integer :: il, im, io, ip, iq, jl, jm, jo, jp, jq
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     Calculates velocity on axisymmetric streamfunction grid
        !     by finte difference of streamfunction values.
        !
        !   Input:
        !     IX        first grid-array dimension
        !     II,JJ     grid i,j size
        !     X(i,j)    grid z coordinates
        !     Y(i,j)    grid r coordinates
        !     XPOS(i)   grid s values
        !     YPOS(j)   grid e values  (axisymmetric streamfunction)
        !
        !   Output:
        !     UG(i,j)   u velocity at grid points
        !     VG(i,j)   v velocity at grid points
        !
        !.............................................................
        !
        !   Uses solution to Thompson's grid-generation equations with
        !   e(x,r) satisfying the axisymmetric streamfunction equation.
        !   Hence,  e = constant  lines are streamlines.
        !
        !       s_zz + s_rr = 0
        !       e_zz + e_rr = e_r / r
        !
        !   The equations for u(s,e), v(s,e) are
        !
        !       u = z_e / (rho r J)
        !
        !       v = r_e / (rho r J)
        !
        !   where
        !
        !         J  =  z_s r_e  -  z_e r_s
        !
        !-------------------------------------------------------------
        !
        !------ go over all interior streamlines
        do jo = 2, jj - 1
            jm = jo - 1
            jp = jo + 1
            !
            !-------- for all interior points on this streamline
            do io = 2, ii - 1
                im = io - 1
                ip = io + 1
                !
                xmm = x(im, jm)
                xom = x(io, jm)
                xpm = x(ip, jm)
                xmo = x(im, jo)
                xoo = x(io, jo)
                xpo = x(ip, jo)
                xmp = x(im, jp)
                xop = x(io, jp)
                xpp = x(ip, jp)
                ymm = y(im, jm)
                yom = y(io, jm)
                ypm = y(ip, jm)
                ymo = y(im, jo)
                yoo = y(io, jo)
                ypo = y(ip, jo)
                ymp = y(im, jp)
                yop = y(io, jp)
                ypp = y(ip, jp)
                !
                dxim = xpos(io) - xpos(im)
                dxip = xpos(ip) - xpos(io)
                dxiav = 0.5 * (dxim + dxip)
                !
                detm = ypos(jo) - ypos(jm)
                detp = ypos(jp) - ypos(jo)
                detav = 0.5 * (detm + detp)
                !
                dxdet = 0.5 * (xop - xom) / detav
                dydet = 0.5 * (yop - yom) / detav
                dxdxi = 0.5 * (xpo - xmo) / dxiav
                dydxi = 0.5 * (ypo - ymo) / dxiav
                !
                aja = dydet * dxdxi - dxdet * dydxi
                yav = yoo
                !
                ug(io, jo) = dxdxi / yav / aja
                vg(io, jo) = dydxi / yav / aja
                !
            enddo
            !
        enddo
        !
        !
        !-------- Calculate velocities on grid boundaries (uses backward differences)
        !
        !-------- Velocities on lower streamline
        jo = 1
        do io = 2, ii - 1
            ip = io + 1
            im = io - 1
            !
            jp = jo + 1
            jq = jo + 2
            !
            xoo = x(io, jo)
            xop = x(io, jp)
            xoq = x(io, jq)
            yoo = y(io, jo)
            yop = y(io, jp)
            yoq = y(io, jq)
            !
            detp = ypos(jp) - ypos(jo)
            detq = ypos(jq) - ypos(jp)
            !
            dxde1 = (xop - xoo) / detp
            dxde2 = (xoq - xop) / detq
            dyde1 = (yop - yoo) / detp
            dyde2 = (yoq - yop) / detq
            !--- backwards difference at boundary
            dxdet = dxde1 - detp * (dxde2 - dxde1) / (detp + detq)
            dydet = dyde1 - detp * (dyde2 - dyde1) / (detp + detq)
            !
            xmo = x(im, jo)
            xoo = x(io, jo)
            xpo = x(ip, jo)
            ymo = y(im, jo)
            yoo = y(io, jo)
            ypo = y(ip, jo)
            !
            dxim = xpos(io) - xpos(im)
            dxip = xpos(ip) - xpos(io)
            dxiav = 0.5 * (dxim + dxip)
            dxdxi = 0.5 * (xpo - xmo) / dxiav
            dydxi = 0.5 * (ypo - ymo) / dxiav
            !
            aja = dydet * dxdxi - dxdet * dydxi
            yav = yoo
            !
            ug(io, jo) = dxdxi / yav / aja
            vg(io, jo) = dydxi / yav / aja
        enddo
        !
        !-------- Velocities on upper streamline
        jo = jj
        do io = 2, ii - 1
            ip = io + 1
            im = io - 1
            !
            jm = jo - 1
            jl = jo - 2
            !
            xoo = x(io, jo)
            xom = x(io, jm)
            xol = x(io, jl)
            yoo = y(io, jo)
            yom = y(io, jm)
            yol = y(io, jl)
            !
            detm = ypos(jo) - ypos(jm)
            detl = ypos(jm) - ypos(jl)
            dxde1 = (xoo - xom) / detm
            dyde1 = (yoo - yom) / detm
            dxde2 = (xom - xol) / detl
            dyde2 = (yom - yol) / detl
            !--- backwards difference at boundary
            dxdet = dxde1 + detm * (dxde1 - dxde2) / (detm + detl)
            dydet = dyde1 + detm * (dyde1 - dyde2) / (detm + detl)
            !
            xmo = x(im, jo)
            xoo = x(io, jo)
            xpo = x(ip, jo)
            ymo = y(im, jo)
            yoo = y(io, jo)
            ypo = y(ip, jo)
            !
            dxim = xpos(io) - xpos(im)
            dxip = xpos(ip) - xpos(io)
            dxiav = 0.5 * (dxim + dxip)
            dxdxi = 0.5 * (xpo - xmo) / dxiav
            dydxi = 0.5 * (ypo - ymo) / dxiav
            !
            aja = dydet * dxdxi - dxdet * dydxi
            yav = yoo
            !
            ug(io, jo) = dxdxi / yav / aja
            vg(io, jo) = dydxi / yav / aja
        enddo
        !
        !-------- Velocities on(i=1) inlet plane (includes corner points J=1,J=JJ)
        io = 1
        ip = io + 1
        iq = io + 2
        !
        do jo = 1, jj
            jl = jo - 2
            jm = jo - 1
            jp = jo + 1
            jq = jo + 2
            !
            if (jo==1) then
                !---- lower streamline inlet corner
                xoo = x(io, jo)
                xop = x(io, jp)
                xoq = x(io, jq)
                yoo = y(io, jo)
                yop = y(io, jp)
                yoq = y(io, jq)
                !
                detp = ypos(jp) - ypos(jo)
                detq = ypos(jq) - ypos(jp)
                dxde1 = (xop - xoo) / detp
                dxde2 = (xoq - xop) / detq
                dyde1 = (yop - yoo) / detp
                dyde2 = (yoq - yop) / detq
                !--------- backwards difference at boundary
                dxdet = dxde1 - detp * (dxde2 - dxde1) / (detp + detq)
                dydet = dyde1 - detp * (dyde2 - dyde1) / (detp + detq)
            elseif (jo==jj) then
                !---- upper streamline inlet corner
                xoo = x(io, jo)
                xom = x(io, jm)
                xol = x(io, jl)
                yoo = y(io, jo)
                yom = y(io, jm)
                yol = y(io, jl)
                !
                detm = ypos(jo) - ypos(jm)
                detl = ypos(jm) - ypos(jl)
                dxde1 = (xoo - xom) / detm
                dyde1 = (yoo - yom) / detm
                dxde2 = (xom - xol) / detl
                dyde2 = (yom - yol) / detl
                !--------- backwards difference at boundary
                dxdet = dxde1 + detm * (dxde1 - dxde2) / (detm + detl)
                dydet = dyde1 + detm * (dyde1 - dyde2) / (detm + detl)

            else
                xoo = x(io, jo)
                xpo = x(ip, jo)
                xqo = x(iq, jo)
                yoo = y(io, jo)
                ypo = y(ip, jo)
                yqo = y(iq, jo)
                !
                dxip = xpos(ip) - xpos(io)
                dxiq = xpos(iq) - xpos(ip)
                dxdx1 = (xpo - xoo) / dxip
                dxdx2 = (xqo - xpo) / dxiq
                dydx1 = (ypo - yoo) / dxip
                dydx2 = (yqo - ypo) / dxiq
                !---------- 2nd-order backward 3-point difference at boundary
                dxdxi = dxdx1 - dxip * (dxdx2 - dxdx1) / (dxip + dxiq)
                dydxi = dydx1 - dxip * (dydx2 - dydx1) / (dxip + dxiq)
            endif
            !
            xom = x(io, jm)
            xoo = x(io, jo)
            xop = x(io, jp)
            yom = y(io, jm)
            yoo = y(io, jo)
            yop = y(io, jp)
            !
            detm = ypos(jo) - ypos(jm)
            detp = ypos(jp) - ypos(jo)
            detav = 0.5 * (detm + detp)
            dxdet = 0.5 * (xop - xom) / detav
            dydet = 0.5 * (yop - yom) / detav
            !
            aja = dydet * dxdxi - dxdet * dydxi
            yav = yoo
            !
            ug(io, jo) = dxdxi / yav / aja
            vg(io, jo) = dydxi / yav / aja
        enddo
        !
        !-------- Velocities on (i=II) outlet plane (includes corner points J=1,J=JJ)
        io = ii
        im = io - 1
        il = io - 2
        do jo = 1, jj
            jl = jo - 2
            jm = jo - 1
            jp = jo + 1
            jq = jo + 2
            !
            if (jo==1) then
                !---- lower streamline outlet corner
                xoo = x(io, jo)
                xop = x(io, jp)
                xoq = x(io, jq)
                yoo = y(io, jo)
                yop = y(io, jp)
                yoq = y(io, jq)
                !
                detp = ypos(jp) - ypos(jo)
                detq = ypos(jq) - ypos(jp)
                dxde1 = (xop - xoo) / detp
                dxde2 = (xoq - xop) / detq
                dyde1 = (yop - yoo) / detp
                dyde2 = (yoq - yop) / detq
                !--------- backwards difference at boundary
                dxdet = dxde1 - detp * (dxde2 - dxde1) / (detp + detq)
                dydet = dyde1 - detp * (dyde2 - dyde1) / (detp + detq)
            elseif (jo==jj) then
                !---- upper streamline outlet corner
                xoo = x(io, jo)
                xom = x(io, jm)
                xol = x(io, jl)
                yoo = y(io, jo)
                yom = y(io, jm)
                yol = y(io, jl)
                !
                detm = ypos(jo) - ypos(jm)
                detl = ypos(jm) - ypos(jl)
                dxde1 = (xoo - xom) / detm
                dyde1 = (yoo - yom) / detm
                dxde2 = (xom - xol) / detl
                dyde2 = (yom - yol) / detl
                !--------- backwards difference at boundary
                dxdet = dxde1 + detm * (dxde1 - dxde2) / (detm + detl)
                dydet = dyde1 + detm * (dyde1 - dyde2) / (detm + detl)

            else
                !
                xoo = x(io, jo)
                xmo = x(im, jo)
                xlo = x(il, jo)
                yoo = y(io, jo)
                ymo = y(im, jo)
                ylo = y(il, jo)
                !
                dxim = xpos(io) - xpos(im)
                dxil = xpos(im) - xpos(il)
                dxdx1 = (xoo - xmo) / dxim
                dxdx2 = (xmo - xlo) / dxil
                dydx1 = (yoo - ymo) / dxim
                dydx2 = (ymo - ylo) / dxil
                !---------- 2nd-order 3-point difference for tangential velocity
                dxdxi = dxdx1 + dxim * (dxdx1 - dxdx2) / (dxim + dxil)
                dydxi = dydx1 + dxim * (dydx1 - dydx2) / (dxim + dxil)
            endif
            !
            xom = x(io, jm)
            xoo = x(io, jo)
            xop = x(io, jp)
            yom = y(io, jm)
            yoo = y(io, jo)
            yop = y(io, jp)
            !
            detm = ypos(jo) - ypos(jm)
            detp = ypos(jp) - ypos(jo)
            detav = 0.5 * (detm + detp)
            dxdet = 0.5 * (xop - xom) / detav
            dydet = 0.5 * (yop - yom) / detav
            !
            aja = dydet * dxdxi - dxdet * dydxi
            yav = yoo
            !
            ug(io, jo) = dxdxi / yav / aja
            vg(io, jo) = dydxi / yav / aja
        enddo
        !
    end subroutine uvgrd
end module m_grdutils
