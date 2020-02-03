module m_xutils
    implicit none
contains
    !*==SETEXP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! WAKMOVR


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

    subroutine setexp(s, ds1, smax, nn)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: ds1, smax
        integer :: nn
        real, dimension(nn) :: s
        !
        ! Local variables
        !
        real :: aaa, bbb, ccc, disc, dratio, dresdr, ds, ratio, &
                & res, rnex, rni, sigma, sigman
        integer :: iter, n, nex
        !
        !*** End of declarations rewritten by SPAG
        !
        !........................................................
        !     Sets geometrically stretched array S:
        !
        !       S(i+1) - S(i)  =  r * [S(i) - S(i-1)]
        !
        !       S     (output)  array to be set
        !       DS1   (input)   first S increment:  S(2) - S(1)
        !       SMAX  (input)   final S value:      S(NN)
        !       NN    (input)   number of points
        !........................................................
        !
        sigma = smax / ds1
        nex = nn - 1
        rnex = float(nex)
        rni = 1.0 / rnex
        !
        !---- solve quadratic for initial geometric ratio guess
        aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0
        bbb = rnex * (rnex - 1.0) / 2.0
        ccc = rnex - sigma
        !
        disc = bbb**2 - 4.0 * aaa * ccc
        disc = max(0.0, disc)
        !
        if (nex<=1) then
            stop 'SETEXP: Cannot fill array.  N too small.'
        elseif (nex==2) then
            ratio = -ccc / bbb + 1.0
        else
            ratio = (-bbb + sqrt(disc)) / (2.0 * aaa) + 1.0
        endif
        !
        if (ratio/=1.0) then
            !
            !---- Newton iteration for actual geometric ratio
            do iter = 1, 100
                sigman = (ratio**nex - 1.0) / (ratio - 1.0)
                res = sigman**rni - sigma**rni
                dresdr = rni * sigman**rni * (rnex * ratio**(nex - 1) - sigman)       &
                        & / (ratio**nex - 1.0)
                !
                dratio = -res / dresdr
                ratio = ratio + dratio
                !
                if (abs(dratio)<1.0e-5) goto 11
                !
            enddo
            write (*, *)                                                    &
                    &'SETEXP: Convergence failed.  Continuing anyway ...'
        endif
        !
        !---- set up stretched array using converged geometric ratio
        11   s(1) = 0.0
        ds = ds1
        do n = 2, nn
            s(n) = s(n - 1) + ds
            ds = ds * ratio
        enddo
        !
    end subroutine setexp
    !*==SETEX2.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine setex2(s, ds1, dsn, smax, n)
        !==========================================================
        !     Sets array S stretched so that a prescribed spacing is
        !     obtained at each end.  The interior spacing is a blend
        !     of two geometric stretchings "shot" from each end.
        !
        !       S     (output)  array to be set
        !       DS1   (input)   approximate first S increment:  S(2) - S(1)
        !       DSN   (input)   approximate last  S increment:  S(N) - S(N-1)
        !       SMAX  (input)   final S value:      S(N)
        !       N     (input)   number of points
        !==========================================================
        !
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: ndim = 500
        !
        ! Dummy arguments
        !
        real :: ds1, dsn, smax
        integer :: n
        real, dimension(n) :: s
        !
        ! Local variables
        !
        real, save :: fend
        integer :: i, in
        real, dimension(ndim) :: s1, sn
        real :: sgn, ss1, ssn, wt1, wtn
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !
        !---- with bigger FEND, the actual end increments will get closer
        !-    to DS1 & DSN, but the interior spacing might get screwy.
        data fend/2.0/
        !
        if (n>ndim) stop 'SETEX2:  Array overflow.'
        !
        !---- calculate spacing arrays each having the prescribed end increment
        call setexp(s1, ds1, smax, n)
        call setexp(sn, dsn, smax, n)
        !
        !c      S(1) = 0.
        !c      S(N) = SMAX
        !
        !---- blend spacing arrays with power-function weights
        do i = 1, n
            in = n - i + 1
            ss1 = s1(i)
            ssn = smax - sn(in)
            !
            !------ power function of integer index
            wt1 = float(n - i)**fend
            wtn = float(i - 1)**fend
            !
            !------ power function of coordinate
            !CC     WT1 = (1.0 - SSN/SMAX)**FEND
            !CC     WTN = (      SS1/SMAX)**FEND
            !
            s(i) = (ss1 * wt1 + ssn * wtn) / (wt1 + wtn)
        enddo
        !
        !---- check for monotonicity
        sgn = sign(1.0, s(n) - s(1))
        do i = 2, n
            if (sgn * s(i)<=sgn * s(i - 1)) then
                write (*, *) 'SETEX2: Warning. Returned array not monotonic.'
                return
            endif
        enddo
        !
    end subroutine setex2
    !*==ATANC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETEX2



    function atanc(y, x, thold)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: thold, x, y
        real :: atanc
        !
        ! Local variables
        !
        real :: dtcorr, dthet, thnew
        real, save :: pi, tpi
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------
        !     ATAN2 function with branch cut checking.
        !
        !     Increments position angle of point X,Y from some previous
        !     value THOLD due to a change in position, ensuring that the
        !     position change does not cross the ATAN2 branch cut
        !     (which is in the -x direction).  For example:
        !
        !       ATANC( -1.0 , -1.0 , 0.75*pi )  returns  1.25*pi , whereas
        !       ATAN2( -1.0 , -1.0 )            returns  -.75*pi .
        !
        !     Typically, ATANC is used to fill an array of angles:
        !
        !        THETA(1) = ATAN2( Y(1) , X(1) )
        !        DO i=2, N
        !          THETA(i) = ATANC( Y(i) , X(i) , THETA(i-1) )
        !        END DO
        !
        !     This will prevent the angle array THETA(i) from jumping by
        !     +/- 2 pi when the path X(i),Y(i) crosses the negative x axis.
        !
        !     Input:
        !       X,Y     point position coordinates
        !       THOLD   position angle of nearby point
        !
        !     Output:
        !       ATANC   position angle of X,Y
        !---------------------------------------------------------------
        data pi/3.1415926535897932384/
        data tpi/6.2831853071795864769/
        !
        !---- set new position angle, ignoring branch cut in ATAN2 function for now
        thnew = atan2(y, x)
        dthet = thnew - thold
        !
        !---- angle change cannot exceed +/- pi, so get rid of any multiples of 2 pi
        dtcorr = dthet - tpi * int((dthet + sign(pi, dthet)) / tpi)
        !
        !---- set correct new angle
        atanc = thold + dtcorr
        !
    end function atanc
    !*==HSORT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ATANC




    subroutine hsort(n, a, indx)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(n) :: a
        integer, dimension(n) :: indx
        !
        ! Local variables
        !
        integer :: i, indxt, ir, j, l
        real :: q
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------
        !     Heapsort algorithm.
        !     Returns INDX(.) such that
        !
        !       A(INDX(i)) < A(INDX(i+1))
        !
        !     Stolen from Numerical Recipes.
        !--------------------------------------
        !
        do i = 1, n
            indx(i) = i
        enddo
        !
        if (n<=1) return
        !
        l = n / 2 + 1
        ir = n
        do
            !
            if (l>1) then
                l = l - 1
                indxt = indx(l)
                q = a(indxt)
            else
                indxt = indx(ir)
                q = a(indxt)
                indx(ir) = indx(1)
                !
                ir = ir - 1
                if (ir==1) then
                    indx(1) = indxt
                    return
                endif
            endif
            !
            i = l
            j = l + l
            do
                !
                if (j<=ir) then
                    if (j<ir) then
                        if (a(indx(j))<a(indx(j + 1))) j = j + 1
                    endif
                    if (q<a(indx(j))) then
                        indx(i) = indx(j)
                        !
                        i = j
                        j = j + j
                    else
                        j = ir + 1
                    endif
                    cycle
                endif
                !
                indx(i) = indxt
                goto 100
            enddo
            exit
        100  enddo
    end subroutine hsort
end module m_xutils
