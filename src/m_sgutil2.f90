module m_sgutil2
    implicit none
contains
    !*==ISGFIND.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETCRV
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

    integer function isgfind(ss, s, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: ss
        real, dimension(n) :: s
        !
        ! Local variables
        !
        integer :: i, ilow, imid
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- Find closest value in S1 array to SS and return index
        !     Assumes array S in ascending value order and uses binary search to
        !     locate interval with S(I-1)<=SS<S(I)
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
        !---- Check S(I-1) and S(I) for closest value to SS
        isgfind = i
        if (i>1) then
            if (abs(ss - s(i))>abs(ss - s(i - 1))) isgfind = i - 1
        endif
        !
    end function isgfind
    !*==SGCOPY.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine sgcopy(s1, s2, n)
        !---- Copy N elements from S1 to S2
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(n) :: s1, s2
        !
        ! Local variables
        !
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        do i = 1, n
            s2(i) = s1(i)
        enddo
        !
    end subroutine sgcopy
    !*==SGCOPF.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine sgcopf(s1, s2, n, soff1, swt1, soff2, swt2, foff, f1, fx1, x1, n1, &
            & f2, fx2, x2, n2)
        use m_spline, only : seval, deval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: foff, soff1, soff2, swt1, swt2
        integer :: n, n1, n2
        real, dimension(n1) :: f1, fx1, x1
        real, dimension(n2) :: f2, fx2, x2
        real, dimension(n) :: s1, s2
        !
        ! Local variables
        !
        real :: dels1, dev1, ds1, fev1, fev2, res, res_s1
        integer :: i, iter
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- reset array S1 so that  f1(s1) = f2(s2) + foff ,
        !-    with  s1 = S1*SWT1 + SOFF1  and  s2 = S2*SWT2 + SOFF2
        !
        dels1 = s1(n) - s1(1)
        do i = 2, n - 1
            fev2 = seval((s2(i) * swt2 + soff2), f2, fx2, x2, n2) + foff
            !
            do iter = 1, 10
                fev1 = seval((s1(i) * swt1 + soff1), f1, fx1, x1, n1)
                dev1 = deval((s1(i) * swt1 + soff1), f1, fx1, x1, n1)
                res = fev1 - fev2
                res_s1 = dev1 * swt1
                ds1 = -res / res_s1
                s1(i) = s1(i) + ds1
                if (abs(ds1 / dels1)<1.0e-5) goto 20
            enddo
            write (*, *) 'SGCOPF: Convergence failed.  dS/Smax =', &
                    & ds1 / dels1
            !
        20   enddo
        !
    end subroutine sgcopf
    !*==SGAVG.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine sgavg(s1, s2, n, c1)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: c1
        integer :: n
        real, dimension(n) :: s1, s2
        !
        ! Local variables
        !
        real :: f1, f2
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- average arrays S1 S2,
        !-    preserving S1 spacing at right point and S2 spacing at left endpoint
        do i = 1, n
            f1 = float(i - 1) * c1
            f2 = float(n - i)
            s1(i) = (s1(i) * f1 + s2(i) * f2) / (f1 + f2)
            s2(i) = s1(i)
        enddo
        !
    end subroutine sgavg
    !*==SGAVG1.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine sgavg1(s1, s2, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(n) :: s1, s2
        !
        ! Local variables
        !
        real :: f1, f2
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- impose spacing of array S1 on S2 at right endpoint
        do i = 1, n
            f1 = float(i - 1)
            f2 = float(n - i)
            s2(i) = (s1(i) * f1 + s2(i) * f2) / (f1 + f2)
        enddo
        !
    end subroutine sgavg1
    !*==SGAVG2.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine sgavg2(s1, s2, n)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(n) :: s1, s2
        !
        ! Local variables
        !
        real :: f1, f2
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- impose spacing of array S2 on S1 at left endpoint
        do i = 1, n
            f1 = float(i - 1)
            f2 = float(n - i)
            s1(i) = (s1(i) * f1 + s2(i) * f2) / (f1 + f2)
        enddo
        !
    end subroutine sgavg2
    !*==SGSHFT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine sgshft(sbeg, send, sold, snew, s1, s2, n)
        !
        !---- Shift S1 values between limits SBEG and SEND such that
        !     value at SOLD is shifted to SNEW in S2.  S1 values
        !     between SBEG and SOLD are scaled to correspond to
        !     the range SBEG to SNEW. S1 values between SOLD and
        !     SEND are scaled to the range SNEW to SEND.
        !
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real :: sbeg, send, snew, sold
        real, dimension(n) :: s1, s2
        !
        ! Local variables
        !
        real :: fac1, fac2
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        if (sbeg>=send) then
            write (*, *) 'SGSHFT: Bad input range SBEG,SEND ', sbeg, send
            stop
        endif
        !
        if (sbeg>sold .or. send<sold) then
            write (*, *) 'SGSHFT: Bad input range SBEG/SOLD/SEND ', sbeg, &
                    & sold, send
            stop
        endif
        !
        if (sbeg>snew .or. send<snew) then
            write (*, *) 'SGSHFT: Bad input range SBEG/SNEW/SEND ', sbeg, &
                    & snew, send
            stop
        endif
        !
        !---- Shift S1 values in SBEG-SEND range to new values such
        !     that SOLD becomes SNEW (linearly interpolated from both
        !     ends of range).
        fac1 = (snew - sbeg) / (sold - sbeg)
        fac2 = (send - snew) / (send - sold)
        do i = 1, n
            if (s1(i)<=sbeg) then
                s2(i) = s1(i)
            elseif (s1(i)>sbeg .and. s1(i)<=sold) then
                s2(i) = sbeg + (s1(i) - sbeg) * fac1
            elseif (s1(i)>sold .and. s1(i)<send) then
                s2(i) = send - (send - s1(i)) * fac2
            else
                s2(i) = s1(i)
            endif
        enddo
        !
    end subroutine sgshft
    !*==SGRENUM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine sgrenum(s, n1, n2)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: nmax = 2000
        !
        ! Dummy arguments
        !
        integer :: n1, n2
        real, dimension(*) :: s
        !
        ! Local variables
        !
        real, dimension(nmax) :: c, x, xi
        real :: cx1, cx2, frac, ri, t
        integer :: i, j
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------------
        !     Interpolates input array S(1:N1) into new number of points 1:N2.
        !     The interpolation is a spline in the array index.
        !     The first and last elements in S will remain the same.
        !----------------------------------------------------------------------
        !
        if (n1>nmax) stop 'NEWNUM:  Array overflow'
        !
        !---- save old input array
        do i = 1, n1
            x(i) = s(i)
        enddo
        !
        !---- spline X(i)  (set up and solve tridiagonal system on the fly)
        c(1) = 0.5
        xi(1) = 1.5 * (x(2) - x(1))
        do i = 2, n1 - 1
            c(i) = 1.0 / (4.0 - c(i - 1))
            xi(i) = (3.0 * (x(i + 1) - x(i - 1)) - xi(i - 1)) * c(i)
        enddo
        i = n1
        xi(i) = (3.0 * (x(i) - x(i - 1)) - xi(i - 1)) / (2.0 - c(i - 1))
        !
        do i = n1 - 1, 1, -1
            xi(i) = xi(i) - c(i) * xi(i + 1)
        enddo
        !
        !---- evaluate s(i) spline at new points
        do j = 1, n2
            frac = float(j - 1) / float(n2 - 1)
            ri = 1.0 + frac * float(n1 - 1)
            i = min(int(ri), n1 - 1)
            !
            t = ri - float(i)
            cx1 = xi(i) - x(i + 1) + x(i)
            cx2 = xi(i + 1) - x(i + 1) + x(i)
            s(j) = t * x(i + 1) + (1.0 - t) * x(i) + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)
        enddo
        !
        !---- make sure new endpoints are exactly the same as old ones
        s(1) = x(1)
        s(n2) = x(n1)
        !
    end subroutine sgrenum
end module m_sgutil2
