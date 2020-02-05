module m_aero_old
    implicit none
contains
    !*==SETIAERO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
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
    !--- Aero data stored for one or more radial aerodynamic sections
    !
    !-- aero data quantities for each defined radial aerodynamic section
    !  NAERO     Number of aerodynamic datasets defined (NAERO>=1)
    !  XIAERO    Radial station r/R where aero dataset is defined
    !  AERODATA  Aerodynamic definition of the blade section at XIAERO
    !            AERODATA( 1,x) = A0 (angle of zero lift)
    !            AERODATA( 2,x) = CLMAX (Max CL)
    !            AERODATA( 3,x) = CLMIN (Min CL)
    !            AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
    !            AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
    !            AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
    !            AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
    !            AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
    !            AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
    !            AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
    !            AERODATA(11,x) = REREF (reference Reynold's number)
    !            AERODATA(12,x) = REXP (Reynold's number exponent cd~re^rexp)
    !            AERODATA(13,x) = MCRIT (critical Mach #)
    !            AERODATA(14,x) = TOC (thickness/chord)
    !            AERODATA(15,x) = DCDCL2S (Scndary, annulus drag param d(Cd)/dCL^2)
    !=========================================================================
    !
    !     Version 070-ES1
    !     Philip Carter, Esotec Developments, February 2009
    !     philip (at) esotec (dot) org
    !
    !     Changes from 0.70:
    !
    !     CL/CD plotting fixed (missing argument).
    !     READ and WRIT commands fixed (LU).
    !     FLIP command to toggle to/from neg-BGam parameters.
    !     DISP formats repaired. Message when set for neg-BGam.
    !     Interpolated A0 (AZERO) stored by GETCLCDCM.
    !     Multi-Re plotting with constant Mach or vice-versa (AEROPLT2).
    !     Modified plotting interface.
    !     All plotting functionality duplicated in EDIT.
    !     PLOTMACH and PLOTREYN subroutines to avoid code duplication.
    !     HARD and ANNO fixed.
    !     Disk index displayed with data for multi-disk cases.
    !     Disk index added to plot titles for multi-disk cases.
    !     Mach constant or Re constant written to plot legends.
    !     Various fixes to control structure and cosmetics.
    !
    !=========================================================================
    !
    !

    subroutine setiaero
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: i, n, nr
        real :: xi
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Sets up indices referring to aero section for
        !     each radial station
        !--------------------------------------------------
        !
        !--- Find lower index of aero data sections XIAERO(N) bounding XI=YRC/RTIP
        do nr = 1, nrotor
            do i = 1, nrc
                iaero(i, nr) = 1
                do n = 1, naero(nr)
                    xi = yrc(i, nr) / rtip(nr)
                    if (xiaero(n, nr)<=xi) iaero(i, nr) = n
                enddo
            enddo
        enddo
    end subroutine setiaero
    !*==GETAERO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine getaero(nr, n, xisect, a0, clmax, clmin, dclda, dclda_stall, &
            & dcl_stall, cdmin, cldmin, dcdcl2, dcdcl2s, cmcon, &
            & mcrit, toc, reref, rexp)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: a0, cdmin, cldmin, clmax, clmin, cmcon, dcdcl2, &
                & dcdcl2s, dclda, dclda_stall, dcl_stall, mcrit, &
                & reref, rexp, toc, xisect
        integer :: n, nr
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------
        !     Gets aero data from stored section array
        !
        !   AERODATA    Aerodynamic definition of the blade section at XIAERO
        !               AERODATA( 1,x) = A0 (angle of zero lift)
        !               AERODATA( 2,x) = CLMAX (Max CL)
        !               AERODATA( 3,x) = CLMIN (Min CL)
        !               AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
        !               AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
        !               AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
        !               AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
        !               AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
        !               AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
        !               AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
        !               AERODATA(11,x) = REREF (reference Reynold's number)
        !               AERODATA(12,x) = REXP (Reynold's number exponent cd~re^rexp)
        !               AERODATA(13,x) = MCRIT (critical Mach #)
        !               AERODATA(14,x) = TOC (thickness/chord)
        !               AERODATA(15,x) = DCDCL2S (Secondary, annulus drag param d(Cd)/dCL^2)
        !---------------------------------------------
        !
        if (nr<1 .or. nr>nrotor) then
            write (*, *) 'Error: blade index of aero section out of bounds'
            return
        endif
        if (n<1 .or. n>naero(nr)) then
            write (*, *) 'Error: index of aero section out of bounds'
            return
        endif
        !
        a0 = aerodata(1, n, nr)
        clmax = aerodata(2, n, nr)
        clmin = aerodata(3, n, nr)
        dclda = aerodata(4, n, nr)
        dclda_stall = aerodata(5, n, nr)
        dcl_stall = aerodata(6, n, nr)
        cdmin = aerodata(7, n, nr)
        cldmin = aerodata(8, n, nr)
        dcdcl2 = aerodata(9, n, nr)
        cmcon = aerodata(10, n, nr)
        reref = aerodata(11, n, nr)
        rexp = aerodata(12, n, nr)
        mcrit = aerodata(13, n, nr)
        toc = aerodata(14, n, nr)
        dcdcl2s = aerodata(15, n, nr)
        xisect = xiaero(n, nr)
        !
    end subroutine getaero
    !*==PUTAERO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine putaero(nr, n, xisect, a0, clmax, clmin, dclda, dclda_stall, &
            & dcl_stall, cdmin, cldmin, dcdcl2, dcdcl2s, cmcon, &
            & mcrit, toc, reref, rexp)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: a0, cdmin, cldmin, clmax, clmin, cmcon, dcdcl2, &
                & dcdcl2s, dclda, dclda_stall, dcl_stall, mcrit, &
                & reref, rexp, toc, xisect
        integer :: n, nr
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------
        !     Puts aero data into stored section array at index N
        !
        !   AERODATA    Aerodynamic definition of the blade section at XIAERO
        !               AERODATA( 1,x) = A0 (angle of zero lift)
        !               AERODATA( 2,x) = CLMAX (Max CL)
        !               AERODATA( 3,x) = CLMIN (Min CL)
        !               AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
        !               AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
        !               AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
        !               AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
        !               AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
        !               AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
        !               AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
        !               AERODATA(11,x) = REREF (reference Reynold's number)
        !               AERODATA(12,x) = REXP (Reynold's number exponent cd~re^rexp)
        !               AERODATA(13,x) = MCRIT (critical Mach #)
        !               AERODATA(14,x) = TOC (thickness/chord)
        !               AERODATA(15,x) = DCDCL2S (Secondary, annulus drag param d(Cd)/dCL^2)
        !--------------------------------------------------------
        !
        if (nr<1 .or. nr>nrx) then
            write (*, *) 'Error: blade index of aero section out of bounds'
            return
        endif
        if (n<1) then
            write (*, *) 'Error: index of aero section out of bounds'
            return
        endif
        if (n>nax) then
            write (*, *) 'Too many aero sections defined...'
            return
        endif
        !
        aerodata(1, n, nr) = a0
        aerodata(2, n, nr) = clmax
        aerodata(3, n, nr) = clmin
        aerodata(4, n, nr) = dclda
        aerodata(5, n, nr) = dclda_stall
        aerodata(6, n, nr) = dcl_stall
        aerodata(7, n, nr) = cdmin
        aerodata(8, n, nr) = cldmin
        aerodata(9, n, nr) = dcdcl2
        aerodata(10, n, nr) = cmcon
        aerodata(11, n, nr) = reref
        aerodata(12, n, nr) = rexp
        aerodata(13, n, nr) = mcrit
        aerodata(14, n, nr) = toc
        aerodata(15, n, nr) = dcdcl2s
        xiaero(n, nr) = xisect
        !
    end subroutine putaero
    !*==SORTAR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine sortar(ns, s, w, ndim)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ndim, ns
        real, dimension(ns) :: s
        real, dimension(ndim, ns) :: w
        !
        ! Local variables
        !
        logical :: done
        integer :: ipass, l, n, np
        real :: temp
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !---- sort arrays by S values
        !     Orders data monotonically increasing in S(i)
        !----------------------------------------------------
        !
        do ipass = 1, 500
            done = .true.
            do n = 1, ns - 1
                np = n + 1
                if (s(np)<s(n)) then
                    temp = s(np)
                    s(np) = s(n)
                    s(n) = temp
                    do l = 1, ndim
                        temp = w(l, np)
                        w(l, np) = w(l, n)
                        w(l, n) = temp
                    enddo
                    done = .false.
                endif
            enddo
            if (done) goto 10
        enddo
        stop 'SORTAR failed'
        !
    10   end subroutine sortar
    !*==GETCLCDCM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SORTAR


    !*************************************************************************
    !  Interpolated aero section properties functions
    !  These routines implement a functional representation of the
    !  blade aero properties (CL,CD,CM) vs ALFA
    !*************************************************************************


    subroutine getclcdcm(nr, is, xi, alf, w, rey, secsig, secstagr, clift, &
            & cl_alf, cl_w, clmax, clmin, dcl_stall, stallf, &
            & cdrag, cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: alf, cdrag, cd_alf, cd_rey, cd_w, clift, clmax, &
                & clmin, cl_alf, cl_w, cmom, cm_al, cm_w, &
                & dcl_stall, rey, secsig, secstagr, w, xi
        integer :: is, nr
        logical :: stallf
        !
        ! Local variables
        !
        real :: a0, a02, cdmin, cdrag2, cd_alf2, cd_rey2, cd_w2, &
                & cldmin, clift2, clmax2, clmin2, cl_alf2, cl_w2, &
                & cmcon, cmom2, cm_al2, cm_w2, dcdcl2, dcdcl2s, &
                & dclda, dclda_stall, dcl_stall2, frac, mcrit, &
                & reref, rexp, toc, xisect1, xisect2
        integer :: n
        logical :: stallf2
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------
        !     CL(alpha),
        !      CD(alpha),
        !       CM(alpha) interpolation function for blade at station IS at XI=r/R
        !-------------------------------------------------------------
        !
        if (xi<0.0 .or. xi>1.0) write (*, *)                             &
                &'Undefined section XI in GETCLCDCM '&
                &, xi
        !
        !--- Check for installed aero data section index
        n = iaero(is, nr)
        if (n<1 .or. n>naero(nr)) then
            !
            if (naero(nr)>1) then
                !--- Find lower index of aero data sections XIAERO(N) bounding XI
                do n = 1, naero(nr)
                    if (xiaero(n, nr)>xi) goto 10
                    !c          write(*,*) 'getcl iaero= ',N,' is= ',is,xiaero(N),xi
                    iaero(is, nr) = n
                enddo
                write (*, *) 'Aero section not found for station ', xi
            endif
            !
            n = 1
            iaero(is, nr) = n
        endif
        !
        !--- Get section aero data from stored section array
        10   a0 = aerodata(1, n, nr)
        clmax = aerodata(2, n, nr)
        clmin = aerodata(3, n, nr)
        dclda = aerodata(4, n, nr)
        dclda_stall = aerodata(5, n, nr)
        dcl_stall = aerodata(6, n, nr)
        cdmin = aerodata(7, n, nr)
        cldmin = aerodata(8, n, nr)
        dcdcl2 = aerodata(9, n, nr)
        cmcon = aerodata(10, n, nr)
        reref = aerodata(11, n, nr)
        rexp = aerodata(12, n, nr)
        mcrit = aerodata(13, n, nr)
        toc = aerodata(14, n, nr)
        dcdcl2s = aerodata(15, n, nr)
        xisect1 = xiaero(n, nr)
        !
        !--- Get data for inner bounding aero section
        call clcdcm(alf, w, rey, vso, secsig, secstagr, clift, cl_alf, cl_w, &
                & stallf, cdrag, cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w, a0, &
                & clmax, clmin, dclda, dclda_stall, dcl_stall, cdmin, cldmin, &
                & dcdcl2, cmcon, mcrit, reref, rexp, toc, dcdcl2s)
        !
        !--- Check for another bounding section, if not we are done,
        !    if we have another section linearly interpolate data to station IS
        if (n<naero(nr)) then
            xisect2 = xiaero(n + 1, nr)
            frac = (xi - xisect1) / (xisect2 - xisect1)
            if (frac<=0.0 .or. frac>1.0) then
                !c         write(*,*) 'CL n,is,xi,frac = ',n,is,xi(is),frac
            endif
            !
            !--- A02 sustituted for A0 in the following (2 places),
            !    to permit A0 interpolation for storage in AZERO
            !
            a02 = aerodata(1, n + 1, nr)
            clmax2 = aerodata(2, n + 1, nr)
            clmin2 = aerodata(3, n + 1, nr)
            dclda = aerodata(4, n + 1, nr)
            dclda_stall = aerodata(5, n + 1, nr)
            dcl_stall2 = aerodata(6, n + 1, nr)
            cdmin = aerodata(7, n + 1, nr)
            cldmin = aerodata(8, n + 1, nr)
            dcdcl2 = aerodata(9, n + 1, nr)
            cmcon = aerodata(10, n + 1, nr)
            reref = aerodata(11, n + 1, nr)
            rexp = aerodata(12, n + 1, nr)
            mcrit = aerodata(13, n + 1, nr)
            toc = aerodata(14, n + 1, nr)
            dcdcl2s = aerodata(15, n + 1, nr)
            !
            !--- Get data for outer bounding aero section
            call clcdcm(alf, w, rey, vso, secsig, secstagr, clift2, cl_alf2, cl_w2, &
                    & stallf2, cdrag2, cd_alf2, cd_w2, cd_rey2, cmom2, cm_al2, &
                    & cm_w2, a02, clmax2, clmin2, dclda, dclda_stall, &
                    & dcl_stall2, cdmin, cldmin, dcdcl2, cmcon, mcrit, reref, &
                    & rexp, toc, dcdcl2s)
            !--- Interpolate aero data to blade station
            stallf = stallf .or. stallf2
            clift = (1.0 - frac) * clift + frac * clift2
            cl_alf = (1.0 - frac) * cl_alf + frac * cl_alf2
            cl_w = (1.0 - frac) * cl_w + frac * cl_w2
            clmax = (1.0 - frac) * clmax + frac * clmax2
            clmin = (1.0 - frac) * clmin + frac * clmin2
            dcl_stall = (1.0 - frac) * dcl_stall + frac * dcl_stall2
            !
            cmom = (1.0 - frac) * cmom + frac * cmom2
            cm_al = (1.0 - frac) * cm_al + frac * cm_al2
            cm_w = (1.0 - frac) * cm_w + frac * cm_w2
            !
            cdrag = (1.0 - frac) * cdrag + frac * cdrag2
            cd_alf = (1.0 - frac) * cd_alf + frac * cd_alf2
            cd_w = (1.0 - frac) * cd_w + frac * cd_w2
            cd_rey = (1.0 - frac) * cd_rey + frac * cd_rey2
            a0 = (1.0 - frac) * a0 + frac * a02
        endif
        !
        azero(is, nr) = a0
        !
    end subroutine getclcdcm
    !*==GETALF.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    !GETCLCDCM



    subroutine getalf(nr, is, xi, secsig, secstagr, clift, w, alf, alf_cl, &
            & alf_w, stallf)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: alf, alf_cl, alf_w, clift, secsig, secstagr, w, &
                & xi
        integer :: is, nr
        logical :: stallf
        !
        ! Local variables
        !
        real :: a0, cdrag, cd_alf, cd_rey, cd_w, clmax, clmin, &
                & cltemp, cl_alf, cl_w, cmom, cm_al, cm_w, dalf, &
                & dcl_stall, rey
        real, save :: eps
        integer :: iter
        integer, save :: niter
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------
        !     Inverse alpha(CL) function
        !     Uses Newton-Raphson iteration to get ALF from CL function
        !------------------------------------------------------------
        data niter/10/
        data eps/1.0e-5/
        !
        stallf = .false.
        !
        !---HHY A0 is now an aero section property
        a0 = aerodata(1, is, nr)
        rey = 0.0
        !
        alf = a0
        do iter = 1, niter
            call getclcdcm(nr, is, xi, alf, w, rey, secsig, secstagr, cltemp, &
                    & cl_alf, cl_w, clmax, clmin, dcl_stall, stallf, cdrag, &
                    & cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w)
            !c      IF(STALLF) GO TO 20
            dalf = -(cltemp - clift) / cl_alf
            alf = alf + dalf
            alf_cl = 1.0 / cl_alf
            alf_w = -cl_w / cl_alf
            if (abs(dalf)<eps) return
        enddo
        !
        write (*, *) 'GETALF: alpha(CL) function inversion failed'
        !      write(*,*) 'is,clift  ',is,clift
        !      write(*,*) 'abs(dalf) ',abs(dalf)
        !      write(*,*) 'cl_alf    ',cl_alf
        !
    end subroutine getalf
    !*==CLCDCM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GETALF



    !*************************************************************************
    !  Basic aero section properties functions
    !  These routines implement a functional representation of the
    !  blade section aero properties (CL,CD,CM) vs ALF
    !*************************************************************************

    subroutine clcdcm(alf, w, rey, vso, secsig, secstagr, clift, cl_alf, cl_w, &
            & stallf, cdrag, cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w, &
            & a0, clmax, clmin, dclda, dclda_stall, dcl_stall, &
            & cdmin, cldmin, dcdcl2, cmcon, mcrit, reref, rexp, toc, &
            & dcdcl2s)
        !------------------------------------------------------------
        !     CL(alpha) function
        !     Note that in addition to setting CLIFT and its derivatives
        !     CLMAX and CLMIN (+ and - stall CL's) are set in this routine
        !     In the compressible range the stall CL is reduced by a factor
        !     proportional to Mcrit-Mach.  Stall limiting for compressible
        !     cases begins when the compressible drag added CDC > CDMstall
        !------------------------------------------------------------
        !     CD(alpha) function - presently CD is assumed to be a sum
        !     of profile drag + stall drag + compressibility drag
        !     In the linear lift range drag is CD0 + quadratic function of CL-CLDMIN
        !     In + or - stall an additional drag is added that is proportional
        !     to the extent of lift reduction from the linear lift value.
        !     Compressible drag is based on adding drag proportional to
        !     (Mach-Mcrit_eff)^MEXP
        !------------------------------------------------------------
        !     CM(alpha) function - presently CM is assumed constant,
        !     varying only with Mach by Prandtl-Glauert scaling
        !------------------------------------------------------------
        !
        ! use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: a0, alf, cdmin, cdrag, cd_alf, cd_rey, cd_w, &
                & cldmin, clift, clmax, clmin, cl_alf, cl_w, cmcon, &
                & cmom, cm_al, cm_w, dcdcl2, dcdcl2s, dclda, &
                & dclda_stall, dcl_stall, mcrit, reref, rexp, rey, &
                & secsig, secstagr, toc, vso, w
        logical :: stallf
        !
        ! Local variables
        !
        real :: cdc, cdcl2, cdc_alf, cdc_w, cdmdd, cdmfactor, &
                & cdmstall, cla, cla_alf, cla_w, clfactor, cllim, &
                & cllim_cla, clmaxm, clmfactor, clminm, clmn, clmx, &
                & critmach, critmach_alf, critmach_w, dcd, dcdx, &
                & dcd_alf, dcd_w, dmdd, dmstall, fac, fac_w, &
                & fstall, mach, mach_w, mexp, msq, msq_w, pgrt, &
                & pgrt_w, rcorr, rcorr_rey
        real :: ecmax, ecmin
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- Factors for compressibility drag model, HHY 10/23/00
        !     Mcrit is set by user
        !     Effective Mcrit is Mcrit_eff = Mcrit - CLMFACTOR*(CL-CLDmin) - DMDD
        !     DMDD is the delta Mach to get CD=CDMDD (usually 0.0020)
        !     Compressible drag is CDC = CDMFACTOR*(Mach-Mcrit_eff)^MEXP
        !     CDMstall is the drag at which compressible stall begins
        !
        cdmfactor = 10.0
        clmfactor = 0.25
        mexp = 3.0
        cdmdd = 0.0020
        cdmstall = 0.1000
        !
        !---- Prandtl-Glauert compressibility factor
        msq = w * w / vso**2
        msq_w = 2.0 * w / vso**2
        if (msq>=1.0) then
            write (*, *) 'CLFUNC: Local Mach number limited to 0.99, was ', &
                    & msq
            msq = 0.99
            msq_w = 0.
        endif
        pgrt = 1.0 / sqrt(1.0 - msq)
        pgrt_w = 0.5 * msq_w * pgrt**3
        !
        !---- Mach number and dependence on velocity
        mach = sqrt(msq)
        mach_w = 0.0
        if (mach/=0.0) mach_w = 0.5 * msq_w / mach
        !
        !------------------------------------------------------------
        !--- Generate CLFACTOR for cascade effects from section solidity
        clfactor = 1.0
        if (secsig>0.0) call getclfactor(secsig, secstagr, clfactor)
        !
        !------------------------------------------------------------
        !--- Generate CL from dCL/dAlpha and Prandtl-Glauert scaling
        cla = dclda * pgrt * (alf - a0) * clfactor
        cla_alf = dclda * pgrt * clfactor
        cla_w = dclda * pgrt_w * (alf - a0) * clfactor
        !
        !ccccccccc
        !      WRITE(*,*)'CL Factor   ALF   A0   DCLDA  CLA'
        !      WRITE(*,*) CLFACTOR,ALF,A0,DCLDA,CLA
        !
        !C--- Effective CLmax is limited by Mach effects
        !    reduces CLmax to match the CL of onset of serious compressible drag
        clmx = clmax
        clmn = clmin
        dmstall = (cdmstall / cdmfactor)**(1.0 / mexp)
        clmaxm = max(0.0, (mcrit + dmstall - mach) / clmfactor) + cldmin
        clmax = min(clmax, clmaxm)
        clminm = min(0.0, -(mcrit + dmstall - mach) / clmfactor) + cldmin
        clmin = max(clmin, clminm)
        !
        !--- CL limiter function (turns on after +-stall
        ecmax = dexp(min(200.0d0, dble((cla - clmax) / dcl_stall)))
        ecmin = dexp(min(200.0d0, dble((clmin - cla) / dcl_stall)))
        cllim = dcl_stall * dlog((1.0d0 + ecmax) / (1.0d0 + ecmin))
        cllim_cla = ecmax / (1.0 + ecmax) + ecmin / (1.0 + ecmin)
        !
        !      if(CLLIM.GT.0.001) then
        !      write(*,999) 'cla,cllim,ecmax,ecmin ',cla,cllim,ecmax,ecmin
        !      endif
        ! 999  format(a,2(1x,f10.6),3(1x,d12.6))
        !
        !--- Subtract off a (nearly unity) fraction of the limited CL function
        !    This sets the dCL/dAlpha in the stalled regions to 1-FSTALL of that
        !    in the linear lift range
        fstall = dclda_stall / dclda
        clift = cla - (1.0 - fstall) * cllim
        cl_alf = cla_alf - (1.0 - fstall) * cllim_cla * cla_alf
        cl_w = cla_w - (1.0 - fstall) * cllim_cla * cla_w
        !
        stallf = .false.
        if (clift>clmax) stallf = .true.
        if (clift<clmin) stallf = .true.
        !
        !
        !------------------------------------------------------------
        !--- CM from CMCON and Prandtl-Glauert scaling
        cmom = pgrt * cmcon
        cm_al = 0.0
        cm_w = pgrt_w * cmcon
        !
        !
        !------------------------------------------------------------
        !--- CD from profile drag, stall drag and compressibility drag
        !
        !---- Reynolds number scaling factor
        if (rey<=0) then
            rcorr = 1.0
            rcorr_rey = 0.0
        else
            rcorr = (rey / reref)**rexp
            rcorr_rey = rexp / rey
        endif
        !
        !--- Include quadratic lift drag terms from airfoil and annulus
        !
        !      CDCL2 = DCDCL2 + DCDCL2S
        cdcl2 = dcdcl2
        ! no chance of getting messed up...
        !
        !--- In the basic linear lift range drag is a function of lift
        !    CD = CD0 (constant) + quadratic with CL)
        cdrag = (cdmin + cdcl2 * (clift - cldmin)**2) * rcorr
        cd_alf = (2.0 * cdcl2 * (clift - cldmin) * cl_alf) * rcorr
        cd_w = (2.0 * cdcl2 * (clift - cldmin) * cl_w) * rcorr
        cd_rey = cdrag * rcorr_rey
        !
        !--- Post-stall drag added
        fstall = dclda_stall / dclda
        dcdx = (1.0 - fstall) * cllim / (pgrt * dclda)
        !      write(*,*) 'cla,cllim,fstall,pg,dclda ',cla,cllim,fstall,pg,dclda
        dcd = 2.0 * dcdx**2
        dcd_alf = 4.0 * dcdx * (1.0 - fstall) * cllim_cla * cla_alf / (pgrt * dclda)
        dcd_w = 4.0 * dcdx * ((1.0 - fstall) * cllim_cla * cla_w / (pgrt * dclda)       &
                & - dcd / pgrt * pgrt_w)
        !      write(*,*) 'alf,cl,dcd,dcd_alf,dcd_w ',alf,clift,dcd,dcd_alf,dcd_w
        !
        !--- Compressibility drag (accounts for drag rise above Mcrit with CL effects
        !    CDC is a function of a scaling factor*(M-Mcrit(CL))**MEXP
        !    DMDD is the Mach difference corresponding to CD rise of CDMDD at MCRIT
        dmdd = (cdmdd / cdmfactor)**(1.0 / mexp)
        critmach = mcrit - clmfactor * abs(clift - cldmin) - dmdd
        critmach_alf = -clmfactor * abs(cl_alf)
        critmach_w = -clmfactor * abs(cl_w)
        if (mach<critmach) then
            cdc = 0.0
            cdc_alf = 0.0
            cdc_w = 0.0
        else
            cdc = cdmfactor * (mach - critmach)**mexp
            cdc_w = mexp * mach_w * cdc / mach - mexp * critmach_w * cdc / critmach
            cdc_alf = -mexp * critmach_alf * cdc / critmach
        endif
        !      write(*,*) 'critmach,mach ',critmach,mach
        !      write(*,*) 'cdc,cdc_w,cdc_alf ',cdc,cdc_w,cdc_alf
        !
        fac = 1.0
        fac_w = 0.0
        !--- Although test data does not show profile drag increases due to Mach #
        !    you could use something like this to add increase drag by Prandtl-Glauert
        !    (or any function you choose)
        !c      FAC   = PG
        !c      FAC_W = PG_W
        !--- Total drag terms
        cdrag = fac * cdrag + dcd + cdc
        cd_alf = fac * cd_alf + dcd_alf + cdc_alf
        cd_w = fac * cd_w + fac_w * cdrag + dcd_w + cdc_alf
        cd_rey = fac * cd_rey
        !
    end subroutine clcdcm
    !*==CHKLIM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CLCDCM



    subroutine chklim(n, nstrt, nend, f, fmax)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: fmax
        integer :: n, nend, nstrt
        real, dimension(n) :: f
        !
        ! Local variables
        !
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !--- Get starting and end index for array values F(i) < FMAX
        nstrt = 1
        nend = n
        !--- Look for first point where F(i)<FMAX
        do i = 1, n
            if (f(i)<fmax) exit
        enddo
        nstrt = max(i - 1, 1)
        !--- Look for last point where F(i)<FMAX
        do i = n, 1, -1
            if (f(i)<fmax) exit
        enddo
        nend = min(i + 1, n)
        !
    end subroutine chklim
    !*==OPFILE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine opfile(lu, fname)
        use m_userio, only : askc, asks
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: fname
        integer :: lu
        !
        ! Local variables
        !
        character(1) :: ans, dummy
        character(4) :: comand
        character(128) :: comarg, tmp
        integer :: k, nf
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !
        !---- get filename if it hasn't been already specified
        if (fname==' ') call asks('Enter output filename^', fname)
        !
        !---- try to open file
        open (lu, file = fname, status = 'OLD', err = 50)
        !
        !---- file exists... ask how to proceed
        nf = index(fname, ' ') - 1
        tmp = 'File  ' // fname(1:nf)                                       &
                & // '  exists.  Overwrite / Append / New file ?^'
        call askc(tmp, comand, comarg)
        ans = comand(1:1)
        !
        !---- ask again if reply is invalid
        if (index('OoAaNn', ans)==0) then
            call askc(' O / A / N  ?^', comand, comarg)
            ans = comand(1:1)
            !
            if (index('OoAaNn', ans)==0) then
                !------- Still bad reply. Give up asking and just return
                write (*, *) 'No action taken'
                return
            endif
        endif
        !
        !---- at this point, file is open and reply is valid
        if (index('Oo', ans)/=0) then
            !------ go to beginning of file to overwrite
            rewind (lu)
            goto 60
        elseif (index('Aa', ans)/=0) then
            !------ go to end of file to append
            do k = 1, 12345678
                read (lu, 1000, end = 40) dummy
                1000       format (a)
            enddo
            40      backspace (lu)
            goto 60
        else
            !------ new file... get filename from command argument, or ask if not supplied
            fname = comarg
            if (fname(1:1)==' ') call asks('Enter output filename^', &
                    & fname)
        endif
        !
        !---- at this point, file FNAME is new or is to be overwritten
        50   open (lu, file = fname, status = 'UNKNOWN', err = 90)
        rewind (lu)
        !
        60   return
        !
        90   write (*, *) 'Bad filename.'
    end subroutine opfile
    !*==GETCLFACTOR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! OPFILE


    subroutine getclfactor(sigma, stagger, clfactor)
        use m_spline, only : sevlin
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        real, parameter :: pi = 3.1415926535897932384, dtr = pi / 180.
        !
        ! Dummy arguments
        !
        real :: clfactor, sigma, stagger
        !
        ! Local variables
        !
        real, dimension(11), save :: a0, a1, a2, x
        real :: aa0, aa1, aa2, daa0, daa1, daa2, sigi, stagr
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------------
        !     Calculates multi-plane cascade effect on lift slope as a
        !     function of solidity and stagger angle
        !
        !     Input:  SIGMA      solidity = Bc/(2*pi*r)
        !             STAGGER    stagger angle (from axis to chordline, in rads)
        !
        !     Output:  CLFACTOR  CLmultiplane/CL2D factor
        !
        !     Implements table-driven quadratic fit to Figure 6-29 in
        !     Wallis, Axial Flow Fans and Ducts.
        !------------------------------------------------------------
        !
        !---- Table of quadratic fit coefficients
        data x/0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, &
                & 1.5/
        data a0/0.4755, 0.5255, 0.5722, 0.6142, 0.6647, 0.7016, &
                & 0.7643, 0.8302, 0.8932, 0.9366, 0.9814/
        data a1/ - 0.367495, -0.341941, -0.300058, -0.255883, &
                & -0.200593, -0.114993, -0.118602, -0.130921, -0.133442, &
                & -0.077980, -0.123071/
        data a2/0.489466, 0.477648, 0.453027, 0.430048, 0.381462, &
                & 0.310028, 0.298309, 0.285309, 0.263084, 0.184165, &
                & 0.251594/
        !
        clfactor = 1.0
        if (sigma<=0.6) return
        !
        !---- Interpolate quadratic fit coefficients by 1/solidity
        sigi = 1.0 / sigma
        call sevlin(sigi, a0, x, 11, aa0, daa0)
        call sevlin(sigi, a1, x, 11, aa1, daa1)
        call sevlin(sigi, a2, x, 11, aa2, daa2)
        !
        !---- Only valid for stagger 20deg to 90deg,
        !     Limit low stagger to 20deg value to give constant lift ratio below that
        stagr = stagger
        if (stagr<20.0 * dtr) stagr = 20.0 * dtr
        if (stagr>90.0 * dtr) stagr = 90.0 * dtr
        !
        !---- Quadratic fit for CLFACTOR at this SIGMA as function of STAGGER
        clfactor = aa0 + aa1 * stagr + aa2 * stagr * stagr
        !---- maximum value of lift ratio should be limited to 1.0
        clfactor = min(1.0, clfactor)
        !
    end subroutine getclfactor
end module m_aero_old
