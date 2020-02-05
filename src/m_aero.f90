module m_aero
    implicit none

    integer, private, allocatable :: n_polars(:), n_polar_points(:), i_polars(:)
    real, private, allocatable :: xi_polars(:), polardata(:, :)
contains

    subroutine setiaero
        use i_dfdc
        implicit none
        integer :: i, n, nr
        real :: xi
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
                    if (xi_polars(sum(n_polars(:nr) + n))<=xi) iaero(i, nr) = n
                enddo
            enddo
        enddo
    end subroutine setiaero
    !*==GETAERO.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020

    subroutine getpolar(i_rotor, i_polar, polar)
        integer, intent(in) :: i_rotor, i_polar
        real, allocatable, intent(out) :: polar(:, :)

        integer :: polar_index

        polar_index = i_polar + sum(n_polars(:i_rotor)) - 1
        polar = polardata(i_polars(polar_index):i_polars(polar_index + 1) - 1, :)
    end subroutine getpolar

    subroutine putpolars(m__n_rotors, m__n_polars, m__n_polar_points, m__xi_polars, m__polardata)
        integer, intent(in) :: m__n_rotors, m__n_polars(m__n_rotors)
        integer, intent(in) ::m__n_polar_points(sum(m__n_polars))
        real, intent(in) :: m__xi_polars(sum(m__n_polars)), m__polardata(sum(m__n_polar_points), 4)

        integer, allocatable :: temp_i_polars(:), temp_n_polars(:), temp_n_polar_points(:)
        real, allocatable :: temp_xi_polars(:), temp_polardata(:, :)

        integer :: i

        temp_i_polars = (/1, (1 + sum(m__n_polar_points(:i)), i=1, sum(m__n_polars))/)
        call move_alloc(temp_i_polars, i_polars)

        temp_n_polars = m__n_polars
        call move_alloc(temp_n_polars, n_polars)

        temp_n_polar_points = m__n_polar_points
        call move_alloc(temp_n_polar_points, n_polar_points)

        temp_xi_polars = m__xi_polars
        call move_alloc(temp_xi_polars, xi_polars)

        temp_polardata = m__polardata
        call move_alloc(temp_polardata, polardata)
    end subroutine putpolars

    subroutine getclcdcm(i_rotor, i_station, xi, alf, w, rey, secsig, secstagr, clift, &
            & cl_alf, cl_w, clmax, clmin, dcl_stall, stallf, &
            & cdrag, cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w)
        use i_dfdc
        implicit none

        integer, intent(in) :: i_station, i_rotor
        real, intent(in) :: xi, alf, w, rey, secsig, secstagr

        real, intent(out) :: cdrag, cd_alf, cd_rey, cd_w, clift, clmax, &
                & clmin, cl_alf, cl_w, cmom, cm_al, cm_w, &
                & dcl_stall
        logical, intent(out) :: stallf

        real, allocatable :: polar(:, :)
        logical stallf2
        real :: xi_polar1, xi_polar2
        integer :: i_polar, i_actual
        real :: clift2, cl_alf2, cl_w2, &
                cdrag2, cd_alf2, cd_w2, cd_rey2, &
                cmom2, cm_al2, cm_w2, &
                clmax2, clmin2, frac, dcl_stall2

        if (xi<0.0 .or. xi>1.0) write (*, *) 'Undefined section XI in GETCLCDCM ', xi

        !
        !--- Check for installed aero data section index
        i_polar = iaero(i_station, i_rotor)
        if (i_polar<1 .or. i_polar>n_polars(i_rotor)) then
            !
            if (n_polars(i_rotor)>1) then
                !--- Find lower index of aero data sections XIAERO(N) bounding XI
                do i_polar = 1, n_polars(i_rotor)
                    if (xiaero(i_polar, i_rotor)>xi) goto 10
                    !c          write(*,*) 'getcl iaero= ',N,' is= ',is,xiaero(N),xi
                    iaero(i_station, i_rotor) = i_polar
                enddo
                write (*, *) 'Aero section not found for station ', xi
            endif
            !
            i_polar = 1
            iaero(i_station, i_rotor) = i_polar
        endif

        i_actual = sum(n_polars(:i_rotor)) + i_polar

        !
        !--- Get section aero data from stored section array
        10 call getpolar(i_rotor, i_polar, polar)
        clmax = maxval(polar(:, 2))
        clmin = minval(polar(:, 2))
        xi_polar1 = xi_polars(i_actual)

        !--- Get data for inner bounding aero section
        call clcdcm(alf, w, rey, vso, secsig, secstagr,  &
                clift, cl_alf, cl_w, stallf, &
                cdrag, cd_alf, cd_w, cd_rey, &
                cmom, cm_al, cm_w, &
                polar, .true.)

        !--- Check for another bounding section, if not we are done,
        !    if we have another section linearly interpolate data to station IS
        if (i_polar<n_polars(i_rotor)) then
            xi_polar2 = xi_polars(i_actual + 1)
            frac = (xi - xi_polar1) / (xi_polar2 - xi_polar1)
            if (frac<=0.0 .or. frac>1.0) then
                !c         write(*,*) 'CL n,is,xi,frac = ',n,is,xi(is),frac
            endif

            call getpolar(i_rotor, i_polar + 1, polar)
            clmax2 = maxval(polar(:, 2))
            clmin2 = minval(polar(:, 2))

            !--- Get data for outer bounding aero section
            call clcdcm(alf, w, rey, vso, secsig, secstagr,  &
                    clift2, cl_alf2, cl_w2, stallf2, &
                    cdrag2, cd_alf2, cd_w2, cd_rey2, &
                    cmom2, cm_al2, cm_w2, &
                    polar, .true.)

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
        end if
    end subroutine getclcdcm

    ! TODO: This is literaly copied from original code.
    !       It is only used for some types of runs, where the blade twist is changed.
    !       So long as those run cases are not done, this does not hurt.
    !       However, things will break if you do try to run it like that.
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

    subroutine clcdcm(alf, w, rey, vso, secsig, secstagr, &
            clift, cl_alf, cl_w, stallf, &
            cdrag, cd_alf, cd_w, cd_rey, &
            cmom, cm_al, cm_w, &
            polar, use_corrections, m_crit)
        use i_dfdc, only : pi
        logical :: stallf
        integer :: n_points
        real :: alf, w, rey, clift, cl_alf, cl_w, cdrag, cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w, &
                secsig, secstagr, vso, &
                deltas(4), new_data(3), deriv(3)
        real :: polar(:, :)
        logical, optional :: use_corrections
        real, optional :: m_crit
        real :: cdmfactor, clfactor, clmfactor, mexp, cdmdd, cdmstall, msq, msq_w, pg, pg_w, mach, mach_w, &
                clmax, clmin, dmstall, clmaxm, clminm, dmdd, critmach, critmach_alf, &
                critmach_w, cdc, cdc_alf, cdc_w, fac, fac_w, cdmin, cldmin
        real :: a_max, a_min

        integer :: i_below, i_above
        real :: f, mcrit

        ! Ensure angle of attack is always between -180 and +180 degrees (-pi and +pi radians)
        if (alf < -pi) then
            alf = alf + 2.*pi
        elseif (alf > pi) then
            alf = alf - 2.*pi
        end if

        a_max = maxval(polar(:, 1))
        a_min = minval(polar(:, 1))
        if (a_min < alf .and. alf < a_max) then
            ! Find the indices of the angles of attack in the polar just below and just above the specified one
            i_below = maxloc(polar(:, 1), 1, polar(:, 1) <= alf)
            i_above = minloc(polar(:, 1), 1, polar(:, 1) >= alf)

            if (i_below == i_above) then
                i_above = i_below + 1
            end if

            ! Compute delta values from alpha(i_below) to alpha(i_above)
            deltas = pack(polar(i_above, :), .true.) - pack(polar(i_below, :), .true.)

            ! Interpolation factor
            f = (alf - polar(i_below, 1)) / deltas(1)
        else
            ! Treat cases where request alf falls outside the known range
            ! This is done by assuming cl, cd, and cm are periodic in alpha with period 2*pi
            ! When this assumption is made, alphas outside of the known range can simply be interpolated using the
            ! edge values.
            i_below = size(polar, 1)
            deltas = pack(polar(1, :), .true.) - pack(polar(i_below, :), .true.)
            deltas(1) = deltas(1) + 2.*pi
            f = (alf - polar(i_below, 1)) / deltas(1)
            if (alf <= a_min) then
                f = f + 2.*pi/deltas(1)
            end if
        end if

        ! Compute interpolate data and local slopes
        new_data = pack(polar(i_below, 2:), .true.) + f * deltas(2:)
        deriv = deltas(2:) / deltas(1)

        ! Store results in proper places
        clift = new_data(1)
        cdrag = new_data(2)
        cmom  = new_data(3)

        cl_alf = deriv(1)
        cd_alf = deriv(2)
        cm_al  = deriv(3)

        cd_rey = 0.

        !------------------------------------------------------------
        !--- Generate CLFACTOR for cascade effects from section solidity
        clfactor = 1.0
        if (secsig>0.0) call getclfactor(secsig, secstagr, clfactor)
        !

        clift = clift * clfactor
        cl_alf = cl_alf * clfactor

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Optional Compressibility Corrections !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (.not. present(use_corrections) .or. .not. use_corrections) then
            cl_w = 0.
            cd_w = 0.
            cm_w = 0.

        else
            if (.not. present(m_crit)) then
                mcrit = 0.6
            else
                mcrit = m_crit
            end if

            ! Compute cl at cdmin
            cdmin = minval(polar(:, 3))
            i_above = maxloc(polar(:, 2), 1, polar(:, 3) == cdmin)
            cldmin = polar(i_below, 2)

            cdmfactor = 10.0
            clmfactor = 0.25
            mexp = 3.0
            cdmdd = 0.0020
            cdmstall = 0.1000
            !
            !---- Prandtl-Glauert compressibility factor
            msq = w * w / vso**2
            msq_w = 2.0 * w / vso**2
            if(msq >= 1.0) then
                write(*, *) 'clfunc: Local Mach number limited to 0.99, was ', msq
                msq = 0.99
                msq_w = 0.
            endif
            pg = 1.0 / sqrt(1.0 - msq)
            pg_w = 0.5 * msq_w * pg**3
            !
            !---- Mach number and dependence on velocity
            mach = sqrt(msq)
            mach_w = 0.0
            if (mach/=0.0) mach_w = 0.5 * msq_w / mach
            !
            !------------------------------------------------------------
            !--- Generate cl from dcl/dAlpha and Prandtl-Glauert scaling
            cl_w = clift * pg_w
            clift = clift * pg
            cl_alf = cl_alf * pg
            !
            !--- Effective cLmax is limited by Mach effects
            !    reduces cLmax to match the cl of onset of serious compressible drag
            clmax = maxval(polar(:, 2))
            clmin = minval(polar(:, 2))
            dmstall = (cdmstall / cdmfactor)**(1.0 / mexp)
            clmaxm = max(0.0, (mcrit + dmstall - mach) / clmfactor) + cldmin
            clmax = min(clmax, clmaxm)
            clminm = min(0.0, -(mcrit + dmstall - mach) / clmfactor) + cldmin
            clmin = max(clmin, clminm)

            !------------------------------------------------------------
            !--- cm from cmcon and Prandtl-Glauert scaling
            cmom = pg * cmom
            cm_al = 0.0
            cm_w = pg_w * cmom

            !--- Compressibility drag (accounts for drag rise above Mcrit with cl effects
            !    cdc is a function of a scaling factor*(m-Mcrit(cl))**mexp
            !    dmdd is the Mach difference corresponding to cd rise of cdmdd at mcrit
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
            !c      fac   = pg
            !c      fac_w = pg_w
            !--- Total drag terms
            cdrag = fac * cdrag  + cdc
            cd_alf = fac * cd_alf  + cdc_alf
            cd_w = fac * cd_w + fac_w * cdrag + cdc_w
            cd_rey = fac * cd_rey
        end if


    end subroutine clcdcm

    subroutine getclfactor(sigma, stagger, clfactor)
        use i_dfdc, only : pi, dtr
        use m_spline, only : sevlin
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
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
end module m_aero
