module api
    use, intrinsic :: iso_c_binding, only : c_float, c_double, &
            & c_int, c_bool, c_char
    implicit none
contains

    subroutine init()
        use m_dfdcsubs, only : dfinit
        ! bind(c, name = 'init')
        implicit none
        logical :: ldebug
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        ldebug = .false.
        call dfinit(ldebug)
        !        CALL LOFTINIT1
        !        CALL LOFTINIT2(1,NRX)
    end subroutine init

    subroutine load_case_old(fname)
        use i_dfdc
        use m_dfdcsubs, only : dfload
        ! bind(c, name = 'set_case')
        implicit none
        character(*) :: fname
        logical :: ferror
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        call dfload(fname, ferror)
        !            CALL LOFTINIT2(1,NROTOR)
        if (.not.ferror) loname = name
    end subroutine load_case_old

    subroutine set_case(&
                m__n_rotors, &
                m__vinf, m__vref, m__rpms, &
                m__rho, m__vso, m__rmu, m__alt, &
                m__xdwake, m__nwake, m__lwkrlx, &
                m__n_polars, m__n_polar_points, m__xi_polars, m__polardata, &
                m__xdisk, m__n_blades, m__nrpdef, m__n_stations, m__rotorgeom, &
                m__n_cb, m__cbgeom, &
                m__n_duct, m__ductgeom)
        use i_dfdc
        use m_aero, only : putpolars
        use m_dfdcsubs, only : geproc

        ! inputs
        integer(c_int), intent(in) :: m__n_rotors
        integer(c_int), intent(in) :: m__nwake, &
                m__n_polars(m__n_rotors), m__n_polar_points(sum(m__n_polars)), &
                m__n_blades(m__n_rotors), m__nrpdef(m__n_rotors), m__n_stations(m__n_rotors), &
                m__n_cb, m__n_duct
        real(c_float), intent(in) :: &
                m__vinf, m__vref, m__rpms(m__n_rotors), &
                m__rho, m__vso, m__rmu, m__alt, &
                m__xdwake, &
                m__xi_polars(sum(m__n_polars)), m__polardata(sum(m__n_polar_points), 3), &
                m__xdisk(m__n_rotors), m__rotorgeom(sum(m__n_stations), 3), m__cbgeom(m__n_cb, 2), m__ductgeom(m__n_duct, 2)
        logical(c_bool), intent(in) :: m__lwkrlx

        ! locals
        integer :: i_rotor

        lbldef = .false.
        lconv = .false.
        nrp = 0
        nrc = 0
        do i_rotor = 1, nrx
            irtype(i_rotor) = 0
            nrdef(i_rotor) = 0
            irtypdef(i_rotor) = 0
        end do

        lrspcdef = .false.
        lrspcusr = .false.

        ndobj = 0


        ! basic properties
        qinf = m__vinf
        qref = m__vref
        omega(:m__n_rotors) = pi * m__rpms / 30.0

        rho = m__rho
        vso = m__vso
        rmu = m__rmu
        alth = m__alt

        xdwklen = m__xdwake
        nwake = m__nwake
        lwrlx = m__lwkrlx

        ! Polars
        call putpolars(m__n_rotors, m__n_polars, m__n_polar_points, m__xi_polars, m__polardata)

        ! Rotor geometry
        xdisk(:m__n_rotors) = m__xdisk
        nrbld(:m__n_rotors) = m__n_blades
        nrdef(:m__n_rotors) = m__n_stations
        do i_rotor = 1, m__n_rotors
            yrdef(:m__n_stations(i_rotor), i_rotor) = m__rotorgeom(:, 1)
            chrdef(:m__n_stations(i_rotor), i_rotor) = m__rotorgeom(:, 2)
            betadef(:m__n_stations(i_rotor), i_rotor) = m__rotorgeom(:, 3) * dtr
            irtypdef(i_rotor) = 2
        end do
        nrotor = m__n_rotors
        nrsta = 11 ! this is the default value, but maybe it should be an input?

        ! Centerbody
        ibfrst(1) = 1
        iblast(1) = m__n_cb
        xb(ibfrst(1):iblast(1)) = m__cbgeom(:, 1)
        yb(ibfrst(1):iblast(1)) = m__cbgeom(:, 2)

        ! Duct
        ibfrst(2) = iblast(1) + 1
        iblast(2) = iblast(1) + m__n_duct
        xb(ibfrst(2):iblast(2)) = m__ductgeom(:, 1)
        yb(ibfrst(2):iblast(2)) = m__ductgeom(:, 2)

        nbel = 2
        nbtype(:2) = 0
        nbtot = iblast(2)

        aname = ''

        ! Initialization
        call geproc
        lload = .true.
        do i_rotor = 1, m__n_rotors
            lbbloft(i_rotor) = .false.
        enddo
    end subroutine set_case

    subroutine get_case(fname)
        use m_dfdcsubs, only : dfsave
        ! bind(c, name = 'get_case')
        implicit none
        character(*) :: fname
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        call dfsave(fname)
    end subroutine get_case

    subroutine oper()
        use i_dfdc
        use m_inigrd, only : setgrdflw
        use m_rotoper, only : rotinitbgam, tqcalc, rotrprt, rotinitbld
        use m_wakesubs, only : wakereset
        use s_rotoper, only : convgthbg
        use m_dfdcsubs, only : gengeom
        ! bind(c, name = 'oper')
        implicit none
        real :: rlx, rlxf, wxeps
        integer :: ip, ir, itrmax, ityp, n
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        ! generate paneled geometry for case
        call gengeom

        if (nptot==0) then
            write (*, *)
            write (*, *) '***  No paneling available  ***'
            return
        endif

        ! clear imposed viscous strengths
        do ip = 1, nptot
            gamvsp(ip) = 0.
            sigvsp(ip) = 0.
        enddo

        ! Calculate solution for current actuator or blade using new solver
        rlx = rlxsolv
        wxeps = epssolv
        itrmax = itrmaxsolv

        rlxf = rlx
        if (nrotor>0) then
            n = 1
            if (irtype(n)==1) then
                write (*, *) 'Using actuator disk...'
                if (.not.lvmav) then
                    call rotinitbgam
                    ! Put flow data into wake grid
                    call setgrdflw
                endif
                call convgthbg(itrmax, rlxf, wxeps)
            elseif (irtype(n)==2) then
                write (*, *) 'Using current blade...'
                ! check for uninitialized BGAM for blade
                do ir = 1, nrc
                    if (bgam(ir, n)/=0.0) goto 32
                enddo
                call rotinitbld
                ! Put flow data into wake grid
                call setgrdflw
                32         call convgthbg(itrmax, rlxf, wxeps)
            endif

            rlxf = rlx
            if (lwrlx) then
                call wakereset
                if (irtype(n)==1) then
                    call convgthbg(itrmax, rlxf, wxeps)
                else
                    call convgthbg(itrmax, rlxf, wxeps)
                endif
            endif

            ityp = 1
            call tqcalc(ityp)
            call rotrprt(6)
        endif
    end subroutine oper

end module api
