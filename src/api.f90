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

    subroutine set_case(fname)
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
