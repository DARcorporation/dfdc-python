module api
    use, intrinsic :: iso_c_binding, only : c_float, c_double, c_int, c_bool, c_char
!    implicit none
    private
    public init, set_case, get_case, oper

    integer, parameter :: dp = kind(0.D0)

    LOGICAL LMODG, LDEBUG, FERROR
contains

    subroutine init() ! bind(c, name = 'init')
        INCLUDE 'DFDC.INC'
        LDEBUG = .FALSE.
        CALL DFINIT(LDEBUG)
!        CALL LOFTINIT1
!        CALL LOFTINIT2(1,NRX)
    end subroutine init

    subroutine set_case(fname) ! bind(c, name = 'set_case')
        INCLUDE 'DFDC.INC'
        CHARACTER*128 FNAME
        CALL DFLOAD(FNAME,FERROR)
        IF(.NOT.FERROR) THEN
            LONAME = NAME
!            CALL LOFTINIT2(1,NROTOR)
        ENDIF
    end subroutine set_case

    subroutine get_case(fname) ! bind(c, name = 'get_case')
        INCLUDE 'DFDC.INC'
        CHARACTER*128 FNAME
        CALL DFSAVE(fname)
    end subroutine get_case

    subroutine oper() ! bind(c, name = 'oper')
        INCLUDE 'DFDC.INC'
        ! generate paneled geometry for case
        CALL GENGEOM

        IF(NPTOT == 0) THEN
            WRITE(*,*)
            WRITE(*,*) '***  No paneling available  ***'
            RETURN
        ENDIF

        ! clear imposed viscous strengths
        DO IP = 1, NPTOT
            GAMVSP(IP) = 0.
            SIGVSP(IP) = 0.
        ENDDO

        ! plot Cp(x) or Q(x)
        ICPQXY = 1
        ! drive total thrust (ISPEC=2), drive rotor thrust for ISPEC=1
        ISPEC = 2
        ! default working disk
        NDSK = 1

        ! Calculate solution for current actuator or blade using new solver
        RLX    = RLXSOLV
        WXEPS  = EPSSOLV
        ITRMAX = ITRMAXSOLV

        RLXF = RLX
        IF(NROTOR > 0) THEN
            N = 1
            IF(IRTYPE(N) == 1) THEN
                WRITE(*,*) 'Using actuator disk...'
                IF(.NOT.LVMAV) THEN
                    CALL ROTINITBGAM
                    ! Put flow data into wake grid
                    CALL SETGRDFLW
                ENDIF
                CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
            ELSEIF(IRTYPE(N) == 2) THEN
                WRITE(*,*) 'Using current blade...'
                ! check for uninitialized BGAM for blade
                DO IR = 1, NRC
                    IF(BGAM(IR,N) /= 0.0) GO TO 32
                END DO
                CALL ROTINITBLD
                ! Put flow data into wake grid
                CALL SETGRDFLW
                32        CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
            ENDIF

            RLXF = RLX
            IF(LWRLX) THEN
                CALL WAKERESET
                IF(IRTYPE(N) == 1) THEN
                    CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
                ELSE
                    CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
                ENDIF
            ENDIF

            ITYP = 1
            CALL TQCALC(ITYP)
            CALL ROTRPRT(6)
        ENDIF
    end subroutine oper

end module api