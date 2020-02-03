!*==API.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
! AWRITE





MODULE api
    USE, intrinsic :: iso_c_binding, only : c_float, c_double, &
            & c_int, c_bool, c_char
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    PRIVATE :: KIND
    !
    ! PARAMETER definitions
    !
    INTEGER, PRIVATE, PARAMETER :: DP = kind(0.D0)
    !
    ! Local variables
    !
    LOGICAL, PRIVATE :: FERROR, LDEBUG, LMODG
    REAL, PUBLIC :: GET_CASE, OPER, SET_CASE
    INTEGER, PUBLIC :: INIT
    !
    !*** End of declarations rewritten by SPAG
    !
    !    implicit none


CONTAINS

    SUBROUTINE init()
        ! bind(c, name = 'init')
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        LDEBUG = .FALSE.
        CALL DFINIT(LDEBUG)
        !        CALL LOFTINIT1
        !        CALL LOFTINIT2(1,NRX)
    END SUBROUTINE INIT

    SUBROUTINE set_case(fname)
        ! bind(c, name = 'set_case')
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        CALL DFLOAD(FNAME, FERROR)
        !            CALL LOFTINIT2(1,NROTOR)
        IF (.NOT.FERROR) LONAME = NAME
    END SUBROUTINE SET_CASE

    SUBROUTINE get_case(fname)
        ! bind(c, name = 'get_case')
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        CALL DFSAVE(fname)
    END SUBROUTINE GET_CASE

    SUBROUTINE oper()
        ! bind(c, name = 'oper')
        USE I_DFDC
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        ! generate paneled geometry for case
        CALL GENGEOM

        IF (NPTOT==0) THEN
            WRITE (*, *)
            WRITE (*, *) '***  No paneling available  ***'
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
        RLX = RLXSOLV
        WXEPS = EPSSOLV
        ITRMAX = ITRMAXSOLV

        RLXF = RLX
        IF (NROTOR>0) THEN
            N = 1
            IF (IRTYPE(N)==1) THEN
                WRITE (*, *) 'Using actuator disk...'
                IF (.NOT.LVMAV) THEN
                    CALL ROTINITBGAM
                    ! Put flow data into wake grid
                    CALL SETGRDFLW
                ENDIF
                CALL CONVGTHBG(ITRMAX, RLXF, WXEPS)
            ELSEIF (IRTYPE(N)==2) THEN
                WRITE (*, *) 'Using current blade...'
                ! check for uninitialized BGAM for blade
                DO IR = 1, NRC
                    IF (BGAM(IR, N)/=0.0) GOTO 32
                ENDDO
                CALL ROTINITBLD
                ! Put flow data into wake grid
                CALL SETGRDFLW
                32         CALL CONVGTHBG(ITRMAX, RLXF, WXEPS)
            ENDIF

            RLXF = RLX
            IF (LWRLX) THEN
                CALL WAKERESET
                IF (IRTYPE(N)==1) THEN
                    CALL CONVGTHBG(ITRMAX, RLXF, WXEPS)
                ELSE
                    CALL CONVGTHBG(ITRMAX, RLXF, WXEPS)
                ENDIF
            ENDIF

            ITYP = 1
            CALL TQCALC(ITYP)
            CALL ROTRPRT(6)
        ENDIF
    END SUBROUTINE OPER

END MODULE API
