MODULE api
    USE, intrinsic :: iso_c_binding, only : c_float, c_double, &
            & c_int, c_bool, c_char
    IMPLICIT NONE
CONTAINS

    SUBROUTINE INIT()
        use m_dfdcsubs, only: dfinit
        ! bind(c, name = 'init')
        IMPLICIT NONE
        LOGICAL :: LDEBUG
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

    SUBROUTINE SET_CASE(fname)
        use i_dfdc
        use m_dfdcsubs, only: dfload
        ! bind(c, name = 'set_case')
        IMPLICIT NONE
        CHARACTER(*) :: FNAME
        LOGICAL :: FERROR
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        CALL DFLOAD(FNAME, FERROR)
        !            CALL LOFTINIT2(1,NROTOR)
        IF (.NOT.FERROR) LONAME = NAME
    END SUBROUTINE SET_CASE

    SUBROUTINE GET_CASE(fname)
        use m_dfdcsubs, only: dfsave
        ! bind(c, name = 'get_case')
        IMPLICIT NONE
        CHARACTER(*) :: FNAME
        !
        !*** Start of declarations rewritten by SPAG
        !
        !*** End of declarations rewritten by SPAG
        !
        CALL DFSAVE(fname)
    END SUBROUTINE GET_CASE

    SUBROUTINE OPER()
        use i_dfdc
        use m_inigrd, only: setgrdflw
        use m_rotoper, only: rotinitbgam, tqcalc, rotrprt, rotinitbld
        use m_wakesubs, only: wakereset
        use s_rotoper, only: convgthbg
        use m_dfdcsubs, only: gengeom
        ! bind(c, name = 'oper')
        IMPLICIT NONE
        REAL :: RLX, RLXF, WXEPS
        INTEGER :: IP, IR, ITRMAX, ITYP, N
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
