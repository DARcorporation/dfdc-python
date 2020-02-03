PROGRAM dfdc_test
    USE api
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Local variables
    !
    REAL :: API
    CHARACTER(128) :: FNAME
    !
    !*** End of declarations rewritten by SPAG
    !
    CALL init
    CALL ASKS(' Enter filename to load)^', FNAME)
    CALL set_case(fname)
    CALL oper
    CALL ASKS(' Enter filename to save)^', FNAME)
    CALL get_case(fname)
END PROGRAM DFDC_TEST
