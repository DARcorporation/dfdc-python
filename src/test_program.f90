PROGRAM dfdc_test
    USE api
    USE m_userio, only : asks
    IMPLICIT NONE
    CHARACTER(128) :: FNAME
    CALL init
    CALL ASKS(' Enter filename to load)^', FNAME)
    CALL set_case(fname)
    CALL oper
    CALL ASKS(' Enter filename to save)^', FNAME)
    CALL get_case(fname)
END PROGRAM DFDC_TEST
