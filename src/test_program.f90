module m_test_program
    implicit none
contains
    !*==DFDC_TEST.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SWDOTS
    
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
end module m_test_program
