program dfdc_test
    use api
    CHARACTER*128 FNAME
    call init
    CALL ASKS(' Enter filename to load)^',FNAME)
    call set_case(fname)
    call oper
    CALL ASKS(' Enter filename to save)^',FNAME)
    call get_case(fname)
end program dfdc_test