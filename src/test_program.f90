program dfdc_test
    use api
    use m_userio, only : asks
    implicit none
    character(128) :: fname
    call init
    call asks(' Enter filename to load)^', fname)
    call set_case(fname)
    call oper
    call asks(' Enter filename to save)^', fname)
    call get_case(fname)
end program dfdc_test
