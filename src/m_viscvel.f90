module m_viscvel
    implicit none
contains
    !*==GETVEL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020



    subroutine getvel(fname1)
        use i_dfdc
        use m_spline, only : spline
        use m_userio, only : asks, askr
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: fname1
        !
        ! Local variables
        !
        character(128) :: fname
        integer :: i, lu, lutemp
        real, save :: uwt, vwt
        real, dimension(irx) :: w1, w2, w3
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------
        !     Reads in incoming velocity profile to be superimposed
        !     on freestream and rotational components.
        !----------------------------------------------------------
        !
        !
        data uwt, vwt/1.0, 1.0/
        !
        lutemp = 3
        lu = lutemp
        !
        fname = fname1
        if (fname(1:1)==' ') call asks('Enter input filename^', fname)
        !
        open (lu, file = fname, status = 'OLD', err = 95)
        do i = 1, irx
            read (lu, *, end = 11, err = 95) w1(i), w2(i), w3(i)
        enddo
        i = irx + 1
        !
        11   if (i<3) goto 96
        if (ninfl>=2) then
            write (*, *)
            write (*, *) '*** Current slipstream profiles overwritten'
        endif
        ninfl = i - 1
        close (lu)
        !
        do i = 1, ninfl
            rinfl(i) = w1(i)
            vainfl(i) = w2(i)
            vtinfl(i) = w3(i)
        enddo
        !
        call askr('Enter axial velocity weight (0 ===> 1)^', uwt)
        call askr('Enter tang. velocity weight (0 , +/-1)^', vwt)
        do i = 1, ninfl
            vainfl(i) = uwt * vainfl(i)
            vtinfl(i) = vwt * vtinfl(i)
        enddo
        !
        call spline(vainfl, vainflr, rinfl, ninfl)
        call spline(vtinfl, vtinflr, rinfl, ninfl)
        !
        write (*, 1200)
        write (*, 1250) (rinfl(i), vainfl(i), vtinfl(i), i = 1, ninfl)
        return
        !
        95   write (*, *) 'File read error'
        96   write (*, *) 'New wake velocities not read'

        !....................................................
        1200 format (/' External slipstream velocity profiles:'/               &
                &'      r (m)     Vaxi (m/s)  Vrot (m/s)')
        !CC                 0.12341     0.12324    -0.08922
        1250 format (1x, 3f12.5)
    end subroutine getvel
    !*==SAVVEL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine savvel(fname1)
        use i_dfdc
        use m_userio, only : asks
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        character(*) :: fname1
        !
        ! Local variables
        !
        character(1) :: ans
        character(128) :: fname
        integer :: i, lu, lutemp, nr
        real :: rdim, udim, vdim
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Writes out induced velocities for current
        !     rotor and operating state.
        !--------------------------------------------------
        !
        lutemp = 3
        lu = lutemp
        !
        if (.not.lconv) then
            write (*, *) 'Must define operating point first.'
            return
        endif
        !
        fname = fname1
        if (fname(1:1)==' ') call asks('Enter output filename^', fname)
        !
        open (lu, file = fname, status = 'OLD', err = 5)
        write (*, *)
        write (*, *) 'Output file exists.  Overwrite?  Y'
        read (*, 1000) ans
        1000 format (a)
        if (index('Nn', ans)==0) goto 6
        !
        close (lu)
        write (*, *) 'Velocities not saved.'
        return
        !
        5    open (lu, file = fname, status = 'NEW', err = 90)
        6    rewind (lu)
        !
        nr = 1
        do i = 1, nrc
            rdim = yrc(i, nr)
            udim = vind(1, i, nr)
            vdim = vind(3, i, nr)
            write (lu, *) rdim, udim, vdim
        enddo
        close (lu)
        return
        !
        90   write (*, *) 'Bad filename.'
        write (*, *) 'Velocities not saved.'
    end subroutine savvel
    !*==UVINFL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine uvinfl(rr, wwa, wwt)
        use i_dfdc
        use m_spline, only : seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: rr, wwa, wwt
        !
        ! Local variables
        !
        real :: rdim
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        wwa = 0.0
        wwt = 0.0
        if (ninfl<=1) return
        !
        rdim = rr
        if (rdim>=rinfl(1) .and. rdim<=rinfl(ninfl)) then
            wwa = seval(rdim, vainfl, vainflr, rinfl, ninfl)
            wwt = seval(rdim, vtinfl, vtinflr, rinfl, ninfl)
        endif
        !
    end subroutine uvinfl
end module m_viscvel
