module m_viscvel
    implicit none
contains
    !*==GETVEL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    
    SUBROUTINE GETVEL(FNAME1)
        use m_spline, only: spline
        use m_userio, only: asks, askr
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        CHARACTER(*) :: FNAME1
        !
        ! Local variables
        !
        CHARACTER(128) :: FNAME
        INTEGER :: I, LU, LUTEMP
        REAL, SAVE :: UWT, VWT
        REAL, DIMENSION(IRX) :: W1, W2, W3
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------
        !     Reads in incoming velocity profile to be superimposed
        !     on freestream and rotational components.
        !----------------------------------------------------------
        !
        !
        DATA UWT, VWT/1.0, 1.0/
        !
        LUTEMP = 3
        LU = LUTEMP
        !
        FNAME = FNAME1
        IF (FNAME(1:1)==' ') CALL ASKS('Enter input filename^', FNAME)
        !
        OPEN (LU, FILE = FNAME, STATUS = 'OLD', ERR = 95)
        DO I = 1, IRX
            READ (LU, *, END = 11, ERR = 95) W1(I), W2(I), W3(I)
        ENDDO
        I = IRX + 1
        !
        11   IF (I<3) GOTO 96
        IF (NINFL>=2) THEN
            WRITE (*, *)
            WRITE (*, *) '*** Current slipstream profiles overwritten'
        ENDIF
        NINFL = I - 1
        CLOSE (LU)
        !
        DO I = 1, NINFL
            RINFL(I) = W1(I)
            VAINFL(I) = W2(I)
            VTINFL(I) = W3(I)
        ENDDO
        !
        CALL ASKR('Enter axial velocity weight (0 ===> 1)^', UWT)
        CALL ASKR('Enter tang. velocity weight (0 , +/-1)^', VWT)
        DO I = 1, NINFL
            VAINFL(I) = UWT * VAINFL(I)
            VTINFL(I) = VWT * VTINFL(I)
        ENDDO
        !
        CALL SPLINE(VAINFL, VAINFLR, RINFL, NINFL)
        CALL SPLINE(VTINFL, VTINFLR, RINFL, NINFL)
        !
        WRITE (*, 1200)
        WRITE (*, 1250) (RINFL(I), VAINFL(I), VTINFL(I), I = 1, NINFL)
        RETURN
        !
        95   WRITE (*, *) 'File read error'
        96   WRITE (*, *) 'New wake velocities not read'
    
        !....................................................
        1200 FORMAT (/' External slipstream velocity profiles:'/               &
                &'      r (m)     Vaxi (m/s)  Vrot (m/s)')
        !CC                 0.12341     0.12324    -0.08922
        1250 FORMAT (1X, 3F12.5)
    END SUBROUTINE GETVEL
    !*==SAVVEL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE SAVVEL(FNAME1)
        use m_userio, only: asks
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        CHARACTER(*) :: FNAME1
        !
        ! Local variables
        !
        CHARACTER(1) :: ANS
        CHARACTER(128) :: FNAME
        INTEGER :: I, LU, LUTEMP, NR
        REAL :: RDIM, UDIM, VDIM
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------
        !     Writes out induced velocities for current
        !     rotor and operating state.
        !--------------------------------------------------
        !
        LUTEMP = 3
        LU = LUTEMP
        !
        IF (.NOT.LCONV) THEN
            WRITE (*, *) 'Must define operating point first.'
            RETURN
        ENDIF
        !
        FNAME = FNAME1
        IF (FNAME(1:1)==' ') CALL ASKS('Enter output filename^', FNAME)
        !
        OPEN (LU, FILE = FNAME, STATUS = 'OLD', ERR = 5)
        WRITE (*, *)
        WRITE (*, *) 'Output file exists.  Overwrite?  Y'
        READ (*, 1000) ANS
        1000 FORMAT (A)
        IF (INDEX('Nn', ANS)==0) GOTO 6
        !
        CLOSE (LU)
        WRITE (*, *) 'Velocities not saved.'
        RETURN
        !
        5    OPEN (LU, FILE = FNAME, STATUS = 'NEW', ERR = 90)
        6    REWIND (LU)
        !
        NR = 1
        DO I = 1, NRC
            RDIM = YRC(I, NR)
            UDIM = VIND(1, I, NR)
            VDIM = VIND(3, I, NR)
            WRITE (LU, *) RDIM, UDIM, VDIM
        ENDDO
        CLOSE (LU)
        RETURN
        !
        90   WRITE (*, *) 'Bad filename.'
        WRITE (*, *) 'Velocities not saved.'
    END SUBROUTINE SAVVEL
    !*==UVINFL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE UVINFL(RR, WWA, WWT)
        use m_spline, only: seval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: RR, WWA, WWT
        !
        ! Local variables
        !
        REAL :: RDIM
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        WWA = 0.0
        WWT = 0.0
        IF (NINFL<=1) RETURN
        !
        RDIM = RR
        IF (RDIM>=RINFL(1) .AND. RDIM<=RINFL(NINFL)) THEN
            WWA = SEVAL(RDIM, VAINFL, VAINFLR, RINFL, NINFL)
            WWT = SEVAL(RDIM, VTINFL, VTINFLR, RINFL, NINFL)
        ENDIF
        !
    END SUBROUTINE UVINFL
end module m_viscvel
