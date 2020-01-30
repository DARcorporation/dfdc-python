      SUBROUTINE GETVEL(FNAME1)
C----------------------------------------------------------
C     Reads in incoming velocity profile to be superimposed
C     on freestream and rotational components.
C----------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DIMENSION W1(IRX), W2(IRX), W3(IRX)
      CHARACTER*128 FNAME
      CHARACTER*(*) FNAME1
C
      DATA UWT, VWT / 1.0, 1.0 /
C
      LUTEMP = 3
      LU = LUTEMP
C
      FNAME = FNAME1
      IF(FNAME(1:1) .EQ. ' ') CALL ASKS('Enter input filename^',FNAME)
C
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=95)
      DO 10 I=1, IRX
        READ(LU,*,END=11,ERR=95) W1(I), W2(I), W3(I)
 10   CONTINUE
      I = IRX+1
C
 11   IF(I.LT.3) GO TO 96
      IF(NINFL.GE.2) THEN
       WRITE(*,*)
       WRITE(*,*) '*** Current slipstream profiles overwritten'
      ENDIF
      NINFL = I-1
      CLOSE(LU)
C
      DO 20 I=1, NINFL
        RINFL(I)  = W1(I)
        VAINFL(I) = W2(I)
        VTINFL(I) = W3(I)
 20   CONTINUE
C
      CALL ASKR('Enter axial velocity weight (0 ===> 1)^',UWT)
      CALL ASKR('Enter tang. velocity weight (0 , +/-1)^',VWT)
      DO 30 I=1, NINFL
        VAINFL(I) = UWT*VAINFL(I)
        VTINFL(I) = VWT*VTINFL(I)
 30   CONTINUE
C
      CALL SPLINE(VAINFL,VAINFLR,RINFL,NINFL)
      CALL SPLINE(VTINFL,VTINFLR,RINFL,NINFL)
C
      WRITE(*,1200)
      WRITE(*,1250) (RINFL(I), VAINFL(I), VTINFL(I), I=1, NINFL)
      RETURN
C
 95   WRITE(*,*) 'File read error'
 96   WRITE(*,*) 'New wake velocities not read'
      
C....................................................
 1200 FORMAT(/' External slipstream velocity profiles:'
     &       /'      r (m)     Vaxi (m/s)  Vrot (m/s)')
CCC                 0.12341     0.12324    -0.08922
 1250 FORMAT(1X, 3F12.5)
      END
      

      SUBROUTINE SAVVEL(FNAME1)
      INCLUDE 'DFDC.INC'
      CHARACTER*128 FNAME
      CHARACTER*(*) FNAME1
C--------------------------------------------------
C     Writes out induced velocities for current
C     rotor and operating state.
C--------------------------------------------------
      CHARACTER*1 ANS
C
      LUTEMP = 3
      LU = LUTEMP
C
      IF(.NOT.LCONV) THEN
       WRITE(*,*) 'Must define operating point first.'
       RETURN
      ENDIF
C
      FNAME = FNAME1
      IF(FNAME(1:1) .EQ. ' ') CALL ASKS('Enter output filename^',FNAME)
C     
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=5)
      WRITE(*,*)
      WRITE(*,*) 'Output file exists.  Overwrite?  Y'
      READ (*,1000) ANS
 1000 FORMAT(A)
      IF(INDEX('Nn',ANS).EQ.0) GO TO 6
C
      CLOSE(LU)
      WRITE(*,*) 'Velocities not saved.'
      RETURN
C
 5    OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=90)
 6    REWIND(LU)
C
      NR = 1
      DO 10 I=1, NRC
        RDIM = YRC(I,NR)
        UDIM = VIND(1,I,NR)
        VDIM = VIND(3,I,NR)
        WRITE(LU,*) RDIM, UDIM, VDIM
 10   CONTINUE
      CLOSE(LU)
      RETURN
C
 90   WRITE(*,*) 'Bad filename.'
      WRITE(*,*) 'Velocities not saved.'
      RETURN
      END


      SUBROUTINE UVINFL(RR,WWA,WWT)
      INCLUDE 'DFDC.INC'
C
      WWA = 0.0
      WWT = 0.0
      IF(NINFL.LE.1) RETURN
C
      RDIM = RR
      IF(RDIM.GE.RINFL(1) .AND. RDIM.LE.RINFL(NINFL)) THEN
       WWA = SEVAL(RDIM,VAINFL,VAINFLR,RINFL,NINFL)
       WWT = SEVAL(RDIM,VTINFL,VTINFLR,RINFL,NINFL)
      ENDIF
C
      RETURN
      END


