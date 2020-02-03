module s_rotoper
    implicit none
contains


    SUBROUTINE CONVGTH(NITER, RLXF, WXEPS)
        use i_dfdc
        use m_rotoper, only : vmavgcalc, vmavginit, updrotvel, gthcalc
        use m_solve, only : gamsolv
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: NITER
        REAL :: RLXF, WXEPS
        !
        ! Local variables
        !
        REAL :: DG, DGTHMAX, RLX, RLXG
        REAL, DIMENSION(IPX) :: DGOLD, GAMTH
        INTEGER :: IP, IPMAX, ITR
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Basic solver for GTH,
        !     Uses underrelaxed iteration for fixed BGAM to converge GTH
        !----------------------------------------------------------------
        !
        !
        RLX = RLXF
        IF (RLX<=0.0) RLX = 0.5
        !
        !---- check for valid solution before iterating
        IF (.NOT.LGAMA) THEN
            IF (.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
            !---- Generate GTH solution for current wakes
            CALL GTHCALC(GAMTH)
            !---- Update wake gamma from initial solution
            DO IP = 1, NPTOT
                GTH(IP) = GAMTH(IP)
            ENDDO
            !---- Generate GAM solution for current RHS
            CALL GAMSOLV
        ENDIF
        !
        DO IP = 1, NPTOT
            DGOLD(IP) = 0.0
        ENDDO
        !
        !----- do specified cycles of under-relaxed iteration
        DO ITR = 1, NITER
            !c         IF(LDBG) WRITE(*,110) ITR
            !
            !---- Set VMavg velocities at wake points
            CALL VMAVGCALC
            !---- Generate GTH solution for current wakes
            CALL GTHCALC(GAMTH)
            !---- Update wake gamma using CSOR
            IPMAX = 0
            DGTHMAX = 0.0
            RLX = RLXF
            DO IP = 1, NPTOT
                DG = GAMTH(IP) - GTH(IP)
                IF (ABS(DG)>ABS(DGTHMAX)) THEN
                    DGTHMAX = DG
                    IPMAX = IP
                ENDIF
                RLXG = RLX
                IF (DG * DGOLD(IP)<0.0) RLXG = 0.6 * RLX
                IF (DG * DGOLD(IP)>0.0) RLXG = 1.2 * RLX
                DGOLD(IP) = DG * RLXG
                GTH(IP) = GTH(IP) + RLXG * DG
            ENDDO
            !---- Generate GAM solution for current RHS
            LGAMA = .FALSE.
            CALL GAMSOLV
            !
            !c         IF(LDBG)
            WRITE (*, 100) ITR, DGTHMAX, IPMAX, RLX
            IF (ABS(DGTHMAX)<WXEPS * QREF) THEN
                LCONV = .TRUE.
                GOTO 20
            ENDIF
        ENDDO
        LCONV = .FALSE.
        !
        100  FORMAT (I3, ' dGTHmax=', F9.5, ' @IP=', I5, ' RLX=', F8.5)
        120  FORMAT (1X, 7G11.4)
        !
        !---- Update rotor velocities
        20   CALL UPDROTVEL
        !
    END SUBROUTINE CONVGTH
    !*==CONVGTHT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    SUBROUTINE CONVGTHT(NITER, RLXF, WXEPS, TSPEC, ISPEC)
        use i_dfdc
        use m_rotoper, only : vmavgcalc, updrotvel, tqcalc, vmavginit, gthcalc
        use m_solve, only : gamsolv
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: ISPEC, NITER
        REAL :: RLXF, TSPEC, WXEPS
        !
        ! Local variables
        !
        REAL :: BGAV, BGMAG, BGMAX, BGMIN, BLDS, DA, DBG, &
                & DBGMAX, DG, DGTHMAX, DR, FDBG, RLX, RLXB, RLXG, &
                & THR, TSCL
        REAL, DIMENSION(IPX) :: DGOLD, GAMTH
        INTEGER :: IP, IPMAX, IR, ITR, NR
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Basic solver for B*GAM,GTH for design mode (fixed thrust)
        !     Uses underrelaxed iteration for fixed BGAM to converge GTH
        !     Input:  TSPEC   thrust specification
        !             ISPEC   1 for rotor thrust spec, 2 for total thrust spec
        !----------------------------------------------------------------
        !
        !
        RLX = RLXF
        IF (RLX<=0.0) RLX = 0.5
        !
        !---- check for valid solution before iterating
        IF (.NOT.LGAMA) THEN
            IF (.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
            !---- Generate GTH solution for current wakes
            CALL GTHCALC(GAMTH)
            !---- Update wake gamma from initial solution
            DO IP = 1, NPTOT
                GTH(IP) = GAMTH(IP)
            ENDDO
            !---- Generate GAM solution for current RHS
            CALL GAMSOLV
        ENDIF
        !
        DO IP = 1, NPTOT
            DGOLD(IP) = 0.0
        ENDDO
        !
        !----- Iteration loop for driving thrust using B*GAM
        NR = 1
        BLDS = FLOAT(NRBLD(NR))
        DO ITR = 1, NITER
            !c         IF(LDBG) WRITE(*,110) ITR
            !---- Update rotor velocities
            CALL UPDROTVEL
            !
            BGAV = 0.0
            BGMAX = BGAM(1, NR)
            BGMIN = BGMAX
            DO IR = 1, NRC
                DR = YRP(IR + 1, NR) - YRP(IR, NR)
                DA = PI * (YRP(IR + 1, NR)**2 - YRP(IR, NR)**2)
                BGAV = BGAV + BGAM(IR, NR)
                BGMAX = MAX(BGMAX, BGAM(IR, NR))
                BGMIN = MIN(BGMIN, BGAM(IR, NR))
            ENDDO
            BGAV = BGAV / FLOAT(NRC)
            BGMAG = MAX(ABS(BGAV), BGMIN, BGMAX, 0.1)
            !
            !---- Calculate current rotor thrust
            CALL TQCALC(1)
            !---- Drive thrust from rotor for ISPEC=1, total thrust for ISPEC=2
            IF (ISPEC==1) THEN
                THR = TTOT
            ELSE
                THR = TTOT + TDUCT
            ENDIF
            !---- Scale factor for BGAM from current value to get desired rotor thrust
            TSCL = 1.0
            IF (THR/=0.0) TSCL = TSPEC / THR
            !
            !---- Check for rational relaxation factors based on BGAM changes
            DBGMAX = 0.0
            RLXB = RLX
            DO IR = 1, NRC
                DBG = (TSCL - 1.0) * BGAM(IR, NR)
                DBGMAX = MAX(DBGMAX, ABS(DBG))
                IF (BGMAG/=0.0) THEN
                    FDBG = ABS(DBG) / BGMAG
                    IF (FDBG * RLXB>0.3) RLXB = 0.3 / FDBG
                ENDIF
            ENDDO
            !
            DO IR = 1, NRC
                DBG = (TSCL - 1.0) * BGAM(IR, NR)
                !---- Update BGAM
                BGAM(IR, NR) = BGAM(IR, NR) + RLXB * DBG
            ENDDO
            IF (TGAP>0.0) BGAM(NRC, NR) = 0.0
            !
            !---- Set VMavg velocities at wake points
            CALL VMAVGCALC
            !---- Generate GTH estimate for updated circulations
            CALL GTHCALC(GAMTH)
            IPMAX = 0
            DGTHMAX = 0.0
            DO IP = 1, NPTOT
                DG = GAMTH(IP) - GTH(IP)
                IF (ABS(DG)>ABS(DGTHMAX)) THEN
                    DGTHMAX = DG
                    IPMAX = IP
                ENDIF
                RLXG = RLXB
                IF (DG * DGOLD(IP)<0.0) RLXG = 0.6 * RLXB
                IF (DG * DGOLD(IP)>0.0) RLXG = 1.2 * RLXB
                DGOLD(IP) = DG * RLXG
                !---- Update GTH wake gamma using CSOR
                GTH(IP) = GTH(IP) + RLXG * DG
            ENDDO
            !
            !---- Generate GAM solution for current RHS
            LGAMA = .FALSE.
            CALL GAMSOLV
            !
            !c       IF(LDBG) THEN
            !c         WRITE(*,*) ' '
            WRITE (*, 100) ITR, DGTHMAX, IPMAX, DBGMAX, RLXB
            IF (ABS(DGTHMAX)<WXEPS * QREF .AND. ABS(DBGMAX)<=0.001 * BGMAG)  &
                    & THEN
                LCONV = .TRUE.
                GOTO 20
            ENDIF
            LGAMA = .FALSE.
        ENDDO
        LCONV = .FALSE.
        !
        100  FORMAT (I3, ' dGTHmax=', F9.5, ' @IP=', I5, '  dBGmax=', F9.5, ' RLX=', &
                & F8.5)
        110  FORMAT (/'Blade velocities on iteration ', I5, &
                &/'     r          Wx         Wr', &
                &'         Wt        Phi       CL       BGam')
        120  FORMAT (1X, 7G11.4)
        !
        !---- Update rotor velocities
        20   CALL UPDROTVEL
        !
    END SUBROUTINE CONVGTHT
    !*==CONVGTHBG.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    SUBROUTINE CONVGTHBG(NITER, RLXF, WXEPS)
        use i_dfdc
        use m_inigrd, only : setgrdflw
        use m_solve, only : gamsolv
        use m_rotoper, only : vmavgcalc, vmavginit, updrotvel, gthcalc
        use m_aero, only : getclcdcm
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: NITER
        REAL :: RLXF, WXEPS
        !
        ! Local variables
        !
        REAL :: ALF, BGAV, BGMAG, BGMAX, BGMIN, BLDS, CD_ALF, &
                & CD_REY, CD_W, CLB, CLMAX, CLMIN, CL_ALF, CL_W, &
                & CMOM, CM_AL, CM_W, DBG, DBGMAX, DCL_STALL, DG, &
                & DGTHMAX, FDBG, PHIB, REY, RLX, RLXB, RLXBG, &
                & RLXG, SECSIG, SECSTAGR, VTBG, WTB, WWB, XI
        REAL, DIMENSION(IRX) :: BGX, DBGOLD
        REAL, DIMENSION(IPX) :: DGOLD, GAMTH
        INTEGER :: IG, IP, IPMAX, IR, IRMAX, ITR, NR
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Basic solver for GTH for a defined blade geometry
        !     Uses underrelaxed iteration to converge GTH, BGAM
        !----------------------------------------------------------------
        !
        !
        RLX = RLXF
        IF (RLX<=0.0) RLX = 0.5
        !
        !---- check for valid solution before iterating
        IF (.NOT.LGAMA) THEN
            IF (.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
            !---- Generate GTH solution for current wakes
            CALL GTHCALC(GAMTH)
            !---- Update wake gamma from initial solution
            DO IP = 1, NPTOT
                GTH(IP) = GAMTH(IP)
            ENDDO
            !---- Generate GAM solution for current RHS
            CALL GAMSOLV
        ENDIF
        !
        DO IP = 1, NPTOT
            DGOLD(IP) = 0.0
        ENDDO
        DO IR = 1, NRC
            DBGOLD(IR) = 0.0
        ENDDO
        !
        !c      write(15,*) 'new solution '
        !----- do several cycles of under-relaxed iteration to converge
        !      BGAM and GTH from specified CL and chord
        DO ITR = 1, NITER
            !c           write(15,*) 'iter ',itr
            !c         IF(LDBG) WRITE(*,110) ITR
            !
            !---- Generate GAM solution for current RHS
            CALL GAMSOLV
            !---- Update rotor velocities
            CALL UPDROTVEL
            !
            RLXG = 0.0
            RLXBG = 0.0
            DBGMAX = 0.0
            DGTHMAX = 0.0
            !
            DO NR = 1, NROTOR
                IF (IRTYPE(NR)==2) THEN
                    IG = IGROTOR(NR)
                    !
                    BLDS = FLOAT(NRBLD(NR))
                    !---- convert to blade relative velocities
                    BGAV = 0.0
                    BGMAX = BGAM(1, NR)
                    BGMIN = BGMAX
                    DO IR = 1, NRC
                        !---- theta velocity at blade lifting line
                        VTBG = BGAM(IR, NR) * PI2I / YRC(IR, NR)
                        WTB = VREL(3, IR, NR) - 0.5 * VTBG
                        !
                        !           IF(NR.EQ.1) THEN
                        !             CIRC = 0.0
                        !           ELSE
                        !             CIRC = BGAMG(IG-1,IR)
                        !           ENDIF
                        !           VTBG = BGAM(IR,NR)*PI2I/YRC(IR,NR)
                        !           VTIN = CIRC*PI2I/YRC(IR,NR)
                        !           VROT = - OMEGA(NR)*YRC(IR,NR)
                        !           WTB = VTIN - OMEGA(NR)*YRC(IR,NR) + 0.5*VTBG
                        !
                        WWB = SQRT(VREL(1, IR, NR)**2 + WTB**2)
                        PHIB = ATAN2(VREL(1, IR, NR), -WTB)
                        !
                        !           write(15,89) 'n,i,vin,vbg,vout,wtb,phi ',nr,ir,VTIN,
                        !     &                     VTBG,VREL(3,IR,NR),WTB,phib/dtr,vrot
                        ! 89        format(a,2i4,6F12.5)
                        !
                        XI = YRC(IR, NR) / RTIP(NR)
                        ALF = BETAR(IR, NR) - PHIB
                        REY = WWB * CHR(IR, NR) * RHO / RMU
                        SECSIG = BLDS * CHR(IR, NR) / (2.0 * PI * YRC(IR, NR))
                        SECSTAGR = 0.5 * PI - BETAR(IR, NR)
                        CALL GETCLCDCM(NR, IR, XI, ALF, WWB, REY, SECSIG, SECSTAGR, &
                                & CLB, CL_ALF, CL_W, CLMAX, CLMIN, DCL_STALL, &
                                & LSTALLR(IR, NR), CDR(IR, NR), CD_ALF, CD_W, &
                                & CD_REY, CMOM, CM_AL, CM_W)
                        CLR(IR, NR) = CLB
                        CLALF(IR, NR) = CL_ALF
                        ALFAR(IR, NR) = ALF
                        !
                        BGX(IR) = BLDS * 0.5 * WWB * CHR(IR, NR) * CLR(IR, NR)
                        BGAV = BGAV + BGAM(IR, NR)
                        BGMAX = MAX(BGMAX, BGAM(IR, NR))
                        BGMIN = MIN(BGMIN, BGAM(IR, NR))
                        !
                        !            WRITE(*,120) YRC(IR,NR),
                        !     &                   VREL(1,IR,NR),VREL(2,IR,NR),WTB,
                        !     &                   PHIB/DTR,CLR(IR,NR),BGX(IR),CLALF(IR,NR)
                        IF (LDBG) WRITE (*, 120) YRC(IR, NR), VREL(1, IR, NR), &
                                & VREL(2, IR, NR), WTB, PHIB / DTR, &
                                & CLR(IR, NR), BGX(IR), CLALF(IR, NR)
                    ENDDO
                    BGAV = BGAV / FLOAT(NRC)
                    IF (BGAV>=0.0) THEN
                        BGMAG = MAX(BGAV, BGMAX, 0.1)
                    ELSE
                        BGMAG = MIN(BGAV, BGMIN, -0.1)
                    ENDIF
                    IF (TGAP>0.0 .AND. OMEGA(NR)/=0.0) THEN
                        BGX(NRC) = 0.0
                        CLR(NRC, NR) = 0.0
                        ALFAR(NRC, NR) = 0.0
                    ENDIF
                    !
                    !---- Check for rational relaxation factors based on BGAM changes
                    RLXB = RLX
                    IRMAX = 0
                    DBGMAX = 0.0
                    DO IR = 1, NRC
                        DBG = (BGX(IR) - BGAM(IR, NR))
                        IF (ABS(DBG)>ABS(DBGMAX)) THEN
                            DBGMAX = DBG
                            IRMAX = IR
                        ENDIF
                        IF (BGMAG/=0.0) THEN
                            !             FDBG = ABS(DBG)/BGAM(IR,NR)
                            !             IF(FDBG*RLXB.GT.0.5) RLXB = 0.5/FDBG
                            !             FDBG = DBG/BGAM(IR,NR)
                            FDBG = DBG / BGMAG
                            IF (FDBG * RLXB<-0.2) RLXB = -0.2 / FDBG
                            IF (FDBG * RLXB>0.4) RLXB = 0.4 / FDBG
                        ENDIF
                    ENDDO
                    !
                    !---- Update blade circulation using CSOR
                    DO IR = 1, NRC
                        DBG = (BGX(IR) - BGAM(IR, NR))
                        RLXBG = 0.5 * RLXB
                        IF (DBG * DBGOLD(IR)<0.0) RLXBG = 0.6 * RLXB
                        !            IF(DBG*DBGOLD(IR).GT.0.0) RLXBG = 1.2*RLXB
                        DBGOLD(IR) = DBG * RLXBG
                        BGAM(IR, NR) = BGAM(IR, NR) + RLXBG * DBG
                    ENDDO
                    IF (TGAP>0.0 .AND. OMEGA(NR)/=0.0) BGAM(NRC, NR) = 0.0
                    !---- Update grid flowfield
                    CALL SETGRDFLW
                ENDIF
                !
            ENDDO ! loop over NROTOR
            !
            !---- Set VMavg velocities at wake points
            CALL VMAVGCALC
            !---- Generate GTH estimate for updated circulations
            CALL GTHCALC(GAMTH)
            IPMAX = 0
            DGTHMAX = 0.0
            DO IP = 1, NPTOT
                DG = GAMTH(IP) - GTH(IP)
                IF (ABS(DG)>ABS(DGTHMAX)) THEN
                    DGTHMAX = DG
                    IPMAX = IP
                ENDIF
                RLXG = RLX
                IF (DG * DGOLD(IP)<0.0) RLXG = 0.6 * RLX
                IF (DG * DGOLD(IP)>0.0) RLXG = 1.2 * RLX
                DGOLD(IP) = DG * RLXG
                !---- Update GTH wake gamma using CSOR
                GTH(IP) = GTH(IP) + RLXG * DG
            ENDDO
            !
            !---- Generate GAM solution for current RHS
            LGAMA = .FALSE.
            !cc         CALL GAMSOLV
            !
            !c       IF(LDBG) THEN
            !c         WRITE(*,*) ' '
            IF (RLXBG/=0.0) THEN
                WRITE (*, 100) ITR, DGTHMAX, IPMAX, DBGMAX, IRMAX, RLXBG
                !c           WRITE(*,100) ITR,DGTHMAX,IPMAX,DBGMAX,IRMAX,RLXBG,BGMAG
            ELSE
                WRITE (*, 105) ITR, DGTHMAX, IPMAX, RLXG
            ENDIF
            IF (ABS(DGTHMAX)<WXEPS * QREF .AND. ABS(DBGMAX)                 &
                    & <=0.001 * ABS(BGMAG)) THEN
                LCONV = .TRUE.
                GOTO 20
            ENDIF
            LGAMA = .FALSE.
        ENDDO
        LCONV = .FALSE.
        !
        100  FORMAT (I3, ' dGTHmax=', F10.5, ' @IP=', I4, '  dBGmax=', F10.5, ' @IR=', &
                & I4, ' RLX=', F8.5)
        !c     &          '  dBGmax=',F10.5,' @IR=',I4,' RLX=',F8.5,' BGMG=',F8.5)
        105  FORMAT (I3, ' dGTHmax=', F10.5, ' @IP=', I4, ' RLX=', F8.5)
        110  FORMAT (/'Disk velocities on iteration ', I4, &
                &/'     r          Wx         Wr', &
                &'         Wt        Phi       CL       BGam      CLalf')
        120  FORMAT (1X, 8G10.4)
        !
        !---- Update rotor velocities
        20   CALL UPDROTVEL
        !
    END SUBROUTINE CONVGTHBG
    !*==CONVGTHBGT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    SUBROUTINE CONVGTHBGT(NITER, RLXF, WXEPS, TSPEC, ISPEC)
        use i_dfdc
        use m_rotoper, only : tqcalc, updrotvel
        use m_aero, only : getclcdcm
        use m_solve, only : gamsolv
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: ISPEC, NITER
        REAL :: RLXF, TSPEC, WXEPS
        !
        ! Local variables
        !
        REAL :: ALF, BGAV, BGEPS, BGFAC, BGMAG, BGMAX, BGMIN, &
                & BLDS, CD_ALF, CD_REY, CD_W, CLB, CLMAX, CLMIN, &
                & CL_ALF, CL_W, CMOM, CM_AL, CM_W, DA, DBETA, DBG, &
                & DBGMAX, DCL_STALL, DR, DTDALF, FDBG, PHIB, REY, &
                & RLX, RLXB, SECSIG, SECSTAGR, THR, TSCL, VTBG, &
                & WTB, WWB, XI
        INTEGER :: IR, ITR, NITER0, NR
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Basic solver for BETA, GTH for analysis mode (fixed thrust)
        !     Uses underrelaxed iteration for blade pitch to converge BGAM,GTH
        !     Input:  TSPEC   thrust specification
        !             ISPEC   1 for rotor thrust spec, 2 for total thrust spec
        !----------------------------------------------------------------
        !
        !
        RLX = RLXF
        IF (RLX<=0.0) RLX = 0.5
        !
        !---- check for valid solution before iterating
        !---- Generate GAM solution for current RHS
        IF (.NOT.LGAMA) CALL GAMSOLV
        !
        !----- do several cycles of under-relaxed iteration to converge
        !      BGAM and GTH from specified CL and chord
        NR = 1
        BLDS = FLOAT(NRBLD(NR))
        DO ITR = 1, NITER
            !c         IF(LDBG) WRITE(*,110) ITR
            !
            BGEPS = 20.0 * WXEPS
            RLXB = RLXF
            NITER0 = MIN(3, NITER / 2)
            CALL CONVGTHBG(NITER, RLXB, BGEPS)
            !
            CALL TQCALC(1)
            !---- Drive thrust from rotor for ISPEC=1, total thrust for ISPEC=2
            IF (ISPEC==1) THEN
                THR = TTOT
            ELSE
                THR = TTOT + TDUCT
            ENDIF
            !
            !---- Check for disk type (actuator disk or bladed)
            IF (IRTYPE(NR)==1) THEN
                !---- Actuator disk
                BGAV = 0.0
                BGMAX = BGAM(1, NR)
                BGMIN = BGMAX
                DO IR = 1, NRC
                    DR = YRP(IR + 1, NR) - YRP(IR, NR)
                    DA = PI * (YRP(IR + 1, NR)**2 - YRP(IR, NR)**2)
                    BGAV = BGAV + BGAM(IR, NR)
                    BGMAX = MAX(BGMAX, BGAM(IR, NR))
                    BGMIN = MIN(BGMIN, BGAM(IR, NR))
                ENDDO
                BGAV = BGAV / FLOAT(NRC)
                BGMAG = MAX(ABS(BGAV), BGMIN, BGMAX, 0.1)
                !---- Scale factor for BGAM from current value to get desired rotor thrust
                TSCL = 1.0
                IF (THR/=0.0) TSCL = TSPEC / THR
                !---- Check for rational relaxation factors based on BGAM changes
                DBGMAX = 0.0
                RLXB = RLX
                DO IR = 1, NRC
                    DBG = (TSCL - 1.0) * BGAM(IR, NR)
                    DBGMAX = MAX(DBGMAX, ABS(DBG))
                    IF (BGMAG/=0.0) THEN
                        FDBG = ABS(DBG) / BGMAG
                        IF (FDBG * RLXB>0.3) RLXB = 0.3 / FDBG
                    ENDIF
                ENDDO
                !---- Update BGAM
                DO IR = 1, NRC
                    BGFAC = (TSCL - 1.0) * BGAM(IR, NR)
                    BGAM(IR, NR) = BGAM(IR, NR) + RLXB * BGFAC
                ENDDO
                IF (TGAP>0.0 .AND. OMEGA(NR)/=0.0) BGAM(NRC, NR) = 0.0
                !
                WRITE (*, 110) ITR, RLXB, DBGMAX * RLXB
                IF (ABS(DBGMAX * RLXB)<=0.001 * BGMAG) THEN
                    LCONV = .TRUE.
                    GOTO 20
                ENDIF
                LGAMA = .FALSE.
                !
            ELSEIF (IRTYPE(NR)==2) THEN
                !---- Bladed disk
                DTDALF = 0.0
                DO IR = 1, NRC
                    !---- theta velocity at blade (use only 1/2 of induced Vt from circulation)
                    VTBG = BGAM(IR, NR) * PI2I / YRC(IR, NR)
                    WTB = VREL(3, IR, NR) - 0.5 * VTBG
                    WWB = SQRT(VREL(1, IR, NR)**2 + WTB**2)
                    PHIB = ATAN2(VREL(1, IR, NR), -WTB)
                    !
                    XI = YRC(IR, NR) / RTIP(NR)
                    ALF = BETAR(IR, NR) - PHIB
                    REY = WWB * CHR(IR, NR) * RHO / RMU
                    SECSIG = BLDS * CHR(IR, NR) / (2.0 * PI * YRC(IR, NR))
                    SECSTAGR = 0.5 * PI - BETAR(IR, NR)
                    CALL GETCLCDCM(NR, IR, XI, ALF, WWB, REY, SECSIG, SECSTAGR, CLB, &
                            & CL_ALF, CL_W, CLMAX, CLMIN, DCL_STALL, &
                            & LSTALLR(IR, NR), CDR(IR, NR), CD_ALF, CD_W, &
                            & CD_REY, CMOM, CM_AL, CM_W)
                    CLR(IR, NR) = CLB
                    CLALF(IR, NR) = CL_ALF
                    ALFAR(IR, NR) = ALF
                    !
                    IF (IR==NRC .AND. TGAP>0.0 .AND. OMEGA(NR)/=0.0) THEN
                        CLR(IR, NR) = 0.0
                        CLALF(IR, NR) = 0.0
                        ALFAR(IR, NR) = 0.0
                    ENDIF
                    !
                    DTDALF = DTDALF + TI_GAM(IR, NR) * 0.5 * WWB * CHR(IR, NR) * CL_ALF
                ENDDO
                !---- Change blade pitch to drive thrust
                RLXB = RLXF
                DBETA = 0.0
                IF (DTDALF/=0.0) DBETA = (TSPEC - THR) / DTDALF
                !---- Limit DBETA changes in iteration to get desired rotor thrust
                IF (ABS(DBETA) * RLXB>0.1) RLXB = 0.1 / ABS(DBETA)
                !---- update BETA by estimate
                DO IR = 1, NRC
                    BETAR(IR, NR) = BETAR(IR, NR) + RLXB * DBETA
                ENDDO
                !
                WRITE (*, 100) ITR, RLXB, DBETA / DTR
                IF (ABS(DBETA)<0.001) THEN
                    LCONV = .TRUE.
                    GOTO 20
                ENDIF
                LGAMA = .FALSE.
                !
            ENDIF
            !
        ENDDO
        LCONV = .FALSE.
        !
        100  FORMAT (I3, ' RLX=', F8.5, ' dBeta=', F9.5)
        110  FORMAT (I3, ' RLX=', F8.5, ' dBGmax=', F9.5)
        !
        !---- Final iterations to converge case
        20   CALL CONVGTHBG(NITER, RLXF, WXEPS)
        !---- Update rotor velocities
        CALL UPDROTVEL
        !
    END SUBROUTINE CONVGTHBGT
    !*==UPDROTVEL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020

end module s_rotoper