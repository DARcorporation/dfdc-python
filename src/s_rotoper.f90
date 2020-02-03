module s_rotoper
    implicit none
contains


    subroutine convgth(niter, rlxf, wxeps)
        use i_dfdc
        use m_rotoper, only : vmavgcalc, vmavginit, updrotvel, gthcalc
        use m_solve, only : gamsolv
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: niter
        real :: rlxf, wxeps
        !
        ! Local variables
        !
        real :: dg, dgthmax, rlx, rlxg
        real, dimension(ipx) :: dgold, gamth
        integer :: ip, ipmax, itr
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Basic solver for GTH,
        !     Uses underrelaxed iteration for fixed BGAM to converge GTH
        !----------------------------------------------------------------
        !
        !
        rlx = rlxf
        if (rlx<=0.0) rlx = 0.5
        !
        !---- check for valid solution before iterating
        if (.not.lgama) then
            if (.not.lvmav) call vmavginit(vavginit)
            !---- Generate GTH solution for current wakes
            call gthcalc(gamth)
            !---- Update wake gamma from initial solution
            do ip = 1, nptot
                gth(ip) = gamth(ip)
            enddo
            !---- Generate GAM solution for current RHS
            call gamsolv
        endif
        !
        do ip = 1, nptot
            dgold(ip) = 0.0
        enddo
        !
        !----- do specified cycles of under-relaxed iteration
        do itr = 1, niter
            !c         IF(LDBG) WRITE(*,110) ITR
            !
            !---- Set VMavg velocities at wake points
            call vmavgcalc
            !---- Generate GTH solution for current wakes
            call gthcalc(gamth)
            !---- Update wake gamma using CSOR
            ipmax = 0
            dgthmax = 0.0
            rlx = rlxf
            do ip = 1, nptot
                dg = gamth(ip) - gth(ip)
                if (abs(dg)>abs(dgthmax)) then
                    dgthmax = dg
                    ipmax = ip
                endif
                rlxg = rlx
                if (dg * dgold(ip)<0.0) rlxg = 0.6 * rlx
                if (dg * dgold(ip)>0.0) rlxg = 1.2 * rlx
                dgold(ip) = dg * rlxg
                gth(ip) = gth(ip) + rlxg * dg
            enddo
            !---- Generate GAM solution for current RHS
            lgama = .false.
            call gamsolv
            !
            !c         IF(LDBG)
            write (*, 100) itr, dgthmax, ipmax, rlx
            if (abs(dgthmax)<wxeps * qref) then
                lconv = .true.
                goto 20
            endif
        enddo
        lconv = .false.
        !
        100  format (i3, ' dGTHmax=', f9.5, ' @IP=', i5, ' RLX=', f8.5)
        120  format (1x, 7g11.4)
        !
        !---- Update rotor velocities
        20   call updrotvel
        !
    end subroutine convgth
    !*==CONVGTHT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine convgtht(niter, rlxf, wxeps, tspec, ispec)
        use i_dfdc
        use m_rotoper, only : vmavgcalc, updrotvel, tqcalc, vmavginit, gthcalc
        use m_solve, only : gamsolv
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ispec, niter
        real :: rlxf, tspec, wxeps
        !
        ! Local variables
        !
        real :: bgav, bgmag, bgmax, bgmin, blds, da, dbg, &
                & dbgmax, dg, dgthmax, dr, fdbg, rlx, rlxb, rlxg, &
                & thr, tscl
        real, dimension(ipx) :: dgold, gamth
        integer :: ip, ipmax, ir, itr, nr
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
        rlx = rlxf
        if (rlx<=0.0) rlx = 0.5
        !
        !---- check for valid solution before iterating
        if (.not.lgama) then
            if (.not.lvmav) call vmavginit(vavginit)
            !---- Generate GTH solution for current wakes
            call gthcalc(gamth)
            !---- Update wake gamma from initial solution
            do ip = 1, nptot
                gth(ip) = gamth(ip)
            enddo
            !---- Generate GAM solution for current RHS
            call gamsolv
        endif
        !
        do ip = 1, nptot
            dgold(ip) = 0.0
        enddo
        !
        !----- Iteration loop for driving thrust using B*GAM
        nr = 1
        blds = float(nrbld(nr))
        do itr = 1, niter
            !c         IF(LDBG) WRITE(*,110) ITR
            !---- Update rotor velocities
            call updrotvel
            !
            bgav = 0.0
            bgmax = bgam(1, nr)
            bgmin = bgmax
            do ir = 1, nrc
                dr = yrp(ir + 1, nr) - yrp(ir, nr)
                da = pi * (yrp(ir + 1, nr)**2 - yrp(ir, nr)**2)
                bgav = bgav + bgam(ir, nr)
                bgmax = max(bgmax, bgam(ir, nr))
                bgmin = min(bgmin, bgam(ir, nr))
            enddo
            bgav = bgav / float(nrc)
            bgmag = max(abs(bgav), bgmin, bgmax, 0.1)
            !
            !---- Calculate current rotor thrust
            call tqcalc(1)
            !---- Drive thrust from rotor for ISPEC=1, total thrust for ISPEC=2
            if (ispec==1) then
                thr = ttot
            else
                thr = ttot + tduct
            endif
            !---- Scale factor for BGAM from current value to get desired rotor thrust
            tscl = 1.0
            if (thr/=0.0) tscl = tspec / thr
            !
            !---- Check for rational relaxation factors based on BGAM changes
            dbgmax = 0.0
            rlxb = rlx
            do ir = 1, nrc
                dbg = (tscl - 1.0) * bgam(ir, nr)
                dbgmax = max(dbgmax, abs(dbg))
                if (bgmag/=0.0) then
                    fdbg = abs(dbg) / bgmag
                    if (fdbg * rlxb>0.3) rlxb = 0.3 / fdbg
                endif
            enddo
            !
            do ir = 1, nrc
                dbg = (tscl - 1.0) * bgam(ir, nr)
                !---- Update BGAM
                bgam(ir, nr) = bgam(ir, nr) + rlxb * dbg
            enddo
            if (tgap>0.0) bgam(nrc, nr) = 0.0
            !
            !---- Set VMavg velocities at wake points
            call vmavgcalc
            !---- Generate GTH estimate for updated circulations
            call gthcalc(gamth)
            ipmax = 0
            dgthmax = 0.0
            do ip = 1, nptot
                dg = gamth(ip) - gth(ip)
                if (abs(dg)>abs(dgthmax)) then
                    dgthmax = dg
                    ipmax = ip
                endif
                rlxg = rlxb
                if (dg * dgold(ip)<0.0) rlxg = 0.6 * rlxb
                if (dg * dgold(ip)>0.0) rlxg = 1.2 * rlxb
                dgold(ip) = dg * rlxg
                !---- Update GTH wake gamma using CSOR
                gth(ip) = gth(ip) + rlxg * dg
            enddo
            !
            !---- Generate GAM solution for current RHS
            lgama = .false.
            call gamsolv
            !
            !c       IF(LDBG) THEN
            !c         WRITE(*,*) ' '
            write (*, 100) itr, dgthmax, ipmax, dbgmax, rlxb
            if (abs(dgthmax)<wxeps * qref .and. abs(dbgmax)<=0.001 * bgmag)  &
                    & then
                lconv = .true.
                goto 20
            endif
            lgama = .false.
        enddo
        lconv = .false.
        !
        100  format (i3, ' dGTHmax=', f9.5, ' @IP=', i5, '  dBGmax=', f9.5, ' RLX=', &
                & f8.5)
        110  format (/'Blade velocities on iteration ', i5, &
                &/'     r          Wx         Wr', &
                &'         Wt        Phi       CL       BGam')
        120  format (1x, 7g11.4)
        !
        !---- Update rotor velocities
        20   call updrotvel
        !
    end subroutine convgtht
    !*==CONVGTHBG.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine convgthbg(niter, rlxf, wxeps)
        use i_dfdc
        use m_inigrd, only : setgrdflw
        use m_solve, only : gamsolv
        use m_rotoper, only : vmavgcalc, vmavginit, updrotvel, gthcalc
        use m_aero, only : getclcdcm
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: niter
        real :: rlxf, wxeps
        !
        ! Local variables
        !
        real :: alf, bgav, bgmag, bgmax, bgmin, blds, cd_alf, &
                & cd_rey, cd_w, clb, clmax, clmin, cl_alf, cl_w, &
                & cmom, cm_al, cm_w, dbg, dbgmax, dcl_stall, dg, &
                & dgthmax, fdbg, phib, rey, rlx, rlxb, rlxbg, &
                & rlxg, secsig, secstagr, vtbg, wtb, wwb, xi
        real, dimension(irx) :: bgx, dbgold
        real, dimension(ipx) :: dgold, gamth
        integer :: ig, ip, ipmax, ir, irmax, itr, nr
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Basic solver for GTH for a defined blade geometry
        !     Uses underrelaxed iteration to converge GTH, BGAM
        !----------------------------------------------------------------
        !
        !
        rlx = rlxf
        if (rlx<=0.0) rlx = 0.5
        !
        !---- check for valid solution before iterating
        if (.not.lgama) then
            if (.not.lvmav) call vmavginit(vavginit)
            !---- Generate GTH solution for current wakes
            call gthcalc(gamth)
            !---- Update wake gamma from initial solution
            do ip = 1, nptot
                gth(ip) = gamth(ip)
            enddo
            !---- Generate GAM solution for current RHS
            call gamsolv
        endif
        !
        do ip = 1, nptot
            dgold(ip) = 0.0
        enddo
        do ir = 1, nrc
            dbgold(ir) = 0.0
        enddo
        !
        !c      write(15,*) 'new solution '
        !----- do several cycles of under-relaxed iteration to converge
        !      BGAM and GTH from specified CL and chord
        do itr = 1, niter
            !c           write(15,*) 'iter ',itr
            !c         IF(LDBG) WRITE(*,110) ITR
            !
            !---- Generate GAM solution for current RHS
            call gamsolv
            !---- Update rotor velocities
            call updrotvel
            !
            rlxg = 0.0
            rlxbg = 0.0
            dbgmax = 0.0
            dgthmax = 0.0
            !
            do nr = 1, nrotor
                if (irtype(nr)==2) then
                    ig = igrotor(nr)
                    !
                    blds = float(nrbld(nr))
                    !---- convert to blade relative velocities
                    bgav = 0.0
                    bgmax = bgam(1, nr)
                    bgmin = bgmax
                    do ir = 1, nrc
                        !---- theta velocity at blade lifting line
                        vtbg = bgam(ir, nr) * pi2i / yrc(ir, nr)
                        wtb = vrel(3, ir, nr) - 0.5 * vtbg
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
                        wwb = sqrt(vrel(1, ir, nr)**2 + wtb**2)
                        phib = atan2(vrel(1, ir, nr), -wtb)
                        !
                        !           write(15,89) 'n,i,vin,vbg,vout,wtb,phi ',nr,ir,VTIN,
                        !     &                     VTBG,VREL(3,IR,NR),WTB,phib/dtr,vrot
                        ! 89        format(a,2i4,6F12.5)
                        !
                        xi = yrc(ir, nr) / rtip(nr)
                        alf = betar(ir, nr) - phib
                        rey = wwb * chr(ir, nr) * rho / rmu
                        secsig = blds * chr(ir, nr) / (2.0 * pi * yrc(ir, nr))
                        secstagr = 0.5 * pi - betar(ir, nr)
                        call getclcdcm(nr, ir, xi, alf, wwb, rey, secsig, secstagr, &
                                & clb, cl_alf, cl_w, clmax, clmin, dcl_stall, &
                                & lstallr(ir, nr), cdr(ir, nr), cd_alf, cd_w, &
                                & cd_rey, cmom, cm_al, cm_w)
                        clr(ir, nr) = clb
                        clalf(ir, nr) = cl_alf
                        alfar(ir, nr) = alf
                        !
                        bgx(ir) = blds * 0.5 * wwb * chr(ir, nr) * clr(ir, nr)
                        bgav = bgav + bgam(ir, nr)
                        bgmax = max(bgmax, bgam(ir, nr))
                        bgmin = min(bgmin, bgam(ir, nr))
                        !
                        !            WRITE(*,120) YRC(IR,NR),
                        !     &                   VREL(1,IR,NR),VREL(2,IR,NR),WTB,
                        !     &                   PHIB/DTR,CLR(IR,NR),BGX(IR),CLALF(IR,NR)
                        if (ldbg) write (*, 120) yrc(ir, nr), vrel(1, ir, nr), &
                                & vrel(2, ir, nr), wtb, phib / dtr, &
                                & clr(ir, nr), bgx(ir), clalf(ir, nr)
                    enddo
                    bgav = bgav / float(nrc)
                    if (bgav>=0.0) then
                        bgmag = max(bgav, bgmax, 0.1)
                    else
                        bgmag = min(bgav, bgmin, -0.1)
                    endif
                    if (tgap>0.0 .and. omega(nr)/=0.0) then
                        bgx(nrc) = 0.0
                        clr(nrc, nr) = 0.0
                        alfar(nrc, nr) = 0.0
                    endif
                    !
                    !---- Check for rational relaxation factors based on BGAM changes
                    rlxb = rlx
                    irmax = 0
                    dbgmax = 0.0
                    do ir = 1, nrc
                        dbg = (bgx(ir) - bgam(ir, nr))
                        if (abs(dbg)>abs(dbgmax)) then
                            dbgmax = dbg
                            irmax = ir
                        endif
                        if (bgmag/=0.0) then
                            !             FDBG = ABS(DBG)/BGAM(IR,NR)
                            !             IF(FDBG*RLXB.GT.0.5) RLXB = 0.5/FDBG
                            !             FDBG = DBG/BGAM(IR,NR)
                            fdbg = dbg / bgmag
                            if (fdbg * rlxb<-0.2) rlxb = -0.2 / fdbg
                            if (fdbg * rlxb>0.4) rlxb = 0.4 / fdbg
                        endif
                    enddo
                    !
                    !---- Update blade circulation using CSOR
                    do ir = 1, nrc
                        dbg = (bgx(ir) - bgam(ir, nr))
                        rlxbg = 0.5 * rlxb
                        if (dbg * dbgold(ir)<0.0) rlxbg = 0.6 * rlxb
                        !            IF(DBG*DBGOLD(IR).GT.0.0) RLXBG = 1.2*RLXB
                        dbgold(ir) = dbg * rlxbg
                        bgam(ir, nr) = bgam(ir, nr) + rlxbg * dbg
                    enddo
                    if (tgap>0.0 .and. omega(nr)/=0.0) bgam(nrc, nr) = 0.0
                    !---- Update grid flowfield
                    call setgrdflw
                endif
                !
            enddo ! loop over NROTOR
            !
            !---- Set VMavg velocities at wake points
            call vmavgcalc
            !---- Generate GTH estimate for updated circulations
            call gthcalc(gamth)
            ipmax = 0
            dgthmax = 0.0
            do ip = 1, nptot
                dg = gamth(ip) - gth(ip)
                if (abs(dg)>abs(dgthmax)) then
                    dgthmax = dg
                    ipmax = ip
                endif
                rlxg = rlx
                if (dg * dgold(ip)<0.0) rlxg = 0.6 * rlx
                if (dg * dgold(ip)>0.0) rlxg = 1.2 * rlx
                dgold(ip) = dg * rlxg
                !---- Update GTH wake gamma using CSOR
                gth(ip) = gth(ip) + rlxg * dg
            enddo
            !
            !---- Generate GAM solution for current RHS
            lgama = .false.
            !cc         CALL GAMSOLV
            !
            !c       IF(LDBG) THEN
            !c         WRITE(*,*) ' '
            if (rlxbg/=0.0) then
                write (*, 100) itr, dgthmax, ipmax, dbgmax, irmax, rlxbg
                !c           WRITE(*,100) ITR,DGTHMAX,IPMAX,DBGMAX,IRMAX,RLXBG,BGMAG
            else
                write (*, 105) itr, dgthmax, ipmax, rlxg
            endif
            if (abs(dgthmax)<wxeps * qref .and. abs(dbgmax)                 &
                    & <=0.001 * abs(bgmag)) then
                lconv = .true.
                goto 20
            endif
            lgama = .false.
        enddo
        lconv = .false.
        !
        100  format (i3, ' dGTHmax=', f10.5, ' @IP=', i4, '  dBGmax=', f10.5, ' @IR=', &
                & i4, ' RLX=', f8.5)
        !c     &          '  dBGmax=',F10.5,' @IR=',I4,' RLX=',F8.5,' BGMG=',F8.5)
        105  format (i3, ' dGTHmax=', f10.5, ' @IP=', i4, ' RLX=', f8.5)
        110  format (/'Disk velocities on iteration ', i4, &
                &/'     r          Wx         Wr', &
                &'         Wt        Phi       CL       BGam      CLalf')
        120  format (1x, 8g10.4)
        !
        !---- Update rotor velocities
        20   call updrotvel
        !
    end subroutine convgthbg
    !*==CONVGTHBGT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine convgthbgt(niter, rlxf, wxeps, tspec, ispec)
        use i_dfdc
        use m_rotoper, only : tqcalc, updrotvel
        use m_aero, only : getclcdcm
        use m_solve, only : gamsolv
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: ispec, niter
        real :: rlxf, tspec, wxeps
        !
        ! Local variables
        !
        real :: alf, bgav, bgeps, bgfac, bgmag, bgmax, bgmin, &
                & blds, cd_alf, cd_rey, cd_w, clb, clmax, clmin, &
                & cl_alf, cl_w, cmom, cm_al, cm_w, da, dbeta, dbg, &
                & dbgmax, dcl_stall, dr, dtdalf, fdbg, phib, rey, &
                & rlx, rlxb, secsig, secstagr, thr, tscl, vtbg, &
                & wtb, wwb, xi
        integer :: ir, itr, niter0, nr
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
        rlx = rlxf
        if (rlx<=0.0) rlx = 0.5
        !
        !---- check for valid solution before iterating
        !---- Generate GAM solution for current RHS
        if (.not.lgama) call gamsolv
        !
        !----- do several cycles of under-relaxed iteration to converge
        !      BGAM and GTH from specified CL and chord
        nr = 1
        blds = float(nrbld(nr))
        do itr = 1, niter
            !c         IF(LDBG) WRITE(*,110) ITR
            !
            bgeps = 20.0 * wxeps
            rlxb = rlxf
            niter0 = min(3, niter / 2)
            call convgthbg(niter, rlxb, bgeps)
            !
            call tqcalc(1)
            !---- Drive thrust from rotor for ISPEC=1, total thrust for ISPEC=2
            if (ispec==1) then
                thr = ttot
            else
                thr = ttot + tduct
            endif
            !
            !---- Check for disk type (actuator disk or bladed)
            if (irtype(nr)==1) then
                !---- Actuator disk
                bgav = 0.0
                bgmax = bgam(1, nr)
                bgmin = bgmax
                do ir = 1, nrc
                    dr = yrp(ir + 1, nr) - yrp(ir, nr)
                    da = pi * (yrp(ir + 1, nr)**2 - yrp(ir, nr)**2)
                    bgav = bgav + bgam(ir, nr)
                    bgmax = max(bgmax, bgam(ir, nr))
                    bgmin = min(bgmin, bgam(ir, nr))
                enddo
                bgav = bgav / float(nrc)
                bgmag = max(abs(bgav), bgmin, bgmax, 0.1)
                !---- Scale factor for BGAM from current value to get desired rotor thrust
                tscl = 1.0
                if (thr/=0.0) tscl = tspec / thr
                !---- Check for rational relaxation factors based on BGAM changes
                dbgmax = 0.0
                rlxb = rlx
                do ir = 1, nrc
                    dbg = (tscl - 1.0) * bgam(ir, nr)
                    dbgmax = max(dbgmax, abs(dbg))
                    if (bgmag/=0.0) then
                        fdbg = abs(dbg) / bgmag
                        if (fdbg * rlxb>0.3) rlxb = 0.3 / fdbg
                    endif
                enddo
                !---- Update BGAM
                do ir = 1, nrc
                    bgfac = (tscl - 1.0) * bgam(ir, nr)
                    bgam(ir, nr) = bgam(ir, nr) + rlxb * bgfac
                enddo
                if (tgap>0.0 .and. omega(nr)/=0.0) bgam(nrc, nr) = 0.0
                !
                write (*, 110) itr, rlxb, dbgmax * rlxb
                if (abs(dbgmax * rlxb)<=0.001 * bgmag) then
                    lconv = .true.
                    goto 20
                endif
                lgama = .false.
                !
            elseif (irtype(nr)==2) then
                !---- Bladed disk
                dtdalf = 0.0
                do ir = 1, nrc
                    !---- theta velocity at blade (use only 1/2 of induced Vt from circulation)
                    vtbg = bgam(ir, nr) * pi2i / yrc(ir, nr)
                    wtb = vrel(3, ir, nr) - 0.5 * vtbg
                    wwb = sqrt(vrel(1, ir, nr)**2 + wtb**2)
                    phib = atan2(vrel(1, ir, nr), -wtb)
                    !
                    xi = yrc(ir, nr) / rtip(nr)
                    alf = betar(ir, nr) - phib
                    rey = wwb * chr(ir, nr) * rho / rmu
                    secsig = blds * chr(ir, nr) / (2.0 * pi * yrc(ir, nr))
                    secstagr = 0.5 * pi - betar(ir, nr)
                    call getclcdcm(nr, ir, xi, alf, wwb, rey, secsig, secstagr, clb, &
                            & cl_alf, cl_w, clmax, clmin, dcl_stall, &
                            & lstallr(ir, nr), cdr(ir, nr), cd_alf, cd_w, &
                            & cd_rey, cmom, cm_al, cm_w)
                    clr(ir, nr) = clb
                    clalf(ir, nr) = cl_alf
                    alfar(ir, nr) = alf
                    !
                    if (ir==nrc .and. tgap>0.0 .and. omega(nr)/=0.0) then
                        clr(ir, nr) = 0.0
                        clalf(ir, nr) = 0.0
                        alfar(ir, nr) = 0.0
                    endif
                    !
                    dtdalf = dtdalf + ti_gam(ir, nr) * 0.5 * wwb * chr(ir, nr) * cl_alf
                enddo
                !---- Change blade pitch to drive thrust
                rlxb = rlxf
                dbeta = 0.0
                if (dtdalf/=0.0) dbeta = (tspec - thr) / dtdalf
                !---- Limit DBETA changes in iteration to get desired rotor thrust
                if (abs(dbeta) * rlxb>0.1) rlxb = 0.1 / abs(dbeta)
                !---- update BETA by estimate
                do ir = 1, nrc
                    betar(ir, nr) = betar(ir, nr) + rlxb * dbeta
                enddo
                !
                write (*, 100) itr, rlxb, dbeta / dtr
                if (abs(dbeta)<0.001) then
                    lconv = .true.
                    goto 20
                endif
                lgama = .false.
                !
            endif
            !
        enddo
        lconv = .false.
        !
        100  format (i3, ' RLX=', f8.5, ' dBeta=', f9.5)
        110  format (i3, ' RLX=', f8.5, ' dBGmax=', f9.5)
        !
        !---- Final iterations to converge case
        20   call convgthbg(niter, rlxf, wxeps)
        !---- Update rotor velocities
        call updrotvel
        !
    end subroutine convgthbgt
    !*==UPDROTVEL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020

end module s_rotoper