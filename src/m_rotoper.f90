module m_rotoper
    implicit none
contains
    !*==ROTINITTHR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! QFCALC
    !=========================================================================
    ! DFDC (Ducted Fan Design Code) is an aerodynamic and aeroacoustic design
    ! and analysis tool for aircraft with propulsors in ducted fan
    ! configurations.
    !
    ! This software was developed under the auspices and sponsorship of the
    ! Tactical Technology Office (TTO) of the Defense Advanced Research
    ! Projects Agency (DARPA).
    !
    ! Copyright (c) 2004, 2005, Booz Allen Hamilton Inc., All Rights Reserved
    !
    ! This program is free software; you can redistribute it and/or modify it
    ! under the terms of the GNU General Public License as published by the
    ! Free Software Foundation; either version 2 of the License, or (at your
    ! option) any later version.
    !
    ! This program is distributed in the hope that it will be useful, but
    ! WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    ! General Public License for more details.
    !
    ! You should have received a copy of the GNU General Public License along
    ! with this program; if not, write to the Free Software Foundation, Inc.,
    ! 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
    !
    ! Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (drela@mit.edu)
    ! Program Management: Brad Tousley, Paul Eremenko (eremenko@alum.mit.edu)
    !
    !
    !=========================================================================
    !
    !     Version 070-ES1
    !     Philip Carter, Esotec Developments, February 2009
    !     philip (at) esotec (dot) org
    !
    !     Changes from 0.70:
    !
    !     Alpha stored in ALPHAR whenever airfoil data is called.
    !     Alpha and CDR zeroed for final station, tip gap cases.
    !     SHOWDUCT, SHOWACTDISK, SHOWBLADE...: argument (LU), formats tweaked.
    !
    !     ROTRPRT: Alpha listing added, formats tweaked.
    !     Beta and Alpha listed in local coords if LCOORD is TRUE and RPM.LE.0.
    !     Global/local message added to output for rotors of RPM.LE.0.
    !
    !     Version 070-ES1a  4 March 2009
    !
    !     ROTRPRT: CP and CT output fixed for neg-rpm disks
    !     Local alfa output bug fixed for neg rpm and tip gap
    !
    !=========================================================================


    subroutine rotinitthr(thr)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: thr
        !
        ! Local variables
        !
        real :: bgama, thrust, vaind, vhsq
        integer :: ir, nr
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Initialize rotor velocities with estimate of
        !     blade circulation and induced velocites derived from thrust
        !
        !     Assumes momentum theory for axial flow calculation
        !----------------------------------------------------------------
        !
        !---- Assume that rotor/act. disk is first disk
        nr = 1
        !---- Initialize rotor with constant circulation derived from thrust
        thrust = thr
        bgama = 2.0 * pi * thrust / (rho * adisk(nr) * omega(nr))
        if (ldbg) write (lundbg, *) 'ROTINITTHR  THRUST  BGAMA ', &
                & thrust, bgama
        !---- Induced velocity from momentum theory
        vhsq = 2.0 * thrust / (rho * adisk(nr))
        vaind = -0.5 * qinf + sqrt((0.5 * qinf)**2 + vhsq)
        !---- Set average velocity in duct
        vavginit = vaind + qinf
        !---- Set blade circulation
        do ir = 1, nrc
            bgam(ir, nr) = bgama
        enddo
        if (tgap>0.0) bgam(nrc, nr) = 0.0
        !---- Initialize using our circulation estimate
        call rotinitbgam
        !
        if (ldbg) then
            write (lundbg, *) 'ROTINITTHR'
            write (lundbg, *) ' Setting circulation from THRUST= ', thrust
            write (lundbg, *) ' Average circulation       B*GAM= ', bgama
            write (lundbg, *) ' Average axial velocity    VAavg= ', &
                    & vavginit
        endif
        !
    end subroutine rotinitthr
    !*==ROTINITBGAM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ROTINITTHR



    subroutine rotinitbgam
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: abgam, area, bgama, da, phir, thrust, vaind, &
                & vhsq, wm, wr, wt, wx
        integer :: ir, nr
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Initialize rotor velocities using current blade circulation
        !     Sets approximate axial induced velocites from thrust
        !
        !     Assumes momentum theory for axial flow calculation
        !----------------------------------------------------------------
        !
        if (ldbg) write (*, *) 'Initializing rotor GAM and induced vel.'
        !
        !---- Find area averaged B*GAMMA for rotor(s) to estimate thrust
        thrust = 0.0
        do nr = 1, nrotor
            abgam = 0.0
            area = 0.0
            do ir = 1, nrc
                da = pi * (yrp(ir + 1, nr)**2 - yrp(ir, nr)**2)
                area = area + da
                abgam = abgam + da * bgam(ir, nr)
            enddo
            bgama = abgam / area
            thrust = thrust + adisk(nr) * rho * omega(nr) * bgama * pi2i
        enddo
        !
        !     Thrust used for estimates for slipstream axial velocities
        if (ldbg) write (*, *) 'Est. rotor thrust (from B*GAM) = ', &
                & thrust
        !---- Induced velocity from momentum theory
        vhsq = 2.0 * thrust / (rho * adisk(1))
        vaind = -0.5 * qinf + sqrt((0.5 * qinf)**2 + vhsq)
        !---- Set average velocity in duct
        vavginit = vaind + qinf
        !c    VAVGINIT = SQRT(THRUST/(RHO*ADISK(1)))  ! old definition
        !
        if (ldbg) then
            write (lundbg, *) 'ROTINITBGAM'
            write (lundbg, *) ' Average circulation       B*GAM= ', bgama
            write (lundbg, *) ' Estimated                THRUST= ', thrust
            write (lundbg, *) ' Average axial velocity    VAavg= ', &
                    & vavginit
        endif
        !
        do nr = 1, nrotor
            do ir = 1, nrc
                !---- Absolute frame induced velocities
                vind(1, ir, nr) = vaind
                vind(2, ir, nr) = 0.0
                vind(3, ir, nr) = bgam(ir, nr) * pi2i / yrc(ir, nr)
            enddo
        enddo
        call setrotvel
        lvmav = .false.
        !c      CALL VMAVGINIT(VAVGINIT)
        !
        if (ldbg) then
            nr = 1
            write (*, 1400) nr
            do ir = 1, nrc
                wx = vrel(1, ir, nr)
                wr = vrel(2, ir, nr)
                wt = vrel(3, ir, nr)
                wm = sqrt(wx * wx + wr * wr)
                if (wt/=0.0) then
                    phir = atan2(wm, -wt)
                else
                    phir = 0.5 * pi
                endif
                write (*, 1401) yrc(ir, nr), wx, wr, wm, wt, phir / dtr, &
                        & bgam(ir, nr)
            enddo
        endif
        !
        1400 format (/'Blade velocities initialized on blade row ', &
                &i3/'     r          Wx         Wr         Wm', &
                &'         Wt        Phi       BGam')
        1401 format (1x, 8g11.4)
        !
    end subroutine rotinitbgam
    !*==ROTINITBLD.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ROTINITBGAM



    subroutine rotinitbld
        use i_dfdc
        use m_viscvel, only : uvinfl
        use m_inigrd, only : rotbg2grd, clrgrdflw
        use m_aero, only : getclcdcm
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: alf, bgamnew, bgamold, blds, cdrag, cd_alf, &
                & cd_rey, cd_w, ci, clb, clmax, clmin, cl_alf, &
                & cl_w, cmom, cm_al, cm_w, dcl_stall, delbgam, dr, &
                & dvaind, phi, rey, rlx, secsig, secstagr, si, &
                & tsum, vaind, vhsq, vtin, w, wsq, wwa, wwt, xi
        integer :: i, ig, iterg, nr
        integer, save :: niterg
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Sets reasonable initial circulation using current
        !     rotor blade geometry (chord, beta).
        !
        !     Initial circulations are set w/o induced effects
        !     An iteration is done using the self-induced velocity
        !     from momentum theory to converge an approximate
        !     induced axial velocity
        !----------------------------------------------------------
        !
        data niterg/10/
        !
        !---- Set up to accumulate blade/disk circulations
        call clrgrdflw
        !
        !---- This initialization assumes blade disks are in streamwise order <- FIX THIS!
        do nr = 1, nrotor
            !
            if (irtype(nr)==1) then
                call rotbg2grd(nr)
                !
            elseif (irtype(nr)==2) then
                ig = igrotor(nr)
                !
                !---- Initialize section circulation neglecting induced swirl velocity
                !---- Start with no circulation and current axial flow estimate
                vaind = vavginit - qinf
                do i = 1, nrc
                    bgam(i, nr) = 0.0
                enddo
                blds = float(nrbld(nr))
                !---- Under-relaxation to reduce transients in CL
                rlx = 0.5
                !
                do iterg = 1, niterg
                    !
                    tsum = 0.0
                    do i = 1, nrc
                        xi = yrc(i, nr) / rtip(nr)
                        dr = yrp(i + 1, nr) - yrp(i, nr)
                        !
                        !---- Use upstream circulation to calculate inflow
                        if (ig<=1) then
                            vtin = 0.0
                        else
                            vtin = bgamg(ig - 1, i) * pi2i / yrc(i, nr)
                        endif
                        si = qinf + vaind
                        ci = vtin - yrc(i, nr) * omega(nr)
                        !
                        wsq = ci * ci + si * si
                        w = sqrt(wsq)
                        phi = atan2(si, -ci)
                        !
                        alf = betar(i, nr) - phi
                        rey = chr(i, nr) * abs(w) * rho / rmu
                        secsig = blds * chr(i, nr) / (2.0 * pi * yrc(i, nr))
                        secstagr = 0.5 * pi - betar(i, nr)
                        call getclcdcm(nr, i, xi, alf, w, rey, secsig, secstagr, clb, &
                                & cl_alf, cl_w, clmax, clmin, dcl_stall, &
                                & lstallr(i, nr), cdrag, cd_alf, cd_w, cd_rey, &
                                & cmom, cm_al, cm_w)
                        clr(i, nr) = clb
                        !c        CLALF(I,NR) = CL_ALF
                        !
                        bgamnew = 0.5 * clr(i, nr) * w * chr(i, nr) * blds
                        !
                        bgamold = bgam(i, nr)
                        delbgam = bgamnew - bgamold
                        bgam(i, nr) = bgamold + rlx * delbgam
                        !
                        tsum = tsum - bgam(i, nr) * rho * ci * dr
                        !
                        !c        write(8,997) 'nr,i,alf,cl,gam,tsum ',nr,i,alf,clr(i,NR),
                        !c     &                                     bgam(i,NR),tsum
                        !
                        call uvinfl(yrc(i, nr), wwa, wwt)
                        !---- Set rotor slipstream velocities from estimates
                        vind(1, i, nr) = vaind
                        vind(2, i, nr) = 0.0
                        vind(3, i, nr) = ci + yrc(i, nr) * omega(nr) + bgam(i, nr)  &
                                & * pi2i / yrc(i, nr)
                    enddo
                    if (tgap>0.0) bgam(nrc, nr) = 0.0
                    !
                    !---- use momentum theory estimate of duct induced axial velocity to set VA
                    vhsq = tsum / (rho * adisk(nr))
                    dvaind = 0.0
                    if (nr==1) then
                        vaind = -0.5 * qinf + sqrt((0.5 * qinf)**2 + vhsq)
                    else
                        dvaind = -0.5 * si + sqrt((0.5 * si)**2 + vhsq)
                    endif
                    !
                    !---- Refine the initial guess with iteration using momentum theory
                    !     to drive the axial velocity
                    !       WRITE(*,*) 'ROTINITBLD noVind TSUM,VA ',TSUM,VAIND,DVAIND
                    !       WRITE(8,*) 'ROTINITBLD noVind TSUM,VA ',TSUM,VAIND,DVAIND
                    !
                enddo
                !
                !---- Set average velocity in duct
                vavginit = vaind + qinf
                !---- Put circulation from disk into grid flow
                call rotbg2grd(nr)
                !
            endif
            !
        enddo
        !
        call setrotvel
        lvmav = .false.
        !c      CALL VMAVGINIT(VAVGINIT)
        !
        997  format (a, ' ', i4, i4, 5(1x, f10.5))
        99   format (i5, 5(1x, f12.6))
        !      WRITE(*,*) 'ROTINITBLD No convergence'
        !
    end subroutine rotinitbld
    !*==SETROTVEL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine setrotvel
        use i_dfdc
        use m_viscvel, only : uvinfl
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: ir, nr
        real :: phia, phir, vfac, vma, vmr, vva, vvr, wm, wr, &
                & wt, wwa, wwt, wx
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Sets absolute and relative frame rotor velocities from
        !     induced velocities
        !     Assumes VIND abs. frame induced velocities (in common) are valid
        !----------------------------------------------------------------
        !     Blade blockage code, Esotec Developments, Sept 2013
        !----------------------------------------------------------------
        !
        do nr = 1, nrotor
            if (lblbl .and. .not.lbbloft(nr)) then
                write (*, 1200) nr
                lblbl = .false.
            endif
        enddo
        !
        1200 format (/, 'No blade blockage factors for Disk', i2)
        !
        do nr = 1, nrotor
            do ir = 1, nrc
                !
                if (lblbl) then
                    vfac = bbvfac(ir, nr)
                else
                    vfac = 1.0
                endif
                !
                call uvinfl(yrc(ir, nr), wwa, wwt)
                !---- Absolute frame velocities
                vabs(3, ir, nr) = vind(3, ir, nr) + wwt
                vabs(1, ir, nr) = (vind(1, ir, nr) + qinf + wwa) * vfac            ! BB entry
                !        VABS(1,IR,NR) = VIND(1,IR,NR) + QINF + WWA      ! v0.70
                vabs(2, ir, nr) = vind(2, ir, nr)
                vma = sqrt(vabs(1, ir, nr)**2 + vabs(2, ir, nr)**2)
                vva = sqrt(vma**2 + vabs(3, ir, nr)**2)
                if (vabs(3, ir, nr)/=0.0) then
                    phia = atan2(vma, -vabs(3, ir, nr))
                else
                    phia = 0.5 * pi
                endif
                !---- Relative frame velocities
                vrel(3, ir, nr) = vabs(3, ir, nr) - omega(nr) * yrc(ir, nr)
                vrel(1, ir, nr) = vabs(1, ir, nr)
                vrel(2, ir, nr) = vabs(2, ir, nr)
                vmr = sqrt(vrel(1, ir, nr)**2 + vrel(2, ir, nr)**2)
                vvr = sqrt(vmr**2 + vrel(3, ir, nr)**2)
                if (vrel(3, ir, nr)/=0.0) then
                    phir = atan2(vmr, -vrel(3, ir, nr))
                else
                    phir = 0.5 * pi
                endif
                !
                if (ldbg) then
                    write (*, 1400)
                    wx = vrel(1, ir, nr)
                    wr = vrel(2, ir, nr)
                    wt = vrel(3, ir, nr)
                    wm = sqrt(wx * wx + wr * wr)
                    if (wt/=0.0) then
                        phir = atan2(wm, -wt)
                    else
                        phir = 0.5 * pi
                    endif
                    write (*, 1401) yrc(ir, nr), wx, wr, wm, wt, &
                            & phir / dtr, bgam(ir, nr)
                endif
            enddo
            !
        enddo
        !
        1400 format (/'Blade slipstream velocities set...'/                    &
                &'     r          Wx         Wr         Wm', &
                &'         Wt        Phi       BGam')
        1401 format (1x, 8g11.4)
        !
    end subroutine setrotvel
    !*==CONVGTH.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine updrotvel
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: aint, da, dr, us, vaint
        integer :: ir, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Update blade or disk velocities based on current solution
        !     Velocities updated include:
        !         induced        velocities  VIND
        !         absolute frame velocities  VABS
        !         relative frame velocities  VREL
        !----------------------------------------------------------------
        !
        !---- get induced velocities upstream and downstream of disk
        do n = 1, nrotor
            call rotorvabs(n, vind(1, 1, n))
        enddo
        !
        !---- set disk downstream velocities in absolute and relative frames
        call setrotvel
        !
        !---- get area-averaged axial velocity over disk
        do n = 1, nrotor
            aint = 0.0
            vaint = 0.0
            do ir = 1, nrc
                dr = yrp(ir + 1, n) - yrp(ir, n)
                da = pi * (yrp(ir + 1, n)**2 - yrp(ir, n)**2)
                us = vabs(1, ir, n)
                aint = aint + da
                vaint = vaint + us * da
            enddo
            vaavg(n) = vaint / aint
            !c       write(*,*) 'n,vaavg ',n,vaavg(n)
        enddo
        !
    end subroutine updrotvel
    !*==VABS2VREL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine vabs2vrel(omeg, ya, vxa, vra, vta, vxr, vrr, vtr, vtr_vta, &
            & vtr_omg)
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: omeg, vra, vrr, vta, vtr, vtr_omg, vtr_vta, vxa, &
                & vxr, ya
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------------
        !     Calculates relative frame induced velocities from
        !     absolute frame velocities at radius YA with rotational
        !     speed OMEG.
        !--------------------------------------------------------------
        !
        vxr = vxa
        vrr = vra
        !---- Blade relative velocity includes rotational speed and swirl effects
        vtr = vta - omeg * ya
        vtr_vta = 1.0
        vtr_omg = -ya
        !
    end subroutine vabs2vrel
    !*==ROTORVABS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine rotorvabs(n, vel)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: n
        real, dimension(3, *) :: vel
        !
        ! Local variables
        !
        real :: bgave, circ
        integer :: ic, ic1, ic2, iel, ig, ir
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------------
        !     Get absolute frame induced velocities downstream of rotor
        !     line center points
        !--------------------------------------------------------------
        !
        !---- get velocities on rotor source line
        iel = ielrotor(n)
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        ig = igrotor(n)
        !
        do ic = ic1, ic2
            ir = ic - ic1 + 1
            !---- Use mean surface velocity on rotor source panel
            vel(1, ir) = qc(1, ic) - qinf
            vel(2, ir) = qc(2, ic)
            !---- Circumferential velocity downstream of rotor due to circulation
            !        IF(IG.EQ.1) THEN
            bgave = bgamg(ig, ir)
            !        ELSE
            !        BGAVE = 0.5*(BGAMG(IG-1,IR)+BGAMG(IG,IR))
            !        ENDIF
            circ = bgave
            vel(3, ir) = circ * pi2i / yc(ic)
        enddo
        !
        98   format (a, i5, 6(1x, f10.6))
        99   format (a, 2i5, 6(1x, f10.6))
        97   format (a, 3i5, 6(1x, f10.6))
        !
    end subroutine rotorvabs
    !*==GETVELABS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine getvelabs(nf, xf, yf, vel)
        use i_dfdc
        use m_qaic, only : qfcalc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! PARAMETER definitions
        !
        integer, parameter :: nfx = irx
        !
        ! Dummy arguments
        !
        integer :: nf
        real, dimension(3, *) :: vel
        real, dimension(*) :: xf, yf
        !
        ! Local variables
        !
        integer :: if, ir
        integer, dimension(nfx) :: iftype, ipfo, ipfp
        real, dimension(2, nfx) :: qf
        real, dimension(2, ipx, nfx) :: qf_gam, qf_gth, qf_sig
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------------------------------
        !     Get absolute frame induced velocities at NF points XF,YF
        !     by potential flow survey
        !--------------------------------------------------------------
        !
        !---- local arrays for calling QFCALC
        !
        !---- get velocities on specified points
        !---- assume the pointa are not on a panel
        !
        do ir = 1, nf
            ipfo(ir) = 0
            ipfp(ir) = 0
            iftype(ir) = 0
        enddo
        !
        !------ evaluate velocity components at points
        call qfcalc(1, nf, xf, yf, ipfo, ipfp, iftype, qf, qf_gam, qf_sig, qf_gth)
        !
        do if = 1, nf
            vel(1, if) = qf(1, if)
            vel(2, if) = qf(2, if)
            vel(3, if) = 0.0
        enddo
        !
    end subroutine getvelabs
    !*==PRTVEL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine prtvel(lu, line, lind, labs, lrel, nr)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        logical :: labs, lind, lrel
        character(*) :: line
        integer :: lu, nr
        !
        ! Local variables
        !
        real :: ang, phib, phir, rpm, vma, vra, vta, vtbg, vva, &
                & vxa, wmr, wrr, wtr, wwr, wxr, yy
        integer :: i
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------
        !     Print out velocities just downstream of rotor
        !     Prints either:
        !        absolute frame induced velocities
        !        absolute frame total velocities
        !        blade relative frame total velocities
        !
        !     LU is logical unit for output
        !     LINE is a title message for table
        !     Print flags LIND,LREL,LABS control what gets printed
        !     NR is disk # for rotor
        !
        !     Takes VIND from current solution for absolute frame
        !     Takes VABS from current solution for absolute frame
        !     Takes VREL from current solution for relative frame
        !----------------------------------------------------------
        !
        !---- Print velocity data
        write (lu, 20) line
        !
        rpm = 30.0 * omega(nr) / pi
        write (lu, 25) qinf, qref, omega(nr), rpm
        !
        !---- Absolute frame induced velocities
        if (lind) then
            write (lu, 30)
            do i = 1, nrc
                yy = yrc(i, nr)
                vxa = vind(1, i, nr)
                vra = vind(2, i, nr)
                vma = sqrt(vxa**2 + vra**2)
                vta = vind(3, i, nr)
                vva = sqrt(vma**2 + vta**2)
                !c          ANG = ATAN2(VTA,VMA)
                ang = atan2(vta, vxa)
                write (lu, 50) yy, vxa, vra, vma, vta, vva, ang / dtr
            enddo
        endif
        !
        !---- Absolute frame velocities
        if (labs) then
            write (lu, 35)
            do i = 1, nrc
                yy = yrc(i, nr)
                vxa = vabs(1, i, nr)
                vra = vabs(2, i, nr)
                vma = sqrt(vxa**2 + vra**2)
                vta = vabs(3, i, nr)
                vva = sqrt(vma**2 + vta**2)
                !c          ANG = ATAN2(VTA,VMA)
                ang = atan2(vta, vxa)
                write (lu, 50) yy, vxa, vra, vma, vta, vva, ang / dtr
            enddo
        endif
        !
        !---- Relative frame velocities
        if (lrel) then
            write (lu, 40)
            do i = 1, nrc
                yy = yrc(i, nr)
                wxr = vrel(1, i, nr)
                wrr = vrel(2, i, nr)
                wmr = sqrt(wxr**2 + wrr**2)
                !---- Blade relative velocity includes rotational speed
                wtr = vrel(3, i, nr)
                wwr = sqrt(wmr**2 + wtr**2)
                if (wtr/=0.0) then
                    !c            PHIR = ATAN2(WMR,-WTR)
                    phir = atan2(wxr, -wtr)
                else
                    phir = 0.5 * pi
                endif
                ang = phir
                write (lu, 50) yy, wxr, wrr, wmr, wtr, wwr, ang / dtr
            enddo
        endif
        !
        !---- Relative frame velocities on blade lifting line
        !     this assumes that the radial component of velocity is parallel
        !     to the blade span and is ignored in velocity magnitude and angle
        !     for the blade.  Radial components are printed however.
        if (lrel) then
            write (lu, 60)
            do i = 1, nrc
                yy = yrc(i, nr)
                wxr = vrel(1, i, nr)
                wrr = vrel(2, i, nr)
                wmr = sqrt(wxr**2 + wrr**2)
                !---- Blade relative velocity includes rotational speed, 1/2 induced swirl
                !---- Angle measured from plane of rotation
                vtbg = bgam(i, nr) * pi2i / yrc(i, nr)
                wtr = vrel(3, i, nr) - 0.5 * vtbg
                wwr = sqrt(wmr**2 + wtr**2)
                if (wtr/=0.0) then
                    phib = atan2(wxr, -wtr)
                else
                    phib = 0.5 * pi
                endif
                ang = phib
                write (lu, 50) yy, wxr, wrr, wmr, wtr, wwr, ang / dtr
            enddo
        endif
        !
        20   format (/, a)
        25   format (' QINF  =', f12.4, /' QREF  =', f12.4, /' OMEGA =', f12.4, &
                &/' RPM   =', f12.4)
        30   format (/'Induced vel, flow angles in absolute frame', &
                &' (downstream of disk)', /'     r          Vxi        Vri', &
                &'        Vmi        Vti        Vi     Swirl(deg)')
        35   format (/'Velocities, flow angles in absolute frame', &
                &' (downstream of disk)', /'     r          Vx         Vr', &
                &'         Vm         Vt         V     Swirl(deg)')
        40   format (/'Velocities, flow angles relative to blade frame', &
                &' (downstream of disk)', /'     r          Wx         Wr', &
                &'         Wm         Wt         W       Phi(deg)')
        60   format (/'Velocities in blade frame,', ' on blade lifting line', &
                &/'flow angle from plane of rotation', &
                &/'     r          Wx         Wr', &
                &'         Wm         Wt         W       Phi(deg)')
        !                  12345678901123456789011234567890112345678901
        50   format (7g11.4)
        !
    end subroutine prtvel
    !*==SHOWDUCT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine showduct(lu)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: lu
        !
        ! Local variables
        !
        integer :: n, nr
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Displays duct geometric data
        !-------------------------------------------------------
        !
        nr = 1
        write (lu, 100) name
        write (lu, 125) qinf, qref, vso, rho, rmu
        write (lu, 150) aname, nbel
        !
        do n = 1, nbel
            write (lu, 300) xble(n), yble(n), xbte(n), ybte(n), &
                    & xbrefe(n), ybrefe(n), volumv(n), asurfv(n), &
                    & xbcenv(n), ybcenv(n), rgyrxv(n), rgyryv(n)
        enddo
        !
        !     &  AREA2DA(NEX),XBCEN2DA(NEX),YBCEN2DA(NEX),
        !     &  EIXX2DA(NEX),EIYY2DA(NEX), EIXY2DA(NEX),
        !     &  AREA2DT(NEX),XBCEN2DT(NEX),YBCEN2DT(NEX),
        !     &  EIXX2DT(NEX),EIYY2DT(NEX), EIXY2DT(NEX),
        !
        !     &  VOLUMV(NEX), ASURFV(NEX),  XBCENV(NEX), YBCENV(NEX),
        !     &  RGYRXV(NEX), RGYRYV(NEX),
        !     &  VOLUMVT(NEX),ASURFVT(NEX),XBCENVT(NEX),YBCENVT(NEX),
        !     &  RGYRXVT(NEX), RGYRYVT(NEX),
        !
        100  format (/, ' DFDC Case: ', a30, /, 1x, 55('-'))
        125  format ('  Qinf  ', f9.4, '     Qref  ', f9.4, /, &
                &'  Speed of sound (m/s)     ', f10.3, /, &
                &'  Air density   (kg/m^3)   ', f10.5, /, &
                &'  Air viscosity (kg/m-s)    ', e11.4)
        !
        150  format (/, ' Geometry name: ', a30, /, ' Current duct geometry with', &
                & i2, ' elements:')
        !
        300  format (/, '  Xle   ', f12.5, '     Rle   ', f12.5, /, '  Xte   ', f12.5, &
                &'     Rte   ', f12.5, /, '  Xref  ', f12.5, '     Rref  ', &
                & f12.5, /, '  Vol   ', f12.6, '     Asurf ', f12.6, /, '  Xcent ', &
                & f12.6, '     Rcent ', f12.6, /, '  rGYRx ', f12.6, &
                & '     rGYRr ', f12.6)
        !
    end subroutine showduct
    !*==SHOWACTDSK.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SHOWDUCT


    subroutine showactdsk(lu)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: lu
        !
        ! Local variables
        !
        integer :: ir, n
        real :: rpm
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Displays parameters on actuator disk
        !-------------------------------------------------------
        !
        do n = 1, nrotor
            !
            if (irtype(n)==1) then
                !        WRITE(LU,100) NAME
                rpm = 30.0 * omega(n) / pi
                !
                write (lu, 120) n, omega(n), rpm, rhub(n), rtip(n), &
                        & adisk(n)
                write (lu, 200)
                !
                do ir = 1, nrc
                    write (lu, 210) yrc(ir, n), yrc(ir, n) / rtip(n), bgam(ir, n)
                enddo
            endif
            !
        enddo
        !
        100  format (/a)
        !
        120  format (/' Current Actuator Disk at Disk', i3, /'  Omega  ', f12.4, &
                &'    Rpm    ', f11.2, /'  Rhub   ', f12.5, '    Rtip   ', &
                & f11.5, /'  Aswept ', f12.5)

        200  format ('     r        r/R       B*Gamma')
        210  format (1x, 7g11.4)
        !
    end subroutine showactdsk
    !*==SHOWBLADE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SHOWACTDSK


    subroutine showblade(lu)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: lu
        !
        ! Local variables
        !
        real :: blds, rpm, sigrot
        integer :: ir, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Displays chord, blade angle and solidity distribution
        !     on rotor blade
        !-------------------------------------------------------
        !
        do n = 1, nrotor
            !
            if (irtype(n)==2) then
                rpm = 30.0 * omega(n) / pi
                !
                write (lu, 120) n, omega(n), rpm, rhub(n), rtip(n), &
                        & nrbld(n), adisk(n)
                blds = float(nrbld(n))
                write (lu, 200)
                !
                do ir = 1, nrc
                    sigrot = blds * chr(ir, n) / (2.0 * pi * yrc(ir, n))
                    write (lu, 210) yrc(ir, n), yrc(ir, n) / rtip(n), chr(ir, n) &
                            &, betar(ir, n) / dtr, sigrot
                enddo
            endif
            !
        enddo
        !
        100  format (/a)
        !
        120  format (/' Current Rotor at Disk', i3, /'  Omega  ', f12.4, &
                &'    Rpm    ', f11.2, /'  Rhub   ', f12.5, '    Rtip   ', &
                & f11.5, /'  #blades', i6, 6x, '    Aswept ', f11.5)
        200  format ('     r         r/R        Ch      Beta(deg)   Solidity')
        210  format (1x, 7g11.4)
        !

    end subroutine showblade
    !*==SHOWDRAGOBJ.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SHOWBLADE


    subroutine showdragobj(nd, lu)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: lu, nd
        !
        ! Local variables
        !
        integer :: id, n, nd1, nd2
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Displays parameters on drag objects
        !-------------------------------------------------------
        !
        if (ndobj<=0) then
            write (lu, *)
            write (lu, *) 'No drag objects defined'
            return
        endif
        !
        if (nd<=0) then
            nd1 = 1
            nd2 = ndobj
            write (lu, 120) ndobj
        else
            nd1 = nd
            nd2 = nd
        endif
        !
        do n = nd1, nd2
            write (lu, 130) n
            write (lu, 200)
            !
            do id = 1, nddef(n)
                write (lu, 210) id, xddef(id, n), yddef(id, n), cdadef(id, n)
            enddo
            !
        enddo
        !
        100  format (/a)
        110  format (' QINF  =', f12.4, /' QREF  =', f12.4, /' OMEGA =', f12.4, &
                &/' RPM   =', f12.4)
        !
        120  format (/, i2, ' Drag Area objects defined:', /)
        130  format (' Drag Area', i3)
        200  format ('    i        x          r        CDave')
        !             a1234567890112345678901123456789011234567890112345678901
        210  format (i5, f11.4, f11.4, 3x, g11.4)
        !
    end subroutine showdragobj
    !*==SETDRGOBJSRC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SHOWDRAGOBJ


    subroutine setdrgobjsrc
        use i_dfdc
        use m_spline, only : segspl, seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: cdav, vrr, vxx, yy
        integer :: i, ic, ic1, ic2, id, iel, ip, ip1, ip2, n
        real, dimension(irx) :: t1, t1s, vdm
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Sets source strength for rotor drag and drag objects
        !-------------------------------------------------------
        !
        !---- Set sources for drag objects
        if (ndobj<=0) then
            !c        WRITE(*,*) 'Drag objects not defined'
            return
            !
        else
            !---- step through all defined drag objects
            do n = 1, ndobj
                !
                if (ieldrgobj(n)<=0) then
                    write (*, *) 'Source element not defined for drag object '&
                            &, n
                    stop
                endif
                !
                !---- if drag objects are turned off clear the source strengths
                if (.not.ldrgobj) then
                    iel = ieldrgobj(n)
                    ip1 = ipfrst(iel)
                    ip2 = iplast(iel)
                    do ip = ip1, ip2
                        sigvsp(ip) = 0.0
                    enddo
                    !
                    !
                    !---- Set up spline arrays for drag object
                elseif (nddef(n)>=2) then
                    do i = 1, nddef(n)
                        t1(i) = cdadef(i, n)
                    enddo
                    !---- Spline drag definition array
                    call segspl(t1, t1s, yddef(1, n), nddef(n))
                    !
                    iel = ieldrgobj(n)
                    ic1 = icfrst(iel)
                    ic2 = iclast(iel)
                    ip1 = ipfrst(iel)
                    ip2 = iplast(iel)
                    !
                    do ic = ic1, ic2
                        vxx = qc(1, ic)
                        vrr = qc(2, ic)
                        id = ic - ic1 + 1
                        vdm(id) = sqrt(vxx**2 + vrr**2)
                    enddo
                    !
                    do ip = ip1, ip2
                        yy = yp(ip)
                        cdav = seval(yy, t1, t1s, yddef(1, n), nddef(n))
                        id = ip - ip1 + 1
                        if (ip==ip1) then
                            sigvsp(ip) = 0.5 * vdm(id) * cdav
                        elseif (ip==ip2) then
                            sigvsp(ip) = 0.5 * vdm(id - 1) * cdav
                        else
                            sigvsp(ip) = 0.25 * (vdm(id) + vdm(id - 1)) * cdav
                        endif
                        !c          WRITE(*,99) 'IP,CDAV,SIGVSP ',IP,CDAV,SIGVSP(IP)
                    enddo
                    !
                endif
                !
            enddo
        endif
        !
        99   format (a, i4, 5(1x, f11.5))
        !
    end subroutine setdrgobjsrc
    !*==SETROTORSRC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETDRGOBJSRC



    subroutine setrotorsrc
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: blds, cdave, chave, vave, vmm, vrr, vtt, vxx
        integer :: ic, ic1, ic2, iel, ip, ip1, ip2, ir, n
        real, dimension(irx) :: vdm
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------
        !     Sets source strength for rotor profile drag
        !-------------------------------------------------------
        !
        do n = 1, nrotor
            !---- Set sources for rotor drag
            if (ielrotor(n)>0) then
                !c        WRITE(*,*) 'Rotor source line not defined'
                !c         RETURN
                !
                iel = ielrotor(n)
                ic1 = icfrst(iel)
                ic2 = iclast(iel)
                ip1 = ipfrst(iel)
                ip2 = iplast(iel)
                blds = float(nrbld(n))
                !
                do ic = ic1, ic2
                    vxx = qc(1, ic)
                    vrr = qc(2, ic)
                    vmm = sqrt(vxx**2 + vrr**2)
                    ir = ic - ic1 + 1
                    vtt = 0.5 * bgam(ir, n) * pi2i / yc(ic) - omega(n) * yc(ic)
                    vdm(ir) = sqrt(vmm**2 + vtt**2)
                enddo
                !
                do ip = ip1, ip2
                    ir = ip - ip1 + 1
                    if (ip==ip1) then
                        sigvsp(ip) = 0.5 * blds * pi2i * vdm(ir) * chr(ir, n) * cdr(ir, n)
                        !c         WRITE(*,99) 'IP,W,CD,SIGVSP ',IP,VDM(IR),CDR(IR),SIGVSP(IP)
                    elseif (ip==ip2) then
                        !---- NOTE: should set tip source to 0.0 for tip gap (no blade defined here)
                        sigvsp(ip) = 0.5 * blds * pi2i * vdm(ir - 1) * chr(ir - 1, n)      &
                                & * cdr(ir - 1, n)
                        !c         WRITE(*,99) 'IP,W,CD,SIGVSP ',IP,VDM(IR-1),CDR(IR-1),SIGVSP(IP)
                    else
                        vave = 0.5 * (vdm(ir) + vdm(ir - 1))
                        cdave = 0.5 * (cdr(ir, n) + cdr(ir - 1, n))
                        chave = 0.5 * (chr(ir, n) + chr(ir - 1, n))
                        sigvsp(ip) = 0.5 * blds * pi2i * vave * chave * cdave
                        !c         WRITE(*,99) 'IP,W,CD,SIGVSP ',IP,VAVE,CDAVE,SIGVSP(IP)
                    endif
                enddo
            endif
            !
        enddo
        !
        99   format (a, i4, 5(1x, f11.5))
        !
    end subroutine setrotorsrc
    !*==VMAVGINIT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! SETROTORSRC



    subroutine vmavginit(vaxial)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: vaxial
        !
        ! Local variables
        !
        integer :: ic, ic1, ic2, iel, ip, ip1, ip2, ir
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Initializes VMAVG using axial velocity estimate
        !---------------------------------------------------------
        !
        if (ldbg) then
            write (*, *) 'VMAVGINIT called with VA=', vaxial
            write (lundbg, *) 'VMAVGINIT called with VA=', vaxial
        endif
        !
        !---- Set VMAVG on all wake element points
        do ir = 1, nrp
            iel = ir2iel(ir)
            !c      IC1 = ICFRST(IEL)
            !c      IC2 = ICLAST(IEL)
            ip1 = ipfrst(iel)
            ip2 = iplast(iel)
            do ip = ip1, ip2
                if (ir==1) then
                    vmavg(ip) = 0.5 * (vaxial + qinf)
                elseif (ir==nrp) then
                    vmavg(ip) = 0.5 * (vaxial + qinf)
                else
                    vmavg(ip) = vaxial
                endif
            enddo
        enddo
        !
        !---- Set VMAVG on CB body wake points
        !      IEL = 1
        !      IC1 = ICFRST(IEL)
        !      IC2 = ICLAST(IEL)
        !      DO IC = IC1, IC2
        !        IP1 = IPCO(IC)
        !        IP2 = IPCP(IC)
        !        IF(IP2IR(IP1).NE.0 .AND. IP2IR(IP2).NE.0) THEN
        !          VMAVG(IP1) = 0.5*(VAXIAL+QINF)
        !          VMAVG(IP2) = 0.5*(VAXIAL+QINF)
        !        ENDIF
        !      END DO
        !
        !---- Set VMAVG on duct body wake points
        !      IEL = 2
        !      IC1 = ICFRST(IEL)
        !      IC2 = ICLAST(IEL)
        !      DO IC = IC2, IC1, -1
        !        IP1 = IPCO(IC)
        !        IP2 = IPCP(IC)
        !        IF(IP2IR(IP1).NE.0 .AND. IP2IR(IP2).NE.0) THEN
        !          VMAVG(IP1) = 0.5*(VAXIAL+QINF)
        !          VMAVG(IP2) = 0.5*(VAXIAL+QINF)
        !        ENDIF
        !      END DO
        !
        if (ldbg) then
            write (lundbg, *) 'At end of VMAVGINIT'
            do iel = 1, nel
                !         IF(NETYPE(IEL).EQ.7) THEN
                ic1 = icfrst(iel)
                ic2 = iclast(iel)
                write (lundbg, *) iel, ic1, ic2
                do ic = ic1, ic2
                    ip1 = ipco(ic)
                    ip2 = ipcp(ic)
                    write (lundbg, *) 'IP1,VMavg ', ip1, vmavg(ip1)
                    write (lundbg, *) 'IP2,VMavg ', ip2, vmavg(ip2)
                enddo
                !        ENDIF
            enddo
        endif
        20   format (a, i5, 5(1x, f10.4))
        !
        lvmav = .true.
        !
    end subroutine vmavginit
    !*==VMAVGCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! VMAVGINIT


    subroutine vmavgcalc
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        integer :: ic, ic1, ic2, iel, ip1, ip2
        real :: qcx, qcy, vmav
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Calculates VMAVG at rotor wake points using current
        !     center point velocities QC(1,.), QC(2,.)
        !---------------------------------------------------------
        !
        if (ldbg) then
            write (*, *) 'Entering VMAVGCALC'
            write (lundbg, *) 'Entering VMAVGCALC'
        endif
        !
        !---- Set VMAVG on all wake points, first and last set from centers
        !     intermediate points get average from nearby panel centers
        do iel = 1, nel
            if (netype(iel)==7) then
                ic1 = icfrst(iel)
                ic2 = iclast(iel)
                do ic = ic1, ic2
                    ip1 = ipco(ic)
                    ip2 = ipcp(ic)
                    !           IF(IEL.NE.IR2IEL(1)) THEN
                    qcx = qc(1, ic)
                    qcy = qc(2, ic)
                    !           ELSE
                    !             QCX = 0.5*QCL(1,IC)
                    !             QCY = 0.5*QCL(2,IC)
                    !           ENDIF
                    vmav = sqrt(qcx**2 + qcy**2)
                    !---- limiter to keep VMAV positive (to fraction of QREF)
                    vmav = max(0.1 * qref, vmav)
                    if (ic==ic1) then
                        vmavg(ip1) = vmav
                    else
                        vmavg(ip1) = 0.5 * (vmavg(ip1) + vmav)
                    endif
                    vmavg(ip2) = vmav
                enddo
            endif
        enddo
        !
        !---- Set VMAVG on CB body vortex wake points
        !      IEL = 1
        !      IR  = 1
        !      IC1 = ICFRST(IEL)
        !      IC2 = ICLAST(IEL)
        !      DO IC = IC1, IC2
        !        IP1 = IPCO(IC)
        !        IP2 = IPCP(IC)
        !        IF(IP2IR(IP1).EQ.IR .AND. IP2IR(IP2).EQ.IR) THEN
        !c         write(77,*) 'iel,ic,ir ',iel,ic,ir
        !         QCX = QC(1,IC)
        !         QCY = QC(2,IC)
        !         VMAV = SQRT(QCX**2 + QCY**2)
        !         IF(IC.EQ.IC1) THEN
        !          VMAVG(IP1) = VMAV
        !         ELSE
        !          VMAVG(IP1) = 0.5*VMAVG(IP1) + 0.5*VMAV
        !         ENDIF
        !         VMAVG(IP2) = VMAV
        !c         write(77,*) 'ip1,vmavg ',ip1,vmavg(ip1)
        !c         write(77,*) 'ip2,vmavg ',ip2,vmavg(ip2)
        !        ENDIF
        !      END DO
        !---- Average velocities on CB TE and wake upstream point
        !      IPB = IPFRST(IEL)
        !      IPW = IPFRST(IR2IEL(IR))
        !      VMAV = 0.5*(VMAVG(IPB) + VMAVG(IPW))
        !      VMAVG(IPB) = VMAV
        !      VMAVG(IPW) = VMAV
        !
        !---- Set VMAVG on duct body vortex wake points
        !      IEL = 2
        !      IR  = NRP
        !      IC1 = ICFRST(IEL)
        !      IC2 = ICLAST(IEL)
        !      DO IC = IC2, IC1, -1
        !        IP1 = IPCO(IC)
        !        IP2 = IPCP(IC)
        !        IF(IP2IR(IP1).EQ.IR .AND. IP2IR(IP2).EQ.IR) THEN
        !         QCX = QC(1,IC)
        !         QCY = QC(2,IC)
        !         VMAV = SQRT(QCX**2 + QCY**2)
        !         IF(IC.EQ.IC2) THEN
        !          VMAVG(IP2) = VMAV
        !         ELSE
        !          VMAVG(IP2) = 0.5*VMAVG(IP1) + 0.5*VMAV
        !         ENDIF
        !         VMAVG(IP1) = VMAV
        !        ENDIF
        !      END DO
        !---- Average velocities on duct TE and wake upstream point
        !      IPB = IPLAST(IEL)
        !      IPW = IPFRST(IR2IEL(IR))
        !      VMAV = 0.5*(VMAVG(IPB) + VMAVG(IPW))
        !      VMAVG(IPB) = VMAV
        !      VMAVG(IPW) = VMAV
        !

        if (ldbg) then
            write (lundbg, *) 'At end of VMAVGCALC'
            do iel = 1, nel
                !          IF(NETYPE(IEL).EQ.7) THEN
                ic1 = icfrst(iel)
                ic2 = iclast(iel)
                write (lundbg, *) 'IEL,IC1,IC2,IPFRST(IEL),IPLAST(IEL)'
                write (lundbg, *) iel, ic1, ic2, ipfrst(iel), iplast(iel)
                do ic = ic1, ic2
                    ip1 = ipco(ic)
                    ip2 = ipcp(ic)
                    write (lundbg, *) 'IP1,VMavg ', ip1, vmavg(ip1)
                    write (lundbg, *) 'IP2,VMavg ', ip2, vmavg(ip2)
                enddo
                !          ENDIF
            enddo
        endif
        20   format (a, i5, 5(1x, f10.4))
        !
        !c      LVMAV = .TRUE.
        !
    end subroutine vmavgcalc
    !*==GTHCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! VMAVGCALC




    subroutine gthcalc(gamth)
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real, dimension(ipx) :: gamth
        !
        ! Local variables
        !
        real :: dg, dh, frac, gth1, y1
        integer :: ic, ic1, ic2, iel, ig, ip, ip1, ip2, ir
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Calculates GTH at rotor wake points using zero
        !     pressure jump relation
        !   Note: this formulation sets GTH on endpoints using
        !         center velocity (not averaged VM at endpoints)
        !
        !---------------------------------------------------------
        !
        if (ldbg) write (*, *) 'Entering GTHCALC'
        !
        do ip = 1, nptot
            gamth(ip) = 0.0
        enddo
        !
        !---- Set GAMTH on intermediate wakes
        do ir = 2, nrp - 1
            iel = ir2iel(ir)
            ic1 = icfrst(iel)
            ic2 = iclast(iel)
            do ic = ic1, ic2
                ip1 = ipco(ic)
                ip2 = ipcp(ic)
                ig = ic2ig(ic)
                dh = dhg(ig, ir) - dhg(ig, ir - 1)
                dg = 0.5 * pi2i**2 * (bgamg(ig, ir)**2 - bgamg(ig, ir - 1)**2)
                if (vmavg(ip1)/=0.0) then
                    gamth(ip1) = (dh - dg / (yp(ip1)**2)) / vmavg(ip1)
                else
                    !            WRITE(*,*) 'Zero VMavg on IP1 =',IP1,VMAVG(IP1)
                    gamth(ip1) = 0.
                endif
                if (vmavg(ip2)/=0.0) then
                    gamth(ip2) = (dh - dg / (yp(ip2)**2)) / vmavg(ip2)
                else
                    !            WRITE(*,*) 'Zero VMavg on IP2 =',IP2,VMAVG(IP2)
                    gamth(ip2) = 0.
                endif
            enddo
        enddo
        !
        !---- Set GAMTH on CB wake
        ir = 1
        iel = ir2iel(ir)
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        do ic = ic1, ic2
            ip1 = ipco(ic)
            ip2 = ipcp(ic)
            ig = ic2ig(ic)
            dh = 0.0
            !c        DH =                DHG(IG,IR)
            dg = 0.5 * pi2i**2 * (bgamg(ig, ir)**2)
            if (vmavg(ip1)/=0.0) then
                gamth(ip1) = (dh - dg / (yp(ip1)**2)) / vmavg(ip1)
            else
                write (*, *) 'Zero VMavg on CB wake IP1 =', ip1, vmavg(ip1)
                gamth(ip1) = 0.
            endif
            if (vmavg(ip2)/=0.0) then
                gamth(ip2) = (dh - dg / (yp(ip2)**2)) / vmavg(ip2)
            else
                write (*, *) 'Zero VMavg on CB wake IP2 =', ip2, vmavg(ip2)
                gamth(ip2) = 0.
            endif
        enddo
        !
        !---- Set GTH on CB to value at first wake point
        ir = 1
        iel = ir2iel(ir)
        ip1 = ipfrst(iel)
        gth1 = gamth(ip1)
        y1 = yp(ip1)
        iel = 1
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        do ic = ic1, ic2
            ip1 = ipco(ic)
            ip2 = ipcp(ic)
            if (ip2ir(ip1)==ir .and. ip2ir(ip2)==ir) then
                gamth(ip1) = 0.0
                gamth(ip2) = 0.0
                !---- taper off the GTH on CB from GTH at TE to 0.0 at rotor
                frac = 1.0 - float(ip1 - ipfrst(iel)) / float(iprotcb(1))
                gamth(ip1) = gth1 * frac
                frac = 1.0 - float(ip2 - ipfrst(iel)) / float(iprotcb(1))
                gamth(ip2) = gth1 * frac
            endif
            !
        enddo
        !
        !---- Set GAMTH on DUCT wake
        ir = nrp
        iel = ir2iel(ir)
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        do ic = ic1, ic2
            ip1 = ipco(ic)
            ip2 = ipcp(ic)
            ig = ic2ig(ic)
            dh = -dhg(ig, ir - 1)
            dg = 0.5 * pi2i**2 * (-bgamg(ig, ir - 1)**2)
            if (vmavg(ip1)/=0.0) then
                gamth(ip1) = (dh - dg / (yp(ip1)**2)) / vmavg(ip1)
            else
                write (*, *) 'Zero VMavg on duct wake IP1 =', ip1, &
                        & vmavg(ip1)
                gamth(ip1) = 0.
            endif
            if (vmavg(ip2)/=0.0) then
                gamth(ip2) = (dh - dg / (yp(ip2)**2)) / vmavg(ip2)
            else
                write (*, *) 'Zero VMavg on duct wake IP2 =', ip2, &
                        & vmavg(ip2)
                gamth(ip2) = 0.
            endif
        enddo
        !
        !---- Set GAMTH on DUCT from first wake GAMTH
        ir = nrp
        iel = ir2iel(ir)
        ip1 = ipfrst(iel)
        gth1 = gamth(ip1)
        y1 = yp(ip1)
        iel = 2
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        do ic = ic1, ic2
            ip1 = ipco(ic)
            ip2 = ipcp(ic)
            if (ip2ir(ip1)==ir .and. ip2ir(ip2)==ir) then
                !          GAMTH(IP1) = 0.0
                !          GAMTH(IP2) = 0.0
                !          GAMTH(IP1) = GTH1
                !          GAMTH(IP2) = GTH1
                !---- taper off the GTH on duct from GTH at TE to 0.0 at rotor
                frac = 1.0 - float(iplast(iel) - ip1)                         &
                        & / float(iplast(iel) - iprotdw(1) + 1)
                gamth(ip1) = gth1 * frac
                frac = 1.0 - float(iplast(iel) - ip2)                         &
                        & / float(iplast(iel) - iprotdw(1) + 1)
                gamth(ip2) = gth1 * frac
            endif
            !
        enddo
        !
        !
        if (ldbg) then
            write (lundbg, *) 'At end of GTHCALC'
            do iel = 1, nel
                ic1 = icfrst(iel)
                ic2 = iclast(iel)
                write (lundbg, *) 'IEL,IC1,IC2,IPFRST(IEL),IPLAST(IEL)'
                write (lundbg, *) iel, ic1, ic2, ipfrst(iel), iplast(iel)
                do ic = ic1, ic2
                    ip1 = ipco(ic)
                    ip2 = ipcp(ic)
                    write (lundbg, *) 'IP1,GAMTH ', ip1, gamth(ip1)
                    write (lundbg, *) 'IP2,GAMTH ', ip2, gamth(ip2)
                enddo
            enddo
        endif
        20   format (a, i5, 5(1x, f10.4))
        !
    end subroutine gthcalc
    !*==TQCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! GTHCALC





    subroutine tqcalc(itype)
        use i_dfdc
        use m_aero, only : getclcdcm, getalf
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: itype
        !
        ! Local variables
        !
        real :: alfa, al_cl, al_dbe, al_gi, al_omg, al_p, al_qnf, &
                & al_va, al_vt, al_w, bdr, blds, brdr, cd_al, &
                & cd_alf, cd_dbe, cd_gi, cd_omg, cd_qnf, cd_re, &
                & cd_rey, cd_va, cd_vt, cd_w, chnew, ch_gi, ch_omg, &
                & ch_qnf, ch_va, ch_vt, ci, ci_omg, ci_vt, clmax, &
                & clmin, cl_al, cl_dbe, cl_gi, cl_omg, cl_qnf, &
                & cl_va, cl_vt, cl_w, cm_al, cm_w, da, dclstall, &
                & dqi, dqi_gi, dqi_qnf, dqi_si, dqi_va, dqv, &
                & dqv_cd, dqv_ch, dqv_ci, dqv_dbe, dqv_gi, dqv_om, &
                & dqv_omg, dqv_qnf, dqv_va, dqv_vt, dqv_w, dr, dti, &
                & dti_ci, dti_gi, dti_omg, dti_vt, dtv, dtv_cd, &
                & dtv_ch, dtv_dbe, dtv_gi, dtv_omg, dtv_qnf, dtv_si, &
                & dtv_va, dtv_vt
        real :: dtv_w, hrwc, hrwc_ch, hrwc_w, phib, p_omg, p_qnf, &
                & p_va, p_vt, ra, rdr, rey, re_ch, re_gi, re_omg, &
                & re_qnf, re_va, re_vt, re_w, secsig, secstagr, si, &
                & si_qnf, si_va, vaa, vainf, vtinf, vtt, w, w_omg, &
                & w_qnf, w_va, w_vt, xi
        integer :: i, j, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------
        !     Sets Thrust, Torque and their sensitivities
        !     wrt  QINF, OMEGA, BETA, chord(i), VA,VT, BGAM
        !----------------------------------------------------------
        !
        !---- total forces accumulation
        ttot = 0.
        tvis = 0.
        qtot = 0.
        qvis = 0.
        ptot = 0.
        pvis = 0.
        !
        do n = 1, nrotor
            !
            !---- forces on blade row
            tinvr(n) = 0.
            qinvr(n) = 0.
            ti_omg(n) = 0.
            qi_omg(n) = 0.
            ti_qnf(n) = 0.
            qi_qnf(n) = 0.
            !
            tvisr(n) = 0.
            qvisr(n) = 0.
            tv_omg(n) = 0.
            qv_omg(n) = 0.
            tv_qnf(n) = 0.
            qv_qnf(n) = 0.
            tv_dbe(n) = 0.
            qv_dbe(n) = 0.
            !
            do i = 1, nrc
                ti_gam(i, n) = 0.
                ti_va(i, n) = 0.
                ti_vt(i, n) = 0.
                qi_gam(i, n) = 0.
                qi_va(i, n) = 0.
                qi_vt(i, n) = 0.
                tv_gam(i, n) = 0.
                tv_va(i, n) = 0.
                tv_vt(i, n) = 0.
                qv_gam(i, n) = 0.
                qv_va(i, n) = 0.
                qv_vt(i, n) = 0.
            enddo
            !
            !c      write(*,*) 'LBLDEF ',LBLDEF
            call rotorvabs(n, vind(1, 1, n))
            !
            !---- go over radial stations, setting thrust and torque
            !
            !---- NOTE: should only go over blade stations to NRC-1 for tip gap case
            !
            blds = float(nrbld(n))
            do i = 1, nrc
                !
                !---- Skip forces for tip gap (only on rotor)
                if (i==nrc .and. tgap>0.0 .and. omega(n)/=0.0) then
                    clr(nrc, n) = 0.0
                    cdr(nrc, n) = 0.0
                    alfar(nrc, n) = 0.0
                    cycle
                endif
                !
                xi = yrc(i, n) / rtip(n)
                ra = 0.5 * (yrp(i + 1, n) + yrp(i, n))
                dr = yrp(i + 1, n) - yrp(i, n)
                rdr = yrc(i, n) * dr
                bdr = blds * dr
                brdr = blds * rdr
                da = pi * ra * dr
                !
                !------ set  W(Qinf,Omeg,Va,Vt)  and  Phi(Qinf,Omeg,Va,Vt)  sensitivities
                call wcalc(n, i, vainf, vtinf, vtt, vaa, ci, ci_omg, ci_vt, si, &
                        & si_qnf, si_va, w, w_omg, w_qnf, w_vt, w_va, phib, p_omg, &
                        & p_qnf, p_vt, p_va)
                !
                alfa = betar(i, n) - phib
                al_dbe = 1.0
                al_p = -1.0
                !
                !---- Set local Mach number in relative frame
                machr(i, n) = w / vso
                !
                !
                !---- Check for rotor type, actuator disk or blade row
                !     If no blade has been defined an inviscid thrust and torque
                !     calculation will be made for the actuator disk using circulation.
                !     Otherwise blade element theory will be used to calculate forces
                !     on the blades using aerodynamic characteristics for the blade elements.
                !
                if (irtype(n)==2) then
                    !
                    if (itype==1) then
                        !------- analysis case:  fix local Beta (except for pitch change)
                        !------- set alfa(Gi,dBeta,Qinf,Omeg,Va,Vt) sensitivites
                        alfa = betar(i, n) - phib
                        alfar(i, n) = alfa
                        al_gi = 0.
                        al_dbe = 1.0
                        al_omg = -p_omg
                        al_qnf = -p_qnf
                        al_vt = -p_vt
                        al_va = -p_va
                        !
                        !------- set CL(Gi,dBeta,Qinf,Omeg,Va,Vt) sensitivites
                        rey = chr(i, n) * abs(w) * rho / rmu
                        secsig = blds * chr(i, n) / (2.0 * pi * yrc(i, n))
                        secstagr = 0.5 * pi - betar(i, n)
                        call getclcdcm(n, i, xi, alfa, w, rey, secsig, secstagr, &
                                & clr(i, n), cl_al, cl_w, clmax, clmin, &
                                & dclstall, lstallr(i, n), cdr(i, n), cd_alf, &
                                & cd_w, cd_rey, cmr(i, n), cm_al, cm_w)
                        clalf(i, n) = cl_al
                        cl_gi = cl_al * al_gi
                        cl_dbe = cl_al * al_dbe
                        cl_omg = cl_al * al_omg + cl_w * w_omg
                        cl_qnf = cl_al * al_qnf + cl_w * w_qnf
                        cl_vt = cl_al * al_vt + cl_w * w_vt
                        cl_va = cl_al * al_va + cl_w * w_va
                        !
                        !------- set c(Gi,Qinf,Omeg,Va,Vt) sensitivites  (chord is fixed)
                        ch_gi = 0.
                        ch_omg = 0.
                        ch_qnf = 0.
                        ch_vt = 0.
                        ch_va = 0.
                        !
                    elseif (itype==2) then
                        !---- design case:  fix local CL and set chord based on circulation
                        !---- update design CL arrays
                        !
                        if (omega(n)>0.0) then
                            do j = 1, nrc
                                cldes(j) = clpos(j)
                            enddo
                        else
                            do j = 1, nrc
                                cldes(j) = clneg(j)
                            enddo
                        endif
                        !
                        !------- set alfa(Gi,dBeta,Adv,Adw,Vt) sensitivites
                        clr(i, n) = cldes(i)
                        !
                        secsig = blds * chr(i, n) / (2.0 * pi * yrc(i, n))
                        secstagr = 0.5 * pi - betar(i, n)
                        call getalf(n, i, xi, secsig, secstagr, clr(i, n), w, alfa, &
                                & al_cl, al_w, lstallr(i, n))
                        !c         write(*,*) 'tq2 getalf i,cl,w,alf ',i,clr(i),w,alfa/dtr
                        !
                        al_gi = 0.
                        al_dbe = 0.
                        al_omg = al_w * w_omg
                        al_qnf = al_w * w_qnf
                        al_vt = al_w * w_vt
                        al_va = al_w * w_va
                        !
                        !------- set CL(Gi,dBeta,Adv,Adw,Vt) sensitivites
                        cl_gi = 0.
                        cl_dbe = 0.
                        cl_omg = 0.
                        cl_qnf = 0.
                        cl_vt = 0.
                        cl_va = 0.
                        !
                        !------- set c(Gi,Adv,Adw,Vt) sensitivites
                        chnew = 2.0 * bgam(i, n) / (blds * w * clr(i, n))
                        !--- Check for chord going zero or negative and use nearby station data
                        !    for this iteration
                        if (chnew<=0.0) then
                            !c           write(*,*) 'TQCALC negative chord @I = ',I,CHNEW
                            if (i==1) then
                                chr(i, n) = chr(i + 1, n)
                            elseif (i==ii) then
                                chr(i, n) = chr(i - 1, n)
                            else
                                chr(i, n) = 0.5 * (chr(i - 1, n) + chr(i + 1, n))
                            endif
                            ch_gi = 0.0
                            ch_omg = 0.0
                            ch_qnf = 0.0
                            ch_vt = 0.0
                            ch_va = 0.0
                        else
                            chr(i, n) = 2.0 * bgam(i, n) / (blds * w * clr(i, n))
                            !c          write(*,*) 'tq2 bgam,cl,ch ',bgam(i,n),clr(i,n),chr(i,n)
                            ch_gi = 2.0 / (blds * w * clr(i, n))
                            ch_omg = (-chr(i, n) / w) * w_omg
                            ch_qnf = (-chr(i, n) / w) * w_qnf
                            ch_vt = (-chr(i, n) / w) * w_vt
                            ch_va = (-chr(i, n) / w) * w_va
                        endif
                        !
                        betar(i, n) = alfa + phib
                        alfar(i, n) = alfa
                        !
                    elseif (itype==3) then
                        !------- design case:  fix local chord and set angles based on CL
                        !
                        !------- set CL(Gi,dBeta,Adv,Adw,Vt) sensitivites
                        clr(i, n) = 2.0 * bgam(i, n) / (blds * w * chr(i, n))
                        cl_gi = 2.0 / (blds * w * chr(i, n))
                        cl_dbe = 0.
                        cl_omg = (-clr(i, n) / w) * w_omg
                        cl_qnf = (-clr(i, n) / w) * w_qnf
                        cl_vt = (-clr(i, n) / w) * w_vt
                        cl_va = (-clr(i, n) / w) * w_va
                        !
                        !------- set alfa(Gi,dBeta,Adv,Adw,Vt) sensitivites
                        secsig = blds * chr(i, n) / (2.0 * pi * yrc(i, n))
                        secstagr = 0.5 * pi - betar(i, n)
                        call getalf(n, i, xi, secsig, secstagr, clr(i, n), w, alfa, &
                                & al_cl, al_w, lstallr(i, n))
                        !c         write(*,*) 'tq3 i,cl,ch,w,alf ',i,clr(i),chr(i),w,alfa/dtr
                        al_gi = al_cl * cl_gi
                        al_dbe = al_cl * cl_dbe
                        al_omg = al_cl * cl_omg + al_w * w_omg
                        al_qnf = al_cl * cl_qnf + al_w * w_qnf
                        al_vt = al_cl * cl_vt + al_w * w_vt
                        al_va = al_cl * cl_va + al_w * w_va
                        !
                        !------- set c(Gi,Adv,Adw,Vt) sensitivites
                        ch_gi = 0.
                        ch_omg = 0.
                        ch_qnf = 0.
                        ch_vt = 0.
                        ch_va = 0.
                        !
                        betar(i, n) = alfa + phib
                        alfar(i, n) = alfa
                        !
                    endif
                    !
                    !
                    !=================================================================
                    !
                    rer(i, n) = chr(i, n) * abs(w) * rho / rmu
                    re_w = chr(i, n) * rho / rmu
                    re_ch = abs(w) * rho / rmu
                    !
                    !------ set Re(Gi,Adv,Adw,Vt) sensitivites
                    re_gi = re_ch * ch_gi
                    re_omg = re_ch * ch_omg + re_w * w_omg
                    re_qnf = re_ch * ch_qnf + re_w * w_qnf
                    re_vt = re_ch * ch_vt + re_w * w_vt
                    re_va = re_ch * ch_va + re_w * w_va
                    !
                    !------ set CM and (not used at present) sensitivites
                    !------ set CD(Gi,dBeta,Adv,Adw,Vt) sensitivites
                    secsig = blds * chr(i, n) / (2.0 * pi * yrc(i, n))
                    secstagr = 0.5 * pi - betar(i, n)
                    call getclcdcm(n, i, xi, alfa, w, rer(i, n), secsig, secstagr, &
                            & clr(i, n), cl_al, cl_w, clmax, clmin, dclstall, &
                            & lstallr(i, n), cdr(i, n), cd_al, cd_w, cd_re, &
                            & cmr(i, n), cm_al, cm_w)
                    clalf(i, n) = cl_al
                    !c        write(*,97) 'tqcalc alfa,cl,cd,cm ',i,alfa,clr(i),cdr(i)
                    97            format (a, i5, 5(1x, f12.6))
                    cd_gi = cd_al * al_gi + cd_re * re_gi
                    cd_omg = cd_al * al_omg + cd_re * re_omg + cd_w * w_omg
                    cd_qnf = cd_al * al_qnf + cd_re * re_qnf + cd_w * w_qnf
                    cd_vt = cd_al * al_vt + cd_re * re_vt + cd_w * w_vt
                    cd_va = cd_al * al_va + cd_re * re_va + cd_w * w_va
                    cd_dbe = cd_al * al_dbe
                    !
                    !
                    hrwc = 0.5 * rho * w * chr(i, n)
                    hrwc_w = 0.5 * rho * chr(i, n)
                    hrwc_ch = 0.5 * rho * w
                    !
                    !*******************************************************
                    !------ Viscous Thrust & Power contributions on real prop
                    !
                    !------ dTv ( Cd , S , W , c ) sensitivites
                    dtv = -hrwc * cdr(i, n) * si * bdr
                    !
                    dtv_cd = -hrwc * si * bdr
                    dtv_si = -hrwc * cdr(i, n) * bdr
                    dtv_w = -hrwc_w * cdr(i, n) * si * bdr
                    dtv_ch = -hrwc_ch * cdr(i, n) * si * bdr
                    !
                    !------ set Tv(Gi,dBeta,Adv,Vt) sensitivites using chain rule
                    dtv_gi = dtv_cd * cd_gi + dtv_ch * ch_gi
                    dtv_dbe = dtv_cd * cd_dbe
                    dtv_omg = dtv_cd * cd_omg + dtv_ch * ch_omg + dtv_w * w_omg
                    dtv_qnf = dtv_cd * cd_qnf + dtv_ch * ch_qnf + dtv_w * w_qnf
                    dtv_vt = dtv_cd * cd_vt + dtv_ch * ch_vt + dtv_w * w_vt
                    dtv_va = dtv_cd * cd_va + dtv_ch * ch_va + dtv_si * si_va + &
                            & dtv_w * w_va
                    !
                    !------ accumulate viscous Thrust and sensitivities
                    tvisr(n) = tvisr(n) + dtv
                    tv_omg(n) = tv_omg(n) + dtv_omg
                    tv_qnf(n) = tv_qnf(n) + dtv_qnf
                    tv_dbe(n) = tv_dbe(n) + dtv_dbe
                    !
                    tv_gam(i, n) = dtv_gi
                    tv_va(i, n) = dtv_va
                    tv_vt(i, n) = dtv_vt
                    !
                    !
                    !------ dQv( Cd , C , W , c )
                    dqv = -hrwc * cdr(i, n) * ci * brdr
                    !
                    dqv_cd = -hrwc * ci * brdr
                    dqv_ci = -hrwc * cdr(i, n) * brdr
                    dqv_w = -hrwc_w * cdr(i, n) * ci * brdr
                    dqv_ch = -hrwc_ch * cdr(i, n) * ci * brdr
                    dqv_om = -hrwc * cdr(i, n) * ci * brdr
                    !
                    !------ set Pv(Gi,dBeta,Adv,Vt) sensitivites using chain rule
                    dqv_gi = dqv_cd * cd_gi + dqv_ch * ch_gi
                    dqv_dbe = dqv_cd * cd_dbe
                    dqv_omg = dqv_om + dqv_cd * cd_omg + dqv_ch * ch_omg + &
                            & dqv_ci * ci_omg + dqv_w * w_omg
                    dqv_qnf = dqv_cd * cd_qnf + dqv_ch * ch_qnf + dqv_w * w_qnf
                    dqv_vt = dqv_cd * cd_vt + dqv_ch * ch_vt + dqv_ci * ci_vt + &
                            & dqv_w * w_vt
                    dqv_va = dqv_cd * cd_va + dqv_ch * ch_va + dqv_w * w_va
                    !
                    !------ accumulate viscous Power and sensitivities
                    qvisr(n) = qvisr(n) + dqv
                    qv_omg(n) = qv_omg(n) + dqv_omg
                    qv_qnf(n) = qv_qnf(n) + dqv_qnf
                    qv_dbe(n) = qv_dbe(n) + dqv_dbe
                    !
                    qv_gam(i, n) = dqv_gi
                    qv_va(i, n) = dqv_va
                    qv_vt(i, n) = dqv_vt
                    !
                endif ! end of check for IRTYPE=2 (blade defined)
                !
                !
                !*******************************************************
                !------ Inviscid Thrust & Power contributions on rotor
                !
                !------ dTi( Gi , C( Omg Vt ) )
                dti = -rho * bgam(i, n) * ci * dr
                !
                dti_ci = -rho * bgam(i, n) * dr
                dti_gi = -rho * ci * dr
                !
                !------ dTi( Adv , Vt(Adw Gj) )
                dti_vt = dti_ci * ci_vt
                dti_omg = dti_ci * ci_omg
                !
                !------ accumulate inviscid Thrust and sensitivities
                tinvr(n) = tinvr(n) + dti
                ti_omg(n) = ti_omg(n) + dti_omg
                !------ Resolve dTi dependencies ( Vt ) to Gamma
                ti_gam(i, n) = dti_gi
                ti_va(i, n) = 0.0
                ti_vt(i, n) = dti_vt
                !
                !
                !------ dQi( S(V Va) , Gi )
                dqi = rho * bgam(i, n) * si * rdr
                !
                dqi_gi = rho * si * rdr
                dqi_si = rho * bgam(i, n) * rdr
                !
                !------ dQi( Vai , Qinf, Gi )
                dqi_qnf = dqi_si * si_qnf
                dqi_va = dqi_si * si_va
                !
                !------ accumulate inviscid Power and sensitivities
                qinvr(n) = qinvr(n) + dqi
                qi_omg(n) = 0.0
                qi_qnf(n) = qi_qnf(n) + dqi_qnf
                !------ Save dQi dependencies to BGAM,VA,VT
                qi_gam(i, n) = dqi_gi
                qi_va(i, n) = dqi_va
                qi_vt(i, n) = 0.0
                !
                !*******************************************************
                !------ Save blade thrust and torque distributions (per blade, per span)
                if (blds/=0.0) then
                    dtii(i, n) = dti / bdr
                    dqii(i, n) = dqi / bdr
                    dtvi(i, n) = dtv / bdr
                    dqvi(i, n) = dqv / bdr
                else
                    !------ or actuator disk thrust and torque distributions (per span)
                    dtii(i, n) = dti / dr
                    dqii(i, n) = dqi / dr
                    dtvi(i, n) = dtv / dr
                    dqvi(i, n) = dqv / dr
                endif
                !------ static pressure rise at this radial station
                dpsi(i, n) = (dti + dtv) / da
                !
            enddo
            !
            !---- forces for this rotor
            ttotr(n) = tinvr(n) + tvisr(n)
            qtotr(n) = qinvr(n) + qvisr(n)
            !
            !---- total forces (all rotors)
            ttot = ttot + ttotr(n)
            tvis = tvis + tvisr(n)
            qtot = qtot + qtotr(n)
            qvis = qvis + qvisr(n)
            !
            !---- Derive some power quantities from torque
            pvis = pvis + qvisr(n) * omega(n)
            ptot = ptot + qtotr(n) * omega(n)
            !
        enddo
        !
        !      write(*,99) 'TQcalc '
        !      write(*,99) ' Tinv,Ti_qnf,Ti_omg       ',tinv,ti_qnf,ti_omg
        !      write(*,99) ' Tvis,Tv_qnf,Tv_omg,Tv_be ',tvis,tv_qnf,tv_omg,tv_dbe
        !      write(*,99) ' Qinv,Qi_qnf,Qi_omg       ',qinv,qi_qnf,qi_omg
        !      write(*,99) ' Qvis,Qv_qnf,Qv_omg,Pv_be ',qvis,qv_qnf,qv_omg,qv_dbe
        !      write(*,99) ' Pinv,Pvis,Ptot           ',pinv,pvis,ptot
        !
        99   format (a, 4(1x, f12.3))
        !
    end subroutine tqcalc
    !*==WCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! TQCALC


    subroutine wcalc(n, i, vain, vtin, vtt, vaa, ci, ci_omg, ci_vt, si, si_qnf, &
            & si_va, w, w_omg, w_qnf, w_vt, w_va, phib, p_omg, p_qnf, &
            & p_vt, p_va)
        use i_dfdc
        use m_viscvel, only : uvinfl
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        real :: ci, ci_omg, ci_vt, phib, p_omg, p_qnf, p_va, &
                & p_vt, si, si_qnf, si_va, vaa, vain, vtin, vtt, &
                & w, w_omg, w_qnf, w_va, w_vt
        integer :: i, n
        !
        ! Local variables
        !
        real :: vfac, vtbg, wsq
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        !---- Calculate velocity components at radial station I on rotor blade
        !
        !
        if (lblbl) then
            vfac = bbvfac(i, n)
        else
            vfac = 1.0
        endif
        !
        !---- At blade lifting line use averaged circulation for tangential velocity
        vtbg = bgam(i, n) * pi2i / yrc(i, n)
        vtt = vabs(3, i, n) - 0.5 * vtbg
        vaa = vind(1, i, n)
        !
        !---- Freestream, body induced and added inflow velocities
        call uvinfl(yrc(i, n), vain, vtin)
        !
        ci = -yrc(i, n) * omega(n) + vtin + vtt
        ci_omg = -yrc(i, n)
        ci_vt = 1.0
        !
        !      SI     = QINF          + VAIN + VAA         ! v0.70
        si = (qinf + vaa + vain) * vfac              ! BB entry
        si_qnf = 1.0
        si_va = 1.0
        !
        wsq = ci * ci + si * si
        w = sqrt(wsq)
        w_omg = (ci * ci_omg) / w
        w_qnf = (si * si_qnf) / w
        w_vt = (ci * ci_vt) / w
        w_va = (si * si_va) / w
        !
        phib = atan2(si, -ci)
        p_omg = (si * ci_omg) / wsq
        p_qnf = (-ci * si_qnf) / wsq
        p_vt = (si * ci_vt) / wsq
        p_va = (-ci * si_va) / wsq
        !
        if (ldbg) then
            write (lundbg, *) 'WCALC @I= ', i
            write (lundbg, 99) 'QINF YRC  ', qinf, yrc(i, n)
            write (lundbg, 99) 'OMEG      ', omega(n)
            write (lundbg, 99) 'VT,VA     ', vtt, vaa
            write (lundbg, 99) 'VTIN,VAIN ', vtin, vain
            write (lundbg, 99) 'CI,SI,W   ', ci, si, w
            write (lundbg, 99) 'PHI       ', phib / dtr
        endif
        !
        99   format (a, 5(1x, f11.6))
        !
    end subroutine wcalc
    !*==ROTRPRT.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! WCALC


    subroutine rotrprt(lu)
        use i_dfdc
        use m_spline, only : spline, seval
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        integer :: lu
        !
        ! Local variables
        !
        real :: adv, aldeg, aref, bdeg, ch34, chi, cp, cp0, ct, &
                & ct0, ctos, dia, dref, effind, efftot, eideal, &
                & en, omegref, pc, pdim, pinv, pvdim, qdim, qinv, &
                & reexp, remax, rpm, rref, sigma, tc, tclim, tdim, &
                & tinv, tvdim, vtip, xi, xre
        integer :: i, iadd, n
        character(30) :: rtype
        character(1) :: schar
        real, dimension(irx) :: w1
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------
        !     Dumps operating state of case to unit LU
        !---------------------------------------------
        !
        write (lu, 1000)
        write (lu, 1001) name
        if (.not.lconv) write (lu, 1002)
        if (ninfl>0) write (lu, 1005)
        !
        1000 format (/1x, 76('-'))
        1001 format (' DFDC  Case:  ', a32)
        1002 format (/19x, '********** NOT CONVERGED **********', /)
        1003 format (1x, 76('-'))
        1004 format (50x)
        1005 format (' (External slipstream present)')
        !
        !---- dimensional thrust, power, torque, rpm
        tdim = ttot + tduct
        qdim = qtot
        pdim = ptot
        tvdim = tvis
        pvdim = pvis
        tinv = ttot - tvis
        qinv = qtot - qvis
        pinv = ptot - pvis
        !
        !---- Define reference quantities for coefficients
        rref = rtip(1)
        aref = adisk(1)
        omegref = omega(1)
        !
        !---- standard thrust/power coefficients based on rotational speed
        dref = 2.0 * rref
        en = abs(omegref * pi2i)
        if (en==0.0) then
            ct = 0.0
            cp = 0.0
        else
            ct = tdim / (rho * en**2 * dref**4)
            cp = pdim / (rho * en**3 * dref**5)
        endif
        !---- standard thrust/power coefficients based on forward speed
        if (qinf>0.0) then
            tc = tdim / (0.5 * rho * qinf**2 * pi * rref**2)
            pc = pdim / (0.5 * rho * qinf**3 * pi * rref**2)
        else
            tc = 0.0
            pc = 0.0
        endif
        !---- thrust/power coefficients based on tip speed
        !     uses helicopter nomenclature for CT0,CP0,FOM
        vtip = abs(omegref * rref)
        adv = qinf / vtip
        if (vtip/=0.0) then
            ct0 = tdim / (rho * aref * vtip**2)
            cp0 = pdim / (rho * aref * vtip**3)
        else
            ct0 = 0.0
            cp0 = 0.0
        endif
        if (ct0>=0.0 .and. cp0/=0.0) then
            fom = abs(ct0)**1.5 / cp0 / 2.0
        else
            fom = 0.0
        endif
        !
        !---- overall efficiency (all thrust components)
        if (pdim/=0.0) efftot = qinf * tdim / pdim
        !---- induced efficiency (all thrust components)
        if (pinv/=0.0) effind = qinf * (tinv + tduct) / pinv
        !---- ideal (actuator disk) efficiency
        if (tc==0) then
            eideal = 0.0
        else
            tclim = max(-1.0, tc)
            eideal = 2.0 / (1.0 + sqrt(tclim + 1.0))
        endif
        !
        !---- Dump overall case data
        !
        write (lu, 1003)
        if (irtype(1)/=2) then
            write (lu, 1008)
        elseif (lblbl) then
            write (lu, 1010)
        else
            write (lu, 1011)
        endif
        !
        write (lu, 1012) qinf, alth, deltat, rho, vso, rmu, tdim, &
                & pdim, efftot, tvdim, pvdim, effind, tduct, &
                & qdim, eideal
        !
        write (lu, 1004)
        !
        write (lu, 1014) aref, rref, omegref
        !---- Thrust/power coefficients based on rotational speed (propeller syntax)
        write (lu, 1015) ct, cp, adv * pi
        !---- Thrust/power coefficients based on forward speed (propeller syntax)
        write (lu, 1016) tc, pc, adv
        !---- Thrust/power coefficients based on tip speed (helicopter nomenclature)
        write (lu, 1017) ct0, cp0, fom
        !
        1008 format (' Flow Condition and total Forces')
        1010 format (' Flow Condition and total Forces', 17x, &
                &'Corrected for blade blockage')
        1011 format (' Flow Condition and total Forces', 17x, &
                &'No blade blockage correction')
        1012 format (/'  Vinf(m/s) :', f10.3, 4x, 'Alt.(km)   :', f9.3, 5x, &
                &'DeltaT(dgC):', f9.4, /' rho(kg/m3) :', f11.4, 3x, &
                &'Vsound(m/s):', f9.3, 5x, 'mu(kg/m-s) :', e11.4, &
                &/' Thrust(N)  :', g11.3, 3x, 'Power(W)   :', g11.3, 3x, &
                &'Efficiency :', f9.4, /' Tvisc (N)  :', f11.4, 3x, &
                &'Pvisc(W)   :', g11.3, 3x, 'Induced Eff:', f9.4, &
                &/' Tduct(N)   :', f11.4, 3x, 'torQue(N-m):', g11.3, 3x, &
                &'Ideal Eff  :', f9.4)
        !
        1014 format ('  Area:', f11.5, '  Radius:', f11.5, ' Omega:', f11.5, &
                &'  Reference data')
        1015 format ('    Ct:', f11.5, '      Cp:', f11.5, '     J:', f11.5, &
                &'  by(Rho,N,Dia)')
        1016 format ('    Tc:', f11.5, '      Pc:', f11.5, '   adv:', f11.5, &
                &'  by(Rho,Vinf,Area)  ')
        1017 format ('   CT0:', f11.5, '     CP0:', f11.5, '   FOM:', f11.5, &
                &'  by(Rho,R*Omg,Area)')
        !
        !
        !---- Display operating state for each rotor
        !
        do n = 1, nrotor
            !
            iadd = 1
            !c      IF(LU.EQ.LUWRIT) IADD = INCR
            !
            write (lu, 1003)
            !
            !---- dimensional thrust, power, torque, rpm
            tdim = ttotr(n)
            qdim = qtotr(n)
            pdim = qdim * omega(n)
            tvdim = tvisr(n)
            pvdim = qvisr(n) * omega(n)
            tinv = tinvr(n)
            qinv = qinvr(n)
            pinv = qinv * omega(n)
            !
            !
            !---- Define reference quantities for coefficients
            rref = rtip(n)
            aref = adisk(n)
            omegref = omega(n)
            rpm = 30.0 * omegref / pi
            !
            !---- standard thrust/power coefficients based on rotational speed
            dia = 2.0 * rref
            en = abs(omegref * pi2i)
            if (en==0.0) then
                ct = 0.0
                cp = 0.0
            else
                ct = tdim / (rho * en**2 * dia**4)
                cp = pdim / (rho * en**3 * dia**5)
            endif
            !
            !---- standard thrust/power coefficients based on forward speed
            if (qinf>0.0) then
                tc = tdim / (0.5 * rho * qinf**2 * pi * rref**2)
                pc = pdim / (0.5 * rho * qinf**3 * pi * rref**2)
            else
                tc = 0.0
                pc = 0.0
            endif
            !---- thrust/power coefficients based on tip speed
            !     uses helicopter nomenclature for CT0,CP0,FOM
            vtip = abs(omegref * rref)
            if (vtip/=0.0) then
                ct0 = tdim / (rho * aref * vtip**2)
                cp0 = pdim / (rho * aref * vtip**3)
                adv = qinf / vtip
            else
                ct0 = 0.0
                cp0 = 0.0
                adv = 0.0
            endif
            !
            !---- efficiency for rotor
            if (pdim/=0.0) efftot = qinf * tdim / pdim
            !---- induced efficiency for rotor
            if (pinv/=0.0) effind = qinf * tinv / pinv
            !---- ideal (actuator disk) efficiency
            if (tc==0) then
                eideal = 0.0
            else
                tclim = max(-1.0, tc)
                eideal = 2.0 / (1.0 + sqrt(tclim + 1.0))
            endif
            !
            sigma = 0.0
            if (irtype(n)==1) then
                if (omega(n)==0.0) then
                    rtype = 'Actuator Disk Stator'
                else
                    rtype = 'Actuator Disk Rotor '
                endif
            elseif (irtype(n)==2) then
                if (omega(n)==0.0) then
                    rtype = 'Bladed Stator '
                else
                    rtype = 'Bladed Rotor  '
                endif
                !---- Overall blade solidity (measured at 3/4 Rtip)
                call spline(chr(1, n), w1, yrc, nrc)
                ch34 = seval(0.75 * rtip(n), chr(1, n), w1, yrc, nrc)
                sigma = float(nrbld(n)) * ch34 / rtip(n) / pi
            else
                rtype = 'Undefined disk'
            endif
            !
            if (sigma/=0.0) then
                ctos = ct0 / sigma
            else
                ctos = 0.0
            endif
            !
            write (lu, 1110) n, rtype, nrbld(n), rpm, adv, tdim, &
                    & pdim, efftot, tvdim, pvdim, effind, qdim, &
                    & qinv, eideal, rtip(n), rhub(n), vaavg(n)
            !
            write (lu, 1004)
            write (lu, 1114) aref, rref, omegref
            !---- Thrust/power coefficients based on rotational speed (propeller syntax)
            write (lu, 1115) ct, cp, adv * pi
            !---- Thrust/power coefficients based on forward speed (propeller syntax)
            write (lu, 1116) tc, pc, adv
            !---- Thrust/power coefficients based on tip speed (helicopter nomenclature)
            write (lu, 1118) ct0, cp0
            write (lu, 1119) sigma, ctos
            !
            write (lu, 1003)
            !
            if (omega(n)<=0.0 .and. irtype(n)==2) then
                if (lcoord) then
                    write (lu, *) '                    -local coords-'
                else
                    write (lu, *) '                    -global coords-'
                endif
            endif
            !
            !
            if (irtype(n)==2) then
                !---- Blade is defined, print blade data
                !
                !---- Overall blade solidity (measured at 3/4 Rtip)
                call spline(chr(1, n), w1, yrc(1, n), nrc)
                ch34 = seval(0.75 * rtip(n), chr(1, n), w1, yrc, nrc)
                sigma = float(nrbld(n)) * ch34 / rtip(n) / pi
                if (sigma/=0.0) then
                    ctos = ct0 / sigma
                else
                    ctos = 0.0
                endif
                !
                !----- find maximum RE on blade
                remax = 0.0
                do i = 1, nrc
                    remax = max(rer(i, n), remax)
                enddo
                reexp = 1.0
                if (remax>=1.0e6) then
                    reexp = 6.0
                elseif (remax>=1.0e3) then
                    reexp = 3.0
                endif
                if (reexp==1.0) then
                    write (lu, 1120)
                else
                    write (lu, 1122) ifix(reexp)
                endif
                !
                !---- NOTE: should only dump blade data to NRC-1 for tip gap case
                !
                !       LSTALLR(NRC,N)=.FALSE.
                !
                do i = 1, nrc, iadd
                    xi = yrc(i, n) / rtip(n)
                    chi = chr(i, n) / rtip(n)
                    xre = rer(i, n) / (10.0**reexp)
                    !
                    if (lcoord .and. omega(n)<=0.0) then
                        bdeg = (pi - betar(i, n)) / dtr
                        aldeg = -alfar(i, n) / dtr
                    else
                        bdeg = betar(i, n) / dtr
                        aldeg = alfar(i, n) / dtr
                    endif
                    !
                    if (i==nrc .and. tgap>0.0 .and. omega(n)/=0.0)         &
                            & aldeg = 0.0
                    !
                    schar = ' '
                    if (lstallr(i, n)) then
                        if (i==nrc .and. tgap>0.0 .and. omega(n)/=0.0) then
                            schar = ' '
                        else
                            schar = 's'
                        endif
                    endif
                    !
                    write (lu, 1130) i, xi, chi, bdeg, aldeg, clr(i, n), &
                            & schar, cdr(i, n), xre, machr(i, n), &
                            & bgam(i, n)
                    !c     &    EFFI,EFFP(I)
                enddo
                !
            else
                !---- Print actuator disk datal
                write (lu, 1220)
                do i = 1, nrc, iadd
                    xi = yrc(i, n) / rtip(n)
                    write (lu, 1230) i, xi, machr(i, n), bgam(i, n)
                enddo
            endif
            !
        enddo
        !
        !c      WRITE(LU,1000)
        !c      WRITE(LU,*   ) ' '
        !
        return
        !....................................................................
        !
        1110 format (1x, 'Disk #', i3, 4x, a, /' # blades   :', i3, 11x, &
                & 'RPM        :', f11.3, 3x, 'adv. ratio :', f9.4, &
                &/' Thrust(N)  :', g11.3, 3x, 'Power(W)   :', g11.3, 3x, &
                &'Efficiency :', f9.4, /' Tvisc (N)  :', f11.4, 3x, &
                &'Pvisc(W)   :', g11.3, 3x, 'Induced Eff:', f9.4, &
                &/' torQue(N-m):', f11.4, 3x, 'Qvisc(N-m) :', g11.3, 3x, &
                &'Ideal Eff  :', f9.4, /' radius(m)  :', f9.4, 5x, &
                &'hub rad.(m):', f9.4, 5x, 'VAavg (m/s):', f9.4)
        !
        1114 format ('  Area:', f11.5, '  Radius:', f11.5, ' Omega:', f11.5, &
                &'  Reference data')
        1115 format ('    Ct:', f11.5, '      Cp:', f11.5, '     J:', f11.5, &
                &'  by(Rho,N,Dia)')
        1116 format ('    Tc:', f11.5, '      Pc:', f11.5, '   adv:', f11.5, &
                &'  by(Rho,Vinf,Area)  ')
        1118 format ('   CT0:', f11.5, '     CP0:', f11.5, 18x, &
                &'  by(Rho,R*Omg,Area)')
        1119 format (' Sigma:', f11.5, ' CT0/Sig:', f11.5)
        !
        !---- Rotor data
        1120 format ('   i   r/R    c/R     beta deg alfa     CL     CD', &
                &'      RE   ', '   Mach    B*Gam')
        1122 format ('   i   r/R    c/R     beta deg alfa     CL     CD', &
                &'    REx10^', i1, '   Mach    B*Gam')
        1130 format (2x, i2, f7.3, f8.4, f8.2, f7.2, 1x, f8.3, a1, f7.4, f8.2, f8.3, f9.3)
        !
        !---- Actuator disk data
        1220 format ('   i   r/R    c/R     beta deg alfa     CL     CD', &
                &'      RE   ', '   Mach    B*Gam')
        1230 format (2x, i2, f7.3, 8x, 8x, 1x, 7x, 1x, 8x, 1x, 7x, 7x, f8.3, f9.3)
        !
        !
        ! 1120 FORMAT(/'   i    r/R     c/R    beta(deg)',
        !     & '    CL       Cd     RE        Mach        B*Gam')
        ! 1030 FORMAT(2X,I2,F7.3,8X,8X,8X,2X,1X,F8.4,1X,
        !     &       F7.3,F10.3)
        !   i    r/R     c/R    beta(deg)alfa    CL      Cd      RE     Mach     B*Gam\
        !   i    r/R     c/R    beta(deg)alfa    CL      Cd    REx10^I  Mach     B*Gam')
        !xxiiffffff7fffffff8fffffff8xxxfffffff8xSffffff8xffffff7xxxxffffff7xfffffffff9
        !xx2i f7.3   f8.4    f8.2    3x  f8.3  x   f8.4 x  f7.2 4x    f7.3 x  f10.3
        !
    end subroutine rotrprt
    !*==NFCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ROTRPRT




    subroutine nfcalc
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: area, da, dp, dpt, dr, dt, dthr, omeg, powerm, &
                & powermt, thrust, thrustm, us, vs, vsq, vtb, ws, &
                & wtb
        integer :: ir, n
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Calculate near-field forces on rotor, momentum, power at disk
        !     Inviscid thrust is calculated from RPM and circulation
        !     This routine is approximate (inviscid only), superseded by
        !     routine TQCALC.
        !----------------------------------------------------------------
        !
        do n = 1, nrotor
            !
            !---- Calculate rotor inviscid thrust from circulation
            omeg = omega(n)
            thrust = 0.0
            do ir = 1, nrc
                dr = yrp(ir + 1, n) - yrp(ir, n)
                da = pi * (yrp(ir + 1, n)**2 - yrp(ir, n)**2)
                !---- use theta velocity at blade (1/2 of induced Vt in slipstream)
                vtb = 0.5 * bgam(ir, n) * pi2i / yrc(ir, n)
                wtb = vtb - omeg * yrc(ir, n)
                dthr = dr * rho * (-wtb) * bgam(ir, n)
                thrust = thrust + dthr
            enddo
            write (*, *) 'Thrust from rotor GAM and OMEGA', thrust
            if (ldbg) write (lundbg, *) 'Thrust from rotor GAM and OMEGA' &
                    &, thrust
            !
            !---- Near-field thrust from momentum theory using rotor velocities
            area = 0.0
            thrustm = 0.0
            powerm = 0.0
            powermt = 0.0
            do ir = 1, nrc
                dr = yrp(ir + 1, n) - yrp(ir, n)
                da = pi * (yrp(ir + 1, n)**2 - yrp(ir, n)**2)
                !
                us = vabs(1, ir, n)
                vs = vabs(2, ir, n)
                ws = vabs(3, ir, n)
                vsq = us * us + vs * vs + ws * ws
                dt = da * rho * us * (us - qinf)
                dp = da * rho * us * (0.5 * vsq - 0.5 * qinf * qinf)
                dpt = da * rho * us * (0.5 * ws**2)
                !
                area = area + da
                thrustm = thrustm + dt
                powerm = powerm + dp
                powermt = powermt + dpt
                !
            enddo
            !
            write (*, *) 'Momentum integration in near-field for rotor # ' &
                    &, n
            write (*, *) '       Area   = ', area
            write (*, *) '   mom.Thrust = ', thrustm
            write (*, *) '   mom.Power  = ', powerm
            write (*, *) ' swirl Power  = ', powermt
            !
        enddo
        !
    end subroutine nfcalc
    !*==FFCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! NFCALC




    subroutine ffcalc
        use i_dfdc
        use m_vels, only : getuv
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: area, da, dp, dpt, dt, dy, fm, fmr, pint, &
                & pinth, pintt, tint, tinth, us, usq, vs, vsq, &
                & ws, xff, yff
        integer :: ielo, ielp, ipo, ipp, ir
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Integrates forces (thrust and power at FF boundary)
        !---------------------------------------------------------
        !
        tint = 0.0
        tinth = 0.0
        pint = 0.0
        pinth = 0.0
        pintt = 0.0
        !
        do ir = 1, nrp - 1
            !
            fmr = rho * vabs(1, ir, 1) * pi * (yrp(ir + 1, 1)**2 - yrp(ir, 1)**2)
            !
            ielo = ir2iel(ir)
            ielp = ir2iel(ir + 1)
            ipo = iplast(ielo)
            ipp = iplast(ielp)
            xff = 0.5 * (xp(ipo) + xp(ipp)) - 0.5 * (xp(ipp) - xp(ipp - 1))
            yff = 0.5 * (yp(ipo) + yp(ipp))
            dy = yp(ipp) - yp(ipo)
            !
            call getuv(xff, yff, us, vs)
            ws = bgamg(ii - 1, ir) * pi2i / yff
            vsq = us * us + vs * vs + ws * ws
            !
            da = 2.0 * pi * yff * dy
            fm = da * rho * us
            dt = da * rho * us * (us - qinf)
            dp = da * rho * us * (0.5 * vsq - 0.5 * qinf * qinf)
            dpt = da * rho * us * (0.5 * ws**2)
            !
            area = area + da
            tint = tint + dt
            pint = pint + dp
            pintt = pintt + dpt
            !
            write (*, *) 'IR,FM,FMR ', ir, fm, fmr
            !
            !---- Use enthalpy, entropy and circulation at far-field
            usq = qinf * qinf - ws * ws + 2.0 * (dhg(ii - 1, ir) - dsg(ii - 1, ir))
            us = sqrt(usq)
            vsq = us * us + ws * ws
            dt = fmr * (us - qinf)
            dp = fmr * (0.5 * vsq - 0.5 * qinf * qinf)
            tinth = tinth + dt
            pinth = pinth + dp
            !
        enddo
        !
        write (*, *) 'FFCALC Area   = ', area
        write (*, *) '       Thrust = ', tint
        write (*, *) '       Power  = ', pint
        write (*, *) ' Swirl Power  = ', pintt
        write (*, *) '   FF  Thrust = ', tinth
        write (*, *) '   FF  Power  = ', pinth
        !
    end subroutine ffcalc
    !*==STGFIND.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020


    subroutine stgfind
        use i_dfdc
        implicit none
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Local variables
        !
        real :: d1, d12, d2, dfrac, qd, qdold, qfrac
        integer :: ic, ic1, ic2, iel, ip1, ip2
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Find stagnation points on CB and duct
        !     The panel center index and X,Y location are saved
        !---------------------------------------------------------
        !
        !---- Stagnation point on axisymmetric CB is simply the upstream point
        iel = 1
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        icstg(iel) = ic2
        xstg(iel) = xp(ip2)
        ystg(iel) = yp(ip2)
        !
        !---- Stagnation point on duct must be found by search in tangential velocity
        iel = 2
        ic1 = icfrst(iel)
        ic2 = iclast(iel)
        ip1 = ipfrst(iel)
        ip2 = iplast(iel)
        !
        icstg(iel) = ic1
        xstg(iel) = xp(ip1)
        ystg(iel) = yp(ip1)
        !
        do ic = ic1, ic2
            qd = -anc(2, ic) * qcr(1, ic) + anc(1, ic) * qcr(2, ic)
            !c        write(*,*) 'ic,qd ',ic,qd
            if (ic>ic1) then
                if (qd * qdold<0.0) then
                    !c            WRITE(*,*) 'Found stagnation point at IC=',IC
                    !c            WRITE(*,*) ' QD, QDOLD ',QD,QDOLD
                    !c            write(*,*) 'xc,yc ',XC(IC),YC(IC)
                    !
                    d1 = 0.5 * dsc(ic - 1)
                    d2 = 0.5 * dsc(ic)
                    d12 = d1 + d2
                    dfrac = d1 / d12
                    qfrac = qd / (qd - qdold)
                    !c            write(*,*) 'd1,d2,dfrac ',d1,d2,dfrac
                    !c            write(*,*) 'qfrac ',qfrac
                    if (qfrac<dfrac) then
                        icstg(iel) = ic - 1
                        xstg(iel) = xc(ic - 1) + (1.0 - qfrac) * d12 * (-anc(2, ic - 1))
                        ystg(iel) = yc(ic - 1) + (1.0 - qfrac) * d12 * (anc(1, ic - 1))
                    else
                        icstg(iel) = ic
                        xstg(iel) = xc(ic) - (qfrac) * d12 * (-anc(2, ic))
                        ystg(iel) = yc(ic) - (qfrac) * d12 * (anc(1, ic))
                    endif
                    exit
                endif
            endif
            qdold = qd
        enddo
        !
        if (ldbg) then
            iel = 1
            write (*, *) 'Element 1 Stag @ IC,X,Y ', icstg(iel), xstg(iel)&
                    &, ystg(iel)
            iel = 2
            write (*, *) 'Element 2 Stag @ IC,X,Y ', icstg(iel), xstg(iel)&
                    &, ystg(iel)
        endif
        !
    end subroutine stgfind
end module m_rotoper
