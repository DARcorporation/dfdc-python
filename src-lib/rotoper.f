C=========================================================================
C DFDC (Ducted Fan Design Code) is an aerodynamic and aeroacoustic design
C and analysis tool for aircraft with propulsors in ducted fan
C configurations.
C 
C This software was developed under the auspices and sponsorship of the
C Tactical Technology Office (TTO) of the Defense Advanced Research
C Projects Agency (DARPA).
C 
C Copyright (c) 2004, 2005, Booz Allen Hamilton Inc., All Rights Reserved
C
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version.
C 
C This program is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C 
C You should have received a copy of the GNU General Public License along
C with this program; if not, write to the Free Software Foundation, Inc.,
C 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
C
C Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (drela@mit.edu)
C Program Management: Brad Tousley, Paul Eremenko (eremenko@alum.mit.edu)
C
C
C=========================================================================
C
C     Version 070-ES1
C     Philip Carter, Esotec Developments, February 2009
C     philip (at) esotec (dot) org
C
C     Changes from 0.70:
C
C     Alpha stored in ALPHAR whenever airfoil data is called.
C     Alpha and CDR zeroed for final station, tip gap cases.
C     SHOWDUCT, SHOWACTDISK, SHOWBLADE...: argument (LU), formats tweaked.
C
C     ROTRPRT: Alpha listing added, formats tweaked. 
C     Beta and Alpha listed in local coords if LCOORD is TRUE and RPM.LE.0. 
C     Global/local message added to output for rotors of RPM.LE.0.
C
C     Version 070-ES1a  4 March 2009
C
C     ROTRPRT: CP and CT output fixed for neg-rpm disks
C     Local alfa output bug fixed for neg rpm and tip gap
C
C=========================================================================

 
      SUBROUTINE ROTINITTHR(THR)
C----------------------------------------------------------------
C     Initialize rotor velocities with estimate of 
C     blade circulation and induced velocites derived from thrust
C
C     Assumes momentum theory for axial flow calculation
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Assume that rotor/act. disk is first disk
      NR = 1
C---- Initialize rotor with constant circulation derived from thrust 
      THRUST = THR
      BGAMA = 2.0*PI*THRUST/(RHO*ADISK(NR)*OMEGA(NR))
      IF(LDBG) WRITE(LUNDBG,*) 'ROTINITTHR  THRUST  BGAMA ',THRUST,BGAMA
C---- Induced velocity from momentum theory
      VHSQ  = 2.0*THRUST/(RHO*ADISK(NR))
      VAIND = -0.5*QINF + SQRT((0.5*QINF)**2 + VHSQ)
C---- Set average velocity in duct
      VAVGINIT = VAIND + QINF
C---- Set blade circulation
      DO IR = 1, NRC
        BGAM(IR,NR) = BGAMA
      END DO
      IF(TGAP.GT.0.0) BGAM(NRC,NR) = 0.0
C---- Initialize using our circulation estimate
      CALL ROTINITBGAM
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'ROTINITTHR'
        WRITE(LUNDBG,*) ' Setting circulation from THRUST= ',THRUST
        WRITE(LUNDBG,*) ' Average circulation       B*GAM= ',BGAMA
        WRITE(LUNDBG,*) ' Average axial velocity    VAavg= ',VAVGINIT
      ENDIF
C
      RETURN
      END ! ROTINITTHR



      SUBROUTINE ROTINITBGAM
C----------------------------------------------------------------
C     Initialize rotor velocities using current blade circulation 
C     Sets approximate axial induced velocites from thrust
C
C     Assumes momentum theory for axial flow calculation
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) WRITE(*,*) 'Initializing rotor GAM and induced vel.'
C
C---- Find area averaged B*GAMMA for rotor(s) to estimate thrust
      THRUST = 0.0
      DO NR = 1, NROTOR
       ABGAM = 0.0
       AREA  = 0.0
       DO IR = 1, NRC
         DA = PI*(YRP(IR+1,NR)**2-YRP(IR,NR)**2)
         AREA = AREA + DA
         ABGAM = ABGAM + DA*BGAM(IR,NR)
       END DO
       BGAMA = ABGAM/AREA  
       THRUST = THRUST + ADISK(NR)*RHO * OMEGA(NR)*BGAMA*PI2I
      END DO
C
C     Thrust used for estimates for slipstream axial velocities
      IF(LDBG) WRITE(*,*) 'Est. rotor thrust (from B*GAM) = ',THRUST
C---- Induced velocity from momentum theory
      VHSQ  = 2.0*THRUST/(RHO*ADISK(1))
      VAIND = -0.5*QINF + SQRT((0.5*QINF)**2 + VHSQ)
C---- Set average velocity in duct
      VAVGINIT = VAIND + QINF
cc    VAVGINIT = SQRT(THRUST/(RHO*ADISK(1)))  ! old definition
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'ROTINITBGAM'
        WRITE(LUNDBG,*) ' Average circulation       B*GAM= ',BGAMA
        WRITE(LUNDBG,*) ' Estimated                THRUST= ',THRUST
        WRITE(LUNDBG,*) ' Average axial velocity    VAavg= ',VAVGINIT
      ENDIF
C
      DO NR = 1, NROTOR
       DO IR = 1, NRC
C---- Absolute frame induced velocities
        VIND(1,IR,NR) = VAIND 
        VIND(2,IR,NR) = 0.0 
        VIND(3,IR,NR) = BGAM(IR,NR)*PI2I/YRC(IR,NR)
       END DO
      END DO
      CALL SETROTVEL
      LVMAV = .FALSE.
cc      CALL VMAVGINIT(VAVGINIT)
C
      IF(LDBG) THEN
       NR = 1
       WRITE(*,1400) NR
       DO IR = 1, NRC
         WX = VREL(1,IR,NR)
         WR = VREL(2,IR,NR)
         WT = VREL(3,IR,NR)
         WM = SQRT(WX*WX + WR*WR)
         IF(WT.NE.0.0) THEN
           PHIR = ATAN2(WM,-WT)
         ELSE
           PHIR = 0.5*PI
         ENDIF
         WRITE(*,1401) YRC(IR,NR),WX,WR,WM,WT,PHIR/DTR,BGAM(IR,NR)
       END DO
      ENDIF
C     
 1400 FORMAT(/'Blade velocities initialized on blade row ',I3
     &       /'     r          Wx         Wr         Wm',
     &        '         Wt        Phi       BGam')
 1401 FORMAT(1X,8G11.4)
C
      RETURN
      END ! ROTINITBGAM



      SUBROUTINE ROTINITBLD
      INCLUDE 'DFDC.INC'
C---------------------------------------------------------
C     Sets reasonable initial circulation using current
C     rotor blade geometry (chord, beta).
C
C     Initial circulations are set w/o induced effects
C     An iteration is done using the self-induced velocity
C     from momentum theory to converge an approximate
C     induced axial velocity
C----------------------------------------------------------
C
      DATA NITERG / 10 /
C
C---- Set up to accumulate blade/disk circulations
      CALL CLRGRDFLW
C
C---- This initialization assumes blade disks are in streamwise order <- FIX THIS!
      DO NR = 1, NROTOR
C
      IF(IRTYPE(NR).EQ.1) THEN
       CALL ROTBG2GRD(NR)
C
      ELSEIF(IRTYPE(NR).EQ.2) THEN
       IG = IGROTOR(NR)
C
C---- Initialize section circulation neglecting induced swirl velocity
C---- Start with no circulation and current axial flow estimate 
       VAIND = VAVGINIT - QINF
       DO I = 1, NRC
         BGAM(I,NR) = 0.0
       END DO
       BLDS = FLOAT(NRBLD(NR))
C---- Under-relaxation to reduce transients in CL
       RLX  = 0.5
C
      DO 100 ITERG = 1, NITERG
C
       TSUM  = 0.0
       DO I = 1, NRC
        XI = YRC(I,NR)/RTIP(NR)
        DR = YRP(I+1,NR)-YRP(I,NR)
C
C---- Use upstream circulation to calculate inflow
        IF(IG.LE.1) THEN
          VTIN = 0.0
        ELSE
          VTIN = BGAMG(IG-1,I)*PI2I/YRC(I,NR)
        ENDIF
        SI =  QINF + VAIND
        CI = VTIN - YRC(I,NR)*OMEGA(NR)    
C
        WSQ = CI*CI + SI*SI
        W = SQRT(WSQ)
        PHI = ATAN2(SI,-CI)
C
        ALF = BETAR(I,NR) - PHI
        REY  = CHR(I,NR) * ABS(W) * RHO/RMU
        SECSIG = BLDS*CHR(I,NR)/(2.0*PI*YRC(I,NR))
        SECSTAGR = 0.5*PI - BETAR(I,NR)
        CALL GETCLCDCM(NR,I,XI,ALF,W,REY,SECSIG,SECSTAGR,
     &                 CLB,CL_ALF,CL_W,
     &                 CLMAX,CLMIN,DCL_STALL,LSTALLR(I,NR),
     &                 CDRAG,CD_ALF,CD_W,CD_REY,
     &                 CMOM,CM_AL,CM_W)
        CLR(I,NR) = CLB
cc        CLALF(I,NR) = CL_ALF
C
        BGAMNEW = 0.5*CLR(I,NR)*W*CHR(I,NR)*BLDS
C
        BGAMOLD = BGAM(I,NR)
        DELBGAM = BGAMNEW-BGAMOLD
        BGAM(I,NR) = BGAMOLD + RLX*DELBGAM
C
        TSUM = TSUM - BGAM(I,NR)*RHO*CI*DR
C
cc        write(8,997) 'nr,i,alf,cl,gam,tsum ',nr,i,alf,clr(i,NR),
cc     &                                     bgam(i,NR),tsum
C
        CALL UVINFL(YRC(I,NR),WWA,WWT)
C---- Set rotor slipstream velocities from estimates
        VIND(1,I,NR) = VAIND
        VIND(2,I,NR) = 0.0
        VIND(3,I,NR) = CI + YRC(I,NR)*OMEGA(NR) 
     &                    + BGAM(I,NR)*PI2I/YRC(I,NR)
       ENDDO
       IF(TGAP.GT.0.0) THEN
         BGAM(NRC,NR) = 0.0
       ENDIF
C
C---- use momentum theory estimate of duct induced axial velocity to set VA
       VHSQ = TSUM/(RHO*ADISK(NR))
       DVAIND = 0.0
       IF(NR.EQ.1) THEN
         VAIND = -0.5*QINF + SQRT((0.5*QINF)**2 + VHSQ)
       ELSE
         DVAIND = -0.5*SI + SQRT((0.5*SI)**2 + VHSQ)
       ENDIF
C
C---- Refine the initial guess with iteration using momentum theory
C     to drive the axial velocity
c       WRITE(*,*) 'ROTINITBLD noVind TSUM,VA ',TSUM,VAIND,DVAIND
c       WRITE(8,*) 'ROTINITBLD noVind TSUM,VA ',TSUM,VAIND,DVAIND
C
 100  CONTINUE
C
C---- Set average velocity in duct
      VAVGINIT = VAIND + QINF
C---- Put circulation from disk into grid flow
      CALL ROTBG2GRD(NR)
C
      ENDIF
C
      END DO
C
      CALL SETROTVEL
      LVMAV = .FALSE.
cc      CALL VMAVGINIT(VAVGINIT)
C
 997  format(A,' ',i4,i4,5(1x,f10.5))
 99   format(i5,5(1x,f12.6))
c      WRITE(*,*) 'ROTINITBLD No convergence'
C
      RETURN
      END


      SUBROUTINE SETROTVEL
C----------------------------------------------------------------
C     Sets absolute and relative frame rotor velocities from 
C     induced velocities 
C     Assumes VIND abs. frame induced velocities (in common) are valid
C----------------------------------------------------------------
C     Blade blockage code, Esotec Developments, Sept 2013
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DO NR = 1,NROTOR
        IF(LBLBL.AND..NOT.LBBLOFT(NR)) THEN
          WRITE(*,1200) NR
          LBLBL = .FALSE.
        ENDIF
      ENDDO
C
 1200 FORMAT(/,'No blade blockage factors for Disk',I2)
C
      DO NR = 1, NROTOR
        DO IR = 1, NRC
C
          IF(LBLBL) THEN
            VFAC = BBVFAC(IR,NR)
          ELSE
            VFAC = 1.0
          ENDIF
C
          CALL UVINFL(YRC(IR,NR),WWA,WWT)
C---- Absolute frame velocities
        VABS(3,IR,NR) =  VIND(3,IR,NR)       + WWT
        VABS(1,IR,NR) = (VIND(1,IR,NR)+QINF+WWA) * VFAC  ! BB entry
C        VABS(1,IR,NR) = VIND(1,IR,NR) + QINF + WWA      ! v0.70
        VABS(2,IR,NR) =  VIND(2,IR,NR)
        VMA = SQRT(VABS(1,IR,NR)**2 + VABS(2,IR,NR)**2)
        VVA = SQRT(VMA**2 + VABS(3,IR,NR)**2)
        IF(VABS(3,IR,NR).NE.0.0) THEN
         PHIA = ATAN2(VMA,-VABS(3,IR,NR))
        ELSE
         PHIA = 0.5*PI
        ENDIF
C---- Relative frame velocities
        VREL(3,IR,NR) = VABS(3,IR,NR) - OMEGA(NR)*YRC(IR,NR)
        VREL(1,IR,NR) = VABS(1,IR,NR)
        VREL(2,IR,NR) = VABS(2,IR,NR)
        VMR = SQRT(VREL(1,IR,NR)**2 + VREL(2,IR,NR)**2)
        VVR = SQRT(VMR**2 + VREL(3,IR,NR)**2)
        IF(VREL(3,IR,NR).NE.0.0) THEN
         PHIR = ATAN2(VMR,-VREL(3,IR,NR))
        ELSE
         PHIR = 0.5*PI
        ENDIF
C
        IF(LDBG) THEN
         WRITE(*,1400)
         WX = VREL(1,IR,NR)
         WR = VREL(2,IR,NR)
         WT = VREL(3,IR,NR)
         WM = SQRT(WX*WX + WR*WR)
         IF(WT.NE.0.0) THEN
           PHIR = ATAN2(WM,-WT)
         ELSE
           PHIR = 0.5*PI
         ENDIF
         WRITE(*,1401) YRC(IR,NR),WX,WR,WM,WT,PHIR/DTR,BGAM(IR,NR)
        ENDIF
       END DO
C
      END DO
C     
 1400 FORMAT(/'Blade slipstream velocities set...'
     &       /'     r          Wx         Wr         Wm',
     &        '         Wt        Phi       BGam')
 1401 FORMAT(1X,8G11.4)
C
      RETURN
      END





      SUBROUTINE CONVGTH(NITER,RLXF,WXEPS)
C----------------------------------------------------------------
C     Basic solver for GTH,
C     Uses underrelaxed iteration for fixed BGAM to converge GTH
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DIMENSION GAMTH(IPX), DGOLD(IPX)
C
      RLX = RLXF
      IF(RLX.LE.0.0) RLX = 0.5
C
C---- check for valid solution before iterating
      IF(.NOT.LGAMA) THEN
         IF(.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
C---- Generate GTH solution for current wakes
         CALL GTHCALC(GAMTH)
C---- Update wake gamma from initial solution
         DO IP = 1, NPTOT
           GTH(IP) = GAMTH(IP)
         END DO
C---- Generate GAM solution for current RHS
         CALL GAMSOLV
      ENDIF
C
      DO IP = 1, NPTOT
        DGOLD(IP) = 0.0
      END DO
C
C----- do specified cycles of under-relaxed iteration
      DO ITR = 1, NITER
cc         IF(LDBG) WRITE(*,110) ITR
C
C---- Set VMavg velocities at wake points
        CALL VMAVGCALC
C---- Generate GTH solution for current wakes
        CALL GTHCALC(GAMTH)
C---- Update wake gamma using CSOR
        IPMAX = 0
        DGTHMAX = 0.0
        RLX = RLXF
        DO IP = 1, NPTOT
          DG = GAMTH(IP) - GTH(IP)
          IF(ABS(DG).GT.ABS(DGTHMAX)) THEN
            DGTHMAX = DG
            IPMAX = IP
          ENDIF
          RLXG = RLX
          IF(DG*DGOLD(IP).LT.0.0) RLXG = 0.6*RLX
          IF(DG*DGOLD(IP).GT.0.0) RLXG = 1.2*RLX
          DGOLD(IP) = DG*RLXG
          GTH(IP) = GTH(IP) + RLXG*DG
        END DO
C---- Generate GAM solution for current RHS
        LGAMA = .FALSE.
        CALL GAMSOLV
C
cc         IF(LDBG) 
        WRITE(*,100) ITR,DGTHMAX,IPMAX,RLX
        IF(ABS(DGTHMAX).LT.WXEPS*QREF) THEN
          LCONV = .TRUE.
          GO TO 20
        ENDIF
      END DO
        LCONV = .FALSE.
C
 100  FORMAT(I3,' dGTHmax=',F9.5,' @IP=',I5,' RLX=',F8.5)
 120  FORMAT(1X,7G11.4)
C
C---- Update rotor velocities
 20   CALL UPDROTVEL
C
      RETURN
      END


      SUBROUTINE CONVGTHT(NITER,RLXF,WXEPS,TSPEC,ISPEC)
C----------------------------------------------------------------
C     Basic solver for B*GAM,GTH for design mode (fixed thrust)
C     Uses underrelaxed iteration for fixed BGAM to converge GTH
C     Input:  TSPEC   thrust specification
C             ISPEC   1 for rotor thrust spec, 2 for total thrust spec
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DIMENSION GAMTH(IPX), DGOLD(IPX)
      DIMENSION BGX(IRX)
C
      RLX = RLXF
      IF(RLX.LE.0.0) RLX = 0.5
C
C---- check for valid solution before iterating
      IF(.NOT.LGAMA) THEN
         IF(.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
C---- Generate GTH solution for current wakes
         CALL GTHCALC(GAMTH)
C---- Update wake gamma from initial solution
         DO IP = 1, NPTOT
           GTH(IP) = GAMTH(IP)
         END DO
C---- Generate GAM solution for current RHS
         CALL GAMSOLV
      ENDIF
C
      DO IP = 1, NPTOT
        DGOLD(IP) = 0.0
      END DO
C
C----- Iteration loop for driving thrust using B*GAM
      NR = 1
      BLDS = FLOAT(NRBLD(NR))
      DO ITR = 1, NITER
cc         IF(LDBG) WRITE(*,110) ITR
C---- Update rotor velocities
         CALL UPDROTVEL
C
         BGAV  = 0.0
         BGMAX = BGAM(1,NR)
         BGMIN = BGMAX
         DO IR = 1, NRC
           DR = YRP(IR+1,NR)-YRP(IR,NR)
           DA = PI*(YRP(IR+1,NR)**2-YRP(IR,NR)**2)
           BGAV = BGAV + BGAM(IR,NR)
           BGMAX = MAX(BGMAX,BGAM(IR,NR))
           BGMIN = MIN(BGMIN,BGAM(IR,NR))
         END DO
         BGAV = BGAV / FLOAT(NRC)
         BGMAG = MAX(ABS(BGAV),BGMIN,BGMAX,0.1)
C
C---- Calculate current rotor thrust
         CALL TQCALC(1)
C---- Drive thrust from rotor for ISPEC=1, total thrust for ISPEC=2
         IF(ISPEC.EQ.1) THEN
           THR = TTOT 
         ELSE
           THR = TTOT + TDUCT
         ENDIF
C---- Scale factor for BGAM from current value to get desired rotor thrust
         TSCL = 1.0
         IF(THR.NE.0.0)  TSCL = TSPEC/THR
C
C---- Check for rational relaxation factors based on BGAM changes
         DBGMAX = 0.0
         RLXB = RLX
         DO IR = 1, NRC
           DBG = (TSCL-1.0)*BGAM(IR,NR)
           DBGMAX = MAX(DBGMAX,ABS(DBG))
           IF(BGMAG.NE.0.0) THEN 
             FDBG = ABS(DBG)/BGMAG
             IF(FDBG*RLXB.GT.0.3) RLXB = 0.3/FDBG
           ENDIF
         END DO
C
         DO IR = 1, NRC
           DBG = (TSCL-1.0)*BGAM(IR,NR)
C---- Update BGAM
           BGAM(IR,NR) = BGAM(IR,NR) + RLXB * DBG
         END DO
         IF(TGAP.GT.0.0) BGAM(NRC,NR) = 0.0
C
C---- Set VMavg velocities at wake points
         CALL VMAVGCALC
C---- Generate GTH estimate for updated circulations
         CALL GTHCALC(GAMTH)
         IPMAX = 0
         DGTHMAX = 0.0
         DO IP = 1, NPTOT
           DG = GAMTH(IP) - GTH(IP)
           IF(ABS(DG).GT.ABS(DGTHMAX)) THEN
             DGTHMAX = DG
             IPMAX = IP
           ENDIF
           RLXG = RLXB
           IF(DG*DGOLD(IP).LT.0.0) RLXG = 0.6*RLXB
           IF(DG*DGOLD(IP).GT.0.0) RLXG = 1.2*RLXB
           DGOLD(IP) = DG*RLXG
C---- Update GTH wake gamma using CSOR
           GTH(IP) = GTH(IP) + RLXG*DG
         END DO
C
C---- Generate GAM solution for current RHS
         LGAMA = .FALSE.
         CALL GAMSOLV
C
cc       IF(LDBG) THEN 
cc         WRITE(*,*) ' '
         WRITE(*,100) ITR,DGTHMAX,IPMAX,DBGMAX,RLXB
         IF(ABS(DGTHMAX).LT.WXEPS*QREF .AND.
     &      ABS(DBGMAX).LE.0.001*BGMAG) THEN
           LCONV = .TRUE.
           GO TO 20
         ENDIF
         LGAMA = .FALSE.
        END DO
        LCONV = .FALSE.
C
 100  FORMAT(I3,' dGTHmax=',F9.5,' @IP=',I5,
     &          '  dBGmax=',F9.5,' RLX=',F8.5)
 110  FORMAT(/'Blade velocities on iteration ',I5,
     &       /'     r          Wx         Wr',
     &        '         Wt        Phi       CL       BGam')
 120  FORMAT(1X,7G11.4)
C
C---- Update rotor velocities
 20   CALL UPDROTVEL
C
      RETURN
      END




      SUBROUTINE CONVGTHBG(NITER,RLXF,WXEPS)
C----------------------------------------------------------------
C     Basic solver for GTH for a defined blade geometry
C     Uses underrelaxed iteration to converge GTH, BGAM
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DIMENSION GAMTH(IPX), DGOLD(IPX)
      DIMENSION BGX(IRX),   DBGOLD(IRX)
C
      RLX = RLXF
      IF(RLX.LE.0.0) RLX = 0.5
C
C---- check for valid solution before iterating
      IF(.NOT.LGAMA) THEN
         IF(.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
C---- Generate GTH solution for current wakes
         CALL GTHCALC(GAMTH)
C---- Update wake gamma from initial solution
         DO IP = 1, NPTOT
           GTH(IP) = GAMTH(IP)
         END DO
C---- Generate GAM solution for current RHS
         CALL GAMSOLV
      ENDIF
C
      DO IP = 1, NPTOT
        DGOLD(IP) = 0.0
      END DO
      DO IR = 1, NRC
        DBGOLD(IR) = 0.0
      END DO
C
cc      write(15,*) 'new solution '
C----- do several cycles of under-relaxed iteration to converge 
C      BGAM and GTH from specified CL and chord
      DO ITR = 1, NITER
cc           write(15,*) 'iter ',itr
cc         IF(LDBG) WRITE(*,110) ITR
C
C---- Generate GAM solution for current RHS
        CALL GAMSOLV
C---- Update rotor velocities
        CALL UPDROTVEL
C
        RLXG  = 0.0
        RLXBG = 0.0
        DBGMAX  = 0.0
        DGTHMAX = 0.0
C
        DO NR = 1, NROTOR
         IF(IRTYPE(NR).EQ.2) THEN
          IG = IGROTOR(NR)
C
          BLDS = FLOAT(NRBLD(NR))
C---- convert to blade relative velocities
          BGAV  = 0.0
          BGMAX = BGAM(1,NR)
          BGMIN = BGMAX
          DO IR = 1, NRC
C---- theta velocity at blade lifting line
           VTBG = BGAM(IR,NR)*PI2I/YRC(IR,NR)
           WTB  = VREL(3,IR,NR) - 0.5*VTBG
C
c           IF(NR.EQ.1) THEN
c             CIRC = 0.0
c           ELSE
c             CIRC = BGAMG(IG-1,IR)
c           ENDIF
c           VTBG = BGAM(IR,NR)*PI2I/YRC(IR,NR)
c           VTIN = CIRC*PI2I/YRC(IR,NR)
c           VROT = - OMEGA(NR)*YRC(IR,NR)
c           WTB = VTIN - OMEGA(NR)*YRC(IR,NR) + 0.5*VTBG 
C
           WWB  = SQRT(VREL(1,IR,NR)**2 + WTB**2)
           PHIB = ATAN2(VREL(1,IR,NR),-WTB)
C
c           write(15,89) 'n,i,vin,vbg,vout,wtb,phi ',nr,ir,VTIN,
c     &                     VTBG,VREL(3,IR,NR),WTB,phib/dtr,vrot
c 89        format(a,2i4,6F12.5)
C
           XI = YRC(IR,NR)/RTIP(NR)
           ALF = BETAR(IR,NR) - PHIB
           REY = WWB*CHR(IR,NR)*RHO/RMU 
           SECSIG = BLDS*CHR(IR,NR)/(2.0*PI*YRC(IR,NR))
           SECSTAGR = 0.5*PI - BETAR(IR,NR)
           CALL GETCLCDCM(NR,IR,XI,ALF,WWB,REY,SECSIG,SECSTAGR,
     &                    CLB,CL_ALF,CL_W,
     &                    CLMAX,CLMIN,DCL_STALL,LSTALLR(IR,NR),
     &                    CDR(IR,NR),CD_ALF,CD_W,CD_REY,
     &                    CMOM,CM_AL,CM_W)
           CLR(IR,NR) = CLB
           CLALF(IR,NR) = CL_ALF
           ALFAR(IR,NR) = ALF
C
           BGX(IR) = BLDS*0.5*WWB*CHR(IR,NR)*CLR(IR,NR)
           BGAV = BGAV + BGAM(IR,NR)
           BGMAX = MAX(BGMAX,BGAM(IR,NR))
           BGMIN = MIN(BGMIN,BGAM(IR,NR))
C
           IF(LDBG) THEN 
c            WRITE(*,120) YRC(IR,NR),
c     &                   VREL(1,IR,NR),VREL(2,IR,NR),WTB,
c     &                   PHIB/DTR,CLR(IR,NR),BGX(IR),CLALF(IR,NR)
            WRITE(*,120) YRC(IR,NR),
     &                   VREL(1,IR,NR),VREL(2,IR,NR),WTB,
     &                   PHIB/DTR,CLR(IR,NR),BGX(IR),CLALF(IR,NR)
           ENDIF
          END DO
          BGAV = BGAV / FLOAT(NRC)
          IF(BGAV.GE.0.0) THEN
            BGMAG = MAX(BGAV,BGMAX,0.1)
          ELSE
            BGMAG = MIN(BGAV,BGMIN,-0.1)
          ENDIF
          IF(TGAP.GT.0.0 .AND. OMEGA(NR).NE.0.0) THEN
            BGX(NRC)    = 0.0
            CLR(NRC,NR) = 0.0
            ALFAR(NRC,NR)= 0.0
          ENDIF
C
C---- Check for rational relaxation factors based on BGAM changes
          RLXB = RLX
          IRMAX = 0
          DBGMAX = 0.0
          DO IR = 1, NRC
           DBG = (BGX(IR)-BGAM(IR,NR))
           IF(ABS(DBG).GT.ABS(DBGMAX)) THEN
             DBGMAX = DBG
             IRMAX = IR
           ENDIF
           IF(BGMAG.NE.0.0) THEN 
c             FDBG = ABS(DBG)/BGAM(IR,NR)
c             IF(FDBG*RLXB.GT.0.5) RLXB = 0.5/FDBG
c             FDBG = DBG/BGAM(IR,NR)
             FDBG = DBG/BGMAG
             IF(FDBG*RLXB.LT.-0.2) RLXB = -0.2/FDBG
             IF(FDBG*RLXB.GT. 0.4) RLXB =  0.4/FDBG
           ENDIF
          END DO
C
C---- Update blade circulation using CSOR
          DO IR = 1, NRC
            DBG = (BGX(IR)-BGAM(IR,NR))
            RLXBG = 0.5*RLXB
            IF(DBG*DBGOLD(IR).LT.0.0) RLXBG = 0.6*RLXB
c            IF(DBG*DBGOLD(IR).GT.0.0) RLXBG = 1.2*RLXB
            DBGOLD(IR) = DBG*RLXBG
            BGAM(IR,NR) = BGAM(IR,NR) + RLXBG*DBG
          END DO
          IF(TGAP.GT.0.0 .AND. OMEGA(NR).NE.0.0) BGAM(NRC,NR) = 0.0
C---- Update grid flowfield
          CALL SETGRDFLW
         ENDIF
C
         END DO ! loop over NROTOR
C
C---- Set VMavg velocities at wake points
         CALL VMAVGCALC
C---- Generate GTH estimate for updated circulations
         CALL GTHCALC(GAMTH)
         IPMAX = 0
         DGTHMAX = 0.0
         DO IP = 1, NPTOT
           DG = GAMTH(IP) - GTH(IP)
           IF(ABS(DG).GT.ABS(DGTHMAX)) THEN
             DGTHMAX = DG
             IPMAX = IP
           ENDIF
           RLXG = RLX
           IF(DG*DGOLD(IP).LT.0.0) RLXG = 0.6*RLX
           IF(DG*DGOLD(IP).GT.0.0) RLXG = 1.2*RLX
           DGOLD(IP) = DG*RLXG
C---- Update GTH wake gamma using CSOR
           GTH(IP) = GTH(IP) + RLXG*DG
         END DO
C
C---- Generate GAM solution for current RHS
         LGAMA = .FALSE.
ccc         CALL GAMSOLV
C
cc       IF(LDBG) THEN 
cc         WRITE(*,*) ' '
         IF(RLXBG.NE.0.0) THEN
           WRITE(*,100) ITR,DGTHMAX,IPMAX,DBGMAX,IRMAX,RLXBG
cc           WRITE(*,100) ITR,DGTHMAX,IPMAX,DBGMAX,IRMAX,RLXBG,BGMAG
         ELSE
           WRITE(*,105) ITR,DGTHMAX,IPMAX,RLXG
         ENDIF
         IF(ABS(DGTHMAX).LT.WXEPS*QREF .AND.
     &      ABS(DBGMAX).LE.0.001*ABS(BGMAG)) THEN
           LCONV = .TRUE.
           GO TO 20
         ENDIF
         LGAMA = .FALSE.
        END DO
        LCONV = .FALSE.
C
 100  FORMAT(I3,' dGTHmax=',F10.5,' @IP=',I4,
     &          '  dBGmax=',F10.5,' @IR=',I4,' RLX=',F8.5)
cc     &          '  dBGmax=',F10.5,' @IR=',I4,' RLX=',F8.5,' BGMG=',F8.5)
 105  FORMAT(I3,' dGTHmax=',F10.5,' @IP=',I4,' RLX=',F8.5)
 110  FORMAT(/'Disk velocities on iteration ',I4,
     &       /'     r          Wx         Wr',
     &        '         Wt        Phi       CL       BGam      CLalf')
 120  FORMAT(1X,8G10.4)
C
C---- Update rotor velocities
 20   CALL UPDROTVEL
C
      RETURN
      END



      SUBROUTINE CONVGTHBGT(NITER,RLXF,WXEPS,TSPEC,ISPEC)
C----------------------------------------------------------------
C     Basic solver for BETA, GTH for analysis mode (fixed thrust)
C     Uses underrelaxed iteration for blade pitch to converge BGAM,GTH
C     Input:  TSPEC   thrust specification
C             ISPEC   1 for rotor thrust spec, 2 for total thrust spec
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DIMENSION GAMTH(IPX)
      DIMENSION BGX(IRX)
C
      RLX = RLXF
      IF(RLX.LE.0.0) RLX = 0.5
C
C---- check for valid solution before iterating
      IF(.NOT.LGAMA) THEN
C---- Generate GAM solution for current RHS
         CALL GAMSOLV
      ENDIF
C
C----- do several cycles of under-relaxed iteration to converge 
C      BGAM and GTH from specified CL and chord
      NR = 1
      BLDS = FLOAT(NRBLD(NR))
      DO ITR = 1, NITER
cc         IF(LDBG) WRITE(*,110) ITR
C
        BGEPS = 20.0*WXEPS
        RLXB = RLXF
        NITER0 = MIN(3,NITER/2)
        CALL CONVGTHBG(NITER,RLXB,BGEPS)
C
        CALL TQCALC(1)
C---- Drive thrust from rotor for ISPEC=1, total thrust for ISPEC=2
        IF(ISPEC.EQ.1) THEN
          THR = TTOT 
        ELSE
          THR = TTOT + TDUCT
        ENDIF
C
C---- Check for disk type (actuator disk or bladed)
        IF(IRTYPE(NR).EQ.1) THEN
C---- Actuator disk
         BGAV  = 0.0
         BGMAX = BGAM(1,NR)
         BGMIN = BGMAX
         DO IR = 1, NRC
           DR = YRP(IR+1,NR)-YRP(IR,NR)
           DA = PI*(YRP(IR+1,NR)**2-YRP(IR,NR)**2)
           BGAV = BGAV + BGAM(IR,NR)
           BGMAX = MAX(BGMAX,BGAM(IR,NR))
           BGMIN = MIN(BGMIN,BGAM(IR,NR))
         END DO
         BGAV = BGAV / FLOAT(NRC)
         BGMAG = MAX(ABS(BGAV),BGMIN,BGMAX,0.1)
C---- Scale factor for BGAM from current value to get desired rotor thrust
         TSCL = 1.0
         IF(THR.NE.0.0)  TSCL = TSPEC/THR
C---- Check for rational relaxation factors based on BGAM changes
         DBGMAX = 0.0
         RLXB = RLX
         DO IR = 1, NRC
           DBG = (TSCL-1.0)*BGAM(IR,NR)
           DBGMAX = MAX(DBGMAX,ABS(DBG))
           IF(BGMAG.NE.0.0) THEN 
             FDBG = ABS(DBG)/BGMAG
             IF(FDBG*RLXB.GT.0.3) RLXB = 0.3/FDBG
           ENDIF
         END DO
C---- Update BGAM
         DO IR = 1, NRC
           BGFAC = (TSCL-1.0)*BGAM(IR,NR)
           BGAM(IR,NR) = BGAM(IR,NR) + RLXB * BGFAC
         END DO
         IF(TGAP.GT.0.0 .AND. OMEGA(NR).NE.0.0) BGAM(NRC,NR) = 0.0
C
        WRITE(*,110) ITR,RLXB,DBGMAX*RLXB
        IF(ABS(DBGMAX*RLXB).LE.0.001*BGMAG) THEN
          LCONV = .TRUE.
          GO TO 20
        ENDIF
        LGAMA = .FALSE.
C
       ELSEIF(IRTYPE(NR).EQ.2) THEN
C---- Bladed disk
        DTDALF = 0.0
        DO IR = 1, NRC
C---- theta velocity at blade (use only 1/2 of induced Vt from circulation)
         VTBG = BGAM(IR,NR)*PI2I/YRC(IR,NR)
         WTB  = VREL(3,IR,NR) - 0.5*VTBG
         WWB  = SQRT(VREL(1,IR,NR)**2 + WTB**2)
         PHIB = ATAN2(VREL(1,IR,NR),-WTB)
C
         XI = YRC(IR,NR)/RTIP(NR)
         ALF = BETAR(IR,NR) - PHIB
         REY = WWB*CHR(IR,NR)*RHO/RMU 
         SECSIG = BLDS*CHR(IR,NR)/(2.0*PI*YRC(IR,NR))
         SECSTAGR = 0.5*PI - BETAR(IR,NR)
         CALL GETCLCDCM(NR,IR,XI,ALF,WWB,REY,SECSIG,SECSTAGR,
     &                  CLB,CL_ALF,CL_W,
     &                  CLMAX,CLMIN,DCL_STALL,LSTALLR(IR,NR),
     &                  CDR(IR,NR),CD_ALF,CD_W,CD_REY,
     &                  CMOM,CM_AL,CM_W)
         CLR(IR,NR) = CLB
         CLALF(IR,NR)=CL_ALF
         ALFAR(IR,NR)=ALF
C
         IF(IR.EQ.NRC .AND. TGAP.GT.0.0 .AND. OMEGA(NR).NE.0.0) THEN
          CLR(IR,NR)  = 0.0
          CLALF(IR,NR)= 0.0
          ALFAR(IR,NR)= 0.0
         ENDIF
C
         DTDALF = DTDALF + TI_GAM(IR,NR)*0.5*WWB*CHR(IR,NR)*CL_ALF
        END DO
C---- Change blade pitch to drive thrust 
        RLXB = RLXF
        DBETA = 0.0
        IF(DTDALF.NE.0.0) DBETA = (TSPEC-THR) / DTDALF
C---- Limit DBETA changes in iteration to get desired rotor thrust
        IF(ABS(DBETA)*RLXB.GT.0.1) RLXB = 0.1/ABS(DBETA)
C---- update BETA by estimate
        DO IR = 1, NRC
           BETAR(IR,NR) = BETAR(IR,NR) + RLXB * DBETA
        END DO
C
        WRITE(*,100) ITR,RLXB,DBETA/DTR
        IF(ABS(DBETA).LT.0.001) THEN
          LCONV = .TRUE.
          GO TO 20
        ENDIF
        LGAMA = .FALSE.
C
       ENDIF
C
      END DO
      LCONV = .FALSE.
C
 100  FORMAT(I3,' RLX=',F8.5,' dBeta=',F9.5)
 110  FORMAT(I3,' RLX=',F8.5,' dBGmax=',F9.5)
C
C---- Final iterations to converge case
 20   CALL CONVGTHBG(NITER,RLXF,WXEPS)
C---- Update rotor velocities
      CALL UPDROTVEL
C
      RETURN
      END



      SUBROUTINE UPDROTVEL
C----------------------------------------------------------------
C     Update blade or disk velocities based on current solution
C     Velocities updated include:
C         induced        velocities  VIND
C         absolute frame velocities  VABS
C         relative frame velocities  VREL
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- get induced velocities upstream and downstream of disk
      DO N = 1, NROTOR
        CALL ROTORVABS(N,VIND(1,1,N))
      END DO
C
C---- set disk downstream velocities in absolute and relative frames
      CALL SETROTVEL 
C
C---- get area-averaged axial velocity over disk
      DO N = 1, NROTOR
       AINT  = 0.0
       VAINT = 0.0
       DO IR = 1, NRC
        DR = YRP(IR+1,N)-YRP(IR,N)
        DA = PI*(YRP(IR+1,N)**2-YRP(IR,N)**2)
        US = VABS(1,IR,N)
        AINT = AINT + DA
        VAINT = VAINT + US*DA 
       END DO        
       VAAVG(N) = VAINT/AINT
cc       write(*,*) 'n,vaavg ',n,vaavg(n)
      END DO
C
      RETURN
      END



      SUBROUTINE VABS2VREL(OMEG,YA,
     &                     VXA,VRA,VTA,
     &                     VXR,VRR,VTR,VTR_VTA,VTR_OMG)
C--------------------------------------------------------------
C     Calculates relative frame induced velocities from 
C     absolute frame velocities at radius YA with rotational 
C     speed OMEG. 
C--------------------------------------------------------------
C
      VXR = VXA
      VRR = VRA
C---- Blade relative velocity includes rotational speed and swirl effects 
      VTR = VTA - OMEG*YA
      VTR_VTA = 1.0
      VTR_OMG = -YA
C
      RETURN
      END



      SUBROUTINE ROTORVABS(N,VEL)
C--------------------------------------------------------------
C     Get absolute frame induced velocities downstream of rotor 
C     line center points 
C--------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION VEL(3,*)
C
C---- get velocities on rotor source line
      IEL = IELROTOR(N)
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      IG = IGROTOR(N)
C
      DO IC = IC1, IC2
        IR = IC - IC1 + 1
C---- Use mean surface velocity on rotor source panel
        VEL(1,IR) = QC(1,IC) - QINF
        VEL(2,IR) = QC(2,IC)
C---- Circumferential velocity downstream of rotor due to circulation
c        IF(IG.EQ.1) THEN
         BGAVE = BGAMG(IG,IR)
c        ELSE
c        BGAVE = 0.5*(BGAMG(IG-1,IR)+BGAMG(IG,IR))
c        ENDIF
        CIRC      = BGAVE
        VEL(3,IR) = CIRC*PI2I/YC(IC)
      END DO
C
 98   FORMAT(A,I5,6(1X,F10.6))
 99   FORMAT(A,2I5,6(1X,F10.6))
 97   FORMAT(A,3I5,6(1X,F10.6))
C
      RETURN
      END



      SUBROUTINE GETVELABS(NF,XF,YF,VEL)
C--------------------------------------------------------------
C     Get absolute frame induced velocities at NF points XF,YF
C     by potential flow survey 
C--------------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION XF(*),YF(*)
      DIMENSION VEL(3,*)
C
C---- local arrays for calling QFCALC
      PARAMETER (NFX=IRX)
      DIMENSION QF(2,NFX),
     &          QF_GAM(2,IPX,NFX),
     &          QF_SIG(2,IPX,NFX),
     &          QF_GTH(2,IPX,NFX)
      DIMENSION IPFO(NFX), IPFP(NFX), IFTYPE(NFX)
C
C---- get velocities on specified points
C---- assume the pointa are not on a panel
C
      DO IR = 1, NF
        IPFO(IR) = 0
        IPFP(IR) = 0
        IFTYPE(IR) = 0
      END DO
C
C------ evaluate velocity components at points
      CALL QFCALC(1,NF, XF,YF, IPFO,IPFP, IFTYPE,
     &            QF,QF_GAM,QF_SIG,QF_GTH)
C
      DO IF = 1, NF
        VEL(1,IF) = QF(1,IF) 
        VEL(2,IF) = QF(2,IF)
        VEL(3,IF) = 0.0
      END DO
C
      RETURN
      END

 

      SUBROUTINE PRTVEL(LU,LINE,LIND,LABS,LREL,NR)
C----------------------------------------------------------
C     Print out velocities just downstream of rotor 
C     Prints either:
C        absolute frame induced velocities
C        absolute frame total velocities
C        blade relative frame total velocities
C
C     LU is logical unit for output
C     LINE is a title message for table
C     Print flags LIND,LREL,LABS control what gets printed
C     NR is disk # for rotor
C
C     Takes VIND from current solution for absolute frame
C     Takes VABS from current solution for absolute frame
C     Takes VREL from current solution for relative frame
C----------------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*) LINE
      LOGICAL LIND, LABS, LREL
C
C---- Print velocity data
      WRITE(LU,20) LINE
C
      RPM = 30.0*OMEGA(NR)/PI
      WRITE(LU,25) QINF,QREF,OMEGA(NR),RPM
C
C---- Absolute frame induced velocities
      IF(LIND) THEN
        WRITE(LU,30)
        DO I = 1, NRC
          YY = YRC(I,NR)
          VXA = VIND(1,I,NR)
          VRA = VIND(2,I,NR)
          VMA = SQRT(VXA**2 + VRA**2)
          VTA = VIND(3,I,NR)
          VVA = SQRT(VMA**2 + VTA**2)
cc          ANG = ATAN2(VTA,VMA)
          ANG = ATAN2(VTA,VXA)
          WRITE(LU,50) YY,VXA,VRA,VMA,VTA,VVA,ANG/DTR
        END DO
      ENDIF
C
C---- Absolute frame velocities
      IF(LABS) THEN
        WRITE(LU,35)
        DO I = 1, NRC
          YY = YRC(I,NR)
          VXA = VABS(1,I,NR)
          VRA = VABS(2,I,NR)
          VMA = SQRT(VXA**2 + VRA**2)
          VTA = VABS(3,I,NR)
          VVA = SQRT(VMA**2 + VTA**2)
cc          ANG = ATAN2(VTA,VMA)
          ANG = ATAN2(VTA,VXA)
          WRITE(LU,50) YY,VXA,VRA,VMA,VTA,VVA,ANG/DTR
        END DO
      ENDIF
C
C---- Relative frame velocities
      IF(LREL) THEN
        WRITE(LU,40)
        DO I = 1, NRC
          YY = YRC(I,NR)
          WXR = VREL(1,I,NR)
          WRR = VREL(2,I,NR)
          WMR = SQRT(WXR**2 + WRR**2)
C---- Blade relative velocity includes rotational speed
          WTR = VREL(3,I,NR)
          WWR = SQRT(WMR**2 + WTR**2)
          IF(WTR.NE.0.0) THEN
cc            PHIR = ATAN2(WMR,-WTR)
            PHIR = ATAN2(WXR,-WTR)
          ELSE
            PHIR = 0.5*PI
          ENDIF
          ANG = PHIR
          WRITE(LU,50) YY,WXR,WRR,WMR,WTR,WWR,ANG/DTR
        END DO
      ENDIF
C
C---- Relative frame velocities on blade lifting line
C     this assumes that the radial component of velocity is parallel 
C     to the blade span and is ignored in velocity magnitude and angle 
C     for the blade.  Radial components are printed however.
      IF(LREL) THEN
        WRITE(LU,60)
        DO I = 1, NRC
          YY = YRC(I,NR)
          WXR = VREL(1,I,NR)
          WRR = VREL(2,I,NR)
          WMR = SQRT(WXR**2 + WRR**2)
C---- Blade relative velocity includes rotational speed, 1/2 induced swirl
C---- Angle measured from plane of rotation
          VTBG = BGAM(I,NR)*PI2I/YRC(I,NR)
          WTR = VREL(3,I,NR) - 0.5*VTBG
          WWR = SQRT(WMR**2 + WTR**2)
          IF(WTR.NE.0.0) THEN
            PHIB = ATAN2(WXR,-WTR)
          ELSE
            PHIB = 0.5*PI
          ENDIF
          ANG = PHIB
          WRITE(LU,50) YY,WXR,WRR,WMR,WTR,WWR,ANG/DTR
        END DO
      ENDIF
C
 20   FORMAT(/,A)
 25   FORMAT(' QINF  =',F12.4,
     &      /' QREF  =',F12.4,
     &      /' OMEGA =',F12.4,
     &      /' RPM   =',F12.4)
 30   FORMAT(/'Induced vel, flow angles in absolute frame',
     &        ' (downstream of disk)',
     &      /'     r          Vxi        Vri',
     &   '        Vmi        Vti        Vi     Swirl(deg)')
 35   FORMAT(/'Velocities, flow angles in absolute frame',
     &        ' (downstream of disk)',
     &      /'     r          Vx         Vr',
     &   '         Vm         Vt         V     Swirl(deg)')
 40   FORMAT(/'Velocities, flow angles relative to blade frame',
     &        ' (downstream of disk)',
     &      /'     r          Wx         Wr',
     &   '         Wm         Wt         W       Phi(deg)')
 60   FORMAT(/'Velocities in blade frame,',
     &        ' on blade lifting line',
     &       /'flow angle from plane of rotation',
     &      /'     r          Wx         Wr',
     &   '         Wm         Wt         W       Phi(deg)')
C                  12345678901123456789011234567890112345678901
 50   FORMAT(7G11.4)
C
      RETURN
      END



      SUBROUTINE SHOWDUCT(LU)
C-------------------------------------------------------
C     Displays duct geometric data 
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      NR = 1
      WRITE(LU,100) NAME
      WRITE(LU,125) QINF,QREF, VSO, RHO, RMU
      WRITE(LU,150) ANAME,NBEL
C
      DO N = 1, NBEL
        WRITE(LU,300) XBLE(N),YBLE(N),
     &               XBTE(N),YBTE(N),
     &               XBREFE(N), YBREFE(N),
     &               VOLUMV(N), ASURFV(N),
     &               XBCENV(N), YBCENV(N),
     &               RGYRXV(N), RGYRYV(N)
      END DO
C
c     &  AREA2DA(NEX),XBCEN2DA(NEX),YBCEN2DA(NEX),
c     &  EIXX2DA(NEX),EIYY2DA(NEX), EIXY2DA(NEX),
c     &  AREA2DT(NEX),XBCEN2DT(NEX),YBCEN2DT(NEX),
c     &  EIXX2DT(NEX),EIYY2DT(NEX), EIXY2DT(NEX),
c
c     &  VOLUMV(NEX), ASURFV(NEX),  XBCENV(NEX), YBCENV(NEX),
c     &  RGYRXV(NEX), RGYRYV(NEX),
c     &  VOLUMVT(NEX),ASURFVT(NEX),XBCENVT(NEX),YBCENVT(NEX),
c     &  RGYRXVT(NEX), RGYRYVT(NEX),
C     
  100 FORMAT(/,' DFDC Case: ',A30,/,1X,55('-'))
  125 FORMAT(  '  Qinf  ',F9.4,'     Qref  ',F9.4,
     &       /,'  Speed of sound (m/s)     ', F10.3,
     &       /,'  Air density   (kg/m^3)   ', F10.5,
     &       /,'  Air viscosity (kg/m-s)    ', E11.4 )
c
  150 FORMAT(/,' Geometry name: ',A30,
     &       /,' Current duct geometry with', I2, ' elements:')
C
  300 FORMAT(/,'  Xle   ',F12.5,'     Rle   ',F12.5,
     &       /,'  Xte   ',F12.5,'     Rte   ',F12.5,
     &       /,'  Xref  ',F12.5,'     Rref  ',F12.5,
     &       /,'  Vol   ',F12.6,'     Asurf ',F12.6,
     &       /,'  Xcent ',F12.6,'     Rcent ',F12.6,
     &       /,'  rGYRx ',F12.6,'     rGYRr ',F12.6)
C
      RETURN
      END ! SHOWDUCT


      SUBROUTINE SHOWACTDSK(LU)
C-------------------------------------------------------
C     Displays parameters on actuator disk
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DO N = 1, NROTOR
C
       IF(IRTYPE(N).EQ.1) THEN
C        WRITE(LU,100) NAME
        RPM = 30.0*OMEGA(N)/PI
C
        WRITE(LU,120) N,OMEGA(N),RPM,RHUB(N),RTIP(N),ADISK(N)
        WRITE(LU,200)
C 
        DO IR = 1, NRC
          WRITE(LU,210) YRC(IR,N),YRC(IR,N)/RTIP(N),BGAM(IR,N)
        END DO
       ENDIF
C
      END DO
C     
  100 FORMAT(/A)
C
  120 FORMAT(/' Current Actuator Disk at Disk',I3,
     &       / '  Omega  ',F12.4, '    Rpm    ',F11.2, 
     &       / '  Rhub   ',F12.5, '    Rtip   ',F11.5,
     &       / '  Aswept ',F12.5)

  200 FORMAT(  '     r        r/R       B*Gamma')
  210 FORMAT(1X,7G11.4)
C
      RETURN
      END ! SHOWACTDSK


      SUBROUTINE SHOWBLADE(LU)
C-------------------------------------------------------
C     Displays chord, blade angle and solidity distribution 
C     on rotor blade
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DO N = 1, NROTOR
C
       IF(IRTYPE(N).EQ.2) THEN
        RPM = 30.0*OMEGA(N)/PI
C
        WRITE(LU,120) N,OMEGA(N),RPM,RHUB(N),
     &                RTIP(N),NRBLD(N),ADISK(N)
        BLDS = FLOAT(NRBLD(N))
        WRITE(LU,200)
C
        DO IR = 1, NRC
          SIGROT = BLDS*CHR(IR,N)/(2.0*PI*YRC(IR,N))
          WRITE(LU,210) YRC(IR,N),YRC(IR,N)/RTIP(N),
     &                 CHR(IR,N),BETAR(IR,N)/DTR,SIGROT
        END DO
       ENDIF
C
      END DO
C     
  100 FORMAT(/A)
C
  120 FORMAT(/ ' Current Rotor at Disk',I3,
     &       / '  Omega  ',F12.4, '    Rpm    ',F11.2, 
     &       / '  Rhub   ',F12.5, '    Rtip   ',F11.5,
     &       / '  #blades',I6, 6X,'    Aswept ',F11.5)
  200 FORMAT(  '     r         r/R        Ch      Beta(deg)   Solidity')
  210 FORMAT(1X,7G11.4)
C

      RETURN
      END ! SHOWBLADE


      SUBROUTINE SHOWDRAGOBJ(ND,LU)
C-------------------------------------------------------
C     Displays parameters on drag objects
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(NDOBJ.LE.0) THEN
        WRITE(LU,*)
        WRITE(LU,*) 'No drag objects defined'
        RETURN
      ENDIF
C
      IF(ND.LE.0) THEN 
        ND1 = 1
        ND2 = NDOBJ
        WRITE(LU,120) NDOBJ 
       ELSE
        ND1 = ND
        ND2 = ND
      ENDIF
C
      DO N = ND1, ND2
        WRITE(LU,130) N
        WRITE(LU,200)
C
        DO ID = 1, NDDEF(N)
          WRITE(LU,210) ID,XDDEF(ID,N),YDDEF(ID,N),CDADEF(ID,N)
        END DO
C
      END DO
C     
  100 FORMAT(/A)
  110 FORMAT(' QINF  =',F12.4,
     &      /' QREF  =',F12.4,
     &      /' OMEGA =',F12.4,
     &      /' RPM   =',F12.4)
C
  120 FORMAT(/,I2,' Drag Area objects defined:',/)
  130 FORMAT(     ' Drag Area', I3)
  200 FORMAT(    '    i        x          r        CDave')
c             a1234567890112345678901123456789011234567890112345678901    
  210 FORMAT(I5,F11.4,F11.4,3X,G11.4)
C
      RETURN
      END ! SHOWDRAGOBJ



      SUBROUTINE DESBLADE(NDES)
C-------------------------------------------------------
C     Sets chord CH and blade angle distribution for current
C     blade circulation and CLDES values
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      NR = NDES
C---- use circulation and local velocity along with CL to define blade
C     sets chord and blade angle using CLDES values on blade stations
C
      CALL SETROTVEL  ! needed for blade blockage toggle
C
C---- update design CL arrays
       IF(OMEGA(NR).GT.0.0)THEN
          DO I=1,NRC
            CLDES(I) = CLPOS(I)
          ENDDO
        ELSE
          DO I=1,NRC
            CLDES(I) = CLNEG(I)
          ENDDO
        ENDIF 
C
      WRITE(*,*)
      WRITE(*,1200)
      IF(LBLBL) THEN
        WRITE(*,1300) NR, NRBLD(NR)
      ELSE
        WRITE(*,1310) NR, NRBLD(NR)
      ENDIF
      WRITE(*,1200)
      WRITE(*,1400)
      BLDS = FLOAT(NRBLD(NR))
C
      DO IR = 1, NRC
C---- use theta velocity at blade lifting line
        VTBG = BGAM(IR,NR)*PI2I/YRC(IR,NR)
        WTB  = VREL(3,IR,NR) - 0.5*VTBG
        WWB  = SQRT(VREL(1,IR,NR)**2 + WTB**2)
        PHIB = ATAN2(VREL(1,IR,NR),-WTB)
C
        IF(PHIB.GT.0.5*PI) THEN
         IF(OMEGA(NR).GT.0.0) THEN
          WRITE(*,*) 'Warning Wtheta reversed ',IR,YRC(IR,NR),PHIB/DTR
         ENDIF
        ENDIF
C
C---- chord set by circulation, velocity and design CL
        CHDES(IR) = ABS(2.0*BGAM(IR,NR)/WWB/CLDES(IR) / BLDS)
        IF(IR.EQ.NRC .AND. TGAP.GT.0.0) THEN
          CHDES(IR) = CHDES(IR-1)
          PHIB = PHIBOLD
        ENDIF
C
C---- estimate blade alfa using 2*pi lift slope initially
        ALF = CLDES(IR)/(2.0*PI)
C---- For case with design BGAM<0, switch sign of ALF
        BETADES(IR) = PHIB + SIGN(ALF,BGAM(IR,NR))
C---- blade angle set by alfa(CLdes) and flow angle
        XI = YRC(IR,NR)/RTIP(NR)
        SECSIG = BLDS*CHDES(IR)/(2.0*PI*YRC(IR,NR))
        SECSTAGR = 0.5*PI - BETADES(IR)
        CALL GETALF(NR,IR,XI,SECSIG,SECSTAGR,
     &              CLDES(IR),WWB,ALF,ALF_CL,ALF_W,LSTALLR(IR,NR))
C
C---- iterate once to set correct stagger angle for multiplane correction
        BETADES(IR) = PHIB + SIGN(ALF,BGAM(IR,NR))
        SECSTAGR = 0.5*PI - BETADES(IR)
        CALL GETALF(NR,IR,XI,SECSIG,SECSTAGR,
     &              CLDES(IR),WWB,ALF,ALF_CL,ALF_W,LSTALLR(IR,NR))
C
C---- iterate once to set correct stagger angle for multiplane correction
        BETADES(IR) = PHIB + SIGN(ALF,BGAM(IR,NR))
        SECSTAGR = 0.5*PI - BETADES(IR)
        CALL GETALF(NR,IR,XI,SECSIG,SECSTAGR,
     &              CLDES(IR),WWB,ALF,ALF_CL,ALF_W,LSTALLR(IR,NR))
C
        CL_ALF = 1.0/ALF_CL
        BETADES(IR) = PHIB + SIGN(ALF,BGAM(IR,NR))
        PHIBOLD = PHIB
C
C---- NOTE: this should only be printed up to NRC-1 for a tip gap case
        IF(TGAP.EQ.0.0 .OR. (TGAP.GT.0.0 .AND. IR.LT.NRC)) THEN
         WRITE(*,1410) XI,YRC(IR,NR),
     &                 CHDES(IR),BETADES(IR)/DTR,CLDES(IR),
     &                 WWB,BGAM(IR,NR),CL_ALF,ALF/DTR
        ENDIF
C
        CHR(IR,NR)   = CHDES(IR)
        BETAR(IR,NR) = BETADES(IR)
        CLR(IR,NR)   = SIGN(CLDES(IR),BGAM(IR,NR))
        ALFAR(IR,NR) = ALF
C
      END DO
C
      WRITE(*,1200)
C
C---- change type for disk so blades are defined
      IRTYPE(NR) = 2
      LBLDEF = .TRUE.
C
C---- update blade blockage data
      IF(LBBLOFT(NR)) CALL BBUPDATE(NR)
C
 1200 FORMAT(1X,78('-'))
 1300 FORMAT(' Disk',I2,' designed with',I2,' blades',
     &     21X,'Corrected for blade blockage')
 1310 FORMAT(' Disk',I2,' designed with',I2,' blades',
     &     21X,'No blade blockage correction')
 1400 FORMAT('     r/R       r       Ch    Beta0deg',
c              1234567891234578901234567890123456789
     &        '    CL       W     BGam ',
c              123456789123456789123456789
     &        '   dCLdA  Alfadeg')
c              123456789123456789
c 1410 FORMAT(F9.4,2F10.4,6F9.3)
 1410 FORMAT(3F9.4,F9.3,F8.3,F9.3,2F8.3,F9.3)
C
      RETURN
      END ! DESBLADE



      SUBROUTINE SETDRGOBJSRC
C-------------------------------------------------------
C     Sets source strength for rotor drag and drag objects
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION VDM(IRX), T1(IRX), T1S(IRX)
C
C---- Set sources for drag objects
      IF(NDOBJ.LE.0) THEN
cc        WRITE(*,*) 'Drag objects not defined'
        RETURN
C
      ELSE
C---- step through all defined drag objects
       DO N = 1, NDOBJ
C
        IF(IELDRGOBJ(N).LE.0) THEN
          WRITE(*,*) 'Source element not defined for drag object ',N
          STOP
        ENDIF
C
C---- if drag objects are turned off clear the source strengths
        IF(.NOT.LDRGOBJ) THEN
         IEL = IELDRGOBJ(N)
         IP1 = IPFRST(IEL)
         IP2 = IPLAST(IEL)
         DO IP = IP1, IP2
           SIGVSP(IP) = 0.0
         END DO
C
        ELSE
C
C---- Set up spline arrays for drag object
         IF(NDDEF(N).GE.2) THEN
          DO I = 1, NDDEF(N)
            T1(I) = CDADEF(I,N)
          END DO
C---- Spline drag definition array
          CALL SEGSPL(T1,T1S,YDDEF(1,N),NDDEF(N))
C
          IEL = IELDRGOBJ(N)
          IC1 = ICFRST(IEL)
          IC2 = ICLAST(IEL)
          IP1 = IPFRST(IEL)
          IP2 = IPLAST(IEL)
C
          DO IC = IC1, IC2
            VXX = QC(1,IC)
            VRR = QC(2,IC)
            ID = IC - IC1 + 1
            VDM(ID) = SQRT(VXX**2 + VRR**2)
          END DO
C
          DO IP = IP1, IP2
            YY = YP(IP)
            CDAV = SEVAL(YY,T1,T1S,YDDEF(1,N),NDDEF(N))
            ID = IP - IP1 + 1
            IF(IP.EQ.IP1) THEN
             SIGVSP(IP) = 0.5*VDM(ID)*CDAV
            ELSEIF(IP.EQ.IP2) THEN
             SIGVSP(IP) = 0.5*VDM(ID-1)*CDAV
            ELSE
             SIGVSP(IP) = 0.25*(VDM(ID)+VDM(ID-1))*CDAV 
            ENDIF
cc          WRITE(*,99) 'IP,CDAV,SIGVSP ',IP,CDAV,SIGVSP(IP)
          END DO
         ENDIF
C
        ENDIF
C
       END DO
      ENDIF
C     
 99   FORMAT(A,I4,5(1X,F11.5))
C
      RETURN
      END ! SETDRGOBJSRC
  


      SUBROUTINE SETROTORSRC
C-------------------------------------------------------
C     Sets source strength for rotor profile drag 
C-------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION VDM(IRX), T1(IRX), T1S(IRX)
C
      DO N = 1, NROTOR
C---- Set sources for rotor drag 
       IF(IELROTOR(N).LE.0) THEN
cc        WRITE(*,*) 'Rotor source line not defined'
cc         RETURN
       ELSE
C
        IEL = IELROTOR(N)
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        IP1 = IPFRST(IEL)
        IP2 = IPLAST(IEL)
        BLDS = FLOAT(NRBLD(N))
C
        DO IC = IC1, IC2
         VXX = QC(1,IC)
         VRR = QC(2,IC)
         VMM = SQRT(VXX**2 + VRR**2)
         IR = IC - IC1 + 1
         VTT = 0.5*BGAM(IR,N)*PI2I/YC(IC) - OMEGA(N)*YC(IC)
         VDM(IR) = SQRT(VMM**2 + VTT**2)
        END DO
C
        DO IP = IP1, IP2
         IR = IP - IP1 + 1
         IF(IP.EQ.IP1) THEN
          SIGVSP(IP) = 0.5*BLDS*PI2I * VDM(IR)*CHR(IR,N)*CDR(IR,N)
cc         WRITE(*,99) 'IP,W,CD,SIGVSP ',IP,VDM(IR),CDR(IR),SIGVSP(IP)
         ELSEIF(IP.EQ.IP2) THEN
C---- NOTE: should set tip source to 0.0 for tip gap (no blade defined here)
          SIGVSP(IP) = 0.5*BLDS*PI2I * VDM(IR-1)*CHR(IR-1,N)*CDR(IR-1,N)
cc         WRITE(*,99) 'IP,W,CD,SIGVSP ',IP,VDM(IR-1),CDR(IR-1),SIGVSP(IP)
         ELSE
          VAVE  = 0.5*(VDM(IR)+VDM(IR-1))
          CDAVE = 0.5*(CDR(IR,N)+CDR(IR-1,N))
          CHAVE = 0.5*(CHR(IR,N)+CHR(IR-1,N))
          SIGVSP(IP) = 0.5*BLDS*PI2I * VAVE*CHAVE*CDAVE
cc         WRITE(*,99) 'IP,W,CD,SIGVSP ',IP,VAVE,CDAVE,SIGVSP(IP)
         ENDIF
        END DO
       ENDIF
C
      END DO
C     
 99   FORMAT(A,I4,5(1X,F11.5))
C
      RETURN
      END ! SETROTORSRC



      SUBROUTINE VMAVGINIT(VAXIAL)
C---------------------------------------------------------
C     Initializes VMAVG using axial velocity estimate
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) THEN
        WRITE(*,*)      'VMAVGINIT called with VA=',VAXIAL
        WRITE(LUNDBG,*) 'VMAVGINIT called with VA=',VAXIAL
      ENDIF
C
C---- Set VMAVG on all wake element points 
      DO IR = 1, NRP
         IEL = IR2IEL(IR)
cc      IC1 = ICFRST(IEL)
cc      IC2 = ICLAST(IEL)
         IP1 = IPFRST(IEL)
         IP2 = IPLAST(IEL)
         DO IP = IP1, IP2
           IF(IR.EQ.1) THEN
             VMAVG(IP) = 0.5*(VAXIAL+QINF)
           ELSEIF(IR.EQ.NRP) THEN
             VMAVG(IP) = 0.5*(VAXIAL+QINF)
           ELSE
             VMAVG(IP) = VAXIAL
           ENDIF
         END DO 
      END DO
C
C---- Set VMAVG on CB body wake points
c      IEL = 1
c      IC1 = ICFRST(IEL)
c      IC2 = ICLAST(IEL)
c      DO IC = IC1, IC2
c        IP1 = IPCO(IC)
c        IP2 = IPCP(IC)
c        IF(IP2IR(IP1).NE.0 .AND. IP2IR(IP2).NE.0) THEN
c          VMAVG(IP1) = 0.5*(VAXIAL+QINF)
c          VMAVG(IP2) = 0.5*(VAXIAL+QINF)
c        ENDIF
c      END DO
C
C---- Set VMAVG on duct body wake points
c      IEL = 2
c      IC1 = ICFRST(IEL)
c      IC2 = ICLAST(IEL)
c      DO IC = IC2, IC1, -1
c        IP1 = IPCO(IC)
c        IP2 = IPCP(IC)
c        IF(IP2IR(IP1).NE.0 .AND. IP2IR(IP2).NE.0) THEN
c          VMAVG(IP1) = 0.5*(VAXIAL+QINF)
c          VMAVG(IP2) = 0.5*(VAXIAL+QINF)
c        ENDIF
c      END DO
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*)  'At end of VMAVGINIT'
        DO IEL = 1, NEL
c         IF(NETYPE(IEL).EQ.7) THEN
           IC1 = ICFRST(IEL)
           IC2 = ICLAST(IEL)
           WRITE(LUNDBG,*)  IEL,IC1,IC2
           DO IC = IC1, IC2
             IP1 = IPCO(IC)
             IP2 = IPCP(IC)
             WRITE(LUNDBG,*)  'IP1,VMavg ',IP1,VMAVG(IP1)
             WRITE(LUNDBG,*)  'IP2,VMavg ',IP2,VMAVG(IP2)
           END DO
c        ENDIF
        END DO
      ENDIF
 20   FORMAT(A,I5,5(1X,F10.4))
C
      LVMAV = .TRUE.
C
      RETURN
      END ! VMAVGINIT


      SUBROUTINE VMAVGCALC
C---------------------------------------------------------
C     Calculates VMAVG at rotor wake points using current
C     center point velocities QC(1,.), QC(2,.)
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) THEN
        WRITE(*,*)      'Entering VMAVGCALC'
        WRITE(LUNDBG,*) 'Entering VMAVGCALC'
      ENDIF
C
C---- Set VMAVG on all wake points, first and last set from centers
C     intermediate points get average from nearby panel centers
      DO IEL = 1, NEL
        IF(NETYPE(IEL).EQ.7) THEN
         IC1 = ICFRST(IEL)
         IC2 = ICLAST(IEL)
         DO IC = IC1, IC2
           IP1 = IPCO(IC)
           IP2 = IPCP(IC)
c           IF(IEL.NE.IR2IEL(1)) THEN
             QCX = QC(1,IC)
             QCY = QC(2,IC)
c           ELSE
c             QCX = 0.5*QCL(1,IC)
c             QCY = 0.5*QCL(2,IC)
c           ENDIF
           VMAV = SQRT(QCX**2 + QCY**2)
C---- limiter to keep VMAV positive (to fraction of QREF)
           VMAV = MAX(0.1*QREF,VMAV)
           IF(IC.EQ.IC1) THEN
            VMAVG(IP1) = VMAV
           ELSE 
            VMAVG(IP1) = 0.5*(VMAVG(IP1) + VMAV)
           ENDIF
           VMAVG(IP2) = VMAV
         END DO
        ENDIF
      END DO
C
C---- Set VMAVG on CB body vortex wake points
c      IEL = 1
c      IR  = 1
c      IC1 = ICFRST(IEL)
c      IC2 = ICLAST(IEL)
c      DO IC = IC1, IC2
c        IP1 = IPCO(IC)
c        IP2 = IPCP(IC)
c        IF(IP2IR(IP1).EQ.IR .AND. IP2IR(IP2).EQ.IR) THEN
cc         write(77,*) 'iel,ic,ir ',iel,ic,ir
c         QCX = QC(1,IC)
c         QCY = QC(2,IC)
c         VMAV = SQRT(QCX**2 + QCY**2)
c         IF(IC.EQ.IC1) THEN
c          VMAVG(IP1) = VMAV
c         ELSE 
c          VMAVG(IP1) = 0.5*VMAVG(IP1) + 0.5*VMAV
c         ENDIF
c         VMAVG(IP2) = VMAV
cc         write(77,*) 'ip1,vmavg ',ip1,vmavg(ip1)
cc         write(77,*) 'ip2,vmavg ',ip2,vmavg(ip2)
c        ENDIF
c      END DO
C---- Average velocities on CB TE and wake upstream point
c      IPB = IPFRST(IEL)
c      IPW = IPFRST(IR2IEL(IR))
c      VMAV = 0.5*(VMAVG(IPB) + VMAVG(IPW))
c      VMAVG(IPB) = VMAV
c      VMAVG(IPW) = VMAV
C
C---- Set VMAVG on duct body vortex wake points
c      IEL = 2
c      IR  = NRP
c      IC1 = ICFRST(IEL)
c      IC2 = ICLAST(IEL)
c      DO IC = IC2, IC1, -1
c        IP1 = IPCO(IC)
c        IP2 = IPCP(IC)
c        IF(IP2IR(IP1).EQ.IR .AND. IP2IR(IP2).EQ.IR) THEN
c         QCX = QC(1,IC)
c         QCY = QC(2,IC)
c         VMAV = SQRT(QCX**2 + QCY**2)
c         IF(IC.EQ.IC2) THEN
c          VMAVG(IP2) = VMAV
c         ELSE 
c          VMAVG(IP2) = 0.5*VMAVG(IP1) + 0.5*VMAV
c         ENDIF
c         VMAVG(IP1) = VMAV 
c        ENDIF
c      END DO
C---- Average velocities on duct TE and wake upstream point
c      IPB = IPLAST(IEL)
c      IPW = IPFRST(IR2IEL(IR))
c      VMAV = 0.5*(VMAVG(IPB) + VMAVG(IPW))
c      VMAVG(IPB) = VMAV
c      VMAVG(IPW) = VMAV
C

      IF(LDBG) THEN
        WRITE(LUNDBG,*)  'At end of VMAVGCALC'
        DO IEL = 1, NEL
c          IF(NETYPE(IEL).EQ.7) THEN
            IC1 = ICFRST(IEL)
            IC2 = ICLAST(IEL)
            WRITE(LUNDBG,*) 'IEL,IC1,IC2,IPFRST(IEL),IPLAST(IEL)'
            WRITE(LUNDBG,*)  IEL,IC1,IC2,IPFRST(IEL),IPLAST(IEL)
            DO IC = IC1, IC2
             IP1 = IPCO(IC)
             IP2 = IPCP(IC)
             WRITE(LUNDBG,*)  'IP1,VMavg ',IP1,VMAVG(IP1)
             WRITE(LUNDBG,*)  'IP2,VMavg ',IP2,VMAVG(IP2)
            END DO
c          ENDIF
        END DO
      ENDIF
 20   FORMAT(A,I5,5(1X,F10.4))
C
cc      LVMAV = .TRUE.
C
      RETURN
      END ! VMAVGCALC




      SUBROUTINE GTHCALC(GAMTH)
C---------------------------------------------------------
C     Calculates GTH at rotor wake points using zero 
C     pressure jump relation
C   Note: this formulation sets GTH on endpoints using 
C         center velocity (not averaged VM at endpoints)
C
C---------------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION GAMTH(IPX)
C
      IF(LDBG) WRITE(*,*) 'Entering GTHCALC'
C
      DO IP = 1, NPTOT
        GAMTH(IP) = 0.0
      END DO
C
C---- Set GAMTH on intermediate wakes
      DO IR = 2, NRP-1
        IEL  = IR2IEL(IR)
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        DO IC = IC1, IC2
          IP1 = IPCO(IC)
          IP2 = IPCP(IC)
          IG = IC2IG(IC)
          DH =                DHG(IG,IR)      - DHG(IG,IR-1)
          DG = 0.5*PI2I**2 * (BGAMG(IG,IR)**2 - BGAMG(IG,IR-1)**2)
          IF(VMAVG(IP1).NE.0.0) THEN
            GAMTH(IP1) = (DH - DG/(YP(IP1)**2))/VMAVG(IP1)
          ELSE
c            WRITE(*,*) 'Zero VMavg on IP1 =',IP1,VMAVG(IP1)
            GAMTH(IP1) = 0.
          ENDIF
          IF(VMAVG(IP2).NE.0.0) THEN
            GAMTH(IP2) = (DH - DG/(YP(IP2)**2))/VMAVG(IP2)
          ELSE
c            WRITE(*,*) 'Zero VMavg on IP2 =',IP2,VMAVG(IP2)
            GAMTH(IP2) = 0.
          ENDIF
        END DO
      END DO
C
C---- Set GAMTH on CB wake
      IR  = 1
      IEL  = IR2IEL(IR)
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      DO IC = IC1, IC2
        IP1 = IPCO(IC)
        IP2 = IPCP(IC)
        IG = IC2IG(IC)
        DH =                0.0
cc        DH =                DHG(IG,IR)
        DG = 0.5*PI2I**2 * (BGAMG(IG,IR)**2)
        IF(VMAVG(IP1).NE.0.0) THEN
          GAMTH(IP1) = (DH - DG/(YP(IP1)**2))/VMAVG(IP1)
        ELSE
          WRITE(*,*) 'Zero VMavg on CB wake IP1 =',IP1,VMAVG(IP1)
          GAMTH(IP1) = 0.
        ENDIF
        IF(VMAVG(IP2).NE.0.0) THEN
          GAMTH(IP2) = (DH - DG/(YP(IP2)**2))/VMAVG(IP2)
        ELSE
          WRITE(*,*) 'Zero VMavg on CB wake IP2 =',IP2,VMAVG(IP2)
          GAMTH(IP2) = 0.
        ENDIF
      END DO
C
C---- Set GTH on CB to value at first wake point
      IR  = 1
      IEL  = IR2IEL(IR)
      IP1 = IPFRST(IEL)
      GTH1 = GAMTH(IP1)
      Y1 = YP(IP1)
      IEL = 1
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      DO IC = IC1, IC2
        IP1 = IPCO(IC)
        IP2 = IPCP(IC)
        IF(IP2IR(IP1).EQ.IR .AND. IP2IR(IP2).EQ.IR) THEN
          GAMTH(IP1) = 0.0
          GAMTH(IP2) = 0.0
C---- taper off the GTH on CB from GTH at TE to 0.0 at rotor
          FRAC = 1.0 - FLOAT(IP1-IPFRST(IEL))/FLOAT(IPROTCB(1))
          GAMTH(IP1) = GTH1*FRAC
          FRAC = 1.0 - FLOAT(IP2-IPFRST(IEL))/FLOAT(IPROTCB(1))
          GAMTH(IP2) = GTH1*FRAC
        ENDIF
C
      END DO
C
C---- Set GAMTH on DUCT wake
      IR  = NRP
      IEL  = IR2IEL(IR)
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      DO IC = IC1, IC2
        IP1 = IPCO(IC)
        IP2 = IPCP(IC)
        IG = IC2IG(IC)
        DH =                                - DHG(IG,IR-1)
        DG = 0.5*PI2I**2 * (                - BGAMG(IG,IR-1)**2)
        IF(VMAVG(IP1).NE.0.0) THEN
          GAMTH(IP1) = (DH - DG/(YP(IP1)**2))/VMAVG(IP1)
        ELSE
          WRITE(*,*) 'Zero VMavg on duct wake IP1 =',IP1,VMAVG(IP1)
          GAMTH(IP1) = 0.
        ENDIF
        IF(VMAVG(IP2).NE.0.0) THEN
          GAMTH(IP2) = (DH - DG/(YP(IP2)**2))/VMAVG(IP2)
        ELSE
          WRITE(*,*) 'Zero VMavg on duct wake IP2 =',IP2,VMAVG(IP2)
          GAMTH(IP2) = 0.
        ENDIF
      END DO
C
C---- Set GAMTH on DUCT from first wake GAMTH
      IR  = NRP
      IEL  = IR2IEL(IR)
      IP1 = IPFRST(IEL)
      GTH1 = GAMTH(IP1)
      Y1 = YP(IP1)
      IEL = 2
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      DO IC = IC1, IC2
        IP1 = IPCO(IC)
        IP2 = IPCP(IC)
        IF(IP2IR(IP1).EQ.IR .AND. IP2IR(IP2).EQ.IR) THEN
c          GAMTH(IP1) = 0.0
c          GAMTH(IP2) = 0.0
c          GAMTH(IP1) = GTH1
c          GAMTH(IP2) = GTH1
C---- taper off the GTH on duct from GTH at TE to 0.0 at rotor
          FRAC = 1.0 - FLOAT(IPLAST(IEL)-IP1) /
     &                 FLOAT(IPLAST(IEL)-IPROTDW(1)+1)
          GAMTH(IP1) = GTH1*FRAC
          FRAC = 1.0 - FLOAT(IPLAST(IEL)-IP2) /
     &                 FLOAT(IPLAST(IEL)-IPROTDW(1)+1)
          GAMTH(IP2) = GTH1*FRAC
        ENDIF
C
      END DO
C
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*)  'At end of GTHCALC'
        DO IEL = 1, NEL
          IC1 = ICFRST(IEL)
          IC2 = ICLAST(IEL)
          WRITE(LUNDBG,*) 'IEL,IC1,IC2,IPFRST(IEL),IPLAST(IEL)'
          WRITE(LUNDBG,*)  IEL,IC1,IC2,IPFRST(IEL),IPLAST(IEL)
          DO IC = IC1, IC2
            IP1 = IPCO(IC)
            IP2 = IPCP(IC)
            WRITE(LUNDBG,*)  'IP1,GAMTH ',IP1,GAMTH(IP1)
            WRITE(LUNDBG,*)  'IP2,GAMTH ',IP2,GAMTH(IP2)
          END DO
        END DO
      ENDIF
 20   FORMAT(A,I5,5(1X,F10.4))
C
      RETURN
      END ! GTHCALC





      SUBROUTINE TQCALC(ITYPE)
      INCLUDE 'DFDC.INC'
C----------------------------------------------------------
C     Sets Thrust, Torque and their sensitivities 
C     wrt  QINF, OMEGA, BETA, chord(i), VA,VT, BGAM
C----------------------------------------------------------
C
C---- total forces accumulation 
      TTOT = 0.
      TVIS = 0.
      QTOT = 0.
      QVIS = 0.
      PTOT = 0.
      PVIS = 0.
C
      DO 2000 N = 1, NROTOR
C
C---- forces on blade row
      TINVR(N)  = 0.
      QINVR(N)  = 0.
      TI_OMG(N) = 0.
      QI_OMG(N) = 0.
      TI_QNF(N) = 0.
      QI_QNF(N) = 0.
C
      TVISR(N)  = 0.
      QVISR(N)  = 0.
      TV_OMG(N) = 0.
      QV_OMG(N) = 0.
      TV_QNF(N) = 0.
      QV_QNF(N) = 0.
      TV_DBE(N) = 0.
      QV_DBE(N) = 0.
C
      DO I=1, NRC
        TI_GAM(I,N) = 0.
        TI_VA (I,N) = 0.
        TI_VT (I,N) = 0.
        QI_GAM(I,N) = 0.
        QI_VA (I,N) = 0.
        QI_VT (I,N) = 0.
        TV_GAM(I,N) = 0.
        TV_VA (I,N) = 0.
        TV_VT (I,N) = 0.
        QV_GAM(I,N) = 0.
        QV_VA (I,N) = 0.
        QV_VT (I,N) = 0.
      ENDDO
C
cc      write(*,*) 'LBLDEF ',LBLDEF
      CALL ROTORVABS(N,VIND(1,1,N))
C
C---- go over radial stations, setting thrust and torque
C
C---- NOTE: should only go over blade stations to NRC-1 for tip gap case
C
      BLDS = FLOAT(NRBLD(N))
      DO 1000 I = 1, NRC
C
C---- Skip forces for tip gap (only on rotor)
        IF(I.EQ.NRC .AND. TGAP.GT.0.0 .AND. OMEGA(N).NE.0.0) THEN
           CLR(NRC,N)  = 0.0
           CDR(NRC,N)  = 0.0
           ALFAR(NRC,N)= 0.0
           GO TO 1000
        ENDIF
C
        XI = YRC(I,N)/RTIP(N)
        RA   = 0.5*(YRP(I+1,N) + YRP(I,N))
        DR   = YRP(I+1,N) - YRP(I,N)
        RDR  = YRC(I,N) * DR
        BDR  = BLDS*DR
        BRDR = BLDS*RDR
        DA   = PI*RA*DR
C
C------ set  W(Qinf,Omeg,Va,Vt)  and  Phi(Qinf,Omeg,Va,Vt)  sensitivities
        CALL WCALC(N,I,VAINF,VTINF,
     &             VTT,
     &             VAA,
     &             CI, CI_OMG,       CI_VT,
     &             SI,        SI_QNF,      SI_VA,
     &             W,   W_OMG, W_QNF, W_VT, W_VA,
     &             PHIB,P_OMG, P_QNF, P_VT, P_VA)
C
        ALFA   = BETAR(I,N) - PHIB
        AL_DBE =  1.0
        AL_P   = -1.0
C
C---- Set local Mach number in relative frame 
        MACHR(I,N) = W/VSO  
C
C
C---- Check for rotor type, actuator disk or blade row
C     If no blade has been defined an inviscid thrust and torque
C     calculation will be made for the actuator disk using circulation.  
C     Otherwise blade element theory will be used to calculate forces 
C     on the blades using aerodynamic characteristics for the blade elements.
C
        IF(IRTYPE(N).EQ.2) THEN
C
        IF(ITYPE.EQ.1) THEN
C------- analysis case:  fix local Beta (except for pitch change)
C------- set alfa(Gi,dBeta,Qinf,Omeg,Va,Vt) sensitivites
         ALFA   = BETAR(I,N) - PHIB
         ALFAR(I,N) = ALFA
         AL_GI  = 0.
         AL_DBE = 1.0
         AL_OMG = -P_OMG
         AL_QNF = -P_QNF
         AL_VT  = -P_VT
         AL_VA  = -P_VA
C
C------- set CL(Gi,dBeta,Qinf,Omeg,Va,Vt) sensitivites
         REY = CHR(I,N)*ABS(W) * RHO/RMU
         SECSIG = BLDS*CHR(I,N)/(2.0*PI*YRC(I,N))
         SECSTAGR = 0.5*PI - BETAR(I,N)
         CALL GETCLCDCM(N,I,XI,ALFA,W,REY,SECSIG,SECSTAGR,
     &                  CLR(I,N),CL_AL,CL_W,
     &                  CLMAX,CLMIN,DCLSTALL,LSTALLR(I,N),
     &                  CDR(I,N),CD_ALF,CD_W,CD_REY,
     &                  CMR(I,N),CM_AL,CM_W)
         CLALF(I,N) = CL_AL
         CL_GI  = CL_AL*AL_GI
         CL_DBE = CL_AL*AL_DBE
         CL_OMG = CL_AL*AL_OMG + CL_W*W_OMG
         CL_QNF = CL_AL*AL_QNF + CL_W*W_QNF
         CL_VT  = CL_AL*AL_VT  + CL_W*W_VT
         CL_VA  = CL_AL*AL_VA  + CL_W*W_VA
C
C------- set c(Gi,Qinf,Omeg,Va,Vt) sensitivites  (chord is fixed)
         CH_GI  = 0.
         CH_OMG = 0.
         CH_QNF = 0.
         CH_VT  = 0.
         CH_VA  = 0.
C
        ELSE IF(ITYPE.EQ.2) THEN
C---- design case:  fix local CL and set chord based on circulation
C---- update design CL arrays
C
        IF(OMEGA(N).GT.0.0)THEN
          DO J=1,NRC
            CLDES(J) = CLPOS(J)
          ENDDO
        ELSE
          DO J=1,NRC
            CLDES(J) = CLNEG(J)
          ENDDO
        ENDIF 
C
C------- set alfa(Gi,dBeta,Adv,Adw,Vt) sensitivites
         CLR(I,N) = CLDES(I)
C
         SECSIG = BLDS*CHR(I,N)/(2.0*PI*YRC(I,N))
         SECSTAGR = 0.5*PI - BETAR(I,N)
         CALL GETALF(N,I,XI,SECSIG,SECSTAGR,
     &               CLR(I,N),W,ALFA,AL_CL,AL_W,LSTALLR(I,N))
cc         write(*,*) 'tq2 getalf i,cl,w,alf ',i,clr(i),w,alfa/dtr
C
         AL_GI  = 0.
         AL_DBE = 0.
         AL_OMG = AL_W*W_OMG
         AL_QNF = AL_W*W_QNF
         AL_VT  = AL_W*W_VT
         AL_VA  = AL_W*W_VA
C
C------- set CL(Gi,dBeta,Adv,Adw,Vt) sensitivites
         CL_GI = 0.
         CL_DBE = 0.
         CL_OMG = 0.
         CL_QNF = 0.
         CL_VT  = 0.
         CL_VA  = 0.
C
C------- set c(Gi,Adv,Adw,Vt) sensitivites
         CHNEW = 2.0*BGAM(I,N) / (BLDS*W*CLR(I,N))
C--- Check for chord going zero or negative and use nearby station data
C    for this iteration
         IF(CHNEW.LE.0.0) THEN
cc           write(*,*) 'TQCALC negative chord @I = ',I,CHNEW
           IF(I.EQ.1) THEN
            CHR(I,N)  = CHR(I+1,N)
           ELSEIF(I.EQ.II) THEN
            CHR(I,N)  = CHR(I-1,N)
           ELSE
            CHR(I,N)  = 0.5*(CHR(I-1,N)+CHR(I+1,N))
           ENDIF
            CH_GI  = 0.0
            CH_OMG = 0.0
            CH_QNF = 0.0
            CH_VT  = 0.0
            CH_VA  = 0.0
         ELSE
          CHR(I,N) = 2.0*BGAM(I,N) / (BLDS*W*CLR(I,N))
cc          write(*,*) 'tq2 bgam,cl,ch ',bgam(i,n),clr(i,n),chr(i,n)
          CH_GI    = 2.0           / (BLDS*W*CLR(I,N))
          CH_OMG = (-CHR(I,N)/W) * W_OMG
          CH_QNF = (-CHR(I,N)/W) * W_QNF
          CH_VT  = (-CHR(I,N)/W) * W_VT
          CH_VA  = (-CHR(I,N)/W) * W_VA
         ENDIF
C
         BETAR(I,N) = ALFA + PHIB
         ALFAR(I,N) = ALFA
C
        ELSE IF(ITYPE.EQ.3) THEN
C------- design case:  fix local chord and set angles based on CL
C
C------- set CL(Gi,dBeta,Adv,Adw,Vt) sensitivites
         CLR(I,N) = 2.0*BGAM(I,N) / (BLDS*W*CHR(I,N))
         CL_GI  =   2.0           / (BLDS*W*CHR(I,N))
         CL_DBE = 0.
         CL_OMG = (-CLR(I,N)/W) * W_OMG
         CL_QNF = (-CLR(I,N)/W) * W_QNF
         CL_VT  = (-CLR(I,N)/W) * W_VT
         CL_VA  = (-CLR(I,N)/W) * W_VA
C
C------- set alfa(Gi,dBeta,Adv,Adw,Vt) sensitivites
         SECSIG = BLDS*CHR(I,N)/(2.0*PI*YRC(I,N))
         SECSTAGR = 0.5*PI - BETAR(I,N)
         CALL GETALF(N,I,XI,SECSIG,SECSTAGR,
     &               CLR(I,N),W,ALFA,AL_CL,AL_W,LSTALLR(I,N))
cc         write(*,*) 'tq3 i,cl,ch,w,alf ',i,clr(i),chr(i),w,alfa/dtr
         AL_GI  = AL_CL*CL_GI
         AL_DBE = AL_CL*CL_DBE
         AL_OMG = AL_CL*CL_OMG + AL_W*W_OMG
         AL_QNF = AL_CL*CL_QNF + AL_W*W_QNF
         AL_VT  = AL_CL*CL_VT  + AL_W*W_VT
         AL_VA  = AL_CL*CL_VA  + AL_W*W_VA
C
C------- set c(Gi,Adv,Adw,Vt) sensitivites
         CH_GI  = 0.
         CH_OMG = 0.
         CH_QNF = 0.
         CH_VT  = 0.
         CH_VA  = 0.
C
         BETAR(I,N) = ALFA + PHIB
         ALFAR(I,N) = ALFA
C
        ENDIF
C
C
C=================================================================
C
        RER(I,N) = CHR(I,N)*ABS(W) * RHO/RMU
        RE_W     = CHR(I,N)        * RHO/RMU
        RE_CH    =          ABS(W) * RHO/RMU
C
C------ set Re(Gi,Adv,Adw,Vt) sensitivites
        RE_GI  = RE_CH*CH_GI
        RE_OMG = RE_CH*CH_OMG + RE_W*W_OMG
        RE_QNF = RE_CH*CH_QNF + RE_W*W_QNF
        RE_VT  = RE_CH*CH_VT  + RE_W*W_VT 
        RE_VA  = RE_CH*CH_VA  + RE_W*W_VA
C
C------ set CM and (not used at present) sensitivites
C------ set CD(Gi,dBeta,Adv,Adw,Vt) sensitivites
        SECSIG = BLDS*CHR(I,N)/(2.0*PI*YRC(I,N))
        SECSTAGR = 0.5*PI - BETAR(I,N)
        CALL GETCLCDCM(N,I,XI,ALFA,W,RER(I,N),SECSIG,SECSTAGR,
     &                 CLR(I,N),CL_AL,CL_W,
     &                 CLMAX,CLMIN,DCLSTALL,LSTALLR(I,N),
     &                 CDR(I,N),CD_AL,CD_W,CD_RE,
     &                 CMR(I,N),CM_AL,CM_W)
        CLALF(I,N) = CL_AL
cc        write(*,97) 'tqcalc alfa,cl,cd,cm ',i,alfa,clr(i),cdr(i)
 97     format(A,I5,5(1x,F12.6))
        CD_GI  = CD_AL*AL_GI  + CD_RE*RE_GI
        CD_OMG = CD_AL*AL_OMG + CD_RE*RE_OMG + CD_W*W_OMG
        CD_QNF = CD_AL*AL_QNF + CD_RE*RE_QNF + CD_W*W_QNF
        CD_VT  = CD_AL*AL_VT  + CD_RE*RE_VT  + CD_W*W_VT
        CD_VA  = CD_AL*AL_VA  + CD_RE*RE_VA  + CD_W*W_VA
        CD_DBE = CD_AL*AL_DBE 
C
C
        HRWC    = 0.5*RHO*W*CHR(I,N)
        HRWC_W  = 0.5*RHO  *CHR(I,N)
        HRWC_CH = 0.5*RHO*W
C
C*******************************************************
C------ Viscous Thrust & Power contributions on real prop
C
C------ dTv ( Cd , S , W , c ) sensitivites
        DTV    = -HRWC   *CDR(I,N)*SI*BDR
C
        DTV_CD = -HRWC            *SI*BDR
        DTV_SI = -HRWC   *CDR(I,N)   *BDR
        DTV_W  = -HRWC_W *CDR(I,N)*SI*BDR
        DTV_CH = -HRWC_CH*CDR(I,N)*SI*BDR
C
C------ set Tv(Gi,dBeta,Adv,Vt) sensitivites using chain rule
        DTV_GI  = DTV_CD*CD_GI  + DTV_CH*CH_GI
        DTV_DBE = DTV_CD*CD_DBE
        DTV_OMG = DTV_CD*CD_OMG + DTV_CH*CH_OMG
     &                          + DTV_W * W_OMG
        DTV_QNF = DTV_CD*CD_QNF + DTV_CH*CH_QNF
     &                          + DTV_W * W_QNF
        DTV_VT  = DTV_CD*CD_VT  + DTV_CH*CH_VT
     &                          + DTV_W * W_VT
        DTV_VA  = DTV_CD*CD_VA  + DTV_CH*CH_VA
     &          + DTV_SI*SI_VA  + DTV_W * W_VA
C
C------ accumulate viscous Thrust and sensitivities
        TVISR(N)  = TVISR(N)  + DTV
        TV_OMG(N) = TV_OMG(N) + DTV_OMG
        TV_QNF(N) = TV_QNF(N) + DTV_QNF
        TV_DBE(N) = TV_DBE(N) + DTV_DBE
C
        TV_GAM(I,N) = DTV_GI
        TV_VA (I,N) = DTV_VA
        TV_VT (I,N) = DTV_VT
C
C
C------ dQv( Cd , C , W , c ) 
        DQV    = -HRWC   *CDR(I,N)*CI*BRDR
C
        DQV_CD = -HRWC          *CI*BRDR
        DQV_CI = -HRWC   *CDR(I,N)   *BRDR
        DQV_W  = -HRWC_W *CDR(I,N)*CI*BRDR
        DQV_CH = -HRWC_CH*CDR(I,N)*CI*BRDR
        DQV_OM = -HRWC   *CDR(I,N)*CI*BRDR
C
C------ set Pv(Gi,dBeta,Adv,Vt) sensitivites using chain rule
        DQV_GI  = DQV_CD*CD_GI  + DQV_CH*CH_GI
        DQV_DBE = DQV_CD*CD_DBE
        DQV_OMG = DQV_OM 
     &          + DQV_CD*CD_OMG + DQV_CH*CH_OMG
     &          + DQV_CI*CI_OMG + DQV_W * W_OMG
        DQV_QNF = DQV_CD*CD_QNF + DQV_CH*CH_QNF
     &                          + DQV_W * W_QNF
        DQV_VT  = DQV_CD*CD_VT  + DQV_CH*CH_VT
     &          + DQV_CI*CI_VT  + DQV_W * W_VT
        DQV_VA  = DQV_CD*CD_VA  + DQV_CH*CH_VA
     &                          + DQV_W * W_VA
C
C------ accumulate viscous Power and sensitivities
        QVISR(N)  = QVISR(N)  + DQV
        QV_OMG(N) = QV_OMG(N) + DQV_OMG
        QV_QNF(N) = QV_QNF(N) + DQV_QNF
        QV_DBE(N) = QV_DBE(N) + DQV_DBE
C
        QV_GAM(I,N) = DQV_GI
        QV_VA (I,N) = DQV_VA
        QV_VT (I,N) = DQV_VT
C
        ENDIF ! end of check for IRTYPE=2 (blade defined)
C
C
C*******************************************************
C------ Inviscid Thrust & Power contributions on rotor
C
C------ dTi( Gi , C( Omg Vt ) )
        DTI    = -RHO*BGAM(I,N)*CI*DR
C
        DTI_CI = -RHO*BGAM(I,N)   *DR
        DTI_GI = -RHO          *CI*DR
C
C------ dTi( Adv , Vt(Adw Gj) )
        DTI_VT  = DTI_CI*CI_VT
        DTI_OMG = DTI_CI*CI_OMG
C
C------ accumulate inviscid Thrust and sensitivities
        TINVR(N)  = TINVR(N)  + DTI
        TI_OMG(N) = TI_OMG(N) + DTI_OMG
C------ Resolve dTi dependencies ( Vt ) to Gamma
        TI_GAM(I,N) = DTI_GI
        TI_VA (I,N) = 0.0
        TI_VT (I,N) = DTI_VT
C
C
C------ dQi( S(V Va) , Gi )
        DQI    = RHO*BGAM(I,N)*SI*RDR
C
        DQI_GI = RHO          *SI*RDR
        DQI_SI = RHO*BGAM(I,N)   *RDR
C
C------ dQi( Vai , Qinf, Gi )
        DQI_QNF = DQI_SI*SI_QNF
        DQI_VA  = DQI_SI*SI_VA
C
C------ accumulate inviscid Power and sensitivities
        QINVR(N)  = QINVR(N)  + DQI
        QI_OMG(N) = 0.0
        QI_QNF(N) = QI_QNF(N) + DQI_QNF 
C------ Save dQi dependencies to BGAM,VA,VT
        QI_GAM(I,N) = DQI_GI
        QI_VA (I,N) = DQI_VA
        QI_VT (I,N) = 0.0
C
C*******************************************************
C------ Save blade thrust and torque distributions (per blade, per span)
        IF(BLDS.NE.0.0) THEN
         DTII(I,N) = DTI/BDR
         DQII(I,N) = DQI/BDR
         DTVI(I,N) = DTV/BDR
         DQVI(I,N) = DQV/BDR
        ELSE
C------ or actuator disk thrust and torque distributions (per span)
         DTII(I,N) = DTI/DR
         DQII(I,N) = DQI/DR
         DTVI(I,N) = DTV/DR
         DQVI(I,N) = DQV/DR
        ENDIF
C------ static pressure rise at this radial station
        DPSI(I,N) = (DTI+DTV)/DA
C
 1000 CONTINUE
C
C---- forces for this rotor
      TTOTR(N) = TINVR(N) + TVISR(N)
      QTOTR(N) = QINVR(N) + QVISR(N)
C
C---- total forces (all rotors)
      TTOT = TTOT + TTOTR(N)
      TVIS = TVIS + TVISR(N)
      QTOT = QTOT + QTOTR(N)
      QVIS = QVIS + QVISR(N)
C
C---- Derive some power quantities from torque
      PVIS = PVIS + QVISR(N) * OMEGA(N)
      PTOT = PTOT + QTOTR(N) * OMEGA(N)
C
 2000 CONTINUE
C
c      write(*,99) 'TQcalc '
c      write(*,99) ' Tinv,Ti_qnf,Ti_omg       ',tinv,ti_qnf,ti_omg
c      write(*,99) ' Tvis,Tv_qnf,Tv_omg,Tv_be ',tvis,tv_qnf,tv_omg,tv_dbe
c      write(*,99) ' Qinv,Qi_qnf,Qi_omg       ',qinv,qi_qnf,qi_omg
c      write(*,99) ' Qvis,Qv_qnf,Qv_omg,Pv_be ',qvis,qv_qnf,qv_omg,qv_dbe
c      write(*,99) ' Pinv,Pvis,Ptot           ',pinv,pvis,ptot
C
 99   format(A,4(1x,F12.3))
C
      RETURN
      END ! TQCALC


      SUBROUTINE WCALC(N,I,VAIN,VTIN,
     &                 VTT,
     &                 VAA,
     &                 CI, CI_OMG,       CI_VT,
     &                 SI,        SI_QNF,      SI_VA,
     &                 W,   W_OMG, W_QNF, W_VT, W_VA,
     &                 PHIB,P_OMG, P_QNF, P_VT, P_VA)
C
C---- Calculate velocity components at radial station I on rotor blade
C
      INCLUDE 'DFDC.INC'
C
      IF(LBLBL)THEN
        VFAC = BBVFAC(I,N)
      ELSE
        VFAC = 1.0
      ENDIF
C        
C---- At blade lifting line use averaged circulation for tangential velocity 
      VTBG = BGAM(I,N)*PI2I/YRC(I,N)
      VTT  = VABS(3,I,N) - 0.5*VTBG
      VAA  = VIND(1,I,N)
C
C---- Freestream, body induced and added inflow velocities     
      CALL UVINFL(YRC(I,N),VAIN,VTIN)
C     
      CI     = -YRC(I,N)*OMEGA(N) + VTIN + VTT
      CI_OMG = -YRC(I,N)
      CI_VT  =                          1.0
C     
c      SI     = QINF          + VAIN + VAA         ! v0.70
      SI     = (QINF + VAA + VAIN) * VFAC          ! BB entry
      SI_QNF = 1.0
      SI_VA  =                          1.0
C
      WSQ = CI*CI + SI*SI
      W = SQRT(WSQ)
      W_OMG = (CI*CI_OMG            )/W
      W_QNF = (            SI*SI_QNF)/W
      W_VT  = (CI*CI_VT             )/W
      W_VA  = (            SI*SI_VA )/W
C     
      PHIB = ATAN2(SI,-CI)
      P_OMG = (             SI*CI_OMG)/WSQ
      P_QNF = (-CI*SI_QNF            )/WSQ
      P_VT  = (             SI*CI_VT )/WSQ
      P_VA  = (-CI*SI_VA             )/WSQ
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*)  'WCALC @I= ',I
        WRITE(LUNDBG,99) 'QINF YRC  ',QINF,YRC(I,N)
        WRITE(LUNDBG,99) 'OMEG      ',OMEGA(N)
        WRITE(LUNDBG,99) 'VT,VA     ',VTT,VAA
        WRITE(LUNDBG,99) 'VTIN,VAIN ',VTIN,VAIN
        WRITE(LUNDBG,99) 'CI,SI,W   ',CI,SI,W
        WRITE(LUNDBG,99) 'PHI       ',PHIB/DTR
      ENDIF
C
 99   FORMAT(A,5(1X,f11.6))
C
      RETURN
      END ! WCALC


      SUBROUTINE ROTRPRT(LU)
      INCLUDE 'DFDC.INC'
      DIMENSION W1(IRX)
      CHARACTER SCHAR*1, RTYPE*30
C---------------------------------------------
C     Dumps operating state of case to unit LU
C---------------------------------------------
C
      WRITE (LU,1000)
      WRITE(LU,1001) NAME
      IF(.NOT.LCONV) WRITE(LU,1002)
      IF(NINFL.GT.0) WRITE(LU,1005)
C
 1000 FORMAT(/1X,76('-'))
 1001 FORMAT(' DFDC  Case:  ', A32)
 1002 FORMAT(/19X,'********** NOT CONVERGED **********',/ )
 1003 FORMAT(1X,76('-'))
 1004 FORMAT(50X)
 1005 FORMAT(' (External slipstream present)' )
C
C---- dimensional thrust, power, torque, rpm
      TDIM  = TTOT + TDUCT
      QDIM  = QTOT
      PDIM  = PTOT
      TVDIM = TVIS
      PVDIM = PVIS
      TINV  = TTOT - TVIS
      QINV  = QTOT - QVIS
      PINV  = PTOT - PVIS
C
C---- Define reference quantities for coefficients
      RREF = RTIP(1)
      AREF = ADISK(1)
      OMEGREF = OMEGA(1)
C
C---- standard thrust/power coefficients based on rotational speed
      DREF = 2.0*RREF
      EN = ABS(OMEGREF*PI2I)
      IF(EN.EQ.0.0) THEN
        CT = 0.0
        CP = 0.0
      ELSE
        CT = TDIM/(RHO*EN**2*DREF**4)
        CP = PDIM/(RHO*EN**3*DREF**5)
      ENDIF
C---- standard thrust/power coefficients based on forward speed
      IF(QINF.GT.0.0) THEN
        TC = TDIM/(0.5*RHO*QINF**2 * PI*RREF**2)
        PC = PDIM/(0.5*RHO*QINF**3 * PI*RREF**2)
      ELSE
        TC = 0.0
        PC = 0.0
      ENDIF
C---- thrust/power coefficients based on tip speed
C     uses helicopter nomenclature for CT0,CP0,FOM
      VTIP = ABS(OMEGREF*RREF)
      ADV  = QINF/VTIP
      IF(VTIP.NE.0.0) THEN
       CT0  = TDIM/(RHO*AREF*VTIP**2)
       CP0  = PDIM/(RHO*AREF*VTIP**3)
      ELSE
       CT0 = 0.0
       CP0 = 0.0
      ENDIF
      IF(CT0.GE.0.0 .AND. CP0.NE.0.0) THEN
        FOM = ABS(CT0)**1.5 / CP0 / 2.0
      ELSE
        FOM = 0.0
      ENDIF
C
C---- overall efficiency (all thrust components)
      IF(PDIM.NE.0.0) EFFTOT = QINF*TDIM/PDIM
C---- induced efficiency (all thrust components)
      IF(PINV.NE.0.0) EFFIND = QINF*(TINV+TDUCT)/PINV
C---- ideal (actuator disk) efficiency
      IF(TC.EQ.0) THEN
        EIDEAL = 0.0
      ELSE
        TCLIM = MAX( -1.0 , TC )
        EIDEAL = 2.0 / (1.0 + SQRT(TCLIM + 1.0))
      ENDIF
C
C---- Dump overall case data
C
      WRITE(LU,1003) 
      IF(IRTYPE(1).EQ.2) THEN
        IF(LBLBL) THEN
          WRITE(LU,1010)
        ELSE
          WRITE(LU,1011) 
        ENDIF
      ELSE
          WRITE(LU,1008)
      ENDIF
C
      WRITE(LU,1012) QINF,  ALTH, DELTAT,
     &               RHO,   VSO,   RMU,
     &               TDIM,  PDIM,  EFFTOT,
     &               TVDIM, PVDIM, EFFIND,
     &               TDUCT, QDIM,  EIDEAL
C
      WRITE(LU,1004)
C
      WRITE(LU,1014) AREF, RREF, OMEGREF
C---- Thrust/power coefficients based on rotational speed (propeller syntax)
      WRITE(LU,1015) CT, CP, ADV*PI
C---- Thrust/power coefficients based on forward speed (propeller syntax)
      WRITE(LU,1016) TC, PC, ADV
C---- Thrust/power coefficients based on tip speed (helicopter nomenclature)
      WRITE(LU,1017) CT0,CP0,FOM
C
 1008 FORMAT(' Flow Condition and total Forces')
 1010 FORMAT(' Flow Condition and total Forces',17X,
     &        'Corrected for blade blockage')
 1011 FORMAT(' Flow Condition and total Forces',17X,
     &        'No blade blockage correction')
 1012 FORMAT(/'  Vinf(m/s) :',F10.3,4X,'Alt.(km)   :',F9.3, 5X,
     &   'DeltaT(dgC):',F9.4,
     &       /' rho(kg/m3) :',F11.4,3X,'Vsound(m/s):',F9.3, 5X,
     &   'mu(kg/m-s) :',E11.4,
     &       /' Thrust(N)  :',G11.3,3X,'Power(W)   :',G11.3,3X,
     &   'Efficiency :',F9.4,
     &       /' Tvisc (N)  :',F11.4,3X,'Pvisc(W)   :',G11.3,3X,
     &   'Induced Eff:',F9.4,
     &       /' Tduct(N)   :',F11.4,3X,'torQue(N-m):',G11.3,3X,
     &   'Ideal Eff  :',F9.4)
C
 1014 FORMAT('  Area:', F11.5,'  Radius:', F11.5,' Omega:', F11.5,
     &   '  Reference data')
 1015 FORMAT('    Ct:', F11.5,'      Cp:', F11.5,'     J:', F11.5,
     &   '  by(Rho,N,Dia)')
 1016 FORMAT('    Tc:', F11.5,'      Pc:', F11.5,'   adv:', F11.5,
     &   '  by(Rho,Vinf,Area)  ')
 1017 FORMAT('   CT0:',F11.5, '     CP0:', F11.5,'   FOM:', F11.5,
     &   '  by(Rho,R*Omg,Area)')
C
C
C---- Display operating state for each rotor
C
      DO N = 1, NROTOR
C
      IADD = 1
cc      IF(LU.EQ.LUWRIT) IADD = INCR
C
      WRITE(LU,1003)
C
C---- dimensional thrust, power, torque, rpm
      TDIM  = TTOTR(N)
      QDIM  = QTOTR(N)
      PDIM  = QDIM     * OMEGA(N)
      TVDIM = TVISR(N)
      PVDIM = QVISR(N) * OMEGA(N)
      TINV  = TINVR(N)
      QINV  = QINVR(N)
      PINV  = QINV     * OMEGA(N)
C
C
C---- Define reference quantities for coefficients
      RREF = RTIP(N)
      AREF = ADISK(N)
      OMEGREF = OMEGA(N)
      RPM = 30.0*OMEGREF/PI
C
C---- standard thrust/power coefficients based on rotational speed
      DIA = 2.0*RREF
      EN = ABS(OMEGREF*PI2I)
      IF(EN.EQ.0.0) THEN
        CT = 0.0
        CP = 0.0
      ELSE
        CT = TDIM/(RHO*EN**2*DIA**4)
        CP = PDIM/(RHO*EN**3*DIA**5)
      ENDIF
C
C---- standard thrust/power coefficients based on forward speed
      IF(QINF.GT.0.0) THEN
        TC = TDIM/(0.5*RHO*QINF**2 * PI*RREF**2)
        PC = PDIM/(0.5*RHO*QINF**3 * PI*RREF**2)
      ELSE
        TC = 0.0
        PC = 0.0
      ENDIF
C---- thrust/power coefficients based on tip speed
C     uses helicopter nomenclature for CT0,CP0,FOM
      VTIP = ABS(OMEGREF*RREF)
      IF(VTIP.NE.0.0) THEN
       CT0  = TDIM/(RHO*AREF*VTIP**2)
       CP0  = PDIM/(RHO*AREF*VTIP**3)
       ADV  = QINF/VTIP
      ELSE
       CT0  = 0.0
       CP0  = 0.0
       ADV  = 0.0
      ENDIF
C
C---- efficiency for rotor 
      IF(PDIM.NE.0.0) EFFTOT = QINF*TDIM/PDIM
C---- induced efficiency for rotor 
      IF(PINV.NE.0.0) EFFIND = QINF*TINV/PINV
C---- ideal (actuator disk) efficiency
      IF(TC.EQ.0) THEN
        EIDEAL = 0.0
      ELSE
        TCLIM = MAX( -1.0 , TC )
        EIDEAL = 2.0 / (1.0 + SQRT(TCLIM + 1.0))
      ENDIF
C
      SIGMA = 0.0
      IF(IRTYPE(N).EQ.1) THEN
       IF(OMEGA(N).EQ.0.0) THEN
        RTYPE = 'Actuator Disk Stator'
       ELSE
        RTYPE = 'Actuator Disk Rotor '
       ENDIF
      ELSEIF(IRTYPE(N).EQ.2) THEN
       IF(OMEGA(N).EQ.0.0) THEN
        RTYPE = 'Bladed Stator '
       ELSE
        RTYPE = 'Bladed Rotor  '
       ENDIF
C---- Overall blade solidity (measured at 3/4 Rtip)
       CALL SPLINE(CHR(1,N),W1,YRC,NRC)
       CH34 = SEVAL(0.75*RTIP(N),CHR(1,N),W1,YRC,NRC)
       SIGMA = FLOAT(NRBLD(N))*CH34/RTIP(N)/PI
      ELSE
       RTYPE = 'Undefined disk'
      ENDIF
C
       IF(SIGMA.NE.0.0) THEN
        CTOS = CT0 / SIGMA
       ELSE
        CTOS = 0.0
       ENDIF
C
      WRITE(LU,1110) N,     RTYPE,
     &               NRBLD(N), RPM,   ADV,
     &               TDIM,  PDIM,  EFFTOT,
     &               TVDIM, PVDIM, EFFIND,
     &               QDIM,  QINV,  EIDEAL,
     &               RTIP(N),RHUB(N),VAAVG(N)
C
      WRITE(LU,1004)
      WRITE(LU,1114) AREF, RREF, OMEGREF
C---- Thrust/power coefficients based on rotational speed (propeller syntax)
      WRITE(LU,1115) CT, CP, ADV*PI
C---- Thrust/power coefficients based on forward speed (propeller syntax)
      WRITE(LU,1116) TC, PC, ADV
C---- Thrust/power coefficients based on tip speed (helicopter nomenclature)
      WRITE(LU,1118) CT0,CP0
      WRITE(LU,1119) SIGMA,CTOS
C
      WRITE(LU,1003)
C
      IF(OMEGA(N).LE.0.0 .AND. IRTYPE(N).EQ.2) THEN
         IF(LCOORD) THEN
           WRITE(LU,*)'                    -local coords-'
         ELSE
           WRITE(LU,*)'                    -global coords-'
         ENDIF
      ENDIF
C
C
      IF(IRTYPE(N).EQ.2) THEN
C---- Blade is defined, print blade data
C
C---- Overall blade solidity (measured at 3/4 Rtip)
       CALL SPLINE(CHR(1,N),W1,YRC(1,N),NRC)
       CH34 = SEVAL(0.75*RTIP(N),CHR(1,N),W1,YRC,NRC)
       SIGMA = FLOAT(NRBLD(N))*CH34/RTIP(N)/PI
       IF(SIGMA.NE.0.0) THEN
        CTOS = CT0 / SIGMA
       ELSE
        CTOS = 0.0
       ENDIF
C
C----- find maximum RE on blade
       REMAX = 0.0
       DO I=1, NRC
         REMAX = MAX(RER(I,N),REMAX)
       END DO
       REEXP = 1.0
       IF(REMAX.GE.1.0E6) THEN
         REEXP = 6.0
       ELSEIF(REMAX.GE.1.0E3) THEN
         REEXP = 3.0
       ENDIF
       IF(REEXP.EQ.1.0) THEN
         WRITE(LU,1120) 
        ELSE
         WRITE(LU,1122) IFIX(REEXP)
       ENDIF
C
C---- NOTE: should only dump blade data to NRC-1 for tip gap case
C
C       LSTALLR(NRC,N)=.FALSE.
C
       DO I=1, NRC, IADD
         XI  = YRC(I,N)/RTIP(N)
         CHI = CHR(I,N)/RTIP(N)
         XRE = RER(I,N)/(10.0**REEXP)
C
         IF(LCOORD .AND. OMEGA(N).LE.0.0) THEN
           BDEG = (PI-BETAR(I,N))/DTR
           ALDEG= -ALFAR(I,N)/DTR
         ELSE
           BDEG = BETAR(I,N)/DTR
           ALDEG= ALFAR(I,N)/DTR
         ENDIF
C
         IF(I.EQ.NRC.AND.TGAP.GT.0.0.AND.OMEGA(N).NE.0.0) ALDEG=0.0
C
         SCHAR = ' '
         IF(LSTALLR(I,N)) THEN
           IF(I.EQ.NRC .AND. TGAP.GT.0.0.AND.OMEGA(N).NE.0.0) THEN
             SCHAR = ' '
           ELSE
             SCHAR = 's'
           ENDIF
         ENDIF
C
         WRITE(LU,1130) I,XI,CHI,BDEG,ALDEG,CLR(I,N),SCHAR,CDR(I,N),
     &                  XRE,MACHR(I,N),BGAM(I,N)
cc     &    EFFI,EFFP(I)
       END DO
C
      ELSE
C---- Print actuator disk datal
       WRITE(LU,1220) 
       DO I=1, NRC, IADD
         XI  = YRC(I,N)/RTIP(N)
         WRITE(LU,1230)
     &     I,XI,MACHR(I,N),BGAM(I,N)
       END DO
      ENDIF
C
      END DO
C
cc      WRITE(LU,1000)
cc      WRITE(LU,*   ) ' '
C
      RETURN
C....................................................................
C
 1110 FORMAT(1X,'Disk #',I3,4X,A,
     &       /' # blades   :',I3,  11X,'RPM        :',F11.3,3X,
     &   'adv. ratio :',F9.4,
     &       /' Thrust(N)  :',G11.3,3X,'Power(W)   :',G11.3,3X,
     &   'Efficiency :',F9.4,
     &       /' Tvisc (N)  :',F11.4,3X,'Pvisc(W)   :',G11.3,3X,
     &   'Induced Eff:',F9.4,
     &       /' torQue(N-m):',F11.4,3X,'Qvisc(N-m) :',G11.3,3X,
     &   'Ideal Eff  :',F9.4,
     &       /' radius(m)  :',F9.4, 5X,'hub rad.(m):',F9.4, 5X,
     &   'VAavg (m/s):',F9.4)
C
 1114 FORMAT('  Area:', F11.5,'  Radius:',F11.5,' Omega:',F11.5,
     &   '  Reference data')
 1115 FORMAT('    Ct:', F11.5,'      Cp:',F11.5,'     J:',F11.5,
     &   '  by(Rho,N,Dia)')
 1116 FORMAT('    Tc:', F11.5,'      Pc:',F11.5,'   adv:',F11.5,
     &   '  by(Rho,Vinf,Area)  ')
 1118 FORMAT('   CT0:', F11.5,'     CP0:',F11.5,18X,
     &   '  by(Rho,R*Omg,Area)')
 1119 FORMAT(' Sigma:', F11.5,' CT0/Sig:',F11.5)
C
C---- Rotor data
 1120 FORMAT('   i   r/R    c/R     beta deg alfa     CL     CD',
     &       '      RE   ',    '   Mach    B*Gam')
 1122 FORMAT('   i   r/R    c/R     beta deg alfa     CL     CD',
     &       '    REx10^',I1,  '   Mach    B*Gam')
 1130 FORMAT(2X,I2,F7.3,F8.4,F8.2,F7.2,1X,F8.3,A1,F7.4,
     &       F8.2,F8.3,F9.3)
C
C---- Actuator disk data
 1220 FORMAT('   i   r/R    c/R     beta deg alfa     CL     CD',
     &       '      RE   ',    '   Mach    B*Gam')
 1230 FORMAT(2X,I2,F7.3,8X, 8X,1X,  7X, 1X,   8X,     1X,7X,
     &       7X,  F8.3,F9.3)
C
C
c 1120 FORMAT(/'   i    r/R     c/R    beta(deg)',
c     & '    CL       Cd     RE        Mach        B*Gam')
c 1030 FORMAT(2X,I2,F7.3,8X,8X,8X,2X,1X,F8.4,1X,
c     &       F7.3,F10.3)
c   i    r/R     c/R    beta(deg)alfa    CL      Cd      RE     Mach     B*Gam\
c   i    r/R     c/R    beta(deg)alfa    CL      Cd    REx10^I  Mach     B*Gam')
cxxiiffffff7fffffff8fffffff8xxxfffffff8xSffffff8xffffff7xxxxffffff7xfffffffff9
cxx2i f7.3   f8.4    f8.2    3x  f8.3  x   f8.4 x  f7.2 4x    f7.3 x  f10.3
C
      END ! ROTRPRT




      SUBROUTINE NFCALC
C----------------------------------------------------------------
C     Calculate near-field forces on rotor, momentum, power at disk 
C     Inviscid thrust is calculated from RPM and circulation
C     This routine is approximate (inviscid only), superseded by 
C     routine TQCALC.
C----------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DO N = 1, NROTOR
C
C---- Calculate rotor inviscid thrust from circulation
      OMEG = OMEGA(N)
      THRUST = 0.0
      DO IR = 1, NRC
        DR = YRP(IR+1,N)-YRP(IR,N)
        DA = PI*(YRP(IR+1,N)**2-YRP(IR,N)**2)
C---- use theta velocity at blade (1/2 of induced Vt in slipstream)
        VTB = 0.5*BGAM(IR,N)*PI2I/YRC(IR,N)
        WTB = VTB - OMEG*YRC(IR,N)
        DTHR = DR*RHO * (-WTB)*BGAM(IR,N)
        THRUST = THRUST + DTHR
      END DO
      WRITE(*,*) 'Thrust from rotor GAM and OMEGA',THRUST
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'Thrust from rotor GAM and OMEGA',THRUST
      ENDIF
C
C---- Near-field thrust from momentum theory using rotor velocities
      AREA    = 0.0
      THRUSTM = 0.0
      POWERM  = 0.0
      POWERMT = 0.0
      DO IR = 1, NRC
        DR = YRP(IR+1,N)-YRP(IR,N)
        DA = PI*(YRP(IR+1,N)**2-YRP(IR,N)**2)
C
        US = VABS(1,IR,N)
        VS = VABS(2,IR,N)
        WS = VABS(3,IR,N)   
        VSQ = US*US + VS*VS + WS*WS
        DT  = DA*RHO*US*(US-QINF)
        DP  = DA*RHO*US*(0.5*VSQ - 0.5*QINF*QINF)
        DPT = DA*RHO*US*(0.5*WS**2)
C
        AREA    = AREA    + DA
        THRUSTM = THRUSTM + DT
        POWERM  = POWERM  + DP
        POWERMT = POWERMT + DPT
C
      END DO       
C
      WRITE(*,*) 'Momentum integration in near-field for rotor # ',N
      WRITE(*,*) '       Area   = ',AREA
      WRITE(*,*) '   mom.Thrust = ',THRUSTM
      WRITE(*,*) '   mom.Power  = ',POWERM
      WRITE(*,*) ' swirl Power  = ',POWERMT
C
      END DO
C
      RETURN
      END ! NFCALC




      SUBROUTINE FFCALC
      INCLUDE 'DFDC.INC'
C---------------------------------------------------------
C     Integrates forces (thrust and power at FF boundary)
C---------------------------------------------------------
C
      TINT  = 0.0
      TINTH = 0.0
      PINT  = 0.0
      PINTH = 0.0
      PINTT = 0.0
C
      DO IR = 1, NRP-1
C
        FMR = RHO*VABS(1,IR,1)*PI*(YRP(IR+1,1)**2 - YRP(IR,1)**2)
C
        IELO = IR2IEL(IR)
        IELP = IR2IEL(IR+1)
        IPO  = IPLAST(IELO) 
        IPP  = IPLAST(IELP) 
        XFF = 0.5*(XP(IPO) + XP(IPP)) - 0.5*(XP(IPP) - XP(IPP-1))
        YFF = 0.5*(YP(IPO) + YP(IPP)) 
        DY = YP(IPP) - YP(IPO)
C
        CALL GETUV(XFF,YFF,US,VS)
        WS = BGAMG(II-1,IR)*PI2I/YFF
        VSQ = US*US + VS*VS + WS*WS
C
        DA = 2.0*PI*YFF*DY
        FM = DA*RHO*US
        DT  = DA*RHO*US*(US-QINF)
        DP  = DA*RHO*US*(0.5*VSQ - 0.5*QINF*QINF)
        DPT = DA*RHO*US*(0.5*WS**2)
C
        AREA  = AREA  + DA
        TINT  = TINT  + DT
        PINT  = PINT  + DP
        PINTT = PINTT + DPT
C
        WRITE(*,*) 'IR,FM,FMR ',IR,FM,FMR        
C
C---- Use enthalpy, entropy and circulation at far-field
        USQ = QINF*QINF - WS*WS + 2.0*(DHG(II-1,IR)-DSG(II-1,IR))
        US = SQRT(USQ)
        VSQ = US*US + WS*WS
        DT  = FMR*(US-QINF)
        DP  = FMR*(0.5*VSQ - 0.5*QINF*QINF)
        TINTH =TINTH + DT
        PINTH =PINTH + DP
C
      END DO       
C
      WRITE(*,*) 'FFCALC Area   = ',AREA
      WRITE(*,*) '       Thrust = ',TINT
      WRITE(*,*) '       Power  = ',PINT
      WRITE(*,*) ' Swirl Power  = ',PINTT
      WRITE(*,*) '   FF  Thrust = ',TINTH
      WRITE(*,*) '   FF  Power  = ',PINTH
C  
      RETURN
      END



      SUBROUTINE STGFIND
      INCLUDE 'DFDC.INC'
C---------------------------------------------------------
C     Find stagnation points on CB and duct
C     The panel center index and X,Y location are saved
C---------------------------------------------------------
C
C---- Stagnation point on axisymmetric CB is simply the upstream point
      IEL = 1
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      ICSTG(IEL) = IC2
      XSTG(IEL) = XP(IP2)
      YSTG(IEL) = YP(IP2)
C
C---- Stagnation point on duct must be found by search in tangential velocity
      IEL = 2
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
      ICSTG(IEL) = IC1
      XSTG(IEL) = XP(IP1)
      YSTG(IEL) = YP(IP1)
C
      DO IC = IC1, IC2
        QD = -ANC(2,IC)*QCR(1,IC) + ANC(1,IC)*QCR(2,IC)
cc        write(*,*) 'ic,qd ',ic,qd
        IF(IC.GT.IC1) THEN
          IF(QD*QDOLD.LT.0.0) THEN
cc            WRITE(*,*) 'Found stagnation point at IC=',IC
cc            WRITE(*,*) ' QD, QDOLD ',QD,QDOLD
cc            write(*,*) 'xc,yc ',XC(IC),YC(IC)
C
            D1 = 0.5*DSC(IC-1) 
            D2 = 0.5*DSC(IC)
            D12 = D1 + D2
            DFRAC = D1/D12
            QFRAC = QD/(QD - QDOLD) 
cc            write(*,*) 'd1,d2,dfrac ',d1,d2,dfrac
cc            write(*,*) 'qfrac ',qfrac
            IF(QFRAC.LT.DFRAC) THEN
             ICSTG(IEL) = IC-1
             XSTG(IEL) = XC(IC-1) + (1.0-QFRAC)*D12*(-ANC(2,IC-1))
             YSTG(IEL) = YC(IC-1) + (1.0-QFRAC)*D12*( ANC(1,IC-1))
            ELSE
             ICSTG(IEL) = IC
             XSTG(IEL) = XC(IC) - (QFRAC)*D12*(-ANC(2,IC))
             YSTG(IEL) = YC(IC) - (QFRAC)*D12*( ANC(1,IC))
            ENDIF  
            GO TO 10
          ENDIF   
        ENDIF
        QDOLD = QD   
      END DO
C
 10   IF(LDBG) THEN
        IEL = 1
        WRITE(*,*) 'Element 1 Stag @ IC,X,Y ',
     &             ICSTG(IEL),XSTG(IEL),YSTG(IEL)
        IEL = 2
        WRITE(*,*) 'Element 2 Stag @ IC,X,Y ',
     &             ICSTG(IEL),XSTG(IEL),YSTG(IEL)
      ENDIF
C  
      RETURN
      END


