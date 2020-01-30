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
C=========================================================================


      SUBROUTINE GAMSOLV
C------------------------------------------------------------------
C     Generates inviscid panel solution with current flow condition
C------------------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(LDBG) THEN
       WRITE(*,*)  'Entering GAMSOLV'
       WRITE(*,*)  ' LNCVP ',LNCVP
       WRITE(*,*)  ' LQAIC ',LQAIC
       WRITE(*,*)  ' LSYSP ',LSYSP
       WRITE(*,*)  ' LGSYS ',LGSYS
       WRITE(*,*)  ' LGAMA ',LGAMA
       WRITE(*,*)  ' LVMAV ',LVMAV
      ENDIF
C
C---- set control points and assign their pointers
      IF(.NOT.LNCVP) CALL CVPGEN
C
C---- set AIC matrix for velocities at control points
      IF(.NOT.LQAIC) CALL QAIC(.FALSE.)
C
C---- set up pointers for direct problem
      IF(.NOT.LSYSP) CALL SYSP
C
C---- set up and factor system for direct problem
      IF(.NOT.LGSYS) THEN
       CALL GSYS(.TRUE.,.TRUE.,.TRUE., 1,NPTOT, 1,0)
       IF(LDBG) WRITE(*,*) 'Factoring system matrix...'
       CALL LUDCMP(NSYSX,NSYS,SYS(1,1))
      ENDIF
C
C---- Set rotor drag source strengths from current velocities
      CALL SETROTORSRC
C---- Set (known) drag object source strengths from current velocities
      CALL SETDRGOBJSRC
C
C---- Initialize VMAVG for wake gamma relaxatio
      IF(.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
C
C---- Solve for current RHS knowns (Qinf,sources,GTH)
      IF(.NOT.LGAMA) THEN
       CALL GSOLVE
      ENDIF
C
C---- Calculate velocities at control points
      CALL QCSUM
C---- Calculate Cp's on surfaces, forces
      CALL QCPFOR
C---- Put flow data into wake grid
      CALL SETGRDFLW
C---- Find stagnation points
      CALL STGFIND
C
cc      write(*,*) 'Qn =', (QNDOF(IEL), IEL=1, 2)
      
      RETURN
      END ! GAMSOLV



      SUBROUTINE GAMSOL
C------------------------------------------------------------------
C     Generates inviscid panel solution by superposing unit-strength
C     solutions
C--------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LBGAM
C
C---- set control points and assign their pointers
      IF(.NOT.LNCVP) CALL CVPGEN
C
C---- set AIC matrix for velocities at control points
      IF(.NOT.LQAIC) CALL QAIC(.FALSE.)
C
C---- set up pointers for direct problem
      IF(.NOT.LSYSP) CALL SYSP
C
C---- set up and factor system for direct problem
      IF(.NOT.LGSYS) THEN
       CALL GSYS(.TRUE.,.TRUE.,.TRUE., 1,NPTOT, 1,0)
       IF(LDBG) WRITE(*,*) 'Factoring system matrix...'
       CALL LUDCMP(NSYSX,NSYS,SYS(1,1))
      ENDIF
C
      IF(.NOT.LGAMU) THEN
C
C----- all r.h.s. vectors are to be considered invalid...
       IF(LDBG) WRITE(*,*) 'Clearing unit-strength solutions...'
       NU = 0
       IUQINF = 0
       DO IR = 1, NRP
         IUVWK(IR) = 0
       ENDDO
       DO IP = 1, NPTOT
         IUGAM(IP) = 0
         IUSIG(IP) = 0
       ENDDO
       DO IU = 1, NUX
         LUSET(IU) = .FALSE.
       ENDDO
       LGAMU = .TRUE.
C
      ENDIF
C
C---- New unit solution for Qinf
      IF(QINF.NE.0.0) THEN
       IF(LDBG) WRITE(*,*) 'Qinf unit calc'
       CALL GUCALC(.TRUE.,.FALSE., 0,-1)
      ENDIF
C
C---- New unit solutions for B*Circ rotor circulation distribution
c      LBGAM = .FALSE.
c      DO IR = 1, NRP
c        LBGAM = (LBGAM .OR. WBGAM(IR).NE.0.0)
c      END DO   
c      IF(LBGAM) THEN
c        IF(LDBG) WRITE(*,*) 'Rotor circulation unit calc'
c        CALL GUCALC(.FALSE.,.TRUE., 0,-1)
c      ENDIF
C
C---- New unit solutions for RHS singularities
      DO IEL = 1, NEL
        IF(GAMSET(IEL).NE.0.0 .OR.
     &     SIGSET(IEL).NE.0.0 .OR.
     &     NETYPE(IEL).EQ.5   .OR. 
     &     NETYPE(IEL).EQ.6     ) THEN
         IP1 = IPFRST(IEL)
         IP2 = IPLAST(IEL)
         IF(LDBG) WRITE(*,*) 'Element set unit calc IEL ',IEL
         CALL GUCALC(.FALSE.,.FALSE., IP1,IP2)
        ENDIF
      ENDDO
C
C---- no imposed viscous strengths
      DO IP = 1, NPTOT
        GAMVSP(IP) = 0.
        SIGVSP(IP) = 0.
      ENDDO
C
C---- Set rotor drag source strengths from current velocities
      CALL SETROTORSRC
C---- Set (known) drag object source strengths from current velocities
      CALL SETDRGOBJSRC
C
C---- Initialize VMAVG for wake gamma relaxatio
      IF(.NOT.LVMAV) CALL VMAVGINIT(VAVGINIT)
C
C---- combine unit vorticities weighted by current singularities and freestream
      IF(.NOT.LGAMA) THEN
       CALL GUSUM
      ENDIF
cccHHY
C---- Calculate velocities at control points
      CALL QCSUM
cccHHY
C
C---- Calculate Cp's on surfaces, forces
      CALL QCPFOR
C
C---- Put flow data into wake grid
      CALL SETGRDFLW
C---- Find stagnation points
      CALL STGFIND
C
cc      write(*,*) 'Qn =', (QNDOF(IEL), IEL=1, 2)
      
      RETURN
      END ! GAMSOL

