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

C---------------------
C  Pointers for referencing polar force coefficients and parameters
C  First 4 pointers must be main polar plot variables.
C
      PARAMETER (
     &    IPAL =  1,    ! alpha
     &    IPCL =  2,    ! CL
     &    IPCD =  3,    ! CD
     &    IPCM =  4,    ! CM
     &    IPCV =  5,    ! CDvisc
     &    IPMA =  7,    ! Mach
     &    IPQI = 10  )  ! Qinf
      PARAMETER (IPTOT=10)

C---------------------
C  Pointers for referencing polar airfoil-element quantities
C
      PARAMETER (
     &    JPCX = 1,    ! Cx
     &    JPCY = 2,    ! Cy
     &    JPCM = 3,    ! Cm
     &    JPGAM = 4,   ! gamma (prescribed)
     &    JPSIG = 5,   ! sigma (prescribed)
     &    JPSOUR = 6,  ! Sigma 
     &    JPCIRC = 7,  ! Gamma
     &    JPDBLX = 8,  ! Doubx
     &    JPDBLY = 9 ) ! Douby
      PARAMETER (JPTOT=9)

C---------------------
C  Pointers for referencing polar airfoil-side quantities
C
      PARAMETER (
     &    KPTP = 1,     ! trip
     &    KPTN = 2 )    ! transition
      PARAMETER (KPTOT=2)

C------------------------------------------------------------------
C     Plot type
C
C     IPTYPE     Plot type takes values (0 for no plot) 
      PARAMETER (JCPLT = 1)  ! Cp plot
      PARAMETER (JQPLT = 2)  ! Q  plot
      PARAMETER (JAPLT = 3)  ! Airfoil plot (with vectors)
