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

      SUBROUTINE ATMO(ALSPEC,DELTAT,VSOALT,RHOALT,RMUALT)
C---------------------------------------------------------
C     Returns speed of sound (VSO) in m/s, density (RHO)
C     in kg/m^3, and dynamic viscosity (RMU) in kg/m-s
C     of standard atmosphere at specified altitude ALSPEC
C     (in kilometers).  If ALSPEC=-1, water properties
C     at 15 Celsius are returned.
C
C     Reference:  "U.S. Standard Atmosphere", NOAA.
C---------------------------------------------------------
      LOGICAL FIRST
C
      PARAMETER ( N = 44 )
      REAL ALT(N), VSO(N), RHO(N), RMU(N)
C
      DATA FIRST / .TRUE. /
      DATA ALT
     &   / 0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
     &    10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
     &    20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
     &    30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 
     &    40.0, 45.0, 60.0, 75.0 /
      DATA VSO
     & / 340.0,336.0,332.0,329.0,325.0,320.0,316.0,312.0,308.0,304.0,
     &   299.0,295.0,295.0,295.0,295.0,295.0,295.0,295.0,295.0,295.0,
     &   295.0,295.8,296.4,297.1,297.8,298.5,299.1,299.8,300.5,301.1,
     &   301.8,302.5,303.1,305.0,306.8,308.7,310.5,312.3,314.0,316.0,
     &   318.0,355.0,372.0,325.0 /
      DATA RHO
     & / 1.226,1.112,1.007,0.909,0.820,0.737,0.660,0.589,0.526,0.467,
     &   0.413,0.364,0.311,0.265,0.227,0.194,0.163,0.141,0.121,0.103,
     &   .0880,.0749,.0637,.0543,.0463,.0395,.0338,.0288,.0246,.0210,
     &   .0180,.0154,.0132,.0113,.0096,.0082,.0070,.0060,.0052,.0044,
     &   0.004,0.002,3.9E-4,8.0E-5 /
      DATA RMU
     & / 1.780,1.749,1.717,1.684,1.652,1.619,1.586,1.552,1.517,1.482,
     &   1.447,1.418,1.418,1.418,1.418,1.418,1.418,1.418,1.418,1.418,
     &   1.418,1.427,1.433,1.438,1.444,1.449,1.454,1.460,1.465,1.471,
     &   1.476,1.481,1.487,1.502,1.512,1.532,1.546,1.561,1.580,1.600,
     &   1.700,1.912,2.047,1.667 /
C
C---- special case: Water at STP
      IF(ALSPEC.EQ.-1.0) THEN
       VSOALT = 1500.
       RHOALT = 1000.
       RMUALT = 1.15E-3
       WRITE(*,*) '                              o        '
       WRITE(*,*) 'ATMO: You are underwater at 15  Celsius'
       RETURN
      ENDIF
C
C---- linearly interpolate quantities from tabulated values
      DO 10 I=2, N
        IF(ALSPEC.GT.ALT(I)) GO TO 10
C
         DALT = ALT(I) - ALT(I-1)
         DVSO = VSO(I) - VSO(I-1)
         DRHO = RHO(I) - RHO(I-1)
         DRMU = RMU(I) - RMU(I-1)
C
         ALFRAC = (ALSPEC - ALT(I-1)) / DALT
C
         VSOALT = VSO(I-1) + DVSO*ALFRAC
         RHOALT = RHO(I-1) + DRHO*ALFRAC
         RMUALT = RMU(I-1) + DRMU*ALFRAC
         RMUALT = RMUALT * 1.0E-5
C
         RETURN
   10 CONTINUE
C
C
      IF(ALSPEC.GT.ALT(N)) THEN
       WRITE(*,*) ' '
       WRITE(*,*) 'ATMO: You''re in low earth orbit.  Good luck.'
       VSOALT = VSO(N)
       RHOALT = RHO(N)
       RMUALT = RMU(N) * 1.0E-5
       RETURN
      ENDIF
C
c      IF(FIRST) THEN
c       DO 20 I=1, N
c         RHO(I) = ALOG(RHO(I))
c 20    CONTINUE
c       CALL SPLINE(VSO,VSOH,ALT,N)
c       CALL SPLIND(RHO,RHOH,ALT,N,999.0,0.0)
c       CALL SPLINE(RMU,RMUH,ALT,N)
c       FIRST = .FALSE.
c      ENDIF
cC
cC---- interpolate quantities from splines
c      VSOALT = SEVAL(ALSPEC,VSO,VSOH,ALT,N)
c      RHOALT = SEVAL(ALSPEC,RHO,RHOH,ALT,N)
c      RMUALT = SEVAL(ALSPEC,RMU,RMUH,ALT,N) * 1.0E-5
c      RHOALT = EXP(RHOALT)
cC
      RETURN
      END ! ATMO


      SUBROUTINE FLOSHO(LU, VSO, RHO, RMU)
      WRITE(LU,10) VSO, RHO, RMU
   10 FORMAT(/' Speed of sound (m/s)  :',F10.3
     &       /' Air density   (kg/m^3):',F10.5
     &       /' Air viscosity (kg/m-s):',E11.4 )
      RETURN
      END

