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


      SUBROUTINE GUIBOX(K, X1,X2,Y1,Y2, COLOR, LABEL)
      CHARACTER*(*) COLOR, LABEL
C----------------------------------------------------------
C     Plots a GUI-button box with label string.
C     Places the box coordinates into the COM_GUI 
C     arrays associated with the button index K.  
C     FUNCTION LGUI can then determine if a cursor 
C     falls within box K.
C----------------------------------------------------------
      COMMON /COM_GUI/ XGUI(2,20), YGUI(2,20)
C
      IF(K.LT.1 .OR. K.GT.20) RETURN
C
      CALL GETORIGIN(XORG,YORG)
      CALL GETFACTORS(XSCALE,YSCALE)
C
      CALL GETCOLOR(ICOL0)
      CALL NEWCOLORNAME(COLOR)
C
C---- set GUI window
      XGUI(1,K) = (X1 - XORG)/XSCALE
      XGUI(2,K) = (X2 - XORG)/XSCALE
      YGUI(1,K) = (Y1 - YORG)/YSCALE
      YGUI(2,K) = (Y2 - YORG)/YSCALE
C
C---- plot GUI window
      CALL PLOT(XGUI(1,K),YGUI(1,K),3)
      CALL PLOT(XGUI(2,K),YGUI(1,K),2)
      CALL PLOT(XGUI(2,K),YGUI(2,K),2)
      CALL PLOT(XGUI(1,K),YGUI(2,K),2)
      CALL PLOT(XGUI(1,K),YGUI(1,K),2)
C
      NL = LEN(LABEL)
      CHA = MIN( (XGUI(2,K)-XGUI(1,K))/FLOAT(NL+1), 
     &           (YGUI(2,K)-YGUI(1,K))/1.8  )
      XCA = 0.5*(XGUI(2,K)+XGUI(1,K)) - 0.5*CHA*FLOAT(NL) + 0.2*CHA
      YCA = 0.5*(YGUI(2,K)+YGUI(1,K)) - 0.6*CHA
      CALL PLCHAR(XCA,YCA,CHA,LABEL,0.0,NL)
C
      CALL NEWCOLOR(ICOL0)
      RETURN
      END ! GUIBOX



      LOGICAL FUNCTION LGUI(K,XC,YC)
C-----------------------------------------------
C     Returns T if location XC,YC falls within 
C     the GUI(K) window defined in GUIBOX.
C-----------------------------------------------
      COMMON /COM_GUI/ XGUI(2,20), YGUI(2,20)
C
      LGUI = .FALSE.
      IF(K.LT.1 .OR. K.GT.20) RETURN
C
      LGUI = XC .GT. XGUI(1,K) .AND. 
     &       XC .LE. XGUI(2,K) .AND.
     &       YC .GT. YGUI(1,K) .AND.
     &       YC .LE. YGUI(2,K)
C
      RETURN
      END ! LGUI
