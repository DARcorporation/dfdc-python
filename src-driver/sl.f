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


      SUBROUTINE GETSL(XS,YS,DS,DIR,
     &                 XBOX,YBOX,QSTAG,
     &                 LPRNT, LUPRNT,
     &                 NPSLDIM,NPSL,
     &                 XSL,YSL,USL,VSL,CPSL)
C-----------------------------------------------
C     Trace a streamline in direction IDIR from
C     point XS,YS. 
C     Streamline points and flow is added to streamline arrays
C-----------------------------------------------
cc      INCLUDE 'DFDC.INC'
      PARAMETER (PI=3.1415926535897932384)
      PARAMETER (DTR=PI/180.)
C
      DIMENSION XSL(*),YSL(*),USL(*),VSL(*),CPSL(*)
      LOGICAL LPRNT, LBOTH
C
C---- local arrays for calling CPCALC
      DIMENSION QF(2), CPF_QF(2)
C
C---- flag for printout of flow data
cc      LPRNT = .FALSE.
C---- speed below this is treated as stagnation
cc      QSTAG = 0.001*QREF
C
C---- Maximum # of SL points 
      MAXPTS = NPSLDIM
C---- minimum dot product between velocity unit vectors used for step control
      DPMIN = 0.99
C---- direction for integration (1.0 forward, -1.0 backward)
      DIR = 1.0
C---- flag for trace in both directions
      LBOTH = .FALSE.
      IF(DIR.EQ.0.0) LBOTH = .TRUE.
C
      NSL  = 0
      NSPT = 0
      NSL = NSL + 1
C
 10   XX  = XS
      YY  = YS
      SS  = DS*DIR
C
 20   CALL GETUV(XX,YY,US,VS)
      UVMAG = SQRT(US*US + VS*VS)
C      
      IF(UVMAG.GT.QSTAG) THEN 
        UUS = US/UVMAG
        VVS = VS/UVMAG
      ELSE
        WRITE(*,*) 'Stagnation point found, streamline ends'
        GO TO 10
      ENDIF
      NSPT = NSPT + 1
      IF(NSPT.GE.MAXPTS) THEN
        WRITE(*,*) 'Too many points, streamline ends'
        GO TO 10        
      ENDIF      
C
        IF(NSPT.GT.1) THEN
C---- Check dot product for new velocity, if too low subdivide interval 
          DOTP = UUOLD*UUS + VVOLD*VVS
          IF(DOTP .LT. DPMIN) THEN
            SS = 0.5*SS
            XX = XOLD + UUOLD*SS
            YY = YOLD + VVOLD*SS
            GO TO 20
          ENDIF
C---- get flow data at point
          QF(1) = US
          QF(2) = VS
          NF = 1
          CALL CPCALC(NF, QF, QINF,QREF, CPF, CPF_QF,CPF_QINF)
C----- add rotor enthalpy, entropy, pressure
          CALL CPADDHS(XX,YY,CPF,DHF,DSF,VTF)
C
C---- report flow data at point
          NPSL = NPSL + 1
          USL(NPSL)  = US
          VSL(NPSL)  = VS
          CPSL(NPSL) = VS
C
          IF(LPRNT) THEN
            ALF  = ATAN2(VS,US)
            SWRL = ATAN2(VTF,US)
            WRITE(*,1800) NSPT, XX, YY,
     &                    US, VS, VTF, UVMAG,
     &                    US/QREF, VS/QREF, VTF/QREF, UVMAG/QREF,
     &                    CPF, ALF/DTR, SWRL/DTR,
     &                    DHF, DSF 
          ENDIF
C---- reset step distance and try new point
          SS = DS*DIR
        ENDIF
C
        UUOLD = UUS
        VVOLD = VVS
        XOLD = XX
        YOLD = YY
C------ see if SL point is out of integration box
        IF(XX.LT.XBOX(1) .OR.
     &     XX.GT.XBOX(2) .OR.
     &     YY.LT.YBOX(1) .OR.
     &     YY.GT.YBOX(2) ) THEN
C------ trace other direction if LBOTH is set
          IF(LBOTH) THEN
            XX  = XS
            YY  = YS
            DIR = -DIR            
            SS  = DS*DIR
            LBOTH = .FALSE.
            NSPT = 0
            GO TO 20
          ENDIF
C-------- terminate SL trajectory integration
          GO TO 10
        ENDIF
C
        XX = XOLD + UUS*SS
        YY = YOLD + VVS*SS
        GO TO 20
C
 1800   FORMAT(/' #SL points ',I4,
     &         /' u     =', G11.5,'  v     =', G11.5,
     &         '  w     =', G11.5,'  q     =', G11.5,
     &         /' u/Qref=', G11.5,'  v/Qref=', G11.5,
     &         '  w/Qref=', G11.5,'  q/Qref=', G11.5,
     &         /' Cp    =', G11.5,'  alfa  =', G11.5,' deg',
     &                            '  swirl =', G11.5,' deg',
     &         /' delH  =', G11.5,'  delS  =', G11.5/ )
C
      END


      SUBROUTINE GETSLC(XS,YS,DS)
C-----------------------------------------------
C     Trace a streamline in direction IDIR from
C     point X0,Y0. 
C     Streamline is added to streamline list
C-----------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      CHARACTER*1 CHKEY
      LOGICAL LPRNT, LBOTH
C
C---- local arrays for calling CPCALC
      DIMENSION QF(2), CPF_QF(2)
C
c      IF(NSL.GE.NSLX) THEN
c       WRITE(*,*) 'Too many streamlines...'
c       RETURN
c      ENDIF
C
      IF(DS.LE.0.0) THEN
        WRITE(*,*) 'Enter delS for streamline step size'
        READ(*,*) DS
      ENDIF
C
C---- speed below this is treated as stagnation
      QSTAG = 0.001*QREF
C---- Maximum # of SL points 
      MAXPTS = 500
C---- minimum dot product between velocity unit vectors used for step control
      DPMIN = 0.99
C---- direction for integration (1.0 forward, -1.0 backward)
      DIR = 1.0
C---- flag for printout of flow data
      LPRNT = .FALSE.
C---- flag for trace in both directions
      LBOTH = .FALSE.
C
      NSL  = 0
      NSPT = 0
C
C---- Interactive input, set up offset and scale for type of plot 
      IF(IPTYPE.EQ.JCPLT .OR. IPTYPE.EQ.JQPLT) THEN
        XOFF = XOFA
        YOFF = YOFA
        FAC  = FACA
      ELSE
        XOFF = -XOFG
        YOFF = -YOFG
        FAC  = GSF
      ENDIF
      WRITE(*,*) 'Streamline tracing'
      WRITE(*,*) '  type U for upstream'
      WRITE(*,*) '  type D for downstream (default)'
      WRITE(*,*) '  type B for upstream/downstream'
C
 10   WRITE(*,*) 'Click on points (type "Q" to finish)...'
C---- interactive input of points with cursor (on plot)
      CALL GETCURSORXY(XPLT,YPLT,CHKEY)
C      
      IF(INDEX('Qq',CHKEY).NE.0) RETURN
      IF(INDEX('Uu',CHKEY).NE.0) DIR = -1.0
      IF(INDEX('Dd',CHKEY).NE.0) DIR =  1.0
      IF(INDEX('Bb',CHKEY).NE.0) THEN
        DIR =  1.0
        LBOTH = .TRUE.
      ENDIF
C
C------- undo plot offset/scaling
      XS = XPLT/FAC - XOFF
      YS = YPLT/FAC - YOFF
C
C------ add point to current geometry plot
      CALL GETCOLOR(ICOL0)
      CALL NEWCOLORNAME('VIOLET')
C------ put a "+" for start point on selected spot
      XPLT0 = (XS + XOFF)*FAC
      YPLT0 = (YS + YOFF)*FAC
      SSIZ = 0.8*SHGT
      CALL PLSYMB(XPLT0,YPLT0,SSIZ,3,0.0,0)
      XX  = XS
      YY  = YS
      SS  = DS*DIR
      NSPT = 0
      NSL = NSL + 1
C
 20   CALL GETUV(XX,YY,US,VS)
      UVMAG = SQRT(US*US + VS*VS)
C      
      IF(UVMAG.GT.QSTAG) THEN 
        UUS = US/UVMAG
        VVS = VS/UVMAG
      ELSE
        WRITE(*,*) 'Stagnation point found, streamline ends'
        CALL PLFLUSH
        GO TO 10
      ENDIF
      NSPT = NSPT + 1
      IF(NSPT.GE.MAXPTS) THEN
        WRITE(*,*) 'Too many points, streamline ends'
        CALL PLFLUSH
        GO TO 10        
      ENDIF      
C
        IF(NSPT.GT.1) THEN
C---- Check dot product for new velocity, if too low subdivide interval 
          DOTP = UUOLD*UUS + VVOLD*VVS
          IF(DOTP .LT. DPMIN) THEN
            SS = 0.5*SS
            XX = XOLD + UUOLD*SS
            YY = YOLD + VVOLD*SS
            GO TO 20
          ENDIF
C---- OK, plot new streamline point 
          XPLT = (XX + XOFF)*FAC
          YPLT = (YY + YOFF)*FAC
          CALL PLOT(XPLT,YPLT,2)
C---- report flow data at point
          QF(1) = US
          QF(2) = VS
          NF = 1
          CALL CPCALC(NF, QF, QINF,QREF, CPF, CPF_QF,CPF_QINF)
C----- add rotor enthalpy, entropy, pressure
          CALL CPADDHS(XX,YY,CPF,DHF,DSF,VTF)
          ALF = ATAN2(VS,US)/DTR
          IF(LPRNT .OR. NSPT.EQ.1) THEN
            WRITE(*,1800) NSPT, XX, YY,
     &                    US/QREF, VS/QREF, VTF/QREF,
     &                    UVMAG/QREF, CPF, ALF, DHF, DSF 
          ENDIF
C---- reset step distance and try new point
          SS = DS*DIR
        ENDIF
C
        UUOLD = UUS
        VVOLD = VVS
        XOLD = XX
        YOLD = YY
C------ see if SL point is out of integration box
        IF(XX.LT.XWBOX(1) .OR.
     &     XX.GT.XWBOX(2) .OR.
     &     YY.LT.YWBOX(1) .OR.
     &     YY.GT.YWBOX(2) ) THEN
C------ trace other direction if LBOTH is set
          IF(LBOTH) THEN
            XX  = XS
            YY  = YS
            XPLT = (XX + XOFF)*FAC
            YPLT = (YY + YOFF)*FAC
            CALL PLOT(XPLT,YPLT,3)
            DIR = -DIR            
            SS  = DS*DIR
            LBOTH = .FALSE.
            NSPT = 0
            GO TO 20
          ENDIF
C-------- terminate SL trajectory integration
          CALL PLFLUSH
          GO TO 10
        ENDIF
C
        XX = XOLD + UUS*SS
        YY = YOLD + VVS*SS
        GO TO 20
C
 1800   FORMAT(/' #SL points ',I4,
     &         /' x      =', G13.5,'   y      =', G13.5,
     &         /' u/Qref =', G13.5,'   v/Qref =', G13.5,
     &                             '   w/Qref =', G13.5,
     &         /' q/Qref =', G13.5,'   Cp     =', G13.5,
     &                             '   alfa   =', G13.5,' deg',
     &         /'  delH  =', G13.5,'   delS   =', G13.5 )
C
C------ put Cp point symbol on Cp(x) plot
c        CALL XYSYMB(1,XFVELS,CPF,-XOFF,FAC,0.0,-CPLFAC,SHGT,0)
c        CALL NEWCOLOR(ICOL0)
C
      END

