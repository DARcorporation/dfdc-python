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

      SUBROUTINE MINMAX(N,X,XMIN,XMAX)
C--------------------------------------
C     Calculates min,max of array X 
C--------------------------------------
      DIMENSION X(*)
C
      I = 1
      XMIN = X(I)
      XMAX = X(I)
      DO I = 2, N
        XMIN = MIN(XMIN,X(I))
        XMAX = MAX(XMAX,X(I))
      ENDDO
C
      RETURN
      END

 
      SUBROUTINE NCALC(X,Y,S,N,XN,YN)
C---------------------------------------
C     Calculates normal unit vector
C     components at airfoil panel nodes
C---------------------------------------
      DIMENSION X(*), Y(*), S(*), XN(*), YN(*)
C
      IF(N.LE.1) RETURN
C
      CALL SEGSPL(X,XN,S,N)
      CALL SEGSPL(Y,YN,S,N)
      DO 10 I=1, N
        SX =  YN(I)
        SY = -XN(I)
        SMOD = SQRT(SX*SX + SY*SY)
        XN(I) = SX/SMOD
        YN(I) = SY/SMOD
   10 CONTINUE
C
C---- average normal vectors at corner points
      DO 20 I=1, N-1
        IF(S(I) .EQ. S(I+1)) THEN
          SX = 0.5*(XN(I) + XN(I+1))
          SY = 0.5*(YN(I) + YN(I+1))
          SMOD = SQRT(SX*SX + SY*SY)
          XN(I)   = SX/SMOD
          YN(I)   = SY/SMOD
          XN(I+1) = SX/SMOD
          YN(I+1) = SY/SMOD
        ENDIF
 20   CONTINUE
C
      RETURN
      END



      SUBROUTINE LEFIND(SLE,X,XP,Y,YP,S,N)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C------------------------------------------------------
C     Locates leading edge spline-parameter value SLE
C
C     The defining condition is
C         
C      (X-XTE,Y-YTE) . (X',Y') = 0     at  S = SLE
C
C     i.e. the surface tangent is normal to the chord
C     line connecting X(SLE),Y(SLE) and the TE point.
C------------------------------------------------------
      IF(N.LE.1) THEN
       SLE = S(1)
       RETURN
      ENDIF
C
C---- convergence tolerance
      DSEPS = (S(N)-S(1)) * 1.0E-5
C
C---- set trailing edge point coordinates
      XTE = 0.5*(X(1) + X(N))
      YTE = 0.5*(Y(1) + Y(N))
C
C---- get first guess for SLE
      DO 10 I=3, N-2
        DXTE = X(I) - XTE
        DYTE = Y(I) - YTE
        DX = X(I+1) - X(I)
        DY = Y(I+1) - Y(I)
        DOTP = DXTE*DX + DYTE*DY
        IF(DOTP .LT. 0.0) GO TO 11
   10 CONTINUE
C
   11 SLE = S(I)
C
C---- check for sharp LE case
      IF(S(I) .EQ. S(I-1)) RETURN
C
C---- Newton iteration to get exact SLE value
      DO 20 ITER=1, 50
        XLE  = SEVAL(SLE,X,XP,S,N)
        YLE  = SEVAL(SLE,Y,YP,S,N)
        DXDS = DEVAL(SLE,X,XP,S,N)
        DYDS = DEVAL(SLE,Y,YP,S,N)
        DXDD = D2VAL(SLE,X,XP,S,N)
        DYDD = D2VAL(SLE,Y,YP,S,N)
C
        XCHORD = XLE - XTE
        YCHORD = YLE - YTE
C
C------ drive dot product between chord line and LE tangent to zero
        RES  = XCHORD*DXDS + YCHORD*DYDS
        RESS = DXDS  *DXDS + DYDS  *DYDS
     &       + XCHORD*DXDD + YCHORD*DYDD
C
        IF(RESS.EQ.0.0) GO TO 21
C
C------ Newton delta for SLE 
        DSLE = -RES/RESS
C
        DSLE = MAX( DSLE , -0.02*ABS(XCHORD+YCHORD) )
        DSLE = MIN( DSLE ,  0.02*ABS(XCHORD+YCHORD) )
        SLE = SLE + DSLE
        IF(ABS(DSLE) .LT. DSEPS) RETURN
 20   CONTINUE
 21   WRITE(*,*) 'LEFIND:  LE point not found.  Continuing...'
C
 25   SLE = S(I)
      RETURN
      END




      SUBROUTINE CNORM(X,XP,Y,YP,S,N, XTRAN,YTRAN,SCALE,ANGLE)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C----------------------------------------------------------
C     Scales and rotates coordinates to get unit chordline
C----------------------------------------------------------
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
C
      XTE = 0.5*(X(1) + X(N))
      YTE = 0.5*(Y(1) + Y(N))
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
C
      CX = XTE - XLE
      CY = YTE - YLE
      CSSQ = CX*CX + CY*CY
C
      IF(CSSQ.EQ.0.0) THEN
       CSQINV = 1.0
       CINV = 1.0
      ELSE
       CSQINV = 1.0/CSSQ
       CINV = 1.0/SQRT(CSSQ)
      ENDIF
C
      DO I=1, N
        XBAR = X(I) - XLE
        YBAR = Y(I) - YLE
        X(I) = ( XBAR*CX + YBAR*CY)*CSQINV
        Y(I) = (-XBAR*CY + YBAR*CX)*CSQINV
        S(I) = S(I)*CINV
      ENDDO
C
      XTRAN = -XLE
      YTRAN = -YLE
      SCALE = CINV
      ANGLE = ATAN2(CY,CX)
C
      RETURN
      END ! CNORM


      SUBROUTINE GEOPAR(X,XP,Y,YP,S,N, T,
     &             SLE,CHORD,RADLE,ANGTE,
     &             AREA,EIXXA,EIYYA,EIXYA,
     &             ASKN,EIXXT,EIYYT,EIXYT,
     &             THICK,CAMBR,
     &             VOLM,VSKN,ASRF,RGXV,RGYV)
      DIMENSION X(*), XP(*), Y(*), YP(*), S(*), T(*)
C------------------------------------------------------
C     Sets geometric parameters for airfoil shape
C------------------------------------------------------
C
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
C
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
C
      CHSQ = (XTE-XLE)**2 + (YTE-YLE)**2
      CHORD = SQRT(CHSQ)
C
      CURVLE = CURV(SLE,X,XP,Y,YP,S,N)
C
      IF(ABS(CURVLE) .GT. 0.001*(S(N)-S(1))) THEN
       RADLE = 1.0 / CURVLE
      ELSE
       RADLE = 0.
      ENDIF
C
      ANG1 = ATAN2( -YP(1) , -XP(1) )
      ANG2 = ATANC(  YP(N) ,  XP(N) , ANG1 )
      ANGTE = ANG2 - ANG1
C
C---- 2D properties
      CALL AECALC(N,X,Y,T, 0, PERIM,AREA,XCENA,YCENA,EIXXA,EIYYA,EIXYA)
      CALL AECALC(N,X,Y,T, 2, PERIM,ASKN,XCENT,YCENT,EIXXT,EIYYT,EIXYT)
      TFAC = 1.0
      CFAC = 1.0
      CALL TCSET(X,XP,Y,YP,S,N,
     &           THICK,CAMBR, TFAC,CFAC, XNEW,YNEW, .FALSE.)
C
C---- Axisymmetric properties
      CALL AXCALC(N,X,Y,T, 0, VOLM,ASRF,XCENV, YCENV, RGXV,RGYV)
      CALL AXCALC(N,X,Y,T, 2, VSKN,ASRF,XCENVT,YCENVT,RGXT,RGYT)
C
      RETURN
      END


      SUBROUTINE AXCALC(N,X,Y,T, ITYPE, 
     &                  VOLM,AREA,XCEN,YCEN,RGX,RGY)
      DIMENSION X(*),Y(*),T(*)
C---------------------------------------------------------------
C     Calculates geometric properties of axisymmetric shape X,Y
C     where X is axial coordinate and Y is radial coordinate
C
C     Input:
C       N      number of points
C       X(.)   shape coordinate point arrays
C       Y(.)
C       T(.)   skin-thickness array, used only if ITYPE = 2
C       ITYPE  = 0 ...   integration is over volume 
C              = 2 ...   integration is over skin with thickness specified
C
C     Output:
C       VOLM  volume 
C       AREA  surface area
C       XCEN  centroid location
C       YCEN  
C       RGX   radius of gyration about X centroid
C       RGY   radius of gyration about Y centroid (Y=0)
C---------------------------------------------------------------
      LOGICAL LAXIBOD, LOPEN
      DATA PI / 3.141592653589793238 /
C
C--- check for shape symmetric about X axis
C    if first or last point is on centerline it is an axisymmetric body
      LAXIBOD = (Y(1).EQ.0.0 .OR. Y(N).EQ.0.0)
C
C--- check shape for "open" or closed body by comparing distance between
C    first and last points to 0.1 of perimeter length
C
      SINT = 0.0
      DO IO = 2, N
        DX =  X(IO) - X(IO-1)
        DY =  Y(IO) - Y(IO-1)
        DS = SQRT(DX*DX + DY*DY)
        SINT = SINT + DS
      END DO
      DX   = X(N) - X(1)
      DY   = Y(N) - Y(1)
      DSTE = SQRT(DX*DX + DY*DY)
      LOPEN = (DSTE.GT.0.1*SINT)
C
C--- integrate around contour to get geometric quantities
      VINT  = 0.0
      AINT  = 0.0
      XINT  = 0.0
      YINT  = 0.0
      XXINT = 0.0
      YYINT = 0.0
C
      DO 10 IO = 1, N
        IF(IO.LT.N) THEN
          IP = IO + 1
        ELSE
C--- if last point is on centerline or body is open end here...
          IF(LAXIBOD .OR. LOPEN) THEN 
           GO TO 10
          ELSE
C--- otherwise close body with last segment between first and last points
           IP = 1
          ENDIF
        ENDIF
C
c        write(54,199) io,x(io),y(io),x(ip),y(ip)
        DX =  X(IO) - X(IP)
        DY =  Y(IO) - Y(IP)
        XA = (X(IO) + X(IP))*0.5
        YA = (Y(IO) + Y(IP))*0.5
        DS = SQRT(DX*DX + DY*DY)
 199    format(I4,7f10.6)
C---- surface area
        DA = 2.0*PI*YA*DS
        AINT  = AINT + DA

        IF(ITYPE.EQ.0) THEN
C-------- treat as a axisymmetric solid body (YCEN=0,intXY=0)
          DV = PI*DX*YA**2
c        write(55,199) io,dx,dy,xa,ya,ds,da,dv
          VINT  = VINT  +       DV
          XINT  = XINT  + XA   *DV
          YINT  = 0.0
          XXINT = XXINT + XA*XA*DV
          YYINT = YYINT + YA*YA*DV/4.0
        ELSE
C--------- integrate over skin thickness
          TA = (T(IO) + T(IP))*0.50
          DV = 2.0*PI*YA*TA*DS
          VINT  = VINT  +       DV
          XINT  = XINT  + XA   *DV
          YINT  = 0.0
          XXINT = XXINT + XA*XA*DV
          YYINT = YYINT + YA*YA*DV
        ENDIF
 10   CONTINUE
C
C---- volume and surface area
      IF(LOPEN .AND. .NOT.LAXIBOD) VINT = 0.0
      VOLM = VINT
      AREA = AINT
C
C---- calculate centroid location
      IF(VINT .LE. 0.0) THEN
       XCEN = 0.
       YCEN = 0.
       RGX  = 0.
       RGY  = 0.
      ELSE
       XCEN = XINT/VINT
       YCEN = YINT/VINT
C---- calculate radii of gyration
       RGXSQ = XXINT/VINT - XCEN**2
       IF(RGXSQ.GE.0.0) RGX = SQRT(RGXSQ)
       RGY = SQRT(YYINT/VINT)
      ENDIF
C
      RETURN
      END ! AXCALC



      SUBROUTINE AECALC(N,X,Y,T, ITYPE, 
     &                  PERIM,AREA,XCEN,YCEN,EIXX,EIYY,EIXY)
      DIMENSION X(*),Y(*),T(*)
C---------------------------------------------------------------
C     Calculates geometric properties of 2D shape X,Y
C
C     Input:
C       N      number of points
C       X(.)   shape coordinate point arrays
C       Y(.)
C       T(.)   skin-thickness array, used only if ITYPE = 2
C       ITYPE  = 0 ...   integration is over whole area  dA = dx dy
C              = 1 ...   integration is over perimeter   dA =   ds
C              = 2 ...   integration is over skin area   dA = t ds
C
C     Output:
C       PERIM perimeter length
C       AREA  area (assumes 2D)
C       XCEN  centroid location
C       YCEN  
C       EIXX  moments of inertia
C       EIYY  
C       EIXY
C---------------------------------------------------------------
      LOGICAL LAXIS
C
C--- check for shape with symmetry about X axis
C    if first or last point is on centerline it is an axisymmetric body
      LAXIS = (Y(1).EQ.0.0 .OR. Y(N).EQ.0.0)
C
      SINT  = 0.0
      VINT  = 0.0
      AINT  = 0.0
      XINT  = 0.0
      YINT  = 0.0
      XXINT = 0.0
      XYINT = 0.0
      YYINT = 0.0
C
      DO 10 IO = 1, N
C
        IF(IO.LT.N) THEN
          IP = IO + 1
        ELSE
C--- if last point is on axis (symmetrical body)
          IF(LAXIS) THEN 
           GO TO 10
          ELSE
C--- otherwise close body with last segment between first and last points
           IP = 1
          ENDIF
        ENDIF
C
        DX =  X(IO) - X(IP)
        DY =  Y(IO) - Y(IP)
        XA = (X(IO) + X(IP))*0.5
        YA = (Y(IO) + Y(IP))*0.5
C
        DS = SQRT(DX*DX + DY*DY)
        SINT = SINT + DS

        IF(ITYPE.EQ.0) THEN
C-------- integrate over airfoil cross-section
          DA = YA*DX
          AINT  = AINT  +       DA
          XINT  = XINT  + XA   *DA
          YINT  = YINT  + YA   *DA/2.0
          XXINT = XXINT + XA*XA*DA
          XYINT = XYINT + XA*YA*DA/2.0
          YYINT = YYINT + YA*YA*DA/3.0
        ELSE
          IF(ITYPE.EQ.1) THEN
C--------- integrate over perimeter
           DA = DS
          ELSE
C--------- integrate over skin thickness
           TA = (T(IO) + T(IP))*0.50
           DA = TA*DS
          ENDIF
          AINT  = AINT  +       DA
          XINT  = XINT  + XA   *DA
          YINT  = YINT  + YA   *DA
          XXINT = XXINT + XA*XA*DA
          XYINT = XYINT + XA*YA*DA
          YYINT = YYINT + YA*YA*DA
        ENDIF
 10   CONTINUE
C
      PERIM = SINT
      AREA  = AINT
C
C---- calculate centroid location
      IF(AINT .EQ. 0.0) THEN
       XCEN = 0.
       YCEN = 0.
      ELSE
       XCEN = XINT/AINT
       YCEN = YINT/AINT
      ENDIF
C
C---- calculate inertias
      EIXX = YYINT - YCEN*YCEN*AINT
      EIYY = XXINT - XCEN*XCEN*AINT
      EIXY = XYINT - XCEN*YCEN*AINT
C
      RETURN
      END ! AECALC



      SUBROUTINE PRAXIS(FXX,FYY,FXY,F11,F22,AP1,AP2)
C-----------------------------------------------------------------
C     Sets prinicipal components F11,F22 of tensor FXX,FYY,FXY.
C     Also sets corresponding principal direction angls AP1,AP2.
C-----------------------------------------------------------------
C
      DATA PI / 3.141592653589793238 /
C
      FSQ  = 0.25*(FXX - FYY)**2  + FXY**2
      SGN = SIGN( 1.0 , FYY-FXX )
      F11 = 0.5*(FXX + FYY) - SGN*SQRT(FSQ)
      F22 = 0.5*(FXX + FYY) + SGN*SQRT(FSQ)
C
      IF(F11.EQ.0.0 .OR. F22.EQ.0.0) THEN
C----- vanishing tensor limit
       AP1 = 0.0
       AP2 = ATAN2(1.0,0.0)
C
      ELSEIF(FSQ/(F11*F22) .LT. 0.0001) THEN
C----- rotationally-invariant tensor (circle, square, etc.)
       AP1 = 0.0
       AP2 = ATAN2(1.0,0.0)
C
      ELSE
C----- normal general tensor case
       C1 = FXY
       S1 = FXX-F11
C
       C2 = FXY
       S2 = FXX-F22
C
       IF(ABS(S1).GT.ABS(S2)) THEN
         AP1 = ATAN2(S1,C1)
         AP2 = AP1 + 0.5*PI
       ELSE
         AP2 = ATAN2(S2,C2)
         AP1 = AP2 - 0.5*PI
       ENDIF

       IF(AP1.LT.-0.5*PI) AP1 = AP1 + PI
       IF(AP1.GT.+0.5*PI) AP1 = AP1 - PI
       IF(AP2.LT.-0.5*PI) AP2 = AP2 + PI
       IF(AP2.GT.+0.5*PI) AP2 = AP2 - PI
C
      ENDIF

      RETURN
      END



      SUBROUTINE TEGAP(X,XS,Y,YS,S,N, DOC,GAPNEW)
C----------------------------------
C     Used to set buffer airfoil 
C     trailing edge gap
C----------------------------------
      DIMENSION X(*),XS(*),Y(*),YS(*),S(*)
C
      CALL LEFIND(SLE,X,XS,Y,YS,S,N)
      XLE = SEVAL(SLE,X,XS,S,N)
      YLE = SEVAL(SLE,Y,YS,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHSQ = (XTE-XLE)**2 + (YTE-YLE)**2
C
      DXN = X(1) - X(N)
      DYN = Y(1) - Y(N)
      GAP = SQRT(DXN**2 + DYN**2)
C
C---- components of unit vector parallel to TE gap
      IF(GAP.GT.0.0) THEN
       DXU = DXN / GAP
       DYU = DYN / GAP
      ELSE
       DXU = -.5*(YS(N) - YS(1))
       DYU = 0.5*(XS(N) - XS(1))
      ENDIF
C
      DGAP = GAPNEW - GAP
C
C---- go over each point, changing the y-thickness appropriately
      DO I = 1, N
C------ chord-based x/c
        XOC = (  (X(I)-XLE)*(XTE-XLE)
     &         + (Y(I)-YLE)*(YTE-YLE) ) / CHSQ
C
C------ thickness factor tails off exponentially away from trailing edge
        IF(DOC .EQ. 0.0) THEN
          TFAC = 0.0
          IF(I.EQ.1 .OR. I.EQ.N) TFAC = 1.0
        ELSE
          ARG = MIN( (1.0-XOC)*(1.0/DOC-1.0) , 15.0 )
          TFAC = EXP(-ARG)
        ENDIF
C
        IF(S(I).LE.SLE) THEN
         X(I) = X(I) + 0.5*DGAP*XOC*TFAC*DXU
         Y(I) = Y(I) + 0.5*DGAP*XOC*TFAC*DYU
        ELSE
         X(I) = X(I) - 0.5*DGAP*XOC*TFAC*DXU
         Y(I) = Y(I) - 0.5*DGAP*XOC*TFAC*DYU
        ENDIF
      ENDDO
C
      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XS,S,N)
      CALL SEGSPL(Y,YS,S,N)
C
      RETURN
      END ! TEGAP



      SUBROUTINE TCSET(X,XP,Y,YP,S,N,
     &                 TMAX,CMAX, TFAC,CFAC, XNEW,YNEW, LNEWSET)
C-----------------------------------------------------
C     Determines max thickness and camber.
C     Scales thickness and camber by TFAC,CFAC
C-----------------------------------------------------
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
      DIMENSION XNEW(*), YNEW(*)
      LOGICAL LNEWSET
C
      DATA EPS / 1.0E-5 /
C
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHBSQ = (XTE-XLE)**2 + (YTE-YLE)**2
C
      IF(CHBSQ.EQ.0.0) THEN
       TMAX = 0.
       CMAX = 0.
       RETURN
      ENDIF
C
C---- set unit chord-line vector
      DXC = (XTE-XLE) / SQRT(CHBSQ)
      DYC = (YTE-YLE) / SQRT(CHBSQ)
C
      STOT = S(N) - S(1)
C
      TMAX = 0.
      CMAX = 0.
      CMIN = 0.
C
C---- go over each point, finding and/or changing the y-thickness
C-    (defined normal to chord line)
      DO 30 I = 1, N
C
C------ coordinates in chord-line axes
        XBAR = (X(I)-XLE)*DXC + (Y(I)-YLE)*DYC
        YBAR = (Y(I)-YLE)*DXC - (X(I)-XLE)*DYC
C
        SOPP = 2.0*SLE - S(I)
        SOPP = MAX(SOPP,S(1))
        SOPP = MIN(SOPP,S(N))
C
        IF(ABS(SOPP-SLE) .LE. EPS*STOT) THEN
         SOPP = SLE
        ELSE
C------- converge on exact opposite point with same chord x value
         DO 302 ITOPP = 1, 12
           XOPP  = SEVAL(SOPP,X,XP,S,N)
           YOPP  = SEVAL(SOPP,Y,YP,S,N)
           XOPPD = DEVAL(SOPP,X,XP,S,N)
           YOPPD = DEVAL(SOPP,Y,YP,S,N)
C
           RES  = (XOPP -XLE)*DXC + (YOPP -YLE)*DYC - XBAR
           RESD =  XOPPD     *DXC +  YOPPD     *DYC
C
           IF(ABS(RES)/(S(N)-S(1)) .LT. EPS) GO TO 305
           IF(RESD .EQ. 0.0) GO TO 303
C
           DSOPP = -RES/RESD
           SOPP = SOPP + DSOPP
C
           IF(ABS(DSOPP/(S(N)-S(1))) .LT. EPS) GO TO 305
 302     CONTINUE
C
 303     WRITE(*,*)
     &      'TCSET: Opposite-point location failed. Continuing...'
         SOPP = 2.0*SLE - S(I)
C
 305     CONTINUE
        ENDIF
C
C------ set point on the opposite side with the same chord x value
        XOPP = SEVAL(SOPP,X,XP,S,N)
        YOPP = SEVAL(SOPP,Y,YP,S,N)
C
        YBAROP = (YOPP-YLE)*DXC - (XOPP-XLE)*DYC
C
        IF(LNEWSET) THEN
C------- set new chord x,y coordinates by changing camber & thickness
         YBARCT = CFAC * 0.5*(YBAR+YBAROP)
     &          + TFAC * 0.5*(YBAR-YBAROP)
         XNEW(I) = XLE + XBAR  *DXC - YBARCT*DYC
         YNEW(I) = YLE + YBARCT*DXC + XBAR  *DYC
        ENDIF
C
        TMAX = MAX( TMAX ,  ABS(YBAR-YBAROP) )
        CMAX = MAX( CMAX , 0.5*(YBAR+YBAROP) )
        CMIN = MIN( CMIN , 0.5*(YBAR+YBAROP) )
   30 CONTINUE
C
      IF(-CMIN .GT. CMAX) CMAX = CMIN
C
      RETURN
      END


      SUBROUTINE YSYM(X,XP,Y,YP,S,NX,N,ISIDE, XNEW,YNEW,NNEW)
C---------------------------------------------------------
C     Makes passed-in airfoil symmetric about chord line.
C---------------------------------------------------------
      DIMENSION X(NX),XP(NX),Y(NX),YP(NX),S(NX)
      DIMENSION XNEW(NX), YNEW(NX)
C
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHSQ = (XTE-XLE)**2 + (YTE-YLE)**2
C
      IF(CHSQ.EQ.0.0) RETURN
C
C---- set unit chord-line vector
      DXC = (XTE-XLE) / SQRT(CHSQ)
      DYC = (YTE-YLE) / SQRT(CHSQ)
C
C---- find index of node ILE which is just before leading edge point
      DO I = 1, N
        IF(S(I) .GE. SLE) GO TO 6
      ENDDO
 6    CONTINUE
      ILE = I-1
C
      DS = S(ILE+1) - S(ILE)
      IF(SLE-S(ILE-1) .LT. 0.1*DS) THEN
C------ point is just before LE, we will move it ahead to LE
        ILE1 = ILE - 1
        ILE2 = ILE + 1
      ELSE IF(S(ILE+1)-SLE .LT. 0.1*DS) THEN
C------ point is just after LE, we will move it back to LE
        ILE1 = ILE
        ILE2 = ILE + 2
      ELSE
C------ no point is near LE ... we will add new point
        ILE1 = ILE
        ILE2 = ILE + 1
      ENDIF
C
C---- set index limits of side which will set symmetric geometry
      IF(ISIDE.EQ.1) THEN
       IG1 = 1
       IG2 = ILE1
       IGDIR = +1
      ELSE
       IG1 = N
       IG2 = ILE2
       IGDIR = -1
      ENDIF
C
C---- set new number of points, including LE point
      NNEW = 2*(IABS(IG2-IG1) + 1) + 1
      IF(NNEW.GT.NX) STOP 'YSYM:  Array overflow on passed arrays.'
C
C---- set symmetric geometry
      DO I = IG1, IG2, IGDIR
C
C------ coordinates in chord-line axes
        XBAR = (X(I)-XLE)*DXC + (Y(I)-YLE)*DYC
        YBAR = (Y(I)-YLE)*DXC - (X(I)-XLE)*DYC
C
        I1 = 1    + (I - IG1)*IGDIR
        I2 = NNEW - (I - IG1)*IGDIR
C
        XNEW(I1) = XLE + XBAR*DXC - YBAR*DYC
        XNEW(I2) = XLE + XBAR*DXC + YBAR*DYC
C
        YNEW(I1) = YLE + YBAR*DXC + XBAR*DYC
        YNEW(I2) = YLE - YBAR*DXC + XBAR*DYC
      ENDDO
C
C---- set new LE point
      XNEW(NNEW/2+1) = XLE
      YNEW(NNEW/2+1) = YLE
C
      RETURN
      END ! YSYM


      SUBROUTINE LERSCL(X,XP,Y,YP,S,N, DOC,RFAC, XNEW,YNEW)
C---------------------------------------------------------
C     Adjusts airfoil to scale LE radius by factor RFAC.
C     Blending of new shape is done with decay length DOC.
C---------------------------------------------------------
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
      DIMENSION XNEW(*), YNEW(*)
C
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD = SQRT((XTE-XLE)**2 + (YTE-YLE)**2)
C
C---- set unit chord-line vector
      DXC = (XTE-XLE) / CHORD
      DYC = (YTE-YLE) / CHORD
C
      SRFAC = SQRT(ABS(RFAC))
C
C---- go over each point, changing the y-thickness appropriately
      DO 30 I = 1, N
C
C------ coordinates in chord-line axes
        XBAR = (X(I)-XLE)*DXC + (Y(I)-YLE)*DYC
        YBAR = (Y(I)-YLE)*DXC - (X(I)-XLE)*DYC
C
        SOPP = 2.0*SLE - S(I)
        SOPP = MAX(SOPP,S( 1))
        SOPP = MIN(SOPP,S(N))
C
        IF(ABS(SOPP/SLE-1.0) .LE. 1.0E-5) THEN
         SOPP = SLE
        ELSE
C------- converge on exact opposite point with same chord x value
         DO 302 ITOPP = 1, 12
           XOPP  = SEVAL(SOPP,X,XP,S,N)
           YOPP  = SEVAL(SOPP,Y,YP,S,N)
           XOPPD = DEVAL(SOPP,X,XP,S,N)
           YOPPD = DEVAL(SOPP,Y,YP,S,N)
C
           RES  = (XOPP -XLE)*DXC + (YOPP -YLE)*DYC - XBAR
           RESD =  XOPPD     *DXC +  YOPPD     *DYC
C
           IF(ABS(RES)/(S(N)-S(1)) .LT. 1.0E-5) GO TO 305
           IF(RESD .EQ. 0.0) GO TO 303
C
           DSOPP = -RES/RESD
           SOPP = SOPP + DSOPP
C
           IF(ABS(DSOPP/(S(N)-S(1))) .LT. 1.0E-5) GO TO 305
 302     CONTINUE
 303     WRITE(*,*)
     &      'LERSCL: Opposite-point location failed. Continuing...'
         SOPP = 2.0*SLE - S(I)
 305     CONTINUE
        ENDIF
C
C------ set point on the opposite side with the same chord x value
        XOPP = SEVAL(SOPP,X,XP,S,N)
        YOPP = SEVAL(SOPP,Y,YP,S,N)
C
        YBAROP = (YOPP-YLE)*DXC - (XOPP-XLE)*DYC
C
C------ thickness factor tails off exponentially towards trailing edge
        XOC = XBAR/CHORD
        ARG = MIN( XOC/DOC , 15.0 )
        TFAC = 1.0 - (1.0-SRFAC)*EXP(-ARG)
C
C------ set new chord x,y coordinates by changing thickness locally
        YBARCT = 0.5*(YBAR+YBAROP) + TFAC * 0.5*(YBAR-YBAROP)
C
        XNEW(I) = XLE + XBAR  *DXC - YBARCT*DYC
        YNEW(I) = YLE + YBARCT*DXC + XBAR  *DYC
   30 CONTINUE
C
      RETURN
      END




      SUBROUTINE FLAP(X,XS,Y,YS,S,N,
     &                XF,YF,ADEF,INSID, 
     &                XNEW,YNEW,NNEW,
     &                TOPS,ATOP,XTOP,YTOP,
     &                BOTS,ABOT,XBOT,YBOT )
C----------------------------------------------------
C     Sets modified airfoil with a deflected flap.
C     Points may be added/subtracted in the flap
C     break vicinity to clean things up.
C     Note: angle deflection ADEF is in radians
C----------------------------------------------------
      DIMENSION X(*),XS(*),Y(*),YS(*),S(*)
      DIMENSION XNEW(*), YNEW(*)
      LOGICAL INSID
C
C
      LOGICAL LT1NEW,LT2NEW,LB1NEW,LB2NEW
C
      IF(INSID) THEN
        ATOP = MAX( 0.0 , -ADEF )
        ABOT = MAX( 0.0 ,  ADEF )
      ELSE
        CHX = DEVAL(BOTS,X,XS,S,N) - DEVAL(TOPS,X,XS,S,N)
        CHY = DEVAL(BOTS,Y,YS,S,N) - DEVAL(TOPS,Y,YS,S,N)
        FVX = SEVAL(BOTS,X,XS,S,N) + SEVAL(TOPS,X,XS,S,N)
        FVY = SEVAL(BOTS,Y,YS,S,N) + SEVAL(TOPS,Y,YS,S,N)
        CRSP = CHX*(YF-0.5*FVY) - CHY*(XF-0.5*FVX)
        IF(CRSP .GT. 0.0) THEN
C-------- flap hinge is above airfoil
          ATOP = MAX( 0.0 ,  ADEF )
          ABOT = MAX( 0.0 ,  ADEF )
        ELSE
C-------- flap hinge is below airfoil
          ATOP = MAX( 0.0 , -ADEF )
          ABOT = MAX( 0.0 , -ADEF )
        ENDIF
      ENDIF
C
C---- find upper and lower surface break arc length values...
      CALL SSS(TOPS,ST1,ST2,ATOP,XF,YF,X,XS,Y,YS,S,N,1)
      CALL SSS(BOTS,SB1,SB2,ABOT,XF,YF,X,XS,Y,YS,S,N,2)
C
C---- ... and x,y coordinates
      XT1 = SEVAL(ST1,X,XS,S,N)
      YT1 = SEVAL(ST1,Y,YS,S,N)
      XT2 = SEVAL(ST2,X,XS,S,N)
      YT2 = SEVAL(ST2,Y,YS,S,N)
      XB1 = SEVAL(SB1,X,XS,S,N)
      YB1 = SEVAL(SB1,Y,YS,S,N)
      XB2 = SEVAL(SB2,X,XS,S,N)
      YB2 = SEVAL(SB2,Y,YS,S,N)
C
C
      WRITE(*,1100) XT1, YT1, XT2, YT2,
     &              XB1, YB1, XB2, YB2
 1100 FORMAT(/' Top breaks: x,y =  ', 2F9.5, 4X, 2F9.5
     &       /' Bot breaks: x,y =  ', 2F9.5, 4X, 2F9.5)
C
C---- find points adjacent to breaks
      DO I = 1, N-1
        IF(S(I).LE.ST1 .AND. S(I+1).GT.ST1) IT1 = I+1
        IF(S(I).LT.ST2 .AND. S(I+1).GE.ST2) IT2 = I
        IF(S(I).LE.SB1 .AND. S(I+1).GT.SB1) IB1 = I
        IF(S(I).LT.SB2 .AND. S(I+1).GE.SB2) IB2 = I+1
      ENDDO
C
      DSAVG = (S(N)-S(1))/FLOAT(N-1)
C
C---- smallest fraction of s increments i+1 and i+2 away from break point
      SFRAC = 0.33333
C
      IF(ATOP .NE. 0.0) THEN
        ST1P = ST1 + SFRAC*(S(IT1  )-ST1)
        ST1Q = ST1 + SFRAC*(S(IT1+1)-ST1)
        IF(S(IT1) .LT. ST1Q) THEN
C-------- simply move adjacent point to ideal SFRAC location
          XT1NEW = SEVAL(ST1Q,X,XS,S,N)
          YT1NEW = SEVAL(ST1Q,Y,YS,S,N)
          LT1NEW = .FALSE.
        ELSE
C-------- make new point at SFRAC location
          XT1NEW = SEVAL(ST1P,X,XS,S,N)
          YT1NEW = SEVAL(ST1P,Y,YS,S,N)
          LT1NEW = .TRUE.
        ENDIF
C
        ST2P = ST2 + SFRAC*(S(IT2 )-ST2)
        IT2Q = MAX(IT2-1,1)
        ST2Q = ST2 + SFRAC*(S(IT2Q)-ST2)
        IF(S(IT2) .GT. ST2Q) THEN
C-------- simply move adjacent point
          XT2NEW = SEVAL(ST2Q,X,XS,S,N)
          YT2NEW = SEVAL(ST2Q,Y,YS,S,N)
          LT2NEW = .FALSE.
        ELSE
C-------- make new point
          XT2NEW = SEVAL(ST2P,X,XS,S,N)
          YT2NEW = SEVAL(ST2P,Y,YS,S,N)
          LT2NEW = .TRUE.
        ENDIF
      ENDIF
C
      IF(ABOT .NE. 0.0) THEN
        SB1P = SB1 + SFRAC*(S(IB1  )-SB1)
        SB1Q = SB1 + SFRAC*(S(IB1-1)-SB1)
        IF(S(IB1) .GT. SB1Q) THEN
C-------- simply move adjacent point
          XB1NEW = SEVAL(SB1Q,X,XS,S,N)
          YB1NEW = SEVAL(SB1Q,Y,YS,S,N)
          LB1NEW = .FALSE.
        ELSE
C-------- make new point
          XB1NEW = SEVAL(SB1P,X,XS,S,N)
          YB1NEW = SEVAL(SB1P,Y,YS,S,N)
          LB1NEW = .TRUE.
        ENDIF
C
        SB2P = SB2 + SFRAC*(S(IB2 )-SB2)
        IB2Q = MIN(IB2+1,N)
        SB2Q = SB2 + SFRAC*(S(IB2Q)-SB2)
        IF(S(IB2) .LT. SB2Q) THEN
C-------- simply move adjacent point
          XB2NEW = SEVAL(SB2Q,X,XS,S,N)
          YB2NEW = SEVAL(SB2Q,Y,YS,S,N)
          LB2NEW = .FALSE.
        ELSE
C-------- make new point
          XB2NEW = SEVAL(SB2P,X,XS,S,N)
          YB2NEW = SEVAL(SB2P,Y,YS,S,N)
          LB2NEW = .TRUE.
        ENDIF
      ENDIF
C
cc      DSTOP = ABS(S(IT2)-S(IT1))
cc      DSBOT = ABS(S(IB2)-S(IB1))
C
      SIND = SIN(ADEF)
      COSD = COS(ADEF)
C
C-------------------------------------------------------------------
C---- initialize accumulator index for new airfoil
      INEW = 0
C
C---- upper flap surface
      DO I = 1, IT2
        XBAR = X(I) - XF
        YBAR = Y(I) - YF
C
        INEW = INEW + 1
        XNEW(INEW) = XF  +  XBAR*COSD  +  YBAR*SIND
        YNEW(INEW) = YF  -  XBAR*SIND  +  YBAR*COSD
      ENDDO
C
C
      IF(ATOP .EQ. 0.0) THEN
C------ arc length of newly created surface on top of airfoil
        DSNEW = ABS(ADEF)*SQRT((XT1-XF)**2 + (YT1-YF)**2)
C
C------ number of points to be added to define newly created surface
        NPADD = INT(1.5*DSNEW/DSAVG + 1.0)
C
        IF(NPADD.GT.0) THEN
C------- add new points along the new surface circular arc segment
         DANG = ADEF / FLOAT(NPADD)
         XBAR = XT1 - XF
         YBAR = YT1 - YF
         DO IP = 1, NPADD
           ANG = DANG*(FLOAT(IP) - 0.5)
           CA = COS(ANG)
           SA = SIN(ANG)
C
           INEW = INEW + 1
           XNEW(INEW) = XF  +  XBAR*CA + YBAR*SA
           YNEW(INEW) = YF  -  XBAR*SA + YBAR*CA
         ENDDO
        ENDIF
C
      ELSE
C------ set point in the corner and possibly two adjacent points
        IF(LT2NEW) THEN
         XBAR = XT2NEW - XF
         YBAR = YT2NEW - YF
         INEW = INEW + 1
         XNEW(INEW) = XF  +  XBAR*COSD + YBAR*SIND
         YNEW(INEW) = YF  -  XBAR*SIND + YBAR*COSD
        ENDIF
C
        INEW = INEW + 1
        XNEW(INEW) = XT1
        YNEW(INEW) = YT1
C
        IF(LT1NEW) THEN
         INEW = INEW + 1
         XNEW(INEW) = XT1NEW
         YNEW(INEW) = YT1NEW
        ENDIF
C
      ENDIF
C
C
C---- unchanged portion ahead of flap breaks
      DO I = IT1, IB1
        INEW = INEW + 1
        XNEW(INEW) = X(I)
        YNEW(INEW) = Y(I)
      ENDDO
C
C
      IF(ABOT .EQ. 0.0) THEN
C------ arc length of newly created surface on top of airfoil
        DSNEW = ABS(ADEF)*SQRT((XB1-XF)**2 + (YB1-YF)**2)
C
C------ number of points to be added to define newly created surface
        NPADD = INT(1.5*DSNEW/DSAVG + 1.0)
C
        IF(NPADD.GT.0) THEN
C------- add new points along the new surface circular arc segment
         DANG = ADEF / FLOAT(NPADD)
         XBAR = XB1 - XF
         YBAR = YB1 - YF
         DO IP=1, NPADD
           ANG = DANG*(FLOAT(IP) - 0.5)
           CA = COS(ANG)
           SA = SIN(ANG)
C
           INEW = INEW + 1
           XNEW(INEW) = XF  +  XBAR*CA + YBAR*SA
           YNEW(INEW) = YF  -  XBAR*SA + YBAR*CA
         ENDDO
        ENDIF
C
      ELSE
        IF(LB1NEW) THEN
         INEW = INEW + 1
         XNEW(INEW) = XB1NEW
         YNEW(INEW) = YB1NEW
        ENDIF
C
        INEW = INEW + 1
        XNEW(INEW) = XB1
        YNEW(INEW) = YB1
C
        IF(LB2NEW) THEN
         XBAR = XB2NEW - XF
         YBAR = YB2NEW - YF
         INEW = INEW + 1
         XNEW(INEW) = XF  +  XBAR*COSD + YBAR*SIND
         YNEW(INEW) = YF  -  XBAR*SIND + YBAR*COSD
        ENDIF
C
      ENDIF
C
C
C---- rotate flap points about the hinge point (XF,YF)
      DO I = IB2, N
C------ rotated part on flap
        XBAR = X(I) - XF
        YBAR = Y(I) - YF
C
        INEW = INEW + 1
        XNEW(INEW) = XF  +  XBAR*COSD  +  YBAR*SIND
        YNEW(INEW) = YF  -  XBAR*SIND  +  YBAR*COSD
      ENDDO
C
C
C---- total number of new points
      NNEW = INEW
C
      XTOP = XT1
      YTOP = YT1
      XBOT = XB1
      YBOT = YB1
C
      RETURN
      END ! FLAP



      LOGICAL FUNCTION INSIDE(X,Y,N, XF,YF)
      DIMENSION X(*),Y(*)
C-------------------------------------
C     Returns .TRUE. if point XF,YF 
C     is inside contour X(i),Y(i).
C-------------------------------------
C
C---- integrate subtended angle around airfoil perimeter
      ANGLE = 0.0
      DO I = 1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
        XB1 = X(I)  - XF
        YB1 = Y(I)  - YF
        XB2 = X(IP) - XF
        YB2 = Y(IP) - YF
        ANGLE = ANGLE + (XB1*YB2 - YB1*XB2)
     &                   / SQRT((XB1**2 + YB1**2)*(XB2**2 + YB2**2))
      ENDDO
C
C---- angle = 0 if XF,YF is outside, angle = +/- 2 pi  if XF,YF is inside
      INSIDE = ABS(ANGLE) .GT. 1.0
C
      RETURN
      END





      SUBROUTINE GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XF,YF)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C
      IF(XF .EQ. -999.0) THEN
       XF = 0.
       CALL ASKR('Enter flap hinge x location^',XF)
      ENDIF
C
C---- find top and bottom y at hinge x location
      TOPS = S(1) + (X(1) - XF)
      BOTS = S(N) - (X(N) - XF)
      CALL SINVRT(TOPS,XF,X,XP,S,N)      
      CALL SINVRT(BOTS,XF,X,XP,S,N)      
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
C
      WRITE(*,1000) TOPY, BOTY
 1000 FORMAT(/'  Top    surface:  y =', F8.4,'     y/t = 1.0'
     &       /'  Bottom surface:  y =', F8.4,'     y/t = 0.0')
C
      IF(YF .EQ. -999.0) THEN
       YF = 999.0
       CALL ASKR(
     &  'Enter flap hinge y location (or 999 to specify y/t)^',YF)
      ENDIF
C
      IF(YF .EQ. 999.0) THEN
        YREL = 0.5
        CALL ASKR('Enter flap hinge relative y/t location^',YREL)
        YF = TOPY*YREL + BOTY*(1.0-YREL)
      ENDIF
C
      RETURN
      END
C

c      IF(LGPARM) THEN
c        CALL NEWPEN(3)
c        CALL GPARPL(XL0,YL0,0.9*CH,.TRUE.,
c     &              CHORDB,AREAB,RADBLE,ANGBTE,
c     &              EI11BA,EI22BA,APX1BA,APX2BA,
c     &              EI11BT,EI22BT,APX1BT,APX2BT,
c     &              THICKB,CAMBRB)
c      ENDIF
cC
c      CALL PLFLUSH
cC
c      LGEOPL = .TRUE.
c      NOVER = 0
cC
c      RETURN
c      END





      SUBROUTINE SSS(SS,S1,S2,DEL,XBF,YBF,X,XP,Y,YP,S,N,ISIDE)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C----------------------------------------------------------------
C     Returns arc length points S1,S2 at flap surface break
C     locations.  S1 is on fixed airfoil part, S2 is on flap.
C     The points are defined according to two cases:
C
C
C     If DEL > 0:  Surface will be eliminated in S1 < s < S2
C
C     Returns the arc length values S1,S2 of the endpoints
C     of the airfoil surface segment which "disappears" as a
C     result of the flap deflection.  The line segments between
C     these enpoints and the flap hinge point (XBF,YBF) have
C     an included angle of DEL.  DEL is therefore the flap
C     deflection which will join up the points at S1,S2.
C     SS is an approximate arc length value near S1 and S2.
C     It is used as an initial guess for the Newton loop 
C     for S1 and S2.
C
C
C     If DEL = 0:  Surface will be created at s = S1 = S2
C
C     If DEL=0, then S1,S2 will cooincide, and will be located
C     on the airfoil surface where the segment joining the
C     point at S1,S2 and the hinge point is perpendicular to
C     the airfoil surface.  This will be the point where the
C     airfoil surface must be broken to permit a gap to open
C     as a result of the flap deflection.
C----------------------------------------------------------------
C
C---- convergence epsilon
      DATA EPS / 1.0E-5 /
C
      STOT = ABS( S(N) - S(1) )
C
      write(*,*) 'sss del ',del
      SIND = SIN(0.5*ABS(DEL))
C
      SSGN = 1.0
      IF(ISIDE.EQ.1) SSGN = -1.0
C
C---- initial guesses for S1, S2
      RSQ = (SEVAL(SS,X,XP,S,N)-XBF)**2 + (SEVAL(SS,Y,YP,S,N)-YBF)**2
      S1 = SS - (SIND*SQRT(RSQ) + EPS*STOT)*SSGN
      S2 = SS + (SIND*SQRT(RSQ) + EPS*STOT)*SSGN
C
C---- Newton iteration loop
      DO 10 ITER=1, 10
        X1  = SEVAL(S1,X,XP,S,N)
        X1P = DEVAL(S1,X,XP,S,N)
        Y1  = SEVAL(S1,Y,YP,S,N)
        Y1P = DEVAL(S1,Y,YP,S,N)
C
        X2  = SEVAL(S2,X,XP,S,N)
        X2P = DEVAL(S2,X,XP,S,N)
        Y2  = SEVAL(S2,Y,YP,S,N)
        Y2P = DEVAL(S2,Y,YP,S,N)
C
        R1SQ = (X1-XBF)**2 + (Y1-YBF)**2
        R2SQ = (X2-XBF)**2 + (Y2-YBF)**2
        R1 = SQRT(R1SQ)
        R2 = SQRT(R2SQ)
C
        RRSQ = (X1-X2)**2 + (Y1-Y2)**2
        RR = SQRT(RRSQ)
C
        IF(R1.LE.EPS*STOT .OR. R2.LE.EPS*STOT) THEN
         S1 = SS
         S2 = SS
         RETURN
        ENDIF
C
        R1_S1 = (X1P*(X1-XBF) + Y1P*(Y1-YBF))/R1
        R2_S2 = (X2P*(X2-XBF) + Y2P*(Y2-YBF))/R2
C
        IF(SIND.GT.0.01) THEN
C
         IF(RR.EQ.0.0) RETURN
C
         RR_S1 =  (X1P*(X1-X2) + Y1P*(Y1-Y2))/RR
         RR_S2 = -(X2P*(X1-X2) + Y2P*(Y1-Y2))/RR
C
C------- Residual 1: set included angle via dot product
         RS1 = ((XBF-X1)*(X2-X1) + (YBF-Y1)*(Y2-Y1))/RR - SIND*R1
         A11 = ((XBF-X1)*( -X1P) + (YBF-Y1)*( -Y1P))/RR
     &       + ((  -X1P)*(X2-X1) + (  -Y1P)*(Y2-Y1))/RR
     &       - ((XBF-X1)*(X2-X1) + (YBF-Y1)*(Y2-Y1))*RR_S1/RRSQ
     &       - SIND*R1_S1
         A12 = ((XBF-X1)*(X2P  ) + (YBF-Y1)*(Y2P  ))/RR
     &       - ((XBF-X1)*(X2-X1) + (YBF-Y1)*(Y2-Y1))*RR_S2/RRSQ
C
C------- Residual 2: set equal length segments
         RS2 = R1 - R2
         A21 = R1_S1
         A22 =    - R2_S2
        ELSE
C
C------- Residual 1: set included angle via small angle approximation
         RS1 = (R1+R2)*SIND + (S1 - S2)*SSGN
         A11 =  R1_S1 *SIND + SSGN
         A12 =  R2_S2 *SIND - SSGN
C
C------- Residual 2: set vector sum of line segments beteen the 
C-       endpoints and flap hinge to be perpendicular to airfoil surface.
         X1PP = D2VAL(S1,X,XP,S,N)
         Y1PP = D2VAL(S1,Y,YP,S,N)
         X2PP = D2VAL(S2,X,XP,S,N)
         Y2PP = D2VAL(S2,Y,YP,S,N)
C
         XTOT = X1+X2 - 2.0*XBF
         YTOT = Y1+Y2 - 2.0*YBF
C
         RS2 = XTOT*(X1P+X2P) + YTOT*(Y1P+Y2P)
         A21 =  X1P*(X1P+X2P) +  Y1P*(Y1P+Y2P) + XTOT*X1PP + YTOT*Y1PP
         A22 =  X2P*(X1P+X2P) +  Y2P*(Y1P+Y2P) + XTOT*X2PP + YTOT*Y2PP
        ENDIF
C
        DET =   A11*A22 - A12*A21
        DS1 = -(RS1*A22 - A12*RS2) / DET
        DS2 = -(A11*RS2 - RS1*A21) / DET
C
        DS1 = MIN( DS1 , 0.01*STOT )
        DS1 = MAX( DS1 , -.01*STOT )
        DS2 = MIN( DS2 , 0.01*STOT )
        DS2 = MAX( DS2 , -.01*STOT )
C
        S1 = S1 + DS1
        S2 = S2 + DS2
        IF(ABS(DS1)+ABS(DS2) .LT. EPS*STOT ) GO TO 11
   10 CONTINUE
      WRITE(*,*) 'SSS: failed to converge subtending angle points'
      S1 = SS
      S2 = SS
C
   11 CONTINUE
C
C---- make sure points are identical if included angle is zero.
      IF(DEL.EQ.0.0) THEN
       S1 = 0.5*(S1+S2)
       S2 = S1
      ENDIF
C
      RETURN
      END


      SUBROUTINE CLIS(X,XP,Y,YP,S,N)
      DIMENSION X(*), XP(*), Y(*), YP(*), S(*)
C-------------------------------------------------------------------
C     Displays curvatures at panel nodes.
C-------------------------------------------------------------------
C
      CMAX = 0.0
      IMAX = 1
C
C---- go over each point, calculating curvature
      WRITE(*,1050)
      DO 30 I=1, N
        CV = CURV(S(I),X,XP,Y,YP,S,N)
        WRITE(*,1100) I, X(I), Y(I), CV
        IF(ABS(CV) .GT. ABS(CMAX)) THEN
         CMAX = CV
         IMAX = I
        ENDIF
   30 CONTINUE
C
      WRITE(*,1200) CMAX, IMAX, X(IMAX), Y(IMAX)
C
      RETURN
C
 1050 FORMAT(/'  i       x        y           curv')
CCC             120   0.2134  -0.0234      2025.322
 1100 FORMAT(1X,I3, 2F9.4, F14.3)
 1200 FORMAT(/' Maximum curvature =', F14.3,
     &        '   at  i,x,y  = ', I3, 2F9.4 )
      END ! CLIS


      SUBROUTINE CANG(X,Y,N,LPRINT)
      DIMENSION X(*), Y(*)
      LOGICAL LPRINT
C-------------------------------------------------------------------
C     LPRINT=t:   Displays all panel node corner angles
C     LPRINT=f:   Displays max panel node corner angle
C-------------------------------------------------------------------
      DATA  PI /3.1415926535897932384/
C
      AMAX = 0.0
      IMAX = 1
C
C---- go over each point, calculating corner angle
      IF(LPRINT) WRITE(*,1050)
      DO 30 I=2, N-1
        DX1 = X(I) - X(I-1)
        DY1 = Y(I) - Y(I-1)
        DX2 = X(I) - X(I+1)
        DY2 = Y(I) - Y(I+1)
C
C------ allow for doubled points
        IF(DX1.EQ.0.0 .AND. DY1.EQ.0.0) THEN
         DX1 = X(I) - X(I-2)
         DY1 = Y(I) - Y(I-2)
        ENDIF
        IF(DX2.EQ.0.0 .AND. DY2.EQ.0.0) THEN
         DX2 = X(I) - X(I+2)
         DY2 = Y(I) - Y(I+2)
        ENDIF
C
        CROSSP = (DX2*DY1 - DY2*DX1)
     &         / SQRT((DX1**2 + DY1**2) * (DX2**2 + DY2**2))
        ANGL = ASIN(CROSSP)*(180.0/PI)
        IF(LPRINT) WRITE(*,1100) I, X(I), Y(I), ANGL
        IF(ABS(ANGL) .GT. ABS(AMAX)) THEN
         AMAX = ANGL
         IMAX = I
        ENDIF
   30 CONTINUE
C
      WRITE(*,1200) AMAX, IMAX, X(IMAX), Y(IMAX)
C
      RETURN
C
 1050 FORMAT(/'  i       x        y      angle')
CCC             120   0.2134  -0.0234   25.322
 1100 FORMAT(1X,I3, 2F9.4, F9.3)
 1200 FORMAT(/' Maximum panel corner angle =', F7.3,
     &        '   at  i,x,y  = ', I3, 2F9.4 )
      END ! CANG


c      SUBROUTINE ADDXYA(IEL,DX,DY,DF,DA)
c      INCLUDE 'AIRPAN.INC'
cC
c      EFE = EXP(DF)
c      SAE = SIN(-DA)
c      CAE = COS(-DA)
c      DXE = DXEL(IEL)
c      DYE = DYEL(IEL)
cC
c      DXEL(IEL) = DXE*CAE - DYE*SAE + DX
c      DYEL(IEL) = DXE*SAE + DYE*CAE + DY
c      DFEL(IEL) = DFEL(IEL) + DF
c      DAEL(IEL) = DAEL(IEL) - DA
cC
c      RETURN
c      END

