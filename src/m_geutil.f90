module m_geutil
    implicit none
contains
    !*==MINMAX.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! ITPSET
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
    !=========================================================================
    
    SUBROUTINE MINMAX(N, X, XMIN, XMAX)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: XMAX, XMIN
        REAL, DIMENSION(*) :: X
        !
        ! Local variables
        !
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !--------------------------------------
        !     Calculates min,max of array X
        !--------------------------------------
        !
        I = 1
        XMIN = X(I)
        XMAX = X(I)
        DO I = 2, N
            XMIN = MIN(XMIN, X(I))
            XMAX = MAX(XMAX, X(I))
        ENDDO
        !
    END SUBROUTINE MINMAX
    !*==NCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE NCALC(X, Y, S, N, XN, YN)
        use m_spline, only: segspl
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(*) :: S, X, XN, Y, YN
        !
        ! Local variables
        !
        INTEGER :: I
        REAL :: SMOD, SX, SY
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------
        !     Calculates normal unit vector
        !     components at airfoil panel nodes
        !---------------------------------------
        !
        IF (N<=1) RETURN
        !
        CALL SEGSPL(X, XN, S, N)
        CALL SEGSPL(Y, YN, S, N)
        DO I = 1, N
            SX = YN(I)
            SY = -XN(I)
            SMOD = SQRT(SX * SX + SY * SY)
            XN(I) = SX / SMOD
            YN(I) = SY / SMOD
        ENDDO
        !
        !---- average normal vectors at corner points
        DO I = 1, N - 1
            IF (S(I)==S(I + 1)) THEN
                SX = 0.5 * (XN(I) + XN(I + 1))
                SY = 0.5 * (YN(I) + YN(I + 1))
                SMOD = SQRT(SX * SX + SY * SY)
                XN(I) = SX / SMOD
                YN(I) = SY / SMOD
                XN(I + 1) = SX / SMOD
                YN(I + 1) = SY / SMOD
            ENDIF
        ENDDO
        !
    END SUBROUTINE NCALC
    !*==LEFIND.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE LEFIND(SLE, X, XP, Y, YP, S, N)
        use m_spline, only: d2val, seval, deval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: SLE
        REAL, DIMENSION(*) :: S, X, XP, Y, YP
        !
        ! Local variables
        !
        REAL :: DOTP, DSEPS, DSLE, DX, DXDD, DXDS, DXTE, DY, &
                & DYDD, DYDS, DYTE, RES, RESS, XCHORD, XLE, XTE, &
                & YCHORD, YLE, YTE
        INTEGER :: I, ITER
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Locates leading edge spline-parameter value SLE
        !
        !     The defining condition is
        !
        !      (X-XTE,Y-YTE) . (X',Y') = 0     at  S = SLE
        !
        !     i.e. the surface tangent is normal to the chord
        !     line connecting X(SLE),Y(SLE) and the TE point.
        !------------------------------------------------------
        IF (N<=1) THEN
            SLE = S(1)
            RETURN
        ENDIF
        !
        !---- convergence tolerance
        DSEPS = (S(N) - S(1)) * 1.0E-5
        !
        !---- set trailing edge point coordinates
        XTE = 0.5 * (X(1) + X(N))
        YTE = 0.5 * (Y(1) + Y(N))
        !
        !---- get first guess for SLE
        DO I = 3, N - 2
            DXTE = X(I) - XTE
            DYTE = Y(I) - YTE
            DX = X(I + 1) - X(I)
            DY = Y(I + 1) - Y(I)
            DOTP = DXTE * DX + DYTE * DY
            IF (DOTP<0.0) EXIT
        ENDDO
        !
        SLE = S(I)
        !
        !---- check for sharp LE case
        IF (S(I)==S(I - 1)) RETURN
        !
        !---- Newton iteration to get exact SLE value
        DO ITER = 1, 50
            XLE = SEVAL(SLE, X, XP, S, N)
            YLE = SEVAL(SLE, Y, YP, S, N)
            DXDS = DEVAL(SLE, X, XP, S, N)
            DYDS = DEVAL(SLE, Y, YP, S, N)
            DXDD = D2VAL(SLE, X, XP, S, N)
            DYDD = D2VAL(SLE, Y, YP, S, N)
            !
            XCHORD = XLE - XTE
            YCHORD = YLE - YTE
            !
            !------ drive dot product between chord line and LE tangent to zero
            RES = XCHORD * DXDS + YCHORD * DYDS
            RESS = DXDS * DXDS + DYDS * DYDS + XCHORD * DXDD + YCHORD * DYDD
            !
            IF (RESS==0.0) EXIT
            !
            !------ Newton delta for SLE
            DSLE = -RES / RESS
            !
            DSLE = MAX(DSLE, -0.02 * ABS(XCHORD + YCHORD))
            DSLE = MIN(DSLE, 0.02 * ABS(XCHORD + YCHORD))
            SLE = SLE + DSLE
            IF (ABS(DSLE)<DSEPS) RETURN
        ENDDO
        WRITE (*, *) 'LEFIND:  LE point not found.  Continuing...'
        !
        SLE = S(I)
    END SUBROUTINE LEFIND
    !*==CNORM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE CNORM(X, XP, Y, YP, S, N, XTRAN, YTRAN, SCALE, ANGLE)
        use m_spline, only: seval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: ANGLE, SCALE, XTRAN, YTRAN
        INTEGER :: N
        REAL, DIMENSION(*) :: S, X, XP, Y, YP
        !
        ! Local variables
        !
        REAL :: CINV, CSQINV, CSSQ, CX, CY, SLE, XBAR, XLE, &
                & XTE, YBAR, YLE, YTE
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------
        !     Scales and rotates coordinates to get unit chordline
        !----------------------------------------------------------
        CALL LEFIND(SLE, X, XP, Y, YP, S, N)
        !
        XTE = 0.5 * (X(1) + X(N))
        YTE = 0.5 * (Y(1) + Y(N))
        XLE = SEVAL(SLE, X, XP, S, N)
        YLE = SEVAL(SLE, Y, YP, S, N)
        !
        CX = XTE - XLE
        CY = YTE - YLE
        CSSQ = CX * CX + CY * CY
        !
        IF (CSSQ==0.0) THEN
            CSQINV = 1.0
            CINV = 1.0
        ELSE
            CSQINV = 1.0 / CSSQ
            CINV = 1.0 / SQRT(CSSQ)
        ENDIF
        !
        DO I = 1, N
            XBAR = X(I) - XLE
            YBAR = Y(I) - YLE
            X(I) = (XBAR * CX + YBAR * CY) * CSQINV
            Y(I) = (-XBAR * CY + YBAR * CX) * CSQINV
            S(I) = S(I) * CINV
        ENDDO
        !
        XTRAN = -XLE
        YTRAN = -YLE
        SCALE = CINV
        ANGLE = ATAN2(CY, CX)
        !
    END SUBROUTINE CNORM
    !*==GEOPAR.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CNORM
    
    
    SUBROUTINE GEOPAR(X, XP, Y, YP, S, N, T, SLE, CHORD, RADLE, ANGTE, AREA, &
            & EIXXA, EIYYA, EIXYA, ASKN, EIXXT, EIYYT, EIXYT, THICK, &
            & CAMBR, VOLM, VSKN, ASRF, RGXV, RGYV)
        use m_spline, only: seval, curv
        use m_xutils, only: atanc
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: ANGTE, AREA, ASKN, ASRF, CAMBR, CHORD, EIXXA, &
                & EIXXT, EIXYA, EIXYT, EIYYA, EIYYT, RADLE, RGXV, &
                & RGYV, SLE, THICK, VOLM, VSKN
        INTEGER :: N
        REAL, DIMENSION(*) :: S, T, X, XP, Y, YP
        !
        ! Local variables
        !
        REAL :: ANG1, ANG2, CFAC, CHSQ, CURVLE, PERIM, RGXT, &
                & RGYT, TFAC, XCENA, XCENT, XCENV, XCENVT, XLE, &
                & XNEW, XTE, YCENA, YCENT, YCENV, YCENVT, YLE, &
                & YNEW, YTE
        REAL, DIMENSION(1) :: XNEW_TEMP, YNEW_TEMP
        !
        !*** End of declarations rewritten by SPAG
        !
        !------------------------------------------------------
        !     Sets geometric parameters for airfoil shape
        !------------------------------------------------------
        !
        CALL LEFIND(SLE, X, XP, Y, YP, S, N)
        !
        XLE = SEVAL(SLE, X, XP, S, N)
        YLE = SEVAL(SLE, Y, YP, S, N)
        XTE = 0.5 * (X(1) + X(N))
        YTE = 0.5 * (Y(1) + Y(N))
        !
        CHSQ = (XTE - XLE)**2 + (YTE - YLE)**2
        CHORD = SQRT(CHSQ)
        !
        CURVLE = CURV(SLE, X, XP, Y, YP, S, N)
        !
        IF (ABS(CURVLE)>0.001 * (S(N) - S(1))) THEN
            RADLE = 1.0 / CURVLE
        ELSE
            RADLE = 0.
        ENDIF
        !
        ANG1 = ATAN2(-YP(1), -XP(1))
        ANG2 = ATANC(YP(N), XP(N), ANG1)
        ANGTE = ANG2 - ANG1
        !
        !---- 2D properties
        CALL AECALC(N, X, Y, T, 0, PERIM, AREA, XCENA, YCENA, EIXXA, EIYYA, EIXYA)
        CALL AECALC(N, X, Y, T, 2, PERIM, ASKN, XCENT, YCENT, EIXXT, EIYYT, EIXYT)
        TFAC = 1.0
        CFAC = 1.0
        CALL TCSET(X, XP, Y, YP, S, N, THICK, CAMBR, TFAC, CFAC, XNEW_TEMP, YNEW_TEMP, .FALSE.)
        XNEW = XNEW_TEMP(1)
        YNEW = YNEW_TEMP(1)
        !
        !---- Axisymmetric properties
        CALL AXCALC(N, X, Y, T, 0, VOLM, ASRF, XCENV, YCENV, RGXV, RGYV)
        CALL AXCALC(N, X, Y, T, 2, VSKN, ASRF, XCENVT, YCENVT, RGXT, RGYT)
        !
    END SUBROUTINE GEOPAR
    !*==AXCALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE AXCALC(N, X, Y, T, ITYPE, VOLM, AREA, XCEN, YCEN, RGX, RGY)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: AREA, RGX, RGY, VOLM, XCEN, YCEN
        INTEGER :: ITYPE, N
        REAL, DIMENSION(*) :: T, X, Y
        !
        ! Local variables
        !
        REAL :: AINT, DA, DS, DSTE, DV, DX, DY, RGXSQ, SINT, &
                & TA, VINT, XA, XINT, XXINT, YA, YINT, YYINT
        INTEGER :: IO, IP
        LOGICAL :: LAXIBOD, LOPEN
        REAL, SAVE :: PI
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------
        !     Calculates geometric properties of axisymmetric shape X,Y
        !     where X is axial coordinate and Y is radial coordinate
        !
        !     Input:
        !       N      number of points
        !       X(.)   shape coordinate point arrays
        !       Y(.)
        !       T(.)   skin-thickness array, used only if ITYPE = 2
        !       ITYPE  = 0 ...   integration is over volume
        !              = 2 ...   integration is over skin with thickness specified
        !
        !     Output:
        !       VOLM  volume
        !       AREA  surface area
        !       XCEN  centroid location
        !       YCEN
        !       RGX   radius of gyration about X centroid
        !       RGY   radius of gyration about Y centroid (Y=0)
        !---------------------------------------------------------------
        DATA PI/3.141592653589793238/
        !
        !--- check for shape symmetric about X axis
        !    if first or last point is on centerline it is an axisymmetric body
        LAXIBOD = (Y(1)==0.0 .OR. Y(N)==0.0)
        !
        !--- check shape for "open" or closed body by comparing distance between
        !    first and last points to 0.1 of perimeter length
        !
        SINT = 0.0
        DO IO = 2, N
            DX = X(IO) - X(IO - 1)
            DY = Y(IO) - Y(IO - 1)
            DS = SQRT(DX * DX + DY * DY)
            SINT = SINT + DS
        ENDDO
        DX = X(N) - X(1)
        DY = Y(N) - Y(1)
        DSTE = SQRT(DX * DX + DY * DY)
        LOPEN = (DSTE>0.1 * SINT)
        !
        !--- integrate around contour to get geometric quantities
        VINT = 0.0
        AINT = 0.0
        XINT = 0.0
        YINT = 0.0
        XXINT = 0.0
        YYINT = 0.0
        !
        DO IO = 1, N
            IF (IO<N) THEN
                IP = IO + 1
            ELSE
                !--- if last point is on centerline or body is open end here...
                IF (LAXIBOD .OR. LOPEN) CYCLE
                !--- otherwise close body with last segment between first and last points
                IP = 1
            ENDIF
            !
            !        write(54,199) io,x(io),y(io),x(ip),y(ip)
            DX = X(IO) - X(IP)
            DY = Y(IO) - Y(IP)
            XA = (X(IO) + X(IP)) * 0.5
            YA = (Y(IO) + Y(IP)) * 0.5
            DS = SQRT(DX * DX + DY * DY)
            199     FORMAT (I4, 7F10.6)
            !---- surface area
            DA = 2.0 * PI * YA * DS
            AINT = AINT + DA
    
            IF (ITYPE==0) THEN
                !-------- treat as a axisymmetric solid body (YCEN=0,intXY=0)
                DV = PI * DX * YA**2
                !        write(55,199) io,dx,dy,xa,ya,ds,da,dv
                VINT = VINT + DV
                XINT = XINT + XA * DV
                YINT = 0.0
                XXINT = XXINT + XA * XA * DV
                YYINT = YYINT + YA * YA * DV / 4.0
            ELSE
                !--------- integrate over skin thickness
                TA = (T(IO) + T(IP)) * 0.50
                DV = 2.0 * PI * YA * TA * DS
                VINT = VINT + DV
                XINT = XINT + XA * DV
                YINT = 0.0
                XXINT = XXINT + XA * XA * DV
                YYINT = YYINT + YA * YA * DV
            ENDIF
        ENDDO
        !
        !---- volume and surface area
        IF (LOPEN .AND. .NOT.LAXIBOD) VINT = 0.0
        VOLM = VINT
        AREA = AINT
        !
        !---- calculate centroid location
        IF (VINT<=0.0) THEN
            XCEN = 0.
            YCEN = 0.
            RGX = 0.
            RGY = 0.
        ELSE
            XCEN = XINT / VINT
            YCEN = YINT / VINT
            !---- calculate radii of gyration
            RGXSQ = XXINT / VINT - XCEN**2
            IF (RGXSQ>=0.0) RGX = SQRT(RGXSQ)
            RGY = SQRT(YYINT / VINT)
        ENDIF
        !
    END SUBROUTINE AXCALC
    !*==AECALC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! AXCALC
    
    
    
    SUBROUTINE AECALC(N, X, Y, T, ITYPE, PERIM, AREA, XCEN, YCEN, EIXX, EIYY, &
            & EIXY)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: AREA, EIXX, EIXY, EIYY, PERIM, XCEN, YCEN
        INTEGER :: ITYPE, N
        REAL, DIMENSION(*) :: T, X, Y
        !
        ! Local variables
        !
        REAL :: AINT, DA, DS, DX, DY, SINT, TA, VINT, XA, &
                & XINT, XXINT, XYINT, YA, YINT, YYINT
        INTEGER :: IO, IP
        LOGICAL :: LAXIS
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------------
        !     Calculates geometric properties of 2D shape X,Y
        !
        !     Input:
        !       N      number of points
        !       X(.)   shape coordinate point arrays
        !       Y(.)
        !       T(.)   skin-thickness array, used only if ITYPE = 2
        !       ITYPE  = 0 ...   integration is over whole area  dA = dx dy
        !              = 1 ...   integration is over perimeter   dA =   ds
        !              = 2 ...   integration is over skin area   dA = t ds
        !
        !     Output:
        !       PERIM perimeter length
        !       AREA  area (assumes 2D)
        !       XCEN  centroid location
        !       YCEN
        !       EIXX  moments of inertia
        !       EIYY
        !       EIXY
        !---------------------------------------------------------------
        !
        !--- check for shape with symmetry about X axis
        !    if first or last point is on centerline it is an axisymmetric body
        LAXIS = (Y(1)==0.0 .OR. Y(N)==0.0)
        !
        SINT = 0.0
        VINT = 0.0
        AINT = 0.0
        XINT = 0.0
        YINT = 0.0
        XXINT = 0.0
        XYINT = 0.0
        YYINT = 0.0
        !
        DO IO = 1, N
            !
            IF (IO<N) THEN
                IP = IO + 1
            ELSE
                !--- if last point is on axis (symmetrical body)
                IF (LAXIS) CYCLE
                !--- otherwise close body with last segment between first and last points
                IP = 1
            ENDIF
            !
            DX = X(IO) - X(IP)
            DY = Y(IO) - Y(IP)
            XA = (X(IO) + X(IP)) * 0.5
            YA = (Y(IO) + Y(IP)) * 0.5
            !
            DS = SQRT(DX * DX + DY * DY)
            SINT = SINT + DS
    
            IF (ITYPE==0) THEN
                !-------- integrate over airfoil cross-section
                DA = YA * DX
                AINT = AINT + DA
                XINT = XINT + XA * DA
                YINT = YINT + YA * DA / 2.0
                XXINT = XXINT + XA * XA * DA
                XYINT = XYINT + XA * YA * DA / 2.0
                YYINT = YYINT + YA * YA * DA / 3.0
            ELSE
                IF (ITYPE==1) THEN
                    !--------- integrate over perimeter
                    DA = DS
                ELSE
                    !--------- integrate over skin thickness
                    TA = (T(IO) + T(IP)) * 0.50
                    DA = TA * DS
                ENDIF
                AINT = AINT + DA
                XINT = XINT + XA * DA
                YINT = YINT + YA * DA
                XXINT = XXINT + XA * XA * DA
                XYINT = XYINT + XA * YA * DA
                YYINT = YYINT + YA * YA * DA
            ENDIF
        ENDDO
        !
        PERIM = SINT
        AREA = AINT
        !
        !---- calculate centroid location
        IF (AINT==0.0) THEN
            XCEN = 0.
            YCEN = 0.
        ELSE
            XCEN = XINT / AINT
            YCEN = YINT / AINT
        ENDIF
        !
        !---- calculate inertias
        EIXX = YYINT - YCEN * YCEN * AINT
        EIYY = XXINT - XCEN * XCEN * AINT
        EIXY = XYINT - XCEN * YCEN * AINT
        !
    END SUBROUTINE AECALC
    !*==PRAXIS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! AECALC
    
    
    
    SUBROUTINE PRAXIS(FXX, FYY, FXY, F11, F22, AP1, AP2)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: AP1, AP2, F11, F22, FXX, FXY, FYY
        !
        ! Local variables
        !
        REAL :: C1, C2, FSQ, S1, S2, SGN
        REAL, SAVE :: PI
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------------------
        !     Sets prinicipal components F11,F22 of tensor FXX,FYY,FXY.
        !     Also sets corresponding principal direction angls AP1,AP2.
        !-----------------------------------------------------------------
        !
        DATA PI/3.141592653589793238/
        !
        FSQ = 0.25 * (FXX - FYY)**2 + FXY**2
        SGN = SIGN(1.0, FYY - FXX)
        F11 = 0.5 * (FXX + FYY) - SGN * SQRT(FSQ)
        F22 = 0.5 * (FXX + FYY) + SGN * SQRT(FSQ)
        !
        IF (F11==0.0 .OR. F22==0.0) THEN
            !----- vanishing tensor limit
            AP1 = 0.0
            AP2 = ATAN2(1.0, 0.0)
            !
        ELSEIF (FSQ / (F11 * F22)<0.0001) THEN
            !----- rotationally-invariant tensor (circle, square, etc.)
            AP1 = 0.0
            AP2 = ATAN2(1.0, 0.0)
            !
        ELSE
            !----- normal general tensor case
            C1 = FXY
            S1 = FXX - F11
            !
            C2 = FXY
            S2 = FXX - F22
            !
            IF (ABS(S1)>ABS(S2)) THEN
                AP1 = ATAN2(S1, C1)
                AP2 = AP1 + 0.5 * PI
            ELSE
                AP2 = ATAN2(S2, C2)
                AP1 = AP2 - 0.5 * PI
            ENDIF
    
            IF (AP1<-0.5 * PI) AP1 = AP1 + PI
            IF (AP1>+0.5 * PI) AP1 = AP1 - PI
            IF (AP2<-0.5 * PI) AP2 = AP2 + PI
            IF (AP2>+0.5 * PI) AP2 = AP2 - PI
            !
        ENDIF
    
    END SUBROUTINE PRAXIS
    !*==TEGAP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE TEGAP(X, XS, Y, YS, S, N, DOC, GAPNEW)
        use m_spline, only: scalc, seval, segspl
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: DOC, GAPNEW
        INTEGER :: N
        REAL, DIMENSION(*) :: S, X, XS, Y, YS
        !
        ! Local variables
        !
        REAL :: ARG, CHSQ, DGAP, DXN, DXU, DYN, DYU, GAP, SLE, &
                & TFAC, XLE, XOC, XTE, YLE, YTE
        INTEGER :: I
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------
        !     Used to set buffer airfoil
        !     trailing edge gap
        !----------------------------------
        !
        CALL LEFIND(SLE, X, XS, Y, YS, S, N)
        XLE = SEVAL(SLE, X, XS, S, N)
        YLE = SEVAL(SLE, Y, YS, S, N)
        XTE = 0.5 * (X(1) + X(N))
        YTE = 0.5 * (Y(1) + Y(N))
        CHSQ = (XTE - XLE)**2 + (YTE - YLE)**2
        !
        DXN = X(1) - X(N)
        DYN = Y(1) - Y(N)
        GAP = SQRT(DXN**2 + DYN**2)
        !
        !---- components of unit vector parallel to TE gap
        IF (GAP>0.0) THEN
            DXU = DXN / GAP
            DYU = DYN / GAP
        ELSE
            DXU = -.5 * (YS(N) - YS(1))
            DYU = 0.5 * (XS(N) - XS(1))
        ENDIF
        !
        DGAP = GAPNEW - GAP
        !
        !---- go over each point, changing the y-thickness appropriately
        DO I = 1, N
            !------ chord-based x/c
            XOC = ((X(I) - XLE) * (XTE - XLE) + (Y(I) - YLE) * (YTE - YLE)) / CHSQ
            !
            !------ thickness factor tails off exponentially away from trailing edge
            IF (DOC==0.0) THEN
                TFAC = 0.0
                IF (I==1 .OR. I==N) TFAC = 1.0
            ELSE
                ARG = MIN((1.0 - XOC) * (1.0 / DOC - 1.0), 15.0)
                TFAC = EXP(-ARG)
            ENDIF
            !
            IF (S(I)<=SLE) THEN
                X(I) = X(I) + 0.5 * DGAP * XOC * TFAC * DXU
                Y(I) = Y(I) + 0.5 * DGAP * XOC * TFAC * DYU
            ELSE
                X(I) = X(I) - 0.5 * DGAP * XOC * TFAC * DXU
                Y(I) = Y(I) - 0.5 * DGAP * XOC * TFAC * DYU
            ENDIF
        ENDDO
        !
        CALL SCALC(X, Y, S, N)
        CALL SEGSPL(X, XS, S, N)
        CALL SEGSPL(Y, YS, S, N)
        !
    END SUBROUTINE TEGAP
    !*==TCSET.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! TEGAP
    
    
    
    SUBROUTINE TCSET(X, XP, Y, YP, S, N, TMAX, CMAX, TFAC, CFAC, XNEW, YNEW, &
            & LNEWSET)
        use m_spline, only: deval, seval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: CFAC, CMAX, TFAC, TMAX
        LOGICAL :: LNEWSET
        INTEGER :: N
        REAL, DIMENSION(*) :: S, X, XNEW, XP, Y, YNEW, YP
        !
        ! Local variables
        !
        REAL :: CHBSQ, CMIN, DSOPP, DXC, DYC, RES, RESD, SLE, &
                & SOPP, STOT, XBAR, XLE, XOPP, XOPPD, XTE, YBAR, &
                & YBARCT, YBAROP, YLE, YOPP, YOPPD, YTE
        REAL, SAVE :: EPS
        INTEGER :: I, ITOPP
        !
        !*** End of declarations rewritten by SPAG
        !
        !-----------------------------------------------------
        !     Determines max thickness and camber.
        !     Scales thickness and camber by TFAC,CFAC
        !-----------------------------------------------------
        !
        DATA EPS/1.0E-5/
        !
        CALL LEFIND(SLE, X, XP, Y, YP, S, N)
        XLE = SEVAL(SLE, X, XP, S, N)
        YLE = SEVAL(SLE, Y, YP, S, N)
        XTE = 0.5 * (X(1) + X(N))
        YTE = 0.5 * (Y(1) + Y(N))
        CHBSQ = (XTE - XLE)**2 + (YTE - YLE)**2
        !
        IF (CHBSQ==0.0) THEN
            TMAX = 0.
            CMAX = 0.
            RETURN
        ENDIF
        !
        !---- set unit chord-line vector
        DXC = (XTE - XLE) / SQRT(CHBSQ)
        DYC = (YTE - YLE) / SQRT(CHBSQ)
        !
        STOT = S(N) - S(1)
        !
        TMAX = 0.
        CMAX = 0.
        CMIN = 0.
        !
        !---- go over each point, finding and/or changing the y-thickness
        !-    (defined normal to chord line)
        DO I = 1, N
            !
            !------ coordinates in chord-line axes
            XBAR = (X(I) - XLE) * DXC + (Y(I) - YLE) * DYC
            YBAR = (Y(I) - YLE) * DXC - (X(I) - XLE) * DYC
            !
            SOPP = 2.0 * SLE - S(I)
            SOPP = MAX(SOPP, S(1))
            SOPP = MIN(SOPP, S(N))
            !
            IF (ABS(SOPP - SLE)<=EPS * STOT) THEN
                SOPP = SLE
            ELSE
                !------- converge on exact opposite point with same chord x value
                DO ITOPP = 1, 12
                    XOPP = SEVAL(SOPP, X, XP, S, N)
                    YOPP = SEVAL(SOPP, Y, YP, S, N)
                    XOPPD = DEVAL(SOPP, X, XP, S, N)
                    YOPPD = DEVAL(SOPP, Y, YP, S, N)
                    !
                    RES = (XOPP - XLE) * DXC + (YOPP - YLE) * DYC - XBAR
                    RESD = XOPPD * DXC + YOPPD * DYC
                    !
                    IF (ABS(RES) / (S(N) - S(1))<EPS) GOTO 305
                    IF (RESD==0.0) EXIT
                    !
                    DSOPP = -RES / RESD
                    SOPP = SOPP + DSOPP
                    !
                    IF (ABS(DSOPP / (S(N) - S(1)))<EPS) GOTO 305
                ENDDO
                !
                WRITE (*, *)                                                 &
                        &'TCSET: Opposite-point location failed. Continuing...'
                SOPP = 2.0 * SLE - S(I)
                !
            ENDIF
            !
            !------ set point on the opposite side with the same chord x value
            305     XOPP = SEVAL(SOPP, X, XP, S, N)
            YOPP = SEVAL(SOPP, Y, YP, S, N)
            !
            YBAROP = (YOPP - YLE) * DXC - (XOPP - XLE) * DYC
            !
            IF (LNEWSET) THEN
                !------- set new chord x,y coordinates by changing camber & thickness
                YBARCT = CFAC * 0.5 * (YBAR + YBAROP) + TFAC * 0.5 * (YBAR - YBAROP)
                XNEW(I) = XLE + XBAR * DXC - YBARCT * DYC
                YNEW(I) = YLE + YBARCT * DXC + XBAR * DYC
            ENDIF
            !
            TMAX = MAX(TMAX, ABS(YBAR - YBAROP))
            CMAX = MAX(CMAX, 0.5 * (YBAR + YBAROP))
            CMIN = MIN(CMIN, 0.5 * (YBAR + YBAROP))
        ENDDO
        !
        IF (-CMIN>CMAX) CMAX = CMIN
        !
    END SUBROUTINE TCSET
    !*==YSYM.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE YSYM(X, XP, Y, YP, S, NX, N, ISIDE, XNEW, YNEW, NNEW)
        use m_spline, only: seval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: ISIDE, N, NNEW, NX
        REAL, DIMENSION(NX) :: S, X, XNEW, XP, Y, YNEW, YP
        !
        ! Local variables
        !
        REAL :: CHSQ, DS, DXC, DYC, SLE, XBAR, XLE, XTE, YBAR, &
                & YLE, YTE
        INTEGER :: I, I1, I2, IG1, IG2, IGDIR, ILE, ILE1, ILE2
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Makes passed-in airfoil symmetric about chord line.
        !---------------------------------------------------------
        !
        CALL LEFIND(SLE, X, XP, Y, YP, S, N)
        XLE = SEVAL(SLE, X, XP, S, N)
        YLE = SEVAL(SLE, Y, YP, S, N)
        XTE = 0.5 * (X(1) + X(N))
        YTE = 0.5 * (Y(1) + Y(N))
        CHSQ = (XTE - XLE)**2 + (YTE - YLE)**2
        !
        IF (CHSQ==0.0) RETURN
        !
        !---- set unit chord-line vector
        DXC = (XTE - XLE) / SQRT(CHSQ)
        DYC = (YTE - YLE) / SQRT(CHSQ)
        !
        !---- find index of node ILE which is just before leading edge point
        DO I = 1, N
            IF (S(I)>=SLE) EXIT
        ENDDO
        ILE = I - 1
        !
        DS = S(ILE + 1) - S(ILE)
        IF (SLE - S(ILE - 1)<0.1 * DS) THEN
            !------ point is just before LE, we will move it ahead to LE
            ILE1 = ILE - 1
            ILE2 = ILE + 1
        ELSEIF (S(ILE + 1) - SLE<0.1 * DS) THEN
            !------ point is just after LE, we will move it back to LE
            ILE1 = ILE
            ILE2 = ILE + 2
        ELSE
            !------ no point is near LE ... we will add new point
            ILE1 = ILE
            ILE2 = ILE + 1
        ENDIF
        !
        !---- set index limits of side which will set symmetric geometry
        IF (ISIDE==1) THEN
            IG1 = 1
            IG2 = ILE1
            IGDIR = +1
        ELSE
            IG1 = N
            IG2 = ILE2
            IGDIR = -1
        ENDIF
        !
        !---- set new number of points, including LE point
        NNEW = 2 * (IABS(IG2 - IG1) + 1) + 1
        IF (NNEW>NX) STOP 'YSYM:  Array overflow on passed arrays.'
        !
        !---- set symmetric geometry
        DO I = IG1, IG2, IGDIR
            !
            !------ coordinates in chord-line axes
            XBAR = (X(I) - XLE) * DXC + (Y(I) - YLE) * DYC
            YBAR = (Y(I) - YLE) * DXC - (X(I) - XLE) * DYC
            !
            I1 = 1 + (I - IG1) * IGDIR
            I2 = NNEW - (I - IG1) * IGDIR
            !
            XNEW(I1) = XLE + XBAR * DXC - YBAR * DYC
            XNEW(I2) = XLE + XBAR * DXC + YBAR * DYC
            !
            YNEW(I1) = YLE + YBAR * DXC + XBAR * DYC
            YNEW(I2) = YLE - YBAR * DXC + XBAR * DYC
        ENDDO
        !
        !---- set new LE point
        XNEW(NNEW / 2 + 1) = XLE
        YNEW(NNEW / 2 + 1) = YLE
        !
    END SUBROUTINE YSYM
    !*==LERSCL.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! YSYM
    
    
    SUBROUTINE LERSCL(X, XP, Y, YP, S, N, DOC, RFAC, XNEW, YNEW)
        use m_spline, only: deval, seval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: DOC, RFAC
        INTEGER :: N
        REAL, DIMENSION(*) :: S, X, XNEW, XP, Y, YNEW, YP
        !
        ! Local variables
        !
        REAL :: ARG, CHORD, DSOPP, DXC, DYC, RES, RESD, SLE, &
                & SOPP, SRFAC, TFAC, XBAR, XLE, XOC, XOPP, XOPPD, &
                & XTE, YBAR, YBARCT, YBAROP, YLE, YOPP, YOPPD, YTE
        INTEGER :: I, ITOPP
        !
        !*** End of declarations rewritten by SPAG
        !
        !---------------------------------------------------------
        !     Adjusts airfoil to scale LE radius by factor RFAC.
        !     Blending of new shape is done with decay length DOC.
        !---------------------------------------------------------
        !
        CALL LEFIND(SLE, X, XP, Y, YP, S, N)
        XLE = SEVAL(SLE, X, XP, S, N)
        YLE = SEVAL(SLE, Y, YP, S, N)
        XTE = 0.5 * (X(1) + X(N))
        YTE = 0.5 * (Y(1) + Y(N))
        CHORD = SQRT((XTE - XLE)**2 + (YTE - YLE)**2)
        !
        !---- set unit chord-line vector
        DXC = (XTE - XLE) / CHORD
        DYC = (YTE - YLE) / CHORD
        !
        SRFAC = SQRT(ABS(RFAC))
        !
        !---- go over each point, changing the y-thickness appropriately
        DO I = 1, N
            !
            !------ coordinates in chord-line axes
            XBAR = (X(I) - XLE) * DXC + (Y(I) - YLE) * DYC
            YBAR = (Y(I) - YLE) * DXC - (X(I) - XLE) * DYC
            !
            SOPP = 2.0 * SLE - S(I)
            SOPP = MAX(SOPP, S(1))
            SOPP = MIN(SOPP, S(N))
            !
            IF (ABS(SOPP / SLE - 1.0)<=1.0E-5) THEN
                SOPP = SLE
            ELSE
                !------- converge on exact opposite point with same chord x value
                DO ITOPP = 1, 12
                    XOPP = SEVAL(SOPP, X, XP, S, N)
                    YOPP = SEVAL(SOPP, Y, YP, S, N)
                    XOPPD = DEVAL(SOPP, X, XP, S, N)
                    YOPPD = DEVAL(SOPP, Y, YP, S, N)
                    !
                    RES = (XOPP - XLE) * DXC + (YOPP - YLE) * DYC - XBAR
                    RESD = XOPPD * DXC + YOPPD * DYC
                    !
                    IF (ABS(RES) / (S(N) - S(1))<1.0E-5) GOTO 305
                    IF (RESD==0.0) EXIT
                    !
                    DSOPP = -RES / RESD
                    SOPP = SOPP + DSOPP
                    !
                    IF (ABS(DSOPP / (S(N) - S(1)))<1.0E-5) GOTO 305
                ENDDO
                WRITE (*, *)                                                 &
                        &'LERSCL: Opposite-point location failed. Continuing...'
                SOPP = 2.0 * SLE - S(I)
            ENDIF
            !
            !------ set point on the opposite side with the same chord x value
            305     XOPP = SEVAL(SOPP, X, XP, S, N)
            YOPP = SEVAL(SOPP, Y, YP, S, N)
            !
            YBAROP = (YOPP - YLE) * DXC - (XOPP - XLE) * DYC
            !
            !------ thickness factor tails off exponentially towards trailing edge
            XOC = XBAR / CHORD
            ARG = MIN(XOC / DOC, 15.0)
            TFAC = 1.0 - (1.0 - SRFAC) * EXP(-ARG)
            !
            !------ set new chord x,y coordinates by changing thickness locally
            YBARCT = 0.5 * (YBAR + YBAROP) + TFAC * 0.5 * (YBAR - YBAROP)
            !
            XNEW(I) = XLE + XBAR * DXC - YBARCT * DYC
            YNEW(I) = YLE + YBARCT * DXC + XBAR * DYC
        ENDDO
        !
    END SUBROUTINE LERSCL
    !*==FLAP.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE FLAP(X, XS, Y, YS, S, N, XF, YF, ADEF, INSID, XNEW, YNEW, NNEW, &
            & TOPS, ATOP, XTOP, YTOP, BOTS, ABOT, XBOT, YBOT)
        use m_spline, only: deval, seval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: ABOT, ADEF, ATOP, BOTS, TOPS, XBOT, XF, XTOP, &
                & YBOT, YF, YTOP
        LOGICAL :: INSID
        INTEGER :: N, NNEW
        REAL, DIMENSION(*) :: S, X, XNEW, XS, Y, YNEW, YS
        !
        ! Local variables
        !
        REAL :: ANG, CA, CHX, CHY, COSD, CRSP, DANG, DSAVG, &
                & DSNEW, FVX, FVY, SA, SB1, SB1P, SB1Q, SB2, &
                & SB2P, SB2Q, SFRAC, SIND, ST1, ST1P, ST1Q, ST2, &
                & ST2P, ST2Q, XB1, XB1NEW, XB2, XB2NEW, XBAR, XT1, &
                & XT1NEW, XT2, XT2NEW, YB1, YB1NEW, YB2, YB2NEW, &
                & YBAR, YT1, YT1NEW, YT2, YT2NEW
        INTEGER :: I, IB1, IB2, IB2Q, INEW, IP, IT1, IT2, IT2Q, &
                & NPADD
        LOGICAL :: LB1NEW, LB2NEW, LT1NEW, LT2NEW
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------
        !     Sets modified airfoil with a deflected flap.
        !     Points may be added/subtracted in the flap
        !     break vicinity to clean things up.
        !     Note: angle deflection ADEF is in radians
        !----------------------------------------------------
        !
        !
        !
        IF (INSID) THEN
            ATOP = MAX(0.0, -ADEF)
            ABOT = MAX(0.0, ADEF)
        ELSE
            CHX = DEVAL(BOTS, X, XS, S, N) - DEVAL(TOPS, X, XS, S, N)
            CHY = DEVAL(BOTS, Y, YS, S, N) - DEVAL(TOPS, Y, YS, S, N)
            FVX = SEVAL(BOTS, X, XS, S, N) + SEVAL(TOPS, X, XS, S, N)
            FVY = SEVAL(BOTS, Y, YS, S, N) + SEVAL(TOPS, Y, YS, S, N)
            CRSP = CHX * (YF - 0.5 * FVY) - CHY * (XF - 0.5 * FVX)
            IF (CRSP>0.0) THEN
                !-------- flap hinge is above airfoil
                ATOP = MAX(0.0, ADEF)
                ABOT = MAX(0.0, ADEF)
            ELSE
                !-------- flap hinge is below airfoil
                ATOP = MAX(0.0, -ADEF)
                ABOT = MAX(0.0, -ADEF)
            ENDIF
        ENDIF
        !
        !---- find upper and lower surface break arc length values...
        CALL SSS(TOPS, ST1, ST2, ATOP, XF, YF, X, XS, Y, YS, S, N, 1)
        CALL SSS(BOTS, SB1, SB2, ABOT, XF, YF, X, XS, Y, YS, S, N, 2)
        !
        !---- ... and x,y coordinates
        XT1 = SEVAL(ST1, X, XS, S, N)
        YT1 = SEVAL(ST1, Y, YS, S, N)
        XT2 = SEVAL(ST2, X, XS, S, N)
        YT2 = SEVAL(ST2, Y, YS, S, N)
        XB1 = SEVAL(SB1, X, XS, S, N)
        YB1 = SEVAL(SB1, Y, YS, S, N)
        XB2 = SEVAL(SB2, X, XS, S, N)
        YB2 = SEVAL(SB2, Y, YS, S, N)
        !
        !
        WRITE (*, 1100) XT1, YT1, XT2, YT2, XB1, YB1, XB2, YB2
        1100 FORMAT (/' Top breaks: x,y =  ', 2F9.5, 4X, &
                &2F9.5/' Bot breaks: x,y =  ', 2F9.5, 4X, 2F9.5)
        !
        !---- find points adjacent to breaks
        DO I = 1, N - 1
            IF (S(I)<=ST1 .AND. S(I + 1)>ST1) IT1 = I + 1
            IF (S(I)<ST2 .AND. S(I + 1)>=ST2) IT2 = I
            IF (S(I)<=SB1 .AND. S(I + 1)>SB1) IB1 = I
            IF (S(I)<SB2 .AND. S(I + 1)>=SB2) IB2 = I + 1
        ENDDO
        !
        DSAVG = (S(N) - S(1)) / FLOAT(N - 1)
        !
        !---- smallest fraction of s increments i+1 and i+2 away from break point
        SFRAC = 0.33333
        !
        IF (ATOP/=0.0) THEN
            ST1P = ST1 + SFRAC * (S(IT1) - ST1)
            ST1Q = ST1 + SFRAC * (S(IT1 + 1) - ST1)
            IF (S(IT1)<ST1Q) THEN
                !-------- simply move adjacent point to ideal SFRAC location
                XT1NEW = SEVAL(ST1Q, X, XS, S, N)
                YT1NEW = SEVAL(ST1Q, Y, YS, S, N)
                LT1NEW = .FALSE.
            ELSE
                !-------- make new point at SFRAC location
                XT1NEW = SEVAL(ST1P, X, XS, S, N)
                YT1NEW = SEVAL(ST1P, Y, YS, S, N)
                LT1NEW = .TRUE.
            ENDIF
            !
            ST2P = ST2 + SFRAC * (S(IT2) - ST2)
            IT2Q = MAX(IT2 - 1, 1)
            ST2Q = ST2 + SFRAC * (S(IT2Q) - ST2)
            IF (S(IT2)>ST2Q) THEN
                !-------- simply move adjacent point
                XT2NEW = SEVAL(ST2Q, X, XS, S, N)
                YT2NEW = SEVAL(ST2Q, Y, YS, S, N)
                LT2NEW = .FALSE.
            ELSE
                !-------- make new point
                XT2NEW = SEVAL(ST2P, X, XS, S, N)
                YT2NEW = SEVAL(ST2P, Y, YS, S, N)
                LT2NEW = .TRUE.
            ENDIF
        ENDIF
        !
        IF (ABOT/=0.0) THEN
            SB1P = SB1 + SFRAC * (S(IB1) - SB1)
            SB1Q = SB1 + SFRAC * (S(IB1 - 1) - SB1)
            IF (S(IB1)>SB1Q) THEN
                !-------- simply move adjacent point
                XB1NEW = SEVAL(SB1Q, X, XS, S, N)
                YB1NEW = SEVAL(SB1Q, Y, YS, S, N)
                LB1NEW = .FALSE.
            ELSE
                !-------- make new point
                XB1NEW = SEVAL(SB1P, X, XS, S, N)
                YB1NEW = SEVAL(SB1P, Y, YS, S, N)
                LB1NEW = .TRUE.
            ENDIF
            !
            SB2P = SB2 + SFRAC * (S(IB2) - SB2)
            IB2Q = MIN(IB2 + 1, N)
            SB2Q = SB2 + SFRAC * (S(IB2Q) - SB2)
            IF (S(IB2)<SB2Q) THEN
                !-------- simply move adjacent point
                XB2NEW = SEVAL(SB2Q, X, XS, S, N)
                YB2NEW = SEVAL(SB2Q, Y, YS, S, N)
                LB2NEW = .FALSE.
            ELSE
                !-------- make new point
                XB2NEW = SEVAL(SB2P, X, XS, S, N)
                YB2NEW = SEVAL(SB2P, Y, YS, S, N)
                LB2NEW = .TRUE.
            ENDIF
        ENDIF
        !
        !c      DSTOP = ABS(S(IT2)-S(IT1))
        !c      DSBOT = ABS(S(IB2)-S(IB1))
        !
        SIND = SIN(ADEF)
        COSD = COS(ADEF)
        !
        !-------------------------------------------------------------------
        !---- initialize accumulator index for new airfoil
        INEW = 0
        !
        !---- upper flap surface
        DO I = 1, IT2
            XBAR = X(I) - XF
            YBAR = Y(I) - YF
            !
            INEW = INEW + 1
            XNEW(INEW) = XF + XBAR * COSD + YBAR * SIND
            YNEW(INEW) = YF - XBAR * SIND + YBAR * COSD
        ENDDO
        !
        !
        IF (ATOP==0.0) THEN
            !------ arc length of newly created surface on top of airfoil
            DSNEW = ABS(ADEF) * SQRT((XT1 - XF)**2 + (YT1 - YF)**2)
            !
            !------ number of points to be added to define newly created surface
            NPADD = INT(1.5 * DSNEW / DSAVG + 1.0)
            !
            IF (NPADD>0) THEN
                !------- add new points along the new surface circular arc segment
                DANG = ADEF / FLOAT(NPADD)
                XBAR = XT1 - XF
                YBAR = YT1 - YF
                DO IP = 1, NPADD
                    ANG = DANG * (FLOAT(IP) - 0.5)
                    CA = COS(ANG)
                    SA = SIN(ANG)
                    !
                    INEW = INEW + 1
                    XNEW(INEW) = XF + XBAR * CA + YBAR * SA
                    YNEW(INEW) = YF - XBAR * SA + YBAR * CA
                ENDDO
            ENDIF
            !
        ELSE
            !------ set point in the corner and possibly two adjacent points
            IF (LT2NEW) THEN
                XBAR = XT2NEW - XF
                YBAR = YT2NEW - YF
                INEW = INEW + 1
                XNEW(INEW) = XF + XBAR * COSD + YBAR * SIND
                YNEW(INEW) = YF - XBAR * SIND + YBAR * COSD
            ENDIF
            !
            INEW = INEW + 1
            XNEW(INEW) = XT1
            YNEW(INEW) = YT1
            !
            IF (LT1NEW) THEN
                INEW = INEW + 1
                XNEW(INEW) = XT1NEW
                YNEW(INEW) = YT1NEW
            ENDIF
            !
        ENDIF
        !
        !
        !---- unchanged portion ahead of flap breaks
        DO I = IT1, IB1
            INEW = INEW + 1
            XNEW(INEW) = X(I)
            YNEW(INEW) = Y(I)
        ENDDO
        !
        !
        IF (ABOT==0.0) THEN
            !------ arc length of newly created surface on top of airfoil
            DSNEW = ABS(ADEF) * SQRT((XB1 - XF)**2 + (YB1 - YF)**2)
            !
            !------ number of points to be added to define newly created surface
            NPADD = INT(1.5 * DSNEW / DSAVG + 1.0)
            !
            IF (NPADD>0) THEN
                !------- add new points along the new surface circular arc segment
                DANG = ADEF / FLOAT(NPADD)
                XBAR = XB1 - XF
                YBAR = YB1 - YF
                DO IP = 1, NPADD
                    ANG = DANG * (FLOAT(IP) - 0.5)
                    CA = COS(ANG)
                    SA = SIN(ANG)
                    !
                    INEW = INEW + 1
                    XNEW(INEW) = XF + XBAR * CA + YBAR * SA
                    YNEW(INEW) = YF - XBAR * SA + YBAR * CA
                ENDDO
            ENDIF
            !
        ELSE
            IF (LB1NEW) THEN
                INEW = INEW + 1
                XNEW(INEW) = XB1NEW
                YNEW(INEW) = YB1NEW
            ENDIF
            !
            INEW = INEW + 1
            XNEW(INEW) = XB1
            YNEW(INEW) = YB1
            !
            IF (LB2NEW) THEN
                XBAR = XB2NEW - XF
                YBAR = YB2NEW - YF
                INEW = INEW + 1
                XNEW(INEW) = XF + XBAR * COSD + YBAR * SIND
                YNEW(INEW) = YF - XBAR * SIND + YBAR * COSD
            ENDIF
            !
        ENDIF
        !
        !
        !---- rotate flap points about the hinge point (XF,YF)
        DO I = IB2, N
            !------ rotated part on flap
            XBAR = X(I) - XF
            YBAR = Y(I) - YF
            !
            INEW = INEW + 1
            XNEW(INEW) = XF + XBAR * COSD + YBAR * SIND
            YNEW(INEW) = YF - XBAR * SIND + YBAR * COSD
        ENDDO
        !
        !
        !---- total number of new points
        NNEW = INEW
        !
        XTOP = XT1
        YTOP = YT1
        XBOT = XB1
        YBOT = YB1
        !
    END SUBROUTINE FLAP
    !*==INSIDE.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! FLAP
    
    
    
    FUNCTION INSIDE(X, Y, N, XF, YF)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL :: XF, YF
        LOGICAL :: INSIDE
        REAL, DIMENSION(*) :: X, Y
        !
        ! Local variables
        !
        REAL :: ANGLE, XB1, XB2, YB1, YB2
        INTEGER :: I, IP
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------
        !     Returns .TRUE. if point XF,YF
        !     is inside contour X(i),Y(i).
        !-------------------------------------
        !
        !---- integrate subtended angle around airfoil perimeter
        ANGLE = 0.0
        DO I = 1, N
            IP = I + 1
            IF (I==N) IP = 1
            XB1 = X(I) - XF
            YB1 = Y(I) - YF
            XB2 = X(IP) - XF
            YB2 = Y(IP) - YF
            ANGLE = ANGLE + (XB1 * YB2 - YB1 * XB2)                              &
                    & / SQRT((XB1**2 + YB1**2) * (XB2**2 + YB2**2))
        ENDDO
        !
        !---- angle = 0 if XF,YF is outside, angle = +/- 2 pi  if XF,YF is inside
        INSIDE = ABS(ANGLE)>1.0
        !
    END FUNCTION INSIDE
    !*==GETXYF.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE GETXYF(X, XP, Y, YP, S, N, TOPS, BOTS, XF, YF)
        use m_spline, only: seval, sinvrt
        use m_userio, only: askr
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: BOTS, TOPS, XF, YF
        INTEGER :: N
        REAL, DIMENSION(*) :: S, X, XP, Y, YP
        !
        ! Local variables
        !
        REAL :: BOTY, TOPY, YREL
        !
        !*** End of declarations rewritten by SPAG
        !
        !
        IF (XF==-999.0) THEN
            XF = 0.
            CALL ASKR('Enter flap hinge x location^', XF)
        ENDIF
        !
        !---- find top and bottom y at hinge x location
        TOPS = S(1) + (X(1) - XF)
        BOTS = S(N) - (X(N) - XF)
        CALL SINVRT(TOPS, XF, X, XP, S, N)
        CALL SINVRT(BOTS, XF, X, XP, S, N)
        TOPY = SEVAL(TOPS, Y, YP, S, N)
        BOTY = SEVAL(BOTS, Y, YP, S, N)
        !
        WRITE (*, 1000) TOPY, BOTY
        1000 FORMAT (/'  Top    surface:  y =', F8.4, &
                &'     y/t = 1.0'/'  Bottom surface:  y =', F8.4, &
                &'     y/t = 0.0')
        !
        IF (YF==-999.0) THEN
            YF = 999.0
            CALL ASKR(&
                    &'Enter flap hinge y location (or 999 to specify y/t)^'&
                    &, YF)
        ENDIF
        !
        IF (YF==999.0) THEN
            YREL = 0.5
            CALL ASKR('Enter flap hinge relative y/t location^', YREL)
            YF = TOPY * YREL + BOTY * (1.0 - YREL)
        ENDIF
        !
    END SUBROUTINE GETXYF
    !*==SSS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    !
    
    !      IF(LGPARM) THEN
    !        CALL NEWPEN(3)
    !        CALL GPARPL(XL0,YL0,0.9*CH,.TRUE.,
    !     &              CHORDB,AREAB,RADBLE,ANGBTE,
    !     &              EI11BA,EI22BA,APX1BA,APX2BA,
    !     &              EI11BT,EI22BT,APX1BT,APX2BT,
    !     &              THICKB,CAMBRB)
    !      ENDIF
    !C
    !      CALL PLFLUSH
    !C
    !      LGEOPL = .TRUE.
    !      NOVER = 0
    !C
    !      RETURN
    !      END
    
    
    
    
    
    SUBROUTINE SSS(SS, S1, S2, DEL, XBF, YBF, X, XP, Y, YP, S, N, ISIDE)
        use m_spline, only: d2val, seval, deval
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        REAL :: DEL, S1, S2, SS, XBF, YBF
        INTEGER :: ISIDE, N
        REAL, DIMENSION(*) :: S, X, XP, Y, YP
        !
        ! Local variables
        !
        REAL :: A11, A12, A21, A22, DET, DS1, DS2, R1, R1SQ, &
                & R1_S1, R2, R2SQ, R2_S2, RR, RRSQ, RR_S1, RR_S2, &
                & RS1, RS2, RSQ, SIND, SSGN, STOT, X1, X1P, X1PP, &
                & X2, X2P, X2PP, XTOT, Y1, Y1P, Y1PP, Y2, Y2P, &
                & Y2PP, YTOT
        REAL, SAVE :: EPS
        INTEGER :: ITER
        !
        !*** End of declarations rewritten by SPAG
        !
        !----------------------------------------------------------------
        !     Returns arc length points S1,S2 at flap surface break
        !     locations.  S1 is on fixed airfoil part, S2 is on flap.
        !     The points are defined according to two cases:
        !
        !
        !     If DEL > 0:  Surface will be eliminated in S1 < s < S2
        !
        !     Returns the arc length values S1,S2 of the endpoints
        !     of the airfoil surface segment which "disappears" as a
        !     result of the flap deflection.  The line segments between
        !     these enpoints and the flap hinge point (XBF,YBF) have
        !     an included angle of DEL.  DEL is therefore the flap
        !     deflection which will join up the points at S1,S2.
        !     SS is an approximate arc length value near S1 and S2.
        !     It is used as an initial guess for the Newton loop
        !     for S1 and S2.
        !
        !
        !     If DEL = 0:  Surface will be created at s = S1 = S2
        !
        !     If DEL=0, then S1,S2 will cooincide, and will be located
        !     on the airfoil surface where the segment joining the
        !     point at S1,S2 and the hinge point is perpendicular to
        !     the airfoil surface.  This will be the point where the
        !     airfoil surface must be broken to permit a gap to open
        !     as a result of the flap deflection.
        !----------------------------------------------------------------
        !
        !---- convergence epsilon
        DATA EPS/1.0E-5/
        !
        STOT = ABS(S(N) - S(1))
        !
        WRITE (*, *) 'sss del ', del
        SIND = SIN(0.5 * ABS(DEL))
        !
        SSGN = 1.0
        IF (ISIDE==1) SSGN = -1.0
        !
        !---- initial guesses for S1, S2
        RSQ = (SEVAL(SS, X, XP, S, N) - XBF)**2 + (SEVAL(SS, Y, YP, S, N) - YBF)**2
        S1 = SS - (SIND * SQRT(RSQ) + EPS * STOT) * SSGN
        S2 = SS + (SIND * SQRT(RSQ) + EPS * STOT) * SSGN
        !
        !---- Newton iteration loop
        DO ITER = 1, 10
            X1 = SEVAL(S1, X, XP, S, N)
            X1P = DEVAL(S1, X, XP, S, N)
            Y1 = SEVAL(S1, Y, YP, S, N)
            Y1P = DEVAL(S1, Y, YP, S, N)
            !
            X2 = SEVAL(S2, X, XP, S, N)
            X2P = DEVAL(S2, X, XP, S, N)
            Y2 = SEVAL(S2, Y, YP, S, N)
            Y2P = DEVAL(S2, Y, YP, S, N)
            !
            R1SQ = (X1 - XBF)**2 + (Y1 - YBF)**2
            R2SQ = (X2 - XBF)**2 + (Y2 - YBF)**2
            R1 = SQRT(R1SQ)
            R2 = SQRT(R2SQ)
            !
            RRSQ = (X1 - X2)**2 + (Y1 - Y2)**2
            RR = SQRT(RRSQ)
            !
            IF (R1<=EPS * STOT .OR. R2<=EPS * STOT) THEN
                S1 = SS
                S2 = SS
                RETURN
            ENDIF
            !
            R1_S1 = (X1P * (X1 - XBF) + Y1P * (Y1 - YBF)) / R1
            R2_S2 = (X2P * (X2 - XBF) + Y2P * (Y2 - YBF)) / R2
            !
            IF (SIND>0.01) THEN
                !
                IF (RR==0.0) RETURN
                !
                RR_S1 = (X1P * (X1 - X2) + Y1P * (Y1 - Y2)) / RR
                RR_S2 = -(X2P * (X1 - X2) + Y2P * (Y1 - Y2)) / RR
                !
                !------- Residual 1: set included angle via dot product
                RS1 = ((XBF - X1) * (X2 - X1) + (YBF - Y1) * (Y2 - Y1)) / RR - SIND * R1
                A11 = ((XBF - X1) * (-X1P) + (YBF - Y1) * (-Y1P))                     &
                        & / RR + ((-X1P) * (X2 - X1) + (-Y1P) * (Y2 - Y1))                 &
                        & / RR - ((XBF - X1) * (X2 - X1) + (YBF - Y1) * (Y2 - Y1)) * RR_S1 / RRSQ - &
                        & SIND * R1_S1
                A12 = ((XBF - X1) * (X2P) + (YBF - Y1) * (Y2P))                       &
                        & / RR - ((XBF - X1) * (X2 - X1) + (YBF - Y1) * (Y2 - Y1)) * RR_S2 / RRSQ
                !
                !------- Residual 2: set equal length segments
                RS2 = R1 - R2
                A21 = R1_S1
                A22 = -R2_S2
            ELSE
                !
                !------- Residual 1: set included angle via small angle approximation
                RS1 = (R1 + R2) * SIND + (S1 - S2) * SSGN
                A11 = R1_S1 * SIND + SSGN
                A12 = R2_S2 * SIND - SSGN
                !
                !------- Residual 2: set vector sum of line segments beteen the
                !-       endpoints and flap hinge to be perpendicular to airfoil surface.
                X1PP = D2VAL(S1, X, XP, S, N)
                Y1PP = D2VAL(S1, Y, YP, S, N)
                X2PP = D2VAL(S2, X, XP, S, N)
                Y2PP = D2VAL(S2, Y, YP, S, N)
                !
                XTOT = X1 + X2 - 2.0 * XBF
                YTOT = Y1 + Y2 - 2.0 * YBF
                !
                RS2 = XTOT * (X1P + X2P) + YTOT * (Y1P + Y2P)
                A21 = X1P * (X1P + X2P) + Y1P * (Y1P + Y2P) + XTOT * X1PP + YTOT * Y1PP
                A22 = X2P * (X1P + X2P) + Y2P * (Y1P + Y2P) + XTOT * X2PP + YTOT * Y2PP
            ENDIF
            !
            DET = A11 * A22 - A12 * A21
            DS1 = -(RS1 * A22 - A12 * RS2) / DET
            DS2 = -(A11 * RS2 - RS1 * A21) / DET
            !
            DS1 = MIN(DS1, 0.01 * STOT)
            DS1 = MAX(DS1, -.01 * STOT)
            DS2 = MIN(DS2, 0.01 * STOT)
            DS2 = MAX(DS2, -.01 * STOT)
            !
            S1 = S1 + DS1
            S2 = S2 + DS2
            IF (ABS(DS1) + ABS(DS2)<EPS * STOT) GOTO 11
        ENDDO
        WRITE (*, *) 'SSS: failed to converge subtending angle points'
        S1 = SS
        S2 = SS
        !
        !
        !---- make sure points are identical if included angle is zero.
        11   IF (DEL==0.0) THEN
            S1 = 0.5 * (S1 + S2)
            S2 = S1
        ENDIF
        !
    END SUBROUTINE SSS
    !*==CLIS.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    
    
    SUBROUTINE CLIS(X, XP, Y, YP, S, N)
        use m_spline, only: curv
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        INTEGER :: N
        REAL, DIMENSION(*) :: S, X, XP, Y, YP
        !
        ! Local variables
        !
        REAL :: CMAX, CV
        INTEGER :: I, IMAX
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     Displays curvatures at panel nodes.
        !-------------------------------------------------------------------
        !
        CMAX = 0.0
        IMAX = 1
        !
        !---- go over each point, calculating curvature
        WRITE (*, 1050)
        DO I = 1, N
            CV = CURV(S(I), X, XP, Y, YP, S, N)
            WRITE (*, 1100) I, X(I), Y(I), CV
            IF (ABS(CV)>ABS(CMAX)) THEN
                CMAX = CV
                IMAX = I
            ENDIF
        ENDDO
        !
        WRITE (*, 1200) CMAX, IMAX, X(IMAX), Y(IMAX)
        !
        RETURN
        !
        1050 FORMAT (/'  i       x        y           curv')
        !CC             120   0.2134  -0.0234      2025.322
        1100 FORMAT (1X, I3, 2F9.4, F14.3)
        1200 FORMAT (/' Maximum curvature =', F14.3, '   at  i,x,y  = ', I3, 2F9.4)
    END SUBROUTINE CLIS
    !*==CANG.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020
    ! CLIS
    
    
    SUBROUTINE CANG(X, Y, N, LPRINT)
        IMPLICIT NONE
        !
        !*** Start of declarations rewritten by SPAG
        !
        ! Dummy arguments
        !
        LOGICAL :: LPRINT
        INTEGER :: N
        REAL, DIMENSION(*) :: X, Y
        !
        ! Local variables
        !
        REAL :: AMAX, ANGL, CROSSP, DX1, DX2, DY1, DY2
        INTEGER :: I, IMAX
        REAL, SAVE :: PI
        !
        !*** End of declarations rewritten by SPAG
        !
        !-------------------------------------------------------------------
        !     LPRINT=t:   Displays all panel node corner angles
        !     LPRINT=f:   Displays max panel node corner angle
        !-------------------------------------------------------------------
        DATA PI/3.1415926535897932384/
        !
        AMAX = 0.0
        IMAX = 1
        !
        !---- go over each point, calculating corner angle
        IF (LPRINT) WRITE (*, 1050)
        DO I = 2, N - 1
            DX1 = X(I) - X(I - 1)
            DY1 = Y(I) - Y(I - 1)
            DX2 = X(I) - X(I + 1)
            DY2 = Y(I) - Y(I + 1)
            !
            !------ allow for doubled points
            IF (DX1==0.0 .AND. DY1==0.0) THEN
                DX1 = X(I) - X(I - 2)
                DY1 = Y(I) - Y(I - 2)
            ENDIF
            IF (DX2==0.0 .AND. DY2==0.0) THEN
                DX2 = X(I) - X(I + 2)
                DY2 = Y(I) - Y(I + 2)
            ENDIF
            !
            CROSSP = (DX2 * DY1 - DY2 * DX1)                                     &
                    & / SQRT((DX1**2 + DY1**2) * (DX2**2 + DY2**2))
            ANGL = ASIN(CROSSP) * (180.0 / PI)
            IF (LPRINT) WRITE (*, 1100) I, X(I), Y(I), ANGL
            IF (ABS(ANGL)>ABS(AMAX)) THEN
                AMAX = ANGL
                IMAX = I
            ENDIF
        ENDDO
        !
        WRITE (*, 1200) AMAX, IMAX, X(IMAX), Y(IMAX)
        !
        RETURN
        !
        1050 FORMAT (/'  i       x        y      angle')
        !CC             120   0.2134  -0.0234   25.322
        1100 FORMAT (1X, I3, 2F9.4, F9.3)
        1200 FORMAT (/' Maximum panel corner angle =', F7.3, '   at  i,x,y  = ', &
                & I3, 2F9.4)
    END SUBROUTINE CANG
end module m_geutil
