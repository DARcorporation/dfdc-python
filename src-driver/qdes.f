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

      SUBROUTINE DFQDES
C------------------------------------------------------
C     Mixed-Inverse design routine. Based on the 
C     same panel formulation as basic analysis method.
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
      CHARACTER*4 COMAND, COMOLD, CIDOTS
      LOGICAL LRECALC, LOFINI
C
      CHARACTER*16 PROMPT
      CHARACTER*128 COMARG, ARGOLD
C
      LOGICAL LPLTEL(NEX)
C
      DIMENSION XQBOX(2), YQBOX(2)
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
      LOGICAL LSAME
C
      SAVE COMOLD, ARGOLD
C
C
      COMAND = '****'
      COMARG = ' '
      LRECALC = .FALSE.
C
      IF(NPTOT.EQ.0) THEN
       WRITE(*,*)
       WRITE(*,*) '***  No paneling available  ***'
       RETURN
      ENDIF
C
C---- startup default: plot all elements
      DO IEL=1, NEL
        LPLTEL(IEL) = .TRUE.
      ENDDO
C
C---- number of sub-intervals for Qspec(s) plotting
      IF(LQSLOP) THEN
       NTQSPL = 10
      ELSE
       NTQSPL = 1
      ENDIF
C
C---- make sure we have a valid target element
      IF(IELQDES.LT.1 .OR. IELQDES.GT.NEL) THEN
       IF(NEL.EQ.1) THEN
        IELQDES = 1
       ELSE
 5      CALL ASKI('Specify target element index^',IELQDES)
        IF(IELQDES.LT.1 .OR. IELQDES.GT.NEL) THEN
         WRITE(*,*) 'Number of elements present:', NEL
         IELQDES = 1
         GO TO 5
        ENDIF
       ENDIF
      ENDIF
C
C---- plot all elements in geometry plot as default case
      IELPLT = 0
C
C---- make sure a current solution exists
      IF(.NOT.LGAMA) THEN
       CALL GAMSOLV
      ENDIF
C
C---- set "current" speed distribution Gamma, arc length, and x/c,y/c arrays
      CALL QSPINI
C
C---- initialize plot limits for all elements
      DO IEL = 1, NEL
        CALL QSPLIM(IEL)
      ENDDO
C
C---- initialize blowup parameters and plot Qspec(s)
      CALL QPLINI(.TRUE.,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
      CALL QSPLOT
C
C---- begin user-interaction loop
C................................................................
C
 500  CONTINUE
      COMOLD = COMAND
      ARGOLD = COMARG
C
 501  CONTINUE
      PROMPT = '.QDES^'
C
      CALL ASKC(PROMPT,COMAND,COMARG)
C
      IF(COMAND.EQ.'    ') THEN
       IF(LPLOT) THEN
        CALL PLEND
        CALL CLRZOOM
       ENDIF
       LPLOT = .FALSE.
C
       LQSPPL = .FALSE.
       LGSPPL = .FALSE.
       RETURN
      ENDIF
C
C---- process previous command ?
      IF(COMAND(1:1).EQ.'!') THEN
        IF(COMOLD.EQ.'****') THEN
          WRITE(*,*) 'Previous .QDES command not valid'
          GO TO 501
        ELSE
          COMAND = COMOLD
          COMARG = ARGOLD
          LRECALC = .TRUE.
        ENDIF
      ELSE
        LRECALC = .FALSE.
      ENDIF
C
C---- extract numbers (if any) from command argument string
      DO I=1, 20
        IINPUT(I) = 0
        RINPUT(I) = 0.0
      ENDDO
      NINPUT = 0
      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
      NINPUT = 0
      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C
C
C--------------------------------------------------------------------
      IF(    COMAND.EQ.'HELP'
     &  .OR. COMAND.EQ.'?   ') THEN
C
        IF(NEL.GT.1) THEN
         CIDOTS = 'ii. '
        ELSE
         CIDOTS = '    '
        ENDIF
C
        WRITE(*,1100) IELQDES, CIDOTS
 1100   FORMAT(
     &    /'   Elem i   Change target element ( currently =', I3, ' )'
     &   //'   QSET     Reset  Qspec <== gamma'
     &   //'   Modi     Modify Qspec'
     &    /'   MARK i   Mark   inverse segment'
     &    /'   UNMA i   Unmark inverse segment'
     &    /'   SLOP     Modified-Qspec endpoint slope matching toggle'
     &    /'   SREV     Reverse plotted s direction'
     &    /'   QREV     Reverse plotted Q direction'
     &   //'   EXEC i   Execute mixed-inverse calculation'
     &    /'   RLX  r   Set maximum under-relaxation factor for update'
     &    /'   QSS      Gamma endpoint curvature constraint toggle'
     &    /'   GMOV     Geometry endpoint movement toggle'
     &    /'   REST     Restore geometry from buffer geometry'
     &    /'   REFL     Reflected-Qspec overlay toggle'
     &   //'   Qplo     Plot Qspec(s) and gamma(s)'
     &    /'   Gplo ',A,'Plot geometry elements'
     &   //'   Blow     Blowup'
     &    /'   Rese     Reset to original plot scale'
     &   //'  .ANNO     Annotate plot'
     &    /'   HARD     Hardcopy current plot'
     &    /'   SIZE r   Change absolute plot size'
     &   //'   Z        Zoom  '
     &    /'   U        Unzoom')
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ELEM'
     &  .OR. COMAND.EQ.'E   ') THEN
    8   IF(NINPUT.GE.1) THEN
         IELQDES = IINPUT(1)
        ELSE
         CALL ASKI('Enter target element index^',IELQDES)
        ENDIF
C
        IF(IELQDES.LT.1 .OR. IELQDES.GT.NEL) THEN
         WRITE(*,*) 'Number of elements present:', NEL
         NINPUT = 0
         GO TO 8
        ENDIF
C
        CALL QPLINI(.FALSE.,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
        CALL QSPLOT
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'QSET') THEN
        CALL GAMQSP(KQTARG,IELQDES)
        CALL QPLINI(.FALSE.,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
        CALL QSPLOT
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'REFL') THEN
        LQREFL = .NOT.LQREFL
        IF(LQREFL) THEN
         WRITE(*,*) 'Reflected Qspec will be plotted'
        ELSE
         WRITE(*,*) 'Reflected Qspec will not be plotted'
         CALL QPLINI(.FALSE.,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
        ENDIF
        CALL QSPLOT
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'MODI'
     &  .OR. COMAND.EQ.'M   ') THEN
        CALL GETCOLOR(ICOL0)
C
C------ put up fresh Qspec(s) plot
        CALL QPLINI(.FALSE.,NEL,IELQDES, LQSREV,LQQREV,
     &               SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &               XQOFF,XQSF,YQOFF,YQSF)
        CALL QSPLOT
C
C------ box in which modification cursor inputs are recognized
        CALL GETZOOMABS(XZOFF,YZOFF,XZFAC,YZFAC)
        XQBOX(1) = MIN(0.0  ,       1.0001*XMARG)/XZFAC - XZOFF
        XQBOX(2) = MIN(XWIND, XPAGE-1.0001*XMARG)/XZFAC - XZOFF
        YQBOX(1) = MAX(0.0  ,       1.0001*YMARG)/YZFAC - YZOFF
        YQBOX(2) = MAX(YWIND, YPAGE-1.0001*YMARG)/YZFAC - YZOFF
C
C------ modify Qspec on target element
        IEL = IELQDES
        NSIDE = 1
        CALL MODIFY(IPX,IPFRST(IEL),IPLAST(IEL),NSIDE,NQSP,
     &              SSPEC,QSPEC,QSPECS, LQSLOP,
     &              IPMOD1,IPMOD2,KSIDE,KQSP,
     &              XQBOX,YQBOX, XQBOX,YQBOX,
     &              XQOFF,YQOFF,XQSF,YQSF, .FALSE.)
C
        IF(IPMOD1.LT.IPMOD2) THEN
C------- spline modified Qspec
         CALL SPLQSP(KQSP,IELQDES)
C
C------- replot just the modified Qspec piece in magenta
         CALL NEWCOLORNAME('MAGENTA')
         CALL QSPPLT(IPMOD1,IPMOD2,KQSP,NTQSPL)
         CALL NEWCOLOR(ICOL0)
C
         CALL PLFLUSH
        ENDIF
C
        COMOLD = COMAND
        ARGOLD = COMARG
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'MARK') THEN
        IF    (NSEG.EQ.0) THEN
         NSEG = 1
         ISEG = 1
         IPSEG1(ISEG) = 1
         IPSEG2(ISEG) = 1
         WRITE(*,*)
         WRITE(*,*) 'New inverse segment declared'
C
        ELSE
         IF(NINPUT.GE.1) THEN
C-------- use number from command argument
          ISEG = IINPUT(1)
C
         ELSE
C-------- no argument was given... show current segments and get segment number
          CALL SEGLIS
C
          ISEG = MIN( NSEG , NSEGX )
          CALL ASKI('Enter number of segment to mark^', ISEG)
         ENDIF
        ENDIF
C
        IF(ISEG.LT.1) GO TO 500
C
C------ segment numbers are consecutive... additional segment must be N+1
        IF(ISEG.GT.NSEG) ISEG = MIN( NSEG+1 , NSEGX )
C
C------ set new number of segments
        NSEG = MAX( NSEG , ISEG )
C
        IF(LQSPPL) THEN
C------- get segment from current Qspec(s) plot
         CALL IQSGET(IELQDES,ISEG)
        ELSE
C------- try geometry plot
         IF(.NOT.LGSPPL) THEN
C-------- but first we need a geometry plot on the screen...
          CALL GPLINI(.TRUE. ,NEL,LPLTEL,
     &                 XPMINE,XPMAXE,YPMINE,YPMAXE,
     &                 XPOFF,XPSF,YPOFF,YPSF)
          CALL GVPLOT(LPLTEL,.TRUE.)
         ENDIF
         CALL IGSGET(LPLTEL,ISEG)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'UNMA') THEN
        IF    (NSEG.EQ.0) THEN
         WRITE(*,*) 'No inverse segments present'
         GO TO 500
C
        ELSEIF(NSEG.EQ.1 .AND. NINPUT.EQ.0) THEN
C------- one segment is present and no segment specified... obvious choice
         ISEG = 1
C
        ELSE
         IF(NINPUT.GE.1) THEN
C-------- unmark segment specified in command argument
          ISEG = IINPUT(1)
C
         ELSE
C-------- no argument was given... show current segments and get segment number
          CALL SEGLIS
C
          ISEG = MIN( NSEG , NSEGX )
          CALL ASKI('Enter number of segment to unmark^', ISEG)
         ENDIF
        ENDIF
C
        IF(ISEG.LT.1 .OR. ISEG.GT.NSEG) THEN
         WRITE(*,*) 'Segment not present.  No action taken.'
        ELSE
C
C------- remove segment ISEG (drop down all segments above it in arrays)
         DO KSEG = ISEG, NSEG-1
           IPSEG1(KSEG) = IPSEG1(KSEG+1)
           IPSEG2(KSEG) = IPSEG2(KSEG+1)
           QSPDOF(1,KSEG) = QSPDOF(1,KSEG+1)
           QSPDOF(2,KSEG) = QSPDOF(2,KSEG+1)
           QSPDOF(3,KSEG) = QSPDOF(3,KSEG+1)
           QSPDOF(4,KSEG) = QSPDOF(4,KSEG+1)
         ENDDO
         NSEG = NSEG-1
C
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SLOP') THEN
        LQSLOP = .NOT.LQSLOP
        IF(LQSLOP) THEN
          WRITE(*,*)
     &      'Modified Qspec piece will be made tangent at endpoints'
        ELSE
          WRITE(*,*)
     &     'Modified Qspec piece will not be made tangent at endpoints'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SREV') THEN
        LQSREV = .NOT.LQSREV
        CALL QPLINI(.TRUE.,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
        IF(LQSPPL) CALL QSPLOT
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'QREV') THEN
        LQQREV = .NOT.LQQREV
        CALL QPLINI(.TRUE.,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
        IF(LQSPPL) CALL QSPLOT
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'RLX ') THEN
        IF(NINPUT.GE.1) THEN
         RLXMAX = RINPUT(1)
        ELSE
         CALL ASKR('Enter max underrelaxation factor^',RLXMAX)
        ENDIF
        IF(RLXMAX.GT.1.0) THEN
         WRITE(*,*) 'RLX > 1  is not recommended'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'QSS ') THEN
        LQCURV = .NOT.LQCURV
        IF(LQCURV) THEN
         WRITE(*,*) 'Gamma curvature constraint used'
        ELSE
         WRITE(*,*) 'Gamma curvature constraint not used'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'GMOV') THEN
        LGNFIX = .NOT.LGNFIX
        IF(LGNFIX) THEN
         WRITE(*,*) 'Geometry endpoints are fixed'
        ELSE
         WRITE(*,*) 'Geometry endpoints are free to move'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'REST') THEN
        CALL PANCOP
        CALL CVPGEN
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'QPLO'
     &  .OR. COMAND.EQ.'Q   ') THEN
        LOFINI = .FALSE.
        CALL QPLINI(LOFINI,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
        CALL QSPLOT
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'GPLO'
     &  .OR. COMAND.EQ.'G   ') THEN
        LOFINI = .NOT.LGSPPL
C
        IF(NINPUT.EQ.0) THEN
C------- no element index specified... plot all elements
         DO IEL=1, NEL
           IF(.NOT.LPLTEL(IEL)) LOFINI = .TRUE.
           LPLTEL(IEL) = .TRUE.
         ENDDO
        ELSE
C------- mark only specified elements for plotting
         DO IEL=1, NEL
           IF(LPLTEL(IEL)) LOFINI = .TRUE.
           LPLTEL(IEL) = .FALSE.
         ENDDO
         DO K=1, NINPUT
           IEL = IINPUT(K)
           IF(IEL.GT.0 .AND. IEL.LE.NEL) THEN
            IF(.NOT.LPLTEL(IEL)) LOFINI = .TRUE.
            LPLTEL(IEL) = .TRUE.
           ENDIF
         ENDDO
        ENDIF
C
        CALL GPLINI(LOFINI,NEL,LPLTEL,
     &              XPMINE,XPMAXE,YPMINE,YPMAXE,
     &              XPOFF,XPSF,YPOFF,YPSF)
        CALL GVPLOT(LPLTEL,.TRUE.)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'BLOW'
     &  .OR. COMAND.EQ.'B   ') THEN
        XWS = XWIND/SIZE
        YWS = YWIND/SIZE
        IF    (LQSPPL) THEN
         LSAME = .FALSE.
         LCURS = .TRUE.
         SHO = 2.0*XWS
         CALL OFFGET(XQOFF,YQOFF,XQSF,YQSF,XWS,YWS, LSAME,LCURS,SHO)
         CALL QPLINI(.FALSE.,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
         CALL QSPLOT
        ELSEIF(LGSPPL) THEN
         LSAME = .TRUE.
         LCURS = .TRUE.
         SHO = 2.0*XWS
         CALL OFFGET(XPOFF,YPOFF,XPSF,YPSF,XWS,YWS, LSAME,LCURS,SHO)
         CALL GPLINI(.FALSE. ,NEL,LPLTEL,
     &                XPMINE,XPMAXE,YPMINE,YPMAXE,
     &                XPOFF,XPSF,YPOFF,YPSF)
         CALL GVPLOT(LPLTEL,.TRUE.)
        ELSE
         WRITE(*,*) 'No active plot to blow up'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'RESE'
     &  .OR. COMAND.EQ.'R   ') THEN
        IF(LGSPPL) THEN
         IF(IELPLT.EQ.0) THEN
          IEL1 = 1
          IEL2 = NEL
         ELSE
          IEL1 = IELPLT
          IEL2 = IELPLT
         ENDIF
         DO IEL = IEL1, IEL2
           I = IPFRST(IEL)
           N = IPLAST(IEL) - IPFRST(IEL) + 1
           CALL MINMAX(N,XP(I),XPMINE(IEL),XPMAXE(IEL))
         ENDDO
         CALL GPLINI(.TRUE. ,NEL,LPLTEL,
     &                XPMINE,XPMAXE,YPMINE,YPMAXE,
     &                XPOFF,XPSF,YPOFF,YPSF)
         CALL GVPLOT(LPLTEL,.TRUE.)
        ELSE
C------- initialize plot limits and replot current target element
         IF(IELQDES.EQ.0) THEN
          IEL1 = 1
          IEL2 = NEL
         ELSE
          IEL1 = IELQDES
          IEL2 = IELQDES
         ENDIF
         DO IEL = IEL1, IEL2
           CALL QSPLIM(IEL)
         ENDDO
         CALL QPLINI(.TRUE.,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
         CALL QSPLOT
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'ANNO') THEN
        IF(LPLOT) THEN
         CALL ANNOT(CH)
        ELSE
         WRITE(*,*) 'No active plot to annotate'
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'HARD') THEN
        IF(LPLOT) CALL PLEND
        LPLOT = .FALSE.
        CALL REPLOT(IDEVRP)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'SIZE') THEN
        IF(NINPUT.GE.1) THEN
         SIZE = RINPUT(1)
        ELSE
         CALL ASKR('Enter plot size^',SIZE)
        ENDIF
C
        IF(LGSPPL) THEN
         CALL GPLINI(.FALSE. ,NEL,LPLTEL,
     &                XPMINE,XPMAXE,YPMINE,YPMAXE,
     &                XPOFF,XPSF,YPOFF,YPSF)
         CALL GVPLOT(LPLTEL,.TRUE.)
        ELSE
         CALL QPLINI(.FALSE.,NEL,IELQDES, LQSREV,LQQREV,
     &            SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &            XQOFF,XQSF,YQOFF,YQSF)
         CALL QSPLOT
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'EXEC') THEN
        IF(NSEG.EQ.0) THEN
         WRITE(*,*) '***  Must mark off target segment first  ***'
         GO TO 500
        ENDIF
C
C------ check if target segment includes stagnation point
        DO ISEG = 1, NSEG
          IPST = 0
          DO IP = IPSEG1(ISEG), IPSEG2(ISEG)-1
            IF(QSGAM(IP).GE.0.0 .AND. QSGAM(IP+1).LT.0.0) IPST = IP
          ENDDO
        ENDDO
        IF(IPST.NE.0) THEN 
         WRITE(*,*)
         WRITE(*,*) 'Target segment cannot include ',
     &              'stagnation point in mixed-inverse.'
         GO TO 500
        ENDIF
C
C------ save current coordinates for restoration if requested
        DO IP=1, NPTOT
          XB(IP) = XP(IP)
          YB(IP) = YP(IP)
          SB(IP) = SP(IP)
          XBS(IP) = XPS(IP)
          YBS(IP) = YPS(IP)
        ENDDO
        NBTOT = NPTOT
CHHY- change to only save first #NBEL elements defined (CB and duct)
        DO IEL = 1, NBEL
          IBFRST(IEL) = IPFRST(IEL)
          IBLAST(IEL) = IPLAST(IEL)
        ENDDO
        LGSAME = .TRUE.
        WRITE(*,*)
        WRITE(*,*) 'Current airfoil saved in buffer airfoil'
C
C------ get number of Newton iterations
        IF(NINPUT.GE.1) THEN
         NITERQ = IINPUT(1)
        ELSE
         NITERQ = 3
         CALL ASKI('Enter max number of iterations^',NITERQ)
        ENDIF
C
C------ set endpoint-fixing flags
        DO ISEG = 1, NSEG
          IPS1 = IPSEG1(ISEG)
          IPS2 = IPSEG2(ISEG)
C
          IEL = IELSEG(ISEG)
C
C-------- fix endoints if requested, or if in interior, or if a closed body
          LNFIX1(ISEG) = LGNFIX .OR. IPS1.GT.IPFRST(IEL) .OR. LBODY(IEL)
          LNFIX2(ISEG) = LGNFIX .OR. IPS2.LT.IPLAST(IEL) .OR. LBODY(IEL)
C
C-------- fix gamma curvatures if requested, but only for interior points
          LDFIX1(ISEG) = LQCURV .AND. IPS1.GT.IPFRST(IEL)
          LDFIX2(ISEG) = LQCURV .AND. IPS2.LT.IPLAST(IEL)
C
C-------- must have at least one endpoint fixed
          IF(.NOT.(LNFIX1(ISEG).OR.LNFIX2(ISEG))) THEN
           LNFIX2(ISEG) = .TRUE.
          ENDIF
        ENDDO
C
        WRITE(*,3050)
 3050   FORMAT(/' segment   endpoints   position_fix   curv_match'
     &         /' -------  ----------   ------------  -----------')
CCC                 2( 2)  123 .. 159      T .. F        T .. F
        DO ISEG = 1, NSEG
          WRITE(*,3100) IELSEG(ISEG), IEL,
     &                  IPSEG1(ISEG), IPSEG2(ISEG),
     &                  LNFIX1(ISEG), LNFIX2(ISEG), 
     &                  LDFIX1(ISEG), LDFIX2(ISEG)
 3100     FORMAT(1X,I3,'(',I2,') ', I4,' ..',I4, 2(6X,L1,' .. ',L1, 2X))
        ENDDO
        WRITE(*,*)
C
C------ do the mixed-inverse iteration
        CALL MIXINV(NITERQ,LPLTEL)
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'Z   ') THEN
        IF(LPLOT) THEN
         CALL USETZOOM(.TRUE.,.TRUE.)
         CALL REPLOT(IDEV)
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'U   ') THEN
        IF(LPLOT) THEN
         CALL CLRZOOM
         CALL REPLOT(IDEV)
        ENDIF
C
C--------------------------------------------------------------------
      ELSE
        WRITE(*,1050) COMAND
 1050   FORMAT(' Command ',A,' not recognized.  Type a " ? " for list.')
        COMAND = '****'
      ENDIF
C
      GO TO 500
C........................................................................
C
      END



      SUBROUTINE SEGLIS
      INCLUDE 'DFDC.INC'
C
      WRITE(*,2020)
      DO I = 1, NSEG
        IP1 = IPSEG1(I)
        IP2 = IPSEG2(I)
        WRITE(*,2021) I, IP1,IP2, XP(IP1),XP(IP2), YP(IP1),YP(IP2)
      ENDDO
      IF(NSEG.LT.NSEGX) THEN
        I = NSEG+1
        WRITE(*,2022) I
      ENDIF
 2020 FORMAT(
     & /'  Current inverse segments:'
     &//'         i1 .. i2       x1   ..   x2          y1   ..   y2   '
     & /'        ---------   -------------------   -------------------')
CCC         3     24   132    -10.1234  -10.1234    -10.1234  -10.1234
CCC         4   (available)
 2021 FORMAT(1X,I3, 1X,2I6, 2X,2F10.4, 2X,2F10.4)
 2022 FORMAT(1X,I3,'   (available)')
C
      RETURN
      END ! SEGLIS



      SUBROUTINE SSPSET(IEL)
      INCLUDE 'DFDC.INC'
C------------------------------------------------------
C     Sets Sspec array for parameterizing Qspec(Sspec)
C------------------------------------------------------
C
C---- don't process non-body elements
      IF( NETYPE(IEL).EQ.3 .OR.
     &    NETYPE(IEL).EQ.4 .OR.
     &    NETYPE(IEL).EQ.5 .OR.
     &    NETYPE(IEL).EQ.6 .OR.
     &    NETYPE(IEL).EQ.7     ) RETURN
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
c      DO IP = IP1, IP2
c        SSPEC(IP) = 0.
c        DO JP = 1, NPTOT
c          SSP_XP(IP,JP) = 0.
c          SSP_YP(IP,JP) = 0.
c        ENDDO
c      ENDDO
C
      DO IP = IP1+1, IP2
        DXP = XP(IP) - XP(IP-1)
        DYP = YP(IP) - YP(IP-1)
        DSP = SQRT(DXP**2 + DYP**2)
C
        SSPEC(IP) = SSPEC(IP-1) + DSP
c        DO JP = IP1, IP2
c          SSP_XP(IP,JP) = SSP_XP(IP-1,JP)
c          SSP_YP(IP,JP) = SSP_YP(IP-1,JP)
c        ENDDO
C
c        IF(DSP.GT.0.0) THEN
c         SSP_XP(IP,IP)   = SSP_XP(IP,IP)   + DXP/DSP
c         SSP_YP(IP,IP)   = SSP_YP(IP,IP)   + DYP/DSP
c         SSP_XP(IP,IP-1) = SSP_XP(IP,IP-1) - DXP/DSP
c         SSP_YP(IP,IP-1) = SSP_YP(IP,IP-1) - DYP/DSP
c        ENDIF
      ENDDO
C
      DSP = SSPEC(IP2)
      DO IP = IP1, IP2
        SSPEC(IP) = SSPEC(IP)/DSP
c        DO JP = IP1, IP2
c          SSP_XP(IP,JP) = (SSP_XP(IP,JP) - SSPEC(IP)*SSP_XP(IP2,JP))/DSP
c          SSP_YP(IP,JP) = (SSP_YP(IP,JP) - SSPEC(IP)*SSP_YP(IP2,JP))/DSP
c        ENDDO
      ENDDO
C
C
C
      CHX = XPTE(IEL) - XPLE(IEL)
      CHY = YPTE(IEL) - YPLE(IEL)
      CHSQ = CHX**2 + CHY**2
      IF(CHSQ.EQ.0.0) THEN
       CHX = 1.0
       CHY = 0.0
       CHSQ = 1.0
      ENDIF
C
      DO IP = IP1, IP2
        XSPOC(IP) = (   (XP(IP)-XPLE(IEL))*CHX
     &                + (YP(IP)-YPLE(IEL))*CHY )/CHSQ
        YSPOC(IP) = (   (YP(IP)-YPLE(IEL))*CHX
     &                - (XP(IP)-XPLE(IEL))*CHY )/CHSQ
      ENDDO
      SSPLE(IEL) = (SPLE(IEL)-SP(IP1)) / (SP(IP2)-SP(IP1))
C
      I = IP1
      N = IP2 - IP1 + 1
      CALL SEGSPL(XSPOC(I),XSPOCS(I),SSPEC(I),N)
      CALL SEGSPL(YSPOC(I),YSPOCS(I),SSPEC(I),N)
C
      RETURN
      END ! SSPSET
C


      SUBROUTINE QSPINI
C--------------------------------------------------------
C     Sets up stuff for mixed-inverse problem.
C     Uses current panel solution as the baseline.
C--------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- set data corresponding to Gamma
      QIGAM = QINF
      CXGAM = CX(0)
      CYGAM = CY(0)
      CMGAM = CM(0)
C
C---- set "old" speed distribution Gamma, arc length, and x/c,y/c arrays
      DO IEL = 1, NEL
        CALL SSPSET(IEL)
      ENDDO
C 
      DO IP = 1, NPTOT
        QSGAM(IP) = GAM(IP) + GTH(IP)
      ENDDO
C
      IF(.NOT.LQSPEC) THEN
C----- initialize Qspec to "old" solution and notify user
       NQSP = 1
       KQTARG = 1
       CALL GAMQSP(KQTARG,IELQDES)
       WRITE(*,1155)
 1155  FORMAT(/' Qspec initialized to current Gamma.'/ )
       LQSPEC = .TRUE.
      ENDIF
C
      RETURN
      END ! QSPINI


      SUBROUTINE SPLQSP(KQSP,IEL)
C------------------------------------------------------
C     Splines Qspec(s).  The end intervals are treated
C     specially to avoid Gibbs-type problems from 
C     blindly splining to the stagnation point.
C------------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
C---- usual spline with natural end BCs
      I = IP1 + 1
      N = IP2 - IP1 - 1
      CALL SPLIND(QSPEC(I,KQSP),QSPECS(I,KQSP),SSPEC(I),N,
     &            -999.0,-999.0)
C
ccC---- pseudo-monotonic spline with simple secant slope calculation
cc      CALL SPLINA(QSPEC(I,KQSP),QSPECS(I,KQSP),SSPEC(I),NI)
C
C---- end intervals are splined separately with natural BCs at
C     the trailing edge and matching slopes at the interior points
C
      I = IP1
      N = 2
      CALL SPLIND(QSPEC(I,KQSP),QSPECS(I,KQSP),SSPEC(I),N,
     &            -999.0,QSPECS(I+1,KQSP))
C
      I = IP2 - 1
      N = 2
      CALL SPLIND(QSPEC(I,KQSP),QSPECS(I,KQSP),SSPEC(I),N,
     &            QSPECS(I,KQSP),-999.0)
C
      RETURN
      END ! SPLQSP

 

      FUNCTION QINCOM(QC,QINF,TKLAM)
C-------------------------------------
C     Sets incompressible speed from
C     Karman-Tsien compressible speed
C-------------------------------------
C
      IF(TKLAM.LT.1.0E-4 .OR. ABS(QC).LT.1.0E-4) THEN
C----- for nearly incompressible case or very small speed, use asymptotic
C      expansion of singular quadratic formula to avoid numerical problems
       QINCOM = QC/(1.0 - TKLAM)
      ELSE
C----- use quadratic formula for typical case
       TMP = 0.5*(1.0 - TKLAM)*QINF/(QC*TKLAM)
       QINCOM = QINF*TMP*(SQRT(1.0 + 1.0/(TKLAM*TMP**2)) - 1.0)
      ENDIF
      RETURN
      END ! QINCOM

  
 
      SUBROUTINE GAMQSP(KQSP,IEL)
C------------------------------------------------
C     Sets Qspec(s,k) from current speed Gamma(s).
C------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      ALQSP(KQSP) = ALGAM
      QIQSP(KQSP) = QIGAM
      CXQSP(KQSP) = CXGAM
      CYQSP(KQSP) = CYGAM
      CLQSP(KQSP) = CLGAM
      CMQSP(KQSP) = CMGAM
C
      DO IP = IPFRST(IEL), IPLAST(IEL)
        QSPEC(IP,KQSP) = QSGAM(IP)
      ENDDO
C
C---- zero out Qspec DOFs
      DO ISEG = 1, NSEGX
        DO IDOF = 1, 4
          QSPDOF(IDOF,ISEG) = 0.0
        ENDDO
      ENDDO
C
      CALL SPLQSP(KQSP,IEL)
C
      RETURN
      END ! GAMQSP


      SUBROUTINE SYMQSP(KQSP,IEL)
C-----------------------------------------
C     Forces symmetry of Qspec(KQSP) array
C-----------------------------------------
      INCLUDE 'DFDC.INC'
C
      ALQSP(KQSP) = ALGAM
      QIQSP(KQSP) = QIGAM
      CXQSP(KQSP) = CXGAM
      CYQSP(KQSP) = CYGAM
      CLQSP(KQSP) = CLGAM
      CMQSP(KQSP) = CMGAM
C
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
C
      SSPMID = 0.5*(SSPEC(IP2) - SSPEC(IP1))
      DO I = IP1, IP1+(IP2+IP1)/2-1
        SSPEC(I) = SSPMID + 0.5*(SSPEC(I)      - SSPEC(IP2-I+IP1     ))
        QSPEC(I,KQSP) =     0.5*(QSPEC(I,KQSP) - QSPEC(IP2-I+IP1,KQSP))
      ENDDO
C
      DO I = (IP2+IP1)/2+1, IP2
        SSPEC(I)      = -SSPEC(IP2-I+IP1)      + 2.0*SSPMID
        QSPEC(I,KQSP) = -QSPEC(IP2-I+IP1,KQSP)
      ENDDO
C
C---- zero out Qspec DOFs
      DO ISEG = 1, NSEGX
        DO IDOF = 1, 4
          QSPDOF(IDOF,ISEG) = 0.0
        ENDDO
      ENDDO
C
      CALL SPLQSP(KQSP,IEL)
C
      WRITE(*,1000) KQSP
 1000 FORMAT(/' Qspec',I2,'  made symmetric')
C
      RETURN
      END ! SYMQSP




      SUBROUTINE MIXINV(NITERQ,LPLTEL)
      INCLUDE 'DFDC.INC'
      LOGICAL LPLTEL(*)
C-----------------------------------------------------
C     Performs NITERQ mixed-inverse Newton iterations
C-----------------------------------------------------
      DIMENSION DGAM(IPX), DNXY(IPX), DSXY(IPX)
      DIMENSION DQSPDOF(4,NSEGX), DQNDOF(NEX)
      DIMENSION XPNEW(IPX), YPNEW(IPX)
C
C---- max Newton changes for dn, dgamma
      DATA DNRLX, DGRLX / 0.02 , 0.5 / 
C
C
C---- make sure direct-problem pointers are set up
      IF(.NOT.LSYSP) CALL SYSP
C
C---- augment system with Qspec DOF variables and geometry-continuity eqn's
      CALL SYSPQ
C
C---- compute, normal vectors along which nodes move
      DO IEL=1, NEL
        CALL ANPSET(IEL)
      ENDDO
C
      KQSP = KQTARG
      DO ISEG = 1, NSEG
        IPS1 = IPSEG1(ISEG)
        IPS2 = IPSEG2(ISEG)
C
        DSPSEG = SSPEC(IPS2) - SSPEC(IPS1)
        DO IP = IPS1, IPS2
          FR1 = (SSPEC(IPS2)-SSPEC(IP  )) / DSPSEG
          FR2 = (SSPEC(IP  )-SSPEC(IPS1)) / DSPSEG
C
          FSPEC(1,IP) =  1.0
          FSPEC(2,IP) = -1.0 + 2.0*FR1
          FSPEC(3,IP) = FR1**3
          FSPEC(4,IP) = FR2**3
C
          GAM(IP) = QSPEC(IP,KQSP) - GTH(IP)
     &            + FSPEC(1,IP)*QSPDOF(1,ISEG)
     &            + FSPEC(2,IP)*QSPDOF(2,ISEG)
     &            + FSPEC(3,IP)*QSPDOF(3,ISEG)
     &            + FSPEC(4,IP)*QSPDOF(4,ISEG)
        ENDDO
      ENDDO
C
C---- Newton iteration loop
      DO 190 ITERQ = 1, NITERQ
        LQSYS = .FALSE.
C
C------ set u,v and their Jacobians w.r.t, gamma,x,y  at all control points
        CALL QAIC(.TRUE.)
C
C------ set up Newton system for Mixed-Inverse problem
        CALL QSYS
C
        WRITE(*,*) 'Solving Newton system...'
C
C------ solve for Newton updates
        CALL LUDCMP(NSYSX,NSYS,SYS(1,1))
        CALL BAKSUB(NSYSX,NSYS,SYS(1,1),RES(1))
C
C------ inactive variables will refer to this dummy Newton delta
        RES(0) = 0.
C
        DGMAX = 0.
        DNMAX = 0.
        IGMAX = 0
        INMAX = 0
        KGMAX = 0
        KNMAX = 0
C
C------ Qspec DOF parameters
        DO ISEG = 1, NSEG
          DO IDOF = 1, 4
            J = JSYSQSP(IDOF,ISEG)
            DQSPDOF(IDOF,ISEG) = -RES(J)
          ENDDO
        ENDDO
C
C------ update node variables
        DO IEL = 1, NEL
          DO IP = IPFRST(IEL), IPLAST(IEL)
            ISEG = KSEGP(IP)
            IF(ISEG.GT.0) THEN
C----------- IP is inside inverse segment... DEL = dn
             J = JSYSDXY(IP)
             DNXY(IP) = -RES(J)
             DSXY(IP) = 0.
             DGAM(IP) = FSPEC(1,IP)*DQSPDOF(1,ISEG)
     &                + FSPEC(2,IP)*DQSPDOF(2,ISEG)
     &                + FSPEC(3,IP)*DQSPDOF(3,ISEG)
     &                + FSPEC(4,IP)*DQSPDOF(4,ISEG)
            ELSE
C----------- IP is outside inverse segment... DEL = dgamma
             J = JSYSGAM(IP)
             DGAM(IP) = -RES(J)
             DNXY(IP) = 0.
             DSXY(IP) = 0.
            ENDIF
C
C---------- set max changes and their indices
            IF(ABS(DNXY(IP)).GT.ABS(DNMAX)) THEN
             DNMAX = DNXY(IP)
             INMAX = IP
             KNMAX = IEL
            ENDIF
            IF(ABS(DGAM(IP)).GT.ABS(DGMAX)) THEN
             DGMAX = DGAM(IP)
             IGMAX = IP
             KGMAX = IEL
            ENDIF
          ENDDO
C
C-------- normal-velocity DOF for closed bodies
          J = JSYSQNC(IEL)
          DQNDOF(IEL) = -RES(J)
        ENDDO
C
C
C------ set underrelaxation factor
        RLX = RLXMAX
C
        IEL = KNMAX
        IF(IEL.GT.0) THEN
         DNLIM = DNRLX * SP(IPLAST(IEL)) - SP(IPFRST(IEL))
         IF(RLX*ABS(DNMAX) .GT. DNLIM) RLX = DNLIM/ABS(DNMAX)
        ENDIF
C
        DGLIM = DGRLX * QREF
        IF(RLX*ABS(DGMAX) .GT. DGLIM) RLX = DGLIM/ABS(DGMAX)
C
C
C------ add Newton update for each element in solution vector
        DO ISEG = 1, NSEG
          DO IDOF = 1, 4
            QSPDOF(IDOF,ISEG) = QSPDOF(IDOF,ISEG)
     &                   + RLX*DQSPDOF(IDOF,ISEG)
          ENDDO
        ENDDO
C
        DO IEL = 1, NEL
          QNDOF(IEL) = QNDOF(IEL) + RLX*DQNDOF(IEL)
        ENDDO
C
        DO IP = 1, NPTOT
          GAM(IP) = GAM(IP) + RLX*DGAM(IP)
          QSGAM(IP) = GAM(IP) + GTH(IP)
c           XP(IP)  = XP(IP)  + RLX*DNXY(IP)*ANP(1,IP)
c           YP(IP)  = YP(IP)  + RLX*DNXY(IP)*ANP(2,IP)
        ENDDO
C
C
C------ move the panel nodes in NMOVE sub-steps
        NMOVE = 20
        DFRAC = RLX/FLOAT(NMOVE)
        DO IMOVE = 1, NMOVE
C-------- apply partial Newton change
          DO IP = 1, NPTOT
            XP(IP) = XP(IP) + DFRAC*DNXY(IP)*ANP(1,IP)
            YP(IP) = YP(IP) + DFRAC*DNXY(IP)*ANP(2,IP)
          ENDDO
C
C-------- spline new shape
          DO IEL = 1, NEL
            CALL XYPSPL(IEL)
          ENDDO
C
C-------- redistribute points along contour to match SSPEC array
          DO IP = 1, NPTOT
            XPNEW(IP) = XP(IP)
            YPNEW(IP) = YP(IP)
          ENDDO
C
          DO ISEG = 1, NSEG
            IEL = IELSEG(ISEG)
            IPS1 = IPSEG1(ISEG)
            IPS2 = IPSEG2(ISEG)
            DSPDSP = (   SP(IPS2) -    SP(IPS1))
     &             / (SSPEC(IPS2) - SSPEC(IPS1))
C
            I = IPFRST(IEL)
            N = IPLAST(IEL) - IPFRST(IEL) + 1
            DO IP = IPS1+1, IPS2-1
              SPI = SP(IPS1) + DSPDSP*(SSPEC(IP) - SSPEC(IPS1))
              DSPI = SPI - SP(IP)
              SPI = SP(IP) + 0.7*DSPI
              XPNEW(IP) = SEVAL(SPI,XP(I),XPS(I),SP(I),N)
              YPNEW(IP) = SEVAL(SPI,YP(I),YPS(I),SP(I),N)
            ENDDO
          ENDDO
C
          DO IP = 1, NPTOT
            XP(IP) = XPNEW(IP)
            YP(IP) = YPNEW(IP)
          ENDDO
C
C---------------------------------------------------------------------
C
C-------- set new node-movement direction vectors
          DO IEL = 1, NEL
            CALL ANPSET(IEL)
          ENDDO
        ENDDO
C
C---------------------------------------------------------------------
C----- A new geometry is available, repanel rotor and regrid wakes
C
C---- Adjust paneling on airfoils for rotor and wake grid
cc          CALL ADJPANL
C---- Set up rotor wake grid from repaneled airfoils 
        CALL INIGRD
C---- Reset points on rotor line for new foil geometry
        CALL ROTORINIT
C---- Set up rotor and wake elements
cc        CALL SETROTWK
C---- Initialize rotor/wake pointers 
cc        CALL ROTPINIT

C
C------ recompute control points, normal vectors
        DO IEL=1, NEL
C-------- but don't do non-body singularities
          IF(NETYPE(IEL).EQ.0 .OR.
     &       NETYPE(IEL).EQ.1 .OR.
     &       NETYPE(IEL).EQ.2 .OR.
     &       NETYPE(IEL).EQ.5 .OR.
     &       NETYPE(IEL).EQ.6 .OR.
     &       NETYPE(IEL).EQ.7     )  THEN
           CALL XYCSET(IEL)
           CALL ANPSET(IEL)
          ENDIF
        ENDDO
C
C
        WRITE(*,4010) DNMAX, INMAX, RLX,
     &                DGMAX, IGMAX
        WRITE(*,4022) (DQNDOF(IEL), IEL=1, NEL)
        DO ISEG = 1, NSEG
          WRITE(*,4025) (DQSPDOF(IDOF,ISEG), IDOF=1, 4), 
     &                   IPSEG1(ISEG), IPSEG2(ISEG)
        ENDDO
        WRITE(*,4032) (QNDOF(IEL), IEL=1, NEL)
        DO ISEG = 1, NSEG
          WRITE(*,4035) (QSPDOF(IDOF,ISEG), IDOF=1, 4), 
     &                   IPSEG1(ISEG), IPSEG2(ISEG)
        ENDDO
C
 4010   FORMAT(/' dNmax =', E11.3, '   (', I4,' )', '   RLX =', F8.4
     &         /' dGmax =', E11.3, '   (', I4,' )' )
 4022   FORMAT( ' dQnorm  =', 9E11.3 )
 4025   FORMAT( ' dQspdof =', 4E11.3, '   (', I4,' ..', I4,' )' )
C
 4032   FORMAT(/'  Qnorm  =', 9E11.3 )
 4035   FORMAT( '  Qspdof =', 4E11.3, '   (', I4,' ..', I4,' )' )
C
C
        IF    (LQSPPL) THEN
         CALL QPLINI(.FALSE.,NEL,IELQDES, LQSREV,LQQREV,
     &              SSPMIN,SSPMAX,QSPMIN,QSPMAX,
     &              XQOFF,XQSF,YQOFF,YQSF)
         CALL QSPLOT
        ELSEIF(LGSPPL) THEN
         CALL GPLINI(.FALSE.,NEL,LPLTEL,
     &              XPMINE,XPMAXE,YPMINE,YPMAXE,
     &              XPOFF,XPSF,YPOFF,YPSF)
         CALL GVPLOT(LPLTEL,.TRUE.)
        ENDIF
C
C------ current gamma is now valid (up to Newton residual error)
        LGAMA = .TRUE.
C
C------ velocity AIC matrix is no longer valid since geometry changed
        LQAIC = .FALSE.
        LQGIC = .FALSE.
        LGSYS = .FALSE.
        LGAMU = .FALSE.
 190  CONTINUE
C
C---- set AIC matrix for final geometry
      CALL QAIC(.FALSE.)
C
C---- set velocities on both sides of panels at all control points
      CALL QQCALC(NCTOT,ANC, IPCO,IPCP, GAM,SIG,GTH, QC, QCL,QCR )
C
C---- set velocities and Cp on both sides of panels at all control points
      CALL CPCALC(NCTOT, QCL, QINF,QREF, CPL,CPL_QCL,CPL_QINF)
      CALL CPCALC(NCTOT, QCR, QINF,QREF, CPR,CPR_QCR,CPR_QINF)
C
C---- integrate Cp to get forces and moment
      COSA = 1.0
      SINA = 0.0
      CX(0) = 0.
      CY(0) = 0.
      CM(0) = 0.
      CD(0) = 0.
      DO IEL = 1, NEL
        I = ICFRST(IEL)
        N = ICLAST(IEL) - ICFRST(IEL) + 1
        IF(N.GE.1) THEN
         CALL FCOEFF(N,CPL(I),CPR(I),XC(I),YC(I),ANC(1,I),DSC(I), 
     &               CX(IEL),CY(IEL),CM(IEL))
         CD(IEL) = CX(IEL)*COSA + CY(IEL)*SINA
C     
         CX(0) = CX(0) + CX(IEL)
         CY(0) = CY(0) + CY(IEL)
         CM(0) = CM(0) + CM(IEL)
         CD(0) = CD(0) + CD(IEL)
        ENDIF
      ENDDO
C
      QIGAM = QINF
      CXGAM = CX(0)
      CYGAM = CY(0)
      CMGAM = CM(0)
C
      LGSAME = .FALSE.
C
      RETURN
      END ! MIXINV
