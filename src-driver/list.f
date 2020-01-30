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
C     Cleaned up output formats in PRFORC
C
C=========================================================================


      SUBROUTINE PRFORC(IELPRT,LU)
      INCLUDE 'DFDC.INC'
C
      IF(IELPRT.EQ.0) THEN
       IEL1 = 1
       IEL2 = NEL
      ELSE
       IF(IELPRT.LT.0 .OR. IELPRT.GT.NEL) RETURN
       IEL1 = IELPRT
       IEL2 = IELPRT
      ENDIF
C
 1010 FORMAT(//,1X,'Force Summary',/,1X,47('-'))
 1050 FORMAT(1X, A,A,A,A,A)
 1100 FORMAT(1X, A, 5G13.5)
 1200 FORMAT(1X, 5(A, G13.5))
C
      WRITE(LU,1010) 
      WRITE(LU,1050) NAME(1:48), '  '
C
      WRITE(LU,*)
      QUREF = 0.5*RHO*QREF**2
C
      WRITE(LU,1200) 'Qinf   =', QINF, '    Rho  =', RHO 
      WRITE(LU,1200) 'Qref   =', QREF, '    Mach =', MACH
      IF(LVISC) THEN
       WRITE(LU,1200) 'unitRN =', UREN, '   q_ref =',QUREF
      ELSE
       WRITE(LU,1200) 'q_ref  =',QUREF
      ENDIF
C
      DO IEL = 0, NEL
C
        FX = CX(IEL)*QUREF
        FY = CY(IEL)*QUREF
        MZ = CM(IEL)*QUREF
C
        FXV   = CXVIS(IEL)*QUREF
        CXINV = CX(IEL) - CXVIS(IEL)
        FXI   = CXINV*QUREF
C
C---- Exclude wake and vortex wake elements from force reporting
        IF(IEL.EQ.0) THEN 
         WRITE(LU,1300) 
     &                FX,  CX(IEL),
     &                FXI, CXINV,         
     &                FXV, CXVIS(IEL)         
C
        ELSEIF(NETYPE(IEL).EQ.0 .OR.
     &         NETYPE(IEL).EQ.5 .OR.
     &         NETYPE(IEL).EQ.6) THEN
         WRITE(LU,1310) IEL, 
     &                 FX,  CX(IEL),
     &                 FXI, CXINV,         
     &                 FXV, CXVIS(IEL)         
         WRITE(LU,1200) 'QnDOF=', QNDOF(IEL)
C
        ENDIF
C
      ENDDO
C
 1300 FORMAT(1X,
     & '-----------------------------------------------',
     &/' Total Forces',
     &/' Fx    =',G13.5,'     Cx    =',G13.5,
     &/' Fxinv =',G13.5,'     Cxinv =',G13.5,
     &/' Fxvis =',G13.5,'     Cxvis =',G13.5)
cc     &/' Fy    =',G13.5,'     Cy    =',G13.5,
cc     &/' Mz    =',G13.5,'     Cm    =',G13.5)
C
 1310 FORMAT(1X,
     & '-----------------------------------------------'
     &/' Element', I3,
     &/' Fx    =',G13.5,'     Cx    =',G13.5,
     &/' Fxinv =',G13.5,'     Cxinv =',G13.5,
     &/' Fxvis =',G13.5,'     Cxvis =',G13.5)
C
      RETURN
      END ! PRFORC



      SUBROUTINE PRDUMP
      INCLUDE 'DFDC.INC'
C
      CHARACTER*128 FNAME
      CHARACTER*1 ANS
      LOGICAL OK
C
 1000 FORMAT(A)
C
      WRITE(*,1010)
 1010 FORMAT(/'    ---------------------------------------'
     &       /'     I ndices and sizes'
     &       /'     P anel vertices'
     &       /'     C ontrol points'
     &       /'     T railing edge panels')
      WRITE(*,2100)
 2100 FORMAT(/'     Select type of output: ',$)
      READ(*,1000) ANS
C
      IF(INDEX('IPCTipct',ANS).EQ.0) RETURN
C
      CALL ASKS('Enter filename (<return> for screen output)^',FNAME)
      IF(FNAME(1:1).EQ.' ') THEN
       LU = 6
      ELSE
       LU = 8
       OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=5)
       CALL ASKL('Output file exists.  Overwrite?^',OK)
       IF(OK) GO TO 6
C
       CLOSE(LU)
       WRITE(*,*) 'Output not written'
       RETURN
C
 5     OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=90)
 6     REWIND(LU)
      ENDIF
C
C-------------------------------------------------
      IF    (INDEX('Ii',ANS).GT.0) THEN
C------ list element sizes and indices
        WRITE(LU,3070) 
     &    '#    iel netype    ic1    ic2    ip1    ip2'
cc        '#123456712345671234567123456712345671234567
        DO IEL = 1, NEL
          WRITE(LU,3074) IEL,NETYPE(IEL),
     &                   ICFRST(IEL),ICLAST(IEL),
     &                   IPFRST(IEL),IPLAST(IEL)
        ENDDO
        WRITE(LU,1000)
C     
C
C-------------------------------------------------
      ELSEIF(INDEX('Pp',ANS).GT.0) THEN
C------ list panel singularity strengths
        DO IEL = 1, NEL
          WRITE(LU,*) '# element ',IEL,
     &                ' ip1 ',IPFRST(IEL),' ip2 ',IPLAST(IEL)
          WRITE(LU,3070) 
     &    '#      x           y           s       ',
     &    '     gam         sig         gth       '
          DO IP=IPFRST(IEL), IPLAST(IEL)
            WRITE(LU,3072)  XP(IP),  YP(IP),  SP(IP),
     &                     GAM(IP), SIG(IP), GTH(IP)
          ENDDO
          WRITE(LU,1000)
        ENDDO
C     
C-------------------------------------------------
      ELSEIF(INDEX('Cc',ANS).GT.0) THEN
C------ list normal and tangential velocities on left and right sides
        DO IEL = 1, NEL
          WRITE(LU,*) '# element ',IEL,
     &                ' ic1 ',ICFRST(IEL),' ic2 ',ICLAST(IEL)
          WRITE(LU,3070) 
     &  '#      x           y       ',
     &   '     Qnl         Qnr      ',
     &   '     Qtl         Qtr      ',
     &   '     Cpl         Cpr      '
          DO IC=ICFRST(IEL), ICLAST(IEL)
            QNR =  ANC(1,IC)*QCR(1,IC)  + ANC(2,IC)*QCR(2,IC)
            QNL =  ANC(1,IC)*QCL(1,IC)  + ANC(2,IC)*QCL(2,IC)
            QTR = -ANC(2,IC)*QCR(1,IC)  + ANC(1,IC)*QCR(2,IC)
            QTL = -ANC(2,IC)*QCL(1,IC)  + ANC(1,IC)*QCL(2,IC)
C     
            WRITE(LU,3072)  XC(IC),  YC(IC), 
     &                     QNL    , QNR    ,
     &                     QTL    , QTR    ,
     &                     CPL(IC), CPR(IC)
          ENDDO
          WRITE(LU,1000)
        ENDDO
C     
C-------------------------------------------------
      ELSEIF(INDEX('Tt',ANS).GT.0) THEN
C------ list TE panel quantities
        WRITE(LU,3070) 
     &  '#     x1           x2           y1           y2      ',
     &   '    gam1         gam2         sig1         sig2     '
        DO IEL = 1, NEL
          IF(LTPAN(IEL)) THEN
           WRITE(LU,3072) 
     &         XPT (1,IEL),XPT (2,IEL), YPT (1,IEL),YPT (2,IEL),
     &         GAMT(1,IEL),GAMT(2,IEL), SIGT(1,IEL),SIGT(2,IEL)
          ENDIF
        ENDDO
C-------------------------------------------------
      ENDIF
C
      RETURN
C
 3070 FORMAT(A,A,A,A)
 3071 FORMAT(A,I5,A,I5,A,I5,A,I5,A)
 3072 FORMAT(1X,20G13.5)
 3074 FORMAT(1X,20I7)
C     

 90   WRITE(*,*) 'Bad filename.'
      WRITE(*,*) 'Output not written.'
C
      RETURN
      END ! PRDUMP
