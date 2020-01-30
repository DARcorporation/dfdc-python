      SUBROUTINE XYBDMP(LU)
C------------------------------------------------
C     Dumps buffer geometry
C------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DO IEL = 1, NBEL
C
      WRITE(LU,*) 'IEL ',IEL
      IB1 = IBFRST(IEL)
      IB2 = IBLAST(IEL)
      WRITE(LU,*) 'IB1,IB2 ',IB1,IB2
C
C---- Dump nodes
      WRITE(LU,*) 'IB, XB, YB'
      DO IB = IB1, IB2
        XX = XB(IB)
        YY = YB(IB)
        WRITE(LU,*) IB,XX,YY
      ENDDO
C
C---- Dump other associated points
      XX = XBLE(IEL)
      YY = YBLE(IEL)
      WRITE(LU,*) 'XBLE, YBLE', XX,YY
C
      XX = XBTE(IEL)
      YY = YBTE(IEL)
      WRITE(LU,*) 'XBTE, YBTE', XX,YY
C
      ENDDO
C
      RETURN
      END ! XYBDMP

      SUBROUTINE XYPDMP(LU)
C------------------------------------------------
C     Dumps buffer geometry
C------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      DO IEL = 1, NBEL
C
      WRITE(LU,*) 'IEL ',IEL
      IP1 = IPFRST(IEL)
      IP2 = IPLAST(IEL)
      WRITE(LU,*) 'IP1,IP2 ',IP1,IP2
C
C---- Dump nodes
      WRITE(LU,*) 'IP, XP, YP'
      DO IP = IP1, IP2
        XX = XP(IP)
        YY = YP(IP)
        WRITE(LU,*) IP,XX,YY
      ENDDO
C
C---- Dump other associated points
      XX = XPLE(IEL)
      YY = YPLE(IEL)
      WRITE(LU,*) 'XPLE, YPLE', XX,YY
C
      XX = XPTE(IEL)
      YY = YPTE(IEL)
      WRITE(LU,*) 'XPTE, YPTE', XX,YY
C
      ENDDO
C
      RETURN
      END ! XYPDMP


      SUBROUTINE CHKVNC
C--------------------------------------
C     Calculates V.n at all control points
C--------------------------------------
      INCLUDE 'DFDC.INC'
C
      WRITE(*,*) 'Calculating V.n at all control points'
C
C----- for each control point find normal flow for this solution
      WRITE(12,*) 'Normal flow at control points ',NCTOT
      DO IEL = 1, NEL
        IC1 = ICFRST(IEL)        
        IC2 = ICLAST(IEL)        
        WRITE(12,*) '  IEL ',IEL, 'NETYPE ',NETYPE(IEL)
        WRITE(12,*) '    IC  ICTYPE       V.n'
        DO IC = IC1, IC2
          VN = ANC(1,IC)*QC(1,IC)
     &       + ANC(2,IC)*QC(2,IC)
          WRITE(12,100) IC,ICTYPE(IC), VN
        END DO
      END DO
C
 100  FORMAT(1X,I5,3X,I5,3X,G12.6)
C
      END ! CHKVNC


