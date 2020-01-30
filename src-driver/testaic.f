C
      IMPLICIT REAL (A-H,M,O-Z)
C
C------------------------------------------------------------------
C     Primary array limits
C
C       IBX    points in buffer geometry (all elements)
C       NEX    number of elements
C       NRPNX  number of panel local-refinement stations
C
C       IPX    number of panel vertices
C       ICX    number of control points
C
C       NQSPX  number of specified-strength distributions
C       NSEGX  number of mixed-inverse target segments
C
C       NMODX  number of geometry parameters
C
C       NPNTX  number of operating points in a polar
C       NPOLX  number of polars saved
C
C       NRSTX  number of rotor radial stations
C
      PARAMETER (IBX=1000, NEX=30, NRPNX=4)
      PARAMETER (IPX=1000, ICX=600)
      PARAMETER (NQSPX=5, NSEGX=10)
      PARAMETER (NMODX=10)
      PARAMETER (NPNTX=200,NPOLX=10)
      PARAMETER (NRSTX=30)

C
C------------------------------------------------------------------
C     Derived array limits
C
C       NCX    number of control points excluding internal points
C       NSX    number of element BL sides
C       NSYSX  size of inviscid linear system for gamma
C       NUX    number of unit-driver gamma distributions
C       NGLVX  number of global variables
C       NGLPX  number of global parameters
C       NRHSX  number of righthand sides
C
      PARAMETER (NCX = ICX-NEX)
      PARAMETER (NSX = 2*NEX)
      PARAMETER (NSYSX = ICX + NEX)
      PARAMETER (NUX   = IPX + NEX + 2)
      PARAMETER (NGLVX = 4*NSEGX)
      PARAMETER (NGLPX = NMODX)
      PARAMETER (NRHSX = NGLVX + NGLPX + 1)
C
      LOGICAL OK
      CHARACTER*80 LINE
C
C---- local arrays for calling QFCALC
c      PARAMETER (NFX=1)
c      DIMENSION QF(2,NFX),
c     &          QF_GAM(2,IPX,NFX),
c     &          QF_SIG(2,IPX,NFX),
c     &          QF_GAMT(2,2,NEX,NFX),
c     &          QF_SIGT(2,2,NEX,NFX),
c     &          CPF(NFX), CPF_QF(2,NFX), CPF_QINF(NFX)


      DIMENSION XP(IPX), YP(IPX)
      DIMENSION IPFRST(IPX), IPLAST(IPX)
      DIMENSION XF(IPX), YF(IPX)
      DIMENSION JPAN1(IPX), JPAN2(IPX), IFTYPE(IPX)
      DIMENSION IPFO(IPX), IPFP(IPX)
      DIMENSION GAM(IPX), SIG(IPX)
      DIMENSION QF(2,IPX)
      DIMENSION QF_GAM(2,IPX,IPX),
     &          QF_SIG(2,IPX,IPX)

C
      write(*,*) 'Enter coordinates of panel'
      read(*,*) xp(1),yp(1),xp(2),yp(2)
      ipfrst(1) = 1
      iplast(1) = 2
C
C------ we will evaluate only this point
      IP1 = 1
      IP2 = 1
      NF = IP2-IP1+1
C
C------ assume the point is not on a panel
      IPFO(1) = 0
      IPFP(1) = 0
      IFTYPE(1) = 0
C
C------ evaluate velocity components
C        CALL QFCALC(IP1,IP2, XFVELS,YFVELS, IPFO,IPFP, IFTYPE,
C     &              QF,QF_GAM,QF_SIG, QF_GAMT,QF_SIGT)
C
 10   DO I = 1, 10
        GAM(I) = 0.0
        SIG(I) = 1.0
        DO K = 1, 2
          QF(K,I) = 0.0
          QF_GAM(K,1,I) = 0.0
          QF_GAM(K,2,I) = 0.0
          QF_SIG(K,1,I) = 0.0
          QF_SIG(K,2,I) = 0.0
        END DO
      END DO
C
      JFRST = 1
      JLAST = 2
      NF = 1
C
      JPAN1(1) = 0
      JPAN2(1) = 0
      IFTYPE(1) = 0
      JDIM = IPX
C
      write(*,*) 'Enter field point X,Y '
      read(*,*,end=100) XF(1),YF(1)
C
      CALL PANAIC(JFRST,JLAST,XP,YP, 
     &            NF,XF,YF, 
     &            JPAN1,JPAN2, IFTYPE,
     &            GAM, SIG,
     &            JDIM, QF, QF_GAM, QF_SIG )
C
      write(*,*) ' '
      write(*,*) 'Influence velocity' 
      write(*,*) 'QF      ',QF(1,1),QF(2,1)
      write(*,*) 'QF_GAM1 ',QF_GAM(1,1,1),QF_GAM(2,1,1)
      write(*,*) 'QF_GAM2 ',QF_GAM(1,2,1),QF_GAM(2,2,1)
      write(*,*) 'QF_SIG1 ',QF_SIG(1,1,1),QF_SIG(2,1,1)
      write(*,*) 'QF_SIG2 ',QF_SIG(1,2,1),QF_SIG(2,2,1)
      write(*,*) ' '
C
      go to 10
 100  STOP
      END
