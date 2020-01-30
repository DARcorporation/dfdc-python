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
C
C-----Version 070-ES1
C
C	Philip Carter, Esotec Developments, February 2009
C	philip (at) esotec (dot) org
C
C       Changes from O.70:
C	ANAME character variable for geometry name (lines 95, 177)
C	ALFAR(IRX,NRX) to store alfa data (line 578)
C       AZERO(IRX,NRX), interpolated A0 for each station (line 646)
C       LCOORD, true if Beta and Alfa output in local coordinates (119,178)
C       LLOAD, true if a case has been loaded from disk (119,178)
C       SVERSION, character string version name (101,126) 
C       (real number VERSION remains for file compatibility)
C
C----Version 070-ES2
C
C	CHRP(IRX,NRX), BETARP (IRX,NRX) (line 597)
C	Chord and Beta at rotor panel boundaries (for geometry output)
C
C       ESLOFT PARAMETERS, declarations and COMMON at end.
C
C=========================================================================
C
      IMPLICIT REAL (A-H,M,O-Z)
C
      PARAMETER (PI=3.1415926535897932384)
      PARAMETER (PI2I=0.5/PI)
      PARAMETER (DTR=PI/180.)
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
      PARAMETER (IBX=1000, NEX=30, NRPNX=4)
      PARAMETER (IPX=2000, ICX=2000)  ! modified from 1000 in 0.70
      PARAMETER (NQSPX=5, NSEGX=10)
      PARAMETER (NMODX=10)
      PARAMETER (NPNTX=200,NPOLX=10)
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
C------------------------------------------------------------------
C
C==================================================================
C=============THE FOLLOWING USED TO BE IN PINDEX.INC===============
C==================================================================
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

C==================================================================
C==================================================================
C==================================================================

C
      CHARACTER*2 PNUM
      CHARACTER*80 NAME,ISPARS,ARGP1,ARGP2,PFILE,CPFILE,FRFILE,ANAME
      CHARACTER*64 PREFIX,SVERSION
      LOGICAL LDBG
      LOGICAL LREVEL, LPGRID, LXYMOV, 
     &        LAGRID, LGGRID, LGTICK, LGPARM, LPSTAG,
     &        LPSHO, LWSHO, LPPLT, LVPLT, LSPLT, LEPLT, LDPLT, LDINT,
     &        LNORM, LRSPC, LRSPCDEF, LRSPCUSR, LPACC,
     &        LQAIC, LQGIC, LVAIC, LVGIC,
     &        LQCNT, LSYSP, LGSYS, LQSYS,
     &        LVMAV, LNCVP, 
     &        LGAMU, LGAMA, 
     &        LSIGP, LSIGM,
     &        LPREF, LFREF,
     &        LNDOF, LTPAN,
     &        LBODY, LXBOD,
     &        LV1ZR, LV2ZR,
     &        LQNZR, LUSET
      LOGICAL LQSPEC,LQSPPL,LGSPPL,LQSVIS,LQREFL,
     &        LQSLOP,LQCURV,LQSREV,LQQREV,
     &        LGNFIX,LGSAME,LGSLOP,LGSYMM,
     &        LGPLOT,LHOMPL,
     &        LNFIX1,LNFIX2,
     &        LDFIX1,LDFIX2,
     &        LCOORD,LLOAD
C
C
      COMMON /CR_GEN/ VERSION, ONE(IBX),ZER(IBX),SVERSION
      COMMON /CI_GEN/ 
     &  NNAME,        ! number of characters in case name
     &  NPREFIX,      ! number of characters in name prefix string
     &  ITMAX,        ! max number of viscous Newton iterations
     &  IMATYP,       ! type of Ma(CL) variation
     &  LUNDBG,       ! logical unit # for debugging output
     &  IEWAKE(NEX),  ! index of wake attached to element
     &  IZERO(IPX)    ! handy dummy zero array
      COMMON /CL_GEN/
     &  LDBG,    ! Debugging flag (set to T for debug output)
     &  LREVEL(NEX), ! element is to be written out in reverse order ?
     &  LPGRID,  ! show grid on Cp(x) plots ?
     &  LXYMOV,  ! x,y panel nodes will be moved ?  (optimization, etc)
     &  LAGRID,  ! show grid on aero plots ?
     &  LGGRID,  ! show grid on geometry plots ?
     &  LGTICK,  ! show tic marks on geometry plots ?
     &  LGPARM,  ! show geometry parameters on geometry plots ?
     &  LPSTAG,  ! plot stagnation points ?
     &  LPSHO,   ! plot point strengths ?
     &  LWSHO,   ! plot wake strengths ?
     &  LPPLT,   ! plot point-singularities ?
     &  LVPLT,   ! plot panel vertices ?
     &  LSPLT,   ! plot symbols on line plots ?
     &  LEPLT,   ! plot element numbers ?
     &  LDINT,   ! compute plotted displacement surface by integration ?
     &  LNORM,   ! normalize airfoil when read in ?
     &  LRSPC,   ! respace input geometry to default or user specified distrib.
     &  LRSPCDEF,  ! respacing data defined?
     &  LRSPCUSR,  ! user input of respacing data?
     &  LPACC,   ! auto accumulate polar sequence ?
     &  LNCVP,   ! plot normal vectors at panel nodes ?
     &  LQAIC,   ! velocity aero-influence matrices dQ/dgam,dQ/dsig available ?
     &  LQGIC,   ! velocity geom-influence matrices dQ/dx  ,dQ/dy   available ?
     &  LVAIC,   ! surface-gamma aero-influence matrices available ?
     &  LVGIC,   ! surface-gamma geom-influence matrices available ?
     &  LQCNT,   ! control points computed ?
     &  LSYSP,   ! system pointers set up ?
     &  LGSYS,   ! system for gamma set up ?
     &  LQSYS,   ! system for specified q(s) (mixed-inverse case) set up ?
     &  LVMAV,   ! if F, VMavg is not initialized 
     &  LGAMU,   ! if F, unit-driver gamma solutions are all invalid
     &  LGAMA,   ! gamma for specified driving Qinf,etc  available ?
     &  LSIGP,   ! pointers for source-influence matrix set up ?
     &  LSIGM,   ! source-influence matrix available ?
     &  LPREF,   ! plot Cp data ?
     &  LFREF,   ! plot force data ?
     &  LNDOF(NEX),  ! does element have a normal-velocity DOF ?
     &  LTPAN(NEX),  ! does element have a TE panel ?
     &  LBODY(NEX),  ! is element a closed body ?
     &  LXBOD(NEX),  ! is element an axisymmetric body on axis?
     &  LV1ZR(NEX),  ! first panel node on element gets zero gamma
     &  LV2ZR(NEX),  ! last  panel node on element gets zero gamma
     &  LQNZR(NEX),  ! normal-velocity DOF of element is set to zero explicitly
     &  LUSET(0:NUX),
     &  LCOORD,      ! Print alfa and beta listings in local coordinates?
     &  LLOAD        ! A case has been loaded from disk
C
      COMMON /CC_GEN/ NAME, ISPARS, ARGP1, ARGP2, PFILE, PNUM(99),
     &                CPFILE,FRFILE,
     &                PREFIX, ANAME
C
      COMMON /CR_PAN/
     &  SRPN1(NRPNX,NEX),
     &  SRPN2(NRPNX,NEX),
     &  CRRAT(NRPNX,NEX),
     &  CVEX(NEX), 
     &  SMOF(NEX), 
     &  FSLE(NEX),
     &  FSTE(NEX),
     &  FPANDEF
      COMMON /CI_PAN/ NPAN(NEX), NPANDEF
C
      COMMON /CR_PLT/
     &  CPLMAX,CPLMIN,CPLDEL,CPLFAC,
     &  QPLMAX,QPLMIN,QPLDEL,QPLFAC,
     &  XPLMAX,XPLMIN,XPLDEL,XPLFAC,
     &  YPLMAX,YPLMIN,YPLDEL,YPLFAC,
     &  XOFA,YOFA,FACA, CPXAR,CPYAR,
     &  XOFG,YOFG,GSF,
     &  PVFAC,QVFAC,
     &  SHF(0:7),SSIZEL(0:7),
     &  XELNUM(NEX),
     &  YELNUM(NEX),
     &  XQOFF,XQSF,
     &  YQOFF,YQSF,
     &  XPOFF,XPSF,
     &  YPOFF,YPSF,
     &  XBOFF,XBSF,
     &  YBOFF,YBSF
      COMMON /CI_PLT/ ICOLEL(NEX), ISPLOT(NEX), IPTYPE
C
      COMMON /CR_GEO/
     &  XB (IBX), YB (IBX), SB(IBX), 
     &  XBS(IBX), YBS(IBX),
     &  XBLE(NEX),YBLE(NEX),SBLE(NEX),
     &  XBTE(NEX),YBTE(NEX),
     &  AREA2DA(NEX),XBCEN2DA(NEX),YBCEN2DA(NEX),
     &  EIXX2DA(NEX),EIYY2DA(NEX), EIXY2DA(NEX),
     &  AREA2DT(NEX),XBCEN2DT(NEX),YBCEN2DT(NEX),
     &  EIXX2DT(NEX),EIYY2DT(NEX), EIXY2DT(NEX),
     &  VOLUMV(NEX), ASURFV(NEX),  XBCENV(NEX), YBCENV(NEX),
     &  RGYRXV(NEX), RGYRYV(NEX),
     &  VOLUMVT(NEX),ASURFVT(NEX),XBCENVT(NEX),YBCENVT(NEX),
     &  RGYRXVT(NEX), RGYRYVT(NEX),
     &  XBREFE(NEX),YBREFE(NEX),
     &  XBMINE(NEX),XBMAXE(NEX),
     &  YBMINE(NEX),YBMAXE(NEX),
     &  XBMIN,XBMAX,
     &  YBMIN,YBMAX,
     &  DXBSUM(NEX),DYBSUM(NEX),AGBSUM(NEX),
     &  XFBSUM(NEX),YFBSUM(NEX),
     &  DTEBOD,DTEPAN,
     &  XWBOX(2),YWBOX(2),
     &  XGBOX(2),YGBOX(2)
      COMMON /CI_GEO/ 
     &  NEL,          ! number of current geometry elements
     &  NBEL,         ! number of buffer  geometry elements
     &  NBTOT,        ! number of buffer  geometry nodes
     &  IBFRST(NEX),  ! index of first buffer geometry node in element
     &  IBLAST(NEX),  ! index of last  buffer geometry node in element
     &  NBTYPE(NEX),  ! type of element in buffer geometry (0 or 5)
     &  NETYPE(NEX),  ! type of element
     &  IELGDES,      ! current target element in GDES
     &  NMOD          ! number of geometry modes defined
C
C  NETYPE  Types          Singularities
C   0     solid surface      unknown vorticity
C   1     wake               source+vorticity specified
C   2     axis line          specified source and doublet
C   3     line or ring       specified source or doublet
C   4     point              specified source or doublet
C   5     source line        specified source only
C   6     rotor source line  specified source only
C   7     vortex wake line   specified vorticity only
C
      COMMON /CR_FLO/ 
     &  QREF, QINF,                ! reference speed, freestream speed
     &  MACH, MACH1,               ! reference Mach, freestream Mach
     &  XTR1(NEX),XTR2(NEX),       ! top, bottom trip locations
     &  CX(0:NEX),CXU(0:NEX,NUX),  ! x-force coeff., dCx/dU
     &  CY(0:NEX),CYU(0:NEX,NUX),  ! y-force coeff., dCy/dU
     &  CM(0:NEX),CMU(0:NEX,NUX),  ! moment  coeff., dCm/dU
     &  CD(0:NEX),CDU(0:NEX,NUX),  ! drag coeff., dCD/dU
cc     &  SOUR(0:NEX),SOURU(0:NEX,NUX),   ! prescribed source
cc     &  CIRC(0:NEX),CIRCU(0:NEX,NUX),   ! prescribed circulation
cc     &  DBLX(0:NEX),DBLXU(0:NEX,NUX),   ! prescribed x-doublet
cc     &  DBLY(0:NEX),DBLYU(0:NEX,NUX),   ! prescribed y-doublet
     &  GAM(IPX), GAMU(IPX,0:NUX),      ! panel-node vortex strength, dGAM/dU
     &  SIG(IPX), SIGU(IPX,0:NUX),      ! panel-node source strength, dSIG/dU
     &  GTH(IPX), GTHU(IPX,0:NUX),      ! node wake  vortex strength, dGTH/dU
     &  QNDOF(NEX), QNDOFU(NEX,0:NUX),  ! elment normal-velocity DOF, dQn/dU
     &  GAMVSP(IPX),                    ! vortex strength from BL
     &  SIGVSP(IPX),                    ! source strength from BL
     &  GAMSET(NEX), GAMP1(NEX),GAMP2(NEX),GAMPD(NEX), ! prescribed GAM
     &  SIGSET(NEX), SIGP1(NEX),SIGP2(NEX),SIGPD(NEX), ! prescribed SIG
     &  QINFP1,QINFP2,QINFPD            ! Qinf sequence first, last, delta
      COMMON /CI_FLO/ 
     &  NU,           ! number of driving variables U
     &  IUQINF,       ! prescribed freestream comp. column index
     &  IUVWK(NEX),   ! prescribed blade circulation B*Circ at rotor stations
     &  IUGAM(IPX),   ! prescribed panel node gamma   "      "
     &  IUSIG(IPX),   ! prescribed panel node sigma   "      "
     &  NPSEQ         ! number of points in a operating sequence
C
      COMMON /CR_CVP/ 
     &  XP(IPX), YP(IPX),                ! panel nodes
     &  SP(IPX),                         ! arc length at panel nodes
     &  XPS(IPX), YPS(IPX),              ! spline derivatives dXP/dSP, ..
     &  XC(ICX), YC(ICX),                ! control points
     &  XPLE(NEX),YPLE(NEX),SPLE(NEX),   ! LE panel nodes
     &  XPTE(NEX),YPTE(NEX),             ! TE panel nodes
     &  ANP(2,IPX), ANP_XYM(2,2,IPX),    ! unit normal at XP,YP
     &              ANP_XYO(2,2,IPX),    
     &              ANP_XYP(2,2,IPX),    
     &  ANC(2,ICX), ANC_DXY(2,2,ICX),    ! unit normal at XC,YC
     &  DSC(ICX)  , DSC_DXY(2,ICX),      ! length of panel
     &  QCR(2,ICX),QCRU(2,ICX,NUX),      ! control point vel. to left  of panel
     &  QCL(2,ICX),QCLU(2,ICX,NUX),      ! control point vel. to right of panel
     &  CPR(ICX),CPRU(ICX,NUX),          ! control point Cp to left  of panel
     &  CPL(ICX),CPLU(ICX,NUX),          ! control point Cp to right of panel
     &  CPR_QCR(2,ICX),CPR_QINF(ICX),
     &  CPL_QCL(2,ICX),CPL_QINF(ICX),
     &  QC(2,ICX),                       ! control point velocity
     &  QC_GAM(2,IPX,ICX), 
c     &  QC_GAM_XP(2,IPX,-1:1,ICX), 
c     &  QC_GAM_YP(2,IPX,-1:1,ICX), 
     &  QC_SIG(2,IPX,ICX), 
c     &  QC_SIG_XP(2,IPX,-1:1,ICX), 
c     &  QC_SIG_YP(2,IPX,-1:1,ICX), 
     &  QC_XP(2,IPX,ICX), 
     &  QC_YP(2,IPX,ICX),
     &  QC_XC(2,ICX), 
     &  QC_YC(2,ICX), 
     &  QCU(2,ICX,0:NUX),
     &  QC_GTH(2,IPX,ICX),               ! velocity influence from vortex wake vorticity 
C
     &  VMAVG(IPX),                      ! mean surface Vm velocity on vortex wake points 
C
     &  XPREFE(NEX),YPREFE(NEX),         ! element reference point
     &  XPMINE(NEX),XPMAXE(NEX),
     &  YPMINE(NEX),YPMAXE(NEX),
     &  XPMIN,XPMAX,
     &  YPMIN,YPMAX,
     &  XPCENT(NEX),YPCENT(NEX),
     &  XSTG(NEX),YSTG(NEX)           ! stagnation point location 
      COMMON /CI_CVP/
     &  NPTOT, IPFRST(NEX), IPLAST(NEX), ! # panel nodes, first,last of element
     &  NCTOT, ICFRST(NEX), ICLAST(NEX), ! # control points
     &  IPCO(ICX), IPCP(ICX),         ! panel nodes on each side of c.p. IC
     &  ICTYPE(ICX),                  ! 0 = flow tangency point, 1 = other
     &  IPTE1(NEX), IPTE2(NEX),       ! TE nodes of element
     &  ICSTG(NEX)                    ! index of panel containing stagnation pt
C
      COMMON /CR_TEP/
     &  XPT(2,NEX), XPT_XP(2,2,NEX),  ! TE panel nodes
     &  YPT(2,NEX), YPT_YP(2,2,NEX),  
     &  GAMT(2,NEX), GAMT_GAM(2,2,NEX), GAMT_SIG(2,2,NEX), ! TE panel strengths
     &  SIGT(2,NEX), SIGT_GAM(2,2,NEX), SIGT_SIG(2,2,NEX),
     &  GAMT_XP(2,2,NEX), GAMT_YP(2,2,NEX),
     &  SIGT_XP(2,2,NEX), SIGT_YP(2,2,NEX),
     &  GAMT_DX(2,NEX), GAMT_DY(2,NEX),
     &  SIGT_DX(2,NEX), SIGT_DY(2,NEX)
C
cc     &  QC_GAMT(2,2,NEX,ICX), QC_SIGT(2,2,NEX,ICX)
C
      COMMON /CR_QSP/ 
     &  RLXMAX,
     &  ALGAM, QIGAM,
     &  CXGAM, CYGAM,
     &  CLGAM, CMGAM,
     &  ALQSP(NQSPX), QIQSP(NQSPX),
     &  CXQSP(NQSPX), CYQSP(NQSPX), 
     &  CLQSP(NQSPX), CMQSP(NQSPX),
     &  SSPEC(IPX), QSPEC(IPX,NQSPX), QSPECS(IPX,NQSPX),
     &  QSGAM(IPX),
     &  XSPOC(IPX), XSPOCS(IPX),
     &  YSPOC(IPX), YSPOCS(IPX),
     &  SSPLE(NEX), 
     &  SSPMIN(NEX),QSPMIN(NEX), 
     &  SSPMAX(NEX),QSPMAX(NEX), 
     &  SSPDEL(NEX),QSPDEL(NEX),
     &  GAMSP(IPX),
     &  SIGSP(IPX),
     &  FSPEC(4,IPX),
     &  QSPDOF(4,NSEGX)
      COMMON /CI_QSP/
     &  NQSP, KQTARG, IELQDES,
     &  KSEGP(IPX), NSEG,
     &  IPSEG1(NSEGX), IPSEG2(NSEGX), IELSEG(NSEGX)
      COMMON /CL_QSP/
     &  LQSPEC,
     &  LQSPPL,
     &  LGSPPL,
     &  LQSVIS,
     &  LQREFL,
     &  LQSLOP,
     &  LQCURV,
     &  LQSREV,
     &  LQQREV,
     &  LGNFIX,
     &  LGSAME,
     &  LGSLOP,
     &  LGSYMM,
     &  LGPLOT,
     &  LHOMPL,
     &  LNFIX1(NSEGX),
     &  LNFIX2(NSEGX),
     &  LDFIX1(NSEGX),
     &  LDFIX2(NSEGX)
C
      COMMON /CR_SYS/ SYS(NSYSX,0:NSYSX), RES(0:NSYSX), 
     &                AICGAM(NSYSX,0:IPX),
     &                AICSIG(NSYSX,0:IPX),
     &                AICQFF(NSYSX,0:1),
     &                AICXYP(NSYSX,2,0:IPX),
     &                AICGTH(NSYSX,0:IPX),
     &                AICVWK(NSYSX,0:NEX)
      COMMON /CI_SYS/ 
     &  NSYS,               ! inviscid system size
     &  NAICGAM,            ! number of forced-gamma nodes
     &  NAICSIG,            ! number of forced-sigma nodes
     &  NAICXYP,            ! number of forced-geometry nodes
     &  NAICGTH,            ! number of wake gamma nodes
     &  NAICVWK,            ! number of rotor vortex wakes 
     &  IELDOFC(ICX),       ! system index of element's normal-velocity DOF
     &  KSYSVNC(ICX),       ! flow-tangency constraint
     &  KSYSGMG(ICX),       ! gamma-continuity constraint
     &  KSYSGAM(IPX),       ! explicit gamma specification
     &  KSYSKUT(NEX),       ! Kutta condition
     &  KSYSQNC(NEX),       ! element normal-velocity DOF constraint
     &  KSYSGSS(NEX),       ! element gamma_ss constraint
     &  KSYSQTT(NEX),       ! element internal q_tan constraint
     &  KSYSDNS(4,0:NSEGX), ! node-movement constraint at inverse endpoints
     &  JSYSGAM(IPX),       ! node gamma at wall nodes
     &  JSYSDXY(IPX),       ! node normal movement along ANV(.v
     &  JSYSQNC(NEX),       ! element normal-velocity DOF
     &  JSYSQSP(4,0:NSEGX), ! Qspec DOF variable
     &  JPSYS(NSYSX),       ! node corresponding to solution vector index
     &  JAICGAM(IPX),       ! node gamma (for specified gamma)
     &  JAICSIG(IPX),       ! node sigma
     &  JAICXYP(IPX),       ! node x,y
     &  JAICGTH(IPX),       ! node gamma at wake nodes
     &  JAICVWK(NEX)        ! rotor vortex wake (const Wthet/Wm*B*CIRC)
C
C
C===========================================================================
C--- Solver data
      LOGICAL LCONV
C
      COMMON /CL_SOLV/ 
     &  LCONV        ! flag for converged solution 
      COMMON /CI_SOLV/ 
     &  ITRMAXSOLV   ! maximum number of iterations for solver
      COMMON /CR_SOLV/ 
     &  VAVGINIT,    ! average axial velocity (for initialization of VMAVG)
     &  RLXSOLV,     ! relaxation factor 
     &  EPSSOLV      ! convergence criteria 
C
C===========================================================================
C--- BL data
C    NBLVX  number of BL variables in BL data array
C           currently X,Y,S,Q,CP,DIV,DS,TH,H,CF,N,FISID
      PARAMETER (NBLVX=12)
C
      LOGICAL LVISC
C
      COMMON /CL_BL/ 
     &  LVISC        ! flag for BL solution
      COMMON /CI_BL/ 
     &  ICBLBEG(2,NEX), ICBLEND(2,NEX)  ! indices for BL side IC start/end 
      COMMON /CR_BL/ 
     &  AMPMAX, UREN,  
     &  SFIX(2,NEX), STRN(2,NEX), SSEP(2,NEX),
     &  CXVIS(0:NEX),          ! CX from viscous forces
     &  BLDATA(NBLVX,ICX)      ! BL data at all panel centers
C
C===========================================================================
C--- Flow condition section
      COMMON /CR_FLO/ 
     &  RHO,    RMU,    VSO,   ALTH,   DELTAT,   GEE 
C
C  RHO      density
C  RMU      viscosity
C  VSO      speed of sound
C  ALTH     altitude (Std Atmosphere model)
C  DELTAT   delta T (from Std Atmosphere model)
C  GEE      gravity acceleration
C
C===========================================================================
C--- Drag area data
C    NDRGX  number of drag areas
C    IDX  number of points defining drag area element
      PARAMETER (NDRGX=10, IDX=30)
C
      LOGICAL LDRGOBJ
C
      COMMON /CL_DRG/ 
     &  LDRGOBJ        ! flag to use drag area(s) 
C
      COMMON /CR_DRG/ 
     &  XDDEF(IDX,NDRGX),  YDDEF(IDX,NDRGX),   CDADEF(IDX,NDRGX)
C
      COMMON /CI_DRG/
     &  NDOBJ,   
     &  NDDEF(NDRGX),
     &  IELDRGOBJ(NDRGX)
C
C  LDRGOBJ   flag for use of drag objects
C
C  XDDEF(.)  axial  coordinate for drag area definition
C  YDDEF(.)  radial coordinate for drag area definition
C  CDADEF(.) drag area at drag area points
C
C  NDOBJ        number of drag objects
C  NDDEF(.)     number of points in drag area definition
C  IELDRGOBJ(.) index of source line element associated with each drag area
C             or 0 if element not associated with a drag area
C
C===========================================================================
C---  Data for vortex wakes and grid setup
      LOGICAL LWRLX
C
      COMMON /CL_WAK/ 
     &  LWRLX    ! flag for automatic vortex wake relaxation

      COMMON /CR_WAK/ 
     &  XWAKE,     XDWKLEN
C
      COMMON /CI_WAK/
     &  NWAKE
C
C  XWAKE     axial location for downstream end of wakes
C  XDWKLEN   dXwake/diam relative length of wake downstream of TE 
C
C  NWAKE     # of points in wake from TE to downstream end
C
C===========================================================================
C    NRX  number of rotors 
C    IRX  number of rotor radial stations
      PARAMETER (NRX=3, IRX=30)
C
C--- Rotor data
      LOGICAL LROTOR, LBLDEF, LCHANGE,
     &        LSTALLR
C
      COMMON /CL_ROT/ 
     &  LBLDEF,  ! flag for rotor blade defined (chord, beta)
     &  LCHANGE, ! flag for geometric or aero changes in operating point
     &  LSTALLR(IRX,NRX)  ! flag for stall in blade section
C
      COMMON /CR_ROT/ 
C--- Total values for thrust, torque, power
     &  TDUCT,  FOM,
     &  TTOT,   TVIS,
     &  QTOT,   QVIS,
     &  PTOT,   PVIS,
C--- Thrust, torque, power on each rotor/stator disk
     &  TINVR(NRX),  TVISR(NRX),  TTOTR(NRX),
     &  QINVR(NRX),  QVISR(NRX),  QTOTR(NRX),
     &  PINVR(NRX),  PVISR(NRX),  PTOTR(NRX),
C--- Rotor/stator thrust, torque sensitivities
     &  TI_OMG(NRX), QI_OMG(NRX),
     &  TI_QNF(NRX), QI_QNF(NRX),
     &  TV_OMG(NRX), QV_OMG(NRX),
     &  TV_QNF(NRX), QV_QNF(NRX),
     &  TV_DBE(NRX), QV_DBE(NRX),
C
     &  TI_GAM(IRX,NRX), TI_VA(IRX,NRX), TI_VT(IRX,NRX),
     &  QI_GAM(IRX,NRX), QI_VA(IRX,NRX), QI_VT(IRX,NRX),
     &  TV_GAM(IRX,NRX), TV_VA(IRX,NRX), TV_VT(IRX,NRX),
     &  QV_GAM(IRX,NRX), QV_VA(IRX,NRX), QV_VT(IRX,NRX),
C
     &  DTII(IRX,NRX), DQII(IRX,NRX), DPII(IRX,NRX), 
     &  DTVI(IRX,NRX), DQVI(IRX,NRX), DPVI(IRX,NRX), 
C
C--- Disk flow properties
     &  OMEGA(NRX),  
     &  VAAVG(NRX),     ! area averaged axial velocity on disk 
     &  DPSI(IRX,NRX),  ! static pressure rise at radial stations
C--- Velocities are stored as arrays VXXX(3,i) where
C      Vxxx(1,i) is axial velocity
C      Vxxx(2,i) is radial velocity
C      Vxxx(3,i) is circumferential velocity, positive in RH sense about X
C--- Absolute frame induced velocities (just downstream of rotor centers)
     &  VIND(3,IRX,NRX),
C--- Absolute frame velocity and sensitivities downstream of rotor centers
     &  VABS(3,IRX,NRX),   
C--- Blade relative frame total velocities downstream of rotor centers
     &  VREL(3,IRX,NRX),     
C
C--- Rotor geometry
     &  XDISK(NRX),    XPAXIS, 
     &  TGAP,          TGAPZL,
     &  RHUB(NRX),     RTIP(NRX),     
     &  ADISK(NRX),    ATIP(NRX), 
     &  XRP(IRX,NRX),  YRP(IRX,NRX), 
     &  XRC(IRX,NRX),  YRC(IRX,NRX),
C
C--- Rotor blade data
     &  BGAM(IRX,NRX),
     &  CHR(IRX,NRX),   BETAR(IRX,NRX), ALFAR(IRX,NRX),
     &  CLR(IRX,NRX),   CDR(IRX,NRX),    CMR(IRX,NRX),     
     &  RER(IRX,NRX),   MACHR(IRX,NRX),  CLALF(IRX,NRX),
C--- Rotor blade design target data
     &  CHDES(IRX), BETADES(IRX), CLDES(IRX),
C--- Rotor blade and actuator disk definition data 
     &  YRDEF(IRX,NRX),   BGAMDEF(IRX,NRX),
     &  CHRDEF(IRX,NRX),  BETADEF(IRX,NRX),
C--- Rotor geometry at panel boundaries (Esotec)
     &  CHRP(IRX,NRX),    BETARP (IRX,NRX)
C
      COMMON /CI_ROT/
     &  NROTOR,          NRUPSTRM,
     &  NRSTA,           NRP,           NRC,
     &  NRDEF(NRX),
     &  IRTYPDEF(NRX),              ! IRTYPE for defined rotor/disk
     &  NRBLD(NRX),         
     &  IRTYPE(NRX),                ! IRTYPE 0 undefined, 1 act.disk, 2 rotor
     &  IELROTOR(NRX),
     &  IPROTCB(NRX),     IPROTDW(NRX), 
     &  IR2IEL(IRX),      IEL2IR(NEX),
     &  IP2IR(IPX),       IAERO(IRX,NRX)
C
C  XDISK    axial coordinate of rotor radial line
C  RHUB     hub radius
C  RTIP     tip radius
C  TGAP     rotor tip gap (physical length)
C  TGAPZL   tip gap for zero loss (physical length)
C           (normally TGAPZL might be ~1% of blade span, default is zero)
C  ADISK    swept disk area of rotor from hub to tip
C  ATIP     disk area using tip radius only
C  OMEGA    rotor angular speed
C
C  XRP(.)   rotor line points axial coordinate array
C  YRP(.)   rotor line points radial coordinate array
C  XRC(.)   rotor line centerpoint axial coordinate array
C  YRC(.)   rotor line centerpoint radial coordinate array
C  BGAM(.)  rotor circulation B*CIRC at radial stations
C  WX(.)    rotor axial velocity at centers
C  WR(.)    rotor radial velocity at centers
C  WM(.)    rotor streamline m velocity at centers
C  WT(.)    rotor tangential velocity at centers (blade relative)
C  W(.)     rotor velocity at centers (blade relative)
C  PHI(.)   rotor flow angle at centers (blade relative)
C
C  IELROTOR element index of source line defined for rotor
C  IPROTCB  index of panel node at rotor on CB wall
C  IPROTDW  index of panel node at rotor on Duct wall
C
C  NRP      # of rotor radial points
C  NRC      # of rotor radial centerpoints (NR-1)
C  NRBLD    # of rotor blades
C
C  NRDEF    # of rotor defining radial stations (YRDEF,CHRDEF,BETADEF,BGAMDEF)
C  NRSTA    # of rotor stations
C
C  IR2IEL(.) index of element for each rotor vortex wake,
C            or 0 if element not a rotor vortex wake
C
C  IP2IR(..) index of rotor station corresponding to each corner point, 
C            or 0 if point is not associated with a rotor wake 
C  IAERO(.) index of rotor radial station for section
C
C===========================================================================
C--- Blade element aero data section
C    NAX - max number of aerodynamic sections defined
C    NDX - number of aerodynamic parameters defined for each section
      PARAMETER (NAX=20, NDX=15)
C
      COMMON /CR_AERO/ 
     &  XIAERO(NAX,NRX),AERODATA(NDX,NAX,NRX),AZERO(IRX,NRX)
C
      COMMON /CI_AERO/
     &  NAERO(NRX)
C
C  XIAERO(..)    r/R for aero section
C  AERODATA(...) aero data for blade element (aero section)
C                AERODATA( 1,x) = A0 (angle of zero lift)
C                AERODATA( 2,x) = CLMAX (Max CL)
C                AERODATA( 3,x) = CLMIN (Min CL)
C                AERODATA( 4,x) = DCLDA (Incompressible 2-D lift curve slope)
C                AERODATA( 5,x) = DCLDA_STALL (2-D lift curve slope at stall)
C                AERODATA( 6,x) = DCL_STALL (CL increment, onset to full stall)
C                AERODATA( 7,x) = CDMIN (Minimum drag coefficient value)
C                AERODATA( 8,x) = CLDMIN (Lift at minimum drag value)
C                AERODATA( 9,x) = DCDCL2 (Parabolic drag param d(Cd)/dCL^2)
C                AERODATA(10,x) = CMCON (Incompressible 2-D pitching moment)
C                AERODATA(11,x) = REREF (reference Reynold's number)
C                AERODATA(12,x) = REXP (Reynold's number exponent Cd~Re^REXP)
C                AERODATA(13,x) = MCRIT (critical Mach #)
C                AERODATA(14,x) = TOC (thickness/chord)
C                AERODATA(15,x) = DCDCL2S (Secondary, annulus drag param d(Cd)/dCL^2)
C
C  NAERO(.)  # of aero sections defined
C
C===========================================================================
C  Grid parameters (used for rotor wake grid)
C
C      IX  number of streamwise points in grid
C      JX  number of streamlines in grid
C
      PARAMETER (IX=200, JX=60)
C
C------------------------------------------------------------------
C     Derived array limits for grid data
C
C--- Vortex wake grid data
      COMMON /CR_GRD/ 
     &  XG(IX,JX),    YG(IX,JX),
     &  QG(IX,JX),    QXG(IX,JX),   QYG(IX,JX),  QTG(IX,JX),  
     &  RG(IX,JX),    PG(IX,JX),    POG(IX,JX),
     &  BGAMG(IX,JX), DSG(IX,JX),   DHG(IX,JX),
     &  XPOS(IX),     YPOS(JX),     XGMIN
      COMMON /CI_GRD/  
     &  II, JJ, IGTECB, IGTEDW, IGROTOR(NRX),
     &  IC2IG(ICX), IP2IG(IPX)
C 
C  XG(..)    X coordinate of rotor wake grid
C  YG(..)    Y coordinate of rotor wake grid
C  QXG(..)   axial  velocity at grid centers (absolute frame)
C  QYG(..)   radial velocity at grid centers (absolute frame)
C  QTG(..)   tangential velocity at grid centers (absolute frame)
C  QG(..)    overall velocity at grid centers (absolute frame)
C  RG(..)    density at grid centers
C  PG(..)    static pressure at grid centers (absolute frame)
C  P0G(..)   total pressure at grid centers (absolute frame)
C  BGAMG(..) B*GAMMA circumferential circulation at grid centers
C  DSG(..)   entropy deltaS at grid centers
C  DHG(..)   enthalpy deltaH at grid centers
C  XPOS(.)   streamwise spacing array for grid 
C  YPOS(.)   streamfunction array for grid 
C  XGMIN     axial location of most upstream grid points
C
C  II        # of streamwise grid points
C  JJ        # of streamlines 
C  IGTECB     streamwise index of CB   TE in grid
C  IGTEDW     streamwise index of duct TE in grid
C  IGROTOR(.) streamwise index of rotor in grid
C  IC2IG(.)   index IG of grid center for panel center (0 for not on grid)
C  IP2IG(.)   index IG of grid line for panel point (0 for point not on grid)
C
      COMMON/CI_VELS/ NINFL
      COMMON/CR_VELS/ RINFL(IRX),
     &                VAINFL(IRX), VTINFL(IRX),
     &                VAINFLR(IRX),VTINFLR(IRX)
C
C
C===========================================================================
C     ESLOFT stuff
C---------------------------------------------------------------------------
C     NAFX       number of parent airfoils
C     NPX        number of points per parent airfoil
C     NLSX       number of lofting stations (includes 2 possible overshoot)
C     MTI        meters to inches
C
C---- xfoil plot parameters
C     IQXX       following xfoil, could be smaller
C     IBXX       same
C
C---- not many loft parameters
      PARAMETER (NAFX= 12)	  ! max number of parent airfoils
      PARAMETER (NPX = 200)       ! points per airfoil
      PARAMETER (NLSX= 36)        ! loft stations - 2 reserved for osho
      PARAMETER (MTI = 39.370079) ! meters to inches
      PARAMETER (DTRX=PI/180.)    ! degrees to radians
      PARAMETER (NTCX= 10)        ! max t/c locations per airfoil
C
C---- xfoil plot parameters
      PARAMETER (IQXX= 360)
C
C---- for work arrays aka xrotor
      PARAMETER (IWX = 200)
      COMMON/TEMP_X/ WE1(IWX),WE2(IWX),WE3(IWX)
C
C-------------------------------------------------------------------------
C---- esloft declarations
C
      DIMENSION XXE(NPX,NAFX),YYE(NPX,NAFX)
      DIMENSION TOCE(NAFX,NRX),TLOC(NAFX,NRX),ALFZ(NAFX,NRX),
     &          CAMLO(NAFX,NRX),CAMLOC(NAFX,NRX),PAREA(NAFX,NRX),
     &          PANGTE(NAFX,NRX),PLERAD(NAFX,NRX),TETH(NAFX,NRX)
      CHARACTER*80 ENAME(NAFX,NRX),SNAME,LONAME
      DIMENSION YLOFT(NLSX),CDLOFT(NLSX),BELOFT(NLSX),TCLOFT(NLSX),
     &          THLOFT(NLSX),AZDFDC(NLSX),AZLOFT(NLSX),BECORR(NLSX),
     &          BLXX(NPX,NLSX),BLYY(NPX,NLSX),
     &          TRXX(NPX,NLSX),TRYY(NPX,NLSX)
      LOGICAL   LL2D,LCALC,LBLEN,LTRAN,LLOFT,LROTATE,LTDEF,LSMOD
      DIMENSION TT1(NLSX),TT1S(NLSX),
     &          TT2(NLSX),TT2S(NLSX),
     &          TT3(NLSX), NTC(NAFX),TCLOC(NAFX,NTCX)
C
C---- multi-disk extensions
C
      DIMENSION XXEQ(NPX,NAFX,NRX),YYEQ(NPX,NAFX,NRX)
      CHARACTER*80 SNAMEQ(NRX)
      DIMENSION NAFQ(NRX),NPP2Q(NRX),PAFQ(NRX),
     &          TKHUBQ(NRX),TKTIPQ(NRX),AXHUBQ(NRX),AXTIPQ(NRX),
     &          AZTIPQ(NRX),NLOFTQ(NRX),NLOFT2Q(NRX),PLOFTQ(NRX),
     &          OSHUBQ(NRX),OSTIPQ(NRX),ITTYPEQ(NRX),PAXISQ(NRX),
     &          TCLOFTQ(NLSX,NRX),THLOFTQ(NLSX,NRX)
      LOGICAL LTDEFQ(NRX)
C
C---- xfoil plot declarations
C---- generally modified symbols to stay clear of DFDC
C
      DIMENSION XBB(NPX),YBB(NPX),SBB(NPX),XBBP(NPX),YBBP(NPX),
     &          WX1(6*IQXX),WX2(6*IQXX),WX3(6*IQXX),IXPLOT(NLSX)
      LOGICAL   LGGRIDX, LGTICKX, LGPARMX, LGEOPLX
      CHARACTER*80 NAMEX
      DIMENSION BLENDATA(8,NLSX)
C
C---- esloft
C
      COMMON /DF_ESLOFT/ NAF,XXE,YYE,TOCE,TLOC,ALFZ,ENAME,SNAME,
     &                   TETH,NPP,NPP2,PAF,TKHUB,TKTIP,AXHUB,AXTIP,
     &	                 OSHUB,OSTIP,NLOFT,NLOFT2,PLOFT,YLOFT,
     &                   CDLOFT,BELOFT,TCLOFT,THLOFT,BECORR,
     &                   IHUB,ITIP,AZDFDC,AZLOFT,ITTYPE,OUTFAC,
     &                   LL2D,LCALC,LBLEN,LTRAN,LONAME,PAXIS,APARA,
     &                   LLOFT,BLXX,BLYY,TRXX,TRYY,CAMLO,CAMLOC,
     &                   LROTATE,AZTIP,NLOLD,LTDEF,LSMOD,
     &                   TT1,TT1S,TT2,TT2S,TT3,
     &                   NTC,TCLOC,PAREA,PLERAD,PANGTE,
     &                   XXEQ,YYEQ,SNAMEQ,NAFQ,NPP2Q,PAFQ,
     &                   TKHUBQ,TKTIPQ,AXHUBQ,AXTIPQ,
     &                   AZTIPQ,NLOFTQ,NLOFT2Q,PLOFTQ,
     &                   OSHUBQ,OSTIPQ,ITTYPEQ,PAXISQ,
     &                   TCLOFTQ,THLOFTQ,LTDEFQ,ICBS
C
C---- xfoil plots - mostly modified symbols
C
      COMMON/XF_PLOTS/XBB,YBB,SBB,NBB,XBBP,YBBP,
     &                NAMEX,SBLEX,CHORDBX,AREABX,
     &                RADBLEX,ANGBTEX,THICKBX,CAMBRBX,
     &                EI11BA,EI22BA,APX1BA,APX2BA,
     &                EI11BT,EI22BT,APX1BT,APX2BT,       
     &                XSFX,YSFX,XOFFX,YOFFX, 
     &                XGMINX,XGMAXX,YGMINX,YGMAXX,DXYGX,
     &                XBMINX,XBMAXX,YBMINX,YBMAXX,
     &                LGGRIDX,LGTICKX,LGPARMX,LGEOPLX,
     &                IXTYPE,NXPLOT,IXPLOT,CHGX,GTICKX,NOVER,
     &                BLENDATA,XTHICKBX,XCAMBRBX,TETHX
C
C--------------------------------------------------------------------------
C
C     XXE	X coordinate of paneled parent airfoils
C     YYE       Y coordinate of paneled parent airfoils
C     TOCE      Max thickness/chord of        "
C     TLOC      Max thickness x/c of          "
C     CAMLO     Max camber of                 "
C     CAMLOC    Max camber x/c of             "
C     TETH      Trailing edge thickness/chord "
C     ALFZ      Zero lift alpha input for     "
C     PAREA     Section area of               "
C     PLERAD    Leading edge radius of        "
C     PANGTE    Trailing edge angle of        "
C     ENAME     Individual airfoil names (imported)
C     SNAME     Esloft airfoil set name
C     NPP       No. of points per paneled airfoil
C     NPP2      No. of points per side
C     PAF       Airfoil panel spacing parameter (del_TE/del_LE)
C
C     TKHUB     Specified blade t at root (mm). 0 to set airfoil.
C     TKTIP     Specified blade t at tip (mm). 0 to set airfoil.
C     AXHUB     Specified pitch axis at root x/c
C     AXTIP     Specified pitch axis at tip x/c
C     AZTIP     Specified dihedral at tip z (mm)
C     OSHUB     Lofted overshoot at root (beyond design dimensions)
C     OSTIP     Lofted overshoot at tip
C     NLOFT     Specified number of lofting stations
C     NLOFT2    Total lofting stations with overshoot stations
C     PLOFT     Lofting station spacing parameter (del_tip/del_root)
C     YLOFT     Radial location of lofting stations (m)
C     CDLOFT    Chords at lofting stations (m)
C     BELOFT    Beta at lofting stations (from DFDC data) (rad)
C     TCLOFT    Thickness/chord at lofting stations
C     THLOFT    Thickness at lofting stations (m)
C     IHUB      Hub lofting station index
C     ITIP      Tip lofting station index (these accommodate overshoot)
C     AZDFDC    A0 at lofting stations - interp. from AERO data (deg)
C     AZLOFT    A0 at lofting stations - interp. from ESLOFT data (deg)
C     BECORR    Beta at lofting stations corrected using A0 data (rad)
C     ITTYPE    Thickness distribution type:
C               1  linear    t/c
C               2  parabolic t/c
C               3  splined   t/c
C               4  linear     t
C               5  parabolic  t
C               6  splined    t
C     PAXIS     Parabola axis location (units of Rtip-Rhub inbd of hub)
C     APARA     Calculated parabolic coefficient A ( t ~ SQRT(Rad/A) )
C     OUTFAC    Output factor for unit conversions
C     LL2D      True to export 2D points files (false is 3D)
C     LCALC     True if current loft input has been calculated
C     LBLEN	True if current blended sections are calculated
C     LTRAN	True if current transformed sections are calculated
C     LLOFT     When true, GENGEOM ignores TGAPZL when paneling blades
C     LROTATE   True if left hand rotation selected for output
C     BLXX      X coordinates of blended sections
C     BLYY      Y coordinates of blended sections
C     TRXX      X coordinates of transformed sections
C     TRYY      Y coordinates of transformed sections
C
C     NLOLD     Number of loft stns for respline
C     LTDEF     True if thickness data defined for spline
C     LSMOD     True if user input required for spline
C     NTC       Number of occurrences of parent t/c
C     TCLOC     Location of parents on blade
C
C     IXTYPE    Type of xfoil plot
C               1  Parent airfoils
C	        2  Blended airfoils
C               3  Transformed airfoils
C     NXPLOT    Number of airfoils to be plotted
C     IXPLOT    Array of airfoil indices to be plotted
C     LGEOPLX   true if there is a current xfoil plot
C     NOVER     number of plot overlays
C
C==========================================================================
C---- MODI declarations
C
      DIMENSION BUFCHD(IRX,NRX), BUFBET(IRX,NRX),
     &          CURCHD(IRX,NRX), CURBET(IRX,NRX),
     &          WB1(IRX),WB2(IRX),WB3(IRX),
     &                   WC2(IRX)
C
      LOGICAL LBA,LBC
C
      COMMON/DF_MODI/BUFCHD,BUFBET,CURCHD,CURBET,
     &               WB1,WB2,WB3,WC2,LBA,LBC,BLFAC
C
C==========================================================================
C---- Blade blockage stuff
C
      DIMENSION BBVFAC(IRX,NRX), CLPOS(IRX),CLNEG(IRX)
      LOGICAL LBLBL, LBBLOFT(NRX)
C
      COMMON/DF_BBLOC/BBVFAC,BBPFAC,CLPOS,CLNEG,LBLBL,LBBLOFT,BPLAST
C
C     BBVFAC   Blade blockage velocity factor
C     BBPFAC   Blade blockage profile factor
C     LBLBL    True if blade blockage corrections are to be applied
C     LBBLOFT  True if blade blockage factors are stored
C     CLPOS    CL array for positive BGAM disks
C     CLNEG    CL array for negative BGAM disks
C     BPLAST   Profile factor of last stored BB data








