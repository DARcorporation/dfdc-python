!*==I_DFDC.f90  processed by SPAG 7.25DB at 08:52 on  3 Feb 2020



module i_dfdc
    implicit none
    !
    ! PARAMETER definitions
    !
    real, parameter :: pi = 3.1415926535897932384, pi2i = 0.5 / pi, &
            & dtr = pi / 180.
    integer, parameter :: ibx = 1000, nex = 30, nrpnx = 4, &
            & ipx = 2000, icx = 2000, nqspx = 5, &
            & nsegx = 10, nmodx = 10, npntx = 200, &
            & npolx = 10, ncx = icx - nex, &
            & nsx = 2 * nex, nsysx = icx + nex, &
            & nux = ipx + nex + 2, nglvx = 4 * nsegx, &
            & nglpx = nmodx, nrhsx = nglvx + nglpx + &
            & 1, ipal = 1, ipcl = 2, ipcd = 3, &
            & ipcm = 4, ipcv = 5, ipma = 7, &
            & ipqi = 10, iptot = 10, jpcx = 1, &
            & jpcy = 2, jpcm = 3, jpgam = 4, &
            & jpsig = 5, jpsour = 6, jpcirc = 7, &
            & jpdblx = 8, jpdbly = 9, jptot = 9, &
            & kptp = 1, kptn = 2, kptot = 2, &
            & jcplt = 1, jqplt = 2, japlt = 3, &
            & nblvx = 12, ndrgx = 10, idx = 30, &
            & nrx = 3, irx = 30, nax = 20, &
            & ndx = 15, ix = 200, jx = 60, &
            & nafx = 12, npx = 200, nlsx = 36
    real, parameter :: mti = 39.370079, dtrx = pi / 180.
    integer, parameter :: ntcx = 10, iqxx = 360, iwx = 200
    !
    ! COMMON /CC_GEN/
    !
    character(80) :: aname, argp1, argp2, cpfile, frfile, &
            & ispars, name, pfile
    character(2), dimension(99) :: pnum
    character(64) :: prefix
    !
    ! COMMON /CI_AERO/
    !
    integer, dimension(nrx) :: naero
    !
    ! COMMON /CI_BL/
    !
    integer, dimension(2, nex) :: icblbeg, icblend
    !
    ! COMMON /CI_CVP/
    !
    integer, dimension(nex) :: icfrst, iclast, icstg, ipfrst, &
            & iplast, ipte1, ipte2
    integer, dimension(icx) :: ictype, ipco, ipcp
    integer :: nctot, nptot
    !
    ! COMMON /CI_DRG/
    !
    integer, dimension(ndrgx) :: ieldrgobj, nddef
    integer :: ndobj
    !
    ! COMMON /CI_FLO/
    !
    integer, dimension(ipx) :: iugam, iusig
    integer :: iuqinf, npseq, nu
    integer, dimension(nex) :: iuvwk
    !
    ! COMMON /CI_GEN/
    !
    integer, dimension(nex) :: iewake
    integer :: imatyp, itmax, lundbg, nname, nprefix
    integer, dimension(ipx) :: izero
    !
    ! COMMON /CI_GEO/
    !
    integer, dimension(nex) :: ibfrst, iblast, nbtype, netype
    integer :: ielgdes, nbel, nbtot, nel, nmod
    !
    ! COMMON /CI_GRD/
    !
    integer, dimension(icx) :: ic2ig
    integer, dimension(nrx) :: igrotor
    integer :: igtecb, igtedw, ii, jj
    integer, dimension(ipx) :: ip2ig
    !
    ! COMMON /CI_PAN/
    !
    integer, dimension(nex) :: npan
    integer :: npandef
    !
    ! COMMON /CI_PLT/
    !
    integer, dimension(nex) :: icolel, isplot
    integer :: iptype
    !
    ! COMMON /CI_QSP/
    !
    integer :: ielqdes, kqtarg, nqsp, nseg
    integer, dimension(nsegx) :: ielseg, ipseg1, ipseg2
    integer, dimension(ipx) :: ksegp
    !
    ! COMMON /CI_ROT/
    !
    integer, dimension(irx, nrx) :: iaero
    integer, dimension(nex) :: iel2ir
    integer, dimension(nrx) :: ielrotor, iprotcb, iprotdw, &
            & irtypdef, irtype, nrbld, nrdef
    integer, dimension(ipx) :: ip2ir
    integer, dimension(irx) :: ir2iel
    integer :: nrc, nrotor, nrp, nrsta, nrupstrm
    !
    ! COMMON /CI_SOLV/
    !
    integer :: itrmaxsolv
    !
    ! COMMON /CI_SYS/
    !
    integer, dimension(icx) :: ieldofc, ksysgmg, ksysvnc
    integer, dimension(ipx) :: jaicgam, jaicgth, jaicsig, &
            & jaicxyp, jsysdxy, jsysgam, &
            & ksysgam
    integer, dimension(nex) :: jaicvwk, jsysqnc, ksysgss, &
            & ksyskut, ksysqnc, ksysqtt
    integer, dimension(nsysx) :: jpsys
    integer, dimension(4, 0:nsegx) :: jsysqsp, ksysdns
    integer :: naicgam, naicgth, naicsig, naicvwk, naicxyp, &
            & nsys
    !
    ! COMMON /CI_VELS/
    !
    integer :: ninfl
    !
    ! COMMON /CI_WAK/
    !
    integer :: nwake
    !
    ! COMMON /CL_BL/
    !
    logical :: lvisc
    !
    ! COMMON /CL_DRG/
    !
    logical :: ldrgobj
    !
    ! COMMON /CL_GEN/
    !
    logical :: lagrid, lcoord, ldbg, ldint, leplt, lfref, &
            & lgama, lgamu, lggrid, lgparm, lgsys, lgtick, &
            & lload, lncvp, lnorm, lpacc, lpgrid, lpplt, &
            & lpref, lpsho, lpstag, lqaic, lqcnt, lqgic, &
            & lqsys, lrspc, lrspcdef, lrspcusr, lsigm, lsigp, &
            & lsplt, lsysp, lvaic, lvgic, lvmav, lvplt, &
            & lwsho, lxymov
    logical, dimension(nex) :: lbody, lndof, lqnzr, lrevel, &
            & ltpan, lv1zr, lv2zr, lxbod
    logical, dimension(0:nux) :: luset
    !
    ! COMMON /CL_QSP/
    !
    logical, dimension(nsegx) :: ldfix1, ldfix2, lnfix1, lnfix2
    logical :: lgnfix, lgplot, lgsame, lgslop, lgsppl, lgsymm, &
            & lhompl, lqcurv, lqqrev, lqrefl, lqslop, lqspec, &
            & lqsppl, lqsrev, lqsvis
    !
    ! COMMON /CL_ROT/
    !
    logical :: lbldef, lchange
    logical, dimension(irx, nrx) :: lstallr
    !
    ! COMMON /CL_SOLV/
    !
    logical :: lconv
    !
    ! COMMON /CL_WAK/
    !
    logical :: lwrlx
    !
    ! COMMON /CR_AERO/
    !
    real, dimension(ndx, nax, nrx) :: aerodata
    real, dimension(irx, nrx) :: azero
    real, dimension(nax, nrx) :: xiaero
    !
    ! COMMON /CR_BL/
    !
    real :: ampmax, uren
    real, dimension(nblvx, icx) :: bldata
    real, dimension(0:nex) :: cxvis
    real, dimension(2, nex) :: sfix, ssep, strn
    !
    ! COMMON /CR_CVP/
    !
    real, dimension(2, icx) :: anc, cpl_qcl, cpr_qcr, dsc_dxy, &
            & qc, qcl, qcr, qc_xc, qc_yc
    real, dimension(2, 2, icx) :: anc_dxy
    real, dimension(2, ipx) :: anp
    real, dimension(2, 2, ipx) :: anp_xym, anp_xyo, anp_xyp
    real, dimension(icx) :: cpl, cpl_qinf, cpr, cpr_qinf, dsc, &
            & xc, yc
    real, dimension(icx, nux) :: cplu, cpru
    real, dimension(2, icx, nux) :: qclu, qcru
    real, dimension(2, icx, 0:nux) :: qcu
    real, dimension(2, ipx, icx) :: qc_gam, qc_gth, qc_sig, &
            & qc_xp, qc_yp
    real, dimension(ipx) :: sp, vmavg, xp, xps, yp, yps
    real, dimension(nex) :: sple, xpcent, xple, xpmaxe, &
            & xpmine, xprefe, xpte, xstg, &
            & ypcent, yple, ypmaxe, ypmine, &
            & yprefe, ypte, ystg
    real :: xpmax, xpmin, ypmax, ypmin
    !
    ! COMMON /CR_DRG/
    !
    real, dimension(idx, ndrgx) :: cdadef, xddef, yddef
    !
    ! COMMON /CR_FLO/
    !
    real :: alth, deltat, gee, mach, mach1, qinf, qinfp1, &
            & qinfp2, qinfpd, qref, rho, rmu, vso
    real, dimension(0:nex) :: cd, cm, cx, cy
    real, dimension(0:nex, nux) :: cdu, cmu, cxu, cyu
    real, dimension(ipx) :: gam, gamvsp, gth, sig, sigvsp
    real, dimension(nex) :: gamp1, gamp2, gampd, gamset, &
            & qndof, sigp1, sigp2, sigpd, &
            & sigset, xtr1, xtr2
    real, dimension(ipx, 0:nux) :: gamu, gthu, sigu
    real, dimension(nex, 0:nux) :: qndofu
    !
    ! COMMON /CR_GEN/
    !
    real, dimension(ibx) :: one, zer
    character(64) :: sversion
    real :: version
    !
    ! COMMON /CR_GEO/
    !
    real, dimension(nex) :: agbsum, area2da, area2dt, asurfv, &
            & asurfvt, dxbsum, dybsum, eixx2da, &
            & eixx2dt, eixy2da, eixy2dt, eiyy2da, &
            & eiyy2dt, rgyrxv, rgyrxvt, rgyryv, &
            & rgyryvt, sble, volumv, volumvt, &
            & xbcen2da, xbcen2dt, xbcenv, &
            & xbcenvt, xble, xbmaxe, xbmine, &
            & xbrefe, xbte, xfbsum, ybcen2da, &
            & ybcen2dt, ybcenv, ybcenvt, yble, &
            & ybmaxe, ybmine, ybrefe, ybte, &
            & yfbsum
    real :: dtebod, dtepan, xbmax, xbmin, ybmax, ybmin
    real, dimension(ibx) :: sb, xb, xbs, yb, ybs
    real, dimension(2) :: xgbox, xwbox, ygbox, ywbox
    !
    ! COMMON /CR_GRD/
    !
    real, dimension(ix, jx) :: bgamg, dhg, dsg, pg, pog, qg, &
            & qtg, qxg, qyg, rg, xg, yg
    real :: xgmin
    real, dimension(ix) :: xpos
    real, dimension(jx) :: ypos
    !
    ! COMMON /CR_PAN/
    !
    real, dimension(nrpnx, nex) :: crrat, srpn1, srpn2
    real, dimension(nex) :: cvex, fsle, fste, smof
    real :: fpandef
    !
    ! COMMON /CR_PLT/
    !
    real :: cpldel, cplfac, cplmax, cplmin, cpxar, cpyar, &
            & faca, gsf, pvfac, qpldel, qplfac, qplmax, qplmin, &
            & qvfac, xboff, xbsf, xofa, xofg, xpldel, xplfac, &
            & xplmax, xplmin, xpoff, xpsf, xqoff, xqsf, yboff, &
            & ybsf, yofa, yofg, ypldel, yplfac, yplmax, yplmin, &
            & ypoff, ypsf, yqoff, yqsf
    real, dimension(0:7) :: shf, ssizel
    real, dimension(nex) :: xelnum, yelnum
    !
    ! COMMON /CR_QSP/
    !
    real :: algam, clgam, cmgam, cxgam, cygam, qigam, rlxmax
    real, dimension(nqspx) :: alqsp, clqsp, cmqsp, cxqsp, &
            & cyqsp, qiqsp
    real, dimension(4, ipx) :: fspec
    real, dimension(ipx) :: gamsp, qsgam, sigsp, sspec, xspoc, &
            & xspocs, yspoc, yspocs
    real, dimension(nex) :: qspdel, qspmax, qspmin, sspdel, &
            & ssple, sspmax, sspmin
    real, dimension(4, nsegx) :: qspdof
    real, dimension(ipx, nqspx) :: qspec, qspecs
    !
    ! COMMON /CR_ROT/
    !
    real, dimension(nrx) :: adisk, atip, omega, pinvr, ptotr, &
            & pvisr, qinvr, qi_omg, qi_qnf, &
            & qtotr, qvisr, qv_dbe, qv_omg, &
            & qv_qnf, rhub, rtip, tinvr, ti_omg, &
            & ti_qnf, ttotr, tvisr, tv_dbe, &
            & tv_omg, tv_qnf, vaavg, xdisk
    real, dimension(irx, nrx) :: alfar, betadef, betar, betarp, &
            & bgam, bgamdef, cdr, chr, &
            & chrdef, chrp, clalf, clr, cmr, &
            & dpii, dpsi, dpvi, dqii, dqvi, &
            & dtii, dtvi, machr, qi_gam, &
            & qi_va, qi_vt, qv_gam, qv_va, &
            & qv_vt, rer, ti_gam, ti_va, &
            & ti_vt, tv_gam, tv_va, tv_vt, &
            & xrc, xrp, yrc, yrdef, yrp
    real, dimension(irx) :: betades, chdes, cldes
    real :: fom, ptot, pvis, qtot, qvis, tduct, tgap, &
            & tgapzl, ttot, tvis, xpaxis
    real, dimension(3, irx, nrx) :: vabs, vind, vrel
    !
    ! COMMON /CR_SOLV/
    !
    real :: epssolv, rlxsolv, vavginit
    !
    ! COMMON /CR_SYS/
    !
    real, dimension(nsysx, 0:ipx) :: aicgam, aicgth, aicsig
    real, dimension(nsysx, 0:1) :: aicqff
    real, dimension(nsysx, 0:nex) :: aicvwk
    real, dimension(nsysx, 2, 0:ipx) :: aicxyp
    real, dimension(0:nsysx) :: res
    real, dimension(nsysx, 0:nsysx) :: sys
    !
    ! COMMON /CR_TEP/
    !
    real, dimension(2, nex) :: gamt, gamt_dx, gamt_dy, sigt, &
            & sigt_dx, sigt_dy, xpt, ypt
    real, dimension(2, 2, nex) :: gamt_gam, gamt_sig, gamt_xp, &
            & gamt_yp, sigt_gam, sigt_sig, &
            & sigt_xp, sigt_yp, xpt_xp, ypt_yp
    !
    ! COMMON /CR_VELS/
    !
    real, dimension(irx) :: rinfl, vainfl, vainflr, vtinfl, &
            & vtinflr
    !
    ! COMMON /CR_WAK/
    !
    real :: xdwklen, xwake
    !
    ! COMMON /DF_BBLOC/
    !
    real :: bbpfac, bplast
    real, dimension(irx, nrx) :: bbvfac
    real, dimension(irx) :: clneg, clpos
    logical, dimension(nrx) :: lbbloft
    logical :: lblbl
    !
    ! COMMON /DF_ESLOFT/
    !
    real, dimension(nafx, nrx) :: alfz, camlo, camloc, pangte, &
            & parea, plerad, teth, tloc, &
            & toce
    real :: apara, axhub, axtip, aztip, oshub, ostip, outfac, &
            & paf, paxis, ploft, tkhub, tktip
    real, dimension(nrx) :: axhubq, axtipq, aztipq, oshubq, &
            & ostipq, pafq, paxisq, ploftq, &
            & tkhubq, tktipq
    real, dimension(nlsx) :: azdfdc, azloft, becorr, beloft, &
            & cdloft, tcloft, thloft, tt1, &
            & tt1s, tt2, tt2s, tt3, yloft
    real, dimension(npx, nlsx) :: blxx, blyy, trxx, tryy
    character(80), dimension(nafx, nrx) :: ename
    integer :: icbs, ihub, itip, ittype, naf, nloft, nloft2, &
            & nlold, npp, npp2
    integer, dimension(nrx) :: ittypeq, nafq, nloft2q, nloftq, &
            & npp2q
    logical :: lblen, lcalc, ll2d, lloft, lrotate, lsmod, &
            & ltdef, ltran
    character(80) :: loname, sname
    logical, dimension(nrx) :: ltdefq
    integer, dimension(nafx) :: ntc
    character(80), dimension(nrx) :: snameq
    real, dimension(nafx, ntcx) :: tcloc
    real, dimension(nlsx, nrx) :: tcloftq, thloftq
    real, dimension(npx, nafx) :: xxe, yye
    real, dimension(npx, nafx, nrx) :: xxeq, yyeq
    !
    ! COMMON /DF_MODI/
    !
    real :: blfac
    real, dimension(irx, nrx) :: bufbet, bufchd, curbet, curchd
    logical :: lba, lbc
    real, dimension(irx) :: wb1, wb2, wb3, wc2
    !
    ! COMMON /TEMP_X/
    !
    real, dimension(iwx) :: we1, we2, we3
    !
    ! COMMON /XF_PLOTS/
    !
    real :: angbtex, apx1ba, apx1bt, apx2ba, apx2bt, areabx, &
            & cambrbx, chgx, chordbx, dxygx, ei11ba, ei11bt, &
            & ei22ba, ei22bt, gtickx, radblex, sblex, tethx, &
            & thickbx, xbmaxx, xbminx, xcambrbx, xgmaxx, xgminx, &
            & xoffx, xsfx, xthickbx, ybmaxx, ybminx, ygmaxx, &
            & ygminx, yoffx, ysfx
    real, dimension(8, nlsx) :: blendata
    integer, dimension(nlsx) :: ixplot
    integer :: ixtype, nbb, nover, nxplot
    logical :: lgeoplx, lggridx, lgparmx, lgtickx
    character(80) :: namex
    real, dimension(npx) :: sbb, xbb, xbbp, ybb, ybbp
end module i_dfdc
