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

      PROGRAM DFDUCT
C-------------------------------------------------------------------
C     Development version 0.70
C
C     This version is based on the new theoretical formulation which
C     sets wake vorticity to satisfy a zero-pressure condition
C
C     This version can analyze and design multiple blade row 
C     fans (rotor, rotor+stator or multiple rotors) 
C====================================================================
C
C     Version ES3a January 2014
C     Blade blockage code
C     ESLOFT upgrades
C     MODI menu
C
C     Version 070-ES2a October 2009
C     DEL fixed in AERO
C     DIH added to ESLOFT
C
C     Version 070-ES2  May 2009
C     LOFT menu item added
C     LOFT initialization added (lines 97,98)
C
C     Version 070-ES1  Feb 2009
C     Changes from 0.70:
C     SVERSION character string version name
C     LLOAD prevents execution of OPER unless a case file is loaded
C
C     Philip Carter, Esotec Developments
C     philip (at) esotec (dot) org
C     
C--------------------------------------------------------------------
C
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      CHARACTER*128 FNAME, LINE
      CHARACTER*24 PROMPT
      LOGICAL LMODG, LDEBUG, FERROR
C
C---- command input strings
      CHARACTER*4 COMAND
      CHARACTER*128 COMARG
C
C---- local arrays for getting keyboard inputs
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
C
 1000 FORMAT(A)
C
C---- print out a GNU license to start off
      CALL WRGNULIC
C
C---- details of this version
C
      WRITE(*,1005)
 1005 FORMAT(
     &    ' ES versions: Esotec Developments, 2009-2015',
     &  /,' philip at esotec.org',
     &  /, 78('='),/)
C
C---- identify this version of DFDC
      VERSION  = 0.70e03           ! This gets saved to file
      SVERSION = '070-ES3.3'      ! A more flexible format
C
      WRITE(*,1010) SVERSION
 1010 FORMAT( '  ---------------------------'
     &       /'    DFDC version ', A20             
     &       /'  ---------------------------')
C
C---- set up everything
      LDEBUG = .FALSE.
      CALL DFINIT(LDEBUG)
      CALL LOFTINIT1
      CALL LOFTINIT2(1,NRX)
C
C---- try to read input dataset from 1st Unix argument filename
      CALL GETARG0(1,ARGP1)
C
C---- try to get 2nd Unix argument string
      CALL GETARG0(2,ARGP2)
C
      FNAME = ARGP1
      IF(FNAME(1:1).NE.' ') CALL DFLOAD(FNAME,FERROR)
C
C---- print top-level menu
      WRITE(*,1500)
C
C---- begin user-interaction loop
C.....................................................................
C
  500 CONTINUE
      PROMPT = ' DFDC^'
C
      CALL ASKC(PROMPT,COMAND,COMARG)
      IF(COMAND.EQ.'    ') GO TO 500
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
C---------------------------------------------------------
      IF(    COMAND.EQ.'HELP'
     &  .OR. COMAND.EQ.'?   ') THEN
       WRITE(*,1500)
 1500 FORMAT(
     &  //'  QUIT   Exit program',
C
     &  //'  LOAD f Read ducted fan case from file',
     &   /'  SAVE f Write current ducted fan case to file',
C     &   /'  GLOA f Read geometry only from file',
C     &   /'  GSAV f Write geometry only to file',
C
     &  //' .OPER   Direct operating point(s)',
     &   /' .AERO   Set or edit blade section aero data',
     &   /' .DRAG   Set or edit drag objects',
     &   /' .MODI   Rotor/stator blade editing facility',
     &   /' .LOFT   Rotor/stator blade lofting facility',
C
     &  //' .GDES   Geometry design facility',
     &   /' .QDES   Mixed-inverse design facility',
     &   /' .PPAR   Show/change paneling',
     &   /'  PANE   Regenerate paneled geometry from buffer geometry',
     &   /'  PCOP   Copy buffer geometry directly into paneling',
     &   /'  PWRT f Write paneling parameters to disk file',
     &   /'  PGET f Get paneling parameters from disk file',
C
     &  //' .PLOP   Plotting options',
     &   /'  NAME s Change case name',
     &   /'  INIT   Reinitialize DFDC')
c     &   /'  DEBU   toggle debugging flag')
C
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'QUIT') THEN
        CALL PLCLOSE
        STOP
C
C---------------------------------------------------------
C--- Load case from file
      ELSEIF(COMAND.EQ.'LOAD') THEN
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') THEN
          CALL ASKS(' Enter filename to load)^',FNAME)
        ENDIF
        CALL DFLOAD(FNAME,FERROR)
        IF(.NOT.FERROR) THEN
          LONAME = NAME
          CALL LOFTINIT2(1,NROTOR)
        ENDIF
C
C---------------------------------------------------------
C--- Save case to file
      ELSEIF(COMAND.EQ.'SAVE') THEN
        IF(COMARG(1:1).NE.' ') PFILE = COMARG
        CALL DFSAVE(COMARG)
C
C---------------------------------------------------------
C--- Load airfoil geometry only
      ELSEIF(COMAND.EQ.'GLOA') THEN
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') THEN
          CALL ASKS(' Enter filename to load)^',FNAME)
        ENDIF
        CALL DFLOAD0(FNAME)
C
C---------------------------------------------------------
C--- Save airfoil geometry only
      ELSEIF(COMAND.EQ.'GSAV') THEN
        CALL DFSAVE0(COMARG)
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'OPER') THEN
        IF(.NOT.LLOAD) THEN
          WRITE(*,*)
          WRITE(*,*)'Load case file first'
        ELSE
          CALL DFOPER
        ENDIF
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'AERO') THEN
        IF(.NOT.LLOAD) THEN
          WRITE(*,*)
          WRITE(*,*)'Load case file first'
        ELSE
          CALL DFAERO
        ENDIF
C
C--------------------------------------------------------------------
      ELSEIF(COMAND.EQ.'DRAG') THEN
C---- Drag area menu
       CALL DRGOPER(LMODG)
       IF(LMODG) CALL GENGEOM
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'GDES') THEN
        CALL DFGDES(LMODG)
        IF(LMODG) CALL GEPROC
c        CALL XYBDMP  
c        CALL GENGEOM
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'QDES') THEN
        CALL DFQDES
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'LOFT') THEN
        CALL ESLOFT 
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'MODI') THEN
        CALL DFMODI 
C
C---------------------------------------------------------
C---- Command to generate panel geometry from buffer (REMOVE LATER)
      ELSEIF(COMAND.EQ.'NGEO') THEN
        CALL GENGEOM
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'PANE') THEN
        CALL PANGEN(.FALSE.)
        CALL GEPROC
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'PPAR') THEN
        CALL PANGEN(.TRUE.)
        CALL GEPROC
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'PCOP') THEN
        CALL PANCOP
        NPAN(1) = IPLAST(1) - IPFRST(1) + 1
        CALL GEPROC
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'PWRT') THEN
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') THEN
          CALL ASKS(' Enter filename for panel parameters)^',FNAME)
        ENDIF
        PFILE = FNAME
        CALL PANWRT
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'PGET') THEN
        FNAME = COMARG
        IF(FNAME(1:1).EQ.' ') THEN
          CALL ASKS(' Enter filename for panel parameters)^',FNAME)
        ENDIF
        PFILE = FNAME
        CALL PANGET(PFILE,ERROR)
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'PLOP') THEN
        CALL OPLSET
        SSIZEL(0) = SHGT*SHF(0)
        SSIZEL(1) = SHGT*SHF(1)
        SSIZEL(2) = SHGT*SHF(2)
        SSIZEL(3) = SHGT*SHF(3)
        SSIZEL(4) = SHGT*SHF(4)
        SSIZEL(5) = SHGT*SHF(5)
        SSIZEL(6) = SHGT*SHF(6)
        SSIZEL(7) = SHGT*SHF(7)
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'NAME') THEN
        WRITE(*,1020) NAME
 1020   FORMAT(/' Current case name is:  ', A)
C
        CALL ASKS(' Enter new case name^',NAME)
        CALL STRIP(NAME,NNAME)
        LONAME=NAME
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'INIT') THEN
        CALL DFINIT(LDEBUG)
        CALL LOFTINIT1
        CALL LOFTINIT2(1,NRX)
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'DEBU') THEN
        LDBG = .NOT.LDBG
        IF(LDBG) THEN
          WRITE(*,*) 'DEBUG output is enabled'
        ELSE
          WRITE(*,*) 'DEBUG output is disabled'
        ENDIF
        LDEBUG = LDBG
C
C---------------------------------------------------------
c      ELSEIF(COMAND.EQ.'DUMP') THEN
C------ dump all panel nodes to fort.11
c        DO IP = 1, NPTOT
c          WRITE(11,*) XP(IP), YP(IP)
c        ENDDO
C
C------ dump all control points to fort.12
c        DO IC = 1, NCTOT
c          WRITE(12,*) XC(IC), YC(IC)
c        ENDDO
C---------------------------------------------------------
      ELSE 
        WRITE(*,1050) COMAND
 1050   FORMAT(1X,A4,' command not recognized.  Type a "?" for list')
C
      ENDIF
C
C---- finished with last command... go back to menu prompt
      GO TO 500
      END




      SUBROUTINE DFLOAD0(FNAME)
C----------------------------------------------------
C     Reads geometry coordinates in XFOIL multi-element
C     format from file FNAME
C     and starts geometry processing for duct analysis
C----------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*)   FNAME
      CHARACTER*128   IFILE, LINE, FILEBL
      LOGICAL LOPEN, LREPANEL
C
C---- local arrays for calling AREAD
      PARAMETER (ITX=300)
      DIMENSION XT(ITX,NEX), YT(ITX,NEX)
      DIMENSION NT(NEX)
C
      LOGICAL ERROR
C
      LU = 1
      CALL STRIP(FNAME,NF)
      IF(FNAME.NE.' ') THEN
        OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=98)
        LOPEN = .TRUE.
       ELSE
        RETURN
      ENDIF
C
      XWAKE = -999.
      XDWKLEN = 1.5
      NWAKE = 20
      XDISK(1) = 999.
C
      NROTOR = 0
      LBLDEF = .FALSE.
      IRTYPE(1) = 0
      IRTYPDEF(1) = 0
C
C---- Read the geometry airfoil file
      FILEBL = ' '
      CALL AREAD(FILEBL,LU, ITX,NEX, XT,YT,
     &           NT,NTEL,
     &           ANAME,ISPARS,IFTYPE)
      IF(IFTYPE.EQ.0 .OR. NTEL.NE.2) THEN
C----- read error occurred
       WRITE(*,*) 'Read error on airfoil geometry read'
       RETURN
      ENDIF
C
      NBEL = NTEL
C
C---- Process airfoil inputs for CB and duct
C---- put element-indexed coordinates XT,YT into linear arrays XB,YB for 
C     buffer airfoil geometry
      IBNEXT = 1
      DO IEL=1, NBEL
C------ first point in element IEL
        IBFRST(IEL) = IBNEXT
C------ initalize accumulation counter
        IB = IBFRST(IEL) - 1
C------ go over all points on element IEL
        DO IT = 1, NT(IEL)
          IB = IB+1
          IF(IB.GT.IBX) STOP 'LOAD: Array overflow on IBX'
          XB(IB) = XT(IT,IEL)
          YB(IB) = YT(IT,IEL)
cc        IF(LDBG) WRITE(20,*) IT,XT(IT,IEL),YT(IT,IEL)
        ENDDO
C---- set buffer airfoils to NBTYPE=0
        NBTYPE(IEL) = 0
        IBLAST(IEL) = IB
        IBNEXT = IB + 1
      END DO
C
C---- set total number of points in buffer geometry
      NBTOT = IBLAST(NBEL)
C
C
C=======================================================================
C---- process buffer foil geometry 
      CALL GEPROC
C
      RETURN
C
C...............................................................
   98 CONTINUE
      WRITE(*,1050) FNAME(1:NF)
      RETURN
C
   99 CONTINUE
      WRITE(*,1100) FNAME(1:NF)
      IF(LOPEN) CLOSE(LU)
      RETURN
C
 1000 FORMAT(A)
 1050 FORMAT(/' File OPEN error:  ', A)
 1100 FORMAT(/' File READ error:  ', A)
C
      END ! DFLOAD0


      SUBROUTINE DFSAVE0(FNAMEIN)
C--------------------------------------
C     Writes out current geometry in XFOIL
C     multi-element airfoil format
C--------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*) FNAMEIN
      LOGICAL OK
C
      CHARACTER*80 FNAME
C
      LU = 2
      IFTYPE = 2
C
C---- save only # elements in buffer geometry
      NELSAV = NBEL
C
C---- get output filename if it was not supplied
      IF(FNAMEIN(1:1) .NE. ' ') THEN
       FNAME = FNAMEIN
      ELSE
       CALL ASKS('Enter output filename^',FNAME)
      ENDIF
C
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=5)
      CALL ASKL('Output file exists.  Overwrite?^',OK)
      IF(OK) GO TO 6
C
      CLOSE(LU)
      WRITE(*,*) 'Current airfoil not saved.'
      RETURN
C
 5    OPEN(LU,FILE=FNAME,STATUS='NEW',ERR=90)
 6    REWIND(LU)
C      
      CALL AWRITE(FNAME,LU,
     &            NELSAV,IPFRST,IPLAST,XP,YP,
     &            ANAME, ISPARS, IFTYPE)
C
      CLOSE(LU)
      RETURN
C
 90   WRITE(*,*) 'Bad filename.'
      WRITE(*,*) 'Current airfoil not saved.'
      RETURN
C
 1000 FORMAT(A)
 1100 FORMAT(1X,2F12.6)
      END ! DFSAVE0







