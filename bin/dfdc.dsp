# Microsoft Developer Studio Project File - Name="dfdc" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) External Target" 0x0106

CFG=dfdc - Win32 Release
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "dfdc.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "dfdc.mak" CFG="dfdc - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "dfdc - Win32 Release" (based on "Win32 (x86) External Target")
!MESSAGE "dfdc - Win32 Debug" (based on "Win32 (x86) External Target")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""

!IF  "$(CFG)" == "dfdc - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir ""
# PROP BASE Intermediate_Dir ""
# PROP BASE Cmd_Line "MAKE /f dfdc.mak"
# PROP BASE Rebuild_Opt "/a"
# PROP BASE Target_File "dfdc.exe"
# PROP BASE Bsc_Name "dfdc.bsc"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ""
# PROP Intermediate_Dir ""
# PROP Cmd_Line "make -f dfdc.mak"
# PROP Rebuild_Opt "/a"
# PROP Target_File ".\dfdc.exe"
# PROP Bsc_Name ""
# PROP Target_Dir ""

!ELSEIF  "$(CFG)" == "dfdc - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Cmd_Line "NMAKE /f dfdc.mak"
# PROP BASE Rebuild_Opt "/a"
# PROP BASE Target_File "dfdc.exe"
# PROP BASE Bsc_Name "dfdc.bsc"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Cmd_Line "make -f   dfdc.mak"
# PROP Rebuild_Opt "/a"
# PROP Target_File "gtest.exe"
# PROP Bsc_Name ""
# PROP Target_Dir ""

!ENDIF 

# Begin Target

# Name "dfdc - Win32 Release"
# Name "dfdc - Win32 Debug"

!IF  "$(CFG)" == "dfdc - Win32 Release"

!ELSEIF  "$(CFG)" == "dfdc - Win32 Debug"

!ENDIF 

# Begin Group "src-driver"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE="..\src-driver\debug.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\dfdc.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\dfdc2.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\esloft.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\gdes.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\gplots.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\gpoint.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\gui.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\list.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\lplots.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\modi.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\modify.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\oper.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\oper1.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\oplots.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\oplset.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\qdes.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\qplots.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\rplots.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\sl.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\testaic.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\wake.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\xfsubs.f"
# End Source File
# Begin Source File

SOURCE="..\src-driver\xmodify.f"
# End Source File
# End Group
# Begin Group "src-lib"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE="..\src-lib\adjpanl.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\aero.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\aic.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\airio.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\atmo.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\blcalc.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\dfdcsubs.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\forces.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\gauss.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\geom.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\geutil.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\grdutils.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\inigrd.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\lamp.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\plutil.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\pnsubs.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\qaic.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\rotoper.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\sgutil.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\sgutil2.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\solve.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\spline.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\system.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\userio.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\userio1.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\userio2.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\vels.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\viscvel.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\wakesubs.f"
# End Source File
# Begin Source File

SOURCE="..\src-lib\xutils.f"
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
