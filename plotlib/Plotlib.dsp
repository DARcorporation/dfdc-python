# Microsoft Developer Studio Project File - Name="plotlib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) External Target" 0x0106

CFG=plotlib - Win32 Release
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Plotlib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Plotlib.mak" CFG="plotlib - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "plotlib - Win32 Release" (based on "Win32 (x86) External Target")
!MESSAGE "plotlib - Win32 Debug" (based on "Win32 (x86) External Target")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""

!IF  "$(CFG)" == "plotlib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Cmd_Line "MAKE /f plotlib.mak"
# PROP BASE Rebuild_Opt "/a"
# PROP BASE Target_File "plotlib.exe"
# PROP BASE Bsc_Name "plotlib.bsc"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Cmd_Line "make -f plotlib.mak"
# PROP Rebuild_Opt "/a"
# PROP Target_File ".\plotlib.exe"
# PROP Bsc_Name ""
# PROP Target_Dir ""

!ELSEIF  "$(CFG)" == "plotlib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Cmd_Line "NMAKE /f plotlib.mak"
# PROP BASE Rebuild_Opt "/a"
# PROP BASE Target_File "plotlib.exe"
# PROP BASE Bsc_Name "plotlib.bsc"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Cmd_Line "make -f   plotlib.mak"
# PROP Rebuild_Opt "/a"
# PROP Target_File "gtest.exe"
# PROP Bsc_Name ""
# PROP Target_Dir ""

!ENDIF 

# Begin Target

# Name "plotlib - Win32 Release"
# Name "plotlib - Win32 Debug"

!IF  "$(CFG)" == "plotlib - Win32 Release"

!ELSEIF  "$(CFG)" == "plotlib - Win32 Debug"

!ENDIF 

# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\gw_subs.f
# End Source File
# Begin Source File

SOURCE=.\plt_3D.f
# End Source File
# Begin Source File

SOURCE=.\plt_base.f
# End Source File
# Begin Source File

SOURCE=.\plt_color.f
# End Source File
# Begin Source File

SOURCE=.\plt_font.f
# End Source File
# Begin Source File

SOURCE=.\plt_old.f
# End Source File
# Begin Source File

SOURCE=.\plt_util.f
# End Source File
# Begin Source File

SOURCE=.\ps_subs.f
# End Source File
# Begin Source File

SOURCE=.\ps_subs_old.f
# End Source File
# Begin Source File

SOURCE=.\set_subs.f
# End Source File
# Begin Source File

SOURCE=".\util-ops.f"
# End Source File
# Begin Source File

SOURCE=.\W32win.c
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
