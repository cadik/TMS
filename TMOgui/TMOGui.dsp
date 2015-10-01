# Microsoft Developer Studio Project File - Name="TMOGui" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=TMOGui - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "TMOGui.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "TMOGui.mak" CFG="TMOGui - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "TMOGui - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "TMOGui - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "$(QTDIR)\include" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "QT_DLL" /D "QT_THREAD_SUPPORT" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x405 /d "NDEBUG"
# ADD RSC /l 0x405 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib imm32.lib wsock32.lib $(QTDIR)\lib\qt-mt230nc.lib $(QTDIR)\lib\qtmain.lib c:\tmo\lib\tmolib.lib c:\tmo\lib\TMOW32.lib /nologo /subsystem:windows /machine:I386 /out:"c:/TMO/bin/TMOGui.exe"

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /GZ /c
# ADD CPP /nologo /MD /W3 /Gm /GX /ZI /Od /I "$(QTDIR)\include" /D "_WINDOWS" /D "QT_DLL" /D "QT_THREAD_SUPPORT" /D "UNICODE" /D "WIN32" /D "_DEBUG" /D "TMOPLUGIN_EXPORTS" /YX /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x405 /d "_DEBUG"
# ADD RSC /l 0x405 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib imm32.lib wsock32.lib $(QTDIR)\lib\qt-mteval333.lib $(QTDIR)\lib\qtmain.lib ../tmolib/Debug/tmolib.lib ../TMOW32/Debug/TMOW32.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "TMOGui - Win32 Release"
# Name "TMOGui - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;for;f90"
# Begin Source File

SOURCE=.\main.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_resources.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_resources1.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_resources2.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_resources3.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_resources4.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_resources5.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIAdjust.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIAdjustValues.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIBitmap.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIFilters.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIHisto.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIImage.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIInfo.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIMenu.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIOutput.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIParameters.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIParametersItem.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIProgressBar.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIRightBar.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIStatistics.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIStatus.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIToneMapping.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_TMOGUIToneSlider.cpp
# End Source File
# Begin Source File

SOURCE=.\moc_tmoguiwindow.cpp
# End Source File
# Begin Source File

SOURCE=.\resources.cpp
# End Source File
# Begin Source File

SOURCE=.\resources1.cpp
# End Source File
# Begin Source File

SOURCE=.\resources2.cpp
# End Source File
# Begin Source File

SOURCE=.\resources3.cpp
# End Source File
# Begin Source File

SOURCE=.\resources4.cpp
# End Source File
# Begin Source File

SOURCE=.\resources5.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIAdjust.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIAdjustValues.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIBitmap.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIFilters.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIHisto.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIImage.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIInfo.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIMenu.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIOutput.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIParameters.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIParametersItem.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIProgressBar.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIRightBar.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIStatistics.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIStatus.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIToneMapping.cpp
# End Source File
# Begin Source File

SOURCE=.\TMOGUIToneSlider.cpp
# End Source File
# Begin Source File

SOURCE=.\tmoguiwindow.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\resources.h
# End Source File
# Begin Source File

SOURCE=.\resources1.h
# End Source File
# Begin Source File

SOURCE=.\resources2.h
# End Source File
# Begin Source File

SOURCE=.\resources3.h
# End Source File
# Begin Source File

SOURCE=.\resources4.h
# End Source File
# Begin Source File

SOURCE=.\resources5.h
# End Source File
# Begin Source File

SOURCE=.\TMOGUIAdjust.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIAdjust.h
InputName=TMOGUIAdjust

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIAdjust.h
InputName=TMOGUIAdjust

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIAdjustValues.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIAdjustValues.h
InputName=TMOGUIAdjustValues

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIAdjustValues.h
InputName=TMOGUIAdjustValues

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIBitmap.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIBitmap.h
InputName=TMOGUIBitmap

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIBitmap.h
InputName=TMOGUIBitmap

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIFilters.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIFilters.h
InputName=TMOGUIFilters

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIFilters.h
InputName=TMOGUIFilters

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIHisto.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIHisto.h
InputName=TMOGUIHisto

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIHisto.h
InputName=TMOGUIHisto

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIImage.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIImage.h
InputName=TMOGUIImage

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIImage.h
InputName=TMOGUIImage

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIInfo.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIInfo.h
InputName=TMOGUIInfo

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIInfo.h
InputName=TMOGUIInfo

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIMenu.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIMenu.h
InputName=TMOGUIMenu

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIMenu.h
InputName=TMOGUIMenu

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIOutput.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIOutput.h
InputName=TMOGUIOutput

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIOutput.h
InputName=TMOGUIOutput

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIParameters.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIParameters.h
InputName=TMOGUIParameters

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIParameters.h
InputName=TMOGUIParameters

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIParametersItem.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIParametersItem.h
InputName=TMOGUIParametersItem

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIParametersItem.h
InputName=TMOGUIParametersItem

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIProgressBar.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIProgressBar.h
InputName=TMOGUIProgressBar

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIProgressBar.h
InputName=TMOGUIProgressBar

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIRightBar.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIRightBar.h
InputName=TMOGUIRightBar

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIRightBar.h
InputName=TMOGUIRightBar

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIStatistics.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIStatistics.h
InputName=TMOGUIStatistics

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIStatistics.h
InputName=TMOGUIStatistics

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIStatus.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIStatus.h
InputName=TMOGUIStatus

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIStatus.h
InputName=TMOGUIStatus

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIToneMapping.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIToneMapping.h
InputName=TMOGUIToneMapping

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIToneMapping.h
InputName=TMOGUIToneMapping

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIToneSlider.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIToneSlider.h
InputName=TMOGUIToneSlider

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\TMOGUIToneSlider.h
InputName=TMOGUIToneSlider

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\tmoguiwindow.h

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\tmoguiwindow.h
InputName=tmoguiwindow

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).h ...
InputDir=.
InputPath=.\tmoguiwindow.h
InputName=tmoguiwindow

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp

# End Custom Build

!ENDIF 

# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;cnt;rtf;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\readme.txt
# End Source File
# Begin Source File

SOURCE=.\resources.ui

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources.ui
InputName=resources

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources.ui
InputName=resources

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\resources1.ui

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources1.ui
InputName=resources1

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources1.ui
InputName=resources1

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\resources2.ui

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources2.ui
InputName=resources2

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources2.ui
InputName=resources2

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\resources3.ui

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources3.ui
InputName=resources3

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources3.ui
InputName=resources3

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\resources4.ui

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources4.ui
InputName=resources4

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources4.ui
InputName=resources4

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\resources5.ui

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources5.ui
InputName=resources5

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Uic'ing $(InputName).ui ...
InputDir=.
InputPath=.\resources5.ui
InputName=resources5

BuildCmds= \
	%qtdir%\bin\uic.exe $(InputPath) -o $(InputDir)\$(InputName).h \
	%qtdir%\bin\uic.exe $(InputPath) -i $(InputName).h -o $(InputDir)\$(InputName).cpp \
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).h -o $(InputDir)\moc_$(InputName).cpp \
	

"$(InputDir)\$(InputName).h" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)

"$(InputDir)\moc_$(InputName).cpp" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
   $(BuildCmds)
# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIHisto.moc

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).cpp ...
InputDir=.
InputPath=.\TMOGUIHisto.moc
InputName=TMOGUIHisto

"$(InputDir)\$(InputName).moc" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).cpp -o $(InputDir)\$(InputName).moc

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).cpp ...
InputDir=.
InputPath=.\TMOGUIHisto.moc
InputName=TMOGUIHisto

"$(InputDir)\$(InputName).moc" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).cpp -o $(InputDir)\$(InputName).moc

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TMOGUIMenu.moc

!IF  "$(CFG)" == "TMOGui - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).cpp ...
InputDir=.
InputPath=.\TMOGUIMenu.moc
InputName=TMOGUIMenu

"$(InputDir)\$(InputName).moc" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).cpp -o $(InputDir)\$(InputName).moc

# End Custom Build

!ELSEIF  "$(CFG)" == "TMOGui - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build - Moc'ing $(InputName).cpp ...
InputDir=.
InputPath=.\TMOGUIMenu.moc
InputName=TMOGUIMenu

"$(InputDir)\$(InputName).moc" : $(SOURCE) "$(INTDIR)$(OUTDIR)"
	%qtdir%\bin\moc.exe $(InputDir)\$(InputName).cpp -o $(InputDir)\$(InputName).moc

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\TODO.txt
# End Source File
# End Target
# End Project
