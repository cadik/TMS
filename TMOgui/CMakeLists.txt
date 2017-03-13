CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)


IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET (QT_MT_REQUIRED true)
SET (QT_MIN_VERSION "3.3.8")
FIND_PACKAGE(Qt3 REQUIRED)
IF (QT_FOUND)
	#SET(QT_UIC_EXECUTABLE "$ENV{QTDIR}/bin/uic") #did not find correct uic on Fedora

# Qt 3 should prefer QTDIR over the PATH
find_program(MY_QT_UIC_EXECUTABLE
  NAMES uic-qt3 uic uic3 uic3-mt 
  HINTS
    ENV QTDIR
  PATHS
  "[HKEY_CURRENT_USER\\Software\\Trolltech\\Qt3Versions\\3.2.1;InstallDir]/include/Qt"
  "[HKEY_CURRENT_USER\\Software\\Trolltech\\Qt3Versions\\3.2.0;InstallDir]/include/Qt"
  "[HKEY_CURRENT_USER\\Software\\Trolltech\\Qt3Versions\\3.1.0;InstallDir]/include/Qt"
  ${GLOB_PATHS_BIN}
    /usr/local/qt
    /usr/lib/qt
    /usr/lib/qt3
    /usr/share/qt3
    C:/Progra~1/qt
    /usr/X11R6
  PATH_SUFFIXES
    bin
  )

  SET(QT_UIC_EXECUTABLE ${MY_QT_UIC_EXECUTABLE})
	
  MESSAGE("Qt3 Found OK ${QT_VERSION_STR}")
   MESSAGE("$ENV{QTDIR}")
   MESSAGE("${QT_MOC_EXECUTABLE}")
   MESSAGE("${QT_UIC_EXECUTABLE}")
   MESSAGE("${QT_QT_LIBRARY}")
   MESSAGE("${QT_INCLUDE_DIR}")
  MESSAGE("${qt_version_str}")
   MESSAGE(${qt_version_str_lib})
  ADD_DEFINITIONS(${QT_DEFINITIONS})
ELSE(QT_FOUND)
  MESSAGE(FATAL_ERROR "No Qt3")
ENDIF(QT_FOUND)


INCLUDE_DIRECTORIES(
${CMAKE_BINARY_DIR}
${QT_INCLUDE_DIR}
../tmolib/
.
)

#Set up libraries to link with all found packages
LINK_LIBRARIES (
${QT_QT_LIBRARY}
${QT_QASSISTANTCLIENT_LIBRARY}
pthread
tmo
)

SET(TMO_SOURCES
lqstring.cpp 
main.cpp 
TMOFilters.cpp 
TMOGUIAdjust.cpp 
TMOGUIAdjustValues.cpp 
TMOGUIBitmap.cpp 
TMOGUICenterView.cpp 
TMOGUIFileToolBar.cpp 
TMOGUIFilters.cpp 
TMOGUIHisto.cpp 
TMOGUIImage.cpp 
TMOGUIInfo.cpp 
TMOGUIInfoPanel.cpp 
TMOGUIInfoTool.cpp 
TMOGUIInfoToolBar.cpp 
TMOGUILineResizer.cpp 
TMOGUIMenu.cpp 
TMOGUIOutput.cpp 
TMOGUIParameters.cpp 
TMOGUIParametersItem.cpp 
TMOGUIProgressBar.cpp 
TMOGUIRightBar.cpp 
TMOGUISaveDialog.cpp 
TMOGUIStatistics.cpp 
TMOGUIStatus.cpp 
TMOGUIToneMapping.cpp 
TMOGUIToneSlider.cpp 
TMOGUITransformation.cpp 
tmoguiwindow.cpp 
TMOGUIZoomTool.cpp 
)

SET(TMO_HEADERS
gamma.h 
lqstring.h 
resources.ui.h 
resources6.ui.h 
TMOFilters.h 
TMOGUIAdjust.h 
TMOGUIAdjustValues.h 
TMOGUIBitmap.h 
TMOGUICenterView.h 
TMOGUIFileToolBar.h 
TMOGUIFilters.h 
TMOGUIHisto.h 
TMOGUIImage.h 
TMOGUIInfo.h 
TMOGUIInfoPanel.h 
TMOGUIInfoTool.h 
TMOGUIInfoToolBar.h 
TMOGUILineResizer.h 
TMOGUIMenu.h 
TMOGUIOutput.h 
TMOGUIParameters.h 
TMOGUIParametersItem.h 
TMOGUIProgressBar.h 
TMOGUIResource.h 
TMOGUIRightBar.h 
TMOGUISaveDialog.h 
TMOGUIStatistics.h 
TMOGUIStatus.h 
TMOGUIToneMapping.h 
TMOGUIToneSlider.h 
TMOGUITransformation.h 
tmoguiwindow.h 
TMOGUIZoomTool.h
)

SET(TMO_FORMS 
resources.ui
resources1.ui
resources2.ui
resources3.ui
resources4.ui
resources5.ui
resources6.ui
)

qt_wrap_cpp(tmogui TMO_SOURCES ${TMO_HEADERS})
qt_wrap_ui(tmogui TMO_HEADERS TMO_SOURCES ${TMO_FORMS})


SET(EXE_NAME "tmogui${TAG_VERSION}" )
ADD_EXECUTABLE(${EXE_NAME}
  ${TMO_SOURCES}
)

add_custom_command( TARGET ${EXE_NAME} POST_BUILD
	COMMAND cp -r tmogui ../
	)
