CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)

SET(CMAKE_AUTOMOC ON)
SET(CMAKE_AUTOUIC ON)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)



IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-DLINUX)
ENDIF(UNIX AND NOT LINUX_SET)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET (QT_MIN_VERSION "4.8.0")
FIND_PACKAGE(Qt4 REQUIRED)

INCLUDE( ${QT_USE_FILE} )
add_definitions(
${QT_DEFINITIONS}
)


INCLUDE_DIRECTORIES(
${CMAKE_BINARY_DIR}
${CMAKE_CURRENT_BINARY_DIR}
${QT_INCLUDES}
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
TMOGUICustomEvents.h
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


QT4_WRAP_CPP(tmogui TMO_SOURCES ${TMO_HEADERS})
#QT4_WRAP_UI(tmogui TMO_HEADERS TMO_SOURCES ${TMO_FORMS})



SET(EXE_NAME "tmogui${TAG_VERSION}" )
ADD_EXECUTABLE(${EXE_NAME}
  ${TMO_SOURCES}
)

target_link_libraries(tmogui Qt4::QtGui Qt4::QtCore ${QT_LIBRARIES})
add_custom_command( TARGET ${EXE_NAME} POST_BUILD
  COMMAND cp -r tmogui ../
)
