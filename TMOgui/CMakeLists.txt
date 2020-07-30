CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)

project(tmogui)

SET(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
SET(CMAKE_AUTOUIC ON)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

cmake_minimum_required(VERSION 3.1.0)

if(CMAKE_VERSION VERSION_LESS "3.7.0")
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

find_package(Qt5 COMPONENTS Widgets PrintSupport REQUIRED)

INCLUDE_DIRECTORIES(
${CMAKE_BINARY_DIR}
${CMAKE_CURRENT_BINARY_DIR}
../tmolib/
.
)

#Set up libraries to link with all found packages
LINK_LIBRARIES (
#${QT_QT_LIBRARY}
#${QT_QASSISTANTCLIENT_LIBRARY}
pthread
tmo
)

SET(TMO_SOURCES
lqstring.cpp
main.cpp
switch.cpp
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
TMOGUIToneMappingChooser.cpp
TMOGUIToneSlider.cpp
TMOGUITransformation.cpp
tmoguiwindow.cpp
TMOGUIZoomTool.cpp
)

SET(TMO_HEADERS
gamma.h
lqstring.h
style.h
switch.h
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
TMOGUIToneMappingChooser.h
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

SET(RESOURCE
icons.qrc
)

qt5_add_resources(RESOURCE_ADDED ${RESOURCE})

SET(EXE_NAME "tmogui${TAG_VERSION}" )
ADD_EXECUTABLE(${EXE_NAME}
  ${TMO_SOURCES}
  ${RESOURCE_ADDED}
)

target_link_libraries(tmogui Qt5::Widgets Qt5::PrintSupport)
add_custom_command( TARGET ${EXE_NAME} POST_BUILD
  COMMAND cp -r tmogui ../
)
