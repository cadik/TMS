TEMPLATE = app
INCLUDEPATH += .

unix {
	DEFINES = LINUX
}
          
isEmpty(OPENEXRPATH){
	OPENEXRPATH = /include/OpenEXR
}


!exists($${OPENEXRPATH}/ImfRgbaFile.h){
 error("File $${OPENEXRPATH}/ImfRgbaFile.h not found. Probably wrong path to OpenEXR headers. Try: qmake "OPENEXRPATH=path/to/OpenExr" ")
}

!exists($${OPENEXRPATH}/ImfStringAttribute.h){
 error("File $${OPENEXRPATH}/ImfStringAttribute.h not found. Probably wrong path to OpenEXR headers. Try: qmake "OPENEXRPATH=path/to/OpenExr" ")
}

!exists($${OPENEXRPATH}/ImfMatrixAttribute.h){
 error("File $${OPENEXRPATH}/ImfMatrixAttribute.h not found. Probably wrong path to OpenEXR headers. Try: qmake "OPENEXRPATH=path/to/OpenExr" ")
}

!exists($${OPENEXRPATH}/ImfArray.h){
 error("File $${OPENEXRPATH}/ImfArray.h not found. Probably wrong path to OpenEXR headers. Try: qmake "OPENEXRPATH=path/to/OpenExr" ")
}

isEmpty(LIBJPEGPATH){
	LIBJPEGPATH = /usr/lib
}

!exists($${LIBJPEGPATH}/libjpeg.so*){
 error("File $${LIBJPEGPATH}/libjpeg.so not found. Probably wrong path to libjpeg files. Try: qmake "LIBJPEGPATH=path/to/libjpeg" ")
}


INCLUDEPATH += $$OPENEXRPATH
LIBS += LIBJPEGPATH


LIBS = -ltiff -lIlmImf -ldl -lqassistantclient

# Input
HEADERS += gamma.h \
           lqstring.h \
           resources.ui.h \
           resources6.ui.h \
           TMOFilters.h \
           TMOGUIAdjust.h \
           TMOGUIAdjustValues.h \
           TMOGUIBitmap.h \
           TMOGUICenterView.h \
           TMOGUIFileToolBar.h \
           TMOGUIFilters.h \
           TMOGUIHisto.h \
           TMOGUIImage.h \
           TMOGUIInfo.h \
           TMOGUIInfoPanel.h \
           TMOGUIInfoTool.h \
           TMOGUIInfoToolBar.h \
           TMOGUILineResizer.h \
           TMOGUIMenu.h \
           TMOGUIOutput.h \
           TMOGUIParameters.h \
           TMOGUIParametersItem.h \
           TMOGUIProgressBar.h \
           TMOGUIResource.h \
           TMOGUIRightBar.h \
           TMOGUISaveDialog.h \
           TMOGUIStatistics.h \
           TMOGUIStatus.h \
           TMOGUIToneMapping.h \
           TMOGUIToneSlider.h \
           TMOGUITransformation.h \
           tmoguiwindow.h \
           TMOGUIZoomTool.h
INTERFACES += resources.ui \
              resources1.ui \
              resources2.ui \
              resources3.ui \
              resources4.ui \
              resources5.ui \
              resources6.ui
SOURCES += lqstring.cpp \
           main.cpp \
           TMOFilters.cpp \
           TMOGUIAdjust.cpp \
           TMOGUIAdjustValues.cpp \
           TMOGUIBitmap.cpp \
           TMOGUICenterView.cpp \
           TMOGUIFileToolBar.cpp \
           TMOGUIFilters.cpp \
           TMOGUIHisto.cpp \
           TMOGUIImage.cpp \
           TMOGUIInfo.cpp \
           TMOGUIInfoPanel.cpp \
           TMOGUIInfoTool.cpp \
           TMOGUIInfoToolBar.cpp \
           TMOGUILineResizer.cpp \
           TMOGUIMenu.cpp \
           TMOGUIOutput.cpp \
           TMOGUIParameters.cpp \
           TMOGUIParametersItem.cpp \
           TMOGUIProgressBar.cpp \
           TMOGUIRightBar.cpp \
           TMOGUISaveDialog.cpp \
           TMOGUIStatistics.cpp \
           TMOGUIStatus.cpp \
           TMOGUIToneMapping.cpp \
           TMOGUIToneSlider.cpp \
           TMOGUITransformation.cpp \
           tmoguiwindow.cpp \
           TMOGUIZoomTool.cpp \
           ../tmolib/TMO.cpp \
           ../tmolib/TMOImage.cpp \
           ../tmolib/TMOParameter.cpp \
           ../TMOLinux/TMOLinux.cpp \
           ../tmolib/TMORadiance.cpp \
           ../tmolib/freejpeghdr.cpp

win32 {	   
	SOURCES+=../TMOW32/TMOW32.cpp
}


