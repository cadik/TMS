CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

FIND_PACKAGE(OpenCV REQUIRED)
IF(OpenCV_VERSION VERSION_LESS "2.4.9")
  message(FATAL_ERROR "ERROR: Can't find OpenCV 2.4.9+")
ELSE(OpenCV_VERSION VERSION_LESS "2.4.9") 
	MESSAGE("OpenCV found OK")
ENDIF()

INCLUDE_DIRECTORIES(
../tmolib/
.
)

LINK_LIBRARIES (
${OpenCV_LIBS}
tmo
)

SET(TMO_SOURCES
 TMOGastal11.cpp
 TMOPlugin.cpp
 ColorspaceConverter.cpp
 MatUtil.cpp
)


SET(TMOWARD_HEADERS
 TMOGastal11.h
 TMOPlugin.h
 ColorspaceConverter.h
 MatUtil.h
)


add_library( gastal11 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET gastal11 POST_BUILD 
	COMMAND cp -f libgastal11.so ../TMOCmd/gastal11.tml
	COMMAND cp -f libgastal11.so ../TMOgui/gastal11.tml 
	COMMAND cp -f libgastal11.so ../gastal11.tml 
#	COMMAND rm -f libgastal11.so 
	)


