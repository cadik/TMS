CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

INCLUDE_DIRECTORIES(
../tmolib/
.
)

LINK_LIBRARIES (
tmo
)

SET(TMO_SOURCES
 TMO2scaleToneManagement.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMO2scaleToneManagement.h
 TMOPlugin.h
)


add_library( TMO2scaleToneManagement SHARED ${TMO_SOURCES} )
add_custom_command( TARGET TMO2scaleToneManagement POST_BUILD 
	COMMAND cp -f libTMO2scaleToneManagement.so ../TMOCmd/TMO2scaleToneManagement.tml
	COMMAND cp -f libTMO2scaleToneManagement.so ../TMOgui/TMO2scaleToneManagement.tml 
	COMMAND cp -f libTMO2scaleToneManagement.so ../TMO2scaleToneManagement.tml 
#	COMMAND rm -f libTMO2scaleToneManagement.so 
	)


