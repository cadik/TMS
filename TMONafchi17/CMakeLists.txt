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
 TMONafchi17.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMONafchi17.h
 TMOPlugin.h
)


add_library( nafchi17 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET nafchi17 POST_BUILD 
	COMMAND cp -f libnafchi17.so ../TMOCmd/nafchi17.tml
	COMMAND cp -f libnafchi17.so ../TMOgui/nafchi17.tml 
	COMMAND cp -f libnafchi17.so ../nafchi17.tml 
#	COMMAND rm -f libnafchi17.so 
	)


