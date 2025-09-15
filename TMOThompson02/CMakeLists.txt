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
 TMOThompson02.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOThompson02.h
 TMOPlugin.h
)


add_library( Thompson02 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET Thompson02 POST_BUILD 
	COMMAND cp -f libThompson02.so ../TMOCmd/Thompson02.tml
	COMMAND cp -f libThompson02.so ../TMOgui/Thompson02.tml 
	COMMAND cp -f libThompson02.so ../Thompson02.tml 
#	COMMAND rm -f libThompson02.so 
	)


