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
 TMOMikamo14.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOMikamo14.h
 TMOPlugin.h
)


add_library( mikamo14 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET mikamo14 POST_BUILD 
	COMMAND cp -f libmikamo14.so ../TMOCmd/mikamo14.tml
	COMMAND cp -f libmikamo14.so ../TMOgui/mikamo14.tml 
	COMMAND cp -f libmikamo14.so ../mikamo14.tml 
#	COMMAND rm -f libmikamo14.so 
	)


