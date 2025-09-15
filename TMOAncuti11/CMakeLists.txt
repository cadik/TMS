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
 TMOAncuti11.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOAncuti11.h
 TMOPlugin.h
)


add_library( ancuti11 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET ancuti11 POST_BUILD 
	COMMAND cp -f libancuti11.so ../TMOCmd/ancuti11.tml
	COMMAND cp -f libancuti11.so ../TMOgui/ancuti11.tml 
	COMMAND cp -f libancuti11.so ../ancuti11.tml 
#	COMMAND rm -f libancuti11.so 
	)


