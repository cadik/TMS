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
 TMOAncuti19.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOAncuti19.h
 TMOPlugin.h
)


add_library( ancuti19 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET ancuti19 POST_BUILD 
	COMMAND cp -f libancuti19.so ../TMOCmd/ancuti19.tml
	COMMAND cp -f libancuti19.so ../TMOgui/ancuti19.tml 
	COMMAND cp -f libancuti19.so ../ancuti19.tml 
#	COMMAND rm -f libancuti19.so 
	)


