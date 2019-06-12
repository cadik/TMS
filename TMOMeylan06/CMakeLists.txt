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
 TMOMeylan06.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOMeylan06.h
 TMOPlugin.h
)


add_library( meylan06 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET meylan06 POST_BUILD 
	COMMAND cp -f libmeylan06.so ../TMOCmd/meylan06.tml
	COMMAND cp -f libmeylan06.so ../TMOgui/meylan06.tml 
	COMMAND cp -f libmeylan06.so ../meylan06.tml 
#	COMMAND rm -f libmeylan06.so 
	)


