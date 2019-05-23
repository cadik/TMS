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
 TMOFarbman08.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOFarbman08.h
 TMOPlugin.h
)


add_library( farbman08 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET farbman08 POST_BUILD 
	COMMAND cp -f libfarbman08.so ../TMOCmd/farbman08.tml
	COMMAND cp -f libfarbman08.so ../TMOgui/farbman08.tml 
	COMMAND cp -f libfarbman08.so ../farbman08.tml 
#	COMMAND rm -f libfarbman08.so 
	)


