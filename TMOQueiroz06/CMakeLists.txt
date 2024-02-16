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
 TMOQueiroz06.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOQueiroz06.h
 TMOPlugin.h
)


add_library( queiroz06 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET queiroz06 POST_BUILD 
	COMMAND cp -f libqueiroz06.so ../TMOCmd/queiroz06.tml
	COMMAND cp -f libqueiroz06.so ../TMOgui/queiroz06.tml 
	COMMAND cp -f libqueiroz06.so ../queiroz06.tml 
#	COMMAND rm -f libqueiroz06.so 
	)


