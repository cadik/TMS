CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

INCLUDE_DIRECTORIES(
../tmolib/
.
)

LINK_LIBRARIES (
tmo
)

SET(TMO_SOURCES
 TMOEisemann04.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOEisemann04.h
 TMOPlugin.h
)


add_library( eisemann04 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET eisemann04 POST_BUILD 
	COMMAND cp -f libeisemann04.so ../TMOCmd/eisemann04.tml
	COMMAND cp -f libeisemann04.so ../TMOgui/eisemann04.tml 
	COMMAND cp -f libeisemann04.so ../eisemann04.tml 
#	COMMAND rm -f libeisemann04.so 
	)



