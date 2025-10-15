CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
    ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

INCLUDE_DIRECTORIES(
    ../tmolib
    .
)

LINK_LIBRARIES (
    tmo
)

SET(TMO_SOURCES
    TMOEilertsen15.cpp
    TMOPlugin.cpp
)

SET(TMO_HEADERS
    TMOEilertsen15.h
    TMOPlugin.h
)

ADD_LIBRARY(Eilertsen15 SHARED ${TMO_SOURCES})

add_custom_command( TARGET Eilertsen15 POST_BUILD 
	COMMAND cp -f libEilertsen15.so ../TMOCmd/Eilertsen15.tml
	COMMAND cp -f libEilertsen15.so ../TMOgui/Eilertsen15.tml 
	COMMAND cp -f libEilertsen15.so ../Eilertsen15.tml 
#	COMMAND rm -f libEilertsen15.so 
)
