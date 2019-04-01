CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_BUILD_TYPE Debug)

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
 TMOCheryl11.cpp
 TMOPlugin.cpp
 Cluster.cpp
 Graph.cpp
)


SET(TMOWARD_HEADERS
 TMOCheryl11.h
 TMOPlugin.h
 Cluster.h
 Graph.h
)


add_library( cheryl11 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET cheryl11 POST_BUILD 
	COMMAND cp -f libcheryl11.so ../TMOCmd/cheryl11.tml
	#COMMAND cp -f libcheryl11.so ../TMOgui/cheryl11.tml 
	COMMAND cp -f libcheryl11.so ../cheryl11.tml 
#	COMMAND rm -f libcheryl11.so 
	)


