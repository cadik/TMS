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

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

SET(TMO_SOURCES
 TMOHsin11.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOHsin11.h
 TMOPlugin.h
)


add_library( hsin11 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET hsin11 POST_BUILD 
	COMMAND cp -f libhsin11.so ../TMOCmd/hsin11.tml
	COMMAND cp -f libhsin11.so ../TMOgui/hsin11.tml 
	COMMAND cp -f libhsin11.so ../hsin11.tml 
#	COMMAND rm -f libhsin11.so 
	)


