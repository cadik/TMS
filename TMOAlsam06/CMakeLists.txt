CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

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

##########
#Eigen
##########
#set(EIGEN_DIR /usr/local/include/eigen3/)
#include_directories(${EIGEN_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})

SET(TMO_SOURCES
 TMOAlsam06.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOAlsam06.h
 TMOPlugin.h
)


add_library( alsam06 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET alsam06 POST_BUILD 
	COMMAND cp -f libalsam06.so ../TMOCmd/alsam06.tml
	COMMAND cp -f libalsam06.so ../TMOgui/alsam06.tml 
	COMMAND cp -f libalsam06.so ../alsam06.tml 
#	COMMAND rm -f libalsam06.so 
	)


