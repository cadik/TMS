CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
# optimization flag is important, program runs a lot faster
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -fPIC -Wall -pedantic -Wfloat-equal -Wshadow -DLINUX")
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

##########
#Eigen
##########
#set(EIGEN_DIR /usr/local/include/eigen3/)
#include_directories(${EIGEN_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})

FIND_PACKAGE(OpenCV REQUIRED)
IF(OpenCV_VERSION VERSION_LESS "2.4.9")
  message(FATAL_ERROR "ERROR: Can't find OpenCV 2.4.9+")
ELSE(OpenCV_VERSION VERSION_LESS "2.4.9") 
	MESSAGE("OpenCV found OK")
ENDIF()

INCLUDE_DIRECTORIES(
${OpenCV_INCLUDE_DIRS}
../tmolib/
.
)

LINK_LIBRARIES (
${OpenCV_LIBS}
tmo
)

SET(TMO_SOURCES
 TMOKou15.cpp
 TMOPlugin.cpp
 L0minimization.cpp
)


SET(TMOWARD_HEADERS
 TMOKou15.h
 TMOPlugin.h
 L0minimization.h
)


add_library( kou15 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET kou15 POST_BUILD 
	COMMAND cp -f libkou15.so ../TMOCmd/kou15.tml
	COMMAND cp -f libkou15.so ../TMOgui/kou15.tml 
	COMMAND cp -f libkou15.so ../kou15.tml 
#	COMMAND rm -f libkou15.so 
	)
