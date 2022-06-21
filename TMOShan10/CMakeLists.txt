CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)
# TODO smazat
SET(CMAKE_BUILD_TYPE Debug)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

##########
#Eigen
##########
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
../tmolib/
./imresize/
.
)

LINK_LIBRARIES (
${OpenCV_LIBS}
tmo
)

SET(TMO_SOURCES
 TMOShan10.cpp
 TMOPlugin.cpp
 MatUtil.cpp 
)


SET(TMOWARD_HEADERS
 TMOShan10.h
 TMOPlugin.h
 MatUtil.h 
)


add_library( shan10 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET shan10 POST_BUILD 
	COMMAND cp -f libshan10.so ../TMOCmd/shan10.tml
	COMMAND cp -f libshan10.so ../TMOgui/shan10.tml 
	COMMAND cp -f libshan10.so ../shan10.tml 
#	COMMAND rm -f libshan10.so 
	)


