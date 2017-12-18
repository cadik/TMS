CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

FIND_PACKAGE(OpenCV REQUIRED)
IF(OpenCV_VERSION VERSION_LESS "2.4.9")
  message(FATAL_ERROR "ERROR: Can't find OpenCV 2.4.9+")
ELSE(OpenCV_VERSION VERSION_LESS "2.4.9") 
	MESSAGE("OpenCV found OK")
ENDIF()

##########################
# opencv
##########################
set(OpenCV_DIR "/usr/local/opencv/share/OpenCV/OpenCVConfig.cmake"
    CACHE PATH "The path where OpenCVConfig.cmake is placed")
find_package(OpenCV REQUIRED)

if(OpenCV_FOUND)
  message(STATUS "Found 'opencv library' ${OpenCV_VERSION}")
else()
  message(SEND_ERROR "Couldn't find OpenCV")
endif()

INCLUDE_DIRECTORIES(
${OpenCV_INCLUDE_DIRS}
../tmolib/
.
)

LINK_LIBRARIES (
${OpenCV_LIBS}
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
 TMOSon14.cpp
 TMOPlugin.cpp
 L0minimization.cpp
 baseAndDetailDecomposition.cpp
)


SET(TMOWARD_HEADERS
 TMOSon14.h
 TMOPlugin.h
 L0minimization.h
 baseAndDetailDecomposition.h
)


add_library( Son14 SHARED ${TMO_SOURCES} )
 add_custom_command( TARGET Son14 POST_BUILD 
 	COMMAND cp -f libSon14.so ../TMOCmd/Son14.tml
 	COMMAND cp -f libSon14.so ../TMOgui/Son14.tml 
 	COMMAND cp -f libSon14.so ../Son14.tml 
 #	COMMAND rm -f libSon14.so 
  	)


