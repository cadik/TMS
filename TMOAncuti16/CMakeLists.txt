CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)


FIND_PACKAGE(OpenCV REQUIRED)
IF (OpenCV_FOUND)
  MESSAGE("OpenCV found OK")
ELSE(OpenCV_FOUND)
  MESSAGE(FATAL_ERROR "OpenCV library not found.")
ENDIF(OpenCV_FOUND)
message(STATUS "opencv DIRECTORY: " ${OpenCV_LIBS})
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
 TMOAncuti16.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOAncuti16.h
 TMOPlugin.h
)


add_library( yourOperator SHARED ${TMO_SOURCES} )
add_custom_command( TARGET yourOperator POST_BUILD 
	COMMAND cp -f libyourOperator.so ../TMOCmd/Ancuti16.tml
	COMMAND cp -f libyourOperator.so ../TMOgui/Ancuti16.tml 
	COMMAND cp -f libyourOperator.so ../Ancuti16.tml 
#	COMMAND rm -f libyourOperator.so 
	)


