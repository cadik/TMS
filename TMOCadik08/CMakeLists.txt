CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
	ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

FIND_PACKAGE(OpenCL REQUIRED)
IF (OpenCL_FOUND)
	MESSAGE("OpenCL found OK")
ELSE(OpenCL_FOUND)
	MESSAGE(FATAL_ERROR "OpenCL library not found.")
ENDIF(OpenCL_FOUND)

INCLUDE_DIRECTORIES(
	../tmolib/
	.
)

LINK_LIBRARIES (
	${OpenCL_LIBRARY}
	${OpenCL_LIBRARIES}
	${OpenCL_3RD_PARTY_LIBRARIES}
	tmo
)

SET(TMO_SOURCES
	TMOCadik08.cpp
	TMOPlugin.cpp
	coloroid.cpp
	spacefill.cpp
	quadtree.cpp
	morton.cpp
)

SET(TMOWARD_HEADERS
	TMOCadik08.h
	TMOPlugin.h
	coloroid.h
	spacefill.h
	quadtree.h
	morton.h
)

add_library( Cadik08 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET Cadik08 POST_BUILD
	COMMAND cp -f libCadik08.so ../TMOCmd/Cadik08.tml
	COMMAND cp -f libCadik08.so ../TMOgui/Cadik08.tml 
	COMMAND cp -f libCadik08.so ../Cadik08.tml
#	COMMAND rm -f libCadik08.so
)

