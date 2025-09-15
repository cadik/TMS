CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

INCLUDE_DIRECTORIES(
../tmolib/
.
)

LINK_LIBRARIES (
gsl
gslcblas
tmo
)

SET(TMO_SOURCES
 TMOMantiuk08.cpp
 TMOPlugin.cpp
 cqp/cqpminimizer
 cqp/initial_point.cpp
 cqp/mg_pdip.cpp
)


SET(TMOWARD_HEADERS
 TMOMantiuk08.h
 TMOPlugin.h
 cqp/initial_point.h
 cqp/gsl_cqp.h
)
#set(TRG TMOMantiuk08)
#add_executable(${TRG} ${TRG}.cpp cqp/cqpminimizer.cpp cqp/initial_point.cpp cqp/mg_pdip.cpp "${GETOPT_OBJECT}")

add_library( mantiuk08 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET mantiuk08 POST_BUILD 
	COMMAND cp -f libmantiuk08.so ../TMOCmd/mantiuk08.tml
	COMMAND cp -f libmantiuk08.so ../TMOgui/mantiuk08.tml 
	COMMAND cp -f libmantiuk08.so ../mantiuk08.tml 
#	COMMAND rm -f libmantiuk08.so 
	)


