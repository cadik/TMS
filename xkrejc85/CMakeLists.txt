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

SET(SOURCES
 xkrejc85.cpp
 TMOPlugin.cpp
)


SET(HEADERS
 xkrejc85.h
 TMOPlugin.h
)

add_library( xkrejc85 SHARED ${SOURCES} ${HEADERS})
add_custom_command( TARGET xkrejc85 POST_BUILD 
	COMMAND cp -f libxkrejc85.so ../TMOCmd/xkrejc85.tml
	COMMAND cp -f libxkrejc85.so ../TMOgui/xkrejc85.tml 
	COMMAND cp -f libxkrejc85.so ../xkrejc85.tml 
#	COMMAND rm -f libqueiroz06.so 
	)
