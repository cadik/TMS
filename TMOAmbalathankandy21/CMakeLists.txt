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
 TMOAmbalathankandy21.cpp
 TMOPlugin.cpp
)


SET(HEADERS
 TMOAmbalathankandy21.h
 TMOPlugin.h
)

add_library( Ambalathankandy21 SHARED ${SOURCES} ${HEADERS})
add_custom_command( TARGET Ambalathankandy21 POST_BUILD 
	COMMAND cp -f libAmbalathankandy21.so ../TMOCmd/Ambalathankandy21.tml
	COMMAND cp -f libAmbalathankandy21.so ../TMOgui/Ambalathankandy21.tml 
	COMMAND cp -f libAmbalathankandy21.so ../Ambalathankandy21.tml 
#	COMMAND rm -f libqueiroz06.so 
	)
