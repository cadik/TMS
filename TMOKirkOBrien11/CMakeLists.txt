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
fftw3
)


SET(TMO_SOURCES
 TMOKirkOBrien11.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOKirkOBrien11.h
 TMOPlugin.h
)


add_library( kirkOBrien11 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET kirkOBrien11 POST_BUILD 
	COMMAND cp -f libkirkOBrien11.so ../TMOCmd/kirkOBrien11.tml
	COMMAND cp -f libkirkOBrien11.so ../TMOgui/kirkOBrien11.tml 
	COMMAND cp -f libkirkOBrien11.so ../kirkOBrien11.tml 
#	COMMAND rm -f libkirkOBrien11.so 
	)


