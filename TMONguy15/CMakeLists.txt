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
-fopenmp
-lpng
-lnlopt
)

SET(TMO_SOURCES
 TMONguy15.cpp
 TMOPlugin.cpp
 mainprepare.cpp
 optimalizer.cpp
 transform.cpp
 color2grey.cpp
)


SET(TMOWARD_HEADERS
 TMONguy15.h
 TMOPlugin.h
 mainprepare.h
 optimalizer.h
 transform.h
 am_sysdep.h
 am_types.h
 math_const.h
 color2grey.h
)


add_library( nguy15 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET nguy15 POST_BUILD 
	COMMAND cp -f libnguy15.so ../TMOCmd/nguy15.tml
	COMMAND cp -f libnguy15.so ../TMOgui/nguy15.tml 
	COMMAND cp -f libnguy15.so ../nguy15.tml 
#	COMMAND rm -f libnguy15.so 
	)


