CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -fPIC")
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)


FIND_PACKAGE(OpenCV REQUIRED)
IF(OpenCV_VERSION VERSION_LESS "2.4.9")
  message(FATAL_ERROR "ERROR: Can't find OpenCV 2.4.9+")
ELSE(OpenCV_VERSION VERSION_LESS "2.4.9") 
	MESSAGE("OpenCV found OK")
ENDIF()

ADD_DEFINITIONS(-fPIC)
ADD_LIBRARY(libpallas  STATIC IMPORTED GLOBAL)
SET_TARGET_PROPERTIES(libpallas PROPERTIES
   IMPORTED_LOCATION /usr/local/lib/libpallas.a)
ADD_LIBRARY(libgtest_main STATIC IMPORTED GLOBAL)
SET_TARGET_PROPERTIES(libgtest_main PROPERTIES
    IMPORTED_LOCATION /usr/local/lib/libgtest_main.a)
ADD_LIBRARY(libgtest STATIC IMPORTED GLOBAL)
SET_TARGET_PROPERTIES(libgtest PROPERTIES
    IMPORTED_LOCATION /usr/local/lib/libgtest.a)

#add_library(A STATIC ../../pallas-solver-master/build/lib/libpallas.a)
#add_library(B STATIC ../../pallas-solver-master/build/lib/libgtest_main.a)
#add_library(C STATIC ../../pallas-solver-master/build/lib/libgtest.a)
INCLUDE_DIRECTORIES(
/include
${OpenCV_INCLUDE_DIRS}
${Dlib_INCLUDE_DIRS}
../tmolib/
/usr/local/include
.
)

LINK_LIBRARIES (
${OpenCV_LIBS}
tmo
ceres
dlib
libpallas
#libgtest_main
#usr/local/lib/libpallas.a
libgtest_main
libgtest


)



SET(TMO_SOURCES
 TMOTao17.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOTao17.h
 TMOPlugin.h
)


add_library( Tao17 SHARED ${TMO_SOURCES} )

 add_custom_command( TARGET Tao17 POST_BUILD 
 	COMMAND cp -f libTao17.so ../TMOCmd/Tao17.tml
 	COMMAND cp -f libTao17.so ../TMOgui/Tao17.tml 
 	COMMAND cp -f libTao17.so ../Tao17.tml 
 #	COMMAND rm -f libHu14.so 
  	)


