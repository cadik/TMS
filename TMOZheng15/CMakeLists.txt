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
fftw3         # added
pthread	      # added
fftw3_threads # added
)

SET(TMO_SOURCES
 TMOZheng15.cpp
 TMOPlugin.cpp  
 ../tmolib/poisson_pde/laplace.cpp
# ../pftools/pde.cpp
# ../pftools/pde_fft.cpp
)


SET(TMOEFFICIENT_HEADERS
 TMOZheng15.h
 TMOPlugin.h 
 ../tmolib/poisson_pde/laplace.h
# ../pftools/pde.h
# ../pftools/array2d.h
# ../pftools/config.h
# ../pftools/pfstmo.h
)

add_library( Zheng15 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET Zheng15 POST_BUILD 
	COMMAND cp -f libZheng15.so ../TMOCmd/Zheng15.tml
	COMMAND cp -f libZheng15.so ../TMOgui/Zheng15.tml 
	COMMAND cp -f libZheng15.so ../Zheng15.tml 
#	COMMAND rm -f libZheng15.so 
	)


