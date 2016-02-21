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
 TMOEfficientC2GConversionForDIInGD.cpp
 TMOPlugin.cpp  
 ../poisson_pde/laplace.cpp
 ../pftools/pde.cpp
 ../pftools/pde_fft.cpp
)


SET(TMOEFFICIENT_HEADERS
 TMOEfficientC2GConversionForDIInGD.h
 TMOPlugin.h 
 ../poisson_pde/laplace.h
 ../pftools/pde.h
 ../pftools/array2d.h
 ../pftools/config.h
 ../pftools/pfstmo.h
)

add_library( TMOEfficientC2GConversionForDIInGD SHARED ${TMO_SOURCES} )
add_custom_command( TARGET TMOEfficientC2GConversionForDIInGD POST_BUILD 
	COMMAND cp -f libTMOEfficientC2GConversionForDIInGD.so ../TMOCmd/TMOEfficientC2GConversionForDIInGD.tml
	COMMAND cp -f libTMOEfficientC2GConversionForDIInGD.so ../TMOgui/TMOEfficientC2GConversionForDIInGD.tml 
	COMMAND cp -f libTMOEfficientC2GConversionForDIInGD.so ../TMOEfficientC2GConversionForDIInGD.tml 
#	COMMAND rm -f libTMOEfficientC2GConversionForDIInGD.so 
	)


