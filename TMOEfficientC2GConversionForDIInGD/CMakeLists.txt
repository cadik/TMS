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
)

SET(TMO_SOURCES
 TMOEfficientC2GConversionForDIInGD.cpp
 TMOPlugin.cpp  
 ../poisson_pde/laplace.cpp
)


SET(TMOEFFICIENT_HEADERS
 TMOEfficientC2GConversionForDIInGD.h
 TMOPlugin.h 
 ../poisson_pde/laplace.h
)

#SET(GCC_COVERAGE_LINK_FLAGS "-lm -lfftw3 -pthread -lfftw3_threads")
#SET (CMAKE_CXX_FLAGS "-lm -lfftw3 -pthread -lfftw3_threads")
#SET (CMAKE_EXE_LINKER_FLAGS "-lm -lfftw3 -pthread -lfftw3_threads")
#target_compile_options(tmo PRIVATE -lm -lfftw3 -pthread -lfftw3_threads)
#target_link_libraries(-lm -lfftw3 -pthread -lfftw3_threads)

add_library( TMOEfficientC2GConversionForDIInGD SHARED ${TMO_SOURCES} )
add_custom_command( TARGET TMOEfficientC2GConversionForDIInGD POST_BUILD 
	COMMAND cp -f libTMOEfficientC2GConversionForDIInGD.so ../TMOCmd/TMOEfficientC2GConversionForDIInGD.tml
	COMMAND cp -f libTMOEfficientC2GConversionForDIInGD.so ../TMOgui/TMOEfficientC2GConversionForDIInGD.tml 
	COMMAND cp -f libTMOEfficientC2GConversionForDIInGD.so ../TMOEfficientC2GConversionForDIInGD.tml 
#	COMMAND rm -f libTMOEfficientC2GConversionForDIInGD.so 
	)


