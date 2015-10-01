CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

FIND_PACKAGE(TIFF REQUIRED)

#OpenEXR
FIND_PATH(OPENEXR_INCLUDE_DIR ImfRgbaFile.h 
	/usr/include/OpenEXR 
	/usr/local/include)

IF(NOT OPENEXR_INCLUDE_DIR)
 MESSAGE(FATAL_ERROR " OpenEXR headers not found. (try setting OPENEXR_INCLUDE_DIR variable to the path to OpenEXR headers)")
ENDIF(NOT OPENEXR_INCLUDE_DIR)

FIND_LIBRARY(OPENEXR_LIBRARY 
	NAMES IlmImf
	PATHS 
	/usr/lib
	/lib
	/usr/local/lib
	${GLOB_PATHS_LIB}
#	NO_CMAKE_ENVIRONMENT_PATH
#	NO_CMAKE_PATH
#	NO_SYSTEM_ENVIRONMENT_PATH
#	NO_CMAKE_SYSTEM_PATH
	)

IF(NOT OPENEXR_LIBRARY)
 MESSAGE(FATAL_ERROR " OpenEXR library not found. (try setting OPENEXR_LIBRARY variable to the path to libIlmImf.so)")
ENDIF(NOT OPENEXR_LIBRARY)

IF(OPENEXR_INCLUDE_DIR AND OPENEXR_LIBRARY)
 MESSAGE("OpenEXR found OK")
ENDIF(OPENEXR_INCLUDE_DIR AND OPENEXR_LIBRARY)
#OpenEXR

INCLUDE_DIRECTORIES(
${CMAKE_BINARY_DIR}
${OPENEXR_INCLUDE_DIR}
../tmolib/
.
)

#Set up libraries to link with all found packages
LINK_LIBRARIES (
${TIFF_LIBRARY}
${OPENEXR_LIBRARY}
)


SET(TMOLIB_SOURCES
../tmolib/TMO.cpp 
../tmolib/TMOImage.cpp 
../tmolib/TMOParameter.cpp 
../TMOLinux/TMOLinux.cpp 
../tmolib/TMORadiance.cpp 
../tmolib/freejpeghdr.cpp
../tmolib/TMO_Tonemap.cpp
)

SET(TMOLIB_HEADERS
../tmolib/TMO.cpp 
../tmolib/TMOImage.cpp 
../tmolib/TMOParameter.cpp 
../TMOLinux/TMOLinux.cpp 
../tmolib/TMORadiance.cpp 
../tmolib/freejpeghdr.cpp
../tmolib/TMO_Tonemap.cpp
)


add_library( tmo SHARED ${TMOLIB_SOURCES} )
#add_library( tmo STATIC ${TMOLIB_SOURCES} )

