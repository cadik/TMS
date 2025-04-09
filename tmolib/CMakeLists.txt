CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

FIND_PACKAGE(TIFF REQUIRED)
FIND_PACKAGE(PNG REQUIRED)

#OpenEXR
FIND_PATH(OPENEXR_INCLUDE_DIR ImfRgbaFile.h
	/usr/include/OpenEXR
	/usr/local/include)
IF(NOT OPENEXR_INCLUDE_DIR)
	MESSAGE(FATAL_ERROR " OpenEXR headers not found. (try setting OPENEXR_INCLUDE_DIR variable to the path to OpenEXR headers)")
ENDIF(NOT OPENEXR_INCLUDE_DIR)

# First, try to find just the right config files
find_package(Imath CONFIG)
if (NOT TARGET Imath::Imath)
    # Couldn't find Imath::Imath, maybe it's older and has IlmBase?
    find_package(IlmBase CONFIG)
endif ()
find_package(OpenEXR CONFIG)
#OpenEXR

FIND_PACKAGE(OpenCV REQUIRED)
IF(OpenCV_VERSION VERSION_LESS "2.4.9")
  message(FATAL_ERROR "ERROR: Can't find OpenCV 2.4.9+")
ELSE(OpenCV_VERSION VERSION_LESS "2.4.9")
	MESSAGE("OpenCV found OK")
ENDIF()

INCLUDE_DIRECTORIES(
${CMAKE_BINARY_DIR}
${OPENEXR_INCLUDE_DIR}
${PNG_INCLUDE_DIR}
${OpenCV_INCLUDE_DIRS}
../tmolib/
.
)

#Set up libraries to link with all found packages
LINK_LIBRARIES (
${TIFF_LIBRARY}
# For OpenEXR/Imath 3.x:
	$<$<TARGET_EXISTS:OpenEXR::OpenEXR>:OpenEXR::OpenEXR>
	$<$<TARGET_EXISTS:Imath::Imath>:Imath::Imath>
	$<$<TARGET_EXISTS:Imath::Half>:Imath::Half>
# For OpenEXR 2.4/2.5:
	$<$<TARGET_EXISTS:OpenEXR::IlmImf>:OpenEXR::IlmImf>
	$<$<TARGET_EXISTS:IlmBase::Imath>:IlmBase::Imath>
	$<$<TARGET_EXISTS:IlmBase::Half>:IlmBase::Half>
	$<$<TARGET_EXISTS:IlmBase::IlmThread>:IlmBase::IlmThread>
	$<$<TARGET_EXISTS:IlmBase::Iex>:IlmBase::Iex>
${PNG_LIBRARY}
${OpenCV_LIBS}
)


SET(TMOLIB_SOURCES
../tmolib/TMO.cpp
../tmolib/TMOv.cpp
../tmolib/TMOImage.cpp 
../tmolib/TMOParameter.cpp 
../TMOLinux/TMOLinux.cpp 
../tmolib/TMORadiance.cpp 
../tmolib/freejpeghdr.cpp
../tmolib/TMO_Tonemap.cpp
../tmolib/TMOVideo.cpp 
)

SET(TMOLIB_HEADERS
../tmolib/TMO.cpp 
../tmolib/TMOv.cpp
../tmolib/TMOImage.cpp 
../tmolib/TMOParameter.cpp 
../TMOLinux/TMOLinux.cpp 
../tmolib/TMORadiance.cpp 
../tmolib/freejpeghdr.cpp
../tmolib/TMO_Tonemap.cpp
../tmolib/TMOVideo.h
)


add_library( tmo SHARED ${TMOLIB_SOURCES} )
#add_library( tmo STATIC ${TMOLIB_SOURCES} )

