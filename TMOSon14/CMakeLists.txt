CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -fPIC")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -fPIC -Wall -pedantic -Wfloat-equal -Wshadow -DLINUX")
# flags from qpOASES, make_linux.mk
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion -Wsign-conversion -D__USE_LONG_INTEGERS__ -D__USE_LONG_FINTS__ -D__NO_COPYRIGHT__")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} -O3 -finline-functions")
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)
	

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

# policy for new style of shared library linking
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


##########
#Eigen
##########
#set(EIGEN_DIR /usr/local/include/eigen3/)
#include_directories(${EIGEN_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})

find_package(CGAL QUIET)
if ( CGAL_FOUND )
	message("cgal founded	")
	include( ${CGAL_USE_FILE} )

	include( CGAL_CreateSingleSourceCGALProgram )
endif()

FIND_PACKAGE(OpenCV REQUIRED)
IF(OpenCV_VERSION VERSION_LESS "2.4.9")
  message(FATAL_ERROR "ERROR: Can't find OpenCV 2.4.9+")
ELSE(OpenCV_VERSION VERSION_LESS "2.4.9") 
	MESSAGE("OpenCV found OK")
ENDIF()

INCLUDE_DIRECTORIES(
${OpenCV_INCLUDE_DIRS}
../tmolib/
# /usr/local/lib
# /cgal_includes
./qpOASES/include
.
)

#target_link_libraries(pallas)
##########
#Ceres
##########
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

find_library(QPOASES_LIB qpOASES ${CMAKE_CURRENT_LIST_DIR}/qpOASES/bin)

LINK_LIBRARIES (
${OpenCV_LIBS}
tmo
# ${CGAL_LIBRARIES}
# ${CGAL_3RD_PARTY_LIBRARIES}
#ceres
${QPOASES_LIB}
)

SET(TMO_SOURCES
 TMOSon14.cpp
 TMOPlugin.cpp
 L0minimization.cpp
 baseAndDetailDecomposition.cpp
 # cgalOptimization.cpp
 #pallasOptimization.cpp
 qpOasesOptimalization.cpp
)

SET(TMOWARD_HEADERS
 TMOSon14.h
 TMOPlugin.h
 L0minimization.h
 baseAndDetailDecomposition.h
 #pallasOptimization.h
 # cgalOptimization.h
 qpOasesOptimalization.h
)

add_library( Son14 SHARED ${TMO_SOURCES} )
 add_custom_command( TARGET Son14 POST_BUILD 
 	COMMAND cp -f libSon14.so ../TMOCmd/Son14.tml
 	COMMAND cp -f libSon14.so ../TMOgui/Son14.tml 
 	COMMAND cp -f libSon14.so ../Son14.tml 
 #	COMMAND rm -f libSon14.so 
  	)



