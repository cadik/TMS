CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

# Find OpenCV
find_package(OpenCV REQUIRED)

cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0079 NEW)

# Find Torch
find_package(Torch 2.0.0 REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
target_link_libraries(tmocmd "${TORCH_LIBRARIES}")
set_property(TARGET tmocmd PROPERTY CXX_STANDARD 14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

INCLUDE_DIRECTORIES(
${OpenCV_INCLUDE_DIRS}
${TORCH_INCLUDE_DIRS}
../tmolib/
.
)

LINK_LIBRARIES (
${OpenCV_LIBS}
${TORCH_LIBRARIES}
tmo
)

SET(TMO_SOURCES
 TMOJah06.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOJah06.h
 TMOPlugin.h
)


add_library( jah06 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET jah06 POST_BUILD 
	COMMAND cp -f libjah06.so ../TMOCmd/jah06.tml
	COMMAND cp -f libjah06.so ../TMOgui/jah06.tml 
	COMMAND cp -f libjah06.so ../jah06.tml
	COMMAND cp -f ../../TMOJah06/model.pt ../model.pt
#	COMMAND rm -f libjah06.so 
	)


