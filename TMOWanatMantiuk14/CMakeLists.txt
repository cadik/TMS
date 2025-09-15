CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
 TMOWanatMantiuk14.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOWanatMantiuk14.h
 TMOPlugin.h
)

find_package(NLopt REQUIRED)

# Ak je nlopt nájdené, pridá jeho cesty k include directories
if(NLopt_FOUND)
  include_directories(${NLOPT_INCLUDE_DIRS})
endif()

# Vytvorenie zdieľanej knižnice
add_library(wanatMantiuk14 SHARED ${TMO_SOURCES})

# Linkovanie s nlopt a ostatnými potrebnými knižnicami
if(NLopt_FOUND)
  target_link_libraries(wanatMantiuk14 ${NLOPT_LIBRARIES})
endif()


add_library( wanatMantiuk14 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET wanatMantiuk14 POST_BUILD 
	COMMAND cp -f libwanatMantiuk14.so ../TMOCmd/wanatMantiuk14.tml
	COMMAND cp -f libwanatMantiuk14.so ../TMOgui/wanatMantiuk14.tml 
	COMMAND cp -f libwanatMantiuk14.so ../wanatMantiuk14.tml 
#	COMMAND rm -f libwanatMantiuk14.so 
	)


