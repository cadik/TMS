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
 TMOWanat14.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOWanat14.h
 TMOPlugin.h
)

find_package(NLopt REQUIRED)

# Ak je nlopt nájdené, pridá jeho cesty k include directories
if(NLopt_FOUND)
  include_directories(${NLOPT_INCLUDE_DIRS})
endif()

# Vytvorenie zdieľanej knižnice
add_library(wanat14 SHARED ${TMO_SOURCES})

# Linkovanie s nlopt a ostatnými potrebnými knižnicami
if(NLopt_FOUND)
  target_link_libraries(wanat14 ${NLOPT_LIBRARIES})
endif()


add_library( wanat14 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET wanat14 POST_BUILD 
	COMMAND cp -f libwanat14.so ../TMOCmd/wanat14.tml
	COMMAND cp -f libwanat14.so ../TMOgui/wanat14.tml 
	COMMAND cp -f libwanat14.so ../wanat14.tml 
#	COMMAND rm -f libwanat14.so 
	)


