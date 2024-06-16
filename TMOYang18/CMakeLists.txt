CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
find_package(PythonLibs REQUIRED)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

INCLUDE_DIRECTORIES(
../tmolib/
.
${PYTHON_INCLUDE_DIRS}
)



LINK_LIBRARIES (
tmo

${PYTHON_LIBRARIES}
)

SET(TMO_SOURCES
 TMOYang18.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOYang18.h
 TMOPlugin.h
)


add_library( Yang18 SHARED ${TMO_SOURCES} )

target_link_libraries(Yang18
    tensorflow  
    #tensorflow_framework
)

add_custom_command( TARGET Yang18 POST_BUILD 
	COMMAND cp -f libYang18.so ../TMOCmd/Yang18.tml
	COMMAND cp -f libYang18.so ../TMOgui/Yang18.tml 
	COMMAND cp -f libYang18.so ../Yang18.tml 
#	COMMAND rm -f libYang18.so 
	)


