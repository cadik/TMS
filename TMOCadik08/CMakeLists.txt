CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
	ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

FIND_PACKAGE(OpenCL REQUIRED)
IF (OpenCL_FOUND)
	MESSAGE("OpenCL found OK")
ELSE(OpenCL_FOUND)
	MESSAGE(FATAL_ERROR "OpenCL library not found.")
ENDIF(OpenCL_FOUND)

INCLUDE_DIRECTORIES(
	../tmolib/
	.
)

LINK_LIBRARIES (
	${OpenCL_LIBRARY}
	${OpenCL_LIBRARIES}
	${OpenCL_3RD_PARTY_LIBRARIES}
	tmo
)

SET(TMO_SOURCES
	TMOCadik08.cpp
	TMOPlugin.cpp
	coloroid.cpp
	spacefill.cpp
	quadtree.cpp
	morton.cpp
)

SET(TMOWARD_HEADERS
	TMOCadik08.h
	TMOPlugin.h
	coloroid.h
	spacefill.h
	quadtree.h
	morton.h
)

add_library( Cadik08 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET Cadik08 POST_BUILD
	COMMAND cp -f libCadik08.so ../TMOCmd/Cadik08.tml
	COMMAND cp -f libCadik08.so ../TMOgui/Cadik08.tml 
	COMMAND cp -f libCadik08.so ../Cadik08.tml
#	COMMAND rm -f libCadik08.so
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/*.txt ${CMAKE_CURRENT_BINARY_DIR}/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/*.TXT ${CMAKE_CURRENT_BINARY_DIR}/

)

#configure_file(COLOR_GRAY_UNIF.txt COLOR_GRAY_UNIF.txt COPYONLY)
#configure_file(resources/kernels/color2gray.cl resources/kernels/color2gray.cl COPYONLY)



# Scan through resource folder for updated files and copy if none existing or changed
file (GLOB_RECURSE resources "resources/*.*")
foreach(resource ${resources})
 get_filename_component(filename ${resource} NAME)
 get_filename_component(dir ${resource} DIRECTORY)
 get_filename_component(dirname ${dir} NAME)
 
 set (output "")
 
 while(NOT ${dirname} STREQUAL resources)
  get_filename_component(path_component ${dir} NAME)
  set (output "${path_component}/${output}")
  get_filename_component(dir ${dir} DIRECTORY)
  get_filename_component(dirname ${dir} NAME)
 endwhile()
 
 set(output "${CMAKE_CURRENT_BINARY_DIR}/resources/${output}/${filename}")
 
 add_custom_command(
  COMMENT "Moving updated resource-file '${filename}'"
  OUTPUT ${output}
  DEPENDS ${resource}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  ${resource}
  ${output}
 )
 add_custom_target(${filename} ALL DEPENDS ${resource} ${output})
 
endforeach()

