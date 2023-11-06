CMAKE_MINIMUM_REQUIRED(VERSION 3.18)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

include(FetchContent)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

INCLUDE_DIRECTORIES(
../tmolib/
.
)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

##########
#Tensorflow
##########
find_package(TensorFlow REQUIRED)
include_directories(${TensorFlow_INCLUDE_DIRS})
add_definitions(${TensorFlow_DEFINITIONS})

# download model
function(download_model url)

FetchContent_Declare(download_pretrained_model
URL ${url}
DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/model
DOWNLOAD_NO_EXTRACT true
TLS_VERIFY true
)

FetchContent_Populate(download_pretrained_model)

endfunction(download_model)

LINK_LIBRARIES (
	tmo
	${TensorFlow_LIBRARIES}
)

message(STATUS "Downloading the pre-trained model for the TMOXia18 operator...")

download_model(
https://github.com/cadik/TMS-storage/raw/main/TMOXia18/tmoxia18_model.zip
)

message(STATUS "Extracting the downloaded model ...")

file(ARCHIVE_EXTRACT
INPUT ${CMAKE_CURRENT_SOURCE_DIR}/model/tmoxia18_model.zip
DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/model
)

SET(TMO_SOURCES
 TMOXia18.cpp
 TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
 TMOXia18.h
 TMOPlugin.h
)


add_library( xia18 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET xia18 POST_BUILD 
	COMMAND cp -f libxia18.so ../TMOCmd/xia18.tml
	COMMAND cp -f libxia18.so ../TMOgui/xia18.tml 
	COMMAND cp -f libxia18.so ../xia18.tml 
#	COMMAND rm -f libxia18.so
	)

