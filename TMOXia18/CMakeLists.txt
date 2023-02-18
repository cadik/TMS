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
${CMAKE_CURRENT_SOURCE_DIR}/include/
)

# download libtensorflow
function(download_lib url)

FetchContent_Declare(download
URL ${url}
DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib
DOWNLOAD_NO_EXTRACT true
TLS_VERIFY true
)

FetchContent_Populate(download)

endfunction(download_lib)

# download libtensorflow_framework
function(download_framework url)

FetchContent_Declare(download_fw
URL ${url}
DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib
DOWNLOAD_NO_EXTRACT true
TLS_VERIFY true
)

FetchContent_Populate(download_fw)

endfunction(download_framework)

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

IF (WIN32)
  # set stuff for windows
  
  	message(STATUS "Downloading the builded TensorFlow library for the TMOXia18 operator ...")
    
  	download_lib(
	https://github.com/davidchocholaty/TMOXia18-libs/raw/windows/libtensorflow_windows.zip
	)

	message(STATUS "Extracting the downloaded library ...")

	file(ARCHIVE_EXTRACT
	INPUT ${CMAKE_CURRENT_SOURCE_DIR}/lib/libtensorflow_windows.zip
	DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
	)

  	LINK_LIBRARIES (
	tmo
	${CMAKE_CURRENT_SOURCE_DIR}/lib/libtensorflow_cc.dll
	${CMAKE_CURRENT_SOURCE_DIR}/lib/libtensorflow_cc.lib
	)
	
	message(STATUS "Downloading the pre-trained model for the TMOXia18 operator ...")
		
	download_model(
	https://github.com/davidchocholaty/TMOXia18-model/raw/main/tmoxia18_model.zip
	)

	message(STATUS "Extracting the downloaded model ...")
	
	file(ARCHIVE_EXTRACT
	INPUT ${CMAKE_CURRENT_SOURCE_DIR}/model/tmoxia18_model.zip
	DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/model
	)

ELSE()
  	# set stuff for other systems
  	message(STATUS "Downloading the builded TensorFlow library for the TMOXia18 operator ...")
  
  	download_lib(
	https://github.com/davidchocholaty/TMOXia18-libs/raw/linux/libtensorflow_cc_linux.tar.gz
	)

  	message(STATUS "Extracting the downloaded library ...")

	file(ARCHIVE_EXTRACT
	INPUT ${CMAKE_CURRENT_SOURCE_DIR}/lib/libtensorflow_cc_linux.tar.gz
	DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
	)
	
  	message(STATUS "Downloading the builded TensorFlow framework for the TMOXia18 operator ...")
	  
	download_framework(
	https://github.com/davidchocholaty/TMOXia18-libs/raw/linux/libtensorflow_framework_linux.tar.gz
	)

  	message(STATUS "Extracting the downloaded framework ...")

	file(ARCHIVE_EXTRACT
	INPUT ${CMAKE_CURRENT_SOURCE_DIR}/lib/libtensorflow_framework_linux.tar.gz
	DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
	)
	  
  	LINK_LIBRARIES (
	tmo
	${CMAKE_CURRENT_SOURCE_DIR}/lib/libtensorflow_cc.so
	)
	
	message(STATUS "Downloading the pre-trained model for the TMOXia18 operator ...")
	  
	download_model(
	https://github.com/davidchocholaty/TMOXia18-model/raw/main/tmoxia18_model.zip
	)
	
	message(STATUS "Extracting the downloaded model ...")

	file(ARCHIVE_EXTRACT
	INPUT ${CMAKE_CURRENT_SOURCE_DIR}/model/tmoxia18_model.zip
	DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/model
	)
	
ENDIF()

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


