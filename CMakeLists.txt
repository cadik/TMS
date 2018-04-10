CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX)
 SET(LINUX_SET TRUE)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX)

add_subdirectory(tmolib/)
add_subdirectory(TMOgui/)
add_subdirectory(TMOCmd/)
add_subdirectory(TMOWard94/)
add_subdirectory(TMOZhao10/)
add_subdirectory(TMOMM16/)
add_subdirectory(TMOAlsam09/)
add_subdirectory(TMOAncuti16/)
#add_subdirectory(TMOYourOperatorNameLinux/)
#add_subdirectory(TMOJin14/)
add_subdirectory(TMOLee10/)
add_subdirectory(TMOKim09/)							# method 1
#add_subdirectory(TMOZheng15/)							# method 2
#add_subdirectory(TMOBae06/)							# method 3 (problem with gui, gui doesnt support TMOString parametr)
add_subdirectory(TMOCadik08/)
add_subdirectory(TMOSon14/)	

INSTALL(FILES tmolib/libtmo.so DESTINATION lib)
INSTALL(FILES TMOgui/tmogui PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION . )
INSTALL(FILES TMOCmd/tmocmd PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION . )
INSTALL(FILES TMOCmd/ward94.tml DESTINATION . )

#pallas

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EIGEN_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/pallas/ext/eigen/")
set(GTEST_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/pallas/ext/googletest/googletest")
set(RAPIDJSON_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/pallas/ext/rapidjson")
set(FIND_GLOG_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/pallas/cmake")
set(CMAKE_MODULE_PATH "${FIND_GLOG_SUBDIR};${CMAKE_MODULE_PATH}")

#set(EIGEN_INCLUDE_DIR ${EIGEN_SUBDIR})
##########
#Eigen
##########
#set(EIGEN_DIR /usr/local/include/eigen3/)
#include_directories(${EIGEN_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/pallas/cmake)
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})

find_package(Ceres REQUIRED)
find_package(Glog REQUIRED)

include_directories(APPEND
                    ${CMAKE_CURRENT_SOURCE_DIR}/pallas/include/
                    ${CMAKE_CURRENT_SOURCE_DIR}/pallas/src/
                    ${CERES_INCLUDE_DIRS}
                    ${EIGEN_SUBDIR}
                    ${GTEST_SUBDIR}
                    ${GLOG_INCLUDE_DIR}
                    ${RAPIDJSON_SUBDIR}/pallas/include)

add_subdirectory(${GTEST_SUBDIR})
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pallas/src/pallas)
