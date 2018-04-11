CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -fPIC")
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)


#pallas
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EIGEN_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/eigen/")
set(GTEST_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/googletest/googletest")
set(RAPIDJSON_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/rapidjson")
set(FIND_GLOG_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${FIND_GLOG_SUBDIR};${CMAKE_MODULE_PATH}")

set(EIGEN_INCLUDE_DIR ${EIGEN_SUBDIR})
##########
#Eigen
##########
#set(EIGEN_DIR /usr/local/include/eigen3/)
#include_directories(${EIGEN_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})

find_package(Ceres REQUIRED)
find_package(Glog REQUIRED)

include_directories(APPEND
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/
                    ${CERES_INCLUDE_DIRS}
                    ${EIGEN_SUBDIR}
                    ${GTEST_SUBDIR}
                    ${GLOG_INCLUDE_DIR}
                    ${RAPIDJSON_SUBDIR}/include)

add_subdirectory(${GTEST_SUBDIR})
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/pallas)