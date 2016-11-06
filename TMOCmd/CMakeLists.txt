CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

INCLUDE_DIRECTORIES(
../tmolib/
.
)

#Set up libraries to link with all found packages
LINK_LIBRARIES (
dl
tmo
pthread
${CMAKE_DL_LIBS}
)

SET(TMOCMD_SOURCES
 TMOCmd.cpp
 main.cpp
)


SET(TMOCMD_HEADERS
 TMOCmd.h
)




SET(EXE_NAME "tmocmd${TAG_VERSION}" )
ADD_EXECUTABLE(${EXE_NAME}
  ${TMOCMD_SOURCES}
)

add_custom_command( TARGET ${EXE_NAME} POST_BUILD
        COMMAND cp -r tmocmd ../
        )

