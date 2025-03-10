CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

INCLUDE_DIRECTORIES(
        ../tmolib/
        .
)

SET(TMO_SOURCES
        TMOKhan20.cpp
        TMOPlugin.cpp
)


SET(TMOWARD_HEADERS
        TMOKhan20.h
        TMOPlugin.h
)



add_library(Khan20 SHARED ${TMO_SOURCES})
target_link_libraries(
        Khan20 PRIVATE tmo
)

add_custom_command( TARGET Khan20 POST_BUILD
        COMMAND cp -f libKhan20.so ../TMOgui/Khan20.tml
        COMMAND cp -f libKhan20.so ../Khan20.tml
)
