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
add_subdirectory(TMOAlsam09/)							#Martin Molek
add_subdirectory(TMOAncuti16/)							#Vladimir Vlkovic
#add_subdirectory(TMOYourOperatorNameLinux/)
#add_subdirectory(TMOJin14/)
add_subdirectory(TMOLee10/)							#Martin Molek
add_subdirectory(TMOKim09/)							#Petr Pospisil: method 1
#add_subdirectory(TMOZheng15/)							#Petr Pospisil: method 2
#add_subdirectory(TMOBae06/)							#Petr Pospisil: method 3 (problem with gui, gui doesnt support TMOString parameter)
add_subdirectory(TMOBiswas05/)							#Tomas Chlubna
add_subdirectory(TMOCadik08/)							#Martin Cadik/Jan Brida: OpenCL implementation (CPU and GPU versions)
#add_subdirectory(TMOXiong17/)							#Filip Brezna: not yet ready 
add_subdirectory(TMOHu14/)							#Vladimir Vlkovic
#add_subdirectory(TMOTao17/)							#Vladimir Vlkovic not ready

INSTALL(FILES tmolib/libtmo.so DESTINATION lib)
INSTALL(FILES TMOgui/tmogui PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION . )
INSTALL(FILES TMOCmd/tmocmd PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION . )
INSTALL(FILES TMOCmd/ward94.tml DESTINATION . )

