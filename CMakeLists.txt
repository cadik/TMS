CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)

project (TMS DESCRIPTION "Tone Mapping Studio")

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX)
 SET(LINUX_SET TRUE)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX)


#Dokumentace se vygeneruje prikazem make doc
find_package(Doxygen)
if(DOXYGEN_FOUND)
	add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
endif(DOXYGEN_FOUND)

add_subdirectory(tmolib/)
add_subdirectory(TMOgui/)
add_subdirectory(TMOCmd/)
add_subdirectory(TMOWard94/)
add_subdirectory(TMOZhao10/)
#add_subdirectory(TMOMM16/)               #Martin Molek
#add_subdirectory(TMOAlsam06/)            #Jiri Zverina
add_subdirectory(TMOAlsam09/)             #Martin Molek
add_subdirectory(TMOAncuti16/)            #Vladimir Vlkovic
#add_subdirectory(TMOCheryl11/)           #Vaclav Chvila
#add_subdirectory(TMOFarbman08/)          #Tomas Hudziec
add_subdirectory(TMOAubry14/)             #Tomas Hudziec
add_subdirectory(TMOKou15/)               #Tomas Hudziec
add_subdirectory(TMOMeylan06/)            
add_subdirectory(TMODrago03/)             #David Chocholaty
add_subdirectory(TMOParis11/)             #David Chocholaty
add_subdirectory(TMOGastal11/)            #David Chocholaty
add_subdirectory(TMOShan10/)              #David Chocholaty
#add_subdirectory(TMOXia18/)              #David Chocholaty
add_subdirectory(TMOYee03/)               #Matus Bicanovsky
#add_subdirectory(TMOYourOperatorNameLinux/)

#add_subdirectory(TMOJin14/)               #
add_subdirectory(TMOLee10/)               #Martin Molek
add_subdirectory(TMOKim09/)               #Petr Pospisil: method 1/Vladimir Vlkovic
add_subdirectory(TMOZheng15/)             #Petr Pospisil: method 2
#add_subdirectory(TMOBae06/)              #Petr Pospisil: method 3 (problem with gui, gui doesnt support TMOString parameter)
add_subdirectory(TMOBiswas05/)            #Tomas Chlubna
#add_subdirectory(TMOCadik08/)            #Martin Cadik/Jan Brida: OpenCL implementation (CPU and GPU versions)
add_subdirectory(TMOXiong17/)             #Filip Brezna: method 1
add_subdirectory(TMOZhongping15/)         #Filip Brezna: method 2
#add_subdirectory(TMOJin17/)              #Filip Brezna: method 3
#add_subdirectory(TMOHu14/)               #Vladimir Vlkovic
#add_subdirectory(TMOTao17/)              #Vladimir Vlkovic not ready
#add_subdirectory(TMOSon14/)              #Sedlar Pavel, Tomas Hudziec
add_subdirectory(TMOCLAHE/)               #Sedlar Pavel
#add_subdirectory(TMOKuhn08/)             #Pavel Melcer
add_subdirectory(TMOEisemann04/)          #Roman Jaska
add_subdirectory(TMOTai08/)               #Michal Vlnas
#add_subdirectory(TMONguy15/)             #Matej Valek
add_subdirectory(TMOWu24/)                #Filip Sapak
add_subdirectory(TMOWu24_2/)							#Branislav Dubec
add_subdirectory(TMONafchi17)							#Tomas Krivanek

INSTALL(FILES tmolib/libtmo.so DESTINATION lib)
INSTALL(FILES TMOgui/tmogui PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION . )
INSTALL(FILES TMOCmd/tmocmd PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION . )
INSTALL(FILES TMOCmd/ward94.tml DESTINATION . )
