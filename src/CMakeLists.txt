
cmake_minimum_required (VERSION 3.16)

project (NanoLive CXX C)

#add_compile_options(-DDEBUG)
#add_compile_options(-DBUILD_SHARED_LIBS=ON)
#add_compile_options(-Wno-format)

add_compile_options(-D_WIN32_WINNT=0x600)
add_compile_options("/MT")
add_compile_options("/std:c++17")
#add_compile_options(-DMSVC_COMPILER)
#add_compile_options("/std:c++latest")
#add_definitions(-D_WIN32_WINNT=0x600)


#set(CMAKE_GENERATORS "MinGW Makefiles")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /MT")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
set(CMAKE_ASM_NASM_OBJECT_FORMAT win32)

if (WIN32)
    set(CMAKE_INSTALL_PREFIX "C:/Apps/${PROJECT_NAME}")
    find_file(libbz2_DLL "libbz2-1__.dll")
    install(FILES ${libbz2_DLL} DESTINATION bin COMPONENT applications)
else()
    set(CMAKE_INSTALL_PREFIX "/usr/local/${PROJECT_NAME}")
endif()

# add all subdirectories
add_subdirectory(spdlog)
add_subdirectory(seqan)
add_subdirectory(libIBF)
add_subdirectory(libBaseCaller)
add_subdirectory(libReadUntil)
add_subdirectory(lyra)
add_subdirectory(main)
#add_subdirectory(test)

# now pack everything into an installer

set(CPACK_SOURCE_GENERATOR "NSIS")
set(CPACK_NSIS_EXECUTION_LEVEL "user")
set(CPACK_PACKAGE_NAME "NanoLive")
set(CPACK_PACKAGE_VENDOR "Jens-Uwe Ulrich")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NanoLive - Installer")
set(CPACK_PACKAGE_VERSION "0.3.0")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "3")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${PROJECT_NAME})

set(CPACK_COMPONENTS_ALL applications libraries headers addfiles)

set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "NanoLive App")
set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "Libraries")
set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ Headers")
set(CPACK_COMPONENT_ADDFILES_DISPLAY_NAME "Additional Files")

set(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION
  "Application for real-time Oxford Nanopore Sequencing data analysis using MinKNOW API and ReadUntil")
set(CPACK_COMPONENT_LIBRARIES_DESCRIPTION
  "Static libraries used to build NanoLive")
set(CPACK_COMPONENT_HEADERS_DESCRIPTION
  "C/C++ header files for use with NanoLive")
set(CPACK_COMPONENT_ADDFILES_DESCRIPTION
  "Additional Files for Basecalling")

set(CPACK_COMPONENT_HEADERS_DEPENDS libraries)

set(CPACK_COMPONENT_APPLICATIONS_GROUP "Runtime")
set(CPACK_COMPONENT_ADDFILES_GROUP "Runtime")
set(CPACK_COMPONENT_LIBRARIES_GROUP "Development")
set(CPACK_COMPONENT_HEADERS_GROUP "Development")

set(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "All of the tools you'll ever need to develop software")

set(CPACK_ALL_INSTALL_TYPES Full Developer)

set(CPACK_COMPONENT_LIBRARIES_INSTALL_TYPES Developer Full)
set(CPACK_COMPONENT_HEADERS_INSTALL_TYPES Developer Full)
set(CPACK_COMPONENT_APPLICATIONS_INSTALL_TYPES Full)
set(CPACK_COMPONENT_ADDFILES_INSTALL_TYPES Full)



# This must always be last!
include(CPack)
