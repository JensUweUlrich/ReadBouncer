
cmake_minimum_required (VERSION 3.16)

project (NanoLiveTk CXX C)

add_compile_options(-DDEBUG)
add_compile_options(-DBUILD_SHARED_LIBS=ON)
add_compile_options(-Wno-format)
add_compile_options(-D_WIN32_WINNT=0x600)
add_definitions(-D_WIN32_WINNT=0x600)

set(CMAKE_GENERATORS "MinGW Makefiles")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
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
add_subdirectory(seqan)
#add_subdirectory(libBaseCaller)
add_subdirectory(libReadUntil)
add_subdirectory(libMinimizer)
add_subdirectory(libBloomFilter)
add_subdirectory(main)
add_subdirectory(test)

# now pack everything into an installer

set(CPACK_SOURCE_GENERATOR "NSIS")
set(CPACK_NSIS_EXECUTION_LEVEL "user")
set(CPACK_PACKAGE_NAME "NanoLiveTk")
set(CPACK_PACKAGE_VENDOR "Jens-Uwe Ulrich")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NanoLiveTk - Installer")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "2")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${PROJECT_NAME})

set(CPACK_COMPONENTS_ALL applications libraries headers)

set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "NanoLiveTk App")
set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "Libraries")
set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ Headers")

set(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION
  "Application for real-time Oxford Nanopore Sequencing data analysis using MinKNOW API and ReadUntil")
set(CPACK_COMPONENT_LIBRARIES_DESCRIPTION
  "Static libraries used to build NanoLiveTk")
set(CPACK_COMPONENT_HEADERS_DESCRIPTION
  "C/C++ header files for use with NanoLiveTk")

set(CPACK_COMPONENT_HEADERS_DEPENDS libraries)

set(CPACK_COMPONENT_APPLICATIONS_GROUP "Runtime")
set(CPACK_COMPONENT_LIBRARIES_GROUP "Development")
set(CPACK_COMPONENT_HEADERS_GROUP "Development")

set(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "All of the tools you'll ever need to develop software")

set(CPACK_ALL_INSTALL_TYPES Full Developer)

set(CPACK_COMPONENT_LIBRARIES_INSTALL_TYPES Developer Full)
set(CPACK_COMPONENT_HEADERS_INSTALL_TYPES Developer Full)
set(CPACK_COMPONENT_APPLICATIONS_INSTALL_TYPES Full)



# This must always be last!
include(CPack)
