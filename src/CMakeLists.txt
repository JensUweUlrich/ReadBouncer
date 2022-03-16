
cmake_minimum_required (VERSION 3.16)

project (ReadBouncer CXX C)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)

set(ARM_BUILD "OFF")
SET(ARM_TARGETS "arm;aarch64")
list (FIND ARM_TARGETS ${CMAKE_SYSTEM_PROCESSOR} _index)
if (${_index} GREATER -1)
    set(ARM_BUILD "ON")
    add_compile_options(-DARM_BUILD)
endif()


if (WIN32)
    add_compile_options(-D_WIN32_WINNT=0x600)
    add_compile_options("/MT")
    add_compile_options("/std:c++17")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /MT")
    set(CMAKE_ASM_NASM_OBJECT_FORMAT win64)
    set(CMAKE_INSTALL_PREFIX "C:/${PROJECT_NAME}")
    find_file(libbz2_DLL "libbz2-1__.dll")
    install(FILES ${libbz2_DLL} DESTINATION bin COMPONENT applications)
else()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    set(CMAKE_INSTALL_PREFIX "/usr/local/${PROJECT_NAME}")
    add_compile_options("-Ofast")
#    add_link_options("stdc++fs")
    add_compile_options("-std=c++17")
endif()

# add all subdirectories
add_subdirectory(interfaces)
add_subdirectory(spdlog)
add_subdirectory(seqan)
add_subdirectory(IBF)
add_subdirectory(basecaller)
add_subdirectory(minknow)
add_subdirectory(lyra)
add_subdirectory(config)
add_subdirectory(main)
#add_subdirectory(test)

# now pack everything into an installer

set(CPACK_SOURCE_GENERATOR "NSIS")
set(CPACK_NSIS_EXECUTION_LEVEL "user")
set(CPACK_PACKAGE_NAME "ReadBouncer")
set(CPACK_PACKAGE_VENDOR "Jens-Uwe Ulrich")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ReadBouncer - Installer")
set(CPACK_PACKAGE_VERSION "1.1.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${PROJECT_NAME})

#set(CPACK_COMPONENTS_ALL applications libraries headers addfiles certificates)
set(CPACK_COMPONENTS_ALL applications addfiles certificates)

set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "ReadBouncer App")
#set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "Libraries")
#set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ Headers")
set(CPACK_COMPONENT_ADDFILES_DISPLAY_NAME "Additional Files")
set(CPACK_COMPONENT_CERTIFICATES_DISPLAY_NAME "Certificates")

set(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION
  "Application for real-time Oxford Nanopore Sequencing data analysis using MinKNOW API and ReadUntil")
#set(CPACK_COMPONENT_LIBRARIES_DESCRIPTION
#  "Static libraries used to build ReadBouncer")
#set(CPACK_COMPONENT_HEADERS_DESCRIPTION
#  "C/C++ header files for use with ReadBouncer")
set(CPACK_COMPONENT_ADDFILES_DESCRIPTION
  "Additional Files for Basecalling")
set(CPACK_COMPONENT_CERTIFICATES_DESCRIPTION
  "SSL/TLS Certificates")

#set(CPACK_COMPONENT_HEADERS_DEPENDS libraries)

#set(CPACK_COMPONENT_APPLICATIONS_GROUP "Runtime")
#set(CPACK_COMPONENT_ADDFILES_GROUP "Runtime")
#set(CPACK_COMPONENT_CERTIFICATES_GROUP "Runtime")
#set(CPACK_COMPONENT_LIBRARIES_GROUP "Development")
#set(CPACK_COMPONENT_HEADERS_GROUP "Development")

#set(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "All of the tools you'll ever need to develop software")

#set(CPACK_ALL_INSTALL_TYPES Full Developer)

#set(CPACK_COMPONENT_LIBRARIES_INSTALL_TYPES Developer Full)
#set(CPACK_COMPONENT_HEADERS_INSTALL_TYPES Developer Full)
#set(CPACK_COMPONENT_APPLICATIONS_INSTALL_TYPES Full)
#set(CPACK_COMPONENT_ADDFILES_INSTALL_TYPES Full)
#set(CPACK_COMPONENT_CERTIFICATES_INSTALL_TYPES Full)


# This must always be last!
include(CPack)
