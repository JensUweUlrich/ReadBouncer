
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(QT NAMES Qt6 Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)

set(Main_SRC
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    ibf_mainwindow.cpp
    ibf_mainwindow.h
    ibf_mainwindow.ui
    classify_mainwindow.cpp
    classify_mainwindow.h
    classify_mainwindow.ui
    live_deplete_mainwindow.cpp
    live_deplete_mainwindow.h
    live_deplete_mainwindow.ui
    connection_test_mainwindow.cpp
    connection_test_mainwindow.h
    connection_test_mainwindow.ui
    QDebugStream.h
    )


qt_add_resources(Main_SRC pics.qrc)

#WIN32 is to delete the second momitor from the APP
file(GLOB MAIN_HDRS "*.hpp")
add_executable (${PROJECT_NAME} WIN32 ${MAIN_HDRS} ${Main_SRC})
#find_package (Threads)

FIND_PATH(DEEPNANO_INCLUDE_DIR DeepNano2.h "${CMAKE_CURRENT_SOURCE_DIR}/../libBaseCaller" DOC "Deepnano2 includes")

# include spd fast logging
set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../spdlog/spdlog/lib/cmake/spdlog)
find_package(spdlog REQUIRED)

include_directories(${DEEPNANO_INCLUDE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/../libBaseCaller" "${CMAKE_CURRENT_SOURCE_DIR}/../util")
#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../weights")

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../lyra/lyra/lyra-src/include ${CMAKE_CURRENT_SOURCE_DIR}/../libBaseCaller ${CMAKE_CURRENT_SOURCE_DIR}/../util)
# using seqan
add_definitions(-DMSVC_COMPILER)
target_link_libraries (${PROJECT_NAME} Qt${QT_VERSION_MAJOR}::Widgets)
target_link_libraries (${PROJECT_NAME} IBF ReadUntil deepnano2 spdlog::spdlog_header_only)
install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT applications)


#Qt

find_package(Qt${QT_VERSION_MAJOR}Widgets)
