cmake_minimum_required(VERSION 3.5)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(QT NAMES Qt6 Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)


set(gui_SRC
    #callGUI.cpp
    mainwindow.cpp
    mainwindow.ui
    ibf_mainwindow.cpp
    ibf_mainwindow.ui
    classify_mainwindow.cpp
    classify_mainwindow.ui
    live_deplete_mainwindow.cpp
    live_deplete_mainwindow.ui
    connection_test_mainwindow.cpp
    connection_test_mainwindow.ui
    )


qt_add_resources(gui_SRC pics.qrc)

set(gui_HDRS
    mainwindow.h
    ibf_mainwindow.h
    classify_mainwindow.h
    live_deplete_mainwindow.h
    connection_test_mainwindow.h
    QDebugStream.h
    )

file(GLOB gui_SRC "*.cpp")
file(GLOB gui_HDRS "*.hpp")

add_library (gui STATIC ${gui_SRC} ${gui_HDRS})

target_include_directories (gui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries (gui PUBLIC Qt${QT_VERSION_MAJOR}::Widgets)


install(TARGETS gui DESTINATION lib COMPONENT libraries)

install(FILES ${gui_HDRS} DESTINATION include COMPONENT headers)


find_package(Qt${QT_VERSION_MAJOR}Widgets)
