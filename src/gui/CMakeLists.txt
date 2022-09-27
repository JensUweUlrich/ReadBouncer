
cmake_minimum_required(VERSION 3.5)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(QT NAMES Qt6 Qt5 COMPONENTS Widgets Qml Quick Gui PrintSupport REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets Qml Quick Gui PrintSupport REQUIRED)

IF (${QT_FOUND})
   MESSAGE(STATUS "Found " Qt${QT_VERSION_MAJOR})
ELSE (${QT_FOUND})
   MESSAGE(FATAL_ERROR "Could not locate Qt.")
ENDIF (${QT_FOUND})

set(gui_SRC )


#qt_add_resources(gui_SRC pics.qrc)
#qt5_add_resources(gui_SRC pics.qrc)
# add the qml.qrc file
qt5_add_resources(qml_QRC qml.qrc)


set(gui_SRC
    #callGUI.cpp
    mainwindow.cpp
    mainwindow.ui
    ibf_mainwindow.cpp
    ibf_mainwindow.ui
    classify_mainwindow.cpp
    classify_mainwindow.ui
    qcustomplot.cpp
    #live_deplete_mainwindow.cpp
    #live_deplete_mainwindow.ui
    #connection_test_mainwindow.cpp
    #connection_test_mainwindow.ui
 #   pics.qrc
    )


#qt_add_resources(gui_SRC pics.qrc)
#qt5_add_resources(gui_SRC pics.qrc)

set(gui_HDRS
    mainwindow.h
    ibf_mainwindow.h
    classify_mainwindow.h
    qcustomplot.h
    #live_deplete_mainwindow.h
    #connection_test_mainwindow.h
    #QDebugStream.h
    )
#qt_add_resources(gui_QRC pics.qrc)
file(GLOB gui_SRC "*.cpp")
file(GLOB gui_HDRS "*.h")
file(GLOB gui_UI "*.ui")
file(GLOB gui_QRC "*.qrc")

add_library (gui STATIC ${gui_SRC} ${gui_HDRS} ${gui_UI} ${gui_QRC})

target_include_directories (gui  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries (gui PUBLIC Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Gui Qt${QT_VERSION_MAJOR}::Quick Qt${QT_VERSION_MAJOR}::Qml Qt${QT_VERSION_MAJOR}::PrintSupport)

install(TARGETS gui DESTINATION lib COMPONENT libraries)

install(FILES ${gui_HDRS}  DESTINATION include COMPONENT headers)


find_package(Qt${QT_VERSION_MAJOR}Widgets)
