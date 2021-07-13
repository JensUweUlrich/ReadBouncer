


set(Main_SRC "main.cpp")



file(GLOB MAIN_HDRS "*.hpp")
add_executable (${PROJECT_NAME} ${MAIN_HDRS} ${Main_SRC})

FIND_PATH(DEEPNANO_INCLUDE_DIR DeepNano2.h "${CMAKE_CURRENT_SOURCE_DIR}/../libBaseCaller" DOC "Deepnano2 includes")

# include spd fast logging
set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../spdlog/spdlog/lib/cmake/spdlog)
find_package(spdlog REQUIRED)

include_directories(${DEEPNANO_INCLUDE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/../libBaseCaller" "${CMAKE_CURRENT_SOURCE_DIR}/../util")

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../lyra/lyra/lyra-src/include ${CMAKE_CURRENT_SOURCE_DIR}/../libBaseCaller ${CMAKE_CURRENT_SOURCE_DIR}/../util)

# using seqan
if (WIN32)
        add_definitions(-DMSVC_COMPILER)
        target_link_libraries (${PROJECT_NAME} gui IBF  ReadUntil  mklcore  mklseq  mkllp64  deepnano2  spdlog::spdlog_header_only)


else()
        target_link_libraries (${PROJECT_NAME} gui IBF   ReadUntil  deepnano2  mkllp64  mklseq  mklcore  spdlog::spdlog_header_only)


endif()


#target_link_libraries (${PROJECT_NAME}  Qt6::Widgets)#caused error ???????????????

install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT applications)


#find_package(Qt${QT_VERSION_MAJOR}Widgets)
#target_link_libraries (${PROJECT_NAME}  Qt6::Widgets)
