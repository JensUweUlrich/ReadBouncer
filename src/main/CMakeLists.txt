set(Main_SRC "main.cpp")
file(GLOB MAIN_HDRS "*.hpp")
add_executable (${PROJECT_NAME} ${MAIN_HDRS} ${Main_SRC})
find_package (Threads)

#FIND_PATH(DEEPNANO_INCLUDE_DIR DeepNano2.h "${CMAKE_CURRENT_SOURCE_DIR}/../libBaseCaller" DOC "Deepnano2 includes")

# include spd fast logging
set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../spdlog/spdlog/lib/cmake/spdlog)
find_package(spdlog REQUIRED)

#include_directories(${DEEPNANO_INCLUDE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/../libBaseCaller" "${CMAKE_CURRENT_SOURCE_DIR}/../util")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../util" "${CMAKE_CURRENT_BINARY_DIR}/../libBaseCaller/osprey/src/src"
                    "${CMAKE_CURRENT_BINARY_DIR}/../libBaseCaller/include")
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../lyra/lyra/lyra-src/include)

#${CMAKE_CURRENT_BINARY_DIR}/../libBaseCaller/osprey/src/src 
#												  ${CMAKE_CURRENT_SOURCE_DIR}/../util ${CMAKE_CURRENT_BINARY_DIR}/../libBaseCaller/bc-build/basecaller-prefix/src/basecaller-build/Library/include 
#												  ${CMAKE_CURRENT_BINARY_DIR}/../libBaseCaller/bc-build/basecaller-prefix/src/basecaller-build/NumCpp/include
#												  ${CMAKE_CURRENT_BINARY_DIR}/../libBaseCaller/bc-build/basecaller-prefix/src/basecaller-build/boost_1_76_0)

# using seqan
if (WIN32)
	add_definitions(-DMSVC_COMPILER)
	#target_link_libraries (${PROJECT_NAME} IBF ReadUntil mklcore mklseq mkllp64 deepnano2 spdlog::spdlog_header_only)
	target_link_libraries (${PROJECT_NAME} IBF ReadUntil osprey spdlog::spdlog_header_only)
else()
	target_link_libraries (${PROJECT_NAME} IBF ReadUntil osprey Threads::Threads decoder mkllp64 mklcore mklseq mklcore spdlog::spdlog_header_only)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT applications)
