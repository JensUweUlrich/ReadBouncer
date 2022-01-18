file(GLOB config_SRC "*.cpp")
file(GLOB config_HDRS "*.hpp")

add_library (config STATIC ${config_SRC} ${config_HDRS})

if (WIN32)
	add_definitions(-DMSVC_COMPILER)
endif()

# include spd fast logging
set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../spdlog/spdlog/lib/cmake/spdlog)
find_package(spdlog REQUIRED)

# include SafeQueue and StopClock
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../util")

# using seqan
target_include_directories (config PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries (config IBF sdsl seqan spdlog::spdlog_header_only)

install(TARGETS config DESTINATION lib COMPONENT libraries)

install(FILES ${config_HDRS} DESTINATION include COMPONENT headers)
