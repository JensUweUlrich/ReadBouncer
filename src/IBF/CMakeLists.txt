file(GLOB IBF_SRC "*.cpp")
file(GLOB IBF_HDRS "*.hpp")
file(GLOB WEIGHTS "${CMAKE_CURRENT_SOURCE_DIR}/../weights/*")

add_library (IBF STATIC ${IBF_SRC} ${IBF_HDRS})

if (WIN32)
	add_definitions(-DMSVC_COMPILER)
endif()

# include spd fast logging
set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../spdlog/spdlog/lib/cmake/spdlog)
find_package(spdlog REQUIRED)

# include SafeQueue and StopClock
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../util" "${CMAKE_CURRENT_SOURCE_DIR}/../interfaces")

# using seqan
target_include_directories (IBF PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/../interfaces")
target_link_libraries (IBF sdsl seqan spdlog::spdlog_header_only interfaces)

install(TARGETS IBF DESTINATION lib COMPONENT libraries)

install(FILES ${IBF_HDRS} DESTINATION include COMPONENT headers)
install(FILES ${WEIGHTS} DESTINATION data COMPONENT addfiles)
