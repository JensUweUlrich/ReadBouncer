# Download and unpack seqan at configure time
configure_file(CMakeLists.txt.in spdlog-build/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spdlog-build)
if(result)
  message(FATAL_ERROR "CMake step for spdlog failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spdlog-build)
if(result)
  message(FATAL_ERROR "Build step for spdlog failed: ${result}")
endif()

set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/spdlog/lib/cmake/spdlog)
find_package(spdlog REQUIRED)

