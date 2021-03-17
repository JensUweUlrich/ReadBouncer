# Download and unpack seqan at configure time
configure_file(CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/lyra-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lyra-download)
if(result)
  message(FATAL_ERROR "CMake step for lyra failed: ${result}")
endif()

set(lyra_DIR ${CMAKE_CURRENT_BINARY_DIR}/lyra/share/lyra/cmake)

execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lyra-download )
if(result)
  message(FATAL_ERROR "Build step for lyra failed: ${result}")
endif()


find_package (lyra REQUIRED)

