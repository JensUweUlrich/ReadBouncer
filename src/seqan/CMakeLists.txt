# Download and unpack seqan at configure time
configure_file(CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download)
if(result)
  message(FATAL_ERROR "CMake step for seqan3 failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download )
if(result)
  message(FATAL_ERROR "Build step for seqan3 failed: ${result}")
endif()

set(SeqAn3_DIR ${CMAKE_CURRENT_BINARY_DIR}/seqan3/seqan3-src/build_system)
find_package (SeqAn3 REQUIRED)

