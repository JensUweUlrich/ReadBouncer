
configure_file(CMakeLists.txt.in flappie-build/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flappie-build)
if(result)
  message(FATAL_ERROR "Generator step for flappie failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flappie-build)
if(result)
  message(FATAL_ERROR "Build step for flappie failed: ${result}")
endif()



