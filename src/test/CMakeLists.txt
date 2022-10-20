# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()


# Add googletest directly to our build. This defines
# the gtest, gtest_main, gmock and gmock_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_toml_file.toml
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

 # Use this command only if you want to change the toml file during the tests, otherwise keep the previous one with "file()"
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_toml_file.toml
#    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# include SafeQueue and StopClock
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../util" "${CMAKE_CURRENT_SOURCE_DIR}/../interfaces")
# add test directories for own tests
add_subdirectory(libConfigReaderTests)
#add_subdirectory(libIBFTests)
#add_subdirectory(classifyTests)
#add_subdirectory(MinKNOWTests)


