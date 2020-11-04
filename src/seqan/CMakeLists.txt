# Download and unpack seqan at configure time
#configure_file(CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download/CMakeLists.txt)
configure_file(CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/seqan2-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download)
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan2-download)
if(result)
  message(FATAL_ERROR "CMake step for seqan failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download )
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan2-download)
if(result)
  message(FATAL_ERROR "Build step for seqan failed: ${result}")
endif()

#set(SeqAn3_DIR ${CMAKE_CURRENT_BINARY_DIR}/seqan3/seqan3-src/build_system)
#find_package (SeqAn3 REQUIRED)

set(SeqAn_DIR ${CMAKE_CURRENT_BINARY_DIR}/seqan2/seqan2-src/util/cmake)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
#add_executable (my_project my_project.cpp)
#target_link_libraries (my_project ${SEQAN_LIBRARIES})

