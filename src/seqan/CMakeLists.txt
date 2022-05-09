# Download and unpack seqan at configure time

if (WIN32)
  set(sdsl_cache_args  "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite
                        -DMSVC_COMPILER")
else()
  set(sdsl_cache_args  "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite")
endif()

message(STATUS "${CMAKE_GENERATOR}")

#configure_file(CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download/CMakeLists.txt)
configure_file(CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/seqan2-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download)
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan2-download)
if(result)
  message(FATAL_ERROR "CMake step for seqan failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release
  RESULT_VARIABLE result
  #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download )
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan2-download)
#if(result)
#  message(FATAL_ERROR "Build step for seqan failed: ${result}")
#endif()

add_library(sdsl STATIC IMPORTED GLOBAL)
if (WIN32)
  set_target_properties(sdsl PROPERTIES 
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/seqan2-download/sdsl-lite-prefix/src/sdsl-lite-build/lib/Release/sdsl.lib")
else()
  set_target_properties(sdsl PROPERTIES 
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/seqan2-download/sdsl-lite-prefix/src/sdsl-lite-build/lib/libsdsl.a")
endif()
target_include_directories( sdsl SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/include )

#set(SeqAn3_DIR ${CMAKE_CURRENT_BINARY_DIR}/seqan3/seqan3-src/build_system)
#find_package (SeqAn3 REQUIRED)
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR}/seqan2/seqan2-src/util/cmake)
set(SEQAN_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/seqan2/seqan2-src/include)
set(SeqAn_DIR ${CMAKE_CURRENT_BINARY_DIR}/seqan2/seqan2-src/util/cmake)
#find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------


find_package( SeqAn CONFIG REQUIRED PATHS ${SeqAn_DIR} )

add_library( seqan INTERFACE )
target_include_directories( seqan SYSTEM INTERFACE ${SEQAN_INCLUDE_DIRS} )
target_link_libraries( seqan INTERFACE ${SEQAN_LIBRARIES} gui)

# Add include directories.

#message(FATAL_ERROR ${SEQAN_INCLUDE_DIRS})
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
#add_executable (my_project my_project.cpp)
#target_link_libraries (my_project ${SEQAN_LIBRARIES})

