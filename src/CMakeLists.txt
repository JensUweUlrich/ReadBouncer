
cmake_minimum_required (VERSION 3.12)

add_compile_options(-DDEBUG)
set(CMAKE_GENERATORS "Unix Makefiles")
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")


project (RUHEL CXX C)
# add all subdirectories
add_subdirectory(seqan)
add_subdirectory(libReadUntil)
add_subdirectory(libMinimizer)
add_subdirectory(libBloomFilter)
add_subdirectory(main)
add_subdirectory(test)

