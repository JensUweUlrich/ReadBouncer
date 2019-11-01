
cmake_minimum_required (VERSION 3.4)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

project (RUHEL CXX)

# add all subdirectories
add_subdirectory(seqan)
add_subdirectory(libReadUntil)
add_subdirectory(libMinimizer)
add_subdirectory(libBloomFilter)
add_subdirectory(main)
add_subdirectory(test)

