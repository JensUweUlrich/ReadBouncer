
cmake_minimum_required (VERSION 3.4)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

project (ReadUntil CXX)

# Seqan specifics
set(SeqAn3_DIR "/home/jens/Development/seqan3/build_system")
find_package (SeqAn3 REQUIRED)

# add all subdirectories
add_subdirectory(libMinimizer)
add_subdirectory(libBloomFilter)
add_subdirectory(main)
add_subdirectory(test)
