
cmake_minimum_required (VERSION 3.4)

project (ReadUntil)

# Seqan specifics
set(SeqAn3_DIR "/home/jens/Development/seqan3/build_system")
find_package (SeqAn3 REQUIRED)

#include_directories(libMinimizer libBloomFilter)
#link_directories(libMinimizer libBloomFilter)

# add all subdirectories
add_subdirectory(libMinimizer)
add_subdirectory(libBloomFilter)
add_subdirectory(main)
#ADD_SUBDIRECTORY(test)
