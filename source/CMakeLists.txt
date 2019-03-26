cmake_minimum_required (VERSION 3.4)
project (MinHash_Containment CXX)
set(SeqAn3_DIR "/home/jens/Development/tutorial/seqan3/build_system")
find_package (SeqAn3 REQUIRED)
add_executable (mhc
                mhc.cpp)
target_link_libraries (mhc seqan3::seqan3)
