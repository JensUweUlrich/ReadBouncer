set(BLOOMSOURCES "bloomFilter.cpp" "bloomFilter.hpp")

add_library (BloomFilter STATIC ${BLOOMSOURCES})

target_link_libraries (BloomFilter seqan3::seqan3)

target_include_directories (BloomFilter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

