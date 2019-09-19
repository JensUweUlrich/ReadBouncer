set(BLOOMSOURCES "customBloomFilter.cpp" "customBloomFilter.hpp" "bloomFilterException.hpp" "bloom_filter.hpp")

add_library (BloomFilter STATIC ${BLOOMSOURCES})

target_link_libraries (BloomFilter seqan3::seqan3)

target_include_directories (BloomFilter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

