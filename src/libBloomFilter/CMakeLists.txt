
file(GLOB BLOOMSRCS "*.cpp")
file(GLOB BLOOMHDRS "*.hpp")

add_library (BloomFilter STATIC ${BLOOMSRCS} ${BLOOMHDRS})

target_link_libraries (BloomFilter seqan3::seqan3)

target_include_directories (BloomFilter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS BloomFilter DESTINATION lib COMPONENT libraries)

install(FILES ${BLOOMHDRS} DESTINATION include COMPONENT headers)