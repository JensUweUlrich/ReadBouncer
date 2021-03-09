file(GLOB IBF_SRC "*.cpp")
file(GLOB IBF_HDRS "*.hpp")

add_library (IBF STATIC ${IBF_SRC} ${IBF_HDRS})

add_definitions(-DMSVC_COMPILER)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../util")

# using seqan
target_include_directories (IBF PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries (IBF sdsl seqan)

install(TARGETS IBF DESTINATION lib COMPONENT libraries)

install(FILES ${IBF_HDRS} DESTINATION include COMPONENT headers)
