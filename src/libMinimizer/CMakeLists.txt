file(GLOB MINIMIZERSRCS "*.cpp")
file(GLOB MINIMIZERHDRS "*.hpp" "*.h")

add_library (Minimizer STATIC ${MINIMIZERSRCS} ${MINIMIZERHDRS})

target_link_libraries (Minimizer seqan3::seqan3)

target_include_directories (Minimizer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS Minimizer DESTINATION lib COMPONENT libraries)

install(FILES ${MINIMIZERHDRS} DESTINATION include COMPONENT headers)