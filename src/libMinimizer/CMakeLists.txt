set(MINISOURCES "minimizer3.cpp" "minimizer3.hpp")

add_library (Minimizer STATIC ${MINISOURCES})

target_link_libraries (Minimizer seqan3::seqan3)

target_include_directories (Minimizer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

