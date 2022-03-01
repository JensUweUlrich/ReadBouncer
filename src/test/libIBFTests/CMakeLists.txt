
add_executable(IBFTest IBFTest.cpp)
target_link_libraries(IBFTest sdsl seqan IBF interfaces gmock_main)
add_test(NAME IBFTest COMMAND IBFTest)
