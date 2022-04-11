
add_executable(classifyTest classifytest.cpp)
target_link_libraries(classifyTest sdsl seqan IBF interfaces gmock_main)
install(TARGETS classifyTest DESTINATION bin COMPONENT applications)