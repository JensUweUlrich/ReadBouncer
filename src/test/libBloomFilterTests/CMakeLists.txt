# Now simply link against gtest or gtest_main as needed. Eg
add_executable(BloomFilterTest BloomFilterTest.cpp)
target_link_libraries(BloomFilterTest gmock_main)
add_test(NAME BloomFilterTest COMMAND BloomFilterTest)
