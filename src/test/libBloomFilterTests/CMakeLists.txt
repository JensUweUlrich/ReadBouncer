# Now simply link against gtest or gtest_main as needed. Eg
add_executable(BloomFilterTest bloomFilterTest.cpp)
target_link_libraries(BloomFilterTest seqan3::seqan3 BloomFilter gmock_main)
add_test(NAME BloomFilterTest COMMAND BloomFilterTest)
