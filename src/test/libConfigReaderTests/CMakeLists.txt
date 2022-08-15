# Now simply link against gtest or gtest_main as needed. Eg
add_executable(ConfigReaderTest ConfigReaderTest.cpp)
target_link_libraries(ConfigReaderTest config gmock_main)
add_test(NAME ConfigReaderTest COMMAND ConfigReaderTest)
