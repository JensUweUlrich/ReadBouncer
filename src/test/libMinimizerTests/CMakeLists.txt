# Now simply link against gtest or gtest_main as needed. Eg
add_executable(MinimizerTest minimizerTest.cpp)
target_link_libraries(MinimizerTest seqan3::seqan3 Minimizer gmock_main)
add_test(NAME MinimizerTest COMMAND MinimizerTest)
