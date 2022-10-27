# Now simply link against gtest or gtest_main as needed.
add_executable(ConfigReaderTest ConfigReaderTest.cpp)
target_link_libraries(ConfigReaderTest config gmock_main)
add_test(NAME ConfigReaderTest COMMAND ConfigReaderTest)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_toml_file.toml ${CMAKE_CURRENT_SOURCE_DIR}/deplete_test1.fasta
          ${CMAKE_CURRENT_SOURCE_DIR}/deplete_test2.fasta ${CMAKE_CURRENT_SOURCE_DIR}/target_test1.fasta
          ${CMAKE_CURRENT_SOURCE_DIR}/target_test2.fasta ${CMAKE_CURRENT_SOURCE_DIR}/reads_test.fastq
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

 # Use this command only if you want to change the toml file during the tests, otherwise keep the previous one with "file()"
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_toml_file.toml
#    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)