# Now simply link against gtest or gtest_main as needed.

add_executable(ConfigReaderTest ConfigReaderTest.cpp)
target_link_libraries(ConfigReaderTest config gmock_main)
add_test(NAME ConfigReaderTest COMMAND ConfigReaderTest)


if (WIN32)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/testData/test_toml_file.toml ${CMAKE_CURRENT_SOURCE_DIR}/testData/deplete_test1.fasta
          ${CMAKE_CURRENT_SOURCE_DIR}/testData/deplete_test2.fasta ${CMAKE_CURRENT_SOURCE_DIR}/testData/target_test1.fasta
          ${CMAKE_CURRENT_SOURCE_DIR}/testData/target_test2.fasta ${CMAKE_CURRENT_SOURCE_DIR}/testData/reads_test.fastq
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Release)
else()
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/testData/test_toml_file.toml ${CMAKE_CURRENT_SOURCE_DIR}/testData/deplete_test1.fasta
          ${CMAKE_CURRENT_SOURCE_DIR}/testData/deplete_test2.fasta ${CMAKE_CURRENT_SOURCE_DIR}/testData/target_test1.fasta
          ${CMAKE_CURRENT_SOURCE_DIR}/testData/target_test2.fasta ${CMAKE_CURRENT_SOURCE_DIR}/testData/reads_test.fastq
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
