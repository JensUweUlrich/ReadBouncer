add_executable(IBFTest IBFTest.cpp)
#target_link_libraries(IBFTest IBF gmock_main)

# using seqan
if (WIN32)
	add_definitions(-DMSVC_COMPILER)
	target_link_libraries (IBFTest IBF gmock_main config ReadUntil mklcore mklseq mkllp64 deepnano2 DeepNano Guppy)
else()
	target_link_libraries (IBFTest IBF gmock_main config ReadUntil Guppy)
endif()

add_test(NAME IBFTest COMMAND IBFTest)


if (WIN32)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/testDataIBF/ref_1.fasta
         ${CMAKE_CURRENT_SOURCE_DIR}/testDataIBF/ref_2.fasta
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Release)
else()
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/testDataIBF/ref_1.fasta
         ${CMAKE_CURRENT_SOURCE_DIR}/testDataIBF/ref_2.fasta
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
