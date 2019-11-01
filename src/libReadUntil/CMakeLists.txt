
configure_file(CMakeLists.txt.in gRPC-build/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gRPC-build )
if(result)
  message(FATAL_ERROR "Generator step for gRPC failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gRPC-build)
if(result)
  message(FATAL_ERROR "Build step for gRPC failed: ${result}")
endif()



set(READUNTILSOURCES "ReadUntilClient.cpp" "ReadUntilClient.hpp")


# Include the gRPC's cmake build (normally grpc source code would live
  # in a git submodule called "third_party/grpc", but this example lives in
  # the same repository as gRPC sources, so we just look a few directories up)
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/grpc ${CMAKE_CURRENT_BINARY_DIR}/grpc EXCLUDE_FROM_ALL)
  message(STATUS "Using gRPC via add_subdirectory.")
  
  # After using add_subdirectory, we can now use the grpc targets directly from
  # this build.
  set(_PROTOBUF_LIBPROTOBUF libprotobuf)
#  set(_PROTOBUF_PROTOC $<TARGET_FILE:protoc>)
  set(_GRPC_GRPCPP_UNSECURE grpc++_unsecure)
#  set(_GRPC_CPP_PLUGIN_EXECUTABLE $<TARGET_FILE:grpc_cpp_plugin>)


add_library (ReadUntil STATIC ${READUNTILSOURCES})

target_link_libraries (ReadUntil  ${_GRPC_GRPCPP_UNSECURE}
    ${_PROTOBUF_LIBPROTOBUF})

target_include_directories (ReadUntil PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

