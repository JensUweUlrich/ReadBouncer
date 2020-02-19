
configure_file(CMakeLists.txt.in gRPC-build/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gRPC-build)
if(result)
  message(FATAL_ERROR "Generator step for gRPC failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -D_WIN32_WINNT=0x600 --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gRPC-build)
if(result)
  message(FATAL_ERROR "Build step for gRPC failed: ${result}")
endif()

file(GLOB READUNTILSOURCES "*.cpp" "*.hpp")


# Include the gRPC's cmake build (normally grpc source code would live
  # in a git submodule called "third_party/grpc", but this example lives in
  # the same repository as gRPC sources, so we just look a few directories up)
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/grpc ${CMAKE_CURRENT_BINARY_DIR}/grpc EXCLUDE_FROM_ALL)
  message(STATUS "Using gRPC via add_subdirectory.")

  #message(STATUS "Using gRPC via add_subdirectory (FetchContent).")
 # include(FetchContent)
 # FetchContent_Declare(
 #   grpc
 #   GIT_REPOSITORY https://github.com/grpc/grpc.git
    # when using gRPC, you will actually set this to an existing tag, such as
    # v1.25.0, v1.26.0 etc..
    # For the purpose of testing, we override the tag used to the commit
    # that's currently under test.
 #   GIT_TAG v1.27.1)

#  FetchContent_MakeAvailable(grpc)
  
  # After using add_subdirectory, we can now use the grpc targets directly from
  # this build.
  set(_PROTOBUF_LIBPROTOBUF libprotobuf)
  set(_PROTOBUF_PROTOC $<TARGET_FILE:protoc>)
  set(_GRPC_GRPCPP_UNSECURE grpc++_unsecure)
  #find_program(_GRPC_CPP_PLUGIN_EXECUTABLE grpc_cpp_plugin)
  set(_GRPC_CPP_PLUGIN_EXECUTABLE $<TARGET_FILE:grpc_cpp_plugin>)

set(Protobuf_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/protobuf/include")
set(Protobuf_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/protobuf/lib")

# Proto file
get_filename_component(reads_proto "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.proto" ABSOLUTE)
get_filename_component(reads_proto_path "${reads_proto}" PATH)

# Generated sources
set(reads_proto_srcs "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.pb.cc")
set(reads_proto_hdrs "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.pb.h")
set(reads_grpc_srcs "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.grpc.pb.cc")
set(reads_grpc_hdrs "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.grpc.pb.h")
add_custom_command(
      OUTPUT "${reads_proto_srcs}" "${reads_proto_hdrs}" "${reads_grpc_srcs}" "${reads_grpc_hdrs}"
      COMMAND ${_PROTOBUF_PROTOC}
      ARGS --grpc_out "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow"
        --cpp_out "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow"
        -I "${reads_proto_path}"
        --plugin=protoc-gen-grpc=${_GRPC_CPP_PLUGIN_EXECUTABLE}
        "${reads_proto}"
      DEPENDS "${reads_proto}")
 
# Include generated *.pb.h files
list(APPEND minknow_srcs ${reads_proto_srcs} ${reads_grpc_srcs})
file(GLOB protofiles "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/rpc/*.proto")

foreach(proto ${protofiles})
  get_filename_component(protofilepath ${proto} PATH)
  get_filename_component(fname ${proto} NAME_WE)
  add_custom_command(
      OUTPUT "${protofilepath}/${fname}.pb.cc" "${protofilepath}/${fname}.pb.h" "${protofilepath}/${fname}.grpc.pb.cc" "${protofilepath}/${fname}.grpc.pb.h"
      COMMAND ${_PROTOBUF_PROTOC}
      ARGS --grpc_out "${CMAKE_CURRENT_BINARY_DIR}/minknow"
        --cpp_out "${CMAKE_CURRENT_BINARY_DIR}/minknow"
        -I "${CMAKE_CURRENT_BINARY_DIR}/minknow"
        -I "${CMAKE_CURRENT_BINARY_DIR}/protobuf/src"
        --plugin=protoc-gen-grpc=${_GRPC_CPP_PLUGIN_EXECUTABLE}
        "${proto}"
      DEPENDS "${proto}")
  message(STATUS "building target for ${fname}.proto")
  list(APPEND minknow_srcs "${protofilepath}/${fname}.pb.cc" "${protofilepath}/${fname}.grpc.pb.cc")
endforeach()

# flappie build parameters
# include hdf5 source directory for ubuntu linux
#include_directories("/usr/include/hdf5/serial/")

#set(flappie_dir "${CMAKE_CURRENT_BINARY_DIR}/flappie")

#configure_file (
#    "${flappie_dir}/src/version.h.in"
#    "${flappie_dir}/include/version.h"
#)




#file(GLOB flappie_srcs "${CMAKE_CURRENT_BINARY_DIR}/flappie/src/*.c")
#file(GLOB flappie_hdrs "${CMAKE_CURRENT_BINARY_DIR}/flappie/src/*.h")

include_directories("${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow" "${CMAKE_CURRENT_BINARY_DIR}/minknow")


#add_library (ReadUntil STATIC ${READUNTILSOURCES} ${minknow_srcs} ${flappie_srcs})
add_library (ReadUntil STATIC ${READUNTILSOURCES} ${minknow_srcs})
#target_include_directories (ReadUntil PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_BINARY_DIR}/minknow" "${CMAKE_CURRENT_BINARY_DIR}/flappie")
target_include_directories (ReadUntil PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_BINARY_DIR}/minknow")
target_link_libraries (ReadUntil  ${_GRPC_GRPCPP_UNSECURE} ${_PROTOBUF_LIBPROTOBUF})



