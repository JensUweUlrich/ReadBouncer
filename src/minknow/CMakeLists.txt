
# the location where protobuf-config.cmake will be installed varies by platform
if (WIN32)
  set(_FINDPACKAGE_PROTOBUF_CONFIG_DIR "${CMAKE_CURRENT_BINARY_DIR}/protobuf/cmake")
else()
  set(_FINDPACKAGE_PROTOBUF_CONFIG_DIR "${CMAKE_CURRENT_BINARY_DIR}/gRPC-build/protobuf/src/protobuf-build/lib/cmake/protobuf")
endif()
set(Protobuf_DIR ${_FINDPACKAGE_PROTOBUF_CONFIG_DIR})

if (WIN32)
  set(grpc_cache_args  "-DgRPC_INSTALL:BOOL=ON
        -DgRPC_BUILD_TESTS:BOOL=OFF
        -D_WIN32_WINNT:STRING=0x600
        -DgRPC_PROTOBUF_PROVIDER:STRING=package
        -DgRPC_PROTOBUF_PACKAGE_TYPE:STRING=CONFIG
        -DProtobuf_DIR:PATH=${Protobuf_DIR}
        -DgRPC_CARES_PROVIDER:STRING=package
        -Dc-ares_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/c-ares/lib/cmake/c-ares
        -DgRPC_ZLIB_PROVIDER:STRING=package
        -DZLIB_ROOT:STRING=${CMAKE_CURRENT_BINARY_DIR}/zlib
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/grpc")
else()
  set(grpc_cache_args  "-DgRPC_INSTALL:BOOL=ON
        -DgRPC_BUILD_TESTS:BOOL=OFF
        -DgRPC_PROTOBUF_PROVIDER:STRING=package
        -DgRPC_PROTOBUF_PACKAGE_TYPE:STRING=CONFIG
        -DProtobuf_DIR:PATH=${Protobuf_DIR}
        -DgRPC_CARES_PROVIDER:STRING=package
        -Dc-ares_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/c-ares/lib/cmake/c-ares
        -DgRPC_ZLIB_PROVIDER:STRING=package
        -DZLIB_ROOT:STRING=${CMAKE_CURRENT_BINARY_DIR}/zlib
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/grpc")
endif()


configure_file(CMakeLists.txt.in gRPC-build/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gRPC-build)
if(result)
  message(FATAL_ERROR "Generator step for gRPC failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gRPC-build)
if(result)
  message(FATAL_ERROR "Build step for gRPC failed: ${result}")
endif()

file(GLOB READUNTILSRCS "*.cpp")
file(GLOB READUNTILHDRS "*.hpp")

# Include the gRPC's cmake build (normally grpc source code would live
# in a git submodule called "third_party/grpc", but this example lives in
# the same repository as gRPC sources, so we just look a few directories up)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/grpc ${CMAKE_CURRENT_BINARY_DIR}/grpc EXCLUDE_FROM_ALL)
message(STATUS "Using gRPC via add_subdirectory.")
  
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
#get_filename_component(reads_proto "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.proto" ABSOLUTE)
#get_filename_component(reads_proto_path "${reads_proto}" PATH)

# Generated sources
#set(reads_proto_srcs "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.pb.cc")
#set(reads_proto_hdrs "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.pb.h")
#set(reads_grpc_srcs "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.grpc.pb.cc")
#set(reads_grpc_hdrs "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow/reads.grpc.pb.h")
#add_custom_command(
#      OUTPUT "${reads_proto_srcs}" "${reads_proto_hdrs}" "${reads_grpc_srcs}" "${reads_grpc_hdrs}"
#      COMMAND ${_PROTOBUF_PROTOC}
#      ARGS --grpc_out "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow"
#        --cpp_out "${CMAKE_CURRENT_BINARY_DIR}/minknow/minknow"
#        -I "${reads_proto_path}"
#        --plugin=protoc-gen-grpc=${_GRPC_CPP_PLUGIN_EXECUTABLE}
#        "${reads_proto}"
#      DEPENDS "${reads_proto}")
 
# Include generated *.pb.h files
#list(APPEND minknow_srcs ${reads_proto_srcs} ${reads_grpc_srcs})
file(GLOB protofiles "${CMAKE_CURRENT_BINARY_DIR}/minknow/proto/minknow_api/*.proto")

foreach(proto ${protofiles})
  get_filename_component(protofilepath ${proto} PATH)
  get_filename_component(fname ${proto} NAME_WE)
  add_custom_command(
      OUTPUT "${protofilepath}/${fname}.pb.cc" "${protofilepath}/${fname}.pb.h" "${protofilepath}/${fname}.grpc.pb.cc" "${protofilepath}/${fname}.grpc.pb.h"
      COMMAND ${_PROTOBUF_PROTOC}
      ARGS --grpc_out "${CMAKE_CURRENT_BINARY_DIR}/minknow/proto"
        --cpp_out "${CMAKE_CURRENT_BINARY_DIR}/minknow/proto"
        -I "${CMAKE_CURRENT_BINARY_DIR}/minknow/proto"
        -I "${CMAKE_CURRENT_BINARY_DIR}/protobuf/src"
        --plugin=protoc-gen-grpc=${_GRPC_CPP_PLUGIN_EXECUTABLE}
        --experimental_allow_proto3_optional
        "${proto}"
      DEPENDS "${proto}")
  message(STATUS "building target for ${fname}.proto")
  list(APPEND minknow_srcs "${protofilepath}/${fname}.pb.cc" "${protofilepath}/${fname}.grpc.pb.cc")
endforeach()


set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../spdlog/spdlog/lib/cmake/spdlog)
find_package(spdlog REQUIRED)

# add certificates
file(GLOB WEIGHTS "${CMAKE_CURRENT_SOURCE_DIR}/../rpc-certs/*")
install(FILES ${WEIGHTS} DESTINATION rpc-certs COMPONENT certificates)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/minknow/proto" "${CMAKE_CURRENT_SOURCE_DIR}/../util" "${CMAKE_CURRENT_BINARY_DIR}/stduuid/include" 
                    "${CMAKE_CURRENT_SOURCE_DIR}/../interfaces")

add_library (ReadUntil STATIC ${READUNTILSRCS} ${READUNTILHDRS} ${minknow_srcs})
target_include_directories (ReadUntil PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_BINARY_DIR}/minknow/proto" 
                                "${CMAKE_CURRENT_BINARY_DIR}/stduuid/include" "${CMAKE_CURRENT_SOURCE_DIR}/../interfaces")


if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
	target_link_libraries(ReadUntil grpc++ ${_PROTOBUF_LIBPROTOBUF} spdlog::spdlog_header_only interfaces stdc++fs)
else()
	target_link_libraries (ReadUntil grpc++ ${_PROTOBUF_LIBPROTOBUF} spdlog::spdlog_header_only interfaces)
endif()
install(TARGETS ReadUntil DESTINATION lib COMPONENT libraries)

install(FILES  ${READUNTILHDRS} DESTINATION include COMPONENT headers)



