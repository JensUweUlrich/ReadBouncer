
cmake_minimum_required (VERSION 3.16)

add_compile_options(-DDEBUG)
#add_compile_options(-DBUILD_SHARED_LIBS=ON)
add_compile_options(-Wno-format)
set(CMAKE_GENERATORS "MinGW Makefiles")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
set(CMAKE_ASM_NASM_OBJECT_FORMAT win32)
add_compile_options(-D_WIN32_WINNT=0x600)



project (RUHEL CXX C)
# add all subdirectories
add_subdirectory(seqan)
add_subdirectory(libReadUntil)
add_subdirectory(libMinimizer)
add_subdirectory(libBloomFilter)
add_subdirectory(main)
add_subdirectory(test)

