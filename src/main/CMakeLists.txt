set(Main_SRC "main.cpp")
file(GLOB MAIN_HDRS "*.hpp")
add_executable (${PROJECT_NAME} ${MAIN_HDRS} ${Main_SRC})
#find_package (Threads)

# include spd fast logging
set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../spdlog/spdlog/lib/cmake/spdlog)
find_package(spdlog REQUIRED)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../util" 
		    "${CMAKE_CURRENT_SOURCE_DIR}/../interfaces" 
		    "${CMAKE_CURRENT_SOURCE_DIR}/../basecaller/Guppy" 
                    "${CMAKE_CURRENT_SOURCE_DIR}/../basecaller")

target_include_directories(${PROJECT_NAME} PUBLIC 
				${CMAKE_CURRENT_BINARY_DIR}/../lyra/lyra/lyra-src/include 
			   	${CMAKE_CURRENT_SOURCE_DIR}/../util 
				${CMAKE_CURRENT_SOURCE_DIR}/../basecaller/Guppy)


if(NOT ARM_BUILD)
#	FIND_PATH(DEEPNANO_INCLUDE_DIR DeepNano2.h "${CMAKE_CURRENT_SOURCE_DIR}/../basecaller/DeepNano" DOC "Deepnano2 includes")
	include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../basecaller/DeepNano")
	target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../basecaller/DeepNano)
endif()

# using seqan
if (WIN32)
	add_definitions(-DMSVC_COMPILER)
	target_link_libraries (${PROJECT_NAME} IBF ReadUntil mklcore mklseq mkllp64 deepnano2 DeepNano Guppy spdlog::spdlog_header_only)
else()
	target_link_libraries (${PROJECT_NAME} IBF ReadUntil Guppy spdlog::spdlog_header_only)
	if (NOT ARM_BUILD)
		target_link_libraries(${PROJECT_NAME} DeepNano deepnano2 mkllp64 mklseq mklcore)
	endif()
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT applications)
