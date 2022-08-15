
add_executable(MinKNOWTests minknowtest.cpp)

#Find spdlog and print message if not/found 
set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR}/../../spdlog/spdlog/lib/cmake/spdlog)
find_package(spdlog REQUIRED)

IF (${SPDLOG_FOUND})
   MESSAGE(STATUS "Found spdlog.")
ELSE (${SPDLOG_FOUND})
   MESSAGE(STATUS "Could not locate spdlog.")
ENDIF (${SPDLOG_FOUND})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../util" 
		    "${CMAKE_CURRENT_SOURCE_DIR}/../interfaces" 
		    "${CMAKE_CURRENT_SOURCE_DIR}/../basecaller/Guppy" 
                    "${CMAKE_CURRENT_SOURCE_DIR}/../basecaller")

target_include_directories(${PROJECT_NAME} PUBLIC 
				${CMAKE_CURRENT_BINARY_DIR}/../lyra/lyra/lyra-src/include 
			   	${CMAKE_CURRENT_SOURCE_DIR}/../util 
				${CMAKE_CURRENT_SOURCE_DIR}/../basecaller/Guppy)
                
target_link_libraries(MinKNOWTests IBF config ReadUntil Guppy spdlog::spdlog_header_only gmock_main)
install(TARGETS MinKNOWTests DESTINATION bin COMPONENT applications)