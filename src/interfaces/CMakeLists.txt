file(GLOB INT_SRC "*.cpp")
file(GLOB INT_HDRS "*.hpp")

set(LIBNAME interfaces)

add_library(${LIBNAME} ${INT_SRC} ${INT_HDRS})

if (WIN32)
	add_definitions(-DMSVC_COMPILER)
endif()

# include spd fast logging
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../util" ${CMAKE_CURRENT_SOURCE_DIR})

# using seqan
target_include_directories (${LIBNAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#install(TARGETS IBF DESTINATION lib COMPONENT libraries)

#install(FILES ${IBF_HDRS} DESTINATION include COMPONENT headers)
#install(FILES ${WEIGHTS} DESTINATION data COMPONENT addfiles)
