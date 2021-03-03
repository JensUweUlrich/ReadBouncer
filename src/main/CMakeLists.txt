set(Main_SRC "main.cpp")
file(GLOB MAIN_HDRS "*.hpp")
add_executable (${PROJECT_NAME} ${MAIN_HDRS} ${Main_SRC})
#find_package (Threads)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../lyra/lyra/lyra-src/include)
# using seqan
add_definitions(-DMSVC_COMPILER)
target_link_libraries (${PROJECT_NAME} IBF ReadUntil deepnano2)
#target_link_libraries (${PROJECT_NAME} ReadUntil)
#target_link_libraries (${PROJECT_NAME} pthread ReadUntil deepnano2)
install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT applications)
