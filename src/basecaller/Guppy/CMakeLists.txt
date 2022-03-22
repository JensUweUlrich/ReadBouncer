
#configure_file(CMakeLists.txt.in Guppy-build/CMakeLists.txt)

#execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#  RESULT_VARIABLE result
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Guppy-build)
#if(result)
#  message(FATAL_ERROR "Generator step for Guppy failed: ${result}")
#endif()

#execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release
#  RESULT_VARIABLE result
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Guppy-build)
#if(result)
#  message(FATAL_ERROR "Build step for Guppy failed: ${result}")
#endif()

file(GLOB GUPPYSRCS "*.cpp")
file(GLOB GUPPYHDRS "*.hpp" "*.h")

add_library(GuppyClient STATIC IMPORTED GLOBAL)
if (WIN32)
  execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "tar" "xvzf" "${CMAKE_CURRENT_SOURCE_DIR}/GuppyCPPClient.tar.gz" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})  
  set_target_properties(GuppyClient PROPERTIES 
      IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/GuppyCPPClient.lib"
   )
elseif(ARM_BUILD)
  execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "tar" "xvzf" "${CMAKE_CURRENT_SOURCE_DIR}/libGuppyCppClient_linux_arm_64.tar.gz" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_target_properties(GuppyClient PROPERTIES 
      IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libGuppyCppClient.a"
   )
else()
  
  execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "tar" "xvzf" "${CMAKE_CURRENT_SOURCE_DIR}/libGuppyCppClient_linux_x86_64.tar.gz" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
      set_target_properties(GuppyClient PROPERTIES 
        IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libGuppyCppClient_gcc-8.a"
      )
  elseif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
      set_target_properties(GuppyClient PROPERTIES 
        IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libGuppyCppClient_gcc-9.a"
      )
  else()
      set_target_properties(GuppyClient PROPERTIES 
        IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libGuppyCppClient_gcc-10.a"
      )
  endif()
endif()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../util" "${CMAKE_CURRENT_SOURCE_DIR}/../../interfaces" ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/..")

add_library (Guppy STATIC ${GUPPYSRCS} ${GUPPYHDRS})
target_include_directories (Guppy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


if (WIN32)
   target_link_libraries(Guppy
        PUBLIC
        GuppyClient
        ws2_32
        Iphlpapi
        )
else()
   target_link_libraries(Guppy GuppyClient)
endif() 

install(TARGETS Guppy DESTINATION lib COMPONENT libraries)

install(FILES ${GUPPYHDRS} DESTINATION include COMPONENT headers)



