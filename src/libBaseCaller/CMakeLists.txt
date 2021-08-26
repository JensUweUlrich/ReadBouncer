if (WIN32)

    configure_file(CMakeLists.txt.in bc-build/CMakeLists.txt)

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bc-build)
    if(result)
        message(FATAL_ERROR "Generator step for basecaller failed: ${result}")
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bc-build)
    if(result)
        message(FATAL_ERROR "Build step for basecaller failed: ${result}")
    endif()

     #download and extract mkl static library
     if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/Library/lib/mkl_core.lib")
        file(DOWNLOAD "https://anaconda.org/intel/mkl-static/2021.3.0/download/win-64/mkl-static-2021.3.0-intel_524.tar.bz2" "${CMAKE_CURRENT_BINARY_DIR}/mkl-static.tar.bz2")
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "tar" "xvz" "${CMAKE_CURRENT_BINARY_DIR}/mkl-static.tar.bz2" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/Library/include/mkl.h")
        file(DOWNLOAD "https://anaconda.org/intel/mkl-include/2021.3.0/download/win-64/mkl-include-2021.3.0-intel_524.tar.bz2" "${CMAKE_CURRENT_BINARY_DIR}/mkl-include.tar.bz2")
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "tar" "xvz" "${CMAKE_CURRENT_BINARY_DIR}/mkl-include.tar.bz2" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    #if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/boost_1_76_0/boost/version.hpp")
    #    file(DOWNLOAD "https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.bz2" "${CMAKE_CURRENT_BINARY_DIR}/boost_1_76_0.tar.bz2")
    #    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "tar" "xvz" "${CMAKE_CURRENT_BINARY_DIR}/boost_1_76_0.tar.bz2" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    #endif()
  
    file(GLOB OSPREYSRCS "${CMAKE_CURRENT_BINARY_DIR}/osprey/src/src/*.cpp")
    file(GLOB OSPREYHDRS "${CMAKE_CURRENT_BINARY_DIR}/osprey/src/src/*.hpp")

    #include_directories("${CMAKE_CURRENT_BINARY_DIR}/NumCpp/include" "${CMAKE_CURRENT_BINARY_DIR}/Library/include" "${CMAKE_CURRENT_BINARY_DIR}/boost_1_76_0")
    include_directories("${CMAKE_CURRENT_BINARY_DIR}/Library/include")

    add_library(decoder STATIC IMPORTED GLOBAL)
        set_target_properties(decoder PROPERTIES 
            IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/osprey/src/src/decoder.lib"
    )

    add_library(mklcore STATIC IMPORTED GLOBAL)
        set_target_properties(mklcore PROPERTIES 
            IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/Library/lib/mkl_core.lib"
    )

    add_library(mklseq STATIC IMPORTED GLOBAL)
        set_target_properties(mklseq PROPERTIES 
            IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/Library/lib/mkl_sequential.lib"
    )
    add_library(mkllp64 STATIC IMPORTED GLOBAL)
        set_target_properties(mkllp64 PROPERTIES 
            IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/Library/lib/mkl_intel_lp64.lib"
    )

    add_library (osprey STATIC ${OSPREYSRCS} ${OSPREYHDRS})

    include_directories("${CMAKE_CURRENT_BINARY_DIR}/osprey/src/src" "${CMAKE_CURRENT_BINARY_DIR}/bc-build/basecaller-prefix/src/basecaller-build/Library/include" 
    #                    "${CMAKE_CURRENT_BINARY_DIR}/bc-build/basecaller-prefix/src/basecaller-build/NumCpp/include" "${CMAKE_CURRENT_BINARY_DIR}/bc-build/basecaller-prefix/src/basecaller-build/boost_1_76_0" 
                        "${CMAKE_CURRENT_SOURCE_DIR}/../util")
    #target_include_directories (osprey PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/osprey/src/src" "${CMAKE_CURRENT_BINARY_DIR}/NumCpp/include" "${CMAKE_CURRENT_BINARY_DIR}/Library/include" "${CMAKE_CURRENT_BINARY_DIR}/boost_1_76_0")
    target_include_directories (osprey PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/osprey/src/src" "${CMAKE_CURRENT_BINARY_DIR}/Library/include")
    target_link_libraries (osprey mklcore mklseq mkllp64 decoder)

    file(GLOB WEIGHTS "${CMAKE_CURRENT_BINARY_DIR}/osprey/src/osprey/weights/*")
    install(FILES ${WEIGHTS} DESTINATION data COMPONENT addfiles)

    install(TARGETS osprey DESTINATION lib COMPONENT libraries)

    #install(FILES  ${READUNTILHDRS} DESTINATION include COMPONENT headers)

    

else()

    file(DOWNLOAD "https://anaconda.org/intel/mkl-static/2020.3/download/linux-64/mkl-static-2020.3-intel_279.tar.bz2" "${CMAKE_CURRENT_BINARY_DIR}/mkl-static.tar.bz2")
    #file(ARCHIVE_EXTRACT INPUT mkl-static.tar.bz2)
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "tar" "xvz" "${CMAKE_CURRENT_BINARY_DIR}/mkl-static.tar.bz2" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    add_library(mklcore STATIC IMPORTED GLOBAL)
    set_target_properties(mklcore PROPERTIES 
        IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libmkl_core.a"
    )
    add_library(mklseq STATIC IMPORTED GLOBAL)
    set_target_properties(mklseq PROPERTIES 
        IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libmkl_sequential.a"
    )
    add_library(mkllp64 STATIC IMPORTED GLOBAL)
    set_target_properties(mkllp64 PROPERTIES 
        IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libmkl_intel_ilp64.a"
    )
    add_library(deepnano2 STATIC IMPORTED GLOBAL)
    set_target_properties(deepnano2 PROPERTIES 
        IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/libdeepnano2.a"
    )
endif()