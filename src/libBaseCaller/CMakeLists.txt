

if (WIN32)
    file(DOWNLOAD "https://anaconda.org/intel/mkl-static/2020.3/download/win-64/mkl-static-2020.3-intel_279.tar.bz2" "${CMAKE_CURRENT_BINARY_DIR}/mkl-static.tar.bz2")
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "tar" "xvz" "${CMAKE_CURRENT_BINARY_DIR}/mkl-static.tar.bz2" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

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
    add_library(deepnano2 STATIC IMPORTED GLOBAL)
    set_target_properties(deepnano2 PROPERTIES 
        IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/deepnano2.lib"
    )
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