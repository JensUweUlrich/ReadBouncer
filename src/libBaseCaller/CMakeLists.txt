
add_library(mklcore STATIC IMPORTED GLOBAL)
set_target_properties(mklcore PROPERTIES 
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/mkl_core.lib"
)
add_library(mklseq STATIC IMPORTED GLOBAL)
set_target_properties(mklseq PROPERTIES 
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/mkl_sequential.lib"
)
add_library(mkllp64 STATIC IMPORTED GLOBAL)
set_target_properties(mkllp64 PROPERTIES 
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/mkl_intel_lp64.lib"
)

add_library(deepnano2 STATIC IMPORTED GLOBAL)
set_target_properties(deepnano2 PROPERTIES 
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/deepnano2.lib"
    )

