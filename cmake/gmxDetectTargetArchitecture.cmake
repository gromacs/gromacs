# - Define function to detect whether the compiler's target
# - architecture is one for which GROMACS has special treatment
# - (e.g. kernel acceleration)
#
# Sets GMX_TARGET_X86 or GMX_TARGET_BGQ if targetting that
# architecture. May set other such variables if/when there is future
# need.

function(gmx_detect_target_architecture)
    if (NOT DEFINED GMX_TARGET_X86)
        try_compile(GMX_TARGET_X86 ${CMAKE_BINARY_DIR}
            "${CMAKE_SOURCE_DIR}/cmake/TestX86.c")
    endif()
    if (NOT DEFINED GMX_TARGET_BGQ)
        try_compile(GMX_TARGET_BGQ ${CMAKE_BINARY_DIR}
            "${CMAKE_SOURCE_DIR}/cmake/TestBlueGeneQ.c")
    endif()
endfunction()
