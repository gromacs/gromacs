# - Define function to detect whether the compiler's target
# - architecture is one for which GROMACS has special treatment
# - (e.g. kernel acceleration)
#
# Sets GMX_IS_X86 or GMX_IS_BGQ if targetting that architecture

function(gmx_detect_target_architecture)
    try_compile(GMX_IS_X86 ${CMAKE_BINARY_DIR}
        "${CMAKE_SOURCE_DIR}/cmake/TestX86.c")
    try_compile(GMX_IS_BGQ ${CMAKE_BINARY_DIR}
        "${CMAKE_SOURCE_DIR}/cmake/TestBlueGeneQ.c")
endfunction()
