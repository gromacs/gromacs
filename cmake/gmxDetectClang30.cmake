macro(gmx_detect_clang_3_0 OUT_VAR)
    try_compile(IS_CLANG_3_0 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/TestClangVersion.c)
    set(${OUT_VAR} ${IS_CLANG_3_0})
endmacro()

