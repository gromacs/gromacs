# Sources that should always be built
file(GLOB NONBONDED_SOURCES *.c nb_kernel_c/*.c)

# For now we enable the (existing) SSE2 kernels for all accelerations
if(((GMX_ACCELERATION STREQUAL "SSE2") OR
    (GMX_ACCELERATION STREQUAL "SSE4.1") OR
    (GMX_ACCELERATION STREQUAL "AVX_128_FMA") OR
    (GMX_ACCELERATION STREQUAL "AVX_256")) AND NOT GMX_DOUBLE)
    file(GLOB NONBONDED_SSE2_SINGLE_SOURCES nb_kernel_sse2_single/*.c)
endif()

# These sources will be used in the parent directory's CMakeLists.txt
set(NONBONDED_SOURCES ${NONBONDED_SOURCES} ${NONBONDED_SSE2_SINGLE_SOURCES} PARENT_SCOPE)



