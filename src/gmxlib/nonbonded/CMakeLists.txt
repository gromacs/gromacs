# Sources that should always be built
file(GLOB NONBONDED_SOURCES *.c nb_kernel_c/*.c)

if(GMX_CPU_ACCELERATION STREQUAL "SSE2" AND NOT GMX_DOUBLE)
    file(GLOB NONBONDED_SSE2_SINGLE_SOURCES nb_kernel_sse2_single/*.c)
endif()

# For now we enable the (existing) SSE4.1 kernels for all higher  accelerations
if(((GMX_CPU_ACCELERATION STREQUAL "SSE4.1") OR
    (GMX_CPU_ACCELERATION STREQUAL "AVX_256")) AND NOT GMX_DOUBLE)
    file(GLOB NONBONDED_SSE4_1_SINGLE_SOURCES nb_kernel_sse4_1_single/*.c)
endif()

if(GMX_CPU_ACCELERATION STREQUAL "AVX_128_FMA" AND NOT GMX_DOUBLE)
    file(GLOB NONBONDED_AVX_128_FMA_SINGLE_SOURCES nb_kernel_avx_128_fma_single/*.c)
endif()

# These sources will be used in the parent directory's CMakeLists.txt
set(NONBONDED_SOURCES ${NONBONDED_SOURCES} ${NONBONDED_SSE2_SINGLE_SOURCES} ${NONBONDED_SSE4_1_SINGLE_SOURCES} ${NONBONDED_AVX_128_FMA_SINGLE_SOURCES} PARENT_SCOPE)



