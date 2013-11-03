######################################
# Output compiler and CFLAGS used
######################################
include(GetCompilerInfo.cmake)
get_compiler_info(C BUILD_C_COMPILER BUILD_CFLAGS)
get_compiler_info(CXX BUILD_CXX_COMPILER BUILD_CXXFLAGS)
if(GMX_GPU)
    get_cuda_compiler_info(CUDA_NVCC_COMPILER_INFO CUDA_NVCC_COMPILER_FLAGS)
endif(GMX_GPU)

configure_file(config.h.cmakein config.h)
configure_file(buildinfo.h.cmakein buildinfo.h)

if (BUILD_TESTING)
    add_custom_target(tests)
    if (GMX_BUILD_UNITTESTS)
        add_subdirectory(external/gmock-1.6.0)
    endif (GMX_BUILD_UNITTESTS)
    include(testutils/TestMacros.cmake)
    if (GMX_BUILD_UNITTESTS)
        add_subdirectory(testutils)
    endif (GMX_BUILD_UNITTESTS)
endif (BUILD_TESTING)

add_subdirectory(gromacs)
add_subdirectory(programs)

if (NOT GMX_FAHCORE)
    add_subdirectory(contrib)
endif()
