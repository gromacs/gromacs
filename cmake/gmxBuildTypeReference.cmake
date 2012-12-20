# Custom build type "Reference", to be used for creating new
# reference values in the Gromacs regression tests.
set( CMAKE_CXX_FLAGS_REFERENCE "-O0 -g" CACHE STRING "C++ flags for regressiontests reference runs." FORCE)
set( CMAKE_C_FLAGS_REFERENCE "-O0 -g" CACHE STRING "C flags for regressiontests reference runs." FORCE)
mark_as_advanced( CMAKE_CXX_FLAGS_REFERENCE CMAKE_C_FLAGS_REFERENCE)

# turn off all fancy options for the regressiontests reference build
if("${CMAKE_BUILD_TYPE}" STREQUAL "Reference")
    set(GMX_GPU OFF CACHE BOOL "Disabled for regressiontests reference builds" FORCE)
    set(GMX_OPENMP OFF CACHE BOOL "Disabled for regressiontests reference builds" FORCE)
    set(GMX_CPU_ACCELERATION "None" CACHE STRING "Disabled for regressiontests reference builds" FORCE)
    set(GMX_FFT_LIBRARY "fftpack" CACHE STRING "Use fftpack for regressiontests reference builds" FORCE)
    set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Disabled for regressiontests reference builds" FORCE)
    set(GMX_THREAD_MPI OFF CACHE BOOL "Disabled for regressiontests reference builds" FORCE)

    # C_COMPILER_VERSION is not defined automatically for CMake below 2.8.8,
    # so we call the GROMACS work-around for that
    include(gmxGetCompilerInfo)
    get_compiler_version()
    if(NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" OR NOT "${C_COMPILER_VERSION}" MATCHES "4.7")
        message(WARNING "Reference values for regressiontests should use Gromacs compiled with "
            "gcc 4.7, but your configuration is using ${CMAKE_C_COMPILER_ID}-${C_COMPILER_VERSION}.")
    endif()
endif()
