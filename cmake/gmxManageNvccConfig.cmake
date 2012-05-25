# Manage CUDA nvcc compilation configuration, try to be smart to ease the users' 
# pain as much as possible: 
# - use the CUDA_NVCC_HOST_COMPILER if defined by the user, otherwise
# - auto-detect compatible nvcc host compiler and set nvcc -ccbin (if not MPI wrapper)
# - set icc compatiblity mode to gcc 4.4 (CUDA 4.0 is not compatible with gcc >v4.4.x)
# - (advanced) variables set:
#   * CUDA_NVCC_HOST_COMPILER       - the compier nvcc is forced to use (via -ccbin)
#   * CUDA_NVCC_HOST_COMPILER_OPTIONS   - the full host-compiler related option list passed to nvcc
if (NOT DEFINED CUDA_NVCC_FLAGS_SET)
    set(CUDA_NVCC_FLAGS_SET TRUE CACHE INTERNAL "True if NVCC flags have been set" FORCE)

    # Set the host compiler for nvcc explicitly if the current cimpiler is
    # supported, otherwise warn if the host compiler is not supported.
    # Note that with MSVC nvcc sets the -compiler-bindir option behind the
    # scenes; to avoid conflicts we shouldn't set -ccbin automatically.
    if (NOT DEFINED CUDA_NVCC_HOST_COMPILER AND NOT MSVC)
        if (NOT CMAKE_COMPILER_IS_GNUCC AND
            NOT (CMAKE_C_COMPILER_ID MATCHES "Intel" AND UNIX AND NOT APPLE))
            message(WARNING "
            Will not set the nvcc host compiler because the current C compiler (ID: ${CMAKE_C_COMPILER_ID}): 
            ${CMAKE_C_COMPILER}
            is not compatible with nvcc. Compatible compilers are: gcc on Linux and Mac OS X,
            Intel Compiler on 64-bit Linux and MSVC on Windows. nvcc will pick the platform
            default; however, note that mixing compilers might lead to errors. 
            To set the nvcc host compiler, edit CUDA_NVCC_FLAGS or re-configure with
            CUDA_NVCC_HOST_COMPILER variable.")
        else()
            # the MPI wrappers might not work for compilation
            if (GMX_MPI AND NOT GMX_THREAD_MPI)
                message(WARNING "
            Will not set the nvcc host compiler because the current C compiler is an MPI 
            compiler wrapper: ${CMAKE_C_COMPILER}
            which is prone to not work with nvcc, but you might get lucky.
            To set the nvcc host compiler, edit CUDA_NVCC_FLAGS or re-configure with 
            CUDA_NVCC_HOST_COMPILER variable.")
            else()
                set(CUDA_NVCC_HOST_COMPILER "${CMAKE_C_COMPILER}")
                set(CUDA_NVCC_HOST_COMPILER_AUTOSET TRUE CACHE INTERNAL
                    "True if CUDA_NVCC_HOST_COMPILER is automatically set" FORCE)
            endif()
        endif()
    endif()

    if(DEFINED CUDA_NVCC_HOST_COMPILER)
        message(STATUS "Setting the nvcc host compiler to: ${CUDA_NVCC_HOST_COMPILER}")
        set(CUDA_NVCC_HOST_COMPILER ${CUDA_NVCC_HOST_COMPILER}
            CACHE PATH "Host compiler for nvcc (do not edit!)" FORCE)
        
        set(CUDA_NVCC_HOST_COMPILER_OPTIONS "-ccbin=${CUDA_NVCC_HOST_COMPILER}")
        # force icc in gcc 4.4 compatiblity mode on *NIX to make nvcc 3.2/4.0 happy
        if (UNIX AND CMAKE_C_COMPILER_ID MATCHES "Intel" AND
            CUDA_NVCC_HOST_COMPILER_AUTOSET)
            message(STATUS "Setting Intel Compiler compatibity mode to gcc 4.4 for nvcc host compilation")
            set(CUDA_NVCC_HOST_COMPILER_OPTIONS "${CUDA_NVCC_HOST_COMPILER_OPTIONS};-Xcompiler;-gcc-version=440;")
        endif()
        set(CUDA_NVCC_HOST_COMPILER_OPTIONS "${CUDA_NVCC_HOST_COMPILER_OPTIONS}"
            CACHE STRING "Host-side compiler and options for nvcc (do not edit!)." FORCE)

        mark_as_advanced(CUDA_NVCC_HOST_COMPILER CUDA_NVCC_HOST_COMPILER_OPTIONS)
    endif()

    # on Linux we need to add -fPIC when building shared gmx libs
    # Note: will add -fPIC for any compiler that supports it as it shouldn't hurt
    if(BUILD_SHARED_LIBS)
        GMX_TEST_CXXFLAG(CXXFLAG_FPIC "-fPIC" _FPIC_NVCC_FLAG)
        if(_FPIC_NVCC_FLAG)
            set(_FPIC_NVCC_FLAG "-Xcompiler;${_FPIC_NVCC_FLAG};")
        endif()
    endif()

    # finally set the damn flags
    set(CUDA_NVCC_FLAGS
        "${_FPIC_NVCC_FLAG}-arch=sm_20;-use_fast_math;${CUDA_NVCC_HOST_COMPILER_OPTIONS}"
        CACHE STRING "Compiler flags for nvcc." FORCE)
endif()
