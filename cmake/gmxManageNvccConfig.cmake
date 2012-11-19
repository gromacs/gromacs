# Manage CUDA nvcc compilation configuration, try to be smart to ease the users'
# pain as much as possible:
# - use the CUDA_NVCC_HOST_COMPILER if defined by the user, otherwise
# - auto-detect compatible nvcc host compiler and set nvcc -ccbin (if not MPI wrapper)
# - set icc compatibility mode to gcc 4.4/4.5 (CUDA 4.0 is not compatible with gcc >v4.4)
# - (advanced) variables set:
#   * CUDA_NVCC_HOST_COMPILER       - the compiler nvcc is forced to use (via -ccbin)
#   * CUDA_NVCC_HOST_COMPILER_OPTIONS   - the full host-compiler related option list passed to nvcc
if (NOT DEFINED CUDA_NVCC_FLAGS_SET)
    set(CUDA_NVCC_FLAGS_SET TRUE CACHE INTERNAL "True if NVCC flags have been set" FORCE)

    # Explicitly set the host compiler for nvcc if the current compiler is
    # supported and it's not an MPI compiler wrapper, otherwise warn the user.
    #
    # Note that even though nvcc compiles host code as C++, we use the
    # CMAKE_C_COMPILER as host compiler. We do this because CUDA versions
    # preceding 5.0 only recognize icc, but not icpc. However, both gcc and icc
    # (i.e. all supported compilers) happily compile C++ code.
    #
    # Also note that with MSVC nvcc sets the -compiler-bindir option behind the
    # scenes; to avoid conflicts we don't set -ccbin automatically.

    if (NOT DEFINED CUDA_NVCC_HOST_COMPILER AND NOT MSVC)
        if (NOT CMAKE_COMPILER_IS_GNUCC AND
            NOT (CMAKE_C_COMPILER_ID MATCHES "Intel" AND UNIX AND NOT APPLE))
            message(WARNING "
            Will not set the nvcc host compiler because the current C compiler is not
            compatible with nvcc:
            ${CMAKE_C_COMPILER} (ID: ${CMAKE_C_COMPILER_ID})
            Compatible compilers are: gcc on Linux and Mac OS X, the Intel Compiler on
            64-bit Linux and MSVC on Windows. If nothing specified, nvcc will automatically
            pick the platform-default compiler; However, as mixing compilers can cause errors.
            To manually set the nvcc host compiler, edit CUDA_NVCC_FLAGS or re-configure
            setting CUDA_NVCC_HOST_COMPILER to the full path of a compatible compiler.
            ")
        else()
            # do not use MPI compiler wrappers, as these are prone to brake nvcc
            if (GMX_MPI AND
                NOT "${${MPI_PREFIX}_FOUND}" AND # FindMPI-based detection
                NOT GMX_THREAD_MPI)
                message(WARNING "
            Will not set the nvcc host compiler because the current C compiler is an MPI
            compiler wrapper: ${CMAKE_C_COMPILER}
            MPI compiler wrappers are prone to not work with nvcc. You might get lucky,
            but the safest is to use the C compiler that the MPI compiler wrapper uses
            (if this is compatible).
            To manually set the nvcc host compiler, edit CUDA_NVCC_FLAGS or re-configure
            setting CUDA_NVCC_HOST_COMPILER to the full path of a compatible compiler.
            ")
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
        # On *nix force icc in gcc 4.4 compatibility with CUDA 3.2/4.0 and
        # gcc 4.5 compatibility mode with later CUDA versions. This is needed
        # as even with icc use as host compiler, when icc's gcc compatibility
        # mode is higher than the max gcc version officially supported by CUDA,
        # nvcc will freak out.
        if (UNIX AND CMAKE_C_COMPILER_ID MATCHES "Intel" AND
            CUDA_NVCC_HOST_COMPILER_AUTOSET)
            if (CUDA_VERSION VERSION_LESS "4.1")
                message(STATUS "Setting Intel Compiler compatibity mode to gcc 4.4 for nvcc host compilation")
                set(CUDA_NVCC_HOST_COMPILER_OPTIONS "${CUDA_NVCC_HOST_COMPILER_OPTIONS};-Xcompiler;-gcc-version=440;")
            else()
                message(STATUS "Setting Intel Compiler compatibity mode to gcc 4.5 for nvcc host compilation")
                set(CUDA_NVCC_HOST_COMPILER_OPTIONS "${CUDA_NVCC_HOST_COMPILER_OPTIONS};-Xcompiler;-gcc-version=450;")
            endif()
        endif()
        set(CUDA_NVCC_HOST_COMPILER_OPTIONS "${CUDA_NVCC_HOST_COMPILER_OPTIONS}"
            CACHE STRING "Host-side compiler and options for it (do not edit!)." FORCE)

        mark_as_advanced(CUDA_NVCC_HOST_COMPILER CUDA_NVCC_HOST_COMPILER_OPTIONS)
    endif()

    # on Linux we need to add -fPIC when building shared gmx libs
    # Note: will add -fPIC for any compiler that supports it as it shouldn't hurt
    if(BUILD_SHARED_LIBS)
        GMX_TEST_CXXFLAG(CXXFLAG_FPIC "-fPIC" _FPIC_NVCC_FLAG)
        if(_FPIC_NVCC_FLAG)
            set(_FPIC_NVCC_FLAG "-Xcompiler;${_FPIC_NVCC_FLAG}")
        endif()
    endif()

    # Set the CUDA GPU architectures to compile for:
    # - with CUDA >v4.2 compute capability 2.0, 2.1 is, but 3.0 is not supported:
    #     => compile sm_20, sm_21 cubin, and compute_20 PTX
    # - with CUDA >=4.2 compute capabity <=3.0 is supported:
    #     => compile sm_20, sm_21, sm_30 cubin, and compute_30 PTX
    # - with CUDA 5.0 compute capabity 3.5 is supported, but generating code
    #   optimized for sm_35 results in lower performance than with sm_30.
    if(CUDA_VERSION VERSION_LESS "4.2")
        set(_CUDA_ARCH_STR "-gencode;arch=compute_20,code=sm_20;-gencode;arch=compute_20,code=sm_21;-gencode;arch=compute_20,code=compute_20")
    else()
        set(_CUDA_ARCH_STR "-gencode;arch=compute_20,code=sm_20;-gencode;arch=compute_20,code=sm_21;-gencode;arch=compute_30,code=sm_30;-gencode;arch=compute_30,code=compute_30")
    endif()

    # finally set the damn flags
    set(CUDA_NVCC_FLAGS
        "${_CUDA_ARCH_STR};-use_fast_math;${CUDA_NVCC_HOST_COMPILER_OPTIONS};${_FPIC_NVCC_FLAG}"
        CACHE STRING "Compiler flags for nvcc." FORCE)
endif()
