# If the user did not set GMX_GPU we'll consider this option to be
# in "auto" mode meaning that we will:
# - search for CUDA and set GMX_GPU=ON we find it
# - check whether GPUs are present
# - if CUDA is not found but GPUs were detected issue a warning
if (NOT DEFINED GMX_GPU)
    set(GMX_GPU_AUTO TRUE CACHE INTERNAL "")
endif()
option(GMX_GPU "Enable GPU acceleration" OFF)

if(GMX_GPU AND GMX_DOUBLE)
    message(FATAL_ERROR "GPU acceleration is not available in double precision, disabled!")
endif()
if(GMX_GPU_AUTO AND GMX_DOUBLE)
    message(WARNING "GPU acceleration is not available in double precision, disabled!")
    set_property(CACHE GMX_GPU PROPERTY VALUE OFF)
    set_property(CACHE GMX_GPU_AUTO PROPERTY VALUE OFF)
endif()

# detect GPUs in the build host machine
if (GMX_GPU OR GMX_GPU_AUTO)
    include(gmxDetectGpu)
    gmx_detect_gpu()
endif()

# Depending on the current vale of GMX_GPU and GMX_GPU_AUTO:
# - OFF, FALSE: Will skip this detection.
# - OFF, TRUE : Will keep GMX_GPU=OFF if no CUDA is detected, but will assemble
#               a warning message which will be issued at the end of the
#               configuration if GPU(s) were found in the build system.
# - ON , FALSE: The user requested GPU builds, will require CUDA and will fail
#               if it is not available.
# - ON , TRUE : Can't happen (GMX_GPU=ON can only be user-set at this point)
if(GMX_GPU OR GMX_GPU_AUTO)
    # We support CUDA >=v3.2 on *nix, but <= v4.1 doesn't work with MSVC
    if(MSVC)
        find_package(CUDA 4.1)
    else()
        find_package(CUDA 3.2)
        #                     PATHS "/usr/local/cuda" "/usr/local/cuda-5.0" "/opt/cuda")
    endif()

    if (EXISTS ${CUDA_TOOLKIT_ROOT_DIR})
        set(CUDA_FOUND TRUE CACHE INTERNAL "Whether the CUDA toolkit was found" FORCE)
    else()
        set(CUDA_FOUND FALSE CACHE INTERNAL "Whether the CUDA toolkit was found" FORCE)
    endif()

    # assemble warning/error message
    if (GMX_DETECT_GPU_AVAILABLE)
        set(_msg "
    ${GMX_DETECT_GPU_COUNT} NVIDIA GPU(s) found in the system")

        # append GPU names
        if (NOT GMX_DETECT_GPU_INFO STREQUAL "")
            set(_msg "${_msg}:")
            foreach(gpu ${GMX_DETECT_GPU_INFO})
                set(_msg "${_msg}
                ${gpu}")
            endforeach()
        endif()

        # TODO remove the second part of the message when we'll have compute
        # capability information from the detection.
        set(_msg "${_msg}
    Compute capability information not available, consult the NVIDIA website:
    https://developer.nvidia.com/cuda-gpus
            ")
    endif()

        set(CUDA_NOTFOUND_MESSAGE "
    mdrun supports native GPU acceleration on NVIDIA hardware with compute
    capability >=2.0. This requires the NVIDIA CUDA library, which was not
    found; its location can be hinted by setting CUDA_TOOLKIT_ROOT_DIR.
    Note that CPU or GPU acceleration can be selected at runtime!
    ${_msg}")
        unset(_msg)

    if (NOT CUDA_FOUND)
        if (GMX_GPU_AUTO)
            # Disable GPU acceleration in auto mode
            message(STATUS "Disabling native GPU acceleration")
            set_property(CACHE GMX_GPU PROPERTY VALUE OFF)
            set(CUDA_NOTFOUND_AUTO ON)
        else ()
            # the user requested CUDA
            message(FATAL_ERROR "${CUDA_NOTFOUND_MESSAGE}")
        endif()
    else()
        if (GMX_GPU_AUTO)
            message(STATUS "Enabling native GPU acceleration")
            set_property(CACHE GMX_GPU PROPERTY VALUE ON)
        endif()
    endif() # NOT CUDA_FOUND
endif()
# Annoyingly enough, FindCUDA leaves a few variables behind as non-advanced.
# We need to mark these advanced outside the conditional, otherwise, if the
# user turns GMX_GPU=OFF after a failed cmake pass, these variables will be
# left behind in the cache.
mark_as_advanced(CUDA_BUILD_CUBIN CUDA_BUILD_EMULATION CUDA_SDK_ROOT_DIR CUDA_VERBOSE_BUILD)

macro(gmx_gpu_setup)
    # set up nvcc options
    include(gmxManageNvccConfig)

    # Check whether we can use atomic operations needed for polling wait for GPU
    # (to avoid the cudaStreamSynchronize + ECC bug).
    # With thread-MPI testing atomics has already been carried out, but without
    # thread-MPI we need to invoke the atomics test independently.
    if (NOT GMX_THREAD_MPI)
        set(TEST_TMPI_ATOMICS_ONLY ON CACHE INTERNAL
            "Test only the atomic operations of thread-MPI.")
        include(ThreadMPI)
    endif()

    # we need this linker flag in case if we have ld >= 2.22 (typically with gcc 4.5+),
    # but it's too cumbersome to check the ld version and the flag should not hurt
    if(CMAKE_COMPILER_IS_GNUCC)
        set(GROMACS_LINKER_FLAGS "-Wl,--add-needed ${GROMACS_LINKER_FLAGS}")
    endif()

    # no OpenMP is no good!
    if(NOT GMX_OPENMP)
        message(WARNING "
    To use GPU acceleration efficiently, mdrun requires OpenMP multi-threading.
    With no OpenMP a single CPU core can be used with a GPU which is not optimal.
    Note that with MPI multiple processes can be forced to use a single GPU, but this
    typically inefficient.")
    endif()
endmacro()
