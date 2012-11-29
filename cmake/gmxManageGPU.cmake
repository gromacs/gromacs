# If the user did not set GMX_GPU we'll consider this option to be
# in "auto" mode meaning that we will:
# - search for CUDA and set GMX_GPU=ON we find it
# - check whether GPUs are present
# - if CUDA is not found but GPUs were detected issue a warning
if (NOT DEFINED GMX_GPU)
    set(GMX_GPU_AUTO TRUE CACHE INTERNAL "GPU acceleration will be selected automatically")
endif()
option(GMX_GPU "Enable GPU acceleration" OFF)

if(GMX_GPU AND GMX_DOUBLE)
    message(FATAL_ERROR "GPU acceleration is not available in double precision!")
endif()
if(GMX_GPU_AUTO AND GMX_DOUBLE)
    message(WARNING "GPU acceleration is not available in double precision, disabled!")
    set_property(CACHE GMX_GPU PROPERTY VALUE OFF)
    set_property(CACHE GMX_GPU_AUTO PROPERTY VALUE OFF)
endif()

# detect GPUs in the build host machine
if (GMX_GPU OR GMX_GPU_AUTO AND NOT GMX_GPU_DETECTION_DONE)
    include(gmxDetectGpu)
    gmx_detect_gpu()
endif()

# We need to call find_package even when we've already done the detection/setup
if(GMX_GPU OR GMX_GPU_AUTO)
    if(GMX_GPU_DETECTION_DONE)
        set(FIND_CUDA_QUIETLY QUIET)
    endif()

    # We support CUDA >=v3.2 on *nix, but <= v4.1 doesn't work with MSVC
    if(MSVC)
        set(FIND_CUDA_VERSION_REQUIRED 4.1)
    else()
        set(FIND_CUDA_VERSION_REQUIRED 3.2)
    endif()
    find_package(CUDA ${FIND_CUDA_VERSION_REQUIRED} ${FIND_CUDA_QUIETLY})
endif()

# Depending on the current vale of GMX_GPU and GMX_GPU_AUTO:
# - OFF, FALSE: Will skip this detection/setup.
# - OFF, TRUE : Will keep GMX_GPU=OFF if no CUDA is detected, but will assemble
#               a warning message which will be issued at the end of the
#               configuration if GPU(s) were found in the build system.
# - ON , FALSE: The user requested GPU builds, will require CUDA and will fail
#               if it is not available.
# - ON , TRUE : Can't happen (GMX_GPU=ON can only be user-set at this point)
if((GMX_GPU OR GMX_GPU_AUTO) AND NOT GMX_GPU_DETECTION_DONE)
    if (EXISTS ${CUDA_TOOLKIT_ROOT_DIR})
        set(CUDA_FOUND TRUE CACHE INTERNAL "Whether the CUDA toolkit was found" FORCE)
    else()
        set(CUDA_FOUND FALSE CACHE INTERNAL "Whether the CUDA toolkit was found" FORCE)
    endif()

     set(CUDA_NOTFOUND_MESSAGE "
    mdrun supports native GPU acceleration on NVIDIA hardware with compute
    capability >=2.0. This requires the NVIDIA CUDA library, which was not
    found; the location can be hinted by setting CUDA_TOOLKIT_ROOT_DIR as
    a CMake option (It does not work as an environment variable).
    The typical location would be /usr/local/cuda[-version].
    Note that CPU or GPU acceleration can be selected at runtime!")

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

        set(CUDA_NOTFOUND_MESSAGE "${_msg}${CUDA_NOTFOUND_MESSAGE}")
        unset(_msg)
    endif()

    if (NOT CUDA_FOUND)
        if (GMX_GPU_AUTO)
            # Disable GPU acceleration in auto mode
            message(STATUS "Disabling native GPU acceleration")
            set_property(CACHE GMX_GPU PROPERTY VALUE OFF)
            set(CUDA_NOTFOUND_AUTO ON)
        else ()
            # the user requested CUDA, but it wasn't found
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
mark_as_advanced(CUDA_BUILD_CUBIN CUDA_BUILD_EMULATION CUDA_SDK_ROOT_DIR CUDA_TOOLKIT_ROOT_DIR CUDA_VERBOSE_BUILD)

macro(gmx_gpu_setup)
    # set up nvcc options
    include(gmxManageNvccConfig)

    # Version info (semicolon used as line separator) for nvcc.
    get_nvcc_version_info()

    # Atomic operations used for polling wait for GPU
    # (to avoid the cudaStreamSynchronize + ECC bug).
    # ThreadMPI is now always included. Thus, we don't check for Atomics anymore here.

    # no OpenMP is no good!
    if(NOT GMX_OPENMP)
        message(WARNING "
    To use GPU acceleration efficiently, mdrun requires OpenMP multi-threading.
    With no OpenMP a single CPU core can be used with a GPU which is not optimal.
    Note that with MPI multiple processes can be forced to use a single GPU, but this
    typically inefficient. Note that you need to set both C and C++ compilers that
    support OpenMP (CC and CXX environment variables, respectively) when using GPUs.")
    endif()
endmacro()
