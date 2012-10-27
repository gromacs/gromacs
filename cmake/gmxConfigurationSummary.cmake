macro (_append_lines TEXT_VAR)
    foreach (line ${ARGN})
        set(${TEXT_VAR} "${${TEXT_VAR}}\n   ${line}")
    endforeach()
endmacro()

macro (gmx_add_to_summary CACHE_VAR)
    set(_outputvar _description)
    set(_option FALSE)
    set(_description)
    set(_warning)
    set(_warning_on)
    set(_warning_off)
    foreach (arg ${ARGN})
        if ("${arg}" STREQUAL "OPTION")
            set(_option TRUE)
        elseif ("${arg}" STREQUAL "DESCRIPTION")
            set(_outputvar _description)
        elseif ("${arg}" STREQUAL "WARNING")
            set(_outputvar _warning)
        elseif ("${arg}" STREQUAL "WARNING_ON")
            set(_outputvar _warning_on)
        elseif ("${arg}" STREQUAL "WARNING_OFF")
            set(_outputvar _warning_off)
        elseif (_outputvar)
            list(APPEND ${_outputvar} ${arg})
        endif()
    endforeach()
    set(_text " * ${CACHE_VAR}=${${CACHE_VAR}}: ")
    if (_option)
        if (${CACHE_VAR})
            set(_text "${_text}Build with ${_description}")
        else()
            set(_text "${_text}Build without ${_description}")
        endif()
    else()
        set(_text "${_text}${_description}")
    endif()
    _append_lines(_text ${_warning})
    if (${CACHE_VAR})
        _append_lines(_text ${_warning_on})
    else()
        _append_lines(_text ${_warning_off})
    endif()
    message(${_text})
endmacro()

option(GMX_CMAKE_SUMMARY
       "Print summary of configuration options at each CMake invocation)"
       ON)
mark_as_advanced(GMX_CMAKE_SUMMARY)
if (GMX_CMAKE_SUMMARY)
    message("Summary of configuration options:")

    if (GMX_DOUBLE)
        set(_precision "double")
    else()
        set(_precision "single")
    endif()
    gmx_add_to_summary(GMX_DOUBLE "Build in ${_precision} precision"
        WARNING_ON "Note: Double precision is significantly slower.")

    if (GMX_OPENMM)
        gmx_add_to_summary(GMX_OPENMM OPTION "GPU support using OpenMM"
            WARNING_ON "Note: GMX_GPU provides native GPU acceleration.")
    endif ()

    if (NOT GMX_OPENMM)
        gmx_add_to_summary(GMX_GPU
            OPTION "native GPU acceleration for NVIDIA hardware"
            WARNING_OFF
            "Note: mdrun supports native GPU acceleration on NVIDIA hardware with"
            "compute capability >=2.0.  This requires the NVIDIA CUDA library."
            "CPU or GPU acceleration can be selected at runtime.  Unless you are sure"
            "you cannot make use of GPU acceleration, you should compile with"
            "GMX_GPU=ON.")

        set(_openmp_warning
            "Compiling without OpenMP support might hurt your performance a lot,"
            "in particular with GPUs.")
        if (GMX_GPU)
            list(APPEND _openmp_warning
                "In order to use GPU acceleration efficiently, mdrun requires OpenMP"
                "multithreading.  Without OpenMP only a single CPU core per GPU can"
                "be used which is suboptimal.  Note that with MPI multiple processes"
                "can be forced to use a single GPU, but this typically inefficient.")
        endif()
        gmx_add_to_summary(GMX_OPENMP OPTION "OpenMP parallelization support"
            WARNING_OFF ${_openmp_warning})

        if (GMX_THREAD_MPI)
            gmx_add_to_summary(GMX_THREAD_MPI
                OPTION "thread MPI for intra-node parallelization")
        else()
            gmx_add_to_summary(GMX_MPI OPTION "MPI parallelization support")
        endif()

        set(_fft_description)
        set(_fft_warning)
        if (${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
            set(_fft_description "Build with FFTW3 FFT library")
        elseif (${GMX_FFT_LIBRARY} STREQUAL "MKL")
            set(_fft_description "Build with Intel MKL FFT library")
        elseif (${GMX_FFT_LIBRARY} STREQUAL "FFTPACK")
            set(_fft_description "Build with built-in fftpack FFT library")
            set(_fft_warning
                "Note: Performance of simulations that use PME is poor with"
                "fftpack.  Consider installing FFTW3 or using Intel MKL if you"
                "are using Intel compiler.")
        endif()
        gmx_add_to_summary(GMX_FFT_LIBRARY ${_fft_description}
            WARNING ${_fft_warning})

        set(_acceleration_x86list_internal
            "NONE" "SSE2" "SSE4.1" "AVX_128_FMA" "AVX_256")
        set(_acceleration_x86list_user
            "None" "SSE2" "SSE4.1" "AVX128+FMA" "AVX256")
        list(FIND _acceleration_x86list_internal ${GMX_ACCELERATION}
             _acceleration_x86_index)
        set(_acceleration_description)
        set(_acceleration_warning)
        if (${GMX_ACCELERATION} STREQUAL "None")
            set(_acceleration_description
                "Build only generic kernels and no platform-specifc optimization")
        elseif (NOT _acceleration_x86_index EQUAL -1)
            list(GET _acceleration_x86list_user ${_acceleration_x86_index}
                 _acceleration_x86_user)
            set(_acceleration_description
                "Build with ${_acceleration_x86_user}-optimized kernels and code")
            if (${GMX_SUGGESTED_ACCELERATION} STREQUAL "None")
                set(_acceleration_warning
                    "Note: Could not detect best acceleration setting for the build host")
            else()
                list(FIND _acceleration_x86list_internal ${GMX_SUGGESTED_ACCELERATION}
                     _acceleration_x86_suggested_index)
                if (_acceleration_x86_index LESS _acceleration_x86_suggested_index)
                    set(_acceleration_warning
                        "Note: The build CPU supports ${GMX_SUGGESTED_ACCELERATION}."
                        "The binary will not take advantage of this.")
                elseif (_acceleration_x86_index GREATER _acceleration_x86_suggested_index)
                    set(_acceleration_warning
                        "Note: The build CPU does not support this level of acceleration."
                        "The binaries will not run on the build host.")
                endif()
            endif()
        elseif (${GMX_ACCELERATION} STREQUAL "BlueGene")
            set(_acceleration_description
                "Build with BlueGene-optimized kernels and code")
        elseif (${GMX_ACCELERATION} STREQUAL "Power6")
            set(_acceleration_description
                "Build with Power6-optimized kernels")
        elseif (${GMX_ACCELERATION} STREQUAL "Fortran")
            set(_acceleration_description
                "Build with Fortran kernels")
        endif()
        gmx_add_to_summary(GMX_ACCELERATION ${_acceleration_description}
            WARNING ${_acceleration_warning})
    endif()
endif()
