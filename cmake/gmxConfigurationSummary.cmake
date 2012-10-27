macro (_gmx_summary_append_messages TEXT_VAR)
    foreach (_line ${ARGN})
        set(${TEXT_VAR} "${${TEXT_VAR}}\n   ${_line}")
    endforeach()
endmacro()

macro (_gmx_summary_check_variable_added CACHE_VAR)
    list(FIND _summary_variables ${CACHE_VAR} _index)
    if (${_index} EQUAL -1)
        message(AUTHOR_WARNING "gmx_summary_add_variable(${CACHE_VAR}) missing")
        return()
    endif()
endmacro()

macro (gmx_summary_add_variable CACHE_VAR)
    list(APPEND _summary_variables ${CACHE_VAR})
endmacro()

macro (gmx_summary_remove_variable CACHE_VAR)
    list(REMOVE_ITEM _summary_variables ${CACHE_VAR})
endmacro()

function (gmx_summary_add_description CACHE_VAR)
    _gmx_summary_check_variable_added(${CACHE_VAR})
    set(_outputvar _description)
    set(_option FALSE)
    set(_description)
    foreach (arg ${ARGN})
        if ("${arg}" STREQUAL "OPTION")
            set(_option TRUE)
        elseif ("${arg}" STREQUAL "DESCRIPTION")
            set(_outputvar _description)
        elseif (_outputvar)
            list(APPEND ${_outputvar} ${arg})
        endif()
    endforeach()
    if (_option)
        if (${CACHE_VAR})
            set(_description "Build with ${_description}")
        else()
            set(_description "Build without ${_description}")
        endif()
    endif()
    set(_summary_${CACHE_VAR}_description ${_description} PARENT_SCOPE)
endfunction()

function (gmx_summary_add_message CACHE_VAR TYPE MSG1)
    _gmx_summary_check_variable_added(${CACHE_VAR})
    # TODO: The parameters need some massaging to make the warning/error
    # look nice.
    # RFC: Do we need to suppress these on subsequent runs? Under which conditions?
    if (${TYPE} STREQUAL "WARNING")
        message(WARNING ${MSG1} ${ARGN})
    elseif (${TYPE} STREQUAL "ERROR")
        message(SEND_ERROR ${MSG1} ${ARGN})
    endif()
    string(TOUPPER ${TYPE} _prefix)
    set(_text "${_prefix}: ${MSG1}")
    list(APPEND _summary_${CACHE_VAR}_messages ${_text} ${ARGN})
    set(_summary_${CACHE_VAR}_messages ${_summary_${CACHE_VAR}_messages} PARENT_SCOPE)
endfunction()

function (gmx_summary_print)
    list(REMOVE_DUPLICATES _summary_variables)
    message("Summary of configuration options:")
    foreach (CACHE_VAR ${_summary_variables})
        set(_description ${_summary_${CACHE_VAR}_description})
        set(_messages    ${_summary_${CACHE_VAR}_messages})
        set(_text " * ${CACHE_VAR}=${${CACHE_VAR}}: ${_description}")
        # RFC: Do we need to ensure that the same messages are reshown on the
        # next invocation of CMake?  Currently, the first run may show a
        # different message than subsequent runs (e.g., if GMX_OPENMP gets
        # disabled because the compiler doesn't support OpenMP).
        _gmx_summary_append_messages(_text ${_messages})
        message(${_text})
    endforeach()
endfunction()

macro (gmx_summary_initialize)
    option(GMX_CMAKE_SUMMARY
           "Print summary of configuration options at each CMake invocation"
           ON)
    mark_as_advanced(GMX_CMAKE_SUMMARY)
    gmx_summary_add_variable(GMX_OPENMM)
    gmx_summary_add_variable(GMX_DOUBLE)
    gmx_summary_add_variable(GMX_CPU_ACCELERATION)
    gmx_summary_add_variable(GMX_GPU)
    gmx_summary_add_variable(GMX_OPENMP)
    gmx_summary_add_variable(GMX_THREAD_MPI)
    gmx_summary_add_variable(GMX_MPI)
    gmx_summary_add_variable(GMX_FFT_LIBRARY)
endmacro()

macro (gmx_summary_process)
    if (GMX_CMAKE_SUMMARY)
        if (GMX_OPENMM)
            gmx_summary_add_description(GMX_OPENMM OPTION "GPU support using OpenMM")
            gmx_summary_add_message(GMX_OPENMM NOTE "GMX_GPU provides native GPU acceleration.")
        else()
            gmx_summary_remove_variable(GMX_OPENMM)
        endif()

        if (GMX_DOUBLE)
            gmx_summary_add_description(GMX_DOUBLE "Build in double precision")
            gmx_summary_add_message(GMX_DOUBLE NOTE "Double precision is significantly slower.")
        else()
            gmx_summary_add_description(GMX_DOUBLE "Build in single precision")
        endif()
        gmx_summary_add_description(GMX_GPU
            OPTION "native GPU acceleration for NVIDIA hardware")
        if (NOT GMX_GPU)
            gmx_summary_add_message(GMX_GPU WARNING
                "mdrun supports native GPU acceleration on NVIDIA hardware with"
                "compute capability >=2.0.  This requires the NVIDIA CUDA library."
                "CPU or GPU acceleration can be selected at runtime.  Unless you are sure"
                "you cannot make use of GPU acceleration, you should compile with"
                "GMX_GPU=ON.")
        endif()
        gmx_summary_add_description(GMX_OPENMP OPTION "OpenMP parallelization support")
        if (NOT GMX_OPENMP)
            gmx_summary_add_message(GMX_OPENMP NOTE
                "Compiling without OpenMP support might hurt your performance a lot,"
                "in particular with GPUs.")
            if (GMX_GPU)
                gmx_summary_add_message(GMX_OPENMP WARNING
                    "To use GPU acceleration efficiently, mdrun requires OpenMP multi-threading."
                    "With no OpenMP a single CPU core can be used with a GPU which is not optimal."
                    "Note that with MPI multiple processes can be forced to use a single GPU, but this"
                    "is typically inefficient. Note that you need to set both C and C++ compilers that"
                    "support OpenMP (CC and CXX environment variables, respectively) when using GPUs.")
            endif()
        endif()

        if (GMX_THREAD_MPI)
            gmx_summary_add_description(GMX_THREAD_MPI
                OPTION "thread MPI for intra-node parallelization")
            gmx_summary_remove_variable(GMX_MPI)
        else()
            gmx_summary_add_description(GMX_MPI
                OPTION "MPI parallelization support")
            gmx_summary_remove_variable(GMX_THREAD_MPI)
        endif()

        if (${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
            gmx_summary_add_description(GMX_FFT_LIBRARY
                "Build with FFTW3 FFT library")
        elseif (${GMX_FFT_LIBRARY} STREQUAL "MKL")
            gmx_summary_add_description(GMX_FFT_LIBRARY
                "Build with Intel MKL FFT library")
        elseif (${GMX_FFT_LIBRARY} STREQUAL "FFTPACK")
            gmx_summary_add_description(GMX_FFT_LIBRARY
                "Build with built-in fftpack FFT library")
            gmx_summary_add_message(GMX_FFT_LIBRARY WARNING
                "Note: Performance of simulations that use PME is poor with"
                "fftpack.  Consider installing FFTW3 or using Intel MKL if you"
                "are using Intel compiler.")
        endif()

        set(_acceleration_x86list_internal
            "NONE" "SSE2" "SSE4.1" "AVX_128_FMA" "AVX_256")
        set(_acceleration_x86list_user
            "None" "SSE2" "SSE4.1" "AVX128+FMA" "AVX256")
        list(FIND _acceleration_x86list_internal ${GMX_CPU_ACCELERATION}
             _acceleration_x86_index)
        if (${GMX_CPU_ACCELERATION} STREQUAL "None")
            gmx_summary_add_description(GMX_CPU_ACCELERATION
                "Build only generic kernels and no platform-specifc optimization")
        elseif (NOT _acceleration_x86_index EQUAL -1)
            list(GET _acceleration_x86list_user ${_acceleration_x86_index}
                 _acceleration_x86_user)
            gmx_summary_add_description(GMX_CPU_ACCELERATION
                "Build with ${_acceleration_x86_user}-optimized kernels and code")
            if (${GMX_SUGGESTED_CPU_ACCELERATION} STREQUAL "None")
                gmx_summary_add_message(GMX_CPU_ACCELERATION NOTE
                    "Could not detect best acceleration setting for the build host")
            else()
                list(FIND _acceleration_x86list_internal ${GMX_SUGGESTED_CPU_ACCELERATION}
                     _acceleration_x86_suggested_index)
                if (_acceleration_x86_index LESS _acceleration_x86_suggested_index)
                    gmx_summary_add_message(GMX_CPU_ACCELERATION NOTE
                        "The build CPU supports ${GMX_SUGGESTED_CPU_ACCELERATION}."
                        "The binary will not take advantage of this.")
                elseif (_acceleration_x86_index GREATER _acceleration_x86_suggested_index)
                    gmx_summary_add_message(GMX_CPU_ACCELERATION NOTE
                        "The build CPU does not support this level of acceleration."
                        "The binaries will not run on the build host.")
                endif()
            endif()
        endif()

        gmx_summary_print()
    endif()
endmacro()
