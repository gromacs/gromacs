set(CONTRIB_PROGRAMS 
     #add here any programs you want to compile

)

# Uncomment the next line to build OpenMM:
#option(GMX_OPENMM "Accelerated execution on GPUs through the OpenMM library (no longer supported" ON)

# At run time, you may need to set the environment variable
# OPENMM_PLUGIN_DIR=PATH_TO_GROMACS/openmm/lib/plugins
# to make things work

if(GMX_OPENMM)
    if(GMX_GPU)
        message(FATAL_ERROR "The OpenMM build is not compatible with the native GPU build")
    endif()

    enable_language(CXX)
    set (GMX_BINARY_SUFFIX "-openmm")
    set (GMX_LIBS_SUFFIX "_openmm")

#######################################################################
# Check for options incompatible with OpenMM build                    #
#######################################################################
    # we'll use the built-in fft to avoid unnecessary dependencies
    string(TOUPPER ${GMX_FFT_LIBRARY} GMX_FFT_LIBRARY)
    if(NOT ${GMX_FFT_LIBRARY} STREQUAL "FFTPACK")
        message(STATUS "No external FFT libraries needed for the OpenMM build, switching to fftpack!")
        set(GMX_FFT_LIBRARY "fftpack" CACHE STRING 
		"No external FFT libraries needed for the OpenMM build, switching to  fftpack!" FORCE)
    endif()
    if(GMX_MPI)
        message(FATAL_ERROR "The OpenMM build is not compatible with MPI!")
    endif(GMX_MPI)
    if(GMX_THREAD_MPI)
        message(STATUS "Thread-MPI not compatible with OpenMM, disabled!")
        set(GMX_THREAD_MPI OFF CACHE BOOL
		"Thread-MPI not compatible with OpenMM build, disabled!" FORCE)
    endif(GMX_THREAD_MPI)
    if(GMX_OPENMP)
        message(STATUS "OpenMP multithreading not compatible with OpenMM, disabled")
        set(GMX_OPENMP OFF CACHE BOOL
            "OpenMP multithreading not compatible with OpenMM, disabled!" FORCE)
    endif()
    if(GMX_SOFTWARE_INVSQRT)
        set(GMX_SOFTWARE_INVSQRT OFF CACHE STRING 
                "The OpenMM build does not need GROMACS software 1/sqrt!" FORCE)
    endif(GMX_SOFTWARE_INVSQRT)
    string(TOUPPER ${GMX_CPU_ACCELERATION} GMX_CPU_ACCELERATION)
    if(NOT GMX_CPU_ACCELERATION STREQUAL "NONE")
        message(STATUS "Switching off CPU-based acceleration, the OpenMM build does not support/need any!")
        set(GMX_CPU_ACCELERATION "None" CACHE STRING
            "Switching off CPU-based acceleration, the OpenMM build does not support/need any!" FORCE)
    endif()
    if(GMX_FAHCORE)
        message(FATAL_ERROR "The OpenMM build does not support FAH build!")
    endif(GMX_FAHCORE)
    if(GMX_DOUBLE)
        message(FATAL_ERROR  "The OpenMM-build does not support double precision calculations!")
    endif()

    set(CUDA_BUILD_EMULATION OFF)
    find_package(CUDA 3.1 REQUIRED)
    add_definitions(-DGMX_OPENMM)
    if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")    
        set(CUDA_VERBOSE_BUILD ON)
    endif()
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src/contrib)
    find_package(OpenMM) 

    # mark as advanced the unused variables
    mark_as_advanced(FORCE GMX_CPU_ACCELERATION GMX_MPI GMX_FFT_LIBRARY 
        GMX_QMMM_PROGRAM GMX_THREAD_MPI GMX_DOUBLE)

else(GMX_OPENMM)

     mark_as_advanced(CLEAR GMX_CPU_ACCELERATION GMX_MPI GMX_FFT_LIBRARY 
        GMX_QMMM_PROGRAM GMX_THREAD_MPI GMX_DOUBLE)

endif(GMX_OPENMM)

foreach(PROG ${CONTRIB_PROGRAMS})
        add_executable(${PROG} ${PROG}.c ${NGMX_COMMON_SOURCE})
	set_target_properties(${PROG} PROPERTIES OUTPUT_NAME "${PROG}${GMX_BINARY_SUFFIX}")
endforeach(PROG)
