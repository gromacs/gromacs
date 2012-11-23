# Manage the MPI setup, assuming that CMAKE_C_COMPILER is an MPI
# (wrapper) compiler.
if(GMX_MPI)
  if(GMX_THREAD_MPI)
    message(STATUS "MPI is not compatible with thread-MPI. Disabling thread-MPI.")
    set(GMX_THREAD_MPI OFF CACHE BOOL
        "Build a thread-MPI-based multithreaded version of GROMACS (not compatible with MPI)" FORCE)
  endif(GMX_THREAD_MPI)

  # Test the CMAKE_C_COMPILER for being an MPI (wrapper) compiler
  TRY_COMPILE(MPI_FOUND ${CMAKE_BINARY_DIR}
    "${CMAKE_SOURCE_DIR}/cmake/TestMPI.c"
    COMPILE_DEFINITIONS )

  # If CMAKE_C_COMPILER is not a MPI wrapper. Try to find MPI using cmake module as fall-back.
  # cmake < 2.8.5 not recommended for fall-back because of unreliability (redmine #851)
  if(NOT MPI_FOUND)
      if(CMAKE_VERSION VERSION_LESS "2.8.5")
          set(MPI_PREFIX MPI)
      else()
          set(MPI_PREFIX MPI_C)
      endif()
      find_package(MPI)
      if(${${MPI_PREFIX}_FOUND})
        set(GROMACS_C_FLAGS ${GROMACS_FLAGS} ${${MPI_PREFIX}_COMPILE_FLAGS})
        set(GROMACS_LINKER_FLAGS ${GROMACS_LINKER_FLAGS} ${${MPI_PREFIX}_LINK_FLAGS})
        include_directories(${${MPI_PREFIX}_INCLUDE_PATH})
        list(APPEND GMX_EXTRA_LIBRARIES ${${MPI_PREFIX}_LIBRARIES})
      endif()
      set(MPI_FOUND ${${MPI_PREFIX}_FOUND})
  endif()

  if(MPI_FOUND)
    include(gmxTestMPI_IN_PLACE)
    if (GMX_MPI_IN_PLACE)
      gmx_test_mpi_in_place(MPI_IN_PLACE_EXISTS)
    endif()

    # Find path of the mpi compilers
    if (${${MPI_PREFIX}_FOUND})
        get_filename_component(_mpi_c_compiler_path "${MPI_C_COMPILER}" PATH)
        get_filename_component(_mpiexec_path "${MPIEXEC}" PATH)
    else()
        get_filename_component(_cmake_c_compiler_path "${CMAKE_C_COMPILER}" PATH)
        get_filename_component(_cmake_cxx_compiler_path "${CMAKE_CXX_COMPILER}" PATH)
    endif()

    # Test for and warn about unsuitable MPI versions
    #
    # Execute the ompi_info binary with the full path of the compiler wrapper
    # found, otherwise we run the risk of false positives.
    find_file(MPI_INFO_BIN ompi_info
              HINTS ${_mpi_c_compiler_path} ${_mpiexec_path}
                    ${_cmake_c_compiler_path} ${_cmake_cxx_compiler_path}
              NO_DEFAULT_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)
    if (MPI_INFO_BIN)
      exec_program(${MPI_INFO_BIN}
        ARGS -v ompi full
        OUTPUT_VARIABLE OPENMPI_TYPE
        RETURN_VALUE OPENMPI_EXEC_RETURN)
      if(OPENMPI_EXEC_RETURN EQUAL 0)
        string(REGEX REPLACE ".*Open MPI: \([0-9]+\\.[0-9]*\\.?[0-9]*\).*" "\\1" OPENMPI_VERSION ${OPENMPI_TYPE})
        if(OPENMPI_VERSION VERSION_LESS "1.4.1")
          MESSAGE(WARNING
             "CMake found OpenMPI version ${OPENMPI_VERSION} on your system. "
             "There are known problems with GROMACS and OpenMPI version < 1.4.1. "
             "Please consider updating your OpenMPI if your MPI wrapper compilers "
             "are using the above OpenMPI version.")
        endif()
        unset(OPENMPI_VERSION)
        unset(OPENMPI_TYPE)
        unset(OPENMPI_EXEC_RETURN)
      endif()
    endif()
    unset(MPI_INFO_BIN CACHE)

    # Execute the mpiname binary with the full path of the compiler wrapper
    # found, otherwise we run the risk of false positives.
    find_file(MPINAME_BIN mpiname
              HINTS ${_mpi_c_compiler_path}
                    ${_cmake_c_compiler_path} ${_cmake_cxx_compiler_path}
              NO_DEFAULT_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)
    if (MPINAME_BIN)
      exec_program(${MPINAME_BIN}
        ARGS -n -v
        OUTPUT_VARIABLE MVAPICH2_TYPE
        RETURN_VALUE MVAPICH2_EXEC_RETURN)
      if(MVAPICH2_EXEC_RETURN EQUAL 0)
        string(REGEX MATCH "MVAPICH2" MVAPICH2_NAME ${MVAPICH2_TYPE})
        # Want to check for MVAPICH2 in case some other library supplies mpiname
        string(REGEX REPLACE "MVAPICH2 \([0-9]+\\.[0-9]*[a-z]?\\.?[0-9]*\)" "\\1" MVAPICH2_VERSION ${MVAPICH2_TYPE})
        if(${MVAPICH2_NAME} STREQUAL "MVAPICH2" AND MVAPICH2_VERSION VERSION_LESS "1.5")
          # This test works correctly even with 1.5a1
          MESSAGE(WARNING
             "CMake found MVAPICH2 version ${MVAPICH2_VERSION} on your system. "
             "There are known problems with GROMACS and MVAPICH2 version < 1.5. "
             "Please consider updating your MVAPICH2 if your MPI wrapper compilers "
             "are using the above MVAPICH2 version.")
       endif()
       unset(MVAPICH2_VERSION)
       unset(MVAPICH2_NAME)
       unset(MVAPICH2_TYPE)
       unset(MVAPICH2_EXEC_RETURN)
      endif()
    endif()
    unset(MPINAME_BIN CACHE)

  else(MPI_FOUND)
    if (CMAKE_VERSION VERSION_LESS "2.8.5")
      message(FATAL_ERROR
        "MPI support requested, but no MPI compiler found. Either set the "
        "C-compiler (CMAKE_C_COMPILER) to the MPI compiler (often called mpicc), "
        "or use a newer cmake version (>=2.8.5) which has improved MPI detection.")
    else()
      message(FATAL_ERROR
        "MPI support requested, but no MPI compiler found. Either set the "
        "C-compiler (CMAKE_C_COMPILER) to the MPI compiler (often called mpicc), "
        "or set the variables reported missing for MPI_C above.")
    endif()

  endif(MPI_FOUND)

  include(gmxTestCatamount)
  gmx_test_catamount(GMX_CRAY_XT3)
  if(GMX_CRAY_XT3)
    set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_CRAY_XT3")
    set(GMX_NO_SYSTEM 1)
    set(GMX_NO_NICE 1)
  endif(GMX_CRAY_XT3)

  set(GMX_LIB_MPI 1)
  set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_LIB_MPI")
endif(GMX_MPI)
