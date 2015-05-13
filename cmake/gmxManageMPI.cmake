#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

# Manage the MPI setup, assuming that CMAKE_C_COMPILER is an MPI
# (wrapper) compiler.
if(GMX_MPI)
  if(GMX_THREAD_MPI)
    message(STATUS "MPI is not compatible with thread-MPI. Disabling thread-MPI.")
    set(GMX_THREAD_MPI OFF CACHE BOOL
        "Build a thread-MPI-based multithreaded version of GROMACS (not compatible with MPI)" FORCE)
  endif()

  # Test the CMAKE_C_COMPILER for being an MPI (wrapper) compiler
  TRY_COMPILE(MPI_FOUND ${CMAKE_BINARY_DIR}
    "${CMAKE_SOURCE_DIR}/cmake/TestMPI.c"
    COMPILE_DEFINITIONS )

  # If CMAKE_C_COMPILER is not a MPI wrapper. Try to find MPI using cmake module as fall-back.
  if(NOT MPI_FOUND)
      find_package(MPI)
      if(MPI_C_FOUND)
        set(MPI_COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
        set(MPI_LINKER_FLAGS ${MPI_C_LINK_FLAGS})
        include_directories(${MPI_C_INCLUDE_PATH})
        list(APPEND GMX_EXTRA_LIBRARIES ${MPI_C_LIBRARIES})
      endif()
      set(MPI_FOUND ${MPI_C_FOUND})
  else()
      # The following defaults are based on FindMPI.cmake in cmake
      # 3.1.2. (That package does not actually do any detection of the
      # flags, but if it ever does then we should re-visit how we use
      # the package.) If we are compiling with an MPI wrapper
      # compiler, then MPI_FOUND will be set above, and will mean that
      # none of these cache variables are populated by the package. We
      # need to do it manually so that test drivers can work using the
      # standard machinery for CMake + FindMPI.cmake.  Users will need
      # to set these to suit their MPI setup in order for tests to
      # work.

      find_program(MPIEXEC
          NAMES mpiexec mpirun lamexec srun aprun poe
          HINTS ${MPI_HOME} $ENV{MPI_HOME}
          PATH_SUFFIXES bin
          DOC "Executable for running MPI programs.")

      set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
      set(MPIEXEC_PREFLAGS     ""    CACHE STRING "These flags will be directly before the executable that is being run by MPIEXEC.")
      set(MPIEXEC_POSTFLAGS    ""    CACHE STRING "These flags will come after all flags given to MPIEXEC.")
      set(MPIEXEC_MAX_NUMPROCS "2"   CACHE STRING "Maximum number of processors available to run MPI applications.")
      mark_as_advanced(MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS MPIEXEC_POSTFLAGS MPIEXEC_MAX_NUMPROCS)
  endif()

  if(MPI_FOUND)
    include(gmxTestMPI_IN_PLACE)
    if (GMX_MPI_IN_PLACE)
      gmx_test_mpi_in_place(MPI_IN_PLACE_EXISTS)
    endif()

    # Find path of the mpi compilers
    if (${MPI_C_FOUND})
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

    # Using find_file() runs the CMake standard module
    # GetPrerequisites.cmake, which adds the file_cmd
    # variable to the top-level CMake namespace. This is
    # fixed in CMake 2.8.10. Meanwhile, clean up for it.
    if(CMAKE_VERSION VERSION_LESS "2.8.10")
        mark_as_advanced(file_cmd)
    endif()

  else()
      message(FATAL_ERROR
        "MPI support requested, but no MPI compiler found. Either set the "
        "C-compiler (CMAKE_C_COMPILER) to the MPI compiler (often called mpicc), "
        "or set the variables reported missing for MPI_C above.")
  endif()

  set(GMX_LIB_MPI 1)
endif()
