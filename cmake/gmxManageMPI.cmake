#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

if (GMX_MPI)
    if (GMX_THREAD_MPI)
        message(STATUS "MPI is not compatible with thread-MPI. Disabling thread-MPI.")
        set(GMX_THREAD_MPI OFF CACHE BOOL
            "Build a thread-MPI-based multithreaded version of GROMACS (not compatible with MPI)" FORCE)
    endif ()
    set(GMX_LIB_MPI 1)
else ()
    set(GMX_LIB_MPI 0)
endif ()

# If we aren't going to use an MPI library, then don't search for one
if (NOT GMXAPI AND NOT GMX_LIB_MPI)
    return()
endif()

# CMake's FindMPI.cmake is not robust enough to cope with a broken MPI
# installation that a GROMACS user might not explicitly want to use,
# but which are searched for because an MPI library is an optional
# dependency of gmxapi even when GROMACS is not built with an MPI
# library.
if (NOT GMX_LIB_MPI AND NOT CMAKE_DISABLE_FIND_PACKAGE_MPI AND NOT MPI_ALREADY_SEARCHED)
    message(STATUS "GROMACS is being built without library MPI support (-DGMX_MPI=no). However "
        "MPI is potentially useful for the gmxapi Python API, so we will search for MPI anyway.  "
        "If this causes problems, disable the check with -DCMAKE_DISABLE_FIND_PACKAGE_MPI=on.")
endif()

# Manage the MPI setup.
# Note that we may want to execute tests or Python with MPI,
# even if we are not using an MPI-enabled GROMACS build.
set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
set(GMX_REQUIRED_MPI_COMPONENTS)
if (GMX_LIB_MPI OR GMXAPI)
    # If we are building GROMACS against an MPI library, we need the CXX component.
    # If the gmxapi interfaces are to be installed, we want to try to help client
    # software to find a compatible MPI toolchain, regardless of the libgromacs configuration.
    list(APPEND GMX_REQUIRED_MPI_COMPONENTS "CXX")
endif ()
if (GMX_LIB_MPI AND GMX_CP2K)
    list(APPEND GMX_REQUIRED_MPI_COMPONENTS "Fortran")
endif ()
# We don't require MPI components here because we report errors elsewhere
# when we can't find a required component, and the MPI target is optional
# in some build configurations (e.g. thread-MPI gmxapi installations).
if (MPI_ALREADY_SEARCHED)
    set(MPI_FIND_QUIETLY ON)
endif()
find_package(MPI COMPONENTS ${GMX_REQUIRED_MPI_COMPONENTS})
set(MPI_ALREADY_SEARCHED TRUE CACHE BOOL "True if a search for MPI has already been done")
mark_as_advanced(MPI_ALREADY_SEARCHED)

if (GMX_LIB_MPI)
    if (NOT MPI_CXX_FOUND)
        message(FATAL_ERROR
                "MPI support requested, but no suitable MPI compiler found. Either set the "
                "MPI_CXX_COMPILER to the MPI compiler wrapper (often called mpicxx or mpic++), "
                "set CMAKE_CXX_COMPILER to a default-MPI-enabled compiler, "
                "or set the variables reported missing for MPI_CXX above.")
    elseif (MPI_CXX_VERSION VERSION_LESS 2.0)
        message(FATAL_ERROR "MPI version 2.0 or higher is required. Please update your MPI library.")
    endif ()
    #TODO(#3672, #3776): These should be acquired through the MPI::MPI_CXX target.
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
    list(APPEND GMX_COMMON_LIBRARIES ${MPI_CXX_LIBRARIES})
endif ()

# Identify particular MPI implementations of interest (for compatibility checks).
if (MPI_CXX_FOUND)
    string(REGEX MATCH ".*Open MPI[:]? [v]?\([0-9]+\\.[0-9]*\\.?[0-9]*\).*" _openmpi_version ${MPI_CXX_LIBRARY_VERSION_STRING})
    if (_openmpi_version)
        string(REGEX REPLACE ".*Open MPI[:]? [v]?\([0-9]+\\.[0-9]*\\.?[0-9]*\).*" "\\1" OPENMPI_VERSION
               ${_openmpi_version})
    endif ()
    string(REGEX MATCH ".*MVAPICH2[:]? [v]?\([0-9]+\\.[0-9]*[a-z]?\\.?[0-9]*\).*" _mvapich2_version ${MPI_CXX_LIBRARY_VERSION_STRING})
    if (_mvapich2_version)
        string(REGEX REPLACE ".*MVAPICH2[:]? [v]?\([0-9]+\\.[0-9]*[a-z]?\\.?[0-9]*\).*" "\\1" MVAPICH2_VERSION
               ${_mvapich2_version})
    endif ()
    unset(_mvapich2_version)
    unset(_openmpi_version)
endif ()

# Test for and warn about unsuitable OpenMPI versions.
# TODO(#4093): Update tests with respect to required (compatible) OpenMPI versions.
if (GMX_LIB_MPI AND OPENMPI_VERSION)
    if (OPENMPI_VERSION VERSION_LESS "1.4.1")
        MESSAGE(WARNING
                "CMake found OpenMPI version ${OPENMPI_VERSION} on your system. "
                "There are known problems with GROMACS and OpenMPI version < 1.4.1. "
                "Please consider updating your OpenMPI if your MPI wrapper compilers "
                "are using the above OpenMPI version.")
    endif ()
    if (OPENMPI_VERSION VERSION_EQUAL "1.8.6")
        MESSAGE(WARNING
                "CMake found OpenMPI version ${OPENMPI_VERSION} on your system. "
                "This OpenMPI version is known to leak memory with GROMACS,"
                "please update to a more recent version. ")
    endif ()
    if (NOT MPI_FIND_QUIETLY)
        MESSAGE(STATUS "GROMACS library will use OpenMPI ${OPENMPI_VERSION}")
    endif ()
endif ()

# Test for and warn about unsuitable MPVAPICH2 versions
# TODO(#4093): Update tests with respect to required (compatible) MVAPICH2 versions.
if (GMX_LIB_MPI AND MVAPICH2_VERSION)
    if (MVAPICH2_VERSION VERSION_LESS "1.5")
        # This test works correctly even with 1.5a1
        MESSAGE(WARNING
                "CMake found MVAPICH2 version ${MVAPICH2_VERSION} on your system. "
                "There are known problems with GROMACS and MVAPICH2 version < 1.5. "
                "Please consider updating your MVAPICH2 if your MPI wrapper compilers "
                "are using the above MVAPICH2 version.")
    endif ()
    if (NOT MPI_FIND_QUIETLY)
        MESSAGE(STATUS "GROMACS library will use MVAPICH2 ${MVAPICH2_VERSION}")
    endif ()
endif ()

# Look for MPI process launchers that may be missed, especially if we didn't
# need to find the full MPI library build system support.
if (NOT MPIEXEC_EXECUTABLE)
    find_program(MPIEXEC
                 NAMES mpiexec mpirun lamexec srun aprun poe
                 HINTS ${MPI_HOME} $ENV{MPI_HOME}
                 PATH_SUFFIXES bin
                 DOC "Executable for running MPI programs.")

    set(MPIEXEC_EXECUTABLE "$MPIEXEC")
    set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
    set(MPIEXEC_PREFLAGS "" CACHE STRING "These flags will be directly before the executable that is being run by MPIEXEC.")
    set(MPIEXEC_POSTFLAGS "" CACHE STRING "These flags will come after all flags given to MPIEXEC.")
    set(MPIEXEC_MAX_NUMPROCS "2" CACHE STRING "Maximum number of processors available to run MPI applications.")
    mark_as_advanced(MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS MPIEXEC_POSTFLAGS MPIEXEC_MAX_NUMPROCS)
endif ()
