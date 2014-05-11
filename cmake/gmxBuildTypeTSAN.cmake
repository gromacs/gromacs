#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014, by the GROMACS development team, led by
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

# Custom build type "TSAN", to be used for compiling GROMACS 
# with clang 3.4 or gcc 4.8 (currently pre-release) with ThreadSanitizer
# (aka "TSan") turned on, so that the tests can be run to detect data races.
#
# The main advantage of the clang version is that there can be a
# suppressions file that acts at compile time, though there is no use of 
# that yet. The main advantage of the gcc version is that it can be used
# with a OpenMP (if gomp is compiled with --disable-linux-futex).
#
# Unfortunately, out of the box Thread-MPI provokes several false
# positives. One example is that tMPI_Event_t contains an atomic int
# field "sync" that is initialized before thread spawn. During a
# collective, this is atomically updated by the source thread, and the
# update is observed by sink threads in tMPI_Event_wait, which do a
# yield wait when no change has occured. This means the read can
# happen before the write (by design, whether or not the read is
# atomic), but the surrounding logic prevents action until the write
# has happened. There is no way for the sink thread(s) to avoid
# reading until the write has happened - that is the point of the
# implementation.
#
# This ought to be able to be suppressed, but my attempts to apply
# suppressions on individual functions don't suppress reporting of the
# race event. Applying the suppression to the whole thread-MPI library
# might work, but seems to defeat the point. We want to be able to
# detect mis-use of the primitives provided by thread-MPI.
#
# This means there needs to be a way for this build type to trigger
# the use of the generic mutex-based fallback implementation within
# thread-MPI.
#
# Later, if a blacklist is needed, use something like
# "-fsanitize-blacklist=${CMAKE_SOURCE_DIR}/cmake/thread-sanitizer.supp"
# TODO find a better home for this and other suppression files
set(_flags "-O1 -g -fsanitize=thread")

foreach(_language C CXX)

    string(REPLACE "X" "+" _human_readable_language ${_language})

    if (CMAKE_${_language}_COMPILER_ID MATCHES "GNU")
        set(CMAKE_${_language}_FLAGS_TSAN "${_flags} -pie -fPIE" CACHE STRING "${_human_readable_language} flags for thread sanitizer")
    else()
        set(CMAKE_${_language}_FLAGS_TSAN ${_flags} CACHE STRING "${_human_readable_language} flags for thread sanitizer")
    endif()
    mark_as_advanced(CMAKE_${_language}_FLAGS_TSAN)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    if (_cmake_build_type STREQUAL TSAN)
        set(TMPI_ATOMICS_DISABLED 1)
        set(TMPI_ATOMICS 0)
        if (NOT((CMAKE_${_language}_COMPILER_ID MATCHES "Clang" AND
                    CMAKE_${_language}_COMPILER_VERSION VERSION_GREATER 3.2.999)
             OR (CMAKE_${_language}_COMPILER_ID MATCHES "GNU" AND
                    CMAKE_${_language}_COMPILER_VERSION VERSION_GREATER 4.7.999)))
            message(FATAL_ERROR "The ThreadSanitizer build is only available with clang ${_human_readable_language} >=3.3 and gnu ${_human_readable_language} >=4.8.")
        endif()
    endif()

endforeach()
