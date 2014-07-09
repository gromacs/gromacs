#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015, by the GROMACS development team, led by
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

# This file implements the custom build type "MSAN", to be used to run
# the Memory Sanitizer checker in clang. The code here follows the
# approaches documented at
# https://code.google.com/p/memory-sanitizer/wiki/LibcxxHowTo and
# https://code.google.com/p/memory-sanitizer/wiki/MemorySanitizer. In
# particular, the dependency stack must also be built with
# MemorySanitizer - notably including the C++ standard library, and
# notably excluding the C and math libraries.
#
# CMake cache variable GMX_MSAN_PATH should be set (e.g. on the
# command line) to point to the location where the dependencies that
# have been built with MemorySanitizer support can be found. This will
# be preprended to preprocessor and linker search paths. The MSAN
# build type is highly unlikely to work usefully without this variable
# set appropriately, and the dependencies installed there. In
# particular, suitable zlib and libxml2 need to be installed there.
#
# Where GROMACS can optionally use bundled dependency code, such as
# fftpack and XDR, only the bundled code is supported. MSAN build type
# does not support static linking at this time.
#
# This build type needs to set linker flags in order to function, and
# CMake does not directly support linker flags that are specific to a
# build type (except on Windows). Thus, it is not clear whether
# changing the build type to/from MSAN, or any other modification to
# CMAKE_*_LINKER_FLAGS will work well. If in doubt, clean the build
# tree and start again.
function(gmxManageMsanBuild)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)

    if (NOT "${_cmake_build_type}" STREQUAL "MSAN")
        # Nothing needs to be done now
        return()
    endif()

    set(_flags "-O2 -g -fsanitize=memory -fno-omit-frame-pointer")

    foreach(_language C CXX)
        string(REPLACE "X" "+" _human_readable_language ${_language})
        if (CMAKE_${_language}_COMPILER_ID MATCHES "Clang" AND
            NOT CMAKE_${_language}_COMPILER_VERSION VERSION_LESS 3.4)
            # Can't do Memory Sanitizer build with this compiler, so don't
            # set up flags for it
	    if(${_language} MATCHES CXX)
                set(_language_flags "${_flags} -stdlib=libc++")
		if(GMX_MSAN_PATH)
                    set(_language_flags "${_language_flags} -I${GMX_MSAN_PATH}/include -I${GMX_MSAN_PATH}/include/c++/v1")
		endif()
	    else()
                set(_language_flags "${_flags}")
            endif()
            set(CMAKE_${_language}_FLAGS_MSAN ${_language_flags} CACHE STRING "${_human_readable_language} flags for Memory Sanitizer")
            mark_as_advanced(CMAKE_${_language}_FLAGS_MSAN)
	else()
            message(FATAL_ERROR "The Memory Sanitizer build is only available with clang ${_human_readable_language} compiler >= 3.4, but it was ${CMAKE_${_language}_COMPILER_ID} and ${CMAKE_${_language}_COMPILER_VERSION}.")
        endif()
    endforeach()

    # Per-build-type linker flags like CMAKE_EXE_LINKER_FLAGS_MSAN
    # only seem to be supported by CMake on Windows, so we can't use
    # that.
    set(_linker_flags "-fsanitize=memory -stdlib=libc++")
    if(GMX_MSAN_PATH)
        set(_linker_flags "${_linker_flags} -L${GMX_MSAN_PATH}/lib -lc++abi -I${GMX_MSAN_PATH}/include -I${GMX_MSAN_PATH}/include/c++/v1 -Wl,-rpath,${GMX_MSAN_PATH}/lib")
    endif()
    set(CMAKE_EXE_LINKER_FLAGS "${_linker_flags} ${CMAKE_EXE_LINKER_FLAGS}" PARENT_SCOPE)
    set(CMAKE_SHARED_LINKER_FLAGS "${_linker_flags} ${CMAKE_SHARED_LINKER_FLAGS}" PARENT_SCOPE)
#   Static linking is not supported for this build type
#    set(CMAKE_STATIC_LINKER_FLAGS ${_linker_flags} CACHE STRING "Linker flags for Memory Sanitizer" FORCE)

    set(GMX_FFT_LIBRARY "fftpack" CACHE STRING "Only fftpack is supported with MSAN")
    set(GMX_SYSTEM_XDR OFF CACHE STRING "Only internal XDR is supported with MSAN")
endfunction()

gmxManageMsanBuild()
