#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

# Manage setup of the different FFT libraries we can use in Gromacs.
set(PKG_FFT "")
set(PKG_FFT_LIBS "")
# Intel 11 and up makes life somewhat easy if you just want to use
# all their stuff. It's not easy if you only want some of their
# stuff...
set(MKL_MANUALLY FALSE)
if (GMX_FFT_LIBRARY STREQUAL "MKL" AND
    NOT (CMAKE_C_COMPILER_ID MATCHES "Intel" AND CMAKE_C_COMPILER_VERSION VERSION_GREATER "11"))
    # The user will have to provide the set of magic libraries in
    # MKL_LIBRARIES (see below), which we cache (non-advanced), so that they
    # don't have to keep specifying it, and can easily see that
    # CMake is still using that information.
    set(MKL_MANUALLY TRUE)
endif()
set(MKL_LIBRARIES_FORMAT_DESCRIPTION "Use full paths to library files, in the right order, and separated by semicolons.")
gmx_dependent_cache_variable(
    MKL_LIBRARIES
    "List of libraries for linking to MKL. ${MKL_LIBRARIES_FORMAT_DESCRIPTION}"
    STRING ""
    MKL_MANUALLY)
gmx_dependent_cache_variable(
    MKL_INCLUDE_DIR
    "Path to mkl.h (non-inclusive)."
    PATH ""
    MKL_MANUALLY)
if(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
    # ${FFTW} must be in upper case
    if(GMX_DOUBLE)
        set(FFTW "FFTW")
    else()
        set(FFTW "FFTWF")
    endif()

    if(GMX_BUILD_OWN_FFTW)

        if(WIN32)
            message(FATAL_ERROR "Cannot build FFTW3 automatically (GMX_BUILD_OWN_FFTW=ON) on Windows")
        endif()

        add_subdirectory(src/contrib/fftw)
        include_directories(BEFORE ${${FFTW}_INCLUDE_DIRS})
        # libgmxfftw is always built static, so libgromacs does not
        # have a dependency on anything, so PKG_FFT should be empty
        set(PKG_FFT "")
        set(FFT_STATUS_MESSAGE "Using external FFT library - FFTW3 build managed by GROMACS")
    else()
        string(TOLOWER "${FFTW}" LOWERFFTW)
        find_package(FFTW COMPONENTS ${LOWERFFTW})

        if(NOT ${FFTW}_FOUND)
            MESSAGE(FATAL_ERROR "Cannot find FFTW 3 (with correct precision - libfftw3f for mixed-precision GROMACS or libfftw3 for double-precision GROMACS). Either choose the right precision, choose another FFT(W) library (-DGMX_FFT_LIBRARY), enable the advanced option to let GROMACS build FFTW 3 for you (-DGMX_BUILD_OWN_FFTW=ON), or use the really slow GROMACS built-in fftpack library (-DGMX_FFT_LIBRARY=fftpack).")
        endif()

        set(PKG_FFT "${${FFTW}_PKG}")
        include_directories(SYSTEM ${${FFTW}_INCLUDE_DIRS})

        if ((${GMX_SIMD_ACTIVE} MATCHES "SSE" OR ${GMX_SIMD_ACTIVE} MATCHES "AVX") AND NOT ${FFTW}_HAVE_SIMD)
            message(WARNING "The fftw library found is compiled without SIMD support, which makes it slow. Consider recompiling it or contact your admin")
        else()
            if(${GMX_SIMD_ACTIVE} MATCHES "AVX" AND NOT (${FFTW}_HAVE_SSE OR ${FFTW}_HAVE_SSE2))
                # If we end up here we have an AVX Gromacs build, and
                # FFTW with SIMD.
                message(WARNING "The FFTW library was compiled with neither --enable-sse nor --enable-sse2; those would have enabled SSE(2) SIMD instructions. This will give suboptimal performance. You should (re)compile the FFTW library with --enable-sse2 and --enable-avx (and --enable-avx2 or --enable-avx512 if supported).")
            endif()
        endif()
        set(FFT_STATUS_MESSAGE "Using external FFT library - FFTW3")
    endif()

    set(FFT_LIBRARIES ${${FFTW}_LIBRARIES})
    set(GMX_FFT_FFTW3 1)
elseif(${GMX_FFT_LIBRARY} STREQUAL "MKL")
    # Intel 11 and up makes life somewhat easy if you just want to use
    # all their stuff. It's not easy if you only want some of their
    # stuff...
    if (NOT MKL_MANUALLY)
        # The next line takes care of everything for MKL
        if (WIN32)
            # This works according to the Intel MKL 10.3 for Windows
            # docs, but on Jenkins Win2k8, icl tries to interpret it
            # as a file. Shrug.
            set(FFT_LINKER_FLAGS "/Qmkl:sequential")
        else()
            set(FFT_LINKER_FLAGS "-mkl=sequential")
        endif()
        # Some versions of icc require this in order that mkl.h can be
        # found at compile time.
        set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} ${FFT_LINKER_FLAGS}")
        set(EXTRA_CXX_FLAGS "${EXTRA_CXX_FLAGS} ${FFT_LINKER_FLAGS}")

        set(MKL_ERROR_MESSAGE "Make sure you have configured your compiler so that ${FFT_LINKER_FLAGS} will work.")
    else()
        include_directories(SYSTEM ${MKL_INCLUDE_DIR})
        set(FFT_LIBRARIES "${MKL_LIBRARIES}")
        set(MKL_ERROR_MESSAGE "The include path to mkl.h in MKL_INCLUDE_DIR, and the link libraries in MKL_LIBRARIES=${MKL_LIBRARIES} need to match what the MKL documentation says you need for your system: ${MKL_LIBRARIES_FORMAT_DESCRIPTION}")
        # Convert the semi-colon separated list to a list of
        # command-line linker arguments so that code using our
        # pkgconfig setup can use it.
        string(REGEX REPLACE ";" " " PKG_FFT_LIBS "${MKL_LIBRARIES}")
    endif()

    # Check MKL works. If we were in a non-global scope, we wouldn't
    # have to play nicely.
    set(old_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FFT_LINKER_FLAGS}")
    set(old_CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_LIBRARIES "${FFT_LIBRARIES}")

    check_function_exists(DftiCreateDescriptor TEST_MKL)

    set(CMAKE_REQUIRED_FLAGS "${old_CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_LIBRARIES "${old_CMAKE_REQUIRED_LIBRARIES}")

    if(NOT TEST_MKL)
        # Hack to help the user vary MKL settings until they work.
        # TODO Make this logic more useful.
        unset(TEST_MKL CACHE)
        message(FATAL_ERROR "Linking with MKL was requested, but was not successful: ${MKL_ERROR_MESSAGE}")
    endif()

    # Set variables to signal that we have MKL available and should use it for FFTs.
    set(GMX_FFT_MKL 1)
    set(HAVE_LIBMKL 1)

    set(FFT_STATUS_MESSAGE "Using external FFT library - Intel MKL")
elseif(${GMX_FFT_LIBRARY} STREQUAL "FFTPACK")
    set(GMX_FFT_FFTPACK 1)
    set(FFT_STATUS_MESSAGE "Using internal FFT library - fftpack")
else()
    gmx_invalid_option_value(GMX_FFT_LIBRARY)
endif()
gmx_check_if_changed(FFT_CHANGED GMX_FFT_LIBRARY)
if (FFT_CHANGED)
    message(STATUS "${FFT_STATUS_MESSAGE}")
endif()

# enable threaded fftw3 if we've found it
if(FFTW3_THREADS OR FFTW3F_THREADS)
    add_definitions(-DFFT5D_FFTW_THREADS)
endif()

