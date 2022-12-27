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

include_guard()

# Manage setup of the different FFT libraries we can use in Gromacs.
set(PKG_FFT "")
set(PKG_FFT_LIBS "")
if(${GMX_FFT_LIBRARY} STREQUAL "FFTW3")
    # ${FFTW} must be in upper case
    if(GMX_DOUBLE)
        set(FFTW "FFTW")
    else()
        set(FFTW "FFTWF")
    endif()

    if(GMX_BUILD_OWN_FFTW)

        if(MSVC)
            message(FATAL_ERROR "Cannot build FFTW3 automatically (GMX_BUILD_OWN_FFTW=ON) in Visual Studio")
        endif()
        if(CMAKE_GENERATOR STREQUAL "Ninja")
            message(FATAL_ERROR "Cannot build FFTW3 automatically (GMX_BUILD_OWN_FFTW=ON) with ninja")
        endif()

        add_subdirectory(src/external/build-fftw)
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

        if(NOT WIN32) # Detection doesn't work on Windows
          if ((${GMX_SIMD_ACTIVE} MATCHES "SSE" OR ${GMX_SIMD_ACTIVE} MATCHES "AVX") AND NOT ${FFTW}_HAVE_SIMD)
              set(FFT_WARNING_MESSAGE "The fftw library found is compiled without SIMD support, which makes it slow. Consider recompiling it or contact your admin")
          else()
              if(${GMX_SIMD_ACTIVE} MATCHES "AVX" AND NOT (${FFTW}_HAVE_SSE OR ${FFTW}_HAVE_SSE2))
                  # If we end up here we have an AVX Gromacs build, and
                  # FFTW with SIMD.
                  set(FFT_WARNING_MESSAGE "The FFTW library was compiled with neither --enable-sse nor --enable-sse2; those would have enabled SSE(2) SIMD instructions. This will give suboptimal performance. You should (re)compile the FFTW library with --enable-sse2 and --enable-avx (and --enable-avx2 or --enable-avx512 if supported).")
              endif()
          endif()
        endif()

        find_path(ARMPL_INCLUDE_DIR "armpl.h" HINTS ${${FFTW}_INCLUDE_DIRS}
            NO_DEFAULT_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH)
        mark_as_advanced(ARMPL_INCLUDE_DIR)
        if (ARMPL_INCLUDE_DIR)
            set(GMX_FFT_ARMPL_FFTW3 1)
            set(FFT_STATUS_MESSAGE "Using external FFT library - ARM Performance Library (FFTW3 compatibility mode)")
        else()
            set(FFT_STATUS_MESSAGE "Using external FFT library - FFTW3")
        endif()
    endif()
    if (NOT GMX_FFT_ARMPL_FFTW3)
        set(GMX_FFT_FFTW3 1)
    endif()

    set(FFT_LIBRARIES ${${FFTW}_LIBRARIES})
  elseif(${GMX_FFT_LIBRARY} STREQUAL "MKL")
    if (TEST_MKL)
        set(FIND_MKL_QUIETLY "QUIET")
    endif()
    set(MKL_THREADING sequential)
    if (GMX_PREFER_STATIC_LIBS)
        set(MKL_LINK static)
    endif()
    find_package(MKL CONFIG REQUIRED ${FIND_MKL_QUIETLY} PATHS $ENV{MKLROOT})
    set(FFT_LIBRARIES MKL::MKL)
    #work-around for https://gitlab.kitware.com/cmake/cmake/-/issues/23844
    get_target_property(_mkl_path MKL::MKL INTERFACE_INCLUDE_DIRECTORIES)
    list(REMOVE_ITEM CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "${_mkl_path}")
    # Check MKL works. If we were in a non-global scope, we wouldn't
    # have to play nicely.
    set(old_CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_LIBRARIES "${FFT_LIBRARIES}")

    check_function_exists(DftiCreateDescriptor TEST_MKL)

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

    # Hide internal MKL options
    get_cmake_property(_VARS VARIABLES)
    foreach (_VARNAME ${_VARS})
        if (_VARNAME MATCHES "^[Mm][Kk][Ll]_")
            mark_as_advanced(${_VARNAME})
        endif()
    endforeach()
    mark_as_advanced(FORCE ENABLE_BLAS95)
    mark_as_advanced(FORCE ENABLE_LAPACK95)
    mark_as_advanced(FORCE ENABLE_BLACS)
    mark_as_advanced(FORCE ENABLE_CDFT)
    mark_as_advanced(FORCE ENABLE_CPARDISO)
    mark_as_advanced(FORCE ENABLE_SCALAPACK)
    mark_as_advanced(FORCE ENABLE_OMP_OFFLOAD)

    set(FFT_STATUS_MESSAGE "Using external FFT library - Intel MKL")
elseif(${GMX_FFT_LIBRARY} STREQUAL "FFTPACK")
    set(GMX_FFT_FFTPACK 1)
    set(FFT_STATUS_MESSAGE "Using internal FFT library - fftpack")
else()
    gmx_invalid_option_value(GMX_FFT_LIBRARY)
endif()
gmx_check_if_changed(FFT_CHANGED GMX_FFT_LIBRARY)
if (FFT_CHANGED)
    if(FFT_WARNING_MESSAGE)
        message(WARNING "${FFT_WARNING_MESSAGE}")
    endif()
    message(STATUS "${FFT_STATUS_MESSAGE}")
endif()

# enable threaded fftw3 if we've found it
if(FFTW3_THREADS OR FFTW3F_THREADS)
    add_definitions(-DFFT5D_FFTW_THREADS)
endif()

