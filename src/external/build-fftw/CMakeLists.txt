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

# Download and build a suitable copy of FFTW.
# The GROMACS team won't distribute source or binaries linked to FFTW
# because we are choosing to be very clear about distributing only
# LGPL-licensed code, to suit requirements from our funding source.
#
# Input: FFTW variable contains the FFTW component to build,
#        either fftw or fftwf for double or single precision

string(TOUPPER "${FFTW}" UPPERFFTW)
string(TOLOWER "${FFTW}" LOWERFFTW)

set(GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION "" CACHE INTERNAL "Optimization flags for FFTW compilation")
if(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[[:space:]].*")
    message(FATAL_ERROR "An internal limitation of FFTW means GROMACS cannot build FFTW in a directory with whitespace in its name. Either use a system FFTW, build it yourself, or build GROMACS in a different location.")
endif()

if(NOT GMX_DOUBLE)
    set(GMX_BUILD_OWN_FFTW_PREC --enable-float)
endif()

# Always build a static lib, so it gets added to libmd and doesn't need to be installed
set(GMX_BUILD_OWN_FFTW_SHARED_FLAG --disable-shared --enable-static)
if ((CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" OR CMAKE_SYSTEM_PROCESSOR STREQUAL "aarch64") AND BUILD_SHARED_LIBS) # FFTW doesn't use -DPIC by default
    set(GMX_BUILD_OWN_FFTW_SHARED_FLAG ${GMX_BUILD_OWN_FFTW_SHARED_FLAG} --with-pic)
endif()

# Set library optimizations
set(_fftw_simd_support_level "")
if(${GMX_SIMD_ACTIVE} MATCHES "^(SSE|AVX)" AND APPLE)
    # OS X and --enable-avx causes compilation issues (fftw always picks gcc by default). It's
    # not an important enough performance loss to bother warning the
    # user about.
    set(_fftw_simd_support_level "--enable-sse2")
elseif(${GMX_SIMD_ACTIVE} MATCHES "^(SSE)")
    set(_fftw_simd_support_level "--enable-sse2")
elseif(${GMX_SIMD_ACTIVE} MATCHES "^(AVX)" AND NOT ${GMX_SIMD_ACTIVE} MATCHES "^(AVX_512)")
    # Testing shows FFTW configured with --enable-sse2 --enable-avx is
    # slightly faster on most architectures than --enable-sse2 alone.
    # Support for --enable-avx2 was only added in 3.3.5, but
    # configuring with it is at worst a warning, even on an earlier
    # version.
    # On platforms capable of AVX512 where we are building with AVX2,
    # enabling AVX512 risks clock-throttling the entire mdrun if
    # fftw happens to pick up an AVX512 kernel (which is not unlikely
    # as fftw tuning is known to produce highly varying results).
    set(_fftw_simd_support_level --enable-sse2;--enable-avx;--enable-avx2)
elseif(${GMX_SIMD_ACTIVE} MATCHES "^(AVX_512)")
    # MSVC, GCC < 4.9, Clang < 3.9 do not support AVX-512, so
    # we should not enable it there. FFTW does not support clang with
    # AVX-512, so we should not enable that either.
    if(MSVC OR (CMAKE_COMPILER_IS_GNUCC AND CMAKE_C_COMPILER_VERSION VERSION_LESS 4.9.0) OR
        (CMAKE_C_COMPILER_ID MATCHES "Clang" AND CMAKE_C_COMPILER_VERSION VERSION_LESS 3.9.0) OR
        (CMAKE_C_COMPILER_ID MATCHES "Clang" AND ${GMX_SIMD_ACTIVE} MATCHES "^(AVX_512)"))
        set(_fftw_simd_support_level --enable-sse2;--enable-avx;--enable-avx2)
    else()
        set(_fftw_simd_support_level --enable-sse2;--enable-avx;--enable-avx2;--enable-avx512)
    endif()
elseif(${GMX_SIMD_ACTIVE} MATCHES "^(VSX)")
    set(_fftw_simd_support_level --enable-vsx)
endif()
set(GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION ${_fftw_simd_support_level} CACHE INTERNAL "Optimization flags for FFTW compilation")

# Allow cross-compiles
if (TARGET_HOST)
    set(GMX_BUILD_OWN_FFTW_TARGET_HOST --host=${TARGET_HOST})
endif()

# Machinery for running the external project
set(EXTERNAL_FFTW_VERSION 3.3.8)
set(GMX_BUILD_OWN_FFTW_URL
    "http://www.fftw.org/fftw-${EXTERNAL_FFTW_VERSION}.tar.gz" CACHE STRING
    "URL from where to download fftw (use an absolute path when offline, adjust GMX_BUILD_OWN_FFTW_MD5 if downloading other version than ${EXTERNAL_FFTW_VERSION})")
set(GMX_BUILD_OWN_FFTW_MD5 8aac833c943d8e90d51b697b27d4384d CACHE STRING
    "Expected MD5 hash for the file at GMX_BUILD_OWN_FFTW_URL")
mark_as_advanced(GMX_BUILD_OWN_FFTW_URL GMX_BUILD_OWN_FFTW_MD5)

# The actual build target.
set(EXTERNAL_FFTW_BUILD_TARGET fftwBuild)
include(ExternalProject)
ExternalProject_add(${EXTERNAL_FFTW_BUILD_TARGET}
        URL "${GMX_BUILD_OWN_FFTW_URL}"
        URL_MD5 ${GMX_BUILD_OWN_FFTW_MD5}
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --libdir=<INSTALL_DIR>/lib --disable-fortran
        ${GMX_BUILD_OWN_FFTW_SHARED_FLAG} ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}
        ${GMX_BUILD_OWN_FFTW_PREC}
        ${GMX_BUILD_OWN_FFTW_TARGET_HOST})
ExternalProject_get_property(${EXTERNAL_FFTW_BUILD_TARGET} INSTALL_DIR)

string(REGEX REPLACE "fftw" "fftw3" FFTW_LIBNAME ${LOWERFFTW})
set(${UPPERFFTW}_LIBRARIES ${INSTALL_DIR}/lib/lib${FFTW_LIBNAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
set(${UPPERFFTW}_INCLUDE_DIRS ${INSTALL_DIR}/include PARENT_SCOPE)

add_library(gmxfftw STATIC IMPORTED GLOBAL)
set_target_properties(gmxfftw PROPERTIES IMPORTED_LOCATION ${${UPPERFFTW}_LIBRARIES})
set(${UPPERFFTW}_LIBRARIES gmxfftw PARENT_SCOPE)
add_dependencies(gmxfftw ${EXTERNAL_FFTW_BUILD_TARGET})

message(STATUS "The GROMACS-managed build of FFTW 3 will configure with the following optimizations: ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}")
