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
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND BUILD_SHARED_LIBS) # FFTW doesn't use -DPIC by default
    set(GMX_BUILD_OWN_FFTW_SHARED_FLAG ${GMX_BUILD_OWN_FFTW_SHARED_FLAG} --with-pic)
endif()

# Set library optimizations
set(_fftw_simd_support_level "")
if(${GMX_SIMD} MATCHES "^(SSE|AVX)" AND APPLE)
    # OS X and --enable-avx causes compilation issues (fftw always picks gcc by default). It's
    # not an important enough performance loss to bother warning the
    # user about.
    set(_fftw_simd_support_level "--enable-sse2")
elseif(${GMX_SIMD} MATCHES "^(SSE)")
    set(_fftw_simd_support_level "--enable-sse2")
elseif(${GMX_SIMD} MATCHES "^(AVX)")
    # Testing shows FFTW configured with --enable-sse2 --enable-avx is
    # slightly faster on most architectures than --enable-sse2 alone
    set(_fftw_simd_support_level --enable-sse2;--enable-avx)
endif()
set(GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION ${_fftw_simd_support_level} CACHE INTERNAL "Optimization flags for FFTW compilation")

# Allow cross-compiles
if (TARGET_HOST)
    set(GMX_BUILD_OWN_FFTW_TARGET_HOST --host=${TARGET_HOST})
endif()

# Machinery for running the external project
set(EXTERNAL_FFTW_VERSION 3.3.4)
# cmake make eats slashes //// -> //
set(GMX_BUILD_OWN_FFTW_URL
    "http:////www.fftw.org/fftw-${EXTERNAL_FFTW_VERSION}.tar.gz" CACHE PATH
    "URL from where to download fftw (use an absolute path when offline, adjust GMX_BUILD_OWN_FFTW_MD5 if downloading other version than ${EXTERNAL_FFTW_VERSION})")
set(GMX_BUILD_OWN_FFTW_MD5 2edab8c06b24feeb3b82bbb3ebf3e7b3 CACHE STRING
    "Expected MD5 hash for the file at GMX_BUILD_OWN_FFTW_URL")
mark_as_advanced(GMX_BUILD_OWN_FFTW_URL GMX_BUILD_OWN_FFTW_MD5)

# ExternalProject at least up to CMake 3.0 prints a confusing error message if
# download fails when MD5 verification is enabled.  So we manage the download
# ourselves so that MD5 sum is not verified there, and then pass a local file
# as the URL to ExternalProject.  This way, ExternalProject still verifies the
# MD5 sum with a proper message if that fails.
set(url "${GMX_BUILD_OWN_FFTW_URL}")
# Determine whether we are actually downloading (this matches the conditions in
# ExternalProject).  ExternalProject works as expected if passed a local file.
set(is_download TRUE)
if (IS_DIRECTORY "${url}" OR "${url}" MATCHES "^file://" OR NOT "${url}" MATCHES "^[a-z]+://")
    set(is_download FALSE)
endif()
if (is_download)
    # For simplicity, don't try to extract the file name from the URL, but use
    # a hard-coded value.
    set(remote_url "${GMX_BUILD_OWN_FFTW_URL}")
    set(local_path "${CMAKE_CURRENT_BINARY_DIR}/fftw.tar.gz")
    set(url ${local_path})
    # Write a script to do our own download step (mimics what ExternalProject
    # would do, but without MD5 sum verification at this step).
    set(download_script ${CMAKE_CURRENT_BINARY_DIR}/fftw-download.cmake)
    configure_file(fftw-download.cmake.cmakein ${download_script} @ONLY)
endif()

# The actual build target.
set(EXTERNAL_FFTW_BUILD_TARGET fftwBuild)
include(ExternalProject)
ExternalProject_add(${EXTERNAL_FFTW_BUILD_TARGET}
        URL "${url}" URL_MD5 ${GMX_BUILD_OWN_FFTW_MD5}
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --libdir=<INSTALL_DIR>/lib --disable-fortran
        ${GMX_BUILD_OWN_FFTW_SHARED_FLAG} ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}
        ${GMX_BUILD_OWN_FFTW_PREC}
        ${GMX_BUILD_OWN_FFTW_TARGET_HOST})
# Add a custom step to do our own download if that is necessary.
if (is_download)
    ExternalProject_add_step(${EXTERNAL_FFTW_BUILD_TARGET} pre-download
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/fftw-download.cmake
            DEPENDERS download)
endif()

ExternalProject_get_property(${EXTERNAL_FFTW_BUILD_TARGET} INSTALL_DIR)

string(REGEX REPLACE "fftw" "fftw3" FFTW_LIBNAME ${LOWERFFTW})
set(${UPPERFFTW}_LIBRARIES ${INSTALL_DIR}/lib/lib${FFTW_LIBNAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
set(${UPPERFFTW}_INCLUDE_DIRS ${INSTALL_DIR}/include PARENT_SCOPE)

add_library(gmxfftw STATIC IMPORTED GLOBAL)
set_target_properties(gmxfftw PROPERTIES IMPORTED_LOCATION ${${UPPERFFTW}_LIBRARIES})
set(${UPPERFFTW}_LIBRARIES gmxfftw PARENT_SCOPE)
add_dependencies(gmxfftw ${EXTERNAL_FFTW_BUILD_TARGET})

message(STATUS "The GROMACS-managed build of FFTW 3 will configure with the following optimizations: ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}")
