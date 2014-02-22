#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

# Testing shows FFTW configured with --enable-avx is never better than --enable-sse2, so we do the latter always.
if(${GMX_SIMD} MATCHES "^(SSE|AVX)")
    set(GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION --enable-sse2 CACHE INTERNAL "Optimization flags for FFTW compilation")
endif()

# Allow cross-compiles
if (TARGET_HOST)
    set(GMX_BUILD_OWN_FFTW_TARGET_HOST --host=${TARGET_HOST})
endif()

# Machinery for running the external project
set(EXTERNAL_FFTW_VERSION 3.3.3)
# cmake make eats slashes //// -> //
set(GMX_BUILD_OWN_FFTW_URL "http:////www.fftw.org/fftw-${EXTERNAL_FFTW_VERSION}.tar.gz" CACHE PATH "URL from where to download fftw, (use an absolute path when offline)")
mark_as_advanced(GMX_BUILD_OWN_FFTW_URL)
set(EXTERNAL_FFTW_MD5SUM 0a05ca9c7b3bfddc8278e7c40791a1c2)
set (EXTERNAL_FFTW_BUILD_TARGET fftwBuild)
include(ExternalProject)
# TODO in master branch - show this warning only on the first run
# by using gmx_check_if_changed_result from I21b791ab8e4f3 when
# that becomes available
message(WARNING "The GROMACS build will download FFTW ${EXTERNAL_FFTW_VERSION} as requested, but it will not know the file it receives is correct. If you now use\nmake\n GROMACS will build and link to FFTW anyway, but there is a possible security risk if you execute a GROMACS tool that calls this library. Instead, you can use\nmake ${EXTERNAL_FFTW_BUILD_TARGET}\n to do just the download and build of FFTW, and then run\nmd5sum src/contrib/fftw/${EXTERNAL_FFTW_BUILD_TARGET}-prefix/src/fftw-${EXTERNAL_FFTW_VERSION}.tar.gz\nto see if it matches ${EXTERNAL_FFTW_MD5SUM}. If so, everything is OK and you should use \nmake\n to proceed with the rest of the GROMACS build. Alternatively, you could stop using GMX_BUILD_OWN_FFTW, and instead follow the GROMACS installation instructions to build FFTW yourself.")
# TODO if/when CMake fixes http://www.cmake.org/Bug/view.php?id=14330
# (ie. at least version > 2.8.11.2), consider reverting to using an
# md5sum check to avoid needing the above warning
    ExternalProject_add(${EXTERNAL_FFTW_BUILD_TARGET}
        URL "${GMX_BUILD_OWN_FFTW_URL}"
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --libdir=<INSTALL_DIR>/lib --disable-fortran
        ${GMX_BUILD_OWN_FFTW_SHARED_FLAG} ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}
        ${GMX_BUILD_OWN_FFTW_PREC}
        ${GMX_BUILD_OWN_FFTW_TARGET_HOST})
externalproject_get_property(${EXTERNAL_FFTW_BUILD_TARGET} INSTALL_DIR)

string(TOUPPER "${FFTW}" UPPERFFTW)
string(TOLOWER "${FFTW}" LOWERFFTW)
string(REGEX REPLACE "fftw" "fftw3" FFTW_LIBNAME ${LOWERFFTW})
set(${UPPERFFTW}_LIBRARIES ${CMAKE_BINARY_PREFIX}/${INSTALL_DIR}/lib/lib${FFTW_LIBNAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
set(${UPPERFFTW}_INCLUDE_DIRS ${CMAKE_BINARY_PREFIX}/${INSTALL_DIR}/include PARENT_SCOPE)
set(${UPPERFFTW}_FOUND TRUE PARENT_SCOPE)
set(${UPPERFFTW}_HAVE_SIMD TRUE PARENT_SCOPE)

add_library(gmxfftw STATIC IMPORTED GLOBAL)
set_target_properties(gmxfftw PROPERTIES IMPORTED_LOCATION ${${UPPERFFTW}_LIBRARIES})
set(${UPPERFFTW}_LIBRARIES gmxfftw PARENT_SCOPE)
add_dependencies(gmxfftw ${EXTERNAL_FFTW_BUILD_TARGET})

message(STATUS "The GROMACS-managed build of FFTW 3 will configure with the following optimizations: ${GMX_BUILD_OWN_FFTW_OPTIMIZATION_CONFIGURATION}")
