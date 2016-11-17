#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2016, by the GROMACS development team, led by
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

set(GMX_TNG_MINIMUM_REQUIRED_VERSION "1.7.10")

gmx_dependent_option(
    GMX_EXTERNAL_TNG
    "Use external TNG instead of compiling the version shipped with GROMACS."
    OFF
    GMX_USE_TNG)
gmx_dependent_option(
    GMX_EXTERNAL_ZLIB
    "Use external ZLIB instead of compiling the version shipped with GROMACS as part of TNG."
    OFF
    "NOT GMX_EXTERNAL_TNG")

if(GMX_USE_TNG)
    # Detect TNG if GMX_EXTERNAL_TNG is explicitly ON
    if(GMX_EXTERNAL_TNG)
        find_package(TNG_IO ${GMX_TNG_MINIMUM_REQUIRED_VERSION})
        if(NOT TNG_IO_FOUND)
            message(FATAL_ERROR "TNG >= ${GMX_TNG_MINIMUM_REQUIRED_VERSION} not found. You can set GMX_EXTERNAL_TNG=OFF to compile the TNG bundled with GROMACS.")
        endif()
    else()
        # Detect zlib if the user requires us to use an external
        # version. If found, it can be used by TNG.
        if(GMX_EXTERNAL_ZLIB)
            find_package(ZLIB)
            if(NOT ZLIB_FOUND)
                message(FATAL_ERROR "External zlib compression library was required but could not be found. Set GMX_EXTERNAL_ZLIB=OFF to compile zlib as part of GROMACS.")
            endif()
            include(gmxTestZLib)
            gmx_test_zlib(HAVE_ZLIB)
            if(NOT HAVE_ZLIB)
                message(FATAL_ERROR "External zlib compression library was required but could not compile and link. Set GMX_EXTERNAL_ZLIB=OFF to compile zlib as part of GROMACS.")
            endif()
        endif()
    endif()
endif()

function(gmx_setup_tng_for_libgromacs)
    set(BUNDLED_TNG_LOCATION "${CMAKE_SOURCE_DIR}/src/external/tng_io")
    if (GMX_USE_TNG)
        if (GMX_EXTERNAL_TNG)
            target_link_libraries(libgromacs PRIVATE tng_io::tng_io)
        else()
            set(_zlib_arg)
            if (NOT GMX_EXTERNAL_ZLIB)
                set(_zlib_arg OWN_ZLIB)
            endif()
            include(${BUNDLED_TNG_LOCATION}/BuildTNG.cmake)
            add_tng_io_library(tng_io OBJECT ${_zlib_arg})
            add_library(tng_io::tng_io ALIAS tng_io)
            target_link_libraries(libgromacs PRIVATE $<BUILD_INTERFACE:tng_io::tng_io>)
        endif()
    endif()
endfunction()
