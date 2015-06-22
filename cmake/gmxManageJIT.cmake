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

# Requires gmxVersionInfo.cmake to have been run already so we have
# GMX_VERSION_STRING_FULL

function(gmx_manage_jit)
    # Make directory structure for GROMACS to store JIT artifacts
    if(APPLE OR UNIX)
        set(GMX_JIT_CACHE_BASE_DIR "$ENV{HOME}/.gromacs")
    elseif(WINDOWS)
        set(GMX_JIT_CACHE_BASE_DIR "$ENV{USERPROFILE}/.gromacs")
    elseif()
        message(FATAL_ERROR "Cannot manage creating a JIT cache on this unknown OS")
    endif()

    # Set up a cache variable so that mdrun knows where to read and
    # write the JIT artifacts.
    #
    # TODO GMX_VERSION_STRING is good enough for managing caching of
    # installed versions of GROMACS, but we should use something like
    # GMX_VERSION_STRING_FULL when we expand this machinery for
    # development-time use. Right now, it is unclear how best to get
    # access to that.
    set(GMX_JIT_CACHE_DIR "${GMX_JIT_CACHE_BASE_DIR}/${GMX_VERSION_STRING}" CACHE INTERNAL "Version-specific location for mdrun to read and write JIT-compiled artifacts")
    mark_as_advanced(GMX_JIT_CACHE_DIR)

    # Configure a CMake script that we use at build time to make
    # the cache directory
    #
    # TODO This requires that the cache directory used at run time is
    # accessible at install time, which might not be possible. It
    # might be more robust in the future for mdrun to be able to make
    # its own cache folder.
    configure_file(cmake/MakeGromacsJitCacheFolder.cmake.cmakein MakeGromacsJitCacheFolder.cmake)
    install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/MakeGromacsJitCacheFolder.cmake)
endfunction()

if(GMX_USE_OPENCL)
    gmx_manage_jit()
endif()