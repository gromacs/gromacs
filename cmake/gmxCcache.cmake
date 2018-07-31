#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018, by the GROMACS development team, led by
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

# Reference https://crascit.com/2016/04/09/using-ccache-with-cmake/
#
# Here we try to make sure that ccache is invoked as `/.../ccache compiler args` to best handle more than one local
# compiler or a compiler wrapper, whereas it is otherwise common to replace the default compilers with symbolic links
# to the ccache binary.
#
# ccache only works for gcc compatible compilers. We should test with anything other than CMAKE_<LANG>_COMPILER_ID==GNU
# Clang is reported to work with some caveats. See https://pspdfkit.com/blog/2015/ccache-for-fun-and-profit/
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    # Check whether C compiler wrapper has been set up.
    if(NOT DEFINED GMX_CACHE_C_COMPILER)
        # Determine whether we have a cacheable compiler.
        set(_cacheable OFF)
        if (CMAKE_C_COMPILER_ID MATCHES "GNU"
            OR CMAKE_C_COMPILER_ID MATCHES "AppleClang"
            OR CMAKE_C_COMPILER_ID MATCHES "Clang")
            message(STATUS "Setting up ccache wrapper for ${CMAKE_C_COMPILER_ID} C compiler ${CMAKE_C_COMPILER}")
            set(_c_launcher "${CCACHE_PROGRAM}")
            configure_file(launch-c.in ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/launch-c)
            unset(_c_launcher)
            file(COPY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/launch-c
                 DESTINATION ${CMAKE_BINARY_DIR}
                 FILE_PERMISSIONS
                 OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                 )
            set(_cacheable ON)
        else()
            message(STATUS "Disabling ccache set up. Not confirmed to work with compiler ID ${CMAKE_C_COMPILER_ID}.")
        endif() # GNU C compiler
        set(GMX_CACHE_C_COMPILER ${_cacheable} CACHE INTERNAL "Whether the C compiler will be wrapped for caching.")
        unset(_cacheable)
    endif() # defined
    # Check whether we should use the wrapper. If so, set CMAKE variables.
    if(GMX_CACHE_C_COMPILER)
        if(CMAKE_GENERATOR STREQUAL "Xcode")
            # Set Xcode project attributes to route compilation and linking
            # through our scripts
            set(CMAKE_XCODE_ATTRIBUTE_CC "${CMAKE_BINARY_DIR}/launch-c")
            set(CMAKE_XCODE_ATTRIBUTE_LD "${CMAKE_BINARY_DIR}/launch-c")
        else()
            # Support Unix Makefiles and Ninja
            set(CMAKE_C_COMPILER_LAUNCHER "${CMAKE_BINARY_DIR}/launch-c")
        endif()
    endif(GMX_CACHE_C_COMPILER)

    # Check whether CXX compiler wrapper has been set up
    if(NOT DEFINED GMX_CACHE_CXX_COMPILER)
        # Determine whether we have a cacheable compiler.
        set(_cacheable OFF)
        if (CMAKE_CXX_COMPILER_ID MATCHES "GNU"
            OR CMAKE_CXX_COMPILER_ID MATCHES "AppleClang"
            OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            message(STATUS "Setting up ccache wrapper for ${CMAKE_CXX_COMPILER_ID} CXX compiler ${CMAKE_CXX_COMPILER}")
            set(_cxx_launcher "${CCACHE_PROGRAM}")
            configure_file(launch-cxx.in ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/launch-cxx)
            unset(_cxx_launcher)
            file(COPY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/launch-cxx
                 DESTINATION ${CMAKE_BINARY_DIR}
                 FILE_PERMISSIONS
                 OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                )
            set(_cacheable ON)
        else()
            message(STATUS "Skipping ccache set up. Not confirmed to work with compiler ID ${CMAKE_CXX_COMPILER_ID}.")
        endif() # GNU C++ compiler
        set(GMX_CACHE_CXX_COMPILER ${_cacheable} CACHE INTERNAL "Whether the C++ compiler will be wrapped for caching.")
        unset(_cacheable)
    endif() # defined
    # Check whether we should use the wrapper. If so, set CMAKE variables.
    if(GMX_CACHE_CXX_COMPILER)
        if(CMAKE_GENERATOR STREQUAL "Xcode")
            # Set Xcode project attributes to route compilation and linking
            # through our scripts
            set(CMAKE_XCODE_ATTRIBUTE_CXX "${CMAKE_BINARY_DIR}/launch-cxx")
            set(CMAKE_XCODE_ATTRIBUTE_LDPLUSPLUS "${CMAKE_BINARY_DIR}/launch-cxx")
        else()
            # Support Unix Makefiles and Ninja
            set(CMAKE_CXX_COMPILER_LAUNCHER "${CMAKE_BINARY_DIR}/launch-cxx")
        endif()
    endif(GMX_CACHE_CXX_COMPILER)
endif(CCACHE_PROGRAM) # ccache program
