#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2016, by the GROMACS development team, led by
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

# This implementation of FindTBB.cmake was adapted from the one present in the
# OGRE repository and modifications by Robert Maynard. To stay in that spirit,
# all the contents of this file is also placed in the public domain and can be
# redistributed freely.
#

# - Find Intel TBB
#
# Usage: find_package(TBB [QUIET] [REQUIRED])
#
# In addition to CMAKE_PREFIX_PATH and default system locations, we also
# look in the usual locations for the commercial TBB version in /opt.
# It is also possible to set TBB_ROOT separately to indicate the location.
#
# For default usage, you probably want to add both TBB_LIBRARIES and
# TBB_MALLOC_LIBRARIES to your link lines, but they are separate variables
# in case you have an even better allocator. Read the TBB documentation for
# information about when/how to use the proxy library.
#
# FindTBB defines the following variables:
#
#  TBB_FOUND                          True if headers and libraries were located
#  TBB_INCLUDE_DIRS                   Where to find TBB headers
#  TBB_LIBRARIES                      The main TBB libraries
#  TBB_MALLOC_LIBRARIES               TBB malloc libraries
#  TBB_MALLOC_PROXY_LIBRARIES         TBB malloc proxy libraries
#  TBB_VERSION                        Overall tbb version (major.minor)
#  TBB_VERSION_MAJOR                  Major Product Version Number
#  TBB_VERSION_MINOR                  Minor Product Version Number
#  TBB_INTERFACE_VERSION              Engineering Focused Version Number
#  TBB_COMPATIBLE_INTERFACE_VERSION   The oldest major interface version
#                                     still supported. This uses the engineering
#                                     focused interface version numbers.
#
#  Note: You can inspect the library variables to check of the malloc or
#  malloc proxy library was found, but this should typically not be critical.
#  All that happens without them is that you use the default allocator,
#  and if you don't have the libraries that is what you will have to do anyway.
#

##############################################
# Start of TBB library search code
##############################################

#
# Locations to search (in addition to CMAKE_PREFIX_PATH)
#
set(TBB_PREFIX_PATH "/opt/intel/tbb" "/usr/local/tbb")

#
# Create the lists with hints about locations for headers and libraries
#
foreach(dir ${TBB_PREFIX_PATH})
    list(APPEND TBB_INCLUDE_DIRS_HINTS ${dir}/include ${dir}/Include ${dir}/include/tbb)
    list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/lib ${dir}/Lib ${dir}/lib/tbb ${dir}/Libs)
endforeach()

##############################################
# Do a bunch of special stuff for windows
##############################################
#
# get the arch, only used by windows to select intel64/ia32
#
if($ENV{TBB_ARCH_PLATFORM})
    set(TBB_ARCH_PLATFORM $ENV{TBB_ARCH_PLATFORM})
endif()

# For Windows, let's assume that the user might be using the precompiled
# TBB packages from the main website. These use a rather awkward directory
# structure (at least for automatically finding the right files) depending
# on platform and compiler, but we'll do our best to accommodate it.
# Not adding the same effort for the precompiled linux builds, though. Those
# have different versions for CC compiler versions and linux kernels which
# will never adequately match the user's setup, so there is no feasible way
# to detect the "best" version to use. The user will have to manually
# select the right files. (Chances are the distributions are shipping their
# custom version of tbb, anyway, so the problem is probably nonexistant.)
if (WIN32 AND MSVC)
    set(COMPILER_PREFIX "vc7.1")
    if (MSVC_VERSION EQUAL 1400)
        set(COMPILER_PREFIX "vc8")
    elseif(MSVC_VERSION EQUAL 1500)
        set(COMPILER_PREFIX "vc9")
    elseif(MSVC_VERSION EQUAL 1600)
        set(COMPILER_PREFIX "vc10")
    elseif(MSVC_VERSION EQUAL 1700)
        set(COMPILER_PREFIX "vc11")
    elseif(MSVC_VERSION EQUAL 1800)
        set(COMPILER_PREFIX "vc12")
    endif()

    # for each prefix path, add ia32/64\${COMPILER_PREFIX}\lib to the
    # list of hints for library locations
    foreach (dir ${TBB_PREFIX_PATH})
        if(CMAKE_CL_64)
            list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/ia64/${COMPILER_PREFIX}/lib)
            list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/lib/ia64/${COMPILER_PREFIX})
            list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/intel64/${COMPILER_PREFIX}/lib)
            list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/lib/intel64/${COMPILER_PREFIX})
        else()
            list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/ia32/${COMPILER_PREFIX}/lib)
            list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/lib/ia32/${COMPILER_PREFIX})
        endif()
    endforeach()
endif()

foreach (dir ${TBB_PREFIX_PATH})
  list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/${TBB_ARCH_PLATFORM}/lib)
  list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/lib/${TBB_ARCH_PLATFORM})
  list(APPEND TBB_LIBRARY_DIRS_HINTS ${dir}/lib)
endforeach ()

##############################
# Locate the main TBB library
##############################
#
# Find the include headers
#
find_path(TBB_INCLUDE_DIRS "tbb/tbb.h" HINTS ${TBB_INCLUDE_DIRS_HINTS})

#
# Find the release & debug versions of the base TBB library
#
find_library(TBB_LIBRARY_RELEASE NAMES "tbb"       HINTS ${TBB_LIBRARY_DIRS_HINTS})
find_library(TBB_LIBRARY_DEBUG   NAMES "tbb_debug" HINTS ${TBB_LIBRARY_DIRS_HINTS})

# The release version library is mandatory - don't set TBB_LIBRARIES if it was not found
if(TBB_LIBRARY_RELEASE AND TBB_LIBRARY_DEBUG)
    set(TBB_LIBRARIES optimized ${TBB_LIBRARY_RELEASE} debug ${TBB_LIBRARY_DEBUG} CACHE STRING "Intel TBB libraries")
elseif(TBB_LIBRARY_RELEASE)
    set(TBB_LIBRARIES ${TBB_LIBRARY_RELEASE} CACHE STRING "Intel TBB libraries")
endif()
# Unset the cached find_library results (since edits must go in TBB_LIBRARIES)
unset(TBB_LIBRARY_RELEASE CACHE)
unset(TBB_LIBRARY_DEBUG CACHE)


if(TBB_INCLUDE_DIRS AND NOT TBB_VERSION)
    #parse all the version numbers from tbb
    # Only read the start of the file. Here we rely on the directory
    # variable only containing a single directory, but that should be fine.
    if("${TBB_INCLUDE_DIRS}" MATCHES "Frameworks")
        file(READ "${TBB_INCLUDE_DIRS}/Headers/tbb_stddef.h"
             TBB_VERSION_CONTENTS LIMIT 2048)
    else()
        file(READ "${TBB_INCLUDE_DIRS}/tbb/tbb_stddef.h"
             TBB_VERSION_CONTENTS LIMIT 2048)
    endif()

    string(REGEX REPLACE ".*#define TBB_VERSION_MAJOR ([0-9]+).*" "\\1"
           VERSION_MAJOR "${TBB_VERSION_CONTENTS}")

    string(REGEX REPLACE ".*#define TBB_VERSION_MINOR ([0-9]+).*" "\\1"
           VERSION_MINOR "${TBB_VERSION_CONTENTS}")

    string(REGEX REPLACE ".*#define TBB_INTERFACE_VERSION ([0-9]+).*" "\\1"
           INTERFACE_VERSION "${TBB_VERSION_CONTENTS}")

    string(REGEX REPLACE ".*#define TBB_COMPATIBLE_INTERFACE_VERSION ([0-9]+).*" "\\1"
           COMPATIBLE_INTERFACE_VERSION "${TBB_VERSION_CONTENTS}")

    set(TBB_VERSION_MAJOR ${VERSION_MAJOR} CACHE INTERNAL "Intel TBB major version")
    set(TBB_VERSION_MINOR ${VERSION_MINOR} CACHE INTERNAL "Intel TBB minor version")
    set(TBB_INTERFACE_VERSION ${INTERFACE_VERSION}
        CACHE INTERNAL "Intel TBB interface version")
    set(TBB_COMPATIBLE_INTERFACE_VERSION ${COMPATIBLE_INTERFACE_VERSION}
        CACHE INTERNAL "Intel TBB compatible interface version")

    # We only expose the overall version to the user
    set(TBB_VERSION "${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR}" CACHE STRING "Intel TBB version")

endif()

find_package_handle_standard_args(TBB FOUND_VAR TBB_FOUND
                                  REQUIRED_VARS TBB_LIBRARIES TBB_INCLUDE_DIRS
                                  VERSION_VAR TBB_VERSION)

if(TBB_FOUND)
    #on unix we need to also link to rt
    if(UNIX AND NOT APPLE)
        list(APPEND TBB_LIBRARIES rt)
    endif()
endif()



################################
# Locate the TBB malloc library
################################
#
# Find the release & debug versions of the TBB malloc library
#
find_library(TBB_MALLOC_LIBRARY_RELEASE NAMES "tbbmalloc"       HINTS ${TBB_LIBRARY_DIRS_HINTS})
find_library(TBB_MALLOC_LIBRARY_DEBUG   NAMES "tbbmalloc_debug" HINTS ${TBB_LIBRARY_DIRS_HINTS})

if(TBB_MALLOC_LIBRARY_RELEASE AND TBB_MALLOC_LIBRARY_DEBUG)
    set(TBB_MALLOC_LIBRARIES
        optimized ${TBB_MALLOC_LIBRARY_RELEASE} debug ${TBB_MALLOC_LIBRARY_DEBUG}
        CACHE STRING "Intel TBB malloc libraries")
elseif(TBB_MALLOC_LIBRARY_RELEASE)
    set(TBB_MALLOC_LIBRARIES ${TBB_MALLOC_LIBRARY_RELEASE}
        CACHE STRING "Intel TBB malloc libraries")
endif()

# Unset the cached find_library results (since edits must go in TBB_MALLOC_LIBRARIES)
unset(TBB_MALLOC_LIBRARY_RELEASE CACHE)
unset(TBB_MALLOC_LIBRARY_DEBUG CACHE)


######################################
# Locate the TBB malloc proxy library
######################################
#
# Find the release & debug versions of the TBB malloc proxy library
#
find_library(TBB_MALLOC_PROXY_LIBRARY_RELEASE NAMES "tbbmalloc_proxy"       HINTS ${TBB_LIBRARY_DIRS_HINTS})
find_library(TBB_MALLOC_PROXY_LIBRARY_DEBUG   NAMES "tbbmalloc_proxy_debug" HINTS ${TBB_LIBRARY_DIRS_HINTS})

if(TBB_MALLOC_PROXY_LIBRARY_RELEASE AND TBB_MALLOC_PROXY_LIBRARY_DEBUG)
    set(TBB_MALLOC_PROXY_LIBRARIES
        optimized ${TBB_MALLOC_PROXY_LIBRARY_RELEASE}
        debug ${TBB_MALLOC_PROXY_LIBRARY_DEBUG}
        CACHE STRING "Intel TBB malloc proxy libraries")
elseif(TBB_MALLOC_PROXY_LIBRARY_RELEASE)
    set(TBB_MALLOC_PROXY_LIBRARIES ${TBB_MALLOC_PROXY_LIBRARY_RELEASE}
        CACHE STRING "Intel TBB malloc proxy libraries")
endif()

# Unset the cached find_library results (since edits must go in TBB_MALLOC_PROXY_LIBRARIES)
unset(TBB_MALLOC_PROXY_LIBRARY_RELEASE CACHE)
unset(TBB_MALLOC_PROXY_LIBRARY_DEBUG CACHE)


mark_as_advanced(TBB_INCLUDE_DIRS TBB_LIBRARIES TBB_MALLOC_LIBRARIES TBB_MALLOC_PROXY_LIBRARIES)







