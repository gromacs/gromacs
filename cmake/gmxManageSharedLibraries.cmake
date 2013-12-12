#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
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

# Manage the Gromacs shared librar setup. 

########################################################################
# Shared/static library settings
########################################################################
# Determine the defaults (this block has no effect if the variables have
# already been set)
if(APPLE OR CYGWIN OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux|.*BSD")
    # Maybe Solaris should be here? Patch this if you know!
    SET(SHARED_LIBS_DEFAULT ON)
elseif(WIN32 OR ${CMAKE_SYSTEM_NAME} MATCHES "BlueGene")
    # Support for shared libs on native Windows is a bit new. Its
    # default might change later if/when we sort things out. Also,
    # Cray should go here. What variable value can detect it?
    SET(SHARED_LIBS_DEFAULT OFF)
else()
    if (NOT DEFINED BUILD_SHARED_LIBS)
        message(STATUS "Defaulting to building static libraries")
    endif()
    SET(SHARED_LIBS_DEFAULT OFF)
endif()
if (GMX_PREFER_STATIC_LIBS)
    if (NOT DEFINED BUILD_SHARED_LIBS AND SHARED_LIBS_DEFAULT)
        message("Searching for static libraries requested, so the GROMACS "
                "libraries will also be static (BUILD_SHARED_LIBS=OFF)")
    endif()
    set(SHARED_LIBS_DEFAULT OFF)
endif()
set(GMX_PREFER_STATIC_LIBS_DEFAULT OFF)
if (WIN32 AND NOT CYGWIN AND NOT BUILD_SHARED_LIBS)
    set(GMX_PREFER_STATIC_LIBS_DEFAULT ON)
endif()

# Declare the user-visible options
option(BUILD_SHARED_LIBS "Enable shared libraries (can be problematic e.g. with MPI, or on some HPC systems)" ${SHARED_LIBS_DEFAULT})
if (UNIX)
    set(GMX_PREFER_STATIC_LIBS_DESCRIPTION
        "When finding libraries prefer static archives (it will only work if static versions of external dependencies are available and found)")
elseif (WIN32 AND NOT CYGWIN)
    set(GMX_PREFER_STATIC_LIBS_DESCRIPTION
        "When finding libraries prefer static system libraries (MT instead of MD)")
endif()
option(GMX_PREFER_STATIC_LIBS "${GMX_PREFER_STATIC_LIBS_DESCRIPTION}"
       ${GMX_PREFER_STATIC_LIBS_DEFAULT})
mark_as_advanced(GMX_PREFER_STATIC_LIBS)

# Act on the set values
if (UNIX AND GMX_PREFER_STATIC_LIBS)
    if(BUILD_SHARED_LIBS)
        # Warn the user about the combination. But don't overwrite the request.
        message(WARNING "Searching for static libraries requested, and "
                        "building shared Gromacs libraries requested. This "
                        "might cause problems linking later.")
    endif()
    # On Linux .a is the static library suffix, on Mac OS X .lib can also
    # be used, so we'll add both to the preference list.
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib;.a" ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif()
IF( WIN32 AND NOT CYGWIN)
  if (NOT BUILD_SHARED_LIBS)
      if(NOT GMX_PREFER_STATIC_LIBS)
          message(WARNING "Shared system libraries requested, and static Gromacs libraries requested.")
      endif()
  else()
      message(FATAL_ERROR "BUILD_SHARED_LIBS=ON not yet working for Windows in the master branch")
      if(GMX_PREFER_STATIC_LIBS)
          #this combination segfaults (illegal passing of file handles)
          message(FATAL_ERROR "Static system libraries requested, and shared Gromacs libraries requested.")
      endif()
      add_definitions(-DUSE_VISIBILITY -DTMPI_USE_VISIBILITY)
      set(PKG_CFLAGS "$PKG_CFLAGS -DUSE_VISIBILITY -DTMPI_USE_VISIBILITY")
  endif()

  IF (GMX_PREFER_STATIC_LIBS)
      #Only setting Debug and Release flags. Others configurations are current not used.
      STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
      STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
      STRING(REPLACE /MD /MT CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
      STRING(REPLACE /MD /MT CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
  ENDIF()
  IF( CMAKE_C_COMPILER_ID MATCHES "Intel" )
    if(BUILD_SHARED_LIBS) #not sure why incremental building with shared libs doesn't work
        STRING(REPLACE "/INCREMENTAL:YES" "" CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS})
    endif()
  ENDIF()
ENDIF()
