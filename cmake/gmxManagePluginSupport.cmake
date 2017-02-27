#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2016,2017, by the GROMACS development team, led by
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

include(gmxOptionUtilities)

# Sets GMX_USE_PLUGINS=ON in the parent scope if the toolchain and
# user selections permit the build to support plugin loading.
function(gmx_manage_plugin_support)
    gmx_option_trivalue(GMX_LOAD_PLUGINS "Compile with plugin support, needed to read VMD supported file formats" AUTO)
    mark_as_advanced(GMX_LOAD_PLUGINS)

    # Find out if non-Windows builds can support plugins. Native Windows
    # neither needs nor has library support.
    if (NOT WIN32)
        # TODO Make a proper find_package for dlopen to find
        # dlfcn.h. The CMake variable CMAKE_DL_LIBS works magically
        # for the library, however.
        include(gmxTestdlopen)
        gmx_test_dlopen(HAVE_DLOPEN)
    endif()

    # Keep the status line quiet unless something relevant has
    # changed.
    gmx_check_if_changed(EMIT_STATUS_MESSAGES GMX_LOAD_PLUGINS BUILD_SHARED_LIBS HAVE_DLOPEN)

    # Whether GROMACS will really try to compile support for VMD
    # plugins.
    set(GMX_USE_PLUGINS OFF)

    # Plugins are supported natively on Windows
    if (WIN32 OR (BUILD_SHARED_LIBS AND HAVE_DLOPEN))
        set(GMX_USE_PLUGINS ${GMX_LOAD_PLUGINS})
    elseif(GMX_LOAD_PLUGINS)
        # Can't support plugins for some reason. If the user required
        # plugins, emit fatal errors. Otherwise, emit status messages
        # for AUTO and be silent for OFF.
        set(message "")
        if (NOT HAVE_DLOPEN)
            set(message "${message}dlopen() support for using dynamic plugins for VMD-supported file formats is missing. ")
        endif()
        if(NOT BUILD_SHARED_LIBS)
            set(message "${message}GROMACS only supports plugins in a build that uses shared libraries, which can be disabled for various reasons. BUILD_SHARED_LIBS=on and a toolchain that supports dynamic linking is required. (Hint: GMX_PREFER_STATIC_LIBS and GMX_BUILD_MDRUN_ONLY can influence the default BUILD_SHARED_LIBS, so if you need plugins, reconsider those choices.) ")
        endif()
        if (GMX_LOAD_PLUGINS_FORCE)
            message(FATAL_ERROR "${message}Cannot build with GMX_LOAD_PLUGINS=${GMX_LOAD_PLUGINS}.")
        endif()
        if (GMX_LOAD_PLUGINS_AUTO AND EMIT_STATUS_MESSAGES)
            message(STATUS "${message}")
        endif()
    endif()

    if(EMIT_STATUS_MESSAGES)
        if(GMX_USE_PLUGINS)
            MESSAGE(STATUS "Using dynamic plugins (e.g VMD-supported file formats)")
        else()
            MESSAGE(STATUS "Not using dynamic plugins (e.g VMD-supported file formats)")
        endif()
    endif()
    set(GMX_USE_PLUGINS ${GMX_USE_PLUGINS} PARENT_SCOPE)
endfunction()

gmx_manage_plugin_support()

if(GMX_USE_PLUGINS)
    list(APPEND GMX_EXTRA_LIBRARIES ${CMAKE_DL_LIBS}) # magic cross-platform pre-set variable for dlopen library
    set(PKG_DL_LIBS "-l${CMAKE_DL_LIBS}")
endif()
