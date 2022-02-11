#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2016- The GROMACS Authors
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

# Executing this macro will add the option GMX_USE_PLUGINS, which
# makes it possible to dynamically load modules at runtime.
# This is a very neat feature, and should virtually always work
# on Linux, but for now it will not work without shared libraries.
# For this reason we disable it by default, to avoid triggering
# errors here when dynamic libraries are disabled elsewhere.
macro(gmx_manage_plugin_support)
    option(GMX_USE_PLUGINS "Enable support for dynamic plugins (e.g. VMD-supported file formats)" OFF)
    mark_as_advanced(GMX_USE_PLUGINS)

    if(GMX_USE_PLUGINS)

        message(STATUS "Checking build environment for dynamic plugins")

        if(NOT BUILD_SHARED_LIBS)
            message(FATAL_ERROR "Shared libraries not built - required for dynamic plugins")
        endif()

        # Plugins are supported natively on Windows, so nothing to check if WIN32 is set

        if (NOT WIN32)
            include(gmxTestdlopen)
            gmx_test_dlopen(HAVE_DLOPEN)
            if(NOT HAVE_DLOPEN)
                message(FATAL_ERROR "dlopen() support missing - required for dynamic plugins")
            endif()
        endif()

        message(STATUS "Checking build environment for dynamic plugins - supported")

        list(APPEND GMX_EXTRA_LIBRARIES ${CMAKE_DL_LIBS}) # magic cross-platform pre-set variable for dlopen library
        set(PKG_DL_LIBS "-l${CMAKE_DL_LIBS}")

    endif()
endmacro()
