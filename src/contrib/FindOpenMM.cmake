#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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
#
# Find OpenMM library.
#
# Looks for the OpenMM libraries at the default (/usr/local) location 
# or custom location found in the OPENMM_ROOT_DIR environment variable. 
#
# The script defines defines: 
#  OpenMM_FOUND     
#  OpenMM_ROOT_DIR
#  OpenMM_INCLUDE_DIR
#  OpenMM_LIBRARY_DIR
#  OpenMM_LIBRARIES      
#  OpenMM_LIBRARIES_D   - debug version of libraries 
#  OpenMM_PLUGIN_DIR
#

# Author: Szilard Pall (pszilard@cbr.su.se)

if(OpenMM_INCLUDE_DIR AND OpenMM_LIBRARY_DIR AND OpenMM_PLUGIN_DIR)
    set(OpenMM_FIND_QUIETLY)
endif()

file(TO_CMAKE_PATH "$ENV{OPENMM_ROOT_DIR}" _env_OPENMM_ROOT_DIR)

if(IS_DIRECTORY ${_env_OPENMM_ROOT_DIR})
    set(OpenMM_ROOT_DIR "${_env_OPENMM_ROOT_DIR}" CACHE PATH "OpenMM installation directory" FORCE)
else()
    if(IS_DIRECTORY "/usr/local/openmm")
        set(OpenMM_ROOT_DIR "/usr/local/openmm" CACHE PATH "OpenMM installation directory" FORCE)
    endif()
endif()

find_library(OpenMM_LIBRARIES
    NAMES OpenMM
    PATHS "${OpenMM_ROOT_DIR}/lib"
    CACHE STRING "OpenMM libraries")

find_library(OpenMM_LIBRARIES_D
    NAMES OpenMM_d
    PATHS "${OpenMM_ROOT_DIR}/lib"
    CACHE STRING "OpenMM debug libraries")

if(OpenMM_LIBRARIES_D AND NOT OpenMM_LIBRARIES)
    set(OpenMM_LIBRARIES ${OpenMM_LIBRARIES_D}
        CACHE STRING "OpenMM libraries" FORCE)
    message(WARNING " Only found debug versions of the OpenMM libraries!")
endif()

get_filename_component(OpenMM_LIBRARY_DIR 
    ${OpenMM_LIBRARIES} 
    PATH
    CACHE STRING "OpenMM library path")

find_path(OpenMM_INCLUDE_DIR 
    NAMES OpenMM.h 
    PATHS "${OpenMM_ROOT_DIR}/include" "${OpenMM_LIBRARY_DIR}/../include"
    CACHE STRING "OpenMM include directory")    

# if we did not manage to set the root dir at the beginning but found the 
# libs then set the ${OpenMM_LIBRARY_DIR}/.. as root
if(NOT IS_DIRECTORY ${OpenMM_ROOT_DIR})
    if (IS_DIRECTORY "${OpenMM_LIBRARY_DIR}/..") # just double-checking
        get_filename_component(OpenMM_ROOT_DIR 
            "${OpenMM_LIBRARY_DIR}/.." 
            ABSOLUTE)
    endif()   
endif()

if(NOT IS_DIRECTORY ${OpenMM_ROOT_DIR})
    message(FATAL_ERROR "Could not find OpenMM! Set the OPENMM_ROOT_DIR environment "
    "variable to contain the path of the OpenMM installation.")
endif()

if(NOT IS_DIRECTORY ${OpenMM_LIBRARY_DIR})
    message(FATAL_ERROR "Can't find OpenMM libraries. Check your OpenMM installation!")
endif()

# now we can be sure that we have the library dir
if(IS_DIRECTORY "${OpenMM_LIBRARY_DIR}/plugins")
    get_filename_component(OpenMM_PLUGIN_DIR
        "${OpenMM_LIBRARY_DIR}/plugins"
        ABSOLUTE)
    set(OpenMM_PLUGIN_DIR ${OpenMM_PLUGIN_DIR} CACHE PATH "OpenMM plugins directory")
else()
    message(WARNING "Could not detect the OpenMM plugin directory at the default location (${OpenMM_LIBRARY_DIR}/plugins)."
            "Check your OpenMM installation or set the OPENMM_PLUGIN_DIR environment variable!")
endif()

if(NOT OpenMM_INCLUDE_DIR)
    message(FATAL_ERROR "Can't find OpenMM includes. Check your OpenMM installation!")
endif()

set(OpenMM_ROOT_DIR ${OpenMM_ROOT_DIR} CACHE PATH "OpenMM installation directory")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMM DEFAULT_MSG 
                                    OpenMM_ROOT_DIR
                                    OpenMM_LIBRARIES 
                                    OpenMM_LIBRARY_DIR 
                                    OpenMM_INCLUDE_DIR)

mark_as_advanced(OpenMM_INCLUDE_DIR
    OpenMM_LIBRARIES
    OpenMM_LIBRARIES_D
    OpenMM_LIBRARY_DIR)
