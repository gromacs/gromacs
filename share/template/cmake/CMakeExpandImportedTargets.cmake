#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014, by the GROMACS development team, led by
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

# CMAKE_EXPAND_IMPORTED_TARGETS(<var> LIBRARIES lib1 lib2...libN
#                                     [CONFIGURATION <config>] )
#
# CMAKE_EXPAND_IMPORTED_TARGETS() takes a list of libraries and replaces
# all imported targets contained in this list with their actual file paths
# of the referenced libraries on disk, including the libraries from their
# link interfaces.
# If a CONFIGURATION is given, it uses the respective configuration of the
# imported targets if it exists. If no CONFIGURATION is given, it uses
# the first configuration from ${CMAKE_CONFIGURATION_TYPES} if set, otherwise
# ${CMAKE_BUILD_TYPE}.
# This macro is used by all Check*.cmake files which use
# TRY_COMPILE() or TRY_RUN() and support CMAKE_REQUIRED_LIBRARIES , so that
# these checks support imported targets in CMAKE_REQUIRED_LIBRARIES:
#    cmake_expand_imported_targets(expandedLibs LIBRARIES ${CMAKE_REQUIRED_LIBRARIES}
#                                               CONFIGURATION "${CMAKE_TRY_COMPILE_CONFIGURATION}" )


#=============================================================================
# Copyright 2012 Kitware, Inc.
# Copyright 2009-2012 Alexander Neundorf <neundorf@kde.org>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:

# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.

# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# ------------------------------------------------------------------------------

# The above copyright and license notice applies to distributions of
# CMake in source and binary form.  Some source files contain additional
# notices of original copyright by their contributors; see each source
# for details.  Third-party software packages supplied with CMake under
# compatible licenses provide their own copyright notices documented in
# corresponding subdirectories.

# ------------------------------------------------------------------------------

# CMake was initially developed by Kitware with the following sponsorship:

#  * National Library of Medicine at the National Institutes of Health
#    as part of the Insight Segmentation and Registration Toolkit (ITK).

#  * US National Labs (Los Alamos, Livermore, Sandia) ASC Parallel
#    Visualization Initiative.

#  * National Alliance for Medical Image Computing (NAMIC) is funded by the
#    National Institutes of Health through the NIH Roadmap for Medical Research,
#    Grant U54 EB005149.

#  * Kitware, Inc.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# Imported by Mark Abraham into GROMACS, without modification, from
# cmake 2.8.8. See comments in CheckCXXLibraryExists.cmake for when
# this file might become superfluous.

include(CMakeParseArguments)

function(CMAKE_EXPAND_IMPORTED_TARGETS _RESULT )

   set(options )
   set(oneValueArgs CONFIGURATION )
   set(multiValueArgs LIBRARIES )

   cmake_parse_arguments(CEIT "${options}" "${oneValueArgs}" "${multiValueArgs}"  ${ARGN})

   if(CEIT_UNPARSED_ARGUMENTS)
      message(FATAL_ERROR "Unknown keywords given to CMAKE_EXPAND_IMPORTED_TARGETS(): \"${CEIT_UNPARSED_ARGUMENTS}\"")
   endif()

   if(NOT CEIT_CONFIGURATION)
      if(CMAKE_CONFIGURATION_TYPES)
         list(GET CMAKE_CONFIGURATION_TYPES 0 CEIT_CONFIGURATION)
      else()
         set(CEIT_CONFIGURATION ${CMAKE_BUILD_TYPE})
      endif()
   endif()

   # handle imported library targets

   set(_CCSR_REQ_LIBS ${CEIT_LIBRARIES})

   set(_CHECK_FOR_IMPORTED_TARGETS TRUE)
   set(_CCSR_LOOP_COUNTER 0)
   while(_CHECK_FOR_IMPORTED_TARGETS)
      math(EXPR _CCSR_LOOP_COUNTER "${_CCSR_LOOP_COUNTER} + 1 ")
      set(_CCSR_NEW_REQ_LIBS )
      set(_CHECK_FOR_IMPORTED_TARGETS FALSE)
      foreach(_CURRENT_LIB ${_CCSR_REQ_LIBS})
         get_target_property(_importedConfigs "${_CURRENT_LIB}" IMPORTED_CONFIGURATIONS)
         if (_importedConfigs)
#            message(STATUS "Detected imported target ${_CURRENT_LIB}")
            # Ok, so this is an imported target.
            # First we get the imported configurations.
            # Then we get the location of the actual library on disk of the first configuration.
            # then we'll get its link interface libraries property,
            # iterate through it and replace all imported targets we find there
            # with there actual location.

            # guard against infinite loop: abort after 100 iterations ( 100 is arbitrary chosen)
            if ("${_CCSR_LOOP_COUNTER}" LESS 100)
               set(_CHECK_FOR_IMPORTED_TARGETS TRUE)
#                else ("${_CCSR_LOOP_COUNTER}" LESS 1)
#                   message(STATUS "********* aborting loop, counter : ${_CCSR_LOOP_COUNTER}")
            endif ("${_CCSR_LOOP_COUNTER}" LESS 100)

            # if one of the imported configurations equals ${CMAKE_TRY_COMPILE_CONFIGURATION},
            # use it, otherwise simply use the first one:
            list(FIND _importedConfigs "${CEIT_CONFIGURATION}" _configIndexToUse)
            if("${_configIndexToUse}" EQUAL -1)
              set(_configIndexToUse 0)
            endif("${_configIndexToUse}" EQUAL -1)
            list(GET _importedConfigs ${_configIndexToUse} _importedConfigToUse)

            get_target_property(_importedLocation "${_CURRENT_LIB}" IMPORTED_LOCATION_${_importedConfigToUse})
            get_target_property(_linkInterfaceLibs "${_CURRENT_LIB}" IMPORTED_LINK_INTERFACE_LIBRARIES_${_importedConfigToUse} )

            list(APPEND _CCSR_NEW_REQ_LIBS  "${_importedLocation}")
#            message(STATUS "Appending lib ${_CURRENT_LIB} as ${_importedLocation}")
            if(_linkInterfaceLibs)
               foreach(_currentLinkInterfaceLib ${_linkInterfaceLibs})
#                  message(STATUS "Appending link interface lib ${_currentLinkInterfaceLib}")
                  if(_currentLinkInterfaceLib)
                     list(APPEND _CCSR_NEW_REQ_LIBS "${_currentLinkInterfaceLib}" )
                  endif(_currentLinkInterfaceLib)
               endforeach(_currentLinkInterfaceLib "${_linkInterfaceLibs}")
            endif(_linkInterfaceLibs)
         else(_importedConfigs)
            # "Normal" libraries are just used as they are.
            list(APPEND _CCSR_NEW_REQ_LIBS "${_CURRENT_LIB}" )
#            message(STATUS "Appending lib directly: ${_CURRENT_LIB}")
         endif(_importedConfigs)
      endforeach(_CURRENT_LIB ${_CCSR_REQ_LIBS})

      set(_CCSR_REQ_LIBS ${_CCSR_NEW_REQ_LIBS} )
   endwhile(_CHECK_FOR_IMPORTED_TARGETS)

   # Finally we iterate once more over all libraries. This loop only removes
   # all remaining imported target names (there shouldn't be any left anyway).
   set(_CCSR_NEW_REQ_LIBS )
   foreach(_CURRENT_LIB ${_CCSR_REQ_LIBS})
      get_target_property(_importedConfigs "${_CURRENT_LIB}" IMPORTED_CONFIGURATIONS)
      if (NOT _importedConfigs)
         list(APPEND _CCSR_NEW_REQ_LIBS "${_CURRENT_LIB}" )
#         message(STATUS "final: appending ${_CURRENT_LIB}")
      else (NOT _importedConfigs)
#             message(STATUS "final: skipping ${_CURRENT_LIB}")
      endif (NOT _importedConfigs)
   endforeach(_CURRENT_LIB ${_CCSR_REQ_LIBS})
#   message(STATUS "setting -${_RESULT}- to -${_CCSR_NEW_REQ_LIBS}-")
   set(${_RESULT} "${_CCSR_NEW_REQ_LIBS}" PARENT_SCOPE)

endfunction()
