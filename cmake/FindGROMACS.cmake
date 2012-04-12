# - Finds parts of gromacs
# Find the native gromacs compents headers and libraries.
#
#  GROMACS_INCLUDE_DIRS - where to find gromacs headers.
#  GROMACS_LIBRARIES    - List of libraries when used by gromacs.
#  GROMACS_FOUND        - True if all gromacs componets were found.
#  GROMACS_DEFINITIONS  - Extra definies needed by gromacs
#  GROMACS_PKG          - The name of the pkg-config package needed
#  GROMACS_VERSION      - Gromacs lib interface version
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

find_package(PkgConfig)
list(LENGTH GROMACS_FIND_COMPONENTS GROMACS_NUM_COMPONENTS_WANTED)
if(${GROMACS_NUM_COMPONENTS_WANTED} LESS 1)
  message(FATAL_ERROR "No gromacs component to search given")
elseif(${GROMACS_NUM_COMPONENTS_WANTED} GREATER 1)
  message(FATAL_ERROR "We only support finding one gromacs component at this point, go and implement it ;-)")
elseif(${GROMACS_FIND_COMPONENTS} MATCHES "^lib(gmx|gromacs)(_d)?$")
  set(GROMACS_PKG "${GROMACS_FIND_COMPONENTS}")
  string(REGEX REPLACE "^lib(.*)" "\\1" GROMACS_LIBRARY_NAME "${GROMACS_PKG}")
else()
  message(FATAL_ERROR "We do not support finding ${GROMACS_FIND_COMPONENTS}, go and implement it ;-)")
endif()

if(GMX_DOUBLE AND NOT "${GROMACS_PKG}" MATCHES "_d$")
  message(FATAL_ERROR "GMX_DOUBLE was true, but I was asked to find ${GROMACS_PKG} (without _d at the end) - unlogical!")
endif(GMX_DOUBLE AND NOT "${GROMACS_PKG}" MATCHES "_d$")
if(NOT GMX_DOUBLE AND "${GROMACS_PKG}" MATCHES "_d$")
  message(FATAL_ERROR "GMX_DOUBLE was false, but I was asked to find ${GROMACS_PKG} (with _d at the end) - unlogical!")
endif(NOT GMX_DOUBLE AND "${GROMACS_PKG}" MATCHES "_d$")

pkg_check_modules(PC_GROMACS ${GROMACS_PKG})
if (GMX_DOUBLE)
  list(APPEND GMX_DEFS "-DGMX_DOUBLE")
endif(GMX_DOUBLE)
if (PC_GROMACS_CFLAGS_OTHER)
  foreach(DEF ${PC_GROMACS_CFLAGS_OTHER})
    if (${DEF} MATCHES "^-D")
      list(APPEND GMX_DEFS ${DEF})
    endif (${DEF} MATCHES "^-D")
  endforeach(DEF)
  list(REMOVE_DUPLICATES GMX_DEFS)
endif (PC_GROMACS_CFLAGS_OTHER)
set(GROMACS_DEFINITIONS "${GMX_DEFS}" CACHE STRING "extra GROMACS definitions")

find_library(GROMACS_LIBRARY NAMES ${GROMACS_LIBRARY_NAME} HINTS ${PC_GROMACS_LIBRARY_DIRS} )
if (GROMACS_LIBRARY)
  if("${GROMACS_LIBRARY}" MATCHES "lib(gmx|gromacs)[^;]*\\.a")
    if(PC_GROMACS_LIBRARIES)
      list(REMOVE_ITEM PC_GROMACS_LIBRARIES ${GROMACS_LIBRARY_NAME})
      foreach (LIB ${PC_GROMACS_LIBRARIES})
        find_library(GROMACS_${LIB} NAMES ${LIB} HINTS ${PC_GROMACS_LIBRARY_DIRS} )
        list(APPEND GMX_DEP_LIBRARIES ${GROMACS_${LIB}})
        unset(GROMACS_${LIB} CACHE)
      endforeach(LIB)
    endif(PC_GROMACS_LIBRARIES)
    if(PC_GROMACS_CFLAGS_OTHER)
      foreach(LIB ${PC_GROMACS_CFLAGS_OTHER})
        if (${LIB} MATCHES "thread")
  	find_package(Threads REQUIRED)
  	list(APPEND GMX_DEP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
        endif (${LIB} MATCHES "thread")
      endforeach(LIB)
    endif(PC_GROMACS_CFLAGS_OTHER)
    set(GROMACS_DEP_LIBRARIES "${GMX_DEP_LIBRARIES}" CACHE FILEPATH "gromacs depency libs (only needed for static (.a) libgmx")
  endif("${GROMACS_LIBRARY}" MATCHES "lib(gmx|gromacs)[^;]*\\.a")
  include(CheckLibraryExists)
  check_library_exists("${GROMACS_LIBRARY};${GROMACS_DEP_LIBRARIES}" GromacsVersion "" FOUND_GROMACS_VERSION)
  if(NOT FOUND_GROMACS_VERSION)
    message(FATAL_ERROR "Could not find GromacsVersion in ${GROMACS_LIBRARY};${GROMACS_DEP_LIBRARIES}, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you don't have pkg-config installed you will most likely have to set GROMACS_LIBRARY and GROMACS_DEP_LIBRARY by hand which sets the gromacs lib and it's depencies (i.e. -DGROMACS_LIBRARY='/path/to/libgmx.so' -DGROMACS_DEP_LIBRARIES='/path/to/libblas.so;/path/to/libm.so') !")
  endif(NOT FOUND_GROMACS_VERSION)
  check_library_exists("${GROMACS_LIBRARY};${GROMACS_DEP_LIBRARIES}" init_mtop "" FOUND_GROMACS_INIT_MTOP)
  if(NOT FOUND_GROMACS_INIT_MTOP)
    message(FATAL_ERROR "Could not find init_mtop in the gromacs library, take look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. This most likely means that your gromacs version is too old, we need at least gromacs 4.0.7.") 
  endif(NOT FOUND_GROMACS_INIT_MTOP)
  set(GMX_VERSION 40)
  check_library_exists("${GROMACS_LIBRARY};${GROMACS_DEP_LIBRARIES}" output_env_init "" FOUND_GROMACS_OUTPUT_ENV_INIT)
  if(FOUND_GROMACS_OUTPUT_ENV_INIT)
    set(GMX_VERSION 45)
  endif(FOUND_GROMACS_OUTPUT_ENV_INIT)
  check_library_exists("${GROMACS_LIBRARY};${GROMACS_DEP_LIBRARIES}" init_domdec_vsites "" FOUND_GROMACS_INIT_DOMDEC_VSITES)
  if(FOUND_GROMACS_INIT_DOMDEC_VSITES)
    set(GMX_VERSION 50)
  endif(FOUND_GROMACS_INIT_DOMDEC_VSITES)
  set(GROMACS_VERSION ${GMX_VERSION} CACHE STRING "Gromacs lib interface version")
else(GROMACS_LIBRARY)
  set(GMX_VERSION 45)
endif (GROMACS_LIBRARY)

if ("${GROMACS_PKG}" MATCHES "libgmx")
  if (${GMX_VERSION} EQUAL 40)
    find_path(GROMACS_INCLUDE_DIR tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
  else(${GMX_VERSION} EQUAL 40)
   find_path(GROMACS_INCLUDE_DIR gromacs/tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
  endif(${GMX_VERSION} EQUAL 40)
elseif("${GROMACS_PKG}" MATCHES "libgromacs")
  find_path(GROMACS_INCLUDE_DIR gromacs/legacyheaders/tpxio.h HINTS ${PC_GROMACS_INCLUDE_DIRS})
endif("${GROMACS_PKG}" MATCHES "libgmx")

set(GROMACS_LIBRARIES "${GROMACS_LIBRARY};${GROMACS_DEP_LIBRARIES}" )
set(GROMACS_INCLUDE_DIRS ${GROMACS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS DEFAULT_MSG GROMACS_LIBRARY GROMACS_INCLUDE_DIR)

mark_as_advanced(GROMACS_INCLUDE_DIR GROMACS_LIBRARY GROMACS_DEFINITIONS GROMACS_PKG GROMACS_VERSION GROMACS_DEP_LIBRARIES)
