# - Find gsl
# Find the native GSL headers and libraries.
#
#  GSL_INCLUDE_DIRS - where to find gsl/gsl_linalg.h, etc.
#  GSL_LIBRARIES    - List of libraries when using gsl.
#  GSL_FOUND        - True if gsl found.
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
# Copyright 2012 The Gromacs Team
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

pkg_check_modules(PC_GSL gsl)
find_path(GSL_INCLUDE_DIR gsl/gsl_linalg.h HINTS ${PC_GSL_INCLUDE_DIRS})
find_library(GSL_LIBRARY NAMES gsl HINTS ${PC_GSL_LIBRARY_DIRS} )
find_library(GSLCBLAS_LIBRARY NAMES gslcblas cblas HINTS ${PC_GSL_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GSL DEFAULT_MSG GSL_LIBRARY GSL_INCLUDE_DIR GSLCBLAS_LIBRARY)

set(GSL_LIBRARIES "${GSL_LIBRARY};${GSLCBLAS_LIBRARY}")
set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR} )

if (GSL_FOUND)
  include(CheckLibraryExists)
  check_library_exists("${GSLCBLAS_LIBRARY};${MATH_LIBRARIES}" cblas_dsyrk "" FOUND_DSYRK)
  if(NOT FOUND_DSYRK)
    message(FATAL_ERROR "Could not find cblas_dsyrk in ${GSLCBLAS_LIBRARY}, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you are using a static lib (.a) make sure you have specified all dependencies of libcblas in GSLCBLAS_LIBRARY by hand (e.g. -DGSLCBLAS_LIBRARY='/path/to/libcblas.so;/path/to/libm.so')! If your gsl was build against an different version of cblas, specify it in GSLCBLAS_LIBRARY")
  endif(NOT FOUND_DSYRK)
  #adding MATH_LIBRARIES here to allow static libs, this does not harm us as we are anyway using it
  check_library_exists("${GSL_LIBRARIES};${MATH_LIBRARIES}" gsl_linalg_QR_decomp "" FOUND_QR_DECOMP)
  if(NOT FOUND_QR_DECOMP)
    message(FATAL_ERROR "Could not find gsl_linalg_QR_decompx in ${GSL_LIBRARY}, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong.  If you are using a static lib (.a) make sure you have specified all dependencies of libgsl in GSL_LIBRARY by hand (e.g. -DGSL_LIBRARY='/path/to/libgsl.so;/path/to/libm.so') !")
  endif(NOT FOUND_QR_DECOMP)
endif (GSL_FOUND)

mark_as_advanced(GSL_INCLUDE_DIR GSL_LIBRARY)
