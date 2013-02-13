#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009-2011, by the VOTCA Development Team (http://www.votca.org)
# Copyright (c) 2012,2013 by the GROMACS development team, led by
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
# - Find gsl
# Find the native GSL headers and libraries.
#
#  GSL_INCLUDE_DIRS - where to find gsl/gsl_linalg.h, etc.
#  GSL_LIBRARIES    - List of libraries when using gsl.
#  GSL_FOUND        - True if gsl found.
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
