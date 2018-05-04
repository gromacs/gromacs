#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

# Build a suitable copy of clFFT.

#FIXME should not be needed in master

set (OpenCL_INCLUDE_DIR ${OPENCL_INCLUDE_DIR})
set (OpenCL_LIBRARY ${OPENCL_LIBRARY})
message (STATUS "clFFT is building with " ${OPENCL_LIBRARY} " " ${OPENCL_INCLUDE_DIR})

#TODO disable boost/fftw detection

set(EXTERNAL_CLFFT_BUILD_TARGET clFFTBuild)
include(ExternalProject)
ExternalProject_add(${EXTERNAL_CLFFT_BUILD_TARGET}
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/external/clFFT/src
  INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR} #???
  
  CMAKE_ARGS
  -DBUILD_EXAMPLES=OFF -DBUILD_TEST=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_CALLBACK_CLIENT=OFF -DBUILD_CLIENT=OFF
  # TODO BUILD_LOADLIBRARIES?
  -DOpenCL_INCLUDE_DIR=${OpenCL_INCLUDE_DIR} -DOpenCL_LIBRARY=${OpenCL_LIBRARY}
  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC
  INSTALL_COMMAND ""
  #cp ./library/libclFFT.a ../../../../../../lib/
  #FIXME
  # ./src/external/build-clFFT/clFFTBuild-prefix/src/clFFTBuild-build/library/libclFFT.a
  )

ExternalProject_get_property(${EXTERNAL_CLFFT_BUILD_TARGET} INSTALL_DIR)
message("clFFT Install dir: ${INSTALL_DIR}")
