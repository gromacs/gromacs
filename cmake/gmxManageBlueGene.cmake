#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

# Managing configuration for all kinds of BlueGene systems
# BlueGene/L is probably obsolete, but does no harm
# BlueGene/P needs testing, but hasn't changed
# BlueGene/Q works
message(STATUS "Configuring for BlueGene")

if (${CMAKE_SYSTEM_NAME} STREQUAL "BlueGeneL")
    # BlueGene/L never had shared lib support.
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Shared libraries not compatible with BlueGene/L, disabled!" FORCE)
endif()
if (${CMAKE_SYSTEM_NAME} MATCHES "BlueGene.*static")
    # BlueGene/P claims shared library support, but Mark Abraham never
    # got it to work. BlueGene/Q claims it, but discourages it for
    # performance reasons. So unless information to the contrary ever
    # comes to light, we should not mess about giving the user options
    # that are useless when they've already selected a static toolchain.
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Static BlueGene build toolchain selected, so shared libraries are disabled" FORCE)
endif()

set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Do not use software reciprocal square root on BlueGene")
set(GMX_X11 OFF CACHE BOOL "X11 not compatible with BlueGene, disabled!" FORCE)
set(GMX_GPU OFF CACHE BOOL "Cannot do GPU acceleration on BlueGene" FORCE)

# It is conceivable you could use ThreadMPI on BlueGene/Q by using its
# facility to run lots of jobs on small chunks of the machine. You
# certainly need proper MPI to use a whole chunk of the machine that
# the scheduler will allocate.
set(GMX_THREAD_MPI OFF CACHE BOOL "GROMACS bundled thread-MPI is not supported on BlueGene" FORCE)
set(GMX_MPI ON CACHE BOOL "MPI is required on BlueGene" FORCE)

# Access to /etc/passwd is not available on the back end of BlueGeneP
# (at least), despite being detected by CMake. This can cause linker
# warnings about harmless things in src/gromacs/utility/cstringutil.h.
set(HAVE_PWD_H OFF)

# The automatic testing for endianness does not work for the BlueGene cross-compiler
set(GMX_FLOAT_FORMAT_IEEE754 1 CACHE INTERNAL "" FORCE)
set(GMX_IEEE754_BIG_ENDIAN_BYTE_ORDER 1 CACHE INTERNAL "BlueGene has big-endian floating-point byte order (by default)" FORCE)
set(GMX_IEEE754_BIG_ENDIAN_WORD_ORDER 1 CACHE INTERNAL "BlueGene has big-endian floating-point word order (by default)" FORCE)
