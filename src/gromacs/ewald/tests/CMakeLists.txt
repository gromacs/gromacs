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

if (GMX_GPU_FFT_VKFFT)
  # Suppress warning that first appeared when compiling tests on macOS with x86.
  check_cxx_compiler_flag("-Wno-unused-parameter" HAS_NO_UNUSED_PARAMETER)
  if(HAS_NO_UNUSED_PARAMETER)
    set_source_files_properties(pmetestcommon.cpp PROPERTIES COMPILE_OPTIONS "-Wno-unused-parameter")
  endif()
endif()

gmx_add_unit_test(EwaldUnitTests ewald-test
    HARDWARE_DETECTION
    DYNAMIC_REGISTRATION
    CPP_SOURCE_FILES
        pmebsplinetest.cpp
        pmegathertest.cpp
        pmesolvetest.cpp
        pmesplinespreadtest.cpp
        pme.cpp
    GPU_CPP_SOURCE_FILES
        pmetestcommon.cpp
)
target_link_libraries(
    ewald-test PRIVATE
        ewald
        domdec
        fft
        gpu_utils
        math
        mdtypes
        pbcutil
        testutils
        topology
        utility
)
