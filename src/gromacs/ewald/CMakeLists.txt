#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

gmx_add_libgromacs_sources(
    calculate-spline-moduli.cpp
    ewald.cpp
    ewald-utils.cpp
    long-range-correction.cpp
    pme.cpp
    pme-gather.cpp
    pme-grid.cpp
    pme-load-balancing.cpp
    pme-only.cpp
    pme-pp.cpp
    pme-redistribute.cpp
    pme-solve.cpp
    pme-spline-work.cpp
    pme-spread.cpp
    # Files that implement stubs
    pme-gpu-program.cpp
    )
if (GMX_USE_CUDA)
    gmx_add_libgromacs_sources(
        # CUDA-specific sources
        pme-gather.cu
        pme-gpu-3dfft.cu
        pme-solve.cu
        pme-spread.cu
        pme-gpu-program-impl.cu
        # GPU-specific sources
        pme-gpu.cpp
        pme-gpu-internal.cpp
        pme-gpu-timings.cpp
        )
    gmx_compile_cpp_as_cuda(
        pme-gpu-internal.cpp
        pme-gpu-program.cpp
        pme-gpu-timings.cpp
        )
elseif (GMX_USE_OPENCL)
    gmx_add_libgromacs_sources(
        # OpenCL-specific sources
        pme-gpu-3dfft-ocl.cpp
        pme-gpu-program-impl-ocl.cpp
        # GPU-specific sources
        pme-gpu.cpp
        pme-gpu-internal.cpp
        pme-gpu-timings.cpp
        )
    file(GLOB EWALD_OPENCL_RUNTIME_SOURCES pme-spread.clh pme-solve.clh pme-gather.clh pme-gpu-utils.clh pme-program.cl pme-gpu-types.h)
    set(MDLIB_OPENCL_KERNELS ${MDLIB_OPENCL_KERNELS} ${EWALD_OPENCL_RUNTIME_SOURCES} PARENT_SCOPE)
else()
    gmx_add_libgromacs_sources(
        # Files that implement stubs
        pme-gpu-program-impl.cpp
        )
endif()

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
