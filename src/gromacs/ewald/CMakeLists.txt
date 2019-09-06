#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
    calculate_spline_moduli.cpp
    ewald.cpp
    ewald_utils.cpp
    long_range_correction.cpp
    pme.cpp
    pme_gather.cpp
    pme_grid.cpp
    pme_load_balancing.cpp
    pme_only.cpp
    pme_pp.cpp
    pme_redistribute.cpp
    pme_solve.cpp
    pme_spline_work.cpp
    pme_spread.cpp
    # Files that implement stubs
    pme_gpu_program.cpp
    pme_pp_comm_gpu_impl.cpp
    pme_coordinate_receiver_gpu_impl.cpp
    pme_force_sender_gpu_impl.cpp
    )
if (GMX_USE_CUDA)
    gmx_add_libgromacs_sources(
        # CUDA-specific sources
        pme_gather.cu
        pme_gpu_3dfft.cu
        pme_solve.cu
        pme_spread.cu
        pme_gpu_program_impl.cu
        pme_pp_comm_gpu_impl.cu
        pme_force_sender_gpu_impl.cu
        pme_coordinate_receiver_gpu_impl.cu
        # GPU-specific sources
        pme_gpu.cpp
        pme_gpu_internal.cpp
        pme_gpu_timings.cpp
        )
    gmx_compile_cpp_as_cuda(
        pme_gpu_internal.cpp
        pme_gpu_program.cpp
        pme_gpu_timings.cpp
        )
elseif (GMX_USE_OPENCL)
    gmx_add_libgromacs_sources(
        # OpenCL-specific sources
        pme_gpu_3dfft_ocl.cpp
        pme_gpu_program_impl_ocl.cpp
        # GPU-specific sources
        pme_gpu.cpp
        pme_gpu_internal.cpp
        pme_gpu_timings.cpp
        )
else()
    gmx_add_libgromacs_sources(
        # Files that implement stubs
        pme_gpu_program_impl.cpp
        )
endif()

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
