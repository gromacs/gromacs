/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/*! \internal \file
 *  \brief Top-level file for generating PME OpenCL kernels.
 * This file includes all the PME OpenCL kernel files multiple times, with additional defines.
 * Plenty of common constants/macros are also expected to be already defined by the compiler
 * (such as "order" which is PME interpolation order).
 * For details, please see how pme-program.cl is compiled in pme-gpu-program-impl-ocl.cpp.
 *
 * This file specifically expects the following work size definitions:
 *
 * - c_spreadWorkGroupSize is a spline/spread kernels work group size.
 * - c_solveMaxWorkGroupSize is a solve kernels maximum work group size.
 * - c_gatherWorkGroupSize is a gather kernels work group size.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

// Assert placeholders, to not rip them out from OpenCL implementation - hopefully they come in handy some day with OpenCL 2
#define static_assert(a, b)
#define assert(a)

#define PmeOpenCLKernelParams PmeGpuKernelParamsBase

/*
 * A section of common helpers, which can also be implemented for CUDA to extract obvious similarities of OpenCL/CUDA kernels.
 * This is non PME-specific and can ideally be moved into a separate utilities file.
 */

//! Thread index in any dimension of the block/workgroup. Corresponds to CUDA threadIdx.x/y/z
inline size_t getThreadLocalIndex(size_t dimIndex)
{
    return get_local_id(dimIndex);
}

//! Thread index in a 2D block/workgroup. Corresponds to CUDA (threadIdx.y * blockDim.x) + threadIdx.x
inline size_t getThreadLocalIndex2d()
{
    return get_local_id(1) * get_local_size(0) + get_local_id(0);
}

//! Thread index in a 3D block/workgroup. Corresponds to CUDA ((threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x) + threadIdx.x
inline size_t getThreadLocalIndex3d()
{
    return (get_local_id(2) * get_local_size(1) + get_local_id(1)) * get_local_size(0) + get_local_id(0);
}

//! Dimension of the block/workgroup - corresponds to CUDA blockDim.x/y/z
inline size_t getBlockSize(size_t dimIndex)
{
    return get_local_size(dimIndex);
}

//! Index of the block/workgroup - corresponds to CUDA blockIdx.x/y/z
inline size_t getBlockIndex(size_t dimIndex)
{
    return get_group_id(dimIndex);
}

/* Kernel instantiations */


/* SPREAD/SPLINE */

#define atomsPerBlock  (c_spreadWorkGroupSize / PME_SPREADGATHER_THREADS_PER_ATOM)

// spline/spread fused
#define computeSplines 1
#define spreadCharges 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineAndSpreadKernel
#include "pme-spread.clh"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME

// spline
#define computeSplines 1
#define spreadCharges 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineKernel
#include "pme-spread.clh"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME

// spread
#define computeSplines 0
#define spreadCharges 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSpreadKernel
#include "pme-spread.clh"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME


#undef atomsPerBlock


/* GATHER */

#define atomsPerBlock (c_gatherWorkGroupSize / PME_SPREADGATHER_THREADS_PER_ATOM)

// gather
#define overwriteForces 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeGatherKernel
#include "pme-gather.clh"
#undef overwriteForces
#undef CUSTOMIZED_KERNEL_NAME

// gather with reduction
#define overwriteForces 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeGatherReduceWithInputKernel
#include "pme-gather.clh"
#undef overwriteForces
#undef CUSTOMIZED_KERNEL_NAME


#undef atomsPerBlock

/* SOLVE */

// GridOrdering enum replacement
#define YZX 1
#define XYZ 2

// solve, YZX dimension order
#define gridOrdering YZX
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXKernel
#include "pme-solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, YZX dimension order
#define gridOrdering YZX
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXEnergyKernel
#include "pme-solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve, XYZ dimension order
#define gridOrdering XYZ
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZKernel
#include "pme-solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, XYZ dimension order
#define gridOrdering XYZ
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZEnergyKernel
#include "pme-solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME
