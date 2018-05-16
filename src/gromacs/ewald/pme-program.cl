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

//! A helper for synchronizing writes to block-visible memory. Corresponds roughly to CUDA __syncthreads().
inline void localMemoryBarrier()
{
    barrier(CLK_LOCAL_MEM_FENCE);
}

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

/* A couple of PME-specific utility functions */

/*! \brief
 * A function for checking the global atom data indices against the atom data array sizes.
 *
 * \param[in] atomDataIndexGlobal  The atom data index.
 * \param[in] nAtomData            The atom data array element count.
 * \returns                        Non-0 if index is within bounds (or PME data padding is enabled), 0 otherwise.
 *
 * This is called from the spline_and_spread and gather PME kernels.
 * The goal is to isolate the global range checks, and allow avoiding them with c_usePadding being true.
 */
int inline pme_gpu_check_atom_data_index(const size_t atomDataIndex, const size_t nAtomData)
{
    return c_usePadding ? 1 : (atomDataIndex < nAtomData);
}

/*! \brief
 * A function for optionally skipping neutral charges, depending on c_skipNeutralAtoms.
 *
 * \param[in] coefficient     The atom charge/coefficient.
 * \returns                   Non-0 if atom should be processed, 0 otherwise.
 */
int inline pme_gpu_check_atom_charge(const float coefficient)
{
    return c_skipNeutralAtoms ? (coefficient != 0.0f) : 1;
}


/* Kernel instantiations */


/* SPREAD/SPLINE */

#define atomsPerBlock  (c_spreadMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM)

// spline/spread fused
#define computeSplines 1
#define spreadCharges 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineAndSpreadKernel
#include "pme-spread.cl"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME

// spline
#define computeSplines 1
#define spreadCharges 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineKernel
#include "pme-spread.cl"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME

// spread
#define computeSplines 0
#define spreadCharges 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSpreadKernel
#include "pme-spread.cl"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME


#undef atomsPerBlock


/* GATHER */

#define atomsPerBlock (c_gatherMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM)

// gather
#define overwriteForces 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeGatherKernel
#include "pme-gather.cl"
#undef overwriteForces
#undef CUSTOMIZED_KERNEL_NAME

// gather with reduction
#define overwriteForces 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeGatherReduceWithInputKernel
#include "pme-gather.cl"
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
#include "pme-solve.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, YZX dimension order
#define gridOrdering YZX
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXEnergyKernel
#include "pme-solve.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve, XYZ dimension order
#define gridOrdering XYZ
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZKernel
#include "pme-solve.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, XYZ dimension order
#define gridOrdering XYZ
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZEnergyKernel
#include "pme-solve.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME
