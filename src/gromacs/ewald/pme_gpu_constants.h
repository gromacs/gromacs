/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief This file defines the PME GPU compile-time constants/macros,
 * used both in device and host code.
 *
 * As OpenCL C is not aware of constexpr, most of this file is
 * forwarded to the OpenCL kernel compilation as defines with same
 * names, for the sake of code similarity.
 *
 * \todo The values are currently common to both CUDA and OpenCL
 * implementations, but should be reconsidered when we tune the OpenCL
 * implementation. See Redmine #2528.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_CONSTANTS_H
#define GMX_EWALD_PME_GPU_CONSTANTS_H

#include "config.h"

#if GMX_GPU == GMX_GPU_CUDA
#    include "gromacs/gpu_utils/cuda_arch_utils.cuh" // for warp_size
#endif

/* General settings for PME GPU behaviour */

/*! \brief
 * false: The atom data GPU buffers are sized precisely according to the number of atoms.
 *        (Except GPU spline data layout which is regardless intertwined for 2 atoms per warp).
 *        The atom index checks in the spread/gather code potentially hinder the performance.
 * true:  The atom data GPU buffers are padded with zeroes so that the possible number of atoms
 *        fitting in is divisible by c_pmeAtomDataAlignment.
 *        The atom index checks are not performed. There should be a performance win, but how big is it, remains to be seen.
 *        Additional cudaMemsetAsync calls are done occasionally (only charges/coordinates; spline data is always recalculated now).
 * \todo Estimate performance differences
 */
constexpr bool c_usePadding = true;

/*! \brief
 * false: Atoms with zero charges are processed by PME. Could introduce some overhead.
 * true:  Atoms with zero charges are not processed by PME. Adds branching to the spread/gather.
 *        Could be good for performance in specific systems with lots of neutral atoms.
 * \todo Estimate performance differences.
 */
constexpr bool c_skipNeutralAtoms = false;

/*! \brief
 * Number of PME solve output floating point numbers.
 * 6 for symmetric virial matrix + 1 for reciprocal energy.
 */
constexpr int c_virialAndEnergyCount = 7;


/* Macros concerning the data layout */

/*
    Here is a current memory layout for the theta/dtheta B-spline float parameter arrays.
    This is the data in global memory used both by spreading and gathering kernels (with same scheduling).
    This example has PME order 4 and 2 particles per warp/data chunk.
    Each particle has 16 threads assigned to it, each thread works on 4 non-sequential global grid contributions.

    ----------------------------------------------------------------------------
    particles 0, 1                                        | particles 2, 3     | ...
    ----------------------------------------------------------------------------
    order index 0           | index 1 | index 2 | index 3 | order index 0 .....
    ----------------------------------------------------------------------------
    tx0 tx1 ty0 ty1 tz0 tz1 | ..........
    ----------------------------------------------------------------------------

    Each data chunk for a single warp is 24 floats. This goes both for theta and dtheta.
    24 = 2 particles per warp * order 4 * 3 dimensions. 48 floats (1.5 warp size) per warp in total.
    I have also tried intertwining theta and theta in a single array (they are used in pairs in gathering stage anyway)
    and it didn't seem to make a performance difference.

    The spline indexing is isolated in the 2 inline functions:
    getSplineParamIndexBase() return a base shared memory index corresponding to the atom in the block;
    getSplineParamIndex() consumes its results and adds offsets for dimension and spline value index.

    The corresponding defines follow.
 */

/*! \brief PME order parameter
 *
 *  Note that the GPU code, unlike the CPU, only supports order 4.
 */
constexpr int c_pmeGpuOrder = 4;

/*! \brief
 * The number of GPU threads used for computing spread/gather contributions of a single atom as function of the PME order.
 * The assumption is currently that any thread processes only a single atom's contributions.
 * TODO: this assumption leads to minimum execution width of 16. See Redmine #2516
 */
constexpr int c_pmeSpreadGatherThreadsPerAtom = c_pmeGpuOrder * c_pmeGpuOrder;

//! Number of threads per atom when order threads are used
constexpr int c_pmeSpreadGatherThreadsPerAtom4ThPerAtom = c_pmeGpuOrder;

/*! \brief Minimum execution width of the PME spread and gather kernels.
 *
 * Due to the one thread per atom and order=4 implementation constraints, order^2 threads
 * should execute without synchronization needed. See c_pmeSpreadGatherThreadsPerAtom
 */
constexpr int c_pmeSpreadGatherMinWarpSize = c_pmeSpreadGatherThreadsPerAtom;

//! Minimum warp size if order threads pera atom are used instead of order^2
constexpr int c_pmeSpreadGatherMinWarpSize4ThPerAtom = c_pmeSpreadGatherThreadsPerAtom4ThPerAtom;

/*! \brief
 * Atom data alignment (in terms of number of atoms).
 * This is the least common multiple of number of atoms processed by
 * a single block/workgroup of the spread and gather kernels.
 * If the GPU atom data buffers are padded (c_usePadding == true),
 * Then the numbers of atoms which would fit in the padded GPU buffers have to be divisible by this.
 * There are debug asserts for this divisibility in pme_gpu_spread() and pme_gpu_gather().
 */
constexpr int c_pmeAtomDataAlignment = 64;

/*
 * The execution widths for PME GPU kernels, used both on host and device for correct scheduling.
 * TODO: those were tuned for CUDA with assumption of warp size 32; specialize those for OpenCL
 * (Redmine #2528).
 * As noted below, these are very approximate maximum sizes; in run time we might have to use
 * smaller block/workgroup sizes, depending on device capabilities.
 */

//! Spreading max block width in warps picked among powers of 2 (2, 4, 8, 16) for max. occupancy and min. runtime in most cases
constexpr int c_spreadMaxWarpsPerBlock = 8;

//! Solving kernel max block width in warps picked among powers of 2 (2, 4, 8, 16) for max.
//! occupancy and min. runtime (560Ti (CC2.1), 660Ti (CC3.0) and 750 (CC5.0)))
constexpr int c_solveMaxWarpsPerBlock = 8;

//! Gathering max block width in warps - picked empirically among 2, 4, 8, 16 for max. occupancy and min. runtime
constexpr int c_gatherMaxWarpsPerBlock = 4;


#if GMX_GPU == GMX_GPU_CUDA

/* All the guys below are dependent on warp_size and should ideally be removed from the host-side code,
 * as we have to do that for OpenCL already.
 * They also express maximum desired block/workgroup sizes, while both with CUDA and OpenCL we have to treat
 * the device runtime limitations gracefully as well.
 */

/*! \brief
 * The number of atoms processed by a single warp in spread/gather.
 * This macro depends on the templated order parameter (2 atoms per warp for order 4 and warp_size
 * of 32). It is mostly used for spline data layout tweaked for coalesced access.
 */
constexpr int c_pmeSpreadGatherAtomsPerWarp = (warp_size / c_pmeSpreadGatherThreadsPerAtom);

//! number of atoms per warp when order threads are used per atom
constexpr int c_pmeSpreadGatherAtomsPerWarp4ThPerAtom =
        (warp_size / c_pmeSpreadGatherThreadsPerAtom4ThPerAtom);

//! Spreading max block size in threads
constexpr int c_spreadMaxThreadsPerBlock = c_spreadMaxWarpsPerBlock * warp_size;

//! Solving kernel max block size in threads
constexpr int c_solveMaxThreadsPerBlock = (c_solveMaxWarpsPerBlock * warp_size);

//! Gathering max block size in threads
constexpr int c_gatherMaxThreadsPerBlock = c_gatherMaxWarpsPerBlock * warp_size;
//! Gathering min blocks per CUDA multiprocessor
constexpr int c_gatherMinBlocksPerMP = GMX_CUDA_MAX_THREADS_PER_MP / c_gatherMaxThreadsPerBlock;

#endif // GMX_GPU == GMX_GPU_CUDA

#endif
