/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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
 * implementation. See Issue #2528.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_CONSTANTS_H
#define GMX_EWALD_PME_GPU_CONSTANTS_H

#include "config.h"

#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/cuda_arch_utils.cuh" // for warp_size
#endif

/* General settings for PME GPU behaviour */

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

/*! \brief The number of GPU threads used for computing spread/gather
 * contributions of a single atom, which relates to the PME order.
 *
 * TODO: this assumption leads to minimum execution width of 16. See Issue #2516
 */
enum class ThreadsPerAtom : int
{
    /*! \brief Use a number of threads equal to the PME order (ie. 4)
     *
     * Only CUDA implements this. See Issue #2516 */
    Order,
    //! Use a number of threads equal to the square of the PME order (ie. 16)
    OrderSquared,
    //! Size of the enumeration
    Count
};

/*
 * The execution widths for PME GPU kernels, used both on host and device for correct scheduling.
 * TODO: those were tuned for CUDA with assumption of warp size 32; specialize those for OpenCL
 * (Issue #2528).
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

#if GMX_GPU_CUDA
/* All the fields below are dependent on warp_size and should
 * ideally be removed from the device-side code, as we have to
 * do that for OpenCL already.
 *
 * They also express maximum desired block/workgroup sizes,
 * while both with CUDA and OpenCL we have to treat the device
 * runtime limitations gracefully as well.
 */

//! Spreading max block size in threads
static constexpr int c_spreadMaxThreadsPerBlock = c_spreadMaxWarpsPerBlock * warp_size;

//! Solving kernel max block size in threads
static constexpr int c_solveMaxThreadsPerBlock = c_solveMaxWarpsPerBlock * warp_size;

//! Gathering max block size in threads
static constexpr int c_gatherMaxThreadsPerBlock = c_gatherMaxWarpsPerBlock * warp_size;
//! Gathering min blocks per CUDA multiprocessor (determined empirically to give best performance)
#    if GMX_PTX_ARCH >= 800
static constexpr int c_gatherMinBlocksPerMP = 12;
#    else
static constexpr int c_gatherMinBlocksPerMP = GMX_CUDA_MAX_THREADS_PER_MP / c_gatherMaxThreadsPerBlock;
#    endif

#endif // GMX_GPU_CUDA

#endif
