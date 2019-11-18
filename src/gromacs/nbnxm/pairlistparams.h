/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 *
 * \brief
 * Declares the PairlistType enum and PairlistParams class
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_PAIRLISTPARAMS_H
#define GMX_NBNXM_PAIRLISTPARAMS_H

#include "config.h"

#include "gromacs/mdtypes/locality.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

namespace Nbnxm
{
enum class KernelType;
}

//! The i-cluster size for CPU kernels, always 4 atoms
static constexpr int c_nbnxnCpuIClusterSize = 4;

//! The i- and j-cluster size for GPU lists, 8 atoms for CUDA, set at compile time for OpenCL
#if GMX_GPU == GMX_GPU_OPENCL
static constexpr int c_nbnxnGpuClusterSize = GMX_OPENCL_NB_CLUSTER_SIZE;
#else
static constexpr int c_nbnxnGpuClusterSize = 8;
#endif

//! The number of clusters along Z in a pair-search grid cell for GPU lists
static constexpr int c_gpuNumClusterPerCellZ = 2;
//! The number of clusters along Y in a pair-search grid cell for GPU lists
static constexpr int c_gpuNumClusterPerCellY = 2;
//! The number of clusters along X in a pair-search grid cell for GPU lists
static constexpr int c_gpuNumClusterPerCellX = 2;
//! The number of clusters in a pair-search grid cell for GPU lists
static constexpr int c_gpuNumClusterPerCell =
        c_gpuNumClusterPerCellZ * c_gpuNumClusterPerCellY * c_gpuNumClusterPerCellX;


/*! \brief The number of sub-parts used for data storage for a GPU cluster pair
 *
 * In CUDA the number of threads in a warp is 32 and we have cluster pairs
 * of 8*8=64 atoms, so it's convenient to store data for cluster pair halves.
 */
static constexpr int c_nbnxnGpuClusterpairSplit = 2;

//! The fixed size of the exclusion mask array for a half GPU cluster pair
static constexpr int c_nbnxnGpuExclSize =
        c_nbnxnGpuClusterSize * c_nbnxnGpuClusterSize / c_nbnxnGpuClusterpairSplit;

//! The available pair list types
enum class PairlistType : int
{
    Simple4x2,
    Simple4x4,
    Simple4x8,
    HierarchicalNxN,
    Count
};

//! Gives the i-cluster size for each pairlist type
static constexpr gmx::EnumerationArray<PairlistType, int> IClusterSizePerListType = {
    { c_nbnxnCpuIClusterSize, c_nbnxnCpuIClusterSize, c_nbnxnCpuIClusterSize, c_nbnxnGpuClusterSize }
};
//! Gives the j-cluster size for each pairlist type
static constexpr gmx::EnumerationArray<PairlistType, int> JClusterSizePerListType = {
    { 2, 4, 8, c_nbnxnGpuClusterSize }
};

/*! \internal
 * \brief The setup for generating and pruning the nbnxn pair list.
 *
 * Without dynamic pruning rlistOuter=rlistInner.
 */
struct PairlistParams
{
    /*! \brief Constructor producing a struct with dynamic pruning disabled
     */
    PairlistParams(Nbnxm::KernelType kernelType, bool haveFep, real rlist, bool haveMultipleDomains);

    PairlistType pairlistType;           //!< The type of cluster-pair list
    bool         haveFep;                //!< Tells whether we have perturbed interactions
    real         rlistOuter;             //!< Cut-off of the larger, outer pair-list
    real         rlistInner;             //!< Cut-off of the smaller, inner pair-list
    bool         haveMultipleDomains;    //!< True when using DD with multiple domains
    bool         useDynamicPruning;      //!< Are we using dynamic pair-list pruning
    int          nstlistPrune;           //!< Pair-list dynamic pruning interval
    int          numRollingPruningParts; //!< The number parts to divide the pair-list into for rolling pruning, a value of 1 gives no rolling pruning
    int          lifetime;               //!< Lifetime in steps of the pair-list
};

#endif
