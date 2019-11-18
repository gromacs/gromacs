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
 * \brief Declares working data structures for the CPU and GPU pairlists
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_PAIRLISTWORK_H
#define GMX_NBNXM_PAIRLISTWORK_H

#include <memory>
#include <vector>

#include "gromacs/simd/simd.h"

#include "boundingboxes.h"
#include "grid.h"
#include "pairlist.h"

/* Working data for the actual i-supercell during pair search */
struct NbnxnPairlistCpuWork
{
    // Struct for storing coordinats and bounding box for an i-entry during search
    struct IClusterData
    {
        IClusterData() :
            bb(1),
            x(c_nbnxnCpuIClusterSize * DIM),
            xSimd(c_nbnxnCpuIClusterSize * DIM * GMX_REAL_MAX_SIMD_WIDTH)
        {
        }

        // The bounding boxes, pbc shifted, for each cluster
        AlignedVector<Nbnxm::BoundingBox> bb;
        // The coordinates, pbc shifted, for each atom
        std::vector<real> x;
        // Aligned list for storing 4*DIM*GMX_SIMD_REAL_WIDTH reals
        AlignedVector<real> xSimd;
    };

    // Protect data from cache pollution between threads
    gmx_cache_protect_t cp0;

    // Work data for generating an IEntry in the pairlist
    IClusterData iClusterData;
    // The current cj_ind index for the current list
    int cj_ind;
    // Temporary j-cluster list, used for sorting on exclusions
    std::vector<nbnxn_cj_t> cj;

    // Nr. of cluster pairs without Coulomb for flop counting
    int ncj_noq;
    // Nr. of cluster pairs with 1/2 LJ for flop count
    int ncj_hlj;

    // Protect data from cache pollution between threads
    gmx_cache_protect_t cp1;
};

/* Working data for the actual i-supercell during pair search */
struct NbnxnPairlistGpuWork
{
    struct ISuperClusterData
    {
        ISuperClusterData() :
            bb(c_gpuNumClusterPerCell),
#if NBNXN_SEARCH_BB_SIMD4
            bbPacked(c_gpuNumClusterPerCell / c_packedBoundingBoxesDimSize * c_packedBoundingBoxesSize),
#endif
            x(c_gpuNumClusterPerCell * c_nbnxnGpuClusterSize * DIM),
            xSimd(c_gpuNumClusterPerCell * c_nbnxnGpuClusterSize * DIM)
        {
        }

        // The bounding boxes, pbc shifted, for each cluster
        AlignedVector<Nbnxm::BoundingBox> bb;
        // As bb, but in packed xxxx format
        AlignedVector<float> bbPacked;
        // The coordinates, pbc shifted, for each atom
        AlignedVector<real> x;
        // Aligned coordinate list used for 4*DIM*GMX_SIMD_REAL_WIDTH floats
        AlignedVector<real> xSimd;
    };

    NbnxnPairlistGpuWork() :
        distanceBuffer(c_gpuNumClusterPerCell),
        sci_sort({}, { gmx::PinningPolicy::PinnedIfSupported })
    {
    }

    // Protect data from cache pollution between threads
    gmx_cache_protect_t cp0;

    // Work data for generating an i-entry in the pairlist
    ISuperClusterData iSuperClusterData;
    // The current j-cluster index for the current list
    int cj_ind;
    // Bounding box distance work array
    AlignedVector<float> distanceBuffer;

    // Buffer for sorting list entries
    std::vector<int> sortBuffer;

    // Second sci array, for sorting
    gmx::HostVector<nbnxn_sci_t> sci_sort;

    // Protect data from cache pollution between threads
    gmx_cache_protect_t cp1;
};

#endif
