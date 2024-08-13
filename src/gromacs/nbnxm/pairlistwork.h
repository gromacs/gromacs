/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include "boundingbox.h"
#include "pairlist.h"

namespace gmx
{

//! Working data for the actual i-supercell during pair search \internal
struct NbnxmPairlistCpuWork
{
    NbnxmPairlistCpuWork(const int iClusterSize) : iClusterData(iClusterSize) {}

    //! Struct for storing coordinates and bounding box for an i-entry during search \internal
    struct IClusterData
    {
        IClusterData(const int iClusterSize) :
            bb(1), x(iClusterSize * DIM), xSimd(iClusterSize * DIM * GMX_REAL_MAX_SIMD_WIDTH)
        {
        }

        //! The bounding boxes, pbc shifted, for each cluster
        AlignedVector<BoundingBox> bb;
        //! The coordinates, pbc shifted, for each atom
        std::vector<real> x;
        //! Aligned list for storing 4*DIM*GMX_SIMD_REAL_WIDTH reals
        AlignedVector<real> xSimd;
    };

    //! Protect data from cache pollution between threads
    gmx_cache_protect_t cp0;

    //! Work data for generating an IEntry in the pairlist
    IClusterData iClusterData;
    //! Temporary j-cluster list, used for sorting on exclusions
    std::vector<nbnxn_cj_t> cj;

    //! Nr. of cluster pairs without Coulomb for flop counting
    int ncj_noq = 0;
    //! Nr. of cluster pairs with 1/2 LJ for flop count
    int ncj_hlj = 0;

    //! Protect data from cache pollution between threads
    gmx_cache_protect_t cp1;
};

/* Working data for the actual i-supercell during pair search */
struct NbnxmPairlistGpuWork
{
    struct ISuperClusterData
    {
        ISuperClusterData();

        //! The bounding boxes, pbc shifted, for each cluster
        AlignedVector<BoundingBox> bb;
        //! As bb, but in packed xxxx format
        AlignedVector<float> bbPacked;
        //! The coordinates, pbc shifted, for each atom
        AlignedVector<real> x;
        //! Aligned coordinate list used for 4*DIM*GMX_SIMD_REAL_WIDTH floats
        AlignedVector<real> xSimd;
    };

    NbnxmPairlistGpuWork();

    //! Protect data from cache pollution between threads
    gmx_cache_protect_t cp0;

    //! Work data for generating an i-entry in the pairlist
    ISuperClusterData iSuperClusterData;
    //! The current j-cluster index for the current list
    int cj_ind;
    //! Bounding box distance work array
    AlignedVector<float> distanceBuffer;

    //! Buffer for sorting list entries
    std::vector<int> sortBuffer;

    //! Second sci array, for sorting
    HostVector<nbnxn_sci_t> sci_sort;

    //! Protect data from cache pollution between threads
    gmx_cache_protect_t cp1;
};

} // namespace gmx

#endif
