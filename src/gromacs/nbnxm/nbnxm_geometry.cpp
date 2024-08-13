/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

#include "gmxpre.h"

#include "nbnxm_geometry.h"

#include <cmath>

#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

#include "pairlist.h"

namespace gmx
{

/* Clusters at the cut-off only increase rlist by 60% of their size */
static constexpr real c_nbnxnRlistIncreaseOutsideFactor = 0.6;

real nbnxmPairlistVolumeRadiusIncrease(const bool useGpu, const real atomDensity)
{
    /* Note that this routine assumes the most common cluster layout.
     * The results should only be used for cost estimates and therefore
     * do not need to be exact.
     */

    int iClusterSize;
    int jClusterSize;
    if (useGpu)
    {
        iClusterSize = c_nbnxnGpuClusterSize;
        // The most common value is 4 (8 atoms split in 2 halves)
        jClusterSize = 4;
    }
    else
    {
        // The i-cluster size is currently always 4
        iClusterSize = 4;
#if GMX_SIMD_HAVE_REAL && ((GMX_SIMD_REAL_WIDTH == 8 && GMX_SIMD_HAVE_FMA) || GMX_SIMD_REAL_WIDTH > 8)
        jClusterSize = 8;
#else
        // The most common value is 4, GMX_SIMD_REAL_WIDTH==4 or plain-C 4x4
        jClusterSize = 4;
#endif
    }

    const real iVolumeIncrease = (iClusterSize - 1) / atomDensity;
    const real jVolumeIncrease = (jClusterSize - 1) / atomDensity;

    return c_nbnxnRlistIncreaseOutsideFactor * std::cbrt(iVolumeIncrease + jVolumeIncrease);
}

real nbnxn_get_rlist_effective_inc(const int clusterSize, const RVec& averageClusterBoundingBox)
{
    /* The average length of the diagonal of a sub cell */
    const real diagonal = std::sqrt(norm2(averageClusterBoundingBox));

    const real volumeRatio = (clusterSize - 1.0_real) / clusterSize;

    return c_nbnxnRlistIncreaseOutsideFactor * square(volumeRatio) * 0.5_real * diagonal;
}

} // namespace gmx
