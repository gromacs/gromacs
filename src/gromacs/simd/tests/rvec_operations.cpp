/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/simd/rvec_operations.h"

#include <math.h>

#include "gromacs/math/vec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

#include "simd.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#if defined GMX_SIMD_HAVE_REAL

//! Convenience typedef
typedef SimdTest Simd4RvecOperationsTest;

TEST_F(Simd4RvecOperationsTest, gmxSimdGatherRvecDistPairIndex)
{
    rvec            xPairsBuffer[GMX_SIMD_REAL_WIDTH*4];
    ArrayRef<rvec>  xPairs(xPairsBuffer);
    int             pairIndexBuffer[GMX_SIMD_REAL_WIDTH*2];
    ArrayRef<int>   pairIndex(pairIndexBuffer);
    real            alignedBuffer[GMX_SIMD_REAL_WIDTH*4];
    real           *alignedPtr = gmx_simd_align_r(alignedBuffer);
    gmx_simd_real_t rSimd[DIM], rSimdExpected[DIM];
    real            expectedResultsBuffer[GMX_SIMD_REAL_WIDTH * (DIM + 1)];
    real           *alignedExpectedResultsPtr = gmx_simd_align_r(expectedResultsBuffer);

    // Create some coordinates
    for (unsigned int i = 0; i < xPairs.size(); ++i)
    {
        for (int d = 0; d < DIM; ++d)
        {
            xPairs[i][d] = i*2 + d;
        }
    }

    // Scatter some pair indices around the buffer
    for (unsigned int i = 0; i < pairIndex.size(); i += 2)
    {
        int size = xPairs.size();
        pairIndex[i]   = size-2 - i;
        pairIndex[i+1] = 2*i+1;
        GMX_RELEASE_ASSERT(pairIndex[i] < size, "Invalid entry in pairIndex");
        GMX_RELEASE_ASSERT(pairIndex[i+1] < size, "Invalid entry in pairIndex");

        // Calculate the expected results without SIMD
        rvec difference, a, b;

        copy_rvec(xPairs[pairIndex[i]], a);
        copy_rvec(xPairs[pairIndex[i+1]], b);
        rvec_sub(a, b, difference);
        for (int d = 0; d < DIM; ++d)
        {
            alignedExpectedResultsPtr[d*GMX_SIMD_REAL_WIDTH + i/2] = difference[d];
        }
    }

    /* Get the expected results into SIMD registers */
    for (int d = 0; d < DIM; ++d)
    {
        rSimdExpected[d] = gmx_simd_load_r(alignedExpectedResultsPtr + GMX_SIMD_REAL_WIDTH*d);
    }

    // Call the function being tested
    gmx_simd_gather_rvec_dist_pair_index(xPairs.data(),
                                         pairIndex.data(),
                                         alignedPtr,
                                         &rSimd[0],
                                         &rSimd[1],
                                         &rSimd[2]);

    for (int d = 0; d < DIM; ++d)
    {
        GMX_EXPECT_SIMD_REAL_EQ(rSimdExpected[d], rSimd[d]);
    }
}

#endif      // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace
