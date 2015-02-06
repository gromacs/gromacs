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

//! Test fixture
class Simd4RvecOperationsTest : public SimdTest
{
    public:
        Simd4RvecOperationsTest()
            : xPairs_(xPairsBuffer_),
              pairIndex_(pairIndexBuffer_)
        {
            alignedPtr_ = gmx_simd_align_r(alignedBuffer_);
            alignedExpectedResultsPtr_ = gmx_simd_align_r(expectedResultsBuffer_);

            // Create some coordinates
            for (unsigned int i = 0; i < xPairs_.size(); ++i)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    xPairs_[i][d] = i*2 + d;
                }
            }

            // Scatter some pair indices around the buffer
            for (unsigned int i = 0; i < pairIndex_.size(); i += 2)
            {
                int size = xPairs_.size();
                pairIndex_[i]   = size-2 - i;
                pairIndex_[i+1] = 2*i+1;
                GMX_RELEASE_ASSERT(pairIndex_[i] < size, "Invalid entry in pairIndex_");
                GMX_RELEASE_ASSERT(pairIndex_[i+1] < size, "Invalid entry in pairIndex_");

                // Calculate the expected results without SIMD
                rvec difference, a, b;

                copy_rvec(xPairs_[pairIndex_[i]], a);
                copy_rvec(xPairs_[pairIndex_[i+1]], b);
                rvec_sub(a, b, difference);
                for (int d = 0; d < DIM; ++d)
                {
                    alignedExpectedResultsPtr_[d*GMX_SIMD_REAL_WIDTH + i/2] = difference[d];
                }
            }

        }

        rvec            xPairsBuffer_[GMX_SIMD_REAL_WIDTH*4];
        ArrayRef<rvec>  xPairs_;
        int             pairIndexBuffer_[GMX_SIMD_REAL_WIDTH*2];
        ArrayRef<int>   pairIndex_;
        real            alignedBuffer_[GMX_SIMD_REAL_WIDTH*4];
        real           *alignedPtr_;
        real            expectedResultsBuffer_[GMX_SIMD_REAL_WIDTH * (DIM + 1)];
        real           *alignedExpectedResultsPtr_;
};

TEST_F(Simd4RvecOperationsTest, gmxSimdGatherRvecDistPairIndex)
{
    gmx_simd_real_t rSimd[DIM], rSimdExpected[DIM];

    // Get the expected results into SIMD registers
    for (int d = 0; d < DIM; ++d)
    {
        rSimdExpected[d] = gmx_simd_load_r(alignedExpectedResultsPtr_ + GMX_SIMD_REAL_WIDTH*d);
    }

    // Call the function being tested
    gmx_simd_gather_rvec_dist_pair_index(xPairs_.data(),
                                         pairIndex_.data(),
                                         alignedPtr_,
                                         &rSimd[0],
                                         &rSimd[1],
                                         &rSimd[2]);

    for (int d = 0; d < DIM; ++d)
    {
        GMX_EXPECT_SIMD_REAL_EQ(rSimdExpected[d], rSimd[d]);
    }
}

TEST_F(Simd4RvecOperationsTest, gmxSimdGatherRvecDistTwoIndex)
{
    int firstBuffer[GMX_SIMD_REAL_WIDTH], secondBuffer[GMX_SIMD_REAL_WIDTH];
    ArrayRef<int> first(firstBuffer), second(secondBuffer);
    gmx_simd_real_t rSimd[DIM], rSimdExpected[DIM];

    // Convert the indices into the form that this function expects
    for (unsigned int i = 0; i != first.size() && i != second.size(); ++i)
    {
        first[i] = pairIndex_[i*2+0];
        second[i] = pairIndex_[i*2+1];
    }

    // Get the expected results into SIMD registers
    for (int d = 0; d < DIM; ++d)
    {
        rSimdExpected[d] = gmx_simd_load_r(alignedExpectedResultsPtr_ + GMX_SIMD_REAL_WIDTH*d);
    }

    // Call the function being tested
    gmx_simd_gather_rvec_dist_two_index(xPairs_.data(),
                                        first.data(),
                                        second.data(),
                                        alignedPtr_,
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
