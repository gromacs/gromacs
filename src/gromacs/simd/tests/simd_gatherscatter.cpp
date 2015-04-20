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

#include "simd.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#if defined(GMX_SIMD_HAVE_REAL) && defined(GMX_SIMD_EXTENDED)

/*! \brief Test fixture for floating-point gather and scatter tests.
 *
 * Inherit from main SimdTest, add code to generate spread-out triplet data.
 */
class SimdGatherScatterTest : public SimdTest
{
    public:
        SimdGatherScatterTest()
        {
            int i;
            // Fill memory we should not read/write
            // We have a total of GMX_SIMD_REAL_WIDTH triplets, padded with as
            // many unused ones (i.e, 2*3*width elements)
            for (i = 0; i < 6*GMX_SIMD_REAL_WIDTH; i++)
            {
                mem_[i] = background_;
            }
            // Create aligned pointers
            m0_ = gmx_simd_align_r(work_);
            m1_ = m0_ + GMX_SIMD_REAL_WIDTH;
            m2_ = m0_ + 2*GMX_SIMD_REAL_WIDTH;

            for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
            {
                // Use every second element so we really check noncontiguous loads
                offset_[i] = 2*i;
                // Integer part is the index, decimal indicates 1st/2nd/3rd component
                // Put the same data in the three aligned memory locations (meaning, transposed)
                mem_[3*offset_[i]]   = m0_[i] = i;
                mem_[3*offset_[i]+1] = m1_[i] = i + 0.1;
                mem_[3*offset_[i]+2] = m2_[i] = i + 0.2;
            }
        }

    protected:
        // Room for GMX_SIMD_REAL_WIDTH triplets, and single real before/between/after.
        real                    mem_[6*GMX_SIMD_REAL_WIDTH];    //!< Mem for checking gather/scatter (3*2*width)
        int                     offset_[GMX_SIMD_REAL_WIDTH];   //!< Array with offsets for gather/scatter ops
        real *                  m0_;                            //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  m1_;                            //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  m2_;                            //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        const real              background_ = -9999.99;         //!< Value for memory we should not touch
    private:
        // Use enough work mem so the aligned-pair-alignment can be 32*GMX_SIMD_REAL_WIDTH.
        real                    work_[33*GMX_SIMD_REAL_WIDTH]; //!< Raw work memory for aligned load/store
};


TEST_F(SimdGatherScatterTest, gmxSimdGatherLoad3TransposeR)
{
    gmx_simd_real_t v0, v1, v2;
    gmx_simd_real_t ref0, ref1, ref2;

    ref0 = gmx_simd_load_r(m0_);
    ref1 = gmx_simd_load_r(m1_);
    ref2 = gmx_simd_load_r(m2_);
    gmx_simd_gather_load_3_transpose_r(mem_, offset_, &v0, &v1, &v2);

    GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
    GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    GMX_EXPECT_SIMD_REAL_EQ(ref2, v2);
}

TEST_F(SimdGatherScatterTest, gmxSimdTransposeScatterStore3R)
{
    gmx_simd_real_t v0, v1, v2;
    real            testmem[6*GMX_SIMD_REAL_WIDTH];
    int             i;

    for (i = 0; i < 6*GMX_SIMD_REAL_WIDTH; i++)
    {
        testmem[i] = background_;
    }
    v0 = gmx_simd_load_r(m0_);
    v1 = gmx_simd_load_r(m1_);
    v2 = gmx_simd_load_r(m2_);
    gmx_simd_transpose_scatter_store_3_r(testmem, offset_, v0, v1, v2);

    for (i = 0; i < 4*GMX_SIMD_REAL_WIDTH+1; i++)
    {
        EXPECT_EQ(mem_[i], testmem[i]);
    }
}

TEST_F(SimdGatherScatterTest, gmxSimdTransposeScatterIncr3R)
{
    gmx_simd_real_t                   v0, v1, v2;
    real                              testmem[6*GMX_SIMD_REAL_WIDTH];
    int                               i;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    for (i = 0; i < 6*GMX_SIMD_REAL_WIDTH; i++)
    {
        testmem[i] = mem_[i];
    }
    v0 = gmx_simd_load_r(m0_);
    v1 = gmx_simd_load_r(m1_);
    v2 = gmx_simd_load_r(m2_);

    // Increment the vectors with themselves
    gmx_simd_transpose_scatter_incr_3_r(testmem, offset_, v0, v1, v2);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        // Check that vectors are now double value
        EXPECT_REAL_EQ_TOL(2*mem_[3*offset_[i]], testmem[3*offset_[i]], tolerance);
        EXPECT_REAL_EQ_TOL(2*mem_[3*offset_[i]+1], testmem[3*offset_[i]+1], tolerance);
        EXPECT_REAL_EQ_TOL(2*mem_[3*offset_[i]+2], testmem[3*offset_[i]+2], tolerance);
        // Each triplet should be followed by three exactly untouched positions
        EXPECT_EQ(background_, testmem[3*offset_[i]+3]);
        EXPECT_EQ(background_, testmem[3*offset_[i]+4]);
        EXPECT_EQ(background_, testmem[3*offset_[i]+5]);
    }
}

TEST_F(SimdGatherScatterTest, gmxSimdTransposeScatterDecr3R)
{
    gmx_simd_real_t                   v0, v1, v2;
    real                              testmem[6*GMX_SIMD_REAL_WIDTH];
    int                               i;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    // Fill test memory with background
    for (i = 0; i < 6*GMX_SIMD_REAL_WIDTH; i++)
    {
        testmem[i] = background_;
    }
    // Set the live vector elements to 4*reference data
    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        testmem[3*offset_[i]]   = 4*mem_[3*offset_[i]];
        testmem[3*offset_[i]+1] = 4*mem_[3*offset_[i]+1];
        testmem[3*offset_[i]+2] = 4*mem_[3*offset_[i]+2];
    }

    v0 = gmx_simd_load_r(m0_);
    v1 = gmx_simd_load_r(m1_);
    v2 = gmx_simd_load_r(m2_);

    // Subtract the vector data from testmem where we have 4*data, to get 3*data
    gmx_simd_transpose_scatter_decr_3_r(testmem, offset_, v0, v1, v2);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        // Check that vectors are now 4*data-data=3*data
        EXPECT_REAL_EQ_TOL(3*mem_[3*offset_[i]], testmem[3*offset_[i]], tolerance);
        EXPECT_REAL_EQ_TOL(3*mem_[3*offset_[i]+1], testmem[3*offset_[i]+1], tolerance);
        EXPECT_REAL_EQ_TOL(3*mem_[3*offset_[i]+2], testmem[3*offset_[i]+2], tolerance);
        // Each triplet should be followed by three exactly untouched positions
        EXPECT_EQ(background_, testmem[3*offset_[i]+3]);
        EXPECT_EQ(background_, testmem[3*offset_[i]+4]);
        EXPECT_EQ(background_, testmem[3*offset_[i]+5]);
    }
}


#endif      // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace
