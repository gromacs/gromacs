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

#include "testutils/testasserts.h"

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

#if defined(GMX_SIMD_HAVE_REAL) && defined(GMX_SIMD_EXTENDED)

/*! \brief Test fixture for higher-level floating-point utility functions.
 *
 * Inherit from main SimdTest, add code to generate aligned memory and data.
 */
class SimdFloatingpointUtilTest : public SimdTest
{
    public:
        SimdFloatingpointUtilTest()
        {
            int   i;
            int   iwork[2*GMX_SIMD_REAL_WIDTH];
            int * tmpoffset = gmx_simd_align_i(iwork);

            // Create aligned pointers
            m0_  = gmx_simd_align_r(work_);
            m1_  = m0_ +     GMX_SIMD_REAL_WIDTH;
            m2_  = m0_ + 2 * GMX_SIMD_REAL_WIDTH;
            m3_  = m0_ + 3 * GMX_SIMD_REAL_WIDTH;
            mem_ = m0_ + 4 * GMX_SIMD_REAL_WIDTH;

            for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
            {
                tmpoffset[i]               = i;
                mem_[4 * tmpoffset[i]]     = m0_[i] = i;
                mem_[4 * tmpoffset[i] + 1] = m1_[i] = i + 0.1;
                mem_[4 * tmpoffset[i] + 2] = m2_[i] = i + 0.2;
                mem_[4 * tmpoffset[i] + 3] = m3_[i] = i + 0.3;
            }
            simdoffset_ = gmx_simd_load_i(tmpoffset);
        }

    protected:
        // Room for GMX_SIMD_REAL_WIDTH triplets, and single real before/between/after.
        real *                  m0_;                             //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  m1_;                             //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  m2_;                             //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  m3_;                             //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  mem_;                            //!< Aligned pointer to 32*GMX_SIMD_REAL_WIDTH mem
        gmx_simd_int32_t        simdoffset_;                     //!< SIMD integer type with offsets
    private:
        real                    work_[37 * GMX_SIMD_REAL_WIDTH]; //!< Raw work memory
};



TEST_F(SimdFloatingpointUtilTest, gmxSimdExpandScalarsToTriplets)
{
    gmx_simd_real_t vs, v0, v1, v2;
    int             i;
    // We borrow the existing aligned memory to store data
    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        m0_[i] = i;
    }

    vs = gmx_simd_load_r(m0_);

    gmx_simd_expand_scalars_to_triplets_r(vs, &v0, &v1, &v2);

    gmx_simd_store_r(m0_, v0);
    gmx_simd_store_r(m1_, v1);
    gmx_simd_store_r(m2_, v2);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_EQ(i / 3, m0_[i]);
        EXPECT_EQ((i + GMX_SIMD_REAL_WIDTH) / 3, m1_[i]);
        EXPECT_EQ((i + 2 * GMX_SIMD_REAL_WIDTH) / 3, m2_[i]);
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdLoadPairs)
{
    int                               i;
    int                               offset[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t                   v0, v1;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    EXPECT_TRUE(gmx_simd_pairs_storage_alignment_r <= 32 * GMX_SIMD_REAL_WIDTH);

    // We cannot use the pre-constructed storage since the alignment is
    // specified by the SIMD implementation, and it might not be 4.
    // Borrow the existing aligned memory to store data. Pair i is (i, i+0.1)
    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        mem_[gmx_simd_pairs_storage_alignment_r * i]     = i;
        mem_[gmx_simd_pairs_storage_alignment_r * i + 1] = i + 0.1;
        offset[i] = i;
    }

    gmx_simd_load_pairs_transpose_r(mem_, offset, &v0, &v1);

    gmx_simd_store_r(m0_, v0);
    gmx_simd_store_r(m1_, v1);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_REAL_EQ_TOL(i, m0_[i], tolerance);
        EXPECT_REAL_EQ_TOL(i + 0.1, m1_[i], tolerance);
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdLoad4Transpose)
{
    gmx_simd_real_t v0, v1, v2, v3;
    gmx_simd_real_t ref0, ref1, ref2, ref3;

    ref0 = gmx_simd_load_r(m0_);
    ref1 = gmx_simd_load_r(m1_);
    ref2 = gmx_simd_load_r(m2_);
    ref3 = gmx_simd_load_r(m3_);

    gmx_simd_load_4_transpose_r(mem_, simdoffset_, &v0, &v1, &v2, &v3);

    GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
    GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    GMX_EXPECT_SIMD_REAL_EQ(ref2, v2);
    GMX_EXPECT_SIMD_REAL_EQ(ref3, v3);
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdLoad2of4Transpose)
{
    gmx_simd_real_t v0, v1;
    gmx_simd_real_t ref0, ref1;

    ref0 = gmx_simd_load_r(m0_);
    ref1 = gmx_simd_load_r(m1_);

    gmx_simd_load_2of4_transpose_r(mem_, simdoffset_, &v0, &v1);

    GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
    GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdReduceIncr4Sum)
{
    int                               i;
    gmx_simd_real_t                   v0, v1, v2, v3;
    real                              sum0, sum1, sum2, sum3, tstsum;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    v0 = gmx_simd_load_r(m0_);
    v1 = gmx_simd_load_r(m1_);
    v2 = gmx_simd_load_r(m2_);
    v3 = gmx_simd_load_r(m3_);

    sum0 = sum1 = sum2 = sum3 = 0;
    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        sum0 += m0_[i];
        sum1 += m1_[i];
        sum2 += m2_[i];
        sum3 += m3_[i];
    }

    // Just put some numbers in memory so we check the addition is correct
    m0_[0] = 5.0;
    m0_[1] = 15.0;
    m0_[2] = 25.0;
    m0_[3] = 35.0;

    tstsum = gmx_simd_reduce_incr_4_return_sum_r(m0_, v0, v1, v2, v3);

    EXPECT_REAL_EQ_TOL( 5.0+sum0, m0_[0], tolerance);
    EXPECT_REAL_EQ_TOL(15.0+sum1, m0_[1], tolerance);
    EXPECT_REAL_EQ_TOL(25.0+sum2, m0_[2], tolerance);
    EXPECT_REAL_EQ_TOL(35.0+sum3, m0_[3], tolerance);

    EXPECT_REAL_EQ_TOL(sum0 + sum1 + sum2 + sum3, tstsum, tolerance);
}


#ifdef GMX_SIMD_HAVE_HSIMD_UTIL

TEST_F(SimdFloatingpointUtilTest, gmxSimdLoadDualHsimd)
{
    gmx_simd_real_t v0, v1;

    // Point m1_ to the upper half of m0_
    m1_ = m0_ + GMX_SIMD_REAL_WIDTH/2;

    v0 = gmx_simd_load_r(m0_);
    v1 = gmx_simd_load_dual_hsimd_r(m0_, m1_);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdLoadDupHsimd)
{
    gmx_simd_real_t v0, v1;
    int             i;
    // Point m1_ to the upper half of m0_
    m1_ = m0_ + GMX_SIMD_REAL_WIDTH/2;
    // Copy data so upper half is identical to lower
    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        m1_[i] = m0_[i];
    }

    v0 = gmx_simd_load_r(m0_);
    v1 = gmx_simd_loaddup_hsimd_r(m0_);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdLoad1DualHsimd)
{
    gmx_simd_real_t v0, v1;
    int             i;
    real            data[2] = { 1, 2 };

    // Point m1_ to the upper half of m0_
    m1_ = m0_ + GMX_SIMD_REAL_WIDTH/2;
    // Set all low elements to data[0], an high to data[1]
    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        m0_[i] = data[0];
        m1_[i] = data[1];
    }

    v0 = gmx_simd_load_r(m0_);
    v1 = gmx_simd_load1_dual_hsimd_r(data);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdStoreDualHsimd)
{
    gmx_simd_real_t v0;
    int             i;

    // Point m1_ to the upper half of m0_
    m1_ = m0_ + GMX_SIMD_REAL_WIDTH/2;

    v0 = gmx_simd_load_r(m2_);
    gmx_simd_store_dual_hsimd_r(m0_, m1_, v0);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_EQ(m2_[i], m0_[i]);
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdDecrHsimd)
{
    gmx_simd_real_t                   v0;
    real                              ref[GMX_SIMD_REAL_WIDTH/2];
    int                               i;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    // Point m2_ to the upper half of m1_
    m2_ = m1_ + GMX_SIMD_REAL_WIDTH/2;
    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        ref[i] = m0_[i] - (m1_[i] + m2_[i]);
    }

    v0 = gmx_simd_load_r(m1_);
    gmx_simd_decr_hsimd_r(m0_, v0);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        EXPECT_REAL_EQ_TOL(ref[i], m0_[i], tolerance);
    }
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdLoadPairsHsimd)
{
    int                               i;
    int                               offset[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t                   v0, v1;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    EXPECT_TRUE(gmx_simd_pairs_storage_alignment_r <= 32 * GMX_SIMD_REAL_WIDTH);

    // We cannot use the pre-constructed storage since the alignment is
    // specified by the SIMD implementation, and it might not be 4.
    // Borrow the existing aligned memory to store data. Pair i is (i, i+0.1)
    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        mem_[gmx_simd_pairs_storage_alignment_r * i]     = i;
        mem_[gmx_simd_pairs_storage_alignment_r * i + 1] = i + 0.1;
        offset[i] = i;
    }

    // Point m0_ to the second half of the pairs, with the same offsets
    m0_ = mem_ + GMX_SIMD_REAL_WIDTH/2*gmx_simd_pairs_storage_alignment_r;

    gmx_simd_load_pairs_transpose_hsimd_r(mem_, m0_, offset, &v0, &v1);

    gmx_simd_store_r(m0_, v0);
    gmx_simd_store_r(m1_, v1);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_REAL_EQ_TOL(i, m0_[i], tolerance);
        EXPECT_REAL_EQ_TOL(i + 0.1, m1_[i], tolerance);
    }
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdReduceIncr4SumHsimd)
{
    int                               i;
    gmx_simd_real_t                   v0, v1;
    real                              sum0, sum1, sum2, sum3, tstsum;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    // Use the half-SIMD storage in memory m0_ and m1_.
    v0 = gmx_simd_load_r(m0_);
    v1 = gmx_simd_load_r(m1_);

    sum0 = sum1 = sum2 = sum3 = 0;
    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        sum0 += m0_[i];
        sum1 += m0_[GMX_SIMD_REAL_WIDTH/2 + i];
        sum2 += m1_[i];
        sum3 += m1_[GMX_SIMD_REAL_WIDTH/2 + i];
    }

    // Just put some numbers in memory so we check the addition is correct
    m0_[0] = 5.0;
    m0_[1] = 15.0;
    m0_[2] = 25.0;
    m0_[3] = 35.0;

    tstsum = gmx_simd_reduce_incr_4_return_sum_hsimd_r(m0_, v0, v1);

    EXPECT_REAL_EQ_TOL( 5.0 + sum0, m0_[0], tolerance);
    EXPECT_REAL_EQ_TOL(15.0 + sum1, m0_[1], tolerance);
    EXPECT_REAL_EQ_TOL(25.0 + sum2, m0_[2], tolerance);
    EXPECT_REAL_EQ_TOL(35.0 + sum3, m0_[3], tolerance);

    EXPECT_REAL_EQ_TOL(sum0 + sum1 + sum2 + sum3, tstsum, tolerance);
}


#endif

#endif      // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace
