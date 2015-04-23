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

#if defined(GMX_SIMD_HAVE_REAL) && defined(GMX_SIMD_V2)

/*! \brief Test fixture for higher-level floating-point utility functions.
 *
 * Inherit from main SimdTest, add code to generate aligned memory and data.
 */
class SimdFloatingpointUtilTest : public SimdTest
{
    public:
        SimdFloatingpointUtilTest()
        {
            int i;
            // Create aligned pointers
            offset_ = gmx_simd_align_i(iwork_);
            val0_   = gmx_simd_align_r(work_);
            val1_   = val0_ +     GMX_SIMD_REAL_WIDTH;
            val2_   = val0_ + 2 * GMX_SIMD_REAL_WIDTH;
            val3_   = val0_ + 3 * GMX_SIMD_REAL_WIDTH;
            mem0_   = val0_ + 4 * GMX_SIMD_REAL_WIDTH;
            mem1_   = mem0_ + 8 * GMX_SIMD_REAL_WIDTH;
            // Set default values for offset and variables m0_ through m3_
            // We cannot fill mem_ here since those values depend on the test.
            for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
            {
                // Use every third point to avoid a continguous access pattern
                offset_[i] = 3 * i;
                val0_[i]   = i;
                val1_[i]   = i + 0.1;
                val2_[i]   = i + 0.2;
                val3_[i]   = i + 0.3;
            }

        }

    protected:
        int *                   offset_;                         //!< Aligned pointer to offset array
        real *                  val0_;                           //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  val1_;                           //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  val2_;                           //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  val3_;                           //!< Aligned pointer to GMX_SIMD_REAL_WIDTH mem
        real *                  mem0_;                           //!< Aligned pointer to 16*GMX_SIMD_REAL_WIDTH mem
        real *                  mem1_;                           //!< Aligned pointer halfway through mem0_
    private:
        // To have a somewhat odd access pattern, we use every
        // third entry, so the largest value of offset[i] is 3*GMX_SIMD_REAL_WIDTH.
        // Then we also allow alignments up to 16, which means the largest index in mem0_[]
        // that we might access is 16*3*GMX_SIMD_REAL_WIDTH+3.
        // On top of this we also have the four m0-m3 locations, and the alignment padding.
        // This makes the work size 53*GMX_SIMD_REAL_WIDTH+3.
        real                    work_[53 * GMX_SIMD_REAL_WIDTH + 3]; //!< Raw work memory (64+4+pad=69)
        int                     iwork_[2 * GMX_SIMD_REAL_WIDTH];     //!< Raw work memory for offsets (1+pad=2)
};



TEST_F(SimdFloatingpointUtilTest, gmxSimdGatherLoadTranspose4)
{
    gmx_simd_real_t v0, v1, v2, v3;
    gmx_simd_real_t ref0, ref1, ref2, ref3;
    const int       nalign                 = 3;
    int             alignment_list[nalign] = { 4, 8, 12 };
    int             i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
            mem0_[align * offset_[j] + 2] = val2_[j];
            mem0_[align * offset_[j] + 3] = val3_[j];
        }

        ref0 = gmx_simd_load_r(val0_);
        ref1 = gmx_simd_load_r(val1_);
        ref2 = gmx_simd_load_r(val2_);
        ref3 = gmx_simd_load_r(val3_);

        if (align == 4)
        {
            gmx_simd_gather_load_transpose_r<4>(mem0_, offset_, v0, v1, v2, v3);
        }
        else if (align == 8)
        {
            gmx_simd_gather_load_transpose_r<8>(mem0_, offset_, v0, v1, v2, v3);
        }
        else if (align == 12)
        {
            gmx_simd_gather_load_transpose_r<12>(mem0_, offset_, v0, v1, v2, v3);
        }
        else
        {
            FAIL();
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
        GMX_EXPECT_SIMD_REAL_EQ(ref2, v2);
        GMX_EXPECT_SIMD_REAL_EQ(ref3, v3);
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdGatherLoadTranspose2)
{
    gmx_simd_real_t v0, v1;
    gmx_simd_real_t ref0, ref1;
    const int       nalign                 = 3;
    int             alignment_list[nalign] = { 2, 4, gmx_simd_best_pair_alignment_r };
    int             i, j, align;

    EXPECT_TRUE(gmx_simd_best_pair_alignment_r <= GMX_SIMD_REAL_WIDTH);

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
        }

        ref0 = gmx_simd_load_r(val0_);
        ref1 = gmx_simd_load_r(val1_);

        if (align == 2)
        {
            gmx_simd_gather_load_transpose_r<2>(mem0_, offset_, v0, v1);
        }
        else if (align == 4)
        {
            gmx_simd_gather_load_transpose_r<4>(mem0_, offset_, v0, v1);
        }
        else if (align == gmx_simd_best_pair_alignment_r)
        {
            gmx_simd_gather_load_transpose_r<gmx_simd_best_pair_alignment_r>(mem0_, offset_, v0, v1);
        }
        else
        {
            FAIL();
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdGatherLoadUTranspose3)
{
    gmx_simd_real_t v0, v1, v2;
    gmx_simd_real_t ref0, ref1, ref2;
    const int       nalign                 = 2;
    int             alignment_list[nalign] = { 3, 4 };
    int             i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
            mem0_[align * offset_[j] + 2] = val2_[j];
        }

        ref0 = gmx_simd_load_r(val0_);
        ref1 = gmx_simd_load_r(val1_);
        ref2 = gmx_simd_load_r(val2_);

        if (align == 3)
        {
            gmx_simd_gather_loadu_transpose_r<3>(mem0_, offset_, v0, v1, v2);
        }
        else if (align == 4)
        {
            gmx_simd_gather_loadu_transpose_r<4>(mem0_, offset_, v0, v1, v2);
        }
        else
        {
            FAIL();
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
        GMX_EXPECT_SIMD_REAL_EQ(ref2, v2);
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdTransposeScatterStoreU3)
{
    gmx_simd_real_t                   v0, v1, v2;
    real                              refmem[12 * GMX_SIMD_REAL_WIDTH]; // Same amount (4*3) as mem0_ in class
    const int                         nalign                 = 2;
    int                               alignment_list[nalign] = { 3, 4 };
    int                               i, j, align;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];

        // Set test and reference memory to background value
        for (j = 0; j < 12 * GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[j] = refmem[j] = 1000.0 + j;
        }

        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            // set values in _reference_ memory (we will then test with mem0_, and compare)
            refmem[align * offset_[j]    ] = val0_[j];
            refmem[align * offset_[j] + 1] = val1_[j];
            refmem[align * offset_[j] + 2] = val2_[j];
        }

        v0 = gmx_simd_load_r(val0_);
        v1 = gmx_simd_load_r(val1_);
        v2 = gmx_simd_load_r(val2_);

        if (align == 3)
        {
            gmx_simd_transpose_scatter_storeu_r<3>(mem0_, offset_, v0, v1, v2);
        }
        else if (align == 4)
        {
            gmx_simd_transpose_scatter_storeu_r<4>(mem0_, offset_, v0, v1, v2);
        }
        else
        {
            FAIL();
        }

        for (j = 0; j < 12 * GMX_SIMD_REAL_WIDTH; j++)
        {
            EXPECT_REAL_EQ_TOL(refmem[j], mem0_[j], tolerance);
        }
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdTransposeScatterIncrU3)
{
    gmx_simd_real_t                   v0, v1, v2;
    real                              refmem[12 * GMX_SIMD_REAL_WIDTH]; // Same amount (4*3) as mem0_ in class
    const int                         nalign                 = 2;
    int                               alignment_list[nalign] = { 3, 4 };
    int                               i, j, align;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];

        // Set test and reference memory to background value
        for (j = 0; j < 12 * GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[j] = refmem[j] = 1000.0 + j;
        }

        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            // Add values to _reference_ memory (we will then test with mem0_, and compare)
            refmem[align * offset_[j]    ] += val0_[j];
            refmem[align * offset_[j] + 1] += val1_[j];
            refmem[align * offset_[j] + 2] += val2_[j];
        }

        v0 = gmx_simd_load_r(val0_);
        v1 = gmx_simd_load_r(val1_);
        v2 = gmx_simd_load_r(val2_);

        if (align == 3)
        {
            gmx_simd_transpose_scatter_incru_r<3>(mem0_, offset_, v0, v1, v2);
        }
        else if (align == 4)
        {
            gmx_simd_transpose_scatter_incru_r<4>(mem0_, offset_, v0, v1, v2);
        }
        else
        {
            FAIL();
        }

        for (j = 0; j < 12 * GMX_SIMD_REAL_WIDTH; j++)
        {
            EXPECT_REAL_EQ_TOL(refmem[j], mem0_[j], tolerance);
        }
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdTransposeScatterDecrU3)
{
    gmx_simd_real_t                   v0, v1, v2;
    real                              refmem[12*GMX_SIMD_REAL_WIDTH]; // Same amount (4*3) as mem0_ in class
    const int                         nalign                 = 2;
    int                               alignment_list[nalign] = { 3, 4 };
    int                               i, j, align;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];

        // Set test and reference memory to background value
        for (j = 0; j < 12 * GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[j] = refmem[j] = 1000.0 + j;
        }

        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            // Subtract values from _reference_ memory (we will then test with mem0_, and compare)
            refmem[align * offset_[j]    ] -= val0_[j];
            refmem[align * offset_[j] + 1] -= val1_[j];
            refmem[align * offset_[j] + 2] -= val2_[j];
        }

        v0 = gmx_simd_load_r(val0_);
        v1 = gmx_simd_load_r(val1_);
        v2 = gmx_simd_load_r(val2_);

        if (align == 3)
        {
            gmx_simd_transpose_scatter_decru_r<3>(mem0_, offset_, v0, v1, v2);
        }
        else if (align == 4)
        {
            gmx_simd_transpose_scatter_decru_r<4>(mem0_, offset_, v0, v1, v2);
        }
        else
        {
            FAIL();
        }

        for (j = 0; j < 12*GMX_SIMD_REAL_WIDTH; j++)
        {
            EXPECT_REAL_EQ_TOL(refmem[j], mem0_[j], tolerance);
        }
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdExpandScalarsToTriplets)
{
    gmx_simd_real_t vs, v0, v1, v2;
    int             i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        mem0_[i] = i;
    }

    vs = gmx_simd_load_r(mem0_);

    gmx_simd_expand_scalars_to_triplets_r(vs, v0, v1, v2);

    gmx_simd_store_r(val0_, v0);
    gmx_simd_store_r(val1_, v1);
    gmx_simd_store_r(val2_, v2);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_EQ(i / 3, val0_[i]);
        EXPECT_EQ((i + GMX_SIMD_REAL_WIDTH) / 3, val1_[i]);
        EXPECT_EQ((i + 2 * GMX_SIMD_REAL_WIDTH) / 3, val2_[i]);
    }
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdGatherLoadBySimdIntTranspose4)
{
    gmx_simd_real_t  v0, v1, v2, v3;
    gmx_simd_real_t  ref0, ref1, ref2, ref3;
    gmx_simd_int32_t simdoffset;
    const int        nalign                 = 3;
    int              alignment_list[nalign] = { 4, 8, 12 };
    int              i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
            mem0_[align * offset_[j] + 2] = val2_[j];
            mem0_[align * offset_[j] + 3] = val3_[j];
        }

        simdoffset = gmx_simd_load_i(offset_);
        ref0       = gmx_simd_load_r(val0_);
        ref1       = gmx_simd_load_r(val1_);
        ref2       = gmx_simd_load_r(val2_);
        ref3       = gmx_simd_load_r(val3_);

        if (align == 4)
        {
            gmx_simd_gather_load_bysimdint_transpose_r<4>(mem0_, simdoffset, v0, v1, v2, v3);
        }
        else if (align == 8)
        {
            gmx_simd_gather_load_bysimdint_transpose_r<8>(mem0_, simdoffset, v0, v1, v2, v3);
        }
        else if (align == 12)
        {
            gmx_simd_gather_load_bysimdint_transpose_r<12>(mem0_, simdoffset, v0, v1, v2, v3);
        }
        else
        {
            FAIL();
        }


        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
        GMX_EXPECT_SIMD_REAL_EQ(ref2, v2);
        GMX_EXPECT_SIMD_REAL_EQ(ref3, v3);
    }
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdGatherLoadBySimdIntTranspose2)
{
    gmx_simd_real_t  v0, v1;
    gmx_simd_real_t  ref0, ref1;
    gmx_simd_int32_t simdoffset;
    const int        nalign                 = 3;
    int              alignment_list[nalign] = { 4, 8, 12 };
    int              i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
        }

        simdoffset = gmx_simd_load_i(offset_);
        ref0       = gmx_simd_load_r(val0_);
        ref1       = gmx_simd_load_r(val1_);

        if (align == 4)
        {
            gmx_simd_gather_load_bysimdint_transpose_r<4>(mem0_, simdoffset, v0, v1);
        }
        else if (align == 8)
        {
            gmx_simd_gather_load_bysimdint_transpose_r<8>(mem0_, simdoffset, v0, v1);
        }
        else if (align == 12)
        {
            gmx_simd_gather_load_bysimdint_transpose_r<12>(mem0_, simdoffset, v0, v1);
        }
        else
        {
            FAIL();
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    }
}

#ifdef GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_REAL
TEST_F(SimdFloatingpointUtilTest, gmxSimdGatherLoadUBySimdIntTranspose2)
{
    gmx_simd_real_t  v0, v1;
    gmx_simd_real_t  ref0, ref1;
    gmx_simd_int32_t simdoffset;
    const int        nalign                 = 3;
    int              alignment_list[nalign] = { 1, 3, 5 };
    int              i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
        }

        simdoffset = gmx_simd_load_i(offset_);
        ref0       = gmx_simd_load_r(val0_);
        ref1       = gmx_simd_load_r(val1_);

        if (align == 1)
        {
            gmx_simd_gather_loadu_bysimdint_transpose_r<1>(mem0_, simdoffset, v0, v1);
        }
        else if (align == 3)
        {
            gmx_simd_gather_loadu_bysimdint_transpose_r<3>(mem0_, simdoffset, v0, v1);
        }
        else if (align == 5)
        {
            gmx_simd_gather_loadu_bysimdint_transpose_r<5>(mem0_, simdoffset, v0, v1);
        }
        else
        {
            FAIL();
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    }
}
#endif

TEST_F(SimdFloatingpointUtilTest, gmxSimdReduceIncr4Sum)
{
    int                               i;
    gmx_simd_real_t                   v0, v1, v2, v3;
    real                              sum0, sum1, sum2, sum3, tstsum;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    v0 = gmx_simd_load_r(val0_);
    v1 = gmx_simd_load_r(val1_);
    v2 = gmx_simd_load_r(val2_);
    v3 = gmx_simd_load_r(val3_);

    sum0 = sum1 = sum2 = sum3 = 0;
    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        sum0 += val0_[i];
        sum1 += val1_[i];
        sum2 += val2_[i];
        sum3 += val3_[i];
    }

    // Just put some numbers in memory so we check the addition is correct
    mem0_[0] = 5.0;
    mem0_[1] = 15.0;
    mem0_[2] = 25.0;
    mem0_[3] = 35.0;

    tstsum = gmx_simd_reduce_incr_4_return_sum_r(mem0_, v0, v1, v2, v3);

    EXPECT_REAL_EQ_TOL( 5.0 + sum0, mem0_[0], tolerance);
    EXPECT_REAL_EQ_TOL(15.0 + sum1, mem0_[1], tolerance);
    EXPECT_REAL_EQ_TOL(25.0 + sum2, mem0_[2], tolerance);
    EXPECT_REAL_EQ_TOL(35.0 + sum3, mem0_[3], tolerance);

    EXPECT_REAL_EQ_TOL(sum0 + sum1 + sum2 + sum3, tstsum, tolerance);
}

#ifdef GMX_SIMD_HAVE_HSIMD_UTIL_REAL

TEST_F(SimdFloatingpointUtilTest, gmxSimdLoadDualHsimd)
{
    gmx_simd_real_t v0, v1;

    // Point val1_ to the upper half of val0_
    val1_ = val0_ + GMX_SIMD_REAL_WIDTH / 2;

    v0 = gmx_simd_load_r(val0_);
    v1 = gmx_simd_load_dual_hsimd_r(val0_, val1_);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdLoadDupHsimd)
{
    gmx_simd_real_t v0, v1;
    int             i;
    // Point val1_ to the upper half of val0_
    val1_ = val0_ + GMX_SIMD_REAL_WIDTH / 2;
    // Copy data so upper half is identical to lower
    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        val1_[i] = val0_[i];
    }

    v0 = gmx_simd_load_r(val0_);
    v1 = gmx_simd_loaddup_hsimd_r(val0_);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdLoad1DualHsimd)
{
    gmx_simd_real_t v0, v1;
    int             i;
    real            data[2] = { 1, 2 };

    // Point val1_ to the upper half of val0_
    val1_ = val0_ + GMX_SIMD_REAL_WIDTH / 2;
    // Set all low elements to data[0], an high to data[1]
    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        val0_[i] = data[0];
        val1_[i] = data[1];
    }

    v0 = gmx_simd_load_r(val0_);
    v1 = gmx_simd_load1_dual_hsimd_r(data);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdStoreDualHsimd)
{
    gmx_simd_real_t v0;
    int             i;

    // Point val1_ to the upper half of val0_
    val1_ = val0_ + GMX_SIMD_REAL_WIDTH / 2;

    v0 = gmx_simd_load_r(val2_);
    gmx_simd_store_dual_hsimd_r(val0_, val1_, v0);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_EQ(val2_[i], val0_[i]);
    }
}

TEST_F(SimdFloatingpointUtilTest, gmxSimdDecrHsimd)
{
    gmx_simd_real_t                   v0;
    real                              ref[GMX_SIMD_REAL_WIDTH / 2];
    int                               i;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    // Point val2_ to the upper half of val1_
    val2_ = val1_ + GMX_SIMD_REAL_WIDTH / 2;
    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        ref[i] = val0_[i] - ( val1_[i] + val2_[i] );
    }

    v0 = gmx_simd_load_r(val1_);
    gmx_simd_decr_hsimd_r(val0_, v0);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        EXPECT_REAL_EQ_TOL(ref[i], val0_[i], tolerance);
    }
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdGatherLoadTranspose2Hsimd)
{
    gmx_simd_real_t v0, v1;
    gmx_simd_real_t ref0, ref1;

    const int       nalign                 = 3;
    int             alignment_list[nalign] = { 2, 4, gmx_simd_best_pair_alignment_r };
    int             i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignment_list[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH / 2; j++)
        {
            // Use mem0_ as base for lower half
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
            // Use mem1_ as base for upper half
            mem1_[align * offset_[j]    ] = val0_[GMX_SIMD_REAL_WIDTH / 2 + j];
            mem1_[align * offset_[j] + 1] = val1_[GMX_SIMD_REAL_WIDTH / 2 + j];

        }

        ref0 = gmx_simd_load_r(val0_);
        ref1 = gmx_simd_load_r(val1_);

        if (align == 2)
        {
            gmx_simd_gather_load_transpose_hsimd_r<2>(mem0_, mem1_, offset_, v0, v1);
        }
        else if (align == 4)
        {
            gmx_simd_gather_load_transpose_hsimd_r<4>(mem0_, mem1_, offset_, v0, v1);
        }
        else if (align == gmx_simd_best_pair_alignment_r)
        {
            gmx_simd_gather_load_transpose_hsimd_r<gmx_simd_best_pair_alignment_r>(mem0_, mem1_, offset_, v0, v1);
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    }
}


TEST_F(SimdFloatingpointUtilTest, gmxSimdReduceIncr4SumHsimd)
{
    int                               i;
    gmx_simd_real_t                   v0, v1;
    real                              sum0, sum1, sum2, sum3, tstsum;
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());

    // Use the half-SIMD storage in memory val0_ and val1_.
    v0 = gmx_simd_load_r(val0_);
    v1 = gmx_simd_load_r(val1_);

    sum0 = sum1 = sum2 = sum3 = 0;
    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        sum0 += val0_[i];
        sum1 += val0_[GMX_SIMD_REAL_WIDTH / 2 + i];
        sum2 += val1_[i];
        sum3 += val1_[GMX_SIMD_REAL_WIDTH / 2 + i];
    }

    // Just put some numbers in memory so we check the addition is correct
    mem0_[0] = 5.0;
    mem0_[1] = 15.0;
    mem0_[2] = 25.0;
    mem0_[3] = 35.0;

    tstsum = gmx_simd_reduce_incr_4_return_sum_hsimd_r(mem0_, v0, v1);

    EXPECT_REAL_EQ_TOL( 5.0 + sum0, mem0_[0], tolerance);
    EXPECT_REAL_EQ_TOL(15.0 + sum1, mem0_[1], tolerance);
    EXPECT_REAL_EQ_TOL(25.0 + sum2, mem0_[2], tolerance);
    EXPECT_REAL_EQ_TOL(35.0 + sum3, mem0_[3], tolerance);

    EXPECT_REAL_EQ_TOL(sum0 + sum1 + sum2 + sum3, tstsum, tolerance);
}


#endif

#endif      // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace
