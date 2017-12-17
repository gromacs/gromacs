/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2017, by the GROMACS development team, led by
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

#include <numeric>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/basedefinitions.h"

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

#if GMX_SIMD_HAVE_REAL

/*! \brief Test fixture for higher-level floating-point utility functions.
 *
 * Inherit from main SimdTest, add code to generate aligned memory and data.
 */
class SimdFloatingpointUtilTest : public SimdTest
{
    public:
        SimdFloatingpointUtilTest()
        {
            // Resize vectors to get the amount of memory we need
            integerMemory_.resize(GMX_SIMD_REAL_WIDTH);

            // The total memory we allocate corresponds to two work arrays
            // and 4 values each of GMX_SIMD_REAL_WIDTH.
            realMemory_.resize(2*s_workMemSize_+4*GMX_SIMD_REAL_WIDTH);

            offset_ = integerMemory_.data();
            val0_   = realMemory_.data();
            val1_   = val0_ + GMX_SIMD_REAL_WIDTH;
            val2_   = val1_ + GMX_SIMD_REAL_WIDTH;
            val3_   = val2_ + GMX_SIMD_REAL_WIDTH;
            mem0_   = val3_ + GMX_SIMD_REAL_WIDTH;
            mem1_   = mem0_ + s_workMemSize_;

            // Set default values for offset and variables val0_ through val3_
            // We cannot fill mem_ here since those values depend on the test.
            for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
            {
                // Use every third point to avoid a continguous access pattern
                offset_[i] = 3 * i;
                // Multiply numbers by 1+100*GMX_REAL_EPS ensures some low bits are
                // set too, so the tests make sure we read all bits correctly.
                val0_[i]   = (i      ) * (1.0 + 100*GMX_REAL_EPS);
                val1_[i]   = (i + 0.1) * (1.0 + 100*GMX_REAL_EPS);
                val2_[i]   = (i + 0.2) * (1.0 + 100*GMX_REAL_EPS);
                val3_[i]   = (i + 0.3) * (1.0 + 100*GMX_REAL_EPS);
            }
        }

    protected:
        //! \brief Size of memory work buffers
        //
        // To have a somewhat odd access pattern, we use every
        // third entry, so the largest value of offset_[i] is 3*GMX_SIMD_REAL_WIDTH.
        // Then we also allow alignments up to 16, which means the largest index in mem0_[]
        // that we might access is 16*3*GMX_SIMD_REAL_WIDTH+3.
        static const std::size_t                    s_workMemSize_ = 16*3*GMX_SIMD_REAL_WIDTH+4;

        std::vector<int, AlignedAllocator<int> >    integerMemory_; //!< Aligned integer memory
        std::vector<real, AlignedAllocator<real> >  realMemory_;    //!< Aligned real memory

        int *   offset_;                                            //!< Pointer to offset indices, aligned memory
        real *  val0_;                                              //!< Pointer to GMX_SIMD_REAL_WIDTH values, aligned
        real *  val1_;                                              //!< Pointer to GMX_SIMD_REAL_WIDTH values, aligned
        real *  val2_;                                              //!< Pointer to GMX_SIMD_REAL_WIDTH values, aligned
        real *  val3_;                                              //!< Pointer to GMX_SIMD_REAL_WIDTH values, aligned

        real *  mem0_;                                              //!< Pointer to aligned memory, s_workMemSize real values
        real *  mem1_;                                              //!< Pointer to aligned memory, s_workMemSize real values
};



TEST_F(SimdFloatingpointUtilTest, gatherLoadTranspose4)
{
    SimdReal        v0, v1, v2, v3;
    SimdReal        ref0, ref1, ref2, ref3;
    const int       nalign                 = 3;
    int             alignmentList[nalign]  = { 4, 8, 12 };
    int             i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
            mem0_[align * offset_[j] + 2] = val2_[j];
            mem0_[align * offset_[j] + 3] = val3_[j];
        }

        ref0 = load<SimdReal>(val0_);
        ref1 = load<SimdReal>(val1_);
        ref2 = load<SimdReal>(val2_);
        ref3 = load<SimdReal>(val3_);

        if (align == 4)
        {
            gatherLoadTranspose<4>(mem0_, offset_, &v0, &v1, &v2, &v3);
        }
        else if (align == 8)
        {
            gatherLoadTranspose<8>(mem0_, offset_, &v0, &v1, &v2, &v3);
        }
        else if (align == 12)
        {
            gatherLoadTranspose<12>(mem0_, offset_, &v0, &v1, &v2, &v3);
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

TEST_F(SimdFloatingpointUtilTest, gatherLoadTranspose2)
{
    SimdReal        v0, v1;
    SimdReal        ref0, ref1;
    const int       nalign                 = 3;
    int             alignmentList[nalign]  = { 2, 4, c_simdBestPairAlignment };
    int             i, j, align;

    EXPECT_TRUE(c_simdBestPairAlignment <= GMX_SIMD_REAL_WIDTH);

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
        }

        ref0 = load<SimdReal>(val0_);
        ref1 = load<SimdReal>(val1_);

        if (align == 2)
        {
            gatherLoadTranspose<2>(mem0_, offset_, &v0, &v1);
        }
        else if (align == 4)
        {
            gatherLoadTranspose<4>(mem0_, offset_, &v0, &v1);
        }
        else if (align == c_simdBestPairAlignment)
        {
            gatherLoadTranspose<c_simdBestPairAlignment>(mem0_, offset_, &v0, &v1);
        }
        else
        {
            FAIL();
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    }
}

TEST_F(SimdFloatingpointUtilTest, gatherLoadUTranspose3)
{
    SimdReal        v0, v1, v2;
    SimdReal        ref0, ref1, ref2;
    const int       nalign                 = 2;
    int             alignmentList[nalign]  = { 3, 4 };
    int             i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
            mem0_[align * offset_[j] + 2] = val2_[j];
        }

        ref0 = load<SimdReal>(val0_);
        ref1 = load<SimdReal>(val1_);
        ref2 = load<SimdReal>(val2_);

        if (align == 3)
        {
            gatherLoadUTranspose<3>(mem0_, offset_, &v0, &v1, &v2);
        }
        else if (align == 4)
        {
            gatherLoadUTranspose<4>(mem0_, offset_, &v0, &v1, &v2);
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

TEST_F(SimdFloatingpointUtilTest, transposeScatterStoreU3)
{
    SimdReal                          v0, v1, v2;
    real                              refmem[s_workMemSize_];
    const int                         nalign                 = 2;
    int                               alignmentList[nalign]  = { 3, 4 };
    int                               i, align;
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];

        // Set test and reference memory to background value
        for (std::size_t j = 0; j < s_workMemSize_; j++)
        {
            // Multiply by 1+100*eps to make sure low bits are also used
            mem0_[j] = refmem[j] = (1000.0 + j) * (1.0 + 100*GMX_REAL_EPS);
        }

        for (std::size_t j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            // set values in _reference_ memory (we will then test with mem0_, and compare)
            refmem[align * offset_[j]    ] = val0_[j];
            refmem[align * offset_[j] + 1] = val1_[j];
            refmem[align * offset_[j] + 2] = val2_[j];
        }

        v0 = load<SimdReal>(val0_);
        v1 = load<SimdReal>(val1_);
        v2 = load<SimdReal>(val2_);

        if (align == 3)
        {
            transposeScatterStoreU<3>(mem0_, offset_, v0, v1, v2);
        }
        else if (align == 4)
        {
            transposeScatterStoreU<4>(mem0_, offset_, v0, v1, v2);
        }
        else
        {
            FAIL();
        }

        for (std::size_t j = 0; j < s_workMemSize_; j++)
        {
            EXPECT_REAL_EQ_TOL(refmem[j], mem0_[j], tolerance);
        }
    }
}

TEST_F(SimdFloatingpointUtilTest, transposeScatterIncrU3)
{
    SimdReal                          v0, v1, v2;
    real                              refmem[s_workMemSize_];
    const int                         nalign                 = 2;
    int                               alignmentList[nalign]  = { 3, 4 };
    int                               i, align;
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];

        // Set test and reference memory to background value
        for (std::size_t j = 0; j < s_workMemSize_; j++)
        {
            // Multiply by 1+100*eps to make sure low bits are also used
            mem0_[j] = refmem[j] = (1000.0 + j) * (1.0 + 100*GMX_REAL_EPS);
        }

        for (std::size_t j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            // Add values to _reference_ memory (we will then test with mem0_, and compare)
            refmem[align * offset_[j]    ] += val0_[j];
            refmem[align * offset_[j] + 1] += val1_[j];
            refmem[align * offset_[j] + 2] += val2_[j];
        }

        v0 = load<SimdReal>(val0_);
        v1 = load<SimdReal>(val1_);
        v2 = load<SimdReal>(val2_);

        if (align == 3)
        {
            transposeScatterIncrU<3>(mem0_, offset_, v0, v1, v2);
        }
        else if (align == 4)
        {
            transposeScatterIncrU<4>(mem0_, offset_, v0, v1, v2);
        }
        else
        {
            FAIL();
        }

        for (std::size_t j = 0; j < s_workMemSize_; j++)
        {
            EXPECT_REAL_EQ_TOL(refmem[j], mem0_[j], tolerance);
        }
    }
}

TEST_F(SimdFloatingpointUtilTest, transposeScatterIncrU3Overlapping)
{
    SimdReal                          v0, v1, v2;
    real                              refmem[s_workMemSize_];
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    // Alter offset_ to make all entries point to the same (first) value, so all entries will overlap
    for (std::size_t j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
    {
        offset_[j] = 0;
    }

    // Set test and reference memory to background value
    for (std::size_t j = 0; j < s_workMemSize_; j++)
    {
        // Multiply by 1+100*eps to make sure low bits are also used
        mem0_[j] = refmem[j] = (1000.0 + j) * (1.0 + 100*GMX_REAL_EPS);
    }

    for (std::size_t j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
    {
        // Add values to _reference_ memory (we will then test with mem0_, and compare)
        refmem[3 * offset_[j]    ] += val0_[j];
        refmem[3 * offset_[j] + 1] += val1_[j];
        refmem[3 * offset_[j] + 2] += val2_[j];
    }

    v0 = load<SimdReal>(val0_);
    v1 = load<SimdReal>(val1_);
    v2 = load<SimdReal>(val2_);

    transposeScatterIncrU<3>(mem0_, offset_, v0, v1, v2);

    for (std::size_t j = 0; j < s_workMemSize_; j++)
    {
        EXPECT_REAL_EQ_TOL(refmem[j], mem0_[j], tolerance);
    }
}

TEST_F(SimdFloatingpointUtilTest, transposeScatterDecrU3)
{
    SimdReal                          v0, v1, v2;
    real                              refmem[s_workMemSize_];
    const int                         nalign                 = 2;
    int                               alignmentList[nalign]  = { 3, 4 };
    int                               i, align;
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];

        // Set test and reference memory to background value
        for (std::size_t j = 0; j < s_workMemSize_; j++)
        {
            // Multiply by 1+100*eps to make sure low bits are also used
            mem0_[j] = refmem[j] = (1000.0 + j) * (1.0 + 100*GMX_REAL_EPS);
        }

        for (std::size_t j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            // Subtract values from _reference_ memory (we will then test with mem0_, and compare)
            refmem[align * offset_[j]    ] -= val0_[j];
            refmem[align * offset_[j] + 1] -= val1_[j];
            refmem[align * offset_[j] + 2] -= val2_[j];
        }

        v0 = load<SimdReal>(val0_);
        v1 = load<SimdReal>(val1_);
        v2 = load<SimdReal>(val2_);

        if (align == 3)
        {
            transposeScatterDecrU<3>(mem0_, offset_, v0, v1, v2);
        }
        else if (align == 4)
        {
            transposeScatterDecrU<4>(mem0_, offset_, v0, v1, v2);
        }
        else
        {
            FAIL();
        }

        for (std::size_t j = 0; j < s_workMemSize_; j++)
        {
            EXPECT_REAL_EQ_TOL(refmem[j], mem0_[j], tolerance);
        }
    }
}

TEST_F(SimdFloatingpointUtilTest, transposeScatterDecrU3Overlapping)
{
    SimdReal                          v0, v1, v2;
    real                              refmem[s_workMemSize_];
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    // Alter offset_ to make all entries point to the same (first) value, so all entries will overlap
    for (std::size_t j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
    {
        offset_[j] = 0;
    }

    // Set test and reference memory to background value
    for (std::size_t j = 0; j < s_workMemSize_; j++)
    {
        // Multiply by 1+100*eps to make sure low bits are also used
        mem0_[j] = refmem[j] = (1000.0 + j) * (1.0 + 100*GMX_REAL_EPS);
    }

    for (std::size_t j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
    {
        // Subtract values from _reference_ memory (we will then test with mem0_, and compare)
        refmem[3 * offset_[j]    ] -= val0_[j];
        refmem[3 * offset_[j] + 1] -= val1_[j];
        refmem[3 * offset_[j] + 2] -= val2_[j];
    }

    v0 = load<SimdReal>(val0_);
    v1 = load<SimdReal>(val1_);
    v2 = load<SimdReal>(val2_);

    transposeScatterDecrU<3>(mem0_, offset_, v0, v1, v2);

    for (std::size_t j = 0; j < s_workMemSize_; j++)
    {
        EXPECT_REAL_EQ_TOL(refmem[j], mem0_[j], tolerance);
    }
}

TEST_F(SimdFloatingpointUtilTest, expandScalarsToTriplets)
{
    SimdReal        vs, v0, v1, v2;
    int             i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        mem0_[i] = i;
    }

    vs = load<SimdReal>(mem0_);

    expandScalarsToTriplets(vs, &v0, &v1, &v2);

    store(val0_, v0);
    store(val1_, v1);
    store(val2_, v2);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_EQ(i / 3, val0_[i]);
        EXPECT_EQ((i + GMX_SIMD_REAL_WIDTH) / 3, val1_[i]);
        EXPECT_EQ((i + 2 * GMX_SIMD_REAL_WIDTH) / 3, val2_[i]);
    }
}


TEST_F(SimdFloatingpointUtilTest, gatherLoadBySimdIntTranspose4)
{
    SimdReal         v0, v1, v2, v3;
    SimdReal         ref0, ref1, ref2, ref3;
    SimdInt32        simdoffset;
    const int        nalign                 = 3;
    int              alignmentList[nalign]  = { 4, 8, 12 };
    int              i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
            mem0_[align * offset_[j] + 2] = val2_[j];
            mem0_[align * offset_[j] + 3] = val3_[j];
        }

        simdoffset = load<SimdInt32>(offset_);
        ref0       = load<SimdReal>(val0_);
        ref1       = load<SimdReal>(val1_);
        ref2       = load<SimdReal>(val2_);
        ref3       = load<SimdReal>(val3_);

        if (align == 4)
        {
            gatherLoadBySimdIntTranspose<4>(mem0_, simdoffset, &v0, &v1, &v2, &v3);
        }
        else if (align == 8)
        {
            gatherLoadBySimdIntTranspose<8>(mem0_, simdoffset, &v0, &v1, &v2, &v3);
        }
        else if (align == 12)
        {
            gatherLoadBySimdIntTranspose<12>(mem0_, simdoffset, &v0, &v1, &v2, &v3);
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


TEST_F(SimdFloatingpointUtilTest, gatherLoadBySimdIntTranspose2)
{
    SimdReal         v0, v1;
    SimdReal         ref0, ref1;
    SimdInt32        simdoffset;
    const int        nalign                 = 3;
    int              alignmentList[nalign]  = { 4, 8, 12 };
    int              i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
        }

        simdoffset = load<SimdInt32>(offset_);
        ref0       = load<SimdReal>(val0_);
        ref1       = load<SimdReal>(val1_);

        if (align == 4)
        {
            gatherLoadBySimdIntTranspose<4>(mem0_, simdoffset, &v0, &v1);
        }
        else if (align == 8)
        {
            gatherLoadBySimdIntTranspose<8>(mem0_, simdoffset, &v0, &v1);
        }
        else if (align == 12)
        {
            gatherLoadBySimdIntTranspose<12>(mem0_, simdoffset, &v0, &v1);
        }
        else
        {
            FAIL();
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    }
}

#if GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_REAL
TEST_F(SimdFloatingpointUtilTest, gatherLoadUBySimdIntTranspose2)
{
    SimdReal         v0, v1;
    SimdReal         ref0, ref1;
    SimdInt32        simdoffset;
    const int        nalign                 = 3;
    int              alignmentList[nalign]  = { 1, 3, 5 };
    int              i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH; j++)
        {
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
        }

        simdoffset = load<SimdInt32>(offset_);
        ref0       = load<SimdReal>(val0_);
        ref1       = load<SimdReal>(val1_);

        if (align == 1)
        {
            gatherLoadUBySimdIntTranspose<1>(mem0_, simdoffset, &v0, &v1);
        }
        else if (align == 3)
        {
            gatherLoadUBySimdIntTranspose<3>(mem0_, simdoffset, &v0, &v1);
        }
        else if (align == 5)
        {
            gatherLoadUBySimdIntTranspose<5>(mem0_, simdoffset, &v0, &v1);
        }
        else
        {
            FAIL();
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    }
}
#endif      // GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_REAL

TEST_F(SimdFloatingpointUtilTest, reduceIncr4Sum)
{
    int                               i;
    SimdReal                          v0, v1, v2, v3;
    real                              sum0, sum1, sum2, sum3, tstsum;
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    v0 = load<SimdReal>(val0_);
    v1 = load<SimdReal>(val1_);
    v2 = load<SimdReal>(val2_);
    v3 = load<SimdReal>(val3_);

    sum0 = sum1 = sum2 = sum3 = 0;
    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        sum0 += val0_[i];
        sum1 += val1_[i];
        sum2 += val2_[i];
        sum3 += val3_[i];
    }

    // Just put some numbers in memory so we check the addition is correct
    mem0_[0] = c0;
    mem0_[1] = c1;
    mem0_[2] = c2;
    mem0_[3] = c3;

    tstsum = reduceIncr4ReturnSum(mem0_, v0, v1, v2, v3);

    EXPECT_REAL_EQ_TOL(c0 + sum0, mem0_[0], tolerance);
    EXPECT_REAL_EQ_TOL(c1 + sum1, mem0_[1], tolerance);
    EXPECT_REAL_EQ_TOL(c2 + sum2, mem0_[2], tolerance);
    EXPECT_REAL_EQ_TOL(c3 + sum3, mem0_[3], tolerance);

    EXPECT_REAL_EQ_TOL(sum0 + sum1 + sum2 + sum3, tstsum, tolerance);
}

#if GMX_SIMD_HAVE_HSIMD_UTIL_REAL

TEST_F(SimdFloatingpointUtilTest, loadDualHsimd)
{
    SimdReal v0, v1;

    // Point p to the upper half of val0_
    real * p = val0_ + GMX_SIMD_REAL_WIDTH / 2;

    v0 = load<SimdReal>(val0_);
    v1 = loadDualHsimd(val0_, p);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}

TEST_F(SimdFloatingpointUtilTest, loadDuplicateHsimd)
{
    SimdReal        v0, v1;
    int             i;
    // Point p to the upper half of val0_
    real          * p = val0_ + GMX_SIMD_REAL_WIDTH / 2;
    // Copy data so upper half is identical to lower
    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        p[i] = val0_[i];
    }

    v0 = load<SimdReal>(val0_);
    v1 = loadDuplicateHsimd(val0_);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}


TEST_F(SimdFloatingpointUtilTest, loadU1DualHsimd)
{
    SimdReal        v0, v1;
    int             i;
    real            data[2] = { 1, 2 };

    // Point p to the upper half of val0_
    real * p = val0_ + GMX_SIMD_REAL_WIDTH / 2;
    // Set all low elements to data[0], an high to data[1]
    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        val0_[i] = data[0];
        p[i]     = data[1];
    }

    v0 = load<SimdReal>(val0_);
    v1 = loadU1DualHsimd(data);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}


TEST_F(SimdFloatingpointUtilTest, storeDualHsimd)
{
    SimdReal        v0;
    int             i;

    // Point p to the upper half of val0_
    real * p = val0_ + GMX_SIMD_REAL_WIDTH / 2;

    v0 = load<SimdReal>(val2_);
    storeDualHsimd(val0_, p, v0);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_EQ(val2_[i], val0_[i]);
    }
}

TEST_F(SimdFloatingpointUtilTest, incrDualHsimd)
{
    real            reference[GMX_SIMD_REAL_WIDTH];
    SimdReal        v0;

    // Create reference values
    for (std::size_t i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        reference[i] = val0_[i] + val2_[i];
    }

    // Point p to the upper half of val0_
    real * p = val0_ + GMX_SIMD_REAL_WIDTH / 2;

    v0 = load<SimdReal>(val2_);
    incrDualHsimd(val0_, p, v0);

    for (std::size_t i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        EXPECT_EQ(reference[i], val0_[i]);
    }
}

TEST_F(SimdFloatingpointUtilTest, incrDualHsimdOverlapping)
{
    real            reference[GMX_SIMD_REAL_WIDTH/2];
    SimdReal        v0;

    // Create reference values
    for (std::size_t i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        reference[i] = val0_[i] + val2_[i] + val2_[GMX_SIMD_REAL_WIDTH/2+i];
    }

    v0 = load<SimdReal>(val2_);
    incrDualHsimd(val0_, val0_, v0);

    for (std::size_t i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        EXPECT_EQ(reference[i], val0_[i]);
    }
}

TEST_F(SimdFloatingpointUtilTest, decrHsimd)
{
    SimdReal                          v0;
    real                              ref[GMX_SIMD_REAL_WIDTH / 2];
    int                               i;
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    // Point p to the upper half of val1_
    real * p = val1_ + GMX_SIMD_REAL_WIDTH / 2;
    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        ref[i] = val0_[i] - ( val1_[i] + p[i] );
    }

    v0 = load<SimdReal>(val1_);
    decrHsimd(val0_, v0);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        EXPECT_REAL_EQ_TOL(ref[i], val0_[i], tolerance);
    }
}


TEST_F(SimdFloatingpointUtilTest, gatherLoadTranspose2Hsimd)
{
    SimdReal        v0, v1;
    SimdReal        ref0, ref1;

    const int       nalign                 = 3;
    int             alignmentList[nalign]  = { 2, 4, c_simdBestPairAlignment };
    int             i, j, align;

    for (i = 0; i < nalign; i++)
    {
        align = alignmentList[i];
        for (j = 0; j < GMX_SIMD_REAL_WIDTH / 2; j++)
        {
            // Use mem0_ as base for lower half
            mem0_[align * offset_[j]    ] = val0_[j];
            mem0_[align * offset_[j] + 1] = val1_[j];
            // Use mem1_ as base for upper half
            mem1_[align * offset_[j]    ] = val0_[GMX_SIMD_REAL_WIDTH / 2 + j];
            mem1_[align * offset_[j] + 1] = val1_[GMX_SIMD_REAL_WIDTH / 2 + j];

        }

        ref0 = load<SimdReal>(val0_);
        ref1 = load<SimdReal>(val1_);

        if (align == 2)
        {
            gatherLoadTransposeHsimd<2>(mem0_, mem1_, offset_, &v0, &v1);
        }
        else if (align == 4)
        {
            gatherLoadTransposeHsimd<4>(mem0_, mem1_, offset_, &v0, &v1);
        }
        else if (align == c_simdBestPairAlignment)
        {
            gatherLoadTransposeHsimd<c_simdBestPairAlignment>(mem0_, mem1_, offset_, &v0, &v1);
        }
        else
        {
            FAIL();
        }

        GMX_EXPECT_SIMD_REAL_EQ(ref0, v0);
        GMX_EXPECT_SIMD_REAL_EQ(ref1, v1);
    }
}


TEST_F(SimdFloatingpointUtilTest, reduceIncr4SumHsimd)
{
    int                               i;
    SimdReal                          v0, v1;
    real                              sum0, sum1, sum2, sum3, tstsum;
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    // Use the half-SIMD storage in memory val0_ and val1_.
    v0 = load<SimdReal>(val0_);
    v1 = load<SimdReal>(val1_);

    sum0 = sum1 = sum2 = sum3 = 0;
    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 2; i++)
    {
        sum0 += val0_[i];
        sum1 += val0_[GMX_SIMD_REAL_WIDTH / 2 + i];
        sum2 += val1_[i];
        sum3 += val1_[GMX_SIMD_REAL_WIDTH / 2 + i];
    }

    // Just put some numbers in memory so we check the addition is correct
    mem0_[0] = c0;
    mem0_[1] = c1;
    mem0_[2] = c2;
    mem0_[3] = c3;

    tstsum = reduceIncr4ReturnSumHsimd(mem0_, v0, v1);

    EXPECT_REAL_EQ_TOL(c0 + sum0, mem0_[0], tolerance);
    EXPECT_REAL_EQ_TOL(c1 + sum1, mem0_[1], tolerance);
    EXPECT_REAL_EQ_TOL(c2 + sum2, mem0_[2], tolerance);
    EXPECT_REAL_EQ_TOL(c3 + sum3, mem0_[3], tolerance);

    EXPECT_REAL_EQ_TOL(sum0 + sum1 + sum2 + sum3, tstsum, tolerance);
}

#endif      // GMX_SIMD_HAVE_HSIMD_UTIL_REAL

//Test Currently doesn't work for GMX_SIMD_REAL_WIDTH<4. Should be fixed by having GMX_EXPECT_SIMD_REAL_EQ which works for both Simd and Simd4
#if GMX_SIMD_HAVE_4NSIMD_UTIL_REAL && GMX_SIMD_REAL_WIDTH >= 4

TEST_F(SimdFloatingpointUtilTest, loadUNDuplicate4)
{
    Simd4NReal      v0, v1;
    int             i;
    real            data[GMX_SIMD_REAL_WIDTH/4];
    std::iota(data, data+GMX_SIMD_REAL_WIDTH/4, 1);

#if defined _ICC && __ICC == 1800 || defined __ICL && __ICL == 1800
#pragma novector /* Work-around for incorrect vectorization for AVX_512(_KNL) */
#endif
    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 4; i++)
    {
        val0_[i*4] = val0_[i*4+1] = val0_[i*4+2] = val0_[i*4+3] = data[i];
    }

    v0 = load<Simd4NReal>(val0_);
    v1 = loadUNDuplicate4(data);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}

TEST_F(SimdFloatingpointUtilTest, load4DuplicateN)
{
    Simd4NReal        v0, v1;
    int               i;
    real              data[4] = { 1, 2, 3, 4};

    for (i = 0; i < GMX_SIMD_REAL_WIDTH / 4; i++)
    {
        val0_[i*4]   = data[0];
        val0_[i*4+1] = data[1];
        val0_[i*4+2] = data[2];
        val0_[i*4+3] = data[3];
    }

    v0 = load<Simd4NReal>(val0_);
    v1 = load4DuplicateN(val0_);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}

TEST_F(SimdFloatingpointUtilTest, loadU4NOffset)
{
    constexpr int   offset  = 6; //non power of 2
    constexpr int   dataLen = 4+offset*(GMX_SIMD_REAL_WIDTH/4-1);
    real            data[dataLen];
    std::iota(data, data+dataLen, 1);

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH / 4; i++)
    {
        val0_[i*4]   = data[0+offset*i];
        val0_[i*4+1] = data[1+offset*i];
        val0_[i*4+2] = data[2+offset*i];
        val0_[i*4+3] = data[3+offset*i];
    }

    const Simd4NReal v0 = load<Simd4NReal>(val0_);
    const Simd4NReal v1 = loadU4NOffset(data, offset);

    GMX_EXPECT_SIMD_REAL_EQ(v0, v1);
}

#endif      // GMX_SIMD_HAVE_4NSIMD_UTIL_REAL

#endif      // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace
