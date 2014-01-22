/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

//! \internal \file
// \brief
// Declares test fixture for SimdTests.
//
// The SIMD tests are both simple and complicated. The actual testing logic
// is _very_ straightforward since we just need to test single values against
// the math library, and for some math functions we need to do it in a loop.
// This could have been achieved in minutes with the default Google Test tools,
// if it wasn't for the problem that we cannot access or compare SIMD contents
// directly without using lots of other SIMD functionality. For this reason
// we have separate the basic testing of load/store operations into a separate
// bootstrapping test. Once this woreks, we can use the utility routines to
// convert SIMD contents to/from stl:vector<> and perform the rest of the tests.
//
//
// \author Erik Lindahl <erik.lindahl@scilifelab.se>
// \ingroup module_simd
//
#ifndef GMX_SIMD_TESTS_UTIL_H
#define GMX_SIMD_TESTS_UTIL_H

#include <vector>
#include <gtest/gtest.h>
#include "gromacs/simd/simd.h"

namespace simdTest
{

// The SIMD module is special in that it is just a (big) collection of simple
// macros that get compiled directly (typically 1-to-1) to assembly.
// This means that the testing needs to assess a large number of very simple
// macros that look like functions, although since they aren't proper functions
// we cannot rely on C++ features such as overloading, which also means the
// class below is mostly a utility container for initializing SIMD data.
//
// Another complication is that the width of the SIMD implementation will
// depend on the hardware and precision. For some simple operations it is
// sufficient to set all SIMD elements to the same value, and check that the
// result is present in all elements. However, for a few more complex
// instructions that might rely on shuffling under-the-hood it is important
// that we can test operations with different elements. We achieve this by
// having test code that can initialize a SIMD variable from an stl::vector
// of arbitrary length; the vector is simply repeated to fill all elements in
// the SIMD variable. We also have similar routines to compare a SIMD result
// with values in a vector, which returns true iff all elements match.
//
// This way we can write simple tests that use different values for all SIMD
// elements. Personally I like using vectors of length _3_, since this means
// there are no simple repeated patterns in low/high halves of SIMD variables
// that are 2,4,8,or 16 elements wide, and we still don't have to care about
// the exact SIMD width of the underlying implementation.
//
// Note that this utility uses a few SIMD load/store instructions internally -
// those have been tested separately in the bootstrap_loadstore.cpp file.
//

// Note that the SIMD constants cannot be part of the class - some x86 compilers
// have problems with the alignments of class-internal data required for AVX.
#ifdef GMX_SIMD_HAVE_REAL
extern const gmx_simd_real_t rSimd_1_2_3;
extern const gmx_simd_real_t rSimd_4_5_6;
extern const gmx_simd_real_t rSimd_7_8_9;
extern const gmx_simd_real_t rSimd_5_7_9;
extern const gmx_simd_real_t rSimd_m1_m2_m3;
extern const gmx_simd_real_t rSimd_3_1_4;
extern const gmx_simd_real_t rSimd_m3_m1_m4;
extern const gmx_simd_real_t rSimd_2p25;
extern const gmx_simd_real_t rSimd_3p75;
extern const gmx_simd_real_t rSimd_m2p25;
extern const gmx_simd_real_t rSimd_m3p75;
extern const gmx_simd_real_t rSimd_Exp;
// Magic FP numbers corresponding to specific bit patterns
extern const gmx_simd_real_t rSimd_Bits1;
extern const gmx_simd_real_t rSimd_Bits2;
extern const gmx_simd_real_t rSimd_Bits3;
extern const gmx_simd_real_t rSimd_Bits4;
extern const gmx_simd_real_t rSimd_Bits5;
extern const gmx_simd_real_t rSimd_Bits6;
// number of bits in default floating-point precision mantissa
extern const int             mantissaBits;
#endif
#ifdef GMX_SIMD_HAVE_INT32
extern const gmx_simd_int32_t iSimd_1_2_3;
extern const gmx_simd_int32_t iSimd_4_5_6;
extern const gmx_simd_int32_t iSimd_7_8_9;
extern const gmx_simd_int32_t iSimd_5_7_9;
extern const gmx_simd_int32_t iSimd_1M_2M_3M;
extern const gmx_simd_int32_t iSimd_4M_5M_6M;
extern const gmx_simd_int32_t iSimd_5M_7M_9M;
extern const gmx_simd_int32_t iSimd_0xF0F0F0F0;
extern const gmx_simd_int32_t iSimd_0xCCCCCCCC;
#endif

//! \brief Test fixture for SIMD tests - contains test settings.
//
// This is a very simple test fixture. It is only used to set default
// tolerances and ranges, mainly for floating-point and math tests, so we
// save a few operations in the test routines and can override the settings
// from the command line.
class SimdTest : public ::testing::Test
{
    public:
        SimdTest()
        {
#ifdef GMX_DOUBLE
            ulpTol_ = 255LL; //!< Aim for ~twice the precision we have in single
#else
            ulpTol_ = 10LL;  //!< Be a bit liberal so compiler optimization doesn't bite us
#endif
            absTol_ = 0;
        }
        void setUlpTol(gmx_int64_t newTol)   { ulpTol_ = newTol; }
        void setAbsTol(real newTol)          { absTol_ = newTol; }

#ifdef GMX_SIMD_HAVE_REAL
        //! \brief Compare two SIMD real.
        //
        // The comparison is applied to each element, and it returns true if all values
        // in the SIMD are within the class tolerances of the corresponding vector values.
        ::testing::AssertionResult
        compareSimdRealUlp(const char * refExpr, const char * tstExpr, gmx_simd_real_t ref, gmx_simd_real_t tst);

        //! \brief Compare scalar against SIMD real.
        //
        // Convenience routine: expands the scalar into vector and calls the vector-flavor
        // of compareRealUlp.
        ::testing::AssertionResult
        compareSimdRealUlp(const char * refExpr, const char * tstExpr, real ref, gmx_simd_real_t tst);
#endif

#ifdef GMX_SIMD_HAVE_INT32
        //! \brief Compare two int SIMD vectors exactly
        ::testing::AssertionResult
        compareSimdInt32(const char * refExpr, const char *  tstExpr,
                         gmx_simd_int32_t ref, gmx_simd_int32_t tst);

        //! \brief Compare a scalar int against all values in SIMD vector
        ::testing::AssertionResult
        compareSimdInt32(const char * refExpr, const char *  tstExpr,
                         int ref,              gmx_simd_int32_t tst);
#endif

    protected:
        gmx_int64_t            ulpTol_; //!< Default tolerance in units-in-last-position
        real                   absTol_; //!< Default absolute tolerance
};

#ifdef GMX_SIMD_HAVE_REAL
//! \brief Convert simd real to std::vector<real>
std::vector<real> simd2Vector(const gmx_simd_real_t simd);

//! \brief Convert std::vector<real> to simd real
gmx_simd_real_t   vector2Simd(const std::vector<real> &v);

//! \brief Create a simd from 3 reals
gmx_simd_real_t   setSimd3R(real r0, real r1, real r2);

//! Print simd real variable - convenient for debugging failing tests
std::string       printSimdReal(const gmx_simd_real_t simd);

#endif

//! \brief Test if a SIMD real is bitwise identical to reference SIMD or real scalar.
//
// If the reference value is a scalar, test is performed against all SIMD elements.
#define GMX_EXPECT_SIMD_REAL_EQ(ref, tst) EXPECT_PRED_FORMAT2(::simdTest::SimdTest::compareSimdRealUlp, ref, tst)

//! \brief Test if a SIMD real is within tolerance of reference SIMD or real scalar.
//
// If the reference value is a scalar, test is performed against all SIMD elements.
#define GMX_EXPECT_SIMD_REAL_NEAR(ref, tst) EXPECT_PRED_FORMAT2(::simdTest::SimdTest::compareSimdRealUlp, ref, tst)

#ifdef GMX_SIMD_HAVE_INT32
//! \brief Convert simd int to std::vector<int>
std::vector<int>  simd2Vector(const gmx_simd_int32_t simd);

//! \brief Convert std::vector<int> to simd int
gmx_simd_int32_t  vector2Simd(const std::vector<int> &v);

//! \brief Create a simd from 3 int
gmx_simd_int32_t  setSimd3I(int i0, int i1, int i2);

//! \brief Print simd int variable - convenient for debugging failing tests
std::string       printSimdInt32(const gmx_simd_int32_t simd);

//! \brief Macro that checks SIMD integer expression against SIMD or reference int.
//
// If the reference argument is a scalar integer it will be expanded into
// the width of the simd register and tested against all elements.
#define GMX_EXPECT_SIMD_INT_EQ(ref, tst)  EXPECT_PRED_FORMAT2(::simdTest::SimdTest::compareSimdInt32, ref, tst)

#endif // GMX_SIMD_HAVE_INT32

}      // namespace SimdTest

#endif // GMX_SIMD_TESTS_UTIL_H
