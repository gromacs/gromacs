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

/*! \internal \file
 * \brief
 * Declares test fixture for SimdTests.
 *
 * The SIMD tests are both simple and complicated. The actual testing logic
 * is \a very straightforward since we just need to test single values against
 * the math library, and for some math functions we need to do it in a loop.
 * This could have been achieved in minutes with the default Google Test tools,
 * if it wasn't for the problem that we cannot access or compare SIMD contents
 * directly without using lots of other SIMD functionality. For this reason
 * we have separate the basic testing of load/store operations into a separate
 * bootstrapping test. Once this works, we can use the utility routines to
 * convert SIMD contents to/from stl:vector<> and perform the rest of the tests.
 *
 * Be careful when editing this file. Since SIMD variables need to be aligned
 * in memory compilers appear to use various non-standard tricks to achieve
 * this, which means they might not behave like normal datatypes, for instance
 * when it comes to overloading or the ability to store static data in classes.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 * \ingroup module_simd
 */
#ifndef GMX_SIMD_TESTS_UTIL_H
#define GMX_SIMD_TESTS_UTIL_H

#include <vector>
#include <gtest/gtest.h>
#include "gromacs/simd/simd.h"

namespace simdTest
{

//! \cond internal
//! \addtogroup module_simd
//! \{

/* The SIMD module is special in that it is just a (big) collection of simple
 * macros that get compiled directly (typically 1-to-1) to assembly.
 * This means that the testing needs to assess a large number of very simple
 * macros that look like functions, although since they aren't proper functions
 * we cannot rely on C++ features such as overloading, which also means the
 * class below is mostly a utility container for initializing SIMD data.
 *
 * Another complication is that the width of the SIMD implementation will
 * depend on the hardware and precision. For some simple operations it is
 * sufficient to set all SIMD elements to the same value, and check that the
 * result is present in all elements. However, for a few more complex
 * instructions that might rely on shuffling under-the-hood it is important
 * that we can test operations with different elements. We achieve this by
 * having test code that can initialize a SIMD variable from an stl::vector
 * of arbitrary length; the vector is simply repeated to fill all elements in
 * the SIMD variable. We also have similar routines to compare a SIMD result
 * with values in a vector, which returns true iff all elements match.
 *
 * This way we can write simple tests that use different values for all SIMD
 * elements. Personally I like using vectors of length 3, since this means
 * there are no simple repeated patterns in low/high halves of SIMD variables
 * that are 2,4,8,or 16 elements wide, and we still don't have to care about
 * the exact SIMD width of the underlying implementation.
 *
 * Note that this utility uses a few SIMD load/store instructions internally -
 * those have been tested separately in the bootstrap_loadstore.cpp file.
 */

/* Note that the SIMD constants cannot be part of the class - some x86 compilers
 * have problems with the alignments of class-internal data required for AVX.
 */
#ifdef GMX_SIMD_HAVE_REAL
extern const gmx_simd_real_t rSimd_1_2_3;     //!< Generic (different) fp values.
extern const gmx_simd_real_t rSimd_4_5_6;     //!< Generic (different) fp values.
extern const gmx_simd_real_t rSimd_7_8_9;     //!< Generic (different) fp values.
extern const gmx_simd_real_t rSimd_5_7_9;     //!< rSimd_1_2_3 + rSimd_4_5_6.
extern const gmx_simd_real_t rSimd_m1_m2_m3;  //!< Generic negative floating-point values.
extern const gmx_simd_real_t rSimd_3_1_4;     //!< Used to test min/max.
extern const gmx_simd_real_t rSimd_m3_m1_m4;  //!< negative rSimd_3_1_4.
extern const gmx_simd_real_t rSimd_2p25;      //!< Value that rounds down.
extern const gmx_simd_real_t rSimd_3p75;      //!< Value that rounds up.
extern const gmx_simd_real_t rSimd_m2p25;     //!< Negative value that rounds up.
extern const gmx_simd_real_t rSimd_m3p75;     //!< Negative value that rounds down.
//! \brief Three large floating-point values whose exponents are >32.
extern const gmx_simd_real_t rSimd_Exp;
// Magic FP numbers corresponding to specific bit patterns
extern const gmx_simd_real_t rSimd_Bits1;       //!< Pattern F0 repeated to fill single/double.
extern const gmx_simd_real_t rSimd_Bits2;       //!< Pattern CC repeated to fill single/double.
extern const gmx_simd_real_t rSimd_Bits3;       //!< Pattern C0 repeated to fill single/double.
extern const gmx_simd_real_t rSimd_Bits4;       //!< Pattern 0C repeated to fill single/double.
extern const gmx_simd_real_t rSimd_Bits5;       //!< Pattern FC repeated to fill single/double.
extern const gmx_simd_real_t rSimd_Bits6;       //!< Pattern 3C repeated to fill single/double.
#endif                                          // GMX_SIMD_HAVE_REAL
#ifdef GMX_SIMD_HAVE_INT32
extern const gmx_simd_int32_t iSimd_1_2_3;      //!< Three generic ints.
extern const gmx_simd_int32_t iSimd_4_5_6;      //!< Three generic ints.
extern const gmx_simd_int32_t iSimd_7_8_9;      //!< Three generic ints.
extern const gmx_simd_int32_t iSimd_5_7_9;      //!< iSimd_1_2_3 + iSimd_4_5_6.
extern const gmx_simd_int32_t iSimd_1M_2M_3M;   //!< Term1 for 32bit add/sub.
extern const gmx_simd_int32_t iSimd_4M_5M_6M;   //!< Term2 for 32bit add/sub.
extern const gmx_simd_int32_t iSimd_5M_7M_9M;   //!< iSimd_1M_2M_3M + iSimd_4M_5M_6M.
extern const gmx_simd_int32_t iSimd_0xF0F0F0F0; //!< Bitpattern to test integer logical operations.
extern const gmx_simd_int32_t iSimd_0xCCCCCCCC; //!< Bitpattern to test integer logical operations.
#endif                                          // GMX_SIMD_HAVE_INT32

/*! \brief Test fixture for SIMD tests - contains test settings.
 *
 * This is a very simple test fixture. It is only used to set default
 * tolerances, mainly for floating-point tests, so we
 * save a few operations in the test routines.
 */
class SimdTest : public ::testing::Test
{
    public:
        /*! \brief Initialize new SIMD test fixture with default tolerances.
         *
         * The default absolute tolerance is set to 0, which means the we always
         * check the ulp tolerance by default (passing the absolute tolerance
         * test would otherwise mean we approve the test instantly).
         * The default ulp tolerance is set to 10 units in single, and 255 units
         * in double precision.
         * Most SIMD math functions achieve 2-3 ulp accuracy in single, but by
         * being a bit liberal we avoid tests failing on aggressive compilers.
         *
         * For double precision we only aim to achieve twice the accuracy of
         * single. This way we can make do with a single extra iteration
         * in some algorithms, in particular 1/sqrt(x).
         */
        SimdTest()
        {
#ifdef GMX_DOUBLE
            ulpTol_       = 255LL; // Aim for roughly twice the precision we have in single.
            mantissaBits_ = 52;
#else
            ulpTol_       = 10LL;  // Be a bit liberal so compiler optimization doesn't bite us.
            mantissaBits_ = 23;
#endif
            absTol_ = 0;
        }
        /*! \brief Adjust ulp tolerance from the default 10 (float) or 255 (double). */
        void setUlpTol(gmx_int64_t newTol)   { ulpTol_ = newTol; }
        /*! \brief Adjust the absolute tolerance from the default 0.
         *
         * If values are closer than the absolute tolerance, the test will pass
         * no matter what their ulp difference is.
         */
        void setAbsTol(real newTol)          { absTol_ = newTol; }

#ifdef GMX_SIMD_HAVE_REAL
        /*! \brief Compare two real SIMD variables.
         *
         * This is an internal implementation routine. YOu should always use
         * GMX_EXPECT_SIMD_REAL_EQ() instead.
         *
         * This routine is designed according to the Google test specs, so the char
         * strings will describe the arguments to the macro.
         *
         * The comparison is applied to each element, and it returns true if each element
         * in the SIMD test variable is within the class tolerances of the corresponding
         * reference element.
         */
        ::testing::AssertionResult
        compareSimdRealUlp(const char * refExpr, const char * tstExpr,
                           const gmx_simd_real_t ref, const gmx_simd_real_t tst);

        /*! \brief Compare scalar against SIMD real.
         *
         * This is an internal implementation routine. YOu should always use
         * GMX_EXPECT_SIMD_REAL_EQ() instead.
         *
         * This version expands the scalar into a SIMD vector and calls the
         * all-SIMD-flavor of compareRealUlp. It will thus return success if
         * all elements are equal to the scalar.
         */
        ::testing::AssertionResult
        compareSimdRealUlp(const char * refExpr, const char * tstExpr,
                           real ref, const gmx_simd_real_t tst);
#endif

#ifdef GMX_SIMD_HAVE_INT32
        /*! \brief Compare two 32-bit integer SIMD variables.
         *
         * This is an internal implementation routine. YOu should always use
         * GMX_EXPECT_SIMD_INT_EQ() instead.
         *
         * This routine is designed according to the Google test specs, so the char
         * strings will describe the arguments to the macro, while the SIMD and
         * tolerance arguments are used to decide if the values are approximately equal.
         *
         * The comparison is applied to each element, and it returns true if each element
         * in the SIMD variable tst is identical to the corresponding reference element.
         */
        ::testing::AssertionResult
        compareSimdInt32(const char * refExpr, const char *  tstExpr,
                         const gmx_simd_int32_t ref, const gmx_simd_int32_t tst);

        /*! \brief Compare scalar integer against SIMD 32-bit integer.
         *
         * This is an internal implementation routine. YOu should always use
         * GMX_EXPECT_SIMD_INT_EQ() instead.
         *
         * This version expands the scalar into a SIMD vector and calls the
         * all-SIMD-flavor of compareSimdInt32. It will thus return success if
         * all elements are equal to the scalar reference.
         */
        ::testing::AssertionResult
        compareSimdInt32(const char * refExpr, const char *  tstExpr,
                         int ref, const gmx_simd_int32_t tst);
#endif

    protected:
        gmx_int64_t            ulpTol_;       //!< Current tolerance in units-in-last-position.
        real                   absTol_;       //!< Current absolute tolerance.
        int                    mantissaBits_; //!< Number of bits in real SIMD mantissa.
};

#ifdef GMX_SIMD_HAVE_REAL
/*! \brief Convert SIMD real to std::vector<real>.
 *
 * The returned vector will have the same length as the SIMD width.
 */
std::vector<real> simdReal2Vector(const gmx_simd_real_t simd);

/*! \brief Return floating-point SIMD value from std::vector<real>.
 *
 * If the vector is longer than SIMD width, only the first elements will be used.
 * If it is shorter, the contents will be repeated to fill the SIMD register.
 */
gmx_simd_real_t   vector2SimdReal(const std::vector<real> &v);

/*! \brief Set SIMD register contents from three real values.
 *
 * Our reason for using three values is that 3 is not a factor in any known
 * SIMD width, so this way there will not be any simple repeated patterns e.g.
 * between the low/high 64/128/256 bits in the SIMD register, which could hide bugs.
 */
gmx_simd_real_t   setSimd3R(real r0, real r1, real r2);

/*! \brief Print SIMD real datatype.
 *
 * Simd datatypes are typically implemented as quite low-level
 * compiler-specific constructs, and at least for intel we have not
 * managed to apply the operator <<, but use a routine instead.
 */
std::string       printSimdReal(const gmx_simd_real_t simd);

#endif // GMX_SIMD_HAVE_REAL

/*! \brief Test if a SIMD real is bitwise identical to reference SIMD or real scalar.
 *
 * If the reference value is a scalar, test is performed against all SIMD elements.
 */
#define GMX_EXPECT_SIMD_REAL_EQ(ref, tst) EXPECT_PRED_FORMAT2(::simdTest::SimdTest::compareSimdRealUlp, ref, tst)

/*! \brief Test if a SIMD real is within tolerance of reference SIMD or real scalar.
 *
 * If the reference value is a scalar, test is performed against all SIMD elements.
 */
#define GMX_EXPECT_SIMD_REAL_NEAR(ref, tst) EXPECT_PRED_FORMAT2(::simdTest::SimdTest::compareSimdRealUlp, ref, tst)

#ifdef GMX_SIMD_HAVE_INT32
/*! \brief Convert SIMD integer to std::vector<int>.
 *
 * The returned vector will have the same length as the SIMD width.
 */
std::vector<int>  simdInt2Vector(const gmx_simd_int32_t simd);

/*! \brief Return 32-bit integer SIMD value from std::vector<int>.
 *
 * If the vector is longer than SIMD width, only the first elements will be used.
 * If it is shorter, the contents will be repeated to fill the SIMD register.
 */
gmx_simd_int32_t  vector2SimdInt(const std::vector<int> &v);

/*! \brief Set SIMD register contents from three int values.
 *
 * Our reason for using three values is that 3 is not a factor in any known
 * SIMD width, so this way there will not be any simple repeated patterns e.g.
 * between the low/high 64/128/256 bits in the SIMD register, which could hide bugs.
 */
gmx_simd_int32_t  setSimd3I(int i0, int i1, int i2);

/*! \brief Print SIMD 32-bit integer datatype.
 *
 * Simd datatypes are typically implemented as quite low-level
 * compiler-specific constructs, and at least for intel we have not
 * managed to apply the operator <<, but use a routine instead.
 */
std::string       printSimdInt32(const gmx_simd_int32_t simd);

/*! \brief Macro that checks SIMD integer expression against SIMD or reference int.
 *
 * If the reference argument is a scalar integer it will be expanded into
 * the width of the SIMD register and tested against all elements.
 */
#define GMX_EXPECT_SIMD_INT_EQ(ref, tst)  EXPECT_PRED_FORMAT2(::simdTest::SimdTest::compareSimdInt32, ref, tst)

#endif // GMX_SIMD_HAVE_INT32

//! \}
//! \endcond

}      // namespace SimdTest

#endif // GMX_SIMD_TESTS_UTIL_H
