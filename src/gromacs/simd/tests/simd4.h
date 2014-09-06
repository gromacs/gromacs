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
#ifndef GMX_SIMD_TESTS_SIMD4_H
#define GMX_SIMD_TESTS_SIMD4_H

/*! \internal \file
 * \brief
 * Declares fixture for testing of SIMD4 functionality.
 *
 * This files specializes the common base test utilities to be used
 * for SIMD4 variables. For detailed documentation, check out the normal
 * SIMD test classes and files.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 * \ingroup module_simd
 */

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/simd/simd.h"

#include "base.h"

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#ifdef GMX_SIMD4_HAVE_REAL
extern const gmx_simd4_real_t rSimd4_1_2_3;     //!< Generic (different) fp values.
extern const gmx_simd4_real_t rSimd4_4_5_6;     //!< Generic (different) fp values.
extern const gmx_simd4_real_t rSimd4_7_8_9;     //!< Generic (different) fp values.
extern const gmx_simd4_real_t rSimd4_5_7_9;     //!< rSimd_1_2_3 + rSimd_4_5_6.
extern const gmx_simd4_real_t rSimd4_m1_m2_m3;  //!< Generic negative fp values.
extern const gmx_simd4_real_t rSimd4_3_1_4;     //!< Used to test min/max.
extern const gmx_simd4_real_t rSimd4_m3_m1_m4;  //!< negative rSimd_3_1_4.
extern const gmx_simd4_real_t rSimd4_2p25;      //!< Value that rounds down.
extern const gmx_simd4_real_t rSimd4_3p75;      //!< Value that rounds up.
extern const gmx_simd4_real_t rSimd4_m2p25;     //!< Negative value that rounds up.
extern const gmx_simd4_real_t rSimd4_m3p75;     //!< Negative value that rounds down.
//! Three large floating-point values whose exponents are >32.
extern const gmx_simd4_real_t rSimd4_Exp;
#    if (defined GMX_SIMD_HAVE_DOUBLE) && (defined GMX_DOUBLE)
// Make sure we also test exponents outside single precision when we use double
extern const gmx_simd4_real_t rSimd4_ExpDouble;
#    endif
extern const gmx_simd4_real_t rSimd4_Bits1; //!< Pattern F0 repeated to fill single/double.
extern const gmx_simd4_real_t rSimd4_Bits2; //!< Pattern CC repeated to fill single/double.
extern const gmx_simd4_real_t rSimd4_Bits3; //!< Pattern C0 repeated to fill single/double.
extern const gmx_simd4_real_t rSimd4_Bits4; //!< Pattern 0C repeated to fill single/double.
extern const gmx_simd4_real_t rSimd4_Bits5; //!< Pattern FC repeated to fill single/double.
extern const gmx_simd4_real_t rSimd4_Bits6; //!< Pattern 3C repeated to fill single/double.

/*! \internal
 * \brief
 * Test fixture for SIMD4 tests - contains test settings.
 *
 * This is a very simple test fixture that basically just takes the common
 * SIMD/SIMD4 functionality from SimdBaseTest and creates wrapper routines
 * specific for SIMD4 functionality.
 */
class Simd4Test : public SimdBaseTest
{
    public:
        /*! \brief Compare two real SIMD4 variables for approximate equality.
         *
         * This is an internal implementation routine. YOu should always use
         * GMX_EXPECT_SIMD4_REAL_NEAR() instead.
         *
         * This routine is designed according to the Google test specs, so the char
         * strings will describe the arguments to the macro.
         *
         * The comparison is applied to each element, and it returns true if each element
         * in the SIMD4 test variable is within the class tolerances of the corresponding
         * reference element.
         */
        ::testing::AssertionResult
        compareSimd4RealUlp(const char * refExpr, const char * tstExpr,
                            const gmx_simd4_real_t ref, const gmx_simd4_real_t tst);

        /*! \brief Compare two real SIMD4 variables for exact equality.
         *
         * This is an internal implementation routine. YOu should always use
         * GMX_EXPECT_SIMD4_REAL_NEAR() instead.
         *
         * This routine is designed according to the Google test specs, so the char
         * strings will describe the arguments to the macro.
         *
         * The comparison is applied to each element, and it returns true if each element
         * in the SIMD4 test variable is within the class tolerances of the corresponding
         * reference element.
         */
        ::testing::AssertionResult
        compareSimd4RealEq(const char * refExpr, const char * tstExpr,
                           const gmx_simd4_real_t ref, const gmx_simd4_real_t tst);
};

/*! \brief Convert SIMD4 real to std::vector<real>.
 *
 * The returned vector will have the same length as the SIMD4 width.
 */
std::vector<real> simd4Real2Vector(const gmx_simd4_real_t simd4);

/*! \brief Return floating-point SIMD4 value from std::vector<real>.
 *
 * If the vector is longer than SIMD4 width, only the first elements will be used.
 * If it is shorter, the contents will be repeated to fill the SIMD4 register.
 */
gmx_simd4_real_t   vector2Simd4Real(const std::vector<real> &v);

/*! \brief Set SIMD4 register contents from three real values.
 *
 * It might seem stupid to use three values when we know that the SIMD4 width
 * is 4, but it simplifies the test organization when the SIMD and SIMD4 tests
 * are completely symmetric.
 */
gmx_simd4_real_t   setSimd4RealFrom3R(real r0, real r1, real r2);

/*! \brief Set SIMD4 register contents from single real value.
 *
 * All elements is set from the given value. This is effectively the same
 * operation as gmx_simd4_set1_r(), but is implemented using only load/store
 * operations that have been tested separately in the bootstrapping tests.
 */
gmx_simd4_real_t   setSimd4RealFrom1R(real value);

/*! \brief Test if a SIMD4 real is bitwise identical to reference SIMD4 value. */
#define GMX_EXPECT_SIMD4_REAL_EQ(ref, tst)   EXPECT_PRED_FORMAT2(compareSimd4RealEq, ref, tst)

/*! \brief Test if a SIMD4 real is within tolerance of reference SIMD4 value. */
#define GMX_EXPECT_SIMD4_REAL_NEAR(ref, tst) EXPECT_PRED_FORMAT2(compareSimd4RealUlp, ref, tst)

#endif  // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace

#endif // GMX_SIMD_TESTS_SIMD4_H
