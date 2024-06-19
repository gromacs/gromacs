/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
#include "gromacs/utility/real.h"

#include "base.h"
#include "data.h"

#if GMX_SIMD

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#    if GMX_SIMD4_HAVE_REAL
extern const Simd4Real rSimd4_c0c1c2; //!< c0,c1,c2 repeated
extern const Simd4Real rSimd4_c3c4c5; //!< c3,c4,c5 repeated
extern const Simd4Real rSimd4_c6c7c8; //!< c6,c7,c8 repeated
extern const Simd4Real rSimd4_c3c0c4; //!< c3,c0,c4 repeated
extern const Simd4Real rSimd4_c4c6c8; //!< c4,c6,c8 repeated
extern const Simd4Real rSimd4_c7c2c3; //!< c7,c2,c3 repeated
extern const Simd4Real rSimd4_m0m1m2; //!< -c0,-c1,-c2 repeated
extern const Simd4Real rSimd4_m3m0m4; //!< -c3,-c0,-c4 repeated
extern const Simd4Real rSimd4_2p25;   //!< Value that rounds down.
extern const Simd4Real rSimd4_3p75;   //!< Value that rounds up.
extern const Simd4Real rSimd4_m2p25;  //!< Negative value that rounds up.
extern const Simd4Real rSimd4_m3p75;  //!< Negative value that rounds down.
//! Three large floating-point values whose exponents are >32.
extern const Simd4Real rSimd4_Exp;

#        if GMX_SIMD_HAVE_LOGICAL
extern const Simd4Real rSimd4_logicalA;         //!< Bit pattern to test logical ops
extern const Simd4Real rSimd4_logicalB;         //!< Bit pattern to test logical ops
extern const Simd4Real rSimd4_logicalResultOr;  //!< Result or bitwise 'or' of A and B
extern const Simd4Real rSimd4_logicalResultAnd; //!< Result or bitwise 'and' of A and B
#        endif                                  // GMX_SIMD_HAVE_LOGICAL

extern const Simd4Real rSimd4_Bits1; //!< Pattern F0 repeated to fill single/double.
extern const Simd4Real rSimd4_Bits2; //!< Pattern CC repeated to fill single/double.
extern const Simd4Real rSimd4_Bits3; //!< Pattern C0 repeated to fill single/double.
extern const Simd4Real rSimd4_Bits4; //!< Pattern 0C repeated to fill single/double.
extern const Simd4Real rSimd4_Bits5; //!< Pattern FC repeated to fill single/double.
extern const Simd4Real rSimd4_Bits6; //!< Pattern 3C repeated to fill single/double.

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
    compareSimd4RealUlp(const char* refExpr, const char* tstExpr, Simd4Real ref, Simd4Real tst);

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
    compareSimd4RealEq(const char* refExpr, const char* tstExpr, Simd4Real ref, Simd4Real tst);
};

/*! \brief Convert SIMD4 real to std::vector<real>.
 *
 * The returned vector will have the same length as the SIMD4 width.
 */
std::vector<real> simd4Real2Vector(Simd4Real simd4);

/*! \brief Return floating-point SIMD4 value from std::vector<real>.
 *
 * If the vector is longer than SIMD4 width, only the first elements will be used.
 * If it is shorter, the contents will be repeated to fill the SIMD4 register.
 */
Simd4Real vector2Simd4Real(const std::vector<real>& v);

/*! \brief Set SIMD4 register contents from three real values.
 *
 * It might seem stupid to use three values when we know that the SIMD4 width
 * is 4, but it simplifies the test organization when the SIMD and SIMD4 tests
 * are completely symmetric.
 */
Simd4Real setSimd4RealFrom3R(real r0, real r1, real r2);

/*! \brief Set SIMD4 register contents from single real value.
 *
 * All elements is set from the given value. This is effectively the same
 * operation as simd4Set1(), but is implemented using only load/store
 * operations that have been tested separately in the bootstrapping tests.
 */
Simd4Real setSimd4RealFrom1R(real value);

/*! \brief Test if a SIMD4 real is bitwise identical to reference SIMD4 value. */
#        define GMX_EXPECT_SIMD4_REAL_EQ(ref, tst) EXPECT_PRED_FORMAT2(compareSimd4RealEq, ref, tst)

/*! \brief Test if a SIMD4 real is within tolerance of reference SIMD4 value. */
#        define GMX_EXPECT_SIMD4_REAL_NEAR(ref, tst) \
            EXPECT_PRED_FORMAT2(compareSimd4RealUlp, ref, tst)

#    endif // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

} // namespace test
} // namespace gmx

#endif // GMX_SIMD

#endif // GMX_SIMD_TESTS_SIMD4_H
