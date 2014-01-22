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
 * Declares test fixture for SIMD4 math tests.
 *
 * This extends the SimdTestBase class with operations to test
 * SIMD math functions over a specified range.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 * \ingroup module_simd
 */
#ifndef GMX_SIMD_TESTS_SIMD4_MATH_H
#define GMX_SIMD_TESTS_SIMD4_MATH_H

#include <vector>
#include <gtest/gtest.h>
#include "gromacs/simd/simd.h"

#include "simd_math.h"
#include "simd4_util.h"

namespace simd4Test
{

/*! \brief Test fixture to use for SIMD4 math funcations.
 *
 * This class simply augments the SimdTestMath fixture with a routine
 * that specifically compares a SIMD4 math function against a reference,
 * since we cannot rely on overloading working for SIMD datatypes.
 */
class Simd4TestMath : public ::simdTest::SimdTestMath
{
    public:
#ifdef GMX_SIMD4_HAVE_REAL
        /*! \brief Compare SIMD4 test math function against reference.
         *
         * This is an internal implementation routine. Do not use it directly,
         * but apply the macro GMX_EXPECT_SIMD4_FUNC_NEAR() instead.
         *
         * The test is performed over the input value range and tolerances specified in the class.
         */
            ::testing::AssertionResult
                                         compareSimd4MathFunction(const char * refFuncExpr, const char *tstFuncExpr,
                                                                  real refFunc(real x), gmx_simd4_real_t tstFunc(gmx_simd4_real_t x));
#endif
};

/*! \brief Test approximate equality of SIMD4 vs reference version of a function.
 *
 * This macro takes vanilla C and SIMD4 flavors of a function and tests it with
 * the number of points, range, and tolerances specified by the test fixture class.
 */
#define GMX_EXPECT_SIMD4_FUNC_NEAR(refFunc, tstFunc) \
    EXPECT_PRED_FORMAT2(::simd4Test::Simd4TestMath::compareSimd4MathFunction, refFunc, tstFunc)

}      // namespace SimdTest

#endif // GMX_SIMD_TESTS_SIMD4_MATH_H
