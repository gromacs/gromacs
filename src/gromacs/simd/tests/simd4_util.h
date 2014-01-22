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
// Declares test fixture for SimdTests
//
// \author Erik Lindahl <erik.lindahl@scilifelab.se>
// \ingroup module_simd
//
#ifndef GMX_SIMD_TESTS_SIMD4_UTIL_H
#define GMX_SIMD_TESTS_SIMD4_UTIL_H

#include <vector>
#include <gtest/gtest.h>
#include "gromacs/simd/simd.h"
#include "util.h"

namespace simd4Test
{

#ifdef GMX_SIMD4_HAVE_REAL
extern const gmx_simd4_real_t rSimd_1_2_3;
extern const gmx_simd4_real_t rSimd_4_5_6;
extern const gmx_simd4_real_t rSimd_7_8_9;
extern const gmx_simd4_real_t rSimd_5_7_9;
extern const gmx_simd4_real_t rSimd_m1_m2_m3;
extern const gmx_simd4_real_t rSimd_3_1_4;
extern const gmx_simd4_real_t rSimd_m3_m1_m4;
extern const gmx_simd4_real_t rSimd_2p25;
extern const gmx_simd4_real_t rSimd_3p75;
extern const gmx_simd4_real_t rSimd_m2p25;
extern const gmx_simd4_real_t rSimd_m3p75;
extern const gmx_simd4_real_t rSimd_Exp;
// Magic FP numbers corresponding to specific bit patterns
extern const gmx_simd4_real_t rSimd_Bits1;
extern const gmx_simd4_real_t rSimd_Bits2;
extern const gmx_simd4_real_t rSimd_Bits3;
extern const gmx_simd4_real_t rSimd_Bits4;
extern const gmx_simd4_real_t rSimd_Bits5;
extern const gmx_simd4_real_t rSimd_Bits6;
// number of bits in default floating-point precision mantissa
extern const int              mantissaBits;
extern const int              defaultTol;

//! \brief Test fixture for SIMD tests - contains test settings.
//
// This is a very simple test fixture. It is only used to set default
// tolerances and ranges, mainly for floating-point and math tests, so we
// save a few operations in the test routines and can override the settings
// from the command line.
class Simd4Test : public ::simdTest::SimdTest
{
    public:

#ifdef GMX_SIMD4_HAVE_REAL
        //! \brief Compare two SIMD real.
        //
        // The comparison is applied to each element, and it returns true if all values
        // in the SIMD are within the class tolerances of the corresponding vector values.
            ::testing::AssertionResult
        compareSimd4RealUlp(const char * refExpr, const char * tstExpr, gmx_simd4_real_t ref, gmx_simd4_real_t tst);

        //! \brief Compare scalar against SIMD real.
        //
        // Convenience routine: expands the scalar into vector and calls the vector-flavor
        // of compareRealUlp.
        ::testing::AssertionResult
        compareSimd4RealUlp(const char * refExpr, const char * tstExpr, real ref, gmx_simd4_real_t tst);
#endif
};


//! \brief Convert simd real to std::vector<real>
std::vector<real> simd2Vector(const gmx_simd4_real_t simd);

//! \brief Convert std::vector<real> to simd real
gmx_simd4_real_t   vector2Simd(const std::vector<real> &v);

//! \brief Create a simd from 3 reals
gmx_simd4_real_t   setSimd3R(real r0, real r1, real r2);

//! \brief Print simd real variable - convenient for debugging failing tests
std::string       printSimd4Real(const gmx_simd4_real_t simd);

//! \brief Test if a SIMD4 real is bitwise identical to reference SIMD4 or real scalar.
//
// If the reference value is a scalar, test is performed against all SIMD4 elements.
#define GMX_EXPECT_SIMD4_REAL_EQ(ref, tst) EXPECT_PRED_FORMAT2(::simd4Test::Simd4Test::compareSimd4RealUlp, ref, tst)

//! \brief Test if a SIMD4 real is within tolerance of reference SIMD4 or real scalar.
//
// If the reference value is a scalar, test is performed against all SIMD4 elements.
#define GMX_EXPECT_SIMD4_REAL_NEAR(ref, tst) EXPECT_PRED_FORMAT2(::simd4Test::Simd4Test::compareSimd4RealUlp, ref, tst)

#endif // GMX_SIMD4_HAVE_REAL

}      // namespace simd4Test

#endif // GMX_SIMD_TESTS_SIMD4_UTIL_H
