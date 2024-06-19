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
#include "gmxpre.h"

#include <cmath>
#include <cstdint>

#include <string>

#include <gtest/gtest.h>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

#include "data.h"

namespace gmx
{
namespace test
{

namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

/************************
 * Floating-point tests *
 ************************/

TEST(SimdScalarTest, load)
{
    real val = load<real>(&c1);

    EXPECT_EQ(c1, val);
}

TEST(SimdScalarTest, loadU)
{
    real val = loadU<real>(&c1);

    EXPECT_EQ(c1, val);
}

TEST(SimdScalarTest, store)
{
    real val;

    store(&val, c1);

    EXPECT_EQ(c1, val);
}

TEST(SimdScalarTest, storeU)
{
    real val;

    store(&val, c1);

    EXPECT_EQ(c1, val);
}

TEST(SimdScalarTest, setZero)
{
    real val = setZero();

    EXPECT_EQ(0, val);
}

TEST(SimdScalarTest, andNot)
{
    real signBit = GMX_REAL_NEGZERO;

    EXPECT_EQ(c1, andNot(signBit, -c1));
    EXPECT_EQ(c2, andNot(signBit, c2));
}

TEST(SimdScalarTest, fma)
{
    EXPECT_REAL_EQ_TOL(c1 * c2 + c3, fma(real(c1), real(c2), real(c3)), defaultRealTolerance());
}

TEST(SimdScalarTest, fms)
{
    EXPECT_REAL_EQ_TOL(c1 * c2 - c3, fms(c1, c2, c3), defaultRealTolerance());
}

TEST(SimdScalarTest, fnma)
{
    EXPECT_REAL_EQ_TOL(-c1 * c2 + c3, fnma(c1, c2, c3), defaultRealTolerance());
}

TEST(SimdScalarTest, fnms)
{
    EXPECT_REAL_EQ_TOL(-c1 * c2 - c3, fnms(c1, c2, c3), defaultRealTolerance());
}

TEST(SimdScalarTest, maskAdd)
{
    EXPECT_REAL_EQ_TOL(c1, maskAdd(c1, c2, false), defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(c1 + c2, maskAdd(c1, c2, true), defaultRealTolerance());
}

TEST(SimdScalarTest, maskzMul)
{
    EXPECT_REAL_EQ_TOL(czero, maskzMul(c1, c2, false), defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(c1 * c2, maskzMul(c1, c2, true), defaultRealTolerance());
}

TEST(SimdScalarTest, maskzFma)
{
    EXPECT_REAL_EQ_TOL(czero, maskzFma(c1, c2, c3, false), defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(c1 * c2 + c3, maskzFma(c1, c2, c3, true), defaultRealTolerance());
}

TEST(SimdScalarTest, abs)
{
    EXPECT_EQ(c1, abs(-c1));
}

TEST(SimdScalarTest, max)
{
    EXPECT_EQ(c3, max(c1, c3));
}

TEST(SimdScalarTest, min)
{
    EXPECT_EQ(c1, min(c1, c3));
}

TEST(SimdScalarTest, round)
{
    EXPECT_EQ(real(2), round(real(2.25)));
    EXPECT_EQ(real(4), round(real(3.75)));
    EXPECT_EQ(real(-2), round(real(-2.25)));
    EXPECT_EQ(real(-4), round(real(-3.75)));
}

TEST(SimdScalarTest, trunc)
{
    EXPECT_EQ(real(2), trunc(real(2.25)));
    EXPECT_EQ(real(3), trunc(real(3.75)));
    EXPECT_EQ(real(-2), trunc(real(-2.25)));
    EXPECT_EQ(real(-3), trunc(real(-3.75)));
}

TEST(SimdScalarTest, reduce)
{
    EXPECT_EQ(c1, reduce(c1));
}

TEST(SimdScalarTest, testBits)
{
    EXPECT_TRUE(testBits(c1));
    EXPECT_TRUE(testBits(GMX_REAL_NEGZERO));
    EXPECT_FALSE(testBits(czero));
}

TEST(SimdScalarTest, anyTrue)
{
    EXPECT_TRUE(anyTrue(true));
    EXPECT_FALSE(anyTrue(false));
}

TEST(SimdScalarTest, selectByMask)
{
    EXPECT_EQ(c1, selectByMask(c1, true));
    EXPECT_EQ(czero, selectByMask(c1, false));
}

TEST(SimdScalarTest, selectByNotMask)
{
    EXPECT_EQ(czero, selectByNotMask(c1, true));
    EXPECT_EQ(c1, selectByNotMask(c1, false));
}

TEST(SimdScalarTest, blend)
{
    EXPECT_EQ(c2, blend(c1, c2, true));
    EXPECT_EQ(c1, blend(c1, c2, false));
}

TEST(SimdScalarTest, cvtR2I)
{
    EXPECT_EQ(std::int32_t(4), cvtR2I(3.75));
    EXPECT_EQ(std::int32_t(-4), cvtR2I(-3.75));
}

TEST(SimdScalarTest, cvttR2I)
{
    EXPECT_EQ(std::int32_t(3), cvttR2I(3.75));
    EXPECT_EQ(std::int32_t(-3), cvttR2I(-3.75));
}

TEST(SimdScalarTest, cvtI2R)
{
    EXPECT_EQ(real(2.0), cvtI2R(2));
    EXPECT_EQ(real(-2.0), cvtI2R(-2));
}

TEST(SimdScalarTest, cvtF2D)
{
    float f = 1.23456789;

    EXPECT_EQ(double(f), cvtF2D(f));
}

TEST(SimdScalarTest, cvtD2D)
{
    double d = 1.23456789;

    EXPECT_EQ(float(d), cvtD2F(d));
}


/*****************
 * Integer tests *
 *****************/


TEST(SimdScalarTest, loadI)
{
    std::int32_t ref = 42;
    std::int32_t val = load<int32_t>(&ref);

    EXPECT_EQ(ref, val);
}

TEST(SimdScalarTest, loadUI)
{
    std::int32_t ref = 42;
    std::int32_t val = loadU<int32_t>(&ref);

    EXPECT_EQ(ref, val);
}

TEST(SimdScalarTest, storeI)
{
    std::int32_t ref = 42;
    std::int32_t val;

    storeU(&val, ref);

    EXPECT_EQ(ref, val);
}

TEST(SimdScalarTest, storeUI)
{
    std::int32_t ref = 42;
    std::int32_t val;

    storeU(&val, ref);

    EXPECT_EQ(ref, val);
}

TEST(SimdScalarTest, andNotI)
{
    EXPECT_EQ(std::int32_t(0x0C0C0C0C), andNot(std::int32_t(0xF0F0F0F0), std::int32_t(0xCCCCCCCC)));
}

TEST(SimdScalarTest, testBitsI)
{
    EXPECT_TRUE(testBits(1));
    EXPECT_FALSE(testBits(0));
}

TEST(SimdScalarTest, selectByMaskI)
{
    EXPECT_EQ(5, selectByMask(5, true));
    EXPECT_EQ(0, selectByMask(5, false));
}

TEST(SimdScalarTest, selectByNotMaskI)
{
    EXPECT_EQ(0, selectByNotMask(5, true));
    EXPECT_EQ(5, selectByNotMask(5, false));
}

TEST(SimdScalarTest, blendI)
{
    EXPECT_EQ(5, blend(3, 5, true));
    EXPECT_EQ(3, blend(3, 5, false));
}

TEST(SimdScalarTest, cvtB2IB)
{
    EXPECT_TRUE(cvtB2IB(true));
    EXPECT_FALSE(cvtB2IB(false));
}

TEST(SimdScalarTest, cvtIB2B)
{
    EXPECT_TRUE(cvtIB2B(true));
    EXPECT_FALSE(cvtIB2B(false));
}

/*! \} */
/*! \endcond internal */

} // namespace
} // namespace test
} // namespace gmx
