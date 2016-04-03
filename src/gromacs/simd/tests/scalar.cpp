/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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

#include <cmath>

#include <gtest/gtest.h>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

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
    real ref = 1.234;
    real val = load(&ref);

    EXPECT_EQ(ref, val);
}

TEST(SimdScalarTest, loadU)
{
    real ref = 1.234;
    real val = loadU(&ref);

    EXPECT_EQ(ref, val);
}

TEST(SimdScalarTest, store)
{
    real ref = 1.234;
    real val;

    store(&val, ref);

    EXPECT_EQ(ref, val);
}

TEST(SimdScalarTest, storeU)
{
    real ref = 1.234;
    real val;

    store(&val, ref);

    EXPECT_EQ(ref, val);
}

TEST(SimdScalarTest, setZero)
{
    real val = setZero();

    EXPECT_EQ(0, val);
}

TEST(SimdScalarTest, andNot)
{
    real signBit = GMX_REAL_NEGZERO;

    EXPECT_EQ(real(1), andNot(signBit, real(-1)));
    EXPECT_EQ(real(2), andNot(signBit, real(2)));
}

TEST(SimdScalarTest, fma)
{
    EXPECT_EQ(real(27), fma(real(3), real(6), real(9))); // 3*6+9=27
}

TEST(SimdScalarTest, fms)
{
    EXPECT_EQ(real(9), fms(real(3), real(6), real(9))); // 3*6-9=9
}

TEST(SimdScalarTest, fnma)
{
    EXPECT_EQ(real(-9), fnma(real(3), real(6), real(9))); // -3*6+9=-9
}

TEST(SimdScalarTest, fnms)
{
    EXPECT_EQ(real(-27), fnms(real(3), real(6), real(9))); // -3*6-9=-27
}

TEST(SimdScalarTest, maskAdd)
{
    EXPECT_EQ(real(3), maskAdd(real(3), real(6), false));
    EXPECT_EQ(real(9), maskAdd(real(3), real(6), true));
}

TEST(SimdScalarTest, maskzMul)
{
    EXPECT_EQ(real(0), maskzMul(real(3), real(6), false));
    EXPECT_EQ(real(18), maskzMul(real(3), real(6), true));
}

TEST(SimdScalarTest, maskzFma)
{
    EXPECT_EQ(real(0), maskzFma(real(3), real(6), real(9), false));
    EXPECT_EQ(real(27), maskzFma(real(3), real(6), real(9), true));
}

TEST(SimdScalarTest, abs)
{
    EXPECT_EQ(real(3), abs(real(-3)));
}

TEST(SimdScalarTest, max)
{
    EXPECT_EQ(real(3), max(real(1), real(3)));
}

TEST(SimdScalarTest, min)
{
    EXPECT_EQ(real(1), min(real(1), real(3)));
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
    EXPECT_EQ(real(5), reduce(real(5)));
}

TEST(SimdScalarTest, testBits)
{
    EXPECT_TRUE(testBits(real(1)));
    EXPECT_TRUE(testBits(GMX_REAL_NEGZERO));
    EXPECT_FALSE(testBits(real(0)));
}

TEST(SimdScalarTest, anyTrue)
{
    EXPECT_TRUE(anyTrue(true));
    EXPECT_FALSE(anyTrue(false));
}

TEST(SimdScalarTest, selectByMask)
{
    EXPECT_EQ(real(5.0), selectByMask(real(5.0), true));
    EXPECT_EQ(real(0.0), selectByMask(real(5.0), false));
}

TEST(SimdScalarTest, selectByNotMask)
{
    EXPECT_EQ(real(0.0), selectByNotMask(real(5.0), true));
    EXPECT_EQ(real(5.0), selectByNotMask(real(5.0), false));
}

TEST(SimdScalarTest, blend)
{
    EXPECT_EQ(real(5.0), blend(real(3.0), real(5.0), true));
    EXPECT_EQ(real(3.0), blend(real(3.0), real(5.0), false));
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
    std::int32_t val = load(&ref);

    EXPECT_EQ(ref, val);
}

TEST(SimdScalarTest, loadUI)
{
    std::int32_t ref = 42;
    std::int32_t val = load(&ref);

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
    EXPECT_EQ(std::int32_t(0x0C0C0C0C),
              andNot(std::int32_t(0xF0F0F0F0), std::int32_t(0xCCCCCCCC)));
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

}      // namespace anonymous
}      // namespace test
}      // namespace gmx
