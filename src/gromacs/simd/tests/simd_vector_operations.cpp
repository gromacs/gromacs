/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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

#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"

#include "data.h"
#include "simd.h"

#if GMX_SIMD

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

/*! \internal \brief Test fixture for vector operations tests (identical to the generic \ref SimdTest) */
typedef SimdTest SimdVectorOperationsTest;

TEST_F(SimdVectorOperationsTest, iprod)
{
    SimdReal aX       = rSimd_c0c1c2;
    SimdReal aY       = rSimd_c3c4c5;
    SimdReal aZ       = rSimd_c6c7c8;
    SimdReal bX       = rSimd_c3c0c4;
    SimdReal bY       = rSimd_c4c6c8;
    SimdReal bZ       = rSimd_c7c2c3;
    SimdReal iprodRef = setSimdRealFrom3R(c0*c3 + c3*c4 + c6*c7,
                                          c1*c0 + c4*c6 + c7*c2,
                                          c2*c4 + c5*c8 + c8*c3);

    setUlpTol(2);
    GMX_EXPECT_SIMD_REAL_NEAR(iprodRef, iprod(aX, aY, aZ, bX, bY, bZ));
}

TEST_F(SimdVectorOperationsTest, norm2)
{
    SimdReal simdX    = rSimd_c0c1c2;
    SimdReal simdY    = rSimd_c3c4c5;
    SimdReal simdZ    = rSimd_c6c7c8;
    SimdReal norm2Ref = setSimdRealFrom3R(c0*c0 + c3*c3 + c6*c6,
                                          c1*c1 + c4*c4 + c7*c7,
                                          c2*c2 + c5*c5 + c8*c8);

    setUlpTol(2);
    GMX_EXPECT_SIMD_REAL_NEAR(norm2Ref, norm2(simdX, simdY, simdZ));
}

TEST_F(SimdVectorOperationsTest, cprod)
{
    SimdReal aX    = rSimd_c0c1c2;
    SimdReal aY    = rSimd_c3c4c5;
    SimdReal aZ    = rSimd_c6c7c8;
    SimdReal bX    = rSimd_c3c0c4;
    SimdReal bY    = rSimd_c4c6c8;
    SimdReal bZ    = rSimd_c7c2c3;
    //The SIMD version might use FMA. If we don't force FMA for the reference value, the compiler is free to use FMA
    //for either product. If the compiler uses FMA for one product and the SIMD version uses FMA for the other, the
    //rounding error of each product adds up and the total possible ulp-error is 12.
    SimdReal refcX = setSimdRealFrom3R( std::fma(-c6, c4, c3*c7), std::fma(-c7, c6, c4*c2), std::fma(-c8, c8, c5*c3));
    SimdReal refcY = setSimdRealFrom3R( std::fma(-c0, c7, c6*c3), std::fma(-c1, c2, c7*c0), std::fma(-c2, c3, c8*c4));
    SimdReal refcZ = setSimdRealFrom3R( std::fma(-c3, c3, c0*c4), std::fma(-c4, c0, c1*c6), std::fma(-c5, c4, c2*c8));
    SimdReal cX, cY, cZ;

    //The test assumes that cprod uses FMA on architectures which have FMA so that the compiler can't choose which
    //product is computed with FMA.
    cprod(aX, aY, aZ, bX, bY, bZ, &cX, &cY, &cZ);

    //The test values cannot be computed without FMA for the case that SIMD has no FMA. Even if no explicit FMA were
    //used, the compiler could choose to use FMA. This causes up to 6upl error because of the product is up to 6 times
    //larger than the final result after the difference.
    setUlpTol(GMX_SIMD_HAVE_FMA ? ulpTol_ : 6);

    GMX_EXPECT_SIMD_REAL_NEAR(refcX, cX);
    GMX_EXPECT_SIMD_REAL_NEAR(refcY, cY);
    GMX_EXPECT_SIMD_REAL_NEAR(refcZ, cZ);
}

#endif      // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace

#endif // GMX_SIMD
