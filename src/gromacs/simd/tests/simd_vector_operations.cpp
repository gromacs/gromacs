/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
    SimdReal aX       = setSimdRealFrom3R(1, 2, 3);
    SimdReal aY       = setSimdRealFrom3R(3, 0, 5);
    SimdReal aZ       = setSimdRealFrom3R(4, 1, 8);
    SimdReal bX       = setSimdRealFrom3R(8, 3, 6);
    SimdReal bY       = setSimdRealFrom3R(2, 3, 1);
    SimdReal bZ       = setSimdRealFrom3R(5, 7, 9);
    SimdReal iprodRef = setSimdRealFrom3R(34, 13, 95);

    setUlpTol(2);
    GMX_EXPECT_SIMD_REAL_NEAR(iprodRef, iprod(aX, aY, aZ, bX, bY, bZ));
}

TEST_F(SimdVectorOperationsTest, norm2)
{
    SimdReal simdX    = setSimdRealFrom3R(1, 2, 3);
    SimdReal simdY    = setSimdRealFrom3R(3, 0, 5);
    SimdReal simdZ    = setSimdRealFrom3R(4, 1, 8);
    SimdReal norm2Ref = setSimdRealFrom3R(26, 5, 98);

    setUlpTol(2);
    GMX_EXPECT_SIMD_REAL_NEAR(norm2Ref, norm2(simdX, simdY, simdZ));
}

TEST_F(SimdVectorOperationsTest, cprod)
{
    SimdReal aX    = setSimdRealFrom3R(1, 2, 3);
    SimdReal aY    = setSimdRealFrom3R(3, 0, 5);
    SimdReal aZ    = setSimdRealFrom3R(4, 1, 8);
    SimdReal bX    = setSimdRealFrom3R(8, 3, 6);
    SimdReal bY    = setSimdRealFrom3R(2, 3, 1);
    SimdReal bZ    = setSimdRealFrom3R(5, 7, 9);
    SimdReal refcX = setSimdRealFrom3R(7, -3, 37);
    SimdReal refcY = setSimdRealFrom3R(27, -11, 21);
    SimdReal refcZ = setSimdRealFrom3R(-22, 6, -27);
    SimdReal cX, cY, cZ;

    cprod(aX, aY, aZ, bX, bY, bZ, &cX, &cY, &cZ);

    setUlpTol(2);
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
