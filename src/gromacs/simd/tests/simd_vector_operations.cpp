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
#include "gmxpre.h"

#include <math.h>

#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"

#include "simd.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#ifdef GMX_SIMD_HAVE_REAL

/*! \internal \brief Test fixture for vector operations tests (identical to the generic \ref SimdTest) */
typedef SimdTest SimdVectorOperationsTest;

TEST_F(SimdVectorOperationsTest, gmxSimdCalcRsqR)
{
    gmx_simd_real_t simdX  = setSimdRealFrom3R(1, 2, 3);
    gmx_simd_real_t simdY  = setSimdRealFrom3R(3, 0, 5);
    gmx_simd_real_t simdZ  = setSimdRealFrom3R(4, 1, 8);
    gmx_simd_real_t simdR2 = setSimdRealFrom3R(26, 5, 98);

    setUlpTol(2);
    GMX_EXPECT_SIMD_REAL_NEAR(simdR2, gmx_simd_calc_rsq_r(simdX, simdY, simdZ));
}

TEST_F(SimdVectorOperationsTest, gmxSimdIprodR)
{
    gmx_simd_real_t aX    = setSimdRealFrom3R(1, 2, 3);
    gmx_simd_real_t aY    = setSimdRealFrom3R(3, 0, 5);
    gmx_simd_real_t aZ    = setSimdRealFrom3R(4, 1, 8);
    gmx_simd_real_t bX    = setSimdRealFrom3R(8, 3, 6);
    gmx_simd_real_t bY    = setSimdRealFrom3R(2, 3, 1);
    gmx_simd_real_t bZ    = setSimdRealFrom3R(5, 7, 9);
    gmx_simd_real_t iprod = setSimdRealFrom3R(34, 13, 95);

    setUlpTol(2);
    GMX_EXPECT_SIMD_REAL_NEAR(iprod, gmx_simd_iprod_r(aX, aY, aZ, bX, bY, bZ));
}

TEST_F(SimdVectorOperationsTest, gmxSimdNorm2R)
{
    gmx_simd_real_t simdX     = setSimdRealFrom3R(1, 2, 3);
    gmx_simd_real_t simdY     = setSimdRealFrom3R(3, 0, 5);
    gmx_simd_real_t simdZ     = setSimdRealFrom3R(4, 1, 8);
    gmx_simd_real_t simdNorm2 = setSimdRealFrom3R(26, 5, 98);

    setUlpTol(2);
    GMX_EXPECT_SIMD_REAL_NEAR(simdNorm2, gmx_simd_norm2_r(simdX, simdY, simdZ));
}

TEST_F(SimdVectorOperationsTest, gmxSimdCprodR)
{
    gmx_simd_real_t aX    = setSimdRealFrom3R(1, 2, 3);
    gmx_simd_real_t aY    = setSimdRealFrom3R(3, 0, 5);
    gmx_simd_real_t aZ    = setSimdRealFrom3R(4, 1, 8);
    gmx_simd_real_t bX    = setSimdRealFrom3R(8, 3, 6);
    gmx_simd_real_t bY    = setSimdRealFrom3R(2, 3, 1);
    gmx_simd_real_t bZ    = setSimdRealFrom3R(5, 7, 9);
    gmx_simd_real_t refcX = setSimdRealFrom3R(7, -3, 37);
    gmx_simd_real_t refcY = setSimdRealFrom3R(27, -11, 21);
    gmx_simd_real_t refcZ = setSimdRealFrom3R(-22, 6, -27);
    gmx_simd_real_t cX, cY, cZ;
    gmx_simd_cprod_r(aX, aY, aZ, bX, bY, bZ, &cX, &cY, &cZ);

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
