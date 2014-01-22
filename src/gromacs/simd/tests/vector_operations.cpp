/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"

#include "util.h"

namespace simdTest
{

#ifdef GMX_SIMD_HAVE_REAL

TEST(SimdTestVectorOps, gmxSimdCalcRsqR)
{
    gmx_simd_real_t simdX  = setSimd3R(1, 2, 3);
    gmx_simd_real_t simdY  = setSimd3R(3, 0, 5);
    gmx_simd_real_t simdZ  = setSimd3R(4, 1, 8);
    gmx_simd_real_t simdR2 = setSimd3R(26, 5, 98);

    GMX_ASSERT_SIMD_REAL_NEAR(simdR2, gmx_simd_calc_rsq_r(simdX, simdY, simdZ), 2, 0); // Allow 2 ulp difference
}

TEST(SimdTestVectorOps, gmxSimdIprodR)
{
    gmx_simd_real_t aX    = setSimd3R(1, 2, 3);
    gmx_simd_real_t aY    = setSimd3R(3, 0, 5);
    gmx_simd_real_t aZ    = setSimd3R(4, 1, 8);
    gmx_simd_real_t bX    = setSimd3R(8, 3, 6);
    gmx_simd_real_t bY    = setSimd3R(2, 3, 1);
    gmx_simd_real_t bZ    = setSimd3R(5, 7, 9);
    gmx_simd_real_t iprod = setSimd3R(34, 13, 95);
    GMX_ASSERT_SIMD_REAL_NEAR(iprod, gmx_simd_iprod_r(aX, aY, aZ, bX, bY, bZ), 2, 0); // Allow 2 ulp difference
}

TEST(SimdTestVectorOps, gmxSimdNorm2R)
{
    gmx_simd_real_t simdX     = setSimd3R(1, 2, 3);
    gmx_simd_real_t simdY     = setSimd3R(3, 0, 5);
    gmx_simd_real_t simdZ     = setSimd3R(4, 1, 8);
    gmx_simd_real_t simdNorm2 = setSimd3R(26, 5, 98);
    GMX_ASSERT_SIMD_REAL_NEAR(simdNorm2, gmx_simd_norm2_r(simdX, simdY, simdZ), 2, 0); // Allow 2 ulp difference
}

TEST(SimdTestBase, gmxSimdCprodR)
{
    gmx_simd_real_t aX    = setSimd3R(1, 2, 3);
    gmx_simd_real_t aY    = setSimd3R(3, 0, 5);
    gmx_simd_real_t aZ    = setSimd3R(4, 1, 8);
    gmx_simd_real_t bX    = setSimd3R(8, 3, 6);
    gmx_simd_real_t bY    = setSimd3R(2, 3, 1);
    gmx_simd_real_t bZ    = setSimd3R(5, 7, 9);
    gmx_simd_real_t refcX = setSimd3R(7, -3, 37);
    gmx_simd_real_t refcY = setSimd3R(27, -11, 21);
    gmx_simd_real_t refcZ = setSimd3R(-22, 6, -27);
    gmx_simd_real_t cX, cY, cZ;
    gmx_simd_cprod_r(aX, aY, aZ, bX, bY, bZ, &cX, &cY, &cZ);
    GMX_ASSERT_SIMD_REAL_NEAR(refcX, cX, 2, 0); // Allow 2 ulp difference
    GMX_ASSERT_SIMD_REAL_NEAR(refcY, cY, 2, 0); // Allow 2 ulp difference
    GMX_ASSERT_SIMD_REAL_NEAR(refcZ, cZ, 2, 0); // Allow 2 ulp difference
}

#endif // GMX_SIMD_HAVE_REAL

}
