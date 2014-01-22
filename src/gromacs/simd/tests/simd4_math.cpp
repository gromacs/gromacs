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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include "gromacs/math/utilities.h"
#include "gromacs/simd/simd_math.h"

#include "simd4_math.h"

namespace simd4Test
{

#ifdef GMX_SIMD4_HAVE_REAL

/// \brief Internal implementation - use the macro GMX_EXPECT_SIMD4_FUNC_NEAR instead.
///
/// The reference and SIMD functions are called for 10,000 points over the
/// specified range. They are considered equal if they are either within
/// the absolute tolerance (with same sign), or within the ULP tolerance.
::testing::AssertionResult
Simd4TestMath::compareSimd4MathFunction(const char * refFuncExpr, const char *simdFuncExpr,
                                        real refFunc(real x),     gmx_simd4_real_t simdFunc(gmx_simd4_real_t x))
{
    std::vector<real>            vsimd;
    const int                    npoints = 100000;
    real                         x, dx, refval;
    gmx_simd4_real_t             simdx;
    gmx_int64_t                  ulpDiff, maxUlpDiff;
    real                         maxUlpDiffPos;
    real                         refValMaxUlpDiff, simdValMaxUlpDiff;
    bool                         eq, signOk;
    int                          i;
#    ifdef GMX_DOUBLE
    union {
        double r; gmx_int64_t i;
    } conv0, conv1;
#    else
    union {
        float  r; int i;
    } conv0, conv1;
#    endif

    maxUlpDiff = 0;
    dx         = (range.second-range.first)/npoints;

    for (x = range.first; x < range.second; x += dx)
    {
        simdx   = gmx_simd4_set1_r(x);
        refval  = refFunc(x);
        vsimd   = simd2Vector(simdFunc(simdx));

        for (i = 0, eq = true, signOk = true; i < GMX_SIMD4_WIDTH && eq == true; i++)
        {
            eq     = eq && ( fabs(refval-vsimd[i]) < absTol );
            signOk = signOk && ( refval*vsimd[i] >= 0 );
        }
        if (eq == true)
        {
            // Go to next point if everything within absolute tolerance
            continue;
        }
        else if (signOk == false)
        {
            return ::testing::AssertionFailure()
                   << "Failing math function comparison due to sign differences." << std::endl
                   << "Reference function: " << refFuncExpr << std::endl
                   << "Simd function:      " << simdFuncExpr << std::endl
                   << "Test range is ( " << range.first << " , " << range.second << " ) " << std::endl
                   << "First sign difference for x=" << std::setprecision(20) << x << std::endl
                   << "Ref value:  " << std::setprecision(20) << refval << std::endl
                   << "SIMD value: " << std::setprecision(20) << vsimd[0] << std::endl;
        }

        for (i = 0; i < GMX_SIMD4_WIDTH; i++)
        {
            conv0.r = refval;
            conv1.r = vsimd[i];
            ulpDiff = llabs(conv0.i-conv1.i);
            if (ulpDiff > maxUlpDiff)
            {
                maxUlpDiff        = ulpDiff;
                maxUlpDiffPos     = x;
                refValMaxUlpDiff  = refval;
                simdValMaxUlpDiff = vsimd[0];
            }
        }
    }

    if (maxUlpDiff <= ulpTol)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing ulp comparison between " << refFuncExpr << " and " << simdFuncExpr << std::endl
               << "Requested ulp tolerance: " << ulpTol << std::endl
               << "Requested abs tolerance: " << absTol << std::endl
               << "Largest Ulp difference occurs for x=" << std::setprecision(20) << maxUlpDiffPos << std::endl
               << "Ref  value: " << std::setprecision(20) << refValMaxUlpDiff << std::endl
               << "SIMD value: " << std::setprecision(20) << simdValMaxUlpDiff << std::endl
               << "Ulp diff.:   " << std::setprecision(20) << maxUlpDiff << std::endl;
    }
}

// Helper function wrapper
real ref_invsqrt(real x)
{
    return 1.0/sqrt(x);
}

TEST_F(Simd4TestMath, gmxSimd4InvsqrtR)
{
    setRange(1e-10,1e10);
    GMX_EXPECT_SIMD4_FUNC_NEAR(ref_invsqrt, gmx_simd4_invsqrt_r);
}

#endif

}
