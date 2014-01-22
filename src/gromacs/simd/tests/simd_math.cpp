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
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/options/basicoptions.h"
#include "simd_math.h"

namespace gmx
{
namespace test
{

int  SimdMathTest::s_nPoints    = 10000;

GMX_TEST_OPTIONS(SimdMathTestOptions, options)
{
    options->addOption(::gmx::IntegerOption("npoints")
                           .store(&SimdMathTest::s_nPoints)
                           .description("Number of points to test for SIMD math functions"));
}

#ifdef GMX_SIMD_HAVE_REAL

::testing::AssertionResult
SimdMathTest::compareSimdMathFunction(const char * refFuncExpr, const char *simdFuncExpr,
                                      real refFunc(real x),     gmx_simd_real_t simdFunc(gmx_simd_real_t x))
{
    std::vector<real>            vsimd;
    real                         x, dx, refval;
    gmx_simd_real_t              simdx;
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
    dx         = (range_.second-range_.first)/SimdMathTest::s_nPoints;

    for (x = range_.first; x < range_.second; x += dx)
    {
        simdx   = gmx_simd_set1_r(x);
        refval  = refFunc(x);
        vsimd   = simdReal2Vector(simdFunc(simdx));

        for (i = 0, eq = true, signOk = true; i < GMX_SIMD_REAL_WIDTH && eq == true; i++)
        {
            eq     = eq && ( fabs(refval-vsimd[i]) < absTol_ );
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
                   << "Test range is ( " << range_.first << " , " << range_.second << " ) " << std::endl
                   << "First sign difference for x=" << std::setprecision(20) << x << std::endl
                   << "Ref value:  " << std::setprecision(20) << refval << std::endl
                   << "SIMD value: " << std::setprecision(20) << vsimd[0] << std::endl;
        }

        for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
        {
            conv0.r = refval;
            conv1.r = vsimd[i];
            ulpDiff = llabs(conv0.i-conv1.i);
            if (ulpDiff > maxUlpDiff)
            {
                maxUlpDiff        = ulpDiff;
                maxUlpDiffPos     = x;
                refValMaxUlpDiff  = refval;
                simdValMaxUlpDiff = vsimd[i];
            }
        }
    }

    if (maxUlpDiff <= ulpTol_)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing ulp comparison between " << refFuncExpr << " and " << simdFuncExpr << std::endl
               << "Requested ulp tolerance: " << ulpTol_ << std::endl
               << "Requested abs tolerance: " << absTol_ << std::endl
               << "Largest Ulp difference occurs for x=" << std::setprecision(20) << maxUlpDiffPos << std::endl
               << "Ref  value: " << std::setprecision(20) << refValMaxUlpDiff << std::endl
               << "SIMD value: " << std::setprecision(20) << simdValMaxUlpDiff << std::endl
               << "Ulp diff.:   " << std::setprecision(20) << maxUlpDiff << std::endl;
    }
}

// Actual math function tests below


namespace
{

TEST_F(SimdMathTest, gmxSimdXorSignR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdFrom3R(-4, 5, 6), gmx_simd_xor_sign_r(setSimdFrom3R(4, 5, 6), setSimdFrom3R(-5, 2, 0)));
    // Make sure that we use the signbit correctly from negative zero too.
    GMX_EXPECT_SIMD_REAL_EQ(setSimdFrom3R(4, -5, 6), gmx_simd_xor_sign_r(setSimdFrom3R(-4, -5, -6), setSimdFrom3R(-5, 0.0, -0.0)));
}

// Helper function wrapper
real ref_invsqrt(real x)
{
    return 1.0/sqrt(x);
}

TEST_F(SimdMathTest, gmxSimdInvsqrtR)
{
    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_invsqrt, gmx_simd_invsqrt_r);
}

// Helper function wrappers to use one argument and return value from invsqrt_pair
gmx_simd_real_t
tst_invsqrt_pair0(gmx_simd_real_t x)
{
    gmx_simd_real_t r0, r1;
    gmx_simd_invsqrt_pair_r(x, x, &r0, &r1);
    return r0;
}
gmx_simd_real_t
tst_invsqrt_pair1(gmx_simd_real_t x)
{
    gmx_simd_real_t r0, r1;
    gmx_simd_invsqrt_pair_r(x, x, &r0, &r1);
    return r1;
}

TEST_F(SimdMathTest, gmxSimdInvsqrtPairR)
{
    setRange(1e-10, 1e10);
    // The accuracy conversions lose a bit of extra accuracy compared to
    // doing the iterations in all-double.
    setUlpTol(4*ulpTol_);

    GMX_EXPECT_SIMD_FUNC_NEAR(ref_invsqrt, tst_invsqrt_pair0);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_invsqrt, tst_invsqrt_pair1);
}

TEST_F(SimdMathTest, gmxSimdSqrtR)
{
    // Just make sure sqrt(0)=0 works and isn't evaluated as 0*1/sqrt(0)=NaN
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdFrom3R(0, 2, 3), gmx_simd_sqrt_r(setSimdFrom3R(0, 4, 9)));
}

// Helper function wrapper
real ref_inv(real x)
{
    return 1.0/x;
}

TEST_F(SimdMathTest, gmxSimdInvR)
{
    // test <0
    setRange(-1e10, -1e-10);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_inv, gmx_simd_inv_r);
    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_inv, gmx_simd_inv_r);
}

// Helper function wrapper
real ref_log(real x)
{
    return log(x);
}

TEST_F(SimdMathTest, gmxSimdLogR)
{
    setRange(1e-30, 1e30);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_log, gmx_simd_log_r);
}

// MSVC does not support exp2(), so we have no reference to test against
#ifndef _MSC_VER
// Helper function wrapper
real ref_exp2(real x)
{
    return exp2(x);
}

TEST_F(SimdMathTest, gmxSimdExp2R)
{
    setRange(-100, 100);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_exp2, gmx_simd_exp2_r);
}
#endif

// Helper function wrapper
real ref_exp(real x)
{
    return exp(x);
}

TEST_F(SimdMathTest, gmxSimdExpR)
{
    setRange(-75, 75);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_exp, gmx_simd_exp_r);
}

// Helper function wrapper. The single-precision gmx_erff() in gmxlib is
// slightly lower precision than the SIMD flavor, so use double for reference.
real ref_erf(real x)
{
    return gmx_erfd(x);
}

TEST_F(SimdMathTest, gmxSimdErfR)
{
    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_erf, gmx_simd_erf_r);
}

// Helper function wrapper. The single-precision gmx_erfcf() in gmxlib is
// slightly lower precision than the SIMD flavor, so use double for reference.
real ref_erfc(real x)
{
    return gmx_erfcd(x);
}

TEST_F(SimdMathTest, gmxSimdErfcR)
{
    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    // Our erfc algorithm has 4 ulp accuracy, so relax defaultTol a bit
    setUlpTol(4*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_erfc, gmx_simd_erfc_r);
}

// Helper function wrapper
real ref_sin(real x)
{
    return sin(x);
}

TEST_F(SimdMathTest, gmxSimdSinR)
{
    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_sin, gmx_simd_sin_r);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    setUlpTol(2*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_sin, gmx_simd_sin_r);
}

// Helper function wrapper
real ref_cos(real x)
{
    return cos(x);
}

TEST_F(SimdMathTest, gmxSimdCosR)
{
    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_cos, gmx_simd_cos_r);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    setUlpTol(2*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_cos, gmx_simd_cos_r);
}

// Helper function wrapper
real ref_tan(real x)
{
    return tan(x);
}

TEST_F(SimdMathTest, gmxSimdTanR)
{
    // Tan(x) is a little sensitive due to the division in the algorithm.
    // Rather than using lots of extra FP operations, we accept the algorithm
    // presently only achieves a ~3 ulp error and use the medium tolerance.
    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_tan, gmx_simd_tan_r);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    setUlpTol(2*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_tan, gmx_simd_tan_r);
}

// Helper function wrapper
real ref_asin(real x)
{
    return asin(x);
}

TEST_F(SimdMathTest, gmxSimdAsinR)
{
    // Our present asin(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_asin, gmx_simd_asin_r);
}

// Helper function wrapper
real ref_acos(real x)
{
    return acos(x);
}

TEST_F(SimdMathTest, gmxSimdAcosR)
{
    // Our present acos(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_acos, gmx_simd_acos_r);
}

// Helper function wrapper
real ref_atan(real x)
{
    return atan(x);
}

TEST_F(SimdMathTest, gmxSimdAtanR)
{
    // Our present atan(x) algorithm achieves 1 ulp accuracy
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_atan, gmx_simd_atan_r);
}

TEST_F(SimdMathTest, gmxSimdAtan2R)
{
    // test each quadrant
    GMX_EXPECT_SIMD_REAL_NEAR(atan2(1.0, 1.0), gmx_simd_atan2_r(rSimd_1_2_3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(atan2(-1.0, 1.0), gmx_simd_atan2_r(rSimd_m1_m2_m3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(atan2(-1.0, -1.0), gmx_simd_atan2_r(rSimd_m1_m2_m3, rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(atan2(1.0, -1.0), gmx_simd_atan2_r(rSimd_1_2_3, rSimd_m1_m2_m3));
    // cases important for calculating angles
    // values on coordinate axes
    GMX_EXPECT_SIMD_REAL_NEAR(atan2(0.0, 1.0), gmx_simd_atan2_r(gmx_simd_setzero_r(), rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(atan2(1.0, 0.0), gmx_simd_atan2_r(rSimd_1_2_3, gmx_simd_setzero_r()));
    GMX_EXPECT_SIMD_REAL_NEAR(atan2(0.0, -1.0), gmx_simd_atan2_r(gmx_simd_setzero_r(), rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(atan2(-1.0, 0.0), gmx_simd_atan2_r(rSimd_m1_m2_m3, gmx_simd_setzero_r()));
    // degenerate value (origin) should return 0.0
    GMX_EXPECT_SIMD_REAL_NEAR(atan2(0.0, 0.0), gmx_simd_atan2_r(setSimdFrom3R(0.0, 0.0, 0.0), gmx_simd_setzero_r()));
}

// Helper function wrapper
real ref_pmecorrF(real x)
{
    real y = sqrt(x);
    return 2*exp(-x)/(sqrt(M_PI)*x) - gmx_erfd(y)/(x*y);
}

// The PME corrections will be added to ~1/r2, so absolute tolerance of EPS is fine.
TEST_F(SimdMathTest, gmxSimdPmecorrForceR)
{
    // Pme correction only needs to be ~1e-6 accuracy single, 1e-10 double
#ifdef GMX_DOUBLE
    setUlpTol((gmx_int64_t)(1e-10/GMX_REAL_EPS));
#else
    setUlpTol((gmx_int64_t)(1e-6/GMX_REAL_EPS));
#endif

    setRange(0.15, 4);
    setAbsTol(GMX_REAL_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_pmecorrF, gmx_simd_pmecorrF_r);
}

// Helper function wrapper
real ref_pmecorrV(real x)
{
    real y = sqrt(x);
    return gmx_erfd(y)/y;
}

// The PME corrections will be added to ~1/r, so absolute tolerance of EPS is fine.
TEST_F(SimdMathTest, gmxSimdPmecorrPotentialR)
{
    // Pme correction only needs to be ~1e-6 accuracy single, 1e-10 double
#ifdef GMX_DOUBLE
    setUlpTol((gmx_int64_t)(1e-10/GMX_REAL_EPS));
#else
    setUlpTol((gmx_int64_t)(1e-6/GMX_REAL_EPS));
#endif
    setRange(0.15, 4);
    setAbsTol(GMX_REAL_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_pmecorrV, gmx_simd_pmecorrV_r);
}

}      // namespace

#endif // GMX_SIMD_HAVE_REAL

}      // namespace
}      // namespace
