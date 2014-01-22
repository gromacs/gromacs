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

#include <vector>

#include "util.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/math/utilities.h"

namespace simdTest
{

#ifdef GMX_SIMD_HAVE_REAL

TEST(SimdTestMath, gmxSimdXorSignR)
{
    GMX_ASSERT_SIMD_REAL_EQ(setSimd3R(-4, 5, 6), gmx_simd_xor_sign_r(setSimd3R(4, 5, 6), setSimd3R(-5, 2, 0)));
    GMX_ASSERT_SIMD_REAL_EQ(setSimd3R(4, -5, -6), gmx_simd_xor_sign_r(setSimd3R(-4, -5, -6), setSimd3R(-5, 2, 0)));
}

// Helper function wrapper
real ref_invsqrt(real x)
{
    return 1.0/sqrt(x);
}

TEST(SimdTestMath, gmxSimdInvsqrtR)
{
    std::pair<real, real>  range(1e-10, 1e10);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_invsqrt, gmx_simd_invsqrt_r, range, defaultTol, 0);
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

TEST(SimdTestMath, gmxSimdInvsqrtPairR)
{
    std::pair<real, real>  range(1e-10, 1e10);
    /* The accuracy conversions lose a bit of extra accuracy compared to
     * doing the iterations in all-double.
     */
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_invsqrt, tst_invsqrt_pair0, range, 4*defaultTol, 0);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_invsqrt, tst_invsqrt_pair1, range, 4*defaultTol, 0);
}

TEST(SimdTestMath, gmxSimdSqrtR)
{
    // Just make sure sqrt(0)=0 works and isn't evaluated as 0*1/sqrt(0)=NaN
    GMX_ASSERT_SIMD_REAL_NEAR(setSimd3R(0, 2, 3), gmx_simd_sqrt_r(setSimd3R(0, 4, 9)), defaultTol, 0);
}

// Helper function wrapper
real ref_inv(real x)
{
    return 1.0/x;
}

TEST(SimdTestMath, gmxSimdInvR)
{
    // test <0
    std::pair<real, real>  range(-1e10, -1e-10);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_inv, gmx_simd_inv_r, range, defaultTol, 0);
    // test >0
    std::pair<real, real>  range2(1e-10, 1e10);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_inv, gmx_simd_inv_r, range2, defaultTol, 0);
}

// Helper function wrapper
real ref_log(real x)
{
    return log(x);
}

TEST(SimdTestMath, gmxSimdLogR)
{
    std::pair<real, real>  range(1e-30, 1e30);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_log, gmx_simd_log_r, range, defaultTol, 0);
}

// MSVC does not support exp2(), so we have no reference to test against
#ifndef _MSC_VER
// Helper function wrapper
real ref_exp2(real x)
{
    return exp2(x);
}

TEST(SimdTestMath, gmxSimdExp2R)
{
    std::pair<real, real>  range(-100, 100);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_exp2, gmx_simd_exp2_r, range, defaultTol, 0);
}
#endif

// Helper function wrapper
real ref_exp(real x)
{
    return exp(x);
}

TEST(SimdTestMath, gmxSimdExpR)
{
    std::pair<real, real>  range(-75, 75);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_exp, gmx_simd_exp_r, range, defaultTol, 0);
}

// Helper function wrapper. The single-precision gmx_erff() in gmxlib is
// slightly lower precision than the SIMD flavor, so use double for reference.
real ref_erf(real x)
{
    return gmx_erfd(x);
}

TEST(SimdTestMath, gmxSimdErfR)
{
    std::pair<real, real>  range(-9, 9);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_erf, gmx_simd_erf_r, range, defaultTol, GMX_REAL_MIN);
}

// Helper function wrapper. The single-precision gmx_erfcf() in gmxlib is
// slightly lower precision than the SIMD flavor, so use double for reference.
real ref_erfc(real x)
{
    return gmx_erfcd(x);
}

TEST(SimdTestMath, gmxSimdErfcR)
{
    // Our erfc algorithm has 4 ulp accuracy, so relax defaultTol a bit
    std::pair<real, real>  range(-9, 9);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_erfc, gmx_simd_erfc_r, range, 4+defaultTol, GMX_REAL_MIN);
}

// Helper function wrapper
real ref_sin(real x)
{
    return sin(x);
}

TEST(SimdTestMath, gmxSimdSinR)
{
    std::pair<real, real>  range(-8*M_PI, 8*M_PI);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_sin, gmx_simd_sin_r, range, defaultTol, 0);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    std::pair<real, real>  range2(-10000, 10000);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_sin, gmx_simd_sin_r, range2, 2*defaultTol, 0);
}

// Helper function wrapper
real ref_cos(real x)
{
    return cos(x);
}

TEST(SimdTestMath, gmxSimdCosR)
{
    std::pair<real, real>  range(-8*M_PI, 8*M_PI);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_cos, gmx_simd_cos_r, range, defaultTol, 0);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    std::pair<real, real>  range2(-10000, 10000);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_cos, gmx_simd_cos_r, range2, 2*defaultTol, 0);
}

// Helper function wrapper
real ref_tan(real x)
{
    return tan(x);
}

TEST(SimdTestMath, gmxSimdTanR)
{
    // Tan(x) is a little sensitive due to the division in the algorithm.
    // Rather than using lots of extra FP operations, we accept the algorithm
    // presently only achieves a ~3 ulp error and use the medium tolerance.
    std::pair<real, real>  range(-8*M_PI, 8*M_PI);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_tan, gmx_simd_tan_r, range, defaultTol, 0);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    std::pair<real, real>  range2(-10000, 10000);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_tan, gmx_simd_tan_r, range2, 2*defaultTol, 0);
}

// Helper function wrapper
real ref_asin(real x)
{
    return asin(x);
}

TEST(SimdTestMath, gmxSimdAsinR)
{
    // Our present asin(x) algorithm achieves 2-3 ulp accuracy
    std::pair<real, real>  range(-1, 1);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_asin, gmx_simd_asin_r, range, defaultTol, 0);
}

// Helper function wrapper
real ref_acos(real x)
{
    return acos(x);
}

TEST(SimdTestMath, gmxSimdAcosR)
{
    // Our present acos(x) algorithm achieves 2-3 ulp accuracy
    std::pair<real, real>  range(-1, 1);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_acos, gmx_simd_acos_r, range, defaultTol, 0);
}

// Helper function wrapper
real ref_atan(real x)
{
    return atan(x);
}

TEST(SimdTestMath, gmxSimdAtanR)
{
    // Our present atan(x) algorithm achieves 1 ulp accuracy
    std::pair<real, real>  range(-10000, 10000);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_atan, gmx_simd_atan_r, range, defaultTol, 0);
}

TEST(SimdTestMath, gmxSimdAtan2R)
{
    // Just test the sign
    gmx_simd_real_t       zero  = gmx_simd_setzero_r();
    gmx_simd_real_t       mzero = gmx_simd_set1_r(-0.0);

    GMX_ASSERT_SIMD_REAL_NEAR(atan2(1.0, 1.0), gmx_simd_atan2_r(rSimd_1_2_3, rSimd_1_2_3), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(-1.0, 1.0), gmx_simd_atan2_r(rSimd_m1_m2_m3, rSimd_1_2_3), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(-1.0, -1.0), gmx_simd_atan2_r(rSimd_m1_m2_m3, rSimd_m1_m2_m3), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(1.0, -1.0), gmx_simd_atan2_r(rSimd_1_2_3, rSimd_m1_m2_m3), defaultTol, 0);
    // special values should be same as for iso standard (see atan2 man page)
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(0.0, -0.0), gmx_simd_atan2_r(zero, mzero), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(-0.0, -0.0), gmx_simd_atan2_r(mzero, mzero), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(0.0, 0.0), gmx_simd_atan2_r(zero, zero), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(-0.0, 0.0), gmx_simd_atan2_r(mzero, zero), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(0.0, -1.0), gmx_simd_atan2_r(zero, rSimd_m1_m2_m3), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(-0.0, -1.0), gmx_simd_atan2_r(mzero, rSimd_m1_m2_m3), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(0.0, -1.0), gmx_simd_atan2_r(zero, rSimd_m1_m2_m3), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(-0.0, 1.0), gmx_simd_atan2_r(mzero, rSimd_1_2_3), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(1.0, 0.0), gmx_simd_atan2_r(rSimd_1_2_3, zero), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(1.0, -0.0), gmx_simd_atan2_r(rSimd_1_2_3, mzero), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(-1.0, 0.0), gmx_simd_atan2_r(rSimd_m1_m2_m3, zero), defaultTol, 0);
    GMX_ASSERT_SIMD_REAL_NEAR(atan2(-1.0, -0.0), gmx_simd_atan2_r(rSimd_m1_m2_m3, mzero), defaultTol, 0);
}

// Helper function wrapper
real ref_pmecorrF(real x)
{
    real y = sqrt(x);
    return 2*exp(-x)/(sqrt(M_PI)*x) - gmx_erfd(y)/(x*y);
}

// The PME corrections will be added to ~1/r2, so absolute tolerance of EPS is fine.
TEST(SimdTestMath, gmxSimdPmecorrForceR)
{
    /* Pme correction only needs to be ~1e-6 accuracy single, 1e-10 double */
#ifdef GMX_DOUBLE
    int                    ulpTol = (int)(1e-10/GMX_REAL_EPS);
#else
    int                    ulpTol = (int)(1e-6/GMX_REAL_EPS);
#endif
    real                   absTol = GMX_REAL_EPS;
    std::pair<real, real>  range(0.15, 3);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_pmecorrF, gmx_simd_pmecorrF_r, range, ulpTol, absTol);
}

// Helper function wrapper
real ref_pmecorrV(real x)
{
    real y = sqrt(x);
    return gmx_erfd(y)/y;
}

// The PME corrections will be added to ~1/r, so absolute tolerance of EPS is fine.
TEST(SimdTestMath, gmxSimdPmecorrPotentialR)
{
    /* Pme correction only needs to be ~1e-6 accuracy single, 1e-10 double */
#ifdef GMX_DOUBLE
    int                    ulpTol = (int)(1e-10/GMX_REAL_EPS);
#else
    int                    ulpTol = (int)(1e-6/GMX_REAL_EPS);
#endif
    real                   absTol = GMX_REAL_EPS;
    std::pair<real, real>  range(0.15, 4);
    GMX_ASSERT_SIMD_FUNC_NEAR(ref_pmecorrV, gmx_simd_pmecorrV_r, range, ulpTol, absTol);
}

#endif

}
