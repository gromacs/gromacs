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

#include "gromacs/simd/simd_math.h"

#include "config.h"

#include <vector>

#include "gromacs/math/utilities.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/simd/simd.h"

#include "simd.h"

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#ifdef GMX_SIMD_HAVE_REAL

class SimdMathTest : public SimdTest
{
    public:
        ::testing::AssertionResult
                            compareSimdMathFunction(const char * refFuncExpr, const char *simdFuncExpr,
                                                    real refFunc(real x),     gmx_simd_real_t gmx_simdcall simdFunc(gmx_simd_real_t x));
};

/*! \brief Test approximate equality of SIMD vs reference version of a function.
 *
 * This macro takes vanilla C and SIMD flavors of a function and tests it with
 * the number of points, range, and tolerances specified by the test fixture class.
 */
#define GMX_EXPECT_SIMD_FUNC_NEAR(refFunc, tstFunc) \
    EXPECT_PRED_FORMAT2(compareSimdMathFunction, refFunc, tstFunc)

/*! \brief Implementation routine to compare SIMD vs reference functions.
 *
 * \param refFuncExpr   Description of reference function expression
 * \param simdFuncExpr  Description of SIMD function expression
 * \param refFunc       Reference math function pointer
 * \param simdFunc      SIMD math function pointer
 *
 * The function will be tested with the range and tolerances specified in
 * the SimdBaseTest class. You should not never call this function directly,
 * but use the macro GMX_EXPECT_SIMD_FUNC_NEAR(refFunc,tstFunc) instead.
 */
::testing::AssertionResult
SimdMathTest::compareSimdMathFunction(const char * refFuncExpr, const char *simdFuncExpr,
                                      real refFunc(real x),     gmx_simd_real_t gmx_simdcall simdFunc(gmx_simd_real_t x))
{
    std::vector<real>            vx(GMX_SIMD_REAL_WIDTH);
    std::vector<real>            vref(GMX_SIMD_REAL_WIDTH);
    std::vector<real>            vtst(GMX_SIMD_REAL_WIDTH);
    real                         dx, absDiff;
    gmx_int64_t                  ulpDiff, maxUlpDiff;
    real                         maxUlpDiffPos;
    real                         refValMaxUlpDiff, simdValMaxUlpDiff;
    bool                         absOk, signOk;
    int                          i, iter;
    int                          niter   = s_nPoints/GMX_SIMD_REAL_WIDTH;
    int                          npoints = niter*GMX_SIMD_REAL_WIDTH;
#    ifdef GMX_DOUBLE
    union {
        double r; gmx_int64_t i;
    } conv0, conv1;
#    else
    union {
        float  r; gmx_int32_t i;
    } conv0, conv1;
#    endif

    maxUlpDiff = 0;
    dx         = (range_.second-range_.first)/npoints;

    for (iter = 0; iter < niter; iter++)
    {
        for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
        {
            vx[i]   = range_.first+dx*(iter*GMX_SIMD_REAL_WIDTH+i);
            vref[i] = refFunc(vx[i]);
        }
        vtst  = simdReal2Vector(simdFunc(vector2SimdReal(vx)));

        for (i = 0, signOk = true, absOk = true; i < GMX_SIMD_REAL_WIDTH; i++)
        {
            absDiff = fabs(vref[i]-vtst[i]);
            absOk   = absOk  && ( absDiff < absTol_ );
            signOk  = signOk && ( (vref[i] >= 0 && vtst[i] >= 0) ||
                                  (vref[i] <= 0 && vtst[i] <= 0));

            if (absDiff >= absTol_)
            {
                /* We replicate the trivial ulp differences comparison here rather than
                 * calling the lower-level routine for comparing them, since this enables
                 * us to run through the entire test range and report the largest deviation
                 * without lots of extra glue routines.
                 */
                conv0.r           = vref[i];
                conv1.r           = vtst[i];
                ulpDiff           = llabs(conv0.i-conv1.i);
                if (ulpDiff > maxUlpDiff)
                {
                    maxUlpDiff        = ulpDiff;
                    maxUlpDiffPos     = vx[i];
                    refValMaxUlpDiff  = vref[i];
                    simdValMaxUlpDiff = vtst[i];
                }
            }
        }
        if ( (absOk == false) && (signOk == false) )
        {
            return ::testing::AssertionFailure()
                   << "Failing SIMD math function comparison due to sign differences." << std::endl
                   << "Reference function: " << refFuncExpr << std::endl
                   << "Simd function:      " << simdFuncExpr << std::endl
                   << "Test range is ( " << range_.first << " , " << range_.second << " ) " << std::endl
                   << "First sign difference around x=" << std::setprecision(20) << ::testing::PrintToString(vx) << std::endl
                   << "Ref values:  " << std::setprecision(20) << ::testing::PrintToString(vref) << std::endl
                   << "SIMD values: " << std::setprecision(20) << ::testing::PrintToString(vtst) << std::endl;
        }
    }

    if (maxUlpDiff <= ulpTol_)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing SIMD math function ulp comparison between " << refFuncExpr << " and " << simdFuncExpr << std::endl
               << "Requested ulp tolerance: " << ulpTol_ << std::endl
               << "Requested abs tolerance: " << absTol_ << std::endl
               << "Largest Ulp difference occurs for x=" << std::setprecision(20) << maxUlpDiffPos << std::endl
               << "Ref  values: " << std::setprecision(20) << refValMaxUlpDiff << std::endl
               << "SIMD values: " << std::setprecision(20) << simdValMaxUlpDiff << std::endl
               << "Ulp diff.:   " << std::setprecision(20) << maxUlpDiff << std::endl;
    }
}

/*! \} */
/*! \endcond */


// Actual math function tests below


namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

TEST_F(SimdMathTest, gmxSimdXorSignR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-4, 5, 6), gmx_simd_xor_sign_r(setSimdRealFrom3R(4, 5, 6), setSimdRealFrom3R(-5, 2, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(4, -5, -6), gmx_simd_xor_sign_r(setSimdRealFrom3R(-4, -5, -6), setSimdRealFrom3R(-5, 2, 0)));
}

/*! \brief Function wrapper to evaluate reference 1/sqrt(x) */
static real
ref_invsqrt(real x)
{
    return 1.0/sqrt(x);
}

TEST_F(SimdMathTest, gmxSimdInvsqrtR)
{
    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_invsqrt, gmx_simd_invsqrt_r);
}

/*! \brief Function wrapper to return first result when testing \ref gmx_simd_invsqrt_pair_r */
gmx_simd_real_t gmx_simdcall
tst_invsqrt_pair0(gmx_simd_real_t x)
{
    gmx_simd_real_t r0, r1;
    gmx_simd_invsqrt_pair_r(x, x, &r0, &r1);
    return r0;
}

/*! \brief Function wrapper to return second result when testing \ref gmx_simd_invsqrt_pair_r */
gmx_simd_real_t gmx_simdcall
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
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(0, 2, 3), gmx_simd_sqrt_r(setSimdRealFrom3R(0, 4, 9)));
}

/*! \brief Function wrapper to evaluate reference 1/x */
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

/*! \brief Function wrapper for log(x), with argument/return in default Gromacs precision */
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
/*! \brief Function wrapper for exp2(x), with argument/return in default Gromacs precision */
real ref_exp2(real x)
{
    return exp2(x);
}

TEST_F(SimdMathTest, gmxSimdExp2R)
{
    setRange(-100, 100);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_exp2, gmx_simd_exp2_r);

    // We do not care about the SIMD implementation getting denormal values right,
    // but they must be clamped to zero rather than producing garbage.
    // Check by setting the absolute tolerance to machine precision.
    setAbsTol(GMX_REAL_EPS);

    // First two values will have denormal results in single, third value in double too.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(ref_exp2(-150.0), ref_exp2(-300.0), ref_exp2(-1050.0)),
                              gmx_simd_exp2_r(setSimdRealFrom3R(-150.0, -300.0, -1050.0)));

    // Reset absolute tolerance to enforce ULP checking
    setAbsTol(0.0);

    // Make sure that underflowing values are set to zero.
    // First two values underflow in single, third value in double too.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(ref_exp2(-200.0), ref_exp2(-600.0), ref_exp2(-1500.0)),
                              gmx_simd_exp2_r(setSimdRealFrom3R(-200.0, -600.0, -1500.0)));
}
#endif

/*! \brief Function wrapper for exp(x), with argument/return in default Gromacs precision */
real ref_exp(real x)
{
    return exp(x);
}

TEST_F(SimdMathTest, gmxSimdExpR)
{
    setRange(-75, 75);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_exp, gmx_simd_exp_r);

    // We do not care about the SIMD implementation getting denormal values right,
    // but they must be clamped to zero rather than producing garbage.
    // Check by setting the absolute tolerance to machine precision.
    setAbsTol(GMX_REAL_EPS);
    // First two values will have denormal results in single, third value in double too.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(ref_exp(-90.0), ref_exp(-100.0), ref_exp(-725.0)),
                              gmx_simd_exp_r(setSimdRealFrom3R(-90.0, -100.0, -725.0)));

    // Reset absolute tolerance to enforce ULP checking
    setAbsTol(0.0);

    // Make sure that underflowing values are set to zero.
    // First two values underflow in single, third value in double too.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(ref_exp(-150.0), ref_exp(-300.0), ref_exp(-800.0)),
                              gmx_simd_exp_r(setSimdRealFrom3R(-150.0, -300.0, -800.0)));
}

/*! \brief Function wrapper for erf(x), with argument/return in default Gromacs precision.
 *
 * \note The single-precision gmx_erff() in gmxlib is slightly lower precision
 * than the SIMD flavor, so we use double for reference.
 */
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

/*! \brief Function wrapper for erfc(x), with argument/return in default Gromacs precision.
 *
 * \note The single-precision gmx_erfcf() in gmxlib is slightly lower precision
 * than the SIMD flavor, so we use double for reference.
 */
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

/*! \brief Function wrapper for sin(x), with argument/return in default Gromacs precision */
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

/*! \brief Function wrapper for cos(x), with argument/return in default Gromacs precision */
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

/*! \brief Function wrapper for tan(x), with argument/return in default Gromacs precision */
real ref_tan(real x)
{
#if (defined(__ibmxl__) || defined(__xlC__)) && !defined(GMX_DOUBLE)
    /* xlC (both version 12 and 13) has some strange behaviour where tan(x) with a float argument incorrectly
     * returns -0.0 for the argument 0.0 at -O3 optimization. tanf() seems to avoid it.
     */
    return tanf(x);
#else
    return tan(x);
#endif
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

/*! \brief Function wrapper for asin(x), with argument/return in default Gromacs precision */
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

/*! \brief Function wrapper for acos(x), with argument/return in default Gromacs precision */
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

/*! \brief Function wrapper for atan(x), with argument/return in default Gromacs precision */
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
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(1.0, 1.0)), gmx_simd_atan2_r(rSimd_1_2_3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(-1.0, 1.0)), gmx_simd_atan2_r(rSimd_m1_m2_m3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(-1.0, -1.0)), gmx_simd_atan2_r(rSimd_m1_m2_m3, rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(1.0, -1.0)), gmx_simd_atan2_r(rSimd_1_2_3, rSimd_m1_m2_m3));
    // cases important for calculating angles
    // values on coordinate axes
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(0.0, 1.0)), gmx_simd_atan2_r(gmx_simd_setzero_r(), rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(1.0, 0.0)), gmx_simd_atan2_r(rSimd_1_2_3, gmx_simd_setzero_r()));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(0.0, -1.0)), gmx_simd_atan2_r(gmx_simd_setzero_r(), rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(-1.0, 0.0)), gmx_simd_atan2_r(rSimd_m1_m2_m3, gmx_simd_setzero_r()));
    // degenerate value (origin) should return 0.0
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(0.0, 0.0)), gmx_simd_atan2_r(setSimdRealFrom3R(0.0, 0.0, 0.0), gmx_simd_setzero_r()));
}

/*! \brief Evaluate reference version of PME force correction. */
real ref_pmecorrF(real x)
{
    if (x != 0)
    {
        real y = sqrt(x);
        return 2*exp(-x)/(sqrt(M_PI)*x) - gmx_erfd(y)/(x*y);
    }
    else
    {
        return -4/(3*sqrt(M_PI));
    }
}

// The PME corrections will be added to ~1/r2, so absolute tolerance of EPS is fine.
TEST_F(SimdMathTest, gmxSimdPmecorrForceR)
{
    // Pme correction only needs to be ~1e-6 accuracy single, 1e-10 double
#ifdef GMX_DOUBLE
    setUlpTol((gmx_int64_t)(5e-10/GMX_REAL_EPS));
#else
    setUlpTol((gmx_int64_t)(5e-6/GMX_REAL_EPS));
#endif

    setRange(0.15, 4);
    setAbsTol(GMX_REAL_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_pmecorrF, gmx_simd_pmecorrF_r);
}

/*! \brief Evaluate reference version of PME potential correction. */
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
    setUlpTol((gmx_int64_t)(5e-10/GMX_REAL_EPS));
#else
    setUlpTol((gmx_int64_t)(5e-6/GMX_REAL_EPS));
#endif
    setRange(0.15, 4);
    setAbsTol(GMX_REAL_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_pmecorrV, gmx_simd_pmecorrV_r);
}






/*
 * Functions that only target single accuracy, even for double SIMD data
 */

TEST_F(SimdMathTest, gmxSimdInvsqrtSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_invsqrt, gmx_simd_invsqrt_singleaccuracy_r);
}

/*! \brief Function wrapper to return first result when testing \ref gmx_simd_invsqrt_pair_singleaccuracy_r */
gmx_simd_real_t gmx_simdcall
tst_invsqrt_singleaccuracy_pair0(gmx_simd_real_t x)
{
    gmx_simd_real_t r0, r1;
    gmx_simd_invsqrt_pair_singleaccuracy_r(x, x, &r0, &r1);
    return r0;
}

/*! \brief Function wrapper to return second result when testing \ref gmx_simd_invsqrt_pair_singleaccuracy_r */
gmx_simd_real_t gmx_simdcall
tst_invsqrt_singleaccuracy_pair1(gmx_simd_real_t x)
{
    gmx_simd_real_t r0, r1;
    gmx_simd_invsqrt_pair_singleaccuracy_r(x, x, &r0, &r1);
    return r1;
}

TEST_F(SimdMathTest, gmxSimdInvsqrtPairSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_invsqrt, tst_invsqrt_singleaccuracy_pair0);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_invsqrt, tst_invsqrt_singleaccuracy_pair1);
}

TEST_F(SimdMathTest, gmxSimdSqrtSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Just make sure sqrt(0)=0 works and isn't evaluated as 0*1/sqrt(0)=NaN
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(0, 2, 3), gmx_simd_sqrt_singleaccuracy_r(setSimdRealFrom3R(0, 4, 9)));
}

TEST_F(SimdMathTest, gmxSimdInvSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // test <0
    setRange(-1e10, -1e-10);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_inv, gmx_simd_inv_singleaccuracy_r);
    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_inv, gmx_simd_inv_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdLogSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(1e-30, 1e30);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_log, gmx_simd_log_singleaccuracy_r);
}

// MSVC does not support exp2(), so we have no reference to test against
#ifndef _MSC_VER
TEST_F(SimdMathTest, gmxSimdExp2SingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-100, 100);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_exp2, gmx_simd_exp2_singleaccuracy_r);
}
#endif

TEST_F(SimdMathTest, gmxSimdExpSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-75, 75);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_exp, gmx_simd_exp_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdErfSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_erf, gmx_simd_erf_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdErfcSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    // Our erfc algorithm has 4 ulp accuracy, so relax tolerance a bit
    setUlpTol(4*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_erfc, gmx_simd_erfc_singleaccuracy_r);
}


TEST_F(SimdMathTest, gmxSimdSinSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_sin, gmx_simd_sin_singleaccuracy_r);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    setUlpTol(2*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_sin, gmx_simd_sin_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdCosSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_cos, gmx_simd_cos_singleaccuracy_r);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_cos, gmx_simd_cos_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdTanSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Tan(x) is a little sensitive due to the division in the algorithm.
    // Rather than using lots of extra FP operations, we accept the algorithm
    // presently only achieves a ~3 ulp error and use the medium tolerance.
    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_tan, gmx_simd_tan_singleaccuracy_r);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_tan, gmx_simd_tan_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdAsinSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Our present asin(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_asin, gmx_simd_asin_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdAcosSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Our present acos(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_acos, gmx_simd_acos_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdAtanSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Our present atan(x) algorithm achieves 1 ulp accuracy
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_atan, gmx_simd_atan_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdAtan2SingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // test each quadrant
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(1.0, 1.0)), gmx_simd_atan2_singleaccuracy_r(rSimd_1_2_3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(-1.0, 1.0)), gmx_simd_atan2_singleaccuracy_r(rSimd_m1_m2_m3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(-1.0, -1.0)), gmx_simd_atan2_singleaccuracy_r(rSimd_m1_m2_m3, rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(1.0, -1.0)), gmx_simd_atan2_singleaccuracy_r(rSimd_1_2_3, rSimd_m1_m2_m3));
    // cases important for calculating angles
    // values on coordinate axes
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(0.0, 1.0)), gmx_simd_atan2_singleaccuracy_r(gmx_simd_setzero_r(), rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(1.0, 0.0)), gmx_simd_atan2_singleaccuracy_r(rSimd_1_2_3, gmx_simd_setzero_r()));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(0.0, -1.0)), gmx_simd_atan2_singleaccuracy_r(gmx_simd_setzero_r(), rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(-1.0, 0.0)), gmx_simd_atan2_singleaccuracy_r(rSimd_m1_m2_m3, gmx_simd_setzero_r()));
    // degenerate value (origin) should return 0.0
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(0.0, 0.0)), gmx_simd_atan2_singleaccuracy_r(setSimdRealFrom3R(0.0, 0.0, 0.0), gmx_simd_setzero_r()));
}

TEST_F(SimdMathTest, gmxSimdPmecorrForceSingleaccuracyR)
{
    // The PME corrections will be added to ~1/r2, so absolute tolerance of EPS is fine.
    // Pme correction only needs to be ~1e-6 accuracy single.
    // Then increase the allowed error by the difference between the actual precision and single.
    setUlpTol( ((gmx_int64_t)(5e-6/GMX_FLOAT_EPS)) * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(0.15, 4);
    setAbsTol(GMX_FLOAT_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_pmecorrF, gmx_simd_pmecorrF_singleaccuracy_r);
}

TEST_F(SimdMathTest, gmxSimdPmecorrPotentialSingleaccuracyR)
{
    // The PME corrections will be added to ~1/r, so absolute tolerance of EPS is fine.
    // Pme correction only needs to be ~1e-6 accuracy single.
    // Then increase the allowed error by the difference between the actual precision and single.
    setUlpTol( ((gmx_int64_t)(5e-6/GMX_FLOAT_EPS)) * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(0.15, 4);
    setAbsTol(GMX_FLOAT_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(ref_pmecorrV, gmx_simd_pmecorrV_singleaccuracy_r);
}

}      // namespace

#endif // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
