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

#include <cmath>
#include <cstdint>

#include <vector>

#include "gromacs/options/basicoptions.h"
#include "gromacs/simd/simd.h"

#include "simd.h"

#if GMX_SIMD

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#if GMX_SIMD_HAVE_REAL

class SimdMathTest : public SimdTest
{
    public:
        ::testing::AssertionResult
                            compareSimdMathFunction(const char * refFuncExpr, const char *simdFuncExpr,
                                                    real refFunc(real x),     SimdReal gmx_simdcall simdFunc(SimdReal x));
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
                                      real refFunc(real x),     SimdReal gmx_simdcall simdFunc(SimdReal x))
{
    std::vector<real>            vx(GMX_SIMD_REAL_WIDTH);
    std::vector<real>            vref(GMX_SIMD_REAL_WIDTH);
    std::vector<real>            vtst(GMX_SIMD_REAL_WIDTH);
    real                         dx, absDiff;
    std::int64_t                 ulpDiff, maxUlpDiff;
    real                         maxUlpDiffPos;
    real                         refValMaxUlpDiff, simdValMaxUlpDiff;
    bool                         absOk, signOk;
    int                          i, iter;
    int                          niter   = s_nPoints/GMX_SIMD_REAL_WIDTH;
    int                          npoints = niter*GMX_SIMD_REAL_WIDTH;
#    ifdef GMX_DOUBLE
    union {
        double r; std::int64_t i;
    } conv0, conv1;
#    else
    union {
        float  r; std::int32_t i;
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

TEST_F(SimdMathTest, copysign)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-4, 5, 6), copysign(setSimdRealFrom3R(4, 5, 6), setSimdRealFrom3R(-5, 2, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-4, 5, 6), copysign(setSimdRealFrom3R(-4, -5, -6), setSimdRealFrom3R(-5, 2, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(4, -5, 6), copysign(setSimdRealFrom3R(4, 5, 6), setSimdRealFrom3R(5, -2, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(4, -5, 6), copysign(setSimdRealFrom3R(-4, -5, -6), setSimdRealFrom3R(5, -2, 0)));
}

/*! \brief Function wrapper to evaluate reference 1/sqrt(x) */
static real
refInvsqrt(real x)
{
    return 1.0/std::sqrt(x);
}

TEST_F(SimdMathTest, invsqrt)
{
    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, invsqrt);
}

TEST_F(SimdMathTest, invsqrtMask)
{
    SimdReal x   = setSimdRealFrom3R(1.0, 0.0, 3.0);
    SimdBool m   = (SimdReal(0.0) < x);
    SimdReal ref = setSimdRealFrom3R(1.0/std::sqrt(1.0), 0.0, 1.0/std::sqrt(3.0));
    GMX_EXPECT_SIMD_REAL_NEAR(ref, invsqrtMask(x, m));
}

/*! \brief Function wrapper to return first result when testing \ref simdInvsqrtPair */
SimdReal gmx_simdcall
tstInvsqrtPair0(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPair(x, x, &r0, &r1);
    return r0;
}

/*! \brief Function wrapper to return second result when testing \ref simdInvsqrtPair */
SimdReal gmx_simdcall
tstInvsqrtPair1(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPair(x, x, &r0, &r1);
    return r1;
}

TEST_F(SimdMathTest, invsqrtPair)
{
    setRange(1e-10, 1e10);
    // The accuracy conversions lose a bit of extra accuracy compared to
    // doing the iterations in all-double.
    setUlpTol(4*ulpTol_);

    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tstInvsqrtPair0);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tstInvsqrtPair1);
}

TEST_F(SimdMathTest, sqrt)
{
    // Just make sure sqrt(0)=0 works and isn't evaluated as 0*1/sqrt(0)=NaN
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(0, 2, 3), sqrt(setSimdRealFrom3R(0, 4, 9)));
}

/*! \brief Function wrapper to evaluate reference 1/x */
real refInv(real x)
{
    return 1.0/x;
}

TEST_F(SimdMathTest, inv)
{
    // test <0
    setRange(-1e10, -1e-10);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv);
    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv);
}

TEST_F(SimdMathTest, invMask)
{
    SimdReal x   = setSimdRealFrom3R(2.0, 0.0, 3.0);
    SimdBool m   = (SimdReal(0.0) < x);
    SimdReal ref = setSimdRealFrom3R(0.5, 0.0, 1.0/3.0);
    GMX_EXPECT_SIMD_REAL_NEAR(ref, invMask(x, m));
}

TEST_F(SimdMathTest, log)
{
    setRange(1e-30, 1e30);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::log, log);
}

TEST_F(SimdMathTest, exp2)
{
    setRange(-100, 100);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2);

    // We do not care about the SIMD implementation getting denormal values right,
    // but they must be clamped to zero rather than producing garbage.
    // Check by setting the absolute tolerance to machine precision.
    setAbsTol(GMX_REAL_EPS);

    // First two values will have denormal results in single, third value in double too.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::exp2(-150.0), std::exp2(-300.0), std::exp2(-1050.0)),
                              exp2(setSimdRealFrom3R(-150.0, -300.0, -1050.0)));

    // Reset absolute tolerance to enforce ULP checking
    setAbsTol(0.0);

    // Make sure that underflowing values are set to zero.
    // First two values underflow in single, third value in double too.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::exp2(-200.0), std::exp2(-600.0), std::exp2(-1500.0)),
                              exp2(setSimdRealFrom3R(-200.0, -600.0, -1500.0)));
}

TEST_F(SimdMathTest, exp)
{
    setRange(-75, 75);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, exp);

    // We do not care about the SIMD implementation getting denormal values right,
    // but they must be clamped to zero rather than producing garbage.
    // Check by setting the absolute tolerance to machine precision.
    setAbsTol(GMX_REAL_EPS);
    // First two values will have denormal results in single, third value in double too.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::exp(-90.0), std::exp(-100.0), std::exp(-725.0)),
                              exp(setSimdRealFrom3R(-90.0, -100.0, -725.0)));

    // Reset absolute tolerance to enforce ULP checking
    setAbsTol(0.0);

    // Make sure that underflowing values are set to zero.
    // First two values underflow in single, third value in double too.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::exp(-150.0), std::exp(-300.0), std::exp(-800.0)),
                              exp(setSimdRealFrom3R(-150.0, -300.0, -800.0)));
}

/*! \brief Function wrapper for erf(x), with argument/return in default Gromacs precision.
 *
 * \note Single-precision erf() in some libraries can be slightly lower precision
 * than the SIMD flavor, so we use a cast to force double precision for reference.
 */
real
refErf(real x)
{
    return std::erf(static_cast<double>(x));
}

TEST_F(SimdMathTest, erf)
{
    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    GMX_EXPECT_SIMD_FUNC_NEAR(refErf, erf);
}

/*! \brief Function wrapper for erfc(x), with argument/return in default Gromacs precision.
 *
 * \note Single-precision erfc() in some libraries can be slightly lower precision
 * than the SIMD flavor, so we use a cast to force double precision for reference.
 */
real
refErfc(real x)
{
    return std::erfc(static_cast<double>(x));
}

TEST_F(SimdMathTest, erfc)
{
    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    // Our erfc algorithm has 4 ulp accuracy, so relax defaultTol a bit
    setUlpTol(4*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(refErfc, erfc);
}

TEST_F(SimdMathTest, sin)
{
    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::sin, sin);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    setUlpTol(2*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::sin, sin);
}

TEST_F(SimdMathTest, cos)
{
    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cos, cos);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    setUlpTol(2*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cos, cos);
}

TEST_F(SimdMathTest, tan)
{
    // Tan(x) is a little sensitive due to the division in the algorithm.
    // Rather than using lots of extra FP operations, we accept the algorithm
    // presently only achieves a ~3 ulp error and use the medium tolerance.
    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::tan, tan);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    setUlpTol(2*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::tan, tan);
}

TEST_F(SimdMathTest, asin)
{
    // Our present asin(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::asin, asin);
}

TEST_F(SimdMathTest, acos)
{
    // Our present acos(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::acos, acos);
}

TEST_F(SimdMathTest, atan)
{
    // Our present atan(x) algorithm achieves 1 ulp accuracy
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::atan, atan);
}

TEST_F(SimdMathTest, atan2)
{
    // test each quadrant
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(std::atan2(1.0, 1.0)), atan2(rSimd_1_2_3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(std::atan2(-1.0, 1.0)), atan2(rSimd_m1_m2_m3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(std::atan2(-1.0, -1.0)), atan2(rSimd_m1_m2_m3, rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(std::atan2(1.0, -1.0)), atan2(rSimd_1_2_3, rSimd_m1_m2_m3));
    // cases important for calculating angles
    // values on coordinate axes
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(std::atan2(0.0, 1.0)), atan2(SimdReal(0.0), rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(std::atan2(1.0, 0.0)), atan2(rSimd_1_2_3, SimdReal(0.0)));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(std::atan2(0.0, -1.0)), atan2(SimdReal(0.0), rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(std::atan2(-1.0, 0.0)), atan2(rSimd_m1_m2_m3, SimdReal(0.0)));
    // degenerate value (origin) should return 0.0
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(std::atan2(0.0, 0.0)), atan2(setSimdRealFrom3R(0.0, 0.0, 0.0), SimdReal(0.0)));
}

/*! \brief Evaluate reference version of PME force correction. */
real
refPmeForceCorrection(real x)
{
    if (x != 0)
    {
        real y = std::sqrt(x);
        return 2*std::exp(-x)/(std::sqrt(M_PI)*x) - std::erf(static_cast<double>(y))/(x*y);
    }
    else
    {
        return -4/(3*std::sqrt(M_PI));
    }
}

// The PME corrections will be added to ~1/r2, so absolute tolerance of EPS is fine.
TEST_F(SimdMathTest, pmeForceCorrection)
{
    // Pme correction only needs to be ~1e-6 accuracy single, 1e-10 double
#ifdef GMX_DOUBLE
    setUlpTol(std::int64_t(5e-10/GMX_REAL_EPS));
#else
    setUlpTol(std::int64_t(5e-6/GMX_REAL_EPS));
#endif

    setRange(0.15, 4);
    setAbsTol(GMX_REAL_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmeForceCorrection, pmeForceCorrection);
}

/*! \brief Evaluate reference version of PME potential correction. */
real
refPmePotentialCorrection(real x)
{
    real y = std::sqrt(x);
    return std::erf(static_cast<double>(y))/y;
}

// The PME corrections will be added to ~1/r, so absolute tolerance of EPS is fine.
TEST_F(SimdMathTest, pmePotentialCorrection)
{
    // Pme correction only needs to be ~1e-6 accuracy single, 1e-10 double
#ifdef GMX_DOUBLE
    setUlpTol(std::int64_t(5e-10/GMX_REAL_EPS));
#else
    setUlpTol(std::int64_t(5e-6/GMX_REAL_EPS));
#endif
    setRange(0.15, 4);
    setAbsTol(GMX_REAL_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmePotentialCorrection, pmePotentialCorrection);
}



#if 0


/*
 * Functions that only target single accuracy, even for double SIMD data
 */

TEST_F(SimdMathTest, invsqrtSingleAccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, invsqrtSingleAccuracy);
}

/*! \brief Function wrapper to return first result when testing \ref simdInvsqrtPairSingleAccuracy */
SimdReal gmx_simdcall
tst_invsqrt_singleaccuracy_pair0(SimdReal x)
{
    SimdReal r0, r1;
    simdInvsqrtPairSingleAccuracy(x, x, &r0, &r1);
    return r0;
}

/*! \brief Function wrapper to return second result when testing \ref simdInvsqrtPairSingleAccuracy */
SimdReal gmx_simdcall
tst_invsqrt_singleaccuracy_pair1(SimdReal x)
{
    SimdReal r0, r1;
    simdInvsqrtPairSingleAccuracy(x, x, &r0, &r1);
    return r1;
}

TEST_F(SimdMathTest, gmxSimdInvsqrtPairSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tst_invsqrt_singleaccuracy_pair0);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tst_invsqrt_singleaccuracy_pair1);
}

TEST_F(SimdMathTest, gmxSimdSqrtSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Just make sure sqrt(0)=0 works and isn't evaluated as 0*1/sqrt(0)=NaN
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(0, 2, 3), simdSqrtSingleAccuracy(setSimdRealFrom3R(0, 4, 9)));
}

TEST_F(SimdMathTest, gmxSimdInvSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // test <0
    setRange(-1e10, -1e-10);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, simdInvSingleAccuracy);
    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, simdInvSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdLogSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(1e-30, 1e30);
    GMX_EXPECT_SIMD_FUNC_NEAR(refLog, simdLogSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdExp2SingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-100, 100);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, simdExp2SingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdExpSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-75, 75);
    GMX_EXPECT_SIMD_FUNC_NEAR(refExp, simdExpSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdErfSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    GMX_EXPECT_SIMD_FUNC_NEAR(refErf, simdErfSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdErfcSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    // Our erfc algorithm has 4 ulp accuracy, so relax tolerance a bit
    setUlpTol(4*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(refErfc, simdErfcSingleAccuracy);
}


TEST_F(SimdMathTest, gmxSimdSinSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(refSin, simdSinSingleAccuracy);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    setUlpTol(2*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(refSin, simdSinSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdCosSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(refCos, simdCosSingleAccuracy);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(refCos, simdCosSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdTanSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Tan(x) is a little sensitive due to the division in the algorithm.
    // Rather than using lots of extra FP operations, we accept the algorithm
    // presently only achieves a ~3 ulp error and use the medium tolerance.
    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(refTan, simdTanSingleAccuracy);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(refTan, simdTanSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdAsinSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Our present asin(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(refAsin, simdAsinSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdAcosSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Our present acos(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(refAcos, simdAcosSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdAtanSingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Our present atan(x) algorithm achieves 1 ulp accuracy
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(refAtan, simdAtanSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdAtan2SingleaccuracyR)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // test each quadrant
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(1.0, 1.0)), simdAtan2SingleAccuracy(rSimd_1_2_3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(-1.0, 1.0)), simdAtan2SingleAccuracy(rSimd_m1_m2_m3, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(-1.0, -1.0)), simdAtan2SingleAccuracy(rSimd_m1_m2_m3, rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(1.0, -1.0)), simdAtan2SingleAccuracy(rSimd_1_2_3, rSimd_m1_m2_m3));
    // cases important for calculating angles
    // values on coordinate axes
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(0.0, 1.0)), simdAtan2SingleAccuracy(SimdReal(0.0), rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(1.0, 0.0)), simdAtan2SingleAccuracy(rSimd_1_2_3, SimdReal(0.0)));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(0.0, -1.0)), simdAtan2SingleAccuracy(SimdReal(0.0), rSimd_m1_m2_m3));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(-1.0, 0.0)), simdAtan2SingleAccuracy(rSimd_m1_m2_m3, SimdReal(0.0)));
    // degenerate value (origin) should return 0.0
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(atan2(0.0, 0.0)), simdAtan2SingleAccuracy(setSimdRealFrom3R(0.0, 0.0, 0.0), SimdReal(0.0)));
}

TEST_F(SimdMathTest, gmxSimdPmecorrForceSingleaccuracyR)
{
    // The PME corrections will be added to ~1/r2, so absolute tolerance of EPS is fine.
    // Pme correction only needs to be ~1e-6 accuracy single.
    // Then increase the allowed error by the difference between the actual precision and single.
    setUlpTol( (std::int64_t(5e-6/GMX_FLOAT_EPS)) * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(0.15, 4);
    setAbsTol(GMX_FLOAT_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmeCorrForce, simdPmeCorrForceSingleAccuracy);
}

TEST_F(SimdMathTest, gmxSimdPmecorrPotentialSingleaccuracyR)
{
    // The PME corrections will be added to ~1/r, so absolute tolerance of EPS is fine.
    // Pme correction only needs to be ~1e-6 accuracy single.
    // Then increase the allowed error by the difference between the actual precision and single.
    setUlpTol( (std::int64_t(5e-6/GMX_FLOAT_EPS)) * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(0.15, 4);
    setAbsTol(GMX_FLOAT_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmeCorrPotential, simdPmeCorrPotentialSingleAccuracy);
}
#endif

}      // namespace

#endif // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace

#endif // GMX_SIMD
