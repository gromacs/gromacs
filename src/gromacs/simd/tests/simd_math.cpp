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

#include "gromacs/simd/simd_math.h"

#include "config.h"

#include <cmath>
#include <cstdint>

#include <vector>

#include "gromacs/math/utilities.h"
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
                            compareSimdMathFunction(const char *  refFuncExpr,
                                                    const char *  simdFuncExpr,
                                                    const char *  denormalsToZeroExpr,
                                                    real          refFunc(real x),
                                                    SimdReal      gmx_simdcall simdFunc(SimdReal x),
                                                    bool          denormalsToZero);
};

/*! \brief Test approximate equality of SIMD vs reference version of a function.
 *
 * This macro takes vanilla C and SIMD flavors of a function and tests it with
 * the number of points, range, and tolerances specified by the test fixture class.
 */
#define GMX_EXPECT_SIMD_FUNC_NEAR(refFunc, tstFunc) \
    EXPECT_PRED_FORMAT3(compareSimdMathFunction, refFunc, tstFunc, false)

/*! \brief Test approximate equality of SIMD vs reference function, denormals can be zero.
 *
 * This macro takes vanilla C and SIMD flavors of a function and tests it with
 * the number of points, range, and tolerances specified by the test fixture class.
 *
 * This version of the function will also return success if the test function
 * returns zero where the reference function returns a denormal value.
 */
#define GMX_EXPECT_SIMD_FUNC_NEAR_DTZ(refFunc, tstFunc) \
    EXPECT_PRED_FORMAT3(compareSimdMathFunction, refFunc, tstFunc, true)

/*! \brief Implementation routine to compare SIMD vs reference functions.
 *
 * \param refFuncExpr         Description of reference function expression
 * \param simdFuncExpr        Description of SIMD function expression
 * \param denormalsToZeroExpr Description of denormal-to-zero setting
 * \param refFunc             Reference math function pointer
 * \param simdFunc            SIMD math function pointer
 * \param denormalsToZero     If true, the function will consider denormal
 *                            values equivalent to 0.0.
 *
 * The function will be tested with the range and tolerances specified in
 * the SimdBaseTest class. You should not never call this function directly,
 * but use the macro GMX_EXPECT_SIMD_FUNC_NEAR(refFunc,tstFunc) instead.
 */
::testing::AssertionResult
SimdMathTest::compareSimdMathFunction(const char              * refFuncExpr,
                                      const char              * simdFuncExpr,
                                      const char              * denormalsToZeroExpr,
                                      real                      refFunc(real x),
                                      SimdReal     gmx_simdcall simdFunc(SimdReal x),
                                      bool                      denormalsToZero)
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
#    if GMX_DOUBLE
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
            if (denormalsToZero)
            {
                // Clamp denormal values to zero if requested
                if (std::abs(vref[i]) <= GMX_REAL_MIN)
                {
                    vref[i] = 0.0;
                }
                if (std::abs(vtst[i]) <= GMX_REAL_MIN)
                {
                    vtst[i] = 0.0;
                }
            }

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
                   << "Denormals can be 0: " << denormalsToZeroExpr << std::endl
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
               << "Denormals can be 0: " << denormalsToZeroExpr << std::endl
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
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-c0,  c1, c2), copysign(setSimdRealFrom3R( c0,  c1,  c2), setSimdRealFrom3R(-c3,  c4, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-c0,  c1, c2), copysign(setSimdRealFrom3R(-c0, -c1, -c2), setSimdRealFrom3R(-c3,  c4, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R( c0, -c1, c2), copysign(setSimdRealFrom3R( c0,  c1,  c2), setSimdRealFrom3R( c3, -c4, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R( c0, -c1, c2), copysign(setSimdRealFrom3R(-c0, -c1, -c2), setSimdRealFrom3R( c3, -c4, 0)));
}

/*! \brief Function wrapper to evaluate reference 1/sqrt(x) */
static real
refInvsqrt(real x)
{
    return 1.0/std::sqrt(x);
}

TEST_F(SimdMathTest, invsqrt)
{
    setRange(1.01*GMX_FLOAT_MIN, GMX_FLOAT_MAX);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, invsqrt);
}

TEST_F(SimdMathTest, maskzInvsqrt)
{
    SimdReal x   = setSimdRealFrom3R(c1, 0.0, c2);
    SimdBool m   = (setZero() < x);
    SimdReal ref = setSimdRealFrom3R(1.0/std::sqrt(c1), 0.0, 1.0/std::sqrt(c2));
    GMX_EXPECT_SIMD_REAL_NEAR(ref, maskzInvsqrt(x, m));
}

/*! \brief Function wrapper to return first result when testing \ref invsqrtPair */
SimdReal gmx_simdcall
tstInvsqrtPair0(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPair(x, x, &r0, &r1);
    return r0;
}

/*! \brief Function wrapper to return second result when testing \ref invsqrtPair */
SimdReal gmx_simdcall
tstInvsqrtPair1(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPair(x, x, &r0, &r1);
    return r1;
}

TEST_F(SimdMathTest, invsqrtPair)
{
    setRange(1.01*GMX_FLOAT_MIN, GMX_FLOAT_MAX);
    // The accuracy conversions lose a bit of extra accuracy compared to
    // doing the iterations in all-double.
    setUlpTol(4*ulpTol_);

    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tstInvsqrtPair0);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tstInvsqrtPair1);
}

/*! \brief Function wrapper to evaluate reference sqrt(x) */
static real
refSqrt(real x)
{
    return std::sqrt(x);
}

/*! \brief Dummy function returning 0.0 to test function ranges that should be zero */
gmx_unused static real
refZero(real gmx_unused x)
{
    return 0.0;
}


TEST_F(SimdMathTest, sqrt)
{
    // The accuracy conversions lose a bit of extra accuracy compared to
    // doing the iterations in all-double.
    setUlpTol(4*ulpTol_);

    // First test that 0.0 and a few other values works
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(0, std::sqrt(c1), std::sqrt(c2)), sqrt(setSimdRealFrom3R(0, c1, c2)));

    // Values smaller-than-or-equal to GMX_FLOAT_MIN will be clamped to 0.0,
    // so only test larger values
    setRange(1.01*GMX_FLOAT_MIN, GMX_FLOAT_MAX);
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrt);

#if GMX_DOUBLE
    // Make sure that values smaller than GMX_FLOAT_MIN lead to result 0.0
    setRange(0.0, 0.99*GMX_FLOAT_MIN);
    GMX_EXPECT_SIMD_FUNC_NEAR(refZero, sqrt);
#endif
}

TEST_F(SimdMathTest, sqrtUnsafe)
{
    // The accuracy conversions lose a bit of extra accuracy compared to
    // doing the iterations in all-double.
    setUlpTol(4*ulpTol_);

    setRange(GMX_FLOAT_MIN, GMX_FLOAT_MAX);
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrt<MathOptimization::Unsafe>);
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

TEST_F(SimdMathTest, maskzInv)
{
    SimdReal x   = setSimdRealFrom3R(c1, 0.0, c2);
    SimdBool m   = (setZero() < x);
    SimdReal ref = setSimdRealFrom3R(1.0/c1, 0.0, 1.0/c2);
    GMX_EXPECT_SIMD_REAL_NEAR(ref, maskzInv(x, m));
}

TEST_F(SimdMathTest, log)
{
    setRange(1e-30, 1e30);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::log, log);
}

TEST_F(SimdMathTest, exp2)
{
    // Test normal/denormal/zero range separately to make errors clearer

    // First test the range where we get normalized (non-denormal) results,
    // since we don't require denormal results to be reproduced correctly.
#if GMX_DOUBLE
    setRange(-1022, 1023);
#else
    setRange(-126, 127);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2);

    // Some implementations might have denormal support, in which case they
    // support an extended range, adding roughly the number of bits in the
    // mantissa to the smallest allowed arg (1023+52 in double, 127+23 single).
    // In this range we allow the value to be either correct (denormal) or 0.0
#if GMX_DOUBLE
    setRange(-1075, -1022);
#else
    setRange(-150, -126);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR_DTZ(std::exp2, exp2);

    // For arguments smaller than the subnormal the result should be zero
    // both in the reference and our implementations.
#if GMX_DOUBLE
    setRange(-1000000.0, -1075.0);
#else
    setRange(-100000.0, -150.0);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2);

    // Test a few very negative values, including values so small that they
    // will start to cause inf values in the polynomial interpolations
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::exp2(-GMX_FLOAT_MAX), std::exp2(-0.1*GMX_REAL_MAX), std::exp2(-GMX_REAL_MAX)),
                              exp2(setSimdRealFrom3R(-GMX_FLOAT_MAX, -0.1*GMX_REAL_MAX, -GMX_REAL_MAX)));
}

TEST_F(SimdMathTest, exp2Unsafe)
{
    // The unsafe version is only defined in this range
#if GMX_DOUBLE
    setRange(-1022, 1023);
#else
    setRange(-126, 127);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2<MathOptimization::Unsafe>);
}

TEST_F(SimdMathTest, exp)
{
    // Test normal/denormal/zero range separately to make errors clearer

    // First test the range where we get normalized (non-denormal) results,
    // since we don't require denormal results to be reproduced correctly.
#if GMX_DOUBLE
    setRange(-708.3, 709.1);
#else
    setRange(-87.3, 88.0);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, exp);

    // Some implementations might have denormal support, in which case they
    // support an extended range, adding roughly the number of bits in the
    // mantissa to the smallest allowed arg (1023+52 in double, 127+23 single).
    // Then multiply with ln(2) to get our limit for exp().
    // In this range we allow the value to be either correct (denormal) or 0.0
#if GMX_DOUBLE
    setRange(-746.0, -708.3);
#else
    setRange(-104.0, -87.3);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR_DTZ(std::exp, exp);

    // For arguments smaller than the subnormal the result should be zero
    // both in the reference and our implementations.
#if GMX_DOUBLE
    setRange(-1000000.0, -746.0);
#else
    setRange(-100000.0, -104.0);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, exp);

    // Test a few very negative values, including values so small that they
    // will start to cause inf values in the polynomial interpolations
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::exp(-GMX_FLOAT_MAX), std::exp(-0.1*GMX_REAL_MAX), std::exp(-GMX_REAL_MAX)),
                              exp(setSimdRealFrom3R(-GMX_FLOAT_MAX, -0.1*GMX_REAL_MAX, -GMX_REAL_MAX)));
}

TEST_F(SimdMathTest, expUnsafe)
{
#if GMX_DOUBLE
    setRange(-708.3, 709.1);
#else
    setRange(-87.3, 88.0);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, exp<MathOptimization::Unsafe>);
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
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(c0, c3), std::atan2(c1, c4), std::atan2(c2, c5)),
                              atan2(rSimd_c0c1c2, rSimd_c3c4c5));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(-c0, c3), std::atan2(-c1, c4), std::atan2(-c2, c5)),
                              atan2(rSimd_m0m1m2, rSimd_c3c4c5));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(-c0, -c3), std::atan2(-c1, -c0), std::atan2(-c2, -c4)),
                              atan2(rSimd_m0m1m2, rSimd_m3m0m4));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(c0, -c3), std::atan2(c1, -c0), std::atan2(c2, -c4)),
                              atan2(rSimd_c0c1c2, rSimd_m3m0m4));

    // cases important for calculating angles
    // values on coordinate axes
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(0, c0), std::atan2(0, c1), std::atan2(0, c2)),
                              atan2(setZero(), rSimd_c0c1c2));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(c0, 0), std::atan2(c1, 0), std::atan2(c2, 0)),
                              atan2(rSimd_c0c1c2, setZero()));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(0, -c0), std::atan2(0, -c1), std::atan2(0, -c2)),
                              atan2(setZero(), rSimd_m0m1m2));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(-c0, 0), std::atan2(-c1, 0), std::atan2(-c2, 0)),
                              atan2(rSimd_m0m1m2, setZero()));
    // degenerate value (origin) should return 0.0. At least IBM xlc 13.1.5 gets the reference
    // value wrong (-nan) at -O3 optimization, so we compare to the correct value (0.0) instead.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(0.0), atan2(setSimdRealFrom3R(0.0, 0.0, 0.0), setZero()));
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
#if GMX_DOUBLE
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
#if GMX_DOUBLE
    setUlpTol(std::int64_t(5e-10/GMX_REAL_EPS));
#else
    setUlpTol(std::int64_t(5e-6/GMX_REAL_EPS));
#endif
    setRange(0.15, 4);
    setAbsTol(GMX_REAL_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmePotentialCorrection, pmePotentialCorrection);
}

// Functions that only target single accuracy, even for double SIMD data

TEST_F(SimdMathTest, invsqrtSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(1.01*GMX_FLOAT_MIN, GMX_FLOAT_MAX);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, invsqrtSingleAccuracy);
}

/*! \brief Function wrapper to return first result when testing \ref invsqrtPairSingleAccuracy */
SimdReal gmx_simdcall
tst_invsqrt_SingleAccuracy_pair0(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPairSingleAccuracy(x, x, &r0, &r1);
    return r0;
}

/*! \brief Function wrapper to return second result when testing \ref invsqrtPairSingleAccuracy */
SimdReal gmx_simdcall
tst_invsqrt_SingleAccuracy_pair1(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPairSingleAccuracy(x, x, &r0, &r1);
    return r1;
}

TEST_F(SimdMathTest, invsqrtPairSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(1.01*GMX_FLOAT_MIN, GMX_FLOAT_MAX);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tst_invsqrt_SingleAccuracy_pair0);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tst_invsqrt_SingleAccuracy_pair1);
}

TEST_F(SimdMathTest, sqrtSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // First test that 0.0 and a few other values works
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(0, std::sqrt(c0), std::sqrt(c1)), sqrtSingleAccuracy(setSimdRealFrom3R(0, c0, c1)));

    // Values smaller-than-or-equal to GMX_FLOAT_MIN will be clamped to 0.0,
    // so only test larger values
    setRange(1.01*GMX_FLOAT_MIN, GMX_FLOAT_MAX);
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrtSingleAccuracy);

#if GMX_DOUBLE
    // Make sure that values smaller than GMX_FLOAT_MIN lead to result 0.0
    setRange(0.0, 0.99*GMX_FLOAT_MIN);
    GMX_EXPECT_SIMD_FUNC_NEAR(refZero, sqrtSingleAccuracy);
#endif
}

TEST_F(SimdMathTest, sqrtSingleAccuracyUnsafe)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Test the full range
    setRange(GMX_FLOAT_MIN, GMX_FLOAT_MAX);
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrtSingleAccuracy<MathOptimization::Unsafe>);
}

TEST_F(SimdMathTest, invSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // test <0
    setRange(-1e10, -1e-10);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, invSingleAccuracy);
    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, invSingleAccuracy);
}

TEST_F(SimdMathTest, logSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(1e-30, 1e30);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::log, logSingleAccuracy);
}

TEST_F(SimdMathTest, exp2SingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

#if GMX_DOUBLE
    setRange(-1022.49, 1023.49);
#else
    setRange(-126, 127.49);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2SingleAccuracy);

    // Test a range that should be zero both for reference and simd versions.
    // Some implementations might have subnormal support, in which case they
    // support an extended range, adding roughly the number of bits in the
    // mantissa to the smallest allowed arg (1023+52 in double, 127+23 single).
#if GMX_DOUBLE
    setRange(-1000000.0, -1075.0);
#else
    setRange(-100000.0, -150.0);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR_DTZ(std::exp2, exp2SingleAccuracy);

    // Test a few very negative values
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::exp2(-GMX_FLOAT_MAX), std::exp2(-0.1*GMX_REAL_MAX), std::exp2(-GMX_REAL_MAX)),
                              exp2SingleAccuracy(setSimdRealFrom3R(-GMX_FLOAT_MAX, -0.1*GMX_REAL_MAX, -GMX_REAL_MAX)));
}

TEST_F(SimdMathTest, exp2SingleAccuracyUnsafe)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

#if GMX_DOUBLE
    setRange(-1022.49, 1023.49);
#else
    setRange(-126, 127.49);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2SingleAccuracy<MathOptimization::Unsafe>);
}

TEST_F(SimdMathTest, expSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

#if GMX_DOUBLE
    setRange(-708.7, 709.4);
#else
    setRange(-87.3, 88.3);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, expSingleAccuracy);

    // Test a range that should be zero both for reference and simd versions.
    // Some implementations might have subnormal support, in which case they
    // support an extended range, adding roughly the number of bits in the
    // mantissa to the smallest result exponent (1023+52 in double, 127+23 single).
    // Then multiply with ln(2) to get our limit for exp().
#if GMX_DOUBLE
    setRange(-1000000.0, -746.0);
#else
    setRange(-100000.0, -104.0);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR_DTZ(std::exp, expSingleAccuracy);

    // Test a few very negative values, including values so small that they
    // will start to cause inf values in the polynomial interpolations
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::exp(-GMX_FLOAT_MAX), std::exp(-0.1*GMX_REAL_MAX), std::exp(-GMX_REAL_MAX)),
                              expSingleAccuracy(setSimdRealFrom3R(-GMX_FLOAT_MAX, -0.1*GMX_REAL_MAX, -GMX_REAL_MAX)));
}

TEST_F(SimdMathTest, expSingleAccuracyUnsafe)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

#if GMX_DOUBLE
    setRange(-708.7, 709.4);
#else
    setRange(-87.3, 88.3);
#endif
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, expSingleAccuracy<MathOptimization::Unsafe>);
}

TEST_F(SimdMathTest, erfSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    GMX_EXPECT_SIMD_FUNC_NEAR(refErf, erfSingleAccuracy);
}

TEST_F(SimdMathTest, erfcSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-9, 9);
    setAbsTol(GMX_REAL_MIN);
    // Our erfc algorithm has 4 ulp accuracy, so relax tolerance a bit
    setUlpTol(4*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(refErfc, erfcSingleAccuracy);
}


TEST_F(SimdMathTest, sinSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::sin, sinSingleAccuracy);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    setUlpTol(2*ulpTol_);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::sin, sinSingleAccuracy);
}

TEST_F(SimdMathTest, cosSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cos, cosSingleAccuracy);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cos, cosSingleAccuracy);
}

TEST_F(SimdMathTest, tanSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Tan(x) is a little sensitive due to the division in the algorithm.
    // Rather than using lots of extra FP operations, we accept the algorithm
    // presently only achieves a ~3 ulp error and use the medium tolerance.
    setRange(-8*M_PI, 8*M_PI);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::tan, tanSingleAccuracy);
    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::tan, tanSingleAccuracy);
}

TEST_F(SimdMathTest, asinSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Our present asin(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::asin, asinSingleAccuracy);
}

TEST_F(SimdMathTest, acosSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Our present acos(x) algorithm achieves 2-3 ulp accuracy
    setRange(-1, 1);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::acos, acosSingleAccuracy);
}

TEST_F(SimdMathTest, atanSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // Our present atan(x) algorithm achieves 1 ulp accuracy
    setRange(-10000, 10000);
    GMX_EXPECT_SIMD_FUNC_NEAR(std::atan, atanSingleAccuracy);
}

TEST_F(SimdMathTest, atan2SingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    // test each quadrant
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(c0, c3), std::atan2(c1, c4), std::atan2(c2, c5)),
                              atan2SingleAccuracy(rSimd_c0c1c2, rSimd_c3c4c5));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(-c0, c3), std::atan2(-c1, c4), std::atan2(-c2, c5)),
                              atan2SingleAccuracy(rSimd_m0m1m2, rSimd_c3c4c5));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(-c0, -c3), std::atan2(-c1, -c0), std::atan2(-c2, -c4)),
                              atan2SingleAccuracy(rSimd_m0m1m2, rSimd_m3m0m4));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(c0, -c3), std::atan2(c1, -c0), std::atan2(c2, -c4)),
                              atan2SingleAccuracy(rSimd_c0c1c2, rSimd_m3m0m4));
    // cases important for calculating angles
    // values on coordinate axes
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(0, c0), std::atan2(0, c1), std::atan2(0, c2)),
                              atan2SingleAccuracy(setZero(), rSimd_c0c1c2));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(c0, 0), std::atan2(c1, 0), std::atan2(c2, 0)),
                              atan2SingleAccuracy(rSimd_c0c1c2, setZero()));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(0, -c0), std::atan2(0, -c1), std::atan2(0, -c2)),
                              atan2SingleAccuracy(setZero(), rSimd_m0m1m2));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(-c0, 0), std::atan2(-c1, 0), std::atan2(-c2, 0)),
                              atan2SingleAccuracy(rSimd_m0m1m2, setZero()));

    // degenerate value (origin) should return 0.0. At least IBM xlc 13.1.5 gets the reference
    // value wrong (-nan) at -O3 optimization, so we compare to the correct value (0.0) instead.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(0.0), atan2SingleAccuracy(setSimdRealFrom3R(0.0, 0.0, 0.0), setZero()));
}

TEST_F(SimdMathTest, pmeForceCorrectionSingleAccuracy)
{
    // The PME corrections will be added to ~1/r2, so absolute tolerance of EPS is fine.
    // Pme correction only needs to be ~1e-6 accuracy single.
    // Then increase the allowed error by the difference between the actual precision and single.
    setUlpTol( (std::int64_t(5e-6/GMX_FLOAT_EPS)) * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(0.15, 4);
    setAbsTol(GMX_FLOAT_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmeForceCorrection, pmeForceCorrectionSingleAccuracy);
}

TEST_F(SimdMathTest, pmePotentialCorrectionSingleAccuracy)
{
    // The PME corrections will be added to ~1/r, so absolute tolerance of EPS is fine.
    // Pme correction only needs to be ~1e-6 accuracy single.
    // Then increase the allowed error by the difference between the actual precision and single.
    setUlpTol( (std::int64_t(5e-6/GMX_FLOAT_EPS)) * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));

    setRange(0.15, 4);
    setAbsTol(GMX_FLOAT_EPS);
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmePotentialCorrection, pmePotentialCorrectionSingleAccuracy);
}

}      // namespace

#endif // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace

#endif // GMX_SIMD
