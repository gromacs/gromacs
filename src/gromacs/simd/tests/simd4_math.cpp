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

#include "config.h"

#include <vector>

#include "gromacs/math/utilities.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"

#include "simd4.h"

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#ifdef GMX_SIMD4_HAVE_REAL

class Simd4MathTest : public Simd4Test
{
    public:
        ::testing::AssertionResult
                             compareSimd4MathFunction(const char * refFuncExpr, const char *simd4FuncExpr,
                                                      real refFunc(real x),     gmx_simd4_real_t gmx_simdcall simd4Func(gmx_simd4_real_t x));
};

/*! \brief Test approximate equality of SIMD4 vs reference version of a function.
 *
 * This macro takes vanilla C and SIMD flavors of a function and tests it with
 * the number of points, range, and tolerances specified by the test fixture class.
 */
#define GMX_EXPECT_SIMD4_FUNC_NEAR(refFunc, tstFunc) \
    EXPECT_PRED_FORMAT2(compareSimd4MathFunction, refFunc, tstFunc)


/*! \brief Implementation routine to compare SIMD4 vs reference functions.
 *
 * \param refFuncExpr    Description of reference function expression
 * \param simd4FuncExpr  Description of SIMD function expression
 * \param refFunc        Reference math function pointer
 * \param simd4Func      SIMD math function pointer
 *
 * The function will be tested with the range and tolerances specified in
 * the SimdBaseTest class. You should not never call this function directly,
 * but use the macro GMX_EXPECT_SIMD4_FUNC_NEAR(refFunc,tstFunc) instead.
 */
::testing::AssertionResult
Simd4MathTest::compareSimd4MathFunction(const char * refFuncExpr, const char *simd4FuncExpr,
                                        real refFunc(real x),     gmx_simd4_real_t gmx_simdcall simd4Func(gmx_simd4_real_t x))
{
    std::vector<real>            vx(GMX_SIMD4_WIDTH);
    std::vector<real>            vref(GMX_SIMD4_WIDTH);
    std::vector<real>            vtst(GMX_SIMD4_WIDTH);
    real                         dx;
    gmx_int64_t                  ulpDiff, maxUlpDiff;
    real                         maxUlpDiffPos;
    real                         refValMaxUlpDiff, simdValMaxUlpDiff;
    bool                         eq, signOk;
    int                          i, iter;
    int                          niter   = s_nPoints/GMX_SIMD4_WIDTH;
    int                          npoints = niter*GMX_SIMD4_WIDTH;
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
        for (i = 0; i < GMX_SIMD4_WIDTH; i++)
        {
            vx[i]   = range_.first+dx*(iter*GMX_SIMD4_WIDTH+i);
            vref[i] = refFunc(vx[i]);
        }
        vtst  = simd4Real2Vector(simd4Func(vector2Simd4Real(vx)));

        for (i = 0, eq = true, signOk = true; i < GMX_SIMD4_WIDTH && eq == true; i++)
        {
            eq     = eq && ( fabs(vref[i]-vtst[i]) < absTol_ );
            signOk = signOk && ( vref[i]*vtst[i] >= 0 );
        }
        if (eq == true)
        {
            // Go to next point if everything within absolute tolerance
            continue;
        }
        else if (signOk == false)
        {
            return ::testing::AssertionFailure()
                   << "Failing SIMD4 math function comparison due to sign differences." << std::endl
                   << "Reference function: " << refFuncExpr << std::endl
                   << "Simd function:      " << simd4FuncExpr << std::endl
                   << "Test range is ( " << range_.first << " , " << range_.second << " ) " << std::endl
                   << "First sign difference around x=" << std::setprecision(20) << ::testing::PrintToString(vx) << std::endl
                   << "Ref values:   " << std::setprecision(20) << ::testing::PrintToString(vref) << std::endl
                   << "SIMD4 values: " << std::setprecision(20) << ::testing::PrintToString(vtst) << std::endl;
        }
        /* We replicate the trivial ulp differences comparison here rather than
         * calling the lower-level routine for comparing them, since this enables
         * us to run through the entire test range and report the largest deviation
         * without lots of extra glue routines.
         */
        for (i = 0; i < GMX_SIMD4_WIDTH; i++)
        {
            conv0.r = vref[i];
            conv1.r = vtst[i];
            ulpDiff = llabs(conv0.i-conv1.i);
            if (ulpDiff > maxUlpDiff)
            {
                maxUlpDiff        = ulpDiff;
                maxUlpDiffPos     = vx[i];
                refValMaxUlpDiff  = vref[i];
                simdValMaxUlpDiff = vtst[i];
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
               << "Failing SIMD4 math function ulp comparison between " << refFuncExpr << " and " << simd4FuncExpr << std::endl
               << "Requested ulp tolerance: " << ulpTol_ << std::endl
               << "Requested abs tolerance: " << absTol_ << std::endl
               << "Largest Ulp difference occurs for x=" << std::setprecision(20) << maxUlpDiffPos << std::endl
               << "Ref  values:  " << std::setprecision(20) << refValMaxUlpDiff << std::endl
               << "SIMD4 values: " << std::setprecision(20) << simdValMaxUlpDiff << std::endl
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

/*! \brief Function wrapper to evaluate reference 1/sqrt(x) */
static real
ref_invsqrt(real x)
{
    return 1.0/sqrt(x);
}

TEST_F(Simd4MathTest, gmxSimd4InvsqrtR)
{
    setRange(1e-10, 1e10);
    GMX_EXPECT_SIMD4_FUNC_NEAR(ref_invsqrt, gmx_simd4_invsqrt_r);
}

TEST_F(Simd4MathTest, gmxSimd4InvsqrtSingleaccuracyR)
{
    setRange(1e-10, 1e10);
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTol(ulpTol_ * (1LL << (std::numeric_limits<real>::digits-std::numeric_limits<float>::digits)));
    GMX_EXPECT_SIMD4_FUNC_NEAR(ref_invsqrt, gmx_simd4_invsqrt_singleaccuracy_r);
}

}      // namespace

#endif // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
