/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \internal \file
 * \brief
 * Tests for arithmetic SIMD functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "general.h"
#include <algorithm>

namespace SIMDTests
{

//! Typedef for the test fixture
typedef SimdFunctionTest<real> SimdFunctionWithOneRealInput;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionWithOneRealInput::callFunction(SimdFunctionSet &simdFunctionSet,
                                           FunctionType     function,
                                           real            *_result)
{
    typename SimdFunctionSet::realType a, result;

    a      = simdFunctionSet.load_pr(inputs_[0]);
    result = function(a);
    simdFunctionSet.store_pr(_result, result);
}

/* SIMD reciprocal square roots are computed by an approximate SIMD
 * operation, and then some Newton-Raphson iterations depending on the
 * accuracy of the underlying hardware operation and the prevailing
 * GROMACS single/double precision.
 *
 * 1) The SIMD rsqrt is never used without subsequent Newton-Raphson
 * iteration, and there is no reference version of approximate inverse
 * square root, so it only makes sense to test the composite
 * operation (i.e. not gmx_invsqrt_pr).
 *
 * 2) GROMACS does not do enough Newton-Raphson iterations for full
 * double precision (because we don't need it), so a wider range for
 * the equality test is required in that case. */
TEST_F(SimdFunctionWithOneRealInput, gmx_rsqrt_pr_Works)
{
#ifdef GMX_DOUBLE
    maxUlps = 80.0; // empirically determined to be enough on x86
#endif
    /* The generation of uniform reals is over [0,1), but we would
       like to test over (0,5]. Ideally, one could subtract 1 and then
       multiply by -5, but that does not work in practice with
       floating-point arithmetic. */
    prepare(1, 1, -1.00001, -5);
    Tester(ReferenceFunctions::gmx_simd_ref_rsqrt_pr,
           TestFunctions::gmx_invsqrt_pr, real(0));
}

/* The SIMD reciprocal is never used without subsequent Newton-Raphson
 * iteration, and there is no reference version of approximate
 * reciprocal, so it only makes sense to test the composite
 * operation. Similarly with rsqrt, the accuracy requirement is too
 * strict in this case.
 *
 * The machinery that generates the random reals over which we test
 * can't easily generate a range with both positive and negative and
 * never zero, so we test positive and negative
 * separately. */
TEST_F(SimdFunctionWithOneRealInput, gmx_inv_pr_WorksOnNegativeValues)
{
#ifdef GMX_DOUBLE
    maxUlps = 40.0; // empirically determined to be enough on x86
#else
    maxUlps = 10.0; // empirically determined to be enough on x86
#endif
    /* Subtracting a simple -1 does not work in practice, as above. */
    prepare(1, 1, -1.00001, -5);
    Tester(ReferenceFunctions::gmx_simd_ref_rcp_pr,
           TestFunctions::gmx_inv_pr, real(0));
}
TEST_F(SimdFunctionWithOneRealInput, gmx_inv_pr_WorksOnPositiveValues)
{
#ifdef GMX_DOUBLE
    maxUlps = 40.0; // empirically determined to be enough on x86
#else
    maxUlps = 10.0; // empirically determined to be enough on x86
#endif
    /* Subtracting a simple -1 does not work in practice, as above. */
    prepare(1, 1, -1.00001, -5);
    Tester(ReferenceFunctions::gmx_simd_ref_rcp_pr,
           TestFunctions::gmx_inv_pr, real(0));
}

TEST_F(SimdFunctionWithOneRealInput, gmx_sqrt_pr_Works)
{
    prepare(1, 1, 0, 5);
    Tester(ReferenceFunctions::gmx_simd_ref_sqrt_pr,
           TestFunctions::gmx_sqrt_pr, real(0));
}

TEST_F(SimdFunctionWithOneRealInput, gmx_exp_pr_Works)
{
    /* gmx_exp_pr uses a rational function approximation, which is not
       accurate enough at double precision. */
#ifdef GMX_DOUBLE
    maxUlps = 8.0; // empirically determined to be enough on x86
#endif
    prepare(1, 1, 0, 5);
    Tester(ReferenceFunctions::gmx_simd_ref_exp_pr,
           TestFunctions::gmx_exp_pr, real(0));
}

TEST_F(SimdFunctionWithOneRealInput, gmx_acos_pr_Works)
{
    /* gmx_acos_pr uses a rational function approximation, which is
       not accurate enough at double precision. */
#ifdef GMX_DOUBLE
    maxUlps = 80.0; // empirically determined to be enough on x86
#endif
    prepare(1, 1, -0.5, 2);
    Tester(ReferenceFunctions::gmx_simd_ref_acos_pr,
           TestFunctions::gmx_acos_pr, real(0));
}

TEST_F(SimdFunctionWithOneRealInput, gmx_round_pr_Works)
{
    /* This function is used for rounding of reals before periodicity
     * checks, so we need to check a realistic range of distances. */
    prepare(1, 1, -0.5, 100);
    Tester(ReferenceFunctions::gmx_simd_ref_round_pr,
           TestFunctions::gmx_round_pr, real(0));
}

/* The set of functions that we'd like to test vary with the kind of
   acceleration being used. */
#ifdef GMX_SIMD_HAVE_FLOOR

TEST_F(SimdFunctionWithOneRealInput, gmx_floor_pr_Works)
{
    /* This function is used for truncation of reals before table
     * lookups, so moderate positive values are those most worth
     * testing. */
    prepare(1, 1, 0, 1024);
    Tester(ReferenceFunctions::gmx_simd_ref_floor_pr,
           TestFunctions::gmx_floor_pr, real(0));
}

#else

namespace ReferenceFunctions
{

/*! Wrapper functions thats test the way gmx_cvtepi32_pr and
   gmx_cvttpr_epi32 are used in practice. */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_truncate_real_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        b.r[i] = (real)(int)(a.r[i]);
    }

    return b;
}

}

namespace TestFunctions
{

/*! Wrapper functions thats test the way gmx_cvtepi32_pr and
   gmx_cvttpr_epi32 are used in practice. */
static gmx_inline gmx_mm_pr gmx_truncate_real_pr(gmx_mm_pr a)
{
    return gmx_cvtepi32_pr(gmx_cvttpr_epi32(a));
}

}

/* This tests both gmx_cvttpr_epi32 and gmx_cvtepi32_pr. These two are
 * inconvenient to test separately. gmx_cvtepi32_pr is only ever used
 * after gmx_cvttpr_epi32, and only when gmx_floor_pr is
 * unavailable. gmx_cvttpr_epi32 is also used before table lookups,
 * which are tested elsewhere. */

TEST_F(SimdFunctionWithOneRealInput, FunctionPairDoingTruncationOfReals_Works)
{
    /* See comment to gmx_round_pr_Works */
    prepare(1, 1, -0.25, 1024);
    Tester(ReferenceFunctions::gmx_simd_ref_truncate_real_pr,
           TestFunctions::gmx_truncate_real_pr, real(0));
}
#endif

} // namespace
