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

typedef SimdFunctionTest<1, 1, real> SimdFunctionWithOneRealInput;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionWithOneRealInput::callFunction(SimdFunctionSet &simdFunctionSet,
                                           FunctionType     function,
                                           real            *_result)
{
    typename SimdFunctionSet::realType a, result;

    a      = simdFunctionSet.load_pr(inputs[0]);
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
    /* The generation of uniform reals is over [0,1), but we would
       like to test over (0,5]. Ideally, one could subtract 1 and then
       multiply by -5, but that does not work in practice with
       floating-point arithmetic. */
    std::fill(shift.begin(), shift.end(), -1.00001);
    std::fill(scale.begin(), scale.end(), -5);
#ifdef GMX_DOUBLE
    maxUlps = 80.0; // empirically determined to be enough on x86
#endif
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
    /* Subtracting a simple -1 does not work in practice, as above. */
    std::fill(shift.begin(), shift.end(), -1.00001);
    std::fill(scale.begin(), scale.end(), 5);
#ifdef GMX_DOUBLE
    maxUlps = 40.0; // empirically determined to be enough on x86
#else
    maxUlps = 10.0; // empirically determined to be enough on x86
#endif
    Tester(ReferenceFunctions::gmx_simd_ref_rcp_pr,
           TestFunctions::gmx_inv_pr, real(0));
}
TEST_F(SimdFunctionWithOneRealInput, gmx_inv_pr_WorksOnPositiveValues)
{
    /* Subtracting a simple -1 does not work in practice, as above. */
    std::fill(shift.begin(), shift.end(), -1.00001);
    std::fill(scale.begin(), scale.end(), -5);
#ifdef GMX_DOUBLE
    maxUlps = 40.0; // empirically determined to be enough on x86
#else
    maxUlps = 10.0; // empirically determined to be enough on x86
#endif
    Tester(ReferenceFunctions::gmx_simd_ref_rcp_pr,
           TestFunctions::gmx_inv_pr, real(0));
}

TEST_F(SimdFunctionWithOneRealInput, gmx_sqrt_pr_Works)
{
    std::fill(shift.begin(), shift.end(), 0);
    std::fill(scale.begin(), scale.end(), 5);
    Tester(ReferenceFunctions::gmx_simd_ref_sqrt_pr,
           TestFunctions::gmx_sqrt_pr, real(0));
}

TEST_F(SimdFunctionWithOneRealInput, gmx_exp_pr_Works)
{
    std::fill(shift.begin(), shift.end(), 0);
    std::fill(scale.begin(), scale.end(), 5);
    /* gmx_exp_pr uses a rational function approximation, which is not
       accurate enough at double precision. */
#ifdef GMX_DOUBLE
    maxUlps = 8.0; // empirically determined to be enough on x86
#endif
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
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 2);
    Tester(ReferenceFunctions::gmx_simd_ref_acos_pr,
           TestFunctions::gmx_acos_pr, real(0));
}

TEST_F(SimdFunctionWithOneRealInput, gmx_round_pr_Works)
{
    /* This function is used for rounding of reals before periodicity
     * checks, so we need to check a realistic range of distances. */
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 100);
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
    std::fill(shift.begin(), shift.end(), 0);
    std::fill(scale.begin(), scale.end(), 1024);
    Tester(ReferenceFunctions::gmx_simd_ref_floor_pr,
           TestFunctions::gmx_floor_pr, real(0));
}

#else

/* First, some functions that help us test the way gmx_cvtepi32_pr and
   gmx_cvttpr_epi32 are used in practice. */
namespace ReferenceFunctions
{
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
    std::fill(shift.begin(), shift.end(), -0.25);
    std::fill(scale.begin(), scale.end(), 1024);
    Tester(ReferenceFunctions::gmx_simd_ref_truncate_real_pr,
           TestFunctions::gmx_truncate_real_pr, real(0));
}
#endif

// ====

typedef SimdFunctionTest<2, 1, real> SimdFunctionWithTwoRealInputs;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionWithTwoRealInputs::callFunction(SimdFunctionSet &simdFunctionSet,
                                            FunctionType     function,
                                            real            *_result)
{
    typename SimdFunctionSet::realType a, b, result;

    a      = simdFunctionSet.load_pr(inputs[0]);
    b      = simdFunctionSet.load_pr(inputs[1]);
    result = function(a, b);
    simdFunctionSet.store_pr(_result, result);
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_add_pr_Works)
{
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 5);
    Tester(ReferenceFunctions::gmx_simd_ref_add_pr,
           TestFunctions::gmx_add_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_mul_pr_Works)
{
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 5);
    Tester(ReferenceFunctions::gmx_simd_ref_mul_pr,
           TestFunctions::gmx_mul_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_sub_pr_Works)
{
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 5);
    Tester(ReferenceFunctions::gmx_simd_ref_sub_pr,
           TestFunctions::gmx_sub_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_max_pr_Works)
{
    std::fill(shift.begin(), shift.end(), 0);
    std::fill(scale.begin(), scale.end(), 5);
    Tester(ReferenceFunctions::gmx_simd_ref_max_pr,
           TestFunctions::gmx_max_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_atan2_pr_Works)
{
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 20);
    /* SIMD implementation uses a table lookup, so exact agreement is not
       expected. */
#ifdef GMX_DOUBLE
    maxUlps = 80.0; // empirically determined to be enough on x86
#else
    maxUlps = 10.0; // empirically determined to be enough on x86
#endif
    Tester(ReferenceFunctions::gmx_simd_ref_atan2_pr,
           TestFunctions::gmx_atan2_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_cpsgn_pr_Works)
{
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 5);
    /* gmx_cpsgn_nonneg_pr combines the sign of a (signed) real with
       the value of a real that is non-negative. The implementation
       assumes the latter, so we have to cater for that when
       testing. */
    shift[1] = 0;
    Tester(ReferenceFunctions::gmx_simd_ref_cpsgn_nonneg_pr,
           TestFunctions::gmx_cpsgn_nonneg_pr, real(0));
}

namespace ReferenceFunctions
{

/* Comparison of SIMD reals for "less than", to mimic x86 treatment
   of true and false for this operation. */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_cmplt_pr_like_simd(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pb c = gmx_simd_ref_cmplt_pr(a, b);
    gmx_simd_ref_pr result;

    BitManipulater  manipulater_d;
    /* When QPX is supported, or GMX_X86 exists, the following lines
     * will need to be versioned. */
    UnsignedIntWithSizeOfReal simdTrue  = ~0;
    UnsignedIntWithSizeOfReal simdFalse = 0;
    for (int i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        /* This might need a change to support QPX SIMD, where all
           less than zero are interpreted as FALSE. */
        manipulater_d.i = c.r[i] ? simdTrue : simdFalse;
        result.r[i]     = manipulater_d.r;
    }

    return result;
}

}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_cmplt_pr_Works)
{
    Tester(ReferenceFunctions::gmx_simd_ref_cmplt_pr_like_simd,
           TestFunctions::gmx_cmplt_pr, real(0));
}

// ===

typedef SimdFunctionTest<3, 1, real> SimdFunctionWithThreeRealInputs;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionWithThreeRealInputs::callFunction(SimdFunctionSet &simdFunctionSet,
                                              FunctionType     function,
                                              real            *_result)
{
    typename SimdFunctionSet::realType a, b, c, result;

    a      = simdFunctionSet.load_pr(inputs[0]);
    b      = simdFunctionSet.load_pr(inputs[1]);
    c      = simdFunctionSet.load_pr(inputs[2]);
    result = function(a, b, c);
    simdFunctionSet.store_pr(_result, result);
}

TEST_F(SimdFunctionWithThreeRealInputs, gmx_madd_pr_Works)
{
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 5);
    Tester(ReferenceFunctions::gmx_simd_ref_madd_pr,
           TestFunctions::gmx_madd_pr, real(0));
}

TEST_F(SimdFunctionWithThreeRealInputs, gmx_nmsub_pr_Works)
{
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 5);
    Tester(ReferenceFunctions::gmx_simd_ref_nmsub_pr,
           TestFunctions::gmx_nmsub_pr, real(0));
}

#ifdef GMX_SIMD_HAVE_BLENDV
TEST_F(SimdFunctionWithThreeRealInputs, gmx_blendv_pr_Works)
{
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 5);
    /* gmx_blendv_pr is only used with a positive first argument and a
       second argument of zero. Only the sign of its third argument is
       actually used. So this test covers a wider range than is
       actually used, but that does no real harm. */
    shift[0] = 0;
    scale[1] = 0;
    Tester(ReferenceFunctions::gmx_simd_ref_blendv_pr,
           TestFunctions::gmx_blendv_pr, real(0));
}
#endif

// ===

typedef SimdFunctionTest<1, 2, real> SimdFunctionSinCos;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionSinCos::callFunction(SimdFunctionSet &simdFunctionSet,
                                 FunctionType     function,
                                 real            *_result)
{
    typename SimdFunctionSet::realType a, result[2];

    a      = simdFunctionSet.load_pr(inputs[0]);
    function(a, &result[0], &result[1]);
    for (size_t i = 0; i != numOutputs; ++i)
    {
        simdFunctionSet.store_pr(_result + i * GMX_SIMD_WIDTH_HERE, result[i]);
    }
}

TEST_F(SimdFunctionSinCos, gmx_sincos_pr_Works)
{
    std::fill(shift.begin(), shift.end(), -0.5);
    std::fill(scale.begin(), scale.end(), 4*M_PI);
    /* SIMD implementation uses a table lookup, so exact agreement is
       not expected. Accuracy is worst near multiples of pi/2, when
       one of the output variables is near zero. */
#ifdef GMX_DOUBLE
    maxUlps = 50000.0;
    /* Empirically determined to be enough on x86. This seems like an
       enormous bound, but it's only at double precision and still
       acceptable for its use cases. */
#else
    maxUlps = 10.0;
#endif
    Tester(ReferenceFunctions::gmx_simd_ref_sincos_pr,
           TestFunctions::gmx_sincos_pr, real(0));
}

} // namespace
