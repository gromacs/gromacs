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
typedef SimdFunctionTest<real> SimdFunctionWithTwoRealInputs;

template<class SimdFunctionSet,
         typename FunctionType>
void
callFunction(SimdFunctionSet     &simdFunctionSet,
             FunctionType         function,
             std::vector<real*>   inputs,
             real                *result)
{
    typename SimdFunctionSet::realType a, b, result_pr;

    a         = simdFunctionSet.load_pr(inputs[0]);
    b         = simdFunctionSet.load_pr(inputs[1]);
    result_pr = function(a, b);
    simdFunctionSet.store_pr(result, result_pr);
}

/* Wrapper function because MSVC can't cope otherwise. */
static gmx_inline gmx_mm_pr gmx_add_pr_wrapper(gmx_mm_pr a, gmx_mm_pr b)
{
    return gmx_add_pr(a, b);
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_add_pr_Works)
{
    prepare(2, 1, -0.5, 5);
    RunTest(ReferenceFunctions::gmx_simd_ref_add_pr,
            gmx_add_pr_wrapper, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_mul_pr_Works)
{
    prepare(2, 1, -0.5, 5);
    RunTest(ReferenceFunctions::gmx_simd_ref_mul_pr,
            gmx_mul_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_sub_pr_Works)
{
    prepare(2, 1, -0.5, 5);
    RunTest(ReferenceFunctions::gmx_simd_ref_sub_pr,
            gmx_sub_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_max_pr_Works)
{
    prepare(2, 1, 0, 5);
    RunTest(ReferenceFunctions::gmx_simd_ref_max_pr,
            gmx_max_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_atan2_pr_Works)
{
    prepare(2, 1, -0.5, 20);
    /* SIMD implementation uses a table lookup, so exact agreement is not
       expected. */
    // TODO Actually, if the ratio is nearly 1 then errors of 1 in 100
    // have been seen. Either need a tunable limit, or a different
    // reference implementation.
#ifdef GMX_DOUBLE
    maxUlps_ = 80.0; // empirically determined to be enough on x86
#else
    maxUlps_ = 10.0; // empirically determined to be enough on x86
#endif
    RunTest(ReferenceFunctions::gmx_simd_ref_atan2_pr,
            gmx_atan2_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_cpsgn_pr_Works)
{
    prepare(2, 1, -0.5, 5);
    /* gmx_cpsgn_nonneg_pr combines the sign of a (signed) real with
       the value of a real that is non-negative. The implementation
       assumes the latter, so we have to cater for that when
       testing. */
    shift_[1] = 0;
    RunTest(ReferenceFunctions::gmx_simd_ref_cpsgn_nonneg_pr,
            gmx_cpsgn_nonneg_pr, real(0));
}

namespace ReferenceFunctions
{

/*! Helper function to do comparison of SIMD reals for "less than",
    that also mimics the SIMD hardware treatment of logical true and
    false. See logic.cpp for more details. */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_cmplt_wrapper_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr result;
    gmx_simd_ref_pb result_pb;

    result_pb = gmx_simd_ref_cmplt_pr(a, b);

    for (int i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        result.r[i] = TestFunctions::HardwareSimdFunctionSet::getSimdBool(result_pb.r[i]);
    }

    return result;
}

}

/* Wrapper function because the macro wrapped is defined with
   parameters. */
static gmx_inline gmx_mm_pr gmx_cmplt_pr_wrapper(gmx_mm_pr a, gmx_mm_pr b)
{
    return gmx_cmplt_pr(a, b);
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_cmplt_pr_Works)
{
    prepare(2, 1);
    RunTest(ReferenceFunctions::gmx_simd_ref_cmplt_wrapper_pr,
            gmx_cmplt_pr_wrapper, real(0));
}

} // namespace
