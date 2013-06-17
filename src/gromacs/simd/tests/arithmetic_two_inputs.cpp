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

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionWithTwoRealInputs::callFunction(SimdFunctionSet &simdFunctionSet,
                                            FunctionType     function,
                                            real            *_result)
{
    typename SimdFunctionSet::realType a, b, result;

    a      = simdFunctionSet.load_pr(inputs_[0]);
    b      = simdFunctionSet.load_pr(inputs_[1]);
    result = function(a, b);
    simdFunctionSet.store_pr(_result, result);
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_add_pr_Works)
{
    prepare(2, 1, -0.5, 5);
    Tester(ReferenceFunctions::gmx_simd_ref_add_pr,
           TestFunctions::gmx_add_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_mul_pr_Works)
{
    prepare(2, 1, -0.5, 5);
    Tester(ReferenceFunctions::gmx_simd_ref_mul_pr,
           TestFunctions::gmx_mul_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_sub_pr_Works)
{
    prepare(2, 1, -0.5, 5);
    Tester(ReferenceFunctions::gmx_simd_ref_sub_pr,
           TestFunctions::gmx_sub_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_max_pr_Works)
{
    prepare(2, 1, 0, 5);
    Tester(ReferenceFunctions::gmx_simd_ref_max_pr,
           TestFunctions::gmx_max_pr, real(0));
}

TEST_F(SimdFunctionWithTwoRealInputs, gmx_atan2_pr_Works)
{
    prepare(2, 1, -0.5, 20);
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
    prepare(2, 1, -0.5, 5);
    /* gmx_cpsgn_nonneg_pr combines the sign of a (signed) real with
       the value of a real that is non-negative. The implementation
       assumes the latter, so we have to cater for that when
       testing. */
    shift_[1] = 0;
    Tester(ReferenceFunctions::gmx_simd_ref_cpsgn_nonneg_pr,
           TestFunctions::gmx_cpsgn_nonneg_pr, real(0));
}

namespace ReferenceFunctions
{

/*! Helper function to do comparison of SIMD reals for "less than", to
    mimic hardware treatment of true and false for this operation. */
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
    prepare(2, 1);
    Tester(ReferenceFunctions::gmx_simd_ref_cmplt_pr_like_simd,
           TestFunctions::gmx_cmplt_pr, real(0));
}

} // namespace
