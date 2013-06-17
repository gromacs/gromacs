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
typedef SimdFunctionTest<real> SimdFunctionWithThreeRealInputs;

template<class SimdFunctionSet,
         typename FunctionType>
void
callFunction(SimdFunctionSet     &simdFunctionSet,
             FunctionType         function,
             std::vector<real*>   inputs,
             real                *result)
{
    typename SimdFunctionSet::realType a, b, c, result_pr;

    a         = simdFunctionSet.load_pr(inputs[0]);
    b         = simdFunctionSet.load_pr(inputs[1]);
    c         = simdFunctionSet.load_pr(inputs[2]);
    result_pr = function(a, b, c);
    simdFunctionSet.store_pr(result, result_pr);
}

/* Wrapper function because the macro wrapped is defined with
   parameters. */
static gmx_inline gmx_mm_pr gmx_madd_pr_wrapper(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c)
{
    return gmx_madd_pr(a, b, c);
}

TEST_F(SimdFunctionWithThreeRealInputs, gmx_madd_pr_Works)
{
    prepare(3, 1, -0.5, 5);
    RunTest(ReferenceFunctions::gmx_simd_ref_madd_pr,
            gmx_madd_pr_wrapper, real(0));
}

/* Wrapper function because the macro wrapped is defined with
   parameters. */
static gmx_inline gmx_mm_pr gmx_nmsub_pr_wrapper(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c)
{
    return gmx_nmsub_pr(a, b, c);
}

TEST_F(SimdFunctionWithThreeRealInputs, gmx_nmsub_pr_Works)
{
    prepare(3, 1, -0.5, 5);
    RunTest(ReferenceFunctions::gmx_simd_ref_nmsub_pr,
            gmx_nmsub_pr_wrapper, real(0));
}

#ifdef GMX_SIMD_HAVE_BLENDV
/* Wrapper function because the macro wrapped is defined with
   parameters. */
static gmx_inline gmx_mm_pr gmx_blendv_pr_wrapper(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c)
{
    return gmx_blendv_pr(a, b, c);
}

TEST_F(SimdFunctionWithThreeRealInputs, gmx_blendv_pr_Works)
{
    prepare(3, 1, -0.5, 5);
    /* gmx_blendv_pr is only used with a positive first argument and a
       second argument of zero. Only the sign of its third argument is
       actually used. So this test covers a wider range than is
       actually used, but that does no real harm. */
    shift_[0] = 0;
    scale_[1] = 0;
    RunTest(ReferenceFunctions::gmx_simd_ref_blendv_pr,
            gmx_blendv_pr_wrapper, real(0));
}
#endif

} // namespace
