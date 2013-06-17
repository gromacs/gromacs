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
 * Tests for SIMD functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "tests.h"

namespace SIMDTests
{

/********************************************************************
 * Tests for SIMD wrapper functions
 */

template<> void
SimdFunctionTest<ReferenceFunction_V_VV, TestFunction_V_VV>::call(ReferenceFunction_V_VV referenceFunction, RealArray theReals)
{
    gmx_simd_ref_pr a, b, result;
    a      = gmx_simd_ref_load_pr(theReals[0]);
    b      = gmx_simd_ref_load_pr(theReals[1]);
    result = referenceFunction(a, b);
    gmx_simd_ref_store_pr(referenceResult[0], result);
}

template<> void
SimdFunctionTest<ReferenceFunction_V_VV, TestFunction_V_VV>::call(TestFunction_V_VV testFunction, RealArray theReals)
{
    gmx_mm_pr a, b, result;
    a      = gmx_load_pr(theReals[0]);
    b      = gmx_load_pr(theReals[1]);
    result = testFunction(a, b);
    gmx_store_pr(testResult[0], result);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_V_VV, TestFunction_V_VV> SimdFunctionWithSignature_V_VV;

TEST_F(SimdFunctionWithSignature_V_VV, gmx_add_pr_Works)
{
    Tester(gmx_simd_ref_add_pr,
           gmx_add_pr,
           reals);
}

TEST_F(SimdFunctionWithSignature_V_VV, gmx_mul_pr_Works)
{
    Tester(gmx_simd_ref_mul_pr,
           gmx_mul_pr,
           reals);
}

TEST_F(SimdFunctionWithSignature_V_VV, gmx_sub_pr_Works)
{
    Tester(gmx_simd_ref_sub_pr,
           gmx_sub_pr,
           reals);
}

TEST_F(SimdFunctionWithSignature_V_VV, gmx_max_pr_Works)
{
    Tester(gmx_simd_ref_max_pr,
           gmx_max_pr,
           reals);
}

TEST_F(SimdFunctionWithSignature_V_VV, gmx_atan2_pr_Works)
{
    Tester(gmx_simd_ref_atan2_pr,
           gmx_atan2_pr,
           reals
#ifdef GMX_DOUBLE
           , 10.0 // empirically determined to be enough on x86
#endif
           );
}

TEST_F(SimdFunctionWithSignature_V_VV, gmx_cpsgn_pr_Works)
{
    /* gmx_cpsgn_nonneg_pr combines the sign of a (signed) real with
       the value of a real that is non-negative, so the test case
       should reproduce that situation. Since positiveReals contains
       non-negative real values, this is easy. */
    for (int i = 0; i < SIMD_TEST_NUM_REAL_VARIABLES; i += SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST)
    {
        for (int j = 0; j < GMX_SIMD_WIDTH_HERE; ++j)
        {
            reals[i+1][j] = positiveReals[i+1][j];
        }
    }
    Tester(gmx_simd_ref_cpsgn_nonneg_pr,
           gmx_cpsgn_nonneg_pr,
           reals);
}

} // namespace
