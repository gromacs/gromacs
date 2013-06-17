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
SimdFunctionTest<ReferenceFunction_V_VVV, TestFunction_V_VVV>::call(ReferenceFunction_V_VVV referenceFunction, RealArray theReals)
{
    gmx_simd_ref_pr a, b, c, result;
    a      = gmx_simd_ref_load_pr(theReals[0]);
    b      = gmx_simd_ref_load_pr(theReals[1]);
    c      = gmx_simd_ref_load_pr(theReals[2]);
    result = referenceFunction(a, b, c);
    gmx_simd_ref_store_pr(referenceResult[0], result);
}

template<> void
SimdFunctionTest<ReferenceFunction_V_VVV, TestFunction_V_VVV>::call(TestFunction_V_VVV testFunction, RealArray theReals)
{
    gmx_mm_pr a, b, c, result;
    a      = gmx_load_pr(theReals[0]);
    b      = gmx_load_pr(theReals[1]);
    c      = gmx_load_pr(theReals[2]);
    result = testFunction(a, b, c);
    gmx_store_pr(testResult[0], result);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_V_VVV, TestFunction_V_VVV> SimdFunctionWithSignature_V_VVV;

TEST_F(SimdFunctionWithSignature_V_VVV, gmx_madd_pr_Works)
{
    Tester(gmx_simd_ref_madd_pr,
           gmx_madd_pr,
           reals);
}

TEST_F(SimdFunctionWithSignature_V_VVV, gmx_nmsub_pr_Works)
{
    Tester(gmx_simd_ref_nmsub_pr,
           gmx_nmsub_pr,
           reals);
}

#ifdef GMX_SIMD_HAVE_BLENDV
TEST_F(SimdFunctionWithSignature_V_VVV, gmx_blendv_pr_Works)
{
    Tester(gmx_simd_ref_blendv_pr,
           gmx_blendv_pr,
           reals);
}
#endif

} // namespace
