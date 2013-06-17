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

/* We can't actually do a SIMD load of a vector of bool, because it is
 * not implemented since there is no need for that in the
 * code. Instead, a vector of bool is always generated from a
 * comparison, so we do that in testing also. gmx_cmplt_pr is tested
 * separately. */

/* Test V_VB functions */

template<> void
SimdFunctionTest<ReferenceFunction_V_VB, TestFunction_V_VB>::call(ReferenceFunction_V_VB referenceFunction, RealArray theReals)
{
    gmx_simd_ref_pr a, b, c, result;
    gmx_simd_ref_pb mask;
    a      = gmx_simd_ref_load_pr(theReals[0]);
    b      = gmx_simd_ref_load_pr(theReals[1]);
    c      = gmx_simd_ref_load_pr(theReals[2]);
    mask   = gmx_simd_ref_cmplt_pr(b, c);
    result = referenceFunction(a, mask);
    gmx_simd_ref_store_pr(referenceResult[0], result);
}

template<> void
SimdFunctionTest<ReferenceFunction_V_VB, TestFunction_V_VB>::call(TestFunction_V_VB testFunction, RealArray theReals)
{
    gmx_mm_pr a, b, c, result;
    gmx_mm_pb mask;
    a         = gmx_load_pr(theReals[0]);
    b         = gmx_load_pr(theReals[1]);
    c         = gmx_load_pr(theReals[2]);
    mask      = gmx_cmplt_pr(b, c);
    result    = testFunction(a, mask);
    gmx_store_pr(testResult[0], result);
}

typedef SimdFunctionTest<ReferenceFunction_V_VB, TestFunction_V_VB> SimdFunctionWithSignature_V_VB;

TEST_F(SimdFunctionWithSignature_V_VB, gmx_blendzero_pr_Works)
{
    Tester(gmx_simd_ref_blendzero_pr,
           gmx_blendzero_pr,
           reals);
}

/* Test V_BVV functions */

template<> void
SimdFunctionTest<ReferenceFunction_V_BVV, TestFunction_V_BVV>::call(ReferenceFunction_V_BVV referenceFunction, RealArray theReals)
{
    gmx_simd_ref_pr a, b, c, d, result;
    gmx_simd_ref_pb mask;
    a      = gmx_simd_ref_load_pr(theReals[0]);
    b      = gmx_simd_ref_load_pr(theReals[1]);
    c      = gmx_simd_ref_load_pr(theReals[2]);
    d      = gmx_simd_ref_load_pr(theReals[3]);
    mask   = gmx_simd_ref_cmplt_pr(c, d);
    result = referenceFunction(mask, a, b);
    gmx_simd_ref_store_pr(referenceResult[0], result);
}

template<> void
SimdFunctionTest<ReferenceFunction_V_BVV, TestFunction_V_BVV>::call(TestFunction_V_BVV testFunction, RealArray theReals)
{
    gmx_mm_pr a, b, c, d, result;
    gmx_mm_pb mask;
    a      = gmx_load_pr(theReals[0]);
    b      = gmx_load_pr(theReals[1]);
    c      = gmx_load_pr(theReals[2]);
    d      = gmx_load_pr(theReals[3]);
    mask   = gmx_cmplt_pr(c, d);
    result = testFunction(mask, a, b);
    gmx_store_pr(testResult[0], result);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_V_BVV, TestFunction_V_BVV> SimdFunctionWithSignature_V_BVV;

TEST_F(SimdFunctionWithSignature_V_BVV, gmx_masknot_add_pr_Works)
{
    Tester(gmx_simd_ref_masknot_add_pr,
           gmx_masknot_add_pr,
           reals);
}

} // namespace
