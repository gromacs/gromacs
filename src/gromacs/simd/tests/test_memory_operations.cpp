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

/* It's hard to test load and store independently, since you need one
   of them to be in a position to test the other. So probably these
   sets of tests will pass/fail in tandem, but if they fail probably
   every other test will be failing also. */

template<> void
SimdFunctionTest<ReferenceFunction_V_C, TestFunction_V_C>::call(ReferenceFunction_V_C referenceFunction, RealArray theReals)
{
    gmx_simd_ref_pr result;
    result = referenceFunction(theReals[0]);
    gmx_simd_ref_store_pr(referenceResult[0], result);
}

template<> void
SimdFunctionTest<ReferenceFunction_V_C, TestFunction_V_C>::call(TestFunction_V_C testFunction, RealArray theReals)
{
    gmx_mm_pr result;
    result = testFunction(theReals[0]);
    gmx_store_pr(testResult[0], result);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_V_C, TestFunction_V_C> SimdFunctionWithSignature_V_C;

TEST_F(SimdFunctionWithSignature_V_C, gmx_load_pr_Works)
{
    Tester(gmx_simd_ref_load_pr,
           gmx_load_pr,
           reals);
}

TEST_F(SimdFunctionWithSignature_V_C, gmx_load1_pr_Works)
{
    Tester(gmx_simd_ref_load1_pr,
           gmx_load1_pr,
           reals);
}

/* Test other "load" functions */

template<> void
SimdFunctionTest<ReferenceFunction_V_R, TestFunction_V_R>::call(ReferenceFunction_V_R referenceFunction, RealArray referenceReals)
{
    gmx_simd_ref_pr result;
    result = referenceFunction(referenceReals[0][0]);
    gmx_simd_ref_store_pr(referenceResult[0], result);
}

template<> void
SimdFunctionTest<ReferenceFunction_V_R, TestFunction_V_R>::call(TestFunction_V_R testFunction, RealArray testReals)
{
    gmx_mm_pr result;
    result = testFunction(testReals[0][0]);
    gmx_store_pr(testResult[0], result);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_V_R, TestFunction_V_R> SimdFunctionWithSignature_V_R;

TEST_F(SimdFunctionWithSignature_V_R, gmx_set1_pr_Works)
{
    Tester(gmx_simd_ref_set1_pr,
           gmx_set1_pr,
           reals);
}

template<> void
SimdFunctionTest<ReferenceFunction_V_none, TestFunction_V_none>::call(ReferenceFunction_V_none referenceFunction, RealArray)
{
    gmx_simd_ref_pr result;
    result = referenceFunction();
    gmx_simd_ref_store_pr(referenceResult[0], result);
}

template<> void
SimdFunctionTest<ReferenceFunction_V_none, TestFunction_V_none>::call(TestFunction_V_none testFunction, RealArray)
{
    gmx_mm_pr result;
    result = testFunction();
    gmx_store_pr(testResult[0], result);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_V_none, TestFunction_V_none> SimdFunctionWithSignature_V_none;

TEST_F(SimdFunctionWithSignature_V_none, gmx_setzero_pr_Works)
{
    Tester(gmx_simd_ref_setzero_pr,
           gmx_setzero_pr,
           reals);
}

/* Test store functions */

template<> void
SimdFunctionTest<ReferenceFunction_none_QV, TestFunction_none_QV>::call(ReferenceFunction_none_QV referenceFunction, RealArray theReals)
{
    gmx_simd_ref_pr a;
    a = gmx_simd_ref_load_pr(theReals[0]);
    referenceFunction(referenceResult[0], a);
}

template<> void
SimdFunctionTest<ReferenceFunction_none_QV, TestFunction_none_QV>::call(TestFunction_none_QV testFunction, RealArray theReals)
{
    gmx_mm_pr a;
    a = gmx_load_pr(theReals[0]);
    testFunction(testResult[0], a);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_none_QV, TestFunction_none_QV> SimdFunctionWithSignature_none_QV;

TEST_F(SimdFunctionWithSignature_none_QV, gmx_store_pr_Works)
{
    Tester(gmx_simd_ref_store_pr,
           gmx_store_pr,
           reals);
}

} // namespace
