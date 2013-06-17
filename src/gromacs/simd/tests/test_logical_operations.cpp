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

/* Ideally tests for SIMD logical operations would load values into
 * SIMD registers and test them. However, there is no SIMD load of a
 * vector of bool, because there is no need for that in the kernels.
 * There is also the issue that all current (i.e. x86) implementations
 * of these functions are actually done as bitwise operations on
 * register types that are identical to the floating-point ones. So we
 * can recycle the V_VV testing machinery. If/when these conditions
 * change, these tests will need updating. In particular, a proper
 * B_BB or B_VV function signature could become necessary. */

/* Test "B_BB" functions */

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_V_VV, TestFunction_V_VV> SimdFunctionWithSignature_V_VV;

/* Bit-wise AND on SIMD reals, to mimic x86. The official
   gmx_simd_ref_and_pb does a logical AND, so we can't use that
   for testing the actual SIMD implementations. */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_bitwise_and_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr result;

    BitManipulater  manipulater_a, manipulater_b, manipulater_c;
    for (int i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        manipulater_a.r = a.r[i];
        manipulater_b.r = b.r[i];
        manipulater_c.i = manipulater_a.i & manipulater_b.i;
        result.r[i]     = manipulater_c.r;
    }

    return result;
}

TEST_F(SimdFunctionWithSignature_V_VV, gmx_and_pb_Works)
{
    Tester(gmx_simd_ref_bitwise_and_pr,
           gmx_and_pb,
           booleanReals);
}

/* Bit-wise OR on SIMD reals, to mimic x86. The official
   gmx_simd_ref_and_pb does a logical OR, so we can't use that for
   testing the actual SIMD implementations. */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_bitwise_or_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr result;

    BitManipulater  manipulater_a, manipulater_b, manipulater_c;
    for (int i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        manipulater_a.r = a.r[i];
        manipulater_b.r = b.r[i];
        manipulater_c.i = manipulater_a.i | manipulater_b.i;
        result.r[i]     = manipulater_c.r;
    }

    return result;
}

TEST_F(SimdFunctionWithSignature_V_VV, gmx_or_pb_Works)
{
    Tester(gmx_simd_ref_bitwise_or_pr,
           gmx_or_pb,
           booleanReals);
}

/* Test "B_VV" functions */

/* Comparison of SIMD reals for "less than", to mimic x86 treatment
   of true and false for this operation. */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_cmplt_pr_like_simd(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pb c = gmx_simd_ref_cmplt_pr(a, b);
    gmx_simd_ref_pr result;

    BitManipulater  manipulater_d;
    for (int i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        /* This might need a change to support QPX SIMD, where all
           less than zero are interpreted as FALSE. */
        manipulater_d.i = c.r[i] ? (UnsignedIntWithSizeOfReal) ~0 : 0;
        result.r[i]     = manipulater_d.r;
    }

    return result;
}

TEST_F(SimdFunctionWithSignature_V_VV, gmx_cmplt_pr_Works)
{
    Tester(gmx_simd_ref_cmplt_pr_like_simd,
           gmx_cmplt_pr,
           reals);
}

/* Test I_V functions */

#ifdef GMX_SIMD_HAVE_ANYTRUE
template<> void
SimdFunctionTest<ReferenceFunction_I_V, TestFunction_I_V>::call(ReferenceFunction_I_V referenceFunction, RealArray theReals)
{
    gmx_simd_ref_pr a;

    a = gmx_simd_ref_load_pr(theReals[0]);
    int result = referenceFunction(a);
    referenceResult[0][0] = (real) result;
}

template<> void
SimdFunctionTest<ReferenceFunction_I_V, TestFunction_I_V>::call(TestFunction_I_V testFunction, RealArray theReals)
{
    gmx_mm_pr a;
    a = gmx_load_pr(theReals[0]);
    int       result = testFunction(a);
    testResult[0][0] = (real) (0 != result);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_I_V, TestFunction_I_V> SimdFunctionWithSignature_I_B;

/* gmx_anytrue_pb reads a gmx_mm_pb, but all the implementations
 * (i.e. x86) just do an implicit conversion of SIMD real to
 * SIMD bool. To test, we have to do the same thing. We could
 * do this by constructing an array of gmx_simd_ref_pb, but so
 * far it is only needed here.
 */
static gmx_inline int
gmx_simd_ref_anytrue_pr(gmx_simd_ref_pr a)
{
    int             i;
    gmx_simd_ref_pb b;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        /* This might need a change to support QPX SIMD, where all
           less than zero are interpreted as FALSE. */
        b.r[i] = (0.0 != a.r[i]);
    }

    return gmx_simd_ref_anytrue_pb(b);
}

TEST_F(SimdFunctionWithSignature_I_B, gmx_anytrue_pb_Works)
{
    Tester(gmx_simd_ref_anytrue_pr,
           gmx_anytrue_pb,
           booleanReals);
}
#endif

} // namespace
