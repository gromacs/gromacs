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
#include <string>

#include "tests.h"

namespace SIMDTests
{

/********************************************************************
 * Tests for SIMD wrapper functions
 */

template<> void
SimdFunctionTest<ReferenceFunction_V_V, TestFunction_V_V>::call(ReferenceFunction_V_V referenceFunction, RealArray referenceReals)
{
    gmx_simd_ref_pr a, result;

    a      = gmx_simd_ref_load_pr(referenceReals[0]);
    result = referenceFunction(a);
    gmx_simd_ref_store_pr(referenceResult[0], result);
}

template<> void
SimdFunctionTest<ReferenceFunction_V_V, TestFunction_V_V>::call(TestFunction_V_V testFunction, RealArray testReals)
{
    gmx_mm_pr a, result;

    a      = gmx_load_pr(testReals[0]);
    result = testFunction(a);
    gmx_store_pr(testResult[0], result);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_V_V, TestFunction_V_V> SimdFunctionWithSignature_V_V;

TEST_F(SimdFunctionWithSignature_V_V, gmx_round_pr_Works)
{
    Tester(gmx_simd_ref_round_pr,
           gmx_round_pr,
           reals);
}

/* SIMD reciprocal square roots are computed by an approximate SIMD
   operation, and then some Newton-Raphson iterations depending on the
   accuracy of the underlying hardware operation and the prevailing
   GROMACS single/double precision.

   1) The SIMD rsqrt is never used without subsequent Newton-Raphson
   iteration, and there is no reference version of approximate inverse
   square root, so it only makes sense to test the composite
   operation.

   2) GROMACS does not do enough Newton-Raphson iterations for full
   double precision (because we don't need it), so the Google Test
   standard criterion does not work in that case. */
TEST_F(SimdFunctionWithSignature_V_V, gmx_rsqrt_pr_Works)
{
    Tester(gmx_simd_ref_rsqrt_pr,
           gmx_invsqrt_pr,
           positiveReals
#ifdef GMX_DOUBLE
           , 8.0 // empirically determined to be enough on x86
#endif
           );
}

/* The SIMD reciprocal is never used without subsequent Newton-Raphson
   iteration, and there is no reference version of approximate
   reciprocal, so it only makes sense to test the composite
   operation. Similarly, the standard Google Test accuracy requirement
   is too strict in this case. */
TEST_F(SimdFunctionWithSignature_V_V, gmx_inv_pr_Works)
{
    Tester(gmx_simd_ref_rcp_pr,
           gmx_inv_pr,
           reals
#ifdef GMX_DOUBLE
           , 2.0 // empirically determined to be enough on x86
#endif
           );
}

TEST_F(SimdFunctionWithSignature_V_V, gmx_sqrt_pr_Works)
{
    Tester(gmx_simd_ref_sqrt_pr,
           gmx_sqrt_pr,
           positiveReals);
}

/* gmx_exp_pr uses a rational function approximation, which is not
   quite accurate enough at double precision. */
TEST_F(SimdFunctionWithSignature_V_V, gmx_exp_pr_Works)
{
    Tester(gmx_simd_ref_exp_pr,
           gmx_exp_pr,
           reals
#ifdef GMX_DOUBLE
           , 2.0
#endif
           );
}

TEST_F(SimdFunctionWithSignature_V_V, gmx_acos_pr_Works)
{
    Tester(gmx_simd_ref_acos_pr,
           gmx_acos_pr,
           reals);
}

/* The set of functions that we'd like to test vary with the kind of
   acceleration being used. */
#ifdef GMX_SIMD_HAVE_FLOOR

TEST_F(SimdFunctionWithSignature_V_V, gmx_floor_pr_Works)
{
    Tester(gmx_simd_ref_floor_pr,
           gmx_floor_pr,
           reals);
}

#else

/* First, some functions that help us test the way gmx_cvtepi32_pr and
   gmx_cvttpr_epi32 are used in practice. */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_truncate_real_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        b.r[i] = (real)(int)(a.r[i]) * 2;
    }

    return b;
}

static gmx_inline gmx_mm_pr gmx_truncate_real_pr(gmx_mm_pr a)
{
    return gmx_cvtepi32_pr(gmx_cvttpr_epi32(a));
}

/* This tests both gmx_cvttpr_epi32 and gmx_cvtepi32_pr. These two are
   inconvenient to test seperately. The latter is only ever used after
   the former, and only when gmx_floor_pr is unavailable. The former
   is also used before table lookups, which are tested next. */

TEST_F(SimdFunctionWithSignature_V_V, PairDoingTruncationOfReals_Works)
{
    Tester(gmx_simd_ref_truncate_real_pr,
           gmx_truncate_real_pr,
           reals);
}
#endif

} // namespace
