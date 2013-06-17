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
SimdFunctionTest<ReferenceFunction_I_VPP, TestFunction_I_VPP>::call(ReferenceFunction_I_VPP referenceFunction, RealArray theReals)
{
    gmx_simd_ref_pr a, b, c;

    a = gmx_simd_ref_load_pr(theReals[0]);
    referenceFunction(a, &b, &c);
    gmx_simd_ref_store_pr(referenceResult[0], b);
    gmx_simd_ref_store_pr(referenceResult[1], c);
}

template<> void
SimdFunctionTest<ReferenceFunction_I_VPP, TestFunction_I_VPP>::call(TestFunction_I_VPP testFunction, RealArray theReals)
{
    gmx_mm_pr a, b, c;
    a = gmx_load_pr(theReals[0]);
    testFunction(a, &b, &c);
    gmx_store_pr(testResult[0], b);
    gmx_store_pr(testResult[1], c);
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_I_VPP, TestFunction_I_VPP> SimdFunctionWithSignature_I_VPP;

TEST_F(SimdFunctionWithSignature_I_VPP, gmx_sincos_pr_Works)
{
    Tester(gmx_simd_ref_sincos_pr,
           gmx_sincos_pr,
           reals,
#ifdef GMX_DOUBLE
           2.0, // empirically determined to be enough on x86
#else
           1.0,
#endif
           2); // Make sure we test both outputs
}

} // namespace
