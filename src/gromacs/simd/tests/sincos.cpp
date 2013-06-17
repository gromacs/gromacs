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
typedef SimdFunctionTest<real> SimdFunctionSinCos;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionSinCos::callFunction(SimdFunctionSet &simdFunctionSet,
                                 FunctionType     function,
                                 real            *_result)
{
    typename SimdFunctionSet::realType a, result[2];

    a      = simdFunctionSet.load_pr(inputs_[0]);
    function(a, &result[0], &result[1]);
    for (size_t i = 0; i != numOutputs_; ++i)
    {
        simdFunctionSet.store_pr(_result + i * GMX_SIMD_WIDTH_HERE, result[i]);
    }
}

TEST_F(SimdFunctionSinCos, gmx_sincos_pr_Works)
{
    prepare(1, 2, -0.5, 4*M_PI);
    /* SIMD implementation uses a table lookup, so exact agreement is
       not expected. Accuracy is worst near multiples of pi/2, when
       one of the output variables is near zero. */
#ifdef GMX_DOUBLE
    maxUlps = 50000.0;
    /* Empirically determined to be enough on x86. This seems like an
       enormous bound, but it's only at double precision and still
       acceptable for its use cases. */
#else
    maxUlps = 10.0;
#endif
    Tester(ReferenceFunctions::gmx_simd_ref_sincos_pr,
           TestFunctions::gmx_sincos_pr, real(0));
}

} // namespace
