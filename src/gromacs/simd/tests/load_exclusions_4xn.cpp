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
 * Tests for 4xn kernel functionality for loading exclusions.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "general.h"
#include "utils.h"

//! TODO this makes more sense as a nbnxn-specific test, rather than a
//! general SIMD test

namespace SIMDTests
{

#ifdef GMX_NBNXN_SIMD_4XN

namespace TestFunctions
{
#include "gromacs/mdlib/nbnxn_kernels/load_4xn_exclusions.h"
#include "gromacs/mdlib/nbnxn_kernels/load_4xn_exclusions_code.h"
}

namespace ReferenceFunctions
{

static gmx_inline void
gmx_simd_ref_load_4xn_exclusions(unsigned int          excl,
                                 gmx_simd_ref_exclmask mask_S0,
                                 gmx_simd_ref_exclmask mask_S1,
                                 gmx_simd_ref_exclmask mask_S2,
                                 gmx_simd_ref_exclmask mask_S3,
                                 gmx_simd_ref_pb      *interact_S0,
                                 gmx_simd_ref_pb      *interact_S1,
                                 gmx_simd_ref_pb      *interact_S2,
                                 gmx_simd_ref_pb      *interact_S3)
{
    /* Load integer interaction mask */
    gmx_simd_ref_exclmask mask_S = gmx_simd_ref_load1_exclmask(excl);
    *interact_S0  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S0);
    *interact_S1  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S1);
    *interact_S2  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S2);
    *interact_S3  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S3);
}

}

typedef SimdFunctionTest<5, 4, unsigned int> SimdFunction_load_4xn_exclusions;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunction_load_4xn_exclusions::callFunction(SimdFunctionSet &simdFunctionSet,
                                               FunctionType     function,
                                               real            *_result)
{
    const size_t numArrayInputs = numInputs - 1;
    typename SimdFunctionSet::exclmaskType masks[numArrayInputs];
    typename SimdFunctionSet::boolType bResults[numArrayInputs];

    for (size_t i = 0; i < numArrayInputs; ++i)
    {
        masks[i] = simdFunctionSet.load_exclmask((unsigned int *) inputs[i]);
    }

    function(inputs[numInputs-1][0],
             masks[0], masks[1],
             masks[2], masks[3],
             &bResults[0], &bResults[1],
             &bResults[2], &bResults[3]);

    for (size_t i = 0; i < numArrayInputs; ++i)
    {
        typename SimdFunctionSet::realType rResult =
            simdFunctionSet.convert_pb_to_pr(bResults[i]);
        simdFunctionSet.store_pr(_result + i * GMX_SIMD_WIDTH_HERE, rResult);
    }
}

TEST_F(SimdFunction_load_4xn_exclusions, gmx_load_4xn_exclusions_Works)
{
    Tester(ReferenceFunctions::gmx_simd_ref_load_4xn_exclusions,
           TestFunctions::gmx_load_4xn_exclusions,
           TestWithSingleBitsSet(0));
}

#endif /* GMX_NBNXN_SIMD_4XN */

}      // namespace
