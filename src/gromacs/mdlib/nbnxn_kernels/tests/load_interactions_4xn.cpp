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
 * Tests for 4xn kernel functionality for loading interaction masks.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "load_interactions.h"

#ifdef GMX_NBNXN_SIMD_4XN

namespace SIMDTests
{

using ReferenceFunctions::gmx_simd_ref_exclfilter;
using ReferenceFunctions::gmx_simd_ref_pb;
using ReferenceFunctions::gmx_simd_ref_load1_exclfilter;
using ReferenceFunctions::gmx_simd_ref_checkbitmask_pb;
using ReferenceFunctions::gmx_simd_ref_load_interaction_mask_pb;

static gmx_inline void
gmx_simd_ref_load_interactions_4xn(unsigned int            excl,
                                   gmx_simd_ref_exclfilter filter_S0,
                                   gmx_simd_ref_exclfilter filter_S1,
                                   gmx_simd_ref_exclfilter filter_S2,
                                   gmx_simd_ref_exclfilter filter_S3,
                                   const char             *interaction_mask_indices,
                                   real                   *simd_interaction_array,
                                   gmx_simd_ref_pb        *interact_S0,
                                   gmx_simd_ref_pb        *interact_S1,
                                   gmx_simd_ref_pb        *interact_S2,
                                   gmx_simd_ref_pb        *interact_S3)
{
    /* Load integer interaction mask */
    gmx_simd_ref_exclfilter mask_S_pr = gmx_simd_ref_load1_exclfilter(excl);
    *interact_S0  = gmx_simd_ref_checkbitmask_pb(mask_S_pr, filter_S0);
    *interact_S1  = gmx_simd_ref_checkbitmask_pb(mask_S_pr, filter_S1);
    *interact_S2  = gmx_simd_ref_checkbitmask_pb(mask_S_pr, filter_S2);
    *interact_S3  = gmx_simd_ref_checkbitmask_pb(mask_S_pr, filter_S3);
}

//! Typedef for the test fixture
typedef SimdFunctionTest_kernel<unsigned> SimdFunction_load_interactions_4xn;

template<class SimdFunctionSet,
         class SimdReal4FunctionSet,
         typename FunctionType>
void
callFunctionLocal(SimdFunctionSet           &simdFunctionSet,
                  SimdReal4FunctionSet      &simdReal4FunctionSet,
                  FunctionType               function,
                  int                        filter_stride,
                  int                        unrollj,
                  std::vector<unsigned int*> inputs,
                  real                      *result)
{
    typename SimdReal4FunctionSet::exclfilterType *masks;
    snew_aligned(masks, UNROLLI, GMX_SIMD_WIDTH_HERE * sizeof(real));
    typename SimdFunctionSet::boolType *bResults;
    snew_aligned(bResults, UNROLLI, GMX_SIMD_WIDTH_HERE * sizeof(real));
    char              interaction_mask_indices[GMX_SIMD_WIDTH_HERE];

    nbnxn_atomdata_t *nbat = NULL;
    prepareToReadExclusions(simdReal4FunctionSet,
                            filter_stride,
                            unrollj,
                            inputs[0][0],
                            nbat,
                            masks,
                            interaction_mask_indices);

    function(inputs[0][0],
             masks[0], masks[1],
             masks[2], masks[3],
             interaction_mask_indices,
#ifdef GMX_CPU_ACCELERATION_IBM_QPX
             nbat->simd_interaction_array,
#else
             NULL, // TODO does this work?
#endif
             &bResults[0], &bResults[1],
             &bResults[2], &bResults[3]);

    for (unsigned i = 0; i < UNROLLI; i++)
    {
        typename SimdFunctionSet::realType rResult =
            simdReal4FunctionSet.convert_pb_to_pr(bResults[i]);
        simdFunctionSet.store_pr(result + i * GMX_SIMD_WIDTH_HERE, rResult);
    }

    cleanUpFromReadingExclusions<SimdReal4FunctionSet>(bResults, nbat, masks);
}

// ===

TEST_F(SimdFunction_load_interactions_4xn, gmx_load_interactions_4xn_Works)
{
    prepare(1, UNROLLI);
    RunTest(gmx_simd_ref_load_interactions_4xn,
            gmx_load_simd_4xn_interactions,
            filter_stride,
            UNROLLJ,
            referenceSimd_,
            referenceSimdReal4_,
            hardwareSimd_,
            hardwareSimdReal4_,
            TestWith32BitInt(0));
}

}

#endif
