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

#if (defined GMX_NBNXN_SIMD_2XNN) || (defined GMX_NBNXN_SIMD_4XN)

#ifdef GMX_NBNXN_SIMD_2XNN
#include "gromacs/mdlib/nbnxn_kernels/load_2xnn_exclusions.h"
#include "gromacs/mdlib/nbnxn_kernels/load_2xnn_exclusions_code.h"
#endif

#ifdef GMX_NBNXN_SIMD_4XN
#include "gromacs/mdlib/nbnxn_kernels/load_4xn_exclusions.h"
#include "gromacs/mdlib/nbnxn_kernels/load_4xn_exclusions_code.h"
#endif

#include "tests.h"

namespace SIMDTests
{

/* Helper function to make it possible to compare the results of the
 * reference code (which produces a boolean value) and the SIMD code
 * (which produces something hardware-specific). */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_convert_pb_to_pr(gmx_simd_ref_pb src)
{
    gmx_simd_ref_pr dest;
    int i;
    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        BitManipulater manipulater;
        /* Might need a new version for QPX */
        UnsignedIntWithSizeOfReal simdTrue = ~0;
        UnsignedIntWithSizeOfReal simdFalse = 0;
        manipulater.i = src.r[i] ? simdTrue : simdFalse;
        dest.r[i] = manipulater.r;
    }

    return dest;
}

/********************************************************************
 * Tests for SIMD wrapper functions
 */

#ifdef GMX_NBNXN_SIMD_2XNN

static inline void
gmx_simd_ref_load_2xnn_exclusions(unsigned int excl,
                                  gmx_simd_ref_exclmask mask_S0,
                                  gmx_simd_ref_exclmask mask_S2,
                                  gmx_simd_ref_pb *interact_S0,
                                  gmx_simd_ref_pb *interact_S2)
{
    /* Load integer interaction mask */
    gmx_simd_ref_exclmask mask_S = gmx_simd_ref_load1_exclmask(excl);
    *interact_S0  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S0);
    *interact_S2  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S2);
}

/* Test "none_IEEDD" function */

const int numMasks_IEEDD = 2;

template<> void
SimdFunctionTest<ReferenceFunction_none_IEEDD, TestFunction_none_IEEDD>::call(ReferenceFunction_none_IEEDD referenceFunction, const UIntArray theUInts)
{
    gmx_simd_ref_exclmask masks[numMasks_IEEDD];
    gmx_simd_ref_pb bResults[numMasks_IEEDD];
    gmx_simd_ref_pr rResults[numMasks_IEEDD];

    for(int i = 0; i < numMasks_IEEDD; ++i)
    {
        masks[i] = gmx_simd_ref_load_exclmask(theUInts[i]);
    }
    referenceFunction(theUInts[numMasks_IEEDD][0],
                      masks[0], masks[1],
                      &bResults[0], &bResults[1]);
    for(int i = 0; i < numMasks_IEEDD; ++i)
    {
        rResults[i] = gmx_simd_ref_convert_pb_to_pr(bResults[i]);
        gmx_simd_ref_store_pr(referenceResult[i], rResults[i]);
    }
}

template<> void
SimdFunctionTest<ReferenceFunction_none_IEEDD, TestFunction_none_IEEDD>::call(TestFunction_none_IEEDD testFunction, const UIntArray theUInts)
{
    gmx_exclmask masks[numMasks_IEEDD];
    gmx_mm_pb results[numMasks_IEEDD];

    for(int i = 0; i < numMasks_IEEDD; ++i)
    {
        masks[i] = gmx_load_exclmask(theUInts[i]);
    }
    testFunction(theUInts[numMasks_IEEDD][0],
                 masks[0], masks[1],
                 &results[0], &results[1]);
    for(int i = 0; i < numMasks_IEEDD; ++i)
    {
        /* On x86 and QPX, gmx_mm_bb is the same type as gmx_mm_pr, so
           this works */
        gmx_store_pr(testResult[i], results[i]);
    }
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_none_IEEDD, TestFunction_none_IEEDD> SimdFunctionWithSignature_none_IEEDD;

TEST_F(SimdFunctionWithSignature_none_IEEDD, gmx_load_2xnn_exclusions_Works)
{
    Tester(gmx_simd_ref_load_2xnn_exclusions,
           gmx_load_2xnn_exclusions,
           bits,
           1.0,
           numMasks_IEEDD);
}

#endif /* GMX_NBNXN_SIMD_2XNN */

#ifdef GMX_NBNXN_SIMD_4XN

static inline void
gmx_simd_ref_load_4xn_exclusions(unsigned int excl,
                                 gmx_simd_ref_exclmask mask_S0,
                                 gmx_simd_ref_exclmask mask_S1,
                                 gmx_simd_ref_exclmask mask_S2,
                                 gmx_simd_ref_exclmask mask_S3,
                                 gmx_simd_ref_pb *interact_S0,
                                 gmx_simd_ref_pb *interact_S1,
                                 gmx_simd_ref_pb *interact_S2,
                                 gmx_simd_ref_pb *interact_S3)
{
    /* Load integer interaction mask */
    gmx_simd_ref_exclmask mask_S = gmx_simd_ref_load1_exclmask(excl);
    *interact_S0  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S0);
    *interact_S1  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S1);
    *interact_S2  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S2);
    *interact_S3  = gmx_simd_ref_checkbitmask_pb(mask_S, mask_S3);
}

/* Test "none_IEEEEDDDD" function */

const int numMasks_IEEEEDDDD = 4;

template<> void
SimdFunctionTest<ReferenceFunction_none_IEEEEDDDD, TestFunction_none_IEEEEDDDD>::call(ReferenceFunction_none_IEEEEDDDD referenceFunction, const UIntArray theUInts)
{
    gmx_simd_ref_exclmask masks[numMasks_IEEEEDDDD];
    gmx_simd_ref_pb bResults[numMasks_IEEEEDDDD];
    gmx_simd_ref_pr rResults[numMasks_IEEEEDDDD];

    for(int i = 0; i < numMasks_IEEEEDDDD; ++i)
    {
        masks[i] = gmx_simd_ref_load_exclmask(theUInts[i]);
    }
    referenceFunction(theUInts[numMasks_IEEEEDDDD][0],
                      masks[0], masks[1],
                      masks[2], masks[3],
                      &bResults[0], &bResults[1],
                      &bResults[2], &bResults[3]);
    for(int i = 0; i < numMasks_IEEEEDDDD; ++i)
    {
        rResults[i] = gmx_simd_ref_convert_pb_to_pr(bResults[i]);
        gmx_simd_ref_store_pr(referenceResult[i], rResults[i]);
    }
}

template<> void
SimdFunctionTest<ReferenceFunction_none_IEEEEDDDD, TestFunction_none_IEEEEDDDD>::call(TestFunction_none_IEEEEDDDD testFunction, const UIntArray theUInts)
{
    gmx_exclmask masks[numMasks_IEEEEDDDD];
    gmx_mm_pb results[numMasks_IEEEEDDDD];

    for(int i = 0; i < numMasks_IEEEEDDDD; ++i)
    {
        masks[i] = gmx_load_exclmask(theUInts[i]);
    }
    testFunction(theUInts[numMasks_IEEEEDDDD][0],
                 masks[0], masks[1],
                 masks[2], masks[3],
                 &results[0], &results[1],
                 &results[2], &results[3]);
    for(int i = 0; i < numMasks_IEEEEDDDD; ++i)
    {
        /* On x86 and QPX, gmx_mm_bb is the same type as gmx_mm_pr, so
           this works */
        gmx_store_pr(testResult[i], results[i]);
    }
}

/* Helper typedef to keep the output pretty */
typedef SimdFunctionTest<ReferenceFunction_none_IEEEEDDDD, TestFunction_none_IEEEEDDDD> SimdFunctionWithSignature_none_IEEEEDDDD;

TEST_F(SimdFunctionWithSignature_none_IEEEEDDDD, gmx_load_4xn_exclusions_Works)
{
    Tester(gmx_simd_ref_load_4xn_exclusions,
           gmx_load_4xn_exclusions,
           bits,
           1.0,
           numMasks_IEEEEDDDD);
}

} // namespace

#endif /* GMX_NBNXN_SIMD_4XN */

#endif
