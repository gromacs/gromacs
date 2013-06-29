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

#include <vector>

namespace SIMDTests
{

namespace TestFunctions
{
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils.h"
}
namespace ReferenceFunctions
{
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_ref.h"
}

/********************************************************************
 * Tests for SIMD wrapper functions
 */

template<class _realType, class _epiType, class _boolType, class _exclmaskType>
class SimdFunctionSet
{
    public:
        typedef _realType realType;
        typedef _epiType epiType;
        typedef _boolType boolType;
        typedef _exclmaskType exclmaskType;

        virtual realType load_pr(const real *a) = 0;
        virtual void store_pr(real *a, realType b) = 0; 
        virtual epiType cvttpr_epi32(realType a) = 0;
        
        virtual ~SimdFunctionSet() {};
};

class x86SimdFunctionSet :
            public SimdFunctionSet<gmx_mm_pr, gmx_epi32, gmx_mm_pb, gmx_exclmask>
{
    public:
        virtual realType load_pr(const real *a)
        {
            return gmx_load_pr(a);
        }

        virtual void store_pr(real *a, realType b)
        {
            return gmx_store_pr(a, b);
        }

        virtual epiType
        cvttpr_epi32(realType a)
        {
            return gmx_cvttpr_epi32(a);
        }
};

class ReferenceSimdFunctionSet :
            public SimdFunctionSet<gmx_simd_ref_pr, gmx_simd_ref_epi32, gmx_simd_ref_pb, gmx_simd_ref_exclmask>
{
    public:
        virtual realType load_pr(const real *a)
        {
            return gmx_simd_ref_load_pr(a);
        }

        virtual void store_pr(real *a, realType b)
        {
            return gmx_simd_ref_store_pr(a, b);
        }

        virtual epiType
        cvttpr_epi32(realType a)
        {
            return gmx_simd_ref_cvttpr_epi32(a);
        }
};

//! TODO this makes more sense as a nbnxn-specific test, rather than a
//! general SIMD test

typedef void
(*ReferenceFunction_load_table_f)(const real *, gmx_simd_ref_epi32, int *,
                                  gmx_simd_ref_pr *, gmx_simd_ref_pr *);

typedef void
(*TestFunction_load_table_f)(const real *, gmx_epi32, int *,
                             gmx_mm_pr *, gmx_mm_pr *);

/* Helper typedef to keep the output pretty */

class SimdFunction_load_table_f :
            public ::testing::Test
{
    public:
        typedef ReferenceFunction_load_table_f refType;
        typedef TestFunction_load_table_f testType;

        // TODO reimplement all the data handling using aligned-memory
        // containers if/when we get such things.
        typedef real *RealArray[SIMD_TEST_NUM_REAL_VARIABLES];

        SimdFunction_load_table_f()
        {
            // Picked an arbitrary table size comparable with the size
            // of tables that are actually used.
            const unsigned int table_size = 1024 * GMX_SIMD_WIDTH_HERE;
            const unsigned int index_size = GMX_SIMD_WIDTH_HERE;

            gmx_rng_t rng;

            rng = gmx_rng_init(SIMD_TEST_RANDOM_SEED);

            // make table with random values
            snew_aligned(table, table_size, GMX_SIMD_WIDTH_HERE * sizeof(real));
    
            snew_aligned(index, index_size, GMX_SIMD_WIDTH_HERE * sizeof(real));

            for(unsigned int i = 0; i != table_size; ++i)
            {
                /* Generate random floating-point numbers. It is not clear
                 * what range makes any particular sense, so pick 0.5 >= x >
                 * -0.5. */
                table[i] = gmx_rng_uniform_real(rng) - 0.5;
            }
            for (unsigned int i = 0; i != index_size; ++i)
            {
                // make random non-overflowing indices into the table
                index[i] = gmx_rng_uniform_real(rng) * (table_size / GMX_SIMD_WIDTH_HERE);
            }

            for (int i = 0; i < SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST; ++i)
            {
                snew_aligned(testResult[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(real));
                snew_aligned(referenceResult[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(real));
            }

            gmx_rng_destroy(rng);
        }

        void call(refType referenceFunction);
        void call(testType testFunction);
        
        template<class SimdFunctionSet,
                 typename FunctionType>
        void callFunction(SimdFunctionSet &support,
                          FunctionType function,
                          real **_result);

        void Tester(refType referenceFunction, testType testFunction, real scaleMaxUlps = 1.0, int numOutputs = 1)
        {
            for (int i = 0; i < SIMD_TEST_NUM_REAL_VARIABLES; i += SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST)
            {
                //                call(referenceFunction);
                //                call(testFunction);

                x86SimdFunctionSet x86Support;
                ReferenceSimdFunctionSet referenceSupport;
                callFunction(referenceSupport, referenceFunction, referenceResult);
                callFunction(x86Support, testFunction, testResult);

                for (int i = 0; i < numOutputs; i++)
                {
                    EXPECT_TRUE(RealArraysAreEqual(referenceResult[i],
                                                   testResult[i],
                                                   GMX_SIMD_WIDTH_HERE,
                                                   scaleMaxUlps));
                }
            }
        }

        // Table of floats to look up
        real *table;

        // SIMD vector of indices to look up in the table. In the kernels,
        // this is computed from reals, so we do the same here.
        real *index;

        gmx_rng_t rng;

        real     *testResult[SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST];
        real     *referenceResult[SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST];
};

template<class SimdFunctionSet,
         typename FunctionType>
void
SimdFunction_load_table_f::callFunction(SimdFunctionSet &support,
                                        FunctionType function,
                                        real **_result)
{
    typename SimdFunctionSet::epiType index_epi;
    typename SimdFunctionSet::realType result[2];

    // Build the temp array some implementations of the function
    // require to be passed to it.
    int ti_array[GMX_SIMD_WIDTH_HERE*2], *ti;
    ti = gmx_simd_align_int(ti_array);

    typename SimdFunctionSet::realType indices = support.load_pr(index);
    index_epi = support.cvttpr_epi32(indices);

    function(table, index_epi, ti,
             &result[0], &result[1]);
    for (unsigned int i = 0; i != 2; ++i)
    {
        support.store_pr(_result[i], result[i]);
    }
}

// void
// SimdFunction_load_table_f::call(refType referenceFunction)
// {
//     gmx_simd_ref_epi32 index_epi32;
//     gmx_simd_ref_pr result[2];

//     // Build the temp array some implementations of the function
//     // require to be passed to it.
//     int ti_array[GMX_SIMD_WIDTH_HERE*2], *ti;
//     ti = gmx_simd_align_int(ti_array);

//     gmx_simd_ref_pr indices = gmx_simd_ref_load_pr(index);
//     index_epi32 = gmx_simd_ref_cvttpr_epi32(indices);

//     referenceFunction(table, index_epi32, ti,
//                       &result[0], &result[1]);
//     for (unsigned int i = 0; i != 2; ++i)
//     {
//         gmx_simd_ref_store_pr(referenceResult[i], result[i]);
//     }
// }

// void
// SimdFunction_load_table_f::call(testType testFunction)
// {
//     gmx_epi32 index_epi32;
//     gmx_mm_pr result[2];

//     // Build the temp array some implementations of the function
//     // require to be passed to it.
//     int ti_array[GMX_SIMD_WIDTH_HERE*2], *ti;
//     ti = gmx_simd_align_int(ti_array);

//     gmx_mm_pr indices = gmx_load_pr(index);
//     index_epi32 = gmx_cvttpr_epi32(indices);

//     testFunction(table, index_epi32, ti,
//                       &result[0], &result[1]);
//     for (unsigned int i = 0; i != 2; ++i)
//     {    
//         gmx_store_pr(testResult[i], result[i]);
//     }
// }

TEST_F(SimdFunction_load_table_f, load_table_f_Works)
{
    Tester(ReferenceFunctions::load_table_f,
           TestFunctions::load_table_f,
           1.0,
           2);
}

} // namespace
