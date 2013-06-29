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

//#include "tests.h"

#include <vector>
#include "typedefs.h"
#include "smalloc.h"
#include "gromacs/legacyheaders/gmx_random.h"
#include "utils.h"

namespace SIMDTests
{

namespace TestFunctions
{
#include "gromacs/simd/types.h"
#include "gromacs/simd/macros.h"
}
namespace ReferenceFunctions
{
#define GMX_SIMD_REFERENCE_PLAIN_C
#include "gromacs/simd/reference_types.h"
#include "gromacs/simd/reference.h"
#undef GMX_SIMD_REFERENCE_PLAIN_C
}

/*! \brief Some helper constants. */
/* TODO make these all configurable from the command line */
// TODO make these static constants in a parent class of SimdFunctionTest?
const int SIMD_TEST_NUM_REAL_VARIABLES = 1024;
// In practice, at most 5 variables are used in any test, but using 8
// is a bit future-proof. This produces 1024/8 = 128 independent sets
// of test variables. That seems like enough to be confident we'll hit
// edge cases decently often.
//const int SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST = 8;
const int SIMD_TEST_RANDOM_SEED                     = 9273;

/********************************************************************
 * Tests for SIMD wrapper functions
 */

/*! \brief Abstract base class to encapsulate all the SIMD support
 * functions needed for testing SIMD wrapper functions.
 *
 * Contains virtual functions that re-wrap the SIMD wrapper functions,
 * so that dynamic dispatch can be used to call an equivalent code
 * path in the test and the reference code, without needing to
 * duplicate code. The production code uses static inline versions of
 * these functions so that there is no overhead, but in test code we
 * prefer convenience to speed.
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
            return TestFunctions::gmx_load_pr(a);
        }

        virtual void store_pr(real *a, realType b)
        {
            return TestFunctions::gmx_store_pr(a, b);
        }

        virtual epiType
        cvttpr_epi32(realType a)
        {
            return gmx_cvttpr_epi32(a);
        }
};

class ReferenceSimdFunctionSet :
            public SimdFunctionSet<ReferenceFunctions::gmx_simd_ref_pr, ReferenceFunctions::gmx_simd_ref_epi32, ReferenceFunctions::gmx_simd_ref_pb, ReferenceFunctions::gmx_simd_ref_exclmask>
{
    public:
        virtual realType load_pr(const real *a)
        {
            return ReferenceFunctions::gmx_simd_ref_load_pr(a);
        }

        virtual void store_pr(real *a, realType b)
        {
            return ReferenceFunctions::gmx_simd_ref_store_pr(a, b);
        }

        virtual epiType
        cvttpr_epi32(realType a)
        {
            return ReferenceFunctions::gmx_simd_ref_cvttpr_epi32(a);
        }
};

//! TODO this makes more sense as a nbnxn-specific test, rather than a
//! general SIMD test

template <size_t N>
class SimdFunctionTest_ :
            public ::testing::Test
{
    public:
        // TODO fix this later
        //#ifdef GMX_X86
        typedef x86SimdFunctionSet HardwareSimdFunctionSet;
        //#endif

        SimdFunctionTest_() : Test()
        {
            rng = gmx_rng_init(SIMD_TEST_RANDOM_SEED);

            snew_aligned(referenceResult,
                         GMX_SIMD_WIDTH_HERE * numOutputs,
                         GMX_SIMD_WIDTH_HERE * sizeof(real));
            snew_aligned(testResult,
                         GMX_SIMD_WIDTH_HERE * numOutputs,
                         GMX_SIMD_WIDTH_HERE * sizeof(real));
        }

        ~SimdFunctionTest_()
        {
            gmx_rng_destroy(rng);
        }

        const unsigned int numOutputs = N;

        ReferenceSimdFunctionSet referenceSimd;
        HardwareSimdFunctionSet hardwareSimd;

        gmx_rng_t rng;

        real     *testResult;
        real     *referenceResult;
};

#if (defined GMX_NBNXN_SIMD_2XNN) || (defined GMX_NBNXN_SIMD_4XN)

#ifdef GMX_NBNXN_SIMD_4XN
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_4xn_outer_header.h"
#endif

#ifdef GMX_NBNXN_SIMD_2XNN
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_2xnn_outer_header.h"
#endif

namespace TestFunctions
{
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils.h"
}

namespace ReferenceFunctions
{
#define GMX_SIMD_REFERENCE_PLAIN_C
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_ref.h"
#undef GMX_SIMD_REFERENCE_PLAIN_C
}

template <size_t N>
class SimdFunction_load_table :
            public SimdFunctionTest_<N>
{
    public:

        SimdFunction_load_table() : SimdFunctionTest_<N>()
        {
            snew_aligned(table, table_size, GMX_SIMD_WIDTH_HERE * sizeof(real));
            for(unsigned int i = 0; i != table_size; ++i)
            {
                /* Generate random floating-point numbers. It is not clear
                 * what range makes any particular sense, so pick 0.5 >= x >
                 * -0.5. */
                table[i] = gmx_rng_uniform_real(rng) - 0.5;
            }

            snew_aligned(index, GMX_SIMD_WIDTH_HERE, GMX_SIMD_WIDTH_HERE * sizeof(real));

            ti = TestFunctions::gmx_simd_align_int(ti_array);
        }

        ~SimdFunction_load_table()
        {
            sfree(table);
            sfree(index);
        }

        //! This function gets defined differently in different
        // specializations of SimdFunction_load_table. Those
        // specializations have different numbers of output values.
        // This trickery supports testing both load_table_f and
        // load_table_f_v with nearly the same code.
        template<class SimdFunctionSet,
                 typename FunctionType>
        void callFunction(SimdFunctionSet &simdFunctionSet,
                          FunctionType function,
                          const real *_table,
                          const real *_index,
                          real *_result);

        template<typename ReferenceFunctionType,
                 typename TestFunctionType>
        void Tester(ReferenceFunctionType referenceFunction,
                    TestFunctionType testFunction)
        {
            for (int i = 0; i < SIMD_TEST_NUM_REAL_VARIABLES; ++i)
            {
                // make random non-overflowing indices into the table
                for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                {
                    index[j] = gmx_rng_uniform_real(rng) * (table_size / GMX_SIMD_WIDTH_HERE);
                }

                callFunction(referenceSimd, referenceFunction, table, index, referenceResult);
                callFunction(hardwareSimd, testFunction, table, index, testResult);

                EXPECT_TRUE(RealArraysAreEqual(referenceResult,
                                               testResult,
                                               GMX_SIMD_WIDTH_HERE * numOutputs,
                                               1));
            }
        }

        using SimdFunctionTest_<N>::numOutputs;
        using SimdFunctionTest_<N>::referenceSimd;
        using SimdFunctionTest_<N>::hardwareSimd;
        using SimdFunctionTest_<N>::rng;
        using SimdFunctionTest_<N>::referenceResult;
        using SimdFunctionTest_<N>::testResult;

        // Picked an arbitrary table size comparable with the size of
        // tables that are actually used.
        const unsigned int table_size = 1024 * GMX_SIMD_WIDTH_HERE;

        // Table of random floats to look up in tests
        real *table;

        // SIMD vector of indices to look up in the table. In the kernels,
        // this is computed from reals, so we do the same here.
        real *index;

        // Temp array required by some implementations of load_table*
        // Allocate twice the necessary size to be sure alignment can
        // be achieved.
        int ti_array[GMX_SIMD_WIDTH_HERE+GMX_SIMD_WIDTH_HERE], *ti;
};

/* Specialization of the class for testing load_table_f */
template<>
template<class SimdFunctionSet,
         typename FunctionType>
void SimdFunction_load_table<2>::callFunction(SimdFunctionSet &simdFunctionSet,
                                              FunctionType function,
                                              const real *_table,
                                              const real *_index,
                                              real *_result)
{
    typename SimdFunctionSet::epiType index_epi;
    typename SimdFunctionSet::realType result[2];
    typename SimdFunctionSet::realType indices;

    indices = simdFunctionSet.load_pr(_index);
    index_epi = simdFunctionSet.cvttpr_epi32(indices);
    
    function(_table, index_epi, ti, &result[0], &result[1]);
    
    for (unsigned int i = 0; i != numOutputs; ++i)
    {
        simdFunctionSet.store_pr(_result + i * GMX_SIMD_WIDTH_HERE, result[i]);
    }
}

/* Specialization of the class for testing load_table_f_v */
template<>
template<class SimdFunctionSet,
         typename FunctionType>
void SimdFunction_load_table<3>::callFunction(SimdFunctionSet &simdFunctionSet,
                                              FunctionType function,
                                              const real *_table,
                                              const real *_index,
                                              real *_result)
{
    typename SimdFunctionSet::epiType index_epi;
    typename SimdFunctionSet::realType result[numOutputs];
    typename SimdFunctionSet::realType indices;

    indices = simdFunctionSet.load_pr(_index);
    index_epi = simdFunctionSet.cvttpr_epi32(indices);
    
#ifdef TAB_FDV0
    function(_table, index_epi, ti, &result[0], &result[1], &result[2]);
#else
    function(_table, _table, index_epi, ti, &result[0], &result[1], &result[2]);
#endif

    for (unsigned int i = 0; i != numOutputs; ++i)
    {
        simdFunctionSet.store_pr(_result + i * GMX_SIMD_WIDTH_HERE, result[i]);
    }
}

/*
typedef void
(*ReferenceFunctionType)(const real *, gmx_simd_ref_epi32, int *,
                         gmx_simd_ref_pr *, gmx_simd_ref_pr *);

typedef void
(*TestFunctionType)(const real *, gmx_epi32, int *,
                    gmx_mm_pr *, gmx_mm_pr *);
*/

typedef SimdFunction_load_table<2> SimdFunction_load_table_f;
TEST_F(SimdFunction_load_table_f, load_table_f_Works)
{
    
    Tester(ReferenceFunctions::load_table_f,
           TestFunctions::load_table_f);
}

typedef SimdFunction_load_table<3> SimdFunction_load_table_f_v;
TEST_F(SimdFunction_load_table_f_v, load_table_f_v_Works)
{
    
    Tester(ReferenceFunctions::load_table_f_v,
           TestFunctions::load_table_f_v);
}

#endif

// ====

typedef SimdFunctionTest_<1> SimdFunctionTestWithOneOutput;

template <size_t M>
class SimdFunctionWithRealInputs :
            public SimdFunctionTestWithOneOutput
{
    public:
        typedef SimdFunctionTest_<1> Parent;

        SimdFunctionWithRealInputs() : Parent(), shift(-0.5), scale(0),
                                       maxUlps(1)
        {
            for (unsigned int i = 0; i < numRealInputs; ++i)
            {
                snew_aligned(reals[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(real));
            }
        }

        ~SimdFunctionWithRealInputs()
        {
            for (unsigned int i = 0; i < numRealInputs; ++i)
            {
                sfree(reals[i]);
            }
        }

        template<class SimdFunctionSet,
                 typename FunctionType>
        void callFunction(SimdFunctionSet &simdFunctionSet,
                          FunctionType function,
                          real *_reals[M],
                          real *_result);

        template<typename ReferenceFunctionType,
                 typename TestFunctionType>
        void Tester(ReferenceFunctionType referenceFunction,
                    TestFunctionType testFunction)
        {
            for (int k = 0; k < SIMD_TEST_NUM_REAL_VARIABLES; ++k)
            {
                for (unsigned int i = 0; i < numRealInputs; ++i)
                {
                    // make random reals
                    for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                    {
                        /* Generate random floating-point numbers in
                         * scale * (shift + [0,1)) */
                        reals[i][j] = scale*(shift+gmx_rng_uniform_real(rng));
                    }
                }

                callFunction(referenceSimd, referenceFunction, reals, referenceResult);
                callFunction(hardwareSimd, testFunction, reals, testResult);

                EXPECT_TRUE(RealArraysAreEqual(referenceResult,
                                               testResult,
                                               GMX_SIMD_WIDTH_HERE * numOutputs,
                                               maxUlps));
            }
        }

        using Parent::numOutputs;
        using Parent::referenceSimd;
        using Parent::hardwareSimd;
        using Parent::rng;
        using Parent::referenceResult;
        using Parent::testResult;

        const unsigned int numRealInputs = M;

        //! TODO reimplement all the data handling using
        // aligned-memory containers if/when we get such things.
        real *reals[M];
        // Variables to control the range of the random reals used for inputs
        real shift, scale;
        // Maximum acceptable difference in "units in last place"
        // (ulps) for equality when comparing reals. Default is 1, and
        // produces standard Google Test behaviour. Larger values are
        // necessary for (e.g.) reduced-precision SIMD followed by
        // Newton-Raphson iterates, particularly in double precision.
        real maxUlps;
};

typedef SimdFunctionWithRealInputs<1> SimdFunctionTest_V_V;

/* Specialization of the class for testing functions with signature V_V */
template<>
template<class SimdFunctionSet,
         typename FunctionType>
void SimdFunctionTest_V_V::callFunction(SimdFunctionSet &simdFunctionSet,
                                        FunctionType function,
                                        real *_reals[1],
                                        real *_result)
{
    typename SimdFunctionSet::realType a, b, result;
    
    a      = simdFunctionSet.load_pr(_reals[0]);
    result = function(a);
    simdFunctionSet.store_pr(_result, result);
}

/* This function is used for truncation of reals before table lookups,
 * so moderate positive values are those most worth testing. And test
 * some negative numbers in case the implementation ever changes. */
TEST_F(SimdFunctionTest_V_V, gmx_round_pr_Works)
{
    shift = -0.25;
    scale = 1024;
    Tester(ReferenceFunctions::gmx_simd_ref_round_pr,
           TestFunctions::gmx_round_pr);
}

/* SIMD reciprocal square roots are computed by an approximate SIMD
 * operation, and then some Newton-Raphson iterations depending on the
 * accuracy of the underlying hardware operation and the prevailing
 * GROMACS single/double precision.
 *
 * 1) The SIMD rsqrt is never used without subsequent Newton-Raphson
 * iteration, and there is no reference version of approximate inverse
 * square root, so it only makes sense to test the composite
 * operation (i.e. not gmx_invsqrt_pr).
 *
 * 2) GROMACS does not do enough Newton-Raphson iterations for full
 * double precision (because we don't need it), so a wider range for
 * the equality test is required in that case. */
TEST_F(SimdFunctionTest_V_V, gmx_rsqrt_pr_Works)
{
    shift = 0;
    scale = 5;
#ifdef GMX_DOUBLE
    maxUlps = 60.0; // empirically determined to be enough on x86
#endif
    Tester(ReferenceFunctions::gmx_simd_ref_rsqrt_pr,
           TestFunctions::gmx_invsqrt_pr);
}

typedef SimdFunctionWithRealInputs<2> SimdFunctionTest_V_VV;

/* Specialization of the class for testing functions with signature V_VV */
template<>
template<class SimdFunctionSet,
         typename FunctionType>
void SimdFunctionTest_V_VV::callFunction(SimdFunctionSet &simdFunctionSet,
                                         FunctionType function,
                                         real *_reals[2],
                                         real *_result)
{
    typename SimdFunctionSet::realType a, b, result;
    
    a      = simdFunctionSet.load_pr(_reals[0]);
    b      = simdFunctionSet.load_pr(_reals[1]);
    result = function(a, b);
    simdFunctionSet.store_pr(_result, result);
}

TEST_F(SimdFunctionTest_V_VV, gmx_mul_pr_Works)
{
    Tester(ReferenceFunctions::gmx_simd_ref_mul_pr,
           TestFunctions::gmx_mul_pr);
}


} // namespace
