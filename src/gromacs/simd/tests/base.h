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

#ifndef _simd_tests_base_h_
#define _simd_tests_base_h_

#include "typedefs.h"
#include "smalloc.h"
#include "gromacs/legacyheaders/gmx_random.h"
#include "utils.h"
#include <maths.h>

namespace SIMDTests
{

/* Note that the test functions should precede reference functions, so
   that types.h #defines GMX_SIMD_WIDTH_HERE, and the reference code
   can know to use that width. */
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

/* Helper function to make it possible to compare the results of some
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
        
        virtual ~SimdFunctionSet() {};

        virtual realType load_pr(const real *a) = 0;
        virtual exclmaskType load_exclmask(const unsigned int *a) = 0;
        virtual realType convert_pb_to_pr(const boolType a) = 0;
        virtual boolType cmplt_pr(realType a, realType b) = 0;
        virtual void store_pr(real *a, realType b) = 0; 
        virtual epiType cvttpr_epi32(realType a) = 0;
};

namespace TestFunctions
{

class x86SimdFunctionSet :
            public SimdFunctionSet<gmx_mm_pr, gmx_epi32, gmx_mm_pb, gmx_exclmask>
{
    public:
        virtual ~x86SimdFunctionSet() {};

        virtual realType load_pr(const real *a)
        {
            return TestFunctions::gmx_load_pr(a);
        }

        virtual exclmaskType load_exclmask(const unsigned int *a)
        {
            /* In double precision, gmx_load_exclmask is intended to
               read a vector that is twice as long and has adjacent
               duplicate values (see nbnxn_atomdata_init()). So we do
               that. */
            exclmaskType return_value;
#ifdef GMX_DOUBLE
            unsigned int *duplicated_ints;
            snew_aligned(duplicated_ints, GMX_SIMD_WIDTH_HERE*2, sizeof(real) * GMX_SIMD_WIDTH_HERE);
            for (unsigned int i = 0; i != GMX_SIMD_WIDTH_HERE; ++i)
            {
                duplicated_ints[2*i+0] = a[i];
                duplicated_ints[2*i+1] = a[i];
            }
            return_value = TestFunctions::gmx_load_exclmask(duplicated_ints);
            sfree(duplicated_ints);
#else
            return_value = TestFunctions::gmx_load_exclmask(a);
#endif
            return return_value;
        }

        virtual realType convert_pb_to_pr(const boolType a)
        {
            return a;
        }

        virtual boolType cmplt_pr(realType a, realType b)
        {
            return TestFunctions::gmx_cmplt_pr(a, b);
        };

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

}

namespace ReferenceFunctions
{

class ReferenceSimdFunctionSet :
            public SimdFunctionSet<ReferenceFunctions::gmx_simd_ref_pr, ReferenceFunctions::gmx_simd_ref_epi32, ReferenceFunctions::gmx_simd_ref_pb, ReferenceFunctions::gmx_simd_ref_exclmask>
{
    public:
        virtual ~ReferenceSimdFunctionSet() {};

        virtual realType load_pr(const real *a)
        {
            return ReferenceFunctions::gmx_simd_ref_load_pr(a);
        }

        virtual exclmaskType load_exclmask(const unsigned int *a)
        {
            return ReferenceFunctions::gmx_simd_ref_load_exclmask(a);
        }

        virtual realType convert_pb_to_pr(const boolType a)
        {
            return ReferenceFunctions::gmx_simd_ref_convert_pb_to_pr(a);
        }

        virtual boolType cmplt_pr(realType a, realType b)
        {
            return ReferenceFunctions::gmx_simd_ref_cmplt_pr(a, b);
        };

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

}

//! TODO this makes more sense as a nbnxn-specific test, rather than a
//! general SIMD test

template <size_t NOutputs>
class SimdFunctionTestBase :
            public ::testing::Test
{
    public:
        // TODO fix this later
        //#ifdef GMX_X86
        typedef TestFunctions::x86SimdFunctionSet HardwareSimdFunctionSet;
        //#endif
        

        SimdFunctionTestBase() : Test()
        {
            rng = gmx_rng_init(SIMD_TEST_RANDOM_SEED);

            snew_aligned(referenceResult,
                         GMX_SIMD_WIDTH_HERE * NOutputs,
                         GMX_SIMD_WIDTH_HERE * sizeof(real));
            snew_aligned(testResult,
                         GMX_SIMD_WIDTH_HERE * NOutputs,
                         GMX_SIMD_WIDTH_HERE * sizeof(real));
        }

        ~SimdFunctionTestBase()
        {
            gmx_rng_destroy(rng);
        }

        static const size_t numOutputs = NOutputs;

        ReferenceFunctions::ReferenceSimdFunctionSet referenceSimd;
        HardwareSimdFunctionSet hardwareSimd;

        gmx_rng_t rng;

        real     *testResult;
        real     *referenceResult;
};

} // namespace

#endif
/* _simd_tests_base_h_ */
