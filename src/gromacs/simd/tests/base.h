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
 * Classes to support tests of SIMD functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#ifndef _simd_tests_base_h_
#define _simd_tests_base_h_

#include "typedefs.h"
#include "smalloc.h"
#include "gromacs/legacyheaders/gmx_random.h"
#include "testutils/testoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/basicoptions.h"
#include "utils.h"
#include <maths.h>

namespace SIMDTests
{

/*! We use some namespaces to be sure the test and reference function
   names cannot collide. This is probably overkill.

   Note that the declaration of the test functions should precede that
   of the reference functions, so that types.h #defines
   GMX_SIMD_WIDTH_HERE, and the reference code can know to use that
   width. */
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

/********************************************************************
 * Tests for SIMD wrapper functions
 */

/*! \brief Base class to encapsulate all the SIMD support functions
 * needed for testing SIMD wrapper functions. Not intended to be
 * instantiated.
 *
 * Contains stubs for functions that re-wrap the SIMD wrapper
 * functions. These are not necessary, but they describe the
 * interface to which the derived classes should adhere if they want
 * to be called by the callFunction() functions in classes derived
 * from SimdFunctionTestBase. These functions here could be made pure
 * virtual, but there's no code that actually uses dynamic dispatch.
 */
template<class _realType, class _epiType, class _boolType, class _exclmaskType>
class SimdFunctionSetBase
{
    public:
        typedef _realType realType;
        typedef _epiType epiType;
        typedef _boolType boolType;
        typedef _exclmaskType exclmaskType;
        
        realType load_pr(const real *a);
        exclmaskType load_exclmask(const unsigned int *a);
        realType convert_pb_to_pr(const boolType a);
        boolType cmplt_pr(realType a, realType b);
        void store_pr(real *a, realType b); 
        epiType cvttpr_epi32(realType a);
};

namespace TestFunctions
{

/*! \brief Class to encapsulate all the support functions needed
 * for testing SIMD wrapper functions on x86. */

class x86SimdFunctionSet :
            public SimdFunctionSetBase<gmx_mm_pr, gmx_epi32, gmx_mm_pb, gmx_exclmask>
{
    public:
        realType load_pr(const real *a)
        {
            return TestFunctions::gmx_load_pr(a);
        }

        exclmaskType load_exclmask(const unsigned int *a)
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

        realType convert_pb_to_pr(const boolType a)
        {
            return a;
        }

        boolType cmplt_pr(realType a, realType b)
        {
            return TestFunctions::gmx_cmplt_pr(a, b);
        };

        void store_pr(real *a, realType b)
        {
            return TestFunctions::gmx_store_pr(a, b);
        }

        epiType
        cvttpr_epi32(realType a)
        {
            return gmx_cvttpr_epi32(a);
        }
};

}

namespace ReferenceFunctions
{

/*! \brief Class to encapsulate all the support functions needed
 * for testing the reference implementation of the SIMD wrapper
 * functions. */

class ReferenceSimdFunctionSet :
            public SimdFunctionSetBase<ReferenceFunctions::gmx_simd_ref_pr,
                                       ReferenceFunctions::gmx_simd_ref_epi32,
                                       ReferenceFunctions::gmx_simd_ref_pb,
                                       ReferenceFunctions::gmx_simd_ref_exclmask>
{
    public:
        realType load_pr(const real *a)
        {
            return ReferenceFunctions::gmx_simd_ref_load_pr(a);
        }

        exclmaskType load_exclmask(const unsigned int *a)
        {
            return ReferenceFunctions::gmx_simd_ref_load_exclmask(a);
        }

        realType convert_pb_to_pr(const boolType a)
        {
            return ReferenceFunctions::gmx_simd_ref_convert_pb_to_pr(a);
        }

        boolType cmplt_pr(realType a, realType b)
        {
            return ReferenceFunctions::gmx_simd_ref_cmplt_pr(a, b);
        };

        void store_pr(real *a, realType b)
        {
            return ReferenceFunctions::gmx_simd_ref_store_pr(a, b);
        }

        epiType
        cvttpr_epi32(realType a)
        {
            return ReferenceFunctions::gmx_simd_ref_cvttpr_epi32(a);
        }
};

}

/*! \brief Base class for test fixtures for SIMD functions
 *
 * \tparam NOutputs The number of output vectors produced by the SIMD function
 *
 * \todo If/when we have containers of SIMD-aligned reals, use them
 * here and in derived classes.
 */
template <size_t NOutputs>
class SimdFunctionTestBase :
            public ::testing::Test
{
    public:
        // TODO fix this later
        //#ifdef GMX_X86
        typedef TestFunctions::x86SimdFunctionSet HardwareSimdFunctionSet;
        //#endif

        SimdFunctionTestBase() : Test(), numberOfRepeats(1000)
        {
            /* This can be set on the command line with
               --gtest_random_seed=SEED. Zero is the GoogleTest
               default that means "seed based on the current time." */
            rng = gmx_rng_init(testing::UnitTest::GetInstance()->random_seed());

            snew_aligned(referenceResult,
                         GMX_SIMD_WIDTH_HERE * NOutputs,
                         GMX_SIMD_WIDTH_HERE * sizeof(real));
            snew_aligned(testResult,
                         GMX_SIMD_WIDTH_HERE * NOutputs,
                         GMX_SIMD_WIDTH_HERE * sizeof(real));

            /* Look for a value from the -numrepeats command line option */
            gmx::Options options(NULL, NULL);
            options.addOption(gmx::IntegerOption("numrepeats").store(&numberOfRepeats));
            gmx::test::parseTestOptions(&options);
            options.finish();
        }

        ~SimdFunctionTestBase()
        {
            gmx_rng_destroy(rng);
        }

        /*! A static constant seems to be found more regularly than a
            template parameter by the the compiler when compiling
            specializations of callFunctions() in derived classes. */
        static const size_t numOutputs = NOutputs;

        ReferenceFunctions::ReferenceSimdFunctionSet referenceSimd;
        HardwareSimdFunctionSet hardwareSimd;

        //! Random number generator
        gmx_rng_t rng;

        //! Controls how many tests on random inputs are run
        int numberOfRepeats;

        //! Vectors to store the outputs for later comparison
        real     *testResult, *referenceResult;
};

} // namespace

#endif /* _simd_tests_base_h_ */
