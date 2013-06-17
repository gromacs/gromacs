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

#ifndef GMX_SIMD_TESTS_BASE_H
#define GMX_SIMD_TESTS_BASE_H

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/gmx_random.h"
#include "testutils/testoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/basicoptions.h"
#include "utils.h"
#include <maths.h>

namespace SIMDTests
{

/*! We use some namespaces to be sure the test and reference versions
   of function names and constants cannot collide.

   Note that the declaration of the test functions should precede that
   of the reference functions, so that macros.h defines
   GMX_SIMD_WIDTH_HERE, and the reference code can know to use that
   width. */
#include "gromacs/simd/macros.h"

namespace ReferenceFunctions
{
#define GMX_SIMD_REFERENCE_PLAIN_C
#include "gromacs/simd/macros_ref.h"
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
 * functions. Declaring these is not necessary (they are not actually
 * inherited), but they describe the interface to which the derived
 * classes should adhere if they want to be called by the
 * callFunction() functions in classes derived from
 * SimdFunctionTestBase. These functions here could be made pure
 * virtual, but there's no code that actually uses dynamic dispatch.
 */
template<class _realType, class _epiType, class _boolType>
class SimdFunctionSetBase
{
    public:
        //! Convenience typedef
        typedef _realType realType;
        //! Convenience typedef
        typedef _epiType epiType;
        //! Convenience typedef
        typedef _boolType boolType;

        //! Wrapper function to load a SIMD real
        realType load_pr(const real *a);
        //! Wrapper function to compare SIMD reals
        boolType cmplt_pr(realType a, realType b);
        //! Wrapper function to store a SIMD real
        void store_pr(real *a, realType b);
        //! Wrapper function to store a SIMD bool
        void store_pb(real *a, boolType b);
        //! Wrapper function to convert a SIMD real to a SIMD integer
        epiType cvttpr_epi32(realType a);
};

namespace TestFunctions
{

/*! \brief Class to encapsulate all the support functions needed
 * for testing SIMD wrapper functions on x86. */
class x86SimdFunctionSet :
    public SimdFunctionSetBase<gmx_mm_pr, gmx_epi32, gmx_mm_pb>
{
    public:
        realType load_pr(const real *a)
        {
            return gmx_load_pr(a);
        }

        boolType cmplt_pr(realType a, realType b)
        {
            return gmx_cmplt_pr(a, b);
        };

        void store_pr(real *a, realType b)
        {
            return gmx_store_pr(a, b);
        }

        void store_pb(real *a, boolType b)
        {
#ifdef __MIC__
            /* On MIC, boolType is a short int, and the low 16 bits
               encode the logical boolean of the 16 vector
               elements. In this function, these are remapped so they
               can be compared with the reference-SIMD bool type in
               the same way as other x86 SIMD types. Also the
               remapping makes clear which bits are wrong when the
               test fails. */
            for (int i = 0; i < GMX_SIMD_WIDTH_HERE; i++)
            {
                a[i] = TestFunctions::HardwareSimdFunctionSet::getSimdBool(b & (1 << i));
            }
#else
            store_pr(a, b);
#endif
        }

        epiType
        cvttpr_epi32(realType a)
        {
            return gmx_cvttpr_epi32(a);
        }

        //! Translate a logical bool to the value that represents it
        //! in the hardware SIMD.
        static real getSimdBool(bool a)
        {
            BitManipulater manipulater;
            manipulater.i = a ? simdTrue : simdFalse;
            return manipulater.r;
        }

        //! Translate a hardware SIMD bool to the logical value that
        //! it represents.
        static bool isSimdTrue(real a)
        {
            /* A direct comparison of the real is impossible, because
               x86 SIMD true is a NaN, and that must compare false. */
            BitManipulater manipulater;
            manipulater.r = a;
            return simdTrue == manipulater.i;
        }

        //! Translate a hardware SIMD bool to the logical value that
        //! it represents.
        static bool isSimdFalse(real a)
        {
            /* A direct comparison of the real is impossible, because
               x86 SIMD true is a NaN, and that must compare false. */
            BitManipulater manipulater;
            manipulater.r = a;
            return simdFalse == manipulater.i;
        }

        //! Constant containing a value used by the SIMD hardware to
        //! represent logical "true."
        static const UnsignedIntWithSizeOfReal simdTrue;
        //! Constant containing a value used by the SIMD hardware to
        //! represent logical "false."
        static const UnsignedIntWithSizeOfReal simdFalse;
};

#ifdef GMX_X86_SSE2

//! Convenience typedef
typedef x86SimdFunctionSet HardwareSimdFunctionSet;

#endif  /* GMX_X86_SSE2 */

/*! \brief Class to encapsulate all the support functions needed
 * for testing SIMD wrapper functions on QPX (i.e. A2 core of BlueGene/Q). */
class QpxSimdFunctionSet :
    public SimdFunctionSetBase<gmx_mm_pr, gmx_epi32, gmx_mm_pb>
{
    public:
        realType load_pr(const real *a)
        {
            return gmx_load_pr(a);
        }

        boolType cmplt_pr(realType a, realType b)
        {
            return gmx_cmplt_pr(a, b);
        };

        void store_pr(real *a, realType b)
        {
            return gmx_store_pr(a, b);
        }

        void store_pb(real *a, boolType b)
        {
            return store_pr(a, b);
        }

        epiType
        cvttpr_epi32(realType a)
        {
            return gmx_cvttpr_epi32(a);
        }

        //! Translate a logical bool to the value that represents it
        //! in the hardware SIMD.
        static real getSimdBool(bool a)
        {
            return (real) (a ? simdTrue : simdFalse);
        }

        //! Translate a hardware SIMD bool to the logical value that
        //! it represents.
        static bool isSimdTrue(real a)
        {
            return simdTrue == a;
        }

        //! Translate a hardware SIMD bool to the logical value that
        //! it represents.
        static bool isSimdFalse(real a)
        {
            return simdFalse == a;
        }

        //! Convenience function to map a value from one form of SIMD
        //! boolean to that of the current hardware
        /*
           boolType mapSimdBooleanToHardwareSimdBoolean(boolType const &a)
           {
            return a;
           };
         */

        //! Constant containing a value used by the SIMD hardware to
        //! represent logical "true."
        static const real simdTrue;
        //! Constant containing a value used by the SIMD hardware to
        //! represent logical "false."
        static const real simdFalse;
};

#ifdef GMX_CPU_ACCELERATION_IBM_QPX

typedef TestFunctions::QpxSimdFunctionSet HardwareSimdFunctionSet;

#endif  /* GMX_CPU_ACCELERATION_IBM_QPX */

}

namespace ReferenceFunctions
{

/*! \brief Class to encapsulate all the support functions needed
 * for testing the reference implementation of the SIMD wrapper
 * functions. */
class ReferenceSimdFunctionSet :
    public SimdFunctionSetBase<ReferenceFunctions::gmx_simd_ref_pr,
                               ReferenceFunctions::gmx_simd_ref_epi32,
                               ReferenceFunctions::gmx_simd_ref_pb>
{
    public:
        realType load_pr(const real *a)
        {
            return ReferenceFunctions::gmx_simd_ref_load_pr(a);
        }

        boolType cmplt_pr(realType a, realType b)
        {
            return ReferenceFunctions::gmx_simd_ref_cmplt_pr(a, b);
        };

        void store_pr(real *a, realType b)
        {
            return ReferenceFunctions::gmx_simd_ref_store_pr(a, b);
        }

        void store_pb(real *a, boolType b)
        {
            realType r;
            for (int i = 0; i < GMX_SIMD_WIDTH_HERE; i++)
            {
                r.r[i] = TestFunctions::HardwareSimdFunctionSet::getSimdBool(b.r[i]);
            }
            return store_pr(a, r);
        }

        epiType
        cvttpr_epi32(realType a)
        {
            return ReferenceFunctions::gmx_simd_ref_cvttpr_epi32(a);
        }
};

}

/*! \brief Class for functor for templated aligned allocation for SIMD
 *  tests
 *
 * \todo If/when we have containers of SIMD-aligned reals, use them
 * here.
 */
template <class InputType>
struct AlignedSimdAllocater
{
    public:
        /*! \brief Constructor
         *
         * \param size Number of array elements to allocate. */
        AlignedSimdAllocater(size_t size) : size_(size) {};

        //! Operator to wrap the aligned allocation
        void operator() (InputType * &pointer)
        {
            snew_aligned(pointer,
                         GMX_SIMD_WIDTH_HERE * size_,
                         GMX_SIMD_WIDTH_HERE * sizeof(InputType));
        }

        //! Number of array elements to allocate.
        size_t size_;
};

/*! \brief Base class for test fixtures for SIMD functions
 *
 * \todo If/when we have containers of SIMD-aligned reals, use them
 * here and in derived classes.
 */
class SimdFunctionTestBase :
    public ::testing::Test
{
    public:
        //! Convenience typedef
        typedef ::testing::Test Parent;
        //! Convenience typedef
        typedef AlignedSimdAllocater<real> OutputAllocater;

        SimdFunctionTestBase() :
            Parent(), numOutputs_(0)
        {
            /* This can be set on the command line with
               --gtest_random_seed=SEED. Zero is the GoogleTest
               default that means "seed based on the current time." */
            rng_ = gmx_rng_init(testing::UnitTest::GetInstance()->random_seed());
        }

        ~SimdFunctionTestBase()
        {
            gmx_rng_destroy(rng_);
            sfree(referenceResult_);
            sfree(testResult_);
        }

        //! Method object for reference SIMD functions
        ReferenceFunctions::ReferenceSimdFunctionSet referenceSimd_;
        //! Method object for hardware SIMD functions
        TestFunctions::HardwareSimdFunctionSet       hardwareSimd_;

        //! Number of output vectors from the SIMD function under test
        size_t    numOutputs_;
        //! Random number generator
        gmx_rng_t rng_;
        /*! Vector to store the concatenated output vectors from the
           test function for later comparison. */
        real     *testResult_;
        /*! Vector to store the concatenated output vectors from the
           reference function for later comparison. */
        real     *referenceResult_;
};

}      // namespace

#endif
