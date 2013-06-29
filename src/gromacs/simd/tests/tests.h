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

#include <string>
#include <ostream>
#include <gtest/gtest.h>
#include "typedefs.h"
#include "smalloc.h"
#include "gromacs/legacyheaders/gmx_random.h"
#include "utils.h"

namespace SIMDTests
{

/********************************************************************
 * Tests for SIMD wrapper functions/macros
 */

/* The real SIMD code must be #included before the reference code,
 * because the latter learns what SIMD width to use from the
 * former. */
#include "gromacs/simd/types.h"
#include "gromacs/simd/macros.h"
#define GMX_SIMD_REFERENCE_PLAIN_C
#include "gromacs/simd/reference_types.h"
#include "gromacs/simd/reference.h"
#undef GMX_SIMD_REFERENCE_PLAIN_C

// TODO Proper doxygen docs

//! TODO tests for hpr stuff?

/* Typedefs for testing both actual SIMD and reference C versions with
 * various type signatures.
 *
 * There are lots of function signatures, some of which are reused, so
 * we will use test classes that are templatized by the function
 * signatures. We need typedefs for those signatures so we don't drown
 * in a swamp of function pointers, and those typedefs need names. So
 * we use a shorthand that lists the return type, an underscore and
 * then each of the argument types. So "I_VPP" means a function
 * returning an int that takes (in order) a vector, a pointer to a
 * vector, and a pointer to a vector.
 *
 * Definitions used in the above scheme:
 *
 * I = int
 * V = a SIMD "vector of real" type
 * J = a SIMD "vector of signed int" type
 * E = a SIMD "vector of exclusion masks" type (an unsigned int on x86)
 * P = a pointer to a SIMD "vector of real" type
 * B = a SIMD "vector of bool" type
 * D = a pointer to a SIMD "vector of bool" type
 * R = a (non-vector) real
 * C = a pointer to a const real (possibly an array whose width is that of the SIMD width)
 * Q = a pointer to a real (possibly an array whose width is that of the SIMD width)
 * none = nothing (i.e. void)
 *
 * Currently, the B_BB and B_VV signatures are not actually
 * required. See test_logical_operations.cpp for details.
 *
 * gmx_store_pb is not tested. It is only used when
 * GMX_SIMD_HAVE_ANYTRUE is undefined, and in practice on x86,
 * GMX_SIMD_HAVE_ANYTRUE is always defined. On QPX, anytrue is not
 * available, but gmx_store_pb will be implemented as vec_st() to some
 * integer type and if vec_st() is broken then none of the testing
 * will work.
 *
 * TODO one of gmx_cvttpr_epi32, gmx_cvtepi32_pr is not yet completely
 * tested
 *
 * gmx_load1_exclmask, gmx_load_exclmask and gmx_checkbitmask_pb are
 * all involved in implementing exclusions. Those implementations are
 * highly architecture-specific, which is a pain to test. The
 * composite operation of loading exclusion masks has thus been
 * bundled into the wrapper functions gmx_load_*_exclusions(), which
 * are tested. */

/*
   typedef gmx_simd_ref_pb
   (*ReferenceFunction_B_BB)(gmx_simd_ref_pb,
                          gmx_simd_ref_pb);

   typedef gmx_mm_pb
   (*TestFunction_B_BB)(gmx_mm_pb,
                     gmx_mm_pb);

   typedef gmx_simd_ref_pb
   (*ReferenceFunction_B_VV)(gmx_simd_ref_pr,
                          gmx_simd_ref_pr);

   typedef gmx_mm_pb
   (*TestFunction_B_VV)(gmx_mm_pr,
                     gmx_mm_pr);
 */

/* Haven't worked out a way to template these typedef pairs. Oh
   well. Not worth much anyway. */
typedef gmx_simd_ref_pr
(*ReferenceFunction_V_VB)(gmx_simd_ref_pr,
                          gmx_simd_ref_pb);

typedef gmx_mm_pr
(*TestFunction_V_VB)(gmx_mm_pr,
                     gmx_mm_pb);

typedef gmx_simd_ref_pr
(*ReferenceFunction_V_none)();

typedef gmx_mm_pr
(*TestFunction_V_none)();

typedef gmx_simd_ref_pr
(*ReferenceFunction_V_R)(real);

typedef gmx_mm_pr
(*TestFunction_V_R)(real);

typedef gmx_simd_ref_pr
(*ReferenceFunction_V_C)(const real *);

typedef gmx_mm_pr
(*TestFunction_V_C)(const real *);

typedef gmx_simd_ref_pr
(*ReferenceFunction_V_V)(gmx_simd_ref_pr);

typedef gmx_mm_pr
(*TestFunction_V_V)(gmx_mm_pr);

typedef gmx_simd_ref_pr
(*ReferenceFunction_V_VV)(gmx_simd_ref_pr,
                          gmx_simd_ref_pr);

typedef gmx_mm_pr
(*TestFunction_V_VV)(gmx_mm_pr,
                     gmx_mm_pr);

typedef gmx_simd_ref_pr
(*ReferenceFunction_V_VVV)(gmx_simd_ref_pr,
                           gmx_simd_ref_pr,
                           gmx_simd_ref_pr);

typedef gmx_mm_pr
(*TestFunction_V_VVV)(gmx_mm_pr,
                      gmx_mm_pr,
                      gmx_mm_pr);

typedef int
(*ReferenceFunction_I_VPP)(gmx_simd_ref_pr,
                           gmx_simd_ref_pr*,
                           gmx_simd_ref_pr*);

typedef int
(*TestFunction_I_VPP)(gmx_mm_pr,
                      gmx_mm_pr*,
                      gmx_mm_pr*);

/* This is a hack - the only functions with signature I_B are only
   implemented on x86, where V and B are equivalent, and it is more
   convenient to use V in the implementation of their tests. */
typedef int
(*ReferenceFunction_I_V)(gmx_simd_ref_pr);

typedef int
(*TestFunction_I_V)(gmx_mm_pr);

typedef void
(*ReferenceFunction_none_QV)(real *,
                             gmx_simd_ref_pr);

typedef void
(*TestFunction_none_QV)(real *,
                        gmx_mm_pr);

typedef gmx_simd_ref_pr
(*ReferenceFunction_V_BVV)(gmx_simd_ref_pb,
                           gmx_simd_ref_pr,
                           gmx_simd_ref_pr);

typedef gmx_mm_pr
(*TestFunction_V_BVV)(gmx_mm_pb,
                      gmx_mm_pr,
                      gmx_mm_pr);

typedef void
(*ReferenceFunction_none_IEEDD)(unsigned int, gmx_simd_ref_exclmask, gmx_simd_ref_exclmask, gmx_simd_ref_pb *, gmx_simd_ref_pb *);

typedef void
(*TestFunction_none_IEEDD)(unsigned int, gmx_exclmask, gmx_exclmask, gmx_mm_pb *, gmx_mm_pb *);

typedef void
(*ReferenceFunction_none_IEEEEDDDD)(unsigned int, gmx_simd_ref_exclmask, gmx_simd_ref_exclmask, gmx_simd_ref_exclmask, gmx_simd_ref_exclmask, gmx_simd_ref_pb *, gmx_simd_ref_pb *, gmx_simd_ref_pb *, gmx_simd_ref_pb *);

typedef void
(*TestFunction_none_IEEEEDDDD)(unsigned int, gmx_exclmask, gmx_exclmask, gmx_exclmask, gmx_exclmask, gmx_mm_pb *, gmx_mm_pb *, gmx_mm_pb *, gmx_mm_pb *);

/*! \brief Some helper constants. */
/* TODO make these all configurable from the command line */
// TODO make these static constants in a parent class of SimdFunctionTest?
const int SIMD_TEST_NUM_REAL_VARIABLES = 1024;
// In practice, at most 5 variables are used in any test, but using 8
// is a bit future-proof. This produces 1024/8 = 128 independent sets
// of test variables. That seems like enough to be confident we'll hit
// edge cases decently often.
const int SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST = 8;
const int SIMD_TEST_RANDOM_SEED                     = 9273;

/*! \brief This class that specifies the kind of SIMD function being
 * tested, i.e. what type signatures the reference and test versions
 * have.
 *
 * It also wraps the calls to the reference and test versions of that
 * function, and provides a function to be called from test bodies.
 *
 * It also prepares data structures for different kinds of tests.
 *
 * TODO Mention why we have booleanReals
 */
template <class ReferenceFunctionType, class TestFunctionType>
class SimdFunctionTest : public ::testing::Test
{
    public:
        typedef ReferenceFunctionType refType;
        typedef TestFunctionType testType;

        // TODO reimplement all the data handling using aligned-memory
        // containers if/when we get such things.
        typedef real *RealArray[SIMD_TEST_NUM_REAL_VARIABLES];
        typedef unsigned int *UIntArray[SIMD_TEST_NUM_REAL_VARIABLES];

        SimdFunctionTest()
        {
            gmx_rng_t rng;

            rng = gmx_rng_init(SIMD_TEST_RANDOM_SEED);

            for (int i = 0; i < SIMD_TEST_NUM_REAL_VARIABLES; ++i)
            {
                snew_aligned(uints[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(unsigned int));
                snew_aligned(bits[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(unsigned int));
                snew_aligned(reals[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(real));
                snew_aligned(positiveReals[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(real));
                snew_aligned(booleanReals[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(real));
                for (int j = 0; j < GMX_SIMD_WIDTH_HERE; ++j)
                {
                    /* Generate random unsigned ints. */
                    uints[i][j] = gmx_rng_uniform_uint32(rng);

                    /* Generate random unsigned bits in a 32-bit
                       integer, i.e. 1^{0,...31}. */
                    bits[i][j] = 1 << (gmx_rng_uniform_uint32(rng) & 31);

                    /* Generate random floating-point numbers. It is
                     * not clear what range makes any particular
                     * sense, so chose -0.5 <= x < 0.5 to exercise the
                     * distance scale of most interest in MD. */
                    reals[i][j] = gmx_rng_uniform_real(rng) - 0.5;

                    /* Likewise, 0 <= x < 1 */
                    positiveReals[i][j] = gmx_rng_uniform_real(rng);

                    /* Generate random booleans according to the SIMD
                       hardware's view of true and false. */
                    BitManipulater manipulater;
                    /* When QPX is supported, or GMX_X86 exists, the
                       following lines will need to be versioned. */
                    UnsignedIntWithSizeOfReal simdTrue = ~0;
                    UnsignedIntWithSizeOfReal simdFalse = 0;
                    manipulater.i      = (gmx_rng_uniform_real(rng) <
                                          gmx_rng_uniform_real(rng)) ?
                        simdTrue : simdFalse;
                    booleanReals[i][j] = manipulater.r;
                }
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

        /* We could template the test and reference versions of these
           functions, but you would have to do so on the function
           type, the SIMD variable type, and the SIMD load/store
           functions. This isn't worth it for four-line
           functions. Normal overloading will do.

           Likewise, we could template on real and unsigned int, but
           there's not much to gain there either. */
        virtual void call(refType referenceFunction, RealArray referenceReals);
        virtual void call(testType testFunction, RealArray testReals);
        void call(refType referenceFunction, const UIntArray referenceUInts);
        void call(testType testFunction, const UIntArray testUInts);

        /*! (regarding scaleMaxUlps) Some SIMD-based computation (like
            rsqrt) is not binary exact (even after Newton-Raphson
            iterations), so a way tune the required tolerance between
            reference and test case results is needed.

            When a test fails, this code shows the observed
            differences in results, because that is useful for seeing
            whether a change in scaleMaxUlps would solve the issue
            (for the tests where that is sensible). Ideally, it would
            show the inputs and outputs for each failing case of each
            function, but I don't think setting up and maintaining
            such a mechanism is worth it. */
        void Tester(refType referenceFunction, testType testFunction, RealArray reals, real scaleMaxUlps = 1.0, int numOutputs = 1)
        {
            for (int i = 0; i < SIMD_TEST_NUM_REAL_VARIABLES; i += SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST)
            {
                call(referenceFunction, reals + i);
                call(testFunction, reals + i);

                for (int i = 0; i < numOutputs; i++)
                {
                    EXPECT_TRUE(RealArraysAreEqual(referenceResult[i],
                                                   testResult[i],
                                                   GMX_SIMD_WIDTH_HERE,
                                                   scaleMaxUlps));
                }
            }
        }

        void Tester(refType referenceFunction, testType testFunction, UIntArray uints, real scaleMaxUlps = 1.0, int numOutputs = 1)
        {
            for (int i = 0; i < SIMD_TEST_NUM_REAL_VARIABLES; i += SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST)
            {
                call(referenceFunction, uints + i);
                call(testFunction, uints + i);

                for (int i = 0; i < numOutputs; i++)
                {
                    //! TODO consider using gmxcheck gear instead
                    EXPECT_TRUE(RealArraysAreEqual(referenceResult[i],
                                                   testResult[i],
                                                   GMX_SIMD_WIDTH_HERE,
                                                   scaleMaxUlps));
                }
            }
        }

        // TODO: recast this with vector<real> (or whatever from Eigen, etc.)
        UIntArray uints;
        UIntArray bits;
        RealArray reals;
        RealArray positiveReals;
        /* Vector containing the SIMD version of truth of whether one
           real was less than another. Different SIMD architectures
           define "true" differently, e.g. on x86 true is 0xffffffff,
           whereas on QPX true is greater than or equal to either kind
           of zero. */
        RealArray booleanReals;
        real     *testResult[SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST];
        real     *referenceResult[SIMD_TEST_MAX_NUM_REAL_VARIABLES_PER_TEST];
};

} // namespace
