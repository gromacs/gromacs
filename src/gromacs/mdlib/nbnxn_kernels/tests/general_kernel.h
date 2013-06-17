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
 * Tests for 4xn kernel functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */
#ifndef _simd_test_general_kernel_h_
#define _simd_test_general_kernel_h_

#include "gromacs/simd/tests/general.h"

#if defined GMX_NBNXN_SIMD_4XN || defined GMX_NBNXN_SIMD_2XNN

#if defined GMX_NBNXN_SIMD_4XN && defined GMX_NBNXN_SIMD_2XNN
//#error "Error: can't test both 4xn and 2xn kernel functionality in the same compilation unit"
#endif

namespace SIMDTests
{

#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_consts.h"

#ifdef GMX_NBNXN_SIMD_4XN
#define GMX_SIMD_J_UNROLL_SIZE 1
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn_common.h"
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
#define GMX_SIMD_J_UNROLL_SIZE 1
#include "gromacs/mdlib/nbnxn_kernels/simd_2xnn/nbnxn_kernel_simd_2xnn_common.h"
#endif
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils.h"

namespace ReferenceFunctions
{
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_ref.h"
}

/* TODO It'd be nice to be able to extend the existing SIMD wrapper
   classes, but I don't know how to do it. Pushing this content there
   does not work because it is too awkward to use all the above
 #includes. */

template<class _realType, class _epiType,
         class _boolType, class _exclfilterType,
         class _real4Type>
class SimdReal4FunctionSetBase :
    public SimdFunctionSetBase<_realType, _epiType, _boolType>
{
    public:
        //! Convenience typedef
        typedef _realType realType;
        //! Convenience typedef
        typedef _boolType boolType;
        //! Convenience typedef
        typedef _real4Type real4Type;
        //! Convenience typedef
        typedef _exclfilterType exclfilterType;

        //! Wrapper function to convert a 4-wide SIMD real to SIMD real
        realType convert_pr4_to_pr(const boolType a);
        //! Wrapper function to convert a SIMD bool to SIMD real
        realType convert_pb_to_pr(const boolType a);
        //! Wrapper function to load a SIMD exclusion filter
        exclfilterType load_exclusion_filter(const unsigned *a);
        //! Wrapper function to load a SIMD interaction mask
        _boolType load_interaction_mask_pb(long a, real *b);
};

namespace TestFunctions
{

#ifdef GMX_X86_SSE2

#if defined GMX_X86_AVX_256 && !defined GMX_USE_HALF_WIDTH_SIMD_HERE && !defined GMX_DOUBLE
//typedef gmx_mm_pr4 real4Type;
#else
typedef gmx_mm_pr gmx_mm_pr4;
#endif

/*! \brief Class to encapsulate all the support functions needed
 * for testing 4-wide SIMD wrapper functions on x86. */
class x86SimdReal4FunctionSet :
    public SimdReal4FunctionSetBase<gmx_mm_pr, gmx_epi32,
                                    gmx_mm_pb, gmx_exclfilter,
                                    gmx_mm_pr4>
{
    public:
        realType convert_pr4_to_pr(real4Type a)
        {
#if defined GMX_X86_AVX_256 && !defined GMX_USE_HALF_WIDTH_SIMD_HERE && !defined GMX_DOUBLE
            return _mm256_castps128_ps256(a);
#else
            return a;
#endif
        }

        realType convert_pb_to_pr(const boolType a)
        {
            return a;
        }

        exclfilterType load_exclusion_filter(const unsigned *a)
        {
            /* In double precision, gmx_load_exclusion_filter is
               intended to read a vector that is twice as long and has
               adjacent duplicate values (see
               nbnxn_atomdata_init()). So we do that. */
            exclfilterType  return_value;
#ifdef GMX_DOUBLE
            unsigned       *duplicated_ints;
            snew_aligned(duplicated_ints, GMX_SIMD_WIDTH_HERE*2, sizeof(real) * GMX_SIMD_WIDTH_HERE);
            for (unsigned int i = 0; i != GMX_SIMD_WIDTH_HERE; ++i)
            {
                duplicated_ints[2*i+0] = a[i];
                duplicated_ints[2*i+1] = a[i];
            }
            return_value = gmx_load_exclusion_filter(duplicated_ints);
            sfree(duplicated_ints);
#else
            return_value = gmx_load_exclusion_filter(a);
#endif
            return return_value;
        }
};

//! Convenience typedef
typedef x86SimdReal4FunctionSet HardwareSimdReal4FunctionSet;

#endif  /* GMX_X86_SSE2 */

#ifdef GMX_CPU_ACCELERATION_IBM_QPX

/*! \brief Class to encapsulate all the support functions needed
 * for testing SIMD wrapper functions on QPX (i.e. A2 core of BlueGene/Q). */
class QpxSimdReal4FunctionSet :
    public SimdReal4FunctionSetBase<gmx_mm_pr, gmx_epi32,
                                    gmx_mm_pb, gmx_exclfilter,
                                    gmx_mm_pr4>
{
    public:
        realType convert_pr4_to_pr(real4Type a)
        {
            return a;
        }

        realType convert_pb_to_pr(const boolType a)
        {
            return a;
        }

        exclfilterType load_exclusion_filter(const unsigned *a)
        {
            /* In double precision, gmx_load_exclusion_filter is
               intended to read a vector that is twice as long and has
               adjacent duplicate values (see
               nbnxn_atomdata_init()). So we do that. */
            exclfilterType  return_value;
#ifdef GMX_DOUBLE
            unsigned       *duplicated_ints;
            snew_aligned(duplicated_ints, GMX_SIMD_WIDTH_HERE*2, sizeof(real) * GMX_SIMD_WIDTH_HERE);
            for (unsigned int i = 0; i != GMX_SIMD_WIDTH_HERE; ++i)
            {
                duplicated_ints[2*i+0] = a[i];
                duplicated_ints[2*i+1] = a[i];
            }
            return_value = gmx_load_exclusion_filter(duplicated_ints);
            sfree(duplicated_ints);
#else
            return_value = gmx_load_exclusion_filter(a);
#endif
            return return_value;
        }

        boolType load_interaction_mask_pb(long a, real *b)
        {
            /* This function is only defined for QPX */
            return TestFunctions::gmx_load_interaction_mask_pb(a, b);
        }
};
typedef TestFunctions::QpxSimdReal4FunctionSet HardwareSimdReal4FunctionSet;

#endif  /* GMX_CPU_ACCELERATION_IBM_QPX */

}

namespace ReferenceFunctions
{

/*! Helper function to make it possible to compare the results of some
    reference code (which produces a boolean value) and the SIMD code
    (which produces something hardware-specific). */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_convert_pb_to_pr(gmx_simd_ref_pb src)
{
    gmx_simd_ref_pr result;
    int             i;
    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        result.r[i] = TestFunctions::HardwareSimdFunctionSet::getSimdBool(src.r[i]);
    }

    return result;
}

/*! Helper function for loading interaction masks */
static gmx_inline gmx_simd_ref_pb
gmx_simd_ref_load_interaction_mask_pb(long a, real *simd_interaction_array)
{
    gmx_simd_ref_pb b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        /* simd_interaction_array contains QPX SIMD boolean values(ie.
           -1 and 1). These need to be mapped to the 0 and 1 used
           elsewhere in the SIMD reference implementations. */
        real interaction = simd_interaction_array[a * GMX_SIMD_REF_WIDTH + i];
        b.r[i] = TestFunctions::QpxSimdFunctionSet::isSimdTrue(interaction) ? 1 : 0;
    }

    return b;
}

/*! \brief Class to encapsulate all the support functions needed
 * for testing the reference implementation of the SIMD wrapper
 * functions. */
class ReferenceSimdReal4FunctionSet :
    public SimdReal4FunctionSetBase<ReferenceFunctions::gmx_simd_ref_pr,
                                    ReferenceFunctions::gmx_simd_ref_epi32,
                                    ReferenceFunctions::gmx_simd_ref_pb,
                                    ReferenceFunctions::gmx_simd_ref_exclfilter,
                                    ReferenceFunctions::gmx_simd_ref_pr4>
{
    public:
        realType convert_pr4_to_pr(real4Type a)
        {
            realType result;
            int      i;

            for (i = 0; i < std::min(4, GMX_SIMD_REF_WIDTH); i++)
            {
                result.r[i] = a.r[i];
            }
            return result;
        }

        realType convert_pb_to_pr(const boolType a)
        {
            return ReferenceFunctions::gmx_simd_ref_convert_pb_to_pr(a);
        }

        exclfilterType load_exclusion_filter(const unsigned *a)
        {
            return ReferenceFunctions::gmx_simd_ref_load_exclusion_filter(a);
        }

        boolType load_interaction_mask_pb(long a, real *b)
        {
            return ReferenceFunctions::gmx_simd_ref_load_interaction_mask_pb(a, b);
        }
};

}   // namespace ReferenceFunctions

template <typename InputType>
class SimdFunctionTest_kernel
    : public SimdFunctionTest<InputType>
{
    public:
        //! Convenience typedef
        typedef SimdFunctionTest<InputType> Parent;

        /*! \brief Do repetitive testing of the SIMD and reference versions of the function.
         *
         * \tparam ReferenceFunctionType The type signature of the reference SIMD function
         * \tparam TestFunctionType The type signature of the reference SIMD function
         * \tparam InputKind The kind of data present in the input vectors, which is distinct from the type of data.
         *
         * Only the compiler needs to worry about the actual type that
         * is ReferenceFunctionType and TestFunctionType. Tests can
         * just pass in the reference and test versions of the
         * function and forget about it.
         *
         * Writes detailed output in failing cases. */
        template<typename ReferenceFunctionType,
                 typename TestFunctionType,
                 typename InputKind>
        void RunTest(ReferenceFunctionType referenceFunction,
                     TestFunctionType      testFunction,
                     int                   outputSimdWidth = GMX_SIMD_WIDTH_HERE,
                     InputKind             inputKind = 0);

        /*! \brief Do repetitive testing of the SIMD and reference versions of the function.
         *
         * \tparam ReferenceFunctionType The type signature of the reference SIMD function
         * \tparam TestFunctionType The type signature of the reference SIMD function
         * \tparam InputKind The kind of data present in the input vectors, which is distinc
           t from the type of data.
         *
         * Only the compiler needs to worry about the actual type that
         * is ReferenceFunctionType and TestFunctionType. Tests can
         * just pass in the reference and test versions of the
         * function and forget about it.
         *
         * Writes detailed output in failing cases. */
        template<typename ReferenceFunctionType,
                 typename TestFunctionType,
                 typename ReferenceSimdFunctionSet,
                 typename ReferenceSimdReal4FunctionSet,
                 typename TestSimdFunctionSet,
                 typename TestSimdReal4FunctionSet,
                 typename InputKind>
        void RunTest(ReferenceFunctionType         referenceFunction,
                     TestFunctionType              testFunction,
                     int                           filter_stride,
                     int                           unrollj,
                     ReferenceSimdFunctionSet      referenceSimdFunctionSet,
                     ReferenceSimdReal4FunctionSet referenceSimdReal4FunctionSet,
                     TestSimdFunctionSet           testSimdFunctionSet,
                     TestSimdReal4FunctionSet      testSimdReal4FunctionSet,
                     InputKind                     inputKind = 0);

        //! Method object for reference SIMD functions
        ReferenceFunctions::ReferenceSimdReal4FunctionSet referenceSimdReal4_;
        //! Method object for hardware SIMD functions
        TestFunctions::HardwareSimdReal4FunctionSet       hardwareSimdReal4_;
};

}      // namespace

#endif /* defined GMX_NBNXN_SIMD_4XN || defined GMX_NBNXN_SIMD_2XNN*/

#endif /* _simd_test_general_kernel_h_ */
