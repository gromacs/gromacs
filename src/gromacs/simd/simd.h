/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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

/* The macros in this file are intended to be used for writing
 * architecture-independent SIMD intrinsics code.
 * To support a new architecture, adding macros here should be (nearly)
 * all that is needed.
 */

#ifndef GMX_SIMD_SIMD_H
#define GMX_SIMD_SIMD_H

#include "gromacs/legacyheaders/types/simple.h"

/* architecture-independent SIMD intrinsics macros.
 *
 * The macros in this file are intended to be used for writing
 * architecture-independent SIMD intrinsics code. Rather than making assumptions
 * based on architecture, use the defines listed below to check what you can
 * do. This way it will be easy to extend GROMACS with support for new
 * simd instruction sets in the future, and speed up your old code too.
 *
 * Unfortunately there is no standard for SIMD architectures. The available
 * features vary a lot, but we still need to use quite a few of them to
 * get the best performance possible. This means some features will only
 * be available on certain platforms, and it is critical that we do NOT make
 * to many assumptions about the storage formats, their size or SIMD width.
 * Just to give a few examples:
 *
 * - On x86, double precision (64-bit) floating-point values always convert
 *   to 32-bit integers, while many other platforms use 64-bit, and some cannot
 *   use 32-bit integers at all. This means we cannot use a mask (boolean)
 *   derived from integer operations to select double-precision floating-point
 *   values, and it could get very complex for higher-level code if all these
 *   decisions were exposed. Instead, we want to keep integers 32-bit since
 *   all algorithms anyway need to work in single precision (w. 32-bit ints).
 * - IBM QPX uses 4-wide SIMD both for single and double precision. Integer
 *   support is highly limited, and the storage format means QPX does not
 *   use x86-style all-ones masks (which have different widths in single/double)
 *   but it uses the sign bit to denote the _false_ value. In particular, this
 *   means we cannot use the bit contents for any fancy mask operations.
 * - AVX1 only supports 4-wide 128-bit integer SIMD arithmetics, but the integer
 *   _conversions_ can still be done 8-wide which corresponds to the single
 *   precision floating-point width. Similarly, with AVX1 conversions between
 *   double-precision and integers use the 32-bit 4-wide 128bit registers where
 *   we can also do integer arithmetics. AVX2 adds proper arithmetics for
 *   8-wide integers. We would severely limit performance if we had to say
 *   that integer support was not present, so instead we stick to 32-bit ints
 *   but limit the operations we expose (and do shuffling internally).
 * - For SSE2 through SSE4.1, double precision is 2-wide, but when we convert
 *   to integers they will be put in the first two elements of a 4-wide integer
 *   type. This means we cannot assume that floating-point SIMD registers and
 *   corresponding integer registers (after conversion) have the same width.
 * - The 2-wide SIMD instructions on BlueGene/L and BlueGene/P cannot do any
 *   floating-point logical operations (and/andnot/or/xor) whatsoever.
 * - Since boolean values can have different width for float/double and the
 *   integers corresponding to float/double, we need to use separate boolean
 *   types for all these values and convert between them if we e.g. want to use
 *   result of an integer compare to select floating-point values.
 *
 * While this might sound complicated, it is actually far easier than writing
 * separate SIMD code for 10 architectures in both single & double. We
 * typcially implement SIMD support for a new architecture in less than a day,
 * and then it will be supported everywhere in Gromacs. For the higer-level
 * code, the only important thing is to never _assume_ anything about the SIMD
 * architecture. This module handles this by introducing a number of datatypes:
 *
 * FLOATING POINT:
 * - Default Gromacs "real" precision: gmx_simd_real_t and gmx_simd_add_r()
 * - Explicit single precision floating-point values: gmx_simd_float_t, and
 *   routines operating on them have names like gmx_simd_add_f()
 * - Explicit double precision floating-point values: gmx_simd_double_t, and
 *   routines operating on them have names like gmx_simd_add_d()
 *
 * INTEGERS CORRESPONDING TO FLOATING-POINT VALUES:
 * - Integers used when converting to/from Gromacs default "real" type:
 *   gmx_simd_int32_t, and gmx_simd_add_i().
 * - Integers used when converting to/from float:
 *   gmx_simd_fint32_t, and gmx_simd_add_fi(). This will also be the widest
 *   integer data type if you want to do pure integer SIMD operations.
 * - Integers used when converting to/from double (narrower type e.g. on AVX):
 *   gmx_simd_dint32_t, and gmx_simd_add_di().
 *
 * Note that all integer load/stores operations defined here load/store 32-bit
 * integers, even when the internal register storage might be 64-bit, and we
 * set the "width" of the SIMD implementation based on how many float/double/
 * integers we load/store - even if the internal width could be larger.
 *
 * BOOLEAN VALUES
 * We need a separate boolean datatype for masks and comparison results, since
 * we cannot assume the are identical either to integers, floats or double -
 * some implementations use specific predicate registers for booleans.
 * - Opeartions on "real": gmx_simd_bool_t, and gmx_simd_or_b()
 * - Operations specifically on float: gmx_simd_fbool_t, and gmx_simd_or_fb()
 * - Operations specifically on double: gmx_simd_dbool_t, and gmx_simd_or_db()
 * - For integers corresponding to "real": gmx_simd_ibool_t, gmx_simd_or_ib()
 * - For integers corresponding to float: gmx_simd_fibool_t, gmx_simd_or_fib()
 * - For integers corresponding to double: gmx_simd_dibool_t, gmx_simd_or_dib()
 *
 * If this seems daunting, in practice you should only need to use these types:
 * - Default Gromacs precision floating-point data: gmx_simd_real_t
 * - Corresponding integer data: gmx_simd_int32_t
 * - Booleans for the floating-point data: gmx_simd_bool_t
 * - Booleans for the integer data: gmx_simd_ibool_t
 *
 * This should be sufficient for code that works with the full SIMD width.
 * Unfortunately reality is not that simple. Some algorithms like lattice
 * summation need quartets of elements, so even when the SIMD width is >4 we
 * need width-4 SIMD if it is supported. These datatypes and operations use the
 * prefix gmx_simd4_, and availability is indicated by GMX_SIMD4_HAVE_FLOAT and
 * amd GMX_SIMD4_HAVE_DOUBLE. For now we only support a small subset of SIMD
 * operations for SIMD4, but that is trivial to extend if we need to.
 *
 * Functionality-wise, we have a small set of core set of features that we
 * require to be present on all platforms, while more avanced features can be
 * used in the code when defines like e.g. GMX_SIMD_HAVE_LOADU are set.
 *
 * This is a summary of the currently available preprocessor defines that
 * you should use to check for support when using the corresponding features.
 * We first list the float/double/int defines set by the _implementation_; in
 * most cases you do not want to check directly for float/double defines, but
 * you should instead use the derived "real" defines set in this file - we list
 * those at the end below.
 *
 * Defines set by the low-level implementation. These are only set if they work
 * for all datatypes; GMX_SIMD_HAVE_LOADU thus means we can load both float,
 * double, and integers from unaligned memory.
 *
 * GMX_SIMD_HAVE_FLOAT              Single-precision instructions available
 * GMX_SIMD_HAVE_DOUBLE             Double-precision instructions available
 * GMX_SIMD_HAVE_HARDWARE           Set when we are NOT emulating SIMD
 * GMX_SIMD_HAVE_LOADU              Load from unaligned memory available
 * GMX_SIMD_HAVE_STOREU             Store to unaligned memory available
 * GMX_SIMD_HAVE_LOGICAL            Support for and/andnot/or/xor on fp regs.
 * GMX_SIMD_HAVE_FMA                Floating-point fused multiply-add.
 *                                  Note: GROMACS provides emulated FMA
 *                                  instructions if you do not have FMA
 *                                  support, but in that case you might be
 *                                  able to code it more efficient w/o FMA.
 * GMX_SIMD_HAVE_FRACTION           Instruction to get decimal fraction.
 *                                  Same as FMA: This denotes hardware support,
 *                                  otherwise instruction will be emulated.
 * GMX_SIMD_HAVE_FINT32             Integer conversions to/from float.
 * GMX_SIMD_HAVE_FINT32_EXTRACT     Support for extracting integer SIMD elements.
 * GMX_SIMD_HAVE_FINT32_LOGICAL     Bitwise shift operations on gmx_simd_fint32_t
 * GMX_SIMD_HAVE_FINT32_ARITHMETICS Arithmetic ops for gmx_simd_fint32_t
 * GMX_SIMD_HAVE_DINT32             Integer conversions to/from double
 * GMX_SIMD_HAVE_DINT32_EXTRACT     Support for extracting integer SIMD elements.
 * GMX_SIMD_HAVE_DINT32_LOGICAL     Bitwise shift operations on gmx_simd_dint32_t
 * GMX_SIMD_HAVE_DINT32_ARITHMETICS Arithmetic ops for gmx_simd_dint32_t
 *
 * Corresponding SIMD4 macros:
 *
 * GMX_SIMD4_HAVE_FLOAT             Support for float width-4-SIMD. This
 *                                  also requires SIMD4 support for all
 *                                  architecture features defined for normal
 *                                  float SIMD (e.g. GMX_SIMD_HAVE_FMA).
 * GMX_SIMD4_HAVE_DOUBLE            Support for double width-4-SIMD. Same
 *                                  comment as for the float version above.

 * Information about the implementation (higher-level code can check these)
 *
 * GMX_SIMD_FLOAT_WIDTH             Number of elements in gmx_simd_float_t,
 *                                  and practical width of gmx_simd_fint32_t.
 * GMX_SIMD_DOUBLE_WIDTH            Number of elements in gmx_simd_double_t,
 *                                  and practical width of gmx_simd_dint32_t.
 * GMX_SIMD_RSQRT_BITS              Accuracy (bits) of 1/sqrt(x) lookup step.
 * GMX_SIMD_RCP_BITS                Accuracy (bits) of 1/x lookup step.
 *
 * After including the low-level architecture-specific implementation, this
 * header sets the following derived defines based on the current precision:
 *
 * GMX_SIMD_HAVE_REAL               GMX_SIMD_HAVE_{FLOAT|DOUBLE}
 * GMX_SIMD4_HAVE_REAL              GMX_SIMD4_HAVE_{FLOAT|DOUBLE}
 * GMX_SIMD_REAL_WIDTH              GMX_SIMD_{FLOAT|DOUBLE}_WIDTH
 * GMX_SIMD_HAVE_INT32              GMX_SIMD_HAVE_{FINT32|DINT32}
 * GMX_SIMD_INT32_WIDTH             GMX_SIMD_{FINT32|DINT32}_WIDTH
 * GMX_SIMD_HAVE_INT32_EXTRACT      GMX_SIMD_HAVE_{FINT32|DINT32}_EXTRACT
 * GMX_SIMD_HAVE_INT32_LOGICAL      GMX_SIMD_HAVE_{FINT32|DINT32}_LOGICAL
 * GMX_SIMD_HAVE_INT32_ARITHMETICS  GMX_SIMD_HAVE_{FINT32|DINT32}_ARITHMETICS
 *
 * GMX_SIMD4_WIDTH                  Always defined to 4. Just a convenience you
 *                                  can use to make it clear that you refer to
 *                                  the SIMD register width, not some other "4".
 *
 * While all these defines are available to specify the features of the
 * hardware, we would strongly recommend that you do NOT sprinkle your code
 * with defines - if nothing else it will be a debug nightmare. Instead you can
 * write a slower generic SIMD function that works everywhere, and then override
 * this with faster architecture-specific versions for some implementations. The
 * recommended way to do that is to e.g. add #ifndef checks for
 * GMX_SIMD_SINCOS_IMPLEMENTED around the gmx_simd_sincos_{d,f}() functions.
 */

/* Forward declarations so we can use memory allocation functions in impl. */
static gmx_inline float *  gmx_simd_align_f(float *p);
static gmx_inline double * gmx_simd_align_d(double *p);
static gmx_inline int *    gmx_simd_align_fi(int *p);
static gmx_inline int *    gmx_simd_align_di(int *p);
static gmx_inline float *  gmx_simd4_align_f(float *p);
static gmx_inline double * gmx_simd4_align_d(double *p);

/* Intel MIC is a bit special since it is a co-processor. This means the rest
 * of GROMACS (which runs on the CPU) should use a default SIMD set like AVX,
 * while the part running on the coprocessor defines __MIC__. All functions in
 * this SIMD module are static, so it will work perfectly fine to include this
 * file with different SIMD definitions for different files.
 */
#if defined __MIC__
#    include "gromacs/simd/impl_intel_mic.h"
#elif defined GMX_SIMD_X86_AVX2_256
#    include "gromacs/simd/impl_x86_avx2_256.h"
#elif defined GMX_SIMD_X86_AVX_256
#    include "gromacs/simd/impl_x86_avx_256.h"
#elif defined GMX_SIMD_X86_AVX_128_FMA
#    include "gromacs/simd/impl_x86_avx_128_fma.h"
#elif defined GMX_SIMD_X86_SSE4_1
#    include "gromacs/simd/impl_x86_sse4_1.h"
#elif defined GMX_SIMD_X86_SSE2
#    include "gromacs/simd/impl_x86_sse2.h"
#elif defined GMX_SIMD_IBM_QPX
#    include "gromacs/simd/impl_ibm_qpx.h"
#elif defined GMX_SIMD_SPARC64_HPC_ACE
#    include "gromacs/simd/impl_sparc64_hpc_ace.h"
#elif (defined GMX_SIMD_REFERENCE) || (defined DOXYGEN)
/* Plain C SIMD reference implementation, also serves as documentation */
#    include "gromacs/simd/impl_reference.h"
#endif

/* Just for convenience to make it clear you refer to SIMD4 width */
#define GMX_SIMD4_WIDTH    4

static gmx_inline float *
gmx_simd_align_f(float *p)
{
#    ifdef GMX_SIMD_HAVE_FLOAT
    return (float *)(((size_t)((p)+GMX_SIMD_FLOAT_WIDTH-1)) & (~((size_t)(GMX_SIMD_FLOAT_WIDTH*sizeof(float)-1))));
#    else
    return p;
#    endif
}

static gmx_inline double *
gmx_simd_align_d(double *p)
{
#    ifdef GMX_SIMD_HAVE_DOUBLE
    return (double *)(((size_t)((p)+GMX_SIMD_DOUBLE_WIDTH-1)) & (~((size_t)(GMX_SIMD_DOUBLE_WIDTH*sizeof(double)-1))));
#    else
    return p;
#    endif
}

static gmx_inline int *
gmx_simd_align_fi(int *p)
{
#    ifdef GMX_SIMD_HAVE_FINT32
    return (int *)(((size_t)((p)+GMX_SIMD_FINT32_WIDTH-1)) & (~((size_t)(GMX_SIMD_FINT32_WIDTH*sizeof(int)-1))));
#    else
    return p;
#    endif
}

static gmx_inline int *
gmx_simd_align_di(int *p)
{
#    ifdef GMX_SIMD_HAVE_DINT32
    return (int *)(((size_t)((p)+GMX_SIMD_DINT32_WIDTH-1)) & (~((size_t)(GMX_SIMD_DINT32_WIDTH*sizeof(int)-1))));
#    else
    return p;
#    endif
}

static gmx_inline float *
gmx_simd4_align_f(float *p)
{
#    ifdef GMX_SIMD4_HAVE_FLOAT
    return (float *)(((size_t)((p)+GMX_SIMD4_WIDTH-1)) & (~((size_t)(GMX_SIMD4_WIDTH*sizeof(float)-1))));
#    else
    return p;
#    endif
}

static gmx_inline double *
gmx_simd4_align_d(double *p)
{
#    ifdef GMX_SIMD4_HAVE_DOUBLE
    return (double *)(((size_t)((p)+GMX_SIMD4_WIDTH-1)) & (~((size_t)(GMX_SIMD4_WIDTH*sizeof(double)-1))));
#    else
    return p;
#    endif
}

/* Define Gromacs "real" precision macros depending on Gromacs config. Note
 * that conversions float-to-double and v.v. are not included here since they
 * are not precision-dependent - find them in the implementation files.
 */
#ifdef GMX_DOUBLE
/* Double floating-point */
#    define gmx_simd_real_t                  gmx_simd_double_t
#    define gmx_simd_load_r                  gmx_simd_load_d
#    define gmx_simd_load1_r                 gmx_simd_load1_d
#    define gmx_simd_set1_r                  gmx_simd_set1_d
#    define gmx_simd_store_r                 gmx_simd_store_d
#    define gmx_simd_loadu_r                 gmx_simd_loadu_d
#    define gmx_simd_storeu_r                gmx_simd_storeu_d
#    define gmx_simd_setzero_r               gmx_simd_setzero_d
#    define gmx_simd_add_r                   gmx_simd_add_d
#    define gmx_simd_sub_r                   gmx_simd_sub_d
#    define gmx_simd_mul_r                   gmx_simd_mul_d
#    define gmx_simd_fmadd_r                 gmx_simd_fmadd_d
#    define gmx_simd_fmsub_r                 gmx_simd_fmsub_d
#    define gmx_simd_fnmadd_r                gmx_simd_fnmadd_d
#    define gmx_simd_fnmsub_r                gmx_simd_fnmsub_d
#    define gmx_simd_and_r                   gmx_simd_and_d
#    define gmx_simd_andnot_r                gmx_simd_andnot_d
#    define gmx_simd_or_r                    gmx_simd_or_d
#    define gmx_simd_xor_r                   gmx_simd_xor_d
#    define gmx_simd_rsqrt_r                 gmx_simd_rsqrt_d
#    define gmx_simd_rcp_r                   gmx_simd_rcp_d
#    define gmx_simd_fabs_r                  gmx_simd_fabs_d
#    define gmx_simd_fneg_r                  gmx_simd_fneg_d
#    define gmx_simd_max_r                   gmx_simd_max_d
#    define gmx_simd_min_r                   gmx_simd_min_d
#    define gmx_simd_round_r                 gmx_simd_round_d
#    define gmx_simd_trunc_r                 gmx_simd_trunc_d
#    define gmx_simd_fraction_r              gmx_simd_fraction_d
#    define gmx_simd_get_exponent_r          gmx_simd_get_exponent_d
#    define gmx_simd_get_mantissa_r          gmx_simd_get_mantissa_d
#    define gmx_simd_set_exponent_r          gmx_simd_set_exponent_d
/* Double integer and conversions */
#    define gmx_simd_int32_t                 gmx_simd_dint32_t
#    define gmx_simd_load_i                  gmx_simd_load_di
#    define gmx_simd_set1_i                  gmx_simd_set1_di
#    define gmx_simd_store_i                 gmx_simd_store_di
#    define gmx_simd_loadu_i                 gmx_simd_loadu_di
#    define gmx_simd_storeu_i                gmx_simd_storeu_di
#    define gmx_simd_setzero_i               gmx_simd_setzero_di
#    define gmx_simd_cvt_r2i                 gmx_simd_cvt_d2i
#    define gmx_simd_cvtt_r2i                gmx_simd_cvtt_d2i
#    define gmx_simd_cvt_i2r                 gmx_simd_cvt_i2d
#    define gmx_simd_extract_i               gmx_simd_extract_di
#    define gmx_simd_slli_i                  gmx_simd_slli_di
#    define gmx_simd_srli_i                  gmx_simd_srli_di
#    define gmx_simd_and_i                   gmx_simd_and_di
#    define gmx_simd_andnot_i                gmx_simd_andnot_di
#    define gmx_simd_or_i                    gmx_simd_or_di
#    define gmx_simd_xor_i                   gmx_simd_xor_di
#    define gmx_simd_add_i                   gmx_simd_add_di
#    define gmx_simd_sub_i                   gmx_simd_sub_di
#    define gmx_simd_mul_i                   gmx_simd_mul_di
/* Double booleans and selection */
#    define gmx_simd_bool_t                  gmx_simd_dbool_t
#    define gmx_simd_cmpeq_r                 gmx_simd_cmpeq_d
#    define gmx_simd_cmplt_r                 gmx_simd_cmplt_d
#    define gmx_simd_cmple_r                 gmx_simd_cmple_d
#    define gmx_simd_and_b                   gmx_simd_and_db
#    define gmx_simd_or_b                    gmx_simd_or_db
#    define gmx_simd_anytrue_b               gmx_simd_anytrue_db
#    define gmx_simd_blendzero_r             gmx_simd_blendzero_d
#    define gmx_simd_blendv_r                gmx_simd_blendv_d
#    define gmx_simd_ibool_t                 gmx_simd_dibool_t
#    define gmx_simd_cmpeq_i                 gmx_simd_cmpeq_di
#    define gmx_simd_cmplt_i                 gmx_simd_cmplt_di
#    define gmx_simd_and_ib                  gmx_simd_and_dib
#    define gmx_simd_or_ib                   gmx_simd_or_dib
#    define gmx_simd_anytrue_ib              gmx_simd_anytrue_dib
#    define gmx_simd_blendzero_i             gmx_simd_blendzero_di
#    define gmx_simd_blendv_i                gmx_simd_blendv_di
/* Conversions between integer and double floating-point booleans */
#    define gmx_simd_cvt_b2ib                gmx_simd_cvt_db2dib
#    define gmx_simd_cvt_ib2b                gmx_simd_cvt_dib2db

/* SIMD4 double fp - we only support a subset of SIMD instructions for SIMD4 */
#    define gmx_simd4_real_t                 gmx_simd4_double_t
#    define gmx_simd4_load_r                 gmx_simd4_load_d
#    define gmx_simd4_load1_r                gmx_simd4_load1_d
#    define gmx_simd4_set1_r                 gmx_simd4_set1_d
#    define gmx_simd4_store_r                gmx_simd4_store_d
#    define gmx_simd4_loadu_r                gmx_simd4_loadu_d
#    define gmx_simd4_storeu_r               gmx_simd4_storeu_d
#    define gmx_simd4_setzero_r              gmx_simd4_setzero_d
#    define gmx_simd4_add_r                  gmx_simd4_add_d
#    define gmx_simd4_sub_r                  gmx_simd4_sub_d
#    define gmx_simd4_mul_r                  gmx_simd4_mul_d
#    define gmx_simd4_fmadd_r                gmx_simd4_fmadd_d
#    define gmx_simd4_fmsub_r                gmx_simd4_fmsub_d
#    define gmx_simd4_fnmadd_r               gmx_simd4_fnmadd_d
#    define gmx_simd4_fnmsub_r               gmx_simd4_fnmsub_d
#    define gmx_simd4_and_r                  gmx_simd4_and_d
#    define gmx_simd4_andnot_r               gmx_simd4_andnot_d
#    define gmx_simd4_or_r                   gmx_simd4_or_d
#    define gmx_simd4_xor_r                  gmx_simd4_xor_d
#    define gmx_simd4_rsqrt_r                gmx_simd4_rsqrt_d
#    define gmx_simd4_fabs_r                 gmx_simd4_fabs_d
#    define gmx_simd4_fneg_r                 gmx_simd4_fneg_d
#    define gmx_simd4_max_r                  gmx_simd4_max_d
#    define gmx_simd4_min_r                  gmx_simd4_min_d
#    define gmx_simd4_round_r                gmx_simd4_round_d
#    define gmx_simd4_trunc_r                gmx_simd4_trunc_d
#    define gmx_simd4_dotproduct3_r          gmx_simd4_dotproduct3_d
#    define gmx_simd4_bool_t                 gmx_simd4_dbool_t
#    define gmx_simd4_cmpeq_r                gmx_simd4_cmpeq_d
#    define gmx_simd4_cmplt_r                gmx_simd4_cmplt_d
#    define gmx_simd4_cmple_r                gmx_simd4_cmple_d
#    define gmx_simd4_and_b                  gmx_simd4_and_db
#    define gmx_simd4_or_b                   gmx_simd4_or_db
#    define gmx_simd4_anytrue_b              gmx_simd4_anytrue_db
#    define gmx_simd4_blendzero_r            gmx_simd4_blendzero_d
#    define gmx_simd4_blendv_r               gmx_simd4_blendv_d

/* Memory allocation */
#    define gmx_simd_align_r                 gmx_simd_align_d
#    define gmx_simd_align_i                 gmx_simd_align_di
#    define gmx_simd4_align_r                gmx_simd4_align_d

#    ifdef GMX_SIMD_HAVE_DOUBLE
#        define GMX_SIMD_HAVE_REAL
#        define GMX_SIMD_REAL_WIDTH          GMX_SIMD_DOUBLE_WIDTH
#    endif
#    ifdef GMX_SIMD_HAVE_DINT32
#        define GMX_SIMD_HAVE_INT32
#        define GMX_SIMD_INT32_WIDTH         GMX_SIMD_DINT32_WIDTH
#    endif
#    ifdef GMX_SIMD_HAVE_DINT32_EXTRACT
#        define GMX_SIMD_HAVE_INT32_EXTRACT
#    endif
#    ifdef GMX_SIMD_HAVE_DINT32_LOGICAL
#        define GMX_SIMD_HAVE_INT32_LOGICAL
#    endif
#    ifdef GMX_SIMD_HAVE_DINT32_ARITHMETICS
#        define GMX_SIMD_HAVE_INT32_ARITHMETICS
#    endif
#    ifdef GMX_SIMD4_HAVE_DOUBLE
#        define GMX_SIMD4_HAVE_REAL
#    endif

#else /* GMX_DOUBLE */

/* Float floating-point */
#    define gmx_simd_real_t                  gmx_simd_float_t
#    define gmx_simd_load_r                  gmx_simd_load_f
#    define gmx_simd_load1_r                 gmx_simd_load1_f
#    define gmx_simd_set1_r                  gmx_simd_set1_f
#    define gmx_simd_store_r                 gmx_simd_store_f
#    define gmx_simd_loadu_r                 gmx_simd_loadu_f
#    define gmx_simd_storeu_r                gmx_simd_storeu_f
#    define gmx_simd_setzero_r               gmx_simd_setzero_f
#    define gmx_simd_add_r                   gmx_simd_add_f
#    define gmx_simd_sub_r                   gmx_simd_sub_f
#    define gmx_simd_mul_r                   gmx_simd_mul_f
#    define gmx_simd_fmadd_r                 gmx_simd_fmadd_f
#    define gmx_simd_fmsub_r                 gmx_simd_fmsub_f
#    define gmx_simd_fnmadd_r                gmx_simd_fnmadd_f
#    define gmx_simd_fnmsub_r                gmx_simd_fnmsub_f
#    define gmx_simd_and_r                   gmx_simd_and_f
#    define gmx_simd_andnot_r                gmx_simd_andnot_f
#    define gmx_simd_or_r                    gmx_simd_or_f
#    define gmx_simd_xor_r                   gmx_simd_xor_f
#    define gmx_simd_rsqrt_r                 gmx_simd_rsqrt_f
#    define gmx_simd_rcp_r                   gmx_simd_rcp_f
#    define gmx_simd_fabs_r                  gmx_simd_fabs_f
#    define gmx_simd_fneg_r                  gmx_simd_fneg_f
#    define gmx_simd_max_r                   gmx_simd_max_f
#    define gmx_simd_min_r                   gmx_simd_min_f
#    define gmx_simd_round_r                 gmx_simd_round_f
#    define gmx_simd_trunc_r                 gmx_simd_trunc_f
#    define gmx_simd_fraction_r              gmx_simd_fraction_f
#    define gmx_simd_get_exponent_r          gmx_simd_get_exponent_f
#    define gmx_simd_get_mantissa_r          gmx_simd_get_mantissa_f
#    define gmx_simd_set_exponent_r          gmx_simd_set_exponent_f
/* Float integer and conversions */
#    define gmx_simd_int32_t                 gmx_simd_fint32_t
#    define gmx_simd_load_i                  gmx_simd_load_fi
#    define gmx_simd_set1_i                  gmx_simd_set1_fi
#    define gmx_simd_store_i                 gmx_simd_store_fi
#    define gmx_simd_loadu_i                 gmx_simd_loadu_fi
#    define gmx_simd_storeu_i                gmx_simd_storeu_fi
#    define gmx_simd_setzero_i               gmx_simd_setzero_fi
#    define gmx_simd_cvt_r2i                 gmx_simd_cvt_f2i
#    define gmx_simd_cvtt_r2i                gmx_simd_cvtt_f2i
#    define gmx_simd_cvt_i2r                 gmx_simd_cvt_i2f
#    define gmx_simd_extract_i               gmx_simd_extract_fi
#    define gmx_simd_slli_i                  gmx_simd_slli_fi
#    define gmx_simd_srli_i                  gmx_simd_srli_fi
#    define gmx_simd_and_i                   gmx_simd_and_fi
#    define gmx_simd_andnot_i                gmx_simd_andnot_fi
#    define gmx_simd_or_i                    gmx_simd_or_fi
#    define gmx_simd_xor_i                   gmx_simd_xor_fi
#    define gmx_simd_add_i                   gmx_simd_add_fi
#    define gmx_simd_sub_i                   gmx_simd_sub_fi
#    define gmx_simd_mul_i                   gmx_simd_mul_fi
/* Float booleans and selection */
#    define gmx_simd_bool_t                  gmx_simd_fbool_t
#    define gmx_simd_cmpeq_r                 gmx_simd_cmpeq_f
#    define gmx_simd_cmplt_r                 gmx_simd_cmplt_f
#    define gmx_simd_cmple_r                 gmx_simd_cmple_f
#    define gmx_simd_and_b                   gmx_simd_and_fb
#    define gmx_simd_or_b                    gmx_simd_or_fb
#    define gmx_simd_anytrue_b               gmx_simd_anytrue_fb
#    define gmx_simd_blendzero_r             gmx_simd_blendzero_f
#    define gmx_simd_blendv_r                gmx_simd_blendv_f
#    define gmx_simd_ibool_t                 gmx_simd_fibool_t
#    define gmx_simd_cmpeq_i                 gmx_simd_cmpeq_fi
#    define gmx_simd_cmplt_i                 gmx_simd_cmplt_fi
#    define gmx_simd_cmple_i                 gmx_simd_cmple_fi
#    define gmx_simd_and_ib                  gmx_simd_and_fib
#    define gmx_simd_or_ib                   gmx_simd_or_fib
#    define gmx_simd_anytrue_ib              gmx_simd_anytrue_fib
#    define gmx_simd_blendzero_i             gmx_simd_blendzero_fi
#    define gmx_simd_blendv_i                gmx_simd_blendv_fi
/* Conversions between integer and single floating-point booleans */
#    define gmx_simd_cvt_b2ib                gmx_simd_cvt_fb2fib
#    define gmx_simd_cvt_ib2b                gmx_simd_cvt_fib2fb
/* SIMD4 float fp - we only support a subset of instructions for SIMD4 */
#    define gmx_simd4_real_t                 gmx_simd4_float_t
#    define gmx_simd4_load_r                 gmx_simd4_load_f
#    define gmx_simd4_load1_r                gmx_simd4_load1_f
#    define gmx_simd4_set1_r                 gmx_simd4_set1_f
#    define gmx_simd4_store_r                gmx_simd4_store_f
#    define gmx_simd4_loadu_r                gmx_simd4_loadu_f
#    define gmx_simd4_storeu_r               gmx_simd4_storeu_f
#    define gmx_simd4_setzero_r              gmx_simd4_setzero_f
#    define gmx_simd4_add_r                  gmx_simd4_add_f
#    define gmx_simd4_sub_r                  gmx_simd4_sub_f
#    define gmx_simd4_mul_r                  gmx_simd4_mul_f
#    define gmx_simd4_fmadd_r                gmx_simd4_fmadd_f
#    define gmx_simd4_fmsub_r                gmx_simd4_fmsub_f
#    define gmx_simd4_fnmadd_r               gmx_simd4_fnmadd_f
#    define gmx_simd4_fnmsub_r               gmx_simd4_fnmsub_f
#    define gmx_simd4_and_r                  gmx_simd4_and_f
#    define gmx_simd4_andnot_r               gmx_simd4_andnot_f
#    define gmx_simd4_or_r                   gmx_simd4_or_f
#    define gmx_simd4_xor_r                  gmx_simd4_xor_f
#    define gmx_simd4_rsqrt_r                gmx_simd4_rsqrt_f
#    define gmx_simd4_fabs_r                 gmx_simd4_fabs_f
#    define gmx_simd4_fneg_r                 gmx_simd4_fneg_f
#    define gmx_simd4_max_r                  gmx_simd4_max_f
#    define gmx_simd4_min_r                  gmx_simd4_min_f
#    define gmx_simd4_round_r                gmx_simd4_round_f
#    define gmx_simd4_trunc_r                gmx_simd4_trunc_f
#    define gmx_simd4_dotproduct3_r          gmx_simd4_dotproduct3_f
#    define gmx_simd4_bool_t                 gmx_simd4_fbool_t
#    define gmx_simd4_cmpeq_r                gmx_simd4_cmpeq_f
#    define gmx_simd4_cmplt_r                gmx_simd4_cmplt_f
#    define gmx_simd4_cmple_r                gmx_simd4_cmple_f
#    define gmx_simd4_and_b                  gmx_simd4_and_fb
#    define gmx_simd4_or_b                   gmx_simd4_or_fb
#    define gmx_simd4_anytrue_b              gmx_simd4_anytrue_fb
#    define gmx_simd4_blendzero_r            gmx_simd4_blendzero_f
#    define gmx_simd4_blendv_r               gmx_simd4_blendv_f

/* Memory allocation */
#    define gmx_simd_align_r                 gmx_simd_align_f
#    define gmx_simd_align_i                 gmx_simd_align_fi
#    define gmx_simd4_align_r                gmx_simd4_align_f

#    ifdef GMX_SIMD_HAVE_FLOAT
#        define GMX_SIMD_HAVE_REAL
#        define GMX_SIMD_REAL_WIDTH          GMX_SIMD_FLOAT_WIDTH
#    endif
#    ifdef GMX_SIMD_HAVE_FINT32
#        define GMX_SIMD_HAVE_INT32
#        define GMX_SIMD_INT32_WIDTH         GMX_SIMD_FINT32_WIDTH
#    endif
#    ifdef GMX_SIMD_HAVE_FINT32_EXTRACT
#        define GMX_SIMD_HAVE_INT32_EXTRACT
#    endif
#    ifdef GMX_SIMD_HAVE_FINT32_LOGICAL
#        define GMX_SIMD_HAVE_INT32_LOGICAL
#    endif
#    ifdef GMX_SIMD_HAVE_FINT32_ARITHMETICS
#        define GMX_SIMD_HAVE_INT32_ARITHMETICS
#    endif
#    ifdef GMX_SIMD4_HAVE_FLOAT
#        define GMX_SIMD4_HAVE_REAL
#    endif

#endif /* GMX_DOUBLE */

#include "simd_math.h"

#endif /* GMX_SIMD_SIMD_H */
