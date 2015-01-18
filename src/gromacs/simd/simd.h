/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

/*! \libinternal
 * \defgroup module_simd SIMD intrinsics interface (simd)
 * \ingroup group_utilitymodules
 *
 * \brief Provides an architecture-independent way of doing SIMD coding.
 *
 * Overview of the SIMD implementation is provided in \ref page_simd.
 * The details are documented in simd.h and the reference implementation
 * impl_reference.h.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 */

#ifndef GMX_SIMD_SIMD_H
#define GMX_SIMD_SIMD_H

/*! \libinternal \file
 *
 * \brief Definitions, capabilities, and wrappers for SIMD module.
 *
 * The macros in this file are intended to be used for writing
 * architecture-independent SIMD intrinsics code.
 * To support a new architecture, adding a new sub-include with macros here
 * should be (nearly) all that is needed.
 *
 * The defines in this top-level file will set default Gromacs real precision
 * operations to either single or double precision based on whether
 * GMX_DOUBLE is defined. The actual implementation - including e.g.
 * conversion operations specifically between single and double - is documented
 * in impl_reference.h.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

#include "config.h"

#include <stddef.h>

#include "gromacs/utility/basedefinitions.h"

/* Forward declarations so memory allocation can be used in implementations */
static gmx_inline float *  gmx_simd_align_f(float *p);
static gmx_inline double * gmx_simd_align_d(double *p);
static gmx_inline int *    gmx_simd_align_fi(int *p);
static gmx_inline int *    gmx_simd_align_di(int *p);
static gmx_inline float *  gmx_simd4_align_f(float *p);
static gmx_inline double * gmx_simd4_align_d(double *p);

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name SIMD predefined macros to describe high-level capabilities
 *
 *  These macros are used to describe the features available in default
 *  Gromacs real precision. They are set from the lower-level implementation
 *  files that have macros describing single and double precision individually,
 *  as well as the implementation details.
 *  \{
 */

/*! \brief
 *  GMX_SIMD indicates that some sort of SIMD support is present in software.
 *
 * It is disabled if no architecture, neither reference SIMD, has been selected.
 */
#define GMX_SIMD


/* Intel MIC is a bit special since it is a co-processor. This means the rest
 * of GROMACS (which runs on the CPU) can use a default SIMD set like AVX.
 * All functions in this SIMD module are static, so it will work perfectly fine
 * to include this file with different SIMD definitions for different files.
 */
#if defined GMX_SIMD_X86_AVX_512ER
#    include "impl_x86_avx_512er/impl_x86_avx_512er.h"
#elif defined GMX_SIMD_X86_AVX_512F
#    include "impl_x86_avx_512f/impl_x86_avx_512f.h"
#elif defined GMX_SIMD_X86_MIC
#    include "impl_intel_mic/impl_intel_mic.h"
#elif defined GMX_SIMD_X86_AVX2_256
#    include "impl_x86_avx2_256/impl_x86_avx2_256.h"
#elif defined GMX_SIMD_X86_AVX_256
#    include "impl_x86_avx_256/impl_x86_avx_256.h"
#elif defined GMX_SIMD_X86_AVX_128_FMA
#    include "impl_x86_avx_128_fma/impl_x86_avx_128_fma.h"
#elif defined GMX_SIMD_X86_SSE4_1
#    include "impl_x86_sse4_1/impl_x86_sse4_1.h"
#elif defined GMX_SIMD_X86_SSE2
#    include "impl_x86_sse2/impl_x86_sse2.h"
#elif defined GMX_SIMD_ARM_NEON
#    include "impl_arm_neon/impl_arm_neon.h"
#elif defined GMX_SIMD_ARM_NEON_ASIMD
#    include "impl_arm_neon_asimd/impl_arm_neon_asimd.h"
#elif defined GMX_SIMD_IBM_QPX
#    include "impl_ibm_qpx/impl_ibm_qpx.h"
#elif defined GMX_SIMD_IBM_VMX
#    include "impl_ibm_vmx/impl_ibm_vmx.h"
#elif defined GMX_SIMD_IBM_VSX
#    include "impl_ibm_vsx/impl_ibm_vsx.h"
#elif defined GMX_SIMD_SPARC64_HPC_ACE
#    include "impl_sparc64_hpc_ace/impl_sparc64_hpc_ace.h"
#elif (defined GMX_SIMD_REFERENCE) || (defined DOXYGEN)
/* Plain C SIMD reference implementation, also serves as documentation.
 * For now this code path will also be taken for Sparc64_HPC_ACE since we have
 * not yet added the verlet kernel extensions there. The group kernels do not
 * depend on this file, so they will still be accelerated with SIMD.
 */
#    include "impl_reference/impl_reference.h"
#else
/* Turn off the GMX_SIMD flag if we do not even have reference support */
#    undef GMX_SIMD
#endif

/*! \brief
 * SIMD4 width is always 4, but use this for clarity in definitions.
 *
 * It improves code readability to allocate e.g. 2*GMX_SIMD4_WIDTH instead of 8.
 */
#define GMX_SIMD4_WIDTH    4

/*! \} */

/*! \name SIMD memory alignment operations
 *  \{
 */

/*! \brief
 * Align a float pointer for usage with SIMD instructions.
 *
 * You should typically \a not call this function directly (unless you explicitly
 * want single precision even when GMX_DOUBLE is set), but use the
 * \ref gmx_simd_align_r macro to align memory in default Gromacs real precision.
 *
 * \param  p Pointer to memory, allocate at least \ref GMX_SIMD_FLOAT_WIDTH extra elements.
 *
 * \return Aligned pointer (>=p) suitable for loading/storing float fp SIMD.
 *         If \ref GMX_SIMD_HAVE_FLOAT is not set, p will be returned unchanged.
 *
 * Start by allocating an extra \ref GMX_SIMD_FLOAT_WIDTH float elements of memory,
 * and then call this function. The returned pointer will be greater or equal
 * to the one you provided, and point to an address inside your provided memory
 * that is aligned to the SIMD width.
 */
static gmx_inline float *
gmx_simd_align_f(float *p)
{
#    ifdef GMX_SIMD_HAVE_FLOAT
    return (float *)(((size_t)((p)+GMX_SIMD_FLOAT_WIDTH-1)) & (~((size_t)(GMX_SIMD_FLOAT_WIDTH*sizeof(float)-1))));
#    else
    return p;
#    endif
}

/*!  \brief
 * Align a double pointer for usage with SIMD instructions.
 *
 * You should typically \a not call this function directly (unless you explicitly
 * want double precision even when GMX_DOUBLE is not set), but use the
 * \ref gmx_simd_align_r macro to align memory in default Gromacs real precision.
 *
 * \param  p Pointer to memory, allocate at least \ref GMX_SIMD_DOUBLE_WIDTH extra elements.
 *
 * \return Aligned pointer (>=p) suitable for loading/storing double fp SIMD.
 *         If \ref GMX_SIMD_HAVE_DOUBLE is not set, p will be returned unchanged.
 *
 * Start by allocating an extra \ref GMX_SIMD_DOUBLE_WIDTH double elements of memory,
 * and then call this function. The returned pointer will be greater or equal
 * to the one you provided, and point to an address inside your provided memory
 * that is aligned to the SIMD width.
 */
static gmx_inline double *
gmx_simd_align_d(double *p)
{
#    ifdef GMX_SIMD_HAVE_DOUBLE
    return (double *)(((size_t)((p)+GMX_SIMD_DOUBLE_WIDTH-1)) & (~((size_t)(GMX_SIMD_DOUBLE_WIDTH*sizeof(double)-1))));
#    else
    return p;
#    endif
}

/*! \brief
 * Align a (float) integer pointer for usage with SIMD instructions.
 *
 * You should typically \a not call this function directly (unless you explicitly
 * want integers corresponding to single precision even when GMX_DOUBLE is
 * set), but use the \ref gmx_simd_align_i macro to align integer memory
 * corresponding to Gromacs default floating-point precision.
 *
 * \param  p Pointer to memory, allocate at least \ref GMX_SIMD_FINT32_WIDTH extra elements.
 *
 * \return Aligned pointer (>=p) suitable for loading/storing float-integer SIMD.
 *         If \ref GMX_SIMD_HAVE_FINT32 is not set, p will be returned unchanged.
 *
 * This routine provides aligned memory for usage with \ref gmx_simd_fint32_t. You
 * should have allocated an extra \ref GMX_SIMD_FINT32_WIDTH * sizeof(int) bytes. The
 * reason why we need to separate float-integer vs. double-integer is that the
 * width of registers after conversions from the floating-point types might not
 * be identical, or even supported, in both cases.
 */
static gmx_inline int *
gmx_simd_align_fi(int *p)
{
#    ifdef GMX_SIMD_HAVE_FINT32
    return (int *)(((size_t)((p)+GMX_SIMD_FINT32_WIDTH-1)) & (~((size_t)(GMX_SIMD_FINT32_WIDTH*sizeof(int)-1))));
#    else
    return p;
#    endif
}

/*! \brief
 * Align a (double) integer pointer for usage with SIMD instructions.
 *
 * You should typically \a not call this function directly (unless you explicitly
 * want integers corresponding to doublele precision even when GMX_DOUBLE is
 * not set), but use the \ref gmx_simd_align_i macro to align integer memory
 * corresponding to Gromacs default floating-point precision.
 *
 * \param  p Pointer to memory, allocate at least \ref GMX_SIMD_DINT32_WIDTH extra elements.
 *
 * \return Aligned pointer (>=p) suitable for loading/storing double-integer SIMD.
 *         If \ref GMX_SIMD_HAVE_DINT32 is not set, p will be returned unchanged.
 *
 * This routine provides aligned memory for usage with \ref gmx_simd_dint32_t. You
 * should have allocated an extra \ref GMX_SIMD_DINT32_WIDTH*sizeof(int) bytes. The
 * reason why we need to separate float-integer vs. double-integer is that the
 * width of registers after conversions from the floating-point types might not
 * be identical, or even supported, in both cases.
 */
static gmx_inline int *
gmx_simd_align_di(int *p)
{
#    ifdef GMX_SIMD_HAVE_DINT32
    return (int *)(((size_t)((p)+GMX_SIMD_DINT32_WIDTH-1)) & (~((size_t)(GMX_SIMD_DINT32_WIDTH*sizeof(int)-1))));
#    else
    return p;
#    endif
}

/*! \brief
 * Align a float pointer for usage with SIMD4 instructions.
 *
 * You should typically \a not call this function directly (unless you explicitly
 * want single precision even when GMX_DOUBLE is set), but use the
 * \ref gmx_simd4_align_r macro to align memory in default Gromacs real precision.
 *
 * \param  p Pointer to memory, allocate at least \ref GMX_SIMD4_WIDTH extra elements.
 *
 * \return Aligned pointer (>=p) suitable for loading/storing float SIMD.
 *         If \ref GMX_SIMD4_HAVE_FLOAT is not set, p will be returned unchanged.
 *
 * This routine provides aligned memory for usage with \ref gmx_simd4_float_t.
 * should have allocated an extra \ref GMX_SIMD4_WIDTH * sizeof(float) bytes.
 */
static gmx_inline float *
gmx_simd4_align_f(float *p)
{
#    ifdef GMX_SIMD4_HAVE_FLOAT
    return (float *)(((size_t)((p)+GMX_SIMD4_WIDTH-1)) & (~((size_t)(GMX_SIMD4_WIDTH*sizeof(float)-1))));
#    else
    return p;
#    endif
}

/*! \brief
 * Align a double pointer for usage with SIMD4 instructions.
 *
 * You should typically \a not call this function directly (unless you explicitly
 * want double precision even when GMX_DOUBLE is not set), but use the
 * \ref gmx_simd4_align_r macro to align memory in default Gromacs real precision.
 *
 * \param  p Pointer to memory, allocate at least \ref GMX_SIMD4_WIDTH extra elements.
 *
 * \return Aligned pointer (>=p) suitable for loading/storing float SIMD.
 *         If \ref GMX_SIMD4_HAVE_DOUBLE is not set, p will be returned unchanged.
 *
 * This routine provides aligned memory for usage with \ref gmx_simd4_double_t.
 * should have allocated an extra \ref GMX_SIMD4_WIDTH * sizeof(double) bytes.
 */
static gmx_inline double *
gmx_simd4_align_d(double *p)
{
#    ifdef GMX_SIMD4_HAVE_DOUBLE
    return (double *)(((size_t)((p)+GMX_SIMD4_WIDTH-1)) & (~((size_t)(GMX_SIMD4_WIDTH*sizeof(double)-1))));
#    else
    return p;
#    endif
}

/*! \} */


/* Define Gromacs "real" precision macros depending on Gromacs config. Note
 * that conversions float-to-double and v.v. are not included here since they
 * are not precision-dependent - find them in the implementation files.
 */
#ifdef GMX_DOUBLE
/* Double floating-point. The documentation is in the float part below */
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
#    define gmx_simd_blendnotzero_r          gmx_simd_blendnotzero_d
#    define gmx_simd_blendv_r                gmx_simd_blendv_d
#    define gmx_simd_reduce_r                gmx_simd_reduce_d
#    define gmx_simd_ibool_t                 gmx_simd_dibool_t
#    define gmx_simd_cmpeq_i                 gmx_simd_cmpeq_di
#    define gmx_simd_cmplt_i                 gmx_simd_cmplt_di
#    define gmx_simd_and_ib                  gmx_simd_and_dib
#    define gmx_simd_or_ib                   gmx_simd_or_dib
#    define gmx_simd_anytrue_ib              gmx_simd_anytrue_dib
#    define gmx_simd_blendzero_i             gmx_simd_blendzero_di
#    define gmx_simd_blendnotzero_i          gmx_simd_blendnotzero_di
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
#    define gmx_simd4_blendnotzero_r         gmx_simd4_blendnotzero_d
#    define gmx_simd4_blendv_r               gmx_simd4_blendv_d
#    define gmx_simd4_reduce_r               gmx_simd4_reduce_d

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

/*! \name SIMD data types
 *
 *  The actual storage of these types is implementation dependent. The
 *  documentation is generated from the reference implementation, but for
 *  normal usage this will likely not be what you are using.
 * \{
 */
/*! \brief Real precision floating-point SIMD datatype.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_REAL is defined.
 *
 * If GMX_DOUBLE is defined, this will be set to \ref gmx_simd_double_t
 * internally, otherwise \ref gmx_simd_float_t.
 */
#    define gmx_simd_real_t                  gmx_simd_float_t

/*! \brief 32-bit integer SIMD type.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_INT32 is defined.
 *
 * If GMX_DOUBLE is defined, this will be set to \ref gmx_simd_dint32_t
 * internally, otherwise \ref gmx_simd_fint32_t. This might seem a strange
 * implementation detail, but it is because some SIMD implementations use
 * different types/widths of integers registers when converting from
 * double vs. single precision floating point. As long as you just use
 * this type you will not have to worry about precision.
 */
#    define gmx_simd_int32_t                 gmx_simd_fint32_t

/*! \brief Boolean SIMD type for usage with \ref gmx_simd_real_t.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_REAL is defined.
 *
 * If GMX_DOUBLE is defined, this will be set to \ref gmx_simd_dbool_t
 * internally, otherwise \ref gmx_simd_fbool_t. This is necessary since some
 * SIMD implementations use bitpatterns for marking truth, so single-
 * vs. double precision booleans are not necessarily exchangable.
 * As long as you just use this type you will not have to worry about precision.
 *
 * See \ref gmx_simd_ibool_t for an explanation of real vs. integer booleans.
 */
#    define gmx_simd_bool_t                  gmx_simd_fbool_t

/*! \brief Boolean SIMD type for usage with \ref gmx_simd_int32_t.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_INT32 is defined.
 *
 * If GMX_DOUBLE is defined, this will be set to \ref gmx_simd_dibool_t
 * internally, otherwise \ref gmx_simd_fibool_t. This is necessary since some
 * SIMD implementations use bitpatterns for marking truth, so single-
 * vs. double precision booleans are not necessarily exchangable, and while
 * a double-precision boolean might be represented with a 64-bit mask, the
 * corresponding integer might only use a 32-bit mask.
 *
 * We provide conversion routines for these cases, so the only thing you need to
 * keep in mind is to use \ref gmx_simd_bool_t when working with
 * \ref gmx_simd_real_t while you pick \ref gmx_simd_ibool_t when working with
 * \ref gmx_simd_int32_t.
 *
 * To convert between them, use \ref gmx_simd_cvt_b2ib and \ref gmx_simd_cvt_ib2b.
 */
#    define gmx_simd_ibool_t                 gmx_simd_fibool_t


/*! \}
 *  \name SIMD load/store operations on gmx_simd_real_t
 *
 *  \note Unaligned load/stores are only available when
 *  \ref GMX_SIMD_HAVE_LOADU and \ref GMX_SIMD_HAVE_STOREU are set, respectively.
 *  \{
 */

/*! \brief Load \ref GMX_SIMD_REAL_WIDTH values from aligned memory to \ref gmx_simd_real_t
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_load_d,
 * otherwise \ref gmx_simd_load_f.
 *
 * \copydetails gmx_simd_load_f
 */
#    define gmx_simd_load_r                  gmx_simd_load_f

/*! \brief Set all elements in \ref gmx_simd_real_t from single value in memory.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_load1_d,
 * otherwise \ref gmx_simd_load1_f.
 *
 * \copydetails gmx_simd_load1_f
 */
#    define gmx_simd_load1_r                 gmx_simd_load1_f

/*! \brief Set all elements in \ref gmx_simd_real_t from a scalar.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_set1_d,
 * otherwise \ref gmx_simd_set1_f.
 *
 * \copydetails gmx_simd_set1_f
 */
#    define gmx_simd_set1_r                  gmx_simd_set1_f

/*! \brief Store \ref GMX_SIMD_REAL_WIDTH values from \ref gmx_simd_real_t to aligned memory.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_store_d,
 * otherwise \ref gmx_simd_store_f.
 *
 * \copydetails gmx_simd_store_f
 */
#    define gmx_simd_store_r                 gmx_simd_store_f

/*! \brief Load \ref GMX_SIMD_REAL_WIDTH values from unaligned memory to \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_loadu_d,
 * otherwise \ref gmx_simd_loadu_f.
 *
 * \copydetails gmx_simd_loadu_f
 */
#    define gmx_simd_loadu_r                 gmx_simd_loadu_f

/*! \brief Store \ref GMX_SIMD_REAL_WIDTH values from \ref gmx_simd_real_t to unaligned memory.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_storeu_d,
 * otherwise \ref gmx_simd_storeu_f.
 *
 * \copydetails gmx_simd_storeu_f
 */
#    define gmx_simd_storeu_r                gmx_simd_storeu_f

/*! \brief Set all elements in \ref gmx_simd_real_t to 0.0.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_setzero_d,
 * otherwise \ref gmx_simd_setzero_f.
 *
 * \copydetails gmx_simd_setzero_f
 */
#    define gmx_simd_setzero_r               gmx_simd_setzero_f

/*! \}
 *  \name SIMD load/store operations on gmx_simd_int32_t
 *
 *  \note Unaligned load/stores are only available when
 *  \ref GMX_SIMD_HAVE_LOADU and \ref GMX_SIMD_HAVE_STOREU are set, respectively.
 *  \{
 */

/*! \brief Load \ref GMX_SIMD_INT32_WIDTH values from aligned memory to \ref gmx_simd_int32_t .
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_load_di ,
 * otherwise \ref gmx_simd_load_fi .
 *
 * \copydetails gmx_simd_load_fi
 */
#    define gmx_simd_load_i                  gmx_simd_load_fi

/*! \brief Set all elements in \ref gmx_simd_int32_t from a single integer.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_set1_di ,
 * otherwise \ref gmx_simd_set1_fi .
 *
 * \copydetails gmx_simd_set1_fi
 */
#    define gmx_simd_set1_i                  gmx_simd_set1_fi

/*! \brief Store \ref GMX_SIMD_REAL_WIDTH values from \ref gmx_simd_int32_t to aligned memory.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_store_di ,
 * otherwise \ref gmx_simd_store_fi .
 *
 * \copydetails gmx_simd_store_fi
 */
#    define gmx_simd_store_i                 gmx_simd_store_fi

/*! \brief Load \ref GMX_SIMD_REAL_WIDTH values from unaligned memory to \ref gmx_simd_int32_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_loadu_di ,
 * otherwise \ref gmx_simd_loadu_fi .
 *
 * \copydetails gmx_simd_loadu_fi
 */
#    define gmx_simd_loadu_i                 gmx_simd_loadu_fi

/*! \brief Store \ref GMX_SIMD_REAL_WIDTH values from \ref gmx_simd_int32_t to unaligned memory.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_storeu_di ,
 * otherwise \ref gmx_simd_storeu_fi .
 *
 * \copydetails gmx_simd_storeu_fi
 */
#    define gmx_simd_storeu_i                gmx_simd_storeu_fi

/*! \brief Extract single integer from \ref gmx_simd_int32_t element.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_extract_di ,
 * otherwise \ref gmx_simd_extract_fi .
 *
 * \copydetails gmx_simd_extract_fi
 */
#    define gmx_simd_extract_i               gmx_simd_extract_fi

/*! \brief Set all elements in \ref gmx_simd_int32_t to 0.
 *
 * If GMX_DOUBLE is defined, it will be aliased to \ref gmx_simd_setzero_di ,
 * otherwise \ref gmx_simd_setzero_fi .
 *
 * \copydetails gmx_simd_setzero_fi
 */
#    define gmx_simd_setzero_i               gmx_simd_setzero_fi


/*! \}
 *  \name SIMD floating-point logical operations on gmx_simd_real_t
 *
 *  These instructions are available if \ref GMX_SIMD_HAVE_LOGICAL is defined.
 *  \{
 */

/*! \brief Bitwise \a and on two \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_and_d,
 * otherwise \ref gmx_simd_and_f.
 *
 * \copydetails gmx_simd_and_f
 */
#    define gmx_simd_and_r                   gmx_simd_and_f

/*! \brief Bitwise \a and-not on two \ref gmx_simd_real_t; 1st arg is complemented.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_andnot_d,
 * otherwise \ref gmx_simd_andnot_f.
 *
 * \copydetails gmx_simd_andnot_f
 */
#    define gmx_simd_andnot_r                gmx_simd_andnot_f

/*! \brief Bitwise \a or on two \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_or_d,
 * otherwise \ref gmx_simd_or_f.
 *
 * \copydetails gmx_simd_or_f
 */
#    define gmx_simd_or_r                    gmx_simd_or_f

/*! \brief Bitwise \a exclusive-or on two \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_xor_d,
 * otherwise \ref gmx_simd_xor_f.
 *
 * \copydetails gmx_simd_xor_f
 */
#    define gmx_simd_xor_r                   gmx_simd_xor_f

/*! \}
 *  \name SIMD floating-point arithmetic operations on gmx_simd_real_t
 *  \{
 */

/*! \brief SIMD a+b for two \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_add_d,
 * otherwise \ref gmx_simd_add_f.
 *
 * \copydetails gmx_simd_add_f
 */
#    define gmx_simd_add_r                   gmx_simd_add_f

/*! \brief SIMD a-b for two \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_sub_d,
 * otherwise \ref gmx_simd_sub_f.
 *
 * \copydetails gmx_simd_sub_f
 */
#    define gmx_simd_sub_r                   gmx_simd_sub_f

/*! \brief SIMD a*b for two \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_mul_d,
 * otherwise \ref gmx_simd_mul_f.
 *
 * \copydetails gmx_simd_mul_f
 */
#    define gmx_simd_mul_r                   gmx_simd_mul_f

/*! \brief SIMD a*b+c for three \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_fmadd_d,
 * otherwise \ref gmx_simd_fmadd_f.
 *
 * \copydetails gmx_simd_fmadd_f
 */
#    define gmx_simd_fmadd_r                 gmx_simd_fmadd_f

/*! \brief SIMD a*b-c for three \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_fmsub_d,
 * otherwise \ref gmx_simd_fmsub_f.
 *
 * \copydetails gmx_simd_fmsub_f
 */
#    define gmx_simd_fmsub_r                 gmx_simd_fmsub_f

/*! \brief SIMD -a*b+c for three \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_fnmadd_d,
 * otherwise \ref gmx_simd_fnmadd_f.
 *
 * \copydetails gmx_simd_fnmadd_f
 */
#    define gmx_simd_fnmadd_r                gmx_simd_fnmadd_f

/*! \brief SIMD -a*b-c for three \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_fnmsub_d,
 * otherwise \ref gmx_simd_fnmsub_f.
 *
 * \copydetails gmx_simd_fnmsub_f
 */
#    define gmx_simd_fnmsub_r                gmx_simd_fnmsub_f

/*! \brief SIMD table lookup for 1/sqrt(x) approximation.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_rsqrt_d,
 * otherwise \ref gmx_simd_rsqrt_f.
 *
 * \copydetails gmx_simd_rsqrt_f
 */
#    define gmx_simd_rsqrt_r                 gmx_simd_rsqrt_f

/*! \brief SIMD table lookup for 1/x approximation.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_rcp_d,
 * otherwise \ref gmx_simd_rcp_f.
 *
 * \copydetails gmx_simd_rcp_f
 */
#    define gmx_simd_rcp_r                   gmx_simd_rcp_f

/*! \brief SIMD fabs(x) for \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_fabs_d,
 * otherwise \ref gmx_simd_fabs_f.
 *
 * \copydetails gmx_simd_fabs_f
 */
#    define gmx_simd_fabs_r                  gmx_simd_fabs_f

/*! \brief SIMD -x for \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_fneg_d,
 * otherwise \ref gmx_simd_fneg_f.
 *
 * \copydetails gmx_simd_fneg_f
 */
#    define gmx_simd_fneg_r                  gmx_simd_fneg_f

/*! \brief SIMD max(a,b) for each element in \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_max_d,
 * otherwise \ref gmx_simd_max_f.
 *
 * \copydetails gmx_simd_max_f
 */
#    define gmx_simd_max_r                   gmx_simd_max_f

/*! \brief SIMD min(a,b) for each element in \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_min_d,
 * otherwise \ref gmx_simd_min_f.
 *
 * \copydetails gmx_simd_min_f
 */
#    define gmx_simd_min_r                   gmx_simd_min_f

/*! \brief Round \ref gmx_simd_real_t to nearest int, return \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_round_d,
 * otherwise \ref gmx_simd_round_f.
 *
 * \copydetails gmx_simd_round_f
 */
#    define gmx_simd_round_r                 gmx_simd_round_f

/*! \brief Truncate \ref gmx_simd_real_t towards 0, return \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_trunc_d,
 * otherwise \ref gmx_simd_trunc_f.
 *
 * \copydetails gmx_simd_trunc_f
 */
#    define gmx_simd_trunc_r                 gmx_simd_trunc_f

/*! \brief SIMD Fraction, i.e. x-trunc(x) for \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_fraction_d,
 * otherwise \ref gmx_simd_fraction_f.
 *
 * \copydetails gmx_simd_fraction_f
 */
#    define gmx_simd_fraction_r              gmx_simd_fraction_f

/*! \brief Return the FP exponent of a SIMD \ref gmx_simd_real_t as a \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_get_exponent_d,
 * otherwise \ref gmx_simd_get_exponent_f.
 *
 * \copydetails gmx_simd_exponent_f
 */
#    define gmx_simd_get_exponent_r          gmx_simd_get_exponent_f

/*! \brief Return the FP mantissa of a SIMD \ref gmx_simd_real_t as a \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_get_mantissa_d,
 * otherwise \ref gmx_simd_get_mantissa_f.
 *
 * \copydetails gmx_simd_mantissa_f
 */
#    define gmx_simd_get_mantissa_r          gmx_simd_get_mantissa_f

/*! \brief Set the exponent of a SIMD \ref gmx_simd_real_t from a \ref gmx_simd_real_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_set_exponent_d,
 * otherwise \ref gmx_simd_set_exponent_f.
 *
 * \copydetails gmx_simd_set_exponent_f
 */
#    define gmx_simd_set_exponent_r          gmx_simd_set_exponent_f

/*! \}
 *  \name SIMD comparison, boolean, and select operations for gmx_simd_real_t
 *  \{
 */

/*! \brief SIMD a==b for \ref gmx_simd_real_t. Returns a \ref gmx_simd_bool_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cmpeq_d,
 * otherwise \ref gmx_simd_cmpeq_f.
 *
 * \copydetails gmx_simd_cmpeq_f
 */
#    define gmx_simd_cmpeq_r                 gmx_simd_cmpeq_f

/*! \brief SIMD a<b for \ref gmx_simd_real_t. Returns a \ref gmx_simd_bool_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cmplt_d,
 * otherwise \ref gmx_simd_cmplt_f.
 *
 * \copydetails gmx_simd_cmplt_f
 */
#    define gmx_simd_cmplt_r                 gmx_simd_cmplt_f

/*! \brief SIMD a<=b for \ref gmx_simd_real_t. Returns a \ref gmx_simd_bool_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cmple_d,
 * otherwise \ref gmx_simd_cmple_f.
 *
 * \copydetails gmx_simd_cmple_f
 */
#    define gmx_simd_cmple_r                 gmx_simd_cmple_f

/*! \brief For each element, the result boolean is true if both arguments are true
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_and_db,
 * otherwise \ref gmx_simd_and_fb.
 *
 * \copydetails gmx_simd_and_fb
 */
#    define gmx_simd_and_b                   gmx_simd_and_fb

/*! \brief For each element, the result boolean is true if either argument is true
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_or_db,
 * otherwise \ref gmx_simd_or_fb.
 *
 * \copydetails gmx_simd_or_fn
 */
#    define gmx_simd_or_b                    gmx_simd_or_fb

/*! \brief Return nonzero if any element in gmx_simd_bool_t is true, otherwise 0.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_anytrue_db,
 * otherwise \ref gmx_simd_anytrue_fb.
 *
 * \copydetails gmx_simd_anytrue_fb
 */
#    define gmx_simd_anytrue_b               gmx_simd_anytrue_fb

/*! \brief Selects elements from \ref gmx_simd_real_t where boolean is true, otherwise 0.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_blendzero_d,
 * otherwise \ref gmx_simd_blendzero_f.
 *
 * \copydetails gmx_simd_blendzero_f
 *
 * \sa gmx_simd_blendzero_i
 */
#    define gmx_simd_blendzero_r             gmx_simd_blendzero_f

/*! \brief Selects elements from \ref gmx_simd_real_t where boolean is false, otherwise 0.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_blendnotzero_d,
 * otherwise \ref gmx_simd_blendnotzero_f.
 *
 * \copydetails gmx_simd_blendnotzero_f
 */
#    define gmx_simd_blendnotzero_r          gmx_simd_blendnotzero_f

/*! \brief Selects from 2nd real SIMD arg where boolean is true, otherwise 1st arg.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_blendv_d,
 * otherwise \ref gmx_simd_blendv_f.
 *
 * \copydetails gmx_simd_blendv_f
 */
#    define gmx_simd_blendv_r                gmx_simd_blendv_f

/*! \brief Return sum of all elements in SIMD floating-point variable.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_reduce_d,
 * otherwise \ref gmx_simd_reduce_f.
 *
 * \copydetails gmx_simd_reduce_f
 */
#    define gmx_simd_reduce_r                gmx_simd_reduce_f

/*! \}
 *  \name SIMD integer logical operations on gmx_simd_int32_t
 *
 *  These instructions are available if \ref GMX_SIMD_HAVE_INT32_LOGICAL is defined.
 *  \{
 */

/*! \brief Shift each element in \ref gmx_simd_int32_t left by immediate
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_slli_di,
 * otherwise \ref gmx_simd_slli_fi.
 *
 * \copydetails gmx_simd_slli_fi
 */
#    define gmx_simd_slli_i                  gmx_simd_slli_fi

/*! \brief Shift each element in \ref gmx_simd_int32_t right by immediate
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_srli_di,
 * otherwise \ref gmx_simd_srli_fi.
 *
 * \copydetails gmx_simd_srli_fi
 */
#    define gmx_simd_srli_i                  gmx_simd_srli_fi

/*! \brief Bitwise \a and on two \ref gmx_simd_int32_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_and_di,
 * otherwise \ref gmx_simd_and_fi.
 *
 * \copydetails gmx_simd_and_fi
 */
#    define gmx_simd_and_i                   gmx_simd_and_fi

/*! \brief Bitwise \a and-not on two \ref gmx_simd_int32_t; 1st arg is complemented.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_andnot_di,
 * otherwise \ref gmx_simd_andnot_fi.
 *
 * \copydetails gmx_simd_andnot_fi
 */
#    define gmx_simd_andnot_i                gmx_simd_andnot_fi

/*! \brief Bitwise \a or on two \ref gmx_simd_int32_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_or_di,
 * otherwise \ref gmx_simd_or_fi.
 *
 * \copydetails gmx_simd_or_fi
 */
#    define gmx_simd_or_i                    gmx_simd_or_fi

/*! \brief Bitwise \a xor on two \ref gmx_simd_int32_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_xor_di,
 * otherwise \ref gmx_simd_xor_fi.
 *
 * \copydetails gmx_simd_xor_fi
 */
#    define gmx_simd_xor_i                   gmx_simd_xor_fi

/*! \}
 *  \name SIMD integer arithmetic operations on gmx_simd_int32_t
 *
 *  These instructions are available if \ref GMX_SIMD_HAVE_INT32_ARITHMETICS is defined.
 *  \{
 */

/*! \brief SIMD a+b for two \ref gmx_simd_int32_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_add_di,
 * otherwise \ref gmx_simd_add_fi.
 *
 * \copydetails gmx_simd_add_fi
 */
#    define gmx_simd_add_i                   gmx_simd_add_fi

/*! \brief SIMD a-b for two \ref gmx_simd_int32_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_sub_di,
 * otherwise \ref gmx_simd_sub_fi.
 *
 * \copydetails gmx_simd_sub_fi
 */
#    define gmx_simd_sub_i                   gmx_simd_sub_fi

/*! \brief SIMD a*b for two \ref gmx_simd_int32_t.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_mul_di,
 * otherwise \ref gmx_simd_mul_fi.
 *
 * \copydetails gmx_simd_mul_fi
 */
#    define gmx_simd_mul_i                   gmx_simd_mul_fi

/*! \}
 *  \name SIMD integer comparison, booleans, and selection on gmx_simd_int32_t
 *
 *  These instructions are available if \ref GMX_SIMD_HAVE_INT32_ARITHMETICS is defined.
 *  \{
 */

/*! \brief Returns boolean describing whether a==b, for \ref gmx_simd_int32_t
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cmpeq_di,
 * otherwise \ref gmx_simd_cmpeq_fi.
 *
 * \copydetails gmx_simd_cmpeq_fi
 */
#    define gmx_simd_cmpeq_i                 gmx_simd_cmpeq_fi

/*! \brief Returns boolean describing whether a<b, for \ref gmx_simd_int32_t
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cmplt_di,
 * otherwise \ref gmx_simd_cmplt_fi.
 *
 * \copydetails gmx_simd_cmplt_fi
 */
#    define gmx_simd_cmplt_i                 gmx_simd_cmplt_fi

/*! \brief For each element, the result boolean is true if both arguments are true
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_and_dib,
 * otherwise \ref gmx_simd_and_fib.
 *
 * \copydetails gmx_simd_and_fib
 */
#    define gmx_simd_and_ib                  gmx_simd_and_fib

/*! \brief For each element, the result boolean is true if either argument is true.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_or_dib,
 * otherwise \ref gmx_simd_or_fib.
 *
 * \copydetails gmx_simd_or_fib
 */
#    define gmx_simd_or_ib                   gmx_simd_or_fib

/*! \brief Return nonzero if any element in gmx_simd_ibool_t is true, otherwise 0.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_anytrue_dib,
 * otherwise \ref gmx_simd_anytrue_fib.
 *
 * \copydetails gmx_simd_anytrue_fib
 */
#    define gmx_simd_anytrue_ib              gmx_simd_anytrue_fib

/*! \brief Selects elements from \ref gmx_simd_int32_t where boolean is true, otherwise 0.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_blendzero_di,
 * otherwise \ref gmx_simd_blendzero_fi.
 *
 * \copydetails gmx_simd_blendzero_fi
 */
#    define gmx_simd_blendzero_i             gmx_simd_blendzero_fi

/*! \brief Selects elements from \ref gmx_simd_int32_t where boolean is false, otherwise 0.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_blendnotzero_di,
 * otherwise \ref gmx_simd_blendnotzero_fi.
 *
 * \copydetails gmx_simd_blendnotzero_fi
 */
#    define gmx_simd_blendnotzero_i          gmx_simd_blendnotzero_fi

/*! \brief Selects from 2nd int SIMD arg where boolean is true, otherwise 1st arg.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_blendv_di,
 * otherwise \ref gmx_simd_blendv_fi.
 *
 * \copydetails gmx_simd_blendv_fi
 */
#    define gmx_simd_blendv_i                gmx_simd_blendv_fi

/*! \}
 *  \name SIMD conversion operations
 *
 *  These instructions are available when both types involved in the conversion
 *  are defined, e.g. \ref GMX_SIMD_HAVE_REAL and \ref GMX_SIMD_HAVE_INT32
 *  for real-to-integer conversion.
 *  \{
 */

/*! \brief Convert gmx_simd_real_t to gmx_simd_int32_t, round to nearest integer.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cvt_d2i,
 * otherwise \ref gmx_simd_cvt_f2i.
 *
 * \copydetails gmx_simd_cvt_f2i
 */
#    define gmx_simd_cvt_r2i                 gmx_simd_cvt_f2i

/*! \brief Convert gmx_simd_real_t to gmx_simd_int32_t, truncate towards zero
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cvtt_d2i,
 * otherwise \ref gmx_simd_cvtt_f2i.
 *
 * \copydetails gmx_simd_cvtt_f2i
 */
#    define gmx_simd_cvtt_r2i                gmx_simd_cvtt_f2i

/*! \brief Convert gmx_simd_int32_t to gmx_simd_real_t
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cvt_i2d,
 * otherwise \ref gmx_simd_cvt_i2f.
 *
 * \copydetails gmx_simd_cvt_i2f
 */
#    define gmx_simd_cvt_i2r                 gmx_simd_cvt_i2f

/*! \brief Convert from gmx_simd_bool_t to gmx_simd_ibool_t
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cvt_db2dib,
 * otherwise \ref gmx_simd_cvt_fb2fib.
 *
 * \copydetails gmx_simd_cvt_fb2fib
 */
#    define gmx_simd_cvt_b2ib                gmx_simd_cvt_fb2fib

/*! \brief Convert from gmx_simd_ibool_t to gmx_simd_bool_t
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_cvt_dib2db,
 * otherwise \ref gmx_simd_cvt_fib2fb.
 *
 * \copydetails gmx_simd_cvt_fib2fb
 */
#    define gmx_simd_cvt_ib2b                gmx_simd_cvt_fib2fb


/*! \}
 *  \name SIMD memory alignment operations
 *  \{
 */

/*! \brief Align real memory for SIMD usage.
 *
 * This routine will only align memory if \ref GMX_SIMD_HAVE_REAL is defined.
 * Otherwise the original pointer will be returned.
 *
 * Start by allocating an extra \ref GMX_SIMD_REAL_WIDTH float elements of memory,
 * and then call this function. The returned pointer will be greater or equal
 * to the one you provided, and point to an address inside your provided memory
 * that is aligned to the SIMD width.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_align_d,
 * otherwise \ref gmx_simd_align_f. For detailed documentation, see the
 * precision-specific implementation routines.
 */
#    define gmx_simd_align_r                 gmx_simd_align_f

/*! \brief Align integer memory for SIMD usage.
 *
 * This routine will only align memory if \ref GMX_SIMD_HAVE_INT32 is defined.
 * Otherwise the original pointer will be returned.
 *
 * Start by allocating an extra \ref GMX_SIMD_INT32_WIDTH elements of memory,
 * and then call this function. The returned pointer will be greater or equal
 * to the one you provided, and point to an address inside your provided memory
 * that is aligned to the SIMD width.
 *
 * If GMX_DOUBLE is defined, this will be aliased to \ref gmx_simd_align_di,
 * otherwise \ref gmx_simd_align_fi. For detailed documentation, see the
 * precision-specific implementation routines.
 */
#    define gmx_simd_align_i                 gmx_simd_align_fi

/*! \} */

/*! \name SIMD4 - constant width-four SIMD datatypes
 *
 * These operations are only meant to be used for a few coordinate
 * manipulation and grid interpolation routines, so we only support a subset
 * of operations for SIMD4. To avoid repeating all the documentation from
 * the generic width SIMD routines, we only provide brief documentation for
 * these operations. Follow the link to the implementation documentation or the
 * reference to the corresponding generic SIMD routine. The format will be
 * exactly the same, but they have SIMD replaced with SIMD4.
 *  \{
 */

/*! \brief SIMD real datatype guaranteed to be 4 elements wide, if available.
 *
 * All the SIMD4 datatypes and operations behave like their counterparts for
 * the generic SIMD implementation, but they might be implemented with different
 * registers, or not supported at all. It is important that you check the
 * define \ref GMX_SIMD4_HAVE_REAL before using it.
 *
 * Just as the normal SIMD operations, all SIMD4 types and routines will
 * be aliased to either single or double precision ones based on whether
 * GMX_DOUBLE is defined.
 *
 * \note There is no support for integer or math operations in SIMD4.
 */
#    define gmx_simd4_real_t                 gmx_simd4_float_t

/*! \brief Boolean for \ref gmx_simd4_real_t comparision/selection */
#    define gmx_simd4_bool_t                 gmx_simd4_fbool_t

/*! \brief Load aligned data to gmx_simd4_real_t.
 *
 * \copydetails gmx_simd4_load_f
 */
#    define gmx_simd4_load_r                 gmx_simd4_load_f

/*! \brief Load single element to gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_load1_f
 */
#    define gmx_simd4_load1_r                gmx_simd4_load1_f

/*! \brief Set gmx_simd4_real_t from scalar value
 *
 * \copydetails gmx_simd4_set1_f
 */
#    define gmx_simd4_set1_r                 gmx_simd4_set1_f

/*! \brief store aligned data from gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_store_f
 */
#    define gmx_simd4_store_r                gmx_simd4_store_f

/*! \brief Load unaligned data to gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_loadu_f
 */
#    define gmx_simd4_loadu_r                gmx_simd4_loadu_f

/*! \brief Store unaligned data from gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_storeu_f
 */
#    define gmx_simd4_storeu_r               gmx_simd4_storeu_f

/*! \brief Set all elements in gmx_simd4_real_t to 0.0
 *
 * \copydetails gmx_simd4_setzero_f
 */
#    define gmx_simd4_setzero_r              gmx_simd4_setzero_f

/*! \brief Bitwise and for two gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_and_f
 */
#    define gmx_simd4_and_r                  gmx_simd4_and_f

/*! \brief Bitwise and-not for two gmx_simd4_real_t. 1st arg is complemented.
 *
 * \copydetails gmx_simd4_andnot_f
 */
#    define gmx_simd4_andnot_r               gmx_simd4_andnot_f

/*! \brief Bitwise or for two gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_or_f
 */
#    define gmx_simd4_or_r                   gmx_simd4_or_f

/*! \brief Bitwise xor for two gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_xor_f
 */
#    define gmx_simd4_xor_r                  gmx_simd4_xor_f

/*! \brief a+b for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_add_f
 */
#    define gmx_simd4_add_r                  gmx_simd4_add_f

/*! \brief a-b for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_sub_f
 */
#    define gmx_simd4_sub_r                  gmx_simd4_sub_f

/*! \brief a*b for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_mul_f
 */
#    define gmx_simd4_mul_r                  gmx_simd4_mul_f

/*! \brief a*b+c for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_fmadd_f
 */
#    define gmx_simd4_fmadd_r                gmx_simd4_fmadd_f

/*! \brief a*b-c for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_fmsub_f
 */
#    define gmx_simd4_fmsub_r                gmx_simd4_fmsub_f

/*! \brief -a*b+c for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_fnmadd_f
 */
#    define gmx_simd4_fnmadd_r               gmx_simd4_fnmadd_f

/*! \brief -a*b-c for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_fnmsub_f
 */
#    define gmx_simd4_fnmsub_r               gmx_simd4_fnmsub_f

/*! \brief 1/sqrt(x) approximate lookup for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_rsqrt_f
 */
#    define gmx_simd4_rsqrt_r                gmx_simd4_rsqrt_f

/*! \brief fabs(x) for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_fabs_f
 */
#    define gmx_simd4_fabs_r                 gmx_simd4_fabs_f

/*! \brief Change sign (-x) for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_fneg_f
 */
#    define gmx_simd4_fneg_r                 gmx_simd4_fneg_f

/*! \brief Select maximum of each pair of elements from args for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_max_f
 */
#    define gmx_simd4_max_r                  gmx_simd4_max_f

/*! \brief Select minimum of each pair of elements from args for \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_min_f
 */
#    define gmx_simd4_min_r                  gmx_simd4_min_f

/*! \brief Round \ref gmx_simd4_real_t to nearest integer, return \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_round_f
 */
#    define gmx_simd4_round_r                gmx_simd4_round_f

/*! \brief Truncate \ref gmx_simd4_real_t towards zero, return \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_trunc_f
 */
#    define gmx_simd4_trunc_r                gmx_simd4_trunc_f

/*! \brief Scalar product of first three elements of two \ref gmx_simd4_real_t *
 *
 * \copydetails gmx_simd4_dotproduct3_f
 */
#    define gmx_simd4_dotproduct3_r          gmx_simd4_dotproduct3_f

/*! \brief Return booleans whether a==b for each element two \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_cmpeq_f
 */
#    define gmx_simd4_cmpeq_r                gmx_simd4_cmpeq_f
/*! \brief Return booleans whether a<b for each element two \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_cmplt_f
 */
#    define gmx_simd4_cmplt_r                gmx_simd4_cmplt_f
/*! \brief Return booleans whether a<=b for each element two \ref gmx_simd4_real_t
 *
 * \copydetails gmx_simd4_cmple_f
 */
#    define gmx_simd4_cmple_r                gmx_simd4_cmple_f

/*! \brief Logical and for two \ref gmx_simd4_bool_t
 *
 * \copydetails gmx_simd4_and_fb
 */
#    define gmx_simd4_and_b                  gmx_simd4_and_fb
/*! \brief Logical or for two \ref gmx_simd4_bool_t
 *
 * \copydetails gmx_simd4_or_fb
 */
#    define gmx_simd4_or_b                   gmx_simd4_or_fb

/*! \brief Return nonzero if any element in \ref gmx_simd4_bool_t is true, otherwise 0
 *
 * \copydetails gmx_simd4_anytrue_fb
 */
#    define gmx_simd4_anytrue_b              gmx_simd4_anytrue_fb

/*! \brief Selects from 2nd real SIMD4 arg where boolean is true, otherwise 1st arg
 *
 * \copydetails gmx_simd4_blendzero_f
 */
#    define gmx_simd4_blendzero_r            gmx_simd4_blendzero_f

/*! \brief Selects from 2nd real SIMD4 arg where boolean is false, otherwise 1st arg
 *
 * \copydetails gmx_simd4_blendnotzero_f
 */
#    define gmx_simd4_blendnotzero_r            gmx_simd4_blendnotzero_f

/*! \brief Selects from 2nd real SIMD4 arg where boolean is true, otherwise 1st arg
 *
 * \copydetails gmx_simd4_blendv_f
 */
#    define gmx_simd4_blendv_r               gmx_simd4_blendv_f

/*! \brief Return sum of all elements in SIMD4 floating-point variable.
 *
 * \copydetails gmx_simd4_reduce_f
 */
#    define gmx_simd4_reduce_r               gmx_simd4_reduce_f

/*! \brief Align real memory for SIMD4 usage.
 *
 * \copydetails gmx_simd4_align_f
 */
#    define gmx_simd4_align_r                gmx_simd4_align_f

/*! \} */

/*! \name SIMD predefined macros to describe high-level capabilities
 *  \{
 */

#    if (defined GMX_SIMD_HAVE_FLOAT) || (defined DOXYGEN)
/*! \brief Defined if gmx_simd_real_t is available.
 *
 *  if GMX_DOUBLE is defined, this will be aliased to
 *  \ref GMX_SIMD_HAVE_DOUBLE, otherwise GMX_SIMD_HAVE_FLOAT.
 */
#        define GMX_SIMD_HAVE_REAL
/*! \brief Width of gmx_simd_real_t.
 *
 *  if GMX_DOUBLE is defined, this will be aliased to
 *  \ref GMX_SIMD_DOUBLE_WIDTH, otherwise GMX_SIMD_FLOAT_WIDTH.
 */
#        define GMX_SIMD_REAL_WIDTH          GMX_SIMD_FLOAT_WIDTH
#    endif
#    if (defined GMX_SIMD_HAVE_FINT32) || (defined DOXYGEN)
/*! \brief Defined if gmx_simd_int32_t is available.
 *
 *  if GMX_DOUBLE is defined, this will be aliased to
 *  \ref GMX_SIMD_HAVE_DINT32, otherwise GMX_SIMD_HAVE_FINT32.
 */
#        define GMX_SIMD_HAVE_INT32
/*! \brief Width of gmx_simd_int32_t.
 *
 *  if GMX_DOUBLE is defined, this will be aliased to
 *  \ref GMX_SIMD_DINT32_WIDTH, otherwise GMX_SIMD_FINT32_WIDTH.
 */
#        define GMX_SIMD_INT32_WIDTH         GMX_SIMD_FINT32_WIDTH
#    endif
#    if (defined GMX_SIMD_HAVE_FINT32_EXTRACT) || (defined DOXYGEN)
/*! \brief Defined if gmx_simd_extract_i() is available.
 *
 *  if GMX_DOUBLE is defined, this will be aliased to
 *  \ref GMX_SIMD_HAVE_DINT32_EXTRACT, otherwise GMX_SIMD_HAVE_FINT32_EXTRACT.
 */
#        define GMX_SIMD_HAVE_INT32_EXTRACT
#    endif
#    if (defined GMX_SIMD_HAVE_FINT32_LOGICAL) || (defined DOXYGEN)
/*! \brief Defined if logical ops are supported on gmx_simd_int32_t.
 *
 *  if GMX_DOUBLE is defined, this will be aliased to
 *  \ref GMX_SIMD_HAVE_DINT32_LOGICAL, otherwise GMX_SIMD_HAVE_FINT32_LOGICAL.
 */
#        define GMX_SIMD_HAVE_INT32_LOGICAL
#    endif
#    if (defined GMX_SIMD_HAVE_FINT32_ARITHMETICS) || (defined DOXYGEN)
/*! \brief Defined if arithmetic ops are supported on gmx_simd_int32_t.
 *
 *  if GMX_DOUBLE is defined, this will be aliased to
 *  \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS, otherwise GMX_SIMD_HAVE_FINT32_ARITHMETICS.
 */
#        define GMX_SIMD_HAVE_INT32_ARITHMETICS
#    endif
#    if (defined GMX_SIMD4_HAVE_FLOAT) || (defined DOXYGEN)
/*! \brief Defined if gmx_simd4_real_t is available.
 *
 *  if GMX_DOUBLE is defined, this will be aliased to
 *  \ref GMX_SIMD4_HAVE_DOUBLE, otherwise GMX_SIMD4_HAVE_FLOAT.
 */
#        define GMX_SIMD4_HAVE_REAL
#    endif

/*! \} */

#endif /* GMX_DOUBLE */

/*! \} */
/*! \endcond */

#if 0
/* Finally, a hack to cover a possible corner case of using an
   explicit GMX_SIMD_HAVE_FLOAT or GMX_SIMD_HAVE_DOUBLE, rather than
   GMX_SIMD_HAVE_REAL.

   Such code is expected to include simd.h to get those symbols
   defined, but the actual definitions are in the implemention headers
   included by simd.h. check-source.py is not a full preprocessor, so
   it does not see the definitions in the implementation headers as
   belonging to simd.h, thus it cannot check that simd.h is being used
   correctly in the above hypothetical corner case. However, the
   checker also does not parse #if 0, so we can fool the checker into
   thinking that definition occurs here, and that will work well
   enough.

   If there's ever other kinds of SIMD code that might have the same
   problem, we might want to add other variables here.
 */
#    define GMX_SIMD_HAVE_FLOAT
#    define GMX_SIMD_HAVE_DOUBLE

#endif /* 0 */

#endif /* GMX_SIMD_SIMD_H */
