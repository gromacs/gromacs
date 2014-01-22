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

/*! \defgroup module_simd SIMD intrinsics interface (simd)
 *  \ingroup group_utilitymodules
 *
 *  \brief Provides an architecture-independent way of doing SIMD coding.
 *
 *  Start by consulting the overview Doxygen SIMD documentation, and then
 *  the details are documented in simd.h and the reference implementation
 *  impl_reference.h
 *
 *  \author Erik Lindahl <erik.lindahl@scilifelab.se>
 */

/*! \libinternal \file
 *
 * \brief Definitions, capabilities, and wrappers for SIMD module.
 *
 * The macros in this file are intended to be used for writing
 * architecture-independent SIMD intrinsics code.
 * To support a new architecture, adding a new sub-include with macros here
 * should be (nearly) all that is needed.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

#ifndef GMX_SIMD_SIMD_H
#define GMX_SIMD_SIMD_H

#include "gromacs/legacyheaders/types/simple.h"

/* Forward declarations so memory allocation can be used in implementations */
static gmx_inline float *  gmx_simd_align_f(float *p);
static gmx_inline double * gmx_simd_align_d(double *p);
static gmx_inline int *    gmx_simd_align_fi(int *p);
static gmx_inline int *    gmx_simd_align_di(int *p);
static gmx_inline float *  gmx_simd4_align_f(float *p);
static gmx_inline double * gmx_simd4_align_d(double *p);

/*! This define indicates that some sort of SIMD support is present.
 *
 * It is disabled if no architecture, neither reference SIMD, has been selected.
 */
#define GMX_SIMD


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
#else
/* Turn off the GMX_SIMD flag if we do not even have reference support */
#    undef GMX_SIMD
#endif

/*! Just for convenience to make it clear when we refer to SIMD4 width */
#define GMX_SIMD4_WIDTH    4

/*! Align a float pointer for usage with SIMD instructions.
 *
 * Start by allocating an extra GMX_SIMD_FLOAT_WIDTH float elements of memory,
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

/*! Align a double pointer for usage with SIMD instructions.
 *
 * Start by allocating an extra GMX_SIMD_DOUBLE_WIDTH double elements of memory,
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

/*! Align a (float) integer pointer for usage with SIMD instructions.
 *
 * This routine provides aligned memory for usage with gmx_simd_fint32_t. You
 * should have allocated an extra GMX_SIMD_FINT32_WIDTH*sizeof(int) bytes.
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

/*! Align a (double) integer pointer for usage with SIMD instructions.
 *
 * This routine provides aligned memory for usage with gmx_simd_dint32_t. You
 * should have allocated an extra GMX_SIMD_DINT32_WIDTH*sizeof(int) bytes.
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

/*! Align a float pointer for usage with SIMD4 instructions.
 *
 * This routine provides aligned memory for usage with gmx_simd4_float_t.
 * should have allocated an extra GMX_SIMD4_WIDTH*sizeof(float) bytes.
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

/*! Align a double pointer for usage with SIMD4 instructions.
 *
 * This routine provides aligned memory for usage with gmx_simd4_double_t.
 * should have allocated an extra GMX_SIMD4_WIDTH*sizeof(double) bytes.
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

/* Single precision floating-point */
/*! Real SIMD datatype */
#    define gmx_simd_real_t                  gmx_simd_float_t
/*! Load aligned real data, width is GMX_SIMD_REAL_WIDTH */
#    define gmx_simd_load_r                  gmx_simd_load_f
/*! Load a single real to a SIMD register */
#    define gmx_simd_load1_r                 gmx_simd_load1_f
/*! Set SIMD real register from FP scalar */
#    define gmx_simd_set1_r                  gmx_simd_set1_f
/*! Store aligned data (width GMX_SIMD_REAL_WIDTH) from real SIMD */
#    define gmx_simd_store_r                 gmx_simd_store_f
/*! Load unaligned data to SIMD real, width GMX_SIMD_REAL_WIDTH */
#    define gmx_simd_loadu_r                 gmx_simd_loadu_f
/*! Store unaligned data from SIMD real, width GMX_SIMD_REAL_WIDTH */
#    define gmx_simd_storeu_r                gmx_simd_storeu_f
/*! Set all real SIMD elements to 0.0 */
#    define gmx_simd_setzero_r               gmx_simd_setzero_f
/*! SIMD a+b */
#    define gmx_simd_add_r                   gmx_simd_add_f
/*! Simd a-b */
#    define gmx_simd_sub_r                   gmx_simd_sub_f
/*! Simd a*b */
#    define gmx_simd_mul_r                   gmx_simd_mul_f
/*! Simd a*b+c. You might improve performance by assigning the result to c. */
#    define gmx_simd_fmadd_r                 gmx_simd_fmadd_f
/*! Simd a*b-c. You might improve performance by assigning the result to c. */
#    define gmx_simd_fmsub_r                 gmx_simd_fmsub_f
/*! Simd -a*b+c. You might improve performance by assigning the result to c. */
#    define gmx_simd_fnmadd_r                gmx_simd_fnmadd_f
/*! Simd -a*b-c. You might improve performance by assigning the result to c. */
#    define gmx_simd_fnmsub_r                gmx_simd_fnmsub_f
/*! Bitwise and on SIMD real */
#    define gmx_simd_and_r                   gmx_simd_and_f
/*! Bitwise and-not (first argument is bitwise complemented) on SIMD vars */
#    define gmx_simd_andnot_r                gmx_simd_andnot_f
/*! Bitwise or for SIMD real */
#    define gmx_simd_or_r                    gmx_simd_or_f
/*! Bitwise xor for SIMD real */
#    define gmx_simd_xor_r                   gmx_simd_xor_f
/*! SIMD 1/sqrt(x) approximation lookup. */
#    define gmx_simd_rsqrt_r                 gmx_simd_rsqrt_f
/*! SIMD 1/x approximation lookup. */
#    define gmx_simd_rcp_r                   gmx_simd_rcp_f
/*! SIMD fabs(x) */
#    define gmx_simd_fabs_r                  gmx_simd_fabs_f
/*! SIMD -x */
#    define gmx_simd_fneg_r                  gmx_simd_fneg_f
/*! SIMD max of each real element */
#    define gmx_simd_max_r                   gmx_simd_max_f
/*! SIMD min of each real element */
#    define gmx_simd_min_r                   gmx_simd_min_f
/*! Round SIMD real to nearest integer, return SIMD real */
#    define gmx_simd_round_r                 gmx_simd_round_f
/*! Truncate SIMD real towards zero, return SIMD real */
#    define gmx_simd_trunc_r                 gmx_simd_trunc_f
/*! Returns SIMD fraction, i.e. x-trunc(x) */
#    define gmx_simd_fraction_r              gmx_simd_fraction_f
/*! Return the exponent of a SIMD real as a SIMD real */
#    define gmx_simd_get_exponent_r          gmx_simd_get_exponent_f
/*! Return the mantissa of a SIMD real as a SIMD real */
#    define gmx_simd_get_mantissa_r          gmx_simd_get_mantissa_f
/*! Set the exponent of a SIMD real from a SIMD real */
#    define gmx_simd_set_exponent_r          gmx_simd_set_exponent_f

/* Float integer and conversions */

/*! SIMD 32-bit integer datatype */

/*! 32-bit integer SIMD type corresponding to gmx_simd_real_t */
#    define gmx_simd_int32_t                 gmx_simd_fint32_t
/*! Load aligned memory (width GMX_SIMD_INT32_WIDTH) to gmx_simd_int32_t */
#    define gmx_simd_load_i                  gmx_simd_load_fi
/*! Set gmx_simd_int32_t from scalar integer */
#    define gmx_simd_set1_i                  gmx_simd_set1_fi
/*! Store aligned memory (width GMX_SIMD_INT32_WIDTH) from gmx_simd_int32_t */
#    define gmx_simd_store_i                 gmx_simd_store_fi
/*! Load unaligned memory (width GMX_SIMD_INT32_WIDTH) to gmx_simd_int32_t */
#    define gmx_simd_loadu_i                 gmx_simd_loadu_fi
/*! Store unaligned memory (width GMX_SIMD_INT32_WIDTH) from gmx_simd_int32_t */
#    define gmx_simd_storeu_i                gmx_simd_storeu_fi
/*! Set all elements in gmx_simd_int32_t to 0 */
#    define gmx_simd_setzero_i               gmx_simd_setzero_fi
/*! Convert gmx_simd_real_t to gmx_simd_int32_t, round to nearest integer */
#    define gmx_simd_cvt_r2i                 gmx_simd_cvt_f2i
/*! Convert gmx_simd_real_t to gmx_simd_int32_t, truncate towards zero */
#    define gmx_simd_cvtt_r2i                gmx_simd_cvtt_f2i
/*! Convert gmx_simd_int32_t to gmx_simd_real_t */
#    define gmx_simd_cvt_i2r                 gmx_simd_cvt_i2f
/*! Extract scalar integer from gmx_simd_int32_t */
#    define gmx_simd_extract_i               gmx_simd_extract_fi
/*! Shift each element in gmx_simd_int32_t left by immediate */
#    define gmx_simd_slli_i                  gmx_simd_slli_fi
/*! Shift each element in gmx_simd_int32_t right by immediate */
#    define gmx_simd_srli_i                  gmx_simd_srli_fi
/*! Bitwise and for two gmx_simd_int32_t */
#    define gmx_simd_and_i                   gmx_simd_and_fi
/*! Bitwise and-not for two gmx_simd_int32_t. First arg is complemented */
#    define gmx_simd_andnot_i                gmx_simd_andnot_fi
/*! Bitwise or for two gmx_simd_int32_t */
#    define gmx_simd_or_i                    gmx_simd_or_fi
/*! Bitwise xord for two gmx_simd_int32_t */
#    define gmx_simd_xor_i                   gmx_simd_xor_fi
/*! a+b for gmx_simd_int32_t */
#    define gmx_simd_add_i                   gmx_simd_add_fi
/*! a-b for gmx_simd_int32_t */
#    define gmx_simd_sub_i                   gmx_simd_sub_fi
/*! a*b for gmx_simd_int32_t */
#    define gmx_simd_mul_i                   gmx_simd_mul_fi

/* Float booleans and selection */

/*! Boolean SIMD corresponding to gmx_simd_real_t */
#    define gmx_simd_bool_t                  gmx_simd_fbool_t
/*! Returns boolean describing whether a==b, for gmx_simd_real_t */
#    define gmx_simd_cmpeq_r                 gmx_simd_cmpeq_f
/*! Returns boolean describing whether a<b, for gmx_simd_real_t */
#    define gmx_simd_cmplt_r                 gmx_simd_cmplt_f
/*! Returns boolean describing whether a<=b, for gmx_simd_real_t */
#    define gmx_simd_cmple_r                 gmx_simd_cmple_f
/*! Logical and between two gmx_simd_bool_t */
#    define gmx_simd_and_b                   gmx_simd_and_fb
/*! Logical or between two gmx_simd_bool_t */
#    define gmx_simd_or_b                    gmx_simd_or_fb
/*! Return nonzero if any element in gmx_simd_bool_t is true, otherwise 0 */
#    define gmx_simd_anytrue_b               gmx_simd_anytrue_fb
/*! Selects elements from 1st real SIMD arg where boolean is true, otherwise 0 */
#    define gmx_simd_blendzero_r             gmx_simd_blendzero_f
/*! Selects from 2nd real SIMD arg where boolean is true, otherwise 1st arg */
#    define gmx_simd_blendv_r                gmx_simd_blendv_f
/*! Boolean simd corresponding to gmx_simd_int32_t */
#    define gmx_simd_ibool_t                 gmx_simd_fibool_t
/*! Returns boolean describing whether a==b, for gmx_simd_int32_t */
#    define gmx_simd_cmpeq_i                 gmx_simd_cmpeq_fi
/*! Returns boolean describing whether a<b, for gmx_simd_int32_t */
#    define gmx_simd_cmplt_i                 gmx_simd_cmplt_fi
/*! Logical and between two gmx_simd_ibool_t */
#    define gmx_simd_and_ib                  gmx_simd_and_fib
/*! Logical or between two gmx_simd_ibool_t */
#    define gmx_simd_or_ib                   gmx_simd_or_fib
/*! Return nonzero if any element in gmx_simd_ibool_t is true, otherwise 0 */
#    define gmx_simd_anytrue_ib              gmx_simd_anytrue_fib
/*! Selects from 2nd int SIMD arg where boolean is true, otherwise 1st arg */
#    define gmx_simd_blendzero_i             gmx_simd_blendzero_fi
/*! Selects elements from 1st int SIMD arg where boolean is true, otherwise 0 */
#    define gmx_simd_blendv_i                gmx_simd_blendv_fi

/* Conversions between integer and single floating-point booleans */

/*! Convert from gmx_simd_bool_t to gmx_simd_ibool_t */
#    define gmx_simd_cvt_b2ib                gmx_simd_cvt_fb2fib
/*! Convert from gmx_simd_ibool_t to gmx_simd_bool_t */
#    define gmx_simd_cvt_ib2b                gmx_simd_cvt_fib2fb

/* SIMD4 float fp - we only support a subset of instructions for SIMD4 */

/*! SIMD real datatype guaranteed to be 4 elements wide, if available. */
#    define gmx_simd4_real_t                 gmx_simd4_float_t
/*! Load aligned data to gmx_simd4_real_t */
#    define gmx_simd4_load_r                 gmx_simd4_load_f
/*! Load single element to gmx_simd4_real_t */
#    define gmx_simd4_load1_r                gmx_simd4_load1_f
/*! Set gmx_simd4_real_t from scalar value */
#    define gmx_simd4_set1_r                 gmx_simd4_set1_f
/*! store aligned data from gmx_simd4_real_t */
#    define gmx_simd4_store_r                gmx_simd4_store_f
/*! Load unaligned data to gmx_simd4_real_t */
#    define gmx_simd4_loadu_r                gmx_simd4_loadu_f
/*! Store unaligned data from gmx_simd4_real_t */
#    define gmx_simd4_storeu_r               gmx_simd4_storeu_f
/*! Set all elements in gmx_simd4_real_t to 0.0 */
#    define gmx_simd4_setzero_r              gmx_simd4_setzero_f
/*! a+b for gmx_simd4_real_t */
#    define gmx_simd4_add_r                  gmx_simd4_add_f
/*! a-b for gmx_simd4_real_t */
#    define gmx_simd4_sub_r                  gmx_simd4_sub_f
/*! a*b for gmx_simd4_real_t */
#    define gmx_simd4_mul_r                  gmx_simd4_mul_f
/*! a*b+c for gmx_simd4_real_t */
#    define gmx_simd4_fmadd_r                gmx_simd4_fmadd_f
/*! a*b-c for gmx_simd4_real_t */
#    define gmx_simd4_fmsub_r                gmx_simd4_fmsub_f
/*! -a*b+c for gmx_simd4_real_t */
#    define gmx_simd4_fnmadd_r               gmx_simd4_fnmadd_f
/*! -a*b-c for gmx_simd4_real_t */
#    define gmx_simd4_fnmsub_r               gmx_simd4_fnmsub_f
/*! Bitwise and for two gmx_simd4_real_t */
#    define gmx_simd4_and_r                  gmx_simd4_and_f
/*! Bitwise and-not for two gmx_simd4_real_t. 1st arg is complemented. */
#    define gmx_simd4_andnot_r               gmx_simd4_andnot_f
/*! Bitwise or for two gmx_simd4_real_t */
#    define gmx_simd4_or_r                   gmx_simd4_or_f
/*! Bitwise xor for two gmx_simd4_real_t */
#    define gmx_simd4_xor_r                  gmx_simd4_xor_f
/*! 1/sqrt(x) approximate lookup for gmx_simd4_real_t */
#    define gmx_simd4_rsqrt_r                gmx_simd4_rsqrt_f
/*! fabs(x) gmx_simd4_real_t */
#    define gmx_simd4_fabs_r                 gmx_simd4_fabs_f
/*! -x gmx_simd4_real_t */
#    define gmx_simd4_fneg_r                 gmx_simd4_fneg_f
/*! Select maximum of each pair of elements from args for gmx_simd4_real_t */
#    define gmx_simd4_max_r                  gmx_simd4_max_f
/*! Select minimum of each pair of elements from args for gmx_simd4_real_t */
#    define gmx_simd4_min_r                  gmx_simd4_min_f
/*! Round gmx_simd4_real_t to nearest integer, return gmx_simd4_real_t */
#    define gmx_simd4_round_r                gmx_simd4_round_f
/*! Truncate gmx_simd4_real_t towards zero, return gmx_simd4_real_t */
#    define gmx_simd4_trunc_r                gmx_simd4_trunc_f
/*! Scalar product of first three elements of two gmx_simd4_real_t */
#    define gmx_simd4_dotproduct3_r          gmx_simd4_dotproduct3_f
/*! Boolean for gmx_simd4_real_t comparision/selection */
#    define gmx_simd4_bool_t                 gmx_simd4_fbool_t
/*! Return booleans whether a==b for each element two gmx_simd4_real_t */
#    define gmx_simd4_cmpeq_r                gmx_simd4_cmpeq_f
/*! Return booleans whether a<b for each element two gmx_simd4_real_t */
#    define gmx_simd4_cmplt_r                gmx_simd4_cmplt_f
/*! Return booleans whether a<=b for each element two gmx_simd4_real_t */
#    define gmx_simd4_cmple_r                gmx_simd4_cmple_f
/*! Logical and for two gmx_simd4_bool_t */
#    define gmx_simd4_and_b                  gmx_simd4_and_fb
/*! Logical or for two gmx_simd4_bool_t */
#    define gmx_simd4_or_b                   gmx_simd4_or_fb
/*! Return nonzero if any element in gmx_simd4_bool_t is true, otherwise 0 */
#    define gmx_simd4_anytrue_b              gmx_simd4_anytrue_fb
/*! Selects from 2nd real SIMD4 arg where boolean is true, otherwise 1st arg */
#    define gmx_simd4_blendzero_r            gmx_simd4_blendzero_f
/*! Selects from 2nd real SIMD4 arg where boolean is true, otherwise 1st arg */
#    define gmx_simd4_blendv_r               gmx_simd4_blendv_f

/* Memory allocation */


/*! Align real memory for SIMD usage. See gmx_simd_align_f for details */
#    define gmx_simd_align_r                 gmx_simd_align_f
/*! Align integer memory for SIMD usage. See gmx_simd_align_fi for details */
#    define gmx_simd_align_i                 gmx_simd_align_fi
/*! Align real memory for SIMD4 usage. See gmx_simd4_align_f for details */
#    define gmx_simd4_align_r                gmx_simd4_align_f

#    if (defined GMX_SIMD_HAVE_FLOAT) || (defined DOXYGEN)
/*! Defined if gmx_simd_real_t is available */
#        define GMX_SIMD_HAVE_REAL
/*! Width of gmx_simd_real_t */
#        define GMX_SIMD_REAL_WIDTH          GMX_SIMD_FLOAT_WIDTH
#    endif
#    if (defined GMX_SIMD_HAVE_FINT32) || (defined DOXYGEN)
/*! Defined if gmx_simd_int32_t is available */
#        define GMX_SIMD_HAVE_INT32
/*! Width of gmx_simd_int32_t */
#        define GMX_SIMD_INT32_WIDTH         GMX_SIMD_FINT32_WIDTH
#    endif
#    if (defined GMX_SIMD_HAVE_FINT32_EXTRACT) || (defined DOXYGEN)
/*! Defined if gmx_simd_extract_i() is available */
#        define GMX_SIMD_HAVE_INT32_EXTRACT
#    endif
#    if (defined GMX_SIMD_HAVE_FINT32_LOGICAL) || (defined DOXYGEN)
/*! Defined if logical ops are supported on gmx_simd_int32_t */
#        define GMX_SIMD_HAVE_INT32_LOGICAL
#    endif
#    if (defined GMX_SIMD_HAVE_FINT32_ARITHMETICS) || (defined DOXYGEN)
/*! Defined if arithmetic ops are supported on gmx_simd_int32_t */
#        define GMX_SIMD_HAVE_INT32_ARITHMETICS
#    endif
#    if (defined GMX_SIMD4_HAVE_FLOAT) || (defined DOXYGEN)
/*! Defined if gmx_simd4_real_t is available */
#        define GMX_SIMD4_HAVE_REAL
#    endif

#endif /* GMX_DOUBLE */

#include "simd_math.h"

#endif /* GMX_SIMD_SIMD_H */
