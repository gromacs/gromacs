/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_INTEL_MIC_H
#define GMX_SIMD_IMPL_INTEL_MIC_H

#ifndef GMX_SIMD_SIMD_H
# error "Never include architecture-specific simd headers directly; use simd.h."
#endif

#include <math.h>
#include <immintrin.h>

/* Intel Xeon Phi, or
 * the-artist-formerly-known-as-Knight's-corner, or
 * the-artist-formerly-formerly-known-as-MIC, or
 * the artist formerly-formerly-formerly-known-as-Larrabee
 * 512-bit SIMD instruction wrappers.
 */

/* Capability definitions for Xeon Phi SIMD */
#define GMX_SIMD_HAVE_FLOAT
#define GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_SIMD_HARDWARE
#define GMX_SIMD_HAVE_LOADU
#define GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#define  GMX_SIMD_HAVE_FINT32_EXTRACT
#define GMX_SIMD_HAVE_FINT32_LOGICAL
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS
#define GMX_SIMD_HAVE_DINT32
#define  GMX_SIMD_HAVE_DINT32_EXTRACT
#define GMX_SIMD_HAVE_DINT32_LOGICAL
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS
#undef  GMX_SIMD4_HAVE_FLOAT
#undef  GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH        16
#define GMX_SIMD_DOUBLE_WIDTH        8
#define GMX_SIMD_FINT32_WIDTH       16
#define GMX_SIMD_DINT32_WIDTH        8
#define GMX_SIMD_RSQRT_BITS         23
#define GMX_SIMD_RCP_BITS           23

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t           __m512
#define gmx_simd_load_f            _mm512_load_ps
#define gmx_simd_load1_f(m)        _mm512_extload_ps(m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE)
#define gmx_simd_set1_f            _mm512_set1_ps
#define gmx_simd_store_f           _mm512_store_ps
#define gmx_simd_loadu_f           gmx_simd_loadu_f_mic
#define gmx_simd_storeu_f          gmx_simd_storeu_f_mic
#define gmx_simd_setzero_f         _mm512_setzero_ps
#define gmx_simd_add_f             _mm512_add_ps
#define gmx_simd_sub_f             _mm512_sub_ps
#define gmx_simd_mul_f             _mm512_mul_ps
#define gmx_simd_fmadd_f           _mm512_fmadd_ps
#define gmx_simd_fmsub_f           _mm512_fmsub_ps
#define gmx_simd_fnmadd_f          _mm512_fnmadd_ps
#define gmx_simd_fnmsub_f          _mm512_fnmsub_ps
#define gmx_simd_and_f(a, b)        _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd_andnot_f(a, b)     _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd_or_f(a, b)         _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd_xor_f(a, b)        _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd_rsqrt_f           _mm512_rsqrt23_ps
#define gmx_simd_rcp_f             _mm512_rcp23_ps
#define gmx_simd_fabs_f(x)         gmx_simd_andnot_f(_mm512_set1_ps(-0.0), x)
#define gmx_simd_fneg_f(x)         _mm512_addn_ps(x, _mm512_setzero_ps())
#define gmx_simd_max_f             _mm512_gmax_ps
#define gmx_simd_min_f             _mm512_gmin_ps
#define gmx_simd_round_f(x)        _mm512_round_ps(x, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd_trunc_f(x)        _mm512_round_ps(x, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd_fraction_f(x)     _mm512_sub_ps(x, gmx_simd_trunc_f(x))
#define gmx_simd_get_exponent_f(x) _mm512_getexp_ps(x)
#define gmx_simd_get_mantissa_f(x) _mm512_getmant_ps(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define gmx_simd_set_exponent_f(x) gmx_simd_set_exponent_f_mic(x)
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t          __m512i
#define gmx_simd_load_fi           _mm512_load_epi32
#define gmx_simd_set1_fi           _mm512_set1_epi32
#define gmx_simd_store_fi          _mm512_store_epi32
#define gmx_simd_loadu_fi          gmx_simd_loadu_fi_mic
#define gmx_simd_storeu_fi         gmx_simd_storeu_fi_mic
#define gmx_simd_extract_fi        gmx_simd_extract_fi_mic
#define gmx_simd_setzero_fi        _mm512_setzero_epi32
#define gmx_simd_cvt_f2i(a)        _mm512_cvtfxpnt_round_adjustps_epi32(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd_cvtt_f2i(a)       _mm512_cvtfxpnt_round_adjustps_epi32(a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd_cvt_i2f(a)        _mm512_cvtfxpnt_round_adjustepi32_ps(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
/* Integer logical ops on gmx_simd_fint32_t */
#define gmx_simd_slli_fi           _mm512_slli_epi32
#define gmx_simd_srli_fi           _mm512_srli_epi32
#define gmx_simd_and_fi            _mm512_and_epi32
#define gmx_simd_andnot_fi         _mm512_andnot_epi32
#define gmx_simd_or_fi             _mm512_or_epi32
#define gmx_simd_xor_fi            _mm512_xor_epi32
/* Integer arithmetic ops on gmx_simd_fint32_t */
#define gmx_simd_add_fi            _mm512_add_epi32
#define gmx_simd_sub_fi            _mm512_sub_epi32
#define gmx_simd_mul_fi            _mm512_mullo_epi32
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t           __mmask16
#define gmx_simd_cmpeq_f(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmplt_f(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LT_OS)
#define gmx_simd_cmple_f(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LE_OS)
#define gmx_simd_and_fb            _mm512_kand
#define gmx_simd_andnot_fb(a, b)   _mm512_knot(_mm512_kor(a, b))
#define gmx_simd_or_fb             _mm512_kor
#define gmx_simd_anytrue_fb        _mm512_mask2int
#define gmx_simd_blendzero_f(a, sel)    _mm512_mask_mov_ps(_mm512_setzero_ps(), sel, a)
#define gmx_simd_blendnotzero_f(a, sel) _mm512_mask_mov_ps(_mm512_setzero_ps(), _mm512_knot(sel), a)
#define gmx_simd_blendv_f(a, b, sel)    _mm512_mask_blend_ps(sel, a, b)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t          __mmask16
#define gmx_simd_cmpeq_fi(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define gmx_simd_cmplt_fi(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define gmx_simd_and_fib           _mm512_kand
#define gmx_simd_or_fib            _mm512_kor
#define gmx_simd_anytrue_fib       _mm512_mask2int
#define gmx_simd_blendzero_fi(a, sel)    _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define gmx_simd_blendnotzero_fi(a, sel) _mm512_mask_mov_epi32(_mm512_setzero_epi32(), _mm512_knot(sel), a)
#define gmx_simd_blendv_fi(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib(x)     (x)
#define gmx_simd_cvt_fib2fb(x)     (x)

/* MIC provides full single precision of some neat functions: */
/* 1/sqrt(x) and 1/x work fine in simd_math.h, and won't use extra iterations */

#define gmx_simd_exp2_f            gmx_simd_exp2_f_mic
//#define gmx_simd_exp_f             gmx_simd_exp_f_mic //only provides 59ulp but tests require 10ulp
#define gmx_simd_log_f             gmx_simd_log_f_mic

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t          __m512d
#define gmx_simd_load_d            _mm512_load_pd
#define gmx_simd_load1_d(m)        _mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE)
#define gmx_simd_set1_d            _mm512_set1_pd
#define gmx_simd_store_d           _mm512_store_pd
#define gmx_simd_loadu_d           gmx_simd_loadu_d_mic
#define gmx_simd_storeu_d          gmx_simd_storeu_d_mic
#define gmx_simd_setzero_d         _mm512_setzero_pd
#define gmx_simd_add_d             _mm512_add_pd
#define gmx_simd_sub_d             _mm512_sub_pd
#define gmx_simd_mul_d             _mm512_mul_pd
#define gmx_simd_fmadd_d           _mm512_fmadd_pd
#define gmx_simd_fmsub_d           _mm512_fmsub_pd
#define gmx_simd_fnmadd_d          _mm512_fnmadd_pd
#define gmx_simd_fnmsub_d          _mm512_fnmsub_pd
#define gmx_simd_and_d(a, b)       _mm512_castsi512_pd(_mm512_and_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_andnot_d(a, b)    _mm512_castsi512_pd(_mm512_andnot_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_or_d(a, b)        _mm512_castsi512_pd(_mm512_or_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_xor_d(a, b)       _mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_rsqrt_d(x)        _mm512_cvtpslo_pd(_mm512_rsqrt23_ps(_mm512_cvtpd_pslo(x)))
#define gmx_simd_rcp_d(x)          _mm512_cvtpslo_pd(_mm512_rcp23_ps(_mm512_cvtpd_pslo(x)))
#define gmx_simd_fabs_d(x)         gmx_simd_andnot_d(_mm512_set1_pd(-0.0), x)
#define gmx_simd_fneg_d(x)         _mm512_addn_pd(x, _mm512_setzero_pd())
#define gmx_simd_max_d             _mm512_gmax_pd
#define gmx_simd_min_d             _mm512_gmin_pd
#define gmx_simd_round_d(a)        _mm512_roundfxpnt_adjust_pd(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd_trunc_d(a)        _mm512_roundfxpnt_adjust_pd(a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd_fraction_d(x)     _mm512_sub_pd(x, gmx_simd_trunc_d(x))
#define gmx_simd_get_exponent_d(x) _mm512_getexp_pd(x)
#define gmx_simd_get_mantissa_d(x) _mm512_getmant_pd(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define gmx_simd_set_exponent_d(x) gmx_simd_set_exponent_d_mic(x)
/* integer datatype corresponding to float: gmx_simd_fint32_t
   Doesn't use mask other than where required. No side effect expected for operating on the (unused) upper 8.
 */
#define gmx_simd_dint32_t          __m512i
#define gmx_simd_load_di(m)        _mm512_mask_loadunpacklo_epi32(_mm512_undefined_epi32(), mask_loh, m)
#define gmx_simd_set1_di           _mm512_set1_epi32
#define gmx_simd_store_di(m, a)    _mm512_mask_packstorelo_epi32(m, mask_loh, a)
#define gmx_simd_loadu_di          gmx_simd_loadu_di_mic
#define gmx_simd_storeu_di         gmx_simd_storeu_di_mic
#define gmx_simd_extract_di        gmx_simd_extract_di_mic
#define gmx_simd_setzero_di        _mm512_setzero_epi32
#define gmx_simd_cvt_d2i(a)        _mm512_cvtfxpnt_roundpd_epi32lo(a, _MM_FROUND_TO_NEAREST_INT)
#define gmx_simd_cvtt_d2i(a)       _mm512_cvtfxpnt_roundpd_epi32lo(a, _MM_FROUND_TO_ZERO)
#define gmx_simd_cvt_i2d           _mm512_cvtepi32lo_pd
/* Integer logical ops on gmx_simd_fint32_t */
#define gmx_simd_slli_di           _mm512_slli_epi32
#define gmx_simd_srli_di           _mm512_srli_epi32
#define gmx_simd_and_di            _mm512_and_epi32
#define gmx_simd_andnot_di         _mm512_andnot_epi32
#define gmx_simd_or_di             _mm512_or_epi32
#define gmx_simd_xor_di            _mm512_xor_epi32
/* Integer arithmetic ops on gmx_simd_fint32_t */
#define gmx_simd_add_di            _mm512_add_epi32
#define gmx_simd_sub_di            _mm512_sub_epi32
#define gmx_simd_mul_di            _mm512_mullo_epi32
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_dbool_t           __mmask8
#define gmx_simd_cmpeq_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmplt_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LT_OS)
#define gmx_simd_cmple_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LE_OS)
#define gmx_simd_and_db            _mm512_kand
#define gmx_simd_or_db             _mm512_kor
#define gmx_simd_anytrue_db(x)     _mm512_mask2int(x)
#define gmx_simd_blendzero_d(a, sel)    _mm512_mask_mov_pd(_mm512_setzero_pd(), sel, a)
#define gmx_simd_blendnotzero_d(a, sel) _mm512_mask_mov_pd(_mm512_setzero_pd(), _mm512_knot(sel), a)
#define gmx_simd_blendv_d(a, b, sel)    _mm512_mask_blend_pd(sel, a, b)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_dibool_t          __mmask16
#define gmx_simd_cmpeq_di(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define gmx_simd_cmplt_di(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define gmx_simd_and_dib           _mm512_kand
#define gmx_simd_or_dib            _mm512_kor
#define gmx_simd_anytrue_dib(x)    (_mm512_mask2int(x)&0xFF)
#define gmx_simd_blendzero_di(a, sel)    _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define gmx_simd_blendnotzero_di(a, sel) _mm512_mask_mov_epi32(_mm512_setzero_epi32(), _mm512_knot(sel), a)
#define gmx_simd_blendv_di(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between booleans. Double & dint stuff is stored in low bits */
#define gmx_simd_cvt_db2dib(x)     (x)
#define gmx_simd_cvt_dib2db(x)     (x)

/* Float/double conversion */
#define gmx_simd_cvt_f2dd          gmx_simd_cvt_f2dd_mic
#define gmx_simd_cvt_dd2f          gmx_simd_cvt_dd2f_mic

#define PERM_LOW2HIGH _MM_PERM_BABA
#define PERM_HIGH2LOW _MM_PERM_DCDC

#define mask_loh _mm512_int2mask(0x00FF) /* would be better a constant - but can't initialize with a function call. */
#define mask_hih _mm512_int2mask(0xFF00)

/* load store float */
static gmx_inline __m512
gmx_simd_loadu_f_mic(float *m)
{
    return _mm512_loadunpackhi_ps(_mm512_loadunpacklo_ps(_mm512_undefined_ps(), m), m+16);
}

static gmx_inline void
gmx_simd_storeu_f_mic(float* m, __m512 s)
{
    _mm512_packstorelo_ps(m, s);
    _mm512_packstorehi_ps(m+16, s);
}

/* load store fint32 */
static gmx_inline __m512i
gmx_simd_loadu_fi_mic(gmx_int32_t *m)
{
    return _mm512_loadunpackhi_epi32(_mm512_loadunpacklo_epi32(_mm512_undefined_epi32(), m), m+16);
}

static gmx_inline void
gmx_simd_storeu_fi_mic(gmx_int32_t* m, __m512i s)
{
    _mm512_packstorelo_epi32(m, s);
    _mm512_packstorehi_epi32(m+16, s);
}

/* load store double */
static gmx_inline __m512d
gmx_simd_loadu_d_mic(double *m)
{
    return _mm512_loadunpackhi_pd(_mm512_loadunpacklo_pd(_mm512_undefined_pd(), m), m+8);
}

static gmx_inline void
gmx_simd_storeu_d_mic(double* m, __m512d s)
{
    _mm512_packstorelo_pd(m, s);
    _mm512_packstorehi_pd(m+8, s);
}

/* load store dint32 */
static gmx_inline __m512i
gmx_simd_loadu_di_mic(gmx_int32_t *m)
{
    return _mm512_mask_loadunpackhi_epi32(_mm512_mask_loadunpacklo_epi32(_mm512_undefined_epi32(), mask_loh, m), mask_loh, m+16);
}

static gmx_inline void
gmx_simd_storeu_di_mic(gmx_int32_t* m, __m512i s)
{
    _mm512_mask_packstorelo_epi32(m, mask_loh, s);
    _mm512_mask_packstorehi_epi32(m+16, mask_loh, s);
}

/* extract */
static gmx_inline gmx_int32_t
gmx_simd_extract_fi_mic(gmx_simd_fint32_t a, int index)
{
    int r;
    _mm512_mask_packstorelo_epi32(&r, _mm512_mask2int(1<<index), a);
    return r;
}

static gmx_inline gmx_int32_t
gmx_simd_extract_di_mic(gmx_simd_dint32_t a, int index)
{
    int r;
    _mm512_mask_packstorelo_epi32(&r, _mm512_mask2int(1<<index), a);
    return r;
}

/* This is likely faster than the built in scale operation (lat 8, t-put 3)
 * since we only work on the integer part and use shifts. TODO: check. given that scale also only does integer
 */
static gmx_inline __m512
gmx_simd_set_exponent_f_mic(__m512 a)
{
    __m512i       iexp         = gmx_simd_cvt_f2i(a);

    const __m512i expbias      = _mm512_set1_epi32(127);
    iexp = _mm512_slli_epi32(_mm512_add_epi32(iexp, expbias), 23);
    return _mm512_castsi512_ps(iexp);

    /* scale alternative:
       return _mm512_scale_ps(_mm512_set1_ps(1), iexp);
     */
}

static gmx_inline __m512d
gmx_simd_set_exponent_d_mic(__m512d a)
{
    /* SIC: We use single storage as an intermediate here, so bias/mantissa _should_ be 127 & 23 */
    const __m512i expbias      = _mm512_set1_epi32(127);
    __m512i       iexp         = _mm512_cvtfxpnt_roundpd_epi32lo(a, _MM_FROUND_TO_NEAREST_INT);

    iexp = _mm512_slli_epi32(_mm512_add_epi32(iexp, expbias), 23);
    return _mm512_cvtpslo_pd(_mm512_castsi512_ps(iexp));
}

static gmx_inline void
gmx_simd_cvt_f2dd_mic(__m512 f, __m512d * d0, __m512d * d1)
{
    /* I might have this permute constant backwards ... */
    __m512i i1 = _mm512_permute4f128_epi32(_mm512_castps_si512(f), _MM_PERM_CDCD);

    *d0 = _mm512_cvtpslo_pd(f);
    *d1 = _mm512_cvtpslo_pd(_mm512_castsi512_ps(i1));
}

static gmx_inline __m512
gmx_simd_cvt_dd2f_mic(__m512d d0, __m512d d1)
{
    __m512 f0 = _mm512_cvtpd_pslo(d0);
    __m512 f1 = _mm512_cvtpd_pslo(d1);
    return _mm512_mask_permute4f128_ps(f0, mask_hih, f1, PERM_LOW2HIGH);
}

static gmx_inline __m512
gmx_simd_exp2_f_mic(__m512 x)
{
    return _mm512_exp223_ps(_mm512_cvtfxpnt_round_adjustps_epi32(x, _MM_ROUND_MODE_NEAREST, _MM_EXPADJ_24));
}

static gmx_inline __m512
gmx_simd_exp_f_mic(__m512 x)
{
    return gmx_simd_exp2_f(_mm512_mul_ps(x, _mm512_set1_ps(1.44269504088896341)));
}

static gmx_inline __m512
gmx_simd_log_f_mic(__m512 x)
{
    return _mm512_mul_ps(_mm512_set1_ps(0.693147180559945286226764), _mm512_log2ae23_ps(x));
}

/* Presently we do not bother with SIMD4 on MIC, since it just does kernels */

/* Function to check whether SIMD operations have resulted in overflow */
static int
gmx_simd_check_and_reset_overflow(void)
{
    int MXCSR;
    int sse_overflow;

    MXCSR = _mm_getcsr();
    /* The overflow flag is bit 3 in the register */
    if (MXCSR & 0x0008)
    {
        sse_overflow = 1;
        /* Set the overflow flag to zero */
        MXCSR = MXCSR & 0xFFF7;
        _mm_setcsr(MXCSR);
    }
    else
    {
        sse_overflow = 0;
    }
    return sse_overflow;
}

#endif /* GMX_SIMD_IMPL_INTEL_MIC_H */
