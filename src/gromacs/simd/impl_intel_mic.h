/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS Development Team.
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
#undef  GMX_SIMD_HAVE_LOADU
#undef  GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#undef  GMX_SIMD_HAVE_FINT32_EXTRACT
#define GMX_SIMD_HAVE_FINT32_LOGICAL     /* AVX1 cannot do 256-bit int shifts */
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS /* AVX1 cannot do 256-bit int +,-,*  */
#define GMX_SIMD_HAVE_DINT32
#undef  GMX_SIMD_HAVE_DINT32_EXTRACT
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
/* No support for unaligned load/store */
#define gmx_simd_setzero_f         _mm512_setzero_ps
#define gmx_simd_add_f             _mm512_add_ps
#define gmx_simd_sub_f             _mm512_sub_ps
#define gmx_simd_mul_f             _mm512_mul_ps
#define gmx_simd_fmadd_f           _mm512_fmadd_ps
#define gmx_simd_fmsub_f           _mm512_fmsub_ps
#define gmx_simd_fnmadd_f          _mm512_fnmadd_ps
#define gmx_simd_fnmsub_f          _mm512_fnmsub_ps
#define gmx_simd_and_f(a, b)        _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_andnot_f(a, b)     _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_or_f(a, b)         _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_xor_f(a, b)        _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_rsqrt_f           _mm512_rsqrt23_ps
#define gmx_simd_rcp_f             _mm512_rcp23_ps
#define gmx_simd_fabs_f(x)         gmx_simd_andnot_f(_mm512_set1_ps(-1.0), x)
#define gmx_simd_fneg_f(x)         _mm512_addn_ps(x, _mm512_setzero_ps())
#define gmx_simd_max_f             _mm512_gmax_ps
#define gmx_simd_min_f             _mm512_gmin_ps
#define gmx_simd_round_f(x)        _mm512_round_ps(x, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd_trunc_f(x)        _mm512_round_ps(x, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd_fraction_f(x)     _mm512_sub_ps(x, gmx_simd_trunc_f(x))
#define gmx_simd_get_exponent_f(x) _mm512_getexp_ps(x)
#define gmx_simd_get_mantissa_f(x) _mm512_getmant_ps(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define gmx_simd_set_exponent_f(x) gmx_simd_set_exponent_f_mic
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t          __m512i
#define gmx_simd_load_fi           _mm512_load_epi32
#define gmx_simd_set1_fi           _mm512_set1_ps
#define gmx_simd_store_fi(m, x)     _mm512_store_epi32
/* No support for unaligned load/store */
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
#define gmx_simd_cmpeq_f(a, b)      _mm512_cmp_ps_mask(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmplt_f(a, b)      _mm512_cmp_ps_mask(a, b, _CMP_LT_OS)
#define gmx_simd_cmple_f(a, b)      _mm512_cmp_ps_mask(a, b, _CMP_LE_OS)
#define gmx_simd_and_fb            _mm512_kand
#define gmx_simd_or_fb             _mm512_kor
#define gmx_simd_anytrue_fb        _mm512_mask2int
#define gmx_simd_blendzero_f(a, sel) _mm512_mask_mov_ps(_mm512_setzero_ps(), sel, a)
#define gmx_simd_blendv_f(a, b, sel)  _mm512_mask_blend_ps(sel, a, b)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t          __mmask16
#define gmx_simd_cmpeq_fi(a, b)     _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define gmx_simd_cmplt_fi(a, b)     _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define gmx_simd_and_fib           _mm512_kand
#define gmx_simd_or_fib            _mm512_kor
#define gmx_simd_anytrue_fib       _mm512_mask2int
#define gmx_simd_blendzero_fi(a, sel) _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define gmx_simd_blendv_fi(a, b, sel)  _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib(x)     (x)
#define gmx_simd_cvt_fib2fb(x)     (x)

/* MIC provides full single precision of some neat functions: */
/* 1/sqrt(x) and 1/x work fine in simd_math.h, and won't use extra iterations */
#define gmx_simd_exp2_f(x)         _mm512_scale_ps(_mm512_set1_ps(1.0), x)
#define gmx_simd_exp_f(x)          gmx_simd_exp2_f(_mm512_mul_ps(x, _mm512_set1_ps(1.44269504088896341)))
#define gmx_simd_log_f(x)          _mm512_mul_ps(_mm512_set1_ps(0.693147180559945286226764g), _mm512_log2ae23_ps(x))

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t          __m512
#define gmx_simd_load_d            _mm512_load_pd
#define gmx_simd_load1_d(m)        _mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE)
#define gmx_simd_set1_d            _mm512_set1_pd
#define gmx_simd_store_d           _mm512_store_pd
/* No support for unaligned load/store */
#define gmx_simd_setzero_d         _mm512_setzero_pd
#define gmx_simd_add_d             _mm512_add_pd
#define gmx_simd_sub_d             _mm512_sub_pd
#define gmx_simd_mul_d             _mm512_mul_pd
#define gmx_simd_fmadd_d           _mm512_fmadd_pd
#define gmx_simd_fmsub_d           _mm512_fmsub_pd
#define gmx_simd_fnmadd_d          _mm512_fnmadd_pd
#define gmx_simd_fnmsub_d          _mm512_fnmsub_pd
#define gmx_simd_and_d(a, b)        _mm512_castsi512_pd(_mm512_and_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_andnot_d(a, b)     _mm512_castsi512_pd(_mm512_andnot_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_or_d(a, b)         _mm512_castsi512_pd(_mm512_or_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_xor_d(a, b)        _mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_rsqrt_d           _mm512_rsqrt23_ps
#define gmx_simd_rcp_d             _mm512_rcp23_ps
#define gmx_simd_fabs_d(x)         gmx_simd_andnot_d(_mm512_set1_pd(-1.0), x)
#define gmx_simd_fneg_d(x)         _mm512_addn_pd(x, _mm512_setzero_pd())
#define gmx_simd_max_d             _mm512_gmax_pd
#define gmx_simd_min_d             _mm512_gmin_pd
#define gmx_simd_round_d(a)        _mm512_roundfxpnt_adjust_pd(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd_trunc_d(a)        _mm512_roundfxpnt_adjust_pd(a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd_fraction_d(x)     _mm512_sub_pd(x, gmx_simd_trunc_d(x))
#define gmx_simd_get_exponent_d(x) _mm512_getexp_pd(x)
#define gmx_simd_get_mantissa_d(x) _mm512_getmant_pd(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define gmx_simd_set_exponent_d(x) gmx_simd_set_exponent_d_mic
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_dint32_t          __m512i
#define gmx_simd_load_di           gmx_simd_load_di_mic
#define gmx_simd_set1_di           _mm512_set1_epi32
#define gmx_simd_store_di          gmx_simd_store_di_mic
/* No support for unaligned load/store */
#define gmx_simd_setzero_di        _mm512_setzero_epi32
#define gmx_simd_cvt_d2i(a)        _mm512_cvt_roundpd_epi32lo(a, _MM_FROUND_TO_NEAREST_INT)
#define gmx_simd_cvtt_d2i(a)       _mm512_cvt_roundpd_epi32lo(a, _MM_FROUND_TO_ZERO)
#define gmx_simd_cvt_i2d(a)        _mm512_cvtepi32lo_pd
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
#define gmx_simd_cmpeq_d(a, b)      _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmplt_d(a, b)      _mm512_cmp_pd_mask(a, b, _CMP_LT_OS)
#define gmx_simd_cmple_d(a, b)      _mm512_cmp_pd_mask(a, b, _CMP_LE_OS)
#define gmx_simd_and_db            _mm512_kand
#define gmx_simd_or_db             _mm512_kor
#define gmx_simd_anytrue_db(x)     _mm512_mask2int(_mm512_kmovlhb(x, x))
#define gmx_simd_blendzero_d(a, sel) _mm512_mask_mov_pd(_mm512_setzero_pd(), sel, a)
#define gmx_simd_blendv_d(a, b, sel)  _mm512_mask_blend_pd(sel, a, b)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_dibool_t          __mmask16
#define gmx_simd_cmpeq_di(a, b)     _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define gmx_simd_cmplt_di(a, b)     _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define gmx_simd_and_dib           _mm512_kand
#define gmx_simd_or_dib            _mm512_kor
#define gmx_simd_anytrue_dib       _mm512_mask2int
#define gmx_simd_blendzero_di(a, sel) _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define gmx_simd_blendv_di(a, b, sel)  _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between booleans. Double & dint stuff is stored in low bits */
#define gmx_simd_cvt_db2dib(x)     (x)
#define gmx_simd_cvt_dib2db(x)     (x)

/* Float/double conversion */
#define gmx_simd_cvt_f2dd          gmx_simd_cvt_f2dd_mic
#define gmx_simd_cvt_dd2f          gmx_simd_cvt_dd2f_mic

/* This is likely faster than the built in scale operation (lat 8, t-put 3)
 * since we only work on the integer part and use shifts.
 */
static gmx_inline __m512
gmx_simd_set_exponent_f_mic(__m512 a)
{
    const __m512i expbias      = _mm512_set1_epi32(127);
    const  int    mantissabits = 23;
    __m512i       iexp         = gmx_simd_cvt_f2i(a);

    iexp = _mm512_slli_epi32(_mm512_add_epi32(iexp, expbias), mantissabits);
    return _mm512_castsi512_ps(iexp);
}

static gmx_inline __m512d
gmx_simd_load_di(int *m)
{
    return _mm512_castsi512_pd(_mm512_extload_epi64(m, _MM_UPCONV_EPI64_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
}

static gmx_inline void
gmx_simd_store_di(int *m, x)
{
    __m512i index = _mm512_set_epi32(16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    _mm512_mask_i32scatter_epi32(m, _mm512_int2mask(0x00FF), index, x, _MM_SCALE_1);
}

static gmx_inline __m512d
gmx_simd_set_exponent_d_mic(__m512d a)
{
    /* SIC: We use single storage as an intermediate here, so bias/mantissa _should_ be 127 & 23 */
    const __m512i expbias      = _mm512_set1_epi32(127);
    const  int    mantissabits = 23;
    __m512i       iexp         = _mm512_cvt_roundpd_epi32lo(a, _MM_FROUND_TO_NEAREST_INT);

    iexp = _mm512_slli_epi32(_mm512_add_epi32(iexp, expbias), mantissabits);
    return _mm512_castpslo_pd(_mm512_castsi512_ps(iexp));
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
gmx_simd_cvt_dd2f_mic(__m512 d0, __m512 d1)
{
    __m512i i0 = _mm512_castps_si512(_mm512_cvtpd_pslo(d0));
    __m512i i1 = _mm512_castps_si512(_mm512_cvtpd_pslo(d1))
        /* I might have both the mask and permute constant backwards ... */
        return _mm512_castsi512_pd(_mm512_mask_permute4f128_epi32(f0, _mm512_int2mask(0xFF00), f1, _MM_PERM_ABAB));
}

/* Presently we do not bother with SIMD4 on MIC, since it just does kernels */

#endif /* GMX_SIMD_IMPL_INTEL_MIC_H */
