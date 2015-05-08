/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_SPARC64_HPC_ACE_H
#define GMX_SIMD_IMPL_SPARC64_HPC_ACE_H

#include <math.h>

/* Fujitsu header borrows the name from SSE2, since some instructions have aliases.
 * Environment/compiler version GM-1.2.0-17 seems to be buggy; when -Xg is
 * defined to enable GNUC extensions, this sets _ISOC99_SOURCE, which in
 * turn causes all intrinsics to be declared inline _instead_ of static. This
 * leads to duplicate symbol errors at link time.
 * To work around this we unset this before including the HPC-ACE header, and
 * reset the value afterwards.
 */
#ifdef _ISOC99_SOURCE
#    undef _ISOC99_SOURCE
#    define SAVE_ISOC99_SOURCE
#endif

#include <emmintrin.h>

#ifdef SAVE_ISOC99_SOURCE
#    define _ISOC99_SOURCE
#    undef SAVE_ISOC99_SOURCE
#endif


/* Sparc64 HPC-ACE SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for defines.
 */

/* Capability definitions for Sparc64 HPC-ACE */
/* HPC-ACE is actually double-only on the register level, but we also implement
 * a single-precision interface where we only offer single-precision accuracy
 * in math functions - this can save quite a few cycles.
 */
#define GMX_SIMD_HAVE_FLOAT
#define GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_HARDWARE
#undef  GMX_SIMD_HAVE_LOADU
#undef  GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#define GMX_SIMD_HAVE_FINT32_EXTRACT
#define GMX_SIMD_HAVE_FINT32_LOGICAL
#undef  GMX_SIMD_HAVE_FINT32_ARITHMETICS
#define GMX_SIMD_HAVE_DINT32
#define GMX_SIMD_HAVE_DINT32_EXTRACT
#define GMX_SIMD_HAVE_DINT32_LOGICAL
#undef  GMX_SIMD_HAVE_DINT32_ARITHMETICS
#undef  GMX_SIMD4_HAVE_FLOAT
#undef  GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH         2
#define GMX_SIMD_DOUBLE_WIDTH        2
#define GMX_SIMD_FINT32_WIDTH        2
#define GMX_SIMD_DINT32_WIDTH        2
#define GMX_SIMD_RSQRT_BITS         10
#define GMX_SIMD_RCP_BITS            9

/* HPC-ACE is a bit strange; some instructions like
 * shifts only work on _integer_ versions of SIMD
 * registers, but there are no intrinsics to load
 * or convert, or even to cast. The only way to use
 * them is to declare unions with the SIMD integer
 * type. However, this will lead to extra load ops,
 * and the normal real-to-int and int-to-real
 * conversions work purely on the v2r8 fp regs.
 * Since our most common usage is to convert and
 * then extract the result for table lookups, we
 * define the gmx_simd_fint32_t datatype to use
 * the v2r8 rather than v2i8 SIMD type.
 */

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t          _fjsp_v2r8
#define gmx_simd_load_f           gmx_simd_load_f_sparc64_hpc_ace
#define gmx_simd_load1_f(m)       _fjsp_set_v2r8((*m), (*m))
#define gmx_simd_set1_f(a)        _fjsp_set_v2r8(a, a)
#define gmx_simd_store_f          gmx_simd_store_f_sparc64_hpc_ace
#define gmx_simd_loadu_f          gmx_simd_load_f
/* No unaligned store of gmx_simd_float_t */
#define gmx_simd_setzero_f        _fjsp_setzero_v2r8
#define gmx_simd_add_f            _fjsp_add_v2r8
#define gmx_simd_sub_f            _fjsp_sub_v2r8
#define gmx_simd_mul_f            _fjsp_mul_v2r8
#define gmx_simd_fmadd_f(a, b, c)   _fjsp_madd_v2r8(a, b, c)
#define gmx_simd_fmsub_f(a, b, c)   _fjsp_msub_v2r8(a, b, c)
#define gmx_simd_fnmadd_f(a, b, c)  _fjsp_nmsub_v2r8(a, b, c)
#define gmx_simd_fnmsub_f(a, b, c)  _fjsp_nmadd_v2r8(a, b, c)
#define gmx_simd_and_f            _fjsp_and_v2r8
#define gmx_simd_andnot_f         _fjsp_andnot1_v2r8
#define gmx_simd_or_f             _fjsp_or_v2r8
#define gmx_simd_xor_f            _fjsp_xor_v2r8
#define gmx_simd_rsqrt_f          _fjsp_rsqrta_v2r8
#define gmx_simd_rcp_f            _fjsp_rcpa_v2r8
#define gmx_simd_fabs_f(x)        _fjsp_abs_v2r8(x)
#define gmx_simd_fneg_f(x)        _fjsp_neg_v2r8(x)
#define gmx_simd_max_f            _fjsp_max_v2r8
#define gmx_simd_min_f            _fjsp_min_v2r8
#define gmx_simd_round_f(x)       gmx_simd_round_d(x)
#define gmx_simd_trunc_f(x)       gmx_simd_trunc_d(x)
#define gmx_simd_fraction_f(x)    gmx_simd_sub_f(x, gmx_simd_trunc_f(x))
#define gmx_simd_get_exponent_f   gmx_simd_get_exponent_d_sparc64_hpc_ace
#define gmx_simd_get_mantissa_f   gmx_simd_get_mantissa_d_sparc64_hpc_ace
#define gmx_simd_set_exponent_f   gmx_simd_set_exponent_d_sparc64_hpc_ace
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t         _fjsp_v2r8
#define gmx_simd_load_fi(m)       gmx_simd_load_di_sparc64_hpc_ace(m)
#define gmx_simd_set1_fi(i)       gmx_simd_set1_di_sparc64_hpc_ace(i)
#define gmx_simd_store_fi(m, x)   gmx_simd_store_di_sparc64_hpc_ace(m, x)
#define gmx_simd_loadu_fi         gmx_simd_load_fi
/* No unaligned store of gmx_simd_fint32_t */
#define gmx_simd_setzero_fi       _fjsp_setzero_v2r8
#define gmx_simd_cvt_f2i          gmx_simd_cvt_d2i
#define gmx_simd_cvtt_f2i         _fjsp_dtox_v2r8
#define gmx_simd_cvt_i2f          _fjsp_xtod_v2r8
#define gmx_simd_extract_fi      gmx_simd_extract_di_sparc64_hpc_ace
/* Integer logical ops on gmx_simd_fint32_t */
/* Shifts are horrible since they require memory re-loads. */
#define gmx_simd_slli_fi          gmx_simd_slli_di_sparc64_hpc_ace
#define gmx_simd_srli_fi          gmx_simd_srli_di_sparc64_hpc_ace
#define gmx_simd_and_fi           _fjsp_and_v2r8
#define gmx_simd_andnot_fi(a, b)   _fjsp_andnot1_v2r8(a, b)
#define gmx_simd_or_fi            _fjsp_or_v2r8
#define gmx_simd_xor_fi           _fjsp_xor_v2r8
/* No integer arithmetic ops on gmx_simd_fint32_t */
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t          _fjsp_v2r8
#define gmx_simd_cmpeq_f          _fjsp_cmpeq_v2r8
#define gmx_simd_cmplt_f          _fjsp_cmplt_v2r8
#define gmx_simd_cmple_f          _fjsp_cmple_v2r8
#define gmx_simd_and_fb           _fjsp_and_v2r8
#define gmx_simd_or_fb            _fjsp_or_v2r8
#define gmx_simd_anytrue_fb       gmx_simd_anytrue_d_sparc64_hpc_ace
#define gmx_simd_blendzero_f      _fjsp_and_v2r8
#define gmx_simd_blendnotzero_f(a, sel) _fjsp_andnot1_v2r8(sel, a)
#define gmx_simd_blendv_f(a, b, s) _fjsp_selmov_v2r8(b, a, s)
#define gmx_simd_reduce_f(a)       gmx_simd_reduce_d_sparc64_hpc_ace(a)
/* No boolean & comparison operations on gmx_simd_fint32_t */
/* No conversions between different booleans */

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t          _fjsp_v2r8
#define gmx_simd_load_d            _fjsp_load_v2r8
#define gmx_simd_load1_d(m)        _fjsp_set_v2r8((*m), (*m))
#define gmx_simd_set1_d(a)         _fjsp_set_v2r8(a, a)
#define gmx_simd_store_d           _fjsp_store_v2r8
#define gmx_simd_loadu_d           gmx_simd_load_d
/* No unaligned store of gmx_simd_double_t */
#define gmx_simd_setzero_d         _fjsp_setzero_v2r8
#define gmx_simd_add_d             _fjsp_add_v2r8
#define gmx_simd_sub_d             _fjsp_sub_v2r8
#define gmx_simd_mul_d             _fjsp_mul_v2r8
#define gmx_simd_fmadd_d(a, b, c)   _fjsp_madd_v2r8(a, b, c)
#define gmx_simd_fmsub_d(a, b, c)   _fjsp_msub_v2r8(a, b, c)
#define gmx_simd_fnmadd_d(a, b, c)  _fjsp_nmsub_v2r8(a, b, c)
#define gmx_simd_fnmsub_d(a, b, c)  _fjsp_nmadd_v2r8(a, b, c)
#define gmx_simd_and_d             _fjsp_and_v2r8
#define gmx_simd_andnot_d          _fjsp_andnot1_v2r8
#define gmx_simd_or_d              _fjsp_or_v2r8
#define gmx_simd_xor_d             _fjsp_xor_v2r8
#define gmx_simd_rsqrt_d(x)        _fjsp_rsqrta_v2r8(x)
#define gmx_simd_rcp_d(x)          _fjsp_rcpa_v2r8(x)
#define gmx_simd_fabs_d(x)         _fjsp_abs_v2r8(x)
#define gmx_simd_fneg_d(x)         _fjsp_neg_v2r8(x)
#define gmx_simd_max_d             _fjsp_max_v2r8
#define gmx_simd_min_d             _fjsp_min_v2r8
#define gmx_simd_round_d(x)        gmx_simd_cvt_i2d(gmx_simd_cvt_d2i(x))
#define gmx_simd_trunc_d(x)        gmx_simd_cvt_i2d(gmx_simd_cvtt_d2i(x))
#define gmx_simd_fraction_d(x)     gmx_simd_sub_d(x, gmx_simd_trunc_d(x))
#define gmx_simd_get_exponent_d    gmx_simd_get_exponent_d_sparc64_hpc_ace
#define gmx_simd_get_mantissa_d    gmx_simd_get_mantissa_d_sparc64_hpc_ace
#define gmx_simd_set_exponent_d    gmx_simd_set_exponent_d_sparc64_hpc_ace
/* integer datatype corresponding to double: gmx_simd_dint32_t */
#define gmx_simd_dint32_t          _fjsp_v2r8
#define gmx_simd_load_di(m)        gmx_simd_load_di_sparc64_hpc_ace(m)
#define gmx_simd_set1_di(i)        gmx_simd_set1_di_sparc64_hpc_ace(i)
#define gmx_simd_store_di(m, x)    gmx_simd_store_di_sparc64_hpc_ace(m, x)
#define gmx_simd_loadu_di          gmx_simd_load_di
/* No unaligned store of gmx_simd_dint32_t */
#define gmx_simd_setzero_di        _fjsp_setzero_v2r8
#define gmx_simd_cvt_d2i           gmx_simd_cvt_d2i_sparc64_hpc_ace
#define gmx_simd_cvtt_d2i          _fjsp_dtox_v2r8
#define gmx_simd_cvt_i2d           _fjsp_xtod_v2r8
#define gmx_simd_extract_di        gmx_simd_extract_di_sparc64_hpc_ace
/* Integer logical ops on gmx_simd_dint32_t */
#define gmx_simd_slli_di           gmx_simd_slli_di_sparc64_hpc_ace
#define gmx_simd_srli_di           gmx_simd_srli_di_sparc64_hpc_ace
#define gmx_simd_and_di            _fjsp_and_v2r8
#define gmx_simd_andnot_di         _fjsp_andnot1_v2r8
#define gmx_simd_or_di             _fjsp_or_v2r8
#define gmx_simd_xor_di            _fjsp_xor_v2r8
/* Integer arithmetic ops on integer datatype corresponding to double */
/* Boolean & comparison operations on gmx_simd_double_t */
#define gmx_simd_dbool_t           _fjsp_v2r8
#define gmx_simd_cmpeq_d           _fjsp_cmpeq_v2r8
#define gmx_simd_cmplt_d           _fjsp_cmplt_v2r8
#define gmx_simd_cmple_d           _fjsp_cmple_v2r8
#define gmx_simd_and_db            _fjsp_and_v2r8
#define gmx_simd_or_db             _fjsp_or_v2r8
#define gmx_simd_anytrue_db         gmx_simd_anytrue_d_sparc64_hpc_ace
#define gmx_simd_blendzero_d        _fjsp_and_v2r8
#define gmx_simd_blendnotzero_d(a, sel)  _fjsp_andnot1_v2r8(sel, a)
#define gmx_simd_blendv_d(a, b, sel) _fjsp_selmov_v2r8(b, a, sel)
#define gmx_simd_reduce_d(a)        gmx_simd_reduce_d_sparc64_hpc_ace(a)

/* No boolean & comparison operations on gmx_simd_dint32_t */
/* Float/double conversion */
#define gmx_simd_cvt_f2d(f)         (f)
#define gmx_simd_cvt_d2f(d)         (d)


/****************************************************
 * SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_load_f_sparc64_hpc_ace(const float *m)
{
    /* We are not allowed to cast single-to-double registers, but we can
     * masquerade the memory location as a variable of type _fjsp_v2r4.
     */
    const _fjsp_v2r4 *p = (const _fjsp_v2r4 *)m;
    _fjsp_v2r4        simd;

    simd = *p;
    return _fjsp_stod_v2r8(simd);
}

static gmx_inline void
gmx_simd_store_f_sparc64_hpc_ace(float *m, gmx_simd_float_t x)
{
    /* We are not allowed to cast single-to-double registers, but we can
     * masquerade the memory location as a variable of type _fjsp_v2r4.
     */
    _fjsp_v2r4 *p = (_fjsp_v2r4 *)m;
    *p = _fjsp_dtos_v2r4(x);
}

static gmx_inline gmx_simd_dint32_t
gmx_simd_load_di_sparc64_hpc_ace(const int *m)
{
    union
    {
        _fjsp_v2r8       simd;
        long long int    i[2];
    }
    conv;

    conv.i[0] = m[0];
    conv.i[1] = m[1];

    return _fjsp_load_v2r8( (double *) &(conv.simd) );
}

static gmx_inline void
gmx_simd_store_di_sparc64_hpc_ace(int *m, gmx_simd_dint32_t x)
{
    union
    {
        _fjsp_v2r8       simd;
        long long int    i[2];
    }
    conv;

    _fjsp_store_v2r8( (double *) &(conv.simd), x );

    m[0] = conv.i[0];
    m[1] = conv.i[1];
}

static gmx_inline gmx_simd_dint32_t
gmx_simd_set1_di_sparc64_hpc_ace(int i)
{
    union
    {
        _fjsp_v2r8       simd;
        long long int    i[2];
    }
    conv;

    conv.i[0] = i;
    conv.i[1] = i;

    return _fjsp_load_v2r8( (double *) &(conv.simd) );
}

static gmx_inline int
gmx_simd_extract_di_sparc64_hpc_ace(gmx_simd_dint32_t x, int i)
{
    long long int res;
    /* This conditional should be optimized away at compile time */
    if (i == 0)
    {
        _fjsp_storel_v2r8((double *)&res, x);
    }
    else
    {
        _fjsp_storeh_v2r8((double *)&res, x);
    }
    return (int)res;
}

static gmx_inline gmx_simd_dint32_t
gmx_simd_slli_di_sparc64_hpc_ace(gmx_simd_dint32_t x, int i)
{
    _fjsp_v2i8 ix = *((_fjsp_v2i8 *)&x);
    ix = _fjsp_slli_v2i8(ix, i);
    x  = *((_fjsp_v2r8 *)&ix);
    return x;
}

static gmx_inline gmx_simd_dint32_t
gmx_simd_srli_di_sparc64_hpc_ace(gmx_simd_dint32_t x, int i)
{
    _fjsp_v2i8 ix = *((_fjsp_v2i8 *)&x);
    ix = _fjsp_srli_v2i8(ix, i);
    x  = *((_fjsp_v2r8 *)&ix);
    return x;
}

static gmx_inline gmx_simd_dint32_t
gmx_simd_cvt_d2i_sparc64_hpc_ace(gmx_simd_double_t x)
{
    _fjsp_v2r8 signbit = _fjsp_set_v2r8(-0.0, -0.0);
    _fjsp_v2r8 half    = _fjsp_set_v2r8(0.5, 0.5);

    x = _fjsp_add_v2r8(x, _fjsp_or_v2r8(_fjsp_and_v2r8(signbit, x), half));
    return _fjsp_dtox_v2r8(x);
}

static gmx_inline int
gmx_simd_anytrue_d_sparc64_hpc_ace(gmx_simd_dbool_t x)
{
    long long int i;
    x = _fjsp_or_v2r8(x, _fjsp_unpackhi_v2r8(x, x));
    _fjsp_storel_v2r8((double *)&i, x);
    return (i != 0LL);
}

static gmx_inline double
gmx_simd_reduce_d_sparc64_hpc_ace(gmx_simd_double_t x)
{
    double d;
    x = _fjsp_add_v2r8(x, _fjsp_unpackhi_v2r8(x, x));
    _fjsp_storel_v2r8(&d, x);
    return d;
}


static gmx_inline gmx_simd_double_t
gmx_simd_get_exponent_d_sparc64_hpc_ace(gmx_simd_double_t x)
{
    /* HPC-ACE cannot cast _fjsp_v2r8 to _fjsp_v4i4, so to perform shifts we
     * would need to store and reload. Since we are only operating on two
     * numbers it is likely more efficient to do the operations directly on
     * normal registers.
     */
    const gmx_int64_t    expmask   = 0x7ff0000000000000LL;
    const gmx_int64_t    expbias   = 1023LL;

    union
    {
        _fjsp_v2r8       simd;
        long long int    i[2];
    }
    conv;

    _fjsp_store_v2r8( (double *) &conv.simd, x);
    conv.i[0] = ((conv.i[0] & expmask) >> 52) - expbias;
    conv.i[1] = ((conv.i[1] & expmask) >> 52) - expbias;
    x         = _fjsp_load_v2r8( (double *) &conv.simd);
    return _fjsp_xtod_v2r8(x);
}

static gmx_inline gmx_simd_double_t
gmx_simd_get_mantissa_d_sparc64_hpc_ace(gmx_simd_double_t x)
{
    gmx_int64_t       mantmask[2] = {0x000fffffffffffffLL, 0x000fffffffffffffLL};
    gmx_simd_double_t one         = _fjsp_set_v2r8(1.0, 1.0);

    x = _fjsp_and_v2r8(x, _fjsp_load_v2r8((double *)mantmask));
    return _fjsp_or_v2r8(x, one);
}

static gmx_inline gmx_simd_double_t
gmx_simd_set_exponent_d_sparc64_hpc_ace(gmx_simd_double_t x)
{
    const gmx_int64_t    expbias   = 1023;
    union
    {
        _fjsp_v2r8       simd;
        long long int    i[2];
    }
    conv;


    _fjsp_store_v2r8( (double *) &conv.simd, gmx_simd_cvt_d2i_sparc64_hpc_ace(x));
    conv.i[0] = (conv.i[0] + expbias) << 52;
    conv.i[1] = (conv.i[1] + expbias) << 52;

    return _fjsp_load_v2r8( (double *) &conv.simd);
}


/* No SIMD4 support, since both single & double are only 2-wide */


#endif /* GMX_SIMD_IMPL_SPARC64_HPC_ACE_H */
