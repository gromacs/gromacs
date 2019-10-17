/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2019, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_SPARC64_HPC_ACE_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_SPARC64_HPC_ACE_SIMD_DOUBLE_H

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

#include <cmath>
#include <cstdint>

#include "impl_sparc64_hpc_ace_common.h"

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdDouble _fjsp_v2r8
#define simdLoadD _fjsp_load_v2r8
#define simdLoad1D(m) _fjsp_set_v2r8((*m), (*m))
#define simdSet1D(a) _fjsp_set_v2r8(a, a)
#define simdStoreD _fjsp_store_v2r8
#define simdLoadUD simdLoadD
/* No unaligned store of SimdDouble */
#define simdSetZeroD _fjsp_setzero_v2r8
#define simdAddD _fjsp_add_v2r8
#define simdSubD _fjsp_sub_v2r8
#define simdMulD _fjsp_mul_v2r8
#define simdFmaddD(a, b, c) _fjsp_madd_v2r8(a, b, c)
#define simdFmsubD(a, b, c) _fjsp_msub_v2r8(a, b, c)
#define simdFnmaddD(a, b, c) _fjsp_nmsub_v2r8(a, b, c)
#define simdFnmsubD(a, b, c) _fjsp_nmadd_v2r8(a, b, c)
#define simdAndD _fjsp_and_v2r8
#define simdAndNotD _fjsp_andnot1_v2r8
#define simdOrD _fjsp_or_v2r8
#define simdXorD _fjsp_xor_v2r8
#define simdRsqrtD(x) _fjsp_rsqrta_v2r8(x)
#define simdRcpD(x) _fjsp_rcpa_v2r8(x)
#define simdAbsD(x) _fjsp_abs_v2r8(x)
#define simdNegD(x) _fjsp_neg_v2r8(x)
#define simdMaxD _fjsp_max_v2r8
#define simdMinD _fjsp_min_v2r8
#define simdRoundD(x) simdCvtI2D(simdCvtD2I(x))
#define simdTruncD(x) simdCvtI2D(simdCvttD2I(x))
#define simdFractionD(x) simdSubD(x, simdTruncD(x))
#define simdGetExponentD simdGetExponentD_sparc64_hpc_ace
#define simdGetMantissaD simdGetMantissaD_sparc64_hpc_ace
#define simdSetExponentD simdSetExponentD_sparc64_hpc_ace
/* integer datatype corresponding to double: SimdDInt32 */
#define SimdDInt32 _fjsp_v2r8
#define simdLoadDI(m) simdLoadDI_sparc64_hpc_ace(m)
#define simdSet1DI(i) simdSet1DI_sparc64_hpc_ace(i)
#define simdStoreDI(m, x) simdStoreDI_sparc64_hpc_ace(m, x)
#define simdLoadUDI simdLoadDI
/* No unaligned store of SimdDInt32 */
#define simdSetZeroDI _fjsp_setzero_v2r8
#define simdCvtD2I simdCvtD2I_sparc64_hpc_ace
#define simdCvttD2I _fjsp_dtox_v2r8
#define simdCvtI2D _fjsp_xtod_v2r8
#define simdExtractDI simdExtractDI_sparc64_hpc_ace
/* Integer logical ops on SimdDInt32 */
#define simdSlliDI simdSlliDI_sparc64_hpc_ace
#define simdSrliDI simdSrliDI_sparc64_hpc_ace
#define simdAndDI _fjsp_and_v2r8
#define simdAndNotDI _fjsp_andnot1_v2r8
#define simdOrDI _fjsp_or_v2r8
#define simdXorDI _fjsp_xor_v2r8
/* Integer arithmetic ops on integer datatype corresponding to double */
/* Boolean & comparison operations on SimdDouble */
#define SimdDBool _fjsp_v2r8
#define simdCmpEqD _fjsp_cmpeq_v2r8
#define simdCmpLtD _fjsp_cmplt_v2r8
#define simdCmpLeD _fjsp_cmple_v2r8
#define simdAndDB _fjsp_and_v2r8
#define simdOrDB _fjsp_or_v2r8
#define simdAnyTrueDB gmx_simd_anytrue_d_sparc64_hpc_ace
#define simdMaskD _fjsp_and_v2r8
#define simdMaskNotD(a, sel) _fjsp_andnot1_v2r8(sel, a)
#define simdBlendD(a, b, sel) _fjsp_selmov_v2r8(b, a, sel)
#define simdReduceD(a) simdReduceD_sparc64_hpc_ace(a)

/* No boolean & comparison operations on SimdDInt32 */
/* Float/double conversion */
#define simdCvtF2D(f) (f)
#define simdCvtD2F(d) (d)


/****************************************************
 * DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static inline SimdDInt32 simdLoadDI_sparc64_hpc_ace(const int* m)
{
    union {
        _fjsp_v2r8    simd;
        long long int i[2];
    } conv;

    conv.i[0] = m[0];
    conv.i[1] = m[1];

    return _fjsp_load_v2r8((double*)&(conv.simd));
}

static inline void simdStoreDI_sparc64_hpc_ace(int* m, SimdDInt32 x)
{
    union {
        _fjsp_v2r8    simd;
        long long int i[2];
    } conv;

    _fjsp_store_v2r8((double*)&(conv.simd), x);

    m[0] = conv.i[0];
    m[1] = conv.i[1];
}

static inline SimdDInt32 simdSet1DI_sparc64_hpc_ace(int i)
{
    union {
        _fjsp_v2r8    simd;
        long long int i[2];
    } conv;

    conv.i[0] = i;
    conv.i[1] = i;

    return _fjsp_load_v2r8((double*)&(conv.simd));
}

static inline int simdExtractDI_sparc64_hpc_ace(SimdDInt32 x, int i)
{
    long long int res;
    /* This conditional should be optimized away at compile time */
    if (i == 0)
    {
        _fjsp_storel_v2r8((double*)&res, x);
    }
    else
    {
        _fjsp_storeh_v2r8((double*)&res, x);
    }
    return (int)res;
}

static inline SimdDInt32 simdSlliDI_sparc64_hpc_ace(SimdDInt32 x, int i)
{
    _fjsp_v2i8 ix = *((_fjsp_v2i8*)&x);
    ix            = _fjsp_slli_v2i8(ix, i);
    x             = *((_fjsp_v2r8*)&ix);
    return x;
}

static inline SimdDInt32 simdSrliDI_sparc64_hpc_ace(SimdDInt32 x, int i)
{
    _fjsp_v2i8 ix = *((_fjsp_v2i8*)&x);
    ix            = _fjsp_srli_v2i8(ix, i);
    x             = *((_fjsp_v2r8*)&ix);
    return x;
}

static inline SimdDInt32 simdCvtD2I_sparc64_hpc_ace(SimdDouble x)
{
    _fjsp_v2r8 signbit = _fjsp_set_v2r8(-0.0, -0.0);
    _fjsp_v2r8 half    = _fjsp_set_v2r8(0.5, 0.5);

    x = _fjsp_add_v2r8(x, _fjsp_or_v2r8(_fjsp_and_v2r8(signbit, x), half));
    return _fjsp_dtox_v2r8(x);
}

static inline int gmx_simd_anytrue_d_sparc64_hpc_ace(SimdDBool x)
{
    long long int i;
    x = _fjsp_or_v2r8(x, _fjsp_unpackhi_v2r8(x, x));
    _fjsp_storel_v2r8((double*)&i, x);
    return (i != 0LL);
}

static inline double simdReduceD_sparc64_hpc_ace(SimdDouble x)
{
    double d;
    x = _fjsp_add_v2r8(x, _fjsp_unpackhi_v2r8(x, x));
    _fjsp_storel_v2r8(&d, x);
    return d;
}


static inline SimdDouble simdGetExponentD_sparc64_hpc_ace(SimdDouble x)
{
    /* HPC-ACE cannot cast _fjsp_v2r8 to _fjsp_v4i4, so to perform shifts we
     * would need to store and reload. Since we are only operating on two
     * numbers it is likely more efficient to do the operations directly on
     * normal registers.
     */
    const std::int64_t expmask = 0x7ff0000000000000LL;
    const std::int64_t expbias = 1023LL;

    union {
        _fjsp_v2r8    simd;
        long long int i[2];
    } conv;

    _fjsp_store_v2r8((double*)&conv.simd, x);
    conv.i[0] = ((conv.i[0] & expmask) >> 52) - expbias;
    conv.i[1] = ((conv.i[1] & expmask) >> 52) - expbias;
    x         = _fjsp_load_v2r8((double*)&conv.simd);
    return _fjsp_xtod_v2r8(x);
}

static inline SimdDouble simdGetMantissaD_sparc64_hpc_ace(SimdDouble x)
{
    std::int64_t mantmask[2] = { 0x000fffffffffffffLL, 0x000fffffffffffffLL };
    SimdDouble   one         = _fjsp_set_v2r8(1.0, 1.0);

    x = _fjsp_and_v2r8(x, _fjsp_load_v2r8((double*)mantmask));
    return _fjsp_or_v2r8(x, one);
}

static inline SimdDouble simdSetExponentD_sparc64_hpc_ace(SimdDouble x)
{
    const std::int64_t expbias = 1023;
    union {
        _fjsp_v2r8    simd;
        long long int i[2];
    } conv;


    _fjsp_store_v2r8((double*)&conv.simd, simdCvtD2I_sparc64_hpc_ace(x));
    conv.i[0] = (conv.i[0] + expbias) << 52;
    conv.i[1] = (conv.i[1] + expbias) << 52;

    return _fjsp_load_v2r8((double*)&conv.simd);
}

#endif /* GMX_SIMD_IMPL_SPARC64_HPC_ACE_SIMD_DOUBLE_H */
