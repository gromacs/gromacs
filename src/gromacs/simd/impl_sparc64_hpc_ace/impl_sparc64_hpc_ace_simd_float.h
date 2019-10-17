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

#ifndef GMX_SIMD_IMPL_SPARC64_HPC_ACE_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_SPARC64_HPC_ACE_SIMD_FLOAT_H

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

#include <math.h>

#include "impl_sparc64_hpc_ace_common.h"

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
 * define the SimdFInt32 datatype to use
 * the v2r8 rather than v2i8 SIMD type.
 */

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdFloat _fjsp_v2r8
#define simdLoadF simdLoadF_sparc64_hpc_ace
#define simdLoad1F(m) _fjsp_set_v2r8((*m), (*m))
#define simdSet1F(a) _fjsp_set_v2r8(a, a)
#define simdStoreF simdStoreF_sparc64_hpc_ace
#define simdLoadUF simdLoadF
/* No unaligned store of SimdFloat */
#define simdSetZeroF _fjsp_setzero_v2r8
#define simdAddF _fjsp_add_v2r8
#define simdSubF _fjsp_sub_v2r8
#define simdMulF _fjsp_mul_v2r8
#define simdFmaddF(a, b, c) _fjsp_madd_v2r8(a, b, c)
#define simdFmsubF(a, b, c) _fjsp_msub_v2r8(a, b, c)
#define simdFnmaddF(a, b, c) _fjsp_nmsub_v2r8(a, b, c)
#define simdFnmsubF(a, b, c) _fjsp_nmadd_v2r8(a, b, c)
#define simdAndF _fjsp_and_v2r8
#define simdAndNotF _fjsp_andnot1_v2r8
#define simdOrF _fjsp_or_v2r8
#define simdXorF _fjsp_xor_v2r8
#define simdRsqrtF _fjsp_rsqrta_v2r8
#define simdRcpF _fjsp_rcpa_v2r8
#define simdAbsF(x) _fjsp_abs_v2r8(x)
#define simdNegF(x) _fjsp_neg_v2r8(x)
#define simdMaxF _fjsp_max_v2r8
#define simdMinF _fjsp_min_v2r8
#define simdRoundF(x) simdRoundD(x)
#define simdTruncF(x) simdTruncD(x)
#define simdFractionF(x) simdSubF(x, simdTruncF(x))
#define simdGetExponentF simdGetExponentD_sparc64_hpc_ace
#define simdGetMantissaF simdGetMantissaD_sparc64_hpc_ace
#define simdSetExponentF simdSetExponentD_sparc64_hpc_ace
/* integer datatype corresponding to float: SimdFInt32 */
#define SimdFInt32 _fjsp_v2r8
#define simdLoadFI(m) simdLoadDI_sparc64_hpc_ace(m)
#define simdSet1FI(i) simdSet1DI_sparc64_hpc_ace(i)
#define simdStoreFI(m, x) simdStoreDI_sparc64_hpc_ace(m, x)
#define simdLoadUFI simdLoadFI
/* No unaligned store of SimdFInt32 */
#define simdSetZeroFI _fjsp_setzero_v2r8
#define simdCvtF2I simdCvtD2I
#define simdCvttF2I _fjsp_dtox_v2r8
#define simdCvtI2F _fjsp_xtod_v2r8
#define simdExtractFI simdExtractDI_sparc64_hpc_ace
/* Integer logical ops on SimdFInt32 */
/* Shifts are horrible since they require memory re-loads. */
#define simdSlliFI simdSlliDI_sparc64_hpc_ace
#define simdSrliFI simdSrliDI_sparc64_hpc_ace
#define simdAndFI _fjsp_and_v2r8
#define simdAndNotFI(a, b) _fjsp_andnot1_v2r8(a, b)
#define simdOrFI _fjsp_or_v2r8
#define simdXorFI _fjsp_xor_v2r8
/* No integer arithmetic ops on SimdFInt32 */
/* Boolean & comparison operations on SimdFloat */
#define SimdFBool _fjsp_v2r8
#define simdCmpEqF _fjsp_cmpeq_v2r8
#define simdCmpLtF _fjsp_cmplt_v2r8
#define simdCmpLeF _fjsp_cmple_v2r8
#define simdAndFB _fjsp_and_v2r8
#define simdOrFB _fjsp_or_v2r8
#define simdAnyTrueFB gmx_simd_anytrue_d_sparc64_hpc_ace
#define simdMaskF _fjsp_and_v2r8
#define simdMaskNotF(a, sel) _fjsp_andnot1_v2r8(sel, a)
#define simdBlendF(a, b, s) _fjsp_selmov_v2r8(b, a, s)
#define simdReduceF(a) simdReduceD_sparc64_hpc_ace(a)
/* No boolean & comparison operations on SimdFInt32 */
/* No conversions between different booleans */

/****************************************************
 * SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static inline SimdFloat simdLoadF_sparc64_hpc_ace(const float* m)
{
    /* We are not allowed to cast single-to-double registers, but we can
     * masquerade the memory location as a variable of type _fjsp_v2r4.
     */
    const _fjsp_v2r4* p = (const _fjsp_v2r4*)m;
    _fjsp_v2r4        simd;

    simd = *p;
    return _fjsp_stod_v2r8(simd);
}

static inline void simdStoreF_sparc64_hpc_ace(float* m, SimdFloat x)
{
    /* We are not allowed to cast single-to-double registers, but we can
     * masquerade the memory location as a variable of type _fjsp_v2r4.
     */
    _fjsp_v2r4* p = (_fjsp_v2r4*)m;
    *p            = _fjsp_dtos_v2r4(x);
}

/* Note that some single precision defines refer to the double precision helpers */

#endif /* GMX_SIMD_IMPL_SPARC64_HPC_ACE_SIMD_FLOAT_H */
