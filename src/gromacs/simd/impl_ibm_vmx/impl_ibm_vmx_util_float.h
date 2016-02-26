/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
#ifndef GMX_SIMD_IMPL_IBM_VMX_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_IBM_VMX_UTIL_FLOAT_H

#include "config.h"

#include <cstddef>
#include <cstdint>

#include "gromacs/utility/basedefinitions.h"

#include "impl_ibm_vmx_definitions.h"
#include "impl_ibm_vmx_simd_float.h"

namespace gmx
{

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float *        base,
                    const std::int32_t   offset[],
                    SimdFloat *          v0,
                    SimdFloat *          v1,
                    SimdFloat *          v2,
                    SimdFloat *          v3)
{
    *v0 = simdLoad( base + align * offset[0] );
    *v1 = simdLoad( base + align * offset[1] );
    *v2 = simdLoad( base + align * offset[2] );
    *v3 = simdLoad( base + align * offset[3] );

    __vector float t0 = vec_mergeh(v0->simdInternal_, v2->simdInternal_);
    __vector float t1 = vec_mergel(v0->simdInternal_, v2->simdInternal_);
    __vector float t2 = vec_mergeh(v1->simdInternal_, v3->simdInternal_);
    __vector float t3 = vec_mergel(v1->simdInternal_, v3->simdInternal_);
    v0->simdInternal_ = vec_mergeh(t0, t2);
    v1->simdInternal_ = vec_mergel(t0, t2);
    v2->simdInternal_ = vec_mergeh(t1, t3);
    v3->simdInternal_ = vec_mergel(t1, t3);
}

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float *        base,
                    const std::int32_t   offset[],
                    SimdFloat *          v0,
                    SimdFloat *          v1)
{
    if (align % 4 == 0)
    {
        SimdFloat t2, t3;

        gatherLoadTranspose<align>(base, offset, v0, v1, &t2, &t3);
    }
    else
    {
        __vector float         t0, t1, t2, t3, t4, t5, t6, t7;
        __vector unsigned char p0, p1, p2, p3;

        // This is REALLY slow, since we have no choice but to load individual
        // elements when we cannot guarantee that we can access beyond the end of
        // the memory. Fortunately, 99% of the usage should be the aligned-to-4
        // case above instead.
        t0                = vec_lde(0, base + align * offset[0]);
        t1                = vec_lde(0, base + align * offset[1]);
        t2                = vec_lde(0, base + align * offset[2]);
        t3                = vec_lde(0, base + align * offset[3]);
        p0                = vec_lvsl(0, base + align * offset[0]);
        p1                = vec_lvsl(0, base + align * offset[1]);
        p2                = vec_lvsl(0, base + align * offset[2]);
        p3                = vec_lvsl(0, base + align * offset[3]);
        t0                = vec_perm(t0, t0, p0);
        t1                = vec_perm(t1, t1, p1);
        t2                = vec_perm(t2, t2, p2);
        t3                = vec_perm(t3, t3, p3);
        t0                = vec_mergeh(t0, t2);
        t1                = vec_mergeh(t1, t3);
        v0->simdInternal_ = vec_mergeh(t0, t1);

        t4                = vec_lde(0, base + align * offset[0] + 1);
        t5                = vec_lde(0, base + align * offset[1] + 1);
        t6                = vec_lde(0, base + align * offset[2] + 1);
        t7                = vec_lde(0, base + align * offset[3] + 1);
        p0                = vec_lvsl(0, base + align * offset[0] + 1);
        p1                = vec_lvsl(0, base + align * offset[1] + 1);
        p2                = vec_lvsl(0, base + align * offset[2] + 1);
        p3                = vec_lvsl(0, base + align * offset[3] + 1);
        t4                = vec_perm(t4, t4, p0);
        t5                = vec_perm(t5, t5, p1);
        t6                = vec_perm(t6, t6, p2);
        t7                = vec_perm(t7, t7, p3);
        t4                = vec_mergeh(t4, t6);
        t5                = vec_mergeh(t5, t7);
        v1->simdInternal_ = vec_mergeh(t4, t5);
    }
}

static const int c_simdBestPairAlignmentFloat = 2;

template <int align>
static inline void gmx_simdcall
gatherLoadUTranspose(const float *        base,
                     const std::int32_t   offset[],
                     SimdFloat *          v0,
                     SimdFloat *          v1,
                     SimdFloat *          v2)
{
    if (align % 4 == 0)
    {
        SimdFloat t3;
        gatherLoadTranspose<align>(base, offset, v0, v1, v2, &t3);
    }
    else
    {
        __vector float         t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
        __vector unsigned char p0, p1, p2, p3;

        // This is REALLY slow, since we have no choice but to load individual
        // elements when we cannot guarantee that we can access beyond the end of
        // the memory. Unfortunately this is likely the most common case.
        t0                = vec_lde(0, base + align * offset[0]);
        t1                = vec_lde(0, base + align * offset[1]);
        t2                = vec_lde(0, base + align * offset[2]);
        t3                = vec_lde(0, base + align * offset[3]);
        p0                = vec_lvsl(0, base + align * offset[0]);
        p1                = vec_lvsl(0, base + align * offset[1]);
        p2                = vec_lvsl(0, base + align * offset[2]);
        p3                = vec_lvsl(0, base + align * offset[3]);
        t0                = vec_perm(t0, t0, p0);
        t1                = vec_perm(t1, t1, p1);
        t2                = vec_perm(t2, t2, p2);
        t3                = vec_perm(t3, t3, p3);
        t0                = vec_mergeh(t0, t2);
        t1                = vec_mergeh(t1, t3);
        v0->simdInternal_ = vec_mergeh(t0, t1);

        t4                = vec_lde(0, base + align * offset[0] + 1);
        t5                = vec_lde(0, base + align * offset[1] + 1);
        t6                = vec_lde(0, base + align * offset[2] + 1);
        t7                = vec_lde(0, base + align * offset[3] + 1);
        p0                = vec_lvsl(0, base + align * offset[0] + 1);
        p1                = vec_lvsl(0, base + align * offset[1] + 1);
        p2                = vec_lvsl(0, base + align * offset[2] + 1);
        p3                = vec_lvsl(0, base + align * offset[3] + 1);
        t4                = vec_perm(t4, t4, p0);
        t5                = vec_perm(t5, t5, p1);
        t6                = vec_perm(t6, t6, p2);
        t7                = vec_perm(t7, t7, p3);
        t4                = vec_mergeh(t4, t6);
        t5                = vec_mergeh(t5, t7);
        v1->simdInternal_ = vec_mergeh(t4, t5);

        t8                = vec_lde(0, base + align * offset[0] + 2);
        t9                = vec_lde(0, base + align * offset[1] + 2);
        t10               = vec_lde(0, base + align * offset[2] + 2);
        t11               = vec_lde(0, base + align * offset[3] + 2);
        p0                = vec_lvsl(0, base + align * offset[0] + 2);
        p1                = vec_lvsl(0, base + align * offset[1] + 2);
        p2                = vec_lvsl(0, base + align * offset[2] + 2);
        p3                = vec_lvsl(0, base + align * offset[3] + 2);
        t8                = vec_perm(t8, t8, p0);
        t9                = vec_perm(t9, t9, p1);
        t10               = vec_perm(t10, t10, p2);
        t11               = vec_perm(t11, t11, p3);
        t8                = vec_mergeh(t8, t10);
        t9                = vec_mergeh(t9, t11);
        v2->simdInternal_ = vec_mergeh(t8, t9);
    }
}


template <int align>
static inline void gmx_simdcall
transposeScatterStoreU(float *              base,
                       const std::int32_t   offset[],
                       SimdFloat            v0,
                       SimdFloat            v1,
                       SimdFloat            v2)
{
    __vector unsigned char p0, p1, p2, p3;

    __vector float         t0 = vec_mergeh(v0.simdInternal_, v2.simdInternal_);
    __vector float         t1 = vec_mergel(v0.simdInternal_, v2.simdInternal_);
    __vector float         t2 = vec_mergeh(v1.simdInternal_, v2.simdInternal_);
    __vector float         t3 = vec_mergel(v1.simdInternal_, v2.simdInternal_);
    __vector float         t4 = vec_mergeh(t0, t2);
    __vector float         t5 = vec_mergel(t0, t2);
    __vector float         t6 = vec_mergeh(t1, t3);
    __vector float         t7 = vec_mergel(t1, t3);

    p0 = vec_lvsr(0, base + align * offset[0]);
    p1 = vec_lvsr(0, base + align * offset[1]);
    p2 = vec_lvsr(0, base + align * offset[2]);
    p3 = vec_lvsr(0, base + align * offset[3]);

    t4 = vec_perm(t4, t4, p0);
    t5 = vec_perm(t5, t5, p1);
    t6 = vec_perm(t6, t6, p2);
    t7 = vec_perm(t7, t7, p3);

    vec_ste(t4, 0, base + align * offset[0]);
    vec_ste(t4, 4, base + align * offset[0]);
    vec_ste(t4, 8, base + align * offset[0]);
    vec_ste(t5, 0, base + align * offset[1]);
    vec_ste(t5, 4, base + align * offset[1]);
    vec_ste(t5, 8, base + align * offset[1]);
    vec_ste(t6, 0, base + align * offset[2]);
    vec_ste(t6, 4, base + align * offset[2]);
    vec_ste(t6, 8, base + align * offset[2]);
    vec_ste(t7, 0, base + align * offset[3]);
    vec_ste(t7, 4, base + align * offset[3]);
    vec_ste(t7, 8, base + align * offset[3]);
}


template <int align>
static inline void gmx_simdcall
transposeScatterIncrU(float *              base,
                      const std::int32_t   offset[],
                      SimdFloat            v0,
                      SimdFloat            v1,
                      SimdFloat            v2)
{
    if (align % 4 == 0)
    {
        __vector float zero = reinterpret_cast<__vector float>(vec_splat_u32(0));
        __vector float t0   = vec_mergeh(v0.simdInternal_, v2.simdInternal_);
        __vector float t1   = vec_mergel(v0.simdInternal_, v2.simdInternal_);
        __vector float t2   = vec_mergeh(v1.simdInternal_, zero);
        __vector float t3   = vec_mergel(v1.simdInternal_, zero);
        __vector float t4   = vec_mergeh(t0, t2);
        __vector float t5   = vec_mergel(t0, t2);
        __vector float t6   = vec_mergeh(t1, t3);
        __vector float t7   = vec_mergel(t1, t3);

        vec_st( vec_add( vec_ld(0, base + align * offset[0]), t4), 0, base + align * offset[0]);
        vec_st( vec_add( vec_ld(0, base + align * offset[1]), t5), 0, base + align * offset[1]);
        vec_st( vec_add( vec_ld(0, base + align * offset[2]), t6), 0, base + align * offset[2]);
        vec_st( vec_add( vec_ld(0, base + align * offset[3]), t7), 0, base + align * offset[3]);
    }
    else
    {
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH) rdata0[GMX_SIMD_FLOAT_WIDTH];
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH) rdata1[GMX_SIMD_FLOAT_WIDTH];
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH) rdata2[GMX_SIMD_FLOAT_WIDTH];

        vec_st(v0.simdInternal_, 0, rdata0);
        vec_st(v1.simdInternal_, 0, rdata1);
        vec_st(v2.simdInternal_, 0, rdata2);

        base[align*offset[0] + 0] += rdata0[0];
        base[align*offset[0] + 1] += rdata1[0];
        base[align*offset[0] + 2] += rdata2[0];
        base[align*offset[1] + 0] += rdata0[1];
        base[align*offset[1] + 1] += rdata1[1];
        base[align*offset[1] + 2] += rdata2[1];
        base[align*offset[2] + 0] += rdata0[2];
        base[align*offset[2] + 1] += rdata1[2];
        base[align*offset[2] + 2] += rdata2[2];
        base[align*offset[3] + 0] += rdata0[3];
        base[align*offset[3] + 1] += rdata1[3];
        base[align*offset[3] + 2] += rdata2[3];
    }
}

template <int align>
static inline void gmx_simdcall
transposeScatterDecrU(float *              base,
                      const std::int32_t   offset[],
                      SimdFloat            v0,
                      SimdFloat            v1,
                      SimdFloat            v2)
{
    if (align % 4 == 0)
    {
        __vector float zero = reinterpret_cast<__vector float>(vec_splat_u32(0));
        __vector float t0   = vec_mergeh(v0.simdInternal_, v2.simdInternal_);
        __vector float t1   = vec_mergel(v0.simdInternal_, v2.simdInternal_);
        __vector float t2   = vec_mergeh(v1.simdInternal_, zero);
        __vector float t3   = vec_mergel(v1.simdInternal_, zero);
        __vector float t4   = vec_mergeh(t0, t2);
        __vector float t5   = vec_mergel(t0, t2);
        __vector float t6   = vec_mergeh(t1, t3);
        __vector float t7   = vec_mergel(t1, t3);

        vec_st( vec_sub( vec_ld(0, base + align * offset[0]), t4), 0, base + align * offset[0]);
        vec_st( vec_sub( vec_ld(0, base + align * offset[1]), t5), 0, base + align * offset[1]);
        vec_st( vec_sub( vec_ld(0, base + align * offset[2]), t6), 0, base + align * offset[2]);
        vec_st( vec_sub( vec_ld(0, base + align * offset[3]), t7), 0, base + align * offset[3]);
    }
    else
    {
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH) rdata0[GMX_SIMD_FLOAT_WIDTH];
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH) rdata1[GMX_SIMD_FLOAT_WIDTH];
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH) rdata2[GMX_SIMD_FLOAT_WIDTH];

        vec_st(v0.simdInternal_, 0, rdata0);
        vec_st(v1.simdInternal_, 0, rdata1);
        vec_st(v2.simdInternal_, 0, rdata2);

        base[align*offset[0] + 0] -= rdata0[0];
        base[align*offset[0] + 1] -= rdata1[0];
        base[align*offset[0] + 2] -= rdata2[0];
        base[align*offset[1] + 0] -= rdata0[1];
        base[align*offset[1] + 1] -= rdata1[1];
        base[align*offset[1] + 2] -= rdata2[1];
        base[align*offset[2] + 0] -= rdata0[2];
        base[align*offset[2] + 1] -= rdata1[2];
        base[align*offset[2] + 2] -= rdata2[2];
        base[align*offset[3] + 0] -= rdata0[3];
        base[align*offset[3] + 1] -= rdata1[3];
        base[align*offset[3] + 2] -= rdata2[3];
    }
}

static inline void gmx_simdcall
expandScalarsToTriplets(SimdFloat    scalar,
                        SimdFloat *  triplets0,
                        SimdFloat *  triplets1,
                        SimdFloat *  triplets2)
{
    const __vector unsigned char perm0 = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7 };
    const __vector unsigned char perm1 = { 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 8, 9, 10, 11 };
    const __vector unsigned char perm2 = { 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15 };

    triplets0->simdInternal_ = vec_perm(scalar.simdInternal_, scalar.simdInternal_, perm0);
    triplets1->simdInternal_ = vec_perm(scalar.simdInternal_, scalar.simdInternal_, perm1);
    triplets2->simdInternal_ = vec_perm(scalar.simdInternal_, scalar.simdInternal_, perm2);
}


template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *  base,
                             SimdFInt32     offset,
                             SimdFloat *    v0,
                             SimdFloat *    v1,
                             SimdFloat *    v2,
                             SimdFloat *    v3)
{
    GMX_ALIGNED(int, GMX_SIMD_FINT32_WIDTH) ioffset[GMX_SIMD_FINT32_WIDTH];

    vec_st( offset.simdInternal_, 0, ioffset);
    gatherLoadTranspose<align>(base, ioffset, v0, v1, v2, v3);
}

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *   base,
                             SimdFInt32      offset,
                             SimdFloat *     v0,
                             SimdFloat *     v1)
{
    GMX_ALIGNED(int, GMX_SIMD_FINT32_WIDTH) ioffset[GMX_SIMD_FINT32_WIDTH];

    vec_st( offset.simdInternal_, 0, ioffset);
    gatherLoadTranspose<align>(base, ioffset, v0, v1);
}



static inline float gmx_simdcall
reduceIncr4ReturnSum(float *    m,
                     SimdFloat  v0,
                     SimdFloat  v1,
                     SimdFloat  v2,
                     SimdFloat  v3)
{
    __vector float t0 = vec_mergeh(v0.simdInternal_, v2.simdInternal_);
    __vector float t1 = vec_mergel(v0.simdInternal_, v2.simdInternal_);
    __vector float t2 = vec_mergeh(v1.simdInternal_, v3.simdInternal_);
    __vector float t3 = vec_mergel(v1.simdInternal_, v3.simdInternal_);
    v0.simdInternal_ = vec_mergeh(t0, t2);
    v1.simdInternal_ = vec_mergel(t0, t2);
    v2.simdInternal_ = vec_mergeh(t1, t3);
    v3.simdInternal_ = vec_mergel(t1, t3);

    v0 = v0 + v1;
    v2 = v2 + v3;
    v0 = v0 + v2;
    v2 = v0 + simdLoad(m);
    store(m, v2);

    return reduce(v0);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_IBM_VMX_UTIL_FLOAT_H
