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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_FLOAT_H

#include "config.h"

#include "gromacs/utility/basedefinitions.h"

#include "impl_ibm_vsx_definitions.h"
#include "impl_ibm_vsx_simd_float.h"

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
    __vector float t0, t1, t2, t3;

    t0                = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<const double *>(base + align * offset[0])));
    t1                = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<const double *>(base + align * offset[1])));
    t2                = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<const double *>(base + align * offset[2])));
    t3                = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<const double *>(base + align * offset[3])));
    t0                = vec_mergeh(t0, t2);
    t1                = vec_mergeh(t1, t3);
    v0->simdInternal_ = vec_mergeh(t0, t1);
    v1->simdInternal_ = vec_mergel(t0, t1);
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
        __vector float               t1, t2, t3, t4, t5, t6, t7, t8;
        const __vector unsigned char perm_lo2hi = { 0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23 };
        const __vector unsigned char perm_hi2lo = { 24, 25, 26, 27, 28, 29, 30, 31, 8, 9, 10, 11, 12, 13, 14, 15 };

        t1  = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<const double *>(base + align * offset[0])));
        t2  = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<const double *>(base + align * offset[1])));
        t3  = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<const double *>(base + align * offset[2])));
        t4  = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<const double *>(base + align * offset[3])));
        t5  = vec_splats( *(base + align * offset[0] + 2) );
        t6  = vec_splats( *(base + align * offset[1] + 2) );
        t7  = vec_splats( *(base + align * offset[2] + 2) );
        t8  = vec_splats( *(base + align * offset[3] + 2) );

        t1                = vec_mergeh(t1, t2);
        t3                = vec_mergeh(t3, t4);
        v0->simdInternal_ = vec_perm(t1, t3, perm_lo2hi);
        v1->simdInternal_ = vec_perm(t3, t1, perm_hi2lo);
        t5                = vec_mergeh(t5, t6);
        t7                = vec_mergeh(t7, t8);
        v2->simdInternal_ = vec_perm(t5, t7, perm_lo2hi);
    }
}


// gcc-4.9 does not recognize that the argument to vec_extract() is used
template <int align>
static inline void gmx_simdcall
transposeScatterStoreU(float *              base,
                       const std::int32_t   offset[],
                       SimdFloat            v0,
                       SimdFloat            v1,
                       SimdFloat gmx_unused v2)
{
    __vector float t1, t2;

    t1   = vec_mergeh(v0.simdInternal_, v1.simdInternal_);
    t2   = vec_mergel(v0.simdInternal_, v1.simdInternal_);
    *reinterpret_cast<double *>( base + align * offset[0] ) = vec_extract(reinterpret_cast<__vector double>(t1), 0);
    base[align*offset[0] + 2]                               = vec_extract(v2.simdInternal_, 0);
    *reinterpret_cast<double *>( base + align * offset[1] ) = vec_extract(reinterpret_cast<__vector double>(t1), 1);
    base[align*offset[1] + 2]                               = vec_extract(v2.simdInternal_, 1);
    *reinterpret_cast<double *>( base + align * offset[2] ) = vec_extract(reinterpret_cast<__vector double>(t2), 0);
    base[align*offset[2] + 2]                               = vec_extract(v2.simdInternal_, 2);
    *reinterpret_cast<double *>( base + align * offset[3] ) = vec_extract(reinterpret_cast<__vector double>(t2), 1);
    base[align*offset[3] + 2]                               = vec_extract(v2.simdInternal_, 3);
}

template <int align>
static inline void gmx_simdcall
transposeScatterIncrU(float *              base,
                      const std::int32_t   offset[],
                      SimdFloat            v0,
                      SimdFloat            v1,
                      SimdFloat            v2)
{
    if (align < 4)
    {
        const __vector unsigned char perm_hi2lo = { 24, 25, 26, 27, 28, 29, 30, 31, 8, 9, 10, 11, 12, 13, 14, 15 };
        __vector float               t0, t1, t2, t3, t4, t5, t6, t7;

        t0 = vec_mergeh(v0.simdInternal_, v1.simdInternal_); // x0 y0 x1 y1
        t1 = vec_perm(t0, t0, perm_hi2lo);                   // x1 y1
        t2 = vec_mergel(v0.simdInternal_, v1.simdInternal_); // x2 y2 x3 y3
        t3 = vec_perm(t2, t2, perm_hi2lo);                   // x3 y3

        t4 = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<double *>(base + align * offset[0])));
        t4 = vec_add(t4, t0);
        *reinterpret_cast<double *>( base + align * offset[0] ) = vec_extract(reinterpret_cast<__vector double>(t4), 0);
        {
            float extracted = vec_extract(v2.simdInternal_, 0);
            base[align*offset[0] + 2] += extracted;
        }

        t5 = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<double *>(base + align * offset[1])));
        t5 = vec_add(t5, t1);
        *reinterpret_cast<double *>( base + align * offset[1] ) = vec_extract(reinterpret_cast<__vector double>(t5), 0);
        {
            float extracted = vec_extract(v2.simdInternal_, 1);
            base[align*offset[1] + 2] += extracted;
        }

        t6 = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<double *>(base + align * offset[2])));
        t6 = vec_add(t6, t2);
        *reinterpret_cast<double *>( base + align * offset[2] ) = vec_extract(reinterpret_cast<__vector double>(t6), 0);
        {
            float extracted = vec_extract(v2.simdInternal_, 2);
            base[align*offset[2] + 2] += extracted;
        }

        t7 = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<double *>(base + align * offset[3])));
        t7 = vec_add(t7, t3);
        *reinterpret_cast<double *>( base + align * offset[3] ) = vec_extract(reinterpret_cast<__vector double>(t7), 0);
        {
            float extracted = vec_extract(v2.simdInternal_, 3);
            base[align*offset[3] + 2] += extracted;
        }
    }
    else
    {
        // Extra elements means we can use full width-4 load/store operations
        SimdFloat      v3;
        __vector float t0 = vec_mergeh(v0.simdInternal_, v2.simdInternal_);
        __vector float t1 = vec_mergel(v0.simdInternal_, v2.simdInternal_);
        __vector float t2 = vec_mergeh(v1.simdInternal_, vec_splats(0.0f));
        __vector float t3 = vec_mergel(v1.simdInternal_, vec_splats(0.0f));
        v0.simdInternal_ = vec_mergeh(t0, t2);
        v1.simdInternal_ = vec_mergel(t0, t2);
        v2.simdInternal_ = vec_mergeh(t1, t3);
        v3.simdInternal_ = vec_mergel(t1, t3);

        store(base + align * offset[0], simdLoad(base + align * offset[0]) + v0 );
        store(base + align * offset[1], simdLoad(base + align * offset[1]) + v1 );
        store(base + align * offset[2], simdLoad(base + align * offset[2]) + v2 );
        store(base + align * offset[3], simdLoad(base + align * offset[3]) + v3 );
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
    if (align < 4)
    {
        const __vector unsigned char perm_hi2lo = { 24, 25, 26, 27, 28, 29, 30, 31, 8, 9, 10, 11, 12, 13, 14, 15 };
        __vector float               t0, t1, t2, t3, t4, t5, t6, t7;

        t0 = vec_mergeh(v0.simdInternal_, v1.simdInternal_); // x0 y0 x1 y1
        t1 = vec_perm(t0, t0, perm_hi2lo);                   // x1 y1
        t2 = vec_mergel(v0.simdInternal_, v1.simdInternal_); // x2 y2 x3 y3
        t3 = vec_perm(t2, t2, perm_hi2lo);                   // x3 y3

        t4 = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<double *>(base + align * offset[0])));
        t4 = vec_sub(t4, t0);
        *reinterpret_cast<double *>( base + align * offset[0] ) = vec_extract(reinterpret_cast<__vector double>(t4), 0);
        {
            float extracted = vec_extract(v2.simdInternal_, 0);
            base[align*offset[0] + 2] -= extracted;
        }

        t5 = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<double *>(base + align * offset[1])));
        t5 = vec_sub(t5, t1);
        *reinterpret_cast<double *>( base + align * offset[1] ) = vec_extract(reinterpret_cast<__vector double>(t5), 0);
        {
            float extracted = vec_extract(v2.simdInternal_, 1);
            base[align*offset[1] + 2] -= extracted;
        }

        t6 = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<double *>(base + align * offset[2])));
        t6 = vec_sub(t6, t2);
        *reinterpret_cast<double *>( base + align * offset[2] ) = vec_extract(reinterpret_cast<__vector double>(t6), 0);
        {
            float extracted = vec_extract(v2.simdInternal_, 2);
            base[align*offset[2] + 2] -= extracted;
        }

        t7 = reinterpret_cast<__vector float>(vec_splats(*reinterpret_cast<double *>(base + align * offset[3])));
        t7 = vec_sub(t7, t3);
        *reinterpret_cast<double *>( base + align * offset[3] ) = vec_extract(reinterpret_cast<__vector double>(t7), 0);
        {
            float extracted = vec_extract(v2.simdInternal_, 3);
            base[align*offset[3] + 2] -= extracted;
        }
    }
    else
    {
        // Extra elements means we can use full width-4 load/store operations
        SimdFloat      v3;
        __vector float t0 = vec_mergeh(v0.simdInternal_, v2.simdInternal_);
        __vector float t1 = vec_mergel(v0.simdInternal_, v2.simdInternal_);
        __vector float t2 = vec_mergeh(v1.simdInternal_, vec_splats(0.0f));
        __vector float t3 = vec_mergel(v1.simdInternal_, vec_splats(0.0f));
        v0.simdInternal_ = vec_mergeh(t0, t2);
        v1.simdInternal_ = vec_mergel(t0, t2);
        v2.simdInternal_ = vec_mergeh(t1, t3);
        v3.simdInternal_ = vec_mergel(t1, t3);

        store(base + align * offset[0], simdLoad(base + align * offset[0]) - v0 );
        store(base + align * offset[1], simdLoad(base + align * offset[1]) - v1 );
        store(base + align * offset[2], simdLoad(base + align * offset[2]) - v2 );
        store(base + align * offset[3], simdLoad(base + align * offset[3]) - v3 );
    }
}

static inline void gmx_simdcall
expandScalarsToTriplets(SimdFloat    scalar,
                        SimdFloat *  triplets0,
                        SimdFloat *  triplets1,
                        SimdFloat *  triplets2)
{
    // These permutes will be translated to immediate permutes (xxpermdi)
    // since they operate on doublewords, which will be faster than loading
    // the constants required for fully flexible permutes.
    // (although the real reason was that the latter was buggy on xlc-13.1).
    __vector unsigned char perm0 = { 0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23 };
    __vector unsigned char perm1 = { 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
    __vector unsigned char perm2 = { 8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31 };
    __vector float         t0, t1;

    t0                       = vec_mergeh(scalar.simdInternal_, scalar.simdInternal_);
    t1                       = vec_mergel(scalar.simdInternal_, scalar.simdInternal_);
    triplets0->simdInternal_ = vec_perm(t0, scalar.simdInternal_, perm0);
    triplets1->simdInternal_ = vec_perm(t0, t1, perm1);
    triplets2->simdInternal_ = vec_perm(scalar.simdInternal_, t1, perm2);
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
    GMX_ALIGNED(std::int32_t, GMX_SIMD_FINT32_WIDTH) ioffset[GMX_SIMD_FINT32_WIDTH];

    store(ioffset, offset );
    gatherLoadTranspose<align>(base, ioffset, v0, v1, v2, v3);
}

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *   base,
                             SimdFInt32      offset,
                             SimdFloat *     v0,
                             SimdFloat *     v1)
{
    GMX_ALIGNED(std::int32_t, GMX_SIMD_FINT32_WIDTH) ioffset[GMX_SIMD_FINT32_WIDTH];

    store(ioffset, offset );
    gatherLoadTranspose<align>(base, ioffset, v0, v1);
}

template <int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const float *  base,
                              SimdFInt32     offset,
                              SimdFloat *    v0,
                              SimdFloat *    v1)
{
    GMX_ALIGNED(std::int32_t, GMX_SIMD_FINT32_WIDTH) ioffset[GMX_SIMD_FINT32_WIDTH];

    store(ioffset, offset );
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

#endif // GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_FLOAT_H
