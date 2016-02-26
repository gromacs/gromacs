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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_UTIL_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_UTIL_FLOAT_H

#include "config.h"

// Assert is buggy on xlc with high optimization, so we skip it for QPX
#include <cstddef>
#include <cstdint>

#ifdef __clang__
#    include <qpxmath.h>
#endif

#include "gromacs/utility/basedefinitions.h"

#include "impl_ibm_qpx_simd_float.h"

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
    v0->simdInternal_ = vec_ld(0, const_cast<float *>(base + align * offset[0]) );
    v1->simdInternal_ = vec_ld(0, const_cast<float *>(base + align * offset[1]) );
    v2->simdInternal_ = vec_ld(0, const_cast<float *>(base + align * offset[2]) );
    v3->simdInternal_ = vec_ld(0, const_cast<float *>(base + align * offset[3]) );

    vector4double t0 = vec_perm(v0->simdInternal_, v2->simdInternal_, vec_gpci(00415));
    vector4double t1 = vec_perm(v0->simdInternal_, v2->simdInternal_, vec_gpci(02637));
    vector4double t2 = vec_perm(v1->simdInternal_, v3->simdInternal_, vec_gpci(00415));
    vector4double t3 = vec_perm(v1->simdInternal_, v3->simdInternal_, vec_gpci(02637));
    v0->simdInternal_ = vec_perm(t0, t2, vec_gpci(00415));
    v1->simdInternal_ = vec_perm(t0, t2, vec_gpci(02637));
    v2->simdInternal_ = vec_perm(t1, t3, vec_gpci(00415));
    v3->simdInternal_ = vec_perm(t1, t3, vec_gpci(02637));
}

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float *        base,
                    const std::int32_t   offset[],
                    SimdFloat *          v0,
                    SimdFloat *          v1)
{
    vector4double t0, t1, t2, t3;

    t0                = vec_ld2(0, const_cast<float *>(base + align * offset[0]) );
    t1                = vec_ld2(0, const_cast<float *>(base + align * offset[1]) );
    t2                = vec_ld2(0, const_cast<float *>(base + align * offset[2]) );
    t3                = vec_ld2(0, const_cast<float *>(base + align * offset[3]) );
    t0                = vec_perm(t0, t2, vec_gpci(00415));
    t1                = vec_perm(t1, t3, vec_gpci(00415));
    v0->simdInternal_ = vec_perm(t0, t1, vec_gpci(00415));
    v1->simdInternal_ = vec_perm(t0, t1, vec_gpci(02637));
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
    vector4double             t1, t2, t3, t4, t5, t6, t7, t8;

    if (align % 4 == 0)
    {
        SimdFloat t3;
        gatherLoadTranspose<align>(base, offset, v0, v1, v2, &t3);
    }
    else
    {
        t1  = vec_perm(vec_splats(base[align * offset[0]]), vec_splats(base[align * offset[0] + 1]), vec_gpci(00415));
        t2  = vec_perm(vec_splats(base[align * offset[1]]), vec_splats(base[align * offset[1] + 1]), vec_gpci(00415));
        t3  = vec_perm(vec_splats(base[align * offset[2]]), vec_splats(base[align * offset[2] + 1]), vec_gpci(00415));
        t4  = vec_perm(vec_splats(base[align * offset[3]]), vec_splats(base[align * offset[3] + 1]), vec_gpci(00415));

        t5  = vec_splats( *(base + align * offset[0] + 2) );
        t6  = vec_splats( *(base + align * offset[1] + 2) );
        t7  = vec_splats( *(base + align * offset[2] + 2) );
        t8  = vec_splats( *(base + align * offset[3] + 2) );

        t1                = vec_perm(t1, t2, vec_gpci(00415));
        t3                = vec_perm(t3, t4, vec_gpci(00415));
        v0->simdInternal_ = vec_perm(t1, t3, vec_gpci(00145));
        v1->simdInternal_ = vec_perm(t3, t1, vec_gpci(06723));
        t5                = vec_perm(t5, t6, vec_gpci(00415));
        t7                = vec_perm(t7, t8, vec_gpci(00415));
        v2->simdInternal_ = vec_perm(t5, t7, vec_gpci(00145));
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
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   m0[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   m1[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   m2[GMX_SIMD_FLOAT_WIDTH];

    store(m0, v0);
    store(m1, v1);
    store(m2, v2);

    base[align * offset[0]    ] = m0[0];
    base[align * offset[0] + 1] = m1[0];
    base[align * offset[0] + 2] = m2[0];
    base[align * offset[1]    ] = m0[1];
    base[align * offset[1] + 1] = m1[1];
    base[align * offset[1] + 2] = m2[1];
    base[align * offset[2]    ] = m0[2];
    base[align * offset[2] + 1] = m1[2];
    base[align * offset[2] + 2] = m2[2];
    base[align * offset[3]    ] = m0[3];
    base[align * offset[3] + 1] = m1[3];
    base[align * offset[3] + 2] = m2[3];
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
        // transpose
        SimdFloat     v3(0.0);
        vector4double t0 = vec_perm(v0.simdInternal_, v2.simdInternal_, vec_gpci(00415));
        vector4double t1 = vec_perm(v0.simdInternal_, v2.simdInternal_, vec_gpci(02637));
        vector4double t2 = vec_perm(v1.simdInternal_, v3.simdInternal_, vec_gpci(00415));
        vector4double t3 = vec_perm(v1.simdInternal_, v3.simdInternal_, vec_gpci(02637));
        v0.simdInternal_ = vec_perm(t0, t2, vec_gpci(00415));
        v1.simdInternal_ = vec_perm(t0, t2, vec_gpci(02637));
        v2.simdInternal_ = vec_perm(t1, t3, vec_gpci(00415));
        v3.simdInternal_ = vec_perm(t1, t3, vec_gpci(02637));
        // increment
        store(base + align*offset[0], simdLoad(base + align*offset[0]) + v0);
        store(base + align*offset[1], simdLoad(base + align*offset[1]) + v1);
        store(base + align*offset[2], simdLoad(base + align*offset[2]) + v2);
        store(base + align*offset[3], simdLoad(base + align*offset[3]) + v3);
    }
    else
    {
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   m0[GMX_SIMD_FLOAT_WIDTH];
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   m1[GMX_SIMD_FLOAT_WIDTH];
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   m2[GMX_SIMD_FLOAT_WIDTH];

        store(m0, v0);
        store(m1, v1);
        store(m2, v2);

        base[align * offset[0]    ] += m0[0];
        base[align * offset[0] + 1] += m1[0];
        base[align * offset[0] + 2] += m2[0];
        base[align * offset[1]    ] += m0[1];
        base[align * offset[1] + 1] += m1[1];
        base[align * offset[1] + 2] += m2[1];
        base[align * offset[2]    ] += m0[2];
        base[align * offset[2] + 1] += m1[2];
        base[align * offset[2] + 2] += m2[2];
        base[align * offset[3]    ] += m0[3];
        base[align * offset[3] + 1] += m1[3];
        base[align * offset[3] + 2] += m2[3];
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
        // transpose
        SimdFloat     v3(0.0);
        vector4double t0 = vec_perm(v0.simdInternal_, v2.simdInternal_, vec_gpci(00415));
        vector4double t1 = vec_perm(v0.simdInternal_, v2.simdInternal_, vec_gpci(02637));
        vector4double t2 = vec_perm(v1.simdInternal_, v3.simdInternal_, vec_gpci(00415));
        vector4double t3 = vec_perm(v1.simdInternal_, v3.simdInternal_, vec_gpci(02637));
        v0.simdInternal_ = vec_perm(t0, t2, vec_gpci(00415));
        v1.simdInternal_ = vec_perm(t0, t2, vec_gpci(02637));
        v2.simdInternal_ = vec_perm(t1, t3, vec_gpci(00415));
        v3.simdInternal_ = vec_perm(t1, t3, vec_gpci(02637));
        // decrement
        store(base + align*offset[0], simdLoad(base + align*offset[0]) - v0);
        store(base + align*offset[1], simdLoad(base + align*offset[1]) - v1);
        store(base + align*offset[2], simdLoad(base + align*offset[2]) - v2);
        store(base + align*offset[3], simdLoad(base + align*offset[3]) - v3);
    }
    else
    {
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   m0[GMX_SIMD_FLOAT_WIDTH];
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   m1[GMX_SIMD_FLOAT_WIDTH];
        GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   m2[GMX_SIMD_FLOAT_WIDTH];

        store(m0, v0);
        store(m1, v1);
        store(m2, v2);

        base[align * offset[0]    ] -= m0[0];
        base[align * offset[0] + 1] -= m1[0];
        base[align * offset[0] + 2] -= m2[0];
        base[align * offset[1]    ] -= m0[1];
        base[align * offset[1] + 1] -= m1[1];
        base[align * offset[1] + 2] -= m2[1];
        base[align * offset[2]    ] -= m0[2];
        base[align * offset[2] + 1] -= m1[2];
        base[align * offset[2] + 2] -= m2[2];
        base[align * offset[3]    ] -= m0[3];
        base[align * offset[3] + 1] -= m1[3];
        base[align * offset[3] + 2] -= m2[3];
    }
}

static inline void gmx_simdcall
expandScalarsToTriplets(SimdFloat    scalar,
                        SimdFloat *  triplets0,
                        SimdFloat *  triplets1,
                        SimdFloat *  triplets2)
{
    triplets0->simdInternal_ = vec_perm(scalar.simdInternal_, scalar.simdInternal_, vec_gpci(00001));
    triplets1->simdInternal_ = vec_perm(scalar.simdInternal_, scalar.simdInternal_, vec_gpci(01122));
    triplets2->simdInternal_ = vec_perm(scalar.simdInternal_, scalar.simdInternal_, vec_gpci(02333));
}

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *  base,
                             SimdFInt32     simdoffset,
                             SimdFloat *    v0,
                             SimdFloat *    v1,
                             SimdFloat *    v2,
                             SimdFloat *    v3)
{
    GMX_ALIGNED(int, GMX_SIMD_FLOAT_WIDTH)   ioffset[GMX_SIMD_FLOAT_WIDTH];

    store(ioffset, simdoffset);
    gatherLoadTranspose<align>(base, ioffset, v0, v1, v2, v3);
}

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *   base,
                             SimdFInt32      simdoffset,
                             SimdFloat *     v0,
                             SimdFloat *     v1)
{
    GMX_ALIGNED(int, GMX_SIMD_FLOAT_WIDTH)   ioffset[GMX_SIMD_FLOAT_WIDTH];

    store(ioffset, simdoffset);
    gatherLoadTranspose<align>(base, ioffset, v0, v1);
}

static inline float gmx_simdcall
reduceIncr4ReturnSum(float *    m,
                     SimdFloat  v0,
                     SimdFloat  v1,
                     SimdFloat  v2,
                     SimdFloat  v3)
{
    vector4double t0 = vec_perm(v0.simdInternal_, v2.simdInternal_, vec_gpci(00415));
    vector4double t1 = vec_perm(v0.simdInternal_, v2.simdInternal_, vec_gpci(02637));
    vector4double t2 = vec_perm(v1.simdInternal_, v3.simdInternal_, vec_gpci(00415));
    vector4double t3 = vec_perm(v1.simdInternal_, v3.simdInternal_, vec_gpci(02637));
    v0.simdInternal_ = vec_perm(t0, t2, vec_gpci(00415));
    v1.simdInternal_ = vec_perm(t0, t2, vec_gpci(02637));
    v2.simdInternal_ = vec_perm(t1, t3, vec_gpci(00415));
    v3.simdInternal_ = vec_perm(t1, t3, vec_gpci(02637));

    v0 = v0 + v1;
    v2 = v2 + v3;
    v0 = v0 + v2;
    v2 = v0 + simdLoad(m);
    store(m, v2);

    return reduce(v0);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPLEMENTATION_IBM_QPX_UTIL_FLOAT_H
