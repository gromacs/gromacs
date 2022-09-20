/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_DOUBLE_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_DOUBLE_H

#include "config.h"

#include "gromacs/simd/impl_ibm_vsx/impl_ibm_vsx_definitions.h"
#include "gromacs/simd/impl_ibm_vsx/impl_ibm_vsx_simd_double.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

template<int align>
static inline void gmx_simdcall gatherLoadTranspose(const double*      base,
                                                    const std::int32_t offset[],
                                                    SimdDouble*        v0,
                                                    SimdDouble*        v1,
                                                    SimdDouble*        v2,
                                                    SimdDouble*        v3)
{
    __vector double t1, t2, t3, t4;

    t1                = *reinterpret_cast<const __vector double*>(base + align * offset[0]);
    t2                = *reinterpret_cast<const __vector double*>(base + align * offset[1]);
    t3                = *reinterpret_cast<const __vector double*>(base + align * offset[0] + 2);
    t4                = *reinterpret_cast<const __vector double*>(base + align * offset[1] + 2);
    v0->simdInternal_ = vec_mergeh(t1, t2);
    v1->simdInternal_ = vec_mergel(t1, t2);
    v2->simdInternal_ = vec_mergeh(t3, t4);
    v3->simdInternal_ = vec_mergel(t3, t4);
}

template<int align>
static inline void gmx_simdcall
gatherLoadTranspose(const double* base, const std::int32_t offset[], SimdDouble* v0, SimdDouble* v1)
{
    __vector double t1, t2;

    t1                = *reinterpret_cast<const __vector double*>(base + align * offset[0]);
    t2                = *reinterpret_cast<const __vector double*>(base + align * offset[1]);
    v0->simdInternal_ = vec_mergeh(t1, t2);
    v1->simdInternal_ = vec_mergel(t1, t2);
}

static const int c_simdBestPairAlignmentDouble = 2;

template<int align>
static inline void gmx_simdcall gatherLoadUTranspose(const double*      base,
                                                     const std::int32_t offset[],
                                                     SimdDouble*        v0,
                                                     SimdDouble*        v1,
                                                     SimdDouble*        v2)
{
    SimdDouble t1, t2;

    t1 = simdLoad(base + align * offset[0]);
    t2 = simdLoad(base + align * offset[1]);

    v0->simdInternal_ = vec_mergeh(t1.simdInternal_, t2.simdInternal_);
    v1->simdInternal_ = vec_mergel(t1.simdInternal_, t2.simdInternal_);
    v2->simdInternal_ = vec_mergeh(vec_splats(*(base + align * offset[0] + 2)),
                                   vec_splats(*(base + align * offset[1] + 2)));
}

template<int align>
static inline void gmx_simdcall transposeScatterStoreU(double*            base,
                                                       const std::int32_t offset[],
                                                       SimdDouble         v0,
                                                       SimdDouble         v1,
                                                       SimdDouble         v2)
{
    SimdDouble t1, t2;

    t1.simdInternal_ = vec_mergeh(v0.simdInternal_, v1.simdInternal_);
    t2.simdInternal_ = vec_mergel(v0.simdInternal_, v1.simdInternal_);

    store(base + align * offset[0], t1);
    base[align * offset[0] + 2] = vec_extract(v2.simdInternal_, 0);
    store(base + align * offset[1], t2);
    base[align * offset[1] + 2] = vec_extract(v2.simdInternal_, 1);
}

template<int align>
static inline void gmx_simdcall
transposeScatterIncrU(double* base, const std::int32_t offset[], SimdDouble v0, SimdDouble v1, SimdDouble v2)
{
    if (align % 4 == 0)
    {
        __vector double t1, t2, t3, t4;
        SimdDouble      t5, t6, t7, t8;

        t1 = vec_mergeh(v0.simdInternal_, v1.simdInternal_);
        t2 = vec_mergel(v0.simdInternal_, v1.simdInternal_);
        t3 = vec_mergeh(v2.simdInternal_, vec_splats(0.0));
        t4 = vec_mergel(v2.simdInternal_, vec_splats(0.0));

        t5               = simdLoad(base + align * offset[0]);
        t6               = simdLoad(base + align * offset[0] + 2);
        t5.simdInternal_ = vec_add(t5.simdInternal_, t1);
        t6.simdInternal_ = vec_add(t6.simdInternal_, t3);
        store(base + align * offset[0], t5);
        store(base + align * offset[0] + 2, t6);

        t5               = simdLoad(base + align * offset[1]);
        t6               = simdLoad(base + align * offset[1] + 2);
        t5.simdInternal_ = vec_add(t5.simdInternal_, t2);
        t6.simdInternal_ = vec_add(t6.simdInternal_, t4);
        store(base + align * offset[1], t5);
        store(base + align * offset[1] + 2, t6);
    }
    else
    {
        __vector double t1, t2;
        SimdDouble      t3, t4;

        t1 = vec_mergeh(v0.simdInternal_, v1.simdInternal_);
        t2 = vec_mergel(v0.simdInternal_, v1.simdInternal_);

        t3               = simdLoad(base + align * offset[0]);
        t3.simdInternal_ = vec_add(t3.simdInternal_, t1);
        store(base + align * offset[0], t3);
        base[align * offset[0] + 2] += vec_extract(v2.simdInternal_, 0);

        t4               = simdLoad(base + align * offset[1]);
        t4.simdInternal_ = vec_add(t4.simdInternal_, t2);
        store(base + align * offset[1], t4);
        base[align * offset[1] + 2] += vec_extract(v2.simdInternal_, 1);
    }
}

template<int align>
static inline void gmx_simdcall
transposeScatterDecrU(double* base, const std::int32_t offset[], SimdDouble v0, SimdDouble v1, SimdDouble v2)
{
    if (align % 4 == 0)
    {
        __vector double t1, t2, t3, t4;
        SimdDouble      t5, t6, t7, t8;

        t1 = vec_mergeh(v0.simdInternal_, v1.simdInternal_);
        t2 = vec_mergel(v0.simdInternal_, v1.simdInternal_);
        t3 = vec_mergeh(v2.simdInternal_, vec_splats(0.0));
        t4 = vec_mergel(v2.simdInternal_, vec_splats(0.0));

        t5               = simdLoad(base + align * offset[0]);
        t6               = simdLoad(base + align * offset[0] + 2);
        t5.simdInternal_ = vec_sub(t5.simdInternal_, t1);
        t6.simdInternal_ = vec_sub(t6.simdInternal_, t3);
        store(base + align * offset[0], t5);
        store(base + align * offset[0] + 2, t6);

        t5               = simdLoad(base + align * offset[1]);
        t6               = simdLoad(base + align * offset[1] + 2);
        t5.simdInternal_ = vec_sub(t5.simdInternal_, t2);
        t6.simdInternal_ = vec_sub(t6.simdInternal_, t4);
        store(base + align * offset[1], t5);
        store(base + align * offset[1] + 2, t6);
    }
    else
    {
        __vector double t1, t2;
        SimdDouble      t3, t4;

        t1 = vec_mergeh(v0.simdInternal_, v1.simdInternal_);
        t2 = vec_mergel(v0.simdInternal_, v1.simdInternal_);

        t3               = simdLoad(base + align * offset[0]);
        t3.simdInternal_ = vec_sub(t3.simdInternal_, t1);
        store(base + align * offset[0], t3);
        base[align * offset[0] + 2] -= vec_extract(v2.simdInternal_, 0);

        t4               = simdLoad(base + align * offset[1]);
        t4.simdInternal_ = vec_sub(t4.simdInternal_, t2);
        store(base + align * offset[1], t4);
        base[align * offset[1] + 2] -= vec_extract(v2.simdInternal_, 1);
    }
}

static inline void gmx_simdcall expandScalarsToTriplets(SimdDouble  scalar,
                                                        SimdDouble* triplets0,
                                                        SimdDouble* triplets1,
                                                        SimdDouble* triplets2)
{
    triplets0->simdInternal_ = vec_mergeh(scalar.simdInternal_, scalar.simdInternal_);
    triplets1->simdInternal_ = scalar.simdInternal_;
    triplets2->simdInternal_ = vec_mergel(scalar.simdInternal_, scalar.simdInternal_);
}

template<int align>
static inline void gmx_simdcall gatherLoadBySimdIntTranspose(const double* base,
                                                             SimdDInt32    offset,
                                                             SimdDouble*   v0,
                                                             SimdDouble*   v1,
                                                             SimdDouble*   v2,
                                                             SimdDouble*   v3)
{
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ioffset[GMX_SIMD_DINT32_WIDTH];

    store(ioffset, offset);
    gatherLoadTranspose<align>(base, ioffset, v0, v1, v2, v3);
}

template<int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const double* base, SimdDInt32 offset, SimdDouble* v0, SimdDouble* v1)
{
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ioffset[GMX_SIMD_DINT32_WIDTH];

    store(ioffset, offset);
    gatherLoadTranspose<align>(base, ioffset, v0, v1);
}


template<int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const double* base, SimdDInt32 offset, SimdDouble* v0, SimdDouble* v1)
{
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ioffset[GMX_SIMD_DINT32_WIDTH];

    store(ioffset, offset);

    SimdDouble t1     = simdLoadU(base + align * ioffset[0]);
    SimdDouble t2     = simdLoadU(base + align * ioffset[1]);
    v0->simdInternal_ = vec_mergeh(t1.simdInternal_, t2.simdInternal_);
    v1->simdInternal_ = vec_mergel(t1.simdInternal_, t2.simdInternal_);
}

static inline double gmx_simdcall
reduceIncr4ReturnSum(double* m, SimdDouble v0, SimdDouble v1, SimdDouble v2, SimdDouble v3)
{
    __vector double t1, t2, t3, t4;

    t1 = vec_mergeh(v0.simdInternal_, v1.simdInternal_);
    t2 = vec_mergel(v0.simdInternal_, v1.simdInternal_);
    t3 = vec_mergeh(v2.simdInternal_, v3.simdInternal_);
    t4 = vec_mergel(v2.simdInternal_, v3.simdInternal_);

    t1 = vec_add(t1, t2);
    t3 = vec_add(t3, t4);

    *reinterpret_cast<__vector double*>(m) += t1;
    *reinterpret_cast<__vector double*>(m + 2) += t3;

    t1 = vec_add(t1, t3);
    return reduce(t1);
}

} // namespace gmx

#endif // GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_DOUBLE_H
