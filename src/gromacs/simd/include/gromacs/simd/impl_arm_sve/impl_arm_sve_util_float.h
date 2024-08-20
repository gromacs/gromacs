/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

/*
 * armv8+sve support to GROMACS was contributed by the Research Organization for
 * Information Science and Technology (RIST).
 * Copyright (c) 2020 Research Organization for Information Science and Technology (RIST).
 */

#ifndef GMX_SIMD_IMPL_ARM_SVE_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_ARM_SVE_UTIL_FLOAT_H

#include "config.h"

#include <arm_sve.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "gromacs/simd/impl_arm_sve/impl_arm_sve_simd_float.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

template<int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float* base, SimdFInt32 offset, SimdFloat* v0, SimdFloat* v1)
{
    // Base pointer must be aligned to the smaller of 2 elements and float SIMD width
    assert(std::size_t(base) % (std::min(GMX_SIMD_FLOAT_WIDTH, 2) * sizeof(float)) == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 2) == 0);

    if (align < 2)
    {
        svbool_t  pg = svptrue_b32();
        svint32_t offsets;
        offsets           = svmul_n_s32_x(pg, offset.simdInternal_, align * 4);
        v0->simdInternal_ = svld1_gather_s32offset_f32(pg, base, offsets);
        offsets           = svadd_n_s32_x(pg, offsets, 4);
        v1->simdInternal_ = svld1_gather_s32offset_f32(pg, base, offsets);
    }
    else if (2 == align)
    {
        svbool_t    pg    = svptrue_b32();
        svfloat32_t t0    = svreinterpret_f32_u64(svld1_gather_s64index_u64(
                svunpklo_b(pg), (uint64_t*)base, svunpklo_s64(offset.simdInternal_)));
        svfloat32_t t1    = svreinterpret_f32_u64(svld1_gather_s64index_u64(
                svunpkhi_b(pg), (uint64_t*)base, svunpkhi_s64(offset.simdInternal_)));
        v0->simdInternal_ = svuzp1(t0, t1);
        v1->simdInternal_ = svuzp2(t0, t1);
    }
    else
    {
        svbool_t    pg      = svptrue_b32();
        svint32_t   offsets = svmul_n_s32_x(pg, offset.simdInternal_, align / 2);
        svfloat32_t t0      = svreinterpret_f32_u64(
                svld1_gather_s64index_u64(svunpklo_b(pg), (uint64_t*)base, svunpklo_s64(offsets)));
        svfloat32_t t1 = svreinterpret_f32_u64(
                svld1_gather_s64index_u64(svunpkhi_b(pg), (uint64_t*)base, svunpkhi_s64(offsets)));
        v0->simdInternal_ = svuzp1(t0, t1);
        v1->simdInternal_ = svuzp2(t0, t1);
    }
}

template<int align>
static inline void gmx_simdcall gatherLoadTranspose(const float*       base,
                                                    const std::int32_t offset[],
                                                    SimdFloat*         v0,
                                                    SimdFloat*         v1,
                                                    SimdFloat*         v2,
                                                    SimdFloat*         v3)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH * sizeof(std::int32_t)) == 0);
    // Base pointer must be aligned to the smaller of 4 elements and float SIMD width
    assert(std::size_t(base) % (std::min(GMX_SIMD_FLOAT_WIDTH, 4) * sizeof(float)) == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 4) == 0);

    svint32_t offsets;
    offsets = svld1_s32(svptrue_b32(), offset);
    gatherLoadBySimdIntTranspose<align>(base, offsets, v0, v1, v2, v3);
}

template<int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float* base, const std::int32_t offset[], SimdFloat* v0, SimdFloat* v1)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH * sizeof(std::int32_t)) == 0);
    // Base pointer must be aligned to the smaller of 2 elements and float SIMD width
    assert(std::size_t(base) % (std::min(GMX_SIMD_FLOAT_WIDTH, 2) * sizeof(float)) == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 2) == 0);

    SimdFInt32 offsets;
    svbool_t   pg         = svptrue_b32();
    offsets.simdInternal_ = svld1(pg, offset);
    gatherLoadBySimdIntTranspose<align>(base, offsets, v0, v1);
}

static const int c_simdBestPairAlignmentFloat = 2;

template<int align>
static inline void gmx_simdcall gatherLoadUTranspose(const float*       base,
                                                     const std::int32_t offset[],
                                                     SimdFloat*         v0,
                                                     SimdFloat*         v1,
                                                     SimdFloat*         v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH * sizeof(std::int32_t)) == 0);

    svint32_t offsets;
    svbool_t  pg      = svptrue_b32();
    offsets           = svmul_n_s32_x(pg, svld1_s32(pg, offset), align * 4);
    v0->simdInternal_ = svld1_gather_s32offset_f32(pg, base, offsets);
    offsets           = svadd_n_s32_x(pg, offsets, 4);
    v1->simdInternal_ = svld1_gather_s32offset_f32(pg, base, offsets);
    offsets           = svadd_n_s32_x(pg, offsets, 4);
    v2->simdInternal_ = svld1_gather_s32offset_f32(pg, base, offsets);
}


template<int align>
static inline void gmx_simdcall
transposeScatterStoreU(float* base, const std::int32_t offset[], SimdFloat v0, SimdFloat v1, SimdFloat v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH * sizeof(std::int32_t)) == 0);

    svint32_t offsets;
    svbool_t  pg = svptrue_b32();
    offsets      = svmul_n_s32_x(pg, svld1_s32(pg, offset), align * 4);
    svst1_scatter_s32offset_f32(pg, base, offsets, v0.simdInternal_);
    offsets = svadd_n_s32_x(pg, offsets, 4);
    svst1_scatter_s32offset_f32(pg, base, offsets, v1.simdInternal_);
    offsets = svadd_n_s32_x(pg, offsets, 4);
    svst1_scatter_s32offset_f32(pg, base, offsets, v2.simdInternal_);
}


template<int align>
static inline void gmx_simdcall
transposeScatterIncrU(float* base, const std::int32_t offset[], SimdFloat v0, SimdFloat v1, SimdFloat v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH * sizeof(std::int32_t)) == 0);

    svbool_t                          pg = svptrue_b32();
    svfloat32x3_t                     v;
    alignas(GMX_SIMD_ALIGNMENT) float tvec[3 * GMX_SIMD_FLOAT_WIDTH];
    v = svcreate3_f32(v0.simdInternal_, v1.simdInternal_, v2.simdInternal_);
    svst3_f32(pg, tvec, v);
    pg = SVE_FLOAT3_MASK;
    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        svfloat32_t t1 = svld1_f32(pg, base + align * offset[i]);
        svfloat32_t t2 = svld1_f32(pg, tvec + 3 * i);
        svfloat32_t t3 = svadd_f32_x(pg, t1, t2);
        svst1_f32(pg, base + align * offset[i], t3);
    }
}

template<int align>
static inline void gmx_simdcall
transposeScatterDecrU(float* base, const std::int32_t offset[], SimdFloat v0, SimdFloat v1, SimdFloat v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH * sizeof(std::int32_t)) == 0);

    svbool_t                          pg = svptrue_b32();
    svfloat32x3_t                     v;
    alignas(GMX_SIMD_ALIGNMENT) float tvec[3 * GMX_SIMD_FLOAT_WIDTH];
    v = svcreate3_f32(v0.simdInternal_, v1.simdInternal_, v2.simdInternal_);
    svst3_f32(pg, tvec, v);
    pg = SVE_FLOAT3_MASK;
    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        svfloat32_t t1 = svld1_f32(pg, base + align * offset[i]);
        svfloat32_t t2 = svld1_f32(pg, tvec + 3 * i);
        svfloat32_t t3 = svsub_f32_x(pg, t1, t2);
        svst1_f32(pg, base + align * offset[i], t3);
    }
}

static inline void gmx_simdcall expandScalarsToTriplets(SimdFloat  scalar,
                                                        SimdFloat* triplets0,
                                                        SimdFloat* triplets1,
                                                        SimdFloat* triplets2)
{
    assert(GMX_SIMD_FLOAT_WIDTH <= 16);
    uint32_t   ind[48] = { 0,  0,  0,  1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,
                         5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9,  9,  10, 10,
                         10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15 };
    svbool_t   pg;
    svuint32_t idx;

    pg                       = svptrue_b32();
    idx                      = svld1_u32(pg, ind);
    triplets0->simdInternal_ = svtbl_f32(scalar.simdInternal_, idx);
    idx                      = svld1_u32(pg, ind + GMX_SIMD_FLOAT_WIDTH);
    triplets1->simdInternal_ = svtbl_f32(scalar.simdInternal_, idx);
    idx                      = svld1_u32(pg, ind + 2 * GMX_SIMD_FLOAT_WIDTH);
    triplets2->simdInternal_ = svtbl_f32(scalar.simdInternal_, idx);
}

template<int align>
static inline void gmx_simdcall gatherLoadBySimdIntTranspose(const float* base,
                                                             SimdFInt32   offset,
                                                             SimdFloat*   v0,
                                                             SimdFloat*   v1,
                                                             SimdFloat*   v2,
                                                             SimdFloat*   v3)
{
    // Base pointer must be aligned to the smaller of 4 elements and float SIMD width
    assert(std::size_t(base) % (std::min(GMX_SIMD_FLOAT_WIDTH, 4) * sizeof(float)) == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 4) == 0);

    svbool_t pg          = svptrue_b32();
    offset.simdInternal_ = svmul_n_s32_x(pg, offset.simdInternal_, align);
    v0->simdInternal_    = svld1_gather_s32index_f32(pg, base, offset.simdInternal_);
    offset.simdInternal_ = svadd_n_s32_x(pg, offset.simdInternal_, 1);
    v1->simdInternal_    = svld1_gather_s32index_f32(pg, base, offset.simdInternal_);
    offset.simdInternal_ = svadd_n_s32_x(pg, offset.simdInternal_, 1);
    v2->simdInternal_    = svld1_gather_s32index_f32(pg, base, offset.simdInternal_);
    offset.simdInternal_ = svadd_n_s32_x(pg, offset.simdInternal_, 1);
    v3->simdInternal_    = svld1_gather_s32index_f32(pg, base, offset.simdInternal_);
}


template<int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const float* base, SimdFInt32 offset, SimdFloat* v0, SimdFloat* v1)
{
    svbool_t  pg      = svptrue_b32();
    svint32_t offsets = svmul_n_s32_x(pg, offset.simdInternal_, align * 4);
    v0->simdInternal_ = svld1_gather_s32offset_f32(pg, base, offsets);
    offsets           = svadd_n_s32_x(pg, offsets, 4);
    v1->simdInternal_ = svld1_gather_s32offset_f32(pg, base, offsets);
}

static inline float gmx_simdcall reduceIncr4ReturnSum(float* m, SimdFloat v0, SimdFloat v1, SimdFloat v2, SimdFloat v3)
{
    // Make sure the memory pointer is aligned to the smaller of 4 elements and float SIMD width
    assert(std::size_t(m) % (std::min(GMX_SIMD_FLOAT_WIDTH, 4) * sizeof(float)) == 0);

    svbool_t    pg = svptrue_b32();
    svfloat32_t _m, _s;
    float32_t   sum[4];
    sum[0] = svadda_f32(pg, 0.0f, v0.simdInternal_);
    sum[1] = svadda_f32(pg, 0.0f, v1.simdInternal_);
    sum[2] = svadda_f32(pg, 0.0f, v2.simdInternal_);
    sum[3] = svadda_f32(pg, 0.0f, v3.simdInternal_);
    pg     = SVE_FLOAT4_MASK;
    _m     = svld1_f32(pg, m);
    _s     = svld1_f32(pg, sum);
    svst1_f32(pg, m, svadd_f32_x(pg, _m, _s));
    return svadda_f32(pg, 0.0f, _s);
}

static inline SimdFloat gmx_simdcall loadDualHsimd(const float* m0, const float* m1)
{
    svfloat32_t v0, v1;
    svbool_t    pg = SVE_FLOAT_HALF_MASK;
    v0             = svld1_f32(pg, m0);
    v1             = svld1_f32(pg, m1);
    return { svsplice_f32(pg, v0, v1) };
}

static inline SimdFloat gmx_simdcall loadDuplicateHsimd(const float* m)
{
    svfloat32_t v;
    svbool_t    pg = SVE_FLOAT_HALF_MASK;
    v              = svld1_f32(pg, m);
    return { svsplice_f32(pg, v, v) };
}

static inline SimdFloat gmx_simdcall loadU1DualHsimd(const float* m)
{
    svfloat32_t v0, v1;
    svbool_t    pg = SVE_FLOAT_HALF_MASK;
    v0             = svdup_n_f32(m[0]);
    v1             = svdup_n_f32(m[1]);
    return { svsplice_f32(pg, v0, v1) };
}

static inline void gmx_simdcall storeDualHsimd(float* m0, float* m1, SimdFloat a)
{
    svbool_t pg = SVE_FLOAT_HALF_MASK;
    svst1_f32(pg, m0, a.simdInternal_);
    svst1_f32(pg, m1, svext_f32(a.simdInternal_, a.simdInternal_, GMX_SIMD_FLOAT_WIDTH / 2));
}

static inline void gmx_simdcall incrDualHsimd(float* m0, float* m1, SimdFloat a)
{
    // Make sure the memory pointer is aligned to half float SIMD width
    assert(std::size_t(m0) % (GMX_SIMD_FLOAT_WIDTH * sizeof(float) / 2) == 0);
    assert(std::size_t(m1) % (GMX_SIMD_FLOAT_WIDTH * sizeof(float) / 2) == 0);

    svbool_t    pg = SVE_FLOAT_HALF_MASK;
    svfloat32_t v0, v2, v3;
    v0 = svld1_f32(pg, m0);
    v2 = svadd_f32_x(pg, v0, a.simdInternal_);
    svst1_f32(pg, m0, v2);
    v0 = svld1_f32(pg, m1);
    v3 = svext_f32(a.simdInternal_, a.simdInternal_, GMX_SIMD_FLOAT_WIDTH / 2);
    v2 = svadd_f32_x(pg, v0, v3);
    svst1_f32(pg, m1, v2);
}

static inline void gmx_simdcall decr3Hsimd(float* m, SimdFloat a0, SimdFloat a1, SimdFloat a2)
{
    svbool_t    pg  = svptrue_b32();
    svbool_t    pg2 = SVE_FLOAT_HALF_MASK;
    svfloat32_t v0, v1, v2, v3;
    v0 = svld1_f32(pg, m);
    v1 = svext_f32(a0.simdInternal_, a1.simdInternal_, GMX_SIMD_FLOAT_WIDTH / 2);
    v2 = svsel_f32(pg2, a0.simdInternal_, a1.simdInternal_);
    v1 = svadd_f32_x(pg, v1, v2);
    v0 = svsub_f32_z(pg, v0, v1);
    svst1_f32(pg, m, v0);
    v0 = svld1_f32(pg2, m + GMX_SIMD_FLOAT_WIDTH);
    v1 = svext_f32(a2.simdInternal_, a0.simdInternal_, GMX_SIMD_FLOAT_WIDTH / 2);
    v2 = svadd_f32_x(pg2, a2.simdInternal_, v1);
    v3 = svsub_f32_x(pg2, v0, v2);
    svst1_f32(pg2, m + GMX_SIMD_FLOAT_WIDTH, v3);
}

static inline float gmx_simdcall reduceIncr4ReturnSumHsimd(float* m, SimdFloat v0, SimdFloat v1)
{
    svbool_t    pg  = SVE_FLOAT_HALF_MASK;
    svbool_t    pg2 = sveor_b_z(svptrue_b32(), pg, svptrue_b32());
    svfloat32_t _m, _s;

    _s = svdup_n_f32(0.0f);
    _s = svinsr_n_f32(_s, svaddv_f32(pg2, v1.simdInternal_));
    _s = svinsr_n_f32(_s, svaddv_f32(pg, v1.simdInternal_));
    _s = svinsr_n_f32(_s, svaddv_f32(pg2, v0.simdInternal_));
    _s = svinsr_n_f32(_s, svaddv_f32(pg, v0.simdInternal_));

    pg = SVE_FLOAT4_MASK;
    _m = svld1_f32(pg, m);
    svst1_f32(pg, m, svadd_f32_x(pg, _m, _s));
    return svaddv_f32(pg, _s);
}

template<int align>
static inline void gmx_simdcall gatherLoadTransposeHsimd(const float*       base0,
                                                         const float*       base1,
                                                         const std::int32_t offset[],
                                                         SimdFloat*         v0,
                                                         SimdFloat*         v1)
{
    svint64_t   offsets = svunpklo_s64(svld1_s32(svptrue_b32(), offset));
    svfloat32_t _v0, _v1;
    if (2 == align)
    {
        _v0 = svreinterpret_f32_f64(svld1_gather_s64index_f64(SVE_DOUBLE_MASK, (double*)base0, offsets));
        _v1 = svreinterpret_f32_f64(svld1_gather_s64index_f64(SVE_DOUBLE_MASK, (double*)base1, offsets));
    }
    else
    {
        offsets = svmul_n_s64_x(svptrue_b64(), offsets, align * 4);
        _v0 = svreinterpret_f32_f64(svld1_gather_s64offset_f64(SVE_DOUBLE_MASK, (double*)base0, offsets));
        _v1 = svreinterpret_f32_f64(svld1_gather_s64offset_f64(SVE_DOUBLE_MASK, (double*)base1, offsets));
    }
    v0->simdInternal_ = svuzp1(_v0, _v1);
    v1->simdInternal_ = svuzp2(_v0, _v1);
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_SVE_UTIL_FLOAT_H
