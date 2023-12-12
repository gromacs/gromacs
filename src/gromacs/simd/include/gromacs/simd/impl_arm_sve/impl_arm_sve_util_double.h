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

#ifndef GMX_SIMD_IMPL_ARM_SVE_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_ARM_SVE_UTIL_DOUBLE_H

#include "config.h"

#include <arm_sve.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "gromacs/simd/impl_arm_sve/impl_arm_sve_simd_double.h"
#include "gromacs/utility/basedefinitions.h"


namespace gmx
{

namespace
{
inline void gmx_simdcall decrHsimd(double* m, SimdDouble a)
{
    // Make sure the memory pointer is aligned to half float SIMD width
    assert(std::size_t(m) % 32 == 0);

    svbool_t    pg = SVE_SIMD_DOUBLE_HALF_MASK;
    svfloat64_t v0, v1, v2, v3;
    v0 = svld1_f64(pg, m);
    v1 = svext_f64(a.simdInternal_, a.simdInternal_, GMX_SIMD_DOUBLE_WIDTH / 2);
    v2 = svadd_f64_x(pg, a.simdInternal_, v1);
    v3 = svsub_f64_x(pg, v0, v2);
    svst1_f64(pg, m, v3);
}
} // namespace

template<int align>
static inline void gmx_simdcall gatherLoadTranspose(const double*      base,
                                                    const std::int32_t offset[],
                                                    SimdDouble*        v0,
                                                    SimdDouble*        v1,
                                                    SimdDouble*        v2,
                                                    SimdDouble*        v3)
{
    assert(std::size_t(offset) % 16 == 0);
    assert(std::size_t(base) % 64 == 0);
    assert(align % 4 == 0);

    svint64_t offsets;
    svbool_t  pg = svptrue_b64();
    offsets      = svmul_n_s64_x(
            pg, svunpklo_s64(svld1_s32(SVE_SIMD_FLOAT_HALF_DOUBLE_MASK, offset)), align * sizeof(double));
    v0->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
    offsets           = svadd_n_s64_x(pg, offsets, sizeof(double));
    v1->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
    offsets           = svadd_n_s64_x(pg, offsets, sizeof(double));
    v2->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
    offsets           = svadd_n_s64_x(pg, offsets, sizeof(double));
    v3->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
}

template<int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const double* base, SimdDInt32 offset, SimdDouble* v0, SimdDouble* v1)
{
    // Base pointer must be aligned to the smaller of 2 elements and float SIMD width
    assert(std::size_t(base) % 8 == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % 2 == 0);

    svbool_t  pg = svptrue_b64();
    svint64_t offsets;
    offsets           = svmul_n_s64_x(pg, offset.simdInternal_, align * sizeof(double));
    v0->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
    offsets           = svadd_n_s64_x(pg, offsets, sizeof(double));
    v1->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
}

template<int align>
static inline void gmx_simdcall
gatherLoadTranspose(const double* base, const std::int32_t offset[], SimdDouble* v0, SimdDouble* v1)
{
    assert(std::size_t(offset) % 64 == 0);
    assert(std::size_t(base) % 8 == 0);
    assert(align % 2 == 0);

    SimdDInt32 offsets;
    svbool_t   pg         = SVE_SIMD_FLOAT_HALF_DOUBLE_MASK;
    offsets.simdInternal_ = svunpklo_s64(svld1_s32(pg, offset));
    gatherLoadBySimdIntTranspose<align>(base, offsets, v0, v1);
}

static const int c_simdBestPairAlignmentDouble = 2;

template<int align>
static inline void gmx_simdcall gatherLoadUTranspose(const double*      base,
                                                     const std::int32_t offset[],
                                                     SimdDouble*        v0,
                                                     SimdDouble*        v1,
                                                     SimdDouble*        v2)
{
    assert(std::size_t(offset) % 16 == 0);

    svint64_t offsets;
    svbool_t  pg = svptrue_b64();
    offsets      = svmul_n_s64_x(
            pg, svunpklo_s64(svld1_s32(SVE_SIMD_FLOAT_HALF_DOUBLE_MASK, offset)), align * sizeof(double));
    v0->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
    offsets           = svadd_n_s64_x(pg, offsets, sizeof(double));
    v1->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
    offsets           = svadd_n_s64_x(pg, offsets, sizeof(double));
    v2->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
}


template<int align>
static inline void gmx_simdcall transposeScatterStoreU(double*            base,
                                                       const std::int32_t offset[],
                                                       SimdDouble         v0,
                                                       SimdDouble         v1,
                                                       SimdDouble         v2)
{
    assert(std::size_t(offset) % 16 == 0);

    svint64_t offsets;
    svbool_t  pg = svptrue_b64();
    offsets      = svmul_n_s64_x(
            pg, svunpklo_s64(svld1_s32(SVE_SIMD_FLOAT_HALF_DOUBLE_MASK, offset)), align * sizeof(double));
    svst1_scatter_s64offset_f64(pg, base, offsets, v0.simdInternal_);
    offsets = svadd_n_s64_x(pg, offsets, sizeof(double));
    svst1_scatter_s64offset_f64(pg, base, offsets, v1.simdInternal_);
    offsets = svadd_n_s64_x(pg, offsets, sizeof(double));
    svst1_scatter_s64offset_f64(pg, base, offsets, v2.simdInternal_);
}


template<int align>
static inline void gmx_simdcall
transposeScatterIncrU(double* base, const std::int32_t offset[], SimdDouble v0, SimdDouble v1, SimdDouble v2)
{
    assert(std::size_t(offset) % 32 == 0);

    svbool_t                           pg = svptrue_b64();
    svfloat64x3_t                      v;
    alignas(GMX_SIMD_ALIGNMENT) double tvec[3 * GMX_SIMD_DOUBLE_WIDTH];
    v = svcreate3_f64(v0.simdInternal_, v1.simdInternal_, v2.simdInternal_);
    svst3_f64(pg, tvec, v);
#if GMX_SIMD_DOUBLE_WIDTH >= 3
    pg = SVE_SIMD4_DOUBLE_MASK;
    for (int i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        svfloat64_t t1 = svld1_f64(pg, base + align * offset[i]);
        svfloat64_t t2 = svld1_f64(pg, tvec + 3 * i);
        svfloat64_t t3 = svadd_f64_x(pg, t1, t2);
        svst1_f64(pg, base + align * offset[i], t3);
    }
#else
    for (std::size_t i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            base[align * offset[i] + j] += tvec[i * 3 + j];
        }
    }
#endif
}

template<int align>
static inline void gmx_simdcall
transposeScatterDecrU(double* base, const std::int32_t offset[], SimdDouble v0, SimdDouble v1, SimdDouble v2)
{
    assert(std::size_t(offset) % 16 == 0);

    svbool_t                           pg = svptrue_b64();
    svfloat64x3_t                      v;
    alignas(GMX_SIMD_ALIGNMENT) double tvec[3 * GMX_SIMD_DOUBLE_WIDTH];
    v = svcreate3_f64(v0.simdInternal_, v1.simdInternal_, v2.simdInternal_);
    svst3_f64(pg, tvec, v);
#if GMX_SIMD_DOUBLE_WIDTH >= 3
    pg = SVE_SIMD4_DOUBLE_MASK;
    for (int i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        svfloat64_t t1 = svld1_f64(pg, base + align * offset[i]);
        svfloat64_t t2 = svld1_f64(pg, tvec + 3 * i);
        svfloat64_t t3 = svsub_f64_x(pg, t1, t2);
        svst1_f64(pg, base + align * offset[i], t3);
    }
#else
    for (std::size_t i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            base[align * offset[i] + j] -= tvec[i * 3 + j];
        }
    }
#endif
}

static inline void gmx_simdcall expandScalarsToTriplets(SimdDouble  scalar,
                                                        SimdDouble* triplets0,
                                                        SimdDouble* triplets1,
                                                        SimdDouble* triplets2)
{
    assert(GMX_SIMD_DOUBLE_WIDTH <= 16);
    uint64_t   ind[48] = { 0,  0,  0,  1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,
                         5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9,  9,  10, 10,
                         10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15 };
    svbool_t   pg;
    svuint64_t idx;

    pg                       = svptrue_b64();
    idx                      = svld1_u64(pg, ind);
    triplets0->simdInternal_ = svtbl_f64(scalar.simdInternal_, idx);
    idx                      = svld1_u64(pg, ind + GMX_SIMD_DOUBLE_WIDTH);
    triplets1->simdInternal_ = svtbl_f64(scalar.simdInternal_, idx);
    idx                      = svld1_u64(pg, ind + 2 * GMX_SIMD_DOUBLE_WIDTH);
    triplets2->simdInternal_ = svtbl_f64(scalar.simdInternal_, idx);
}

template<int align>
static inline void gmx_simdcall gatherLoadBySimdIntTranspose(const double* base,
                                                             SimdDInt32    offset,
                                                             SimdDouble*   v0,
                                                             SimdDouble*   v1,
                                                             SimdDouble*   v2,
                                                             SimdDouble*   v3)
{
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ioffset[GMX_SIMD_FINT32_WIDTH];

    assert(std::size_t(base) % 16 == 0);
    assert(align % 4 == 0);

    store(ioffset, offset);
    gatherLoadTranspose<align>(base, ioffset, v0, v1, v2, v3);
}


template<int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const double* base, SimdDInt32 offset, SimdDouble* v0, SimdDouble* v1)
{
    svbool_t  pg      = svptrue_b64();
    svint64_t offsets = svmul_n_s64_x(pg, offset.simdInternal_, align * sizeof(double));
    v0->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
    offsets           = svadd_n_s64_x(pg, offsets, sizeof(double));
    v1->simdInternal_ = svld1_gather_s64offset_f64(pg, base, offsets);
}

static inline double gmx_simdcall
reduceIncr4ReturnSum(double* m, SimdDouble v0, SimdDouble v1, SimdDouble v2, SimdDouble v3)
{
    assert(std::size_t(m) % 16 == 0);
    svbool_t pg = svptrue_b64();
    double   sum[4];
    sum[0] = svadda_f64(pg, 0.0, v0.simdInternal_);
    sum[1] = svadda_f64(pg, 0.0, v1.simdInternal_);
    sum[2] = svadda_f64(pg, 0.0, v2.simdInternal_);
    sum[3] = svadda_f64(pg, 0.0, v3.simdInternal_);
#if GMX_SIMD_DOUBLE_WIDTH >= 4
    pg             = SVE_SIMD4_DOUBLE_MASK;
    svfloat64_t _m = svld1_f64(pg, m);
    svfloat64_t _s = svld1_f64(pg, sum);
    svst1_f64(pg, m, svadd_f64_x(pg, _m, _s));
    return svadda_f64(pg, 0.0, _s);
#else
    double res = 0;
    for (int i = 0; i < 4; i++)
    {
        m[i] += sum[i];
        res += sum[i];
    }
    return res;
#endif
}

static inline SimdDouble gmx_simdcall loadDualHsimd(const double* m0, const double* m1)
{
    svfloat64_t v0, v1;
    svbool_t    pg = SVE_SIMD_DOUBLE_HALF_MASK;
    v0             = svld1_f64(pg, m0);
    v1             = svld1_f64(pg, m1);
    return { svsplice_f64(pg, v0, v1) };
}

static inline SimdDouble gmx_simdcall loadDuplicateHsimd(const double* m)
{
    svfloat64_t v;
    svbool_t    pg = SVE_SIMD_DOUBLE_HALF_MASK;
    v              = svld1_f64(pg, m);
    return { svsplice_f64(pg, v, v) };
}

static inline SimdDouble gmx_simdcall loadU1DualHsimd(const double* m)
{
    svfloat64_t v0, v1;
    svbool_t    pg = SVE_SIMD_DOUBLE_HALF_MASK;
    v0             = svdup_n_f64(m[0]);
    v1             = svdup_n_f64(m[1]);
    return { svsplice_f64(pg, v0, v1) };
}

static inline void gmx_simdcall storeDualHsimd(double* m0, double* m1, SimdDouble a)
{
    svbool_t pg = SVE_SIMD_DOUBLE_HALF_MASK;
    svst1_f64(pg, m0, a.simdInternal_);
    pg = sveor_b_z(svptrue_b64(), pg, svptrue_b64());
    svst1_f64(pg, m1 - GMX_SIMD_DOUBLE_WIDTH / 2, a.simdInternal_);
}

static inline void gmx_simdcall incrDualHsimd(double* m0, double* m1, SimdDouble a)
{
    // Make sure the memory pointer is aligned to half float SIMD width
    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    svbool_t    pg = SVE_SIMD_DOUBLE_HALF_MASK;
    svfloat64_t v0, v2, v3;
    v0 = svld1_f64(pg, m0);
    v2 = svadd_f64_x(pg, v0, a.simdInternal_);
    svst1_f64(pg, m0, v2);
    v0 = svld1_f64(pg, m1);
    v3 = svext_f64(a.simdInternal_, a.simdInternal_, GMX_SIMD_DOUBLE_WIDTH / 2);
    v2 = svadd_f64_x(pg, v0, v3);
    svst1_f64(pg, m1, v2);
}

static inline void gmx_simdcall decr3Hsimd(double* m, SimdDouble a0, SimdDouble a1, SimdDouble a2)
{
    decrHsimd(m, a0);
    decrHsimd(m + GMX_SIMD_DOUBLE_WIDTH / 2, a1);
    decrHsimd(m + GMX_SIMD_DOUBLE_WIDTH, a2);
}

static inline double gmx_simdcall reduceIncr4ReturnSumHsimd(double* m, SimdDouble v0, SimdDouble v1)
{
    svbool_t pg = SVE_SIMD_DOUBLE_HALF_MASK;
    double   sum[4];
    sum[0] = svadda_f64(pg, 0.0, v0.simdInternal_);
    sum[2] = svadda_f64(pg, 0.0, v1.simdInternal_);
    pg     = sveor_b_z(svptrue_b64(), pg, svptrue_b64());
    sum[1] = svadda_f64(pg, 0.0, v0.simdInternal_);
    sum[3] = svadda_f64(pg, 0.0, v1.simdInternal_);

#if GMX_SIMD_DOUBLE_WIDTH >= 4
    pg             = SVE_SIMD4_DOUBLE_MASK;
    svfloat64_t _m = svld1_f64(pg, m);
    svfloat64_t _s = svld1_f64(pg, sum);
    svst1_f64(pg, m, svadd_f64_x(pg, _m, _s));
    return svadda_f64(pg, 0.0, _s);
#else
    double res = 0;
    for (int i = 0; i < 4; i++)
    {
        m[i] += sum[i];
        res += sum[i];
    }
    return res;
#endif
}

template<int align>
static inline void gmx_simdcall gatherLoadTransposeHsimd(const double*      base0,
                                                         const double*      base1,
                                                         const std::int32_t offset[],
                                                         SimdDouble*        v0,
                                                         SimdDouble*        v1)
{
    svint64_t   offsets;
    svbool_t    pg = SVE_SIMD_DOUBLE_HALF_MASK;
    svfloat64_t _v0, _v1;
    offsets = svmul_n_s64_x(
            pg, svunpklo(svld1_s32(SVE_SIMD_FLOAT_HALF_DOUBLE_MASK, offset)), align * sizeof(double));
    _v0               = svld1_gather_s64offset_f64(pg, base0, offsets);
    _v1               = svld1_gather_s64offset_f64(pg, base1, offsets);
    v0->simdInternal_ = svsplice_f64(pg, _v0, _v1);
    offsets           = svadd_n_s64_x(pg, offsets, sizeof(double));
    _v0               = svld1_gather_s64offset_f64(pg, base0, offsets);
    _v1               = svld1_gather_s64offset_f64(pg, base1, offsets);
    v1->simdInternal_ = svsplice_f64(pg, _v0, _v1);
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_SVE_UTIL_DOUBLE_H
