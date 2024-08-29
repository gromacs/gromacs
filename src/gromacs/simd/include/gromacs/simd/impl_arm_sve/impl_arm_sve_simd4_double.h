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

#ifndef GMX_SIMD_IMPL_ARM_SVE_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_ARM_SVE_SIMD4_DOUBLE_H

#include "config.h"

#include <arm_sve.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "gromacs/math/utilities.h"

namespace gmx
{

class Simd4Double
{
private:
    typedef svfloat64_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    Simd4Double() {}

    Simd4Double(const double d) { this->simdInternal_ = svdup_n_f64(d); }

    Simd4Double(svfloat64_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

class Simd4DBool
{
private:
    typedef svbool_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    Simd4DBool() {}

    Simd4DBool(const bool b) { this->simdInternal_ = svdup_n_b64(b); }

    Simd4DBool(svbool_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

static inline Simd4Double gmx_simdcall load4(const double* m)
{
    assert(std::size_t(m) % (GMX_SIMD4_WIDTH * sizeof(double)) == 0);
    svbool_t pg = SVE_DOUBLE4_MASK;
    return { svld1_f64(pg, m) };
}

static inline void gmx_simdcall store4(double* m, Simd4Double a)
{
    assert(std::size_t(m) % (GMX_SIMD4_WIDTH * sizeof(double)) == 0);
    svbool_t pg = SVE_DOUBLE4_MASK;
    svst1_f64(pg, m, a.simdInternal_);
}

static inline Simd4Double gmx_simdcall load4U(const double* m)
{
    svbool_t pg = SVE_DOUBLE4_MASK;
    return { svld1_f64(pg, m) };
}

static inline void gmx_simdcall store4U(double* m, Simd4Double a)
{
    svbool_t pg = SVE_DOUBLE4_MASK;
    svst1_f64(pg, m, a.simdInternal_);
}

static inline Simd4Double gmx_simdcall simd4SetZeroD()
{
    return { svdup_n_f64(0.0) };
}

static inline Simd4Double gmx_simdcall operator&(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svreinterpret_f64_s64(svand_s64_z(
            pg, svreinterpret_s64_f64(a.simdInternal_), svreinterpret_s64_f64(b.simdInternal_))) };
}

static inline Simd4Double gmx_simdcall andNot(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svreinterpret_f64_s64(svbic_s64_z(
            pg, svreinterpret_s64_f64(b.simdInternal_), svreinterpret_s64_f64(a.simdInternal_))) };
}

static inline Simd4Double gmx_simdcall operator|(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svreinterpret_f64_s64(svorr_s64_z(
            pg, svreinterpret_s64_f64(a.simdInternal_), svreinterpret_s64_f64(b.simdInternal_))) };
}

static inline Simd4Double gmx_simdcall operator^(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svreinterpret_f64_s64(sveor_s64_z(
            pg, svreinterpret_s64_f64(a.simdInternal_), svreinterpret_s64_f64(b.simdInternal_))) };
}

static inline Simd4Double gmx_simdcall operator+(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svadd_f64_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Double gmx_simdcall operator-(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svsub_f64_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Double gmx_simdcall operator-(Simd4Double a)
{
    svbool_t pg = svptrue_b64();
    return { svneg_f64_z(pg, a.simdInternal_) };
}

static inline Simd4Double gmx_simdcall operator*(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svmul_f64_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Double gmx_simdcall fma(Simd4Double a, Simd4Double b, Simd4Double c)
{
    svbool_t pg = svptrue_b64();
    return { svmad_f64_z(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Double gmx_simdcall fms(Simd4Double a, Simd4Double b, Simd4Double c)
{
    svbool_t pg = svptrue_b64();
    return { svnmsb_f64_z(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Double gmx_simdcall fnma(Simd4Double a, Simd4Double b, Simd4Double c)
{
    svbool_t pg = svptrue_b64();
    return { svmsb_f64_z(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Double gmx_simdcall fnms(Simd4Double a, Simd4Double b, Simd4Double c)
{
    svbool_t pg = svptrue_b64();
    return { svnmad_f64_z(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Double gmx_simdcall rsqrt(Simd4Double x)
{
    svbool_t    pg = SVE_DOUBLE4_MASK;
    svfloat64_t f  = svsplice_f64(pg, x.simdInternal_, svdup_n_f64(1.0));
    return { svrsqrte_f64(f) };
}

static inline Simd4Double gmx_simdcall abs(Simd4Double x)
{
    svbool_t pg = svptrue_b64();
    return { svabs_f64_z(pg, x.simdInternal_) };
}

static inline Simd4Double gmx_simdcall max(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svmax_f64_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Double gmx_simdcall min(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svmin_f64_z(pg, a.simdInternal_, b.simdInternal_) };
}

// Round and trunc operations are defined at the end of this file, since they
// need to use float-to-integer and integer-to-float conversions.

static inline double gmx_simdcall reduce(Simd4Double a)
{
    svbool_t pg = SVE_DOUBLE4_MASK;
    return svadda_f64(pg, 0.0, a.simdInternal_);
}

static inline Simd4DBool gmx_simdcall operator==(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svcmpeq_f64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline Simd4DBool gmx_simdcall operator!=(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svcmpne_f64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline Simd4DBool gmx_simdcall operator<(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svcmplt_f64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline Simd4DBool gmx_simdcall operator<=(Simd4Double a, Simd4Double b)
{
    svbool_t pg = svptrue_b64();
    return { svcmple_f64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline Simd4DBool gmx_simdcall operator&&(Simd4DBool a, Simd4DBool b)
{
    svbool_t pg = svptrue_b64();
    return { svand_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline Simd4DBool gmx_simdcall operator||(Simd4DBool a, Simd4DBool b)
{
    svbool_t pg = svptrue_b64();
    return { svorr_b_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(Simd4DBool a)
{
    svbool_t pg = SVE_DOUBLE4_MASK;
    return svptest_any(pg, a.simdInternal_);
}

static inline Simd4Double gmx_simdcall selectByMask(Simd4Double a, Simd4DBool m)
{
    return { svsel_f64(m.simdInternal_, a.simdInternal_, svdup_n_f64(0.0)) };
}

static inline Simd4Double gmx_simdcall selectByNotMask(Simd4Double a, Simd4DBool m)
{
    return { svsel_f64(m.simdInternal_, svdup_n_f64(0.0), a.simdInternal_) };
}

static inline Simd4Double gmx_simdcall blend(Simd4Double a, Simd4Double b, Simd4DBool sel)
{
    return { svsel_f64(sel.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline Simd4Double gmx_simdcall round(Simd4Double x)
{
    svbool_t pg = svptrue_b64();
    return { svrinta_f64_z(pg, x.simdInternal_) };
}

static inline Simd4Double gmx_simdcall trunc(Simd4Double x)
{
    svbool_t pg = SVE_DOUBLE4_MASK;
    return { svcvt_f64_z(pg, svcvt_s64_z(pg, x.simdInternal_)) };
}

static inline double gmx_simdcall dotProduct(Simd4Double a, Simd4Double b)
{
    svbool_t pg = SVE_DOUBLE3_MASK;
    return svadda_f64(pg, 0.0, svmul_f64_z(pg, a.simdInternal_, b.simdInternal_));
}

static inline void gmx_simdcall transpose(Simd4Double* v0, Simd4Double* v1, Simd4Double* v2, Simd4Double* v3)
{
    svbool_t pg = SVE_DOUBLE4_MASK;
    double   tmp[16];
    svst1_f64(pg, tmp, v0->simdInternal_);
    svst1_f64(pg, tmp + 4, v1->simdInternal_);
    svst1_f64(pg, tmp + 8, v2->simdInternal_);
    svst1_f64(pg, tmp + 12, v3->simdInternal_);

    svfloat64x4_t vec = svld4_f64(pg, tmp);

    v0->simdInternal_ = svget4_f64(vec, 0);
    v1->simdInternal_ = svget4_f64(vec, 1);
    v2->simdInternal_ = svget4_f64(vec, 2);
    v3->simdInternal_ = svget4_f64(vec, 3);
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_SVE_SIMD4_DOUBLE_H
