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

#ifndef GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD4_FLOAT_H

#include "config.h"

#include <arm_neon.h>

#include "gromacs/simd/impl_arm_neon/impl_arm_neon_simd4_float.h"

namespace gmx
{

static inline Simd4Float gmx_simdcall
fma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               vnegq_f32(vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
    };
}

static inline Simd4Float gmx_simdcall
fnma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fnms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               vnegq_f32(vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
    };
}

static inline Simd4Float gmx_simdcall
round(Simd4Float x)
{
    return {
               vrndnq_f32(x.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
trunc(Simd4Float x)
{
    return {
               vrndq_f32(x.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(Simd4FBool a)
{
    return (vmaxvq_u32(a.simdInternal_) != 0);
}

static inline float gmx_simdcall
reduce(Simd4Float a)
{
    float32x4_t b = a.simdInternal_;
    b = vpaddq_f32(b, b);
    b = vpaddq_f32(b, b);
    return vgetq_lane_f32(b, 0);
}

static inline float gmx_simdcall
dotProduct(Simd4Float a, Simd4Float b)
{
    Simd4Float c;

    c = a * b;
    /* set 4th element to 0, then add all of them */
    c.simdInternal_ = vsetq_lane_f32(0.0f, c.simdInternal_, 3);
    return reduce(c);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD4_FLOAT_H
