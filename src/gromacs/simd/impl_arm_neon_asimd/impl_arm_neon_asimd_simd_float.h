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

#ifndef GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_FLOAT_H

#include "config.h"

#include <arm_neon.h>

#include "gromacs/simd/impl_arm_neon/impl_arm_neon_simd_float.h"

namespace gmx
{

static inline SimdFloat gmx_simdcall fma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return { vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdFloat gmx_simdcall fms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return { vnegq_f32(vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)) };
}

static inline SimdFloat gmx_simdcall fnma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return { vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdFloat gmx_simdcall fnms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return { vnegq_f32(vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)) };
}

static inline SimdFloat gmx_simdcall round(SimdFloat x)
{
    return { vrndnq_f32(x.simdInternal_) };
}

static inline SimdFloat gmx_simdcall trunc(SimdFloat x)
{
    return { vrndq_f32(x.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall cvtR2I(SimdFloat a)
{
    return { vcvtnq_s32_f32(a.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(SimdFBool a)
{
    return (vmaxvq_u32(a.simdInternal_) != 0);
}

static inline bool gmx_simdcall anyTrue(SimdFIBool a)
{
    return (vmaxvq_u32(a.simdInternal_) != 0);
}

static inline float gmx_simdcall reduce(SimdFloat a)
{
    float32x4_t b = a.simdInternal_;
    b             = vpaddq_f32(b, b);
    b             = vpaddq_f32(b, b);
    return vgetq_lane_f32(b, 0);
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_FLOAT_H
