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

#ifndef GMX_SIMD_IMPL_ARM_SVE_GENERAL_H
#define GMX_SIMD_IMPL_ARM_SVE_GENERAL_H

namespace gmx
{

static inline void simdPrefetch(void* m)
{
#ifdef __GNUC__
    __builtin_prefetch(m);
#endif
}

#define SVE_DOUBLE_MASK svptrue_b64()
#define SVE_DINT32_MASK svptrue_b64()
#define SVE_SIMD_FLOAT_HALF_DOUBLE_MASK svwhilelt_b32(0, (int32_t)GMX_SIMD_DINT32_WIDTH)
#define SVE_SIMD_DOUBLE_HALF_MASK svwhilelt_b64(0, (int32_t)GMX_SIMD_DOUBLE_WIDTH / 2)
#define SVE_FLOAT_HALF_MASK svwhilelt_b32(0, GMX_SIMD_FLOAT_WIDTH / 2)
#define SVE_FINT32_HALF_MASK svwhilelt_b32(0, GMX_SIMD_FLOAT_WIDTH / 2)
#define SVE_FLOAT4_MASK svptrue_pat_b32(SV_VL4)
#define SVE_FLOAT3_MASK svptrue_pat_b32(SV_VL3)
#define SVE_DOUBLE4_MASK svptrue_pat_b64(SV_VL4)
#define SVE_DOUBLE3_MASK svptrue_pat_b64(SV_VL3)
} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_SVE_GENERAL_H
