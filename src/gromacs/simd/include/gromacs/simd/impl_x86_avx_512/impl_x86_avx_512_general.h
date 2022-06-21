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

#ifndef GMX_SIMD_IMPL_X86_AVX_512_GENERAL_H
#define GMX_SIMD_IMPL_X86_AVX_512_GENERAL_H

#include <immintrin.h>

namespace gmx
{

static inline void simdPrefetch(const void* m)
{
    _mm_prefetch(reinterpret_cast<const char*>(m), _MM_HINT_T0);
}

/*! \brief Return integer from AVX-512 mask
 *
 *  \param m  Mask suitable for use with AVX-512 instructions
 *
 *  \return Short integer representation of mask
 */
static inline short avx512Mask2Int(__mmask16 m)
{
    return static_cast<short>(m);
}

/*! \brief Return AVX-512 mask from integer
 *
 *  \param m  Short integer
 *
 *  \return Mask suitable for use with AVX-512 instructions.
 */
static inline __mmask16 avx512Int2Mask(short i)
{
    return static_cast<__mmask16>(i);
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_512_GENERAL_H
