/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

#ifndef GMX_GPU_UTILS_DEVICE_UTILS_HIP_SYCL_H
#define GMX_GPU_UTILS_DEVICE_UTILS_HIP_SYCL_H

/*! \file
 *  \brief Shared device methods for SYCL and HIP targets.
 *
 *  This file contains the following shared methods: *
 *  * Cross lane move operations for AMD targets.
 *
 *  TODO add more shared methods to this file
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 */

// We only want to use the methods in this header when we are actually compiling device code
#if (defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)) || defined(__HIPCC__)

#    include "config.h"

#    include "gromacs/math/functions.h"

// We need to properly define the attributes so that some compilers don't choke on them.
#    if GMX_GPU_HIP
#        define GMX_DEVICE_ATTRIBUTE __device__
#    else
#        define GMX_DEVICE_ATTRIBUTE
#    endif

/* !\brief Cross-lane move operation using AMD DPP (Data-Parallel Primitives).
 *
 * Uses the __builtin_amdgcn_update_dpp intrinsic which expressed the data
 * movement but the compiler will combine these with subsequent instructions
 * if possible.
 *
 * Note that this is a generic implementation for any type T (for current use
 * it could be more simple).
 *
 * Ref: https://gpuopen.com/learn/amd-gcn-assembly-cross-lane-operations
 */
template<class T, int dppCtrl, int rowMask = 0xf, int bankMask = 0xf, bool boundCtrl = true>
GMX_DEVICE_ATTRIBUTE __attribute__((always_inline)) T amdDppUpdateShfl(const T& input)
{
    static constexpr int c_wordCount = gmx::divideRoundUp(sizeof(T), sizeof(int));

    struct V
    {
        int words[c_wordCount];
    };
    V wordList = __builtin_bit_cast(V, input);

#    pragma unroll
    for (int i = 0; i < c_wordCount; i++)
    {
        wordList.words[i] = __builtin_amdgcn_update_dpp(
                0, wordList.words[i], dppCtrl, rowMask, bankMask, boundCtrl);
    }

    return __builtin_bit_cast(T, wordList);
}

#    undef GMX_DEVICE_ATTRIBUTE

#endif /* Device code only */

#endif /* GMX_GPU_UTILS_WAVE_MOVE_DPP_H */
