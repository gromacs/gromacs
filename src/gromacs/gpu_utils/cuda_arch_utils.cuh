/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef CUDA_ARCH_UTILS_CUH_
#define CUDA_ARCH_UTILS_CUH_

#include "config.h"

#include "gromacs/utility/basedefinitions.h"

/*! \file
 *  \brief CUDA arch dependent definitions.
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */

/* GMX_PTX_ARCH is set to the virtual arch (PTX) version targeted by
 * the current compiler pass or zero for the host pass and it is
 * intended to be used instead of __CUDA_ARCH__.
 */
#ifndef __CUDA_ARCH__
    #define GMX_PTX_ARCH 0
#else
    #define GMX_PTX_ARCH __CUDA_ARCH__
#endif

/* Until CC 5.2 and likely for the near future all NVIDIA architectures
   have a warp size of 32, but this could change later. If it does, the
   following constants should depend on the value of GMX_PTX_ARCH.
 */
static const int warp_size      = 32;
static const int warp_size_log2 = 5;
/*! \brief Bitmask corresponding to all threads active in a warp.
 *  NOTE that here too we assume 32-wide warps.
 */
static const unsigned int c_fullWarpMask = 0xffffffff;

/* Below are backward-compatibility wrappers for CUDA 9 warp-wide intrinsics. */

/*! \brief Compatibility wrapper around the CUDA __syncwarp() instrinsic.  */
static __forceinline__ __device__
void gmx_syncwarp(const unsigned int activeMask = c_fullWarpMask)
{
#if GMX_CUDA_VERSION < 9000
    /* no sync needed on pre-Volta. */
    GMX_UNUSED_VALUE(activeMask);
#else
    __syncwarp(activeMask);
#endif
}

/*! \brief Compatibility wrapper around the CUDA __ballot()/__ballot_sync() instrinsic.  */
static __forceinline__ __device__
unsigned int gmx_ballot_sync(const unsigned int activeMask,
                             const int          pred)
{
#if GMX_CUDA_VERSION < 9000
    GMX_UNUSED_VALUE(activeMask);
    return __ballot(pred);
#else
    return __ballot_sync(activeMask, pred);
#endif
}

/*! \brief Compatibility wrapper around the CUDA __any()/__any_sync() instrinsic.  */
static __forceinline__ __device__
int gmx_any_sync(const unsigned int activeMask,
                 const int          pred)
{
#if GMX_CUDA_VERSION < 9000
    GMX_UNUSED_VALUE(activeMask);
    return __any(pred);
#else
    return __any_sync(activeMask, pred);
#endif
}

/*! \brief Compatibility wrapper around the CUDA __shfl_up()/__shfl_up_sync() instrinsic.  */
template <typename T>
static __forceinline__ __device__
T gmx_shfl_up_sync(const unsigned int activeMask,
                   const T            var,
                   unsigned int       offset,
                   int                width = warp_size)
{
#if GMX_CUDA_VERSION < 9000
    GMX_UNUSED_VALUE(activeMask);
    return __shfl_up(var, offset, width);
#else
    return __shfl_up_sync(activeMask, var, offset, width);
#endif
}

/*! \brief Compatibility wrapper around the CUDA __shfl_down()/__shfl_down_sync() instrinsic.  */
template <typename T>
static __forceinline__ __device__
T gmx_shfl_down_sync(const unsigned int activeMask,
                     const T            var,
                     unsigned int       offset,
                     int                width = warp_size)
{
#if GMX_CUDA_VERSION < 9000
    GMX_UNUSED_VALUE(activeMask);
    return __shfl_down(var, offset, width);
#else
    return __shfl_down_sync(activeMask, var, offset, width);
#endif
}

/*! \brief Allow disabling CUDA textures using the GMX_DISABLE_CUDA_TEXTURES macro.
 *
 *  This option will not influence functionality. All features using textures ought
 *  to have fallback for texture-less reads (direct/LDG loads), all new code needs
 *  to provide fallback code.
 */
#if defined GMX_DISABLE_CUDA_TEXTURES
#define DISABLE_CUDA_TEXTURES 1
#else
#define DISABLE_CUDA_TEXTURES 0
#endif

/*! \brief True if the use of texture fetch in the CUDA kernels is disabled. */
static const bool c_disableCudaTextures = DISABLE_CUDA_TEXTURES;


/* CUDA architecture technical characteristics. Needs macros because it is used
 * in the __launch_bounds__ function qualifiers and might need it in preprocessor
 * conditionals.
 *
 */
#if GMX_PTX_ARCH > 0
    #if   GMX_PTX_ARCH <= 210  // CC 2.x
        #define GMX_CUDA_MAX_BLOCKS_PER_MP   8
        #define GMX_CUDA_MAX_THREADS_PER_MP  1536
    #elif GMX_PTX_ARCH <= 370  // CC 3.x
        #define GMX_CUDA_MAX_BLOCKS_PER_MP   16
        #define GMX_CUDA_MAX_THREADS_PER_MP  2048
    #else // CC 5.x, 6.x
          /* Note that this final branch covers all future architectures (current gen
           * is 6.x as of writing), hence assuming that these *currently defined* upper
           * limits will not be lowered.
           */
        #define GMX_CUDA_MAX_BLOCKS_PER_MP   32
        #define GMX_CUDA_MAX_THREADS_PER_MP  2048
    #endif
#else
        #define GMX_CUDA_MAX_BLOCKS_PER_MP   0
        #define GMX_CUDA_MAX_THREADS_PER_MP  0
#endif

#endif /* CUDA_ARCH_UTILS_CUH_ */
