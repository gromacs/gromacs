/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
#ifndef CUDA_ARCH_UTILS_CUH_
#define CUDA_ARCH_UTILS_CUH_

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
#    define GMX_PTX_ARCH 0
#else
#    define GMX_PTX_ARCH __CUDA_ARCH__
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

/*! \brief Allow disabling CUDA textures using the GMX_DISABLE_CUDA_TEXTURES macro.
 *
 *  Only texture objects supported.
 *  Disable texture support missing in clang (all versions up to <=5.0-dev as of writing).
 *  Disable texture support on CC 7.0 and 8.0 for performance reasons (Issue #3845).
 *
 *  This option will not influence functionality. All features using textures ought
 *  to have fallback for texture-less reads (direct/LDG loads), all new code needs
 *  to provide fallback code.
 */
#if defined(GMX_DISABLE_CUDA_TEXTURES) || (defined(__clang__) && defined(__CUDA__)) \
        || (GMX_PTX_ARCH == 700) || (GMX_PTX_ARCH >= 800)
#    define DISABLE_CUDA_TEXTURES 1
#else
#    define DISABLE_CUDA_TEXTURES 0
#endif

/*! \brief True if the use of texture fetch in the CUDA kernels is disabled. */
static const bool c_disableCudaTextures = DISABLE_CUDA_TEXTURES;


/* CUDA architecture technical characteristics. Needs macros because it is used
 * in the __launch_bounds__ function qualifiers and might need it in preprocessor
 * conditionals.
 *
 */
#if GMX_PTX_ARCH > 0
#    if GMX_PTX_ARCH <= 370 // CC 3.x
#        define GMX_CUDA_MAX_BLOCKS_PER_MP 16
#        define GMX_CUDA_MAX_THREADS_PER_MP 2048
#    elif GMX_PTX_ARCH == 750 // CC 7.5, lower limits compared to 7.0
#        define GMX_CUDA_MAX_BLOCKS_PER_MP 16
#        define GMX_CUDA_MAX_THREADS_PER_MP 1024
#    elif (GMX_PTX_ARCH == 860 || GMX_PTX_ARCH == 890) // CC 8.6,8.9 lower limits compared to 8.0
#        define GMX_CUDA_MAX_BLOCKS_PER_MP 16
#        define GMX_CUDA_MAX_THREADS_PER_MP 1536
#    else // CC 5.x, 6.x, 7.0, 8.0
/* Note that this final branch covers all future architectures (current gen
 * is 8.x as of writing), hence assuming that these *currently defined* upper
 * limits will not be lowered.
 */
#        define GMX_CUDA_MAX_BLOCKS_PER_MP 32
#        define GMX_CUDA_MAX_THREADS_PER_MP 2048
#    endif
#else
#    define GMX_CUDA_MAX_BLOCKS_PER_MP 0
#    define GMX_CUDA_MAX_THREADS_PER_MP 0
#endif

// Macro defined for clang CUDA device compilation in the presence of debug symbols
// used to work around codegen bug that breaks some kernels when assertions are on
// at -O1 and higher (tested with clang 6-8).
#if defined(__clang__) && defined(__CUDA__) && defined(__CUDA_ARCH__) && !defined(NDEBUG)
#    define CLANG_DISABLE_OPTIMIZATION_ATTRIBUTE __attribute__((optnone))
#else
#    define CLANG_DISABLE_OPTIMIZATION_ATTRIBUTE
#endif


#endif /* CUDA_ARCH_UTILS_CUH_ */
