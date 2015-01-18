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

#ifndef _nbnxn_simd_h
#define _nbnxn_simd_h

#include "config.h"

#include "gromacs/legacyheaders/typedefs.h"

/* Include SIMD, below we select kernels based on the SIMD width */
#include "gromacs/simd/simd.h"

#ifdef GMX_SIMD_REFERENCE
#define GMX_NBNXN_SIMD
#endif

/* As we modularize the verlet kernels, we should remove stuff like this
 * that checks internal SIMD implementation details.
 */
#if (defined GMX_SIMD_X86_SSE2) || (defined GMX_SIMD_X86_SSE4_1) || \
    (defined GMX_SIMD_X86_AVX_128_FMA) || (defined GMX_SIMD_X86_AVX_256) || \
    (defined GMX_SIMD_X86_AVX2_256) || (defined GMX_SIMD_IBM_QPX)
/* Use SIMD accelerated nbnxn search and kernels */
#define GMX_NBNXN_SIMD
#endif

/* MIC for double is implemented in the SIMD module but so far missing in
   mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_x86_mic.h */
#if defined GMX_SIMD_X86_MIC && !defined GMX_DOUBLE
#define GMX_NBNXN_SIMD
#endif

#ifdef GMX_NBNXN_SIMD
/* The nbnxn SIMD 4xN and 2x(N+N) kernels can be added independently.
 * Currently the 2xNN SIMD kernels only make sense with:
 *  8-way SIMD: 4x4 setup, works with AVX-256 in single precision
 * 16-way SIMD: 4x8 setup, works with Intel MIC in single precision
 */
#if GMX_SIMD_REAL_WIDTH == 2 || GMX_SIMD_REAL_WIDTH == 4 || GMX_SIMD_REAL_WIDTH == 8
#define GMX_NBNXN_SIMD_4XN
#endif
#if GMX_SIMD_REAL_WIDTH == 8 || GMX_SIMD_REAL_WIDTH == 16
#define GMX_NBNXN_SIMD_2XNN
#endif

#if !(defined GMX_NBNXN_SIMD_4XN || defined GMX_NBNXN_SIMD_2XNN)
#error "No SIMD kernel type defined"
#endif

#endif /* GMX_NBNXN_SIMD */

#endif /* _nbnxn_simd_h */
