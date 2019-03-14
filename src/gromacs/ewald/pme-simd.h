/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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
#ifndef GMX_EWALD_PME_SIMD_H
#define GMX_EWALD_PME_SIMD_H

/* Include the SIMD macro file and then check for support */
#include "gromacs/simd/simd.h"

/* Check if we have 4-wide SIMD macro support */
#if GMX_SIMD4_HAVE_REAL
/* Do PME spread and gather with 4-wide SIMD.
 * NOTE: SIMD is only used with PME order 4 and 5 (which are the most common).
 */
#    define PME_SIMD4_SPREAD_GATHER

#    if GMX_SIMD_HAVE_LOADU && GMX_SIMD_HAVE_STOREU
/* With PME-order=4 on x86, unaligned load+store is slightly faster
 * than doubling all SIMD operations when using aligned load+store.
 */
#        define PME_SIMD4_UNALIGNED
#    endif
#endif

#ifdef PME_SIMD4_SPREAD_GATHER
#    define SIMD4_ALIGNMENT  (GMX_SIMD4_WIDTH*sizeof(real))
#else
/* We can use any alignment, apart from 0, so we use 4 reals */
#    define SIMD4_ALIGNMENT  (4*sizeof(real))
#endif

/* Check if we can use SIMD with packs of 4 for gather with order 4 */
#if GMX_SIMD_HAVE_4NSIMD_UTIL_REAL && GMX_SIMD_REAL_WIDTH <= 16
#    define PME_4NSIMD_GATHER  1
#else
#    define PME_4NSIMD_GATHER  0
#endif

#endif
