/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013 by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

/*! \file
 * \brief
 * The types in this file are intended to be used for supporting
 * writing new architecture-independent SIMD intrinsics code. To
 * support a new architecture, defining types equivalent to those here
 * and and implementing functions that use them should be (nearly) all
 * that is needed.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_simd
 */

#ifndef _gmx_simd_reference_types_h_
#define _gmx_simd_reference_types_h_

#include "typedefs.h"


#ifndef GMX_SIMD_REFERENCE_PLAIN_C
#error Reference NBNXN kernel include file used without setting GMX_SIMD_REFERENCE_PLAIN_C
#endif

#define GMX_HAVE_SIMD_MACROS

/* In general the reference SIMD supports any SIMD width, including 1.
 * For the nbnxn 4xn kernels all widths (2, 4 and 8) are supported.
 * The nbnxn 2xnn kernels are currently not supported.
 */
#ifdef GMX_SIMD_WIDTH_HERE
/* This must be test code that has already #included the real SIMD
 * code that is being tested, so copy its SIMD width. */
#define GMX_SIMD_REF_WIDTH GMX_SIMD_WIDTH_HERE
#else
/* This must be someone running the reference SIMD implementation
 * manually. They can modify this to match what they want to do. */
#define GMX_SIMD_REF_WIDTH  4
#define GMX_SIMD_WIDTH_HERE  GMX_SIMD_REF_WIDTH
#endif

/* float/double SIMD register type */
typedef struct {
    real r[GMX_SIMD_REF_WIDTH];
} gmx_simd_ref_pr;

/* boolean SIMD register type */
typedef struct {
    char r[GMX_SIMD_REF_WIDTH];
} gmx_simd_ref_pb;

/* integer SIMD register type, only for table indexing and exclusion masks */
typedef struct {
    int r[GMX_SIMD_REF_WIDTH];
} gmx_simd_ref_epi32;
#define GMX_SIMD_REF_EPI32_WIDTH  GMX_SIMD_REF_WIDTH

/* unsigned integer SIMD register type, only for exclusion masks */
typedef struct {
    unsigned int r[GMX_SIMD_REF_WIDTH];
} gmx_simd_ref_exclmask;

#endif /* _gmx_simd_reference_types_h_ */
