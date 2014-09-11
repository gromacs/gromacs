/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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
#ifndef GMX_MATH_3DTRANSFORMS_H
#define GMX_MATH_3DTRANSFORMS_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Index for the fourth dimension for `vec4`. */
#define WW 3

/*! \brief
 * 4D vector type used in 3D transformations.
 *
 * In \Gromacs, only a limited set of 3D transformations are used, and all of
 * them operate on coordinates, so the fourth element is assumed to be one and
 * ignored in all contexts.
 */
typedef real vec4[4];

/*! \brief
 * 4D matrix type used in 3D transformations.
 */
typedef real mat4[4][4];

void gmx_mat4_copy(mat4 a, mat4 b);

void gmx_mat4_transform_point(mat4 m, rvec x, vec4 v);

/*! \brief
 * Computes the product of two `mat4` matrices as A = B * C.
 *
 * Note that the order of operands is different from mmul() in vec.h!
 */
void gmx_mat4_mmul(mat4 A, mat4 B, mat4 C);

void gmx_mat4_init_unity(mat4 m);

void gmx_mat4_init_rotation(int axis, real angle, mat4 A);

void gmx_mat4_init_translation(real tx, real ty, real tz, mat4 A);

void gmx_mat4_print(FILE *fp, const char *s, mat4 A);

void gmx_vec4_print(FILE *fp, const char *s, vec4 a);

#ifdef __cplusplus
}
#endif

#endif
