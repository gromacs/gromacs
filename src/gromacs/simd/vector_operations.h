/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

/* The macros in this file are intended to be used for writing
 * architecture-independent SIMD intrinsics code.
 * To support a new architecture, adding macros here should be (nearly)
 * all that is needed.
 */

/* This file contains vector operation functions using SIMD intrinsics.
 * gromacs/simd/macros.h should be included before including this file.
 */

#ifndef GMX_SIMD_VECTOR_OPERATIONS_H
#define GMX_SIMD_VECTOR_OPERATIONS_H

#ifndef GMX_SIMD_MACROS_H
#error "gromacs/simd/macros.h was not included before including gromacs/simd/vector_operations.h"
#endif


/* x^2 + y^2 + z^2 */
static gmx_inline gmx_simd_real_t
gmx_simd_calc_rsq_r(gmx_simd_real_t x, gmx_simd_real_t y, gmx_simd_real_t z)
{
    return gmx_simd_fmadd_r(z, z, gmx_simd_fmadd_r(y, y, gmx_simd_mul_r(x, x)));
}

/* inner-product of multiple vectors */
static gmx_inline gmx_simd_real_t
gmx_simd_iprod_r(gmx_simd_real_t ax, gmx_simd_real_t ay, gmx_simd_real_t az,
                 gmx_simd_real_t bx, gmx_simd_real_t by, gmx_simd_real_t bz)
{
    gmx_simd_real_t ret;

    ret = gmx_simd_mul_r(ax, bx);
    ret = gmx_simd_fmadd_r(ay, by, ret);
    ret = gmx_simd_fmadd_r(az, bz, ret);

    return ret;
}

/* norm squared of multiple vectors */
static gmx_inline gmx_simd_real_t
gmx_simd_norm2_r(gmx_simd_real_t ax, gmx_simd_real_t ay, gmx_simd_real_t az)
{
    gmx_simd_real_t ret;

    ret = gmx_simd_mul_r(ax, ax);
    ret = gmx_simd_fmadd_r(ay, ay, ret);
    ret = gmx_simd_fmadd_r(az, az, ret);

    return ret;
}

/* cross-product of multiple vectors */
static gmx_inline void
gmx_simd_cprod_r(gmx_simd_real_t ax, gmx_simd_real_t ay, gmx_simd_real_t az,
                 gmx_simd_real_t bx, gmx_simd_real_t by, gmx_simd_real_t bz,
                 gmx_simd_real_t *cx, gmx_simd_real_t *cy, gmx_simd_real_t *cz)
{
    *cx = gmx_simd_mul_r(ay, bz);
    *cx = gmx_simd_fnmadd_r(az, by, *cx);

    *cy = gmx_simd_mul_r(az, bx);
    *cy = gmx_simd_fnmadd_r(ax, bz, *cy);

    *cz = gmx_simd_mul_r(ax, by);
    *cz = gmx_simd_fnmadd_r(ay, bx, *cz);
}

/* a + b + c + d (not really a vector operation, but where else put this?) */
static gmx_inline gmx_simd_real_t
gmx_simd_sum4_r(gmx_simd_real_t a, gmx_simd_real_t b, gmx_simd_real_t c, gmx_simd_real_t d)
{
    return gmx_simd_add_r(gmx_simd_add_r(a, b), gmx_simd_add_r(c, d));
}


#endif
