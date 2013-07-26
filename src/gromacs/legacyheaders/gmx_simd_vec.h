/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS Development Team
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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

/* The macros in this file are intended to be used for writing
 * architecture-independent SIMD intrinsics code.
 * To support a new architecture, adding macros here should be (nearly)
 * all that is needed.
 */

/* This file contains vector operation functions using SIMD intrinsics.
 * gmx_simd_macros.h should be included before including this file.
 */

#ifndef _gmx_simd_vec_h_
#define _gmx_simd_vec_h_

#ifndef _gmx_simd_macros_h_
#error "gmx_simd_macros.h was not included before including gmx_simd_vec.h"
#endif


/* x^2 + y^2 + z^2 */
static gmx_inline gmx_mm_pr
gmx_calc_rsq_pr(gmx_mm_pr x, gmx_mm_pr y, gmx_mm_pr z)
{
    return gmx_madd_pr(z, z, gmx_madd_pr(y, y, gmx_mul_pr(x, x)));
}

/* inner-product of multiple vectors */
static gmx_inline gmx_mm_pr
gmx_iprod_pr(gmx_mm_pr ax, gmx_mm_pr ay, gmx_mm_pr az,
             gmx_mm_pr bx, gmx_mm_pr by, gmx_mm_pr bz)
{
    gmx_mm_pr ret;

    ret = gmx_mul_pr(ax, bx);
    ret = gmx_madd_pr(ay, by, ret);
    ret = gmx_madd_pr(az, bz, ret);

    return ret;
}

/* norm squared of multiple vectors */
static gmx_inline gmx_mm_pr
gmx_norm2_pr(gmx_mm_pr ax, gmx_mm_pr ay, gmx_mm_pr az)
{
    gmx_mm_pr ret;

    ret = gmx_mul_pr(ax, ax);
    ret = gmx_madd_pr(ay, ay, ret);
    ret = gmx_madd_pr(az, az, ret);

    return ret;
}

/* cross-product of multiple vectors */
static gmx_inline void
gmx_cprod_pr(gmx_mm_pr ax, gmx_mm_pr ay, gmx_mm_pr az,
             gmx_mm_pr bx, gmx_mm_pr by, gmx_mm_pr bz,
             gmx_mm_pr *cx, gmx_mm_pr *cy, gmx_mm_pr *cz)
{
    *cx = gmx_mul_pr(ay, bz);
    *cx = gmx_nmsub_pr(az, by, *cx);

    *cy = gmx_mul_pr(az, bx);
    *cy = gmx_nmsub_pr(ax, bz, *cy);

    *cz = gmx_mul_pr(ax, by);
    *cz = gmx_nmsub_pr(ay, bx, *cz);
}

/* a + b + c + d (not really a vector operation, but where else put this?) */
static gmx_inline gmx_mm_pr
gmx_sum4_pr(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c, gmx_mm_pr d)
{
    return gmx_add_pr(gmx_add_pr(a, b), gmx_add_pr(c, d));
}


#endif /* _gmx_simd_vec_h_ */
