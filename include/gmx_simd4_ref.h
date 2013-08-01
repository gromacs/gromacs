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

#ifndef _gmx_simd4_ref_h_
#define _gmx_simd4_ref_h_

/* This file contains a reference plain-C implementation of 4-wide SIMD.
 * This code is only useful for testing and documentation.
 * Either float or double precision is supported through gmx_simd4_real,
 * which is set in gmx_simd4_macros.h
 */


#include <math.h>

/* float/double SIMD register type */
typedef struct {
    gmx_simd4_real r[GMX_SIMD4_WIDTH];
} gmx_simd4_ref_pr;

/* boolean SIMD register type */
typedef struct {
    char r[GMX_SIMD4_WIDTH];
} gmx_simd4_ref_pb;


/* Load GMX_SIMD4_WIDTH reals for memory starting at r */
static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_load_pr(const gmx_simd4_real *r)
{
    gmx_simd4_ref_pr a;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        a.r[i] = r[i];
    }

    return a;
}

/* Set all SIMD register elements to r */
static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_set1_pr(gmx_simd4_real r)
{
    gmx_simd4_ref_pr a;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        a.r[i] = r;
    }

    return a;
}

/* Set all SIMD register elements to 0 */
static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_setzero_pr()
{
    gmx_simd4_ref_pr a;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        a.r[i] = 0.0;
    }

    return a;
}

static gmx_inline void
gmx_simd4_ref_store_pr(gmx_simd4_real *dest, gmx_simd4_ref_pr src)
{
    int i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        dest[i] = src.r[i];
    }
}

static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_add_pr(gmx_simd4_ref_pr a, gmx_simd4_ref_pr b)
{
    gmx_simd4_ref_pr c;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }

    return c;
}

static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_sub_pr(gmx_simd4_ref_pr a, gmx_simd4_ref_pr b)
{
    gmx_simd4_ref_pr c;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }

    return c;
}

static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_mul_pr(gmx_simd4_ref_pr a, gmx_simd4_ref_pr b)
{
    gmx_simd4_ref_pr c;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = a.r[i]*b.r[i];
    }

    return c;
}

static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_madd_pr(gmx_simd4_ref_pr a, gmx_simd4_ref_pr b, gmx_simd4_ref_pr c)
{
    gmx_simd4_ref_pr d;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        d.r[i] = a.r[i]*b.r[i] + c.r[i];
    }

    return d;
}

static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_nmsub_pr(gmx_simd4_ref_pr a, gmx_simd4_ref_pr b, gmx_simd4_ref_pr c)
{
    gmx_simd4_ref_pr d;
    int             i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        d.r[i] = -a.r[i]*b.r[i] + c.r[i];
    }

    return d;
}

static gmx_inline gmx_simd4_real
gmx_simd4_ref_dotproduct3(gmx_simd4_ref_pr a, gmx_simd4_ref_pr b)
{
    gmx_simd4_real dp;
    int            i;

    dp = 0.0;
    for (i = 0; i < 3; i++)
    {
        dp += a.r[i]*b.r[i];
    }

    return dp;
}

static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_min_pr(gmx_simd4_ref_pr a, gmx_simd4_ref_pr b)
{
    gmx_simd4_ref_pr c;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = (a.r[i] <= b.r[i] ? a.r[i] : b.r[i]);
    }

    return c;
}

static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_max_pr(gmx_simd4_ref_pr a, gmx_simd4_ref_pr b)
{
    gmx_simd4_ref_pr c;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = (a.r[i] >= b.r[i] ? a.r[i] : b.r[i]);
    }

    return c;
}

static gmx_inline gmx_simd4_ref_pr
gmx_simd4_ref_blendzero_pr(gmx_simd4_ref_pr a, gmx_simd4_ref_pb b)
{
    gmx_simd4_ref_pr c;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = (b.r[i] ? a.r[i] : 0.0);
    }

    return c;
}

/* Comparison */
static gmx_inline gmx_simd4_ref_pb
gmx_simd4_ref_cmplt_pr(gmx_simd4_ref_pr a, gmx_simd4_ref_pr b)
{
    gmx_simd4_ref_pb c;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = (a.r[i] < b.r[i]);
    }

    return c;
}

/* Logical AND on SIMD booleans */
static gmx_inline gmx_simd4_ref_pb
gmx_simd4_ref_and_pb(gmx_simd4_ref_pb a, gmx_simd4_ref_pb b)
{
    gmx_simd4_ref_pb c;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = (a.r[i] && b.r[i]);
    }

    return c;
}

/* Logical OR on SIMD booleans */
static gmx_inline gmx_simd4_ref_pb
gmx_simd4_ref_or_pb(gmx_simd4_ref_pb a, gmx_simd4_ref_pb b)
{
    gmx_simd4_ref_pb c;
    int              i;

    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = (a.r[i] || b.r[i]);
    }

    return c;
}

/* gmx_anytrue_pb(x) returns if any of the boolean is x is True */
static gmx_inline int
gmx_simd4_ref_anytrue_pb(gmx_simd4_ref_pb a)
{
    int anytrue;
    int i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        if (a.r[i])
        {
            anytrue = 1;
        }
    }

    return anytrue;
}

#endif /* _gmx_simd4_ref_h_ */
