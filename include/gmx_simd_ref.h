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

#ifndef _gmx_simd_ref_h_
#define _gmx_simd_ref_h_

/* This file contains a reference plain-C implementation of arbitrary width.
 * This code is only useful for testing and documentation.
 * The SIMD width is set by defining GMX_SIMD_REF_WIDTH before including.
 */


#ifndef GMX_SIMD_REF_WIDTH
#error "GMX_SIMD_REF_WIDTH should be defined before including gmx_simd_ref.h"
#endif

#include <math.h>

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

/* Load GMX_SIMD_REF_WIDTH reals for memory starting at r */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_load_pr(const real *r)
{
    gmx_simd_ref_pr a;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        a.r[i] = r[i];
    }

    return a;
}

/* Set all SIMD register elements to *r */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_load1_pr(const real *r)
{
    gmx_simd_ref_pr a;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        a.r[i] = *r;
    }

    return a;
}

/* Set all SIMD register elements to r */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_set1_pr(real r)
{
    gmx_simd_ref_pr a;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        a.r[i] = r;
    }

    return a;
}

/* Set all SIMD register elements to 0 */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_setzero_pr()
{
    gmx_simd_ref_pr a;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        a.r[i] = 0.0;
    }

    return a;
}

static gmx_inline void
gmx_simd_ref_store_pr(real *dest, gmx_simd_ref_pr src)
{
    int i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        dest[i] = src.r[i];
    }
}

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_add_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }

    return c;
}

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_sub_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }

    return c;
}

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_mul_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = a.r[i]*b.r[i];
    }

    return c;
}

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_madd_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b, gmx_simd_ref_pr c)
{
    gmx_simd_ref_pr d;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        d.r[i] = a.r[i]*b.r[i] + c.r[i];
    }

    return d;
}

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_nmsub_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b, gmx_simd_ref_pr c)
{
    gmx_simd_ref_pr d;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        d.r[i] = -a.r[i]*b.r[i] + c.r[i];
    }

    return d;
}

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_max_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = (a.r[i] >= b.r[i] ? a.r[i] : b.r[i]);
    }

    return c;
}

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_blendzero_pr(gmx_simd_ref_pr a, gmx_simd_ref_pb b)
{
    gmx_simd_ref_pr c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = (b.r[i] ? a.r[i] : 0.0);
    }

    return c;
}

/* Note that this reference implementation rounds away from zero,
 * whereas most SIMD intrinsics will round to nearest even. Since this
 * function is only used for periodic image calculations, the rounding
 * of mantissas close to 0.5 is irrelevant, except in testing. This
 * could be fixed by using rint/rintf, but the bigger problem is that
 * MSVC does not support full C99, and none of the round or rint
 * functions are defined. It's much easier to approximately implement
 * round() than rint(), so we do that and hope we never get bitten in
 * testing. (Thanks, Microsoft.)
 */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_round_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
#ifdef _MSC_VER
        int temp = (a.r[i] >= 0.)
            ? (a.r[i] + 0.5)
            : (a.r[i] - 0.5);
        b.r[i] = (real) temp;
#elif defined GMX_DOUBLE
        b.r[i] = round(a.r[i]);
#else
        b.r[i] = roundf(a.r[i]);
#endif
    }

    return b;
}

/* Not required, only used to speed up the nbnxn tabulated PME kernels */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_floor_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
#ifdef GMX_DOUBLE
        b.r[i] = floor(a.r[i]);
#else
        b.r[i] = floorf(a.r[i]);
#endif
    }

    return b;
}

/* Not required, only used when blendv is faster than comparison */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_blendv_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b, gmx_simd_ref_pr c)
{
    gmx_simd_ref_pr d;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        d.r[i] = (c.r[i] >= 0) ? a.r[i] : b.r[i];
    }

    return d;
}

/* Copy the sign of a to b, assumes b >= 0 for efficiency */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_cpsgn_nonneg_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = (a.r[i] >= 0) ? b.r[i] : -b.r[i];
    }

    return c;
}

/* Very specific operation required in the non-bonded kernels */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_masknot_add_pr(gmx_simd_ref_pb a, gmx_simd_ref_pr b, gmx_simd_ref_pr c)
{
    gmx_simd_ref_pr d;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        d.r[i] = a.r[i] ? b.r[i] : b.r[i] + c.r[i];
    }

    return d;
}

/* Comparison */
static gmx_inline gmx_simd_ref_pb
gmx_simd_ref_cmplt_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pb c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = (a.r[i] < b.r[i]);
    }

    return c;
}

/* Logical AND on SIMD booleans. Can't be static or it can't be a
   template parameter (at least on XLC for BlueGene/Q) */
gmx_inline gmx_simd_ref_pb
gmx_simd_ref_and_pb(gmx_simd_ref_pb a, gmx_simd_ref_pb b)
{
    gmx_simd_ref_pb c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = (a.r[i] && b.r[i]);
    }

    return c;
}

/* Logical OR on SIMD booleans. Can't be static or it can't be a
   template parameter (at least on XLC for BlueGene/Q) */
gmx_inline gmx_simd_ref_pb
gmx_simd_ref_or_pb(gmx_simd_ref_pb a, gmx_simd_ref_pb b)
{
    gmx_simd_ref_pb c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = (a.r[i] || b.r[i]);
    }

    return c;
}

/* Returns a single int (0/1) which tells if any of the booleans is True */
static gmx_inline int
gmx_simd_ref_anytrue_pb(gmx_simd_ref_pb a)
{
    int anytrue;
    int i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        if (a.r[i])
        {
            anytrue = 1;
        }
    }

    return anytrue;
}

/* Conversions only used for PME table lookup */
static gmx_inline gmx_simd_ref_epi32
gmx_simd_ref_cvttpr_epi32(gmx_simd_ref_pr a)
{
    gmx_simd_ref_epi32 b;
    int                i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        b.r[i] = (int)a.r[i];
    }

    return b;
};

/* These two function only need to be approximate, Newton-Raphson iteration
 * is used for full accuracy in gmx_invsqrt_pr and gmx_inv_pr.
 */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_rsqrt_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
#ifdef GMX_DOUBLE
        b.r[i] = 1.0/sqrt(a.r[i]);
#else
        b.r[i] = 1.0/sqrtf(a.r[i]);
#endif
    }

    return b;
};

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_rcp_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        b.r[i] = 1.0/a.r[i];
    }

    return b;
};

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_exp_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
#ifdef GMX_DOUBLE
        b.r[i] = exp(a.r[i]);
#else
        b.r[i] = expf(a.r[i]);
#endif
    }

    return b;
};

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_sqrt_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
#ifdef GMX_DOUBLE
        b.r[i] = sqrt(a.r[i]);
#else
        b.r[i] = sqrtf(a.r[i]);
#endif
    }

    return b;
}

static gmx_inline int
gmx_simd_ref_sincos_pr(gmx_simd_ref_pr a,
                       gmx_simd_ref_pr *s, gmx_simd_ref_pr *c)
{
    int i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        s->r[i] = sin(a.r[i]);
        c->r[i] = cos(a.r[i]);
    }

    return 0;
}

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_acos_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        b.r[i] = acos(a.r[i]);
    }

    return b;
}

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_atan2_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = atan2(a.r[i], b.r[i]);
    }

    return c;
}

#endif /* _gmx_simd_ref_h_ */
