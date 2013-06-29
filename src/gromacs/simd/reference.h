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

#include "reference_types.h"
#include "typedefs.h"

/* This file contains a reference plain-C implementation of arbitrary width.
 * This code is only useful for testing and documentation.
 * The SIMD width is set by defining GMX_SIMD_REF_WIDTH before including.
 */


#ifndef GMX_SIMD_REF_WIDTH
#error "GMX_SIMD_REF_WIDTH should be defined before including gromacs/simd/reference.h"
#endif

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

static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_round_pr(gmx_simd_ref_pr a)
{
    gmx_simd_ref_pr b;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        b.r[i] = (real)(int)(a.r[i] + ((a.r[i] >= 0) ? 0.5 : -0.5));
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
        b.r[i] = (real)(int)((a.r[i] >= 0) ? a.r[i] : a.r[i] - 1);
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

/* Logical AND on SIMD booleans */
static gmx_inline gmx_simd_ref_pb
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

/* Logical OR on SIMD booleans */
static gmx_inline gmx_simd_ref_pb
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

/* Not required, gmx_anytrue_pb(x) returns if any of the boolean is x is True.
 * If this is not present, define GMX_SIMD_IS_TRUE(real x),
 * which should return x==True, where True is True as defined in SIMD.
 */
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

/* If we don't have gmx_anytrue_pb, we need to store gmx_mm_pb */
static gmx_inline void
gmx_simd_ref_store_pb(real *dest, gmx_simd_ref_pb src)
{
    int i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        dest[i] = src.r[i];
    }
};

static gmx_inline gmx_simd_ref_exclmask
gmx_simd_ref_load1_exclmask(unsigned int src)
{
    gmx_simd_ref_exclmask a;
    int                   i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        a.r[i] = src;
    }

    return a;
}

static gmx_inline gmx_simd_ref_exclmask
gmx_simd_ref_load_exclmask(const unsigned int *src)
{
    gmx_simd_ref_exclmask a;
    int                   i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        a.r[i] = src[i];
    }

    return a;
}

/* For topology exclusion-pair checking we need: (anytrue(a & b) ?
 * False : True). The x86 implementations use hardware-suitable
 * integer- and/or real-valued SIMD operations and a bit-wise "and" to
 * do this. The reference implementation normally uses logical
 * operations for logic, but in this case the i- and j-atom exclusion
 * masks computed during searching expect to be combined with bit-wise
 * "and," so we do that.
 *
 * If the same bit is set in both input masks, return TRUE, else
 * FALSE. This function is only called with a single bit set in b.
 */
static gmx_inline gmx_simd_ref_pb
gmx_simd_ref_checkbitmask_pb(gmx_simd_ref_exclmask a, gmx_simd_ref_exclmask b)
{
    gmx_simd_ref_pb c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = ((a.r[i] & b.r[i]) != 0);
    }

    return c;
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
        b.r[i] = 1.0/sqrt(a.r[i]);
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
        b.r[i] = exp(a.r[i]);
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
        b.r[i] = sqrt(a.r[i]);
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
