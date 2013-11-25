/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_BITMASK_H /* TODO: fix if not used with Ia2b9030324 */
#define GMX_MDLIB_BITMASK_H

#include <string.h>

#include "gromacs/utility/basedefinitions.h"

#define BITMASK_SIZE 1 /* number of bits is 64*BITMASK_SIZE */

#if BITMASK_SIZE == 1
typedef gmx_uint64_t gmx_bitmask_t;

gmx_inline static void bitmask_clear(gmx_bitmask_t* m)
{
    *m = 0;
}

gmx_inline static void bitmask_set_bit(gmx_bitmask_t* m, int b)
{
    *m |= (1ULL << b);
}

gmx_inline static void bitmask_init_bit(gmx_bitmask_t* m, int b)
{
    *m = (1ULL << b);
}

gmx_inline static void bitmask_init_low_bits(gmx_bitmask_t* m, int b)
{
    *m = (1ULL << b) - 1;
}

gmx_inline static gmx_bool bitmask_is_set(gmx_bitmask_t m, int b)
{
    return (m & (1ULL << b)) != 0;
}

gmx_inline static gmx_bool bitmask_is_disjoint(gmx_bitmask_t a, gmx_bitmask_t b)
{
    return !(a & b);
}

gmx_inline static gmx_bool bitmask_is_equal(gmx_bitmask_t a, gmx_bitmask_t b)
{
    return a == b;
}

gmx_inline static gmx_bool bitmask_is_zero(gmx_bitmask_t m)
{
    return !m;
}

gmx_inline static void bitmask_union(gmx_bitmask_t* a, gmx_bitmask_t b)
{
    *a |= b;
}
#else
typedef gmx_uint64_t gmx_bitmask_t[BITMASK_SIZE];

gmx_inline static void bitmask_clear(gmx_bitmask_t* m)
{
    memset(m, 0, 8*BITMASK_SIZE);
}

gmx_inline static void bitmask_set_bit(gmx_bitmask_t* m, int b)
{
    (*m)[b/64] |= (1ULL << (b%64));
}

gmx_inline static void bitmask_init_bit(gmx_bitmask_t* m, int b)
{
    bitmask_clear(m);
    (*m)[b/64] = (1ULL << (b%64));
}

gmx_inline static void bitmask_init_low_bits(gmx_bitmask_t* m, int b)
{
    memset(m, 255, b/64*8);
    (*m)[b/64] = (1ULL << (b%64)) - 1;
    memset(m+b/64+64, 0, (8-b/64-1)*8);
}

gmx_inline static gmx_bool bitmask_is_set(gmx_bitmask_t m, int b)
{
    return (m[b/64] & (1ULL << (b%64))) != 0;
}

gmx_inline static gmx_bool bitmask_is_disjoint(gmx_bitmask_t a, gmx_bitmask_t b)
{
    int      i;
    gmx_bool r = 1;
    for (i = 0; i < BITMASK_SIZE; i++)
    {
        r = r && !(a[i] & b[i]);
    }
    return r;
}

gmx_inline static gmx_bool bitmask_is_equal(gmx_bitmask_t a, gmx_bitmask_t b)
{
    int      i;
    gmx_bool r = 1;
    for (i = 0; i < BITMASK_SIZE; i++)
    {
        r = r && (a[i] == b[i]);
    }
    return r;
}

gmx_inline static gmx_bool bitmask_is_zero(gmx_bitmask_t m)
{
    int      i;
    gmx_bool r = 1;
    for (i = 0; i < BITMASK_SIZE; i++)
    {
        r = r && !m[i];
    }
    return r;
}

gmx_inline static void bitmask_union(gmx_bitmask_t* a, gmx_bitmask_t b)
{
    int i;
    for (i = 0; i < BITMASK_SIZE; i++)
    {
        (*a)[i] |= b[i];
    }
}
#endif

#endif
