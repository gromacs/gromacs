/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "utilities.h"

#include "config.h"

#include <assert.h>
#include <limits.h>

#include <cmath>

#include <algorithm>

#if !GMX_NATIVE_WINDOWS
// for fp exception control stuff
#include <fenv.h>
#endif

int
gmx_within_tol(double   f1,
               double   f2,
               double   tol)
{
    /* The or-equal is important - otherwise we return false if f1==f2==0 */
    if (fabs(f1-f2) <= tol*0.5*(fabs(f1)+fabs(f2)) )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int
gmx_numzero(double a)
{
    return gmx_within_tol(a, 0.0, GMX_REAL_MIN/GMX_REAL_EPS);
}


gmx_bool
check_int_multiply_for_overflow(gmx_int64_t  a,
                                gmx_int64_t  b,
                                gmx_int64_t *result)
{
    gmx_int64_t sign = 1;
    if ((0 == a) || (0 == b))
    {
        *result = 0;
        return TRUE;
    }
    if (a < 0)
    {
        a    = -a;
        sign = -sign;
    }
    if (b < 0)
    {
        b    = -b;
        sign = -sign;
    }
    if (GMX_INT64_MAX / b < a)
    {
        *result = (sign > 0) ? GMX_INT64_MAX : GMX_INT64_MIN;
        return FALSE;
    }
    *result = sign * a * b;
    return TRUE;
}

int gmx_greatest_common_divisor(int p, int q)
{
    int tmp;
    while (q != 0)
    {
        tmp = q;
        q   = p % q;
        p   = tmp;
    }
    return p;
}

int gmx_feenableexcept()
{
#if HAVE_FEENABLEEXCEPT
    return feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#elif (defined(__i386__) || defined(__x86_64__)) && defined(__APPLE__)
    /* Author:  David N. Williams
     * License:  Public Domain
     *
     * Might also work on non-Apple Unix. But should be tested
     * before enabling.
     */
    unsigned int  excepts = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW;
    static fenv_t fenv;
    unsigned int  new_excepts = excepts & FE_ALL_EXCEPT,
                  old_excepts; // previous masks

    if (fegetenv (&fenv) )
    {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // unmask
    fenv.__control &= ~new_excepts;
    fenv.__mxcsr   &= ~(new_excepts << 7);

    return ( fesetenv (&fenv) ? -1 : old_excepts );
#else
    return -1;
#endif
}

real max_cutoff(real cutoff1, real cutoff2)
{
    if (cutoff1 == 0 || cutoff2 == 0)
    {
        return 0;
    }
    else
    {
        return std::max(cutoff1, cutoff2);
    }
}
