/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "xdrf.h"
#include "gmx_fatal.h"
#include "smalloc.h"

int xdr_real(XDR *xdrs, real *r)
{
#ifdef GMX_DOUBLE
    float f;
    int   ret;

    f   = *r;
    ret = xdr_float(xdrs, &f);
    *r  = f;

    return ret;
#else
    return xdr_float(xdrs, (float *)r);
#endif
}

int xdr3drcoord(XDR *xdrs, real *fp, int *size, real *precision)
{
#ifdef GMX_DOUBLE
    float *ffp;
    float  fprec;
    int    i, ret, isize;

    isize = *size*DIM;
    if (isize <= 0)
    {
        gmx_fatal(FARGS, "Don't know what to malloc for ffp, isize = %d", isize);
    }

    snew(ffp, isize);

    for (i = 0; (i < isize); i++)
    {
        ffp[i] = fp[i];
    }
    fprec = *precision;
    ret   = xdr3dfcoord(xdrs, ffp, size, &fprec);

    *precision = fprec;
    for (i = 0; (i < isize); i++)
    {
        fp[i] = ffp[i];
    }

    sfree(ffp);
    return ret;
#else
    return xdr3dfcoord(xdrs, (float *)fp, size, (float *)precision);
#endif
}

int xdr_gmx_large_int(XDR *xdrs, gmx_large_int_t *i, const char *warn)
{
    /* This routine stores values compatible with xdr_int64_t */

    int imaj, imin;
    int ret;

#if ((defined SIZEOF_GMX_LARGE_INT) && SIZEOF_GMX_LARGE_INT == 8)
    static const gmx_large_int_t two_p32_m1 = 0xFFFFFFFF;
    gmx_large_int_t              imaj64, imin64;

    imaj64 = ((*i)>>32) & two_p32_m1;
    imin64 = (*i) & two_p32_m1;
    imaj   = (int)imaj64;
    imin   = (int)imin64;
#else
    /* Our code has 4 bytes, but we should make sure that this value
     * will be correctly read by 8 byte code.
     */
    if (*i >= 0)
    {
        imaj = 0;
    }
    else
    {
        imaj = -1;
    }
    imin = *i;
#endif
    ret = xdr_int(xdrs, &imaj);
    ret = xdr_int(xdrs, &imin);

#if ((defined SIZEOF_GMX_LARGE_INT) && SIZEOF_GMX_LARGE_INT == 8)
    *i = (((gmx_large_int_t)imaj << 32) | ((gmx_large_int_t)imin & two_p32_m1));
#else
    *i = imin;

    if (warn != NULL && (imaj < -1 || imaj > 0))
    {
        fprintf(stderr, "\nWARNING during %s:\n", warn);
        fprintf(stderr, "a step value written by code supporting 64bit integers is read by code that only supports 32bit integers, out of range step value has been converted to %d\n\n", *i);
    }
#endif

    return ret;
}
