/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "gromacs/fileio/xdrf.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

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

int xdr_int64(XDR *xdrs, gmx_int64_t *i)
{
    /* This routine stores values compatible with xdr_int64_t */

    int imaj, imin;
    int ret;

    static const gmx_int64_t two_p32_m1 = 0xFFFFFFFF;
    gmx_int64_t              imaj64, imin64;

    imaj64 = ((*i)>>32) & two_p32_m1;
    imin64 = (*i) & two_p32_m1;
    imaj   = (int)imaj64;
    imin   = (int)imin64;
    ret    = xdr_int(xdrs, &imaj);
    ret   |= xdr_int(xdrs, &imin);

    *i = (((gmx_int64_t)imaj << 32) | ((gmx_int64_t)imin & two_p32_m1));

    return ret;
}
