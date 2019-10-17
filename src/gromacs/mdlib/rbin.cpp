/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2015,2018,2019, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "rbin.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"

t_bin* mk_bin()
{
    t_bin* b;

    snew(b, 1);

    return b;
}

void destroy_bin(t_bin* b)
{
    if (b->maxreal > 0)
    {
        sfree(b->rbuf);
    }

    sfree(b);
}

void reset_bin(t_bin* b)
{
    b->nreal = 0;
}

int add_binr(t_bin* b, int nr, const real r[])
{
#define MULT 4
    int     i, rest, index;
    double* rbuf;

    if (b->nreal + nr > b->maxreal)
    {
        b->maxreal = b->nreal + nr;
        rest       = b->maxreal % MULT;
        if (rest != 0)
        {
            b->maxreal += MULT - rest;
        }
        srenew(b->rbuf, b->maxreal);
    }
    /* Copy pointer */
    rbuf = b->rbuf + b->nreal;

#if (defined __ICC && __ICC >= 1500 || defined __ICL && __ICL >= 1500) && defined __MIC__
#    pragma novector /* Work-around for incorrect vectorization */
#endif
    for (i = 0; (i < nr); i++)
    {
        rbuf[i] = r[i];
    }

    index = b->nreal;
    b->nreal += nr;

    return index;
}

int add_binr(t_bin* b, gmx::ArrayRef<const real> r)
{
    return add_binr(b, r.size(), r.data());
}

int add_bind(t_bin* b, int nr, const double r[])
{
#define MULT 4
    int     i, rest, index;
    double* rbuf;

    if (b->nreal + nr > b->maxreal)
    {
        b->maxreal = b->nreal + nr;
        rest       = b->maxreal % MULT;
        if (rest != 0)
        {
            b->maxreal += MULT - rest;
        }
        srenew(b->rbuf, b->maxreal);
    }
    /* Copy pointer */
    rbuf = b->rbuf + b->nreal;
    for (i = 0; (i < nr); i++)
    {
        rbuf[i] = r[i];
    }

    index = b->nreal;
    b->nreal += nr;

    return index;
}

void sum_bin(t_bin* b, const t_commrec* cr)
{
    int i;

    for (i = b->nreal; (i < b->maxreal); i++)
    {
        b->rbuf[i] = 0;
    }
    gmx_sumd(b->maxreal, b->rbuf, cr);
}

void extract_binr(t_bin* b, int index, int nr, real r[])
{
    int     i;
    double* rbuf;

    rbuf = b->rbuf + index;
    for (i = 0; (i < nr); i++)
    {
        r[i] = rbuf[i];
    }
}

void extract_binr(t_bin* b, int index, gmx::ArrayRef<real> r)
{
    extract_binr(b, index, r.size(), r.data());
}

void extract_bind(t_bin* b, int index, int nr, double r[])
{
    int     i;
    double* rbuf;

    rbuf = b->rbuf + index;
    for (i = 0; (i < nr); i++)
    {
        r[i] = rbuf[i];
    }
}
