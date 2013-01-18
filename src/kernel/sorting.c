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

#include <limits.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "sorting.h"

/*****************************************************************************
 *                                                                           *
 *                   Block sorting on coordinates                            *
 *                                                                           *
 *****************************************************************************/

static rvec *make_xblock(t_block *block, rvec x[])
{
    int   i, j, k, nr, n;
    rvec *xblock;

    nr = block->nr;
    snew(xblock, nr);
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            xblock[i][j] = 0.0;
        }
        for (j = block->index[i]; j < (int)(block->index[i+1]); j++)
        {
            for (k = 0; k < DIM; k++)
            {
                xblock[i][k] += x[j][k];
            }
        }
        n = block->index[i+1]-block->index[i];
        for (k = 0; k < DIM; k++)
        {
            xblock[i][k] /= n;
        }
    }
    return xblock;
}

static rvec *xblock; /* just global to bcomp1, used in qsort */

static int bomp1(const void *p1, const void *p2)
{
    int i, i1, i2;

    i1 = *(int *)p1;
    i2 = *(int *)p2;
    for (i = 0; i < DIM; i++)
    {
        if (xblock[i1][i] < xblock[i2][i])
        {
            return -1;
        }
        else if (xblock[i1][i] > xblock[i2][i])
        {
            return 1;
        }
    }
    return 0;
}

void sort_xblock(t_block *block, rvec x[], int renum[])
{
    int i, nr, *invnum;

    nr = block->nr;
    snew(invnum, nr);
    xblock = make_xblock(block, x);
    for (i = 0; i < nr; i++)
    {
        invnum[i] = i;
    }
    qsort((void *)invnum, nr, (size_t)sizeof(invnum[0]), bomp1);
    for (i = 0; i < nr; i++)
    {
        renum[invnum[i]] = i;
    }
    sfree(xblock);
    sfree(invnum);
}

/*****************************************************************************
 *                                                                           *
 *                        Bond list sorting                                  *
 *                                                                           *
 *****************************************************************************/

static int bcomp2(const void *p1, const void *p2)
{
    int done;

    if ((((atom_id *)p1)[0]) != (((atom_id *)p2)[0]))
    {
        done = ((((atom_id *)p1)[0])-(((atom_id *)p2)[0]));
    }
    else
    {
        done = ((((atom_id *)p1)[1])-(((atom_id *)p2)[1]));
    }
#ifdef DEBUG
    printf("bcomp2: [%d,%d] with [%d,%d] result %d\n",
           ((atom_id *)p1)[0], ((atom_id *)p1)[1],
           ((atom_id *)p2)[0], ((atom_id *)p2)[1], done);
#endif
    return done;
}

void sort_bond_list(t_bond bonds[], int nr)
{
    qsort((void *)bonds, nr, (size_t)sizeof(bonds[0]), bcomp2);
}
