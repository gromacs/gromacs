/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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

#include "invblock.h"

#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

atom_id *make_invblock(const t_block *block, int nr)
{
    int      i, j;
    atom_id *invblock;

    snew(invblock, nr+1);
    /* Mark unused numbers */
    for (i = 0; i <= nr; i++)
    {
        invblock[i] = NO_ATID;
    }
    for (i = 0; (i < block->nr); i++)
    {
        for (j = block->index[i]; (j < block->index[i+1]); j++)
        {
            if (invblock[j] == NO_ATID)
            {
                invblock[j] = i;
            }
            else
            {
                gmx_fatal(FARGS, "Double entries in block structure. Item %d is in blocks %d and %d\n"
                          " Cannot make an unambiguous inverse block.",
                          j, i, invblock[j]);
            }
        }
    }
    return invblock;
}

atom_id *make_invblocka(const t_blocka *block, int nr)
{
    int      i, j;
    atom_id *invblock;

    snew(invblock, nr+1);
    /* Mark unused numbers */
    for (i = 0; i <= nr; i++)
    {
        invblock[i] = NO_ATID;
    }
    for (i = 0; (i < block->nr); i++)
    {
        for (j = block->index[i]; (j < block->index[i+1]); j++)
        {
            if (invblock[block->a[j]] == NO_ATID)
            {
                invblock[block->a[j]] = i;
            }
            else
            {
                gmx_fatal(FARGS, "Double entries in block structure. Item %d is in blocks %d and %d\n"
                          " Cannot make an unambiguous inverse block.",
                          j, i, invblock[block->a[j]]);
            }
        }
    }
    return invblock;
}
