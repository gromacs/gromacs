/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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

#include "energyhistory.h"

#include <cstring>

#include <algorithm>

#include "gromacs/utility/smalloc.h"

static void done_delta_h_history(delta_h_history_t *dht)
{
    int i;

    for (i = 0; i < dht->nndh; i++)
    {
        sfree(dht->dh[i]);
    }
    sfree(dht->dh);
    sfree(dht->ndh);
}

void init_energyhistory(energyhistory_t * enerhist)
{
    enerhist->nener = 0;

    enerhist->ener_ave     = NULL;
    enerhist->ener_sum     = NULL;
    enerhist->ener_sum_sim = NULL;
    enerhist->dht          = NULL;

    enerhist->nsteps     = 0;
    enerhist->nsum       = 0;
    enerhist->nsteps_sim = 0;
    enerhist->nsum_sim   = 0;

    enerhist->dht = NULL;
}

void done_energyhistory(energyhistory_t * enerhist)
{
    sfree(enerhist->ener_ave);
    sfree(enerhist->ener_sum);
    sfree(enerhist->ener_sum_sim);

    if (enerhist->dht != NULL)
    {
        done_delta_h_history(enerhist->dht);
        sfree(enerhist->dht);
    }
}
