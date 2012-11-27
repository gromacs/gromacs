/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 * Copyright (c) 2012, by the GROMACS development team, led by
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

#include "nbnxn_kernel_common.h"

void
clear_f(const nbnxn_atomdata_t *nbat,real *f)
{
    int i;

    for(i=0; i<nbat->natoms*nbat->fstride; i++)
    {
        f[i] = 0;
    }
}

void
clear_fshift(real *fshift)
{
    int i;

    for(i=0; i<SHIFTS*DIM; i++)
    {
        fshift[i] = 0;
    }
}

void
reduce_energies_over_lists(const nbnxn_atomdata_t     *nbat,
                           int                        nlist,
                           real                       *Vvdw,
                           real                       *Vc)
{
    int nb;
    int i,j,ind,indr;

    for(nb=0; nb<nlist; nb++)
    {
        for(i=0; i<nbat->nenergrp; i++)
        {
            /* Reduce the diagonal terms */
            ind = i*nbat->nenergrp + i;
            Vvdw[ind] += nbat->out[nb].Vvdw[ind];
            Vc[ind]   += nbat->out[nb].Vc[ind];

            /* Reduce the off-diagonal terms */
            for(j=i+1; j<nbat->nenergrp; j++)
            {
                /* The output should contain only one off-diagonal part */
                ind  = i*nbat->nenergrp + j;
                indr = j*nbat->nenergrp + i;
                Vvdw[ind] += nbat->out[nb].Vvdw[ind] + nbat->out[nb].Vvdw[indr];
                Vc[ind]   += nbat->out[nb].Vc[ind]   + nbat->out[nb].Vc[indr];
            }
        }
    }
}
