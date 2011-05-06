/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "typedefs.h"
#include "vec.h"
#include "smalloc.h"
#include "force.h"
#include "nb_cell_kernel.h"


/* Include the force+energy kernel */
#define CALC_SHIFTFORCES
#define CALC_ENERGIES
#include "nb_cell_kernel_sse2_single.h"
#undef CALC_ENERGIES
#undef CALC_SHIFTFORCES

/* Include the force only kernel */
#include "nb_cell_kernel_sse2_single.h"


static void clear_f(const gmx_nb_atomdata_t *nbat,
                    real *f)
{
    int i;

    for(i=0; i<nbat->natoms*nbat->xstride; i++)
    {
        f[i] = 0;
    }
}

static void clear_fshift(real *fshift)
{
    int i;

    for(i=0; i<SHIFTS*DIM; i++)
    {
        fshift[i] = 0;
    }
}

void
nb_cell_kernel(int                        nnbl,
               gmx_nblist_t               **nbl,
               const gmx_nb_atomdata_t    *nbat,
               const t_forcerec *         fr,
               real                       tabscale,  
               const real *               VFtab,
               int                        force_flags,
               gmx_bool                   clearF,
               real *                     Vc,
               real *                     Vvdw)
{
    int i;

#pragma omp parallel for schedule(static)
    for(i=0; i<nnbl; i++)
    {
        if (clearF)
        {
            clear_f(nbat,nbat->out[i].f);

            if (nnbl > 1)
            {
                clear_fshift(nbat->out[i].fshift);
            }
        }

        if (force_flags & GMX_FORCE_ENERGY)
        {
            nbat->out[i].Vc   = 0;
            nbat->out[i].Vvdw = 0;
            nb_cell_kernel_sse2_single_ener(nbl[i],nbat,
                                            fr,tabscale,VFtab,
                                            nbat->out[i].f,
                                            nnbl == 1 ?
                                            fr->fshift[0] :
                                            nbat->out[i].fshift,
                                            &nbat->out[i].Vc,
                                            &nbat->out[i].Vvdw);
        }
        else
        {
            nb_cell_kernel_sse2_single_noener(nbl[i],nbat,
                                              fr,tabscale,VFtab,
                                              nbat->out[i].f);
        }
    }

    if (force_flags & GMX_FORCE_ENERGY)
    {
        /* Reduce the energies */
        for(i=0; i<nnbl; i++)
        {
            *Vc   += nbat->out[i].Vc;
            *Vvdw += nbat->out[i].Vvdw;
        }
    }
}
