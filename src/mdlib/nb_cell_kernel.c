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


/* Include the force+energy kernels */
#define CALC_SHIFTFORCES
#define CALC_ENERGIES
#define LJ_COMB_GEOM
#include "nb_cell_kernel_sse2_single.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nb_cell_kernel_sse2_single.h"
#undef LJ_COMB_LB
#include "nb_cell_kernel_sse2_single.h"
#undef CALC_ENERGIES
#undef CALC_SHIFTFORCES

/* Include the force only kernels */
#define LJ_COMB_GEOM
#include "nb_cell_kernel_sse2_single.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nb_cell_kernel_sse2_single.h"
#undef LJ_COMB_LB
#include "nb_cell_kernel_sse2_single.h"

typedef void (*p_nbk_func_ener)(const gmx_nblist_t         *nbl,
                                const gmx_nb_atomdata_t    *nbat,
                                const interaction_const_t  *ic,
                                real                       tabscale,
                                const real                 *VFtab,
                                rvec                       *shift_vec,
                                real                       *f,
                                real                       *fshift,
                                real                       *Vc,
                                real                       *Vvdw);

typedef void (*p_nbk_func_noener)(const gmx_nblist_t         *nbl,
                                  const gmx_nb_atomdata_t    *nbat,
                                  const interaction_const_t  *ic,
                                  real                       tabscale,  
                                  const real                 *VFtab,
                                  rvec                       *shift_vec,
                                  real                       *f);

p_nbk_func_ener p_nbk_ener[ljcrNR] =
{ nb_cell_kernel_sse2_single_comb_geom_ener,
  nb_cell_kernel_sse2_single_comb_lb_ener,
  nb_cell_kernel_sse2_single_comb_none_ener };

p_nbk_func_noener p_nbk_noener[ljcrNR] =
{ nb_cell_kernel_sse2_single_comb_geom_noener,
  nb_cell_kernel_sse2_single_comb_lb_noener,
  nb_cell_kernel_sse2_single_comb_none_noener };

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
               const interaction_const_t  *ic,
               real                       tabscale,  
               const real                 *VFtab,
               rvec                       *shift_vec, 
               int                        force_flags,
               gmx_bool                   clearF,
               real                       *fshift,
               real                       *Vc,
               real                       *Vvdw)
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

            p_nbk_ener[nbat->comb_rule](nbl[i],nbat,
                                        ic,tabscale,VFtab,
                                        shift_vec,
                                        nbat->out[i].f,
                                        nnbl == 1 ?
                                        fshift :
                                        nbat->out[i].fshift,
                                        &nbat->out[i].Vc,
                                        &nbat->out[i].Vvdw);
        }
        else
        {
            p_nbk_noener[nbat->comb_rule](nbl[i],nbat,
                                          ic,tabscale,VFtab,
                                          shift_vec,
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
