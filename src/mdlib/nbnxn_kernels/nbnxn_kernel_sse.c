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
#include "gmx_omp_nthreads.h"

#if ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) )

#include "nbnxn_kernel_sse.h"

/* Include all flavors of the SSE kernel loops */

/* Analytical reaction-field kernels */
#define CALC_COUL_RF

/* Include the force+energy kernels */
#define CALC_ENERGIES
#define LJ_COMB_GEOM
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef CALC_ENERGIES

/* Include the force+energygroups kernels */
#define CALC_ENERGIES
#define ENERGY_GROUPS
#define LJ_COMB_GEOM
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef ENERGY_GROUPS
#undef CALC_ENERGIES

/* Include the force only kernels */
#define LJ_COMB_GEOM
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"

#undef CALC_COUL_RF


/* Tabulated Coulomb kernels */

/* Include the force+energy kernels */
#define CALC_ENERGIES
#define LJ_COMB_GEOM
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef CALC_ENERGIES

/* Include the force+energygroups kernels */
#define CALC_ENERGIES
#define ENERGY_GROUPS
#define LJ_COMB_GEOM
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef ENERGY_GROUPS
#undef CALC_ENERGIES

/* Include the force only kernels */
#define LJ_COMB_GEOM
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_sse_outer.h"


typedef void (*p_nbk_func_ener)(const nbnxn_pairlist_t     *nbl,
                                const nbnxn_atomdata_t     *nbat,
                                const interaction_const_t  *ic,
                                rvec                       *shift_vec,
                                real                       *f,
                                real                       *fshift,
                                real                       *Vvdw,
                                real                       *Vc);

typedef void (*p_nbk_func_noener)(const nbnxn_pairlist_t     *nbl,
                                  const nbnxn_atomdata_t     *nbat,
                                  const interaction_const_t  *ic,
                                  rvec                       *shift_vec,
                                  real                       *f,
                                  real                       *fshift);

enum { coultRF, coultTAB, coultNR };


p_nbk_func_ener p_nbk_ener[coultNR][ljcrNR] =
{ { nbnxn_kernel_sse_rf_comb_geom_ener,
    nbnxn_kernel_sse_rf_comb_lb_ener,
    nbnxn_kernel_sse_rf_comb_none_ener },
  { nbnxn_kernel_sse_tab_comb_geom_ener,
    nbnxn_kernel_sse_tab_comb_lb_ener,
    nbnxn_kernel_sse_tab_comb_none_ener } };

p_nbk_func_ener p_nbk_energrp[coultNR][ljcrNR] =
{ { nbnxn_kernel_sse_rf_comb_geom_energrp,
    nbnxn_kernel_sse_rf_comb_lb_energrp,
    nbnxn_kernel_sse_rf_comb_none_energrp },
  { nbnxn_kernel_sse_tab_comb_geom_energrp,
    nbnxn_kernel_sse_tab_comb_lb_energrp,
    nbnxn_kernel_sse_tab_comb_none_energrp } };

p_nbk_func_noener p_nbk_noener[coultNR][ljcrNR] =
{ { nbnxn_kernel_sse_rf_comb_geom_noener,
    nbnxn_kernel_sse_rf_comb_lb_noener,
    nbnxn_kernel_sse_rf_comb_none_noener },
  { nbnxn_kernel_sse_tab_comb_geom_noener,
    nbnxn_kernel_sse_tab_comb_lb_noener,
    nbnxn_kernel_sse_tab_comb_none_noener } };

#endif /* SSE */


static void clear_f(const nbnxn_atomdata_t *nbat,
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

static void reduce_group_energies(int ng,int ng_2log,
                                  const real *VSvdw,const real *VSc,
                                  real *Vvdw,real *Vc)
{
    int ng_p2,i,j,j0,j1,c,s;

#ifndef GMX_DOUBLE
#define SIMD_WIDTH   4
#define SIMD_WIDTH_2 2
#else
#define SIMD_WIDTH   2
#define SIMD_WIDTH_2 1
#endif

    ng_p2 = (1<<ng_2log);

    /* The size of the SSE energy group buffer array is:
     * stride^3*SIMD_WIDTH_2*SIMD_WIDTH
     */
    for(i=0; i<ng; i++)
    {
        for(j=0; j<ng; j++)
        {
            Vvdw[i*ng+j] = 0;
            Vc[i*ng+j]   = 0;
        }

        for(j1=0; j1<ng; j1++)
        {
            for(j0=0; j0<ng; j0++)
            {
                c = ((i*ng + j1)*ng_p2 + j0)*SIMD_WIDTH_2*SIMD_WIDTH;
                for(s=0; s<SIMD_WIDTH_2; s++)
                {
                    Vvdw[i*ng+j0] += VSvdw[c+0];
                    Vvdw[i*ng+j1] += VSvdw[c+1];
                    Vc  [i*ng+j0] += VSc  [c+0];
                    Vc  [i*ng+j1] += VSc  [c+1];
                    c += SIMD_WIDTH + 2;
                }
            }
        }
    }
}

void
nbnxn_kernel_sse(nbnxn_pairlist_set_t       *nbl_list,
                 const nbnxn_atomdata_t     *nbat,
                 const interaction_const_t  *ic,
                 rvec                       *shift_vec, 
                 int                        force_flags,
                 gmx_bool                   clearF,
                 real                       *fshift,
                 real                       *Vc,
                 real                       *Vvdw)
#if ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) )
{
    int              nnbl;
    nbnxn_pairlist_t **nbl;
    int coult;
    int nb;

    nnbl = nbl_list->nnbl;
    nbl  = nbl_list->nbl;

    if (EEL_RF(ic->eeltype) || ic->eeltype == eelCUT)
    {
        coult = coultRF;
    }
    else
    {
        coult = coultTAB;
    }

#pragma omp parallel for schedule(static) num_threads(gmx_omp_nthreads_get(emntNonbonded))
    for(nb=0; nb<nnbl; nb++)
    {
        nbnxn_atomdata_output_t *out;
        real *fshift_p;

        out = &nbat->out[nb];

        if (clearF)
        {
            clear_f(nbat,out->f);
        }

        if ((force_flags & GMX_FORCE_VIRIAL) && nnbl == 1)
        {
            fshift_p = fshift;
        }
        else
        {
            fshift_p = out->fshift;

            if (clearF)
            {
                clear_fshift(fshift_p);
            }
        }

        /* With Ewald type electrostatics we the forces for excluded atom pairs
         * should not contribute to the virial sum. The exclusion forces
         * are not calculate in the energy kernels, but are in _noener.
         */
        if (!((force_flags & GMX_FORCE_ENERGY) ||
              (EEL_FULL(ic->eeltype) && (force_flags & GMX_FORCE_VIRIAL))))
        {
            /* Don't calculate energies */
            p_nbk_noener[coult][nbat->comb_rule](nbl[nb],nbat,
                                                 ic,
                                                 shift_vec,
                                                 out->f,
                                                 fshift_p);
        }
        else if (out->nV == 1 || !(force_flags & GMX_FORCE_ENERGY))
        {
            /* No energy groups */
            out->Vvdw[0] = 0;
            out->Vc[0]   = 0;

            p_nbk_ener[coult][nbat->comb_rule](nbl[nb],nbat,
                                               ic,
                                               shift_vec,
                                               out->f,
                                               fshift_p,
                                               out->Vvdw,
                                               out->Vc);
        }
        else
        {
            /* Calculate energy group contributions */
            int i;

            for(i=0; i<out->nVS; i++)
            {
                out->VSvdw[i] = 0;
            }
            for(i=0; i<out->nVS; i++)
            {
                out->VSc[i] = 0;
            }

            p_nbk_energrp[coult][nbat->comb_rule](nbl[nb],nbat,
                                                  ic,
                                                  shift_vec,
                                                  out->f,
                                                  fshift_p,
                                                  out->VSvdw,
                                                  out->VSc);

            reduce_group_energies(nbat->nenergrp,nbat->neg_2log,
                                  out->VSvdw,out->VSc,
                                  out->Vvdw,out->Vc);
        }
    }

    if (force_flags & GMX_FORCE_ENERGY)
    {
        /* Reduce the energies */
        for(nb=0; nb<nnbl; nb++)
        {
            int i;

            for(i=0; i<nbat->out[nb].nV; i++)
            {
                Vvdw[i] += nbat->out[nb].Vvdw[i];
                Vc[i]   += nbat->out[nb].Vc[i];
            }
        }
        /* printf("vdw %f c %f\n",Vvdw[0],Vc[0]); */
    }
}
#else
{
    gmx_incons("nbnxn_kernel_sse called while GROMACS was configured with SSE enabled");
}
#endif
