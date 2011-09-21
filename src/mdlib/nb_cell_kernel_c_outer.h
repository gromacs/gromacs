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


#define UNROLLI    4
#define UNROLLJ    4


/* All functionality defines are set here, except for:
 * CALC_ENERGIES, ENERGY_GROUPS which are defined before.
 * CHECK_EXCLS, which is set just before including the inner loop contents.
 */

/* We always calculate shift forces, because it's cheap anyhow */
#define CALC_SHIFTFORCES

#ifdef CALC_COUL_RF
#define NBK_FUNC_NAME(x,y) x##_rf_##y
#else
#define NBK_FUNC_NAME(x,y) x##_tab_##y
#endif

static void
#ifndef CALC_ENERGIES
NBK_FUNC_NAME(nb_cell_kernel_c,noener)
#else
#ifndef ENERGY_GROUPS
NBK_FUNC_NAME(nb_cell_kernel_c,ener)
#else
NBK_FUNC_NAME(nb_cell_kernel_c,energrp)
#endif
#endif
#undef NBK_FUNC_NAME
                            (const gmx_nblist_t         *nbl,
                             const gmx_nb_atomdata_t    *nbat,
                             const interaction_const_t  *ic,
                             rvec                       *shift_vec, 
                             real                       *f
#ifdef CALC_SHIFTFORCES
                             ,
                             real                       *fshift
#endif
#ifdef CALC_ENERGIES
                             ,
                             real                       *Vvdw,
                             real                       *Vc
#endif
                            )
{
    const gmx_nbl_ci_t *nbln;
    const gmx_nbl_cj_t *cj;
    const int          *type;
    const real         *q;
    const real         *shiftvec;
    const real         *x;
    const real         *nbfp;
    real       rcut2;
    int        ntype2;
    real       facel;
    real       *nbfp_i;
    int        n,si;
    int        ish3;
    gmx_bool   half_LJ,do_coul;
    int        ssi;
    int        sjind0,sjind1,sjind;
    int        ip,jp;

    real       xi[UNROLLI*DIM];
    real       fi[UNROLLI*DIM];
    real       qi[UNROLLI];

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS
    real       Vvdw_ci,Vc_ci;
#else
    int        egp_sh_i[UNROLLI];
#endif
#endif
    
#ifdef CALC_COUL_RF
    real       k_rf2;
#ifdef CALC_ENERGIES
    real       k_rf,c_rf;
#endif
#else
    real       tabscale;
#ifdef CALC_ENERGIES
    real       halfsp;
#endif
    const real *tab_coul_FDV0;
#endif

    int ninner;

#ifdef COUNT_PAIRS
    int npair=0;
#endif

#ifdef CALC_COUL_RF
    k_rf2 = 2*ic->k_rf;
#ifdef CALC_ENERGIES
    k_rf = ic->k_rf;
    c_rf = ic->c_rf;
#endif
#else
    tabscale = ic->tabq_scale;
#ifdef CALC_ENERGIES
    halfsp = 0.5/ic->tabq_scale;
#endif

    tab_coul_FDV0 = ic->tabq_coul_FDV0;
#endif

    rcut2               = ic->rvdw*ic->rvdw;

    ntype2              = nbat->ntype*2;
    nbfp                = nbat->nbfp;
    q                   = nbat->q;
    type                = nbat->type;
    facel               = ic->epsfac;
    shiftvec            = shift_vec[0];
    x                   = nbat->x;

    cj = nbl->cj;

    ninner = 0;
    for(n=0; n<nbl->nci; n++)
    {
        int i,d;

        nbln = &nbl->ci[n];

        ish3             = 3*(nbln->shift & NBL_CI_SHIFT);
        sjind0           = nbln->cj_ind_start;      
        sjind1           = nbln->cj_ind_end;    
        /* Currently only works super-cells equal to sub-cells */
        si               = nbln->ci;

        half_LJ = (nbln->shift & NBL_CI_HALF_LJ(0));
        do_coul = (nbln->shift & NBL_CI_DO_COUL(0));

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS
        Vvdw_ci = 0;
        Vc_ci   = 0;
#else
        for(i=0; i<UNROLLI; i++)
        {
            egp_sh_i[i] = ((nbat->energrp[si]>>(8*i)) & 255)*nbat->nenergrp;
        }
#endif
#endif

        for(i=0; i<UNROLLI; i++)
        {
            for(d=0; d<DIM; d++)
            {
                xi[i*DIM+d] = x[(si*UNROLLI+i)*DIM+d] + shiftvec[ish3+d];
                fi[i*DIM+d] = 0;
            }
        }

        /* With half_LJ we currently always calculate Coulomb interactions */
        if (do_coul || half_LJ)
        {
            for(i=0; i<UNROLLI; i++)
            {
                qi[i] = facel*q[si*UNROLLI+i];
            }
        }

        sjind = sjind0;
        while (sjind < sjind1 && nbl->cj[sjind].excl != 0xffff)
        {
#define CHECK_EXCLS
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "nb_cell_kernel_c_inner.h"
#undef HALF_LJ
            }
            else if (do_coul)
            {
#include "nb_cell_kernel_c_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "nb_cell_kernel_c_inner.h"
            }
#undef CHECK_EXCLS
            sjind++;
        }

        for(; (sjind<sjind1); sjind++)
        {
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "nb_cell_kernel_c_inner.h"
#undef HALF_LJ
            }
            else if (do_coul)
            {
#include "nb_cell_kernel_c_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "nb_cell_kernel_c_inner.h"
            }
        }
        ninner += sjind1 - sjind0;

        /* Add i forces to mem force list */
        for(i=0; i<UNROLLI; i++)
        {
            for(d=0; d<DIM; d++)
            {
                f[(si*UNROLLI+i)*DIM+d] += fi[i*DIM+d];
            }
        }
#ifdef CALC_SHIFTFORCES
        if (fshift != NULL)
        {
            /* Add i forces to shifted force list */
            for(i=0; i<UNROLLI; i++)
            {
                for(d=0; d<DIM; d++)
                {
                    fshift[ish3+d] += fi[i*DIM+d];
                }
            }
        }
#endif

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS
        *Vvdw += Vvdw_ci;
        *Vc   += Vc_ci;
#endif
#endif
	}

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n",npair);
#endif
}

#undef CALC_SHIFTFORCES

#undef UNROLLI   
#undef UNROLLJ   
