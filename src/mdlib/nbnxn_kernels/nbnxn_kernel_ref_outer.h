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


#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    NBNXN_CPU_CLUSTER_I_SIZE


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
NBK_FUNC_NAME(nbnxn_kernel_ref,noener)
#else
#ifndef ENERGY_GROUPS
NBK_FUNC_NAME(nbnxn_kernel_ref,ener)
#else
NBK_FUNC_NAME(nbnxn_kernel_ref,energrp)
#endif
#endif
#undef NBK_FUNC_NAME
                            (const nbnxn_pairlist_t     *nbl,
                             const nbnxn_atomdata_t     *nbat,
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
    const nbnxn_ci_t   *nbln;
    const nbnxn_cj_t   *l_cj;
    const int          *type;
    const real         *q;
    const real         *shiftvec;
    const real         *x;
    const real         *nbfp;
    real       rcut2;
    int        ntype2;
    real       facel;
    real       *nbfp_i;
    int        n,ci,ci_sh;
    int        ish,ish3;
    gmx_bool   half_LJ,do_coul;
    int        cjind0,cjind1,cjind;
    int        ip,jp;

    real       xi[UNROLLI*DIM];
    real       fi[UNROLLI*DIM];
    real       qi[UNROLLI];

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS
    
    real       Vvdw_ci,Vc_ci;
#else
    int        egp_mask;
    int        egp_sh_i[UNROLLI];
#endif
    real       sh_invrc6;
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
#ifndef GMX_DOUBLE
    const real *tab_coul_FDV0;
#else
    const real *tab_coul_F;
    const real *tab_coul_V;
#endif
#endif

    int ninner;

#ifdef COUNT_PAIRS
    int npair=0;
#endif

#ifdef CALC_ENERGIES
    sh_invrc6 = ic->sh_invrc6;
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

#ifndef GMX_DOUBLE
    tab_coul_FDV0 = ic->tabq_coul_FDV0;
#else
    tab_coul_F    = ic->tabq_coul_F;
    tab_coul_V    = ic->tabq_coul_V;
#endif
#endif

#ifdef ENERGY_GROUPS
    egp_mask = (1<<nbat->neg_2log) - 1;
#endif


    rcut2               = ic->rvdw*ic->rvdw;

    ntype2              = nbat->ntype*2;
    nbfp                = nbat->nbfp;
    q                   = nbat->q;
    type                = nbat->type;
    facel               = ic->epsfac;
    shiftvec            = shift_vec[0];
    x                   = nbat->x;

    l_cj = nbl->cj;

    ninner = 0;
    for(n=0; n<nbl->nci; n++)
    {
        int i,d;

        nbln = &nbl->ci[n];

        ish              = (nbln->shift & NBNXN_CI_SHIFT);
        ish3             = ish*3;
        cjind0           = nbln->cj_ind_start;      
        cjind1           = nbln->cj_ind_end;    
        /* Currently only works super-cells equal to sub-cells */
        ci               = nbln->ci;
        ci_sh            = (ish == CENTRAL ? ci : -1);

        half_LJ = (nbln->shift & NBNXN_CI_HALF_LJ(0));
        do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS
        Vvdw_ci = 0;
        Vc_ci   = 0;
#else
        for(i=0; i<UNROLLI; i++)
        {
            egp_sh_i[i] = ((nbat->energrp[ci]>>(i*nbat->neg_2log)) & egp_mask)*nbat->nenergrp;
        }
#endif
#endif

        for(i=0; i<UNROLLI; i++)
        {
            for(d=0; d<DIM; d++)
            {
                xi[i*DIM+d] = x[(ci*UNROLLI+i)*DIM+d] + shiftvec[ish3+d];
                fi[i*DIM+d] = 0;
            }
        }

        /* With half_LJ we currently always calculate Coulomb interactions */
        if (do_coul || half_LJ)
        {
            for(i=0; i<UNROLLI; i++)
            {
                qi[i] = facel*q[ci*UNROLLI+i];

#ifdef CALC_ENERGIES
                if (l_cj[nbln->cj_ind_start].cj == ci_sh)
                {
#ifdef ENERGY_GROUPS
                    Vc[egp_sh_i[i]+((nbat->energrp[ci]>>(i*nbat->neg_2log)) & egp_mask)]
#else
                    Vc[0]
#endif
#ifdef CALC_COUL_RF
                        -= qi[i]*q[ci*UNROLLI+i]*0.5*c_rf;
#else
                        -= qi[i]*q[ci*UNROLLI+i]*ic->ewaldcoeff*0.564189583548;
#endif
                }
#endif
            }
        }

        cjind = cjind0;
        while (cjind < cjind1 && nbl->cj[cjind].excl != 0xffff)
        {
#define CHECK_EXCLS
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
            }
            else if (do_coul)
            {
#include "nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "nbnxn_kernel_ref_inner.h"
            }
#undef CHECK_EXCLS
            cjind++;
        }

        for(; (cjind<cjind1); cjind++)
        {
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
            }
            else if (do_coul)
            {
#include "nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "nbnxn_kernel_ref_inner.h"
            }
        }
        ninner += cjind1 - cjind0;

        /* Add i forces to mem force list */
        for(i=0; i<UNROLLI; i++)
        {
            for(d=0; d<DIM; d++)
            {
                f[(ci*UNROLLI+i)*DIM+d] += fi[i*DIM+d];
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
