/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.6.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2011, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "network.h"
#include "calcgrid.h"
#include "pme.h"
#include "vec.h"
#include "domdec.h"
#include "nbnxn_cuda_data_mgmt.h"
#include "force.h"
#include "pme_loadbal.h"

typedef struct {
    real rcut;
    real rlist;
    real spacing;
    ivec grid;
    real grid_eff;
    real coeff;
    gmx_pme_t pmedata;

    int  count;
    double cycles;
} pme_setup_t;

/* In the initial scan, step by grids that are at least a factor 0.8 coarser */
#define PMES_GRID_SCALE_FAC  0.8
/* In the initial scan, try to skip grids with uneven x/y/z spacing,
 * checking if the "efficiency" is more than 5% worse than the previous grid.
 */
#define PMES_GRID_EFF_FAC  1.05
/* Rerun up till 12% slower setups than the fastest up till now */
#define PMES_SLOW_FAC  1.12
/* If setups get more than 2% faster, do another round to avoid
 * choosing a slower setup due to acceleration or fluctuations.
 */
#define PMES_ACCEL_TOL 1.02

typedef struct pme_loadbal {
    int  nstage;        /* the current maximum number of stages */

    real cut_spacing;   /* the minimum cutoff / PME grid spacing ratio */
    real rbuf;          /* the pairlist buffer size */
    matrix box_start;   /* the initial simulation box */
    int n;              /* the count of setup as well as the allocation size */
    pme_setup_t *setup; /* the PME+cutoff setups */
    int cur;            /* the current setup */
    int fastest;        /* fastest setup up till now */
    int start;          /* start of setup range to consider in stage>0 */
    int end;            /* end   of setup range to consider in stage>0 */

    int stage;          /* the current stage */
} t_pme_loadbal;

void pme_loadbal_init(pme_loadbal_t *pmes_p,
                      const t_inputrec *ir,matrix box,
                      const interaction_const_t *ic,
                      gmx_pme_t pmedata)
{
    pme_loadbal_t pmes;
    real spm,sp;
    int  d;

    snew(pmes,1);

    /* Any number of stages >= 2 is supported */
    pmes->nstage   = 2;

    pmes->rbuf = ic->rlist - ic->rcoulomb;

    copy_mat(box,pmes->box_start);
    if (ir->ePBC==epbcXY && ir->nwall==2)
    {
        svmul(ir->wall_ewald_zfac,pmes->box_start[ZZ],pmes->box_start[ZZ]);
    }

    pmes->n = 1;
    snew(pmes->setup,pmes->n);

    pmes->cur = 0;
    pmes->setup[0].rcut     = ic->rcoulomb;
    pmes->setup[0].rlist    = ic->rlist;
    pmes->setup[0].grid[XX] = ir->nkx;
    pmes->setup[0].grid[YY] = ir->nky;
    pmes->setup[0].grid[ZZ] = ir->nkz;
    pmes->setup[0].coeff    = ic->ewaldcoeff;

    pmes->setup[0].pmedata  = pmedata;
    
    spm = 0;
    for(d=0; d<DIM; d++)
    {
        sp = norm(pmes->box_start[d])/pmes->setup[0].grid[d];
        if (sp > spm)
        {
            spm = sp;
        }
    }
    pmes->setup[0].spacing = spm;

    if (ir->fourier_spacing > 0)
    {
        pmes->cut_spacing = ir->rcoulomb/ir->fourier_spacing;
    }
    else
    {
        pmes->cut_spacing = ir->rcoulomb/pmes->setup[0].spacing;
    }

    pmes->stage = 0;

    pmes->fastest = 0;
    pmes->start   = 0;

    *pmes_p = pmes;
}

static gmx_bool pme_loadbal_increase_cutoff(pme_loadbal_t pmes,int pme_order)
{
    pme_setup_t *set;
    real fac,sp;
    int d;

    /* Try to add a new setup with next larger cut-off to the list */
    pmes->n++;
    srenew(pmes->setup,pmes->n);
    set = &pmes->setup[pmes->n-1];
    set->pmedata = NULL;

    fac = 1;
    do
    {
        fac *= 1.01;
        clear_ivec(set->grid);
        sp = calc_grid(NULL,pmes->box_start,
                       fac*pmes->setup[pmes->cur].spacing,
                       &set->grid[XX],
                       &set->grid[YY],
                       &set->grid[ZZ]);

        /* In parallel we can't have grids smaller than 2*pme_order,
         * and we would anyhow not gain much speed at these grid sizes.
         */
        for(d=0; d<DIM; d++)
        {
            if (set->grid[d] <= 2*pme_order)
            {
                pmes->n--;

                return FALSE;
            }
        }
    }
    while (sp <= 1.001*pmes->setup[pmes->cur].spacing);

    set->rcut    = pmes->cut_spacing*sp;
    set->rlist   = set->rcut + pmes->rbuf;
    set->spacing = sp;
    /* The grid efficiency is the size wrt a grid with uniform x/y/z spacing */
    set->grid_eff = 1;
    for(d=0; d<DIM; d++)
    {
        set->grid_eff *= (set->grid[d]*sp)/norm(pmes->box_start[d]);
    }
    /* The Ewald coefficient is inversly proportional to the cut-off */
    set->coeff   = pmes->setup[0].coeff*pmes->setup[0].rcut/set->rcut;

    set->count   = 0;
    set->cycles  = 0;

    if (debug)
    {
        fprintf(debug,"PME loadbal: grid %d %d %d, cutoff %f\n",
                set->grid[XX],set->grid[YY],set->grid[ZZ],set->rcut);
    }

    return TRUE;
}

static void print_grid(FILE *fp_err,FILE *fp_log,
                       const char *pre,
                       const char *desc,
                       const pme_setup_t *set,
                       double cycles)
{
    char buf[STRLEN],buft[STRLEN];
    
    if (cycles >= 0)
    {
        sprintf(buft,": %.1f M-cycles",cycles*1e-6);
    }
    else
    {
        buft[0] = '\0';
    }
    sprintf(buf,"%-11s%10s pme grid %d %d %d, cutoff %.3f%s",
            pre,
            desc,set->grid[XX],set->grid[YY],set->grid[ZZ],set->rcut,
            buft);
    if (fp_err != NULL)
    {
        fprintf(fp_err,"%s\n",buf);
    }
    if (fp_log != NULL)
    {
        fprintf(fp_log,"%s\n",buf);
    }
}

static void switch_to_stage1(pme_loadbal_t pmes)
{
    pmes->start = 0;
    while (pmes->start+1 < pmes->n &&
           (pmes->setup[pmes->start].count == 0 ||
            pmes->setup[pmes->start].cycles >
            pmes->setup[pmes->fastest].cycles*PMES_SLOW_FAC))
    {
        pmes->start++;
    }
    while (pmes->start > 0 && pmes->setup[pmes->start-1].cycles == 0)
    {
        pmes->start--;
    }

    pmes->end = pmes->n;
    if (pmes->setup[pmes->end-1].count > 0 &&
        pmes->setup[pmes->end-1].cycles >
        pmes->setup[pmes->fastest].cycles*PMES_SLOW_FAC)
    {
        pmes->end--;
    }

    pmes->stage = 1;

    /* Start add start, 1 will be added immediately after returning */
    pmes->cur = pmes->start - 1;
}

gmx_bool pme_loadbalance(pme_loadbal_t pmes,
                         t_commrec *cr,
                         FILE *fp_err,
                         FILE *fp_log,
                         t_inputrec *ir,
                         t_state *state,
                         double cycles,
                         interaction_const_t *ic,
                         nonbonded_verlet_t *nbv,
                         gmx_pme_t *pmedata,
                         int step)
{
    gmx_bool OK;
    pme_setup_t *set;
    double cycles_fast;
    char buf[STRLEN];

    if (pmes->stage == pmes->nstage)
    {
        return FALSE;
    }

    if (PAR(cr))
    {
        gmx_sumd(1,&cycles,cr);
        cycles /= cr->nnodes;
    }

    set = &pmes->setup[pmes->cur];

    set->count++;
    if (set->count % 2 == 1)
    {
        /* Skip the first cycle, because the first step after a switch
         * is much slower due to allocation and/or caching effects.
         */
        return TRUE;
    }

    sprintf(buf, "step %4d: ", step);
    print_grid(fp_err,fp_log,buf,"timed with",set,cycles);

    if (set->count <= 2)
    {
        set->cycles = cycles;
    }
    else
    {
        if (cycles*PMES_ACCEL_TOL < set->cycles &&
            pmes->stage == pmes->nstage - 1)
        {
            /* The performance went up a lot (due to e.g. DD load balancing).
             * Add a stage, keep the minima, but rescan all setups.
             */
            pmes->nstage++;

            if (debug)
            {
                fprintf(debug,"The performance for grid %d %d %d went from %.3f to %.1f M-cycles, this is more than %f\n"
                        "Increased the number stages to %d"
                        " and ignoring the previous performance\n",
                        set->grid[XX],set->grid[YY],set->grid[ZZ],
                        cycles*1e-6,set->cycles*1e-6,PMES_ACCEL_TOL,
                        pmes->nstage);
            }
        }
        set->cycles = min(set->cycles,cycles);
    }

    if (set->cycles < pmes->setup[pmes->fastest].cycles)
    {
        pmes->fastest = pmes->cur;
    }
    cycles_fast = pmes->setup[pmes->fastest].cycles;

    /* Check in stage 0 if we should stop scanning grids.
     * Stop when the time is more than SLOW_FAC longer than the fastest.
     */
    if (pmes->stage == 0 && pmes->cur > 0 &&
        cycles > pmes->setup[pmes->fastest].cycles*PMES_SLOW_FAC)
    {
        pmes->n = pmes->cur + 1;
        /* Done with scanning, go to stage 1 */
        switch_to_stage1(pmes);
    }

    if (pmes->stage == 0)
    {
        int gridsize_start;

        gridsize_start = set->grid[XX]*set->grid[YY]*set->grid[ZZ];

        do
        {
            if (pmes->cur+1 < pmes->n)
            {
                /* We had already generated the next setup */
                OK = TRUE;
            }
            else
            {
                /* Find the next setup */
                OK = pme_loadbal_increase_cutoff(pmes,ir->pme_order);
            }
                
            if (OK && ir->ePBC != epbcNONE)
            {
                OK = (sqr(pmes->setup[pmes->cur+1].rlist)
                      <= max_cutoff2(ir->ePBC,state->box));
            }

            if (OK)
            {
                pmes->cur++;

                if (DOMAINDECOMP(cr))
                {
                    OK = change_dd_cutoff(cr,state,ir,
                                          pmes->setup[pmes->cur].rlist);
                    if (!OK)
                    {
                        /* Failed: do not use this setup */
                        pmes->cur--;
                    }
                }
            }
            if (!OK)
            {
                /* We hit the upper limit for the cut-off,
                 * the setup should not go further than cur.
                 */
                pmes->n = pmes->cur + 1;
                /* Switch to the next stage */
                switch_to_stage1(pmes);
            }
        }
        while (OK &&
               !(pmes->setup[pmes->cur].grid[XX]*
                 pmes->setup[pmes->cur].grid[YY]*
                 pmes->setup[pmes->cur].grid[ZZ] <
                 gridsize_start*PMES_GRID_SCALE_FAC
                 &&
                 pmes->setup[pmes->cur].grid_eff <
                 pmes->setup[pmes->cur-1].grid_eff*PMES_GRID_EFF_FAC));
    }

    if (pmes->stage > 0 && pmes->end == 1)
    {
        pmes->cur = 0;
        pmes->stage = pmes->nstage;
    }
    else if (pmes->stage > 0 && pmes->end > 1)
    {
        /* If stage = nstage-1:
         *   scan over all setups, rerunning only those setups
         *   which are not much slower than the fastest
         * else:
         *   use the next setup
         */
        do
        {
            pmes->cur++;
            if (pmes->cur == pmes->end)
            {
                pmes->stage++;
                pmes->cur = pmes->start;
            }
        }
        while (pmes->stage == pmes->nstage - 1 &&
               pmes->setup[pmes->cur].count > 0 &&
               pmes->setup[pmes->cur].cycles > cycles_fast*PMES_SLOW_FAC);

        if (pmes->stage == pmes->nstage)
        {
            /* We are done optiming, use the fastest setup we found */
            pmes->cur = pmes->fastest;
        }
    }

    if (DOMAINDECOMP(cr) && pmes->stage > 0)
    {
        OK = change_dd_cutoff(cr,state,ir,pmes->setup[pmes->cur].rlist);
        if (!OK)
        {
            /* Failsafe solution */
            if (pmes->cur > 1 && pmes->stage == pmes->nstage)
            {
                pmes->stage--;
            }
            pmes->fastest = 0;
            pmes->start   = 0;
            pmes->end     = pmes->cur;
            pmes->cur     = pmes->start;
        }
    }

    /* Change the Coulomb cut-off and the PME grid */

    set = &pmes->setup[pmes->cur];

    ic->rcoulomb   = set->rcut;
    ic->rlist      = set->rlist;
    ic->ewaldcoeff = set->coeff;

    if (nbv->grp[0].kernel_type == nbk8x8x8_CUDA)
    {
        nbnxn_cuda_pmetune_update_param(nbv->cu_nbv,ic);
    }
    else
    {
        init_interaction_const_tables(NULL,ic,nbv->grp[0].kernel_type);
    }

    if (nbv->ngrp > 1)
    {
        init_interaction_const_tables(NULL,ic,nbv->grp[1].kernel_type);
    }

    if (cr->duty & DUTY_PME)
    {
        if (pmes->setup[pmes->cur].pmedata == NULL)
        {
            /* Generate a new PME data structure,
             * copying part of the old pointers.
             */
            gmx_pme_reinit(&set->pmedata,
                           cr,pmes->setup[0].pmedata,ir,
                           set->grid);
        }
        *pmedata = set->pmedata;
    }
    else
    {
        /* Tell our PME-only node to switch grid */
        gmx_pme_send_switch(cr, set->grid, set->coeff);
    }

    if (debug)
    {
        print_grid(NULL,debug,"","switched to",set,-1);
    }

    if (pmes->stage == pmes->nstage)
    {
        print_grid(fp_err,fp_log,"","optimal",set,-1);
    }

    return TRUE;
}

void restart_pme_loadbal(pme_loadbal_t pmes, int n)
{
    pmes->nstage += n;
}

static int pme_grid_points(const pme_setup_t *setup)
{
    return setup->grid[XX]*setup->grid[YY]*setup->grid[ZZ];
}

static void print_pme_loadbal_setting(FILE *fplog,
                                     char *name,
                                     const pme_setup_t *setup)
{
    fprintf(fplog,
            "   %7s %6.3f nm %6.3f nm     %3d %3d %3d   %5.3f nm  %5.3f nm\n",
            name,
            setup->rcut,setup->rlist,
            setup->grid[XX],setup->grid[YY],setup->grid[ZZ],
            setup->spacing,1/setup->coeff);
}

static void print_pme_loadbal_settings(pme_loadbal_t pmes, FILE *fplog)
{
    double pp_ratio,grid_ratio;

    pp_ratio   = pow(pmes->setup[pmes->cur].rlist/pmes->setup[0].rlist,3.0);
    grid_ratio = pme_grid_points(&pmes->setup[pmes->cur])/
        (double)pme_grid_points(&pmes->setup[0]);

    fprintf(fplog,"\n");
    fprintf(fplog,"       P P   -   P M E   L O A D   B A L A N C I N G\n");
    fprintf(fplog,"\n");
    fprintf(fplog," PP/PME load balancing changed the cut-off and PME settings:\n");
    fprintf(fplog,"           particle-particle                    PME\n");
    fprintf(fplog,"            rcoulomb  rlist            grid      spacing   1/beta\n");
    print_pme_loadbal_setting(fplog,"initial",&pmes->setup[0]);
    print_pme_loadbal_setting(fplog,"optimal",&pmes->setup[pmes->cur]);
    fprintf(fplog," cost-ratio           %4.2f             %4.2f\n",
            pp_ratio,grid_ratio);
    fprintf(fplog," (note that these numbers concern only part of the total PP and PME load)\n");
    fprintf(fplog,"\n");
}

void pme_loadbal_done(pme_loadbal_t pmes, FILE *fplog)
{
    if (fplog != NULL && pmes->cur > 0)
    {
        print_pme_loadbal_settings(pmes,fplog);
    }

    /* TODO: Here we should free all pointers in pmes,
     * but as it contains pme data structures,
     * we need to first make pme.c free all data.
     */
}
