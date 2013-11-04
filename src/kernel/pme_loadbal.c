/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2011, The GROMACS development team,
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

#include "smalloc.h"
#include "network.h"
#include "calcgrid.h"
#include "pme.h"
#include "vec.h"
#include "domdec.h"
#include "nbnxn_cuda_data_mgmt.h"
#include "force.h"
#include "md_logging.h"
#include "pme_loadbal.h"

/* Parameters and setting for one PP-PME setup */
typedef struct {
    real      rcut_coulomb;    /* Coulomb cut-off                              */
    real      rlist;           /* pair-list cut-off                            */
    real      rlistlong;       /* LR pair-list cut-off                         */
    int       nstcalclr;       /* frequency of evaluating long-range forces for group scheme */
    real      spacing;         /* (largest) PME grid spacing                   */
    ivec      grid;            /* the PME grid dimensions                      */
    real      grid_efficiency; /* ineffiency factor for non-uniform grids <= 1 */
    real      ewaldcoeff;      /* the Ewald coefficient                        */
    gmx_pme_t pmedata;         /* the data structure used in the PME code      */

    int       count;           /* number of times this setup has been timed    */
    double    cycles;          /* the fastest time for this setup in cycles    */
} pme_setup_t;

/* In the initial scan, step by grids that are at least a factor 0.8 coarser */
#define PME_LB_GRID_SCALE_FAC  0.8
/* In the initial scan, try to skip grids with uneven x/y/z spacing,
 * checking if the "efficiency" is more than 5% worse than the previous grid.
 */
#define PME_LB_GRID_EFFICIENCY_REL_FAC  1.05
/* Rerun up till 12% slower setups than the fastest up till now */
#define PME_LB_SLOW_FAC  1.12
/* If setups get more than 2% faster, do another round to avoid
 * choosing a slower setup due to acceleration or fluctuations.
 */
#define PME_LB_ACCEL_TOL 1.02

enum {
    epmelblimNO, epmelblimBOX, epmelblimDD, epmelblimPMEGRID, epmelblimNR
};

const char *pmelblim_str[epmelblimNR] =
{ "no", "box size", "domain decompostion", "PME grid restriction" };

struct pme_load_balancing {
    int          nstage;             /* the current maximum number of stages */

    real         cut_spacing;        /* the minimum cutoff / PME grid spacing ratio */
    real         rcut_vdw;           /* Vdw cutoff (does not change) */
    real         rcut_coulomb_start; /* Initial electrostatics cutoff */
    int          nstcalclr_start;    /* Initial electrostatics cutoff */
    real         rbuf_coulomb;       /* the pairlist buffer size */
    real         rbuf_vdw;           /* the pairlist buffer size */
    matrix       box_start;          /* the initial simulation box */
    int          n;                  /* the count of setup as well as the allocation size */
    pme_setup_t *setup;              /* the PME+cutoff setups */
    int          cur;                /* the current setup */
    int          fastest;            /* fastest setup up till now */
    int          start;              /* start of setup range to consider in stage>0 */
    int          end;                /* end   of setup range to consider in stage>0 */
    int          elimited;           /* was the balancing limited, uses enum above */
    int          cutoff_scheme;      /* Verlet or group cut-offs */

    int          stage;              /* the current stage */
};

void pme_loadbal_init(pme_load_balancing_t *pme_lb_p,
                      const t_inputrec *ir, matrix box,
                      const interaction_const_t *ic,
                      gmx_pme_t pmedata)
{
    pme_load_balancing_t pme_lb;
    real                 spm, sp;
    int                  d;

    snew(pme_lb, 1);

    /* Any number of stages >= 2 is supported */
    pme_lb->nstage   = 2;

    pme_lb->cutoff_scheme = ir->cutoff_scheme;

    if (pme_lb->cutoff_scheme == ecutsVERLET)
    {
        pme_lb->rbuf_coulomb = ic->rlist - ic->rcoulomb;
        pme_lb->rbuf_vdw     = pme_lb->rbuf_coulomb;
    }
    else
    {
        if (ic->rcoulomb > ic->rlist)
        {
            pme_lb->rbuf_coulomb = ic->rlistlong - ic->rcoulomb;
        }
        else
        {
            pme_lb->rbuf_coulomb = ic->rlist - ic->rcoulomb;
        }
        if (ic->rvdw > ic->rlist)
        {
            pme_lb->rbuf_vdw = ic->rlistlong - ic->rvdw;
        }
        else
        {
            pme_lb->rbuf_vdw = ic->rlist - ic->rvdw;
        }
    }

    copy_mat(box, pme_lb->box_start);
    if (ir->ePBC == epbcXY && ir->nwall == 2)
    {
        svmul(ir->wall_ewald_zfac, pme_lb->box_start[ZZ], pme_lb->box_start[ZZ]);
    }

    pme_lb->n = 1;
    snew(pme_lb->setup, pme_lb->n);

    pme_lb->rcut_vdw              = ic->rvdw;
    pme_lb->rcut_coulomb_start    = ir->rcoulomb;
    pme_lb->nstcalclr_start       = ir->nstcalclr;

    pme_lb->cur                   = 0;
    pme_lb->setup[0].rcut_coulomb = ic->rcoulomb;
    pme_lb->setup[0].rlist        = ic->rlist;
    pme_lb->setup[0].rlistlong    = ic->rlistlong;
    pme_lb->setup[0].nstcalclr    = ir->nstcalclr;
    pme_lb->setup[0].grid[XX]     = ir->nkx;
    pme_lb->setup[0].grid[YY]     = ir->nky;
    pme_lb->setup[0].grid[ZZ]     = ir->nkz;
    pme_lb->setup[0].ewaldcoeff   = ic->ewaldcoeff;

    pme_lb->setup[0].pmedata  = pmedata;

    spm = 0;
    for (d = 0; d < DIM; d++)
    {
        sp = norm(pme_lb->box_start[d])/pme_lb->setup[0].grid[d];
        if (sp > spm)
        {
            spm = sp;
        }
    }
    pme_lb->setup[0].spacing = spm;

    if (ir->fourier_spacing > 0)
    {
        pme_lb->cut_spacing = ir->rcoulomb/ir->fourier_spacing;
    }
    else
    {
        pme_lb->cut_spacing = ir->rcoulomb/pme_lb->setup[0].spacing;
    }

    pme_lb->stage = 0;

    pme_lb->fastest  = 0;
    pme_lb->start    = 0;
    pme_lb->end      = 0;
    pme_lb->elimited = epmelblimNO;

    *pme_lb_p = pme_lb;
}

static gmx_bool pme_loadbal_increase_cutoff(pme_load_balancing_t pme_lb,
                                            int                  pme_order,
                                            const gmx_domdec_t   *dd)
{
    pme_setup_t *set;
    int          npmenodes_x, npmenodes_y;
    real         fac, sp;
    real         tmpr_coulomb, tmpr_vdw;
    int          d;
    gmx_bool     grid_ok;

    /* Try to add a new setup with next larger cut-off to the list */
    pme_lb->n++;
    srenew(pme_lb->setup, pme_lb->n);
    set          = &pme_lb->setup[pme_lb->n-1];
    set->pmedata = NULL;

    get_pme_nnodes(dd, &npmenodes_x, &npmenodes_y);

    fac = 1;
    do
    {
        /* Avoid infinite while loop, which can occur at the minimum grid size.
         * Note that in practice load balancing will stop before this point.
         * The factor 2.1 allows for the extreme case in which only grids
         * of powers of 2 are allowed (the current code supports more grids).
         */
        if (fac > 2.1)
        {
            pme_lb->n--;

            return FALSE;
        }

        fac *= 1.01;
        clear_ivec(set->grid);
        sp = calc_grid(NULL, pme_lb->box_start,
                       fac*pme_lb->setup[pme_lb->cur].spacing,
                       &set->grid[XX],
                       &set->grid[YY],
                       &set->grid[ZZ]);

        /* As here we can't easily check if one of the PME nodes
         * uses threading, we do a conservative grid check.
         * This means we can't use pme_order or less grid lines
         * per PME node along x, which is not a strong restriction.
         */
        gmx_pme_check_restrictions(pme_order,
                                   set->grid[XX], set->grid[YY], set->grid[ZZ],
                                   npmenodes_x, npmenodes_y,
                                   TRUE,
                                   FALSE,
                                   &grid_ok);
    }
    while (sp <= 1.001*pme_lb->setup[pme_lb->cur].spacing || !grid_ok);

    set->rcut_coulomb = pme_lb->cut_spacing*sp;

    if (pme_lb->cutoff_scheme == ecutsVERLET)
    {
        set->rlist        = set->rcut_coulomb + pme_lb->rbuf_coulomb;
        /* We dont use LR lists with Verlet, but this avoids if-statements in further checks */
        set->rlistlong    = set->rlist;
    }
    else
    {
        tmpr_coulomb          = set->rcut_coulomb + pme_lb->rbuf_coulomb;
        tmpr_vdw              = pme_lb->rcut_vdw + pme_lb->rbuf_vdw;
        set->rlist            = min(tmpr_coulomb, tmpr_vdw);
        set->rlistlong        = max(tmpr_coulomb, tmpr_vdw);

        /* Set the long-range update frequency */
        if (set->rlist == set->rlistlong)
        {
            /* No long-range interactions if the short-/long-range cutoffs are identical */
            set->nstcalclr = 0;
        }
        else if (pme_lb->nstcalclr_start == 0 || pme_lb->nstcalclr_start == 1)
        {
            /* We were not doing long-range before, but now we are since rlist!=rlistlong */
            set->nstcalclr = 1;
        }
        else
        {
            /* We were already doing long-range interactions from the start */
            if (pme_lb->rcut_vdw > pme_lb->rcut_coulomb_start)
            {
                /* We were originally doing long-range VdW-only interactions.
                 * If rvdw is still longer than rcoulomb we keep the original nstcalclr,
                 * but if the coulomb cutoff has become longer we should update the long-range
                 * part every step.
                 */
                set->nstcalclr = (tmpr_vdw > tmpr_coulomb) ? pme_lb->nstcalclr_start : 1;
            }
            else
            {
                /* We were not doing any long-range interaction from the start,
                 * since it is not possible to do twin-range coulomb for the PME interaction.
                 */
                set->nstcalclr = 1;
            }
        }
    }

    set->spacing      = sp;
    /* The grid efficiency is the size wrt a grid with uniform x/y/z spacing */
    set->grid_efficiency = 1;
    for (d = 0; d < DIM; d++)
    {
        set->grid_efficiency *= (set->grid[d]*sp)/norm(pme_lb->box_start[d]);
    }
    /* The Ewald coefficient is inversly proportional to the cut-off */
    set->ewaldcoeff =
        pme_lb->setup[0].ewaldcoeff*pme_lb->setup[0].rcut_coulomb/set->rcut_coulomb;

    set->count   = 0;
    set->cycles  = 0;

    if (debug)
    {
        fprintf(debug, "PME loadbal: grid %d %d %d, coulomb cutoff %f\n",
                set->grid[XX], set->grid[YY], set->grid[ZZ], set->rcut_coulomb);
    }
    return TRUE;
}

static void print_grid(FILE *fp_err, FILE *fp_log,
                       const char *pre,
                       const char *desc,
                       const pme_setup_t *set,
                       double cycles)
{
    char buf[STRLEN], buft[STRLEN];

    if (cycles >= 0)
    {
        sprintf(buft, ": %.1f M-cycles", cycles*1e-6);
    }
    else
    {
        buft[0] = '\0';
    }
    sprintf(buf, "%-11s%10s pme grid %d %d %d, coulomb cutoff %.3f%s",
            pre,
            desc, set->grid[XX], set->grid[YY], set->grid[ZZ], set->rcut_coulomb,
            buft);
    if (fp_err != NULL)
    {
        fprintf(fp_err, "\r%s\n", buf);
    }
    if (fp_log != NULL)
    {
        fprintf(fp_log, "%s\n", buf);
    }
}

static int pme_loadbal_end(pme_load_balancing_t pme_lb)
{
    /* In the initial stage only n is set; end is not set yet */
    if (pme_lb->end > 0)
    {
        return pme_lb->end;
    }
    else
    {
        return pme_lb->n;
    }
}

static void print_loadbal_limited(FILE *fp_err, FILE *fp_log,
                                  gmx_large_int_t step,
                                  pme_load_balancing_t pme_lb)
{
    char buf[STRLEN], sbuf[22];

    sprintf(buf, "step %4s: the %s limits the PME load balancing to a coulomb cut-off of %.3f",
            gmx_step_str(step, sbuf),
            pmelblim_str[pme_lb->elimited],
            pme_lb->setup[pme_loadbal_end(pme_lb)-1].rcut_coulomb);
    if (fp_err != NULL)
    {
        fprintf(fp_err, "\r%s\n", buf);
    }
    if (fp_log != NULL)
    {
        fprintf(fp_log, "%s\n", buf);
    }
}

static void switch_to_stage1(pme_load_balancing_t pme_lb)
{
    pme_lb->start = 0;
    while (pme_lb->start+1 < pme_lb->n &&
           (pme_lb->setup[pme_lb->start].count == 0 ||
            pme_lb->setup[pme_lb->start].cycles >
            pme_lb->setup[pme_lb->fastest].cycles*PME_LB_SLOW_FAC))
    {
        pme_lb->start++;
    }
    while (pme_lb->start > 0 && pme_lb->setup[pme_lb->start-1].cycles == 0)
    {
        pme_lb->start--;
    }

    pme_lb->end = pme_lb->n;
    if (pme_lb->setup[pme_lb->end-1].count > 0 &&
        pme_lb->setup[pme_lb->end-1].cycles >
        pme_lb->setup[pme_lb->fastest].cycles*PME_LB_SLOW_FAC)
    {
        pme_lb->end--;
    }

    pme_lb->stage = 1;

    /* Next we want to choose setup pme_lb->start, but as we will increase
     * pme_ln->cur by one right after returning, we subtract 1 here.
     */
    pme_lb->cur = pme_lb->start - 1;
}

gmx_bool pme_load_balance(pme_load_balancing_t pme_lb,
                          t_commrec           *cr,
                          FILE                *fp_err,
                          FILE                *fp_log,
                          t_inputrec          *ir,
                          t_state             *state,
                          double               cycles,
                          interaction_const_t *ic,
                          nonbonded_verlet_t  *nbv,
                          gmx_pme_t           *pmedata,
                          gmx_large_int_t      step)
{
    gmx_bool     OK;
    pme_setup_t *set;
    double       cycles_fast;
    char         buf[STRLEN], sbuf[22];
    real         rtab;
    gmx_bool     bUsesSimpleTables = TRUE;

    if (pme_lb->stage == pme_lb->nstage)
    {
        return FALSE;
    }

    if (PAR(cr))
    {
        gmx_sumd(1, &cycles, cr);
        cycles /= cr->nnodes;
    }

    set = &pme_lb->setup[pme_lb->cur];
    set->count++;

    rtab = ir->rlistlong + ir->tabext;

    if (set->count % 2 == 1)
    {
        /* Skip the first cycle, because the first step after a switch
         * is much slower due to allocation and/or caching effects.
         */
        return TRUE;
    }

    sprintf(buf, "step %4s: ", gmx_step_str(step, sbuf));
    print_grid(fp_err, fp_log, buf, "timed with", set, cycles);

    if (set->count <= 2)
    {
        set->cycles = cycles;
    }
    else
    {
        if (cycles*PME_LB_ACCEL_TOL < set->cycles &&
            pme_lb->stage == pme_lb->nstage - 1)
        {
            /* The performance went up a lot (due to e.g. DD load balancing).
             * Add a stage, keep the minima, but rescan all setups.
             */
            pme_lb->nstage++;

            if (debug)
            {
                fprintf(debug, "The performance for grid %d %d %d went from %.3f to %.1f M-cycles, this is more than %f\n"
                        "Increased the number stages to %d"
                        " and ignoring the previous performance\n",
                        set->grid[XX], set->grid[YY], set->grid[ZZ],
                        cycles*1e-6, set->cycles*1e-6, PME_LB_ACCEL_TOL,
                        pme_lb->nstage);
            }
        }
        set->cycles = min(set->cycles, cycles);
    }

    if (set->cycles < pme_lb->setup[pme_lb->fastest].cycles)
    {
        pme_lb->fastest = pme_lb->cur;

        if (DOMAINDECOMP(cr))
        {
            /* We found a new fastest setting, ensure that with subsequent
             * shorter cut-off's the dynamic load balancing does not make
             * the use of the current cut-off impossible. This solution is
             * a trade-off, as the PME load balancing and DD domain size
             * load balancing can interact in complex ways.
             * With the Verlet kernels, DD load imbalance will usually be
             * mainly due to bonded interaction imbalance, which will often
             * quickly push the domain boundaries beyond the limit for the
             * optimal, PME load balanced, cut-off. But it could be that
             * better overal performance can be obtained with a slightly
             * shorter cut-off and better DD load balancing.
             */
            change_dd_dlb_cutoff_limit(cr);
        }
    }
    cycles_fast = pme_lb->setup[pme_lb->fastest].cycles;

    /* Check in stage 0 if we should stop scanning grids.
     * Stop when the time is more than SLOW_FAC longer than the fastest.
     */
    if (pme_lb->stage == 0 && pme_lb->cur > 0 &&
        cycles > pme_lb->setup[pme_lb->fastest].cycles*PME_LB_SLOW_FAC)
    {
        pme_lb->n = pme_lb->cur + 1;
        /* Done with scanning, go to stage 1 */
        switch_to_stage1(pme_lb);
    }

    if (pme_lb->stage == 0)
    {
        int gridsize_start;

        gridsize_start = set->grid[XX]*set->grid[YY]*set->grid[ZZ];

        do
        {
            if (pme_lb->cur+1 < pme_lb->n)
            {
                /* We had already generated the next setup */
                OK = TRUE;
            }
            else
            {
                /* Find the next setup */
                OK = pme_loadbal_increase_cutoff(pme_lb, ir->pme_order, cr->dd);

                if (!OK)
                {
                    pme_lb->elimited = epmelblimPMEGRID;
                }
            }

            if (OK && ir->ePBC != epbcNONE)
            {
                OK = (sqr(pme_lb->setup[pme_lb->cur+1].rlistlong)
                      <= max_cutoff2(ir->ePBC, state->box));
                if (!OK)
                {
                    pme_lb->elimited = epmelblimBOX;
                }
            }

            if (OK)
            {
                pme_lb->cur++;

                if (DOMAINDECOMP(cr))
                {
                    OK = change_dd_cutoff(cr, state, ir,
                                          pme_lb->setup[pme_lb->cur].rlistlong);
                    if (!OK)
                    {
                        /* Failed: do not use this setup */
                        pme_lb->cur--;
                        pme_lb->elimited = epmelblimDD;
                    }
                }
            }
            if (!OK)
            {
                /* We hit the upper limit for the cut-off,
                 * the setup should not go further than cur.
                 */
                pme_lb->n = pme_lb->cur + 1;
                print_loadbal_limited(fp_err, fp_log, step, pme_lb);
                /* Switch to the next stage */
                switch_to_stage1(pme_lb);
            }
        }
        while (OK &&
               !(pme_lb->setup[pme_lb->cur].grid[XX]*
                 pme_lb->setup[pme_lb->cur].grid[YY]*
                 pme_lb->setup[pme_lb->cur].grid[ZZ] <
                 gridsize_start*PME_LB_GRID_SCALE_FAC
                 &&
                 pme_lb->setup[pme_lb->cur].grid_efficiency <
                 pme_lb->setup[pme_lb->cur-1].grid_efficiency*PME_LB_GRID_EFFICIENCY_REL_FAC));
    }

    if (pme_lb->stage > 0 && pme_lb->end == 1)
    {
        pme_lb->cur   = 0;
        pme_lb->stage = pme_lb->nstage;
    }
    else if (pme_lb->stage > 0 && pme_lb->end > 1)
    {
        /* If stage = nstage-1:
         *   scan over all setups, rerunning only those setups
         *   which are not much slower than the fastest
         * else:
         *   use the next setup
         */
        do
        {
            pme_lb->cur++;
            if (pme_lb->cur == pme_lb->end)
            {
                pme_lb->stage++;
                pme_lb->cur = pme_lb->start;
            }
        }
        while (pme_lb->stage == pme_lb->nstage - 1 &&
               pme_lb->setup[pme_lb->cur].count > 0 &&
               pme_lb->setup[pme_lb->cur].cycles > cycles_fast*PME_LB_SLOW_FAC);

        if (pme_lb->stage == pme_lb->nstage)
        {
            /* We are done optimizing, use the fastest setup we found */
            pme_lb->cur = pme_lb->fastest;
        }
    }

    if (DOMAINDECOMP(cr) && pme_lb->stage > 0)
    {
        OK = change_dd_cutoff(cr, state, ir, pme_lb->setup[pme_lb->cur].rlistlong);
        if (!OK)
        {
            /* Failsafe solution */
            if (pme_lb->cur > 1 && pme_lb->stage == pme_lb->nstage)
            {
                pme_lb->stage--;
            }
            pme_lb->fastest  = 0;
            pme_lb->start    = 0;
            pme_lb->end      = pme_lb->cur;
            pme_lb->cur      = pme_lb->start;
            pme_lb->elimited = epmelblimDD;
            print_loadbal_limited(fp_err, fp_log, step, pme_lb);
        }
    }

    /* Change the Coulomb cut-off and the PME grid */

    set = &pme_lb->setup[pme_lb->cur];

    ic->rcoulomb   = set->rcut_coulomb;
    ic->rlist      = set->rlist;
    ic->rlistlong  = set->rlistlong;
    ir->nstcalclr  = set->nstcalclr;
    ic->ewaldcoeff = set->ewaldcoeff;

    bUsesSimpleTables = uses_simple_tables(ir->cutoff_scheme, nbv, 0);
    if (pme_lb->cutoff_scheme == ecutsVERLET &&
        nbv->grp[0].kernel_type == nbnxnk8x8x8_CUDA)
    {
        nbnxn_cuda_pme_loadbal_update_param(nbv->cu_nbv, ic);
    }
    else
    {
        init_interaction_const_tables(NULL, ic, bUsesSimpleTables,
                                      rtab);
    }

    if (pme_lb->cutoff_scheme == ecutsVERLET && nbv->ngrp > 1)
    {
        init_interaction_const_tables(NULL, ic, bUsesSimpleTables,
                                      rtab);
    }

    if (cr->duty & DUTY_PME)
    {
        if (pme_lb->setup[pme_lb->cur].pmedata == NULL)
        {
            /* Generate a new PME data structure,
             * copying part of the old pointers.
             */
            gmx_pme_reinit(&set->pmedata,
                           cr, pme_lb->setup[0].pmedata, ir,
                           set->grid);
        }
        *pmedata = set->pmedata;
    }
    else
    {
        /* Tell our PME-only node to switch grid */
        gmx_pme_send_switchgrid(cr, set->grid, set->ewaldcoeff);
    }

    if (debug)
    {
        print_grid(NULL, debug, "", "switched to", set, -1);
    }

    if (pme_lb->stage == pme_lb->nstage)
    {
        print_grid(fp_err, fp_log, "", "optimal", set, -1);
    }

    return TRUE;
}

void restart_pme_loadbal(pme_load_balancing_t pme_lb, int n)
{
    pme_lb->nstage += n;
}

static int pme_grid_points(const pme_setup_t *setup)
{
    return setup->grid[XX]*setup->grid[YY]*setup->grid[ZZ];
}

static real pme_loadbal_rlist(const pme_setup_t *setup)
{
    /* With the group cut-off scheme we can have twin-range either
     * for Coulomb or for VdW, so we need a check here.
     * With the Verlet cut-off scheme rlist=rlistlong.
     */
    if (setup->rcut_coulomb > setup->rlist)
    {
        return setup->rlistlong;
    }
    else
    {
        return setup->rlist;
    }
}

static void print_pme_loadbal_setting(FILE              *fplog,
                                      char              *name,
                                      const pme_setup_t *setup)
{
    fprintf(fplog,
            "   %-7s %6.3f nm %6.3f nm     %3d %3d %3d   %5.3f nm  %5.3f nm\n",
            name,
            setup->rcut_coulomb, pme_loadbal_rlist(setup),
            setup->grid[XX], setup->grid[YY], setup->grid[ZZ],
            setup->spacing, 1/setup->ewaldcoeff);
}

static void print_pme_loadbal_settings(pme_load_balancing_t pme_lb,
                                       t_commrec           *cr,
                                       FILE                *fplog,
                                       gmx_bool             bNonBondedOnGPU)
{
    double pp_ratio, grid_ratio;

    pp_ratio   = pow(pme_loadbal_rlist(&pme_lb->setup[pme_lb->cur])/pme_loadbal_rlist(&pme_lb->setup[0]), 3.0);
    grid_ratio = pme_grid_points(&pme_lb->setup[pme_lb->cur])/
        (double)pme_grid_points(&pme_lb->setup[0]);

    fprintf(fplog, "\n");
    fprintf(fplog, "       P P   -   P M E   L O A D   B A L A N C I N G\n");
    fprintf(fplog, "\n");
    /* Here we only warn when the optimal setting is the last one */
    if (pme_lb->elimited != epmelblimNO &&
        pme_lb->cur == pme_loadbal_end(pme_lb)-1)
    {
        fprintf(fplog, " NOTE: The PP/PME load balancing was limited by the %s,\n",
                pmelblim_str[pme_lb->elimited]);
        fprintf(fplog, "       you might not have reached a good load balance.\n");
        if (pme_lb->elimited == epmelblimDD)
        {
            fprintf(fplog, "       Try different mdrun -dd settings or lower the -dds value.\n");
        }
        fprintf(fplog, "\n");
    }
    fprintf(fplog, " PP/PME load balancing changed the cut-off and PME settings:\n");
    fprintf(fplog, "           particle-particle                    PME\n");
    fprintf(fplog, "            rcoulomb  rlist            grid      spacing   1/beta\n");
    print_pme_loadbal_setting(fplog, "initial", &pme_lb->setup[0]);
    print_pme_loadbal_setting(fplog, "final", &pme_lb->setup[pme_lb->cur]);
    fprintf(fplog, " cost-ratio           %4.2f             %4.2f\n",
            pp_ratio, grid_ratio);
    fprintf(fplog, " (note that these numbers concern only part of the total PP and PME load)\n");

    if (pp_ratio > 1.5 && !bNonBondedOnGPU)
    {
        md_print_warn(cr, fplog,
                      "NOTE: PME load balancing increased the non-bonded workload by more than 50%%.\n"
                      "      For better performance use (more) PME nodes (mdrun -npme),\n"
                      "      or in case you are beyond the scaling limit, use less nodes in total.\n");
    }
    else
    {
        fprintf(fplog, "\n");
    }
}

void pme_loadbal_done(pme_load_balancing_t pme_lb,
                      t_commrec *cr, FILE *fplog,
                      gmx_bool bNonBondedOnGPU)
{
    if (fplog != NULL && (pme_lb->cur > 0 || pme_lb->elimited != epmelblimNO))
    {
        print_pme_loadbal_settings(pme_lb, cr, fplog, bNonBondedOnGPU);
    }

    /* TODO: Here we should free all pointers in pme_lb,
     * but as it contains pme data structures,
     * we need to first make pme.c free all data.
     */
}
