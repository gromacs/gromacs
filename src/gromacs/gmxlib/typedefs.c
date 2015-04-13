/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/typedefs.h"

#include <string.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/random.h"
#include "gromacs/utility/smalloc.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

int gmx_int64_to_int(gmx_int64_t step, const char *warn)
{
    int i;

    i = (int)step;

    if (warn != NULL && (step < INT_MIN || step > INT_MAX))
    {
        fprintf(stderr, "\nWARNING during %s:\n", warn);
        fprintf(stderr, "step value ");
        fprintf(stderr, "%"GMX_PRId64, step);
        fprintf(stderr, " does not fit in int, converted to %d\n\n", i);
    }

    return i;
}

void init_inputrec(t_inputrec *ir)
{
    memset(ir, 0, (size_t)sizeof(*ir));
    snew(ir->fepvals, 1);
    snew(ir->expandedvals, 1);
    snew(ir->simtempvals, 1);
}

static void done_pull_group(t_pull_group *pgrp)
{
    if (pgrp->nat > 0)
    {
        sfree(pgrp->ind);
        sfree(pgrp->weight);
    }
}

static void done_pull_params(pull_params_t *pull)
{
    int i;

    for (i = 0; i < pull->ngroup+1; i++)
    {
        done_pull_group(pull->group);
    }

    sfree(pull->group);
    sfree(pull->coord);
}

void done_inputrec(t_inputrec *ir)
{
    int m;

    for (m = 0; (m < DIM); m++)
    {
        if (ir->ex[m].a)
        {
            sfree(ir->ex[m].a);
        }
        if (ir->ex[m].phi)
        {
            sfree(ir->ex[m].phi);
        }
        if (ir->et[m].a)
        {
            sfree(ir->et[m].a);
        }
        if (ir->et[m].phi)
        {
            sfree(ir->et[m].phi);
        }
    }

    sfree(ir->opts.nrdf);
    sfree(ir->opts.ref_t);
    sfree(ir->opts.annealing);
    sfree(ir->opts.anneal_npoints);
    sfree(ir->opts.anneal_time);
    sfree(ir->opts.anneal_temp);
    sfree(ir->opts.tau_t);
    sfree(ir->opts.acc);
    sfree(ir->opts.nFreeze);
    sfree(ir->opts.QMmethod);
    sfree(ir->opts.QMbasis);
    sfree(ir->opts.QMcharge);
    sfree(ir->opts.QMmult);
    sfree(ir->opts.bSH);
    sfree(ir->opts.CASorbitals);
    sfree(ir->opts.CASelectrons);
    sfree(ir->opts.SAon);
    sfree(ir->opts.SAoff);
    sfree(ir->opts.SAsteps);
    sfree(ir->opts.bOPT);
    sfree(ir->opts.bTS);

    if (ir->pull)
    {
        done_pull_params(ir->pull);
        sfree(ir->pull);
    }
}

static void zero_history(history_t *hist)
{
    hist->disre_initf  = 0;
    hist->ndisrepairs  = 0;
    hist->disre_rm3tav = NULL;
    hist->orire_initf  = 0;
    hist->norire_Dtav  = 0;
    hist->orire_Dtav   = NULL;
}

static void zero_ekinstate(ekinstate_t *eks)
{
    eks->ekin_n         = 0;
    eks->ekinh          = NULL;
    eks->ekinf          = NULL;
    eks->ekinh_old      = NULL;
    eks->ekinscalef_nhc = NULL;
    eks->ekinscaleh_nhc = NULL;
    eks->vscale_nhc     = NULL;
    eks->dekindl        = 0;
    eks->mvcos          = 0;
}

static void init_swapstate(swapstate_t *swapstate)
{
    int ii, ic;

    swapstate->eSwapCoords = 0;
    swapstate->nAverage    = 0;

    /* Ion/water position swapping */
    for (ic = 0; ic < eCompNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            swapstate->nat_req[ic][ii]        = 0;
            swapstate->nat_req_p[ic][ii]      = NULL;
            swapstate->inflow_netto[ic][ii]   = 0;
            swapstate->inflow_netto_p[ic][ii] = NULL;
            swapstate->nat_past[ic][ii]       = NULL;
            swapstate->nat_past_p[ic][ii]     = NULL;
            swapstate->fluxfromAtoB[ic][ii]   = 0;
            swapstate->fluxfromAtoB_p[ic][ii] = NULL;
        }
    }
    swapstate->fluxleak               = NULL;
    swapstate->nions                  = 0;
    swapstate->comp_from              = NULL;
    swapstate->channel_label          = NULL;
    swapstate->bFromCpt               = 0;
    swapstate->nat[eChan0]            = 0;
    swapstate->nat[eChan1]            = 0;
    swapstate->xc_old_whole[eChan0]   = NULL;
    swapstate->xc_old_whole[eChan1]   = NULL;
    swapstate->xc_old_whole_p[eChan0] = NULL;
    swapstate->xc_old_whole_p[eChan1] = NULL;
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

void init_gtc_state(t_state *state, int ngtc, int nnhpres, int nhchainlength)
{
    int i, j;

    state->ngtc          = ngtc;
    state->nnhpres       = nnhpres;
    state->nhchainlength = nhchainlength;
    if (state->ngtc > 0)
    {
        snew(state->nosehoover_xi, state->nhchainlength*state->ngtc);
        snew(state->nosehoover_vxi, state->nhchainlength*state->ngtc);
        snew(state->therm_integral, state->ngtc);
        for (i = 0; i < state->ngtc; i++)
        {
            for (j = 0; j < state->nhchainlength; j++)
            {
                state->nosehoover_xi[i*state->nhchainlength + j]   = 0.0;
                state->nosehoover_vxi[i*state->nhchainlength + j]  = 0.0;
            }
        }
        for (i = 0; i < state->ngtc; i++)
        {
            state->therm_integral[i]  = 0.0;
        }
    }
    else
    {
        state->nosehoover_xi  = NULL;
        state->nosehoover_vxi = NULL;
        state->therm_integral = NULL;
    }

    if (state->nnhpres > 0)
    {
        snew(state->nhpres_xi, state->nhchainlength*nnhpres);
        snew(state->nhpres_vxi, state->nhchainlength*nnhpres);
        for (i = 0; i < nnhpres; i++)
        {
            for (j = 0; j < state->nhchainlength; j++)
            {
                state->nhpres_xi[i*nhchainlength + j]   = 0.0;
                state->nhpres_vxi[i*nhchainlength + j]  = 0.0;
            }
        }
    }
    else
    {
        state->nhpres_xi  = NULL;
        state->nhpres_vxi = NULL;
    }
}


void init_state(t_state *state, int natoms, int ngtc, int nnhpres, int nhchainlength, int nlambda)
{
    int i;

    state->natoms    = natoms;
    state->flags     = 0;
    state->fep_state = 0;
    state->lambda    = 0;
    snew(state->lambda, efptNR);
    for (i = 0; i < efptNR; i++)
    {
        state->lambda[i] = 0;
    }
    state->veta   = 0;
    clear_mat(state->box);
    clear_mat(state->box_rel);
    clear_mat(state->boxv);
    clear_mat(state->pres_prev);
    clear_mat(state->svir_prev);
    clear_mat(state->fvir_prev);
    init_gtc_state(state, ngtc, nnhpres, nhchainlength);
    state->nalloc = state->natoms;
    if (state->nalloc > 0)
    {
        snew(state->x, state->nalloc);
        snew(state->v, state->nalloc);
    }
    else
    {
        state->x = NULL;
        state->v = NULL;
    }
    state->sd_X = NULL;
    state->cg_p = NULL;
    zero_history(&state->hist);
    zero_ekinstate(&state->ekinstate);
    init_energyhistory(&state->enerhist);
    init_df_history(&state->dfhist, nlambda);
    init_swapstate(&state->swapstate);
    state->ddp_count       = 0;
    state->ddp_count_cg_gl = 0;
    state->cg_gl           = NULL;
    state->cg_gl_nalloc    = 0;
}

void done_state(t_state *state)
{
    if (state->x)
    {
        sfree(state->x);
    }
    if (state->v)
    {
        sfree(state->v);
    }
    if (state->sd_X)
    {
        sfree(state->sd_X);
    }
    if (state->cg_p)
    {
        sfree(state->cg_p);
    }
    state->nalloc = 0;
    if (state->cg_gl)
    {
        sfree(state->cg_gl);
    }
    state->cg_gl_nalloc = 0;
    if (state->lambda)
    {
        sfree(state->lambda);
    }
    if (state->ngtc > 0)
    {
        sfree(state->nosehoover_xi);
        sfree(state->nosehoover_vxi);
        sfree(state->therm_integral);
    }
}

t_state *serial_init_local_state(t_state *state_global)
{
    int      i;
    t_state *state_local;

    snew(state_local, 1);

    /* Copy all the contents */
    *state_local = *state_global;
    snew(state_local->lambda, efptNR);
    /* local storage for lambda */
    for (i = 0; i < efptNR; i++)
    {
        state_local->lambda[i] = state_global->lambda[i];
    }

    return state_local;
}

static void do_box_rel(t_inputrec *ir, matrix box_rel, matrix b, gmx_bool bInit)
{
    int d, d2;

    for (d = YY; d <= ZZ; d++)
    {
        for (d2 = XX; d2 <= (ir->epct == epctSEMIISOTROPIC ? YY : ZZ); d2++)
        {
            /* We need to check if this box component is deformed
             * or if deformation of another component might cause
             * changes in this component due to box corrections.
             */
            if (ir->deform[d][d2] == 0 &&
                !(d == ZZ && d2 == XX && ir->deform[d][YY] != 0 &&
                  (b[YY][d2] != 0 || ir->deform[YY][d2] != 0)))
            {
                if (bInit)
                {
                    box_rel[d][d2] = b[d][d2]/b[XX][XX];
                }
                else
                {
                    b[d][d2] = b[XX][XX]*box_rel[d][d2];
                }
            }
        }
    }
}

void set_box_rel(t_inputrec *ir, t_state *state)
{
    /* Make sure the box obeys the restrictions before we fix the ratios */
    correct_box(NULL, 0, state->box, NULL);

    clear_mat(state->box_rel);

    if (PRESERVE_SHAPE(*ir))
    {
        do_box_rel(ir, state->box_rel, state->box, TRUE);
    }
}

void preserve_box_shape(t_inputrec *ir, matrix box_rel, matrix b)
{
    if (PRESERVE_SHAPE(*ir))
    {
        do_box_rel(ir, box_rel, b, FALSE);
    }
}

real max_cutoff(real cutoff1, real cutoff2)
{
    if (cutoff1 == 0 || cutoff2 == 0)
    {
        return 0;
    }
    else
    {
        return max(cutoff1, cutoff2);
    }
}

void init_df_history(df_history_t *dfhist, int nlambda)
{
    int i;

    dfhist->nlambda  = nlambda;
    dfhist->bEquil   = 0;
    dfhist->wl_delta = 0;

    if (nlambda > 0)
    {
        snew(dfhist->sum_weights, dfhist->nlambda);
        snew(dfhist->sum_dg, dfhist->nlambda);
        snew(dfhist->sum_minvar, dfhist->nlambda);
        snew(dfhist->sum_variance, dfhist->nlambda);
        snew(dfhist->n_at_lam, dfhist->nlambda);
        snew(dfhist->wl_histo, dfhist->nlambda);

        /* allocate transition matrices here */
        snew(dfhist->Tij, dfhist->nlambda);
        snew(dfhist->Tij_empirical, dfhist->nlambda);

        /* allocate accumulators for various transition matrix
           free energy methods here */
        snew(dfhist->accum_p, dfhist->nlambda);
        snew(dfhist->accum_m, dfhist->nlambda);
        snew(dfhist->accum_p2, dfhist->nlambda);
        snew(dfhist->accum_m2, dfhist->nlambda);

        for (i = 0; i < dfhist->nlambda; i++)
        {
            snew(dfhist->Tij[i], dfhist->nlambda);
            snew(dfhist->Tij_empirical[i], dfhist->nlambda);
            snew((dfhist->accum_p)[i], dfhist->nlambda);
            snew((dfhist->accum_m)[i], dfhist->nlambda);
            snew((dfhist->accum_p2)[i], dfhist->nlambda);
            snew((dfhist->accum_m2)[i], dfhist->nlambda);
        }
    }
}

extern void copy_df_history(df_history_t *df_dest, df_history_t *df_source)
{
    int i, j;

    /* Currently, there should not be any difference in nlambda between the two,
       but this is included for completeness for potential later functionality */
    df_dest->nlambda  = df_source->nlambda;
    df_dest->bEquil   = df_source->bEquil;
    df_dest->wl_delta = df_source->wl_delta;

    for (i = 0; i < df_dest->nlambda; i++)
    {
        df_dest->sum_weights[i]  = df_source->sum_weights[i];
        df_dest->sum_dg[i]       = df_source->sum_dg[i];
        df_dest->sum_minvar[i]   = df_source->sum_minvar[i];
        df_dest->sum_variance[i] = df_source->sum_variance[i];
        df_dest->n_at_lam[i]     = df_source->n_at_lam[i];
        df_dest->wl_histo[i]     = df_source->wl_histo[i];
    }

    for (i = 0; i < df_dest->nlambda; i++)
    {
        for (j = 0; j < df_dest->nlambda; j++)
        {
            df_dest->accum_p[i][j]        = df_source->accum_p[i][j];
            df_dest->accum_m[i][j]        = df_source->accum_m[i][j];
            df_dest->accum_p2[i][j]       = df_source->accum_p2[i][j];
            df_dest->accum_m2[i][j]       = df_source->accum_m2[i][j];
            df_dest->Tij[i][j]            = df_source->Tij[i][j];
            df_dest->Tij_empirical[i][j]  = df_source->Tij_empirical[i][j];
        }
    }
}

void done_df_history(df_history_t *dfhist)
{
    int i;

    if (dfhist->nlambda > 0)
    {
        sfree(dfhist->n_at_lam);
        sfree(dfhist->wl_histo);
        sfree(dfhist->sum_weights);
        sfree(dfhist->sum_dg);
        sfree(dfhist->sum_minvar);
        sfree(dfhist->sum_variance);

        for (i = 0; i < dfhist->nlambda; i++)
        {
            sfree(dfhist->Tij[i]);
            sfree(dfhist->Tij_empirical[i]);
            sfree(dfhist->accum_p[i]);
            sfree(dfhist->accum_m[i]);
            sfree(dfhist->accum_p2[i]);
            sfree(dfhist->accum_m2[i]);
        }
    }
    dfhist->bEquil   = 0;
    dfhist->nlambda  = 0;
    dfhist->wl_delta = 0;
}
