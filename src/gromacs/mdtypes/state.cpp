/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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

#include "state.h"

#include <cstring>

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/math/veccompare.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/smalloc.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

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
    /* Ion/water position swapping */
    swapstate->eSwapCoords            = 0;
    swapstate->nIonTypes              = 0;
    swapstate->nAverage               = 0;
    swapstate->fluxleak               = 0;
    swapstate->fluxleak_p             = NULL;
    swapstate->bFromCpt               = 0;
    swapstate->nat[eChan0]            = 0;
    swapstate->nat[eChan1]            = 0;
    swapstate->xc_old_whole[eChan0]   = NULL;
    swapstate->xc_old_whole[eChan1]   = NULL;
    swapstate->xc_old_whole_p[eChan0] = NULL;
    swapstate->xc_old_whole_p[eChan1] = NULL;
    swapstate->ionType                = NULL;
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


void init_state(t_state *state, int natoms, int ngtc, int nnhpres, int nhchainlength, int dfhistNumLambda)
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
        /* We need to allocate one element extra, since we might use
         * (unaligned) 4-wide SIMD loads to access rvec entries.
         */
        snew(state->x, state->nalloc + 1);
        snew(state->v, state->nalloc + 1);
    }
    else
    {
        state->x = NULL;
        state->v = NULL;
    }
    state->cg_p = NULL;
    zero_history(&state->hist);
    zero_ekinstate(&state->ekinstate);
    snew(state->enerhist, 1);
    init_energyhistory(state->enerhist);
    if (dfhistNumLambda > 0)
    {
        snew(state->dfhist, 1);
        init_df_history(state->dfhist, dfhistNumLambda);
    }
    else
    {
        state->dfhist = NULL;
    }
    state->swapstate       = NULL;
    state->edsamstate      = NULL;
    state->ddp_count       = 0;
    state->ddp_count_cg_gl = 0;
    state->cg_gl           = NULL;
    state->cg_gl_nalloc    = 0;
}

void done_state(t_state *state)
{
    sfree(state->x);
    sfree(state->v);
    sfree(state->cg_p);
    sfree(state->enerhist);
    state->nalloc = 0;
    sfree(state->cg_gl);
    state->cg_gl_nalloc = 0;
    sfree(state->lambda);
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

void comp_state(const t_state *st1, const t_state *st2,
                gmx_bool bRMSD, real ftol, real abstol)
{
    int i, j, nc;

    fprintf(stdout, "comparing flags\n");
    cmp_int(stdout, "flags", -1, st1->flags, st2->flags);
    fprintf(stdout, "comparing box\n");
    cmp_rvecs(stdout, "box", DIM, st1->box, st2->box, FALSE, ftol, abstol);
    fprintf(stdout, "comparing box_rel\n");
    cmp_rvecs(stdout, "box_rel", DIM, st1->box_rel, st2->box_rel, FALSE, ftol, abstol);
    fprintf(stdout, "comparing boxv\n");
    cmp_rvecs(stdout, "boxv", DIM, st1->boxv, st2->boxv, FALSE, ftol, abstol);
    if (st1->flags & (1<<estSVIR_PREV))
    {
        fprintf(stdout, "comparing shake vir_prev\n");
        cmp_rvecs(stdout, "svir_prev", DIM, st1->svir_prev, st2->svir_prev, FALSE, ftol, abstol);
    }
    if (st1->flags & (1<<estFVIR_PREV))
    {
        fprintf(stdout, "comparing force vir_prev\n");
        cmp_rvecs(stdout, "fvir_prev", DIM, st1->fvir_prev, st2->fvir_prev, FALSE, ftol, abstol);
    }
    if (st1->flags & (1<<estPRES_PREV))
    {
        fprintf(stdout, "comparing prev_pres\n");
        cmp_rvecs(stdout, "pres_prev", DIM, st1->pres_prev, st2->pres_prev, FALSE, ftol, abstol);
    }
    cmp_int(stdout, "ngtc", -1, st1->ngtc, st2->ngtc);
    cmp_int(stdout, "nhchainlength", -1, st1->nhchainlength, st2->nhchainlength);
    if (st1->ngtc == st2->ngtc && st1->nhchainlength == st2->nhchainlength)
    {
        for (i = 0; i < st1->ngtc; i++)
        {
            nc = i*st1->nhchainlength;
            for (j = 0; j < nc; j++)
            {
                cmp_real(stdout, "nosehoover_xi",
                         i, st1->nosehoover_xi[nc+j], st2->nosehoover_xi[nc+j], ftol, abstol);
            }
        }
    }
    cmp_int(stdout, "nnhpres", -1, st1->nnhpres, st2->nnhpres);
    if (st1->nnhpres == st2->nnhpres && st1->nhchainlength == st2->nhchainlength)
    {
        for (i = 0; i < st1->nnhpres; i++)
        {
            nc = i*st1->nhchainlength;
            for (j = 0; j < nc; j++)
            {
                cmp_real(stdout, "nosehoover_xi",
                         i, st1->nhpres_xi[nc+j], st2->nhpres_xi[nc+j], ftol, abstol);
            }
        }
    }

    cmp_int(stdout, "natoms", -1, st1->natoms, st2->natoms);
    if (st1->natoms == st2->natoms)
    {
        if ((st1->flags & (1<<estX)) && (st2->flags & (1<<estX)))
        {
            fprintf(stdout, "comparing x\n");
            cmp_rvecs(stdout, "x", st1->natoms, st1->x, st2->x, bRMSD, ftol, abstol);
        }
        if ((st1->flags & (1<<estV)) && (st2->flags & (1<<estV)))
        {
            fprintf(stdout, "comparing v\n");
            cmp_rvecs(stdout, "v", st1->natoms, st1->v, st2->v, bRMSD, ftol, abstol);
        }
    }
}
