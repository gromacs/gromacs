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

#include "state.h"

#include <cstring>

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
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
