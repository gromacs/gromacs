/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "history.h"

#include <assert.h>

#include <cstring>

#include <algorithm>

#include "gromacs/awh/awh.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "correlation-history.h"
#include "grid.h"
#include "types.h"

/* These are just dummy initializers used when initializing the state because all AWH parameters
   are not (readily) available then/there (if the awh struct hasn't been initialized). TODO. */
static void init_bias_history(awh_bias_history_t *bias_history)
{
    bias_history->in_initial                = 0;
    bias_history->histsize                  = 0;
    bias_history->npoints                   = 0;
    bias_history->coordpoint                = NULL;
    bias_history->coord_refvalue_index      = 0;
    bias_history->ndim                      = 0;
    bias_history->bForce_correlation        = FALSE;
    bias_history->forcecorr_hist            = NULL;
    bias_history->log_relative_sampleweight = 0;
}

static void init_bias_history_from_state(awh_bias_history_t *bias_history, const awh_bias_t *bias)
{
    init_bias_history(bias_history);

    bias_history->npoints            = bias->npoints;
    bias_history->ndim               = bias->ndim;

    snew(bias_history->coordpoint, bias_history->npoints);

    if (bias_history->bForce_correlation)
    {
        bias_history->forcecorr_hist = init_correlation_grid_history_from_state(bias->forcecorr);
    }
}

void init_awh_history(awh_history_t *awh_history)
{
    awh_history->nbias                              = 0;
    awh_history->awh_bias_history                   = NULL;
    awh_history->convolved_bias_shift               = 0;
}

void init_awh_history_from_state(awh_history_t *awh_history, const awh_t *awh)
{
    awh_history->convolved_bias_shift      = awh->convolved_bias_shift;

    awh_history->nbias = awh->nbias;
    snew(awh_history->awh_bias_history, awh_history->nbias);
    for (int k = 0; k < awh_history->nbias; k++)
    {
        init_bias_history_from_state(&awh_history->awh_bias_history[k], &awh->awh_bias[k]);
    }
}

static void update_bias_history(awh_bias_history_t *bias_history, awh_bias_t *bias)
{
    GMX_RELEASE_ASSERT(bias_history->coordpoint != NULL, "AWH history coord point not initialized when updating history");
    GMX_RELEASE_ASSERT(bias_history->npoints > 0, "AWH history number of points not initialized when updating history");
    GMX_RELEASE_ASSERT(bias_history->ndim > 0, "AWH history number of dimensions not initialized when updating history");

    bias_history->coord_refvalue_index = bias->coord_refvalue_index;

    for (int m = 0; m < bias_history->npoints; m++)
    {
        bias_history->coordpoint[m].target                               = bias->coordpoint[m].target;
        bias_history->coordpoint[m].free_energy                          = bias->coordpoint[m].free_energy;
        bias_history->coordpoint[m].bias                                 = bias->coordpoint[m].bias;
        bias_history->coordpoint[m].weightsum_iteration                  = bias->coordpoint[m].weightsum_iteration;
        bias_history->coordpoint[m].weightsum_covering                   = bias->coordpoint[m].weightsum_covering;
        bias_history->coordpoint[m].weightsum_tot                        = bias->coordpoint[m].weightsum_tot;
        bias_history->coordpoint[m].weightsum_ref                        = bias->coordpoint[m].weightsum_ref;
        bias_history->coordpoint[m].last_update_index                    = bias->coordpoint[m].last_update_index;
        bias_history->coordpoint[m].log_pmfsum                           = bias->coordpoint[m].log_pmfsum;
        bias_history->coordpoint[m].visits_iteration                     = bias->coordpoint[m].visits_iteration;
        bias_history->coordpoint[m].visits_tot                           = bias->coordpoint[m].visits_tot;
    }

    bias_history->in_initial             = bias->in_initial;
    bias_history->histsize               = bias->histsize;

    bias_history->origin_index_updatelist = multidim_gridindex_to_linear(bias->grid,
                                                                         bias->origin_updatelist);
    bias_history->end_index_updatelist = multidim_gridindex_to_linear(bias->grid,
                                                                      bias->end_updatelist);

    if (bias->bForce_correlation)
    {
        GMX_RELEASE_ASSERT(bias_history->forcecorr_hist != NULL, "AWH history force correlation not initialized when updating history");
        update_correlation_grid_history(bias_history->forcecorr_hist, bias->forcecorr);
    }

    bias_history->log_relative_sampleweight  = bias->log_relative_sampleweight;
}

static void restore_bias_state_from_history(const awh_bias_history_t *bias_history, awh_bias_t *bias)
{
    bias->coord_refvalue_index = bias_history->coord_refvalue_index;

    for (int m = 0; m < bias->npoints; m++)
    {
        bias->coordpoint[m].target                                = bias_history->coordpoint[m].target;
        bias->coordpoint[m].free_energy                           = bias_history->coordpoint[m].free_energy;
        bias->coordpoint[m].bias                                  = bias_history->coordpoint[m].bias;
        bias->coordpoint[m].weightsum_iteration                   = bias_history->coordpoint[m].weightsum_iteration;
        bias->coordpoint[m].weightsum_covering                    = bias_history->coordpoint[m].weightsum_covering;
        bias->coordpoint[m].weightsum_tot                         = bias_history->coordpoint[m].weightsum_tot;
        bias->coordpoint[m].weightsum_ref                         = bias_history->coordpoint[m].weightsum_ref;
        bias->coordpoint[m].last_update_index                     = bias_history->coordpoint[m].last_update_index;
        bias->coordpoint[m].log_pmfsum                            = bias_history->coordpoint[m].log_pmfsum;
        bias->coordpoint[m].visits_iteration                      = bias_history->coordpoint[m].visits_iteration;
        bias->coordpoint[m].visits_tot                            = bias_history->coordpoint[m].visits_tot;
    }

    bias->in_initial      = bias_history->in_initial;
    bias->histsize        = bias_history->histsize;

    linear_gridindex_to_multidim(bias->grid, bias_history->origin_index_updatelist, bias->origin_updatelist);
    linear_gridindex_to_multidim(bias->grid, bias_history->end_index_updatelist, bias->end_updatelist);

    if (bias->bForce_correlation)
    {
        restore_correlation_grid_state_from_history(bias_history->forcecorr_hist, bias->forcecorr);
    }

    bias->log_relative_sampleweight  = bias_history->log_relative_sampleweight;
}

void update_awh_history(awh_history_t *awh_history, const awh_t *awh)
{
    GMX_RELEASE_ASSERT(awh_history->awh_bias_history != NULL, "AWH history not initialized when updating history");

    awh_history->convolved_bias_shift = awh->convolved_bias_shift;

    for (int k = 0; k < awh_history->nbias; k++)
    {
        update_bias_history(&awh_history->awh_bias_history[k], &awh->awh_bias[k]);
    }
}

static void broadcast_initialize_bias_history(awh_bias_history_t *bias_history, const t_commrec *cr)
{
    gmx_bcast(sizeof(awh_bias_history_t), bias_history, cr);

    if (!MASTER(cr))
    {
        snew(bias_history->coordpoint, bias_history->npoints);
    }
    gmx_bcast(sizeof(awh_coordpoint_history_t)*bias_history->npoints, bias_history->coordpoint, cr);

    if (bias_history->bForce_correlation)
    {
        correlation_grid_history_t *corrgrid_hist = init_correlation_grid_history_from_checkpoint(bias_history->forcecorr_hist, cr);

        if (!MASTER(cr))
        {
            bias_history->forcecorr_hist = corrgrid_hist;
        }
    }
}

static void broadcast_initialize_awh_history(awh_history_t *awh_history, const t_commrec *cr)
{
    gmx_bcast(sizeof(awh_history_t), awh_history, cr);
    GMX_RELEASE_ASSERT(awh_history->nbias > 0, "nbias needs to be > 0");

    if (!MASTER(cr))
    {
        snew(awh_history->awh_bias_history, awh_history->nbias);
    }

    for (int k = 0; k < awh_history->nbias; k++)
    {
        broadcast_initialize_bias_history(&awh_history->awh_bias_history[k], cr);
    }
}

void init_awh_history_from_checkpoint(awh_history_t *awh_history, const t_commrec *cr)
{
    /* Only master has read the checkpoint and only master has an initialized AWH history struct.
       Non-master ranks need to get the their history allocated and initialized by broadcasting from the master rank. */
    if (PAR(cr))
    {
        GMX_RELEASE_ASSERT((MASTER(cr) && awh_history->awh_bias_history != NULL) || (!MASTER(cr) && awh_history->awh_bias_history == NULL),
                           "(Non-)master rank(s) should (not) have been initialized before restoring awh_history");
        broadcast_initialize_awh_history(awh_history, cr);
    }
}

void restore_awh_state_from_history(const awh_history_t *awh_history, awh_t *awh)
{
    GMX_RELEASE_ASSERT(awh_history->awh_bias_history != NULL, "AWH history was not initialized when attempting to restore state from history");
    GMX_RELEASE_ASSERT(awh != NULL, "AWH bias not initialized when attempting to restore state from history");

    /* Restore the history to the current state */
    awh->convolved_bias_shift       = awh_history->convolved_bias_shift;

    for (int k = 0; k < awh_history->nbias; k++)
    {
        restore_bias_state_from_history(&awh_history->awh_bias_history[k], &awh->awh_bias[k]);
    }
}
