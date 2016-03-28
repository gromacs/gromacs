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
static void init_awh_bias_history(awh_bias_history_t *awh_bias_history)
{
    awh_bias_history->in_initial                = 0;
    awh_bias_history->histsize                  = 0;
    awh_bias_history->npoints                   = 0;
    awh_bias_history->coordpoint                = NULL;
    awh_bias_history->coord_refvalue_index      = 0;
    awh_bias_history->ndim                      = 0;
    awh_bias_history->bForce_correlation        = FALSE;
    awh_bias_history->forcecorr_hist            = NULL;
    awh_bias_history->log_relative_sampleweight = 0;
}

static void init_awh_bias_history_from_state(awh_bias_history_t *awh_bias_history, const awh_bias_t *awh_bias)
{
    init_awh_bias_history(awh_bias_history);

    awh_bias_history->npoints            = awh_bias->npoints;
    awh_bias_history->ndim               = awh_bias->ndim;
    awh_bias_history->bForce_correlation = awh_bias->bForce_correlation;

    snew(awh_bias_history->coordpoint, awh_bias_history->npoints);

    if (awh_bias_history->bForce_correlation)
    {
        awh_bias_history->forcecorr_hist = init_correlation_grid_history_from_state(awh_bias->forcecorr);
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
        init_awh_bias_history_from_state(&awh_history->awh_bias_history[k], &awh->awh_bias[k]);
    }
}

static void update_awh_bias_history(awh_bias_history_t *awh_bias_history, awh_bias_t *awh_bias)
{
    GMX_RELEASE_ASSERT(awh_bias_history->coordpoint != NULL, "AWH history coord point not initialized when updating history");
    GMX_RELEASE_ASSERT(awh_bias_history->npoints > 0, "AWH history number of points not initialized when updating history");
    GMX_RELEASE_ASSERT(awh_bias_history->ndim > 0, "AWH history number of dimensions not initialized when updating history");

    awh_bias_history->coord_refvalue_index = awh_bias->coord_refvalue_index;

    for (int m = 0; m < awh_bias_history->npoints; m++)
    {
        awh_bias_history->coordpoint[m].target                               = awh_bias->coordpoint[m].target;
        awh_bias_history->coordpoint[m].free_energy                          = awh_bias->coordpoint[m].free_energy;
        awh_bias_history->coordpoint[m].bias                                 = awh_bias->coordpoint[m].bias;
        awh_bias_history->coordpoint[m].weightsum_iteration                  = awh_bias->coordpoint[m].weightsum_iteration;
        awh_bias_history->coordpoint[m].weightsum_covering                   = awh_bias->coordpoint[m].weightsum_covering;
        awh_bias_history->coordpoint[m].weightsum_tot                        = awh_bias->coordpoint[m].weightsum_tot;
        awh_bias_history->coordpoint[m].weightsum_ref                        = awh_bias->coordpoint[m].weightsum_ref;
        awh_bias_history->coordpoint[m].last_update_index                    = awh_bias->coordpoint[m].last_update_index;
        awh_bias_history->coordpoint[m].log_pmfsum                           = awh_bias->coordpoint[m].log_pmfsum;
        awh_bias_history->coordpoint[m].visits_iteration                     = awh_bias->coordpoint[m].visits_iteration;
        awh_bias_history->coordpoint[m].visits_tot                           = awh_bias->coordpoint[m].visits_tot;
    }

    awh_bias_history->in_initial             = awh_bias->in_initial;
    awh_bias_history->histsize               = awh_bias->histsize;

    awh_bias_history->origin_index_updatelist = multidim_gridindex_to_linear(awh_bias->grid,
                                                                             awh_bias->origin_updatelist);
    awh_bias_history->end_index_updatelist = multidim_gridindex_to_linear(awh_bias->grid,
                                                                          awh_bias->end_updatelist);

    if (awh_bias->bForce_correlation)
    {
        GMX_RELEASE_ASSERT(awh_bias_history->forcecorr_hist != NULL, "AWH history force correlation not initialized when updating history");
        update_correlation_grid_history(awh_bias_history->forcecorr_hist, awh_bias->forcecorr);
    }

    awh_bias_history->log_relative_sampleweight  = awh_bias->log_relative_sampleweight;
}

static void restore_awh_state_from_history(const awh_bias_history_t *awh_bias_history, awh_bias_t *awh_bias)
{
    awh_bias->coord_refvalue_index = awh_bias_history->coord_refvalue_index;

    for (int m = 0; m < awh_bias->npoints; m++)
    {
        awh_bias->coordpoint[m].target                                = awh_bias_history->coordpoint[m].target;
        awh_bias->coordpoint[m].free_energy                           = awh_bias_history->coordpoint[m].free_energy;
        awh_bias->coordpoint[m].bias                                  = awh_bias_history->coordpoint[m].bias;
        awh_bias->coordpoint[m].weightsum_iteration                   = awh_bias_history->coordpoint[m].weightsum_iteration;
        awh_bias->coordpoint[m].weightsum_covering                    = awh_bias_history->coordpoint[m].weightsum_covering;
        awh_bias->coordpoint[m].weightsum_tot                         = awh_bias_history->coordpoint[m].weightsum_tot;
        awh_bias->coordpoint[m].weightsum_ref                         = awh_bias_history->coordpoint[m].weightsum_ref;
        awh_bias->coordpoint[m].last_update_index                     = awh_bias_history->coordpoint[m].last_update_index;
        awh_bias->coordpoint[m].log_pmfsum                            = awh_bias_history->coordpoint[m].log_pmfsum;
        awh_bias->coordpoint[m].visits_iteration                      = awh_bias_history->coordpoint[m].visits_iteration;
        awh_bias->coordpoint[m].visits_tot                            = awh_bias_history->coordpoint[m].visits_tot;
    }

    awh_bias->in_initial      = awh_bias_history->in_initial;
    awh_bias->histsize        = awh_bias_history->histsize;

    linear_gridindex_to_multidim(awh_bias->grid, awh_bias_history->origin_index_updatelist, awh_bias->origin_updatelist);
    linear_gridindex_to_multidim(awh_bias->grid, awh_bias_history->end_index_updatelist, awh_bias->end_updatelist);

    if (awh_bias->bForce_correlation)
    {
        restore_correlation_grid_state_from_history(awh_bias_history->forcecorr_hist, awh_bias->forcecorr);
    }

    awh_bias->log_relative_sampleweight  = awh_bias_history->log_relative_sampleweight;
}

void update_awh_history(awh_history_t *awh_history, const awh_t *awh)
{
    GMX_RELEASE_ASSERT(awh_history->awh_bias_history != NULL, "AWH history not initialized when updating history");

    awh_history->convolved_bias_shift = awh->convolved_bias_shift;

    for (int k = 0; k < awh_history->nbias; k++)
    {
        update_awh_bias_history(&awh_history->awh_bias_history[k], &awh->awh_bias[k]);
    }
}

static void broadcast_initialize_awh_bias_history(awh_bias_history_t *awh_bias_history, const t_commrec *cr)
{
    gmx_bcast(sizeof(awh_bias_history_t), awh_bias_history, cr);

    if (!MASTER(cr))
    {
        snew(awh_bias_history->coordpoint, awh_bias_history->npoints);
    }
    gmx_bcast(sizeof(awh_coordpoint_history_t)*awh_bias_history->npoints, awh_bias_history->coordpoint, cr);

    if (awh_bias_history->bForce_correlation)
    {
        correlation_grid_history_t *corrgrid_hist = init_correlation_grid_history_from_checkpoint(awh_bias_history->forcecorr_hist, cr);

        if (!MASTER(cr))
        {
            awh_bias_history->forcecorr_hist = corrgrid_hist;
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
        broadcast_initialize_awh_bias_history(&awh_history->awh_bias_history[k], cr);
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
        restore_awh_state_from_history(&awh_history->awh_bias_history[k], &awh->awh_bias[k]);
    }
}
