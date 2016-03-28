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

#include "awhhistory.h"

#include <algorithm>
#include <assert.h>
#include <cstring>

#include "gromacs/awh/awh.h"
#include "gromacs/awh/awh-correlation-history.h"
#include "gromacs/awh/awh-types.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/fatalerror.h"

/* These are just dummy initializers used when initializing the state because all AWH parameters
   are not (readily) available then/there (if the awhbias struct hasn't been initialized). TODO. */
static void init_awhhistory(awhhistory_t *awhhist)
{
    awhhist->in_initial                = 0;
    awhhist->histsize                  = 0;
    awhhist->npoints                   = 0;
    awhhist->coord_point               = NULL;
    awhhist->coord_refvalue_index      = 0;
    awhhist->ndim                      = 0;
    awhhist->bForce_correlation        = FALSE;
    awhhist->forcecorr_hist            = NULL;
    awhhist->log_relative_sampleweight = 0;
}

static void init_awhhistory_from_state(awhhistory_t *awhhist, const t_awh *awh)
{
    init_awhhistory(awhhist);

    awhhist->npoints            = awh->npoints;
    awhhist->ndim               = awh->ndim;
    awhhist->bForce_correlation = awh->bForce_correlation;

    snew(awhhist->coord_point, awhhist->npoints);

    if (awhhist->bForce_correlation)
    {
        snew(awhhist->forcecorr_hist, 1);
        init_correlation_grid_history_from_state(awhhist->forcecorr_hist, awh->forcecorr);
    }
}

void init_awhbiashistory(awhbiashistory_t *awhbiashist)
{
    awhbiashist->used                      = FALSE;
    awhbiashist->convolved_bias_shift      = 0;
    awhbiashist->nawhhist                  = 0;
    awhbiashist->awhhist                   = NULL;
}

void init_awhbiashistory_from_state(awhbiashistory_t *awhbiashist, const t_awhbias *awhbias)
{
    awhbiashist->used                      = TRUE;
    awhbiashist->convolved_bias_shift      = awhbias->convolved_bias_shift;

    awhbiashist->nawhhist = awhbias->nawh;
    snew(awhbiashist->awhhist, awhbiashist->nawhhist);
    for (int k = 0; k < awhbiashist->nawhhist; k++)
    {
        init_awhhistory_from_state(&awhbiashist->awhhist[k], &awhbias->awh[k]);
    }
}

static void update_awhhistory(awhhistory_t *awhhist, t_awh *awh)
{
    GMX_RELEASE_ASSERT(awhhist->coord_point != NULL, "AWH history coord point not initialized when updating history");
    GMX_RELEASE_ASSERT(awhhist->npoints > 0, "AWH history number of points not initialized when updating history");
    GMX_RELEASE_ASSERT(awhhist->ndim > 0, "AWH history number of dimensions not initialized when updating history");

    awhhist->coord_refvalue_index = awh->coord_refvalue_index;

    for (int m = 0; m < awhhist->npoints; m++)
    {
        awhhist->coord_point[m].target                               = awh->coord_point[m].target;
        awhhist->coord_point[m].free_energy                          = awh->coord_point[m].free_energy;
        awhhist->coord_point[m].bias                                 = awh->coord_point[m].bias;
        awhhist->coord_point[m].weightsum_iteration                  = awh->coord_point[m].weightsum_iteration;
        awhhist->coord_point[m].weightsum_covering                   = awh->coord_point[m].weightsum_covering;
        awhhist->coord_point[m].weightsum_tot                        = awh->coord_point[m].weightsum_tot;
        awhhist->coord_point[m].weightsum_ref                        = awh->coord_point[m].weightsum_ref;
        awhhist->coord_point[m].last_update_index                    = awh->coord_point[m].last_update_index;
        awhhist->coord_point[m].log_pmfsum                           = awh->coord_point[m].log_pmfsum;
        awhhist->coord_point[m].visits_iteration                     = awh->coord_point[m].visits_iteration;
        awhhist->coord_point[m].visits_tot                           = awh->coord_point[m].visits_tot;
    }

    awhhist->in_initial             = awh->in_initial;
    awhhist->histsize               = awh->histsize;

    for (int d = 0; d < awhhist->ndim; d++)
    {
        awhhist->origin_updatelist[d] = awh->origin_updatelist[d];
        awhhist->end_updatelist[d]    = awh->end_updatelist[d];
    }

    if (awh->bForce_correlation)
    {
        GMX_RELEASE_ASSERT(awhhist->forcecorr_hist != NULL, "AWH history force correlation not initialized when updating history");
        update_correlation_grid_history(awhhist->forcecorr_hist, awh->forcecorr);
    }

    awhhist->log_relative_sampleweight  = awh->log_relative_sampleweight;
}

static void restore_awh_state_from_history(const awhhistory_t *awhhist, t_awh *awh)
{
    awh->coord_refvalue_index = awhhist->coord_refvalue_index;

    for (int m = 0; m < awh->npoints; m++)
    {
        awh->coord_point[m].target                                = awhhist->coord_point[m].target;
        awh->coord_point[m].free_energy                           = awhhist->coord_point[m].free_energy;
        awh->coord_point[m].bias                                  = awhhist->coord_point[m].bias;
        awh->coord_point[m].weightsum_iteration                   = awhhist->coord_point[m].weightsum_iteration;
        awh->coord_point[m].weightsum_covering                    = awhhist->coord_point[m].weightsum_covering;
        awh->coord_point[m].weightsum_tot                         = awhhist->coord_point[m].weightsum_tot;
        awh->coord_point[m].weightsum_ref                         = awhhist->coord_point[m].weightsum_ref;
        awh->coord_point[m].last_update_index                     = awhhist->coord_point[m].last_update_index;
        awh->coord_point[m].log_pmfsum                            = awhhist->coord_point[m].log_pmfsum;
        awh->coord_point[m].visits_iteration                      = awhhist->coord_point[m].visits_iteration;
        awh->coord_point[m].visits_tot                            = awhhist->coord_point[m].visits_tot;
    }

    awh->in_initial      = awhhist->in_initial;
    awh->histsize        = awhhist->histsize;

    for (int d = 0; d < awh->ndim; d++)
    {
        awh->origin_updatelist[d] = awhhist->origin_updatelist[d];
        awh->end_updatelist[d]    = awhhist->end_updatelist[d];
    }

    if (awh->bForce_correlation)
    {
        restore_correlation_grid_state_from_history(awhhist->forcecorr_hist, awh->forcecorr);
    }

    awh->log_relative_sampleweight  = awhhist->log_relative_sampleweight;
}

void update_awhbiashistory(awhbiashistory_t *awhbiashist, const t_awhbias *awhbias)
{
    GMX_RELEASE_ASSERT(awhbiashist->awhhist != NULL, "AWH history not initialized when updating history");

    awhbiashist->convolved_bias_shift = awhbias->convolved_bias_shift;

    for (int k = 0; k < awhbiashist->nawhhist; k++)
    {
        update_awhhistory(&awhbiashist->awhhist[k], &awhbias->awh[k]);
    }
}

static void broadcast_initialize_awhhistory(awhhistory_t *awhhist, const t_commrec *cr)
{
    gmx_bcast(sizeof(awhhistory_t), awhhist, cr);

    if (!MASTER(cr))
    {
        snew(awhhist->coord_point, awhhist->npoints);
    }
    gmx_bcast(sizeof(awhhistory_coord_point_t)*awhhist->npoints, awhhist->coord_point, cr);

    if (awhhist->bForce_correlation)
    {
        if (!MASTER(cr))
        {
            snew(awhhist->forcecorr_hist, 1);
        }
        init_correlation_grid_history_from_checkpoint(awhhist->forcecorr_hist, cr);
    }
}

static void broadcast_initialize_awhbiashistory(awhbiashistory_t *awhbiashist, const t_commrec *cr)
{
    gmx_bcast(sizeof(awhbiashistory_t), awhbiashist, cr);
    GMX_RELEASE_ASSERT(awhbiashist->nawhhist > 0, "nawhhist needs to be > 0");

    if (!MASTER(cr))
    {
        snew(awhbiashist->awhhist, awhbiashist->nawhhist);
    }

    for (int k = 0; k < awhbiashist->nawhhist; k++)
    {
        broadcast_initialize_awhhistory(&awhbiashist->awhhist[k], cr);
    }
}

void init_awhbiashistory_from_checkpoint(awhbiashistory_t *awhbiashist, const t_commrec *cr)
{
    /* Only master has read the checkpoint and only master has an initialized AWH history struct.
       Non-master ranks need to get the their history allocated and initialized by broadcasting from the master rank. */
    if (PAR(cr))
    {
        GMX_RELEASE_ASSERT((MASTER(cr) && awhbiashist->awhhist != NULL) || (!MASTER(cr) && awhbiashist->awhhist == NULL),
                           "(Non-)master rank(s) should (not) have been initialized before restoring awhbiashistory");
        broadcast_initialize_awhbiashistory(awhbiashist, cr);
    }
}

void restore_awhbias_state_from_history(const awhbiashistory_t *awhbiashist, t_awhbias *awhbias)
{
    GMX_RELEASE_ASSERT(awhbiashist->awhhist != NULL, "AWH history was not initialized when attempting to restore state from history");
    GMX_RELEASE_ASSERT(awhbias != NULL, "AWH bias not initialized when attempting to restore state from history");

    /* Restore the history to the current state */
    awhbias->convolved_bias_shift       = awhbiashist->convolved_bias_shift;

    for (int k = 0; k < awhbiashist->nawhhist; k++)
    {
        restore_awh_state_from_history(&awhbiashist->awhhist[k], &awhbias->awh[k]);
    }
}

bool force_correlation_needs_checkpointing(const awhbiashistory_t *awhbiashist)
{
    /* Currently each AWH history has its own variable although they all have the same value */
    return awhbiashist->nawhhist > 0 &&  awhbiashist->awhhist[0].bForce_correlation;
}
