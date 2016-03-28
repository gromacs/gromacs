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

#include "grid.h"
#include "types.h"

/*! \brief
 * Trivially initialize bias history.
 *
 * These are just dummy initializers used when initializing the state because all AWH parameters
 * are not (readily) available then/there (if the awh struct hasn't been initialized).
 *
 * \param[in,out] bias_history   Bias history struct.
 */
static void init_bias_history(awh_bias_history_t *bias_history)
{
    bias_history->in_initial                = 0;
    bias_history->equilibrateHistogram      = 0;
    bias_history->histsize                  = 0;
    bias_history->npoints                   = 0;
    bias_history->coordpoint                = NULL;
    bias_history->refCoordpoint             = 0;
    bias_history->ndim                      = 0;
    bias_history->scaledSampleWeight        = 0;
    bias_history->maxScaledSampleWeight     = 0;
}

/*! \brief
 * Initialize and allocated for bias history using parameters from bias state.
 *
 * \param[in,out] bias_history   Bias history struct.
 * \param[in] bias               Bias state.
 */
static void init_bias_history_from_state(awh_bias_history_t *bias_history, const awh_bias_t *bias)
{
    init_bias_history(bias_history);

    bias_history->npoints            = bias->npoints;
    bias_history->ndim               = bias->ndim;

    snew(bias_history->coordpoint, bias_history->npoints);
}

/* Allocate and initialize an AWH history with the given AWH state. */
void init_awh_history_from_state(awh_history_t *awh_history, const awh_t *awh)
{
    init_awh_history(awh_history);

    awh_history->nbias = awh->nbias;
    snew(awh_history->bias, awh_history->nbias);
    for (int k = 0; k < awh_history->nbias; k++)
    {
        init_bias_history_from_state(&awh_history->bias[k], &awh->awh_bias[k]);
    }
}

/*! \brief
 * Update the bias history with a new state.
 *
 * \param[in,out] bias_history   Bias history struct.
 * \param[in] bias               Bias state.
 */
static void update_bias_history(awh_bias_history_t *bias_history, const awh_bias_t *bias)
{
    GMX_RELEASE_ASSERT(bias_history->coordpoint != NULL, "AWH history coord point not initialized when updating history");
    GMX_RELEASE_ASSERT(bias_history->npoints > 0, "AWH history number of points not initialized when updating history");
    GMX_RELEASE_ASSERT(bias_history->ndim > 0, "AWH history number of dimensions not initialized when updating history");

    bias_history->refCoordpoint = bias->refCoordpoint;

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
    bias_history->equilibrateHistogram   = bias->equilibrateHistogram;
    bias_history->histsize               = bias->histsize;

    bias_history->origin_index_updatelist = multidim_gridindex_to_linear(bias->grid,
                                                                         bias->origin_updatelist);
    bias_history->end_index_updatelist = multidim_gridindex_to_linear(bias->grid,
                                                                      bias->end_updatelist);

    bias_history->scaledSampleWeight     = bias->scaledSampleWeight;
    bias_history->maxScaledSampleWeight  = bias->maxScaledSampleWeight;
}

/*! \brief
 * Restore the bias state from history.
 *
 * \param[in] bias_history   Bias history struct.
 * \param[in,out] bias       Bias state.
 */
static void restore_bias_state_from_history(const awh_bias_history_t *bias_history, awh_bias_t *bias)
{
    bias->refCoordpoint = bias_history->refCoordpoint;

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

    bias->in_initial             = bias_history->in_initial;
    bias->equilibrateHistogram   = bias_history->equilibrateHistogram;
    bias->histsize               = bias_history->histsize;

    linear_gridindex_to_multidim(bias->grid, bias_history->origin_index_updatelist, bias->origin_updatelist);
    linear_gridindex_to_multidim(bias->grid, bias_history->end_index_updatelist, bias->end_updatelist);

    bias->scaledSampleWeight     = bias_history->scaledSampleWeight;
    bias->maxScaledSampleWeight  = bias_history->maxScaledSampleWeight;
}

/* Update the AWH bias history for checkpointing. */
void update_awh_history(awh_history_t *awh_history, const awh_t *awh)
{
    GMX_RELEASE_ASSERT(awh_history->bias != NULL, "AWH history not initialized when updating history");

    awh_history->potential_offset = awh->potential_offset;

    for (int k = 0; k < awh_history->nbias; k++)
    {
        update_bias_history(&awh_history->bias[k], &awh->awh_bias[k]);
    }
}

/*! \brief
 * Initialize and allocate for the bias history when only master rank has been initialized.
 *
 * \param[in,out] bias_history   Bias history struct.
 * \param[in] cr                 Struct for communication.
 */
static void broadcast_initialize_bias_history(awh_bias_history_t *bias_history, const t_commrec *cr)
{
    gmx_bcast(sizeof(awh_bias_history_t), bias_history, cr);

    if (!MASTER(cr))
    {
        snew(bias_history->coordpoint, bias_history->npoints);
    }
    gmx_bcast(sizeof(awh_coordpoint_history_t)*bias_history->npoints, bias_history->coordpoint, cr);
}

/*! \brief
 * Initialize and allocate for the AWH history when only master rank has been initialized.
 *
 * \param[in,out] awh_history   AWH history struct.
 * \param[in] cr                Struct for communication.
 */
static void broadcast_initialize_awh_history(awh_history_t *awh_history, const t_commrec *cr)
{
    gmx_bcast(sizeof(awh_history_t), awh_history, cr);
    GMX_RELEASE_ASSERT(awh_history->nbias > 0, "nbias needs to be > 0");

    if (!MASTER(cr))
    {
        snew(awh_history->bias, awh_history->nbias);
    }

    for (int k = 0; k < awh_history->nbias; k++)
    {
        broadcast_initialize_bias_history(&awh_history->bias[k], cr);
    }
}

/* Initialize an AWH history when restarting a from checkpoint. */
void init_awh_history_from_checkpoint(awh_history_t *awh_history, const t_commrec *cr)
{
    /* Only master has read the checkpoint and only master has an initialized AWH history struct.
       Non-master ranks need to get the their history allocated and initialized by broadcasting from the master rank. */
    if (PAR(cr))
    {
        GMX_RELEASE_ASSERT((MASTER(cr) && awh_history->bias != NULL) || (!MASTER(cr) && awh_history->bias == NULL),
                           "(Non-)master rank(s) should (not) have been initialized before restoring awh_history");
        broadcast_initialize_awh_history(awh_history, cr);
    }
}

/* Restores the AWH state from history. */
void restore_awh_state_from_history(const awh_history_t *awh_history, awh_t *awh)
{
    GMX_RELEASE_ASSERT(awh_history->bias != NULL, "AWH history was not initialized when attempting to restore state from history");
    GMX_RELEASE_ASSERT(awh != NULL, "AWH bias not initialized when attempting to restore state from history");

    /* Restore the history to the current state */
    awh->potential_offset = awh_history->potential_offset;

    for (int k = 0; k < awh_history->nbias; k++)
    {
        restore_bias_state_from_history(&awh_history->bias[k], &awh->awh_bias[k]);
    }
}
