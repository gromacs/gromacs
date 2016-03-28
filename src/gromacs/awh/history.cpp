/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
#include "pointstate.h"

/* Allocate and initialize an AWH history with the given AWH state. */
void AwhBiasCollection::initHistoryFromState(AwhHistory *awhHistory) const
{
    GMX_RELEASE_ASSERT(awhHistory != nullptr, "Should be called with a valid history struct (and only on the master rank)");

    awhHistory->bias.clear();
    awhHistory->bias.resize(bias_.size());

    for (size_t k = 0; k < awhHistory->bias.size(); k++)
    {
        AwhBiasHistory *biasHistory = &awhHistory->bias[k];
        biasHistory->pointState.resize(bias_[k].pointState.size());
    }
}

/*! \brief
 * Update the bias history with a new state.
 *
 * \param[in,out] biasHistory  Bias history struct.
 * \param[in] bias             Bias state.
 */
static void updateBiasHistory(AwhBiasHistory *biasHistory, const Bias &bias)
{
    GMX_RELEASE_ASSERT(biasHistory->pointState.size() == bias.pointState.size(), "The AWH history setup does not match the AWH state.");

    const BiasState     &state            = bias.state;
    AwhBiasStateHistory *stateHistory     = &biasHistory->state;
    stateHistory->refGridpoint            = state.refGridpoint;

    for (size_t m = 0; m < biasHistory->pointState.size(); m++)
    {
        const PointState     &pointState  = bias.pointState[m];
        AwhPointStateHistory *psh         = &biasHistory->pointState[m];

        psh->target                       = pointState.target;
        psh->free_energy                  = pointState.freeEnergy();
        psh->bias                         = pointState.bias();
        psh->weightsum_iteration          = pointState.weightsum_iteration;
        psh->weightsum_covering           = pointState.weightsum_covering;
        psh->weightsum_tot                = pointState.weightsumTot;
        psh->weightsum_ref                = pointState.weightsumRef();
        psh->last_update_index            = pointState.lastUpdateIndex();
        psh->log_pmfsum                   = pointState.logPmfsum();
        psh->visits_iteration             = pointState.visits_iteration;
        psh->visits_tot                   = pointState.visits_tot;
    }

    stateHistory->in_initial              = state.in_initial;
    stateHistory->equilibrateHistogram    = state.equilibrateHistogram;
    stateHistory->histSize                = state.histSize;

    stateHistory->origin_index_updatelist = multidim_gridindex_to_linear(bias.grid.get(),
                                                                         state.origin_updatelist);
    stateHistory->end_index_updatelist    = multidim_gridindex_to_linear(bias.grid.get(),
                                                                         state.end_updatelist);

    stateHistory->scaledSampleWeight      = state.scaledSampleWeight;
    stateHistory->maxScaledSampleWeight   = state.maxScaledSampleWeight;
}

/*! \brief
 * Restore the bias state from history.
 *
 * \param[in] biasHistory  Bias history struct.
 * \param[in,out] bias     Bias state.
 */
static void restoreBiasStateFromHistory(const AwhBiasHistory *biasHistory, Bias *bias)
{
    BiasState                 *state        = &bias->state;
    const AwhBiasStateHistory &stateHistory = biasHistory->state;

    state->refGridpoint = stateHistory.refGridpoint;

    for (size_t m = 0; m < bias->pointState.size(); m++)
    {
        bias->pointState[m].setFromHistory(biasHistory->pointState[m]);
    }

    state->in_initial             = stateHistory.in_initial;
    state->equilibrateHistogram   = stateHistory.equilibrateHistogram;
    state->histSize               = stateHistory.histSize;

    linear_gridindex_to_multidim(bias->grid.get(), stateHistory.origin_index_updatelist, state->origin_updatelist);
    linear_gridindex_to_multidim(bias->grid.get(), stateHistory.end_index_updatelist, state->end_updatelist);

    state->scaledSampleWeight     = stateHistory.scaledSampleWeight;
    state->maxScaledSampleWeight  = stateHistory.maxScaledSampleWeight;
}

/*! \brief
 * Broadcast the bias data over the MPI ranks.
 *
 * \param[in,out] bias  The AWH Bias, with complete data on the master rank.
 * \param[in] cr        Struct for communication.
 */
static void broadcastBias(Bias *bias, const t_commrec *cr)
{
    gmx_bcast(bias->pointState.size()*sizeof(PointState), bias->pointState.data(), cr);

    gmx_bcast(sizeof(BiasState), &bias->state, cr);
}

/* Restores the AWH state from history. */
void AwhBiasCollection::restoreStateFromHistory(const AwhHistory &awhHistory,
                                                const t_commrec  *cr)
{
    /* Restore the history to the current state */
    if (MASTER(cr))
    {
        GMX_RELEASE_ASSERT(awhHistory.bias.size() == bias_.size(), "AWH state and history bias count should match");

        potentialOffset_ = awhHistory.potentialOffset;
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(potentialOffset_), &potentialOffset_, cr);
    }

    for (size_t k = 0; k < bias_.size(); k++)
    {
        if (MASTER(cr))
        {
            restoreBiasStateFromHistory(&awhHistory.bias[k], &bias_[k]);
        }

        if (PAR(cr))
        {
            broadcastBias(&bias_[k], cr);
        }
    }
}

/* Update the AWH bias history for checkpointing. */
void AwhBiasCollection::updateHistory(AwhHistory *awhHistory) const
{
    /* This assert will also catch a non-master rank calling this function. */
    GMX_RELEASE_ASSERT(awhHistory->bias.size() == bias_.size(), "AWH state and history bias count should match");

    awhHistory->potentialOffset = potentialOffset_;

    for (size_t k = 0; k < awhHistory->bias.size(); k++)
    {
        updateBiasHistory(&awhHistory->bias[k], bias_[k]);
    }
}
