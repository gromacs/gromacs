/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \internal \file
 * \brief Defines functions of the AWH history used by modular checkpointing
 *
 * \author Pascal Merz <pascal.merz@me.com>
 */

#include "gmxpre.h"

#include "awh_history.h"

#include <string>

#include "gromacs/utility/strconvert.h"

#include "checkpointdata.h"

namespace gmx
{

namespace
{
/*!
 * \brief Enum describing the contents all AWH history classes implemented in this file
 *        write to modular checkpoint
 *
 * When changing the checkpoint content, add a new element just above Count, and adjust the
 * checkpoint functionality.
 */
enum class CheckpointVersion
{
    Base, //!< First version of modular checkpointing
    Count //!< Number of entries. Add new versions right above this!
};
constexpr auto c_currentVersion = CheckpointVersion(int(CheckpointVersion::Count) - 1);

template<CheckpointDataOperation operation>
void doCheckpoint(CheckpointData<operation> checkpointData, AwhPointStateHistory* awhPointStateHistory)
{
    checkpointData.scalar("bias", &awhPointStateHistory->bias);
    checkpointData.scalar("free_energy", &awhPointStateHistory->free_energy);
    checkpointData.scalar("target", &awhPointStateHistory->target);
    checkpointData.scalar("weightsum_iteration", &awhPointStateHistory->weightsum_iteration);
    checkpointData.scalar("weightsum_covering", &awhPointStateHistory->weightsum_covering);
    checkpointData.scalar("weightsum_tot", &awhPointStateHistory->weightsum_tot);
    checkpointData.scalar("weightsum_ref", &awhPointStateHistory->weightsum_ref);
    checkpointData.scalar("last_update_index", &awhPointStateHistory->last_update_index);
    checkpointData.scalar("log_pmfsum", &awhPointStateHistory->log_pmfsum);
    checkpointData.scalar("visits_iteration", &awhPointStateHistory->visits_iteration);
    checkpointData.scalar("visits_tot", &awhPointStateHistory->visits_tot);
    checkpointData.scalar("localWeightSum", &awhPointStateHistory->localWeightSum);
}

template<CheckpointDataOperation operation>
void doCheckpoint(CheckpointData<operation> checkpointData, AwhBiasStateHistory* awhBiasStateHistory)
{
    checkpointData.scalar("umbrellaGridpoint", &awhBiasStateHistory->umbrellaGridpoint);
    checkpointData.scalar("origin_index_updatelist", &awhBiasStateHistory->origin_index_updatelist);
    checkpointData.scalar("end_index_updatelist", &awhBiasStateHistory->end_index_updatelist);
    checkpointData.scalar("in_initial", &awhBiasStateHistory->in_initial);
    checkpointData.scalar("equilibrateHistogram", &awhBiasStateHistory->equilibrateHistogram);
    checkpointData.scalar("histSize", &awhBiasStateHistory->histSize);
    checkpointData.scalar("logScaledSampleWeight", &awhBiasStateHistory->logScaledSampleWeight);
    checkpointData.scalar("maxLogScaledSampleWeight", &awhBiasStateHistory->maxLogScaledSampleWeight);
    checkpointData.scalar("numUpdates", &awhBiasStateHistory->numUpdates);
}

template<CheckpointDataOperation operation>
void doCheckpoint(CheckpointData<operation>    checkpointData,
                  CorrelationBlockDataHistory* correlationBlockDataHistory)
{
    checkpointData.scalar("blockSumWeight", &correlationBlockDataHistory->blockSumWeight);
    checkpointData.scalar("blockSumSquareWeight", &correlationBlockDataHistory->blockSumSquareWeight);
    checkpointData.scalar("blockSumWeightX", &correlationBlockDataHistory->blockSumWeightX);
    checkpointData.scalar("blockSumWeightY", &correlationBlockDataHistory->blockSumWeightY);
    checkpointData.scalar("sumOverBlocksSquareBlockWeight",
                          &correlationBlockDataHistory->sumOverBlocksSquareBlockWeight);
    checkpointData.scalar("sumOverBlocksBlockSquareWeight",
                          &correlationBlockDataHistory->sumOverBlocksBlockSquareWeight);
    checkpointData.scalar("sumOverBlocksBlockWeightBlockWeightX",
                          &correlationBlockDataHistory->sumOverBlocksBlockWeightBlockWeightX);
    checkpointData.scalar("sumOverBlocksBlockWeightBlockWeightY",
                          &correlationBlockDataHistory->sumOverBlocksBlockWeightBlockWeightY);
    checkpointData.scalar("blockLength", &correlationBlockDataHistory->blockLength);
    checkpointData.scalar("previousBlockIndex", &correlationBlockDataHistory->previousBlockIndex);
    checkpointData.scalar("correlationIntegral", &correlationBlockDataHistory->correlationIntegral);
}

template<CheckpointDataOperation operation>
void doCheckpoint(CheckpointData<operation> checkpointData, CorrelationGridHistory* correlationGridHistory)
{
    checkpointData.scalar("numCorrelationTensors", &correlationGridHistory->numCorrelationTensors);
    checkpointData.scalar("tensorSize", &correlationGridHistory->tensorSize);
    checkpointData.scalar("blockDataListSize", &correlationGridHistory->blockDataListSize);

    int blockDataBufferSize = correlationGridHistory->blockDataBuffer.size();
    checkpointData.scalar("blockDataBufferSize", &blockDataBufferSize);
    if (operation == CheckpointDataOperation::Read)
    {
        correlationGridHistory->blockDataBuffer.resize(blockDataBufferSize);
    }
    int counter = 0;
    for (auto& buffer : correlationGridHistory->blockDataBuffer)
    {
        doCheckpoint(checkpointData.subCheckpointData("blockDataBuffer " + toString(counter)), &buffer);
        counter++;
    }
}

template<CheckpointDataOperation operation>
void doCheckpointData(CheckpointData<operation> checkpointData, AwhBiasHistory* awhBiasHistory)
{
    int pointStateSize = awhBiasHistory->pointState.size();
    checkpointData.scalar("pointStateSize", &pointStateSize);
    if (operation == CheckpointDataOperation::Read)
    {
        awhBiasHistory->pointState.resize(pointStateSize);
    }

    int counter = 0;
    for (auto& point : awhBiasHistory->pointState)
    {
        doCheckpoint(checkpointData.subCheckpointData("pointState " + toString(counter)), &point);
        counter++;
    }
    doCheckpoint(checkpointData.subCheckpointData("state"), &awhBiasHistory->state);
    doCheckpoint(checkpointData.subCheckpointData("forceCorrelationGrid"),
                 &awhBiasHistory->forceCorrelationGrid);
}
} // namespace

// This trips doxygen up
//! \cond
template<CheckpointDataOperation operation>
void AwhHistory::doCheckpoint(CheckpointData<operation> checkpointData)
{
    // Keep one version for all objects in this file
    checkpointVersion(&checkpointData, "AwhHistory version", c_currentVersion);

    checkpointData.scalar("potentialOffset", &potentialOffset);

    int biasSize = bias.size();
    checkpointData.scalar("biasSize", &biasSize);
    if (operation == CheckpointDataOperation::Read)
    {
        bias.resize(biasSize);
    }
    int counter = 0;
    for (auto& b : bias)
    {
        doCheckpointData(checkpointData.subCheckpointData("bias " + toString(counter)), &b);
        counter++;
    }
}

// explicit template instantiation
template void AwhHistory::doCheckpoint(CheckpointData<CheckpointDataOperation::Read>);
template void AwhHistory::doCheckpoint(CheckpointData<CheckpointDataOperation::Write>);
//! \endcond

} // namespace gmx
