/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * \brief
 * Implements bias sharing checking functionality.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "biassharing.h"

#include "config.h"

#include <algorithm>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

//! Determines and returns which of the local biases are shared with who how many other simulations
std::multiset<int> getGlobalShareIndices(ArrayRef<const int> localShareIndices, MPI_Comm simulationMainComm)
{
#if GMX_MPI
    int numSimulations;
    MPI_Comm_size(simulationMainComm, &numSimulations);
    int ourRank;
    MPI_Comm_rank(simulationMainComm, &ourRank);
    std::vector<int> biasCountsIn(numSimulations, 0);
    std::vector<int> biasCounts(numSimulations, 0);
    biasCountsIn[ourRank] = localShareIndices.size();
    MPI_Allreduce(biasCountsIn.data(), biasCounts.data(), numSimulations, MPI_INT, MPI_SUM, simulationMainComm);
    // Now we need to gather the share indices to all (main) ranks.
    // We could use MPI_Allgatherv, but thread-MPI does not support that and using
    // MPI_Allreduce produces simpler code, so we use that.
    int totNumBiases = 0;
    int ourOffset    = 0;
    for (int rank = 0; rank < numSimulations; rank++)
    {
        if (rank == ourRank)
        {
            ourOffset = totNumBiases;
        }
        totNumBiases += biasCounts[rank];
    }
    // Fill a buffer with zeros and our part of sharing indices
    std::vector<int> shareIndicesAllIn(totNumBiases, 0);
    std::copy(localShareIndices.begin(), localShareIndices.end(), shareIndicesAllIn.begin() + ourOffset);
    // Gather all sharing indices to all (main) ranks
    std::vector<int> shareIndicesAll(totNumBiases);
    MPI_Allreduce(shareIndicesAllIn.data(), shareIndicesAll.data(), totNumBiases, MPI_INT, MPI_SUM, simulationMainComm);
#else
    GMX_UNUSED_VALUE(simulationMainComm);

    ArrayRef<const int> shareIndicesAll = localShareIndices;
#endif // GMX_MPI

    std::multiset<int> shareIndicesSet;
    for (int shareIndex : shareIndicesAll)
    {
        if (shareIndex > 0)
        {
            shareIndicesSet.insert(shareIndex);
        }
    }

    return shareIndicesSet;
}

} // namespace

BiasSharing::BiasSharing(const AwhParams& awhParams, const t_commrec& commRecord, MPI_Comm simulationMainComm) :
    commRecord_(commRecord)
{
    if (MAIN(&commRecord))
    {
        std::vector<int> localShareIndices;
        int              shareGroupPrev = 0;
        for (int k = 0; k < awhParams.numBias(); k++)
        {
            const int shareGroup = awhParams.awhBiasParams(k).shareGroup();
            GMX_RELEASE_ASSERT(shareGroup >= 0, "Bias share group values should be >= 0");
            localShareIndices.push_back(shareGroup);
            if (shareGroup > 0)
            {
                if (shareGroup <= shareGroupPrev)
                {
                    GMX_THROW(
                            InvalidInputError("AWH biases that are shared should use increasing "
                                              "share-group values"));
                }
                shareGroupPrev = shareGroup;
            }
        }
        std::multiset<int> globalShareIndices =
                getGlobalShareIndices(localShareIndices, simulationMainComm);

        int numSimulations = 1;
#if GMX_MPI
        MPI_Comm_size(simulationMainComm, &numSimulations);
        int myRank;
        MPI_Comm_rank(simulationMainComm, &myRank);
#endif // GMX_MPI

        numSharingSimulations_.resize(awhParams.numBias(), 1);
        sharingSimulationIndices_.resize(awhParams.numBias(), 0);
        multiSimCommPerBias_.resize(awhParams.numBias(), MPI_COMM_NULL);

        for (int shareIndex : globalShareIndices)
        {
            if (globalShareIndices.count(shareIndex) > 1)
            {
                const auto& findBiasIndex =
                        std::find(localShareIndices.begin(), localShareIndices.end(), shareIndex);
                const Index localBiasIndex = (findBiasIndex == localShareIndices.end()
                                                      ? -1
                                                      : findBiasIndex - localShareIndices.begin());
                MPI_Comm    splitComm;
                if (static_cast<int>(globalShareIndices.count(shareIndex)) == numSimulations)
                {
                    splitComm = simulationMainComm;
                }
                else
                {
#if GMX_MPI
                    const int haveLocally = (localBiasIndex >= 0 ? 1 : 0);
                    MPI_Comm_split(simulationMainComm, haveLocally, myRank, &splitComm);
                    createdCommList_.push_back(splitComm);
#else
                    GMX_RELEASE_ASSERT(false, "Can not have sharing without MPI");
#endif // GMX_MPI
                }
                if (localBiasIndex >= 0)
                {
                    numSharingSimulations_[localBiasIndex] = globalShareIndices.count(shareIndex);
#if GMX_MPI
                    MPI_Comm_rank(splitComm, &sharingSimulationIndices_[localBiasIndex]);
#endif // GMX_MPI
                    multiSimCommPerBias_[localBiasIndex] = splitComm;
                }
            }
        }
    }

#if GMX_MPI
    if (commRecord.nnodes > 1)
    {
        numSharingSimulations_.resize(awhParams.numBias());
        MPI_Bcast(
                numSharingSimulations_.data(), numSharingSimulations_.size(), MPI_INT, 0, commRecord.mpi_comm_mygroup);
    }
#endif // GMX_MPI
}

BiasSharing::~BiasSharing()
{
#if GMX_MPI
    for (MPI_Comm comm : createdCommList_)
    {
        MPI_Comm_free(&comm);
    }
#endif // GMX_MPI
}

namespace
{

#if GMX_MPI

template<typename T>
std::enable_if_t<std::is_same_v<T, int>, MPI_Datatype> mpiType()
{
    return MPI_INT;
}

template<typename T>
std::enable_if_t<std::is_same_v<T, long>, MPI_Datatype> mpiType()
{
    return MPI_LONG;
}

template<typename T>
std::enable_if_t<std::is_same_v<T, double>, MPI_Datatype> mpiType()
{
    return MPI_DOUBLE;
}

#endif // GMX_MPI

} // namespace

/*! \brief
 * Sum an array over all simulations on main ranks or all ranks of each simulation.
 *
 * This assumes the data is identical on all ranks within each simulation.
 *
 * \param[in,out] data          The data to sum.
 * \param[in]     multiSimComm  Communicator for the main ranks of sharing simulations.
 * \param[in]     broadcastWithinSimulation  Broadcast the result to all ranks within the simulation
 * \param[in]     commRecord    Struct for intra-simulation communication.
 */
template<typename T>
void sumOverSimulations(ArrayRef<T>      data,
                        MPI_Comm         multiSimComm,
                        const bool       broadcastWithinSimulation,
                        const t_commrec& commRecord)
{
#if GMX_MPI
    if (MAIN(&commRecord))
    {
        MPI_Allreduce(MPI_IN_PLACE, data.data(), data.size(), mpiType<T>(), MPI_SUM, multiSimComm);
    }
    if (broadcastWithinSimulation && commRecord.nnodes > 1)
    {
        gmx_bcast(data.size() * sizeof(T), data.data(), commRecord.mpi_comm_mygroup);
    }
#else
    GMX_UNUSED_VALUE(data);
    GMX_UNUSED_VALUE(commRecord);
    GMX_UNUSED_VALUE(broadcastWithinSimulation);
    GMX_UNUSED_VALUE(multiSimComm);
#endif // GMX_MPI
}

void BiasSharing::sumOverSharingMainRanks(ArrayRef<int> data, const int biasIndex) const
{
    sumOverSimulations(data, multiSimCommPerBias_[biasIndex], false, commRecord_);
}

void BiasSharing::sumOverSharingMainRanks(ArrayRef<long> data, const int biasIndex) const
{
    sumOverSimulations(data, multiSimCommPerBias_[biasIndex], false, commRecord_);
}

void BiasSharing::sumOverSharingMainRanks(ArrayRef<double> data, const int biasIndex) const
{
    sumOverSimulations(data, multiSimCommPerBias_[biasIndex], false, commRecord_);
}

void BiasSharing::sumOverSharingSimulations(ArrayRef<int> data, const int biasIndex) const
{
    sumOverSimulations(data, multiSimCommPerBias_[biasIndex], true, commRecord_);
}

void BiasSharing::sumOverSharingSimulations(ArrayRef<double> data, const int biasIndex) const
{
    sumOverSimulations(data, multiSimCommPerBias_[biasIndex], true, commRecord_);
}

bool haveBiasSharingWithinSimulation(const AwhParams& awhParams)
{
    bool haveSharing = false;

    for (int k = 0; k < awhParams.numBias(); k++)
    {
        int shareGroup = awhParams.awhBiasParams(k).shareGroup();
        if (shareGroup > 0)
        {
            for (int i = k + 1; i < awhParams.numBias(); i++)
            {
                if (awhParams.awhBiasParams(i).shareGroup() == shareGroup)
                {
                    haveSharing = true;
                }
            }
        }
    }

    return haveSharing;
}

void biasesAreCompatibleForSharingBetweenSimulations(const AwhParams&       awhParams,
                                                     ArrayRef<const size_t> pointSize,
                                                     const BiasSharing&     biasSharing)
{
    /* Check the point sizes. This is a sufficient condition for running
     * as shared multi-sim run. No physics checks are performed here.
     */
    const auto& awhBiasParams = awhParams.awhBiasParams();
    for (int b = 0; b < gmx::ssize(awhBiasParams); b++)
    {
        if (awhBiasParams[b].shareGroup() > 0)
        {
            const int numSim = biasSharing.numSharingSimulations(b);
            if (numSim == 1)
            {
                // This bias is not actually shared
                continue;
            }
            const int        simIndex = biasSharing.sharingSimulationIndex(b);
            std::vector<int> intervals(numSim * 2);
            intervals[numSim * 0 + simIndex] = awhParams.nstSampleCoord();
            intervals[numSim * 1 + simIndex] = awhParams.numSamplesUpdateFreeEnergy();
            biasSharing.sumOverSharingMainRanks(intervals, b);
            for (int sim = 1; sim < numSim; sim++)
            {
                if (intervals[sim] != intervals[0])
                {
                    GMX_THROW(InvalidInputError(
                            "All simulations should have the same AWH sample interval"));
                }
                if (intervals[numSim + sim] != intervals[numSim])
                {
                    GMX_THROW(
                            InvalidInputError("All simulations should have the same AWH "
                                              "free-energy update interval"));
                }
            }

            std::vector<long> pointSizes(numSim);
            pointSizes[simIndex] = pointSize[b];
            biasSharing.sumOverSharingMainRanks(pointSizes, b);
            for (int sim = 1; sim < numSim; sim++)
            {
                if (pointSizes[sim] != pointSizes[0])
                {
                    GMX_THROW(InvalidInputError(
                            gmx::formatString("Shared AWH bias %d has different grid sizes in "
                                              "different simulations\n",
                                              b + 1)));
                }
            }
        }
    }
}

} // namespace gmx
