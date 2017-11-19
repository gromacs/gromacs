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

/*! \internal \file
 * \brief
 * Implements the Awh class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "awh.h"

#include <assert.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/enxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"

#include "bias.h"
#include "biassharing.h"
#include "correlationgrid.h"
#include "pointstate.h"

namespace gmx
{

/*! \internal
 * \brief A bias and its coupling to the system.
 *
 * This struct is used to separate the bias machinery in the Bias class,
 * which should be independent from the reaction coordinate, from the
 * obtaining of the reaction coordinate values and passing the computed forces.
 * Currently the AWH method couples to the system by mapping each
 * AWH bias to a pull coordinate. This can easily be generalized here.
 */
struct BiasCoupledToSystem
{
    /*! \brief Constructor, couple a bias to a set of pull coordinates.
     *
     * \param[in] bias            The bias.
     * \param[in] pullCoordIndex  The pull coordinate indices.
     */
    BiasCoupledToSystem(Bias                    bias,
                        const std::vector<int> &pullCoordIndex);

    Bias                   bias;           /**< The bias. */
    const std::vector<int> pullCoordIndex; /**< The pull coordinates this bias acts on. */

    /* Here AWH can be extended to work on other coordinates than pull. */
};

BiasCoupledToSystem::BiasCoupledToSystem(Bias                    bias,
                                         const std::vector<int> &pullCoordIndex) :
    bias(std::move(bias)),
    pullCoordIndex(pullCoordIndex)
{
    /* We already checked for this in grompp, but check again here. */
    GMX_RELEASE_ASSERT(static_cast<size_t>(bias.ndim()) == pullCoordIndex.size(), "The bias dimensionality should match the number of pull coordinates.");
}

Awh::Awh(FILE              *fplog,
         const t_inputrec  &inputRecord,
         const t_commrec   *commRecord,
         const AwhParams   &awhParams,
         const std::string &biasInitFilename,
         pull_t            *pull_work) :
    seed_(awhParams.seed),
    nstout_(awhParams.nstOut),
    commRecord_(commRecord),
    pull_(pull_work),
    potentialOffset_(0)
{
    /* We already checked for this in grompp, but check again here. */
    GMX_RELEASE_ASSERT(inputRecord.pull != nullptr, "With AWH we should have pull parameters");
    GMX_RELEASE_ASSERT(pull_work != nullptr, "With AWH pull should be initialized before initializing AWH");

    if (fplog != nullptr)
    {
        please_cite(fplog, "Lindahl2014");
    }

    if (haveBiasSharingWithinSimulation(awhParams))
    {
        /* This has likely been checked by grompp, but throw anyhow. */
        GMX_THROW(InvalidInputError("Biases within a simulation are shared, currently sharing of biases is only supported between simulations"));
    }

    int numSharingSimulations = 1;
    if (awhParams.shareBiasMultisim && MULTISIM(commRecord_))
    {
        numSharingSimulations = commRecord_->ms->nsim;
    }

    /* Initialize all the biases */
    const double beta = 1/(BOLTZ*inputRecord.opts.ref_t[0]);
    for (int k = 0; k < awhParams.numBias; k++)
    {
        const AwhBiasParams    &awhBiasParams = awhParams.awhBiasParams[k];

        std::vector<int>        pullCoordIndex;
        std::vector<DimParams>  dimParams;
        for (int d = 0; d < awhBiasParams.ndim; d++)
        {
            const AwhDimParams &awhDimParams      = awhBiasParams.dimParams[d];
            GMX_RELEASE_ASSERT(awhDimParams.eCoordProvider == eawhcoordproviderPULL, "Currently only the pull code is supported as coordinate provider");
            const t_pull_coord &pullCoord         = inputRecord.pull->coord[awhDimParams.coordIndex];
            double              conversionFactor  = pull_coordinate_is_angletype(&pullCoord) ? DEG2RAD : 1;
            dimParams.push_back(DimParams(conversionFactor, awhDimParams.forceConstant, beta));

            pullCoordIndex.push_back(awhDimParams.coordIndex);
        }

        /* Construct the bias and couple it to the system. */
        Bias::ThisRankWillDoIO thisRankWillDoIO = (MASTER(commRecord_) ? Bias::ThisRankWillDoIO::Yes : Bias::ThisRankWillDoIO::No);
        biasCoupledToSystem_.emplace_back(Bias(k, awhParams, awhParams.awhBiasParams[k], dimParams, beta, inputRecord.delta_t, numSharingSimulations, biasInitFilename, thisRankWillDoIO),
                                          pullCoordIndex);

        biasCoupledToSystem_.back().bias.printInitializationToLog(fplog);
    }

    /* Need to register the AWH coordinates to be allowed to apply forces to the pull coordinates. */
    registerAwhWithPull(awhParams, pull_);

    if (numSharingSimulations > 1 && MASTER(commRecord_))
    {
        std::vector<size_t> pointSize;
        for (auto const &biasCts : biasCoupledToSystem_)
        {
            pointSize.push_back(biasCts.bias.state().points().size());
        }
        /* Ensure that the shared biased are compatible between simulations */
        biasesAreCompatibleForSharingBetweenSimulations(awhParams, pointSize, commRecord_->ms);
    }
}

Awh::~Awh() = default;

bool Awh::isOutputStep(gmx_int64_t step) const
{
    return (nstout_ > 0 && step % nstout_ == 0);
}

real Awh::applyBiasForcesAndUpdateBias(int                     ePBC,
                                       const t_mdatoms        &mdatoms,
                                       const matrix            box,
                                       gmx::ForceWithVirial   *forceWithVirial,
                                       double                  t,
                                       gmx_int64_t             step,
                                       gmx_wallcycle          *wallcycle,
                                       FILE                   *fplog)
{
    GMX_ASSERT(forceWithVirial, "Need a valid ForceWithVirial object");

    wallcycle_start(wallcycle, ewcAWH);

    t_pbc  pbc;
    set_pbc(&pbc, ePBC, box);

    /* During the AWH update the potential can instantaneously jump due to either
       an bias update or moving the umbrella. The jumps are kept track of and
       subtracted from the potential in order to get a useful conserved energy quantity. */
    double awhPotential = potentialOffset_;

    for (auto &biasCts : biasCoupledToSystem_)
    {
        /* Update the AWH coordinate values with those of the corresponding
         * pull coordinates.
         */
        awh_dvec coordValue = { 0, 0, 0, 0 };
        for (int d = 0; d < biasCts.bias.ndim(); d++)
        {
            coordValue[d] = get_pull_coord_value(pull_, biasCts.pullCoordIndex[d], &pbc);
        }

        /* Perform an AWH biasing step: this means, at regular intervals,
         * sampling observables based on the input pull coordinate value,
         * setting the bias force and/or updating the AWH bias state.
         */
        double              biasPotential;
        double              biasPotentialJump;
        /* Note: In the near future this call will be split in calls
         *       to supports bias sharing within a single simulation.
         */
        gmx::ArrayRef<const double> biasForce =
            biasCts.bias.calcForceAndUpdateBias(coordValue,
                                                &biasPotential, &biasPotentialJump,
                                                commRecord_->ms,
                                                t, step, seed_, fplog);

        awhPotential += biasPotential;

        /* Keep track of the total potential shift needed to remove the potential jumps. */
        potentialOffset_ -= biasPotentialJump;

        /* Communicate the bias force to the pull struct.
         * The bias potential is returned at the end of this function,
         * so that it can be added externally to the correct energy data block.
         */
        for (int d = 0; d < biasCts.bias.ndim(); d++)
        {
            apply_external_pull_coord_force(pull_, biasCts.pullCoordIndex[d],
                                            biasForce[d], &mdatoms,
                                            forceWithVirial);
        }

        if (isOutputStep(step))
        {
            /* We might have skipped updates for part of the grid points.
             * Ensure all points are updated before writing out their data.
             */
            biasCts.bias.doSkippedUpdatesForAllPoints();
        }
    }

    wallcycle_stop(wallcycle, ewcAWH);

    return MASTER(commRecord_) ? static_cast<real>(awhPotential) : 0;
}

std::shared_ptr<AwhHistory> Awh::initHistoryFromState() const
{
    if (MASTER(commRecord_))
    {
        std::shared_ptr<AwhHistory> awhHistory(new AwhHistory);
        awhHistory->bias.clear();
        awhHistory->bias.resize(biasCoupledToSystem_.size());

        for (size_t k = 0; k < awhHistory->bias.size(); k++)
        {
            biasCoupledToSystem_[k].bias.initHistoryFromState(&awhHistory->bias[k]);
        }

        return awhHistory;
    }
    else
    {
        /* Return an empty pointer */
        return std::shared_ptr<AwhHistory>();
    }
}

void Awh::restoreStateFromHistory(const AwhHistory *awhHistory)
{
    /* Restore the history to the current state */
    if (MASTER(commRecord_))
    {
        GMX_RELEASE_ASSERT(awhHistory != nullptr, "The master rank should have a valid awhHistory when restoring the state from history.");

        if (awhHistory->bias.size() != biasCoupledToSystem_.size())
        {
            GMX_THROW(InvalidInputError("AWH state and history contain different numbers of biases. Likely you provided a checkpoint from a different simulation."));
        }

        potentialOffset_ = awhHistory->potentialOffset;
    }
    if (PAR(commRecord_))
    {
        gmx_bcast(sizeof(potentialOffset_), &potentialOffset_, commRecord_);
    }

    for (size_t k = 0; k < biasCoupledToSystem_.size(); k++)
    {
        biasCoupledToSystem_[k].bias.restoreStateFromHistory(awhHistory ? &awhHistory->bias[k] : nullptr, commRecord_);
    }
}

void Awh::updateHistory(AwhHistory *awhHistory) const
{
    if (!MASTER(commRecord_))
    {
        return;
    }

    /* This assert will also catch a non-master rank calling this function. */
    GMX_RELEASE_ASSERT(awhHistory->bias.size() == biasCoupledToSystem_.size(), "AWH state and history bias count should match");

    awhHistory->potentialOffset = potentialOffset_;

    for (size_t k = 0; k < awhHistory->bias.size(); k++)
    {
        biasCoupledToSystem_[k].bias.updateHistory(&awhHistory->bias[k]);
    }
}

const char * Awh::externalPotentialString()
{
    return "AWH";
}

void Awh::registerAwhWithPull(const AwhParams &awhParams,
                              pull_t          *pull_work)
{
    GMX_RELEASE_ASSERT(pull_work, "Need a valid pull object");

    for (int k = 0; k < awhParams.numBias; k++)
    {
        const AwhBiasParams &biasParams = awhParams.awhBiasParams[k];

        for (int d = 0; d < biasParams.ndim; d++)
        {
            register_external_pull_potential(pull_work, biasParams.dimParams[d].coordIndex, Awh::externalPotentialString());
        }
    }
}

/* Fill the AWH data block of an energy frame with data (if there is any). */
void Awh::writeToEnergyFrame(gmx_int64_t  step,
                             t_enxframe  *frame) const
{
    GMX_ASSERT(MASTER(commRecord_), "writeToEnergyFrame should only be called on the master rank");
    GMX_ASSERT(frame != nullptr, "Need a valid energy frame");

    if (!isOutputStep(step))
    {
        /* This is not an AWH output step, don't write any AWH data */
        return;
    }

    /* Get the total number of energy subblocks that AWH needs */
    int numSubblocks  = 0;
    for (auto &biasCoupledToSystem : biasCoupledToSystem_)
    {
        numSubblocks += biasCoupledToSystem.bias.numEnergySubblocksToWrite();
    }
    GMX_ASSERT(numSubblocks > 0, "We should always have data to write");

    /* Add 1 energy block */
    add_blocks_enxframe(frame, frame->nblock + 1);

    /* Take the block that was just added and set the number of subblocks. */
    t_enxblock *awhEnergyBlock = &(frame->block[frame->nblock - 1]);
    add_subblocks_enxblock(awhEnergyBlock, numSubblocks);

    /* Claim it as an AWH block. */
    awhEnergyBlock->id = enxAWH;

    /* Transfer AWH data blocks to energy sub blocks */
    int energySubblockCount = 0;
    for (auto &biasCoupledToSystem : biasCoupledToSystem_)
    {
        energySubblockCount += biasCoupledToSystem.bias.writeToEnergySubblocks(&(awhEnergyBlock->sub[energySubblockCount]));
    }
}

} // namespace gmx
