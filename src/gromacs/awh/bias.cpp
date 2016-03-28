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
 * Implements the Bias class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "bias.h"

#include <assert.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "correlationgrid.h"
#include "correlationhistory.h"
#include "pointstate.h"

namespace gmx
{

void Bias::warnForHistogramAnomalies(double t, gmx_int64_t step, FILE *fplog)
{
    const int    maxNumWarningsInCheck = 1;   /* The maximum number of warnings to print per check */
    const int    maxNumWarningsInRun   = 10;  /* The maximum number of warnings to print in a run */

    if (fplog == nullptr || numWarningsIssued_ >= maxNumWarningsInRun || state_.inInitialStage() ||
        !params_.isCheckStep(state_.points().size(), step))
    {
        return;
    }

    numWarningsIssued_ +=
        state_.warnForHistogramAnomalies(grid_, biasIndex(), t, fplog,
                                         maxNumWarningsInCheck);

    if (numWarningsIssued_ >= maxNumWarningsInRun)
    {
        fprintf(fplog, "\nawh%d: suppressing future AWH warnings.\n", biasIndex() + 1);
    }
}

void Bias::doSkippedUpdatesForAllPoints()
{
    if (params_.skipUpdates())
    {
        state_.doSkippedUpdatesForAllPoints(params_);
    }
}

gmx::ArrayRef<const double>
Bias::calcForceAndUpdateBias(const awh_dvec        coordValue,
                             double               *awhPotential,
                             double               *potentialJump,
                             const gmx_multisim_t *ms,
                             double                t,
                             gmx_int64_t           step,
                             gmx_int64_t           seed,
                             FILE                 *fplog)
{
    if (step < 0)
    {
        GMX_THROW(InvalidInputError("The step number is negative which is not supported by the AWH code."));
    }

    state_.setCoordValue(grid_, coordValue);

    std::vector < double, AlignedAllocator < double>> &probWeightNeighbor = alignedTempWorkSpace_;

    /* If the convolved force is needed or this is a sampling step,
     * the bias in the current neighborhood needs to be up-to-date
     * and the probablity weights need to be calculated.
     */
    const bool sampleCoord   = params_.isSampleCoordStep(step);
    const bool moveUmbrella  = (sampleCoord || step == 0);
    double     convolvedBias = 0;
    if (params_.convolveForce || moveUmbrella || sampleCoord)
    {
        if (params_.skipUpdates())
        {
            state_.doSkippedUpdatesInNeighborhood(params_, grid_);
        }

        convolvedBias = state_.updateProbabilityWeightsAndConvolvedBias(dimParams_, grid_, &probWeightNeighbor);

        if (sampleCoord)
        {
            updateForceCorrelationGrid(probWeightNeighbor, t);

            state_.sampleCoordAndPmf(grid_, probWeightNeighbor, convolvedBias);
        }
    }

    const CoordState &coordState = state_.coordState();

    /* Set the bias force and get the potential contribution from this bias.
     * The potential jump occurs at different times depending on how
     * the force is applied (and how the potential is normalized).
     * For the convolved force it happens when the bias is updated,
     * for the umbrella when the umbrella is moved.
     */
    *potentialJump = 0;
    double potential;
    if (params_.convolveForce)
    {
        state_.calcConvolvedForce(dimParams_, grid_, probWeightNeighbor,
                                  tempForce_, biasForce_);

        potential = -convolvedBias*params_.invBeta;
    }
    else
    {
        /* Umbrella force */
        GMX_RELEASE_ASSERT(state_.points()[coordState.umbrellaGridpoint()].inTargetRegion(),
                           "AWH bias grid point for the umbrella reference value is outside of the target region.");
        potential =
            state_.calcUmbrellaForceAndPotential(dimParams_, grid_, coordState.umbrellaGridpoint(), biasForce_);

        /* Moving the umbrella results in a force correction and
         * a new potential. The umbrella center is sampled as often as
         * the coordinate so we know the probability weights needed
         * for moving the umbrella are up-to-date.
         */
        if (moveUmbrella)
        {
            double newPotential = state_.moveUmbrella(dimParams_, grid_, probWeightNeighbor, biasForce_, step, seed, params_.biasIndex);
            *potentialJump      = newPotential - potential;
        }
    }

    /* Update the free energy estimates and bias and other history dependent method parameters */
    if (params_.isUpdateFreeEnergyStep(step))
    {
        state_.updateFreeEnergyAndAddSamplesToHistogram(dimParams_, grid_,
                                                        params_,
                                                        ms, t, step, fplog,
                                                        &updateList_);

        if (params_.convolveForce)
        {
            /* The update results in a potential jump, so we need the new convolved potential. */
            double newPotential = -calcConvolvedBias(coordState.coordValue())*params_.invBeta;
            *potentialJump      = newPotential - potential;
        }
    }

    /* Return the potential. */
    *awhPotential = potential;

    /* Check the sampled histograms and potentially warn user if something is suspicious */
    warnForHistogramAnomalies(t, step, fplog);

    return biasForce_;
}

void Bias::restoreStateFromHistory(const AwhBiasHistory *biasHistory,
                                   const t_commrec      *cr)
{
    GMX_RELEASE_ASSERT(thisRankDoesIO_ == MASTER(cr), "The master rank should do I/O, the other ranks should not");

    if (MASTER(cr))
    {
        GMX_RELEASE_ASSERT(biasHistory != nullptr, "On the master rank we need a valid history object to restore from");
        state_.restoreFromHistory(*biasHistory, grid_);

        if (forceCorrelationGrid_ != nullptr)
        {
            forceCorrelationGrid_->restoreStateFromHistory(biasHistory->forceCorrelationGrid);
        }
    }

    if (PAR(cr))
    {
        state_.broadcast(cr);
    }
}

void Bias::initHistoryFromState(AwhBiasHistory *biasHistory) const
{
    GMX_RELEASE_ASSERT(biasHistory != nullptr, "Need a valid biasHistory");

    state_.initHistoryFromState(biasHistory);

    if (forceCorrelationGrid_ != nullptr)
    {
        biasHistory->forceCorrelationGrid = initCorrelationGridHistoryFromState(forceCorrelationGrid());
    }
}

void Bias::updateHistory(AwhBiasHistory *biasHistory) const
{
    GMX_RELEASE_ASSERT(biasHistory != nullptr, "Need a valid biasHistory");

    state_.updateHistory(biasHistory, grid_);

    if (forceCorrelationGrid_ != nullptr)
    {
        updateCorrelationGridHistory(&biasHistory->forceCorrelationGrid, forceCorrelationGrid());
    }
}

Bias::Bias(int                             biasIndexInCollection,
           const AwhParams                &awhParams,
           const AwhBiasParams            &awhBiasParams,
           const std::vector<DimParams>   &dimParamsInit,
           double                          beta,
           double                          mdTimeStep,
           int                             numSharingSimulations,
           const std::string              &biasInitFilename,
           ThisRankWillDoIO                thisRankWillDoIO,
           BiasParams::DisableUpdateSkips  disableUpdateSkips) :
    dimParams_(dimParamsInit),
    grid_(dimParamsInit, awhBiasParams.dimParams),
    params_(awhParams, awhBiasParams, dimParams_, beta, mdTimeStep, disableUpdateSkips, numSharingSimulations, grid_.axis(), biasIndexInCollection),
    state_(awhBiasParams, params_.initialHistogramSize, dimParams_, grid_),
    thisRankDoesIO_(thisRankWillDoIO == ThisRankWillDoIO::Yes),
    biasForce_(ndim()),
    alignedTempWorkSpace_(),
    tempForce_(ndim()),
    numWarningsIssued_(0)
{
    /* For a global update updateList covers all points, so reserve that */
    updateList_.reserve(grid_.numPoints());

    state_.initGridPointState(awhBiasParams, dimParams_, grid_, params_, biasInitFilename, awhParams.numBias);

    if (thisRankDoesIO_)
    {
        /* Set up the force correlation object. */

        /* We let the correlation init function set its parameters
         * to something useful for now.
         */
        double blockLength = 0;
        /* Construct the force correlation object. */
        forceCorrelationGrid_ =
            gmx::compat::make_unique<CorrelationGrid>(state_.points().size(), ndim(),
                                                      blockLength, CorrelationGrid::BlockLengthMeasure::Time,
                                                      awhParams.nstSampleCoord*mdTimeStep);

        writer_ = std::unique_ptr<BiasWriter>(new BiasWriter(*this));
    }
}

void Bias::printInitializationToLog(FILE *fplog) const
{
    if (fplog != nullptr && forceCorrelationGrid_ != nullptr)
    {
        std::string prefix =
            gmx::formatString("\nawh%d:", params_.biasIndex + 1);

        fprintf(fplog,
                "%s initial force correlation block length = %g %s"
                "%s force correlation number of blocks = %d",
                prefix.c_str(), forceCorrelationGrid().getBlockLength(),
                forceCorrelationGrid().blockLengthMeasure == CorrelationGrid::BlockLengthMeasure::Weight ? "" : "ps",
                prefix.c_str(), forceCorrelationGrid().getNumBlocks());
    }
}

void Bias::updateForceCorrelationGrid(gmx::ArrayRef<const double> probWeightNeighbor,
                                      double                      t)
{
    if (forceCorrelationGrid_ == nullptr)
    {
        return;
    }

    const std::vector<int> &neighbor = grid_.point(state_.coordState().gridpointIndex()).neighbor;

    gmx::ArrayRef<double>   forceFromNeighbor = tempForce_;
    for (size_t n = 0; n < neighbor.size(); n++)
    {
        double weightNeighbor = probWeightNeighbor[n];
        int    indexNeighbor  = neighbor[n];

        /* Add the force data of this neighbor point. Note: the sum of these forces is the convolved force.

           We actually add the force normalized by beta which has the units of 1/length. This means that the
           resulting correlation time integral is directly in units of friction time/length^2 which is really what
           we're interested in. */
        state_.calcUmbrellaForceAndPotential(dimParams_, grid_, indexNeighbor, forceFromNeighbor);

        /* Note: we might want to give a whole list of data to add instead and have this loop in the data adding function */
        forceCorrelationGrid_->addData(indexNeighbor, weightNeighbor, forceFromNeighbor, t);
    }
}

/* Return the number of data blocks that have been prepared for writing. */
int Bias::numEnergySubblocksToWrite() const
{
    GMX_RELEASE_ASSERT(writer_ != nullptr, "Should only request data from an initialized writer");

    return writer_->numBlocks();
}

/* Write bias data blocks to energy subblocks. */
int Bias::writeToEnergySubblocks(t_enxsubblock *subblock) const
{
    GMX_RELEASE_ASSERT(writer_ != nullptr, "Should only request data from an initialized writer");

    return writer_->writeToEnergySubblocks(*this, subblock);
}

} // namespace gmx
