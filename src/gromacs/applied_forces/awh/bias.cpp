/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Implements the Bias class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "bias.h"

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <memory>

#include "gromacs/applied_forces/awh/biasgrid.h"
#include "gromacs/applied_forces/awh/biasparams.h"
#include "gromacs/applied_forces/awh/biasstate.h"
#include "gromacs/applied_forces/awh/biaswriter.h"
#include "gromacs/applied_forces/awh/coordstate.h"
#include "gromacs/applied_forces/awh/dimparams.h"
#include "gromacs/applied_forces/awh/histogramsize.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/awh_correlation_history.h"
#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "biassharing.h"
#include "correlationgrid.h"
#include "correlationhistory.h"
#include "pointstate.h"

namespace gmx
{

void Bias::warnForHistogramAnomalies(double t, int64_t step, FILE* fplog)
{
    const int maxNumWarningsInCheck = 1;  /* The maximum number of warnings to print per check */
    const int maxNumWarningsInRun   = 10; /* The maximum number of warnings to print in a run */

    if (fplog == nullptr || numWarningsIssued_ >= maxNumWarningsInRun || state_.inInitialStage()
        || !params_.isCheckHistogramForAnomaliesStep(step))
    {
        return;
    }

    numWarningsIssued_ +=
            state_.warnForHistogramAnomalies(grid_, biasIndex(), t, fplog, maxNumWarningsInCheck);

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

gmx::ArrayRef<const double> Bias::calcForceAndUpdateBias(const awh_dvec         coordValue,
                                                         ArrayRef<const double> neighborLambdaEnergies,
                                                         ArrayRef<const double> neighborLambdaDhdl,
                                                         double*                awhPotential,
                                                         double*                potentialJump,
                                                         double                 t,
                                                         int64_t                step,
                                                         int64_t                seed,
                                                         FILE*                  fplog)
{
    if (step < 0)
    {
        GMX_THROW(InvalidInputError(
                "The step number is negative which is not supported by the AWH code."));
    }

    GMX_RELEASE_ASSERT(!(params_.convolveForce && grid_.hasLambdaAxis()),
                       "When using AWH to sample an FEP lambda dimension the AWH potential cannot "
                       "be convolved.");

    state_.setCoordValue(grid_, coordValue);

    std::vector<double, AlignedAllocator<double>>& probWeightNeighbor = alignedTempWorkSpace_;

    /* If the convolved force is needed or this is a sampling step,
     * the bias in the current neighborhood needs to be up-to-date
     * and the probablity weights need to be calculated.
     */
    const bool        isSampleCoordStep = params_.isSampleCoordStep(step);
    const bool        moveUmbrella      = (isSampleCoordStep || step == 0);
    double            convolvedBias     = 0;
    const CoordState& coordState        = state_.coordState();

    if (params_.convolveForce || moveUmbrella || isSampleCoordStep)
    {
        if (params_.skipUpdates())
        {
            state_.doSkippedUpdatesInNeighborhood(params_, grid_);
        }
        convolvedBias = state_.updateProbabilityWeightsAndConvolvedBias(
                dimParams_, grid_, moveUmbrella ? neighborLambdaEnergies : ArrayRef<const double>{}, &probWeightNeighbor);

        if (isSampleCoordStep)
        {
            updateForceCorrelationGrid(probWeightNeighbor, neighborLambdaDhdl, t);

            state_.sampleCoordAndPmf(dimParams_, grid_, probWeightNeighbor, convolvedBias);
        }
    }

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
        state_.calcConvolvedForce(dimParams_,
                                  grid_,
                                  probWeightNeighbor,
                                  moveUmbrella ? neighborLambdaDhdl : ArrayRef<const double>{},
                                  tempForce_,
                                  biasForce_);

        potential = -convolvedBias * params_.invBeta;
    }
    else
    {
        /* Umbrella force */
        GMX_RELEASE_ASSERT(state_.points()[coordState.umbrellaGridpoint()].inTargetRegion(),
                           "AWH bias grid point for the umbrella reference value is outside of the "
                           "target region.");
        potential = state_.calcUmbrellaForceAndPotential(
                dimParams_,
                grid_,
                coordState.umbrellaGridpoint(),
                moveUmbrella ? neighborLambdaDhdl : ArrayRef<const double>{},
                biasForce_);

        /* Moving the umbrella results in a force correction and
         * a new potential. The umbrella center is sampled as often as
         * the coordinate so we know the probability weights needed
         * for moving the umbrella are up-to-date.
         */
        if (moveUmbrella)
        {
            const bool onlySampleUmbrellaGridpoint = false;
            double     newPotential                = state_.moveUmbrella(dimParams_,
                                                      grid_,
                                                      probWeightNeighbor,
                                                      neighborLambdaDhdl,
                                                      biasForce_,
                                                      step,
                                                      seed,
                                                      params_.biasIndex_,
                                                      onlySampleUmbrellaGridpoint);
            *potentialJump                         = newPotential - potential;
        }
    }

    /* Update the free energy estimates and bias and other history dependent method parameters */
    if (params_.isUpdateFreeEnergyStep(step))
    {
        state_.updateFreeEnergyAndAddSamplesToHistogram(
                dimParams_, grid_, params_, forceCorrelationGrid(), t, step, fplog, &updateList_);

        if (params_.convolveForce)
        {
            /* The update results in a potential jump, so we need the new convolved potential. */
            double newPotential = -calcConvolvedBias(coordState.coordValue()) * params_.invBeta;
            *potentialJump      = newPotential - potential;
        }
    }
    /* If there is a lambda axis it is still controlled using an umbrella even if the force
     * is convolved in the other dimensions. */
    if (moveUmbrella && params_.convolveForce && grid_.hasLambdaAxis())
    {
        const bool onlySampleUmbrellaGridpoint = true;
        state_.moveUmbrella(dimParams_,
                            grid_,
                            probWeightNeighbor,
                            neighborLambdaDhdl,
                            biasForce_,
                            step,
                            seed,
                            params_.biasIndex_,
                            onlySampleUmbrellaGridpoint);
    }

    /* Return the potential. */
    *awhPotential = potential;

    /* Check the sampled histograms and potentially warn user if something is suspicious */
    warnForHistogramAnomalies(t, step, fplog);

    return biasForce_;
}

/*! \brief
 * Count the total number of samples / sample weight over all grid points.
 *
 * \param[in] pointState  The state of the points in a bias.
 * \returns the total sample count.
 */
static int64_t countSamples(ArrayRef<const PointState> pointState)
{
    double numSamples = 0;
    for (const PointState& point : pointState)
    {
        numSamples += point.weightSumTot();
    }

    return gmx::roundToInt64(numSamples);
}

/*! \brief
 * Check if the state (loaded from checkpoint) and the run are consistent.
 *
 * When the state and the run setup are inconsistent, an exception is thrown.
 *
 * \param[in] params  The parameters of the bias.
 * \param[in] state   The state of the bias.
 */
static void ensureStateAndRunConsistency(const BiasParams& params, const BiasState& state)
{
    int64_t numSamples = countSamples(state.points());
    int64_t numUpdatesFromSamples =
            numSamples / (params.numSamplesUpdateFreeEnergy_ * params.numSharedUpdate);
    int64_t numUpdatesExpected = state.histogramSize().numUpdates();
    if (numUpdatesFromSamples != numUpdatesExpected)
    {
        std::string mesg = gmx::formatString(
                "The number of AWH updates in the checkpoint file (%" PRId64
                ") does not match the total number of AWH samples divided by the number of samples "
                "per update for %d sharing AWH bias(es) (%" PRId64 "/%d=%" PRId64 ")",
                numUpdatesExpected,
                params.numSharedUpdate,
                numSamples,
                params.numSamplesUpdateFreeEnergy_ * params.numSharedUpdate,
                numUpdatesFromSamples);
        mesg += " Maybe you changed AWH parameters.";
        /* Unfortunately we currently do not store the number of simulations
         * sharing the bias or the state to checkpoint. But we can hint at
         * a mismatch in the number of sharing simulations.
         */
        if (numUpdatesFromSamples % state.histogramSize().numUpdates() == 0)
        {
            mesg += gmx::formatString(
                    " Or the run you continued from used %" PRId64
                    " sharing simulations, whereas you now specified %d sharing simulations.",
                    numUpdatesFromSamples / state.histogramSize().numUpdates(),
                    params.numSharedUpdate);
        }
        GMX_THROW(InvalidInputError(mesg));
    }
}

void Bias::restoreStateFromHistory(const AwhBiasHistory* biasHistory, const t_commrec* cr)
{
    GMX_RELEASE_ASSERT(thisRankDoesIO_ == MAIN(cr),
                       "The main rank should do I/O, the other ranks should not");

    if (MAIN(cr))
    {
        GMX_RELEASE_ASSERT(biasHistory != nullptr,
                           "On the main rank we need a valid history object to restore from");
        state_.restoreFromHistory(*biasHistory, grid_);

        /* Ensure that the state is consistent with our current run setup,
         * since the user can have changed AWH parameters or the number
         * of simulations sharing the bias.
         */
        ensureStateAndRunConsistency(params_, state_);

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

void Bias::initHistoryFromState(AwhBiasHistory* biasHistory) const
{
    GMX_RELEASE_ASSERT(biasHistory != nullptr, "Need a valid biasHistory");

    state_.initHistoryFromState(biasHistory);

    if (forceCorrelationGrid_ != nullptr)
    {
        biasHistory->forceCorrelationGrid = initCorrelationGridHistoryFromState(forceCorrelationGrid());
    }
}

void Bias::updateHistory(AwhBiasHistory* biasHistory) const
{
    GMX_RELEASE_ASSERT(biasHistory != nullptr, "Need a valid biasHistory");

    state_.updateHistory(biasHistory, grid_);

    if (forceCorrelationGrid_ != nullptr)
    {
        updateCorrelationGridHistory(&biasHistory->forceCorrelationGrid, forceCorrelationGrid());
    }
}

Bias::Bias(int                            biasIndexInCollection,
           const AwhParams&               awhParams,
           const AwhBiasParams&           awhBiasParams,
           ArrayRef<const DimParams>      dimParamsInit,
           double                         beta,
           double                         mdTimeStep,
           const BiasSharing*             biasSharing,
           const std::string&             biasInitFilename,
           ThisRankWillDoIO               thisRankWillDoIO,
           BiasParams::DisableUpdateSkips disableUpdateSkips) :
    dimParams_(dimParamsInit.begin(), dimParamsInit.end()),
    grid_(dimParamsInit, awhBiasParams.dimParams()),
    params_(awhParams,
            awhBiasParams,
            dimParams_,
            beta,
            mdTimeStep,
            disableUpdateSkips,
            biasSharing ? biasSharing->numSharingSimulations(biasIndexInCollection) : 1,
            grid_.axis(),
            biasIndexInCollection),
    state_(awhBiasParams, params_.initialHistogramSize, dimParams_, grid_, biasSharing),
    thisRankDoesIO_(thisRankWillDoIO == ThisRankWillDoIO::Yes),
    biasForce_(ndim()),
    tempForce_(ndim()),
    numWarningsIssued_(0)
{
    /* For a global update updateList covers all points, so reserve that */
    updateList_.reserve(grid_.numPoints());

    /* Set up the force correlation object. */

    /* We let the correlation init function set its parameters
     * to something useful for now.
     */
    double blockLength = 0;
    /* Construct the force correlation object. */
    forceCorrelationGrid_ = std::make_unique<CorrelationGrid>(state_.points().size(),
                                                              ndim(),
                                                              blockLength,
                                                              CorrelationGrid::BlockLengthMeasure::Time,
                                                              awhParams.nstSampleCoord() * mdTimeStep);

    state_.initGridPointState(awhBiasParams,
                              dimParams_,
                              grid_,
                              params_,
                              forceCorrelationGrid(),
                              biasInitFilename,
                              awhParams.numBias());

    if (thisRankDoesIO_)
    {
        writer_ = std::make_unique<BiasWriter>(*this);
    }
}

void Bias::printInitializationToLog(FILE* fplog) const
{
    if (fplog != nullptr && forceCorrelationGrid_ != nullptr)
    {
        std::string prefix = gmx::formatString("\nawh%d:", params_.biasIndex_ + 1);

        fprintf(fplog,
                "%s initial force correlation block length = %g %s"
                "%s force correlation number of blocks = %d",
                prefix.c_str(),
                forceCorrelationGrid().getBlockLength(),
                forceCorrelationGrid().blockLengthMeasure_ == CorrelationGrid::BlockLengthMeasure::Weight
                        ? ""
                        : "ps",
                prefix.c_str(),
                forceCorrelationGrid().getNumBlocks());
    }
}

void Bias::updateForceCorrelationGrid(gmx::ArrayRef<const double> probWeightNeighbor,
                                      ArrayRef<const double>      neighborLambdaDhdl,
                                      double                      t)
{
    if (forceCorrelationGrid_ == nullptr)
    {
        return;
    }

    const std::vector<int>& neighbor = grid_.point(state_.coordState().gridpointIndex()).neighbor;

    gmx::ArrayRef<double> forceFromNeighbor = tempForce_;
    for (size_t n = 0; n < neighbor.size(); n++)
    {
        double weightNeighbor = probWeightNeighbor[n];
        int    indexNeighbor  = neighbor[n];

        /* Add the force data of this neighbor point. Note: the sum of these forces is the convolved force.

           We actually add the force normalized by beta which has the units of 1/length. This means that the
           resulting correlation time integral is directly in units of friction time/length^2 which is really what
           we're interested in. */
        state_.calcUmbrellaForceAndPotential(
                dimParams_, grid_, indexNeighbor, neighborLambdaDhdl, forceFromNeighbor);

        /* Note: we might want to give a whole list of data to add instead and have this loop in the data adding function */
        forceCorrelationGrid_->addData(indexNeighbor, weightNeighbor, forceFromNeighbor, t);
    }
}

void Bias::updateBiasStateSharedCorrelationTensorTimeIntegral()
{
    state_.updateSharedCorrelationTensorTimeIntegral(params_, *forceCorrelationGrid_, false);
}

/* Return the number of data blocks that have been prepared for writing. */
int Bias::numEnergySubblocksToWrite() const
{
    GMX_RELEASE_ASSERT(writer_ != nullptr, "Should only request data from an initialized writer");

    return writer_->numBlocks();
}

/* Write bias data blocks to energy subblocks. */
int Bias::writeToEnergySubblocks(t_enxsubblock* subblock) const
{
    GMX_RELEASE_ASSERT(writer_ != nullptr, "Should only request data from an initialized writer");

    return writer_->writeToEnergySubblocks(*this, subblock);
}

bool Bias::isSampleCoordStep(const int64_t step) const
{
    return params_.isSampleCoordStep(step);
}


} // namespace gmx
