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
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "grid.h"
#include "math.h"
#include "pointstate.h"

namespace gmx
{

void Bias::checkHistograms(double t, gmx_int64_t step, FILE *fplog)
{
    const int    maxNumWarningsInCheck = 1;   /* The maximum number of warnings to print per check */
    const int    maxNumWarningsInRun   = 10;  /* The maximum number of warnings to print in a run */

    if (fplog == nullptr || numWarningsIssued_ >= maxNumWarningsInRun || state_.inInitialStage() ||
        !params_.isCheckStep(state_.points().size(), step))
    {
        return;
    }

    numWarningsIssued_ +=
        state_.checkHistograms(grid(), biasIndex(), t, fplog,
                               maxNumWarningsInCheck);

    if (numWarningsIssued_ >= maxNumWarningsInRun)
    {
        fprintf(fplog, "\nawh%d: suppressing future AWH warnings.\n", biasIndex() + 1);
    }
}

void Bias::doSkippedUpdatesForAllPoints()
{
    state_.doSkippedUpdatesForAllPoints(params_);
}

void Bias::calcForceAndUpdateBias(const awh_dvec coordValue,
                                  awh_dvec biasForce,
                                  double *awhPotential, double *potentialJump,
                                  const gmx_multisim_t *ms,
                                  double t, gmx_int64_t step, gmx_int64_t seed,
                                  FILE *fplog)
{
    if (step < 0)
    {
        gmx_fatal(FARGS, "The step number is negative which is not supported by the AWH code.");
    }

    state_.setCoordValue(grid(), coordValue);

    std::vector<double> *probWeightNeighbor = &tempWorkSpace_;

    /* If the convolved force is needed or this is a sampling step,
     * the bias in the current neighborhood needs to be up-to-date
     * and the probablity weights need to be calculated.
     */
    const bool sampleCoord   = params_.isSampleCoordStep(step);
    double     convolvedBias = 0;
    if (params_.convolveForce || sampleCoord)
    {
        if (params_.skipUpdates())
        {
            state_.doSkippedUpdatesInNeighborhood(params_, grid());
        }

        convolvedBias = state_.updateProbabilityWeightsAndConvolvedBias(dimParams_, grid(), probWeightNeighbor);

        if (step > 0 && sampleCoord)
        {
            state_.sampleCoordAndPmf(grid(), *probWeightNeighbor, convolvedBias);
        }
    }

    const CoordinateState &coordinateState = state_.coordinateState();

    /* Set the bias force and get the potential contribution from this bias.
     * The potential jump occurs at different times depending on how
     * the force is applied (and how the potential is normalized).
     * For the convolved force it happens when the bias is updated,
     * for the umbrella when the umbrella is moved.
     */
    double potential, newPotential;
    if (params_.convolveForce)
    {
        state_.calcConvolvedForce(dimParams_, grid(), *probWeightNeighbor,
                                  biasForce);

        potential    = -convolvedBias*params_.invBeta;
        newPotential = potential; /* Assume no jump */
    }
    else
    {
        /* Umbrella force */
        GMX_RELEASE_ASSERT(state_.points()[coordinateState.umbrellaGridpoint()].inTargetRegion(),
                           "AWH bias grid point for the umbrella reference value is outside of the target region.");
        potential =
            state_.calcUmbrellaForceAndPotential(dimParams_, grid(), coordinateState.umbrellaGridpoint(), biasForce);

        /* Moving the umbrella results in a force correction and
         * a new potential. The umbrella center is sampled as often as
         * the coordinate so we know the probability weights needed
         * for moving the umbrella are up-to-date.
         */
        if (sampleCoord)
        {
            newPotential = state_.moveUmbrella(dimParams_, grid(), *probWeightNeighbor, biasForce, step, seed, params_.biasIndex);
        }
        else
        {
            newPotential = potential;
        }
    }

    /* Update the free energy estimates and bias and other history dependent method parameters */
    if (params_.isUpdateFreeEnergyStep(step))
    {
        state_.updateFreeEnergyAndAddSamplesToHistogram(dimParams_, grid(),
                                                        params_,
                                                        ms, t, step, fplog,
                                                        &updateList_);

        if (params_.convolveForce)
        {
            /* The update results in a potential jump, so we need the new convolved potential. */
            newPotential = -calcConvolvedBias(dimParams_, grid(), state_.points(), coordinateState.coordValue())*params_.invBeta;
        }
    }

    /* Return the new potential. */
    *awhPotential  = potential;
    /* Return the potential jump of this bias. */
    *potentialJump = (newPotential - potential);

    /* Check the sampled histograms and potentially warn user if something is suspicious */
    checkHistograms(t, step, fplog);
}

void Bias::restoreStateFromHistory(const AwhBiasHistory *biasHistory,
                                   const t_commrec      *cr)
{
    if (MASTER(cr))
    {
        GMX_RELEASE_ASSERT(biasHistory != nullptr, "On the master rank we need a valid history object to restore from");
        state_.restoreFromHistory(*biasHistory, grid());
    }

    if (PAR(cr))
    {
        state_.broadcast(cr);
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
           BiasParams::DisableUpdateSkips  disableUpdateSkips) :
    dimParams_(dimParamsInit),
    grid_(new Grid(dimParamsInit, awhBiasParams.dimParams)),
    params_(awhParams, awhBiasParams, dimParams_, beta, mdTimeStep, disableUpdateSkips, numSharingSimulations, grid_->axis(), biasIndexInCollection),
    state_(awhBiasParams, params_.histSizeInitial, dimParams_, grid()),
    tempWorkSpace_(),
    numWarningsIssued_(0)
{
    /* For a global update updateList covers all points, so reserve that */
    updateList_.reserve(grid_->numPoints());

    state_.initGridPointState(awhBiasParams, dimParams_, grid(), params_, biasInitFilename, awhParams.numBias);
}

} // namespace gmx
