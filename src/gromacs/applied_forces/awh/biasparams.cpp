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
 * Implements the initialization of the BiasParams class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "biasparams.h"

#include <cmath>
#include <cstddef>

#include <algorithm>
#include <limits>

#include "gromacs/applied_forces/awh/dimparams.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "biasgrid.h"

namespace gmx
{

namespace
{

/*! \brief Determines the interval for updating the target distribution.
 *
 * The interval value is based on the target distrbution type
 * (this could be made a user-option but there is most likely
 * no big need for tweaking this for most users).
 *
 * \param[in] awhParams      AWH parameters.
 * \param[in] awhBiasParams  Bias parameters.
 * \returns the target update interval in steps.
 */
int64_t calcTargetUpdateInterval(const AwhParams& awhParams, const AwhBiasParams& awhBiasParams)
{
    int64_t numStepsUpdateTarget = 0;
    /* Set the target update frequency based on the target distrbution type
     * (this could be made a user-option but there is most likely no big need
     * for tweaking this for most users).
     */
    switch (awhBiasParams.targetDistribution())
    {
        case AwhTargetType::Constant:
            /* If the target distribution should be scaled by the friction metric,
             * set numStepsUpdateTarget below */
            if (!awhBiasParams.scaleTargetByMetric())
            {
                numStepsUpdateTarget = 0;
                break;
            }
            [[fallthrough]];
        case AwhTargetType::Cutoff:
        case AwhTargetType::Boltzmann:
            /* Updating the target generally requires updating the whole grid so to keep the cost
               down we generally update the target less often than the free energy (unless the free
               energy update step is set to > 100 samples). */
            numStepsUpdateTarget = std::max(100 % awhParams.numSamplesUpdateFreeEnergy(),
                                            awhParams.numSamplesUpdateFreeEnergy())
                                   * awhParams.nstSampleCoord();
            break;
        case AwhTargetType::LocalBoltzmann:
            /* The target distribution is set equal to the reference histogram which is updated every free energy update.
               So the target has to be updated at the same time. This leads to a global update each time because it is
               assumed that a target distribution update can take any form. This is a bit unfortunate for a "local"
               target distribution. One could avoid the global update by making a local target update function (and
               postponing target updates for non-local points as for the free energy update). We avoid such additions
               for now and accept that this target type always does global updates. */
            numStepsUpdateTarget = awhParams.numSamplesUpdateFreeEnergy() * awhParams.nstSampleCoord();
            break;
        default: GMX_RELEASE_ASSERT(false, "Unknown AWH target type"); break;
    }

    return numStepsUpdateTarget;
}

/*! \brief Determines the step interval for checking for covering.
 *
 * \param[in] awhParams      AWH parameters.
 * \param[in] dimParams         Parameters for the dimensions of the coordinate.
 * \param[in] gridAxis          The BiasGrid axes.
 * \returns the check interval in steps.
 */
int64_t calcCheckCoveringInterval(const AwhParams&          awhParams,
                                  ArrayRef<const DimParams> dimParams,
                                  ArrayRef<const GridAxis>  gridAxis)
{
    /* Each sample will have a width of sigma. To cover the axis a
       minimum number of samples of width sigma is required. */
    int minNumSamplesCover = 0;
    for (size_t d = 0; d < gridAxis.size(); d++)
    {
        int numSamplesCover;
        if (dimParams[d].isPullDimension())
        {
            GMX_RELEASE_ASSERT(
                    dimParams[d].pullDimParams().betak > 0,
                    "Inverse temperature (beta) and force constant (k) should be positive.");
            double sigma = 1.0 / std::sqrt(dimParams[d].pullDimParams().betak);

            /* The additional sample is here because to cover a discretized
            axis of length sigma one needs two samples, one for each
            end point. */
            GMX_RELEASE_ASSERT(gridAxis[d].length() / sigma < std::numeric_limits<int>::max(),
                               "The axis length in units of sigma should fit in an int");
            numSamplesCover = static_cast<int>(std::ceil(gridAxis[d].length() / sigma)) + 1;
        }
        else
        {
            numSamplesCover = dimParams[d].fepDimParams().numFepLambdaStates;
        }
        /* The minimum number of samples needed for simultaneously
        covering all axes is limited by the axis requiring most
        samples. */
        minNumSamplesCover = std::max(minNumSamplesCover, numSamplesCover);
    }

    /* Convert to number of steps using the sampling frequency. The
       check interval should be a multiple of the update step
       interval. */
    int numStepsUpdate = awhParams.numSamplesUpdateFreeEnergy() * awhParams.nstSampleCoord();
    GMX_RELEASE_ASSERT(awhParams.numSamplesUpdateFreeEnergy() > 0,
                       "When checking for AWH coverings, the number of samples per AWH update need "
                       "to be > 0.");
    int numUpdatesCheck = std::max(1, minNumSamplesCover / awhParams.numSamplesUpdateFreeEnergy());
    int numStepsCheck   = numUpdatesCheck * numStepsUpdate;

    GMX_RELEASE_ASSERT(numStepsCheck % numStepsUpdate == 0,
                       "Only check covering at free energy update steps");

    return numStepsCheck;
}

/*! \brief
 * Estimate a reasonable initial reference weight histogram size.
 *
 * \param[in] awhBiasParams     Bias parameters.
 * \param[in] gridAxis          The BiasGrid axes.
 * \param[in] beta              1/(k_B T).
 * \param[in] samplingTimestep  Sampling frequency of probability weights.
 * \returns estimate of initial histogram size.
 */
double getInitialHistogramSizeEstimate(const AwhBiasParams&     awhBiasParams,
                                       ArrayRef<const GridAxis> gridAxis,
                                       double                   beta,
                                       double                   samplingTimestep)
{
    /* Get diffusion factor */
    double              maxCrossingTime = 0.;
    std::vector<double> x;
    const auto          awhDimParams = awhBiasParams.dimParams();
    for (size_t d = 0; d < gridAxis.size(); d++)
    {
        GMX_RELEASE_ASSERT(awhDimParams[d].diffusion() > 0, "We need positive diffusion");
        // With diffusion it takes on average T = L^2/2D time to cross length L
        double axisLength   = gridAxis[d].isFepLambdaAxis() ? 1.0 : gridAxis[d].length();
        double crossingTime = (axisLength * axisLength) / (2 * awhDimParams[d].diffusion());
        maxCrossingTime     = std::max(maxCrossingTime, crossingTime);
    }
    GMX_RELEASE_ASSERT(maxCrossingTime > 0, "We need at least one dimension with non-zero length");
    double errorInitialInKT = beta * awhBiasParams.initialErrorEstimate();
    double histogramSize    = maxCrossingTime / (gmx::square(errorInitialInKT) * samplingTimestep);

    return histogramSize;
}

/*! \brief Returns the number of simulations sharing bias updates.
 *
 * \param[in] awhBiasParams          Bias parameters.
 * \param[in] numSharingSimulations  The number of simulations to share the bias across.
 * \returns the number of shared updates.
 */
int getNumSharedUpdate(const AwhBiasParams& awhBiasParams, int numSharingSimulations)
{
    GMX_RELEASE_ASSERT(numSharingSimulations >= 1, "We should ''share'' at least with ourselves");

    int numShared = 1;

    if (awhBiasParams.shareGroup() > 0)
    {
        /* We do not yet support sharing within a simulation */
        int numSharedWithinThisSimulation = 1;
        numShared                         = numSharingSimulations * numSharedWithinThisSimulation;
    }

    return numShared;
}

} // namespace

BiasParams::BiasParams(const AwhParams&          awhParams,
                       const AwhBiasParams&      awhBiasParams,
                       ArrayRef<const DimParams> dimParams,
                       double                    beta,
                       double                    mdTimeStep,
                       DisableUpdateSkips        disableUpdateSkips,
                       int                       numSharingSimulations,
                       ArrayRef<const GridAxis>  gridAxis,
                       int                       biasIndex) :
    invBeta(beta > 0 ? 1 / beta : 0),
    numStepsSampleCoord_(awhParams.nstSampleCoord()),
    numSamplesUpdateFreeEnergy_(awhParams.numSamplesUpdateFreeEnergy()),
    numStepsUpdateTarget_(calcTargetUpdateInterval(awhParams, awhBiasParams)),
    numStepsCheckCovering_(calcCheckCoveringInterval(awhParams, dimParams, gridAxis)),
    eTarget(awhBiasParams.targetDistribution()),
    scaleTargetByMetric(awhBiasParams.scaleTargetByMetric()),
    targetMetricScalingLimit(awhBiasParams.targetMetricScalingLimit()),
    freeEnergyCutoffInKT(beta * awhBiasParams.targetCutoff()),
    temperatureScaleFactor(awhBiasParams.targetBetaScaling()),
    idealWeighthistUpdate(eTarget != AwhTargetType::LocalBoltzmann),
    numSharedUpdate(getNumSharedUpdate(awhBiasParams, numSharingSimulations)),
    updateWeight(numSamplesUpdateFreeEnergy_ * numSharedUpdate),
    localWeightScaling(eTarget == AwhTargetType::LocalBoltzmann ? temperatureScaleFactor : 1),
    initialErrorInKT(beta * awhBiasParams.initialErrorEstimate()),
    initialHistogramSize(
            getInitialHistogramSizeEstimate(awhBiasParams, gridAxis, beta, numStepsSampleCoord_ * mdTimeStep)),
    convolveForce(awhParams.potential() == AwhPotentialType::Convolved),
    biasIndex_(biasIndex),
    disableUpdateSkips_(disableUpdateSkips == DisableUpdateSkips::yes)
{
    if (beta <= 0)
    {
        GMX_THROW(InvalidInputError("To use AWH, the beta=1/(k_B T) should be > 0"));
    }

    const auto& awhDimParams = awhBiasParams.dimParams();
    for (int d = 0; d < gmx::ssize(awhDimParams); d++)
    {
        /* The spacing in FEP dimensions is 1. The calculated coverRadius will be in lambda states
         * (cf points in other dimensions). */
        double coverRadiusInNm =
                0.5 * dimParams[d].scaleUserInputToInternal(awhDimParams[d].coverDiameter());
        double spacing  = gridAxis[d].spacing();
        coverRadius_[d] = spacing > 0 ? static_cast<int>(std::round(coverRadiusInNm / spacing)) : 0;
    }
}

} // namespace gmx
