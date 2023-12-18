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
 * Implements the BiasState class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "biasstate.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <optional>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "biasgrid.h"
#include "biassharing.h"
#include "correlationgrid.h"
#include "correlationtensor.h"
#include "dimparams.h"
#include "pointstate.h"

namespace gmx
{

void BiasState::getPmf(gmx::ArrayRef<float> pmf) const
{
    GMX_ASSERT(pmf.size() == points_.size(), "pmf should have the size of the bias grid");

    /* The PMF is just the negative of the log of the sampled PMF histogram.
     * Points with zero target weight are ignored, they will mostly contain noise.
     */
    for (size_t i = 0; i < points_.size(); i++)
    {
        pmf[i] = points_[i].inTargetRegion() ? -points_[i].logPmfSum() : GMX_FLOAT_MAX;
    }
}

namespace
{

/*! \brief
 * Sum PMF over multiple simulations, when requested.
 *
 * \param[in,out] pointState         The state of the points in the bias.
 * \param[in]     numSharedUpdate    The number of biases sharing the histogram.
 * \param[in]     biasSharing        Object for sharing bias data over multiple simulations
 * \param[in]     biasIndex          Index of this bias in the total list of biases in this simulation
 */
void sumPmf(gmx::ArrayRef<PointState> pointState, int numSharedUpdate, const BiasSharing* biasSharing, const int biasIndex)
{
    if (numSharedUpdate == 1)
    {
        return;
    }
    GMX_ASSERT(biasSharing != nullptr
                       && numSharedUpdate % biasSharing->numSharingSimulations(biasIndex) == 0,
               "numSharedUpdate should be a multiple of multiSimComm->numSimulations_");
    GMX_ASSERT(numSharedUpdate == biasSharing->numSharingSimulations(biasIndex),
               "Sharing within a simulation is not implemented (yet)");

    std::vector<double> buffer(pointState.size());

    /* Need to temporarily exponentiate the log weights to sum over simulations */
    for (size_t i = 0; i < buffer.size(); i++)
    {
        buffer[i] = pointState[i].inTargetRegion() ? std::exp(pointState[i].logPmfSum()) : 0;
    }

    biasSharing->sumOverSharingSimulations(gmx::ArrayRef<double>(buffer), biasIndex);

    /* Take log again to get (non-normalized) PMF */
    double normFac = 1.0 / numSharedUpdate;
    for (gmx::Index i = 0; i < pointState.ssize(); i++)
    {
        if (pointState[i].inTargetRegion())
        {
            pointState[i].setLogPmfSum(std::log(buffer[i] * normFac));
        }
    }
}

/*! \brief
 * Find the minimum free energy value.
 *
 * \param[in] pointState  The state of the points.
 * \returns the minimum free energy value.
 */
double freeEnergyMinimumValue(gmx::ArrayRef<const PointState> pointState)
{
    double fMin = GMX_FLOAT_MAX;

    for (auto const& ps : pointState)
    {
        if (ps.inTargetRegion() && ps.freeEnergy() < fMin)
        {
            fMin = ps.freeEnergy();
        }
    }

    return fMin;
}

/*! \brief
 * Find and return the log of the probability weight of a point given a coordinate value.
 *
 * The unnormalized weight is given by
 * w(point|value) = exp(bias(point) - U(value,point)),
 * where U is a harmonic umbrella potential.
 *
 * \param[in] dimParams              The bias dimensions parameters
 * \param[in] points                 The point state.
 * \param[in] grid                   The grid.
 * \param[in] pointIndex             Point to evaluate probability weight for.
 * \param[in] pointBias              Bias for the point (as a log weight).
 * \param[in] value                  Coordinate value.
 * \param[in] neighborLambdaEnergies The energy of the system in neighboring lambdas states. Can be
 * empty when there are no free energy lambda state dimensions.
 * \param[in] gridpointIndex         The index of the current grid point.
 * \returns the log of the biased probability weight.
 */
double biasedLogWeightFromPoint(ArrayRef<const DimParams>  dimParams,
                                ArrayRef<const PointState> points,
                                const BiasGrid&            grid,
                                int                        pointIndex,
                                double                     pointBias,
                                const awh_dvec             value,
                                ArrayRef<const double>     neighborLambdaEnergies,
                                int                        gridpointIndex)
{
    double logWeight = detail::c_largeNegativeExponent;

    /* Only points in the target region have non-zero weight */
    if (points[pointIndex].inTargetRegion())
    {
        logWeight = pointBias;

        /* Add potential for all parameter dimensions */
        for (size_t d = 0; d < dimParams.size(); d++)
        {
            if (dimParams[d].isFepLambdaDimension())
            {
                /* If this is not a sampling step or if this function is called from
                 * calcConvolvedBias(), when writing energy subblocks, neighborLambdaEnergies will
                 * be empty. No convolution is required along the lambda dimension. */
                if (!neighborLambdaEnergies.empty())
                {
                    const int pointLambdaIndex     = grid.point(pointIndex).coordValue[d];
                    const int gridpointLambdaIndex = grid.point(gridpointIndex).coordValue[d];

                    const double energyDiff = neighborLambdaEnergies[pointLambdaIndex]
                                              - neighborLambdaEnergies[gridpointLambdaIndex];
                    if (dimParams[d].fepDimParams().beta * energyDiff < -0.5 * detail::c_largePositiveExponent)
                    {
                        GMX_THROW(SimulationInstabilityError(gmx::formatString(
                                "AWH lambda dimension encountered a too large negative neighbor "
                                "energy difference of %f kJ/mol",
                                energyDiff)));
                    }
                    logWeight -= dimParams[d].fepDimParams().beta * energyDiff;
                }
            }
            else
            {
                double dev = getDeviationFromPointAlongGridAxis(grid, d, pointIndex, value[d]);
                logWeight -= 0.5 * dimParams[d].pullDimParams().betak * dev * dev;
            }
        }
    }
    return logWeight;
}

/*! \brief
 * Calculates the marginal distribution (marginal probability) for each value along
 * a free energy lambda axis.
 * The marginal distribution of one coordinate dimension value is the sum of the probability
 * distribution of all values (herein all neighbor values) with the same value in the dimension
 * of interest.
 * \param[in] grid               The bias grid.
 * \param[in] neighbors          The points to use for the calculation of the marginal distribution.
 * \param[in] probWeightNeighbor Probability weights of the neighbors.
 * \returns The calculated marginal distribution in a 1D array with
 * as many elements as there are points along the axis of interest.
 */
std::vector<double> calculateFELambdaMarginalDistribution(const BiasGrid&        grid,
                                                          ArrayRef<const int>    neighbors,
                                                          ArrayRef<const double> probWeightNeighbor)
{
    const std::optional<int> lambdaAxisIndex = grid.lambdaAxisIndex();
    GMX_RELEASE_ASSERT(lambdaAxisIndex,
                       "There must be a free energy lambda axis in order to calculate the free "
                       "energy lambda marginal distribution.");
    const int           numFepLambdaStates = grid.numFepLambdaStates();
    std::vector<double> lambdaMarginalDistribution(numFepLambdaStates, 0);

    for (size_t i = 0; i < neighbors.size(); i++)
    {
        const int neighbor    = neighbors[i];
        const int lambdaState = grid.point(neighbor).coordValue[lambdaAxisIndex.value()];
        lambdaMarginalDistribution[lambdaState] += probWeightNeighbor[i];
    }
    return lambdaMarginalDistribution;
}

} // namespace

void BiasState::calcConvolvedPmf(ArrayRef<const DimParams> dimParams,
                                 const BiasGrid&           grid,
                                 std::vector<float>*       convolvedPmf) const
{
    size_t numPoints = grid.numPoints();

    convolvedPmf->resize(numPoints);

    /* Get the PMF to convolve. */
    std::vector<float> pmf(numPoints);
    getPmf(pmf);

    for (size_t m = 0; m < numPoints; m++)
    {
        double           freeEnergyWeights = 0;
        const GridPoint& point             = grid.point(m);
        for (const auto& neighbor : point.neighbor)
        {
            /* Do not convolve the bias along a lambda axis - only use the pmf from the current point */
            if (!pointsHaveDifferentLambda(grid, m, neighbor))
            {
                /* The negative PMF is a positive bias. */
                double biasNeighbor = -pmf[neighbor];

                /* Add the convolved PMF weights for the neighbors of this point.
                Note that this function only adds point within the target > 0 region.
                Sum weights, take the logarithm last to get the free energy. */
                double logWeight = biasedLogWeightFromPoint(
                        dimParams, points_, grid, neighbor, biasNeighbor, point.coordValue, {}, m);
                freeEnergyWeights += std::exp(logWeight);
            }
        }

        GMX_RELEASE_ASSERT(freeEnergyWeights > 0,
                           "Attempting to do log(<= 0) in AWH convolved PMF calculation.");
        // We should cast to float after taking the logarithm to avoid underflows
        (*convolvedPmf)[m] = static_cast<float>(-std::log(freeEnergyWeights));
    }
}

double BiasState::calculateAverageNonZeroMetric()
{
    int    elementCount = 0;
    double sumVolume    = 0;
    for (gmx::Index pointIndex = 0; pointIndex < ssize(points_); pointIndex++)
    {
        std::vector<double> correlationIntegral     = getSharedPointCorrelationIntegral(pointIndex);
        double              correlationTensorVolume = getSqrtDeterminant(correlationIntegral);

        if (correlationTensorVolume > 0)
        {
            sumVolume += correlationTensorVolume;
            elementCount++;
        }
    }
    double averageVolume = 0;
    if (elementCount != 0)
    {
        averageVolume = sumVolume / elementCount;
    }
    return averageVolume;
}

double BiasState::scaleTargetByMetric(double targetMetricScalingLimit)
{
    GMX_RELEASE_ASSERT(targetMetricScalingLimit > 1,
                       "When scaling by friction metric, the upper scaling limit must be > 1.");
    /* Calculate the average of non-zero correlation tensor volume elements before
     * scaling by the friction tensor. The average will be used to scale the target
     * relatively and to avoid scaling (scaleFactor = 1) where there is not enough
     * data (no valid correlationTensorVolume) to use for scaling. */
    double averageVolume = calculateAverageNonZeroMetric();
    if (averageVolume == 0)
    {
        averageVolume = 1;
    }
    double sumTarget = 0;

    /* There is a global scaling factor limit. The lower limit can be modified per
     * point depending on its sampling, but it cannot be lower than lowerScalingLimit. */
    double upperScalingLimit = averageVolume * targetMetricScalingLimit;
    double lowerScalingLimit = averageVolume / targetMetricScalingLimit;
    for (gmx::Index pointIndex = 0; pointIndex < ssize(points_); pointIndex++)
    {
        PointState& ps = points_[pointIndex];

        std::vector<double> correlationIntegral     = getSharedPointCorrelationIntegral(pointIndex);
        double              correlationTensorVolume = getSqrtDeterminant(correlationIntegral);

        /* Points may have a very low correlation tensor from not being sampled enough.
         * This sets a lower limit on the scaling based on the amount of samples */
        const double weightSumTot     = points_[pointIndex].weightSumTot();
        double pointLowerScalingLimit = weightSumTot > 1 ? averageVolume / weightSumTot : averageVolume;
        pointLowerScalingLimit        = std::max(lowerScalingLimit, pointLowerScalingLimit);

        /* If there is no correlation tensor volume from this point use the average
         * volume. This will result in no friction tensor scaling for the target
         * distribution of this point. */
        if (correlationTensorVolume == 0)
        {
            correlationTensorVolume = averageVolume;
        }
        /* Clamp the scaling to the allowed interval */
        else
        {
            correlationTensorVolume =
                    std::clamp(correlationTensorVolume, pointLowerScalingLimit, upperScalingLimit);
        }
        double scaleFactor = correlationTensorVolume / averageVolume;
        ps.scaleTarget(scaleFactor);

        sumTarget += ps.target();
    }
    return sumTarget;
}

namespace
{
/*! \brief
 * Puts together a string describing a grid point.
 *
 * \param[in] grid         The grid.
 * \param[in] point        BiasGrid point index.
 * \returns a string for the point.
 */
std::string gridPointValueString(const BiasGrid& grid, int point)
{
    std::string pointString;

    pointString += "(";

    for (int d = 0; d < grid.numDimensions(); d++)
    {
        pointString += gmx::formatString("%g", grid.point(point).coordValue[d]);
        if (d < grid.numDimensions() - 1)
        {
            pointString += ",";
        }
        else
        {
            pointString += ")";
        }
    }

    return pointString;
}

} // namespace

void BiasState::updateTargetDistribution(const BiasParams& params, const CorrelationGrid& forceCorrelation)
{
    double freeEnergyCutoff = 0;
    if (params.eTarget == AwhTargetType::Cutoff)
    {
        freeEnergyCutoff = freeEnergyMinimumValue(points_) + params.freeEnergyCutoffInKT;
    }

    double sumTarget = 0;
    for (PointState& ps : points_)
    {
        sumTarget += ps.updateTargetWeight(params, freeEnergyCutoff);
    }
    GMX_RELEASE_ASSERT(sumTarget > 0, "We should have a non-zero distribution");

    /* Scale the target distribution by the friction metric - normalize afterwards */
    if (params.scaleTargetByMetric && !inInitialStage())
    {
        updateSharedCorrelationTensorTimeIntegral(params, forceCorrelation, true);
        sumTarget = scaleTargetByMetric(params.targetMetricScalingLimit);
    }

    /* Normalize to 1 */
    double invSum = 1.0 / sumTarget;
    for (PointState& ps : points_)
    {
        ps.scaleTarget(invSum);
    }
}

int BiasState::warnForHistogramAnomalies(const BiasGrid& grid, int biasIndex, double t, FILE* fplog, int maxNumWarnings) const
{
    GMX_ASSERT(fplog != nullptr, "Warnings can only be issued if there is log file.");
    const double maxHistogramRatio = 0.5; /* Tolerance for printing a warning about the histogram ratios */

    /* Sum up the histograms and get their normalization */
    double sumVisits  = 0;
    double sumWeights = 0;
    for (const auto& pointState : points_)
    {
        if (pointState.inTargetRegion())
        {
            sumVisits += pointState.numVisitsTot();
            sumWeights += pointState.weightSumTot();
        }
    }
    GMX_RELEASE_ASSERT(sumVisits > 0, "We should have visits");
    GMX_RELEASE_ASSERT(sumWeights > 0, "We should have weight");
    double invNormVisits = 1.0 / sumVisits;
    double invNormWeight = 1.0 / sumWeights;

    /* Check all points for warnings */
    int    numWarnings = 0;
    size_t numPoints   = grid.numPoints();
    for (size_t m = 0; m < numPoints; m++)
    {
        /* Skip points close to boundary or non-target region */
        const GridPoint& gridPoint = grid.point(m);
        bool             skipPoint = false;
        for (size_t n = 0; (n < gridPoint.neighbor.size()) && !skipPoint; n++)
        {
            int neighbor = gridPoint.neighbor[n];
            skipPoint    = !points_[neighbor].inTargetRegion();
            for (int d = 0; (d < grid.numDimensions()) && !skipPoint; d++)
            {
                const GridPoint& neighborPoint = grid.point(neighbor);
                skipPoint                      = neighborPoint.index[d] == 0
                            || neighborPoint.index[d] == grid.axis(d).numPoints() - 1;
            }
        }

        /* Warn if the coordinate distribution is less than the target distribution with a certain fraction somewhere */
        const double relativeWeight = points_[m].weightSumTot() * invNormWeight;
        const double relativeVisits = points_[m].numVisitsTot() * invNormVisits;
        if (!skipPoint && relativeVisits < relativeWeight * maxHistogramRatio)
        {
            std::string pointValueString = gridPointValueString(grid, m);
            std::string warningMessage   = gmx::formatString(
                    "\nawh%d warning: "
                    "at t = %g ps the obtained coordinate distribution at coordinate value %s "
                    "is less than a fraction %g of the reference distribution at that point. "
                    "If you are not certain about your settings you might want to increase your "
                    "pull force constant or "
                    "modify your sampling region.\n",
                    biasIndex + 1,
                    t,
                    pointValueString.c_str(),
                    maxHistogramRatio);
            gmx::TextLineWrapper wrapper;
            wrapper.settings().setLineLength(c_linewidth);
            fprintf(fplog, "%s", wrapper.wrapToString(warningMessage).c_str());

            numWarnings++;
        }
        if (numWarnings >= maxNumWarnings)
        {
            break;
        }
    }

    return numWarnings;
}

double BiasState::calcUmbrellaForceAndPotential(ArrayRef<const DimParams> dimParams,
                                                const BiasGrid&           grid,
                                                int                       point,
                                                ArrayRef<const double>    neighborLambdaDhdl,
                                                ArrayRef<double>          force) const
{
    double potential = 0;
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        if (dimParams[d].isFepLambdaDimension())
        {
            /* The force we set here is only used for computing the friction metric */
            if (!neighborLambdaDhdl.empty())
            {
                const int coordpointLambdaIndex = grid.point(point).coordValue[d];
                force[d]                        = neighborLambdaDhdl[coordpointLambdaIndex];
                /* The potential should not be affected by the lambda dimension. */
            }
        }
        else
        {
            double deviation =
                    getDeviationFromPointAlongGridAxis(grid, d, point, coordState_.coordValue()[d]);
            double k = dimParams[d].pullDimParams().k;

            /* Force from harmonic potential 0.5*k*dev^2 */
            force[d] = -k * deviation;
            potential += 0.5 * k * deviation * deviation;
        }
    }

    return potential;
}

void BiasState::calcConvolvedForce(ArrayRef<const DimParams> dimParams,
                                   const BiasGrid&           grid,
                                   ArrayRef<const double>    probWeightNeighbor,
                                   ArrayRef<const double>    neighborLambdaDhdl,
                                   ArrayRef<double>          forceWorkBuffer,
                                   ArrayRef<double>          force) const
{
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        force[d] = 0;
    }

    /* Only neighboring points have non-negligible contribution. */
    const std::vector<int>& neighbor          = grid.point(coordState_.gridpointIndex()).neighbor;
    gmx::ArrayRef<double>   forceFromNeighbor = forceWorkBuffer;
    for (size_t n = 0; n < neighbor.size(); n++)
    {
        double weightNeighbor = probWeightNeighbor[n];
        int    indexNeighbor  = neighbor[n];

        /* Get the umbrella force from this point. The returned potential is ignored here. */
        calcUmbrellaForceAndPotential(dimParams, grid, indexNeighbor, neighborLambdaDhdl, forceFromNeighbor);

        /* Add the weighted umbrella force to the convolved force. */
        for (size_t d = 0; d < dimParams.size(); d++)
        {
            force[d] += forceFromNeighbor[d] * weightNeighbor;
        }
    }
}

double BiasState::moveUmbrella(ArrayRef<const DimParams> dimParams,
                               const BiasGrid&           grid,
                               ArrayRef<const double>    probWeightNeighbor,
                               ArrayRef<const double>    neighborLambdaDhdl,
                               ArrayRef<double>          biasForce,
                               int64_t                   step,
                               int64_t                   seed,
                               int                       indexSeed,
                               bool                      onlySampleUmbrellaGridpoint)
{
    /* Generate and set a new coordinate reference value */
    coordState_.sampleUmbrellaGridpoint(
            grid, coordState_.gridpointIndex(), probWeightNeighbor, step, seed, indexSeed);

    if (onlySampleUmbrellaGridpoint)
    {
        return 0;
    }

    std::vector<double> newForce(dimParams.size());
    double              newPotential = calcUmbrellaForceAndPotential(
            dimParams, grid, coordState_.umbrellaGridpoint(), neighborLambdaDhdl, newForce);

    /*  A modification of the reference value at time t will lead to a different
        force over t-dt/2 to t and over t to t+dt/2. For high switching rates
        this means the force and velocity will change signs roughly as often.
        To avoid any issues we take the average of the previous and new force
        at steps when the reference value has been moved. E.g. if the ref. value
        is set every step to (coord dvalue +/- delta) would give zero force.
     */
    for (gmx::Index d = 0; d < biasForce.ssize(); d++)
    {
        /* Average of the current and new force */
        biasForce[d] = 0.5 * (biasForce[d] + newForce[d]);
    }

    return newPotential;
}

namespace
{

/*! \brief
 * Sets the histogram rescaling factors needed to control the histogram size.
 *
 * For sake of robustness, the reference weight histogram can grow at a rate
 * different from the actual sampling rate. Typically this happens for a limited
 * initial time, alternatively growth is scaled down by a constant factor for all
 * times. Since the size of the reference histogram sets the size of the free
 * energy update this should be reflected also in the PMF. Thus the PMF histogram
 * needs to be rescaled too.
 *
 * This function should only be called by the bias update function or wrapped by a function that
 * knows what scale factors should be applied when, e.g,
 * getSkippedUpdateHistogramScaleFactors().
 *
 * \param[in]  params             The bias parameters.
 * \param[in]  newHistogramSize   New reference weight histogram size.
 * \param[in]  oldHistogramSize   Previous reference weight histogram size (before adding new samples).
 * \param[out] weightHistScaling  Scaling factor for the reference weight histogram.
 * \param[out] logPmfSumScaling   Log of the scaling factor for the PMF histogram.
 */
void setHistogramUpdateScaleFactors(const BiasParams& params,
                                    double            newHistogramSize,
                                    double            oldHistogramSize,
                                    double*           weightHistScaling,
                                    double*           logPmfSumScaling)
{

    /* The two scaling factors below are slightly different (ignoring the log factor) because the
       reference and the PMF histogram apply weight scaling differently. The weight histogram
       applies is  locally, i.e. each sample is scaled down meaning all samples get equal weight.
       It is done this way because that is what target type local Boltzmann (for which
       target = weight histogram) needs. In contrast, the PMF histogram is rescaled globally
       by repeatedly scaling down the whole histogram. The reasons for doing it this way are:
       1) empirically this is necessary for converging the PMF; 2) since the extraction of
       the PMF is theoretically only valid for a constant bias, new samples should get more
       weight than old ones for which the bias is fluctuating more. */
    *weightHistScaling =
            newHistogramSize / (oldHistogramSize + params.updateWeight * params.localWeightScaling);
    *logPmfSumScaling = std::log(newHistogramSize / (oldHistogramSize + params.updateWeight));
}

} // namespace

void BiasState::getSkippedUpdateHistogramScaleFactors(const BiasParams& params,
                                                      double*           weightHistScaling,
                                                      double*           logPmfSumScaling) const
{
    GMX_ASSERT(params.skipUpdates(),
               "Calling function for skipped updates when skipping updates is not allowed");

    if (inInitialStage())
    {
        /* In between global updates the reference histogram size is kept constant so we trivially
           know what the histogram size was at the time of the skipped update. */
        double histogramSize = histogramSize_.histogramSize();
        setHistogramUpdateScaleFactors(
                params, histogramSize, histogramSize, weightHistScaling, logPmfSumScaling);
    }
    else
    {
        /* In the final stage, the reference histogram grows at the sampling rate which gives trivial scale factors. */
        *weightHistScaling = 1;
        *logPmfSumScaling  = 0;
    }
}

void BiasState::doSkippedUpdatesForAllPoints(const BiasParams& params)
{
    double weightHistScaling;
    double logPmfsumScaling;

    getSkippedUpdateHistogramScaleFactors(params, &weightHistScaling, &logPmfsumScaling);

    for (auto& pointState : points_)
    {
        bool didUpdate = pointState.performPreviouslySkippedUpdates(
                params, histogramSize_.numUpdates(), weightHistScaling, logPmfsumScaling);

        /* Update the bias for this point only if there were skipped updates in the past to avoid calculating the log unnecessarily */
        if (didUpdate)
        {
            pointState.updateBias();
        }
    }
}

void BiasState::doSkippedUpdatesInNeighborhood(const BiasParams& params, const BiasGrid& grid)
{
    double weightHistScaling;
    double logPmfsumScaling;

    getSkippedUpdateHistogramScaleFactors(params, &weightHistScaling, &logPmfsumScaling);

    /* For each neighbor point of the center point, refresh its state by adding the results of all past, skipped updates. */
    const std::vector<int>& neighbors = grid.point(coordState_.gridpointIndex()).neighbor;
    for (const auto& neighbor : neighbors)
    {
        bool didUpdate = points_[neighbor].performPreviouslySkippedUpdates(
                params, histogramSize_.numUpdates(), weightHistScaling, logPmfsumScaling);

        if (didUpdate)
        {
            points_[neighbor].updateBias();
        }
    }
}

namespace
{

/*! \brief
 * Merge update lists from multiple sharing simulations.
 *
 * \param[in,out] updateList    Update list for this simulation (assumed >= npoints long).
 * \param[in]     numPoints     Total number of points.
 * \param[in]     biasSharing   Object for sharing bias data over multiple simulations
 * \param[in]     biasIndex     Index of this bias in the total list of biases in this simulation
 */
void mergeSharedUpdateLists(std::vector<int>*  updateList,
                            int                numPoints,
                            const BiasSharing& biasSharing,
                            const int          biasIndex)
{
    std::vector<int> numUpdatesOfPoint;

    /* Flag the update points of this sim.
       TODO: we can probably avoid allocating this array and just use the input array. */
    numUpdatesOfPoint.resize(numPoints, 0);
    for (auto& pointIndex : *updateList)
    {
        numUpdatesOfPoint[pointIndex] = 1;
    }

    /* Sum over the sims to get all the flagged points */
    biasSharing.sumOverSharingSimulations(arrayRefFromArray(numUpdatesOfPoint.data(), numPoints), biasIndex);

    /* Collect the indices of the flagged points in place. The resulting array will be the merged update list.*/
    updateList->clear();
    for (int m = 0; m < numPoints; m++)
    {
        if (numUpdatesOfPoint[m] > 0)
        {
            updateList->push_back(m);
        }
    }
}

/*! \brief
 * Generate an update list of points sampled since the last update.
 *
 * \param[in] grid              The AWH bias.
 * \param[in] points            The point state.
 * \param[in] originUpdatelist  The origin of the rectangular region that has been sampled since
 *                              last update.
 * \param[in] endUpdatelist     The end of the rectangular that has been sampled since
 *                              last update.
 * \param[in,out] updateList    Local update list to set (assumed >= npoints long).
 */
void makeLocalUpdateList(const BiasGrid&            grid,
                         ArrayRef<const PointState> points,
                         const awh_ivec             originUpdatelist,
                         const awh_ivec             endUpdatelist,
                         std::vector<int>*          updateList)
{
    awh_ivec origin;
    awh_ivec numPoints;

    /* Define the update search grid */
    for (int d = 0; d < grid.numDimensions(); d++)
    {
        origin[d]    = originUpdatelist[d];
        numPoints[d] = endUpdatelist[d] - originUpdatelist[d] + 1;

        /* Because the end_updatelist is unwrapped it can be > (npoints - 1) so that numPoints can be > npoints in grid.
           This helps for calculating the distance/number of points but should be removed and fixed when the way of
           updating origin/end updatelist is changed (see sampleProbabilityWeights). */
        numPoints[d] = std::min(grid.axis(d).numPoints(), numPoints[d]);
    }

    /* Make the update list */
    updateList->clear();
    int  pointIndex  = -1;
    bool pointExists = true;
    while (pointExists)
    {
        pointExists = advancePointInSubgrid(grid, origin, numPoints, &pointIndex);

        if (pointExists && points[pointIndex].inTargetRegion())
        {
            updateList->push_back(pointIndex);
        }
    }
}

} // namespace

void BiasState::resetLocalUpdateRange(const BiasGrid& grid)
{
    const int gridpointIndex = coordState_.gridpointIndex();
    for (int d = 0; d < grid.numDimensions(); d++)
    {
        /* This gives the  minimum range consisting only of the current closest point. */
        originUpdatelist_[d] = grid.point(gridpointIndex).index[d];
        endUpdatelist_[d]    = grid.point(gridpointIndex).index[d];
    }
}

namespace
{

/*! \brief
 * Add partial histograms (accumulating between updates) to accumulating histograms.
 *
 * \param[in,out] pointState         The state of the points in the bias.
 * \param[in,out] weightSumCovering  The weights for checking covering.
 * \param[in]     numSharedUpdate    The number of biases sharing the histrogram.
 * \param[in]     biasSharing        Object for sharing bias data over multiple simulations
 * \param[in]     biasIndex          Index of this bias in the total list of biases in this
 * simulation \param[in]     localUpdateList    List of points with data.
 */
void sumHistograms(gmx::ArrayRef<PointState> pointState,
                   gmx::ArrayRef<double>     weightSumCovering,
                   int                       numSharedUpdate,
                   const BiasSharing*        biasSharing,
                   const int                 biasIndex,
                   const std::vector<int>&   localUpdateList)
{
    /* The covering checking histograms are added before summing over simulations, so that the
       weights from different simulations are kept distinguishable. */
    for (int globalIndex : localUpdateList)
    {
        weightSumCovering[globalIndex] += pointState[globalIndex].weightSumIteration();
    }

    /* Sum histograms over multiple simulations if needed. */
    if (numSharedUpdate > 1)
    {
        GMX_ASSERT(numSharedUpdate == biasSharing->numSharingSimulations(biasIndex),
                   "Sharing within a simulation is not implemented (yet)");

        /* Collect the weights and counts in linear arrays to be able to use gmx_sumd_sim. */
        std::vector<double> weightSum;
        std::vector<double> coordVisits;

        weightSum.resize(localUpdateList.size());
        coordVisits.resize(localUpdateList.size());

        for (size_t localIndex = 0; localIndex < localUpdateList.size(); localIndex++)
        {
            PointState& ps = pointState[localUpdateList[localIndex]];

            weightSum[localIndex]   = ps.weightSumIteration();
            coordVisits[localIndex] = ps.numVisitsIteration();
            ps.addLocalWeightSum();
        }

        biasSharing->sumOverSharingSimulations(gmx::ArrayRef<double>(weightSum), biasIndex);
        biasSharing->sumOverSharingSimulations(gmx::ArrayRef<double>(coordVisits), biasIndex);

        /* Transfer back the result */
        for (size_t localIndex = 0; localIndex < localUpdateList.size(); localIndex++)
        {
            PointState& ps = pointState[localUpdateList[localIndex]];

            ps.setPartialWeightAndCount(weightSum[localIndex], coordVisits[localIndex]);

            /* Now add the partial counts and weights to the accumulating histograms.
            Note: we still need to use the weights for the update so we wait
            with resetting them until the end of the update. */
            ps.addPartialWeightAndCount();
        }
    }
    else
    {
        for (int globalIndex : localUpdateList)
        {
            PointState& ps = pointState[globalIndex];
            ps.addLocalWeightSum();

            /* Now add the partial counts and weights to the accumulating histograms.
            Note: we still need to use the weights for the update so we wait
            with resetting them until the end of the update. */
            ps.addPartialWeightAndCount();
        }
    }
}

/*! \brief
 * Label points along an axis as covered or not.
 *
 * A point is covered if it is surrounded by visited points up to a radius = coverRadius.
 *
 * \param[in]     visited        Visited? For each point.
 * \param[in]     checkCovering  Check for covering? For each point.
 * \param[in]     numPoints      The number of grid points along this dimension.
 * \param[in]     period         Period in number of points.
 * \param[in]     coverRadius    Cover radius, in points, needed for defining a point as covered.
 * \param[in,out] covered        In this array elements are 1 for covered points and 0 for
 * non-covered points, this routine assumes that \p covered has at least size \p numPoints.
 */
void labelCoveredPoints(const std::vector<bool>& visited,
                        const std::vector<bool>& checkCovering,
                        int                      numPoints,
                        int                      period,
                        int                      coverRadius,
                        gmx::ArrayRef<int>       covered)
{
    GMX_ASSERT(covered.ssize() >= numPoints, "covered should be at least as large as the grid");

    bool haveFirstNotVisited = false;
    int  firstNotVisited     = -1;
    int  notVisitedLow       = -1;
    int  notVisitedHigh      = -1;

    for (int n = 0; n < numPoints; n++)
    {
        if (checkCovering[n] && !visited[n])
        {
            if (!haveFirstNotVisited)
            {
                notVisitedLow       = n;
                firstNotVisited     = n;
                haveFirstNotVisited = true;
            }
            else
            {
                notVisitedHigh = n;

                /* Have now an interval I = [notVisitedLow,notVisitedHigh] of visited points bounded
                   by unvisited points. The unvisted end points affect the coveredness of the
                   visited with a reach equal to the cover radius. */
                int notCoveredLow  = notVisitedLow + coverRadius;
                int notCoveredHigh = notVisitedHigh - coverRadius;
                for (int i = notVisitedLow; i <= notVisitedHigh; i++)
                {
                    covered[i] = static_cast<int>((i > notCoveredLow) && (i < notCoveredHigh));
                }

                /* Find a new interval to set covering for. Make the notVisitedHigh of this interval
                   the notVisitedLow of the next. */
                notVisitedLow = notVisitedHigh;
            }
        }
    }

    /* Have labelled all the internal points. Now take care of the boundary regions. */
    if (!haveFirstNotVisited)
    {
        /* No non-visited points <=> all points visited => all points covered. */

        for (int n = 0; n < numPoints; n++)
        {
            covered[n] = 1;
        }
    }
    else
    {
        int lastNotVisited = notVisitedLow;

        /* For periodic boundaries, non-visited points can influence points
           on the other side of the boundary so we need to wrap around. */

        /* Lower end. For periodic boundaries the last upper end not visited point becomes the low-end not visited point.
           For non-periodic boundaries there is no lower end point so a dummy value is used. */
        int notVisitedHigh = firstNotVisited;
        int notVisitedLow  = period > 0 ? (lastNotVisited - period) : -(coverRadius + 1);

        int notCoveredLow  = notVisitedLow + coverRadius;
        int notCoveredHigh = notVisitedHigh - coverRadius;

        for (int i = 0; i <= notVisitedHigh; i++)
        {
            /* For non-periodic boundaries notCoveredLow = -1 will impose no restriction. */
            covered[i] = static_cast<int>((i > notCoveredLow) && (i < notCoveredHigh));
        }

        /* Upper end. Same as for lower end but in the other direction. */
        notVisitedHigh = period > 0 ? (firstNotVisited + period) : (numPoints + coverRadius);
        notVisitedLow  = lastNotVisited;

        notCoveredLow  = notVisitedLow + coverRadius;
        notCoveredHigh = notVisitedHigh - coverRadius;

        for (int i = notVisitedLow; i <= numPoints - 1; i++)
        {
            /* For non-periodic boundaries notCoveredHigh = numPoints will impose no restriction. */
            covered[i] = static_cast<int>((i > notCoveredLow) && (i < notCoveredHigh));
        }
    }
}

} // namespace

bool BiasState::isSamplingRegionCovered(const BiasParams&         params,
                                        ArrayRef<const DimParams> dimParams,
                                        const BiasGrid&           grid) const
{
    /* Allocate and initialize arrays: one for checking visits along each dimension,
       one for keeping track of which points to check and one for the covered points.
       Possibly these could be kept as AWH variables to avoid these allocations. */
    struct CheckDim
    {
        std::vector<bool> visited;
        std::vector<bool> checkCovering;
        // We use int for the covering array since we might use gmx_sumi_sim.
        std::vector<int> covered;
    };

    std::vector<CheckDim> checkDim;
    checkDim.resize(grid.numDimensions());

    for (int d = 0; d < grid.numDimensions(); d++)
    {
        const size_t numPoints = grid.axis(d).numPoints();
        checkDim[d].visited.resize(numPoints, false);
        checkDim[d].checkCovering.resize(numPoints, false);
        checkDim[d].covered.resize(numPoints, 0);
    }

    /* Set visited points along each dimension and which points should be checked for covering.
       Specifically, points above the free energy cutoff (if there is one) or points outside
       of the target region are ignored. */

    /* Set the free energy cutoff */
    double maxFreeEnergy = GMX_FLOAT_MAX;

    if (params.eTarget == AwhTargetType::Cutoff)
    {
        maxFreeEnergy = freeEnergyMinimumValue(points_) + params.freeEnergyCutoffInKT;
    }

    /* Set the threshold weight for a point to be considered visited. */
    double weightThreshold = 1;
    for (int d = 0; d < grid.numDimensions(); d++)
    {
        if (grid.axis(d).isFepLambdaAxis())
        {
            /* Do not modify the weight threshold based on a FEP lambda axis. The spread
             * of the sampling weights is not depending on a Gaussian distribution (like
             * below). */
            weightThreshold *= 1.0;
        }
        else
        {
            /* The spacing is proportional to 1/sqrt(betak). The weight threshold will be
             * approximately (given that the spacing can be modified if the dimension is periodic)
             * proportional to sqrt(1/(2*pi)). */
            weightThreshold *= grid.axis(d).spacing()
                               * std::sqrt(dimParams[d].pullDimParams().betak * 0.5 * M_1_PI);
        }
    }

    /* Project the sampling weights onto each dimension */
    for (size_t m = 0; m < grid.numPoints(); m++)
    {
        const PointState& pointState = points_[m];

        for (int d = 0; d < grid.numDimensions(); d++)
        {
            int n = grid.point(m).index[d];

            /* Is visited if it was already visited or if there is enough weight at the current point */
            checkDim[d].visited[n] = checkDim[d].visited[n] || (weightSumCovering_[m] > weightThreshold);

            /* Check for covering if there is at least point in this slice that is in the target region and within the cutoff */
            checkDim[d].checkCovering[n] =
                    checkDim[d].checkCovering[n]
                    || (pointState.inTargetRegion() && pointState.freeEnergy() < maxFreeEnergy);
        }
    }

    /* Label each point along each dimension as covered or not. */
    for (int d = 0; d < grid.numDimensions(); d++)
    {
        labelCoveredPoints(checkDim[d].visited,
                           checkDim[d].checkCovering,
                           grid.axis(d).numPoints(),
                           grid.axis(d).numPointsInPeriod(),
                           params.coverRadius()[d],
                           checkDim[d].covered);
    }

    /* Now check for global covering. Each dimension needs to be covered separately.
       A dimension is covered if each point is covered.  Multiple simulations collectively
       cover the points, i.e. a point is covered if any of the simulations covered it.
       However, visited points are not shared, i.e. if a point is covered or not is
       determined by the visits of a single simulation. In general the covering criterion is
       all points covered => all points are surrounded by visited points up to a radius = coverRadius.
       For 1 simulation, all points covered <=> all points visited. For multiple simulations
       however, all points visited collectively !=> all points covered, except for coverRadius = 0.
       In the limit of large coverRadius, all points covered => all points visited by at least one
       simulation (since no point will be covered until all points have been visited by a
       single simulation). Basically coverRadius sets how much "connectedness" (or mixing) a point
       needs with surrounding points before sharing covering information with other simulations. */

    /* Communicate the covered points between sharing simulations if needed. */
    if (params.numSharedUpdate > 1)
    {
        /* For multiple dimensions this may not be the best way to do it. */
        for (int d = 0; d < grid.numDimensions(); d++)
        {
            biasSharing_->sumOverSharingSimulations(
                    gmx::arrayRefFromArray(checkDim[d].covered.data(), grid.axis(d).numPoints()),
                    params.biasIndex_);
        }
    }

    /* Now check if for each dimension all points are covered. Break if not true. */
    bool allPointsCovered = true;
    for (int d = 0; d < grid.numDimensions() && allPointsCovered; d++)
    {
        for (int n = 0; n < grid.axis(d).numPoints() && allPointsCovered; n++)
        {
            allPointsCovered = (checkDim[d].covered[n] != 0);
        }
    }

    return allPointsCovered;
}

/*! \brief
 * Normalizes the free energy and PMF sum.
 *
 * \param[in] pointState  The state of the points.
 */
static void normalizeFreeEnergyAndPmfSum(std::vector<PointState>* pointState)
{
    double minF = freeEnergyMinimumValue(*pointState);

    for (PointState& ps : *pointState)
    {
        ps.normalizeFreeEnergyAndPmfSum(minF);
    }
}

void BiasState::updateFreeEnergyAndAddSamplesToHistogram(ArrayRef<const DimParams> dimParams,
                                                         const BiasGrid&           grid,
                                                         const BiasParams&         params,
                                                         const CorrelationGrid&    forceCorrelation,
                                                         double                    t,
                                                         int64_t                   step,
                                                         FILE*                     fplog,
                                                         std::vector<int>*         updateList)
{
    /* Note hat updateList is only used in this scope and is always
     * re-initialized. We do not use a local vector, because that would
     * cause reallocation every time this funtion is called and the vector
     * can be the size of the whole grid.
     */

    /* Make a list of all local points, i.e. those that could have been touched since
       the last update. These are the points needed for summing histograms below
       (non-local points only add zeros). For local updates, this will also be the
       final update list. */
    makeLocalUpdateList(grid, points_, originUpdatelist_, endUpdatelist_, updateList);
    if (params.numSharedUpdate > 1)
    {
        mergeSharedUpdateLists(updateList, points_.size(), *biasSharing_, params.biasIndex_);
    }

    /* Reset the range for the next update */
    resetLocalUpdateRange(grid);

    /* Add samples to histograms for all local points and sync simulations if needed */
    sumHistograms(points_, weightSumCovering_, params.numSharedUpdate, biasSharing_, params.biasIndex_, *updateList);

    sumPmf(points_, params.numSharedUpdate, biasSharing_, params.biasIndex_);

    /* Renormalize the free energy if values are too large. */
    bool needToNormalizeFreeEnergy = false;
    for (int& globalIndex : *updateList)
    {
        /* We want to keep the absolute value of the free energies to be less
           c_largePositiveExponent to be able to safely pass these values to exp(). The check below
           ensures this as long as the free energy values grow less than 0.5*c_largePositiveExponent
           in a return time to this neighborhood. For reasonable update sizes it's unlikely that
           this requirement would be broken. */
        if (std::abs(points_[globalIndex].freeEnergy()) > 0.5 * detail::c_largePositiveExponent)
        {
            needToNormalizeFreeEnergy = true;
            break;
        }
    }

    /* Update target distribution? */
    bool doScaleTargetByMetric = params.scaleTargetByMetric && !inInitialStage();
    bool needToUpdateTargetDistribution =
            ((params.eTarget != AwhTargetType::Constant || doScaleTargetByMetric)
             && params.isUpdateTargetStep(step));

    /* In the initial stage, the histogram grows dynamically as a function of the number of coverings. */
    bool detectedCovering = false;
    if (inInitialStage())
    {
        detectedCovering =
                (params.isCheckCoveringStep(step) && isSamplingRegionCovered(params, dimParams, grid));
    }

    /* The weighthistogram size after this update. */
    double newHistogramSize = histogramSize_.newHistogramSize(
            params, t, detectedCovering, points_, weightSumCovering_, fplog);

    /* Make the update list. Usually we try to only update local points,
     * but if the update has non-trivial or non-deterministic effects
     * on non-local points a global update is needed. This is the case when:
     * 1) a covering occurred in the initial stage, leading to non-trivial
     *    histogram rescaling factors; or
     * 2) the target distribution will be updated, since we don't make any
     *    assumption on its form; or
     * 3) the AWH parameters are such that we never attempt to skip non-local
     *    updates; or
     * 4) the free energy values have grown so large that a renormalization
     *    is needed.
     */
    if (needToUpdateTargetDistribution || detectedCovering || !params.skipUpdates() || needToNormalizeFreeEnergy)
    {
        /* Global update, just add all points. */
        updateList->clear();
        for (size_t m = 0; m < points_.size(); m++)
        {
            if (points_[m].inTargetRegion())
            {
                updateList->push_back(m);
            }
        }
    }

    /* Set histogram scale factors. */
    double weightHistScalingSkipped = 0;
    double logPmfsumScalingSkipped  = 0;
    if (params.skipUpdates())
    {
        getSkippedUpdateHistogramScaleFactors(params, &weightHistScalingSkipped, &logPmfsumScalingSkipped);
    }
    double weightHistScalingNew;
    double logPmfsumScalingNew;
    setHistogramUpdateScaleFactors(
            params, newHistogramSize, histogramSize_.histogramSize(), &weightHistScalingNew, &logPmfsumScalingNew);

    /* Update free energy and reference weight histogram for points in the update list. */
    for (int pointIndex : *updateList)
    {
        PointState* pointStateToUpdate = &points_[pointIndex];

        /* Do updates from previous update steps that were skipped because this point was at that time non-local. */
        if (params.skipUpdates())
        {
            pointStateToUpdate->performPreviouslySkippedUpdates(
                    params, histogramSize_.numUpdates(), weightHistScalingSkipped, logPmfsumScalingSkipped);
        }

        /* Now do an update with new sampling data. */
        pointStateToUpdate->updateWithNewSampling(
                params, histogramSize_.numUpdates(), weightHistScalingNew, logPmfsumScalingNew);
    }

    /* Only update the histogram size after we are done with the local point updates */
    histogramSize_.setHistogramSize(newHistogramSize, weightHistScalingNew);

    if (needToNormalizeFreeEnergy)
    {
        normalizeFreeEnergyAndPmfSum(&points_);
    }

    if (needToUpdateTargetDistribution)
    {
        /* The target distribution is always updated for all points at once. */
        updateTargetDistribution(params, forceCorrelation);
    }

    /* Update the bias. The bias is updated separately and last since it simply a function of
       the free energy and the target distribution and we want to avoid doing extra work. */
    for (int pointIndex : *updateList)
    {
        points_[pointIndex].updateBias();
    }

    /* Increase the update counter. */
    histogramSize_.incrementNumUpdates();
}

double BiasState::updateProbabilityWeightsAndConvolvedBias(ArrayRef<const DimParams> dimParams,
                                                           const BiasGrid&           grid,
                                                           ArrayRef<const double> neighborLambdaEnergies,
                                                           std::vector<double, AlignedAllocator<double>>* weight) const
{
    /* Only neighbors of the current coordinate value will have a non-negligible chance of getting sampled */
    const std::vector<int>& neighbors = grid.point(coordState_.gridpointIndex()).neighbor;

#if GMX_SIMD_HAVE_DOUBLE
    typedef SimdDouble PackType;
    constexpr int      packSize = GMX_SIMD_DOUBLE_WIDTH;
#else
    typedef double PackType;
    constexpr int  packSize = 1;
#endif
    /* Round the size of the weight array up to packSize */
    const int weightSize = ((neighbors.size() + packSize - 1) / packSize) * packSize;
    weight->resize(weightSize);

    double* gmx_restrict weightData = weight->data();
    PackType             weightSumPack(0.0);
    for (size_t i = 0; i < neighbors.size(); i += packSize)
    {
        for (size_t n = i; n < i + packSize; n++)
        {
            if (n < neighbors.size())
            {
                const int neighbor = neighbors[n];
                (*weight)[n]       = biasedLogWeightFromPoint(dimParams,
                                                        points_,
                                                        grid,
                                                        neighbor,
                                                        points_[neighbor].bias(),
                                                        coordState_.coordValue(),
                                                        neighborLambdaEnergies,
                                                        coordState_.gridpointIndex());
            }
            else
            {
                /* Pad with values that don't affect the result */
                (*weight)[n] = detail::c_largeNegativeExponent;
            }
        }
        PackType weightPack = load<PackType>(weightData + i);
        weightPack          = gmx::exp(weightPack);
        weightSumPack       = weightSumPack + weightPack;
        store(weightData + i, weightPack);
    }
    /* Sum of probability weights */
    double weightSum = reduce(weightSumPack);
    GMX_RELEASE_ASSERT(weightSum > 0,
                       "zero probability weight when updating AWH probability weights.");

    /* Normalize probabilities to sum to 1 */
    double invWeightSum = 1 / weightSum;

    /* When there is a free energy lambda state axis remove the convolved contributions along that
     * axis from the total bias. This must be done after calculating invWeightSum (since weightSum
     * will be modified), but before normalizing the weights (below). */
    if (grid.hasLambdaAxis())
    {
        /* If there is only one axis the bias will not be convolved in any dimension. */
        if (grid.axis().size() == 1)
        {
            weightSum = gmx::exp(points_[coordState_.gridpointIndex()].bias());
        }
        else
        {
            for (size_t i = 0; i < neighbors.size(); i++)
            {
                const int neighbor = neighbors[i];
                if (pointsHaveDifferentLambda(grid, coordState_.gridpointIndex(), neighbor))
                {
                    weightSum -= weightData[i];
                }
            }
            /* Subtracting potentially very large values above may lead to rounding errors, causing
             * negative weightSum, even in double precision. */
            weightSum = std::max(weightSum, GMX_DOUBLE_MIN);
        }
    }

    for (double& w : *weight)
    {
        w *= invWeightSum;
    }

    /* Return the convolved bias */
    return std::log(weightSum);
}

double BiasState::calcConvolvedBias(ArrayRef<const DimParams> dimParams,
                                    const BiasGrid&           grid,
                                    const awh_dvec&           coordValue) const
{
    int              point     = grid.nearestIndex(coordValue);
    const GridPoint& gridPoint = grid.point(point);

    /* Sum the probability weights from the neighborhood of the given point */
    double weightSum = 0;
    for (int neighbor : gridPoint.neighbor)
    {
        /* No convolution is required along the lambda dimension. */
        if (pointsHaveDifferentLambda(grid, point, neighbor))
        {
            continue;
        }
        double logWeight = biasedLogWeightFromPoint(
                dimParams, points_, grid, neighbor, points_[neighbor].bias(), coordValue, {}, point);
        weightSum += std::exp(logWeight);
    }

    /* Returns -GMX_FLOAT_MAX if no neighboring points were in the target region. */
    return (weightSum > 0) ? std::log(weightSum) : -GMX_FLOAT_MAX;
}

void BiasState::sampleProbabilityWeights(const BiasGrid& grid, gmx::ArrayRef<const double> probWeightNeighbor)
{
    const std::vector<int>& neighbor = grid.point(coordState_.gridpointIndex()).neighbor;

    /* Save weights for next update */
    for (size_t n = 0; n < neighbor.size(); n++)
    {
        points_[neighbor[n]].increaseWeightSumIteration(probWeightNeighbor[n]);
    }

    /* Update the local update range. Two corner points define this rectangular
     * domain. We need to choose two new corner points such that the new domain
     * contains both the old update range and the current neighborhood.
     * In the simplest case when an update is performed every sample,
     * the update range would simply equal the current neighborhood.
     */
    int neighborStart = neighbor[0];
    int neighborLast  = neighbor[neighbor.size() - 1];
    for (int d = 0; d < grid.numDimensions(); d++)
    {
        int origin = grid.point(neighborStart).index[d];
        int last   = grid.point(neighborLast).index[d];

        if (origin > last)
        {
            /* Unwrap if wrapped around the boundary (only happens for periodic
             * boundaries). This has been already for the stored index interval.
             */
            /* TODO: what we want to do is to find the smallest the update
             * interval that contains all points that need to be updated.
             * This amounts to combining two intervals, the current
             * [origin, end] update interval and the new touched neighborhood
             * into a new interval that contains all points from both the old
             * intervals.
             *
             * For periodic boundaries it becomes slightly more complicated
             * than for closed boundaries because then it needs not be
             * true that origin < end (so one can't simply relate the origin/end
             * in the min()/max() below). The strategy here is to choose the
             * origin closest to a reference point (index 0) and then unwrap
             * the end index if needed and choose the largest end index.
             * This ensures that both intervals are in the new interval
             * but it's not necessarily the smallest.
             * Currently we solve this by going through each possibility
             * and checking them.
             */
            last += grid.axis(d).numPointsInPeriod();
        }

        originUpdatelist_[d] = std::min(originUpdatelist_[d], origin);
        endUpdatelist_[d]    = std::max(endUpdatelist_[d], last);
    }
}

void BiasState::sampleCoordAndPmf(const std::vector<DimParams>& dimParams,
                                  const BiasGrid&               grid,
                                  gmx::ArrayRef<const double>   probWeightNeighbor,
                                  double                        convolvedBias)
{
    /* Sampling-based deconvolution extracting the PMF.
     * Update the PMF histogram with the current coordinate value.
     *
     * Because of the finite width of the harmonic potential, the free energy
     * defined for each coordinate point does not exactly equal that of the
     * actual coordinate, the PMF. However, the PMF can be estimated by applying
     * the relation exp(-PMF) = exp(-bias_convolved)*P_biased/Z, i.e. by keeping a
     * reweighted histogram of the coordinate value. Strictly, this relies on
     * the unknown normalization constant Z being either constant or known. Here,
     * neither is true except in the long simulation time limit. Empirically however,
     * it works (mainly because how the PMF histogram is rescaled).
     */

    const int                gridPointIndex  = coordState_.gridpointIndex();
    const std::optional<int> lambdaAxisIndex = grid.lambdaAxisIndex();

    /* Update the PMF of points along a lambda axis with their bias. */
    if (lambdaAxisIndex)
    {
        const std::vector<int>& neighbors = grid.point(gridPointIndex).neighbor;

        std::vector<double> lambdaMarginalDistribution =
                calculateFELambdaMarginalDistribution(grid, neighbors, probWeightNeighbor);

        awh_dvec coordValueAlongLambda = { coordState_.coordValue()[0],
                                           coordState_.coordValue()[1],
                                           coordState_.coordValue()[2],
                                           coordState_.coordValue()[3] };
        for (size_t i = 0; i < neighbors.size(); i++)
        {
            const int neighbor = neighbors[i];
            double    bias;
            if (pointsAlongLambdaAxis(grid, gridPointIndex, neighbor))
            {
                const double neighborLambda = grid.point(neighbor).coordValue[lambdaAxisIndex.value()];
                if (neighbor == gridPointIndex)
                {
                    bias = convolvedBias;
                }
                else
                {
                    coordValueAlongLambda[lambdaAxisIndex.value()] = neighborLambda;
                    bias = calcConvolvedBias(dimParams, grid, coordValueAlongLambda);
                }

                const double probWeight = lambdaMarginalDistribution[neighborLambda];
                const double weightedBias = bias - std::log(std::max(probWeight, GMX_DOUBLE_MIN)); // avoid log(0)

                if (neighbor == gridPointIndex && grid.covers(coordState_.coordValue()))
                {
                    points_[neighbor].samplePmf(weightedBias);
                }
                else
                {
                    points_[neighbor].updatePmfUnvisited(weightedBias);
                }
            }
        }
    }
    else
    {
        /* Only save coordinate data that is in range (the given index is always
         * in range even if the coordinate value is not).
         */
        if (grid.covers(coordState_.coordValue()))
        {
            /* Save PMF sum and keep a histogram of the sampled coordinate values */
            points_[gridPointIndex].samplePmf(convolvedBias);
        }
    }

    /* Save probability weights for the update */
    sampleProbabilityWeights(grid, probWeightNeighbor);
}

void BiasState::initHistoryFromState(AwhBiasHistory* biasHistory) const
{
    biasHistory->pointState.resize(points_.size());
}

void BiasState::updateHistory(AwhBiasHistory* biasHistory, const BiasGrid& grid) const
{
    GMX_RELEASE_ASSERT(biasHistory->pointState.size() == points_.size(),
                       "The AWH history setup does not match the AWH state.");

    AwhBiasStateHistory* stateHistory = &biasHistory->state;
    stateHistory->umbrellaGridpoint   = coordState_.umbrellaGridpoint();

    for (size_t m = 0; m < biasHistory->pointState.size(); m++)
    {
        AwhPointStateHistory* psh = &biasHistory->pointState[m];

        points_[m].storeState(psh);

        psh->weightsum_covering = weightSumCovering_[m];
    }

    histogramSize_.storeState(stateHistory);

    stateHistory->origin_index_updatelist = multiDimGridIndexToLinear(grid, originUpdatelist_);
    stateHistory->end_index_updatelist    = multiDimGridIndexToLinear(grid, endUpdatelist_);
}

void BiasState::restoreFromHistory(const AwhBiasHistory& biasHistory, const BiasGrid& grid)
{
    const AwhBiasStateHistory& stateHistory = biasHistory.state;

    coordState_.restoreFromHistory(stateHistory);

    if (biasHistory.pointState.size() != points_.size())
    {
        GMX_THROW(
                InvalidInputError("Bias grid size in checkpoint and simulation do not match. "
                                  "Likely you provided a checkpoint from a different simulation."));
    }
    for (size_t m = 0; m < points_.size(); m++)
    {
        points_[m].setFromHistory(biasHistory.pointState[m]);
    }

    for (size_t m = 0; m < weightSumCovering_.size(); m++)
    {
        weightSumCovering_[m] = biasHistory.pointState[m].weightsum_covering;
    }

    histogramSize_.restoreFromHistory(stateHistory);

    linearGridindexToMultiDim(grid, stateHistory.origin_index_updatelist, originUpdatelist_);
    linearGridindexToMultiDim(grid, stateHistory.end_index_updatelist, endUpdatelist_);
}

void BiasState::broadcast(const t_commrec* commRecord)
{
    gmx_bcast(sizeof(coordState_), &coordState_, commRecord->mpi_comm_mygroup);

    gmx_bcast(points_.size() * sizeof(PointState), points_.data(), commRecord->mpi_comm_mygroup);

    gmx_bcast(weightSumCovering_.size() * sizeof(double), weightSumCovering_.data(), commRecord->mpi_comm_mygroup);

    gmx_bcast(sizeof(histogramSize_), &histogramSize_, commRecord->mpi_comm_mygroup);
}

void BiasState::setFreeEnergyToConvolvedPmf(ArrayRef<const DimParams> dimParams, const BiasGrid& grid)
{
    std::vector<float> convolvedPmf;

    calcConvolvedPmf(dimParams, grid, &convolvedPmf);

    for (size_t m = 0; m < points_.size(); m++)
    {
        points_[m].setFreeEnergy(convolvedPmf[m]);
    }
}

void BiasState::updateSharedCorrelationTensorTimeIntegral(const BiasParams&      biasParams,
                                                          const CorrelationGrid& forceCorrelation,
                                                          bool shareAcrossAllRanks)
{
    const int numCorrelation = forceCorrelation.tensorSize();
    const int numPoints      = points_.size();

    GMX_ASSERT(numPoints > 0, "There must be points in the AWH coordinate grid.");
    GMX_ASSERT(static_cast<int>(sharedCorrelationTensorTimeIntegral_.size()) == numPoints,
               "The number of points in sharedCorrelationTensorTimeIntegral_ does not match the "
               "number of AWH points.");
    GMX_ASSERT(
            static_cast<int>(sharedCorrelationTensorTimeIntegral_[0].size()) == numCorrelation,
            "The number of correlations in sharedCorrelationTensorTimeIntegral_ does not match the "
            "AWH correlation grid tensor size.");

    std::vector<double> buffer(numCorrelation * numPoints, 0);

    if (biasParams.numSharedUpdate > 1)
    {
        GMX_ASSERT(biasSharing_ != nullptr
                           && biasParams.numSharedUpdate
                                              % biasSharing_->numSharingSimulations(biasParams.biasIndex_)
                                      == 0,
                   "numSharedUpdate should be a multiple of multiSimComm->numSimulations_");
        GMX_ASSERT(biasParams.numSharedUpdate == biasSharing_->numSharingSimulations(biasParams.biasIndex_),
                   "Sharing within a simulation is not implemented (yet)");

        for (int gridPointIndex = 0; gridPointIndex < numPoints; gridPointIndex++)
        {
            if (points_[gridPointIndex].inTargetRegion())
            {
                for (int correlationTensorIndex = 0; correlationTensorIndex < numCorrelation;
                     correlationTensorIndex++)
                {
                    int index     = gridPointIndex * numCorrelation + correlationTensorIndex;
                    buffer[index] = forceCorrelation.tensors()[gridPointIndex].getTimeIntegral(
                                            correlationTensorIndex, forceCorrelation.dtSample_)
                                    * points_[gridPointIndex].localWeightSum();
                }
            }
        }

        if (shareAcrossAllRanks)
        {
            biasSharing_->sumOverSharingSimulations(buffer, biasParams.biasIndex_);
        }
        else
        {
            biasSharing_->sumOverSharingMainRanks(buffer, biasParams.biasIndex_);
        }

        for (int gridPointIndex = 0; gridPointIndex < numPoints; gridPointIndex++)
        {
            for (int correlationTensorIndex = 0; correlationTensorIndex < numCorrelation;
                 correlationTensorIndex++)
            {
                int index = gridPointIndex * numCorrelation + correlationTensorIndex;
                if (points_[gridPointIndex].weightSumTot() > 0 && points_[gridPointIndex].inTargetRegion())
                {
                    sharedCorrelationTensorTimeIntegral_[gridPointIndex][correlationTensorIndex] =
                            buffer[index] / points_[gridPointIndex].weightSumTot();
                }
                else
                {
                    sharedCorrelationTensorTimeIntegral_[gridPointIndex][correlationTensorIndex] = 0;
                }
            }
        }
    }
    else
    {
        for (int gridPointIndex = 0; gridPointIndex < numPoints; gridPointIndex++)
        {
            for (int correlationTensorIndex = 0; correlationTensorIndex < numCorrelation;
                 correlationTensorIndex++)
            {
                sharedCorrelationTensorTimeIntegral_[gridPointIndex][correlationTensorIndex] =
                        forceCorrelation.tensors()[gridPointIndex].getTimeIntegral(
                                correlationTensorIndex, forceCorrelation.dtSample_);
            }
        }
    }
}

double BiasState::getSharedCorrelationTensorTimeIntegral(const int gridPointIndex,
                                                         const int correlationTensorIndex) const
{
    if (!points_[gridPointIndex].inTargetRegion() || points_[gridPointIndex].weightSumTot() <= 0)
    {
        return 0;
    }

    return sharedCorrelationTensorTimeIntegral_[gridPointIndex][correlationTensorIndex];
}

const std::vector<double>& BiasState::getSharedPointCorrelationIntegral(const int gridPointIndex) const
{
    if (!points_[gridPointIndex].inTargetRegion() || points_[gridPointIndex].weightSumTot() <= 0)
    {
        for (size_t i = 0; i < sharedCorrelationTensorTimeIntegral_[gridPointIndex].size(); i++)
        {
            GMX_RELEASE_ASSERT(sharedCorrelationTensorTimeIntegral_[gridPointIndex][i] == 0,
                               "Correlation tensor time integral of unvisited points should be 0.");
        }
    }

    return sharedCorrelationTensorTimeIntegral_[gridPointIndex];
}

/*! \brief
 * Count trailing data rows containing only zeros.
 *
 * \param[in] data        2D data array.
 * \param[in] numRows     Number of rows in array.
 * \param[in] numColumns  Number of cols in array.
 * \returns the number of trailing zero rows.
 */
static int countTrailingZeroRows(const MultiDimArray<std::vector<double>, dynamicExtents2D>& data,
                                 int numRows,
                                 int numColumns)
{
    const auto& dataView = data.asConstView();

    int numZeroRows = 0;
    for (int m = numRows - 1; m >= 0; m--)
    {
        bool rowIsZero = true;
        for (int d = 0; d < numColumns; d++)
        {
            if (dataView[d][m] != 0)
            {
                rowIsZero = false;
                break;
            }
        }

        if (!rowIsZero)
        {
            /* At a row with non-zero data */
            break;
        }
        else
        {
            /* Still at a zero data row, keep checking rows higher up. */
            numZeroRows++;
        }
    }

    return numZeroRows;
}

/*! \brief
 * Initializes the PMF and target with data read from an input table.
 *
 * \param[in]     dimParams   The dimension parameters.
 * \param[in]     grid        The grid.
 * \param[in]     filename    The filename to read PMF and target from.
 * \param[in]     numBias     Number of biases.
 * \param[in]     biasIndex   The index of the bias.
 * \param[in,out] pointState  The state of the points in this bias.
 */
static void readUserPmfAndTargetDistribution(ArrayRef<const DimParams> dimParams,
                                             const BiasGrid&           grid,
                                             const std::string&        filename,
                                             int                       numBias,
                                             int                       biasIndex,
                                             std::vector<PointState>*  pointState)
{
    /* Read the PMF and target distribution.
       From the PMF, the convolved PMF, or the reference value free energy, can be calculated
       base on the force constant. The free energy and target together determine the bias.
     */
    std::string filenameModified(filename);
    if (numBias > 1)
    {
        size_t n = filenameModified.rfind('.');
        GMX_RELEASE_ASSERT(n != std::string::npos,
                           "The filename should contain an extension starting with .");
        filenameModified.insert(n, formatString("%d", biasIndex));
    }

    std::string correctFormatMessage = formatString(
            "%s is expected in the following format. "
            "The first ndim column(s) should contain the coordinate values for each point, "
            "each column containing values of one dimension (in ascending order). "
            "For a multidimensional coordinate, points should be listed "
            "in the order obtained by traversing lower dimensions first. "
            "E.g. for two-dimensional grid of size nxn: "
            "(1, 1), (1, 2),..., (1, n), (2, 1), (2, 2), ..., , (n, n - 1), (n, n). "
            "Column ndim +  1 should contain the PMF value for each coordinate value. "
            "The target distribution values should be in column ndim + 2  or column ndim + 5. "
            "Make sure the input file ends with a new line but has no trailing new lines.",
            filename.c_str());
    gmx::TextLineWrapper wrapper;
    wrapper.settings().setLineLength(c_linewidth);
    correctFormatMessage = wrapper.wrapToString(correctFormatMessage);

    gmx::MultiDimArray<std::vector<double>, gmx::dynamicExtents2D> data = readXvgData(filenameModified);
    const int                                                      numColumns = data.extent(0);
    const int                                                      numRows    = data.extent(1);

    /* Check basic data properties here. BiasGrid takes care of more complicated things. */

    if (numRows <= 0)
    {
        std::string mesg = gmx::formatString(
                "%s is empty!.\n\n%s", filename.c_str(), correctFormatMessage.c_str());
        GMX_THROW(InvalidInputError(mesg));
    }

    /* Less than 2 points is not useful for PMF or target. */
    if (numRows < 2)
    {
        std::string mesg = gmx::formatString(
                "%s contains too few data points (%d)."
                "The minimum number of points is 2.",
                filename.c_str(),
                numRows);
        GMX_THROW(InvalidInputError(mesg));
    }

    /* Make sure there are enough columns of data.

       Two formats are allowed. Either with columns  {coords, PMF, target} or
       {coords, PMF, x, y, z, target, ...}. The latter format is allowed since that
       is how AWH output is written (x, y, z being other AWH variables). For this format,
       trailing columns are ignored.
     */
    int columnIndexTarget;
    int numColumnsMin  = dimParams.size() + 2;
    int columnIndexPmf = dimParams.size();
    if (numColumns == numColumnsMin)
    {
        columnIndexTarget = columnIndexPmf + 1;
    }
    else
    {
        columnIndexTarget = columnIndexPmf + 4;
    }

    if (numColumns < numColumnsMin)
    {
        std::string mesg = gmx::formatString(
                "The number of columns in %s should be at least %d."
                "\n\n%s",
                filename.c_str(),
                numColumnsMin,
                correctFormatMessage.c_str());
        GMX_THROW(InvalidInputError(mesg));
    }

    const auto& dataView = data.asView();

    /* read_xvg can give trailing zero data rows for trailing new lines in the input. We allow 1 zero row,
       since this could be real data. But multiple trailing zero rows cannot correspond to valid data. */
    int numZeroRows = countTrailingZeroRows(data, numRows, numColumns);
    if (numZeroRows > 1)
    {
        std::string mesg = gmx::formatString(
                "Found %d trailing zero data rows in %s. Please remove trailing empty lines and "
                "try again.",
                numZeroRows,
                filename.c_str());
        GMX_THROW(InvalidInputError(mesg));
    }

    /* Convert from user units to internal units before sending the data of to grid. */
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        double scalingFactor = dimParams[d].scaleUserInputToInternal(1);
        if (scalingFactor == 1)
        {
            continue;
        }
        for (size_t m = 0; m < pointState->size(); m++)
        {
            dataView[d][m] *= scalingFactor;
        }
    }

    /* Get a data point for each AWH grid point so that they all get data. */
    std::vector<int> gridIndexToDataIndex(grid.numPoints());
    mapGridToDataGrid(&gridIndexToDataIndex, data, numRows, filename, grid, correctFormatMessage);

    // The bounds for the PMF such that exp(pmf) doesn't over/underflow in double precision
    const double c_pmfMax = 700.0;

    /* Extract the data for each grid point.
     * We check if the target distribution is zero for all points.
     */
    bool targetDistributionIsZero = true;
    for (size_t m = 0; m < pointState->size(); m++)
    {
        const double pmf = dataView[columnIndexPmf][gridIndexToDataIndex[m]];
        if (pmf < -c_pmfMax || pmf > c_pmfMax)
        {
            GMX_THROW(InvalidInputError(
                    "A value in the user input PMF is beyond the bounds of +-700 kT"));
        }
        (*pointState)[m].setLogPmfSum(-pmf);
        double target = dataView[columnIndexTarget][gridIndexToDataIndex[m]];

        /* Check if the values are allowed. */
        if (target < 0)
        {
            std::string mesg = gmx::formatString(
                    "Target distribution weight at point %zu (%g) in %s is negative.",
                    m,
                    target,
                    filename.c_str());
            GMX_THROW(InvalidInputError(mesg));
        }
        if (target > 0)
        {
            targetDistributionIsZero = false;
        }
        (*pointState)[m].setTargetConstantWeight(target);
    }

    if (targetDistributionIsZero)
    {
        std::string mesg =
                gmx::formatString("The target weights given in column %d in %s are all 0",
                                  columnIndexTarget,
                                  filename.c_str());
        GMX_THROW(InvalidInputError(mesg));
    }
}

void BiasState::normalizePmf(int numSharingSims)
{
    /* The normalization of the PMF estimate matters because it determines how big effect the next sample has.
       Approximately (for large enough force constant) we should have:
       sum_x(exp(-pmf(x)) = nsamples*sum_xref(exp(-f(xref)).
     */

    /* Calculate the normalization factor, i.e. divide by the pmf sum, multiply by the number of samples and the f sum */
    double expSumPmf = 0;
    double expSumF   = 0;
    for (const PointState& pointState : points_)
    {
        if (pointState.inTargetRegion())
        {
            expSumPmf += std::exp(pointState.logPmfSum());
            expSumF += std::exp(-pointState.freeEnergy());
        }
    }
    double numSamples = histogramSize_.histogramSize() / numSharingSims;

    /* Renormalize */
    double logRenorm = std::log(numSamples * expSumF / expSumPmf);
    for (PointState& pointState : points_)
    {
        if (pointState.inTargetRegion())
        {
            pointState.setLogPmfSum(pointState.logPmfSum() + logRenorm);
        }
    }
}

void BiasState::initGridPointState(const AwhBiasParams&      awhBiasParams,
                                   ArrayRef<const DimParams> dimParams,
                                   const BiasGrid&           grid,
                                   const BiasParams&         params,
                                   const CorrelationGrid&    forceCorrelation,
                                   const std::string&        filename,
                                   int                       numBias)
{
    /* Modify PMF, free energy and the constant target distribution factor
     * to user input values if there is data given.
     */
    if (awhBiasParams.userPMFEstimate())
    {
        readUserPmfAndTargetDistribution(dimParams, grid, filename, numBias, params.biasIndex_, &points_);
        setFreeEnergyToConvolvedPmf(dimParams, grid);
    }

    /* The local Boltzmann distribution is special because the target distribution is updated as a function of the reference weighthistogram. */
    GMX_RELEASE_ASSERT(params.eTarget != AwhTargetType::LocalBoltzmann || points_[0].weightSumRef() != 0,
                       "AWH reference weight histogram not initialized properly with local "
                       "Boltzmann target distribution.");

    updateTargetDistribution(params, forceCorrelation);

    for (PointState& pointState : points_)
    {
        if (pointState.inTargetRegion())
        {
            pointState.updateBias();
        }
        else
        {
            /* Note that for zero target this is a value that represents -infinity but should not be used for biasing. */
            pointState.setTargetToZero();
        }
    }

    /* Set the initial reference weighthistogram. */
    const double histogramSize = histogramSize_.histogramSize();
    for (auto& pointState : points_)
    {
        pointState.setInitialReferenceWeightHistogram(histogramSize);
    }

    /* Make sure the pmf is normalized consistently with the histogram size.
       Note: the target distribution and free energy need to be set here. */
    normalizePmf(params.numSharedUpdate);
}

BiasState::BiasState(const AwhBiasParams&      awhBiasParams,
                     double                    histogramSizeInitial,
                     ArrayRef<const DimParams> dimParams,
                     const BiasGrid&           grid,
                     const BiasSharing*        biasSharing) :
    coordState_(awhBiasParams, dimParams, grid),
    points_(grid.numPoints()),
    weightSumCovering_(grid.numPoints()),
    histogramSize_(awhBiasParams, histogramSizeInitial),
    biasSharing_(biasSharing),
    /* The number of tensor elements are dim*(dim+1)/2 */
    sharedCorrelationTensorTimeIntegral_(
            grid.numPoints(),
            std::vector<double>(dimParams.size() * (dimParams.size() + 1) / 2, 0))
{
    /* The minimum and maximum multidimensional point indices that are affected by the next update */
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        int index            = grid.point(coordState_.gridpointIndex()).index[d];
        originUpdatelist_[d] = index;
        endUpdatelist_[d]    = index;
    }
}

} // namespace gmx
