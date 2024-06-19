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
 * Implements the HistogramSize class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "histogramsize.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>

#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "biasparams.h"
#include "pointstate.h"

namespace gmx
{

HistogramSize::HistogramSize(const AwhBiasParams& awhBiasParams, double histogramSizeInitial) :
    numUpdates_(0),
    histogramSize_(histogramSizeInitial),
    inInitialStage_(awhBiasParams.growthType() == AwhHistogramGrowthType::ExponentialLinear),
    growthFactor_(awhBiasParams.growthFactor()),
    equilibrateHistogram_(awhBiasParams.equilibrateHistogram()),
    logScaledSampleWeight_(0),
    maxLogScaledSampleWeight_(0),
    havePrintedAboutCovering_(false)
{
}

double HistogramSize::newHistogramSizeInitialStage(const BiasParams& params,
                                                   double            t,
                                                   bool              detectedCovering,
                                                   ArrayRef<double>  weightsumCovering,
                                                   FILE*             fplog)
{
    /* The histogram size is kept constant until the sampling region has been covered
       and the current sample weight is large enough and the histogram is ready. */
    if (!detectedCovering || (logScaledSampleWeight_ < maxLogScaledSampleWeight_) || equilibrateHistogram_)
    {
        return histogramSize_;
    }

    /* Reset the covering weight histogram. If we got this far we are either entering a
       new covering stage with a new covering histogram or exiting the initial stage
       altogether. */
    std::fill(weightsumCovering.begin(), weightsumCovering.end(), 0);

    /*  The current sample weigth is now the maximum. */
    double prevMaxLogScaledSampleWeight = maxLogScaledSampleWeight_;
    maxLogScaledSampleWeight_           = logScaledSampleWeight_;

    /* Increase the histogram size by a constant scale factor if we can, i.e. if the sample weight
       resulting from such a scaling is still larger than the previous maximum sample weight
       (ensuring that the sample weights at the end of each covering stage are monotonically
       increasing). If we cannot, exit the initial stage without changing the histogram size. */

    /* The scale factor is in most cases very close to the histogram growth factor. */
    double scaleFactor =
            growthFactor_ / (1. + params.updateWeight * params.localWeightScaling / histogramSize_);

    bool exitInitialStage =
            (logScaledSampleWeight_ - std::log(scaleFactor) <= prevMaxLogScaledSampleWeight);
    double newHistogramSize = exitInitialStage ? histogramSize_ : histogramSize_ * growthFactor_;

    /* Update the AWH bias about the exit. */
    inInitialStage_ = !exitInitialStage;

    /* Print information about coverings and if there was an exit. */
    if (fplog != nullptr)
    {
        std::string prefix = gmx::formatString("\nawh%d:", params.biasIndex_ + 1);
        fprintf(fplog, "%s covering at t = %g ps. Decreased the update size.\n", prefix.c_str(), t);

        if (exitInitialStage)
        {
            fprintf(fplog, "%s out of the initial stage at t = %g.\n", prefix.c_str(), t);
            /* It would be nice to have a way of estimating a minimum time until exit but it
               is difficult because the exit time is determined by how long it takes to cover
               relative to the time it takes to "regaining" enough sample weight. The latter
               is easy to calculate, but how the former depends on the histogram size
               is not known. */
        }
        fflush(fplog);
    }
    return newHistogramSize;
}

namespace
{

/*! \brief
 * Checks if the histogram has equilibrated to the target distribution.
 *
 * The histogram is considered equilibrated if, for a minimum fraction of
 * the target region, the relative error of the sampled weight relative
 * to the target is less than a tolerance value.
 *
 * \param[in] pointStates  The state of the bias points.
 * \returns true if the histogram is equilibrated.
 */
bool histogramIsEquilibrated(ArrayRef<const PointState> pointStates)
{
    /* Get the total weight of the total weight histogram; needed for normalization. */
    double totalWeight     = 0;
    int    numTargetPoints = 0;
    for (const auto& pointState : pointStates)
    {
        if (!pointState.inTargetRegion())
        {
            continue;
        }
        totalWeight += pointState.weightSumTot();
        numTargetPoints++;
    }
    GMX_RELEASE_ASSERT(totalWeight > 0, "No samples when normalizing AWH histogram.");
    double inverseTotalWeight = 1. / totalWeight;

    /* Points with target weight below a certain cutoff are ignored. */
    static const double minTargetCutoff = 0.05;
    double              minTargetWeight = 1. / numTargetPoints * minTargetCutoff;

    /* Points with error less than this tolerance pass the check.*/
    static const double errorTolerance = 0.2;

    /* Sum up weight of points that do or don't pass the check. */
    double equilibratedWeight    = 0;
    double notEquilibratedWeight = 0;
    for (const auto& pointState : pointStates)
    {
        double targetWeight  = pointState.target();
        double sampledWeight = pointState.weightSumTot() * inverseTotalWeight;

        /* Ignore these points. */
        if (!pointState.inTargetRegion() || targetWeight < minTargetWeight)
        {
            continue;
        }

        if (std::abs(sampledWeight / targetWeight - 1) > errorTolerance)
        {
            notEquilibratedWeight += targetWeight;
        }
        else
        {
            equilibratedWeight += targetWeight;
        }
    }

    /* It is enough if sampling in at least a fraction of the target region follows the target
       distribution. Boundaries will in general fail and this should be ignored (to some extent). */
    static const double minFraction = 0.8;

    return equilibratedWeight / (equilibratedWeight + notEquilibratedWeight) > minFraction;
}

} // namespace

double HistogramSize::newHistogramSize(const BiasParams&          params,
                                       double                     t,
                                       bool                       covered,
                                       ArrayRef<const PointState> pointStates,
                                       ArrayRef<double>           weightsumCovering,
                                       FILE*                      fplog)
{
    double newHistogramSize;
    if (inInitialStage_)
    {
        /* Only bother with checking equilibration if we have covered already. */
        if (equilibrateHistogram_ && covered)
        {
            /* The histogram is equilibrated at most once. */
            equilibrateHistogram_ = !histogramIsEquilibrated(pointStates);

            if (fplog != nullptr)
            {
                std::string prefix = gmx::formatString("\nawh%d:", params.biasIndex_ + 1);
                if (!equilibrateHistogram_)
                {
                    fprintf(fplog, "%s equilibrated histogram at t = %g ps.\n", prefix.c_str(), t);
                }
                else if (!havePrintedAboutCovering_)
                {
                    fprintf(fplog,
                            "%s covered but histogram not equilibrated at t = %g ps.\n",
                            prefix.c_str(),
                            t);
                    havePrintedAboutCovering_ = true; /* Just print once. */
                }
            }
        }

        /* In the initial stage, the histogram grows dynamically as a function of the number of coverings. */
        newHistogramSize = newHistogramSizeInitialStage(params, t, covered, weightsumCovering, fplog);
    }
    else
    {
        /* If not in the initial stage, the histogram grows at a linear, possibly scaled down, rate. */
        newHistogramSize = histogramSize_ + params.updateWeight * params.localWeightScaling;
    }

    return newHistogramSize;
}

void HistogramSize::setHistogramSize(double histogramSize, double weightHistogramScalingFactor)
{
    GMX_ASSERT(histogramSize > 0, "The histogram should not be empty");
    GMX_ASSERT(weightHistogramScalingFactor > 0, "The histogram scaling factor should be positive");

    histogramSize_ = histogramSize;

    /* The weight of new samples relative to previous ones change
     * when the histogram is rescaled. We keep the log since this number
     * can become very large.
     */
    logScaledSampleWeight_ -= std::log(weightHistogramScalingFactor);
};

void HistogramSize::restoreFromHistory(const AwhBiasStateHistory& stateHistory)
{
    numUpdates_               = stateHistory.numUpdates;
    histogramSize_            = stateHistory.histSize;
    inInitialStage_           = stateHistory.in_initial;
    equilibrateHistogram_     = stateHistory.equilibrateHistogram;
    logScaledSampleWeight_    = stateHistory.logScaledSampleWeight;
    maxLogScaledSampleWeight_ = stateHistory.maxLogScaledSampleWeight;
    havePrintedAboutCovering_ = false;
}

void HistogramSize::storeState(AwhBiasStateHistory* stateHistory) const
{
    stateHistory->numUpdates               = numUpdates_;
    stateHistory->histSize                 = histogramSize_;
    stateHistory->in_initial               = inInitialStage_;
    stateHistory->equilibrateHistogram     = equilibrateHistogram_;
    stateHistory->logScaledSampleWeight    = logScaledSampleWeight_;
    stateHistory->maxLogScaledSampleWeight = maxLogScaledSampleWeight_;
    /* We'll print again about covering when restoring the state */
}

} // namespace gmx
