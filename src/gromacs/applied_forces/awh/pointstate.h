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
 *
 * \brief
 * Declares and defines the PointState class.
 *
 * Since nearly all operations on PointState objects occur in loops over
 * (parts of) the grid of an AWH bias, all these methods should be inlined.
 * Only samplePmf() is called only once per step and is thus not inlined.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_POINTSTATE_H
#define GMX_AWH_POINTSTATE_H

#include <cmath>

#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "biasparams.h"

namespace gmx
{

namespace detail
{

//! A value that can be passed to exp() with result 0, also with SIMD
constexpr double c_largeNegativeExponent = -10000.0;

//! The largest acceptable positive exponent for variables that are passed to exp().
constexpr double c_largePositiveExponent = 700.0;

} // namespace detail

/*! \internal
 * \brief The state of a coordinate point.
 *
 * This class contains all the state variables of a coordinate point
 * (on the bias grid) and methods to update the state of a point.
 */
class PointState
{
public:
    /*! \brief Constructs a point state with default values. */
    PointState() :
        bias_(0),
        freeEnergy_(0),
        target_(1),
        targetConstantWeight_(1),
        weightSumIteration_(0),
        weightSumTot_(0),
        weightSumRef_(1),
        lastUpdateIndex_(0),
        logPmfSum_(0),
        numVisitsIteration_(0),
        numVisitsTot_(0),
        localWeightSum_(0)
    {
    }

    /*! \brief
     * Set all values in the state to those from a history.
     *
     * \param[in] psh  Coordinate point history to copy from.
     */
    void setFromHistory(const AwhPointStateHistory& psh)
    {
        target_             = psh.target;
        freeEnergy_         = psh.free_energy;
        bias_               = psh.bias;
        weightSumIteration_ = psh.weightsum_iteration;
        weightSumTot_       = psh.weightsum_tot;
        weightSumRef_       = psh.weightsum_ref;
        lastUpdateIndex_    = psh.last_update_index;
        logPmfSum_          = psh.log_pmfsum;
        numVisitsIteration_ = psh.visits_iteration;
        numVisitsTot_       = psh.visits_tot;
        localWeightSum_     = psh.localWeightSum;
    }

    /*! \brief
     * Store the state of a point in a history struct.
     *
     * \param[in,out] psh  Coordinate point history to copy to.
     */
    void storeState(AwhPointStateHistory* psh) const
    {
        psh->target              = target_;
        psh->free_energy         = freeEnergy_;
        psh->bias                = bias_;
        psh->weightsum_iteration = weightSumIteration_;
        psh->weightsum_tot       = weightSumTot_;
        psh->weightsum_ref       = weightSumRef_;
        psh->last_update_index   = lastUpdateIndex_;
        psh->log_pmfsum          = logPmfSum_;
        psh->visits_iteration    = numVisitsIteration_;
        psh->visits_tot          = numVisitsTot_;
        psh->localWeightSum      = localWeightSum_;
    }

    /*! \brief
     * Query if the point is in the target region.
     *
     * \returns true if the point is in the target region.
     */
    bool inTargetRegion() const { return target_ > 0; }

    /*! \brief Return the bias function estimate. */
    double bias() const { return bias_; }

    /*! \brief Set the target to zero and the bias to minus infinity. */
    void setTargetToZero()
    {
        target_ = 0;
        /* the bias = log(target) + const = -infty */
        bias_ = detail::c_largeNegativeExponent;
    }

    /*! \brief Return the free energy. */
    double freeEnergy() const { return freeEnergy_; }

    /*! \brief Set the free energy, only to be used at initialization.
     *
     * \param[in] freeEnergy  The free energy.
     */
    void setFreeEnergy(double freeEnergy) { freeEnergy_ = freeEnergy; }

    /*! \brief Return the target distribution value. */
    double target() const { return target_; }

    /*! \brief Return the weight accumulated since the last update. */
    double weightSumIteration() const { return weightSumIteration_; }

    /*! \brief Increases the weight accumulated since the last update.
     *
     * \param[in] weight  The amount to add to the weight
     */
    void increaseWeightSumIteration(double weight) { weightSumIteration_ += weight; }

    /*! \brief Returns the accumulated weight */
    double weightSumTot() const { return weightSumTot_; }

    /*! \brief Return the reference weight histogram. */
    double weightSumRef() const { return weightSumRef_; }

    /*! \brief Return log(PmfSum). */
    double logPmfSum() const { return logPmfSum_; }

    /*! \brief Set log(PmfSum).
     *
     * TODO: Replace this setter function with a more elegant solution.
     *
     * \param[in] logPmfSum  The log(PmfSum).
     */
    void setLogPmfSum(double logPmfSum) { logPmfSum_ = logPmfSum; }

    /*! \brief Return the number of visits since the last update */
    double numVisitsIteration() const { return numVisitsIteration_; }

    /*! \brief Return the total number of visits */
    double numVisitsTot() const { return numVisitsTot_; }

    /*! \brief Return the local contribution to the accumulated weight */
    double localWeightSum() const { return localWeightSum_; }

    /*! \brief Set the constant target weight factor.
     *
     * \param[in] targetConstantWeight  The target weight factor.
     */
    void setTargetConstantWeight(double targetConstantWeight)
    {
        targetConstantWeight_ = targetConstantWeight;
    }

    /*! \brief Updates the bias of a point. */
    void updateBias()
    {
        GMX_ASSERT(target_ > 0, "AWH target distribution must be > 0 to calculate the point bias.");

        bias_ = freeEnergy() + std::log(target_);
    }

    /*! \brief Set the initial reference weighthistogram.
     *
     * \param[in] histogramSize  The weight histogram size.
     */
    void setInitialReferenceWeightHistogram(double histogramSize)
    {
        weightSumRef_ = histogramSize * target_;
    }

    /*! \brief Correct free energy and PMF sum for the change in minimum.
     *
     * \param[in] minimumFreeEnergy  The free energy at the minimum;
     */
    void normalizeFreeEnergyAndPmfSum(double minimumFreeEnergy)
    {
        if (inTargetRegion())
        {
            /* The sign of the free energy and PMF constants are opposite
             * because the PMF samples are reweighted with the negative
             * bias e^(-bias) ~ e^(-free energy).
             */
            freeEnergy_ -= minimumFreeEnergy;
            logPmfSum_ += minimumFreeEnergy;
        }
    }

    /*! \brief Apply previous updates that were skipped.
     *
     * An update can only be skipped if the parameters needed for the update are constant or
     * deterministic so that the same update can be performed at a later time.
     * Here, the necessary parameters are the sampled weight and scaling factors for the
     * histograms. The scaling factors are provided as arguments only to avoid recalculating
     * them for each point
     *
     * The last update index is also updated here.
     *
     * \param[in] params             The AWH bias parameters.
     * \param[in] numUpdates         The global number of updates.
     * \param[in] weighthistScaling  Scale factor for the reference weight histogram.
     * \param[in] logPmfSumScaling   Scale factor for the reference PMF histogram.
     * \returns true if at least one update was applied.
     */
    bool performPreviouslySkippedUpdates(const BiasParams& params,
                                         int64_t           numUpdates,
                                         double            weighthistScaling,
                                         double            logPmfSumScaling)
    {
        GMX_ASSERT(params.skipUpdates(),
                   "Calling function for skipped updates when skipping updates is not allowed");

        if (!inTargetRegion())
        {
            return false;
        }

        /* The most current past update */
        int64_t lastUpdateIndex   = numUpdates;
        int64_t numUpdatesSkipped = lastUpdateIndex - lastUpdateIndex_;

        if (numUpdatesSkipped == 0)
        {
            /* Was not updated */
            return false;
        }

        for (int64_t i = 0; i < numUpdatesSkipped; i++)
        {
            /* This point was non-local at the time of the update meaning no weight */
            updateFreeEnergyAndWeight(params, 0, weighthistScaling, logPmfSumScaling);
        }

        /* Only past updates are applied here. */
        lastUpdateIndex_ = lastUpdateIndex;

        return true;
    }

    /*! \brief Apply a point update with new sampling.
     *
     * \note The last update index is also updated here.
     * \note The new sampling containers are cleared here.
     *
     * \param[in] params              The AWH bias parameters.
     * \param[in] numUpdates          The global number of updates.
     * \param[in] weighthistScaling   Scaling factor for the reference weight histogram.
     * \param[in] logPmfSumScaling    Log of the scaling factor for the PMF histogram.
     */
    void updateWithNewSampling(const BiasParams& params, int64_t numUpdates, double weighthistScaling, double logPmfSumScaling)
    {
        GMX_RELEASE_ASSERT(lastUpdateIndex_ == numUpdates,
                           "When doing a normal update, the point update index should match the "
                           "global index, otherwise we lost (skipped?) updates.");

        updateFreeEnergyAndWeight(params, weightSumIteration_, weighthistScaling, logPmfSumScaling);
        lastUpdateIndex_ += 1;

        /* Clear the iteration collection data */
        weightSumIteration_ = 0;
        numVisitsIteration_ = 0;
    }


    /*! \brief Update the PMF histogram with the current coordinate value.
     *
     * \param[in] convolvedBias  The convolved bias.
     */
    void samplePmf(double convolvedBias);

    /*! \brief Update the PMF histogram of unvisited coordinate values
     * (along a lambda axis)
     *
     * \param[in] bias  The bias to update with.
     */
    void updatePmfUnvisited(double bias);

private:
    /*! \brief Update the free energy estimate of a point.
     *
     * The free energy update here is inherently local, i.e. it just depends on local sampling and
     * on constant AWH parameters. This assumes that the variables used here are kept constant, at
     * least in between global updates.
     *
     * \param[in] params          The AWH bias parameters.
     * \param[in] weightAtPoint   Sampled probability weight at this point.
     */
    void updateFreeEnergy(const BiasParams& params, double weightAtPoint)
    {
        double weighthistSampled = weightSumRef() + weightAtPoint;
        double weighthistTarget  = weightSumRef() + params.updateWeight * target_;

        double df = -std::log(weighthistSampled / weighthistTarget);
        freeEnergy_ += df;

        if (std::abs(freeEnergy_) > detail::c_largePositiveExponent)
        {
            GMX_THROW(InvalidInputError(
                    "An AWH free energy difference is larger than 700 kT, which is not supported"));
        }
    }

    /*! \brief Update the reference weight histogram of a point.
     *
     * \param[in] params         The AWH bias parameters.
     * \param[in] weightAtPoint  Sampled probability weight at this point.
     * \param[in] scaleFactor    Factor to rescale the histogram with.
     */
    void updateWeightHistogram(const BiasParams& params, double weightAtPoint, double scaleFactor)
    {
        if (params.idealWeighthistUpdate)
        {
            /* Grow histogram using the target distribution. */
            weightSumRef_ += target_ * params.updateWeight * params.localWeightScaling;
        }
        else
        {
            /* Grow using the actual samples (which are distributed ~ as target). */
            weightSumRef_ += weightAtPoint * params.localWeightScaling;
        }

        weightSumRef_ *= scaleFactor;
    }

    /*! \brief Apply a point update.
     *
     * This updates local properties that can be updated without
     * accessing or affecting all points.
     * This excludes updating the size of reference weight histogram and
     * the target distribution. The bias update is excluded only because
     * if updates have been skipped this function will be called multiple
     * times, while the bias only needs to be updated once (last).
     *
     * Since this function only performs the update with the given
     * arguments and does not know anything about the time of the update,
     * the last update index is not updated here. The caller should take
     * care of updating the update index.
     *
     * \param[in] params             The AWH bias parameters.
     * \param[in] weightAtPoint      Sampled probability weight at this point.
     * \param[in] weighthistScaling  Scaling factor for the reference weight histogram.
     * \param[in] logPmfSumScaling   Log of the scaling factor for the PMF histogram.
     */
    void updateFreeEnergyAndWeight(const BiasParams& params,
                                   double            weightAtPoint,
                                   double            weighthistScaling,
                                   double            logPmfSumScaling)
    {
        updateFreeEnergy(params, weightAtPoint);
        updateWeightHistogram(params, weightAtPoint, weighthistScaling);
        logPmfSum_ += logPmfSumScaling;
    }


public:
    /*! \brief Update the target weight of a point.
     *
     * Note that renormalization over all points is needed after the update.
     *
     * \param[in] params            The AWH bias parameters.
     * \param[in] freeEnergyCutoff  The cut-off for the free energy for target type "cutoff".
     * \returns the updated value of the target.
     */
    double updateTargetWeight(const BiasParams& params, double freeEnergyCutoff)
    {
        switch (params.eTarget)
        {
            case AwhTargetType::Constant: target_ = 1; break;
            case AwhTargetType::Cutoff:
            {
                double df = freeEnergy_ - freeEnergyCutoff;
                target_   = 1 / (1 + std::exp(df));
                break;
            }
            case AwhTargetType::Boltzmann:
                target_ = std::exp(-params.temperatureScaleFactor * freeEnergy_);
                break;
            case AwhTargetType::LocalBoltzmann: target_ = weightSumRef_; break;
            default: GMX_RELEASE_ASSERT(false, "Unhandled enum");
        }

        /* All target types can be modulated by a constant factor. */
        target_ *= targetConstantWeight_;

        return target_;
    }

    /*! \brief Set the weight and count accumulated since the last update.
     *
     * \param[in] weightSum  The weight-sum value
     * \param[in] numVisits  The number of visits
     */
    void setPartialWeightAndCount(double weightSum, double numVisits)
    {
        weightSumIteration_ = weightSum;
        numVisitsIteration_ = numVisits;
    }

    /*! \brief Add the weights and counts accumulated between updates. */
    void addPartialWeightAndCount()
    {
        weightSumTot_ += weightSumIteration_;
        numVisitsTot_ += numVisitsIteration_;
    }

    /*! \brief Add the local weight contribution accumulated between updates. */
    void addLocalWeightSum() { localWeightSum_ += weightSumIteration_; }

    /*! \brief Scale the target weight of the point.
     *
     * \param[in] scaleFactor  Factor to scale with.
     */
    void scaleTarget(double scaleFactor) { target_ *= scaleFactor; }

private:
    double bias_;                 /**< Current biasing function estimate */
    double freeEnergy_;           /**< Current estimate of the convolved free energy/PMF. */
    double target_;               /**< Current target distribution, normalized to 1 */
    double targetConstantWeight_; /**< Constant target weight, from user data. */
    double weightSumIteration_; /**< Accumulated weight this iteration; note: only contains data for this Bias, even when sharing biases. */
    double weightSumTot_;       /**< Accumulated weights, never reset or scaled. */
    double weightSumRef_; /**< The reference weight histogram determining the free energy updates */
    int64_t lastUpdateIndex_; /**< The last update that was performed at this point, in units of number of updates. */
    double  logPmfSum_;          /**< Logarithm of the PMF histogram */
    double  numVisitsIteration_; /**< Visits to this bin this iteration; note: only contains data for this Bias, even when sharing biases. */
    double  numVisitsTot_;       /**< Accumulated visits to this bin */
    double  localWeightSum_; /**< The contributed weight sum from the local Bias. This is used for computing the average shared friction metric. Never reset or scaled. */
};

} // namespace gmx

#endif /* GMX_AWH_POINTSTATE_H */
