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
 *
 * \brief
 * Declares and defined the PointState class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_POINTSTATE_H
#define GMX_AWH_POINTSTATE_H

#include <cmath>

#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/utility/gmxassert.h"

#include "biasparams.h"
#include "math.h"

namespace gmx
{

//! A value that can be passed to exp() with result 0, also with SIMD
static const double c_largeNegativeExponent = -10000.0;

//! The largest acceptable positive exponent for variables that are passed to exp().
static const double c_largePositiveExponent =  700.0;

/*! \internal \brief The state of a coordinate point.
 */
class PointState
{
    public:
        /*! \brief Constructs a point state with default values. */
        PointState() : bias_(0),
                       freeEnergy_(0),
                       target_(1),
                       targetConstantWeight_(1),
                       weightsum_iteration(0),
                       weightsum_covering(0),
                       weightsumTot(0),
                       weightsumRef_(1),
                       lastUpdateIndex_(0),
                       logPmfsum_(0),
                       visits_iteration(0),
                       visits_tot(0)
        {
        };

        /*! \brief
         * Set all values in the state to those from a history.
         *
         * \param[in] psh  Coordinate point history to copy from.
         */
        void setFromHistory(const AwhPointStateHistory &psh)
        {
            target_             = psh.target;
            freeEnergy_         = psh.free_energy;
            bias_               = psh.bias;
            weightsum_iteration = psh.weightsum_iteration;
            weightsum_covering  = psh.weightsum_covering;
            weightsumTot        = psh.weightsum_tot;
            weightsumRef_       = psh.weightsum_ref;
            lastUpdateIndex_    = psh.last_update_index;
            logPmfsum_          = psh.log_pmfsum;
            visits_iteration    = psh.visits_iteration;
            visits_tot          = psh.visits_tot;
        }

        /*! \brief
         * Query if the point is in the target region.
         *
         * \returns true if the point is in the target region.
         */
        bool inTargetRegion() const
        {
            return target_ > 0;
        }

        /*! \brief Return the bias function estimate. */
        double bias() const
        {
            return bias_;
        }

        /*! \brief Set the target to zero and the bias to minus infinity. */
        void setTargetToZero()
        {
            target_ = 0;
            /* the bias = log(target) + const = -infty */
            bias_   = c_largeNegativeExponent;
        }

        /*! \brief Return the free energy. */
        double freeEnergy() const
        {
            return freeEnergy_;
        }

        /*! \brief Set the free energy.
         *
         * TODO: Replace this setter function with a more elegant solution.
         *
         * \param[in] freeEnergy  The free energy.
         */
        void setFreeEnergy(double freeEnergy)
        {
            freeEnergy_ = freeEnergy;
        }

        /*! \brief Return the target. */
        double target() const
        {
            return target_;
        }

        /*! \brief Return the reference weight histogram. */
        double weightsumRef() const
        {
            return weightsumRef_;
        }

        /*! \brief Return the last update that was performed (in units of number of updates). */
        int lastUpdateIndex() const
        {
            return lastUpdateIndex_;
        }

        /*! \brief Return log(PMFsum). */
        double logPmfsum() const
        {
            return logPmfsum_;
        }

        /*! \brief Set log(PMFsum).
         *
         * TODO: Replace this setter function with a more elegant solution.
         *
         * \param[in] logPmfsum  The log(PMFsum).
         */
        void setLogPmfsum(double logPmfsum)
        {
            logPmfsum_ = logPmfsum;
        }

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
            weightsumRef_ = histogramSize*target_;
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
                logPmfsum_  += minimumFreeEnergy;
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
         * \param[in] logPmfsumScaling      Scale factor for the reference PMF histogram.
         * \returns true if at least one update was applied.
         */
        bool updateSkipped(const BiasParams &params, int numUpdates,
                           double weighthistScaling, double logPmfsumScaling)
        {
            GMX_ASSERT(params.skipUpdates(), "Calling function for skipped updates when skipping updates is not allowed");

            if (!inTargetRegion())
            {
                return false;
            }

            /* The most current past update */
            int lastUpdateIndex   = numUpdates;
            int numUpdatesSkipped = lastUpdateIndex - lastUpdateIndex_;

            if (numUpdatesSkipped == 0)
            {
                /* Was not updated */
                return false;
            }

            for (int i = 0; i < numUpdatesSkipped; i++)
            {
                /* This point was non-local at the time of the update meaning no weight */
                update(params, 0, weighthistScaling, logPmfsumScaling);
            }

            /* Only past updates are applied here. */
            lastUpdateIndex_ = lastUpdateIndex;

            return true;
        }

        /*! \brief Apply a point update with new sampling data.
         *
         * The last update index is also updated here.
         *
         * \param[in] params              The AWH bias parameters.
         * \param[in] numUpdates          The global number of updates.
         * \param[in] weighthistScaling   Scaling factor for the reference weight histogram.
         * \param[in] logPmfsumScaling    Log of the scaling factor for the PMF histogram.
         */
        void updateNew(const BiasParams &params, int numUpdates,
                       double weighthistScaling, double logPmfsumScaling)
        {
            GMX_RELEASE_ASSERT(lastUpdateIndex_ == numUpdates, "When doing a normal update, the point update index should match the global index, otherwise we lost (skipped?) updates.");

            update(params, weightsum_iteration, weighthistScaling, logPmfsumScaling);
            lastUpdateIndex_ += 1;
        }


        /*! \brief Update the PMF histogram with the current coordinate value.
         *
         * \param[in] convolvedBias  The convolved bias.
         */
        void samplePmf(double convolvedBias)
        {
            if (inTargetRegion())
            {
                logPmfsum_        = expsum(logPmfsum_, -convolvedBias);
                visits_iteration += 1;
            }
        }

    private:
        /*! \brief Update the free energy estimate of a point.
         *
         * The free energy update here is inherently local, i.e. it just depends on local sampling and on constant
         * AWH parameters. This assumes that the variables used here are kept constant, at least in between
         * global updates.
         *
         * \param[in] params          The AWH bias parameters.
         * \param[in] weightAtPoint   Sampled probability weight at this point.
         */
        void updateFreeEnergy(const BiasParams &params, double weightAtPoint)
        {
            double weighthistSampled  = weightsumRef() + weightAtPoint;
            double weighthistTarget   = weightsumRef() + params.updateWeight*target_;

            double df                 = -std::log(weighthistSampled/weighthistTarget);
            freeEnergy_              += df;

            GMX_RELEASE_ASSERT(std::abs(freeEnergy_) < c_largePositiveExponent,
                               "Very large free energy differences or badly normalized free energy in AWH update.");
        }

        /*! \brief Update the reference weight histogram of a point.
         *
         * \param[in] params         The AWH bias parameters.
         * \param[in] weightAtPoint  Sampled probability weight at this point.
         * \param[in] scaleFactor    Factor to rescale the histogram with.
         */
        void updateWeightHistogram(const BiasParams &params, double weightAtPoint, double scaleFactor)
        {
            if (params.idealWeighthistUpdate)
            {
                /* Grow histogram using the target distribution. */
                weightsumRef_ += target_*params.updateWeight*params.localWeightScaling;
            }
            else
            {
                /* Grow using the actual samples (which are distributed ~ as target). */
                weightsumRef_ += weightAtPoint*params.localWeightScaling;
            }

            weightsumRef_ *= scaleFactor;
        }

        /*! \brief Apply a point update.
         *
         * This updates local properties that can be updated without
         * accessing or affecting all points.
         * This excludes updating the size of reference weight histogram and
         * the target distribution. The bias update is excluded only because if updates
         * have been skipped this function will be called multiple times, while the bias
         * only needs to be updated once (last).
         *
         * Since this function only performs the update with the given arguments and does
         * not know anything about the time of the update, the last update index is not
         * updated here. The caller should take care of updating the update index
         *
         * \param[in] params             The AWH bias parameters.
         * \param[in] weightAtPoint      Sampled probability weight at this point.
         * \param[in] weighthistScaling  Scaling factor for the reference weight histogram.
         * \param[in] logPmfsumScaling   Log of the scaling factor for the PMF histogram.
         */
        void update(const BiasParams &params,
                    double weightAtPoint, double weighthistScaling, double logPmfsumScaling)
        {
            updateFreeEnergy(params, weightAtPoint);
            updateWeightHistogram(params, weightAtPoint, weighthistScaling);
            logPmfsum_ += logPmfsumScaling;
        }


    public:
        /*! \brief Update the target wight of a point.
         *
         * Note that renormalization over all points is needed after the update.
         *
         * \param[in] params            The AWH bias parameters.
         * \param[in] freeEnergyCutoff  The cut-off for the free energy for target type "cutoff".
         * \returns the updated value of the target.
         */
        double updateTarget(const BiasParams &params, double freeEnergyCutoff)
        {
            switch (params.eTarget)
            {
                case eawhtargetCONSTANT:
                    target_   = 1;
                    break;
                case eawhtargetCUTOFF:
                {
                    double df = freeEnergy_ - freeEnergyCutoff;
                    target_   = 1/(1 + std::exp(df));
                    break;
                }
                case eawhtargetBOLTZMANN:
                    target_   = std::exp(-params.targetParam*freeEnergy_);
                    break;
                case eawhtargetLOCALBOLTZMANN:
                    target_   = weightsumRef_;
                    break;
            }

            /* All target types can be modulated by a constant factor. */
            target_ *= targetConstantWeight_;

            return target_;
        }

        /*! \brief Scale the target weight of the point.
         *
         * \param[in] scaleFactor  Factor to scale with.
         */
        void scaleTarget(double scaleFactor)
        {
            target_ *= scaleFactor;
        }

    private:
        double bias_;                  /**< Current biasing function estimate */
        double freeEnergy_;            /**< Current estimate of the convolved free energy/PMF. */
        double target_;                /**< Current target distribution, normalized to 1 */
        double targetConstantWeight_;  /**< Constant target weight, from user data. */
    public:
        double weightsum_iteration;    /**< Accumulated weight this iteration. */
        double weightsum_covering;     /**< Accumulated weights for covering checks */
        double weightsumTot;           /**< Accumulated weights, never reset */
    private:
        double weightsumRef_;          /**< The reference weight histogram determining the free energy updates */
        int    lastUpdateIndex_;       /**< The last update that was performed at this point (in units of number of updates). */
        double logPmfsum_;             /**< Logarithm of the PMF histogram (for 1 replica) */
    public:
        double visits_iteration;       /**< Visits to this bin this iteration. */
        double visits_tot;             /**< Accumulated visits to this bin */
};

}      // namespace gmx

#endif /* GMX_AWH_POINTSTATE_H */
