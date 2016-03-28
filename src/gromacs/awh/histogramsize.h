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
 * Declares the HistogramSize class.
 *
 * The data members of this class keep track of global size and update related
 * properties of the bias histogram and the evolution of the histogram size.
 * Initially histogramSize_ (and thus the convergence rate) is controlled
 * heuristically to get good initial estimates,  i.e. increase the robustness
 * of the method.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_HISTOGRAMSIZE_H
#define GMX_AWH_HISTOGRAMSIZE_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

struct AwhBiasStateHistory;
struct AwhBiasParams;
class BiasParams;
class PointState;

/*! \internal
 * \brief Tracks global size related properties of the bias histogram.
 *
 * Tracks the number of updates and the histogram size.
 * Also keep track of the stage (initial/final of the AWH method
 * and printing warnings about covering.
 *
 * \note Histogram sizes are floating-point values, since the histogram uses weighted
 *        entries and we can assign a floating-point scaling factor when changing it.
 */
class HistogramSize
{
    public:
        /*! \brief Constructor.
         *
         * \param[in] awhBiasParams         The Bias parameters from inputrec.
         * \param[in] histogramSizeInitial  The initial histogram size.
         */
        HistogramSize(const AwhBiasParams &awhBiasParams,
                      double               histogramSizeInitial);

    private:
        /*! \brief
         * Returns the new size of the reference weight histogram in the initial stage.
         *
         * This function also takes care resetting the histogram used for covering checks
         * and for exiting the initial stage.
         *
         * \param[in]     params             The bias parameters.
         * \param[in]     t                  Time.
         * \param[in]     detectedCovering   True if we detected that the sampling interval has been sufficiently covered.
         * \param[in,out] weightsumCovering  The weight sum for checking covering.
         * \param[in,out] fplog              Log file.
         * \returns the new histogram size.
         */
        double newHistogramSizeInitialStage(const BiasParams &params,
                                            double            t,
                                            bool              detectedCovering,
                                            ArrayRef<double>  weightsumCovering,
                                            FILE             *fplog);

    public:
        /*! \brief
         * Return the new reference weight histogram size for the current update.
         *
         * This function also takes care of checking for covering in the initial stage.
         *
         * \param[in]     params             The bias parameters.
         * \param[in]     t                  Time.
         * \param[in]     covered            True if the sampling interval has been covered enough.
         * \param[in]     pointStates        The state of the grid points.
         * \param[in,out] weightsumCovering  The weight sum for checking covering.
         * \param[in,out] fplog              Log file.
         * \returns the new histogram size.
         */
        double newHistogramSize(const BiasParams              &params,
                                double                         t,
                                bool                           covered,
                                const std::vector<PointState> &pointStates,
                                ArrayRef<double>               weightsumCovering,
                                FILE                          *fplog);

        /*! \brief Restores the histogram size from history.
         *
         * \param[in] stateHistory  The AWH bias state history.
         */
        void restoreFromHistory(const AwhBiasStateHistory &stateHistory);

        /*! \brief Store the histogram size state in a history struct.
         *
         * \param[in,out] stateHistory  The AWH bias state history.
         */
        void storeState(AwhBiasStateHistory *stateHistory) const;

        /*! \brief Returns the number of updates since the start of the simulation.
         */
        int numUpdates() const
        {
            return numUpdates_;
        };

        /*! \brief Increments the number of updates by 1.
         */
        void incrementNumUpdates()
        {
            numUpdates_ += 1;
        }

        /*! \brief Returns the histogram size.
         */
        double histogramSize() const
        {
            return histogramSize_;
        };

        /*! \brief Sets the histogram size.
         *
         * \param[in] histogramSize                 The new histogram size.
         * \param[in] weightHistogramScalingFactor  The factor to scale the weight by.
         */
        void setHistogramSize(double histogramSize,
                              double weightHistogramScalingFactor);

        /*! \brief Returns true if we are in the initial stage of the AWH method.
         */
        bool inInitialStage() const
        {
            return inInitialStage_;
        };

        /*! \brief Returns The log of the current sample weight, scaled because of the histogram rescaling.
         */
        double logScaledSampleWeight() const
        {
            return logScaledSampleWeight_;
        };

    private:
        gmx_int64_t numUpdates_; /**< The number of updates performed since the start of the simulation. */

        /* The histogram size sets the update size and so controls the convergence rate of the free energy and bias. */
        double      histogramSize_; /**< Size of reference weight histogram. */

        /* Values that control the evolution of the histogram size. */
        bool      inInitialStage_;           /**< True if in the intial stage. */
        bool      equilibrateHistogram_;     /**< True if samples are kept from accumulating until the sampled distribution is close enough to the target. */
        double    logScaledSampleWeight_;    /**< The log of the current sample weight, scaled because of the histogram rescaling. */
        double    maxLogScaledSampleWeight_; /**< Maximum sample weight obtained for previous (smaller) histogram sizes. */

        /* Bool to avoid printing multiple, not so useful, messages to log */
        bool      havePrintedAboutCovering_; /**< True if we have printed about covering to the log while equilibrateHistogram==true */
};

}      // namespace gmx

#endif /* GMX_AWH_HISTOGRAMSIZE_H */
