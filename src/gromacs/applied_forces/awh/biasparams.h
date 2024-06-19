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
 * Declares the BiasParams class.
 *
 * This class holds the parameters for the bias. Most are direct copies
 * of the input that the user provided. Some are a combination of user
 * input and properties of the simulated system.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_BIASPARAMS_H
#define GMX_AWH_BIASPARAMS_H

#include <cstdint>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

#include "dimparams.h"

namespace gmx
{

template<typename>
class ArrayRef;
class AwhBiasParams;
class AwhParams;
struct DimParams;
class GridAxis;
enum class AwhTargetType : int;

/*! \internal \brief Constant parameters for the bias.
 */
class BiasParams
{
public:
    /*! \brief Switch to turn off update skips, useful for testing.
     */
    enum class DisableUpdateSkips
    {
        no, /**< Allow update skips (when supported by the method) */
        yes /**< Disable update skips */
    };

    /*! \brief
     * Check if the parameters permit skipping updates.
     *
     * Generally, we can skip updates of points that are non-local
     * at the time of the update if we for later times, when the points
     * with skipped updates have become local, know exactly how to apply
     * the previous updates. The free energy updates only depend
     * on local sampling, but the histogram rescaling factors
     * generally depend on the histogram size (all samples).
     * If the histogram size is kept constant or the scaling factors
     * are trivial, this is not a problem. However, if the histogram growth
     * is scaled down by some factor the size at the time of the update
     * needs to be known. It would be fairly simple to, for a deterministically
     * growing histogram, backtrack and calculate this value, but currently
     * we just disallow this case. This is not a restriction because it
     * only affects the local Boltzmann target type for which every update
     * is currently anyway global because the target is always updated globally.
     *
     * \returns true when we can skip updates.
     */
    inline bool skipUpdates() const { return (!disableUpdateSkips_ && localWeightScaling == 1); }

    /*! \brief
     * Returns the radius that needs to be sampled around a point before it is considered covered.
     */
    inline const awh_ivec& coverRadius() const { return coverRadius_; }

    /*! \brief
     * Returns whether we should sample the coordinate.
     *
     * \param[in] step  The MD step number.
     */
    inline bool isSampleCoordStep(int64_t step) const
    {
        return (step > 0 && step % numStepsSampleCoord_ == 0);
    }

    /*! \brief
     * Returns whether we should update the free energy.
     *
     * \param[in] step  The MD step number.
     */
    inline bool isUpdateFreeEnergyStep(int64_t step) const
    {
        int stepIntervalUpdateFreeEnergy = numSamplesUpdateFreeEnergy_ * numStepsSampleCoord_;
        return (step > 0 && step % stepIntervalUpdateFreeEnergy == 0);
    }

    /*! \brief
     * Returns whether we should update the target distribution.
     *
     * \param[in] step  The MD step number.
     */
    inline bool isUpdateTargetStep(int64_t step) const { return step % numStepsUpdateTarget_ == 0; }

    /*! \brief
     * Returns if to do checks for covering in the initial stage.
     *
     * To avoid overhead due to expensive checks, we do not check
     * at every free energy update. However, if checks are
     * performed too rarely the detection of coverings will be
     * delayed, ultimately affecting free energy convergence.
     *
     * \param[in] step                  Time step.
     * \returns true at steps where checks should be performed.
     * \note  Only returns true at free energy update steps.
     */
    bool isCheckCoveringStep(int64_t step) const
    {
        return step > 0 && (step % numStepsCheckCovering_ == 0);
    }

    /*! \brief
     * Returns if to perform checks for anomalies in the histogram.
     *
     * To avoid overhead due to expensive checks, we do not check
     * at every free energy update. These checks are only used for
     * warning the user and can be made as infrequently as
     * neccessary without affecting the algorithm itself.
     *
     * \param[in] step                  Time step.
     * \returns true at steps where checks should be performed.
     * \note Only returns true at free energy update steps.
     * \todo Currently this function just calls isCheckCoveringStep but the checks could be done less frequently.
     */
    bool isCheckHistogramForAnomaliesStep(int64_t step) const { return isCheckCoveringStep(step); }

    /*! \brief Constructor.
     *
     * The local Boltzmann target distibution is defined by
     * 1) Adding the sampled weights instead of the target weights to the reference weight histogram.
     * 2) Scaling the weights of these samples by the beta scaling factor.
     * 3) Setting the target distribution equal the reference weight histogram.
     * This requires the following special update settings:
     *   localWeightScaling      = targetParam
     *   idealWeighthistUpdate   = false
     * Note: these variables could in principle be set to something else also for other target distribution types.
     * However, localWeightScaling < 1  is in general expected to give lower efficiency and, except for local Boltzmann,
     * idealWeightHistUpdate = false gives (in my experience) unstable, non-converging results.
     *
     * \param[in] awhParams              AWH parameters.
     * \param[in] awhBiasParams          Bias parameters.
     * \param[in] dimParams              Bias dimension parameters.
     * \param[in] beta                   1/(k_B T) in units of 1/(kJ/mol), should be > 0.
     * \param[in] mdTimeStep             The MD time step.
     * \param[in] numSharingSimulations  The number of simulations to share the bias across.
     * \param[in] gridAxis               The grid axes.
     * \param[in] disableUpdateSkips     If to disable update skips, useful for testing.
     * \param[in] biasIndex              Index of the bias.
     */
    BiasParams(const AwhParams&          awhParams,
               const AwhBiasParams&      awhBiasParams,
               ArrayRef<const DimParams> dimParams,
               double                    beta,
               double                    mdTimeStep,
               DisableUpdateSkips        disableUpdateSkips,
               int                       numSharingSimulations,
               ArrayRef<const GridAxis>  gridAxis,
               int                       biasIndex);

    /* Data members */
    const double invBeta; /**< 1/beta = kT in kJ/mol */
private:
    const int64_t numStepsSampleCoord_; /**< Number of steps per coordinate value sample. */
public:
    const int numSamplesUpdateFreeEnergy_; /**< Number of samples per free energy update. */
private:
    const int64_t numStepsUpdateTarget_; /**< Number of steps per updating the target distribution. */
    const int64_t numStepsCheckCovering_; /**< Number of steps per checking for covering. */
public:
    const AwhTargetType eTarget;    /**< Type of target distribution. */
    const bool scaleTargetByMetric; /**< Scale the target distribution based on the friction metric? */
    const double targetMetricScalingLimit; /**< The upper limit for the scaling factor based on the friction metric. The lower limit is the inverse. */
    const double freeEnergyCutoffInKT; /**< Free energy cut-off in kT for cut-off target distribution. */
    const double temperatureScaleFactor; /**< Temperature scaling factor for temperature scaled targed distributions. */
    const bool   idealWeighthistUpdate; /**< Update reference weighthistogram using the target distribution? Otherwise use the realized distribution. */
    const int    numSharedUpdate; /**< The number of (multi-)simulations sharing the bias update */
    const double updateWeight;    /**< The probability weight accumulated for each update. */
    const double localWeightScaling; /**< Scaling factor applied to a sample before adding it to the reference weight histogram (= 1, usually). */
    const double initialErrorInKT;   /**< Estimated initial free energy error in kT. */
    const double initialHistogramSize; /**< Initial reference weight histogram size. */
private:
    awh_ivec coverRadius_; /**< The radius (in points) that needs to be sampled around a point before it is considered covered. */
public:
    const bool convolveForce; /**< True if we convolve the force, false means use MC between umbrellas. */
    const int biasIndex_; /**< Index of the bias, used as a second random seed and for priting. */
private:
    const bool disableUpdateSkips_; /**< If true, we disallow update skips, even when the method supports it. */
};

} // namespace gmx

#endif /* GMX_AWH_BIASPARAMS_H */
