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

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

#include "dimparams.h"

struct awh_bias_params_t;
struct awh_params_t;
struct DimParams;
struct gmx_multisim_t;
class GridAxis;
struct t_commrec;

/*! \internal \brief Constant parameters for the bias.
 */
class BiasParams
{
    public:
        /*! \brief Switch to turn off update skips, useful for testing.
         */
        enum class DisableUpdateSkips
        {
            no,  /**< Allow update skips (when supported by the method) */
            yes  /**< Disable update skips */
        };

        /*! \brief
         * Check if the parameters permit skipping updates.
         *
         * Generally, we can skip  updates of points that are non-local
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
        bool skipUpdates() const
        {
            return
                disableUpdateSkips_ == DisableUpdateSkips::no &&
                localWeightScaling == 1;
        }

        /*! \brief Constructor.
         *
         * \param[in] awhParams              AWH parameters.
         * \param[in] awhBiasParams          Bias parameters.
         * \param[in] dimParams              Bias dimension parameters.
         * \param[in] beta                   1/(k_B T).
         * \param[in] mdTimeStep             The MD time step.
         * \param[in] cr                     Struct for communication.
         * \param[in] gridAxis               The grid axes.
         * \param[in] disableUpdateSkips     If to disable update skips, useful for testing.
         * \param[in] biasIndex              Index of the bias.
         */
        BiasParams(const awh_params_t           &awhParams,
                   const awh_bias_params_t      &awhBiasParams,
                   const std::vector<DimParams> &dimParams,
                   double                        beta,
                   double                        mdTimeStep,
                   DisableUpdateSkips            disableUpdateSkips,
                   const t_commrec              *cr,
                   const std::vector<GridAxis>  &gridAxis,
                   int                           biasIndex);

        double             invBeta;                    /**< 1/beta = kT */
        int                numStepsSampleCoord;        /**< Number of steps per coordinate value sample. */
        int                numSamplesUpdateFreeEnergy; /**< Number of samples per free energy update. */
        int                nstupdate_target;           /**< Number of steps per updating the target distribution. */
        int                eTarget;                    /**< Type of target distribution. */
        double             targetParam;                /**< Target distribution parameter (meaning depends on eTarget). */
        bool               idealWeighthistUpdate;      /**< Update reference weighthistogram using the target distribution? Otherwise use the realized distribution. */
        double             update_weight;              /**< The probability weight accumulated for each update. */
        double             localWeightScaling;         /**< Scaling factor applied to a sample before adding it to the reference weight histogram (= 1, usually). */
        double             histSizeInitial;            /**< Initial reference weight histogram size. */
        int                numSharedUpdate;            /**< The number of (multi-)simulations sharing the bias update */
        awh_ivec           coverRadius;                /**< The radius (in points) that needs to be sampled around a point before it is considered covered. */
        bool               convolveForce;              /**< True if we convolve the force, false means use MC between umbrellas. */
        int                biasIndex;                  /**< Index of the bias, used as a second random seed and for priting. */
    private:
        DisableUpdateSkips disableUpdateSkips_;        /**< If true, we disallow update skips, even when the method supports it. */
};

#endif  /* GMX_AWH_BIASPARAMS_H */
