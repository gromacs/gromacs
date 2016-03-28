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

/*! \libinternal \file
 *
 * \brief
 * Declares and the BiasParams and DimParams structs.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_BIASPARAMS_H
#define GMX_AWH_BIASPARAMS_H

/*! \cond INTERNAL */

class Grid;

//! Constant parameters for each dimension of the coordinate.
struct DimParams
{
    /*! \brief
     * Constructor.
     *
     * \param[in] awhDimParams  AWH input parameters for this dim.
     * \param[in] pullParams    The pull parameters.
     * \param[in] beta          1/(k_B T).
     */
    DimParams(const awh_dim_params_t &awhDimParams,
              const pull_params_t    &pullParams,
              double                  beta);

    double                  betak;                    /**< Inverse variance (1/nm^2) for each coordinate dimension. */
    double                  k;                        /**< Force constant (kJ/mol/nm^2) for each coordinate dimension. */
    int                     pullCoordIndex;           /**< Indices of the pull coordinates for each coordinate dimension. */
    double                  userCoordUnitsToInternal; /**< Conversion factor coordinate units. */

    /*! \brief Convert internal coordinate units to external, user coordinate units.
     *
     * \param[in] value               Value to convert.
     * \returns the converted value.
     */
    double scaleInternalToUserInput(double value) const
    {
        return value/userCoordUnitsToInternal;
    }

    /*! \brief Convert external, user coordinate units to internal coordinate units.
     *
     * \param[in] value               Value to convert.
     * \returns the converted value.
     */
    double scaleUserInputToInternal(double value) const
    {
        return value*userCoordUnitsToInternal;
    }
};

//! Constant parameters for a Bias.
struct BiasParams
{
    /*! \brief Constructor.
     */
    BiasParams(const awh_params_t           &awhParams, 
               const awh_bias_params_t      &awhBiasParams,
               const std::vector<DimParams> &dimParams,
               double                        beta,
               double                        delta_t,
               const t_commrec              *cr,
               const Grid                   *grid);

    double    invBeta;               /**< 1/beta = kT */
    int       nstsample_coord;       /**< Number of steps per coordinate value sample. */
    int       nstupdate_free_energy; /**< Number of steps per free energy update. */
    int       nstupdate_target;      /**< Number of steps per updating the target distribution. */
    int       eTarget;               /**< Type of target distribution. */
    double    target_param;          /**< Target distribution parameter (meaning depends on eTarget). */
    bool      idealWeighthistUpdate; /**< Update reference weighthistogram using the target distribution? Otherwise use the realized distribution. */
    double    update_weight;         /**< The probability weight accumulated for each update. */
    double    localWeightScaling;    /**< Scaling factor applied to a sample before adding it to the reference weight histogram (= 1, usually). */
    double    histSizeInitial;       /**< Initial reference weight histogram size. */
    int       numSharedUpdate;       /**< The number of (multi-)simulations sharing the bias update */
    awh_ivec  coverRadius;           /**< The radius (in points) that needs to be sampled around a point before it is considered covered. */

    /*! \brief
     * Check if the AWH constant parameters permit skipping updates.
     *
     * Generally, we can skip (postpone) updates of points that are non-local
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
    bool canSkipUpdates() const
    {
        return localWeightScaling == 1;
    }
};

/*! \endcond */

#endif  /* GMX_AWH_BIASPARAMS_H */
