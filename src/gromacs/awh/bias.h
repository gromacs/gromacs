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
 * Declares Bias class and its helpers.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_BIAS_H
#define GMX_AWH_BIAS_H

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

//! The maximum dimensionality of the AWH coordinate.
static const int c_biasMaxNumDim = 4;

//! A real vector in AWH coordinate space.
typedef double awh_dvec[c_biasMaxNumDim];

//! An integer vector in AWH coordinate space.
typedef int awh_ivec[c_biasMaxNumDim];

struct AwhHistory;
class AwhPointStateHistory;
struct awh_params_t;
struct BiasParams;
struct BiasState;
struct gmx_multisim_t;
struct t_commrec;
class Grid;
struct t_inputrec;
struct t_mdatoms;

/*! \cond INTERNAL */

//! A value that can be passed to exp() with result 0, also with SIMD
static const double c_largeNegativeExponent = -10000.0;

//! The largest acceptable positive exponent for variables that are passed to exp().
static const double c_largePositiveExponent =  700.0;

//! The state of a point in the AWH grid.
class PointState
{
    private:
        double bias_;                  /**< Current biasing function estimate */
        double freeEnergy_;            /**< Current estimate of the convolved free energy/PMF. */
    public:
        double target;                 /**< Current target distribution, normalized to 1 */
        double target_constant_weight; /**< Constant target weight, from user data. */
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

        /*! \brief Constructor. */
        PointState();

        /*! \brief
         * Set all values in the state to those from a history.
         *
         * \param[in] cph  Coordinate point history to copy from.
         */
        void setFromHistory(const AwhPointStateHistory &psh);

        /*! \brief
         * Query if the grid point is in the target region.
         *
         * \returns true if the point is in the target region.
         */
        bool inTargetRegion() const
        {
            return target > 0;
        }

        /*! \brief Return the bias function estimate. */
        double bias() const
        {
            return bias_;
        }

        /*! \brief Set the bias to minus infinity. */
        void setBiasToMinusInfinity()
        {
            bias_ = c_largeNegativeExponent;
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
        /*! \brief Updates the bias of a grid point. */
        void updateBias();

        /*! \brief Set the initial reference weighthistogram.
         *
         * \param[in] histogramSize  Size of reference weight histogram.
         */
        void setInitialReferenceWeightHistogram(const BiasState &biasState);

        /*! \brief Correct free energy and PMF sum for the change in minimum.
         *
         * \param[in] minimumFreeEnergy  The free energy at the minimum;
         */
        void normalizeFreeEnergyAndPmfSum(double minimumFreeEnergy);

        /*! \brief Apply previous updates that were skipped.
         *
         * The last update index is also updated here.
         *
         * \param[in] params             The AWH bias parameters.
         * \param[in] step               Current time step.
         * \param[in] weighthistScaling  Scale factor for the reference weight histogram.
         * \param[in] logPmfsumScaling      Scale factor for the reference PMF histogram.
         * \returns true if at least one update was applied.
         */
        bool updateSkipped(const BiasParams &params,
                           gmx_int64_t step, double weighthistScaling, double logPmfsumScaling);

        /*! \brief Apply a point update with new sampling data.
         *
         * The last update index is also updated here.
         *
         * \param[in] params              The AWH bias parameters.
         * \param[in] step                Time step for updating the last update index.
         * \param[in] weighthistScaling   Scaling factor for the reference weight histogram.
         * \param[in] logPmfsumScaling    Log of the scaling factor for the PMF histogram.
         */
        void updateNew(const BiasParams &params,
                       gmx_int64_t step, double weighthistScaling, double logPmfsumScaling);

        /*! \brief Update the PMF histogram with the current coordinate value.
         *
         * \param[in] convolvedBias  The convolved bias.
         */
        void samplePmf(double convolvedBias);

    private:
        /*! \brief Update the free energy estimate of a point.
         *
         * \param[in] params          The AWH bias parameters.
         * \param[in] weightAtPoint   Sampled probability weight at this point.
         */
        void updateFreeEnergy(const BiasParams &params, double weightAtPoint);

        /*! \brief Update the reference weight histogram of a point.
         *
         * \param[in] params         The AWH bias parameters.
         * \param[in] weightAtPoint  Sampled probability weight at this point.
         * \param[in] scaleFactor    Factor to rescale the histogram with.
         */
        void updateWeightHistogram(const BiasParams &params, double weightAtPoint, double scaleFactor);

        /*! \brief Apply a point update.
         *
         * This updates local properties that can be updated without
         * accessing or affecting all points.
         *
         * \param[in] params             The AWH bias parameters.
         * \param[in] weightAtPoint      Sampled probability weight at this point.
         * \param[in] weighthistScaling  Scaling factor for the reference weight histogram.
         * \param[in] logPmfsumScaling   Log of the scaling factor for the PMF histogram.
         */
        void update(const BiasParams &params,
                    double weightAtPoint, double weighthistScaling, double logPmfsumScaling);
};

//! Constant parameters for each dimension of the coordinate.
struct DimParams
{
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

//! Constant parameters for the method.
struct BiasParams
{
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

//! Variables defining the global state of a bias.
struct BiasState
{
    awh_dvec  coordValue;             /**< Current coordinate value in (nm or rad) */
    int       gridpointIndex;         /**< The grid point index for the current coordinate value */
    int       refGridpoint;           /**< Index for the current reference grid point (for umbrella potential type) */
    awh_ivec  origin_updatelist;      /**< The origin of the subgrid that has been touched since last update. */
    awh_ivec  end_updatelist;         /**< The end of the subgrid that has been touched since last update. */
    bool      in_initial;             /**< True if in the intial stage. */
    bool      equilibrateHistogram;   /**< True if samples are kept from accumulating until the sampled distribution is close enough to the target. */
    double    histSize;               /**< Size of reference weight histogram. */
    double    scaledSampleWeight;     /**< The log of the current sample weight, scaled because of the histogram rescaling. */
    double    maxScaledSampleWeight;  /**< Maximum sample weight obtained for previous (smaller) histogram sizes. */
};

//! An AWH bias, i.e. all data associated with biasing reaction coordinates that contribute to one histogram.
struct Bias
{
    int                     biasIndex;                  /**< The index of this bias in AwhBiasCollection */
    int                     ndim;                       /**< Dimensionality of the AWH coordinate. */
    std::unique_ptr<Grid>   grid;                       /**< The multidimensional grid organizing the coordinate point locations. */
    std::vector<PointState> pointState;                 /**< Grid point states. */
    std::vector<int>        updateList;                 /**< List of points for update for temporary use. */

    DimParams               dimParams[c_biasMaxNumDim]; /**< Parameters for each dimension. */

    BiasParams              params;                     /**< Constant parameters for the method. */

    BiasState               state;                      /**< The global state. */
};

/*! \endcond */

#endif  /* GMX_AWH_BIAS_H */
