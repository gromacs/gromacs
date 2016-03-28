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

struct AwhBiasHistory;
struct AwhHistory;
struct AwhPointStateHistory;
struct awh_bias_params_t;
struct awh_dim_params_t;
struct awh_params_t;
struct gmx_multisim_t;
class Grid;
struct GridAxis;
class PointState;
struct t_commrec;

/*! \internal \brief Constant parameters for each dimension of the coordinate.
 */
struct DimParams
{
    /*! \brief
     * Constructor.
     *
     * \param[in] conversionFactor  Conversion factor from user coordinate units to bias internal units (=DEG2RAD for angles).
     * \param[in] forceConstant     The harmonic force constant.
     * \param[in] beta              1/(k_B T).
     */
    DimParams(double conversionFactor,
              double forceConstant,
              double beta);

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

    const double k;                        /**< Force constant (kJ/mol/nm^2) for each coordinate dimension. */
    const double betak;                    /**< Inverse variance (1/nm^2) for each coordinate dimension. */
    const double userCoordUnitsToInternal; /**< Conversion factor coordinate units. */
};

/*! \internal \brief Constant parameters for the method.
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
         */
        BiasParams(const awh_params_t           &awhParams,
                   const awh_bias_params_t      &awhBiasParams,
                   const std::vector<DimParams> &dimParams,
                   double                        beta,
                   double                        mdTimeStep,
                   DisableUpdateSkips            disableUpdateSkips,
                   const t_commrec              *cr,
                   const std::vector<GridAxis>  &gridAxis);

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
    private:
        DisableUpdateSkips disableUpdateSkips_;        /**< If true, we disallow update skips, even when the method supports it. */
};

/*! \internal \brief Variables defining the global state of a bias.
 */
struct BiasState
{
    awh_dvec  coordValue;             /**< Current coordinate value in (nm or rad) */
    int       gridpointIndex;         /**< The grid point index for the current coordinate value */
    int       refGridpoint;           /**< Index for the current reference grid point (for umbrella potential type) */
    awh_ivec  originUpdatelist;       /**< The origin of the subgrid that has been touched since last update. */
    awh_ivec  endUpdatelist;          /**< The end of the subgrid that has been touched since last update. */
    bool      inInitialStage;         /**< True if in the intial stage. */
    bool      equilibrateHistogram;   /**< True if samples are kept from accumulating until the sampled distribution is close enough to the target. */
    double    histSize;               /**< Size of reference weight histogram. */
    double    scaledSampleWeight;     /**< The log of the current sample weight, scaled because of the histogram rescaling. */
    double    maxScaledSampleWeight;  /**< Maximum sample weight obtained for previous (smaller) histogram sizes. */
    int       numUpdates;             /**< The number of updates performed since the start of the simulation. */
};

/*! \internal \brief An AWH bias, i.e. all high-level data and methods associated with biasing reaction coordinates that contribute to one histogram.
 */
class Bias
{
    public:
        /*! \brief
         * Constructor.
         *
         * \param[in,out] fplog              Log file.
         * \param[in] cr                     Struct for communication.
         * \param[in] biasIndexInCollection  Index of the bias in collection.
         * \param[in] awhParams              AWH parameters.
         * \param[in] awhBiasParams          Bias parameters.
         * \param[in] dimParams              Bias dimension parameters.
         * \param[in] beta                   1/(k_B T).
         * \param[in] mdTimeStep             The MD time step.
         * \param[in] disableUpdateSkips     If to disable update skips, useful for testing
         */
        Bias(FILE                          *fplog,
             const t_commrec               *cr,
             int                            biasIndexInCollection,
             const awh_params_t            &awhParams,
             const awh_bias_params_t       &awhBiasParams,
             const std::vector<DimParams>  &dimParams,
             double                         beta,
             double                         mdTimeStep,
             BiasParams::DisableUpdateSkips disableUpdateSkips = BiasParams::DisableUpdateSkips::no);

        /*! \brief
         * Update the coordinate value of dimension \p dim.
         *
         * \param[in] dim         The dimension.
         * \param[in] coordValue  The coordinate value.
         */
        void setCoordValue(int dim, double coordValue);

        /*! \brief
         * Do what the AWH bias needs to do every step.
         *
         * The main objectives of the AWH bias step is to 1) set the bias force
         * and potential; 2) update the free energy and bias; and 3) sample
         * observables (e.g. the PMF).
         *
         * \param[out]    biasForce      The AWH bias force.
         * \param[in,out] awhPotential   Bias potential from all AWH biases.
         * \param[in,out] potential_jump Change in bias potential for all AWH biases.
         * \param[in] ms                 Struct for multi-simulation communication.
         * \param[in] t                  Time.
         * \param[in] step               Time step.
         * \param[in] seed               Random seed.
         * \param[in,out] fplog          Log file.
         */
        void doStep(awh_dvec biasForce,
                    double *awhPotential, double *potential_jump,
                    const gmx_multisim_t *ms,
                    double t, gmx_int64_t step, int seed,
                    FILE *fplog);

        /*! \brief
         * Update the bias history with a new state.
         *
         * \param[in,out] biasHistory  Bias history struct.
         */
        void updateHistory(AwhBiasHistory *biasHistory) const;

        /*! \brief
         * Restore the bias state from history.
         *
         * \param[in] biasHistory  Bias history struct.
         */
        void restoreStateFromHistory(const AwhBiasHistory *biasHistory);

        /*! \brief
         * Broadcast the bias data over the MPI ranks in this simulation.
         *
         * \param[in] cr  Struct for communication.
         */
        void broadcast(const t_commrec *cr);

        //! Returns the dimensionality of the bias.
        int ndim() const
        {
            return dimParams_.size();
        }

        /*! \brief Returns the dimension parameters for a dimension.
         *
         * \param[in] dim  The dimension.
         */
        const DimParams &dimParams(int dim) const
        {
            return dimParams_[dim];
        }

        //! Returns the grid
        const Grid *grid() const
        {
            return grid_.get();
        }

        //! Returns the bias parameters
        const BiasParams &params() const
        {
            return params_;
        }

        //! Returns the point state.
        const std::vector<PointState> &pointState() const
        {
            return pointState_;
        }

        //! Returns the global state of the bias.
        const BiasState &state() const
        {
            return state_;
        }

    private:
        /*! \brief
         * Performs an update of the bias.
         *
         * The objective of the update is to use collected samples (probability weights)
         * to improve the free energy estimate. For sake of efficiency, the update is
         * local whenever possible meaning that only points that have actually been sampled
         * are accessed and updated here. For certain AWH settings or at certain steps
         * however, global need to be performed. Besides the actual free energy update, this
         * function takes care of ensuring future convergence of the free energy. Convergence
         * is obtained by increasing the size of the reference weight histogram in a controlled
         * (sometimes dynamic) manner. Also, there are AWH variables that are direct functions
         * of the free energy or sampling history that need to be updated here, namely the target
         * distribution and the bias function.
         *
         * \param[in] ms         Struct for multi-simulation communication.
         * \param[in] t          Time.
         * \param[in] step       Time step.
         * \param[in,out] fplog  Log file.
         */
        void updateBias(const gmx_multisim_t *ms, double t, gmx_int64_t step, FILE *fplog);

        /*! \brief
         * Move the center point of the umbrella potential.
         *
         * A new umbrella center is sampled from the biased distibution. Also, the bias
         * force is updated and the new potential is return.
         *
         * This function should only be called when the bias force is not being convolved.
         * It is assumed that the probability distribution has been updated.
         *
         * \param[in] probWeightNeighbor  Probability weights of the neighbors.
         * \param[in,out] biasForce       The AWH bias force.
         * \param[in] step                Step number, needed for the random number generator.
         * \param[in] seed                Random seed.
         * \returns the new potential value.
         */
        double moveUmbrella(const std::vector<double> &probWeightNeighbor,
                            awh_dvec biasForce,
                            gmx_int64_t step, int seed);

    public:
        /*! \brief
         * Do all previously skipped updates.
         * Public for use by tests.
         */
        void doSkippedUpdatesForAllPoints();

    private:
        /*! \brief
         * Do previously skipped updates in this neighborhood.
         */
        void doSkippedUpdatesInNeighborhood();

        /*! \brief
         * Reset the range used to make the local update list.
         */
        void resetLocalUpdateRange();

        /*! \brief
         * Return the new reference weight histogram size for the current update.
         *
         * This function also takes care of checking for covering in the initial stage.
         *
         * \param[in] t          Time.
         * \param[in] covered    True if the sampling interval has been covered enough.
         * \param[in,out] fplog  Log file.
         * \returns the new histogram size.
         */
        double newHistSize(double t, bool covered, FILE *fplog);

        /*! \brief
         * Returns the new size of the reference weight histogram in the initial stage.
         *
         * This function also takes care resetting the histogram used for covering checks
         * and for exiting the initial stage.
         *
         * \param[in] t            Time.
         * \param[in] haveCovered  True if the sampling interval has been covered enough.
         * \param[in,out] fplog    Log file.
         * \returns the new histogram size.
         */
        double newHistSizeInitialStage(double t, bool haveCovered, FILE *fplog);

        /*! \brief
         * Sample observables for future updates or analysis.
         *
         * \param[in] probWeightNeighbor  Probability weights of the neighbors.
         * \param[in]     convolvedBias   The convolved bias.
         */
        void doObservableSampling(const std::vector<double> &probWeightNeighbor,
                                  double                     convolvedBias);

        /*! \brief
         * Save the current probability weights for future updates and analysis.
         *
         * Points in the current neighborhood will now have data meaning they
         * need to be included in the local update list of the next update.
         * Therefore, the local update range is also update here.
         *
         * \param[in] probWeightNeighbor  Probability weights of the neighbors.
         */
        void sampleProbabilityWeights(const std::vector<double> &probWeightNeighbor);

        /* Data members. */
    public:
        const int                    biasIndex;    /**< The index of this bias in AwhBiasCollection */
    private:
        const std::vector<DimParams> dimParams_;   /**< Parameters for each dimension. */
        std::unique_ptr<Grid>        grid_;        /**< The multidimensional grid organizing the coordinate point locations. */

        BiasParams                   params_;      /**< Constant parameters for the method. */

        BiasState                    state_;       /**< The global state. */
        std::vector<PointState>      pointState_;  /**< Grid point states. */
        std::vector<int>             updateList_;  /**< List of points for update for temporary use. */

        /* Temporary working vector used during the update.
         * This only here to avoid allocation at every MD step.
         */
        std::vector<double>          tempWorkSpace_;  /**< Working vector of doubles. */
};

#endif  /* GMX_AWH_BIAS_H */
