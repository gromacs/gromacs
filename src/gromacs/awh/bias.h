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
 * Declares the Bias class.
 *
 * This class is essentially a wrapper around the BiasState class.
 * In addition to BiasState, it holds all data that BiasState needs
 * to update the bias. Interaction of the outside world, such as updating
 * BiasState or extracting bias data all happen through Bias.
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
#include "gromacs/utility/gmxassert.h"

#include "biasparams.h"
#include "biasstate.h"
#include "dimparams.h"

struct gmx_multisim_t;
struct t_commrec;
struct t_enxsubblock;

namespace gmx
{

struct AwhBiasHistory;
struct AwhBiasParams;
struct AwhHistory;
struct AwhParams;
struct AwhPointStateHistory;
class BiasParams;
class BiasState;
class BiasWriter;
class CorrelationGrid;
class Grid;
class GridAxis;
class PointState;

/*! \internal
 * \brief A bias acting on a multidimensional coordinate.
 *
 * At each step AWH should provide its biases with updated
 * values of their coordinates. Each bias provides AWH with an updated
 * bias forces and the corresponding potential.
 *
 * The bias is responsible for keeping and updating a free energy estimate
 * along the coordinate. The bias potential is basically a function of the
 * free energy estimate and so also changes by the update.
 * The free energy update is based on information from coordinate samples
 * collected at a contant bias potential, beween updates.
 *
 * The bias keeps a grid with coordinate points that organizes spatial
 * information about the coordinate. The grid has the the same geometry
 * as the coordinate, i.e. they has the same dimensionality and periodicity
 * (if any). The number of points in the grid sets the resolution of
 * the collected data and its extent defines the sampling region of interest.
 *
 * Each coordinate point has further statistical properties and function values
 * which a grid point does not know about. E.g., for the bias each coordinate point
 * is associated with values of the bias, free energy and target distribution,
 * accumulated sampling weight, etc. For this the bias attaches to each grid
 * point a state. The grid + vector of point states are the bias coordinate points.
 *
 * The bias has a fairly complex global state keeping track of where the system
 * (coordinate) currently is, where it has sampled since the last update,
 * and controlling the free energy convergence rate.
 *
 * Partly, the complexity comes from the bias having two convergence stages: an
 * initial stage which in an heuristic, non-deterministic way restricts the early
 * convergence rate for sake of robustness; and a final stage where the convergence rate is constant.
 * The length of the initial stage depends on the sampling and is unknown beforehand.
 *
 * Another complexity comes from the fact that coordinate points, for sake of efficiency
 * in the case of many grid points, are typically only accessed in recently sampled
 * regions even though the free energy update is inherently global and affects all points.
 * The bias allows points thay are non-local at the time the update was issued to postpone
 * ("skip", as it is called in the code) the update. A non-local point is defined
 * as a point which has not been sampled since the last update. Local points are points that have
 * been sampled since the last update. The (current) set of local points are kept track of by
 * the bias state and reset after every update. An update is called local if it only updates
 * local points. Non-local points will temporarily "skip" the update until next time they are
 * local (or when a global update is issued). For this to work, the bias keeps a global "clock"
 * of the number of issued updates. Each point state also has its own local "clock" with the
 * counting the number of updates it has pulled through. When a point updates its state it
 * asserts that its local clock is synched with the global clock.
 */
class Bias
{
    public:
        /*! \brief
         * Constructor.
         *
         * \param[in] cr                     Struct for communication.
         * \param[in] biasIndexInCollection  Index of the bias in collection.
         * \param[in] awhParams              AWH parameters.
         * \param[in] awhBiasParams          Bias parameters.
         * \param[in] dimParams              Bias dimension parameters.
         * \param[in] beta                   1/(k_B T).
         * \param[in] mdTimeStep             The MD time step.
         * \param[in] disableUpdateSkips     If to disable update skips, useful for testing
         */
        Bias(const t_commrec               *cr,
             int                            biasIndexInCollection,
             const AwhParams               &awhParams,
             const AwhBiasParams           &awhBiasParams,
             const std::vector<DimParams>  &dimParams,
             double                         beta,
             double                         mdTimeStep,
             BiasParams::DisableUpdateSkips disableUpdateSkips = BiasParams::DisableUpdateSkips::no);

        /*! \brief Destructor */
        ~Bias();

        /*! \brief
         * Print information about initialization to log file.
         *
         * Prints information about AWH variables that are set internally
         * but might be of interest to the user.
         *
         * \param[in,out] fplog  Log file, can be nullptr.
         */
        void printInitializationToLog(FILE *fplog) const;

        /*! \brief
         * Update the coordinate value of dimension \p dim.
         *
         * Currently public because AWH calls it.
         *
         * \param[in] dim         The dimension.
         * \param[in] coordValue  The coordinate value.
         */
        inline void setCoordValue(int dim, double coordValue)
        {
            GMX_RELEASE_ASSERT(dim >= 0 && dim < ndim(), "The dimension should be in range");

            state_.setCoordValue(grid(), dim, coordValue);
        };

        /*! \brief
         * Evolves the bias at every step.
         *
         * At each step the bias step needs to:
         * - set the bias force and potential;
         * - update the free energy and bias if needed;
         * - reweight samples to extract the PMF.
         *
         * \param[out]    biasForce   The bias force.
         * \param[out] awhPotential   Bias potential.
         * \param[out] potentialJump  Change in bias potential for this bias.
         * \param[in] ms              Struct for multi-simulation communication.
         * \param[in] t               Time.
         * \param[in] step            Time step.
         * \param[in] seed            Random seed.
         * \param[in,out] fplog       Log file.
         */
        void calcForceAndUpdateBias(awh_dvec              biasForce,
                                    double               *awhPotential,
                                    double               *potentialJump,
                                    const gmx_multisim_t *ms,
                                    double                t,
                                    gmx_int64_t           step,
                                    int                   seed,
                                    FILE                 *fplog);

        /*! \brief
         * Restore the bias state from history on the master rank and broadcast it.
         *
         * \param[in] biasHistory  Bias history struct, only allowed to be nullptr on non-master ranks.
         * \param[in] cr           The communication record.
         */
        void restoreStateFromHistory(const AwhBiasHistory *biasHistory,
                                     const t_commrec      *cr);

        /*! \brief
         * Do all previously skipped updates.
         * Public for use by tests.
         */
        void doSkippedUpdatesForAllPoints();

        //! Returns the dimensionality of the bias.
        inline int ndim() const
        {
            return dimParams_.size();
        }

        /*! \brief Returns the dimension parameters.
         */
        inline const std::vector<DimParams> &dimParams() const
        {
            return dimParams_;
        }

        //! Returns the grid
        inline const Grid &grid() const
        {
            return *grid_.get();
        }

        //! Returns the bias parameters
        inline const BiasParams &params() const
        {
            return params_;
        }

        //! Returns the global state of the bias.
        inline const BiasState &state() const
        {
            return state_;
        }

        //! Returns the index of the bias.
        inline int biasIndex() const
        {
            return params_.biasIndex;
        }

    private:
        /*! \brief
         * Makes checks for the collected histograms and warns if issues are detected.
         *
         * \param[in] t            Time.
         * \param[in] step         Time step.
         * \param[in,out] fplog    Output file for warnings.
         */
        void checkHistograms(double t, gmx_int64_t step, FILE *fplog);

        /*! \brief
         * Collect samples for the force correlation analysis.
         *
         * \param[in] probWeightNeighbor  Probability weight of the neighboring points.
         * \param[in] t                   The time.
         */
        void updateForceCorrelation(const std::vector<double>    &probWeightNeighbor,
                                    double                        t);

    public:
        /*! \brief Return a const reference to the force correlation data.
         */
        const CorrelationGrid &forceCorr() const
        {
            GMX_ASSERT(forceCorr_ != nullptr, "forceCorr() should only be called with a valid force correlation object");

            return *forceCorr_.get();
        }

        /*! \brief Prepare data for writing to energy frame.
         *
         * \param[in] ms    Struct for multi-simulation communication.
         */
        void prepareOutput(const gmx_multisim_t *ms);

        /*! \brief Return the number of data blocks that have been prepared for writing.
         */
        int numEnergySubblocksToWrite() const;

        /*! \brief Write bias data blocks to energy subblocks.
         *
         * \param[in,out] subblock  Energy subblocks to write to.
         * \returns the number of subblocks written.
         */
        int writeToEnergySubblocks(t_enxsubblock *subblock) const;

        /* Data members. */
    private:
        const std::vector<DimParams> dimParams_;   /**< Parameters for each dimension. */
        std::unique_ptr<Grid>        grid_;        /**< The multidimensional grid organizing the coordinate point locations. */

        BiasParams                   params_;      /**< Constant parameters for the method. */

        BiasState                    state_;       /**< The state, both global and of the grid points */
        std::vector<int>             updateList_;  /**< List of points for update for temporary use (could be made another tempWorkSpace) */

        /* Force correlation */
        std::unique_ptr<CorrelationGrid> forceCorr_;   /**< Takes care of force correlation statistics. */

        /* I/O */
        std::unique_ptr<BiasWriter>  writer_;      /**< Takes care of AWH data output. */

        /* Temporary working vector used during the update.
         * This only here to avoid allocation at every MD step.
         */
        std::vector<double>          tempWorkSpace_;     /**< Working vector of doubles. */

        int                          numWarningsIssued_; /**< The number of warning issued in the current run. */
};

}      // namespace gmx

#endif /* GMX_AWH_BIAS_H */
