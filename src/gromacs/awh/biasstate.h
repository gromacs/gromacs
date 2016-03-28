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
 * Declares the BiasState class.
 *
 * The data members of this class are the state variables of the bias.
 * All interaction from the outside happens through the Bias class, which
 * holds important helper classes such as DimParams and Grid.
 * This class holds many methods, but more are const methods that compute
 * properties of the state.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_BIASSTATE_H
#define GMX_AWH_BIASSTATE_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

#include "coordinatestate.h"
#include "dimparams.h"
#include "histogramsize.h"

struct gmx_multisim_t;
struct t_commrec;

namespace gmx
{

struct AwhBiasHistory;
struct AwhBiasParams;
class BiasParams;
class Grid;
class GridAxis;
class PointState;

/*! \internal
 * \brief The state of a bias.
 *
 * The bias state has the current coordinate state: its value and the grid point
 * it maps to (the grid point of the umbrella potential if needed). It contains
 * a vector with the state for each point on the grid. It also
 * counts the number of updates issued and tracks which points have been sampled
 * since last update. Finally, the convergence state is a global property set
 * ultimately by the histogram size histSize in the sub-class HistogramSize,
 * since the update sizes are ~ 1/histSize.
 */
class BiasState
{
    public:
        /*! \brief Constructor.
         *
         * Constructs the global state and the point states on a provided
         * geometric grid passed in \p grid.
         *
         * \param[in] awhBiasParams    The Bias parameters from inputrec.
         * \param[in] histSizeInitial  The initial histogram size.
         * \param[in] dimParams        The dimension Parameters.
         * \param[in] grid             The grid.
         */
        BiasState(const AwhBiasParams          &awhBiasParams,
                  double                        histSizeInitial,
                  const std::vector<DimParams> &dimParams,
                  const Grid                   &grid);

        /*! \brief
         * Restore the bias state from history.
         *
         * \param[in] biasHistory  Bias history struct.
         * \param[in] grid         The grid.
         */
        void restoreFromHistory(const AwhBiasHistory &biasHistory,
                                const Grid           &grid);

        /*! \brief
         * Broadcast the bias state over the MPI ranks in this simulation.
         *
         * \param[in] cr  Struct for communication.
         */
        void broadcast(const t_commrec *cr);

        /*! \brief
         * Update the bias history with a new state.
         *
         * \param[out] biasHistory  Bias history struct.
         * \param[in]  grid         The grid.
         */
        void updateHistory(AwhBiasHistory *biasHistory,
                           const Grid     &grid) const;

        /*! \brief
         * Partition the sampling domain.
         *
         * \param[in] grid       The grid setup.
         * \param[in] pointMin   Minimum boundary point index.
         * \param[in] pointMax   Maximum boundary point index.
         * \returns the new histogram size.
         */
        double partitionDomain(const Grid &grid,
                               int         pointMin,
                               int         pointMax);

    private:
        /*! \brief
         * Convolves the PMF and sets the initial free energy to its convolution.
         *
         * \param[in] dimParams  The bias dimensions parameters
         * \param[in] grid       The grid.
         * \param[in] params     The bias parameters.
         * \param[in] ms         Struct for multi-simulation communication.
         */
        void setFreeEnergyToConvolvedPmf(const std::vector<DimParams>  &dimParams,
                                         const Grid                    &grid,
                                         const BiasParams              &params,
                                         const gmx_multisim_t          *ms);

        /*! \brief
         * Normalize the PMF histogram.
         *
         * \param[in] numSharingSims  The number of simulations sharing the bias.
         */
        void normalizePmf(int numSharingSims);

    public:
        /*! \brief
         * Initialize the state of grid coordinate points.
         *
         * \param[in] awhBiasParams   Bias parameters from inputrec.
         * \param[in] dimParams       The dimension parameters.
         * \param[in] grid            The grid.
         * \param[in] params          The bias parameters.
         * \param[in] numBias         The number of biases.
         * \param[in] ms              Struct for multi-simulation communication.
         */
        void initGridPointState(const AwhBiasParams           &awhBiasParams,
                                const std::vector<DimParams>  &dimParams,
                                const Grid                    &grid,
                                const BiasParams              &params,
                                int                            numBias,
                                const gmx_multisim_t          *ms);

        /*! \brief
         * Makes checks for the collected histograms and warns if issues are detected.
         * \param[in]     grid            The grid.
         * \param[in]     biasIndex       The index of the bias we are checking for.
         * \param[in]     t               Time.
         * \param[in,out] fplog           Output file for warnings.
         * \param[in]     maxNumWarnings  Don't issue more than this number of warnings.
         * \returns the number of warnings issued.
         */
        int checkHistograms(const Grid  &grid,
                            int          biasIndex,
                            double       t,
                            FILE        *fplog,
                            int          maxNumWarnings) const;

        /*! \brief
         * Calculates and sets the force the coordinate experiences from an umbrella centered at the given point.
         *
         * The umbrella potential is an harmonic potential given by 0.5k(coord value - point value)^2. This
         * value is also returned.
         *
         * \param[in] dimParams    The bias dimensions parameters.
         * \param[in] grid         The grid.
         * \param[in] point        Point for umbrella center.
         * \param[in,out] force    Force vector to set.
         * Returns the umbrella potential.
         */
        double calcUmbrellaForceAndPotential(const std::vector<DimParams> &dimParams,
                                             const Grid                   &grid,
                                             int                           point,
                                             awh_dvec                      force) const;

        /*! \brief
         * Calculates and sets the convolved force acting on the coordinate.
         *
         * The convolved force is the weighted sum of forces from umbrellas
         * located at each point in the grid.
         *
         * \param[in] dimParams           The bias dimensions parameters.
         * \param[in] grid                The grid.
         * \param[in] probWeightNeighbor  Probability weights of the neighbors.
         * \param[in,out] force           Bias force vector to set.
         */
        void calcConvolvedForce(const std::vector<DimParams> &dimParams,
                                const Grid                   &grid,
                                const std::vector<double>    &probWeightNeighbor,
                                awh_dvec                      force) const;

        /*! \brief
         * Move the center point of the umbrella potential.
         *
         * A new umbrella center is sampled from the biased distibution. Also, the bias
         * force is updated and the new potential is return.
         *
         * This function should only be called when the bias force is not being convolved.
         * It is assumed that the probability distribution has been updated.
         *
         * \param[in] dimParams           Bias dimension parameters.
         * \param[in] grid                The grid.
         * \param[in] probWeightNeighbor  Probability weights of the neighbors.
         * \param[in,out] biasForce       The AWH bias force.
         * \param[in] step                Step number, needed for the random number generator.
         * \param[in] seed                Random seed.
         * \param[in] indexSeed           Second random seed, should be the bias Index.
         * \returns the new potential value.
         */
        double moveUmbrella(const std::vector<DimParams> &dimParams,
                            const Grid                   &grid,
                            const std::vector<double>    &probWeightNeighbor,
                            awh_dvec                      biasForce,
                            gmx_int64_t                   step,
                            int                           seed,
                            int                           indexSeed);

    private:
        /*! \brief
         * Sets the histogram rescaling factors needed for skipped updates.
         *
         * \param[in]  params             The bias parameters.
         * \param[out] weighthistScaling  Scaling factor for the reference weight histogram.
         * \param[out] logPmfsumScaling   Log of the scaling factor for the PMF histogram.
         */
        void setSkippedUpdateHistogramScaleFactors(const BiasParams &params,
                                                   double           *weighthistScaling,
                                                   double           *logPmfsumScaling) const;

    public:
        /*! \brief
         * Do all previously skipped updates.
         * Public for use by tests.
         *
         * \param[in] params  The bias parameters.
         */
        void doSkippedUpdatesForAllPoints(const BiasParams &params);

        /*! \brief
         * Do previously skipped updates in this neighborhood.
         *
         * \param[in] params  The bias parameters.
         * \param[in] grid    The grid.
         */
        void doSkippedUpdatesInNeighborhood(const BiasParams &params,
                                            const Grid       &grid);

    private:
        /*! \brief
         * Reset the range used to make the local update list.
         *
         * \param[in] grid  The grid.
         */
        void resetLocalUpdateRange(const Grid &grid);

        /*! \brief
         * Returns the new size of the reference weight histogram in the initial stage.
         *
         * This function also takes care resetting the histogram used for covering checks
         * and for exiting the initial stage.
         *
         * \param[in] params            The bias parameters.
         * \param[in] t                 Time.
         * \param[in] detectedCovering  True if we detected that the sampling interval has been sufficiently covered.
         * \param[in,out] fplog         Log file.
         * \returns the new histogram size.
         */
        double newHistSizeInitialStage(const BiasParams &params,
                                       double            t,
                                       bool              detectedCovering,
                                       FILE             *fplog);

        /*! \brief
         * Check if the sampling region has been covered "enough" or not.
         *
         * A one-dimensional interval is defined as covered if each point has
         * accumulated the same weight as is in the peak of a discretized normal
         * distribution. For multiple dimensions, the weights are simply projected
         * onto each dimension and the multidimensional space is covered if each
         * dimension is.
         *
         * \note The covering criterion for multiple dimensions could improved, e.g.
         * by using a path finding algorithm.
         *
         * \param[in] params     The bias parameters.
         * \param[in] dimParams  Bias dimension parameters.
         * \param[in] grid       The grid.
         * \param[in] ms         Struct for multi-simulation communication.
         * \returns true if covered.
         */
        bool isCovered(const BiasParams             &params,
                       const std::vector<DimParams> &dimParams,
                       const Grid                   &grid,
                       const gmx_multisim_t         *ms) const;

        /*! \brief
         * Return the new reference weight histogram size for the current update.
         *
         * This function also takes care of checking for covering in the initial stage.
         *
         * \param[in] params     The bias parameters.
         * \param[in] t          Time.
         * \param[in] covered    True if the sampling interval has been covered enough.
         * \param[in,out] fplog  Log file.
         * \returns the new histogram size.
         */
        double newHistSize(const BiasParams &params,
                           double            t,
                           bool              covered,
                           FILE             *fplog);

    public:
        /*! \brief
         * Update the coordinate value of dimension \p dim.
         *
         * Currently public because AWH calls it.
         *
         * \param[in] grid        The grid.
         * \param[in] dim         The dimension.
         * \param[in] coordValue  The coordinate value.
         */
        inline void setCoordValue(const Grid &grid, int dim, double coordValue)
        {
            coordinateState_.setCoordValue(grid, dim, coordValue);
        };

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
         * \param[in]     dimParams   The dimension parameters.
         * \param[in]     grid        The grid.
         * \param[in]     params      The bias parameters.
         * \param[in]     ms          Struct for multi-simulation communication.
         * \param[in]     t           Time.
         * \param[in]     step        Time step.
         * \param[in,out] fplog       Log file.
         * \param[in,out] updateList  Work space to store a temporary list.
         */
        void updateFreeEnergyAndAddSamplesToHistogram(const std::vector<DimParams> &dimParams,
                                                      const Grid                   &grid,
                                                      const BiasParams             &params,
                                                      const gmx_multisim_t         *ms,
                                                      double                        t,
                                                      gmx_int64_t                   step,
                                                      FILE                         *fplog,
                                                      std::vector<int>             *updateList);

        /*! \brief
         * Update the probability weights and the convolved bias.
         *
         * Given a coordinate value, each grid point is assigned a probability
         * weight, w(point|value), that depends on the current bias function. The sum
         * of these weights is needed for normalizing the probability sum to 1 but
         * also equals the effective, or convolved, biasing weight for this coordinate
         * value. The convolved bias is needed e.g. for extracting the PMF, so we save
         * it here since this saves us from doing extra exponential function evaluations
         * later on.
         *
         * \param[in]  dimParams  The bias dimensions parameters
         * \param[in]  grid       The grid.
         * \param[out] weight     Probability weights of the neighbors.
         * \returns the convolved bias.
         */

        double updateProbabilityWeightsAndConvolvedBias(const std::vector<DimParams> &dimParams,
                                                        const Grid                   &grid,
                                                        std::vector<double>          *weight) const;

        /*! \brief
         * Save the current probability weights for future updates and analysis.
         *
         * Points in the current neighborhood will now have data meaning they
         * need to be included in the local update list of the next update.
         * Therefore, the local update range is also update here.
         *
         * \param[in] grid                The grid.
         * \param[in] probWeightNeighbor  Probability weights of the neighbors.
         */
        void sampleProbabilityWeights(const Grid                &grid,
                                      const std::vector<double> &probWeightNeighbor);

        /*! \brief
         * Sample observables for future updates or analysis.
         *
         * \param[in] grid                The grid.
         * \param[in] probWeightNeighbor  Probability weights of the neighbors.
         * \param[in] convolvedBias       The convolved bias.
         */
        void sampleCoordAndPmf(const Grid                &grid,
                               const std::vector<double> &probWeightNeighbor,
                               double                     convolvedBias);

    public:
        /*! \brief Returns the current coordinate state.
         */
        inline const CoordinateState &coordinateState() const
        {
            return coordinateState_;
        };

        /*! \brief Returns a const reference to the point state.
         */
        inline const std::vector<PointState> &points() const
        {
            return points_;
        };

        /*! \brief Returns true if we are in the initial stage.
         */
        inline bool inInitialStage() const
        {
            return histogramSize_.inInitialStage();
        };

        /* Data members */
    private:
        CoordinateState coordinateState_; /**< The Current coordinate state */

        /* The grid point state */
        std::vector<PointState> points_; /**< Vector of state of the grid points */

        /* Covering values for each point on the grid */
        std::vector<double> weightsumCovering_; /**< Accumulated weights for covering checks */

        HistogramSize       histogramSize_;     /**< GLobal histogram size related values. */

        /* Track the part of the grid sampled since the last update. */
        awh_ivec  originUpdatelist_;  /**< The origin of the rectangular region that has been sampled since last update. */
        awh_ivec  endUpdatelist_;     /**< The end of the rectangular that has been sampled since last update. */
};

/* Here follow some utility functions used by multiple files in AWH */

/*! \brief
 * Calculates the convolved bias for a given coordinate value.
 *
 * The convolved bias is the effective bias acting on the coordinate.
 * Since the bias here has arbitrary normalization, this only makes
 * sense as a relative, to other coordinate values, measure of the bias.
 *
 * \note If it turns out to be costly to calculate this pointwise
 * the convolved bias for the whole grid could be returned instead.
 *
 * \param[in] dimParams   The bias dimensions parameters
 * \param[in] grid        The grid.
 * \param[in] points      The point state.
 * \param[in] coordValue  Coordinate value.
 * \returns the convolved bias >= -GMX_DOUBLE_MAX.
 */
double calcConvolvedBias(const std::vector<DimParams>  &dimParams,
                         const Grid                    &grid,
                         const std::vector<PointState> &points,
                         const awh_dvec                &coordValue);

/*! \brief
 * Sets the given array with PMF values.
 *
 * Points outside of the biasing target region will get PMF = GMX_DOUBLE_MAX.
 * In the simplest case the PMF is simply the negative of the PMF histogram.
 * If there are sharing replicas however, histograms need to be summed
 * across multiple simulations. The output PMF is not normalized.
 *
 * \param[in] params      The bias parameters.
 * \param[in] points      The point state.
 * \param[in] ms          Struct for multi-simulation communication, needed for bias sharing replicas.
 * \param[in,out] pmf     Array returned will be of the same length as the AWH grid to store the PMF in.
 */
void calculatePmf(const BiasParams              &params,
                  const std::vector<PointState> &points,
                  const gmx_multisim_t          *ms,
                  std::vector<float>            *pmf);

/*! \brief
 * Query if something should be done at this step.
 *
 * \param[in] stepInterval   Step interval for doing something.
 * \param[in] step           Current step.
 * \returns true if something should be done at this step.
 */
static inline bool doAtStep(int         stepInterval,
                            gmx_int64_t step)
{
    GMX_ASSERT(stepInterval > 0, "All step intervals in AWH should be > 0");

    return (step % stepInterval == 0);
};

/*! \brief
 * Returns if to do checks, only returns true at free-energy update steps.
 *
 * To avoid overhead due to expensive checks, we only do checks when we
 * have taken at least as many samples as we have points.
 *
 * \param[in] params      The AWH bias parameters.
 * \param[in] pointState  The state of the points.
 * \param[in] step        Time step.
 * \returns true at steps where checks should be performed.
 */
bool isCheckStep(const BiasParams              &params,
                 const std::vector<PointState> &pointState,
                 gmx_int64_t                    step);

//! Linewidth used for warning output
static const int c_linewidth = 78;

//! Indent used for warning output
static const int c_indent    = 0;

}      // namespace gmx

#endif /* GMX_AWH_BIASSTATE_H */
