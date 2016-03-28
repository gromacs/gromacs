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
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

#include "coordstate.h"
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
 * ultimately by the histogram size histogramSize in the sub-class HistogramSize,
 * since the update sizes are ~ 1/histogramSize.
 */
class BiasState
{
    public:
        /*! \brief Constructor.
         *
         * Constructs the global state and the point states on a provided
         * geometric grid passed in \p grid.
         *
         * \param[in] awhBiasParams         The Bias parameters from inputrec.
         * \param[in] histogramSizeInitial  The estimated initial histogram size.
         *                                  This is floating-point, since histograms use weighted
         *                                  entries and grow by a floating-point scaling factor.
         * \param[in] dimParams             The dimension parameters.
         * \param[in] grid                  The bias grid.
         */
        BiasState(const AwhBiasParams          &awhBiasParams,
                  double                        histogramSizeInitial,
                  const std::vector<DimParams> &dimParams,
                  const Grid                   &grid);

        /*! \brief
         * Restore the bias state from history.
         *
         * \param[in] biasHistory  Bias history struct.
         * \param[in] grid         The bias grid.
         */
        void restoreFromHistory(const AwhBiasHistory &biasHistory,
                                const Grid           &grid);

        /*! \brief
         * Broadcast the bias state over the MPI ranks in this simulation.
         *
         * \param[in] commRecord  Struct for communication.
         */
        void broadcast(const t_commrec *commRecord);

        /*! \brief
         * Allocate and initialize a bias history with the given bias state.
         *
         * This function will be called at the start of a new simulation.
         * Note that this only sets the correct size and does produce
         * a valid history object, but with all data set to zero.
         * Actual history data is set by \ref updateHistory.
         *
         * \param[in,out] biasHistory  AWH history to initialize.
         */
        void initHistoryFromState(AwhBiasHistory *biasHistory) const;

        /*! \brief
         * Update the bias state history with the current state.
         *
         * \param[out] biasHistory  Bias history struct.
         * \param[in]  grid         The bias grid.
         */
        void updateHistory(AwhBiasHistory *biasHistory,
                           const Grid     &grid) const;

    private:
        /*! \brief Convolves the given PMF using the given AWH bias.
         *
         * \note: The PMF is in single precision, because it is a statistical
         *        quantity and therefore never reaches full float precision.
         *
         * \param[in] dimParams     The bias dimensions parameters
         * \param[in] grid          The grid.
         * \param[in,out] convolvedPmf  Array returned will be of the same length as the AWH grid to store the convolved PMF in.
         */
        void calcConvolvedPmf(const std::vector<DimParams> &dimParams,
                              const Grid                   &grid,
                              std::vector<float>           *convolvedPmf) const;

        /*! \brief
         * Convolves the PMF and sets the initial free energy to its convolution.
         *
         * \param[in] dimParams  The bias dimensions parameters
         * \param[in] grid       The bias grid.
         */
        void setFreeEnergyToConvolvedPmf(const std::vector<DimParams>  &dimParams,
                                         const Grid                    &grid);

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
         * \param[in] filename        Name of file to read PMF and target from.
         * \param[in] numBias         The number of biases.
         */
        void initGridPointState(const AwhBiasParams           &awhBiasParams,
                                const std::vector<DimParams>  &dimParams,
                                const Grid                    &grid,
                                const BiasParams              &params,
                                const std::string             &filename,
                                int                            numBias);

        /*! \brief
         * Performs statistical checks on the collected histograms and warns if issues are detected.
         *
         * \param[in]     grid            The grid.
         * \param[in]     biasIndex       The index of the bias we are checking for.
         * \param[in]     t               Time.
         * \param[in,out] fplog           Output file for warnings.
         * \param[in]     maxNumWarnings  Don't issue more than this number of warnings.
         * \returns the number of warnings issued.
         */
        int warnForHistogramAnomalies(const Grid  &grid,
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
         * \param[in]     dimParams  The bias dimensions parameters.
         * \param[in]     grid       The grid.
         * \param[in]     point      Point for umbrella center.
         * \param[in,out] force      Force vector to set.
         * Returns the umbrella potential.
         */
        double calcUmbrellaForceAndPotential(const std::vector<DimParams> &dimParams,
                                             const Grid                   &grid,
                                             int                           point,
                                             gmx::ArrayRef<double>         force) const;

        /*! \brief
         * Calculates and sets the convolved force acting on the coordinate.
         *
         * The convolved force is the weighted sum of forces from umbrellas
         * located at each point in the grid.
         *
         * \param[in]     dimParams           The bias dimensions parameters.
         * \param[in]     grid                The grid.
         * \param[in]     probWeightNeighbor  Probability weights of the neighbors.
         * \param[in]     forceWorkBuffer     Force work buffer, values only used internally.
         * \param[in,out] force               Bias force vector to set.
         */
        void calcConvolvedForce(const std::vector<DimParams> &dimParams,
                                const Grid                   &grid,
                                gmx::ArrayRef<const double>   probWeightNeighbor,
                                gmx::ArrayRef<double>         forceWorkBuffer,
                                gmx::ArrayRef<double>         force) const;

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
                            gmx::ArrayRef<const double>   probWeightNeighbor,
                            gmx::ArrayRef<double>         biasForce,
                            gmx_int64_t                   step,
                            gmx_int64_t                   seed,
                            int                           indexSeed);

    private:
        /*! \brief
         * Gets the histogram rescaling factors needed for skipped updates.
         *
         * \param[in]  params             The bias parameters.
         * \param[out] weighthistScaling  Scaling factor for the reference weight histogram.
         * \param[out] logPmfsumScaling   Log of the scaling factor for the PMF histogram.
         */
        void getSkippedUpdateHistogramScaleFactors(const BiasParams &params,
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
         * \param[in]     params            The bias parameters.
         * \param[in]     t                 Time.
         * \param[in]     detectedCovering  True if we detected that the sampling interval has been sufficiently covered.
         * \param[in,out] fplog             Log file.
         * \returns the new histogram size.
         */
        double newHistogramSizeInitialStage(const BiasParams &params,
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
         * \param[in] params        The bias parameters.
         * \param[in] dimParams     Bias dimension parameters.
         * \param[in] grid          The grid.
         * \param[in] multiSimComm  Struct for multi-simulation communication.
         * \returns true if covered.
         */
        bool isSamplingRegionCovered(const BiasParams             &params,
                                     const std::vector<DimParams> &dimParams,
                                     const Grid                   &grid,
                                     const gmx_multisim_t         *multiSimComm) const;

        /*! \brief
         * Return the new reference weight histogram size for the current update.
         *
         * This function also takes care of checking for covering in the initial stage.
         *
         * \param[in]     params   The bias parameters.
         * \param[in]     t        Time.
         * \param[in]     covered  True if the sampling interval has been covered enough.
         * \param[in,out] fplog    Log file.
         * \returns the new histogram size.
         */
        double newHistogramSize(const BiasParams &params,
                                double            t,
                                bool              covered,
                                FILE             *fplog);

    public:
        /*! \brief
         * Update the reaction coordinate value.
         *
         * \param[in] grid        The bias grid.
         * \param[in] coordValue  The current reaction coordinate value (there are no limits on allowed values).
         */
        void setCoordValue(const Grid     &grid,
                           const awh_dvec  coordValue)
        {
            coordState_.setCoordValue(grid, coordValue);
        };

        /*! \brief
         * Performs an update of the bias.
         *
         * The objective of the update is to use collected samples (probability weights)
         * to improve the free energy estimate. For sake of efficiency, the update is
         * local whenever possible, meaning that only points that have actually been sampled
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
         * \param[out] weight     Probability weights of the neighbors, SIMD aligned.
         * \returns the convolved bias.
         */

        double updateProbabilityWeightsAndConvolvedBias(const std::vector<DimParams>                  &dimParams,
                                                        const Grid                                    &grid,
                                                        std::vector < double, AlignedAllocator < double>> *weight) const;

        /*! \brief
         * Take samples of the current probability weights for future updates and analysis.
         *
         * Points in the current neighborhood will now have data meaning they
         * need to be included in the local update list of the next update.
         * Therefore, the local update range is also update here.
         *
         * \param[in] grid                The grid.
         * \param[in] probWeightNeighbor  Probability weights of the neighbors.
         */
        void sampleProbabilityWeights(const Grid                  &grid,
                                      gmx::ArrayRef<const double>  probWeightNeighbor);

        /*! \brief
         * Sample the reaction coordinate and PMF for future updates or analysis.
         *
         * These samples do not affect the (future) sampling and are thus
         * pure observables. Statisics of these are stored in the energy file.
         *
         * \param[in] grid                The grid.
         * \param[in] probWeightNeighbor  Probability weights of the neighbors.
         * \param[in] convolvedBias       The convolved bias.
         */
        void sampleCoordAndPmf(const Grid                  &grid,
                               gmx::ArrayRef<const double>  probWeightNeighbor,
                               double                       convolvedBias);
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
         * \param[in] coordValue  Coordinate value.
         * \returns the convolved bias >= -GMX_FLOAT_MAX.
         */
        double calcConvolvedBias(const std::vector<DimParams>  &dimParams,
                                 const Grid                    &grid,
                                 const awh_dvec                &coordValue) const;

        /*! \brief
         * Fills the given array with PMF values.
         *
         * Points outside of the biasing target region will get PMF = GMX_FLOAT_MAX.
         * \note: The PMF is in single precision, because it is a statistical
         *        quantity and therefore never reaches full float precision.
         *
         * \param[out] pmf  Array(ref) to be filled with the PMF values, should have the same size as the bias grid.
         */
        void getPmf(gmx::ArrayRef<float>) const;

        /*! \brief Returns the current coordinate state.
         */
        const CoordState &coordState() const
        {
            return coordState_;
        };

        /*! \brief Returns a const reference to the point state.
         */
        const std::vector<PointState> &points() const
        {
            return points_;
        };

        /*! \brief Returns true if we are in the initial stage.
         */
        bool inInitialStage() const
        {
            return histogramSize_.inInitialStage();
        };

        /*! \brief Returns the current histogram size.
         */
        inline HistogramSize histogramSize() const
        {
            return histogramSize_;
        };

        /* Data members */
    private:
        CoordState             coordState_; /**< The Current coordinate state */

        /* The grid point state */
        std::vector<PointState> points_; /**< Vector of state of the grid points */

        /* Covering values for each point on the grid */
        std::vector<double> weightSumCovering_; /**< Accumulated weights for covering checks */

        HistogramSize       histogramSize_;     /**< Global histogram size related values. */

        /* Track the part of the grid sampled since the last update. */
        awh_ivec            originUpdatelist_;  /**< The origin of the rectangular region that has been sampled since last update. */
        awh_ivec            endUpdatelist_;     /**< The end of the rectangular region that has been sampled since last update. */
};

//! Linewidth used for warning output
static const int c_linewidth = 80 - 2;

//! Indent used for warning output
static const int c_indent    = 0;

}      // namespace gmx

#endif /* GMX_AWH_BIASSTATE_H */
