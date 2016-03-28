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
 * \brief
 * Implements the BiasState class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "biasstate.h"

#include <assert.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "grid.h"
#include "internal.h"
#include "math.h"
#include "pointstate.h"


/* Sets the given array with PMF values. */
void calculatePmf(const BiasParams              &params,
                  const std::vector<PointState> &points,
                  const gmx_multisim_t          *ms,
                  std::vector<float>            *pmf)
{
    /* The PMF is just the negative of the log of the sampled PMF histogram.
     * Points with zero target weight are ignored, they will mostly contain noise.
     */

    pmf->resize(points.size());

    if (params.numSharedUpdate > 1)
    {
        /* Need to temporarily exponentiate the log weights to sum over simulations */
        for (size_t i = 0; i < points.size(); i++)
        {
            (*pmf)[i] = points[i].inTargetRegion() ? std::exp(static_cast<float>(points[i].logPmfsum())) : 0;
        }

        gmx_sumf_sim(points.size(), pmf->data(), ms);

        /* Take log again to get (non-normalized) PMF */
        for (size_t i = 0; i < points.size(); i++)
        {
            (*pmf)[i] = points[i].inTargetRegion() ? -std::log((*pmf)[i]) : GMX_FLOAT_MAX;
        }
    }
    else
    {
        for (size_t i = 0; i < points.size(); i++)
        {
            (*pmf)[i] = points[i].inTargetRegion() ? -points[i].logPmfsum() : GMX_FLOAT_MAX;
        }
    }
}

/*! \brief
 * Find the minimum free energy value.
 *
 * \param[in] pointState  The state of the points.
 * \returns the minimum free energy value.
 */
static double freeEnergyMinimumValue(const std::vector<PointState> &pointState)
{
    double fMin = GMX_DOUBLE_MAX;

    for (auto const &ps : pointState)
    {
        if (ps.inTargetRegion() && ps.freeEnergy() < fMin)
        {
            fMin = ps.freeEnergy();
        }
    }

    return fMin;
}

/*! \brief
 * Find the probability weight of a point given a coordinate value.
 *
 * The unnormalized weight is given by
 * w(point|value) = exp(bias(point) - U(value,point)),
 * where U is an harmonic umbrella potential.
 *
 * \param[in] dimParams     The bias dimensions parameters
 * \param[in] points        The point state.
 * \param[in] grid          The grid.
 * \param[in] pointIndex    Point to evaluate probability weight for.
 * \param[in] pointBias     Bias for the point (as a log weight).
 * \param[in] value         Coordinate value.
 * \returns the biased probability weight.
 */
static double biasedWeightFromPoint(const std::vector<DimParams>  &dimParams,
                                    const std::vector<PointState> &points,
                                    const Grid                    &grid,
                                    int                            pointIndex,
                                    double                         pointBias,
                                    const awh_dvec                 value)
{
    double weight = 0;

    /* Only points in the target reigon have non-zero weight */
    if (points[pointIndex].inTargetRegion())
    {
        double log_weight = pointBias;

        /* Add potential for all parameter dimensions */
        for (size_t d = 0; d < dimParams.size(); d++)
        {
            double dev = getDeviationFromPointAlongGridaxis(grid, d, pointIndex, value[d]);
            log_weight -= 0.5*dimParams[d].betak*dev*dev;
        }

        weight = std::exp(log_weight);
    }

    return weight;
}

/* Convolves the given PMF using the given AWH bias. */
void calculateConvolvedPmf(const std::vector<DimParams>  &dimParams,
                           const Grid                    &grid,
                           const BiasParams              &params,
                           const std::vector<PointState> &points,
                           const gmx_multisim_t          *ms,
                           std::vector<float>            *convolvedPmf)
{
    size_t             numPoints = grid.numPoints();

    std::vector<float> pmf;

    convolvedPmf->resize(numPoints);

    /* Get the PMF to convolve. */
    calculatePmf(params, points, ms, &pmf);

    for (size_t m = 0; m < numPoints; m++)
    {
        double           freeEnergyWeights = 0;
        const GridPoint &point             = grid.point(m);
        for (auto &neighbor : point.neighbor)
        {
            /* The negative PMF is a positive bias. */
            double biasNeighbor = -pmf[neighbor];

            /* Add the convolved PMF weights for the neighbors of this point.
               Note that this function only adds point within the target > 0 region.
               Sum weights, take the logarithm last to get the free energy. */
            freeEnergyWeights += biasedWeightFromPoint(dimParams, points, grid,
                                                       neighbor, biasNeighbor,
                                                       point.coordValue);
        }

        GMX_RELEASE_ASSERT(freeEnergyWeights > 0, "Attempting to do log(<= 0) in AWH convolved PMF calculation.");
        (*convolvedPmf)[m] = -std::log(static_cast<float>(freeEnergyWeights));
    }
}

/* Updates the target distribution for all points. */
void updateTarget(std::vector<PointState> *pointState,
                  const BiasParams        &params)
{
    double freeEnergyCutoff = 0;
    if (params.eTarget == eawhtargetCUTOFF)
    {
        freeEnergyCutoff = freeEnergyMinimumValue(*pointState) + params.targetParam;
    }

    double sumTarget = 0;
    for (auto &ps : *pointState)
    {
        sumTarget += ps.updateTarget(params, freeEnergyCutoff);
    }

    /* Normalize to 1 */
    double invSum = 1./sumTarget;
    for (auto &ps : *pointState)
    {
        ps.scaleTarget(invSum);
    }
}

/*! \brief
 * Puts together a string describing a grid point.
 *
 * \param[in] grid         The grid.
 * \param[in] point        Grid point index.
 * \param[in,out] pointstr String for the point.
 */
static void set_grid_point_value_string(const Grid &grid, int point, char *pointstr)
{
    char buf[STRLEN];

    strcpy(pointstr, "(");

    for (int d = 0; d < grid.ndim(); d++)
    {
        sprintf(buf, "%g", grid.point(point).coordValue[d]);
        strcat(pointstr, buf);
        if (d < grid.ndim() - 1)
        {
            sprintf(buf, ",");
            strcat(pointstr, buf);
        }
        else
        {
            sprintf(buf, ")");
            strcat(pointstr, buf);
        }
    }
}

/* Makes checks for the collected histograms and warns if issues are detected. */
int BiasState::checkHistograms(const Grid  &grid,
                               int          biasIndex,
                               double       t,
                               FILE        *fplog,
                               int          maxNumWarnings) const
{
    const double maxHistogramRatio = 0.5; /* Tolerance for printing a warning about the histogram ratios */

    /* Sum up the histograms and get their normalization */
    double sumVisits = 0;
    double sumW      = 0;
    for (auto &pointState : points_)
    {
        if (pointState.inTargetRegion())
        {
            sumVisits += pointState.visits_tot;
            sumW      += pointState.weightsumTot;
        }
    }
    double invNormVisits = 1./sumVisits;
    double invNormW      = 1./sumW;

    /* Check all points for warnings */
    int         numWarnings  = 0;
    size_t      numPoints    = grid.numPoints();
    for (size_t m = 0; m < numPoints; m++)
    {
        /* Skip points close to boundary or non-target region */
        const GridPoint &gridPoint = grid.point(m);
        bool             skipPoint = false;
        for (size_t n = 0; (n < gridPoint.neighbor.size()) && !skipPoint; n++)
        {
            int ineighbor = gridPoint.neighbor[n];
            skipPoint     = !points_[ineighbor].inTargetRegion();
            for (int d = 0; (d < grid.ndim()) && !skipPoint; d++)
            {
                const GridPoint &neighborPoint = grid.point(ineighbor);
                skipPoint =
                    neighborPoint.index[d] == 0 ||
                    neighborPoint.index[d] == grid.axis(d).numPoints() - 1;
            }
        }

        /* Warn if the coordinate distribution less than the target distribution with a certain fraction somewhere */
        if (!skipPoint &&
            (points_[m].weightsumTot*invNormW*maxHistogramRatio > points_[m].visits_tot*invNormVisits))
        {
            char pointvaluestr[STRLEN], warningmsg[STRLEN];

            set_grid_point_value_string(grid, m, pointvaluestr);
            sprintf(warningmsg, "\nawh%d warning: "
                    "at t = %g ps the obtained coordinate distribution at coordinate value %s "
                    "is less than a fraction %g of the reference distribution at that point. "
                    "If you are not certain about your settings you might want to increase your pull force constant or "
                    "modify your sampling region.\n",
                    biasIndex + 1, t, pointvaluestr, maxHistogramRatio);
            fprintf(fplog, "%s", wrap_lines(warningmsg, linewidth, indent, FALSE));

            numWarnings++;
        }
        if (numWarnings >= maxNumWarnings)
        {
            break;
        }
    }

    return numWarnings;
}

/* Calculates and sets the force the coordinate experiences from an umbrella centered at the given point.
 */
double BiasState::calcUmbrellaForceAndPotential(const std::vector<DimParams> &dimParams,
                                                const Grid                   &grid,
                                                int                           point,
                                                awh_dvec                      force) const
{
    double potential = 0;
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        double dev = getDeviationFromPointAlongGridaxis(grid, d, point, coordValue_[d]);

        double k   = dimParams[d].k;

        /* Force from harmonic potential 0.5*k*dev^2 */
        force[d]   = -k*dev;
        potential += 0.5*k*dev*dev;
    }

    return potential;
}

/* Calculates and sets the convolved force acting on the coordinate. */
void BiasState::calcConvolvedForce(const std::vector<DimParams> &dimParams,
                                   const Grid                   &grid,
                                   const std::vector<double>    &probWeightNeighbor,
                                   awh_dvec                      force) const
{
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        force[d] = 0;
    }

    /* Only neighboring points have non-negligible contribution. */
    const std::vector<int> &neighbor = grid.point(gridpointIndex_).neighbor;
    for (size_t n = 0; n < neighbor.size(); n++)
    {
        double weightNeighbor = probWeightNeighbor[n];
        int    indexNeighbor  = neighbor[n];

        /* Get the umbrella force from this point. The returned potential is ignored here. */
        awh_dvec forceFromNeighbor;
        calcUmbrellaForceAndPotential(dimParams, grid, indexNeighbor,
                                      forceFromNeighbor);

        /* Add the weighted umbrella force to the convolved force. */
        for (size_t d = 0; d < dimParams.size(); d++)
        {
            force[d] += forceFromNeighbor[d]*weightNeighbor;
        }
    }
}

/*! \brief
 * Sample a new reference point given the current coordinate value.
 *
 * It is assumed that the probability distribution has been updated.
 *
 * \param[in] grid                The grid.
 * \param[in] gridpointIndex      The grid point, sets the neighborhood.
 * \param[in] probWeightNeighbor  Probability weights of the neighbors.
 * \param[in] step                Step number, needed for the random number generator.
 * \param[in] seed                Random seed.
 * \param[in] indexSeed           Second random seed, should be the bias Index.
 * \returns the index of the sampled point.
 */
static int sampleReferenceGridpoint(const Grid                &grid,
                                    int                        gridpointIndex,
                                    const std::vector<double> &probWeightNeighbor,
                                    gmx_int64_t                step,
                                    int                        seed,
                                    int                        indexSeed)
{
    /* Sample new reference value from the probability distribution which is defined for the neighboring
       points of the current coordinate value.*/
    const std::vector<int> &neighbor = grid.point(gridpointIndex).neighbor;

    /* In order to use the same seed for all AWH biases and get independent
       samples we use the index of the bias. */
    int n_sampled  = get_sample_from_distribution(probWeightNeighbor, neighbor.size(),
                                                  step, seed, indexSeed);

    return neighbor[n_sampled];
}

/* Move the center point of the umbrella potential. */
double BiasState::moveUmbrella(const std::vector<DimParams> &dimParams,
                               const Grid                   &grid,
                               const std::vector<double>    &probWeightNeighbor,
                               awh_dvec                      biasForce,
                               gmx_int64_t                   step,
                               int                           seed,
                               int                           indexSeed)
{
    /* Generate and set a new coordinate reference value */
    refGridpoint_ = sampleReferenceGridpoint(grid, gridpointIndex_, probWeightNeighbor, step, seed, indexSeed);

    awh_dvec newForce;
    double   newPotential =
        calcUmbrellaForceAndPotential(dimParams, grid, refGridpoint_, newForce);

    /*  A modification of the reference value at time t will lead to a different
        force over t-dt/2 to t and over t to t+dt/2. For high switching rates
        this means the force and velocity will change signs roughly as often.
        To avoid any issues we take the average of the previous and new force
        at steps when the reference value has been moved. E.g. if the ref. value
        is set every step to (coord dvalue +/- delta) would give zero force.
     */
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        /* clang thinks newForce[d] can be garbage */
#ifndef __clang_analyzer__
        /* Average of the current and new force */
        biasForce[d] = 0.5*(biasForce[d] + newForce[d]);
#endif
    }

    return newPotential;
}

/*! \brief
 * Sets the histogram rescaling factors needed to control the histogram size.
 *
 * For sake of robustness, the reference weight histogram can grow at a rate
 * different from the actual sampling rate. Typically this happens for a limited
 * initial time, alternatively growth is scaled down by a constant factor for all
 * times. Since the size of the reference histogram sets the size of the free
 * energy update this should be reflected also in the PMF. Thus the PMF histogram
 * needs to be rescaled too.
 *
 * This function should only be called by the bias update function or wrapped by a function that
 * knows what scale factors should be applied when, e.g,
 * setSkippedUpdateHistogramScaleFactors().
 *
 * \param[in]  params             The bias parameters.
 * \param[in]  histSizeNew        New reference weight histogram size.
 * \param[in]  histSizeOld        Previous reference weight histogram size (before adding new samples).
 * \param[out] weighthistScaling  Scaling factor for the reference weight histogram.
 * \param[out] logPmfsumScaling   Log of the scaling factor for the PMF histogram.
 */
static void setHistogramUpdateScaleFactors(const BiasParams &params, double histSizeNew, double histSizeOld,
                                           double *weighthistScaling, double *logPmfsumScaling)
{

    /* The two scaling factors below are slightly different (ignoring the log factor) because the
       reference and the PMF histogram apply weight scaling differently. The weight histogram
       applies is  locally, i.e. each sample is scaled down meaning all samples get equal weight.
       It is done this way because that is what target type local Boltzmann (for which
       target = weight histogram) needs. In contrast, the PMF histogram is rescaled globally
       by repeatedly scaling down the whole histogram. The reasons for doing it this way are:
       1) empirically this is necessary for converging the PMF; 2) since the extraction of
       the PMF is theoretically only valid for a constant bias, new samples should get more
       weight than old ones for which the bias is fluctuating more. */
    *weighthistScaling = histSizeNew/(histSizeOld + params.update_weight*params.localWeightScaling);
    *logPmfsumScaling  = std::log(histSizeNew/(histSizeOld + params.update_weight));
}

/* Sets the histogram rescaling factors needed for skipped updates. */
void BiasState::setSkippedUpdateHistogramScaleFactors(const BiasParams &params,
                                                      double           *weighthistScaling,
                                                      double           *logPmfsumScaling) const
{
    GMX_ASSERT(params.skipUpdates(), "Calling function for skipped updates when skipping updates is not allowed");

    if (inInitialStage_)
    {
        /* In between global updates the reference histogram size is kept constant so we trivially know what the
            histogram size was at the time of the skipped update. */
        setHistogramUpdateScaleFactors(params, histSize_, histSize_,
                                       weighthistScaling, logPmfsumScaling);
    }
    else
    {
        /* In the final stage, the reference histogram grows at the sampling rate which gives trivial scale factors. */
        *weighthistScaling = 1;
        *logPmfsumScaling  = 0;
    }
}

/* Do all previously skipped updates. */
void BiasState::doSkippedUpdatesForAllPoints(const BiasParams &params)
{
    double weighthistScaling, logPmfsumScaling;

    setSkippedUpdateHistogramScaleFactors(params, &weighthistScaling, &logPmfsumScaling);

    for (auto &pointState : points_)
    {
        bool bUpdated = pointState.updateSkipped(params, numUpdates_, weighthistScaling, logPmfsumScaling);

        /* Update the bias for this point only if there were skipped updates in the past to avoid calculating the log unneccessarily */
        if (bUpdated)
        {
            pointState.updateBias();
        }
    }
}

/* Do previously skipped updates in this neighborhood. */
void BiasState::doSkippedUpdatesInNeighborhood(const BiasParams &params,
                                               const Grid       &grid)
{
    double weighthistScaling, logPmfsumScaling;

    setSkippedUpdateHistogramScaleFactors(params, &weighthistScaling, &logPmfsumScaling);

    /* For each neighbor point of the center point, refresh its state by adding the results of all past, skipped updates. */
    const std::vector<int> &neighbors = grid.point(gridpointIndex_).neighbor;
    for (auto &neighbor : neighbors)
    {
        bool didUpdate = points_[neighbor].updateSkipped(params, numUpdates_, weighthistScaling, logPmfsumScaling);

        if (didUpdate)
        {
            points_[neighbor].updateBias();
        }
    }
}

/*! \brief
 * Merge update lists from multiple sharing simulations.
 *
 * \param[in,out] updateList  Update list for this simulation (assumed >= npoints long).
 * \param[in] numPoints       Total number of points.
 * \param[in] ms              Struct for multi-simulation communication.
 */
static void mergeSharedUpdateLists(std::vector<int> *updateList, int numPoints, const gmx_multisim_t *ms)
{
    std::vector<int> numUpdatesOfPoint;

    /* Flag the update points of this sim.
       TODO: we can probably avoid allocating this array and just use the input array. */
    numUpdatesOfPoint.resize(numPoints, 0);
    for (auto &pointIndex : *updateList)
    {
        numUpdatesOfPoint[pointIndex] = 1;
    }

    /* Sum over the sims to get all the flagged points */
    gmx_sumi_sim(numPoints, numUpdatesOfPoint.data(), ms);

    /* Collect the indices of the flagged points in place. The resulting array will be the merged update list.*/
    updateList->resize(0);
    for (int m = 0; m < numPoints; m++)
    {
        if (numUpdatesOfPoint[m] > 0)
        {
            updateList->push_back(m);
        }
    }
}

/*! \brief
 * Generate an update list of points sampled since the last update.
 *
 * \param[in] grid              The AWH bias.
 * \param[in] points            The point state.
 * \param[in] originUpdatelist  The origin of the rectangular region that has been sampled since last update.
 * \param[in] endUpdatelist     The end of the rectangular that has been sampled since last update.
 * \param[in,out] updateList    Local update list to set (assumed >= npoints long).
 */
static void makeLocalUpdateList(const Grid                    &grid,
                                const std::vector<PointState> &points,
                                const awh_ivec                 originUpdatelist,
                                const awh_ivec                 endUpdatelist,
                                std::vector<int>              *updateList)
{
    awh_ivec origin_dim, npoints_dim;

    /* Define the update search grid */
    for (int d = 0; d < grid.ndim(); d++)
    {
        origin_dim[d]  = originUpdatelist[d];
        npoints_dim[d] = endUpdatelist[d] - originUpdatelist[d] + 1;

        /* Because the end_updatelist is unwrapped it can be > (npoints - 1) so that npoints_dim can be > npoints in grid.
           This helps for calculating the distance/number of points but should be removed and fixed when the way of
           updating origin/end updatelist is changed (see sampleProbabilityWeights). */
        npoints_dim[d] = std::min(grid.axis(d).numPoints(), npoints_dim[d]);
    }

    /* Make the update list */
    updateList->clear();
    int  point_index = -1;
    bool pointExists = true;
    while (pointExists)
    {
        pointExists = getNextPointInSubgrid(grid, origin_dim, npoints_dim, &point_index);

        if (pointExists && points[point_index].inTargetRegion())
        {
            updateList->push_back(point_index);
        }
    }
}

/* Reset the range used to make the local update list. */
void BiasState::resetLocalUpdateRange(const Grid &grid)
{
    for (int d = 0; d < grid.ndim(); d++)
    {
        /* This gives the  minimum range consisting only of the current closest point. */
        originUpdatelist_[d] = grid.point(gridpointIndex_).index[d];
        endUpdatelist_[d]    = grid.point(gridpointIndex_).index[d];
    }
}

/*! \brief
 * Add partial histograms (accumulating between updates) to accumulating histograms.
 *
 * \param[in] pointState         The state of the points in the bias.
 * \param[in] numSharedUpdate    The number of biases sharing the histrogram.
 * \param[in] ms                 Struct for multi-simulation communication.
 * \param[in] localUpdateList    List of points with data.
 */
static void sumHistograms(std::vector<PointState> *pointState,
                          int numSharedUpdate, const gmx_multisim_t *ms,
                          const std::vector<int> &localUpdateList)
{
    /* The covering checking histograms are added before summing over simulations, so that the weights from different
       simulations are kept distinguishable. */
    for (auto &iglobal : localUpdateList)
    {
        (*pointState)[iglobal].weightsum_covering += (*pointState)[iglobal].weightsum_iteration;
    }

    /* Sum histograms over multiple simulations if needed. */
    if (numSharedUpdate > 1)
    {
        /* Collect the weights and counts in linear arrays to be able to use gmx_sumd_sim. */
        std::vector<double> weightDistr, coordVisits;

        weightDistr.resize(localUpdateList.size());
        coordVisits.resize(localUpdateList.size());

        for (size_t ilocal = 0; ilocal < localUpdateList.size(); ilocal++)
        {
            int iglobal         = localUpdateList[ilocal];

            weightDistr[ilocal] = (*pointState)[iglobal].weightsum_iteration;
            coordVisits[ilocal] = (*pointState)[iglobal].visits_iteration;
        }

        gmx_sumd_sim(weightDistr.size(), weightDistr.data(), ms);
        gmx_sumd_sim(coordVisits.size(), coordVisits.data(), ms);

        /* Transfer back the result */
        for (size_t ilocal = 0; ilocal < localUpdateList.size(); ilocal++)
        {
            int iglobal                                = localUpdateList[ilocal];

            (*pointState)[iglobal].weightsum_iteration = weightDistr[ilocal];
            (*pointState)[iglobal].visits_iteration    = coordVisits[ilocal];
        }
    }

    /* Now add the partial counts and weights to the accumulating histograms.
       Note: we still need to use the weights for the update so we wait
       with resetting them until the end of the update. */
    for (auto &iglobal : localUpdateList)
    {
        double weights = (*pointState)[iglobal].weightsum_iteration;
        double visits  = (*pointState)[iglobal].visits_iteration;
        (*pointState)[iglobal].weightsumTot += weights;
        (*pointState)[iglobal].visits_tot   += visits;
    }
}

/* Returns the new size of the reference weight histogram in the initial stage. */
double BiasState::newHistSizeInitialStage(const BiasParams &params,
                                          double            t,
                                          bool              detectedCovering,
                                          FILE             *fplog)
{
    /* The histogram size is kept constant until the sampling region has been covered
       and the the current sample weight is large enough and the histogram is ready. */
    if (!detectedCovering ||
        (scaledSampleWeight_ < maxScaledSampleWeight_) ||
        equilibrateHistogram_)
    {
        return histSize_;
    }

    /* Reset the covering weight histogram. If we got this far we are either entering a
       new covering stage with a new covering histogram or exiting the initial stage
       altogether. */
    for (auto &pointState : points_)
    {
        pointState.weightsum_covering = 0.;
    }

    /*  The current sample weigth is now the maximum. */
    double prevMaxScaledSampleWeight = maxScaledSampleWeight_;
    maxScaledSampleWeight_ = scaledSampleWeight_;

    /* Increase the histogram size by a constant scale factor if we can, i.e. if the sample weight
       resulting from such a scaling is still larger than the previous maximum sample weight
       (ensuring that the sample weights at the end of each covering stage are monotonically
       increasing). If we cannot, exit the initial stage without changing the histogram size. */

    /* The scale factor. The value is not very critical but should obviously be > 1 (or the exit
       will happen very late) and probably < 5 or so (or there will be no initial stage). */
    static const double growthFactor = 3;

    /* The scale factor is in most cases very close to the histogram growth factor. */
    double              scaleFactor = growthFactor/(1. + params.update_weight*params.localWeightScaling/histSize_);

    bool                bExit       = (scaledSampleWeight_ - std::log(scaleFactor) <= prevMaxScaledSampleWeight);
    double              histSizeNew = bExit ? histSize_ : histSize_*growthFactor;

    /* Update the AWH bias about the exit. */
    inInitialStage_ = !bExit;

    /* Print information about coverings and if there was an exit. */
    if (fplog != nullptr)
    {
        char                buf[STRLEN];
        sprintf(buf, "\nawh%d:", params.biasIndex + 1);
        fprintf(fplog, "%s covering at t = %g ps. Decreased the update size.\n", buf, t);

        if (bExit)
        {
            fprintf(fplog, "%s out of the initial stage at t = %g.\n", buf, t);
            /* It would be nice to have a way of estimating a minimum time until exit but it
               is difficult because the exit time is determined by how long it takes to cover
               relative to the time it takes to "regaining" enough sample weight. The latter
               is easy to calculate, but how the former depends on the histogram size
               is not known. */
        }
        fflush(fplog);
    }
    return histSizeNew;
}

/*! \brief
 * Label points along an axis as covered or not.
 *
 * A point is covered if it is surrounded by visited points up to a radius = coverRadius.
 *
 * \param[in] visited           Visited? For each point.
 * \param[in] checkCovering     Check for covering? For each point.
 * \param[in] npoints           Number of points.
 * \param[in] period            Period in number of points.
 * \param[in] coverRadius       Cover radius, in points, needed for defining a point as covered.
 * \param[in,out] covered       Covered? For each point.
 */
static void labelCoveredPoints(const std::vector<bool> &visited, const std::vector<bool> &checkCovering, int npoints, int period,
                               int coverRadius, std::vector<int> *covered)
{

    bool haveFirstNotVisited = false;
    int  firstNotVisited     = -1, notVisitedLow = -1, notVisitedHigh = -1;

    for (int n = 0; n < npoints; n++)
    {
        if (checkCovering[n] && !visited[n])
        {
            if (!haveFirstNotVisited)
            {
                notVisitedLow       = n;
                firstNotVisited     = n;
                haveFirstNotVisited = true;
            }
            else
            {
                notVisitedHigh = n;

                /* Have now an interval I = [notVisitedLow,notVisitedHigh] of visited points bounded
                   by unvisited points. The unvisted end points affect the coveredness of the visited
                   with a reach equal to the cover radius. */
                int notCoveredLow  = notVisitedLow + coverRadius;
                int notCoveredHigh = notVisitedHigh - coverRadius;
                for (int i = notVisitedLow; i <= notVisitedHigh; i++)
                {
                    (*covered)[i] = (i > notCoveredLow) && (i < notCoveredHigh);
                }

                /* Find a new interval to set covering for. Make the notVisitedHigh of this interval the
                   notVisitedLow of the next. */
                notVisitedLow = notVisitedHigh;
            }
        }
    }

    /* Have labelled all the internal points. Now take care of the boundary regions. */
    if (!haveFirstNotVisited)
    {
        /* No non-visited points <=> all points visited => all points covered. */

        for (int n = 0; n < npoints; n++)
        {
            (*covered)[n] = TRUE;
        }
    }
    else
    {
        int lastNotVisited = notVisitedLow;

        /* For periodic boundaries, non-visited points can influence points
           on the other side of the boundary so we need to wrap around. */

        /* Lower end. For periodic boundaries the last upper end not visited point becomes the low-end not visited point.
           For non-periodic boundaries there is no lower end point so a dummy value is used. */
        int notVisitedHigh = firstNotVisited;
        int notVisitedLow  = period > 0 ? (lastNotVisited - period) : -(coverRadius + 1);

        int notCoveredLow  =  notVisitedLow + coverRadius;
        int notCoveredHigh = notVisitedHigh - coverRadius;

        for (int i = 0; i <= notVisitedHigh; i++)
        {
            /* For non-periodic boundaries notCoveredLow = -1 will impose no restriction. */
            (*covered)[i] = (i > notCoveredLow) && (i < notCoveredHigh);
        }

        /* Upper end. Same as for lower end but in the other direction. */
        notVisitedHigh = period > 0 ? (firstNotVisited + period) : (npoints + coverRadius);
        notVisitedLow  = lastNotVisited;

        notCoveredLow  = notVisitedLow + coverRadius;
        notCoveredHigh = notVisitedHigh - coverRadius;

        for (int i = notVisitedLow; i <= npoints - 1; i++)
        {
            /* For non-periodic boundaries notCoveredHigh = npoints will impose no restriction. */
            (*covered)[i] = (i > notCoveredLow) && (i < notCoveredHigh);
        }
    }
}

/* Check if the sampling region has been covered "enough" or not. */
bool BiasState::isCovered(const BiasParams             &params,
                          const std::vector<DimParams> &dimParams,
                          const Grid                   &grid,
                          const gmx_multisim_t         *ms) const
{
    /* Allocate and initialize arrays: one for checking visits along each dimension,
       one for keeping track of which points to check and one for the covered points.
       Possibly these could be kept as AWH variables to avoid these allocations. */
    struct CheckDim
    {
        std::vector<bool> visited;
        std::vector<bool> checkCovering;
        // We use int for the covering array since we might use gmx_sumi_sim.
        std::vector<int>  covered;
    };

    std::vector<CheckDim> checkDim;
    checkDim.resize(grid.ndim());

    for (int d = 0; d < grid.ndim(); d++)
    {
        size_t npoints = grid.axis(d).numPoints();
        checkDim[d].visited.resize(npoints, false);
        checkDim[d].checkCovering.resize(npoints, false);
        checkDim[d].covered.resize(npoints, 0);
    }

    /* Set visited points along each dimension and which points should be checked for covering.
       Specifically, points above the free energy cutoff (if there is one) or points outside
       of the target region are ignored. */

    /* Set the free energy cutoff */
    double maxFreeEnergy = GMX_DOUBLE_MAX;

    if (params.eTarget == eawhtargetCUTOFF)
    {
        maxFreeEnergy = freeEnergyMinimumValue(points_) + params.targetParam;
    }

    /* Set the cutoff weight for a point to be considered visited. */
    double weight_peak = 1;
    for (int d = 0; d < grid.ndim(); d++)
    {
        weight_peak *= grid.axis(d).spacing()*std::sqrt(dimParams[d].betak*0.5*M_1_PI);
    }

    /* Project the sampling weights onto each dimension */
    for (size_t m = 0; m < grid.numPoints(); m++)
    {
        const PointState &pointState = points_[m];

        for (int d = 0; d < grid.ndim(); d++)
        {
            int n = grid.point(m).index[d];

            /* Is visited if it was already visited or if there is enough weight at the current point */
            checkDim[d].visited[n] = checkDim[d].visited[n] || (pointState.weightsum_covering > weight_peak);

            /* Check for covering if there is at least point in this slice that is in the target region and within the cutoff */
            checkDim[d].checkCovering[n] = checkDim[d].checkCovering[n] || (pointState.inTargetRegion() && pointState.freeEnergy() < maxFreeEnergy);
        }
    }

    /* Label each point along each dimension as covered or not. */
    for (int d = 0; d < grid.ndim(); d++)
    {
        labelCoveredPoints(checkDim[d].visited, checkDim[d].checkCovering, grid.axis(d).numPoints(), grid.axis(d).numPointsInPeriod(),
                           params.coverRadius[d], &checkDim[d].covered);
    }

    /* Now check for global covering. Each dimension needs to be covered separately.
       A dimension is covered if each point is covered.  Multiple simulations collectively
       cover the points, i.e. a point is covered if any of the simulations covered it.
       However, visited points are not shared, i.e. if a point is covered or not is
       determined by the visits of a single simulation. In general the covering criterion is
       all points covered => all points are surrounded by visited points up to a radius = coverRadius.
       For 1 simulation, all points covered <=> all points visited. For multiple simulations
       however, all points visited collectively !=> all points covered, except for coverRadius = 0.
       In the limit of large coverRadius, all points covered => all points visited by at least one
       simulation (since no point will be covered until all points have been visited by a
       single simulation). Basically coverRadius sets how much "connectedness" (or mixing) a point
       needs with surrounding points before sharing covering information with other simulations. */

    /* Communicate the covered points between sharing simulations if needed. */
    if (params.numSharedUpdate > 1)
    {
        /* For multiple dimensions this may not be the best way to do it. */
        for (int d = 0; d < grid.ndim(); d++)
        {
            gmx_sumi_sim(grid.axis(d).numPoints(), checkDim[d].covered.data(), ms);
        }
    }

    /* Now check if for each dimension all points are covered. Break if not true. */
    bool allPointsCovered = true;
    for (int d = 0; d < grid.ndim() && allPointsCovered; d++)
    {
        for (int n = 0; n < grid.axis(d).numPoints() && allPointsCovered; n++)
        {
            allPointsCovered = checkDim[d].covered[n];
        }
    }

    return allPointsCovered;
}

/*! \brief
 * Checks if the histogram has equilibrated to the target distribution.
 *
 * The histogram is considered equilibrated if, for a minimum fraction of
 * the target region, the relative error of the sampled weight relative
 * to the target is less than a tolerance value.
 *
 * \param[in] pointStates  The state of the bias points.
 * \returns true if the histogram is equilibrated.
 */
static bool histogramIsEquilibrated(const std::vector<PointState> &pointStates)
{
    /* Get the total weight of the total weight histogram; needed for normalization. */
    double totalWeight     = 0;
    int    numTargetPoints = 0;
    for (auto &pointState : pointStates)
    {
        if (!pointState.inTargetRegion())
        {
            continue;
        }
        totalWeight += pointState.weightsumTot;
        numTargetPoints++;
    }
    GMX_ASSERT(totalWeight > 0, "No samples when normalizing AWH histogram.");
    double              inverseTotalWeight   = 1./totalWeight;

    /* Points with target weight below a certain cutoff are ignored. */
    static const double minTargetCutoff  = 0.05;
    double              minTargetWeight  = 1./numTargetPoints*minTargetCutoff;

    /* Points with error less than this tolerance pass the check.*/
    static const double errorTolerance   = 0.2;

    /* Sum up weight of points that do or don't pass the check. */
    double equilibratedWeight    = 0;
    double notEquilibratedWeight = 0;
    for (auto &pointState : pointStates)
    {
        double targetWeight  = pointState.target();
        double sampledWeight = pointState.weightsumTot*inverseTotalWeight;

        /* Ignore these points. */
        if (!pointState.inTargetRegion() || targetWeight < minTargetWeight)
        {
            continue;
        }

        if (std::fabs(sampledWeight/targetWeight - 1) > errorTolerance)
        {
            notEquilibratedWeight += targetWeight;
        }
        else
        {
            equilibratedWeight += targetWeight;
        }
    }

    /* It is enough if sampling in at least a fraction of the target region follows the target
       distribution. Boundaries will in general fail and this should be ignored (to some extent). */
    static const double minFraction = 0.8;

    return equilibratedWeight/(equilibratedWeight + notEquilibratedWeight) > minFraction;;
}

/* Return the new reference weight histogram size for the current update. */
double BiasState::newHistSize(const BiasParams &params,
                              double            t,
                              bool              covered,
                              FILE             *fplog)
{
    double histSizeNew;
    if (inInitialStage_)
    {
        /* Only bother with checking equilibration if we have covered already. */
        if (equilibrateHistogram_ && covered)
        {
            /* The histogram is equilibrated at most once. */
            equilibrateHistogram_ = !histogramIsEquilibrated(points_);

            static bool printAboutCovering = true;
            if (!equilibrateHistogram_)
            {
                char buf[STRLEN];
                sprintf(buf, "\nawh%d:", params.biasIndex + 1);
                fprintf(fplog, "%s equilibrated histogram at t = %g ps.\n", buf, t);
            }
            else if (printAboutCovering)
            {
                char buf[STRLEN];
                sprintf(buf, "\nawh%d:", params.biasIndex + 1);
                fprintf(fplog, "%s covered but histogram not equilibrated at t = %g ps.\n", buf, t);
                printAboutCovering = false; /* Just print once. */
            }
        }

        /* In the initial stage, the histogram grows dynamically as a function of the number of coverings. */
        histSizeNew = newHistSizeInitialStage(params, t, covered, fplog);
    }
    else
    {
        /* If not in the initial stage, the histogram grows at a linear, possibly scaled down, rate. */
        histSizeNew = histSize_ + params.update_weight*params.localWeightScaling;
    }

    return histSizeNew;
}

/*! \brief
 * Normalizes the free energy and PMF sum.
 *
 * \param[in] pointState  The state of the points.
 */
static void normalizeFreeEnergyAndPmfSum(std::vector<PointState> *pointState)
{
    double minF = freeEnergyMinimumValue(*pointState);

    for (auto &ps : *pointState)
    {
        ps.normalizeFreeEnergyAndPmfSum(minF);
    }
}

/* Performs an update of the bias. */
void BiasState::updateFreeEnergyAndAddSamplesToHistogram(const std::vector<DimParams> &dimParams,
                                                         const Grid                   &grid,
                                                         const BiasParams             &params,
                                                         const gmx_multisim_t         *ms,
                                                         double                        t,
                                                         gmx_int64_t                   step,
                                                         FILE                         *fplog,
                                                         std::vector<int>             *updateList)
{
    /* Note hat updateList is only used in this scope and is always
     * re-initialized. We do not use a local vector, because that would
     * cause reallocation every time this funtion is called and the vector
     * can be the size of the whole grid.
     */

    /* Make a list of all local points, i.e. those that could have been touched since
       the last update. These are the points needed for summing histograms below
       (non-local points only add zeros). For local updates, this will also be the
       final update list. */
    makeLocalUpdateList(grid, points_, originUpdatelist_, endUpdatelist_,
                        updateList);
    if (params.numSharedUpdate > 1)
    {
        mergeSharedUpdateLists(updateList, points_.size(), ms);
    }

    /* Reset the range for the next update */
    resetLocalUpdateRange(grid);

    /* Add samples to histograms for all local points and sync simulations if needed */
    sumHistograms(&points_, params.numSharedUpdate, ms, *updateList);

    /* Renormalize the free energy if values are too large. */
    bool needToNormalizeFreeEnergy = false;
    for (auto &iglobal : *updateList)
    {
        /* We want to keep the absolute value of the free energies to be less c_largePositiveExponent
           to be able to safely pass these values to exp(). The check below ensures this as long as
           the free energy values grow less than 0.5*c_largePositiveExponent in a return time to this
           neighborhood. For reasonable update sizes it's unlikely that this requirement would be
           broken. */
        if (std::fabs(points_[iglobal].freeEnergy()) > 0.5*c_largePositiveExponent)
        {
            needToNormalizeFreeEnergy = true;
            break;
        }
    }

    /* Update target distribution? */
    bool doUpdateTarget = (params.eTarget != eawhtargetCONSTANT &&
                           doAtStep(params.nstupdate_target, step));

    /* In the initial stage, the histogram grows dynamically as a function of the number of coverings. */
    bool detectedCovering = false;
    if (inInitialStage_)
    {
        detectedCovering = (isCheckStep(params, points_, step) &&
                            isCovered(params, dimParams, grid, ms));
    }

    /* The weighthistogram size after this update. */
    double histSizeNew = newHistSize(params, t, detectedCovering, fplog);

    /* Make the update list. Usually we try to only update local points but if the update has
       non-trivial or non-deterministic affect on non-local points a global update is needed.
       This is the case when: 1) a covering occurred in the initial stage, leading to non-trivial
       histogram rescaling factors; or 2) the target distribution will be updated since we don't
       make any assumption on its form; or 3) the AWH parameters are such that we never attempt
       to skip non-local updates; or 4) the free energy values have grown so large that the
       a renormalization is needed. */
    if (doUpdateTarget || detectedCovering || !params.skipUpdates() || needToNormalizeFreeEnergy)
    {
        /* Global update, just add all points. */
        updateList->clear();
        for (size_t m = 0; m < points_.size(); m++)
        {
            if (points_[m].inTargetRegion())
            {
                updateList->push_back(m);
            }
        }
    }

    /* Set histogram scale factors. */
    double   weighthistScalingSkipped = 0, logPmfsumScalingSkipped = 0;
    double   weighthistScalingNew, logPmfsumScalingNew;
    if (params.skipUpdates())
    {
        setSkippedUpdateHistogramScaleFactors(params, &weighthistScalingSkipped, &logPmfsumScalingSkipped);
    }
    setHistogramUpdateScaleFactors(params, histSizeNew, histSize_,
                                   &weighthistScalingNew, &logPmfsumScalingNew);

    /* Update free energy and reference weight histogram for points in the update list. */
    for (auto &pointIndex : *updateList)
    {
        PointState *pointStateToUpdate = &points_[pointIndex];

        /* Do updates from previous update steps that were skipped because this point was at that time non-local. */
        if (params.skipUpdates())
        {
            pointStateToUpdate->updateSkipped(params, numUpdates_, weighthistScalingSkipped, logPmfsumScalingSkipped);
        }

        /* Now do an update with new sampling data. */
        pointStateToUpdate->updateNew(params, numUpdates_, weighthistScalingNew, logPmfsumScalingNew);

        /* Reset histograms collected for this update. */
        pointStateToUpdate->weightsum_iteration = 0;
        pointStateToUpdate->visits_iteration    = 0;
    }

    /* Only update the histogram size after we are done with the local point updates */
    histSize_ = histSizeNew;

    /* The weight of new samples relative to previous ones change when the histogram is
       rescaled. We keep the log since this number can become very large. */
    scaledSampleWeight_ -= std::log(weighthistScalingNew);

    if (needToNormalizeFreeEnergy)
    {
        normalizeFreeEnergyAndPmfSum(&points_);
    }

    if (doUpdateTarget)
    {
        /* The target distribution is always updated for all points at once. */
        updateTarget(&points_, params);
    }

    /* Update the bias. The bias is updated separately and last since it simply a function of
       the free energy and the target distribution and we want to avoid doing extra work. */
    for (auto &pointIndex : *updateList)
    {
        points_[pointIndex].updateBias();
    }

    /* Increase the update counter. */
    numUpdates_ += 1;
}

/* Update the probability weights and the convolved bias. */
double BiasState::updateProbabilityWeightsAndConvolvedBias(const std::vector<DimParams> &dimParams,
                                                           const Grid                   &grid,
                                                           std::vector<double>          *weight) const
{
    /* Only neighbors of the current coordinate value will have a non-negligible chance of getting sampled */
    const std::vector<int> &neighbors = grid.point(gridpointIndex_).neighbor;

    weight->resize(neighbors.size());

    double weightSum = 0;
    for (size_t n = 0; n < neighbors.size(); n++)
    {
        int neighbor = neighbors[n];
        (*weight)[n] = biasedWeightFromPoint(dimParams, points_, grid,
                                             neighbor, points_[neighbor].bias(),
                                             coordValue_);
        weightSum   += (*weight)[n];
    }
    GMX_RELEASE_ASSERT(weightSum > 0, "zero probability weight when updating AWH probability weights.");

    /* Normalize probabilities to sum to 1 */
    double invWeightSum = 1/weightSum;
    for (auto &w : *weight)
    {
        w *= invWeightSum;
    }

    /* Return the convolved bias */
    return std::log(weightSum);
}

/* Calculates the convolved bias for the given coordinate value. */
double calcConvolvedBias(const std::vector<DimParams>  &dimParams,
                         const Grid                    &grid,
                         const std::vector<PointState> &points,
                         const awh_dvec                &coordValue)
{
    int              point      = grid.nearestIndex(coordValue);
    const GridPoint &gridPoint  = grid.point(point);

    /* Sum the probability weights from the neighborhood of the given point */
    double weight_sum = 0;
    for (auto &neighbor : gridPoint.neighbor)
    {
        weight_sum += biasedWeightFromPoint(dimParams, points, grid,
                                            neighbor, points[neighbor].bias(),
                                            coordValue);
    }

    /* Returns -GMX_DOUBLE_MAX if no neighboring points where in the target region. */
    return (weight_sum > 0) ? std::log(weight_sum) : -GMX_DOUBLE_MAX;
}

/* Save the current probability weights for future updates and analysis. */
void BiasState::sampleProbabilityWeights(const Grid                &grid,
                                         const std::vector<double> &probWeightNeighbor)
{
    const std::vector<int> &neighbor = grid.point(gridpointIndex_).neighbor;

    /* Save weights for next update */
    for (size_t n = 0; n < neighbor.size(); n++)
    {
        points_[neighbor[n]].weightsum_iteration += probWeightNeighbor[n];
    }

    /* Update the local update range. Two corner points define this rectangular domain. We need to
       choose two new corner points such that the new domain contains both the old update range and
       the current neighborhood. In the simplest case when an update is performed every sample, the
        update range would simply equal the current neighborhood. */
    int m_neighbor_origin = neighbor[0];
    int m_neighbor_end    = neighbor[neighbor.size() - 1];
    for (int d = 0; d < grid.ndim(); d++)
    {
        int origin_d, end_d;

        origin_d = grid.point(m_neighbor_origin).index[d];
        end_d    = grid.point(m_neighbor_end).index[d];

        if (origin_d > end_d)
        {
            /* Unwrap if wrapped around the boundary (only happens for periodic boundaries).
               This has been already for the stored index interval. */

            /* TODO: what we want to do is to find the smallest the update interval that contains all points that need to be updated.
               This amounts to combining two intervals, the current [origin, end] update interval and the new touched neighborhood
               into a new interval that contains all points from both the old intervals.

               For periodic boundaries it becomes slightly more complicated than for closed boundaries because the it needs not be
               true that origin < end (so one can't simply relate the origin/end in the min()/max() below). The strategy here is to choose the
               origin closest to a reference point (index 0) and then unwrap the end index if needed and choose the largest end index.
               This ensures that both intervals are in the new interval but it's not necessarily the smallest.
               I can't think of a better way of solving this than going through each possibility and checking them.
             */
            end_d += grid.axis(d).numPointsInPeriod();
        }

        originUpdatelist_[d] = std::min(originUpdatelist_[d], origin_d);
        endUpdatelist_[d]    = std::max(endUpdatelist_[d], end_d);
    }
}

/* Sample observables for future updates or analysis. */
void BiasState::sampleCoordAndPmf(const Grid                &grid,
                                  const std::vector<double> &probWeightNeighbor,
                                  double                     convolvedBias)
{
    /* Sampling-based deconvolution extracting the PMF.
     * Update the PMF histogram with the current coordinate value.
     *
     * Because of the finite width of the harmonic potential, the free energy
     * defined for each coordinate point does not exactly equal that of the
     * actual coordinate, the PMF. However, the PMF can be estimated by applying
     * the relation exp(-PMF) = exp(-bias_convolved)*P_biased/Z, i.e. by keeping a
     * reweighted histogram of the coordinate value. Strictly, this relies on
     * the unknown normalization constant Z being either constant or known. Here,
     * neither is true except in the long simulation time limit. Empirically however,
     * it works (mainly because how the PMF histogram is rescaled).
     */

    /* Only save coordinate data that is in range (the given index is always
     * in range even if the coordinate value is not).
     */
    if (grid.covers(coordValue_))
    {
        /* Save PMF sum and keep a histogram of the sampled coordinate values */
        points_[gridpointIndex_].samplePmf(convolvedBias);
    }

    /* Save probability weights for the update */
    sampleProbabilityWeights(grid, probWeightNeighbor);
}

/* Update the coordinate value with coordValue. */
void BiasState::setCoordValue(const Grid &grid, int dim, double coordValue)
{
    coordValue_[dim] = coordValue;

    /* The grid point closest to the coordinate value defines the current
     * neighborhood of points. Besides at steps when global updates and/or
     * checks are performed, only the neighborhood will be touched.
     */
    gridpointIndex_ = grid.nearestIndex(coordValue_);
}

/* Update the bias history with a new state. */
void BiasState::updateHistory(AwhBiasHistory *biasHistory,
                              const Grid     &grid) const
{
    GMX_RELEASE_ASSERT(biasHistory->pointState.size() == points_.size(), "The AWH history setup does not match the AWH state.");

    AwhBiasStateHistory *stateHistory = &biasHistory->state;
    stateHistory->refGridpoint        = refGridpoint_;

    for (size_t m = 0; m < biasHistory->pointState.size(); m++)
    {
        const PointState     &ps  = points_[m];
        AwhPointStateHistory *psh = &biasHistory->pointState[m];

        psh->target               = ps.target();
        psh->free_energy          = ps.freeEnergy();
        psh->bias                 = ps.bias();
        psh->weightsum_iteration  = ps.weightsum_iteration;
        psh->weightsum_covering   = ps.weightsum_covering;
        psh->weightsum_tot        = ps.weightsumTot;
        psh->weightsum_ref        = ps.weightsumRef();
        psh->last_update_index    = ps.lastUpdateIndex();
        psh->log_pmfsum           = ps.logPmfsum();
        psh->visits_iteration     = ps.visits_iteration;
        psh->visits_tot           = ps.visits_tot;
    }

    stateHistory->in_initial              = inInitialStage_;
    stateHistory->equilibrateHistogram    = equilibrateHistogram_;
    stateHistory->histSize                = histSize_;

    stateHistory->origin_index_updatelist = multidim_gridindex_to_linear(grid,
                                                                         originUpdatelist_);
    stateHistory->end_index_updatelist    = multidim_gridindex_to_linear(grid,
                                                                         endUpdatelist_);

    stateHistory->scaledSampleWeight      = scaledSampleWeight_;
    stateHistory->maxScaledSampleWeight   = maxScaledSampleWeight_;
    stateHistory->numUpdates              = numUpdates_;
}

/* Restore the bias state from history. */
void BiasState::restoreFromHistory(const AwhBiasHistory &biasHistory,
                                   const Grid           &grid)
{
    const AwhBiasStateHistory &stateHistory = biasHistory.state;

    refGridpoint_ = stateHistory.refGridpoint;

    for (size_t m = 0; m < points_.size(); m++)
    {
        points_[m].setFromHistory(biasHistory.pointState[m]);
    }

    inInitialStage_         = stateHistory.in_initial;
    equilibrateHistogram_   = stateHistory.equilibrateHistogram;
    histSize_               = stateHistory.histSize;

    linear_gridindex_to_multidim(grid, stateHistory.origin_index_updatelist, originUpdatelist_);
    linear_gridindex_to_multidim(grid, stateHistory.end_index_updatelist, endUpdatelist_);

    scaledSampleWeight_     = stateHistory.scaledSampleWeight;
    maxScaledSampleWeight_  = stateHistory.maxScaledSampleWeight;
    numUpdates_             = stateHistory.numUpdates;
}

/* Broadcast the bias state over the MPI ranks in this simulation. */
void BiasState::broadcast(const t_commrec *cr)
{
    std::vector<PointState> pointsTmp = move(points_);
    GMX_RELEASE_ASSERT(points_.capacity() == 0, "Can only broadcast the state when points has no memory allocated");

    /* Broadcast the state memory with zero capacity std::vector */
    gmx_bcast(sizeof(BiasState), this, cr);

    points_ = move(pointsTmp);

    gmx_bcast(points_.size()*sizeof(PointState), points_.data(), cr);
}


/* Partition sampling domain. */
double BiasState::partitionDomain(const Grid &grid,
                                  int         pointMin,
                                  int         pointMax)
{
    awh_ivec    domain_imin_dim, domain_imax_dim, npoints_dim;

    /* Convert linear index to multidimensional index */
    for (int d = 0; d < grid.ndim(); d++)
    {
        npoints_dim[d] = grid.axis(d).numPoints();
    }
    linear_array_index_to_multidim(pointMin, grid.ndim(), npoints_dim, domain_imin_dim);
    linear_array_index_to_multidim(pointMax, grid.ndim(), npoints_dim, domain_imax_dim);

    double targetSum = 0.;
    for (size_t m = 0; m < points_.size(); m++)
    {
        PointState &pointState = points_[m];
        for (int d = 0; d < grid.ndim(); d++)
        {
            int index_d = grid.point(m).index[d];
            if (index_d < domain_imin_dim[d] ||
                index_d > domain_imax_dim[d])
            {
                pointState.setTargetToZero();
            }
        }
        targetSum += pointState.target();
    }

    /* Renormalize target distribution to 1.
     * NOTE: The total histogram size in state still needs to be scaled down.
     */
    double invTargetSum = 1/targetSum;
    for (auto &pointState : points_)
    {
        pointState.scaleTarget(invTargetSum);
    }

    histSize_ *= targetSum;

    return histSize_;
}

/* Convolves the PMF and sets the initial free energy to its convolution. */
void BiasState::setFreeEnergyToConvolvedPmf(const std::vector<DimParams>  &dimParams,
                                            const Grid                    &grid,
                                            const BiasParams              &params,
                                            const gmx_multisim_t          *ms)
{
    std::vector<float> convolvedPmf;

    calculateConvolvedPmf(dimParams, grid, params, points_, ms, &convolvedPmf);

    for (size_t m = 0; m < points_.size(); m++)
    {
        points_[m].setFreeEnergy(convolvedPmf[m]);
    }
}

/*! \brief
 * Find trailing data rows containing only zeros.
 *
 * \param[in] data    2D data array.
 * \param[in] nrows   Number of rows in array.
 * \param[in] ncols   Number of cols in array.
 * \returns the number of trailing zero rows.
 */
static int findTrailingZeroRows(const double* const *data, int nrows, int ncols)
{
    int nZeroRows = 0;
    for (int m = nrows - 1; m >= 0; m--)
    {
        bool rowIsZero = true;
        for (int d = 0; d < ncols; d++)
        {
            if (data[d][m] != 0)
            {
                rowIsZero = false;
                break;
            }
        }

        if (!rowIsZero)
        {
            /* At a row with non-zero data */
            break;
        }
        else
        {
            /* Still at a zero data row, keep checking rows higher up. */
            nZeroRows++;
        }
    }

    return nZeroRows;
}

/*! \brief
 * Initializes the PMF and target with data read from an input table.
 *
 * \param[in] dimParams       The dimension parameters.
 * \param[in] grid            The grid.
 * \param[in] numBias         Number of biases.
 * \param[in] biasIndex       The index of the bias.
 * \param[in,out] pointState  The state of the points in this bias.
 */
static void readUserPmfAndTargetDistribution(const std::vector<DimParams> &dimParams,
                                             const Grid                   &grid,
                                             int                           numBias,
                                             int                           biasIndex,
                                             std::vector<PointState>      *pointState)
{
    char          filename[STRLEN];

    /* Read the PMF and target distribution.
       From the PMF, the convolved PMF, or the reference value free energy, can be calculated
       base on the force constant. The free energy and target together determine the bias.
     */
    if (numBias == 1)
    {
        sprintf(filename, "awh-init.xvg");
    }
    else
    {
        sprintf(filename, "awh%d-init.xvg", biasIndex + 1);
    }

    char buf[STRLEN];
    sprintf(buf,
            "%s is expected in the following format. "
            "The first ndim column(s) should contain the coordinate values for each point, "
            " each column containing values of one dimension (in ascending order). "
            "For a multidimensional coordinate, points should be listed "
            "in the order obtained by traversing lower dimensions first. "
            "E.g. for two-dimensional grid of size nxn: "
            "(1, 1), (1, 2),..., (1, n), (2, 1), (2, 2), ..., , (n, n - 1), (n, n). "
            "Column ndim +  1 should contain the PMF value for each coordinate value. "
            "The target distribution values should be in column ndim + 2  or column ndim + 5. "
            "Make sure there input file ends with a new line but has no trailing new lines.",
            filename);
    char correctFormatMessage[STRLEN];
    sprintf(correctFormatMessage, "%s", wrap_lines(buf, linewidth, indent, FALSE));

    double  **data;
    int       nrows, ncols;
    nrows = read_xvg(filename, &data, &ncols);

    /* Check basic data properties here. Grid takes care of more complicated things. */

    if (nrows <= 0)
    {
        gmx_fatal(FARGS, "%s is empty!.\n\n%s", filename, correctFormatMessage);
    }

    /* Less than 2 points is not useful for PMF or target. */
    if (nrows <  2)
    {
        gmx_fatal(FARGS, "%s contains too few data points (%d)."
                  "The minimum number of points is 2.",
                  filename, nrows);
    }

    /* Make sure there are enough columns of data.

       Two formats are allowed. Either with columns  {coords, PMF, target} or
       {coords, PMF, x, y, z, target, ...}. The latter format is allowed since that
       is how AWH output is written (x, y, z being other AWH variables). For this format,
       trailing columns are ignored.
     */
    int colindex_target;
    int ncols_min    = dimParams.size() + 2;
    int colindex_pmf = dimParams.size();
    if (ncols == ncols_min)
    {
        colindex_target = colindex_pmf + 1;
    }
    else
    {
        colindex_target = colindex_pmf + 4;
    }

    if (ncols < ncols_min)
    {
        gmx_fatal(FARGS, "The number of columns in %s (%d) should be at least %d."
                  "\n\n%s",
                  filename, correctFormatMessage);
    }

    /* read_xvg can give trailing zero data rows for trailing new lines in the input. We allow 1 zero row,
       since this could be real data. But multiple trailing zero rows cannot correspond to valid data. */
    int nZeroRows = findTrailingZeroRows(data, nrows, ncols);
    if (nZeroRows > 1)
    {
        gmx_fatal(FARGS, "Found %d trailing zero data rows in %s. Please remove trailing empty lines and try again.",
                  nZeroRows, filename);
    }

    /* Convert from user units to internal units before sending the data of to grid. */
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        double scalingFactor = dimParams[d].scaleUserInputToInternal(1);
        if (scalingFactor == 1)
        {
            continue;
        }
        for (size_t m = 0; m < pointState->size(); m++)
        {
            data[d][m] *= scalingFactor;
        }
    }

    /* Get a data point for each AWH grid point so that they all get data. */
    std::vector<int> grid2data_index;
    grid2data_index.resize(grid.numPoints());
    mapGridToDatagrid(&grid2data_index, data, nrows, filename, grid, correctFormatMessage);

    /* Extract the data for each grid point */
    bool target_is_zero = true;
    for (size_t m = 0; m < pointState->size(); m++)
    {
        double target;

        (*pointState)[m].setLogPmfsum(-data[colindex_pmf][grid2data_index[m]]);
        target = data[colindex_target][grid2data_index[m]];

        /* Check if the values are allowed. */
        if (target < 0)
        {
            gmx_fatal(FARGS, "Target distribution weight at point %d (%g) in %s is negative.",
                      m, target, filename);
        }
        if (target > 0)
        {
            target_is_zero = false;
        }
        (*pointState)[m].setTargetConstantWeight(target);
    }

    if (target_is_zero)
    {
        gmx_fatal(FARGS, "The target weights given in column %d in %s are all 0",
                  filename, colindex_target);
    }

    /* Free the arrays. */
    for (int m = 0; m < ncols; m++)
    {
        sfree(data[m]);
    }
    sfree(data);
}

/* Normalize the PMF histogram. */
void BiasState::normalizePmf(int numSharingSims)
{
    /* The normalization of the PMF estimate matters because it determines how big effect the next sample has.
       Approximately (for large enough force constant) we should have:
       sum_x(exp(-pmf(x)) = nsamples*sum_xref(exp(-f(xref)).
     */

    /* Calculate the normalization factor, i.e. divide by the pmf sum, multiply by the number of samples and the f sum */
    double expSumPmf = 0;
    double expSumF   = 0;
    for (auto &pointState : points_)
    {
        if (pointState.inTargetRegion())
        {
            expSumPmf += std::exp( pointState.logPmfsum());
            expSumF   += std::exp(-pointState.freeEnergy());
        }
    }
    double numSamples = histSize_/numSharingSims;

    /* Renormalize */
    double logRenorm = std::log(numSamples*expSumF/expSumPmf);
    for (auto &pointState : points_)
    {
        if (pointState.inTargetRegion())
        {
            pointState.setLogPmfsum(pointState.logPmfsum() + logRenorm);
        }
    }
}

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
void BiasState::initGridPointState(const awh_bias_params_t       &awhBiasParams,
                                   const std::vector<DimParams>  &dimParams,
                                   const Grid                    &grid,
                                   const BiasParams              &params,
                                   int                            numBias,
                                   const gmx_multisim_t          *ms)
{
    /* Modify PMF, free energy and the constant target distribution factor
     * to user input values if there is data given.
     */
    if (awhBiasParams.bUser_data)
    {
        readUserPmfAndTargetDistribution(dimParams, grid, numBias, params.biasIndex, &points_);
        setFreeEnergyToConvolvedPmf(dimParams, grid, params, ms);
    }

    /* The local Boltzmann distribution is special because the target distribution is updated as a function of the reference weighthistogram. */
    GMX_RELEASE_ASSERT((params.eTarget != eawhtargetLOCALBOLTZMANN) ||
                       (params.eTarget == eawhtargetLOCALBOLTZMANN && points_[0].weightsumRef() != 0),
                       "AWH reference weight histogram not initialized properly with local Boltzmann target distribution.");

    updateTarget(&points_, params);

    for (auto &pointState : points_)
    {
        if (pointState.inTargetRegion())
        {
            pointState.updateBias();
        }
        else
        {
            /* Note that for zero target this is a value that represents -infinity but should not be used for biasing. */
            pointState.setTargetToZero();
        }
    }

    /* Set the initial reference weighthistogram. */
    const double histogramSize = histSize_;
    for (auto &pointState : points_)
    {
        pointState.setInitialReferenceWeightHistogram(histogramSize);
    }

    /* Make sure the pmf is normalized consistently with the histogram size.
       Note: the target distribution and free energy need to be set here. */
    normalizePmf(params.numSharedUpdate);
}

BiasState::BiasState(const awh_bias_params_t      &awhBiasParams,
                     double                        histSizeInitial,
                     const std::vector<DimParams> &dimParams,
                     const Grid                   &grid) :
    numUpdates_(0),
    histSize_(histSizeInitial),
    inInitialStage_(awhBiasParams.eGrowth == eawhgrowthEXP_LINEAR),
    equilibrateHistogram_(awhBiasParams.equilibrateHistogram),
    /* The initial sample weight is set to 1 and we keep the logarithm. */
    scaledSampleWeight_(0),
    maxScaledSampleWeight_(0)
{
    awh_dvec             coordValueInit;

    for (size_t d = 0; d < dimParams.size(); d++)
    {
        coordValueInit[d] = dimParams[d].scaleUserInputToInternal(awhBiasParams.dim_params[d].coord_value_init);
    }

    /* Set initial coordinate reference value to the one closest to the initial reference value given in pull.
       More correctly one would sample from the biased distribution, but it doesn't really matter. */
    gridpointIndex_ = grid.nearestIndex(coordValueInit);
    refGridpoint_   = gridpointIndex_;

    /* The minimum and maximum multidimensional point indices that are affected by the next update */
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        int index_d          = grid.point(gridpointIndex_).index[d];
        originUpdatelist_[d] = index_d;
        endUpdatelist_[d]    = index_d;
    }

    /* Initialize free energy functions and biases */
    points_.resize(grid.numPoints());
}

/* Returns if to do checks, only returns true at free-energy update steps. */
bool isCheckStep(const BiasParams &params, const std::vector<PointState> &pointState, gmx_int64_t step)
{
    int numStepsUpdateFreeEnergy = params.numSamplesUpdateFreeEnergy*params.numStepsSampleCoord;
    int numStepsCheck            = (1 + pointState.size()/params.numSamplesUpdateFreeEnergy)*numStepsUpdateFreeEnergy;

    return step % numStepsCheck == 0;
}
