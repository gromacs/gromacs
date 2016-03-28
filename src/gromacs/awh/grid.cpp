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
 * Implements functions in grid.h.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "grid.h"

#include <assert.h>

#include <cmath>
#include <cstring>

#include <algorithm>

#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

/*! \brief
 * Modify x so that it is periodic in [-period/2, +period/2).
 *
 * x is modified by shifting its value by a +/- a period if
 * needed. Thus, it is assumed that x is at most one period
 * away from this interval. For period = 0, x is not modified.
 *
 * \param[in,out] x       Pointer to the value to modify.
 * \param[in]     period  The period, or 0 if not periodic.
 */
void centerPeriodicValueAroundZero(double *x,
                                   double  period)
{
    GMX_ASSERT(period >= 0, "Periodic should not be negative");

    const double halfPeriod = period*0.5;

    if (*x >= halfPeriod)
    {
        *x -= period;
    }
    else if (*x < -halfPeriod)
    {
        *x += period;
    }
}

/*! \brief
 * If period>0, retrun x so that it is periodic in [0, period), else return x.
 *
 * Return x is shifted its value by a +/- a period, if
 * needed. Thus, it is assumed that x is at most one period
 * away from this interval. For this domain and period > 0
 * this is equivalent to x = x % period. For period = 0,
 * x is not modified.
 *
 * \param[in,out] x       Pointer to the value to modify, should be >= 0.
 * \param[in]     period  The period, or 0 if not periodic.
 * \returns for period>0: index value witin [0, period), otherwise: \p x.
 */
int indexWithinPeriod(int x,
                      int period)
{
    GMX_ASSERT(period >= 0, "Periodic should not be negative");

    if (period == 0)
    {
        return x;
    }

    GMX_ASSERT(x > -period && x < 2*period, "x should not be more shifted by more than one period");

    if (x >= period)
    {
        return x - period;
    }
    else if (x < 0)
    {
        return x + period;
    }
    else
    {
        return x;
    }
}

/*! \brief
 * Get the length of the interval (origin, end).
 *
 * This returns the distance obtained by connecting the origin point to
 * the end point in the positive direction. Note that this is generally
 * not the shortest distance. For period > 0, both origin and
 * end are expected to take values in the same periodic interval,
 * ie. |origin - end| < period.
 *
 * \param[in] origin    Start value of the interval.
 * \param[in] end       End value of the interval.
 * \param[in] period    The period, or 0 if not periodic.
 * \returns the interval length from origin to end.
 */
double getIntervalLengthPeriodic(double origin,
                                 double end,
                                 double period)
{
    double length = end - origin;
    if (length < 0)
    {
        /* The interval wraps around the +/- boundary which has a discontinuous jump of -period. */
        length += period;
    }

    GMX_RELEASE_ASSERT(length >= 0, "Negative AWH grid axis length.");
    GMX_RELEASE_ASSERT(period == 0 || length <= period, "Interval length longer than period.");

    return length;
}

/*! \brief
 * Get the deviation x - x0.
 *
 * For period > 0, the deviation with minimum absolute value is returned,
 * i.e. with a value in the interval [-period/2, +period/2).
 * Also for period > 0, it is assumed that |x - x0| < period.
 *
 * \param[in] x        From value.
 * \param[in] x0       To value.
 * \param[in] period   The period, or 0 if not periodic.
 * \returns the deviation from x to x0.
 */
double getDeviationPeriodic(double x,
                            double x0,
                            double period)
{
    double dev = x - x0;

    if (period > 0)
    {
        centerPeriodicValueAroundZero(&dev, period);
    }

    return dev;
}

}   // namespace

double getDeviationFromPointAlongGridAxis(const Grid &grid,
                                          int         dimIndex,
                                          int         pointIndex,
                                          double      value)
{
    double coordValue = grid.point(pointIndex).coordValue[dimIndex];

    return getDeviationPeriodic(value, coordValue, grid.axis(dimIndex).period());
}

void linearArrayIndexToMultiDim(int indexLinear, int numDimensions, const awh_ivec numPointsDim, awh_ivec indexMulti)
{
    for (int d = 0; d < numDimensions; d++)
    {
        int stride = 1;

        /* Workaround for bug in clang */
#ifndef __clang_analyzer__
        for (int k = d + 1; k < numDimensions; k++)
        {
            stride *= numPointsDim[k];
        }
#endif

        indexMulti[d] = indexLinear/stride;
        indexLinear  -= indexMulti[d]*stride;
    }
}

void linearGridindexToMultiDim(const Grid &grid,
                               int         indexLinear,
                               awh_ivec    indexMulti)
{
    awh_ivec numPointsDim;

    for (int d = 0; d < grid.numDimensions(); d++)
    {
        numPointsDim[d] = grid.axis(d).numPoints();
    }

    linearArrayIndexToMultiDim(indexLinear, grid.numDimensions(), numPointsDim, indexMulti);
}


int multiDimArrayIndexToLinear(const awh_ivec indexMulti,
                               int            numDimensions,
                               const awh_ivec numPointsDim)
{
    int stride      = 1;
    int indexLinear = 0;
    for (int d = numDimensions - 1; d >= 0; d--)
    {
        indexLinear += stride*indexMulti[d];
        stride      *= numPointsDim[d];
    }

    return indexLinear;
}

namespace
{

/*! \brief Convert a multidimensional grid point index to a linear one.
 *
 * \param[in] axis       The grid axes.
 * \param[in] indexMulti Multidimensional grid point index to convert to a linear one.
 * \returns the linear index.
 */
int multiDimGridIndexToLinear(const std::vector<GridAxis> &axis,
                              const awh_ivec               indexMulti)
{
    awh_ivec numPointsDim = { 0 };

    for (size_t d = 0; d < axis.size(); d++)
    {
        numPointsDim[d] = axis[d].numPoints();
    }

    return multiDimArrayIndexToLinear(indexMulti, axis.size(), numPointsDim);
}

}   // namespace

int multiDimGridIndexToLinear(const Grid     &grid,
                              const awh_ivec  indexMulti)
{
    return multiDimGridIndexToLinear(grid.axis(), indexMulti);
}

namespace
{

/*! \brief
 * Take a step in a multidimensional array.
 *
 * The multidimensional index gives the starting point to step from. Dimensions are
 * stepped through in order of decreasing dimensional index such that the index is
 * incremented in the highest dimension possible. If the starting point is the end
 * of the array, a step cannot be taken and the index is not modified.
 *
 * \param[in] numDim        Number of dimensions of the array.
 * \param[in] numPoints     Vector with the number of points along each dimension.
 * \param[in,out] indexDim  Multidimensional index, each with values in [0, numPoints[d] - 1].
 * \returns true if a step was taken, false if not.
 */
bool stepInMultiDimArray(int            numDim,
                         const awh_ivec numPoints,
                         awh_ivec       indexDim)
{
    bool haveStepped = false;

    for (int d = numDim - 1; d >= 0 && !haveStepped; d--)
    {
        if (indexDim[d] < numPoints[d] - 1)
        {
            /* Not at a boundary, just increase by 1. */
            indexDim[d]++;
            haveStepped = true;
        }
        else
        {
            /* At a boundary. If we are not at the end of the array,
               reset the index and check if we can step in higher dimensions */
            if (d > 0)
            {
                indexDim[d] = 0;
            }
        }
    }

    return haveStepped;
}

/*! \brief
 * Transforms a grid point index to to the multidimensional index of a subgrid.
 *
 * The subgrid is defined by the location of its origin and the number of points
 * along each dimension. The index transformation thus consists of a projection
 * of the linear index onto each dimension, followed by a translation of the origin.
 * The subgrid may have parts that don't overlap with the grid. E.g. the origin
 * vector can have negative components meaning the origin lies outside of the grid.
 * However, the given point needs to be both a grid and subgrid point.
 *
 * Periodic boundaries are taken care of by wrapping the subgrid around the grid.
 * Thus, for periodic dimensions the number of subgrid points need to be less than
 * the number of points in a period to prevent problems of wrapping around.
 *
 * \param[in]     grid            The grid.
 * \param[in]     subgridOrigin   Vector locating the subgrid origin relative to the grid origin.
 * \param[in]     subgridNpoints  The number of subgrid points in each dimension.
 * \param[in]     point           Grid point to get subgrid index for.
 * \param[in,out] subgridIndex    Subgrid multidimensional index.
 */
void gridToSubgridIndex(const Grid     &grid,
                        const awh_ivec  subgridOrigin,
                        const awh_ivec  subgridNpoints,
                        int             point,
                        awh_ivec        subgridIndex)
{
    /* Get the subgrid index of the given grid point, for each dimension. */
    for (int d = 0; d < grid.numDimensions(); d++)
    {
        /* The multidimensional grid point index relative to the subgrid origin. */
        subgridIndex[d] =
            indexWithinPeriod(grid.point(point).index[d] - subgridOrigin[d],
                              grid.axis(d).numPointsInPeriod());

        /* The given point should be in the subgrid. */
        GMX_RELEASE_ASSERT((subgridIndex[d] >= 0) && (subgridIndex[d] < subgridNpoints[d]),
                           "Attempted to convert an AWH grid point index not in subgrid to out of bounds subgrid index");
    }
}

/*! \brief
 * Transform a multidimensional subgrid index to a grid point index.
 *
 * If the given subgrid point is not a grid point the transformation will not be successful
 * and the grid point index will not be set. Periodic boundaries are taken care of by
 * wrapping the subgrid around the grid.
 *
 * \param[in]     grid           The grid.
 * \param[in]     subgridOrigin  Vector locating the subgrid origin relative to the grid origin.
 * \param[in]     subgridIndex   Subgrid multidimensional index to get grid point index for.
 * \param[in,out] gridIndex      Grid point index.
 * \returns true if the transformation was successful.
 */
bool subgridToGridIndex(const Grid     &grid,
                        const awh_ivec  subgridOrigin,
                        const awh_ivec  subgridIndex,
                        int            *gridIndex)
{
    awh_ivec globalIndexDim;

    /* Check and apply boundary conditions for each dimension */
    for (int d = 0; d < grid.numDimensions(); d++)
    {
        /* Transform to global multidimensional indexing by adding the origin */
        globalIndexDim[d] = subgridOrigin[d] + subgridIndex[d];

        /* The local grid is allowed to stick out on the edges of the global grid. Here the boundary conditions are applied.*/
        if (globalIndexDim[d] < 0 ||  globalIndexDim[d] > grid.axis(d).numPoints() - 1)
        {
            /* Try to wrap around if periodic. Otherwise, the transformation failed so return. */
            if (!grid.axis(d).isPeriodic())
            {
                return false;
            }

            /* The grid might not contain a whole period. Can only wrap around if this gap is not too large. */
            int gap = grid.axis(d).numPointsInPeriod() - grid.axis(d).numPoints();

            int bridge;
            int numWrapped;
            if (globalIndexDim[d] < 0)
            {
                bridge     = -globalIndexDim[d];
                numWrapped = bridge - gap;
                if (numWrapped > 0)
                {
                    globalIndexDim[d] = grid.axis(d).numPoints() - numWrapped;
                }
            }
            else
            {
                bridge     = globalIndexDim[d] - (grid.axis(d).numPoints() - 1);
                numWrapped = bridge - gap;
                if (numWrapped > 0)
                {
                    globalIndexDim[d] = numWrapped - 1;
                }
            }

            if (numWrapped <= 0)
            {
                return false;
            }
        }
    }

    /* Translate from multidimensional to linear indexing and set the return value */
    (*gridIndex) = multiDimGridIndexToLinear(grid, globalIndexDim);

    return true;
}

}   // namespace

bool advancePointInSubgrid(const Grid     &grid,
                           const awh_ivec  subgridOrigin,
                           const awh_ivec  subgridNumPoints,
                           int            *gridPointIndex)
{
    /* Initialize the subgrid index to the subgrid origin. */
    awh_ivec subgridIndex = { 0 };

    /* Get the subgrid index of the given grid point index. */
    if (*gridPointIndex >= 0)
    {
        gridToSubgridIndex(grid, subgridOrigin, subgridNumPoints, *gridPointIndex, subgridIndex);
    }
    else
    {
        /* If no grid point is given we start at the subgrid origin (which subgridIndex is initialized to).
           If this is a valid grid point then we're done, otherwise keep looking below. */
        /* TODO: separate into a separate function (?) */
        if (subgridToGridIndex(grid, subgridOrigin, subgridIndex, gridPointIndex))
        {
            return true;
        }
    }

    /* Traverse the subgrid and look for the first point that is also in the grid. */
    while (stepInMultiDimArray(grid.numDimensions(), subgridNumPoints, subgridIndex))
    {
        /* If this is a valid grid point, the grid point index is updated.*/
        if (subgridToGridIndex(grid, subgridOrigin, subgridIndex, gridPointIndex))
        {
            return true;
        }
    }

    return false;
}

/*! \brief
 * Returns the point distance between from value x to value x0 along the given axis.
 *
 * Note that the returned distance may be negative or larger than the
 * number of points in the axis. For a periodic axis, the distance is chosen
 * to be in [0, period), i.e. always positive but not the shortest one.
 *
 * \param[in]  axis   Grid axis.
 * \param[in]  x      From value.
 * \param[in]  x0     To value.
 * \returns (x - x0) in number of points.
 */
static int pointDistanceAlongAxis(const GridAxis &axis, double x, double x0)
{
    int distance = 0;

    if (axis.spacing() > 0)
    {
        /* Get the real-valued distance. For a periodic axis, the shortest one. */
        double period = axis.period();
        double dx     = getDeviationPeriodic(x, x0, period);

        /* Transform the distance into a point distance.
           Shift by +0.5 so we can use floor or integer casting below to get the integer index */
        distance = static_cast<int>(floor(dx/axis.spacing() + 0.5));

        /* If periodic, shift the point distance to be in [0, period) */
        distance = indexWithinPeriod(distance, axis.numPointsInPeriod());
    }

    return distance;
}

/*! \brief
 * Query if a value is in range of the grid.
 *
 * \param[in] value   Value to check.
 * \param[in] axis    The grid axes.
 * \returns true if the value is in the grid.
 */
static bool valueIsInGrid(const awh_dvec               value,
                          const std::vector<GridAxis> &axis)
{
    /* For each dimension get the one-dimensional index and check if it is in range. */
    for (size_t d = 0; d < axis.size(); d++)
    {
        /* The index is computed as the point distance from the origin. */
        int index = pointDistanceAlongAxis(axis[d], value[d], axis[d].origin());

        if (!(index >= 0 && index < axis[d].numPoints()))
        {
            return false;
        }
    }

    return true;
}

bool Grid::covers(const awh_dvec value) const
{
    return valueIsInGrid(value, axis());
}

int GridAxis::nearestIndex(double value) const
{
    /* Get the point distance to the origin. This may by an out of index range for the axis. */
    int index = pointDistanceAlongAxis(*this, value, origin_);

    if (index < 0 || index >= numPoints_)
    {
        if (isPeriodic())
        {
            GMX_RELEASE_ASSERT(index >= 0 && index < numPointsInPeriod_,
                               "Index not in periodic interval 0 for AWH periodic axis");
            int endDistance    = (index - (numPoints_ - 1));
            int originDistance = (numPointsInPeriod_ - index);
            index = originDistance < endDistance ? 0 : numPoints_ - 1;
        }
        else
        {
            index = (index < 0) ? 0 : (numPoints_ - 1);
        }
    }

    return index;
}

/*! \brief
 * Map a value to the nearest point in the grid.
 *
 * \param[in] value  Value.
 * \param[in] axis   The grid axes.
 * \returns the point index nearest to the value.
 */
static int getNearestIndexInGrid(const awh_dvec               value,
                                 const std::vector<GridAxis> &axis)
{
    awh_ivec indexMulti;

    /* If the index is out of range, modify it so that it is in range by choosing the nearest point on the edge. */
    for (size_t d = 0; d < axis.size(); d++)
    {
        indexMulti[d] = axis[d].nearestIndex(value[d]);
    }

    return multiDimGridIndexToLinear(axis, indexMulti);
}

int Grid::nearestIndex(const awh_dvec value) const
{
    return getNearestIndexInGrid(value, axis());
}

namespace
{

/*! \brief
 * Find and set the neighbors of a grid point.
 *
 * The search space for neighbors is a subgrid with size set by a scope cutoff.
 * In general not all point within scope will be valid grid points.
 *
 * \param[in]     pointIndex           Grid point index.
 * \param[in]     grid                 The grid.
 * \param[in,out] neighborIndexArray   Array to fill with neighbor indices.
 */
void setNeighborsOfGridPoint(int               pointIndex,
                             const Grid       &grid,
                             std::vector<int> *neighborIndexArray)
{
    const int c_maxNeighborsAlongAxis = 1 + 2*static_cast<int>(Grid::c_numPointsPerSigma*Grid::c_scopeCutoff);

    awh_ivec  numCandidates           = {0};
    awh_ivec  subgridOrigin           = {0};
    for (int d = 0; d < grid.numDimensions(); d++)
    {
        /* The number of candidate points along this dimension is given by the scope cutoff. */
        numCandidates[d] = std::min(c_maxNeighborsAlongAxis,
                                    grid.axis(d).numPoints());

        /* The origin of the subgrid to search */
        int centerIndex  = grid.point(pointIndex).index[d];
        subgridOrigin[d] = centerIndex - numCandidates[d]/2;
    }

    /* Find and set the neighbors */
    int  neighborIndex = -1;
    bool aPointExists  = true;

    /* Keep looking for grid points while traversing the subgrid. */
    while (aPointExists)
    {
        /* The point index is updated if a grid point was found. */
        aPointExists = advancePointInSubgrid(grid, subgridOrigin, numCandidates, &neighborIndex);

        if (aPointExists)
        {
            neighborIndexArray->push_back(neighborIndex);
        }
    }
}

}   // namespace

void Grid::initPoints()
{
    awh_ivec     numPointsDimWork = { 0 };
    awh_ivec     indexWork        = { 0 };

    for (size_t d = 0; d < axis_.size(); d++)
    {
        /* Temporarily gather the number of points in each dimension in one array */
        numPointsDimWork[d] = axis_[d].numPoints();
    }

    for (auto &point : point_)
    {
        for (size_t d = 0; d < axis_.size(); d++)
        {
            point.coordValue[d] = axis_[d].origin() + indexWork[d]*axis_[d].spacing();

            if (axis_[d].period() > 0)
            {
                /* Do we always want the values to be centered around 0 ? */
                centerPeriodicValueAroundZero(&point.coordValue[d], axis_[d].period());
            }

            point.index[d] = indexWork[d];
        }

        stepInMultiDimArray(axis_.size(), numPointsDimWork, indexWork);
    }
}

GridAxis::GridAxis(double origin, double end,
                   double period, double pointDensity) :
    origin_(origin),
    period_(period)
{
    length_ = getIntervalLengthPeriodic(origin_, end, period_);

    /* Automatically determine number of points based on the user given endpoints
       and the expected fluctuations in the umbrella. */
    if (length_ == 0)
    {
        numPoints_ = 1;
    }
    else if (pointDensity == 0)
    {
        numPoints_ = 2;
    }
    else
    {
        double lengthInPoints = length_*pointDensity;

        numPoints_            = 1 + static_cast<int>(std::ceil(lengthInPoints));
    }

    /* Set point spacing based on the number of points */
    if (isPeriodic())
    {
        /* Set the grid spacing so that a period is matched exactly by an integer number of points.
           The number of points in a period is equal to the number of grid spacings in a period
           since the endpoints are connected.  */
        numPointsInPeriod_ = length_ > 0 ? static_cast<int>(std::ceil(period/length_*(numPoints_ - 1))) : 1;
        spacing_           = period_/numPointsInPeriod_;

        /* Modify the number of grid axis points to be compatible with the period dependent spacing. */
        numPoints_         = std::min(static_cast<int>(round(length_/spacing_)) + 1,
                                      numPointsInPeriod_);
    }
    else
    {
        numPointsInPeriod_ = 0;
        spacing_           = numPoints_ > 1 ? length_/(numPoints_ - 1) : 0;
    }
}

GridAxis::GridAxis(double origin, double end,
                   double period, int numPoints) :
    origin_(origin),
    period_(period),
    numPoints_(numPoints)
{
    length_            = getIntervalLengthPeriodic(origin_, end, period_);
    spacing_           = numPoints_ > 1 ? length_/(numPoints_ - 1) : period_;
    numPointsInPeriod_ = static_cast<int>(std::round(period_/spacing_));
}

Grid::Grid(const std::vector<DimParams> &dimParams,
           const AwhDimParams           *awhDimParams)
{
    /* Define the discretization along each dimension */
    awh_dvec period;
    int      numPoints = 1;
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        double origin = dimParams[d].scaleUserInputToInternal(awhDimParams[d].origin);
        double end    = dimParams[d].scaleUserInputToInternal(awhDimParams[d].end);
        period[d]     = dimParams[d].scaleUserInputToInternal(awhDimParams[d].period);
        static_assert(c_numPointsPerSigma >= 1.0, "The number of points per sigma should be at least 1.0 to get a uniformly covering the reaction using Gaussians");
        double pointDensity = std::sqrt(dimParams[d].betak)*c_numPointsPerSigma;
        axis_.push_back(GridAxis(origin, end, period[d], pointDensity));
        numPoints *= axis_[d].numPoints();
    }

    point_.resize(numPoints);

    /* Set their values */
    initPoints();

    /* Keep a neighbor list for each point.
     * Note: could also generate neighbor list only when needed
     * instead of storing them for each point.
     */
    for (size_t m = 0; m < point_.size(); m++)
    {
        std::vector<int> *neighbor = &point_[m].neighbor;

        setNeighborsOfGridPoint(m, *this, neighbor);
    }
}

void mapGridToDataGrid(std::vector<int>    *gridpointToDatapoint,
                       const double* const *data,
                       int                  numDataPoints,
                       const std::string   &dataFilename,
                       const Grid          &grid,
                       const std::string   &correctFormatMessage)
{
    /* Transform the data into a grid in order to map each grid point to a data point
       using the grid functions. */
    std::vector<GridAxis> axis_;

    /* Count the number of points for each dimension. Each dimension
       has its own stride. */
    int              stride           = 1;
    int              numPointsCounted = 0;
    std::vector<int> numPoints(grid.numDimensions());
    for (int d = grid.numDimensions() - 1; d >= 0; d--)
    {
        int    numPointsInDim = 0;
        int    pointIndex     = 0;
        double firstValue     = data[d][pointIndex];
        do
        {
            numPointsInDim++;
            pointIndex       += stride;
        }
        while (pointIndex < numDataPoints &&
               !gmx_within_tol(firstValue, data[d][pointIndex], GMX_REAL_EPS));

        /* The stride in dimension dimension d - 1 equals the number of points
           dimension d. */
        stride = numPointsInDim;

        numPointsCounted = (numPointsCounted == 0) ? numPointsInDim : numPointsCounted*numPointsInDim;

        numPoints[d]     = numPointsInDim;
    }

    if (numPointsCounted != numDataPoints)
    {
        std::string mesg = gmx::formatString("Could not extract data properly from %s. Wrong data format?"
                                             "\n\n%s",
                                             dataFilename.c_str(), correctFormatMessage.c_str());
        GMX_THROW(InvalidInputError(mesg));
    }

    /* The data grid has the data that was read and the properties of the AWH grid */
    for (int d = 0; d < grid.numDimensions(); d++)
    {
        axis_.push_back(GridAxis(data[d][0], data[d][numDataPoints - 1],
                                 grid.axis(d).period(), numPoints[d]));
    }

    /* Map each grid point to a data point. No interpolation, just pick the nearest one.
     * It is assumed that the given data is uniformly spaced for each dimension.
     */
    for (size_t m = 0; m < grid.numPoints(); m++)
    {
        /* We only define what we need for the datagrid since it's not needed here which is a bit ugly */

        if (!valueIsInGrid(grid.point(m).coordValue, axis_))
        {
            std::string mesg =
                gmx::formatString("%s does not contain data for all coordinate values. "
                                  "Make sure your input data covers the whole sampling domain "
                                  "and is correctly formatted. \n\n%s",
                                  dataFilename.c_str(), correctFormatMessage.c_str());
            GMX_THROW(InvalidInputError(mesg));
        }
        (*gridpointToDatapoint)[m] = getNearestIndexInGrid(grid.point(m).coordValue, axis_);
    }
}

} // namespace gmx
