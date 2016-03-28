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

#include "gmxpre.h"

#include "grid.h"

#include <assert.h>

#include <cmath>
#include <cstring>

#include <algorithm>

#include "gromacs/awh/awh.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "math.h"
#include "types.h"

/*! \brief
 * The point density per sigma of the Gaussian distribution in an umbrella.
 *
 * This value should be at least 1 to uniformly cover the reaction coordinate
 * range with density and having it larger than 1 does not add information.
 */
static const double c_pointDensity          = 1.0;

//! Cut-off in sigma for considering points, neglects 4e-8 of the density.
static const double c_scopeCutoff           = 5.5;

//! The maximum number of neighboring points to consider along an axis.
static const int    c_maxNeighborsAlongAxis = 1 + 2*static_cast<int>(c_pointDensity*c_scopeCutoff);


/*! \brief
 * Modify x so that it is periodic in [-period/2, +period/2).
 *
 * x is modified by shifting its value by a +/- a period if
 * needed. Thus, it is assumed that x is at most one period
 * away from this interval. For period = 0, x is not modified.
 *
 * \param[in,out] x    Pointer to the value to modify.
 * \param[in] period   The period, or 0 if not periodic.
 */
static void make_value_periodic(double *x, double period)
{
    const double halfperiod = period*0.5;

    if (*x >= halfperiod)
    {
        *x -= period;
    }
    else if (*x < -halfperiod)
    {
        *x += period;
    }
}

/*! \brief
 * Modify x so that it is periodic in [0, period).
 *
 * x is modified by shifting its value by a +/- a period if
 * needed. Thus, it is assumed that x is at most one period
 * away from this interval. For this domain and period > 0
 * this is equivalent to x = x % period. For period = 0,
 * x is not modified.
 *
 * \param[in,out] x    Pointer to the value to modify.
 * \param[in] period   The period, or 0 if not periodic.
 */
static void make_index_periodic(int *x, int period)
{
    if (*x >= period)
    {
        *x -= period;
    }
    else if (*x < 0)
    {
        *x += period;
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
static double get_interval_length_periodic(double origin, double end, double period)
{
    double L;

    L = end - origin;
    if (L < 0)
    {
        /* The interval wraps around the +/- boundary which has a discontinuous jump of -period. */
        L += period;
    }

    GMX_ASSERT(L >= 0, "Negative AWH grid axis length.");
    GMX_ASSERT(period == 0 || L <= period, "Interval length longer than period.");

    return L;
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
static double get_deviation_periodic(double x, double x0, double period)
{
    double dev = x - x0;

    if (period > 0)
    {
        make_value_periodic(&dev, period);
    }

    return dev;
}

/*! \brief
 * Query if the axis has periodic boundaries.
 *
 * \param[in] axis   The axis to query about.
 * \returns true if the axis is periodic.
 */
static bool gridaxis_is_periodic(const GridAxis *axis)
{
    return axis->npoints_period  > 0;
}

/*! \brief
 * Get the period of the grid along the given dimension.
 *
 * \param[in] grid      The axis to query about.
 * \param[in] dimindex  Dimensional index.
 * \returns the period.
 */
static double get_gridaxis_period(const Grid *grid, int dimindex)
{
    return grid->axis[dimindex].npoints_period*grid->axis[dimindex].spacing;
}

/* Get the deviation along one dimension from the given value to a point in the grid. */
double get_deviation_from_point_along_gridaxis(const Grid *grid, int dimindex, int pointindex, double value)
{
    double pointvalue = grid->point[pointindex].value[dimindex];

    return get_deviation_periodic(value, pointvalue, get_gridaxis_period(grid, dimindex));
}

/* Gets the length of the grid for the given dimension. */
double get_gridaxis_length(const Grid *grid, int dimindex)
{
    return get_interval_length_periodic(grid->axis[dimindex].origin,
                                        grid->axis[dimindex].end,
                                        get_gridaxis_period(grid, dimindex));
}

/* Convert a linear index to a multidimensional index. */
void linear_array_index_to_multidim(int index_linear, int ndim, const awh_ivec npoints_dim, awh_ivec index_multi)
{
    for (int d = 0; d < ndim; d++)
    {
        int stride = 1;

        for (int k = d + 1; k < ndim; k++)
        {
            stride *= npoints_dim[k];
        }

        index_multi[d] = index_linear/stride;
        index_linear  -= index_multi[d]*stride;
    }
}

/* Convert a linear grid point index to a multidimensional one. */
void linear_gridindex_to_multidim(const Grid *grid, int index_linear, awh_ivec index_multi)
{
    awh_ivec npoints_dim;

    for (int d = 0; d < grid->ndim; d++)
    {
        npoints_dim[d] = grid->axis[d].npoints;
    }

    linear_array_index_to_multidim(index_linear, grid->ndim, npoints_dim, index_multi);
}


/* Convert multidimensional array index to a linear one. */
int multidim_array_index_to_linear(const awh_ivec index_multi, int ndim, const awh_ivec npoints_dim)
{
    int stride       = 1;
    int index_linear = 0;
    /* Workaround for bug in CLANG */
#ifndef __clang_analyzer__
    for (int d = ndim - 1; d >= 0; d--)
    {
        index_linear += stride*index_multi[d];
        stride       *= npoints_dim[d];
    }
#endif

    return index_linear;
}

/* Convert a multidimensional grid point index to a linear one. */
int multidim_gridindex_to_linear(const Grid *grid, const awh_ivec index_multi)
{
    awh_ivec npoints_dim;

    for (int d = 0; d < grid->ndim; d++)
    {
        npoints_dim[d] = grid->axis[d].npoints;
    }

    return multidim_array_index_to_linear(index_multi, grid->ndim, npoints_dim);
}

/*! \brief
 * Take a step in a multidimensional array.
 *
 * The multidimensional index gives the starting point to step from. Dimensions are
 * stepped through in order of decreasing dimensional index such that the index is
 * incremented in the highest dimension possible. If the starting point is the end
 * of the array, a step cannot be taken and the index is not modified.
 *
 * \param[in] ndim            Number of dimensions of the array.
 * \param[in] npoints_dim     Vector with the number of points along each dimension.
 * \param[in,out] index_dim   Multidimensional index, each with values in [0, npoints_dim[d] - 1].
 * \returns true if a step was taken, false if not.
 */
static bool step_in_multidim_array(int ndim, const awh_ivec npoints_dim, awh_ivec index_dim)
{
    bool bStepped = false;

    for (int d = ndim - 1; d >= 0 && !bStepped; d--)
    {
        if (index_dim[d] < npoints_dim[d] - 1)
        {
            /* Not at a boundary, just increase by 1. */
            index_dim[d]++;
            bStepped = true;
        }
        else
        {
            /* At a boundary. If we are not at the end of the array,
               reset the index and check if we can step in higher dimensions */
            if (d > 0)
            {
                index_dim[d] = 0;
            }
        }
    }

    return bStepped;
}

/*! \brief
 * Transforms a grid point index to to the  multidimensional index of a subgrid.
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
 * \param[in] grid                  The grid.
 * \param[in] subgridOrigin         Vector locating the subgrid origin relative to the grid origin.
 * \param[in] subgridNpoints        The number of subgrid points in each dimension.
 * \param[in] point                 Grid point to get subgrid index for.
 * \param[in,out] subgridIndex      Subgrid multidimensional index.
 */
static void gridToSubgridIndex(const Grid *grid,
                               const awh_ivec subgridOrigin, const awh_ivec subgridNpoints,
                               int point, awh_ivec subgridIndex)
{
    /* Get the subgrid index of the given grid point, for each dimension. */
    for (int d = 0; d < grid->ndim; d++)
    {
        /* The multidimensional grid point index relative to the subgrid origin. */
        subgridIndex[d] = grid->point[point].index[d] - subgridOrigin[d];

        /* The subgrid wraps around periodic boundaries of the grid. */
        make_index_periodic(&subgridIndex[d], grid->axis[d].npoints_period);

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
 * \param[in] grid                  The grid.
 * \param[in] subgridOrigin         Vector locating the subgrid origin relative to the grid origin.
 * \param[in] subgridIndex          Subgrid multidimensional index to get grid point index for.
 * \param[in,out] gridIndex         Grid point index.
 * \returns true if the transformation was successful.
 */
static bool subgridToGridIndex(const Grid     *grid,
                               const awh_ivec  subgridOrigin,
                               const awh_ivec  subgridIndex,
                               int            *gridIndex)
{
    awh_ivec global_index_dim;

    /* Check and apply boundary conditions for each dimension */
    for (int d = 0; d < grid->ndim; d++)
    {
        /* Transform to global multidimensional indexing by adding the origin */
        global_index_dim[d] = subgridOrigin[d] + subgridIndex[d];

        /* The local grid is allowed to stick out on the edges of the global grid. Here the boundary conditions are applied.*/
        if (global_index_dim[d] < 0 ||  global_index_dim[d] > grid->axis[d].npoints - 1)
        {
            /* Try to wrap around if periodic. Otherwise, the transformation failed so return. */
            if (!gridaxis_is_periodic(&grid->axis[d]))
            {
                return false;
            }

            /* The grid might not contain a whole period. Can only wrap around if this gap is not too large. */
            int gap = grid->axis[d].npoints_period - grid->axis[d].npoints;

            int bridge, nwrapped;
            if (global_index_dim[d] < 0)
            {
                bridge   = -global_index_dim[d];
                nwrapped = bridge - gap;
                if (nwrapped > 0)
                {
                    global_index_dim[d] = grid->axis[d].npoints - nwrapped;
                }
            }
            else
            {
                bridge   = global_index_dim[d] - (grid->axis[d].npoints - 1);
                nwrapped = bridge - gap;
                if (nwrapped > 0)
                {
                    global_index_dim[d] = nwrapped - 1;
                }
            }

            if (nwrapped <= 0)
            {
                return false;
            }
        }
    }

    /* Translate from multidimensional to linear indexing and set the return value */
    (*gridIndex) = multidim_gridindex_to_linear(grid, global_index_dim);

    return true;
}

/* Find the next grid point in the subgrid given a starting point. */
bool get_next_point_in_local_grid(const Grid *grid,
                                  const awh_ivec subgridOrigin, const awh_ivec subgridNpoints,
                                  int *gridPointIndex)
{
    /* Initialize the subgrid index to the subgrid origin. */
    awh_ivec  subgridIndex = {0};

    /* Get the subgrid index of the given grid point index. */
    if (*gridPointIndex >= 0)
    {
        gridToSubgridIndex(grid, subgridOrigin, subgridNpoints, *gridPointIndex, subgridIndex);
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
    while (step_in_multidim_array(grid->ndim, subgridNpoints, subgridIndex))
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
static int pointDistanceAlongAxis(const GridAxis *axis, double x, double x0)
{
    int distance = 0;

    if (axis->spacing > 0)
    {
        /* Get the real-valued distance. For a periodic axis, the shortest one. */
        double period = axis->npoints_period*axis->spacing;
        double dx     = get_deviation_periodic(x, x0, period);

        /* Transform the distance into a point distance.
           Shift by +0.5 so we can use floor or integer casting below to get the integer index */
        distance = static_cast<int>(floor(dx/axis->spacing + 0.5));

        /* If periodic, shift the point distance to be in [0, period) */
        make_index_periodic(&distance, axis->npoints_period);
    }

    return distance;
}

/*! \brief
 * Query if a value is in range of the grid.
 *
 * \param[in] value   Value to check.
 * \param[in]  grid   The grid.
 * \returns true if the value is in the grid.
 */
bool value_is_in_grid(const awh_dvec value, const Grid *grid)
{
    /* For each dimension get the one-dimensional index and check if it is in range. */
    for (int d = 0; d < grid->ndim; d++)
    {
        int index = pointDistanceAlongAxis(&grid->axis[d], value[d], grid->axis[d].origin);

        if (!(index >= 0 && index < grid->axis[d].npoints))
        {
            return false;
        }
    }

    return true;
}

/*! \brief
 * Map a value to the closest point index along an axis.
 *
 * \param[in] axis   Grid axis.
 * \param[in] value  Value along the axis.
 * \returns the index closest to the value.
 */
static int get_closest_index_in_grid_along_axis_for_value(const GridAxis *axis, double value)
{
    /* Get the point distance to the origin. This may by an out of index range for the axis. */
    int index = pointDistanceAlongAxis(axis, value, axis->origin);

    if (index < 0 || index >= axis->npoints)
    {
        if (gridaxis_is_periodic(axis))
        {
            GMX_RELEASE_ASSERT(index >= 0 && index < axis->npoints_period,
                               "Index not in periodic interval 0 for AWH periodic axis");
            int endDistance    = (index - (axis->npoints - 1));
            int originDistance = (axis->npoints_period - index);
            index = originDistance < endDistance ? 0 : axis->npoints - 1;
        }
        else
        {
            index = (index < 0) ? 0 : (axis->npoints - 1);
        }
    }

    return index;
}

/*! \brief
 * Map a value to the closest point in the grid.
 *
 * \param[in] value  Value.
 * \param[in] grid   The grid.
 * \returns the point index closest to the value.
 */
int get_closest_index_in_grid(const awh_dvec value, const Grid *grid)
{
    awh_ivec index_multi;

    /* If the index is out of range, modify it so that it is in range by choosing the closest point on the edge. */
    for (int d = 0; d < grid->ndim; d++)
    {
        index_multi[d] = get_closest_index_in_grid_along_axis_for_value(&grid->axis[d], value[d]);
    }

    return multidim_gridindex_to_linear(grid, index_multi);
}

/*! \brief
 * Find and set the neighbors of a grid point.
 *
 * The search space for neighbors is a subgrid with size
 * set by a scope cutoff. In general not all subgrid points
 * will be valid grid points.
 *
 * \param[in] point_index              Grid point index.
 * \param[in] grid                     The grid.
 * \param[in,out] neighborIndexArray   Array to fill with neighbor indices.
 */
static void set_neighbors_of_grid_point(int point_index, const Grid *grid,
                                        std::vector<int> *neighborIndexArray)
{
    awh_ivec       ncandidates_dim      = {0};
    awh_ivec       subgridOrigin        = {0};

    /* Set up a subgrid  array for searching for neighbors in a local region of the given grid point.
       This subgrid generally contain points that are not in the  grid since boundary conditions are not
       taking into account here. */
    for (int d = 0; d < grid->ndim; d++)
    {
        int center_index_d;

        /* The number of candidate points along this dimension is given by the scope cutoff. */
        ncandidates_dim[d] = std::min(c_maxNeighborsAlongAxis, grid->axis[d].npoints);

        /* The origin of the subgrid to search */
        center_index_d       = grid->point[point_index].index[d];
        subgridOrigin[d]     = center_index_d - ncandidates_dim[d]/2;
    }

    /* Find and set the neighbors */
    int  neighbor_index = -1;
    bool aPointExists   = true;

    /* Keep looking for grid points while traversing the subgrid. */
    while (aPointExists)
    {
        /* The point index is updated if a grid point was found. */
        aPointExists = get_next_point_in_local_grid(grid, subgridOrigin, ncandidates_dim, &neighbor_index);

        if (aPointExists)
        {
            neighborIndexArray->push_back(neighbor_index);
        }
    }
}

/*! \brief
 * Allocate and initialize the grid point neighbor arrays.
 *
 * \param[in,out] grid       The grid.
 * returns the number of neighbors allocated for each point.
 */
static size_t init_grid_point_neighbors(Grid *grid)
{
    size_t maxNumNeighbors = 1;

    /* Note: could also generate neighbor list only when needed instead of storing them for each point */
    for (size_t m = 0; m < grid->point.size(); m++)
    {
        set_neighbors_of_grid_point(m, grid, &grid->point[m].neighbor);

        maxNumNeighbors = std::max(maxNumNeighbors,
                                   grid->point[m].neighbor.size());
    }

    return maxNumNeighbors;
}

/*! \brief
 * Allocate and initialize the grid points.
 *
 * \param[in,out] grid   The grid.
 * \param[in] period     The period for each dimension (0 if not periodic).
 */
static void init_grid_points(Grid *grid, const awh_dvec period)
{
    awh_ivec     npoints_dim_work;
    awh_ivec     index_work = {0};

    for (int d = 0; d < grid->ndim; d++)
    {
        /* Temporarily gather the number of points in each dimension in one array */
        npoints_dim_work[d] = grid->axis[d].npoints;
    }

    for (auto &point : grid->point)
    {
        for (int d = 0; d < grid->ndim; d++)
        {
            point.value[d] = grid->axis[d].origin + index_work[d]*grid->axis[d].spacing;

            if (period[d] > 0)
            {
                make_value_periodic(&point.value[d], period[d]);
            }

            point.index[d] = index_work[d];
        }

        step_in_multidim_array(grid->ndim, npoints_dim_work, index_work);
    }

    /* Depending on how the spacing was chosen the input end value might differ slightly from the set point value.
       Make sure they are the same. */
    for (int d = 0; d < grid->ndim; d++)
    {
        grid->axis[d].end = grid->point[grid->point.size() - 1].value[d];
    }
}

/*! \brief
 * Initializes a grid axis.
 *
 * The spacing and number of points are set as a function of a given expected variance
 * of the data that will live on the grid. The number of points scales as 1/sqrt(variance).
 * The inverse variance is given as input so that an infinite number of points is prohibited.
 * Conversely, infinite variance translates to having the minimum number of points (2).
 *
 * \param[in,out] axis           Grid axis to intialize.
 * \param[in] origin             Starting value.
 * \param[in] end                End value.
 * \param[in] inv_pointvariance  Inverse of the expected variance of data on the grid.
 * \param[in] period             Period (0 if not periodic).
 */
static void init_gridaxis(GridAxis *axis, double origin, double end,
                          double inv_pointvariance, double period)
{
    const bool bAxisIsPeriodic = (period > 0);

    axis->origin      = origin;
    axis->end         = end;

    double axisLength = get_interval_length_periodic(axis->origin, axis->end, period);

    /* Automatically determine number of points based on the user given endpoints
       and the expected fluctuations in the umbrella. */
    if (axisLength == 0)
    {
        axis->npoints = 1;
    }
    else if (inv_pointvariance == 0)
    {
        axis->npoints = 2;
    }
    else
    {
        /* Note: nstdevs > 0.  */
        double nstdevs    = axisLength*std::sqrt(inv_pointvariance);

        /* npoints >= 2 since nstdevs > 0. */
        axis->npoints     = 1 + static_cast<int>(std::ceil(c_pointDensity*nstdevs));
    }

    /* Set point spacing based on the number of points */
    if (bAxisIsPeriodic)
    {
        /* Set the grid spacing so that a period is matched exactly by an integer number of points.
           The number of points in a period is equal to the number of grid spacings in a period
           since the endpoints are connected.  */
        axis->npoints_period = axisLength > 0 ? static_cast<int>(std::ceil(period/axisLength*(axis->npoints - 1))) : 1;
        axis->spacing        = period/axis->npoints_period;

        /* Modify the number of grid axis points to be compatible with the period dependent spacing. */
        axis->npoints        = std::min(static_cast<int>(round(axisLength/axis->spacing)) + 1, axis->npoints_period);
    }
    else
    {
        axis->npoints_period = 0;
        axis->spacing        = axis->npoints > 1 ?
            (axis->end - axis->origin)/(axis->npoints - 1) : 0;
    }
}

/* Allocate, initialize and return a grid. */
Grid *init_grid(int ndim, const awh_dvec betak, const awh_dvec origin, const awh_dvec end,
                const awh_dvec period, int *nneighbors_alloc_ptr)
{
    Grid *grid = new Grid;

    grid->ndim = ndim;
    grid->axis.resize(grid->ndim);

    /* Define the discretization along each dimension */
    for (int d = 0; d < grid->ndim; d++)
    {
        init_gridaxis(&grid->axis[d], origin[d], end[d], betak[d], period[d]);
    }

    /* The number of grid points in total */
    int npoints = 1;
    for (int d = 0; d < grid->ndim; d++)
    {
        npoints *= grid->axis[d].npoints;
    }

    grid->point.resize(npoints);

    /* Set their values */
    init_grid_points(grid, period);

    /* Keep a neighbor list for each point */
    *nneighbors_alloc_ptr = init_grid_point_neighbors(grid);

    return grid;
}

/* Maps each point in the grid to a point in the data grid. */
void map_grid_to_datagrid(std::vector<int> *gridpoint_to_datapoint, const double* const *data, int ndata,
                          const char *datafilename, const Grid *grid,
                          const char *correctFormatMessage)
{
    /* Transform the data into a grid in order to map each grid point to a data point
       using the grid functions. */
    Grid dataGrid;
    dataGrid.ndim = grid->ndim;
    dataGrid.point.resize(ndata);
    dataGrid.axis.resize(dataGrid.ndim);

    /* Count the number of points for each dimension. Each dimension
       has its own stride. */
    int stride           = 1;
    int npointsCounted   = 0;
    for (int d = grid->ndim - 1; d >= 0; d--)
    {
        int    npoints_dim   = 0;
        int    pointIndex    = 0;
        double firstValue    = data[d][pointIndex];
        do
        {
            npoints_dim++;
            pointIndex     += stride;
        }
        while ((pointIndex < ndata) && (!gmx_within_tol(firstValue, data[d][pointIndex], GMX_REAL_EPS)));

        /* The stride in dimension dimension d - 1 equals the number of points
           dimension d. */
        stride = npoints_dim;

        npointsCounted = (npointsCounted == 0) ?  npoints_dim : npointsCounted*npoints_dim;

        dataGrid.axis[d].npoints = npoints_dim;
    }

    if (npointsCounted != ndata)
    {
        gmx_fatal(FARGS, "Could not extract data properly from %s. Wrong data format?"
                  "\n\n%s",
                  datafilename, correctFormatMessage);
    }

    /* The data grid has the data that was read and the properties of the AWH grid */
    for (int d = 0; d < dataGrid.ndim; d++)
    {
        /* Note: we don't use init_gridaxis here since that automatically determines the
           number of points and the spacing. Here we already know the grid parameters. */
        int npoints = dataGrid.axis[d].npoints;

        dataGrid.axis[d].origin               = data[d][0];
        dataGrid.axis[d].end                  = data[d][ndata - 1];

        /* The period should be the same for the AWH grid and the data grid. */
        double period     = get_gridaxis_period(grid, d);
        double axisLength = get_interval_length_periodic(dataGrid.axis[d].origin, dataGrid.axis[d].end, period);

        dataGrid.axis[d].spacing              = npoints > 1 ? axisLength/(npoints - 1) : period;
        dataGrid.axis[d].npoints_period       = static_cast<int>(std::round(period/dataGrid.axis[d].spacing));
    }

    /* Map each grid point to a data point. No interpolation, just pick the closest one.
     * It is assumed that the given data is uniformly spaced for each dimension.
     */
    for (size_t m = 0; m < grid->point.size(); m++)
    {
        /* We only define what we need for the datagrid since it's not needed here which is a bit ugly */

        if (!value_is_in_grid(grid->point[m].value, &dataGrid))
        {
            gmx_fatal(FARGS, "%s does not contain data for all coordinate values. "
                      "Make sure your input data covers the whole sampling domain "
                      "and is correctly formatted. \n\n%s",
                      datafilename, correctFormatMessage);
        }
        (*gridpoint_to_datapoint)[m] = get_closest_index_in_grid(grid->point[m].value, &dataGrid);
    }
}
