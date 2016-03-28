/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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

/* Constants for defining and moving in the grid */
/* The point density per sigma of the Gaussian distribution in an Umbrella.
 * This value should be at least 1 to uniformly cover the reaction coordinate
 * range with density and having it larger than 1 does not add information.
 */
static const double c_pointDensity          = 1.0;
/* Cut-off in sigma for considering points, neglects 0.6e-6 of the density */
static const double c_scopeCutoff           = 5.0;
/* The maximum number of neighboring points to consider along an axis */
static const int    c_maxNeighborsAlongAxis = static_cast<int>(c_pointDensity*c_scopeCutoff*2 + 1);


/* Modify x so that it is periodic in [-period/2, +period/2).
   This function assumes that x is at most one period away from this interval. */
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

/* Get the length of the interval (origin, end) where origin, end should take values in (-period/2, +period/2)
 * Note that this does not return the minimum distance between the interval endpoints.
 * For period = 0, end - origin is returned. */
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
    return L;
}

/* Get the smallest deviation of x - x0 in the periodic interval P=(-period/2, period/2).
 * It is assumed that |x - x0| < period. For period = 0, x - x0 is returned.
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

static bool gridaxis_is_periodic(const grid_axis_t *axis)
{
    return axis->npoints_period  > 0;
}

static double get_gridaxis_period(const grid_t *grid, int dimindex)
{
    return grid->axis[dimindex].npoints_period*grid->axis[dimindex].spacing;
}

double get_deviation_from_point_along_gridaxis(const grid_t *grid, int dimindex, int pointindex, double value)
{
    double pointvalue = grid->point[pointindex].value[dimindex];

    return get_deviation_periodic(value, pointvalue, get_gridaxis_period(grid, dimindex));
}

double get_gridaxis_length(const grid_t *grid, int dimindex)
{
    return get_interval_length_periodic(grid->axis[dimindex].origin,
                                        grid->axis[dimindex].end,
                                        get_gridaxis_period(grid, dimindex));
}

/* Convert a linear index to a multidimensional index where:
   index_multi = (i1, i2,.., indim) and index_linear = i1*N2*N3*..*Nndim + i2*N3*..*Nndim + ... + indim */
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

void linear_gridindex_to_multidim(const grid_t *grid, int index_linear, awh_ivec index_multi)
{
    awh_ivec npoints_dim;

    for (int d = 0; d < grid->ndim; d++)
    {
        npoints_dim[d] = grid->axis[d].npoints;
    }

    linear_array_index_to_multidim(index_linear, grid->ndim, npoints_dim, index_multi);
}

/* Convert a multidim array index to a linear index for the given grid */
int multidim_array_index_to_linear(const awh_ivec index_multi, int ndim, const awh_ivec npoints_dim)
{
    int stride       = 1;
    int index_linear = 0;
    for (int d = ndim - 1; d >= 0; d--)
    {
        index_linear += stride*index_multi[d];
        stride       *= npoints_dim[d];
    }

    return index_linear;
}

int multidim_gridindex_to_linear(const grid_t *grid, const awh_ivec index_multi)
{
    awh_ivec npoints_dim;

    for (int d = 0; d < grid->ndim; d++)
    {
        npoints_dim[d] = grid->axis[d].npoints;
    }

    return multidim_array_index_to_linear(index_multi, grid->ndim, npoints_dim);
}

/* Take a step in a multidimensional array and increase the multidimensional index in the highest dimension possible. */
static bool step_in_multidim_array(int ndim, const awh_ivec npoints_dim, awh_ivec index_dim)
{
    bool bStepped = false;

    for (int d = ndim - 1; d >= 0 && !bStepped; d--)
    {
        /* Workaround for bug in CLANG 3.4, should be removed with 3.8 */
#ifndef __clang_analyzer__
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
#endif
    }

    return bStepped;
}

/*! \brief Transform the global linear index to a multidimensional index of the given local array.
 *
 * \note
 * This function is when traversing a subgrid of the grid (get_next_point_in_local_grid) to transform
 * the initial global grid index to to the initial subgrid index.
 */
static void global_to_local_index(const grid_t *grid, const awh_ivec origin_local_dim,
                                  int global_linear_index, awh_ivec local_index_dim)
{
    for (int d = 0; d < grid->ndim; d++)
    {
        /* The grid point index relative to the local origin index */
        local_index_dim[d] = grid->point[global_linear_index].index[d] - origin_local_dim[d];

        /* The local, relative index may not be in the local array.
           If the grid has periodic boundaries, the local array and index could wrap around.
           We unwrap the local index here. */
        make_index_periodic(&local_index_dim[d], grid->axis[d].npoints_period);
    }
}

/*! \brief
 * Transform the multidimensional index of the given local array to the global linear grid index.
 *
 * \note
 * This function is when traversing a subgrid of the grid (get_next_point_in_local_grid) to transform
 *  the updated subgrid index back to the global grid indexing.
 */
static bool local_to_global_index(const grid_t   *grid,
                                  const awh_ivec  origin_local_dim,
                                  const awh_ivec  local_index_dim,
                                  int            *global_linear_index_ptr)
{
    awh_ivec global_index_dim;

    /* Check and apply boundary conditions for each dimension */
    for (int d = 0; d < grid->ndim; d++)
    {
        /* Transform to global multidimensional indexing by adding the origin */
        global_index_dim[d] = origin_local_dim[d] + local_index_dim[d];

        /* The local grid is allowed to stick out on the edges of the global grid. Here the boundary conditions are applied.*/
        if (global_index_dim[d] < 0 ||  global_index_dim[d] > grid->axis[d].npoints - 1)
        {
            /* Try to wrap around if periodic. */
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

    awh_ivec npoints_dim;

    /* Translate from multidimensional to linear indexing and set the return value */
    for (int d = 0; d < grid->ndim; d++)
    {
        npoints_dim[d] = grid->axis[d].npoints;
    }

    (*global_linear_index_ptr) = multidim_array_index_to_linear(global_index_dim, grid->ndim, npoints_dim);

    return true;
}

/*! \brief Find the next point in the local subgrid.
 * \note This function is used by awh.cpp to generate the update list, a subset the total number of points. It's also used
 * for initializing the neighbors of each grid point.
 */
bool get_next_point_in_local_grid(const grid_t *grid, const awh_ivec origin_local_dim, const awh_ivec npoints_local_dim,
                                  int *global_linear_index_ptr)
{
    awh_ivec  local_index_dim = {0};

    if (*global_linear_index_ptr >= 0)
    {
        /* Set the local index corresponding to the given global index.
           Note: the local index is unwrapped any periodic boundaries, i.e. it is not "aware"
           of boundaries of the global grid which simplifies stepping in the local array (below).
           The boundary conditions of the global grid are taken care of when transforming back.
         */
        global_to_local_index(grid, origin_local_dim, *global_linear_index_ptr, local_index_dim);
        for (int d = 0; d < grid->ndim; d++)
        {
            GMX_RELEASE_ASSERT((local_index_dim[d] >= 0) && (local_index_dim[d] < npoints_local_dim[d]),
                               "Attempted to convert the AWH grid point index not in local array to out of bounds local array index");
        }

    }
    else
    {
        /* If no index is given we start at the local origin (which local_index_dim is initialized to).
           If this is a valid point in the global grid we're done, otherwise keep looking below. */
        if (local_to_global_index(grid, origin_local_dim, local_index_dim, global_linear_index_ptr))
        {
            return true;
        }
    }

    while (step_in_multidim_array(grid->ndim, npoints_local_dim, local_index_dim))
    {
        /* If a point is found, the corresponding global index is set */
        if (local_to_global_index(grid, origin_local_dim, local_index_dim, global_linear_index_ptr))
        {
            return true;
        }
    }

    return false;
}

static int get_index_along_axis_for_value(const grid_axis_t *axis, double value)
{
    int index = 0;

    if (axis->spacing > 0)
    {
        /* Shift by +0.5 so we can use floor or integer casting below to get the integer index */
        index = static_cast<int>(floor((value - axis->origin)/axis->spacing + 0.5));
    }

    return index;
}

bool value_is_in_grid(const awh_dvec value, const grid_t *grid)
{
    /* For each dimension get the one-dimensional index and check if it is in range. */
    for (int d = 0; d < grid->ndim; d++)
    {
        int index = get_index_along_axis_for_value(&grid->axis[d], value[d]);
        if (!(index >= 0 && index < grid->axis[d].npoints))
        {
            return false;
        }
    }

    return true;
}

static int get_closest_index_in_grid_along_axis_for_value(const grid_axis_t *axis, double value)
{
    int index = get_index_along_axis_for_value(axis, value);

    if (index < 0 || index >= axis->npoints)
    {
        if (gridaxis_is_periodic(axis))
        {
            int index_shifted, index_distance_end, index_distance_origin;

            /* Since we are out of range, either the index is negative or >= the number of points along this axis.
             * If the index < 0, add a period to get the the get the shortest positive distance to the interval end point.
             * The distance to the other origin is obtained by knowing that the two distances + the interval length = 1 period.
             * Could probably do this cleaner, e.g. using modulo.
             */
            index_shifted = index < 0 ? index + axis->npoints_period : index;

            index_distance_end    = index_shifted - (axis->npoints - 1);
            index_distance_origin = axis->npoints_period - (axis->npoints - 1) - index_distance_end;
            index                 = index_distance_origin < index_distance_end ? 0 : axis->npoints - 1;
        }
        else
        {
            index = index < 0 ? 0 : axis->npoints - 1;
        }
    }

    return index;
}

int get_closest_index_in_grid(const awh_dvec value, const grid_t *grid)
{
    awh_ivec index_multi;

    /* If the index is out of range, modify it so that it is in range by choosing the closest point on the edge. */
    for (int d = 0; d < grid->ndim; d++)
    {
        index_multi[d] = get_closest_index_in_grid_along_axis_for_value(&grid->axis[d], value[d]);
    }

    return multidim_gridindex_to_linear(grid, index_multi);
}

static int set_neighbors_of_grid_point(int point_index, grid_t *grid, int *neighbor_index_array)
{
    awh_ivec       ncandidates_dim      = {0};
    awh_ivec       origin_local_dim     = {0};

    /* Set up an array for searching for neighbors in a local region of the given grid point. This local array may contain points
       that are not in the global grid since boundary conditions are treated in get_next_point_in_local_grid. */
    for (int d = 0; d < grid->ndim; d++)
    {
        int center_index_d;

        /* The number of candidate points along this dimension */
        ncandidates_dim[d] = std::min(c_maxNeighborsAlongAxis, grid->axis[d].npoints);

        /* The origin of the local array to search */
        center_index_d       = grid->point[point_index].index[d];
        origin_local_dim[d]  = center_index_d - ncandidates_dim[d]/2;
    }

    /* Find and store the neighbors */
    int  nneighbors     = 0;
    int  neighbor_index = -1;
    bool aPointExists   = true;
    while (aPointExists)
    {
        aPointExists = get_next_point_in_local_grid(grid, origin_local_dim, ncandidates_dim, &neighbor_index);

        if (aPointExists)
        {
            neighbor_index_array[nneighbors] = neighbor_index;
            nneighbors++;
        }
    }

    return nneighbors;
}

static int init_grid_point_neighbors(grid_t *grid)
{
    /* Allocate assuming the maximum number of neighbors for all points */
    int nneighbors_max = 1;
    for (int d = 0; d < grid->ndim; d++)
    {
        /* At most all points along a dimension are neighbors */
        nneighbors_max *= c_maxNeighborsAlongAxis;
    }

    for (int m = 0; m < grid->npoints; m++)
    {
        snew(grid->point[m].neighbor, nneighbors_max);
    }

    /* Note: could also generate neighbor list only when needed instead of storing them for each point */
    for (int m = 0; m < grid->npoints; m++)
    {
        grid->point[m].nneighbors = set_neighbors_of_grid_point(m, grid, grid->point[m].neighbor);
    }

    return nneighbors_max;
}

static void init_grid_points(grid_t *grid, const awh_dvec period)
{
    awh_ivec     npoints_dim_work;
    awh_ivec     index_work = {0};

    for (int d = 0; d < grid->ndim; d++)
    {
        /* Temporarily gather the number of points in each dimension in one array */
        npoints_dim_work[d] = grid->axis[d].npoints;
    }

    for (int m = 0; m < grid->npoints; m++)
    {
        for (int d = 0; d < grid->ndim; d++)
        {
            grid->point[m].value[d] = grid->axis[d].origin + index_work[d]*grid->axis[d].spacing;

            if (period[d] > 0)
            {
                make_value_periodic(&grid->point[m].value[d], period[d]);
            }

            grid->point[m].index[d] = index_work[d];
        }

        step_in_multidim_array(grid->ndim, npoints_dim_work, index_work);
    }

    /* Depending on how the spacing was chosen the input end value might differ slightly from the set point value.
       Make sure they are the same. */
    for (int d = 0; d < grid->ndim; d++)
    {
        grid->axis[d].end = grid->point[grid->npoints - 1].value[d];
    }
}

static void init_gridaxis(grid_axis_t *axis, double origin, double end,
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

grid_t *init_grid(int ndim, const awh_dvec betak, const awh_dvec origin, const awh_dvec end,
                  const awh_dvec period, int *nneighbors_alloc_ptr)
{
    grid_t     *grid;

    snew(grid, 1);
    grid->ndim = ndim;
    snew(grid->axis, grid->ndim);

    /* Define the discretization along each dimension */
    for (int d = 0; d < grid->ndim; d++)
    {
        init_gridaxis(&grid->axis[d], origin[d], end[d], betak[d], period[d]);
    }

    /* The number of grid points in total */
    grid->npoints = 1;
    for (int d = 0; d < grid->ndim; d++)
    {
        grid->npoints *= grid->axis[d].npoints;
    }

    snew(grid->point, grid->npoints);

    /* Set their values */
    init_grid_points(grid, period);

    /* Keep a neighbor list for each point */
    *nneighbors_alloc_ptr = init_grid_point_neighbors(grid);

    return grid;
}

/* The value of data grid point i along dimension d is given by data[d][i].
   The number of dimensions of the data should equal that of the grid. */
void map_grid_to_datagrid(int *gridpoint_to_datapoint, const double* const *data, int ndata,
                          const char *datafilename, const grid_t *grid)
{
    int           itmp;
    double        stride;
    grid_t       *datagrid;

    /* Transform the data into a grid in order to map each grid point to a data point
       using the grid functions. */
    snew(datagrid, 1);
    datagrid->ndim    = grid->ndim;
    datagrid->npoints = ndata;
    snew(datagrid->axis, datagrid->ndim);

    /* Count the number of points for each dimension.
     * Each dimension has different stride.
     */
    stride = 1;
    itmp   = 1;
    for (int d = grid->ndim-1; d >= 0; d--)
    {
        int ndata_dim = 0;
        int n         = 0, n_prev;

        do
        {
            ndata_dim++;
            n_prev = n;
            n     += stride;
        }
        while ((n < ndata) && (data[d][n] - data[d][n_prev] > 0));

        stride = ndata_dim;

        itmp *= ndata_dim;

        datagrid->axis[d].npoints      =  ndata_dim;
    }

    if (itmp != ndata)
    {
        gmx_fatal(FARGS, "In file %s, could not extract data properly. Wrong data format?",
                  datafilename);
    }

    /* The data grid has the data that was read and the properties of the AWH grid */
    for (int d = 0; d < datagrid->ndim; d++)
    {
        datagrid->axis[d].origin               = data[d][0];
        datagrid->axis[d].end                  = data[d][datagrid->axis[d].npoints - 1];
        datagrid->axis[d].spacing              = (data[d][ndata - 1] - data[d][0])/(datagrid->axis[d].npoints - 1);
        datagrid->axis[d].npoints_period       = grid->axis[d].npoints_period;
        datagrid->axis[d].npoints_period       = static_cast<int>(round(get_gridaxis_period(grid, d)/datagrid->axis[d].spacing));
    }

    /* Map each grid point to a data point. No interpolation, just pick the closest one.
     * It is assumed that the given data is uniformly spaced for each dimension.
     */
    for (int m = 0; m < grid->npoints; m++)
    {
        /* We only define what we need for the datagrid since it's not needed here which is a bit ugly */

        if (!value_is_in_grid(grid->point[m].value, datagrid))
        {
            gmx_fatal(FARGS, "File does not contain data %s for all coordinate grid points. "
                      "Make sure your input data covers the whole AWH domain.",
                      datafilename);
        }
        gridpoint_to_datapoint[m] = get_closest_index_in_grid(grid->point[m].value, datagrid);
    }

    sfree(datagrid->axis);
    sfree(datagrid);
}
