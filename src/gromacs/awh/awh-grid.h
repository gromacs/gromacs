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

/*! \libinternal \file
 *
 *
 * \brief
 * This file contains datatypes and function declarations necessary
   for AWH to interface with the AWH grid code.
 *
 * \author Viveca Lindahl
 * \inlibraryapi
 */

#ifndef GMX_AWH_GRID_H
#define GMX_AWH_GRID_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/awh/awh-types.h" /* This currently needed for awh_dvec */

struct awhdim_params_t;
struct t_awh_grid;

//! Grid data structures
typedef struct t_awh_grid_axis {
    double   origin;                   /* Interval start value */
    double   end;                      /* Interval end value */
    double   spacing;                  /* Point spacing */
    int      npoints;                  /* Number of points in the interval */
    int      npoints_period;           /* Number of points in a period (0 if no periodicity) */
} t_awh_grid_axis;

typedef struct t_awh_grid_point {
    double *value;                      /* Multidimensional value of this point */
    int    *index;                      /* Multidimensional point indices */
    int     nneighbors;                 /* Number of neighboring points (including this point) */
    int    *neighbor;                   /* Linear point indices of the neighboring points */
} t_awh_grid_point;

typedef struct t_awh_grid {
    int               npoints;          /* Number of points in the grid */
    t_awh_grid_point *point;            /* Grid points array of length npoints */
    int               ndim;             /* Number of dimensions of the grid */
    t_awh_grid_axis  *axis;             /* Grid axes array of length ndim */
} t_awh_grid;


/*! \brief Allocate, initialize and return a grid.
 *
 * \param[in] ndim        Number of dimensions of the grid.
 * \param[in] betak       Expected inverse variance of the coordinate living on the grid (determines the grid spacing).
 * \param[in] origin      Start value for each dimension.
 * \param[in] origin      End value for each dimension.
 * \param[in] period      Period of each dimension (0 if no periodicity).
 * \param[out] nneighbors_alloc      The number of neighbors allocated memory for for each gridpoint.
 * \returns the initialized grid.
 */
t_awh_grid *init_grid(int ndim, const awh_dvec betak, const awh_dvec origin, const awh_dvec end,
                      const awh_dvec period, int *nneighbors_alloc_ptr);

/*! \brief Convert a multidimensional grid point index to a linear one.
 *
 * \param[in] grid        The grid.
 * \param[in] index_multi Multidimensional grid point index to convert to a linear one.
 * \returns the linear index.
 */
int multidim_gridindex_to_linear(const t_awh_grid *grid, const awh_ivec index_multi);

/*! \brief Convert multidimensional array index to a linear one.
 *
 * \param[in] index_multi Multidimensional index to convert to a linear one.
 * \param[in] ndim        Number of dimensions of the array.
 * \param[in] npoints_dim Number of points of the array.
 * \returns the linear index.
 * \note This function can be used without having an initialized grid.
 */
int multidim_array_index_to_linear(const awh_ivec index_multi,  int ndim, const awh_ivec npoints_dim);

/*! \brief Convert a linear grid point index to a multidimensional one.
 *
 * \param[in] grid                The grid.
 * \param[in] index_linear        Linear grid point index to convert to a multidimensional one.
 * \param[out] index_multi        The multidimensional index.
 */
void linear_gridindex_to_multidim(const t_awh_grid *grid, int index_linear, awh_ivec index_multi);

/*! \brief Convert a linear array index to a multidimensional one.
 *
 * \param[in] index_linear        Linear array index
 * \param[in] ndim                Number of dimensions of the array.
 * \param[out] index_multi        The multidimensional index.
 */
void linear_array_index_to_multidim(int index_linear, int ndim, const awh_ivec npoints_dim, awh_ivec index_multi);

/*! \brief Find the grid point with value closest to the given value.
 *
 * \param[in] value               Value vector.
 * \param[in] grid                The grid.
 * \returns the grid point index.
 */
int get_closest_index_in_grid(const awh_dvec value, const t_awh_grid *grid);


/*! \brief Query if the value is in the grid.
 *
 * \param[in] value               Value vector.
 * \param[in] grid                The grid.
 * \returns true if the value is in the grid.
 * \note It is assumed that any periodicity of value has already been taken care of.
 */
bool value_is_in_grid(const awh_dvec value, const t_awh_grid *grid);

/*! \brief Finds the next grid point contained in the given local subgrid.
 *
 * The given initial point index is increased until a point index within the
 * local subgrid is found.
 *
 * \param[in] grid                        The grid.
 * \param[in] origin_local_dim           Origin vector of the local grid.
 * \param[in] npoints_local_dim          Number of points of the local grid.
 * \param[in,out] global_linear_index_ptr Pointer to the initial/next local point index.
 * \returns true if the next point was found, false if the initial point was the end point.
 */
bool get_next_point_in_local_grid(const t_awh_grid *grid, const awh_ivec origin_local_dim,
                                  const awh_ivec npoints_local_dim, int *global_linear_index_ptr);

/*! \brief Gets the length of the grid along the given dimension.
 *
 * \param[in] grid                        The grid.
 * \param[in] dimindex                    Dimension index in [0, ndim - 1].
 * \returns the length.
 */
double get_gridaxis_length(const t_awh_grid *grid, int dimindex);

/*! \brief Maps each point in the grid to a point in the data grid.
 *
 * A fatal error is thrown if extracting the data fails or the data does not cover the whole grid.
 *
 * \param[out] gridpoint_to_datapoint     Array mapping each grid point to a data point index.
 * \param[in] datavalues                  2D array in format ndim x ndatapoints with data grid point values.
 * \param[in] ndatapoints                 Number of data points.
 * \param[in] datafilename                The data filename.
 * \param[in] grid                        The grid.
 */
void map_grid_to_datagrid(int *gridpoint_to_datapoint, const double* const *datagrid, int ndatapoints,
                          const char *datafilename, const t_awh_grid *grid);

/*! \brief Gets deviation of a value from a point, along a single dimension.
 * \param[in] grid                        The grid.
 * \param[in] grid                        Dimensional index in [0, ndim - 1].
 * \param[in] pointindex                  Grid point index.
 * \param[in] value                       Coordinate value along the given dimension.
 */
double get_deviation_from_point_along_gridaxis(const t_awh_grid *grid, int dimindex, int pointindex, double value);

#endif
