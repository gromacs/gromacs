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
 * This file contains datatypes and function declarations necessary
 * for AWH to interface with the grid code.
 *
 * The grid organizes spatial properties of the AWH coordinate points.
 * This includes traversing points in a specific order, locating
 * neighboring points and calculating distances. Multiple dimensions
 * as well as periodic dimensions are supported.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_GRID_H
#define GMX_AWH_GRID_H

#include "types.h" /* This currently needed for awh_dvec */

/*! \cond INTERNAL */

struct awh_dim_params_t;

//! An axis, i.e. dimension, of the grid.
struct GridAxis
{
    double   origin;                /**< Interval start value */
    double   end;                   /**< Interval end value */
    double   spacing;               /**< Point spacing */
    int      npoints;               /**< Number of points in the interval */
    int      npoints_period;        /**< Number of points in a period (0 if no periodicity) */
};

//! A point in the grid.
struct GridPoint
{
    awh_dvec         value;         /**< Multidimensional value of this point */
    awh_ivec         index;         /**< Multidimensional point indices */
    std::vector<int> neighbor;      /**< Linear point indices of the neighboring points */
};

//! The grid, generally multidimensional and periodic.
class Grid
{
    std::vector<GridPoint> point_; /**< Grid points array */
    std::vector<GridAxis>  axis_;  /**< Grid axes array */

    /*! \brief Initializes the grid point array.
     *
     * \param[in] period  The period of the grid in each dimension.
     */
    void initPoints(const awh_dvec period);

    /*! \brief Initializes the neighbors of all grid points.
     */
    void initPointNeighbors();

    public:
        /*! \brief Constructor.
         *
         * \param[in] numDims                 Number of dimensions of the grid.
         * \param[in] dimParams               Dimension parameters including the expected inverse variance of the coordinate living on the grid (determines the grid spacing).
         * \param[in] origin                  Start value for each dimension.
         * \param[in] end                     End value for each dimension.
         * \param[in] period                  Period of each dimension (0 if no periodicity).
         */
        Grid(int numDims, const DimParams *dimParams,
             const awh_dvec origin, const awh_dvec end,
             const awh_dvec period);

        /*! \brief Returns the number of points in the grid.
         *
         * \returns the number of points in the grid.
         */
        size_t numPoints() const
        {
            return point_.size();
        }

        /*! \brief Returns a reference to a point on the grid.
         *
         * \returns a constant reference to a point on the grid.
         */
        const GridPoint &point(size_t pointIndex) const
        {
            return point_[pointIndex];
        }

        /*! \brief Returns the dimensionality of the grid.
         *
         * \returns the dimensionality of the grid.
         */
        int ndim() const
        {
            return axis_.size();
        }

        /*! \brief Returns the grid axes.
         *
         * \returns a constant reference to the grid axes.
         */
        const std::vector<GridAxis> &axis() const
        {
            return axis_;
        }

        /*! \brief Returns a grid axis.
         *
         * param[in] dim  Dimension to return the grid axis for.
         * \returns a constant reference to the grid axis.
         */
        const GridAxis &axis(int dim) const
        {
            return axis_[dim];
        }
};

/*! \endcond */

/*! \brief Convert a multidimensional grid point index to a linear one.
 *
 * \param[in] grid        The grid.
 * \param[in] index_multi Multidimensional grid point index to convert to a linear one.
 * \returns the linear index.
 */
int multidim_gridindex_to_linear(const Grid *grid, const awh_ivec index_multi);

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
void linear_gridindex_to_multidim(const Grid *grid, int index_linear, awh_ivec index_multi);

/*! \brief Convert a linear array index to a multidimensional one.
 *
 * \param[in] index_linear        Linear array index
 * \param[in] ndim                Number of dimensions of the array.
 * \param[in] npoints_dim         Number of points for each dimension.
 * \param[out] index_multi        The multidimensional index.
 */
void linear_array_index_to_multidim(int index_linear, int ndim, const awh_ivec npoints_dim, awh_ivec index_multi);

/*! \brief Find the grid point with value closest to the given value.
 *
 * \param[in] value               Value vector.
 * \param[in] grid                The grid.
 * \returns the grid point index.
 */
int get_closest_index_in_grid(const awh_dvec value, const Grid *grid);


/*! \brief Query if the value is in the grid.
 *
 * \param[in] value               Value vector.
 * \param[in] grid                The grid.
 * \returns true if the value is in the grid.
 * \note It is assumed that any periodicity of value has already been taken care of.
 */
bool value_is_in_grid(const awh_dvec value, const Grid *grid);

/*! \brief
 * Find the next grid point in the sub-part of the grid given a starting point.
 *
 * The given grid point index is updated to the next valid grid point index
 * by traversing the sub-part of the grid, here termed the subgrid.
 * Since the subgrid range might extend beyond the actual size of the grid,
 * the subgrid is traversed until a point both in the subgrid and grid is
 * found. If no point is found, the function returns false and the index is
 * not modified. The starting point needs to be inside of the subgrid. However,
 * if this index is not given, meaning < 0, then the search is initialized at
 * the subgrid origin, i.e. in this case the "next" grid point index is
 * defined to be the first common grid/subgrid point.
 *
 * \param[in] grid                The grid.
 * \param[in] subgridOrigin       Vector locating the subgrid origin relative to the grid origin.
 * \param[in] subgridNpoints      Number of points along each subgrid dimension.
 * \param[in,out] gridPointIndex  Pointer to the starting/next grid point index.
 * \returns true if the grid point was updated.
 */
bool get_next_point_in_local_grid(const Grid *grid,
                                  const awh_ivec subgridOrigin, const awh_ivec subgridNpoints,
                                  int *gridPointIndex);

/*! \brief Gets the length of the given grid dimension.
 *
 * \param[in] axis  The grid axis.
 * \returns the length.
 */
double get_gridaxis_length(const GridAxis &axis);

/*! \brief Maps each point in the grid to a point in the data grid.
 *
 * This functions maps an AWH bias grid to a user provided input data grid.
 * The value of data grid point i along dimension d is given by data[d][i].
 * The number of dimensions of the data should equal that of the grid.
 * A fatal error is thrown if extracting the data fails or the data does not cover the whole grid.
 *
 * \param[out] gridpoint_to_datapoint     Array mapping each grid point to a data point index.
 * \param[in] data                        2D array in format ndim x ndatapoints with data grid point values.
 * \param[in] ndatapoints                 Number of data points.
 * \param[in] datafilename                The data filename.
 * \param[in] grid                        The grid.
 * \param[in] correctFormatMessage        String to include in error message if extracting the data fails.
 */
void map_grid_to_datagrid(std::vector<int> *gridpoint_to_datapoint, const double* const *data, int ndatapoints,
                          const char *datafilename, const Grid *grid, const char *correctFormatMessage);

/*! \brief
 * Get the deviation along one dimension from the given value to a point in the grid.
 *
 * \param[in] grid        The grid.
 * \param[in] dimindex    Dimensional index in [0, ndim -1].
 * \param[in] pointindex  Grid point index.
 * \param[in] value       Value along the given dimension.
 * \returns the deviation of the given value to the given point.
 */
double get_deviation_from_point_along_gridaxis(const Grid *grid, int dimindex, int pointindex, double value);

/*! \endcond */

#endif
