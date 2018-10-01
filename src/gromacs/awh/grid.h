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
 * \todo: Replace this by a more generic grid class once that is available.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_GRID_H
#define GMX_AWH_GRID_H

#include <memory>
#include <string>

#include "dimparams.h" /* This is needed for awh_dvec */

namespace gmx
{

struct AwhDimParams;

/*! \internal
 * \brief An axis, i.e. dimension, of the grid.
 */
class GridAxis
{
    public:
        /*! \brief Constructor.
         *
         * The spacing and number of points are set such that we have
         * at least the requested point density.
         * Requesting 0 point density results in the minimum number
         * of points (2).
         *
         * \param[in] origin           Starting value.
         * \param[in] end              End value.
         * \param[in] period           Period, pass 0 if not periodic.
         * \param[in] pointDensity     Requested number of point per unit of axis length.
         */
        GridAxis(double origin, double end,
                 double period, double pointDensity);

        /*! \brief Constructor.
         *
         * \param[in] origin           Starting value.
         * \param[in] end              End value.
         * \param[in] period           Period, pass 0 if not periodic.
         * \param[in] numPoints        The number of points.
         */
        GridAxis(double origin, double end,
                 double period, int numPoints);

        /*! \brief Returns if the axis has periodic boundaries.
         */
        bool isPeriodic() const
        {
            return period_ > 0;
        }

        /*! \brief Returns the period of the grid along the axis.
         */
        double period() const
        {
            return period_;
        }

        /*! \brief Returns the grid origin along the axis.
         */
        double origin() const
        {
            return origin_;
        }

        /*! \brief Returns the grid point spacing along the axis.
         */
        double spacing() const
        {
            return spacing_;
        }

        /*! \brief Returns the number of grid points along the axis.
         */
        int numPoints() const
        {
            return numPoints_;
        }

        /*! \brief Returns the period of the grid in points along the axis.
         *
         * Returns 0 if the axis is not periodic.
         */
        int numPointsInPeriod() const
        {
            return numPointsInPeriod_;
        }

        /*! \brief Returns the length of the interval.
         *
         * This returns the distance obtained by connecting the origin point to
         * the end point in the positive direction. Note that this is generally
         * not the shortest distance. For period > 0, both origin and
         * end are expected to take values in the same periodic interval.
         */
        double length() const
        {
            return length_;
        }

        /*! \brief Map a value to the nearest point index along an axis.
         *
         * \param[in] value  Value along the axis.
         * \returns the index nearest to the value.
         */
        int nearestIndex(double value) const;

    private:
        double   origin_;            /**< Interval start value */
        double   length_;            /**< Interval length */
        double   period_;            /**< The period, 0 if not periodic */
        double   spacing_;           /**< Point spacing */
        int      numPoints_;         /**< Number of points in the interval */
        int      numPointsInPeriod_; /**< Number of points in a period (0 if no periodicity) */
};

/*! \internal
 * \brief A point in the grid.
 *
 * A grid point has a coordinate value and a coordinate index of the same dimensionality as the grid.
 * It knows the the linear indices of its neighboring point (which are useful only when handed up to
 * the grid).
 */
struct GridPoint
{
    awh_dvec         coordValue;    /**< Multidimensional coordinate value of this point */
    awh_ivec         index;         /**< Multidimensional point indices */
    std::vector<int> neighbor;      /**< Linear point indices of the neighboring points */
};

/*! \internal
 * \brief The grid, generally multidimensional and periodic.
 *
 * The grid discretizes a multidimensional space with some given resolution.
 * Each dimension is represented by an axis which sets the spatial extent,
 * point spacing and periodicity of the grid in that direction.
 */
class Grid
{
    private:
        /*! \brief Initializes the grid points.
         */
        void initPoints();

    public:
        /*! \brief
         * The point density per sigma of the Gaussian distribution in an umbrella.
         *
         * This value should be at least 1 to uniformly cover the reaction coordinate
         * range with density and having it larger than 1 does not add information.
         */
        static constexpr double c_numPointsPerSigma = 1.0;

        //! Cut-off in sigma for considering points, neglects 4e-8 of the density.
        static constexpr double c_scopeCutoff       = 5.5;

        /*! \brief Construct a grid using AWH input parameters.
         *
         * \param[in] dimParams     Dimension parameters including the expected inverse variance of the coordinate living on the grid (determines the grid spacing).
         * \param[in] awhDimParams  Dimension params from inputrec.
         */
        Grid(const std::vector<DimParams> &dimParams,
             const AwhDimParams           *awhDimParams);

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
        int numDimensions() const
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

        /*! \brief Find the grid point with value nearest to the given value.
         *
         * \param[in] value  Value vector.
         * \returns the grid point index.
         */
        int nearestIndex(const awh_dvec value) const;

        /*! \brief Query if the value is in the grid.
         *
         * \param[in] value  Value vector.
         * \returns true if the value is in the grid.
         * \note It is assumed that any periodicity of value has already been taken care of.
         */
        bool covers(const awh_dvec value) const;

    private:
        std::vector<GridPoint> point_; /**< Points on the grid */
        std::vector<GridAxis>  axis_;  /**< Axes, one for each dimension. */
};

/*! \endcond */

/*! \brief Convert a multidimensional grid point index to a linear one.
 *
 * \param[in] grid        The grid.
 * \param[in] indexMulti  Multidimensional grid point index to convert to a linear one.
 * \returns the linear index.
 */
int multiDimGridIndexToLinear(const Grid &grid, const awh_ivec indexMulti);

/*! \brief Convert multidimensional array index to a linear one.
 *
 * \param[in] indexMulti    Multidimensional index to convert to a linear one.
 * \param[in] numDim        Number of dimensions of the array.
 * \param[in] numPointsDim  Number of points of the array.
 * \returns the linear index.
 * \note This function can be used without having an initialized grid.
 */
int multiDimArrayIndexToLinear(const awh_ivec indexMulti,
                               int            numDim,
                               const awh_ivec numPointsDim);

/*! \brief Convert a linear grid point index to a multidimensional one.
 *
 * \param[in]  grid         The grid.
 * \param[in]  indexLinear  Linear grid point index to convert to a multidimensional one.
 * \param[out] indexMulti   The multidimensional index.
 */
void linearGridindexToMultiDim(const Grid &grid,
                               int         indexLinear,
                               awh_ivec    indexMulti);

/*! \brief Convert a linear array index to a multidimensional one.
 *
 * \param[in]  indexLinear   Linear array index
 * \param[in]  ndim          Number of dimensions of the array.
 * \param[in]  numPointsDim  Number of points for each dimension.
 * \param[out] indexMulti    The multidimensional index.
 */
void linearArrayIndexToMultiDim(int            indexLinear,
                                int            ndim,
                                const awh_ivec numPointsDim,
                                awh_ivec       indexMulti);

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
 * \param[in]     grid            The grid.
 * \param[in]     subgridOrigin   Vector locating the subgrid origin relative to the grid origin.
 * \param[in]     subgridNpoints  Number of points along each subgrid dimension.
 * \param[in,out] gridPointIndex  Pointer to the starting/next grid point index.
 * \returns true if the grid point was updated.
 */
bool advancePointInSubgrid(const Grid     &grid,
                           const awh_ivec  subgridOrigin,
                           const awh_ivec  subgridNpoints,
                           int            *gridPointIndex);

/*! \brief Maps each point in the grid to a point in the data grid.
 *
 * This functions maps an AWH bias grid to a user provided input data grid.
 * The value of data grid point i along dimension d is given by data[d][i].
 * The number of dimensions of the data should equal that of the grid.
 * A fatal error is thrown if extracting the data fails or the data does not cover the whole grid.
 *
 * \param[out] gridpointToDatapoint  Array mapping each grid point to a data point index.
 * \param[in]  data                  2D array in format ndim x ndatapoints with data grid point values.
 * \param[in]  numDataPoints         Number of data points.
 * \param[in]  dataFilename          The data filename.
 * \param[in]  grid                  The grid.
 * \param[in]  correctFormatMessage  String to include in error message if extracting the data fails.
 */
void mapGridToDataGrid(std::vector<int>    *gridpointToDatapoint,
                       const double* const *data,
                       int                  numDataPoints,
                       const std::string   &dataFilename,
                       const Grid          &grid,
                       const std::string   &correctFormatMessage);

/*! \brief
 * Get the deviation along one dimension from the given value to a point in the grid.
 *
 * \param[in] grid        The grid.
 * \param[in] dimIndex    Dimensional index in [0, ndim -1].
 * \param[in] pointIndex  Grid point index.
 * \param[in] value       Value along the given dimension.
 * \returns the deviation of the given value to the given point.
 */
double getDeviationFromPointAlongGridAxis(const Grid &grid,
                                          int         dimIndex,
                                          int         pointIndex,
                                          double      value);

} // namespace gmx

#endif
