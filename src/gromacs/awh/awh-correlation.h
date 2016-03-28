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
 * \brief
 * Declares data types and function for AWH to do force correlation
 * statistics.
 *
 * \author Viveca Lindahl
 * \inlibraryapi
 */

#ifndef GMX_AWH_CORRELATION_H
#define GMX_AWH_CORRELATION_H

#include "gromacs/utility/basedefinitions.h"

//! Correlation block averaging data.
typedef struct t_correlation_blockdata {
    double correlation_integral;       /**< The time integral of the correlation function of x and y, corr(x(0), y(t)). */
    double blocklength;                /**< The length of each block used for block averaging. */
    double sum_wx;                     /**< Weighted sum of x for current block. */
    double sum_wy;                     /**< Weighted sum of y for current block. */
    double sum_w;                      /**< Sum weights for current block. */
    int    blockindex_prev;            /**< The last block index data was added to (needed only for block length in terms of time). */
} t_correlation_blockdata;

//! Correlation matrix element.
typedef struct t_correlation {
    int                       xy_index[2];         /**< The dimensional indices of the vector components x and y. */
    double                    covariance;          /**< The covariance, covar(x, y) = corr(x(t = 0), y(t = 0)). */
    int                       nblockdata;          /**< To be able to increase the block length later on, data is saved for several block lengths. */
    t_correlation_blockdata  *blockdata;           /**< The block averaging data */
} t_correlation;

//! Correlation matrix.
typedef struct correlation_matrix_t {
    int            ncorr;                /**< The number of stored correlation matrix elements. */
    t_correlation *corr;                 /**< Array with the correlation elements corr(x, y) in the matrix, where x, y are vector components. */
} correlation_matrix_t;

//! Grid of local correlation matrices.
typedef struct correlation_grid_t {
    double                  dt_sample;                     /**< Time in between samples. */
    bool                    bBlocklength_in_weight;        /**< If true, a block is measured in probability weight, otherwise in time. */
    int                     nstsample;                     /**< Number of steps per sample. */
    int                     ncorrmatrix;                   /**< Number correlation matrices in the grid (equal to the number of points). */
    correlation_matrix_t   *corrmatrix;                    /**< Correlation matrices */
} correlation_grid_t;

/*! \brief Allocate, initialize and return a correlation grid struct.
 *
 * \param[in] npoints                     Number of points in the grid.
 * \param[in] ndim                        Number of dimensions of the grid.
 * \param[in] nblocks                     Number of blocks to use for the block averaging.
 * \param[in] blocklength                 Length of the blocks used for block averaging.
 * \param[in] bBlocklength_in_weight      If true, a block is measured in probability weight, otherwise in time.
 * \param[in] dt_step                     Time step.
 * \param[in] nstsample                   Number of steps per sample.
 * \returns a pointer to the initialized correlation grid.
 */
correlation_grid_t *init_correlation_grid(int npoints, int ndim,
                                          int nblocks, double blocklength,
                                          bool bBlocklength_in_weight,
                                          double dt_step,
                                          int nstsample);

/*! \brief Returns if correlation should be sampled at the given step.
 *
 * \param[in] corrgrid      The correlation grid.
 * \param[in] step          Time step.
 * \returns true if the step is a correlation sampling step
 */
bool time_to_sample_correlation(const correlation_grid_t *corrgrid, gmx_int64_t step);

/*! \brief Adds a weighted data vector to one point in the correlation grid.
 *
 * \param[in] corrgrid      The correlation grid.
 * \param[in] pointindex    Index of the point to add data to.
 * \param[in] weight        Weight to assign to the data.
 * \param[in] dimdata       Data vector of the same dimensionality as the grid.
 * \param[in] time          The time when the data was sampled.
 */
void add_data_to_correlation_matrix(correlation_grid_t *corrgrid, int pointindex,
                                    double weight, const double *dimdata, double t);

/*! \brief Returns the number of correlation elements in each correlation matrix in the grid.
 *
 * \param[in] corrgrid      The correlation grid.
 * \returns the number of correlation elements.
 */
int get_ncorrelation(const correlation_grid_t *corrgrid);

/*! \brief Returns an element of the time integrated correlation matrix at a given point in the grid.
 *
 * The units of the integral are time*(units of data)^2. This will be friction units time/length^2
 * if the data unit is 1/length.
 *
 * The correlation index lists the elements of the upper-triangular correlation matrix row-wise.
 * (TODO: this should ideally not have to be known by the caller.)
 *
 * \param[in] corrgrid      The correlation grid.
 * \param[in] pointindex    Point index in the grid.
 * \param[in] corrindex     Correlation element index.
 * \returns the integral.
 */
double get_correlation_timeintegral(const correlation_grid_t *corrgrid, int pointindex, int corrindex);

/*! \brief Returns an element of the correlation time matrix at a given point in the correlation grid.
 *
 *  The correlation time matrix is the time-integrated correlation function matrix in units of time.
 *
 * \param[in] corrgrid      The correlation grid.
 * \param[in] pointindex    Point index in the grid.
 * \param[in] corrindex     Correlation element index.
 * \returns the correlation time.
 */
double get_correlation_time(const correlation_grid_t *corrgrid, int pointindex, int corrindex);

/*! \brief
 * Returns the volume element of the correlation metric.
 *
 * The matrix of the metric equals the time-integrated correlation matrix. The volume element of
 * the metric therefore equals the square-root of the absolute value of its determinant according
 * to the standard formula for a volume element in a metric space.
 *
 * Since the units of the metric matrix elements are time*<units of data>^2, the volume element has
 * units of (sqrt(time)*<units of data>)^<ndim of data>
 *
 * \param[in] corrgrid      The correlation grid.
 * \param[in] pointindex    Point index in the grid.
 * \returns the volume element.
 */
double get_correlation_volelem(const correlation_grid_t *corrgrid, int pointindex);

/* Right now the below functions are only used for an initial log printing. */

//! Get the current blocklength.
double get_blocklength(const correlation_grid_t *corrgrid);

//! Get the current number of blocks.
int get_nblocks(const correlation_grid_t *corrgrid);

#endif
