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
 * Contains datatypes and function declarations needed by AWH to
 * have its force correlation data checkpointed.
 *
 * \author Viveca Lindahl
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_AWH_CORRELATION_HISTORY_H
#define GMX_MDTYPES_AWH_CORRELATION_HISTORY_H

struct correlation_blockdata_history_t;
struct correlation_grid_t;
struct t_commrec;

/*! \cond INTERNAL */

//! Correlation block averaging data.
struct correlation_blockdata_history_t
{
    double correlation_integral;       /**< The time integral of the correlation function of x and y, corr(x(0), y(t)). */
    double blocklength;                /**< The length of each block used for block averaging. */
    double sum_wx;                     /**< Weighted sum of x for current block. */
    double sum_wy;                     /**< Weighted sum of y for current block. */
    double sum_w;                      /**< Sum weights for current block. */
    int    blockindex_prev;            /**< The last block index data was added to (needed only for block length in terms of time). */
};

//! Correlation matrix element
struct correlation_history_t
{
    double                           covariance;   /**< The covariance, covar(x, y) = corr(x(t = 0), y(t = 0)). */
    correlation_blockdata_history_t *blockdata;    /**< The block averaging data */
};

/* The below structs do not contain any history dependence and are just for organizing the data in the structs above. */

//! Correlation matrix.
struct correlation_matrix_history_t
{
    correlation_history_t *corr;                  /**< Array with the correlation elements corr(x, y) in the matrix, where x, y are vector components. */
};

//! Grid of local correlation matrices.
struct correlation_grid_history_t
{
    /* These counts here since we curently need them for initializing the correlation grid when reading a checkpoint */
    int                           ncorrmatrix;        /**< Number correlation matrices in the grid (equal to the number of points). */
    int                           ncorr;              /**< The number of stored correlation matrix elements. */
    int                           nblockdata;         /**< To be able to increase the block length later on, data is saved for several block lengths. */
    correlation_matrix_history_t *corrmatrix;         /**< Correlation matrices. */
};

/*! \endcond */

/*! \brief Allocate a correlation grid history struct.
 *
 * This would be called if there is not an correlation grid initialized but
 * basic parameters are available.
 *
 * \param[in,out] corrgrid_hist  Correlation grid history to initialize.
 * \param[int]    ncorrmatrix    Number of correlation matrices in the grid.
 * \param[int]    ncorr          Number of correlation elements in each matrix.
 * \param[int]    ncorr          Number of block data averages.
 */
void init_correlation_grid_history(correlation_grid_history_t *corrgrid_hist,
                                   int ncorrmatrix, int ncorr, int nblockdata);

#endif
