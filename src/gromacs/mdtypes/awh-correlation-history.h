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

/*! \libinternal \file
 *
 * \brief
 * Contains datatypes and function declarations needed by AWH to
 * have its force correlation data checkpointed.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_AWH_CORRELATION_HISTORY_H
#define GMX_MDTYPES_AWH_CORRELATION_HISTORY_H

/*! \cond INTERNAL */

//! Correlation block averaging data.
struct CorrelationBlockdataHistory
{
    double blockSumW;            /**< Sum weights for current block. */
    double blockSumSqrW;         /**< Sum weights^2 for current block. */
    double blockSumWX;           /**< Weighted sum of x for current block. */
    double blockSumWY;           /**< Weighted sum of y for current block. */
    double simSumSqrBlockW;      /**< Sum over blocks of block weight^2. */
    double simSumBlockSqrW;      /**< Sum over blocks of weight^2. */
    double simSumBlockWBlockWX;  /**< Sum over blocks of block weight times blockSumWX. */
    double simSumBlockWBlockWY;  /**< Sum over blocks of block weight times blockSumWY. */
    double blockLength;          /**< The length of each block used for block averaging. */
    int    previousBlockIndex;   /**< The last block index data was added to (needed only for block length in terms of time). */
    double correlationIntegral;  /**< The time integral of the correlation function of x and y, corr(x(0), y(t)). */
};

//! Correlation matrix element
struct CorrelationHistory
{
    CorrelationBlockdataHistory *blockData; /**< The block averaging data */
};

/* The below structs do not contain any history dependence and are just for organizing the data in the structs above. */

//! Correlation matrix.
struct CorrelationTensorHistory
{
    CorrelationHistory *corr; /**< Array with the correlation elements corr(x, y) in the matrix, where x, y are vector components. */
};

//! Grid of local correlation matrices.
struct CorrelationGridHistory
{
    /* These counts here since we curently need them for initializing the correlation grid when reading a checkpoint */
    int                       numCorrTensor; /**< Number correlation tensors in the grid (equal to the number of points). */
    int                       tensorSize;    /**< The number of stored correlation matrix elements. */
    int                       numBlockData;  /**< To be able to increase the block length later on, data is saved for several block lengths. */
    CorrelationTensorHistory *corrTensor;    /**< Correlation tensors. */
};

/*! \endcond */

/*! \brief
 * Allocate and initialize correlation grid history.
 *
 * \param[in,out] corrGridHist  Correlation grid history for master rank.
 * \param[in] numCorrTensors    Number of correlation tensors in the grid.
 * \param[in] tensorSize        Number of correlation elements in each tensor.
 * \param[in] numBlockData      Number of block data structs needed for each correlation element.
 */
void initCorrelationGridHistory(CorrelationGridHistory *corrGridHist, int numCorrTensors, int tensorSize, int numBlockData);

#endif
