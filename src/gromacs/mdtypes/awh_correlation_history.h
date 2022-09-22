/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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

#include <vector>

namespace gmx
{

/*! \cond INTERNAL */

//! Correlation block averaging data.
struct CorrelationBlockDataHistory
{
    double blockSumWeight;       /**< Sum weights for current block. */
    double blockSumSquareWeight; /**< Sum weights^2 for current block. */
    double blockSumWeightX;      /**< Weighted sum of x for current block. */
    double blockSumWeightY;      /**< Weighted sum of y for current block. */
    double sumOverBlocksSquareBlockWeight; /**< Sum over all blocks in the simulation of block weight^2 over the whole simulation. */
    double sumOverBlocksBlockSquareWeight; /**< Sum over all blocks in the simulation of weight^2 over the whole simulation. */
    double sumOverBlocksBlockWeightBlockWeightX; /**< Sum over all blocks in the simulation of block weight times blockSumWeightX over the whole simulation. */
    double sumOverBlocksBlockWeightBlockWeightY; /**< Sum over all blocks in the simulation of block weight times blockSumWeightY over the whole simulation. */
    double blockLength; /**< The length of each block used for block averaging. */
    int    previousBlockIndex; /**< The last block index data was added to (needed only for block length in terms of time). */
    double correlationIntegral; /**< The time integral of the correlation function of x and y, corr(x(0), y(t)). */
};

//! Grid of local correlation matrices.
struct CorrelationGridHistory
{
    /* These counts here since we curently need them for initializing the correlation grid when reading a checkpoint */
    int numCorrelationTensors =
            0; /**< Number correlation tensors in the grid (equal to the number of points). */
    int tensorSize = 0; /**< The number of stored correlation matrix elements. */
    int blockDataListSize =
            0; /**< To be able to increase the block length later on, data is saved for several block lengths for each element. */

    /* We store all tensor sequentially in a buffer */
    std::vector<CorrelationBlockDataHistory> blockDataBuffer; /**< Buffer that contains the correlation data. */
};

/*! \endcond */

/*! \brief
 * Initialize correlation grid history, sets all sizes.
 *
 * \param[in,out] correlationGridHistory  Correlation grid history for main rank.
 * \param[in] numCorrelationTensors       Number of correlation tensors in the grid.
 * \param[in] tensorSize                  Number of correlation elements in each tensor.
 * \param[in] blockDataListSize           The number of blocks in the list of each tensor element.
 */
void initCorrelationGridHistory(CorrelationGridHistory* correlationGridHistory,
                                int                     numCorrelationTensors,
                                int                     tensorSize,
                                int                     blockDataListSize);

} // namespace gmx

#endif /* GMX_MDTYPES_AWH_CORRELATION_HISTORY_H */
