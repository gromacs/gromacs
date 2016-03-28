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
 * Implements the CorrelationGrid class to collect correlation statistics on a grid, using several block lengths.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "correlationgrid.h"

#include "gromacs/math/utilities.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

namespace
{

/*! \brief
 * Return the number of block data structs needed for keeping a certain number of blocks.
 *
 * The list start with 1 block and doubles, so we need 1 + 2log(numBlocks).
 *
 * \param[in] numBlocks  Number of blocks.
 * \returns the number of block data structs.
 */
int getBlockDataListSize(int numBlocks)
{
    int blockDataListSize = 1;

    while (numBlocks > (1 << (blockDataListSize - 1)))
    {
        blockDataListSize++;
    }

    GMX_RELEASE_ASSERT((1 << (blockDataListSize - 1)) == numBlocks, "numBlocks should be a power of 2");

    return blockDataListSize;
}

}   // namespace

CorrelationGrid::CorrelationGrid(int                numPoints,
                                 int                numDim,
                                 double             blockLengthInit,
                                 BlockLengthMeasure blockLengthMeasure,
                                 double             dtSample) :
    dtSample(dtSample),
    blockLengthMeasure(blockLengthMeasure)
{
    /* Set the initial block length for the block averaging. The length doesn't really matter
       after the block length has been doubled a few times, as long as it's set small enough */
    if (blockLengthMeasure == BlockLengthMeasure::Weight)
    {
        blockLengthInit = blockLengthInit > 0 ? blockLengthInit : 1;
    }
    else
    {
        blockLengthInit = blockLengthInit > 0 ? blockLengthInit : dtSample;
    }

    /* Set the number of blocks. The number of blocks determines the current span of the data
       and how many different block lengths (nblockdata) we need to keep track of to be able to
       increase the block length later */
    int numBlocks         = CorrelationTensor::c_numCorrelationBlocks;
    int BlockDataListSize = getBlockDataListSize(numBlocks);

    tensors_.resize(numPoints, CorrelationTensor(numDim, BlockDataListSize, blockLengthInit));
}

int CorrelationGrid::getNumBlocks() const
{
    const auto &blockDataList  = tensors()[0].blockDataList();
    double      maxBlockLength = blockDataList.back().blockLength();
    double      minBlockLength = blockDataList[0].blockLength();

    /* If we have a finite block span we have a constant number of blocks, otherwise we are always adding more blocks (and we don't keep track of the number) */
    if (maxBlockLength < GMX_DOUBLE_MAX)
    {
        return static_cast<int>(maxBlockLength/minBlockLength);
    }
    else
    {
        return -1;
    }
}

double CorrelationGrid::getBlockLength() const
{
    /* Return the  minimum blocklength */
    return tensors()[0].blockDataList()[0].blockLength();
}

} // namespace gmx
