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
 * Implements helper functions for checkpointing the AWH state and observables history.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "correlationhistory.h"

#include <assert.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/awh-correlation-history.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "correlationgrid.h"

namespace gmx
{

void initCorrelationGridHistory(CorrelationGridHistory *correlationGridHistory,
                                int                     numCorrelationTensors,
                                int                     tensorSize,
                                int                     blockDataListSize)
{
    correlationGridHistory->numCorrelationTensors = numCorrelationTensors;
    correlationGridHistory->tensorSize            = tensorSize;
    correlationGridHistory->blockDataListSize     = blockDataListSize;

    correlationGridHistory->blockDataBuffer.resize(numCorrelationTensors*tensorSize*blockDataListSize);
}

CorrelationGridHistory initCorrelationGridHistoryFromState(const CorrelationGrid &correlationGrid)
{
    CorrelationGridHistory correlationGridHistory;

    initCorrelationGridHistory(&correlationGridHistory, correlationGrid.tensors().size(), correlationGrid.tensorSize(), correlationGrid.blockDataListSize());

    return correlationGridHistory;
}

/* Update the correlation grid history for checkpointing. */
void updateCorrelationGridHistory(CorrelationGridHistory *correlationGridHistory,
                                  const CorrelationGrid  &correlationGrid)
{
    GMX_RELEASE_ASSERT(correlationGridHistory != nullptr, "We need a valid history object");

    gmx::ArrayRef<CorrelationBlockDataHistory> blockDataBuffer = correlationGridHistory->blockDataBuffer;

    /* Store the grid in a linear array */
    size_t bufferIndex = 0;
    for (const CorrelationTensor &tensor : correlationGrid.tensors())
    {
        const int                numDims    = tensor.blockDataList()[0].coordData().size();
        const int                tensorSize = tensor.blockDataList()[0].correlationIntegral().size();

        /* Loop of the tensor elements, ignore the symmetric data */
        int       d1         = 0;
        int       d2         = 0;
        for (int k = 0; k < tensorSize; k++)
        {
            /* BlockData for each correlation element */
            for (const CorrelationBlockData &blockData : tensor.blockDataList())
            {
                const CorrelationBlockData::CoordData &cx  = blockData.coordData()[d1];
                const CorrelationBlockData::CoordData &cy  = blockData.coordData()[d2];

                CorrelationBlockDataHistory           &bdh = blockDataBuffer[bufferIndex];

                bdh.blockSumWeight                       = blockData.blockSumWeight();
                bdh.blockSumSquareWeight                 = blockData.blockSumSquareWeight();
                bdh.blockSumWeightX                      = cx.blockSumWeightX;
                bdh.blockSumWeightY                      = cy.blockSumWeightX;
                bdh.sumOverBlocksSquareBlockWeight       = blockData.sumOverBlocksSquareBlockWeight();
                bdh.sumOverBlocksBlockSquareWeight       = blockData.sumOverBlocksBlockSquareWeight();
                bdh.sumOverBlocksBlockWeightBlockWeightX = cx.sumOverBlocksBlockWeightBlockWeightX;
                bdh.sumOverBlocksBlockWeightBlockWeightY = cy.sumOverBlocksBlockWeightBlockWeightX;
                bdh.previousBlockIndex                   = blockData.previousBlockIndex();
                bdh.blockLength                          = blockData.blockLength();
                bdh.correlationIntegral                  = blockData.correlationIntegral()[k];

                bufferIndex++;
            }

            d1++;
            if (d1 == numDims)
            {
                d1 = 0;
                d2++;
            }
        }
    }

    GMX_RELEASE_ASSERT(bufferIndex == blockDataBuffer.size(), "We should store exactly as many elements as the buffer size");
}

void CorrelationBlockData::restoreFromHistory(const CorrelationBlockDataHistory                  &blockHistory,
                                              const std::vector<CorrelationBlockData::CoordData> &coordData,
                                              const std::vector<double>                          &correlationIntegral)
{
    blockSumWeight_                 = blockHistory.blockSumWeight;
    blockSumSquareWeight_           = blockHistory.blockSumSquareWeight;
    sumOverBlocksSquareBlockWeight_ = blockHistory.sumOverBlocksSquareBlockWeight;
    sumOverBlocksBlockSquareWeight_ = blockHistory.sumOverBlocksBlockSquareWeight;
    previousBlockIndex_             = blockHistory.previousBlockIndex;
    blockLength_                    = blockHistory.blockLength;
    coordData_                      = coordData;
    correlationIntegral_            = correlationIntegral;
}

/* Restore a correlation element from history. */
void CorrelationTensor::restoreFromHistory(const std::vector<CorrelationBlockDataHistory> &blockDataBuffer,
                                           size_t                                         *bufferIndex)
{
    /* Blockdata for each correlation element */
    for (CorrelationBlockData &blockData : blockDataList_)
    {
        /* Correlation elements for each tensor */
        const int numDims    = blockDataList_[0].coordData().size();
        const int tensorSize = blockDataList_[0].correlationIntegral().size();
        int       d1         = 0;
        int       d2         = 0;

        /* Temporary containers to collect data */
        std::vector<CorrelationBlockData::CoordData> coordData(numDims);
        std::vector<double>                          correlationIntegral(tensorSize);
        for (int k = 0; k < tensorSize; k++)
        {
            if (*bufferIndex >= blockDataBuffer.size())
            {
                GMX_THROW(InvalidInputError("Mismatch of the correlation tensor size for the force correlation between checkpoint and simulation. Likely you have provided a checkpoint from a different simulation."));
            }
            const CorrelationBlockDataHistory &blockHistory = blockDataBuffer[*bufferIndex];

            /* To simplify the checkpointing, CorrelationBlockDataHistory
             * duplicates some weight data for all tensor elements.
             * Here we collect the coordinate and tensor data
             * in temporary buffers.
             */
            coordData[d1].blockSumWeightX                      = blockHistory.blockSumWeightX;
            coordData[d2].blockSumWeightX                      = blockHistory.blockSumWeightY;
            coordData[d1].sumOverBlocksBlockWeightBlockWeightX = blockHistory.sumOverBlocksBlockWeightBlockWeightX;
            coordData[d2].sumOverBlocksBlockWeightBlockWeightX = blockHistory.sumOverBlocksBlockWeightBlockWeightY;

            correlationIntegral[k] = blockHistory.correlationIntegral;

            /* Check if we collected all data needed for blockData */
            if (k == tensorSize - 1)
            {
                blockData.restoreFromHistory(blockHistory, coordData, correlationIntegral);
            }

            (*bufferIndex)++;

            d1++;
            if (d1 == numDims)
            {
                d1 = 0;
                d2++;
            }
        }
    }
}

/* Restores the correlation grid state from the correlation grid history. */
void CorrelationGrid::restoreStateFromHistory(const CorrelationGridHistory &correlationGridHistory)
{
    if (tensors_.size() != static_cast<size_t>(correlationGridHistory.numCorrelationTensors))
    {
        GMX_THROW(InvalidInputError("Mismatch of the grid size for the force correlation between checkpoint and simulation. Likely you have provided a checkpoint from a different simulation."));
    }

    /* Extract the state from the linear history array */
    size_t bufferIndex = 0;
    for (CorrelationTensor &tensor : tensors_)
    {
        tensor.restoreFromHistory(correlationGridHistory.blockDataBuffer,
                                  &bufferIndex);
    }

    if (bufferIndex != correlationGridHistory.blockDataBuffer.size())
    {
        GMX_THROW(InvalidInputError("Mismatch of the correlation tensor size for the force correlation between checkpoint and simulation. Likely you have provided a checkpoint from a different simulation."));
    }
}

} // namespace gmx
