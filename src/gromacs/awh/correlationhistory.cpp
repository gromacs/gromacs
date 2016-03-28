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
 * Implements functions needed by AWH to have its force correlation data checkpointed.
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
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "correlationgrid.h"

namespace gmx
{

/*! \brief
 * Packs the correlation history elements from a linear array into the correlation grid struct.
 *
 * \param[in,out] corr_flat           Flattened array of correlation elements.
 * \param[in] corrgrid_hist           Correlation grid to pack correlation elements into.
 */
static void pack_correlation_history(CorrelationHistory *corr_flat, CorrelationGridHistory *corrgrid_hist)
{
    int j = 0;
    for (int i = 0; i < corrgrid_hist->numCorrTensor; i++)
    {
        corrgrid_hist->corrTensor[i].corr = &(corr_flat[j]);
        j += corrgrid_hist->tensorSize;
    }
}

/*! \brief
 * Packs the block data history structs from a linear array into the correlation grid struct.
 *
 * \param[in,out] blockdata_flat      Flattened array of block data structs.
 * \param[in] corrgrid_hist           Correlation grid to pack block data structs into.
 */
static void pack_blockdata_history(CorrelationBlockdataHistory *blockdata_flat, CorrelationGridHistory *corrgrid_hist)
{
    int k = 0;
    for (int i = 0; i < corrgrid_hist->numCorrTensor; i++)
    {
        for (int j = 0; j < corrgrid_hist->tensorSize; j++)
        {
            corrgrid_hist->corrTensor[i].corr[j].blockData = &(blockdata_flat[k]);
            k += corrgrid_hist->numBlockData;
        }
    }
}

/* Allocate and initialize correlation grid history. */
void initCorrelationGridHistory(CorrelationGridHistory *corrGridHist, int numCorrTensors, int tensorSize, int numBlockData)
{
    corrGridHist->numCorrTensor = numCorrTensors;
    corrGridHist->tensorSize    = tensorSize;
    corrGridHist->numBlockData  = numBlockData;

    snew(corrGridHist->corrTensor, numCorrTensors);

    /* Allocate memory for correlation elements and block data in one block and distribute */
    CorrelationHistory          *corr_flat;
    CorrelationBlockdataHistory *blockdata_flat;
    snew(corr_flat, numCorrTensors*tensorSize);
    snew(blockdata_flat, numCorrTensors*tensorSize*numBlockData);
    pack_correlation_history(corr_flat, corrGridHist);
    pack_blockdata_history(blockdata_flat, corrGridHist);
}

/* Allocate a correlation grid history when restarting from a checkpoint. */
CorrelationGridHistory *initCorrelationGridHistoryFromState(const CorrelationGrid &corrGrid)
{
    CorrelationGridHistory *corrGridHist;

    snew(corrGridHist, 1);

    initCorrelationGridHistory(corrGridHist, corrGrid.corr().size(), corrGrid.tensorSize(), corrGrid.blockDataSize());

    return corrGridHist;
}

/* Update the correlation grid history for checkpointing. */
void updateCorrelationGridHistory(CorrelationGridHistory *corrGridHist, const CorrelationGrid &corrGrid)
{
    GMX_RELEASE_ASSERT(corrGridHist->corrTensor != nullptr, "AWH force correlation matrix not initialized when updating history");

    /* Matrix for each grid point */
    for (size_t m = 0; m < corrGrid.corr().size(); m++)
    {
        CorrelationTensorHistory *corrTensorHist = &corrGridHist->corrTensor[m];
        const CorrelationTensor  &corr           = corrGrid.corr()[m];

        /* Correlation elements for each tensor.
         * Currently the history blockdata struct contains all data needed
         * for a tensor element. This leads to data duplication, but this
         * simplifies the storage.
         */
        const int numDims    = corr.blockData()[0].coordData.size();
        const int tensorSize = corr.blockData()[0].correlationIntegral.size();
        int       d1         = 0;
        int       d2         = 0;
        for (int k = 0; k < tensorSize; k++)
        {
            CorrelationHistory *corrHist = &corrTensorHist->corr[k];

            /* Blockdatas for each correlation element */
            for (size_t l = 0; l < corr.blockData().size(); l++)
            {
                CorrelationBlockdataHistory *blockDataHist = &corrHist->blockData[l];
                const BlockData             &blockData     = corr.blockData()[l];
                const CoordData             &cx            = blockData.coordData[d1];
                const CoordData             &cy            = blockData.coordData[d2];

                blockDataHist->blockSumW           = blockData.blockSumW;
                blockDataHist->blockSumSqrW        = blockData.blockSumSqrW;
                blockDataHist->blockSumWX          = cx.blockSumWX;
                blockDataHist->blockSumWY          = cy.blockSumWX;
                blockDataHist->simSumSqrBlockW     = blockData.simSumSqrBlockW;
                blockDataHist->simSumBlockSqrW     = blockData.simSumBlockSqrW;
                blockDataHist->simSumBlockWBlockWX = cx.simSumBlockWBlockWX;
                blockDataHist->simSumBlockWBlockWY = cy.simSumBlockWBlockWX;
                blockDataHist->previousBlockIndex  = blockData.previousBlockIndex;
                blockDataHist->blockLength         = blockData.blockLength;
                blockDataHist->correlationIntegral = blockData.correlationIntegral[l];
            }

            d1++;
            if (d1 == numDims)
            {
                d1 = 0;
                d2++;
            }
        }
    }
}

/* Restore a correlation element from history. */
void CorrelationTensor::restoreFromHistory(const CorrelationTensorHistory &corrTensorHist)
{
    /* Blockdata for each correlation element */
    for (size_t l = 0; l < blockData_.size(); l++)
    {
        BlockData *blockData = &blockData_[l];

        /* Correlation elements for each tensor */
        const int numDims    = blockData_[0].coordData.size();
        const int tensorSize = blockData_[0].correlationIntegral.size();
        int       d1         = 0;
        int       d2         = 0;
        for (int k = 0; k < tensorSize; k++)
        {
            const CorrelationBlockdataHistory &bdh = corrTensorHist.corr[k].blockData[l];

            blockData->blockSumW                         = bdh.blockSumW;
            blockData->blockSumSqrW                      = bdh.blockSumSqrW;
            blockData->coordData[d1].blockSumWX          = bdh.blockSumWX;
            blockData->coordData[d2].blockSumWX          = bdh.blockSumWY;
            blockData->simSumSqrBlockW                   = bdh.simSumSqrBlockW;
            blockData->simSumBlockSqrW                   = bdh.simSumBlockSqrW;
            blockData->coordData[d1].simSumBlockWBlockWX = bdh.simSumBlockWBlockWX;
            blockData->coordData[d2].simSumBlockWBlockWX = bdh.simSumBlockWBlockWY;
            blockData->previousBlockIndex                = bdh.previousBlockIndex;
            blockData->blockLength                       = bdh.blockLength;
            blockData->correlationIntegral[k]            = bdh.correlationIntegral;

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
void CorrelationGrid::restoreStateFromHistory(const CorrelationGridHistory &corrGridHist)
{
    GMX_RELEASE_ASSERT(corr_.size() == static_cast<size_t>(corrGridHist.numCorrTensor), "The size of the grid in the state and history should match");

    /* Matrix for each grid point */
    for (size_t m = 0; m < corr_.size(); m++)
    {
        corr_[m].restoreFromHistory(corrGridHist.corrTensor[m]);
    }
}

} // namespace gmx
