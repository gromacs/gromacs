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

#include "gmxpre.h"

#include "correlation-history.h"

#include <assert.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/awh-correlation-history.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "correlation.h"

/*! \brief
 * Packs the correlation history elements from a linear array into the correlation grid struct.
 *
 * \param[in,out] corr_flat           Flattened array of correlation elements.
 * \param[in] corrgrid_hist           Correlation grid to pack correlation elements into.
 */
static void pack_correlation_history(correlation_history_t *corr_flat, correlation_grid_history_t *corrgrid_hist)
{
    int j = 0;
    for (int i  = 0; i < corrgrid_hist->ncorrmatrix; i++)
    {
        corrgrid_hist->corrmatrix[i].corr = &(corr_flat[j]);
        j += corrgrid_hist->ncorr;
    }
}

/*! \brief
 * Packs the block data history structs from a linear array into the correlation grid struct.
 *
 * \param[in,out] blockdata_flat      Flattened array of block data structs.
 * \param[in] corrgrid_hist           Correlation grid to pack block data structs into.
 */
static void pack_blockdata_history(correlation_blockdata_history_t *blockdata_flat, correlation_grid_history_t *corrgrid_hist)
{
    int k = 0;
    for (int i = 0; i < corrgrid_hist->ncorrmatrix; i++)
    {
        for (int j = 0; j < corrgrid_hist->ncorr; j++)
        {
            corrgrid_hist->corrmatrix[i].corr[j].blockdata = &(blockdata_flat[k]);
            k += corrgrid_hist->nblockdata;
        }
    }
}

/* Allocate and initialize correlation grid history. */
void initCorrelationGridHistory(correlation_grid_history_t *corrGridHist, int numCorrTensors, int tensorSize, int numBlockData)
{
    corrGridHist->ncorrmatrix = numCorrTensors;
    corrGridHist->ncorr       = tensorSize;
    corrGridHist->nblockdata  = numBlockData;

    snew(corrGridHist->corrmatrix, numCorrTensors);

    /* Allocate memory for correlation elements and block data in one block and distribute */
    correlation_history_t           *corr_flat;
    correlation_blockdata_history_t *blockdata_flat;
    snew(corr_flat, numCorrTensors*tensorSize);
    snew(blockdata_flat, numCorrTensors*tensorSize*numBlockData);
    pack_correlation_history(corr_flat, corrGridHist);
    pack_blockdata_history(blockdata_flat, corrGridHist);
}

/* Allocate a correlation grid history when restarting from a checkpoint. */
correlation_grid_history_t *initCorrelationGridHistoryFromState(const CorrelationGrid &corrGrid)
{
    correlation_grid_history_t *corrGridHist;

    snew(corrGridHist, 1);

    initCorrelationGridHistory(corrGridHist, corrGrid.corr.size(), corrGrid.tensorSize(), corrGrid.blockDataSize());

    return corrGridHist;
}

/* Update the correlation grid history for checkpointing. */
void updateCorrelationGridHistory(correlation_grid_history_t *corrGridHist, const CorrelationGrid &corrGrid)
{
    GMX_RELEASE_ASSERT(corrGridHist->corrmatrix != nullptr, "AWH force correlation matrix not initialized when updating history");

    /* Matrix for each grid point */
    for (size_t m = 0; m < corrGrid.corr.size(); m++)
    {
        correlation_matrix_history_t *corrTensorHist = &corrGridHist->corrmatrix[m];
        const CorrelationTensor      &corr           = corrGrid.corr[m];

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
            correlation_history_t *corrHist = &corrTensorHist->corr[k];

            /* Blockdatas for each correlation element */
            for (size_t l = 0; l < corr.blockData().size(); l++)
            {
                correlation_blockdata_history_t *blockDataHist = &corrHist->blockdata[l];
                const BlockData                 &blockData     = corr.blockData()[l];
                const CoordData                 &cx            = blockData.coordData[d1];
                const CoordData                 &cy            = blockData.coordData[d2];

                blockDataHist->blocklength          = blockData.blockLength;
                blockDataHist->correlation_integral = blockData.correlationIntegral[l];
                blockDataHist->sum_wx               = cx.blockSumWX;
                blockDataHist->sum_wy               = cy.blockSumWX;
                blockDataHist->sum_w                = blockData.blockSumW;
                blockDataHist->sum_w2_tot           = blockData.simSumSqrW;
                blockDataHist->sum_wbwx             = cx.simSumBlockWBlockWX;
                blockDataHist->sum_wbwy             = cy.simSumBlockWBlockWX;
                blockDataHist->sum_wb2              = blockData.simSumSqrBlockW;
                blockDataHist->blockindex_prev      = blockData.previousBlockIndex;
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
void CorrelationTensor::restoreFromHistory(const correlation_matrix_history_t &corrTensorHist)
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
            const correlation_blockdata_history_t &bdh = corrTensorHist.corr[k].blockdata[l];

            blockData->coordData[d1].blockSumWX          = bdh.sum_wx;
            blockData->coordData[d2].blockSumWX          = bdh.sum_wy;
            blockData->blockSumW                         = bdh.sum_w;
            blockData->simSumSqrW                        = bdh.sum_w2_tot;
            blockData->coordData[d1].simSumBlockWBlockWX = bdh.sum_wbwx;
            blockData->coordData[d2].simSumBlockWBlockWX = bdh.sum_wbwy;
            blockData->simSumSqrBlockW                   = bdh.sum_wb2;
            blockData->previousBlockIndex                = bdh.blockindex_prev;
            blockData->blockLength                       = bdh.blocklength;
            blockData->correlationIntegral[k]            = bdh.correlation_integral;

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
void restoreCorrelationGridStateFromHistory(const correlation_grid_history_t &corrGridHist, CorrelationGrid *corrGrid)
{
    /* Matrix for each grid point */
    for (size_t m = 0; m < corrGrid->corr.size(); m++)
    {
        corrGrid->corr[m].restoreFromHistory(corrGridHist.corrmatrix[m]);
    }
}
