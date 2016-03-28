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
 * Implements AWH force correlation statistics.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "correlationgrid.h"

#include <assert.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "math.h"

namespace gmx
{

/* Get the current number of blocks. */
int getNumBlocks(const CorrelationGrid &corrGrid)
{
    const std::vector<BlockData> &blockData      = corrGrid.corr()[0].blockData();
    double                        maxBlockLength = blockData[blockData.size() - 1].blockLength;
    double                        minBlockLength = blockData[0].blockLength;

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

/* Get the current blocklength. */
double getBlockLength(const CorrelationGrid &corrGrid)
{
    /* Return the  minimum blocklength */
    return corrGrid.corr()[0].blockData()[0].blockLength;
}

/* Returns an element of the time integrated correlation matrix at a given point in the grid. */
double getCorrelationTimeIntegral(const CorrelationTensor &corr,
                                  int tensorIndex, double dtSample)
{
    const BlockData &blockData = corr.blockData()[0];
    double           weight    = blockData.simSumBlockSqrW;
    double           corrsum   = 0;
    if (weight > 0)
    {
        corrsum = blockData.correlationIntegral[tensorIndex]/weight;
    }

    return 0.5*corrsum*dtSample;
}

/* Returns the volume element of the correlation metric. */
double getCorrelationVolumeElement(const CorrelationTensor &corr, double dtSample)
{
    double det, a, b, c, d, e, f;

    switch (corr.blockData()[0].correlationIntegral.size())
    {
        case 1:
            /* 1-dimensional tensor: [a] */
            det = getCorrelationTimeIntegral(corr, 0, dtSample);
            break;
        case 3:
            /* 2-dimensional tensor: [a b; b c] */
            a   = getCorrelationTimeIntegral(corr, 0, dtSample);
            b   = getCorrelationTimeIntegral(corr, 1, dtSample);
            c   = getCorrelationTimeIntegral(corr, 2, dtSample);

            det = a*c - b*b;
            break;
        case 6:
            /* 3-dimensional tensor: [a b d; b c e; d e f] */
            a   = getCorrelationTimeIntegral(corr, 0, dtSample);
            b   = getCorrelationTimeIntegral(corr, 1, dtSample);
            c   = getCorrelationTimeIntegral(corr, 2, dtSample);
            d   = getCorrelationTimeIntegral(corr, 3, dtSample);
            e   = getCorrelationTimeIntegral(corr, 4, dtSample);
            f   = getCorrelationTimeIntegral(corr, 5, dtSample);

            det = a*c*f + 2*b*d*e - d*c*d - b*b*f - a*e*e;
            break;
        default:
            det = 0;
            /* meh */
    }

    /* Returns 0 if no data, not supported number of dims
       or not enough data to give a positive determinant (as it should be) */
    return det > 0 ? std::sqrt(det) : 0;
}

/* Updates the block length by doubling. */
void CorrelationTensor::doubleBlockLengths()
{
    /* We need to shift the data so that a given blockdata gets the data for double the block length.
       The data for the shortest block length is not needed anymore. */

    for (size_t i = 0; i < blockData_.size() - 1; i++)
    {
        blockData_[i] = blockData_[i + 1];
    }

    /* The blockdata which has 1 block is the same as the old one but with double the block length */
    blockData_[blockData_.size() - 1].blockLength *= 2.;
}

/* Updates the block length such that data fits. */
void CorrelationTensor::updateBlockLengths(double samplingLength)
{
    /* How many times do we need to double the longest block length to fit the data? */
    double longestLength = blockData_[blockData_.size() - 1].blockLength;
    int    numDoublings  = 0;
    while (samplingLength > longestLength)
    {
        numDoublings++;
        longestLength *= 2;
    }

    while (numDoublings > 0)
    {
        doubleBlockLengths();
        numDoublings--;
    }
}

/* Adds a filled data block to correlation time integral. */
void BlockData::addBlockToCorrelationIntegral()
{
    const bool firstBlock = (simSumSqrBlockW == 0);

    if (!firstBlock)
    {
        const int numDims     = coordData.size();
        int       tensorIndex = 0;
        for (int d1 = 0; d1 < numDims; d1++)
        {
            const CoordData &c1 = coordData[d1];

            for (int d2 = 0; d2 <= d1; d2++)
            {
                const CoordData &c2 = coordData[d2];

                /* Need the old average, before the data of this block was added */
                double avg_x_old    = c1.simSumBlockWBlockWX/simSumSqrBlockW;

                double sum_wbwy_new = c2.simSumBlockWBlockWX + blockSumW*c2.blockSumWX;
                double sum_wb2_new  = simSumSqrBlockW + gmx::square(blockSumW);
                double avg_y_new    = sum_wbwy_new/sum_wb2_new;

                double sum_dx_old   = c1.blockSumWX - avg_x_old*blockSumW;
                double sum_dy_new   = c2.blockSumWX - avg_y_new*blockSumW;

                /* Update the correlation integral using the changes in averages. */
                correlationIntegral[tensorIndex] += sum_dx_old*sum_dy_new;
                tensorIndex++;
            }
        }
    }

    /* Add the weights of the block to the block sums and clear the weights */
    simSumSqrBlockW           += gmx::square(blockSumW);
    simSumBlockSqrW           += blockSumSqrW;
    for (auto &c : coordData)
    {
        c.simSumBlockWBlockWX += blockSumW*c.blockSumWX;
        /* Reset */
        c.blockSumWX           = 0;
    }
    /* Reset */
    blockSumW                  = 0;
    blockSumSqrW               = 0;
}

/*! \brief
 * Get the block index for the current block and simulation length.
 *
 * \param[in] blockLength  Block length.
 * \param[in] length       Sampling length of all data, in time or weight.
 * \returns the block index.
 */
static int getBlockIndex(double blockLength, double length)
{
    return static_cast<int>(length/blockLength);
}

/* Adds new data to block. */
void BlockData::addData(double weight, const double *dimData)
{
    blockSumW    += weight;
    blockSumSqrW += weight*weight;

    for (size_t d = 0; d < coordData.size(); d++)
    {
        coordData[d].blockSumWX += weight*dimData[d];
    }
}

/* Get the total weight of the data in the correlation matrix. */
double CorrelationTensor::getWeight() const
{
    /* The last blockdata has only 1 block containing all data */
    return blockData()[blockData().size() - 1].blockSumW;
}

/* Adds a weighted data vector to one point in the correlation grid. */
void CorrelationTensor::addData(double weight, const double *dimData,
                                bool blockLengthInWeight, double t)
{
    if (weight < 1e-6)
    {
        /* Nothing to add */
        return;
    }

    /*  The sampling length is measured either in the total (local) weight or the current time */
    double samplingLength = blockLengthInWeight ? getWeight() + weight : t;

    /* Make sure the blocks are long enough to fit all data */
    updateBlockLengths(samplingLength);

    /* Store the data for each block length considered. First update the longest block which has all data since it's
       used for updating the correlation function for the other block lengths. */
    for (size_t i = 0; i < blockData_.size() - 1; i++)
    {
        BlockData *bd = &blockData_[i];

        /* Find current block index given block length. */
        int blockIndex = getBlockIndex(bd->blockLength, samplingLength);

        if (bd->previousBlockIndex >= 0 && bd->previousBlockIndex != blockIndex)
        {
            /* Changed block. Update correlation with old data before adding to new block. */
            bd->addBlockToCorrelationIntegral();
        }

        /* Keep track of which block index data was last added to */
        bd->previousBlockIndex = blockIndex;

        /* Store the data */
        bd->addData(weight, dimData);
    }

    /* The last blockdata has only 1 block which contains all data so far.
     * Add the data for the largest block length.
     */
    blockData_[blockData_.size() - 1].addData(weight, dimData);
}

/* Adds a weighted data vector to one point in the correlation grid. */
void CorrelationGrid::addData(int pointIndex, double weight, const double *dimData, double t)
{
    corr_[pointIndex].addData(weight, dimData, blockLengthInWeight, t);
}

/* Constructor. */
CorrelationTensor::CorrelationTensor(int numDims, int numBlockData,
                                     double blockLengthInit)
{
    int scaling = 1;

    for (int n = 0; n < numBlockData; n++)
    {
        blockData_.push_back(BlockData(numDims, scaling*blockLengthInit));
        scaling <<= 1; /* Double block length */
    }
}

/*! \brief
 * Return the number of block data structs needed for keeping a certain number of blocks.
 *
 * \param[in] numBlocks  Number of blocks.
 * \returns the number of block data structs.
 */
static int getNumBlockData(int numBlocks)
{
    int numBlockData = 0;

    while (numBlocks > 0)
    {
        numBlocks >>= 1; /* divide by 2 */
        numBlockData++;
    }
    return numBlockData;
}

/* Constructor. */
CorrelationGrid::CorrelationGrid(int numPoints, int numDims,
                                 double blockLengthInit,
                                 bool measureBlockLengthInWeight,
                                 double dtSample) :
    dtSample(dtSample),
    blockLengthInWeight(measureBlockLengthInWeight)
{
    /* Set the initial block length for the block averaging. The length doesn't really matter
       after the block length has been doubled a few times, as long as it's set small enough */
    if (blockLengthInWeight)
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
    int numBlocks    = CorrelationTensor::c_numCorrelationBlocks;
    int numBlockData = getNumBlockData(numBlocks);

    corr_.resize(numPoints, CorrelationTensor(numDims, numBlockData, blockLengthInit));
}

} // namespace gmx
