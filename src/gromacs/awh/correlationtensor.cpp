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
 * Implements the CorrelationTensor class for correlation tensor statistics.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "correlationtensor.h"

#include <assert.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

namespace
{

/*! \brief
 * Get the block index for the current block and simulation length.
 *
 * This is simply how many blocks of length \p blockLength fit completely
 * into \p totalAccumulatedLength, which is either the current time minus
 * the start time, which time weighting, or the total weight.
 *
 * \param[in] blockLength               Block length.
 * \param[in] currentAccumulatedLength  Sampling length of all data collected up to the current step, in time or weight.
 * \returns the block index.
 */
int getBlockIndex(double blockLength,
                  double currentAccumulatedLength)
{
    return static_cast<int>(currentAccumulatedLength/blockLength);
}

}   // namespace

double CorrelationTensor::getTimeIntegral(int    tensorIndex,
                                          double dtSample) const
{
    const CorrelationBlockData &blockData           = blockDataList_[0];
    double                      weight              = blockData.sumOverBlocksBlockSquareWeight();
    double                      correlationIntegral = 0;
    if (weight > 0)
    {
        correlationIntegral = blockData.correlationIntegral()[tensorIndex]/weight;
    }

    return 0.5*correlationIntegral*dtSample;
}

double CorrelationTensor::getVolumeElement(double dtSample) const
{
    double det;

    switch (blockDataList_[0].correlationIntegral().size())
    {
        case 1:
            /* 1-dimensional tensor: [a] */
            det = getTimeIntegral(0, dtSample);
            break;
        case 3:
        {
            /* 2-dimensional tensor: [a b; b c] */
            double a = getTimeIntegral(0, dtSample);
            double b = getTimeIntegral(1, dtSample);
            double c = getTimeIntegral(2, dtSample);

            det      = a*c - b*b;
        }
        break;
        case 6:
        {
            /* 3-dimensional tensor: [a b d; b c e; d e f] */
            double a = getTimeIntegral(0, dtSample);
            double b = getTimeIntegral(1, dtSample);
            double c = getTimeIntegral(2, dtSample);
            double d = getTimeIntegral(3, dtSample);
            double e = getTimeIntegral(4, dtSample);
            double f = getTimeIntegral(5, dtSample);

            det      = a*c*f + 2*b*d*e - d*c*d - b*b*f - a*e*e;
        }
        break;
        default:
            det = 0;
            /* meh */
            break;
    }

    /* Returns 0 if no data, not supported number of dims
       or not enough data to give a positive determinant (as it should be) */
    return det > 0 ? std::sqrt(det) : 0;
}

void CorrelationTensor::doubleBlockLengths()
{
    /* We need to shift the data so that a given blockdata gets the data for double the block length.
       The data for the shortest block length is not needed anymore. */

    for (size_t i = 0; i < blockDataList_.size() - 1; i++)
    {
        blockDataList_[i] = blockDataList_[i + 1];
    }

    /* The blockdata which has 1 block is the same as the old one but with double the block length */
    blockDataList_.back().doubleBlockLength();
}

void CorrelationTensor::updateBlockLengths(double samplingLength)
{
    /* How many times do we need to double the longest block length to fit the data? */
    double longestLength = blockDataList_.back().blockLength();
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
void CorrelationBlockData::addBlockToCorrelationIntegral()
{
    const bool firstBlock = (sumOverBlocksSquareBlockWeight_ == 0);

    if (!firstBlock)
    {
        const int numDim      = coordData_.size();
        int       tensorIndex = 0;
        for (int d1 = 0; d1 < numDim; d1++)
        {
            const CoordData &c1 = coordData_[d1];

            for (int d2 = 0; d2 <= d1; d2++)
            {
                const CoordData &c2 = coordData_[d2];

                /* Compute change in correlaion integral due to adding
                 * the block through computing the difference of the block
                 * average with the old average for one component (we use x)
                 * and with the new component (we use y).
                 */
                /* Need the old average, before the data of this block was added */
                GMX_ASSERT(sumOverBlocksSquareBlockWeight_, "Denominator should be > 0 (should be guaranteed by the conditional above)");
                double oldAverageX              = c1.sumOverBlocksBlockWeightBlockWeightX/sumOverBlocksSquareBlockWeight_;

                double newSumWeightBlockWeightY = c2.sumOverBlocksBlockWeightBlockWeightX + blockSumWeight_*c2.blockSumWeightX;
                double newSumSquareBlockWeight  = sumOverBlocksSquareBlockWeight_ + gmx::square(blockSumWeight_);

                GMX_ASSERT(newSumSquareBlockWeight > 0, "Denominator should be > 0");
                double newAverageY              = newSumWeightBlockWeightY/newSumSquareBlockWeight;

                double diffBlockWithOldAverageX = c1.blockSumWeightX - oldAverageX*blockSumWeight_;
                double diffBlockWithNewAverageY = c2.blockSumWeightX - newAverageY*blockSumWeight_;

                /* Update the correlation integral using the changes in averages. */
                correlationIntegral_[tensorIndex] += diffBlockWithOldAverageX*diffBlockWithNewAverageY;
                tensorIndex++;
            }
        }
    }

    /* Add the weights of the block to the block sums and clear the weights */
    sumOverBlocksSquareBlockWeight_ += gmx::square(blockSumWeight_);
    sumOverBlocksBlockSquareWeight_ += blockSumSquareWeight_;
    for (auto &c : coordData_)
    {
        c.sumOverBlocksBlockWeightBlockWeightX += blockSumWeight_*c.blockSumWeightX;
        /* Reset */
        c.blockSumWeightX                           = 0;
    }
    /* Reset */
    blockSumWeight_       = 0;
    blockSumSquareWeight_ = 0;
}

void CorrelationTensor::addData(double                       weight,
                                gmx::ArrayRef<const double>  data,
                                bool                         blockLengthInWeight,
                                double                       t)
{
    /* We should avoid adding data with very small weight to avoid
     * divergence close to 0/0. The total spread weight for each sample is 1,
     * so 1e-6 is a completely negligible amount.
     */
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
    for (size_t i = 0; i < blockDataList_.size() - 1; i++)
    {
        CorrelationBlockData &bd = blockDataList_[i];

        /* Find current block index given the block length. */
        int blockIndex = getBlockIndex(bd.blockLength(), samplingLength);

        if (bd.previousBlockIndex() >= 0 &&
            bd.previousBlockIndex() != blockIndex)
        {
            /* Changed block. Update correlation with old data before adding to new block. */
            bd.addBlockToCorrelationIntegral();
        }

        /* Keep track of which block index data was last added to */
        bd.setPreviousBlockIndex(blockIndex);

        /* Store the data */
        bd.addData(weight, data);
    }

    /* The last blockdata has only 1 block which contains all data so far.
     * Add the data for the largest block length.
     */
    blockDataList_.back().addData(weight, data);
}

CorrelationTensor::CorrelationTensor(int    numDim,
                                     int    numBlockData,
                                     double blockLengthInit)
{
    unsigned int scaling = 1;

    GMX_RELEASE_ASSERT(numBlockData < static_cast<int>(sizeof(scaling)*8), "numBlockData should we smaller than the number of bits in scaling");

    for (int n = 0; n < numBlockData; n++)
    {
        blockDataList_.push_back(CorrelationBlockData(numDim, scaling*blockLengthInit));
        scaling <<= 1; /* Double block length */
    }
}

} // namespace gmx
