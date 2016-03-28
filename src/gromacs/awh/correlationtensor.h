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
 * Declares the CorrelationTensor class for correlation tensor statistics.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_CORRELATIONTENSOR_H
#define GMX_AWH_CORRELATIONTENSOR_H

#include <cstddef>

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

struct CorrelationBlockDataHistory;

/*! \internal \brief Correlation block averaging weight-only data.
 */
class CorrelationBlockData
{
    public:
        /*! \internal \brief Correlation block averaging data.
         */
        struct CoordData
        {
            /*! \brief Constructor.
             */
            CoordData() :
                blockSumWeightX(0),
                sumOverBlocksBlockWeightBlockWeightX(0)
            {
            };

            double blockSumWeightX;                      /**< Weighted sum of x for current block. */
            double sumOverBlocksBlockWeightBlockWeightX; /**< Sum over all blocks in the simulation of block weight times sum_wx. */
        };

        /*! \brief Constructor.
         *
         * \param[in] numDim           The dimensionality.
         * \param[in] blockLengthInit  The initial block length.
         */
        CorrelationBlockData(int    numDim,
                             double blockLengthInit) :
            blockSumWeight_(0),
            blockSumSquareWeight_(0),
            sumOverBlocksSquareBlockWeight_(0),
            sumOverBlocksBlockSquareWeight_(0),
            blockLength_(blockLengthInit),
            previousBlockIndex_(-1),
            coordData_(numDim),
            correlationIntegral_(numDim*(numDim + 1)/2)
        {
        };

        /*! \brief Restore the state from history.
         *
         * \param[in] blockHistory         The block data history containing the weight sums.
         * \param[in] coordData            The coordinate data.
         * \param[in] correlationIntegral  The correlation integral for all tensor elements.
         */
        void restoreFromHistory(const CorrelationBlockDataHistory &blockHistory,
                                const std::vector<CoordData>      &coordData,
                                const std::vector<double>         &correlationIntegral);

        /*! \brief Adds a weighted data vector to one point in the correlation grid.
         *
         * \param[in] weight  The weight of the data.
         * \param[in] data    One data point for each grid dimension.
         */
        void addData(double                      weight,
                     gmx::ArrayRef<const double> data)
        {
            GMX_ASSERT(data.size() == coordData_.size(), "Size of data should match the size of coordData");

            blockSumWeight_       += weight;
            blockSumSquareWeight_ += weight*weight;

            for (size_t d = 0; d < coordData_.size(); d++)
            {
                coordData_[d].blockSumWeightX += weight*data[d];
            }
        };

        /*! \brief Adds a filled data block to correlation time integral.
         */
        void addBlockToCorrelationIntegral();

        /*! \brief Returns the sum weights for current block. */
        double blockSumWeight() const
        {
            return blockSumWeight_;
        };

        /*! \brief Returns the sum weights^2 for current block. */
        double blockSumSquareWeight() const
        {
            return blockSumSquareWeight_;
        };

        /*! \brief Returns the sum over blocks of block weight^2. */
        double sumOverBlocksSquareBlockWeight() const
        {
            return sumOverBlocksSquareBlockWeight_;
        };

        /*! \brief Returns the sum over blocks of weight^2. */
        double sumOverBlocksBlockSquareWeight() const
        {
            return sumOverBlocksBlockSquareWeight_;
        };

        /*! \brief Returns the length of each block used for block averaging. */
        double blockLength() const
        {
            return blockLength_;
        };

        /*! \brief Double the length of each block used for block averaging. */
        void doubleBlockLength()
        {
            blockLength_ *= 2;
        };

        /*! \brief Return the last block index data was added to (needed only for block length in terms of time). */
        int previousBlockIndex() const
        {
            return previousBlockIndex_;
        }

        /*! \brief Set the last block index data was added to.
         *
         * \param[in] blockIndex  The block index.
         */
        void setPreviousBlockIndex(int blockIndex)
        {
            previousBlockIndex_ = blockIndex;
        }

        /*! \brief Return sums for each coordinate dimension. */
        const std::vector<CoordData> &coordData() const
        {
            return coordData_;
        };

        /*! \brief Return the correlation integral tensor. */
        const std::vector<double> &correlationIntegral() const
        {
            return correlationIntegral_;
        };

    private:
        /* Weight sum data, indentical for all dimensions */
        double blockSumWeight_;                 /**< Sum weights for current block. */
        double blockSumSquareWeight_;           /**< Sum weights^2 for current block. */
        double sumOverBlocksSquareBlockWeight_; /**< Sum over all blocks in the simulation of block weight^2. */
        double sumOverBlocksBlockSquareWeight_; /**< Sum over all blocks in the simulation of weight^2. */
        double blockLength_;                    /**< The length of each block used for block averaging */
        int    previousBlockIndex_;             /**< The last block index data was added to (needed only for block length in terms of time). */

        /* Sums for each coordinate dimension. */
        std::vector<CoordData> coordData_;      /**< Array with sums for each coordinate dimension. */

        /* Correlation tensor. */
        std::vector<double> correlationIntegral_; /**< Array with the correlation elements corr(x, y) in the tensor, where x, y are vector components. */
};

/*! \internal
 * \brief Correlation data for computing the correlation tensor of one grid point.
 *
 * The time integrated autocorrelation of the desired quantity is computed using
 * block averages, which is a computationally efficient and low memory method.
 * Most of the work here goes into computing the block averages for weights
 * and the coordinate quantity. This is done for a number of blocks in
 * the range of \p c_numCorrelationBlocks/2 + 1 to \p c_numCorrelationBlocks,
 * depending on the current simulation length. As the simulation time
 * progresses, the blocks get longer. This is implemented in an efficient
 * way by keeping track of log2(\p c_numCorrelationBlocks) \p BlockData
 * data blocks with block length increasing progressively by a factor of 2.
 * Once \p c_numCorrelationBlocks are reached, all block lengths are doubled.
 */
class CorrelationTensor
{
    public:
        /*! \brief 64 blocks is a good trade-off between signal and noise */
        static constexpr int c_numCorrelationBlocks = 64;

        /*! \brief Constructor.
         *
         * \param[in] numDim          The dimensionality.
         * \param[in] numBlockData     The number of data block data structs.
         * \param[in] blockLengthInit  The initial block length.
         */
        CorrelationTensor(int    numDim,
                          int    numBlockData,
                          double blockLengthInit);

        /*! \brief Get a const reference to the list of block data.
         */
        const std::vector<CorrelationBlockData> &blockDataList() const
        {
            return blockDataList_;
        }

        /*! \brief
         * Get the total weight of the data in the correlation matrix.
         *
         * \returns the weight of the added data.
         */
        double getWeight() const
        {
            /* The last blockdata has only 1 block containing all data */
            return blockDataList().back().blockSumWeight();
        };

        /*! \brief Restore a correlation element from history.
         *
         * \param[in]     blockDataBuffer  The linear correlation grid data history buffer.
         * \param[in,out] bufferIndex      The index in \p blockDataBuffer to start reading, is increased with the number of blocks read.
         */
        void restoreFromHistory(const std::vector<CorrelationBlockDataHistory> &blockDataBuffer,
                                size_t                                         *bufferIndex);

    private:
        /*! \brief Updates the block length by doubling.
         *
         * The length of all blocks is doubled. This is achieved by removing
         * the shortest block, moving all other blocks and duplicating
         * the data of longest block to a nw block of double length (but
         * currenly only half filled with data).
         */
        void doubleBlockLengths();

        /*! \brief Updates the block length such that data fits.
         *
         * \param[in] samplingLength  Sampling length of all data, in time or weight.
         */
        void updateBlockLengths(double samplingLength);

    public:
        /*! \brief Adds a weighted data vector to one point in the correlation grid.
         *
         * \note To avoid rounding noise, data with weight smaller than 1e-6
         *       is ignored.
         *
         * \param[in] weight               The weight of the data.
         * \param[in] data                 One data point for each grid dimension.
         * \param[in] blockLengthInWeight  If true, a block is measured in probability weight, otherwise in time.
         * \param[in] t                    The simulation time.
         */
        void addData(double                       weight,
                     gmx::ArrayRef<const double>  data,
                     bool                         blockLengthInWeight,
                     double                       t);

        /*! \brief Returns an element of the time integrated correlation tensor at a given point in the grid.
         *
         * The units of the integral are time*(units of data)^2. This will be friction units time/length^2
         * if the data unit is 1/length.
         *
         * The correlation index lists the elements of the upper-triangular
         * correlation matrix row-wise, so e.g. in 3D:
         * 0 (0,0), 1 (1,0), 2 (1,1), 3 (2,0), 4 (2,1), 5 (2,2).
         * (TODO: this should ideally not have to be known by the caller.)
         *
         * \param[in] tensorIndex  Index in the tensor.
         * \param[in] dtSample     The sampling interval length.
         * \returns the integral.
         */
        double getTimeIntegral(int    tensorIndex,
                               double dtSample) const;

        /*! \brief
         * Returns the volume element of the correlation metric.
         *
         * The matrix of the metric equals the time-integrated correlation matrix. The volume element of
         * the metric therefore equals the square-root of the absolute value of its determinant according
         * to the standard formula for a volume element in a metric space.
         *
         * Since the units of the metric matrix elements are time*(units of data)^2, the volume element has
         * units of (sqrt(time)*(units of data))^(ndim of data).
         *
         * \param[in] dtSample  The sampling interval length.
         * \returns the volume element.
         */
        double getVolumeElement(double dtSample) const;

    private:
        std::vector<CorrelationBlockData> blockDataList_; /**< The block data for different, consecutively doubling block lengths. */
};

}      // namespace gmx

#endif /* GMX_AWH_CORRELATIONTENSOR_H */
