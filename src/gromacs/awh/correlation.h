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
 * Declares data types and function for AWH to do force correlation
 * statistics.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_CORRELATION_H
#define GMX_AWH_CORRELATION_H

#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

struct correlation_matrix_history_t;

//! Correlation block averaging data.
struct CoordData
{
    /*! \brief Constructor.
     */
    CoordData() :
        blockSumWX(0),
        simSumBlockWBlockWX(0)
    {
    };

    double blockSumWX;          /**< Weighted sum of x for current block. */
    double simSumBlockWBlockWX; /**< Sum over blocks of block weight times sum_wx. */
};

/*! \internal \brief Correlation block averaging weight-only data.
 */
struct BlockData
{
    /*! \brief Constructor.
     *
     * \param[in] numDims          The dimensionality.
     * \param[in] blockLengthInit  The initial block length.
     */
    BlockData(int numDims, double blockLengthInit) :
        blockSumW(0),
        simSumSqrW(0),
        simSumSqrBlockW(0),
        blockLength(blockLengthInit),
        previousBlockIndex(-1),
        coordData(numDims),
        correlationIntegral(numDims*(numDims + 1)/2)
    {
    };

    /*! \brief Adds a weighted data vector to one point in the correlation grid.
     *
     * \param[in] weight   The weight of the data.
     * \param[in] dimData  Data for each dimension.
     */
    void addData(double weight, const double *dimData);

    /*! \brief Adds a filled data block to correlation time integral.
     */
    void addBlockToCorrelationIntegral();

    /* Weight sum data, indentical for all dimensions */
    double blockSumW;          /**< Sum weights for current block. */
    double simSumSqrW;         /**< Sum weights^2 over all samples. */
    double simSumSqrBlockW;    /**< Sum over blocks of block weight^2. */
    double blockLength;        /**< The length of each block used for block averaging */
    int    previousBlockIndex; /**< The last block index data was added to (needed only for block length in terms of time). */

    /* Sums for each coordinate dimension. */
    std::vector<CoordData> coordData; /**< Array with sums for each coordinate dimension. */

    /* Correlation tensor. */
    std::vector<double>    correlationIntegral; /**< Array with the correlation elements corr(x, y) in the tensor, where x, y are vector components. */
};

/*! \internal \brief Correlation data for computing the correlation tensor of one grid point.
 */
class CorrelationTensor
{
    public:
        /*! \brief Constructor.
         *
         * \param[in] numDims          The dimensionality.
         * \param[in] numBlockData     The number of data block data structs.
         * \param[in] blockLengthInit  The initial block length.
         */
        CorrelationTensor(int numDims, int numBlockData, double blockLengthInit);

        /*! Get a const reference to the block data.
         */
        const std::vector<BlockData> &blockData() const
        {
            return blockData_;
        }

        /*! \brief
         * Get the total weight of the data in the correlation matrix.
         *
         * \returns the weight of the added data.
         */
        double getWeight() const;

        /*! \brief Restore a correlation element from history.
         *
         * \param[in] corrTensorHist  Correlation history to restore from.
         */
        void restoreFromHistory(const correlation_matrix_history_t &corrTensorHist);

    private:
        /*! \brief Updates the block length by doubling.
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
         * \param[in] weight               The weight of the data.
         * \param[in] dimData              Data for each dimension.
         * \param[in] blockLengthInWeight  If true, a block is measured in probability weight, otherwise in time.
         * \param[in] t                    The simulation time.
         */
        void addData(double weight, const double *dimData,
                     bool blockLengthInWeight, double t);

    private:
        std::vector<BlockData> blockData_; /**< The block data for different block lenghts. */
};

/*! \internal \brief Grid of local correlation matrices.
 */
struct CorrelationGrid
{
    /*! \brief Constructor.
     *
     * \param[in] numPoints                   Number of points in the grid.
     * \param[in] numDims                     Number of dimensions of the grid.
     * \param[in] blockLengthInit             Initial length of the blocks used for block averaging.
     * \param[in] measureBlockLengthInWeight  If true, a block is measured in probability weight, otherwise in time.
     * \param[in] dtSample                    Time step for sampling correlations.
     */
    CorrelationGrid(int numPoints, int numDims,
                    double blockLengthInit,
                    bool measureBlockLengthInWeight,
                    double dtSample);

    /*! \brief Returns the number of elements in the tensor: dim*(dim+1)/2.
     */
    int tensorSize() const
    {
        GMX_RELEASE_ASSERT(corr.size() > 0, "Should only call tensorSize on a valid grid");

        return corr[0].blockData()[0].correlationIntegral.size();
    }

    /*! \brief Returns the size of the block data arrays.
     */
    int blockDataSize() const
    {
        GMX_RELEASE_ASSERT(corr.size() > 0, "Should only call tensorSize on a valid grid");

        return corr[0].blockData().size();
    }

    const double                   dtSample;            /**< Time in between samples. */
    const bool                     blockLengthInWeight; /**< If true, a block is measured in probability weight, otherwise in time. */
    std::vector<CorrelationTensor> corr;                /**< Correlation tensor grid */
};

/*! \brief Adds a weighted data vector to one point in the correlation grid.
 *
 * \param[in] corrGrid      The correlation grid.
 * \param[in] pointIndex    Index of the point to add data to.
 * \param[in] weight        Weight to assign to the data.
 * \param[in] dimdata       Data vector of the same dimensionality as the grid.
 * \param[in] t             The time when the data was sampled.
 */
void addDataToCorrelationTensor(CorrelationGrid *corrGrid, int pointIndex,
                                double weight, const double *dimdata, double t);

/*! \brief Returns an element of the time integrated correlation tensor at a given point in the grid.
 *
 * The units of the integral are time*(units of data)^2. This will be friction units time/length^2
 * if the data unit is 1/length.
 *
 * The correlation index lists the elements of the upper-triangular correlation matrix row-wise.
 * (TODO: this should ideally not have to be known by the caller.)
 *
 * \param[in] corr         The correlation grid element.
 * \param[in] tensorIndex  Index in the tensor.
 * \param[in] dtSample     The sampling interval length.
 * \returns the integral.
 */
double getCorrelationTimeIntegral(const CorrelationTensor &corr, int tensorIndex, double dtSample);

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
 * \param[in] corr      The correlation grid element.
 * \param[in] dtSample  The sampling interval length.
 * \returns the volume element.
 */
double getCorrelationVolumeElement(const CorrelationTensor &corr, double dtSample);

/* Right now the below functions are only used for an initial log printing. */

//! Get the current blocklength.
double getBlockLength(const CorrelationGrid &corrGrid);

//! Get the current number of blocks.
int getNumBlocks(const CorrelationGrid &corrgrid);

#endif
