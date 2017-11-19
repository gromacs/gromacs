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

#include "biaswriter.h"

#include <assert.h>

#include <cmath>

#include "gromacs/awh/awh.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "bias.h"
#include "correlationgrid.h"
#include "grid.h"
#include "pointstate.h"

namespace gmx
{

namespace
{

/*! \brief
 * Map the output entry type to a normalization type.
 *
 * The data is written to energy file blocks in the order given by
 * the iterator of this map, which is based on the enum value
 * (and matches the order of the lines below).
 */
static const std::map<AwhOutputEntryType, Normalization> outputTypeToNormalization =
{
    { AwhOutputEntryType::MetaData,               Normalization::None },
    { AwhOutputEntryType::CoordValue,             Normalization::Coordinate },
    { AwhOutputEntryType::Pmf,                    Normalization::FreeEnergy },
    { AwhOutputEntryType::Bias,                   Normalization::FreeEnergy },
    { AwhOutputEntryType::Visits,                 Normalization::Distribution },
    { AwhOutputEntryType::Weights,                Normalization::Distribution },
    { AwhOutputEntryType::Target,                 Normalization::Distribution },
    { AwhOutputEntryType::ForceCorrelationVolume, Normalization::Distribution },
    { AwhOutputEntryType::FrictionTensor,         Normalization::None }
};

/*! \brief
 * Gets the coordinate normalization value for the given dimension.
 *
 * \param[in] bias      The AWH bias.
 * \param[in] dimIndex  Dimensional index.
 * \returns the coordinate normalization value.
 */
float getCoordNormalizationValue(const Bias &bias,
                                 int         dimIndex)
{
    /* AWH may use different units internally but here we convert to user units */
    return bias.dimParams()[dimIndex].scaleInternalToUserInput(1);
}

/*! \brief
 * Gets the normalization value for the given output entry type.
 *
 * \param[in] outputType  Output entry type.
 * \param[in] bias        The AWH bias.
 * \param[in] numBlocks   The number of blocks for this output type.
 * \returns the normalization value.
 */
float getNormalizationValue(AwhOutputEntryType  outputType,
                            const Bias         &bias,
                            int                 numBlocks)
{
    float normalizationValue = 0;

    switch (outputType)
    {
        case AwhOutputEntryType::CoordValue:
            normalizationValue = getCoordNormalizationValue(bias, numBlocks);
            break;
        case AwhOutputEntryType::Visits:
        case AwhOutputEntryType::Weights:
        case AwhOutputEntryType::Target:
            normalizationValue = static_cast<float>(bias.state().points().size());
            break;
        case AwhOutputEntryType::ForceCorrelationVolume:
            normalizationValue = static_cast<double>(bias.state().points().size());
            break;
        default:
            break;
    }

    return normalizationValue;
}

}   // namespace

AwhEnergyBlock::AwhEnergyBlock(int            numPoints,
                               Normalization  normalizationType,
                               float          normalizationValue) :
    normalizationType(normalizationType),
    normalizationValue(normalizationValue),
    data_(numPoints)
{
}

BiasWriter::BiasWriter(const Bias &bias)
{
    std::map<AwhOutputEntryType, int> outputTypeNumBlock; /* Number of blocks per output type */

    /* Different output variable types need different number of blocks.
     * We keep track of the starting block for each variable.
     */
    int blockCount = 0;
    for (const auto &pair : outputTypeToNormalization)
    {
        const AwhOutputEntryType outputType = pair.first;
        {
            outputTypeToBlock_[outputType] = blockCount;

            if (outputType == AwhOutputEntryType::CoordValue)
            {
                outputTypeNumBlock[outputType] = bias.ndim();
            }
            else if (outputType == AwhOutputEntryType::FrictionTensor)
            {
                outputTypeNumBlock[outputType] = bias.forceCorrelationGrid().tensorSize();
            }
            else
            {
                /* Most output variable types need one block */
                outputTypeNumBlock[outputType] = 1;
            }
        }
        blockCount += outputTypeNumBlock[outputType];
    }

    /* Initialize the data blocks for each variable */
    for (const auto &pair : outputTypeToNormalization)
    {
        const AwhOutputEntryType outputType = pair.first;
        int                      numPoints;
        if (outputType == AwhOutputEntryType::MetaData)
        {
            numPoints = static_cast<int>(AwhOutputMetaData::Count);
        }
        else
        {
            numPoints = bias.state().points().size();
        }
        for (int b = 0; b < outputTypeNumBlock[outputType]; b++)
        {
            block_.push_back(AwhEnergyBlock(numPoints,
                                            pair.second,
                                            getNormalizationValue(outputType, bias, b)));
        }
    }
}

/*! \brief
 * Normalizes block data for output.
 *
 * \param[in,out] block  The block to normalize.
 * \param[in]     bias   The AWH bias.
 */
static void normalizeBlock(AwhEnergyBlock *block, const Bias &bias)
{
    gmx::ArrayRef<float> data = block->data();

    /* Here we operate on float data (which is accurate enough, since it
     * is statistical data that will never reach full float precision).
     * But since we can have very many data points, we sum into a double.
     */
    double sum       = 0;
    float  minValue  = GMX_FLOAT_MAX;
    float  recipNorm = 0;

    switch (block->normalizationType)
    {
        case Normalization::None:
            break;
        case Normalization::Coordinate:
            /* Normalize coordinate values by a scale factor */
            for (float &point : data)
            {
                point *= block->normalizationValue;
            }
            break;
        case Normalization::FreeEnergy:
            /* Normalize free energy values by subtracting the minimum value */
            for (size_t index = 0; index < data.size(); index++)
            {
                if (bias.state().points()[index].inTargetRegion() && data[index] < minValue)
                {
                    minValue = data[index];
                }
            }
            for (size_t index = 0; index < data.size(); index++)
            {
                if (bias.state().points()[index].inTargetRegion())
                {
                    data[index] -= minValue;
                }
            }

            break;
        case Normalization::Distribution:
            /* Normalize distribution values by normalizing their sum */
            for (float &point : data)
            {
                sum += point;
            }
            if (sum > 0)
            {
                recipNorm = block->normalizationValue/static_cast<float>(sum);
            }
            for (float &point : data)
            {
                point *= recipNorm;
            }
            break;
        default:
            GMX_RELEASE_ASSERT(false, "Unknown AWH normalization type");
            break;
    }
}

void BiasWriter::transferMetaDataToWriter(size_t             metaDataIndex,
                                          AwhOutputMetaData  metaDataType,
                                          const Bias        &bias)
{
    gmx::ArrayRef<float> data = block_[getVarStartBlock(AwhOutputEntryType::MetaData)].data();
    GMX_ASSERT(metaDataIndex < data.size(), "Attempt to transfer AWH meta data to block for index out of range");

    /* Transfer the point data of this variable to the right block(s) */
    switch (metaDataType)
    {
        case AwhOutputMetaData::NumBlock:
            /* The number of subblocks per awh (needed by gmx_energy) */
            data[metaDataIndex] = static_cast<double>(block_.size());
            /* Note: a single subblock takes only a single type and we need doubles. */
            break;
        case AwhOutputMetaData::TargetError:
            /* The theoretical target error */
            data[metaDataIndex] = bias.params().initialErrorInKT*std::sqrt(bias.params().initialHistogramSize/bias.state().histogramSize().histogramSize());
            break;
        case AwhOutputMetaData::ScaledSampleWeight:
            /* The logarithm of the sample weight relative to a sample weight of 1 at the initial time.
               In the normal case: this will increase in the initial stage and then stay at a constant value. */
            data[metaDataIndex] = bias.state().histogramSize().logScaledSampleWeight();
            break;
        case AwhOutputMetaData::Count:
            break;
    }
}

void
BiasWriter::transferPointDataToWriter(AwhOutputEntryType          outputType,
                                      int                         pointIndex,
                                      const Bias                 &bias,
                                      gmx::ArrayRef<const float>  pmf)
{
    /* The starting block index of this output type.
     * Note that some variables need several (contiguous) blocks.
     */
    int blockStart = getVarStartBlock(outputType);
    GMX_ASSERT(pointIndex < static_cast<int>(block_[blockStart].data().size()), "Attempt to transfer AWH data to block for point index out of range");

    const CorrelationGrid &forceCorrelation = bias.forceCorrelationGrid();
    int                    numCorrelation   = forceCorrelation.tensorSize();

    /* Transfer the point data of this variable to the right block(s) */
    int b = blockStart;
    switch (outputType)
    {
        case AwhOutputEntryType::MetaData:
            GMX_RELEASE_ASSERT(false, "MetaData is handled by a different function");
            break;
        case AwhOutputEntryType::CoordValue:
        {
            const awh_dvec &coordValue = bias.getGridCoordValue(pointIndex);
            for (int d = 0; d < bias.ndim(); d++)
            {
                block_[b].data()[pointIndex] = coordValue[d];
                b++;
            }
        }
        break;
        case AwhOutputEntryType::Pmf:
            block_[b].data()[pointIndex] = bias.state().points()[pointIndex].inTargetRegion() ? pmf[pointIndex] : 0;
            break;
        case AwhOutputEntryType::Bias:
        {
            const awh_dvec &coordValue = bias.getGridCoordValue(pointIndex);
            block_[b].data()[pointIndex] = bias.state().points()[pointIndex].inTargetRegion() ? bias.calcConvolvedBias(coordValue) : 0;
        }
        break;
        case AwhOutputEntryType::Visits:
            block_[b].data()[pointIndex] = bias.state().points()[pointIndex].numVisitsTot();
            break;
        case AwhOutputEntryType::Weights:
            block_[b].data()[pointIndex] = bias.state().points()[pointIndex].weightSumTot();
            break;
        case AwhOutputEntryType::Target:
            block_[b].data()[pointIndex] = bias.state().points()[pointIndex].target();
            break;
        case AwhOutputEntryType::ForceCorrelationVolume:
            block_[b].data()[pointIndex] = forceCorrelation.tensors()[pointIndex].getVolumeElement(forceCorrelation.dtSample);
            break;
        case AwhOutputEntryType::FrictionTensor:
            /* Store force correlation in units of friction, i.e. time/length^2 */
            for (int n = 0; n < numCorrelation; n++)
            {
                block_[b].data()[pointIndex] = forceCorrelation.tensors()[pointIndex].getTimeIntegral(n, forceCorrelation.dtSample);
                b++;
            }
            break;
        default:
            GMX_RELEASE_ASSERT(false, "Unknown AWH output variable");
            break;
    }
}

void BiasWriter::prepareBiasOutput(const Bias &bias)
{
    /* Pack the AWH data into the writer data. */

    /* Evaluate the PMF for all points */
    gmx::ArrayRef<float> pmf = block_[getVarStartBlock(AwhOutputEntryType::Pmf)].data();
    bias.state().getPmf(pmf);

    /* Pack the the data point by point.
     * Unfortunately we can not loop over a class enum, so we cast to int.
     * \todo Use strings instead of enum when we port the output to TNG.
     */
    for (int i = 0; i < static_cast<int>(AwhOutputMetaData::Count); i++)
    {
        transferMetaDataToWriter(i, static_cast<AwhOutputMetaData>(i), bias);
    }
    for (const auto &pair : outputTypeToNormalization)
    {
        const AwhOutputEntryType outputType = pair.first;
        /* Skip metadata (transfered above) and unused blocks */
        if (outputType == AwhOutputEntryType::MetaData || !hasVarBlock(outputType))
        {
            continue;
        }
        for (size_t m = 0; m < bias.state().points().size(); m++)
        {
            transferPointDataToWriter(outputType, m, bias, pmf);
        }
    }

    /* For looks of the output, normalize it */
    for (AwhEnergyBlock &block : block_)
    {
        normalizeBlock(&block, bias);
    }
}

int BiasWriter::writeToEnergySubblocks(const Bias    &bias,
                                       t_enxsubblock *sub)
{
    prepareBiasOutput(bias);

    for (size_t b = 0; b < block_.size(); b++)
    {
        sub[b].type = xdr_datatype_float;
        sub[b].nr   = block_[b].data().size();
        sub[b].fval = block_[b].data().data();
    }

    return block_.size();
}

} // namepace gmx
