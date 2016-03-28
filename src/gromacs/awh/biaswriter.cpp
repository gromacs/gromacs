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
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "bias.h"
#include "grid.h"
#include "pointstate.h"

namespace gmx
{

/*! \brief
 * Map the variable type to a normalization type.
 *
 * \param[in] var     Value for variable type enum.
 * \returns enum value for normalization type.
 */
static Normalization getNormType(AwhVar var)
{
    Normalization normType;

    switch (var)
    {
        case AwhVar::CoordValue:
            normType = Normalization::Coordinate;
            break;
        case AwhVar::Pmf:
        case AwhVar::Bias:
            normType = Normalization::FreeEnergy;
            break;
        case AwhVar::Visits:
        case AwhVar::Weights:
        case AwhVar::Target:
            normType = Normalization::Distribution;
            break;
        default:
            normType = Normalization::None;
            break;
    }

    return normType;
}

/*! \brief
 * Gets the coordinate normalization value for the given dimension.
 *
 * \param[in] bias      The AWH bias.
 * \param[in] dimIndex  Dimensional index.
 * \returns the coordinate normalization value.
 */
static double getCoordNormvalue(const Bias &bias, int dimIndex)
{
    /* AWH may use different units internally but here we convert to user units */
    return bias.dimParams()[dimIndex].scaleInternalToUserInput(1);
}

/*! \brief
 * Gets the normalization value for the given variable.
 *
 * \param[in] var    Value for variable type enum.
 * \param[in] bias   The AWH bias.
 * \param[in] count  Index of the variable.
 * \returns the normalization value.
 */
static double getNormvalue(AwhVar var, const Bias &bias, int count)
{
    double normValue = 0;

    switch (var)
    {
        case AwhVar::CoordValue:
            normValue = getCoordNormvalue(bias, count);
            break;
        case AwhVar::Visits:
        case AwhVar::Weights:
        case AwhVar::Target:
            normValue = static_cast<double>(bias.state().points().size());
            break;
        default:
            break;
    }

    return normValue;
}

Block::Block(int            numPoints,
             Normalization  normType,
             double         normValue) :
    normType(normType),
    normValue(static_cast<float>(normValue)),
    data(numPoints)
{
}

BiasWriter::BiasWriter(const Bias &bias)
{
    constexpr int awhVarCount = static_cast<int>(AwhVar::Count);
    int           varNumBlock[awhVarCount]; /* Number of blocks per variable */

    /* Different variables need different number of blocks. We keep track of the starting block for
       each variable. */
    int blockCount = 0;
    for (int varIndex = 0; varIndex < awhVarCount; varIndex++)
    {
        switch (static_cast<AwhVar>(varIndex))
        {
            case AwhVar::CoordValue:
                varToBlock_[varIndex] = blockCount;
                varNumBlock[varIndex] = bias.ndim();
                break;
            default:
                /* Most variables need one block */
                varToBlock_[varIndex] = blockCount;
                varNumBlock[varIndex] = 1;
                break;
        }
        blockCount += varNumBlock[varIndex];
    }

    /* Initialize the data blocks for each variable */
    for (int varIndex = 0; varIndex < awhVarCount; varIndex++)
    {
        const AwhVar var = static_cast<AwhVar>(varIndex);
        int          numPoints;
        if (var == AwhVar::MetaData)
        {
            numPoints = static_cast<int>(MetaData::Count);
        }
        else
        {
            numPoints = bias.state().points().size();
        }
        for (int b = 0; b < varNumBlock[varIndex]; b++)
        {
            block_.push_back(Block(numPoints,
                                   getNormType(var),
                                   getNormvalue(var, bias, b)));
        }
    }

    /* No real data yet */
    haveDataToWrite_ = false;
}

/*! \brief
 * Normalizes block data for output.
 *
 * \param[in,out] block  The block to normalize.
 * \param[in] bias       The AWH bias.
 */
static void normalizeBlock(Block *block, const Bias &bias)
{
    /* Here we operate on float data (which is accurate enough, since it
     * is statistical data that will never reach full float precision).
     * But since we can have very many data points, we sum into a double.
     */
    double sum       = 0;
    float  minValue  = GMX_FLOAT_MAX;
    float  recipNorm = 0;

    switch (block->normType)
    {
        case Normalization::None:
            break;
        case Normalization::Coordinate:
            /* Normalize coordinate values by a scale factor */
            for (auto &point : block->data)
            {
                point *= block->normValue;
            }
            break;
        case Normalization::FreeEnergy:
            /* Normalize free energy values by subtracting the minimum value */
            for (size_t m = 0; m < block->data.size(); m++)
            {
                if (bias.state().points()[m].inTargetRegion() && block->data[m] < minValue)
                {
                    minValue = block->data[m];
                }
            }
            for (size_t m = 0; m < block->data.size(); m++)
            {
                if (bias.state().points()[m].inTargetRegion())
                {
                    block->data[m] -= minValue;
                }
            }

            break;
        case Normalization::Distribution:
            /* Normalize distribution values by normalizing their sum */
            for (auto &point : block->data)
            {
                sum += point;
            }
            if (sum > 0)
            {
                recipNorm = block->normValue/static_cast<float>(sum);
            }
            for (auto &point : block->data)
            {
                point *= recipNorm;
            }
            break;
        default:
            gmx_incons("Unknown AWH norm type");
    }
}

inline void
BiasWriter::transferVariablePointDataToWriter(AwhVar                    var,
                                              int                       m,
                                              const Bias               &bias,
                                              const std::vector<float> &pmf)
{
    /* The starting block index of this variable. Note that some variables need several (contiguous) blocks. */
    int blockStart = getVarStartBlock(var);
    GMX_ASSERT(m < static_cast<int>(block_[blockStart].data.size()), "Attempt to transfer AWH data to block for point index out of range");

    /* Transfer the point data of this variable to the right block(s) */
    int b = blockStart;
    switch (var)
    {
        case AwhVar::MetaData:
            switch (static_cast<MetaData>(m))
            {
                case MetaData::NumBlock:
                    /* The number of subblocks per awh (needed by gmx_energy) */
                    block_[b].data[m] = static_cast<double>(block_.size());
                    /* Note: a single subblock takes only a single type and we need doubles. */
                    break;
                case MetaData::TargetError:
                    /* The theoretical target error */
                    block_[b].data[m] = bias.params().errorInitial*std::sqrt(bias.params().histSizeInitial/bias.state().histogramSize().histSize());
                    break;
                case MetaData::ScaledSampleWeight:
                    /* The logarithm of the sample weight relative to a sample weight of 1 at the initial time.
                       In the normal case: this will increase in the initial stage and then stay at a constant value. */
                    block_[b].data[m] = bias.state().histogramSize().scaledSampleWeight();
                    break;
                case MetaData::Count:
                    break;
            }
            break;
        case AwhVar::CoordValue:
            {
                const awh_dvec &coordValue = bias.getGridCoordValue(m);
                for (int d = 0; d < bias.ndim(); d++)
                {
                    block_[b].data[m] = coordValue[d];
                    b++;
                }
            }
            break;
        case AwhVar::Pmf:
            block_[b].data[m] = bias.state().points()[m].inTargetRegion() ? pmf[m] : 0;
            break;
        case AwhVar::Bias:
            {
                const awh_dvec &coordValue = bias.getGridCoordValue(m);
                block_[b].data[m] = bias.state().points()[m].inTargetRegion() ? bias.calcConvolvedBias(coordValue) : 0;
            }
            break;
        case AwhVar::Visits:
            block_[b].data[m] = bias.state().points()[m].numVisitsTot();
            break;
        case AwhVar::Weights:
            block_[b].data[m] = bias.state().points()[m].weightsumTot();
            break;
        case AwhVar::Target:
            block_[b].data[m] = bias.state().points()[m].target();
            break;
        default:
            gmx_incons("Unknown AWH output variable");
            break;
    }
}

void BiasWriter::prepareBiasOutput(const Bias &bias)
{
    /* Pack the AWH data into the writer data. */

    /* Evaluate the PMF for all points */
    std::vector<float> *pmf = &block_[getVarStartBlock(AwhVar::Pmf)].data;
    bias.state().getPmf(pmf);

    /* Pack the the data point by point. */
    for (int i = 0; i < static_cast<int>(MetaData::Count); i++)
    {
        transferVariablePointDataToWriter(AwhVar::MetaData, i, bias, *pmf);
    }
    for (int varIndex = 0; varIndex < static_cast<int>(AwhVar::Count); varIndex++)
    {
        const AwhVar var = static_cast<AwhVar>(varIndex);
        if (var == AwhVar::MetaData || !hasVarBlock(var))
        {
            continue;
        }
        for (size_t m = 0; m < bias.state().points().size(); m++)
        {
            transferVariablePointDataToWriter(var, m, bias, *pmf);
        }
    }

    /* For looks of the output, normalize it */
    for (auto &block : block_)
    {
        normalizeBlock(&block, bias);
    }

    haveDataToWrite_ = true;
}

int BiasWriter::writeToEnergySubblocks(t_enxsubblock *sub)
{
    for (size_t b = 0; b < block_.size(); b++)
    {
        sub[b].type = xdr_datatype_float;
        sub[b].nr   = block_[b].data.size();
        sub[b].fval = block_[b].data.data();
    }

    GMX_ASSERT(haveDataToWrite_, "All writers should have the same haveDataToWrite value");
    /* Don't write again until a new dataset has been prepared */
    haveDataToWrite_ = false;

    return block_.size();
}

} // namepace gmx
