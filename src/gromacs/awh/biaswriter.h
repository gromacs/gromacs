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
 * This file contains the BiasWriter class that prepares and writes data of a Bias to an energy file.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_BIASWRITER_H
#define GMX_AWH_BIASWRITER_H

#include <map>
#include <vector>

#include "gromacs/fileio/enxio.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_multisim_t;
struct t_enxsubblock;

namespace gmx
{
class Bias;

/* TODO: the post-simulations AWH reader and this AWH writer are totally
 * disconnected although they read/write the same data. I'm not sure how
 * to handle that or if it should be left as it is until the writing is done
 * in a differen format (i.e. TNG) than the current energy file.
 */

//! Enum with the AWH variables to write.
enum class AwhOutputEntryType
{
    MetaData,               //!< Meta data.
    CoordValue,             //!< Coordinate value.
    Pmf,                    //!< The pmf.
    Bias,                   //!< The bias.
    Visits,                 //!< The number of visits.
    Weights,                //!< The weights.
    Target,                 //!< The target distribition.
    ForceCorrelationVolume, //!< The volume of the force correlation tensor.
    FrictionTensor          //!< The full friction tensor.
};

//! Enum with the types of metadata to write.
enum class AwhOutputMetaData
{
    NumBlock,           //!< The number of blocks.
    TargetError,        //!< The target error.
    ScaledSampleWeight, //!< The logarithm of the sample weight relative to a sample weight of 1 at the initial time.
    Count               //!< The number of enum values, not including Count.
};

//! Enum with different ways of normalizing the output.
enum class Normalization
{
    None,         //!< No normalization.
    Coordinate,   //!< Scale using the internal/user input coordinate scaling factor.
    FreeEnergy,   //!< Normalize free energy values by subtracting the minimum value.
    Distribution  //!< Normalize the distribution to 1.
};

/*! \internal \brief AWH output data block that can be written to an energy file block.
 */
class AwhEnergyBlock
{
    public:
        /*! \brief Constructor
         *
         * \param[in] numPoints           Number of points in block.
         * \param[in] normalizationType   Type of normalization.
         * \param[in] normalizationValue  Value to normalization with.
         */
        AwhEnergyBlock(int            numPoints,
                       Normalization  normalizationType,
                       float          normalizationValue);

        /*! \brief Returns an ArrarRef to the data in the block.
         */
        gmx::ArrayRef<float> data()
        {
            return data_;
        }

        const Normalization normalizationType;  /**< How to normalize the output data */
        const float         normalizationValue; /**< The normalization value */
    private:
        std::vector<float>  data_;              /**< The data, always float which is enough since this is statistical data */
};

/*! \internal \brief Class organizing the output data storing and writing of an AWH bias.
 */
class BiasWriter
{
    public:
        /*! \brief Constructor.
         *
         * \param[in] bias  The AWH bias.
         */
        BiasWriter(const Bias &bias);

        /*! \brief Returns the number of data blocks.
         *
         * \returns the number of data blocks.
         */
        int numBlocks() const
        {
            return block_.size();
        }

        /*! \brief Collect AWH bias data in blocks and write to energy subblocks.
         *
         * \param[in]     bias      The AWH Bias.
         * \param[in,out] subblock  Energy subblocks to write to.
         * \returns the number of blocks written.
         */
        int writeToEnergySubblocks(const Bias &bias, t_enxsubblock *subblock);

    private:
        /*! \brief Query if the writer has a block for the given variable.
         *
         * \param[in] outputType  Output entry type.
         */
        bool hasVarBlock(AwhOutputEntryType outputType) const
        {
            return outputTypeToBlock_.find(outputType)->second >= 0;
        }

        /*! \brief* Find the first block containing the given variable.
         *
         * \param[in] outputType  Output entry type.
         * \returns the first block index for the variable, or -1 there is no block.
         */
        int getVarStartBlock(AwhOutputEntryType outputType) const
        {
            return outputTypeToBlock_.find(outputType)->second;
        }

        /*! \brief Transfer AWH point data to writer data blocks.
         *
         * \param[in] metaDataIndex  Meta data block index.
         * \param[in] metaDataType   The type of meta data to write.
         * \param[in] bias           The AWH Bias.
         */
        void transferMetaDataToWriter(size_t             metaDataIndex,
                                      AwhOutputMetaData  metaDataType,
                                      const Bias        &bias);

        /*! \brief Transfer AWH point data to writer data blocks.
         *
         * \param[in] outputType  Output entry type.
         * \param[in] pointIndex  The point index.
         * \param[in] bias        The AWH Bias.
         * \param[in] pmf         PMF values.
         */
        void transferPointDataToWriter(AwhOutputEntryType          outputType,
                                       int                         pointIndex,
                                       const Bias                 &bias,
                                       gmx::ArrayRef<const float>  pmf);

        /*! \brief
         * Prepare the bias output data.
         *
         * \param[in] bias  The AWH Bias.
         */
        void prepareBiasOutput(const Bias &bias);

    private:
        std::vector<AwhEnergyBlock>       block_;             /**< The data blocks */
        std::map<AwhOutputEntryType, int> outputTypeToBlock_; /**< Start block index for each variable, -1 when variable should not be written */
};

}       // namespace gmx

#endif  /* GMX_AWH_BIASWRITER_H */
