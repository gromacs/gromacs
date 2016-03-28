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
 * This file contains function declarations needed internally by AWH
 * for preparing and writing output data to an energy frame.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_ENERGY_WRITER_H
#define GMX_AWH_ENERGY_WRITER_H

#include <vector>

#include "gromacs/fileio/enxio.h"
#include "gromacs/utility/basedefinitions.h"

class AwhBiasCollection;
struct Bias;
struct awh_params_t;
struct gmx_multisim_t;

/* TODO: the post-simulations AWH reader and this AWH writer are totally
 * disconnected although they read/write the same data. I'm not sure how
 * to handle that or if it should be left as it is until the writing is done
 * in a differen format (i.e. TNG) than the current energy file.
 */

//! Enum with the AWH variables to write
enum {
    evarMETA, evarCOORDVALUE, evarPMF, evarBIAS, evarVISITS, evarWEIGHTS,
    evarTARGET, evarNR
};

//! Enum with the types of metadata to write
enum {
    emetadataNBLOCK, emetadataTARGETERROR, emetadataSCALEDSAMPLEWEIGHT, emetadataNR
};

//! Enum with different ways of normalizing the output
enum {
    enormtypeNONE, enormtypeCOORD, enormtypeFREE_ENERGY, enormtypeDISTRIBUTION
};

//! Output data block.
struct Block
{
    int                normtype;          /**< How to normalize the output data */
    float              normvalue;         /**< The normalization value */
    std::vector<float> data;              /**< The data, always float which is enough since this is statistical data */
};

//! Class organizing the output data storing and writing of an AWH bias.
class BiasWriter
{
    std::vector<Block> block_;              /**< The data blocks */
    int                varToBlock_[evarNR]; /**< Start block index for each variable */
    bool               haveDataToWrite_;    /**< Tells if we have data to write to file */

    /*! \brief Query if the writer has a block for the given variable.
     *
     * \param[in] evar     Value for variable type enum.
     */
    bool hasVarBlock(int evar) const
    {
        return varToBlock_[evar] >= 0;
    }

    /*! \brief* Find the first block containing the given variable.
     *
     * \param[in] evar     Value for variable type enum.
     * \returns the first block index for the variable, or -1 there is no block.
     */
    int getVarStartBlock(int evar) const
    {
        return varToBlock_[evar];
    }

    /*! \brief Transfer AWH point data to writer data blocks.
     *
     * \param[in] evar        Value for variable type enum.
     * \param[in] m           Point index.
     * \param[in] bias        The AWH Bias.
     * \param[in] pmf         PMF values.
     */
    void transferVariablePointDataToWriter(int evar, int m,
                                           const Bias &bias,
                                           const std::vector<float> &pmf);

    public:
        /*! \brief Constructor.
         *
         * \param[in] bias  The AWH bias.
         */
        BiasWriter(const Bias &bias);

        /*! \brief Returns if we have data to write.
         *
         * \returns if we have data to write.
         */
        bool haveDataToWrite() const
        {
            return haveDataToWrite_;
        }

        /*! \brief Returns the number of data blocks.
         *
         * \returns the number of data blocks.
         */
        int numBlocks() const
        {
            return block_.size();
        }

        /*! \brief
         * Prepare the bias output data.
         *
         * \param[in] bias  The AWH Bias.
         * \param[in] ms    Struct for multi-simulation communication.
         */
        void prepareBiasOutput(const Bias &bias, const gmx_multisim_t *ms);

        /*! \brief Write AWH bias data blocks to energy subblocks.
         *
         * \param[in,out] sub  Energy subblock to write to.
         * \returns the number of blocks written.
         */
        int writeSubblocks(t_enxsubblock *sub);
};

#endif  /* GMX_AWH_ENERGY_WRITER_H */
