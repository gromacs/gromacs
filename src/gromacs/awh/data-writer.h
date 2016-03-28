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

struct AwhBiasCollection;
struct AwhBias;
struct awh_params_t;
struct gmx_multisim_t;

/* TODO: the post-simulations AWH reader and  and this AWH writer are totally disconnected although they
   read/write the same data. I'm not sure how to handle that or if it should be left as it is until the
   writing is done in a differen format (i.e. TNG) than the current energy file. */

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

//! Struct organizing the output data storing and writing of an AWH bias.
struct BiasWriter
{
    std::vector<Block> block;                /**< The data blocks */
    int                var_to_block[evarNR]; /**< Start block index for each variable */
    bool               bPrint;               /**< Ready to print? */
};

/*! \brief Initializes and returns a bias writer.
 *
 * \param[in] awh_bias     AWH bias.
 * \returns a bias writer.
 */
BiasWriter *init_bias_writer(const AwhBias &awh_bias);

/*! \brief Prepare AWH output data for writing.
 *
 * \param[in,out] awh        AWH working struct.
 * \param[in]     awh_params AWH parameters.
 * \param[in]     ms         Struct for multi-simulation communication, needed for bias sharing replicas.
 */
void prep_awh_output(AwhBiasCollection *awh, const awh_params_t *awh_params,
                     const gmx_multisim_t *ms);

/*! \brief Fill the AWH data block of an energy frame with data (if there is any).
 *
 * \param[in,out] frame     Energy data frame.
 * \param[in,out] awh       AWH working struct.
 */
void write_awh_to_frame(t_enxframe *frame, AwhBiasCollection *awh);

#endif  /* GMX_AWH_ENERGY_WRITER_H */
