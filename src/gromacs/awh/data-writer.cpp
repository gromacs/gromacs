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

#include "data-writer.h"

#include <assert.h>

#include <cmath>

#include "gromacs/fileio/enxio.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "grid.h"
#include "internal.h"
#include "types.h"


/*! \brief
 * Query if the writer has a block for the given variable.
 *
 * \param[in] writer   Bias writer.
 * \param[in] evar     Value for variable type enum.
 */
static bool has_var_block(BiasWriter *writer, int evar)
{
    return writer->var_to_block[evar] >= 0;
}

/*! \brief
 * Find the first block containing the given variable.
 *
 * \param[in] writer   Bias writer.
 * \param[in] evar     Value for variable type enum.
 * \returns the first block index for the variable, or -1 there is no block.
 */
static int get_var_startblock(BiasWriter *writer, int evar)
{
    return writer->var_to_block[evar];
}

/*! \brief
 * Map the variable type to a normalization type.
 *
 * \param[in] evar     Value for variable type enum.
 * \returns enum value for normalization type.
 */
static int get_normtype(int evar)
{
    int normtype;

    switch (evar)
    {
        case evarCOORDVALUE:
            normtype = enormtypeCOORD; break;
        case evarPMF:
        case evarBIAS:
            normtype = enormtypeFREE_ENERGY; break;
        case evarVISITS:
        case evarWEIGHTS:
        case evarTARGET:
            normtype = enormtypeDISTRIBUTION; break;
        default:
            normtype = enormtypeNONE; break;
    }

    return normtype;
}

/*! \brief
 * Gets the coordinate normalization value for the given dimension.
 *
 * \param[in] awh_bias       AWH bias.
 * \param[in] dimIndex       Dimensional index.
 * \returns the coordinate normalization value.
 */
static double get_coord_normvalue(const AwhBias &awh_bias, int dimIndex)
{
    /* AWH may use different units internally but here we convert to user units */
    return awh_bias.dimParams[dimIndex].scaleInternalToUserInput(1);
}

/*! \brief
 * Gets the normalization value for the given variable.
 *
 * \param[in] evar          Value for variable type enum.
 * \param[in] awh_bias      AWH bias.
 * \param[in] count         Index of the variable.
 * \returns the normalization value.
 */
static double get_normvalue(int evar, const AwhBias &awh_bias, int count)
{
    double normvalue = 0;

    switch (evar)
    {
        case evarCOORDVALUE:
            normvalue = get_coord_normvalue(awh_bias, count);
            break;
        case evarVISITS:
        case evarWEIGHTS:
        case evarTARGET:
            normvalue = static_cast<double>(awh_bias.coordPoint.size());
            break;
        default:
            break;
    }

    return normvalue;
}

/*! \brief
 * Initializes a data block.
 *
 * \param[in,out] block    Data block.
 * \param[in] npoints      Number of points in block.
 * \param[in] normtype     Value for normalization type enum.
 * \param[in] normvalue    Normalization value.
 */
static void init_block(Block *block, int npoints, int normtype, double normvalue)
{
    block->data.resize(npoints);
    block->normtype  = normtype;
    block->normvalue = static_cast<float>(normvalue);
}

BiasWriter *init_bias_writer(const AwhBias &awh_bias)
{
    BiasWriter *writer = new BiasWriter;

    int var_nblock[evarNR];                /* Number of blocks per variable   */

    /* Different variables need different number of blocks. We keep track of the starting block for
       each variable. */
    int blockcount = 0;
    for (int evar = 0; evar < evarNR; evar++)
    {
        switch (evar)
        {
            case evarCOORDVALUE:
                writer->var_to_block[evar] = blockcount;
                var_nblock[evar]           = awh_bias.ndim;
                break;
            default:
                /* Most variables need one block */
                writer->var_to_block[evar] = blockcount;
                var_nblock[evar]           = 1;
                break;
        }
        blockcount += var_nblock[evar];
    }

    /* Initialize the data blocks for each variable */
    writer->block.resize(blockcount);

    for (int evar = 0; evar < evarNR; evar++)
    {
        int npoints;
        if (evar == evarMETA)
        {
            npoints        = emetadataNR;
        }
        else
        {
            npoints        = awh_bias.coordPoint.size();
        }
        int bstart         = writer->var_to_block[evar];
        int var_blockcount = 0;
        for (int b = bstart; b < bstart + var_nblock[evar]; b++)
        {
            init_block(&writer->block[b], npoints,
                       get_normtype(evar),
                       get_normvalue(evar, awh_bias, var_blockcount));
            blockcount++;
        }
    }

     /* No real data yet */
    writer->bPrint = FALSE;

    return writer;
}

/*! \brief
 * Normalizes data for output.
 *
 * \param[in,out] writer   Data writer.
 * \param[in] awh_bias      AWH bias.
 */
static void normalize_data(BiasWriter *writer, const AwhBias &awh_bias)
{
    /* Normalize the data */
    for (auto &block : writer->block)
    {
        /* Here we operate on float data (which is accurate enough, since it
         * is statistical data that will never reach full float precision).
         * But since we can have very many data points, we sum into a double.
         */
        double       sum      = 0;
        float        minval   = GMX_FLOAT_MAX;
        float        inv_norm = 0;

        switch (block.normtype)
        {
            case enormtypeNONE:
                break;
            case enormtypeCOORD:
                /* Normalize coordinate values by a scale factor */
                for (auto &point : block.data)
                {
                    point *= block.normvalue;
                }
                break;
            case enormtypeFREE_ENERGY:
                /* Normalize free energy values by subtracting the minimum value */
                for (size_t m = 0; m < block.data.size(); m++)
                {
                    if (awh_bias.coordPoint[m].inTargetRegion() && block.data[m] < minval)
                    {
                        minval = block.data[m];
                    }
                }
                for (size_t m = 0; m < block.data.size(); m++)
                {
                    if (awh_bias.coordPoint[m].inTargetRegion())
                    {
                        block.data[m] -= minval;
                    }
                }

                break;
            case enormtypeDISTRIBUTION:
                /* Normalize distribution values by normalizing their sum */
                for (auto &point : block.data)
                {
                    sum += point;
                }
                if (sum > 0)
                {
                    inv_norm = block.normvalue/static_cast<float>(sum);
                }
                for (auto &point : block.data)
                {
                    point *= inv_norm;
                }
                break;
            default:
                gmx_incons("Unknown enormtype");
        }
    }
}

/*! \brief
 * Transfer AWH point data to writer data blocks.
 *
 * \param[in] evar             Value for variable type enum.
 * \param[in,out] writer       Data writer.
 * \param[in] m                Point index.
 * \param[in] awh_bias_params  Bias parameters.
 * \param[in] awh_bias         Bias.
 * \param[in] pmf              PMF values.
 */
static void transfer_variable_point_data_to_writer(int evar, BiasWriter *writer, int m, const awh_bias_params_t *awh_bias_params,
                                                   const AwhBias &awh_bias, const std::vector<float> &pmf)
{
    /* All variables are generally not written */
    if (!has_var_block(writer, evar))
    {
        return;
    }

    /* The starting block index of this variable. Note that some variables need several (contiguous) blocks. */
    int                 bstart = get_var_startblock(writer, evar);
    std::vector<Block> &block  = writer->block;
    GMX_ASSERT(m < static_cast<int>(block[bstart].data.size()), "Attempt to transfer AWH data to block for point index out of range");

    /* Transfer the point data of this variable to the right block(s) */
    int b = bstart;
    switch (evar)
    {
        case evarMETA:
            switch (m)
            {
                case emetadataNBLOCK:
                    /* The number of subblocks per awh (needed by gmx_energy) */
                    block[b].data[m] = static_cast<double>(writer->block.size());
                    /* Note: a single subblock takes only a single type and we need doubles. */
                    break;
                case emetadataTARGETERROR:
                    /* The theoretical target error */
                    block[b].data[m] = awh_bias_params->error_initial*std::sqrt(awh_bias.params.histsize_initial/awh_bias.state.histsize);
                    break;
                case emetadataSCALEDSAMPLEWEIGHT:
                    /* The logarithm of the sample weight relative to a sample weight of 1 at the initial time.
                       In the normal case: this will increase in the initial stage and then stay at a constant value. */
                    block[b].data[m] = awh_bias.state.scaledSampleWeight;
                    break;
            }
            break;
        case evarCOORDVALUE:
            for (int d = 0; d < awh_bias.ndim; d++)
            {
                block[b].data[m] = awh_bias.grid->point[m].value[d];
                b++;
            }
            break;
        case evarPMF:
            block[b].data[m] = awh_bias.coordPoint[m].inTargetRegion() ? pmf[m] : 0;
            break;
        case evarBIAS:
            block[b].data[m] = awh_bias.coordPoint[m].inTargetRegion() ? calc_convolved_bias(awh_bias, awh_bias.grid->point[m].value) : 0;
            break;
        case evarVISITS:
            block[b].data[m] = awh_bias.coordPoint[m].visits_tot;
            break;
        case evarWEIGHTS:
            block[b].data[m] = awh_bias.coordPoint[m].weightsum_tot;
            break;
        case evarTARGET:
            block[b].data[m] =  awh_bias.coordPoint[m].target;
            break;
        default:
            gmx_incons("Unknown AWH output variable");
            break;
    }
}


/*! \brief
 * Prepare the bias output data.
 *
 * \param[in,out] writer       Bias writer.
 * \param[in] awh_bias_params  Bias parameters.
 * \param[in] awh_bias         Bias.
 * \param[in] ms               Struct for multi-simulation communication.
 */
static void prep_awh_bias_output(BiasWriter *writer, const awh_bias_params_t *awh_bias_params, const AwhBias &awh_bias, const gmx_multisim_t *ms)
{
    /* Pack the AWH data into the writer data. */

    /* Evaluate the PMF for all points */
    std::vector<float> *pmf = &writer->block[writer->var_to_block[evarPMF]].data;
    getPmf(awh_bias, ms, pmf);

    /* Pack the the data point by point. */

    /* The metadata has a fixed number of points so we do this variable separately.
       We could also switch the loop order below (to variables outer, points inner).
       I'm thinking this loop order performs better but it would need testing. */
    for (int i = 0; i < emetadataNR; i++)
    {
        transfer_variable_point_data_to_writer(evarMETA, writer, i, awh_bias_params, awh_bias, *pmf);
    }
    for (size_t m = 0; m < awh_bias.coordPoint.size(); m++)
    {
        for (int evar = 0; evar < evarNR; evar++)
        {
            if (evar == evarMETA)
            {
                continue;
            }
            transfer_variable_point_data_to_writer(evar, writer, m, awh_bias_params, awh_bias, *pmf);
        }
    }

    /* For looks of the output, normalize it */
    normalize_data(writer, awh_bias);
}

/* Prepare AWH output data for writing. */
void prep_awh_output(AwhBiasCollection    *awh,
                     const awh_params_t   *awh_params,
                     const gmx_multisim_t *ms)
{
    for (size_t k = 0; k < awh->awhBias.size(); k++)
    {
        AwhBias &awhBias = awh->awhBias[k];
        prep_awh_bias_output(awhBias.writer.get(), &awh_params->awh_bias_params[k], awhBias, ms);
        awhBias.writer->bPrint = TRUE;
    }
}

/*! \brief
 * Write AWH bias data blocks to energy subblocks.
 *
 * \param[in,out] sub       Energy subblock.
 * \param[in] writer        Bias writer.
 */
static void write_awh_to_subblocks(t_enxsubblock *sub, BiasWriter *writer)
{
    for (size_t b = 0; b < writer->block.size(); b++)
    {
        sub[b].type  = xdr_datatype_float;
        sub[b].nr    = writer->block[b].data.size();
        sub[b].fval  = writer->block[b].data.data();
    }
}

/* Fill the AWH data block of an energy frame with data (if there is any). */
void write_awh_to_frame(t_enxframe *frame, AwhBiasCollection *awh)
{
    if (awh->awhBias[0].writer->bPrint)
    {
        t_enxblock *awhenergyblock;

        /* Get the total number of energy subblocks that AWH needs */
        int nsub  = 0;
        for (auto &awhBias : awh->awhBias)
        {
            nsub  += awhBias.writer->block.size();
        }

        /* Add 1 energy block */
        add_blocks_enxframe(frame, frame->nblock + 1);

        /* Take the block that was just added and set the number of subblocks. */
        awhenergyblock = &(frame->block[frame->nblock - 1]);
        add_subblocks_enxblock(awhenergyblock, nsub);

        /* Claim it as an AWH block. */
        awhenergyblock->id = enxAWH;

        /* Transfer AWH data blocks to energy sub blocks */
        int energy_subblock_count = 0;
        for (auto &awhBias: awh->awhBias)
        {
            BiasWriter *writer = awhBias.writer.get();
            write_awh_to_subblocks(&(awhenergyblock->sub[energy_subblock_count]), writer);
            energy_subblock_count += writer->block.size();

            GMX_ASSERT(writer->bPrint, "All writers should have the same bPrint value");
            /* Don't print again until a new dataset has been prepared */
            writer->bPrint = FALSE;
        }
    }
}
