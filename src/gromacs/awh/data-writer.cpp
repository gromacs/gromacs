/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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

#include "correlation.h"
#include "grid.h"
#include "internal.h"
#include "types.h"

/* TODO: the post-simulations AWH reader and  and this AWH writer are totally disconnected although they
   read/write the same data. I'm not sure how to handle that or if it should be left as it is until the
   writing is done in a differen format (i.e. TNG) than the current energy file. */

//! Enum with the AWH variables to write
enum {
    evarMETA, evarCOORDVALUE, evarPMF, evarBIAS, evarVISITS, evarWEIGHTS,
    evarTARGET, evarFORCECORRVOL, evarFRICTION, evarCORRTIME, evarNR
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
struct block_t {
    int     npoints;                      /**< Number of data points in block */
    int     normtype;                     /**< How to normalize the output data */
    float   normvalue;                    /**< The normalization value */
    float  *data;                         /**< The data, always float which is enough since this is statistical data */
};

//! Struct organizing the output data storing and writing of an AWH bias.
struct bias_writer_t {
    int             nblock;               /**< Number of blocks */
    int             var_to_block[evarNR]; /**< Start block index for each variable */
    block_t        *block;                /**< The data blocks */
};

//!  Struct organizing AWH output data storing and writing.
struct awh_energywriter_t {
    int              nstout;               /**< Printing interval */
    int              nwriter;              /**< Number of AWH biases to write */
    bias_writer_t   *writer;               /**< The writers */
    bool             bPrint;               /**< Ready to print? */
};


/*! \brief
 * Query if the writer has a block for the given variable.
 *
 * \param[in] writer   Bias writer.
 * \param[in] evar     Value for variable type enum.
 */
static bool has_var_block(bias_writer_t *writer, int evar)
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
static int get_var_startblock(bias_writer_t *writer, int evar)
{
    return writer->var_to_block[evar];
}

/* Query if output should be written at the given step. */
bool time_to_write(gmx_int64_t step, const awh_energywriter_t *writer)
{
    return (writer->nstout > 0) && (step % writer->nstout == 0);
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
        case evarFORCECORRVOL:
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
static double get_coord_normvalue(const awh_bias_t *awh_bias, int dimIndex)
{
    /* AWH may use different units internally but here we convert to user units */
    return scaleInternalToUserInput(awh_bias, dimIndex, 1);
}

/*! \brief
 * Gets the normalization value for the given variable.
 *
 * \param[in] evar          Value for variable type enum.
 * \param[in] awh_bias      AWH bias.
 * \param[in] count         Index of the variable.
 * \returns the normalization value.
 */
static double get_normvalue(int evar, const awh_bias_t *awh_bias, int count)
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
        case evarFORCECORRVOL:
            normvalue = static_cast<float>(awh_bias->npoints);
            break;
        default:
            break;
    }

    return normvalue;
}

/*! \brief
 * Distributes allocated data to the data blocks.
 *
 * \param[in,out] writer    Data writer.
 * \param[in,out] data      Data pointer.
 */
static void set_writer_data(bias_writer_t *writer, float *data)
{
    /* Distribute the data to the blocks */
    int idata = 0;
    for (int b = 0; b < writer->nblock; b++)
    {
        writer->block[b].data = &data[idata];
        idata                += writer->block[b].npoints;
    }
}

/*! \brief
 * Initializes a data block.
 *
 * \param[in,out] block    Data block.
 * \param[in] npoints      Number of points in block.
 * \param[in] normtype     Value for normalization type enum.
 * \param[in] normvalue    Normalization value.
 */
static void init_block(block_t *block, int npoints, int normtype, double normvalue)
{
    block->data      = NULL;
    block->npoints   = npoints;
    block->normtype  = normtype;
    block->normvalue = static_cast<float>(normvalue);
}

/*! \brief
 * Initializes a bias writer.
 *
 * \param[in,out] writer   Data writer.
 * \param[in] awh_bias     AWH bias.
 */
static void init_bias_writer(bias_writer_t *writer, const awh_bias_t *awh_bias)
{
    int blockcount;
    int var_nblock[evarNR];                /* Number of blocks per variable   */

    /* Different variables need different number of blocks. We keep track of the starting block for
       each variable. */
    writer->nblock = 0;
    blockcount     = 0;
    for (int evar = 0; evar < evarNR; evar++)
    {
        switch (evar)
        {
            case evarCOORDVALUE:
                writer->var_to_block[evar] = blockcount;
                var_nblock[evar]           = awh_bias->ndim;
                break;
            case evarFORCECORRVOL:
                writer->var_to_block[evar] = blockcount;
                var_nblock[evar]           = 1;
                break;
            case evarFRICTION:
                writer->var_to_block[evar] = blockcount;
                var_nblock[evar]           = awh_bias->forcecorr->corrmatrix[0].ncorr;
                break;
            case evarCORRTIME:
                /* Currently not enabled. */
                writer->var_to_block[evar] = -1;
                var_nblock[evar]           = 0;
                break;
            default:
                /* Most variables need one block */
                writer->var_to_block[evar] = blockcount;
                var_nblock[evar]           = 1;
                break;
        }
        blockcount += var_nblock[evar];
    }
    writer->nblock = blockcount;

    /* Initialize the data blocks for each variable */
    snew(writer->block, writer->nblock);

    for (int evar = 0; evar < evarNR; evar++)
    {
        int npoints, bstart, var_blockcount;

        npoints        = (evar == evarMETA) ? emetadataNR : awh_bias->npoints;
        bstart         = writer->var_to_block[evar];
        var_blockcount = 0;
        for (int b = bstart; b < bstart + var_nblock[evar]; b++)
        {
            init_block(&writer->block[b], npoints,
                       get_normtype(evar),
                       get_normvalue(evar, awh_bias, var_blockcount));
            blockcount++;
        }
    }
}

/* Allocate, initialize and return an AWH writer. */
awh_energywriter_t *init_awh_energywriter(int nstout, const awh_t *awh)
{
    float                  *data;
    awh_energywriter_t     *energywriter;

    snew(energywriter, 1);

    energywriter->nstout = nstout;

    /* Each AWH is mapped to a writer */
    energywriter->nwriter = awh->nbias;
    snew(energywriter->writer, energywriter->nwriter);
    for (int k = 0; k < energywriter->nwriter; k++)
    {
        init_bias_writer(&energywriter->writer[k], &awh->awh_bias[k]);
    }

    /* Allocate space for flattened output data */
    int ndata_alloc = 0;
    for (int k = 0; k < energywriter->nwriter; k++)
    {
        bias_writer_t *writer = &energywriter->writer[k];

        for (int b = 0; b < writer->nblock; b++)
        {
            ndata_alloc += writer->block[b].npoints;
        }
    }
    snew(data, ndata_alloc);

    /* Give the data to the writers */
    int ndata_given = 0;
    for (int k = 0; k < energywriter->nwriter; k++)
    {
        bias_writer_t *writer = &energywriter->writer[k];

        set_writer_data(writer, &data[ndata_given]);

        for (int b = 0; b < writer->nblock; b++)
        {
            ndata_given += writer->block[b].npoints;
        }
    }
    data = NULL;

    /* No real data yet */
    energywriter->bPrint = FALSE;

    return energywriter;
}

/*! \brief
 * Normalizes data for output.
 *
 * \param[in,out] writer   Data writer.
 * \param[in] awh_bias      AWH bias.
 */
static void normalize_data(bias_writer_t *writer, const awh_bias_t *awh_bias)
{
    /* Normalize the data */
    for (int b = 0; b < writer->nblock; b++)
    {
        block_t     *block    = &writer->block[b];
        /* Here we operate on float data (which is accurate enough, since it
         * is statistical data that will never reach full float precision).
         * But since we can have very many data points, we sum into a double.
         */
        double       sum      = 0;
        float        minval   = GMX_FLOAT_MAX;
        float        inv_norm = 0;

        switch (block->normtype)
        {
            case enormtypeNONE:
                break;
            case enormtypeCOORD:
                /* Normalize coordinate values by a scale factor */
                for (int m = 0; m < block->npoints; m++)
                {
                    block->data[m] *= block->normvalue;
                }
                break;
            case enormtypeFREE_ENERGY:
                /* Normalize free energy values by subtracting the minimum value */
                for (int m = 0; m < block->npoints; m++)
                {
                    if (in_target_region(&awh_bias->coordpoint[m]) && block->data[m] < minval)
                    {
                        minval = block->data[m];
                    }
                }
                for (int m = 0; m < block->npoints; m++)
                {
                    if (in_target_region(&awh_bias->coordpoint[m]))
                    {
                        block->data[m] -= minval;
                    }
                }

                break;
            case enormtypeDISTRIBUTION:
                /* Normalize distribution value so that average equals 1. */
                for (int m = 0; m < block->npoints; m++)
                {
                    sum += block->data[m];
                }
                if (sum > 0)
                {
                    inv_norm = block->normvalue/static_cast<float>(sum);
                }
                for (int m = 0; m < block->npoints; m++)
                {
                    block->data[m] *= inv_norm;
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
static void transfer_variable_point_data_to_writer(int evar, bias_writer_t *writer, int m, const awh_bias_params_t *awh_bias_params,
                                                   const awh_bias_t *awh_bias, const float *pmf)
{
    /* All variables are generally not written */
    if (!has_var_block(writer, evar))
    {
        return;
    }

    /* The starting block index of this variable. Note that some variables need several (contiguous) blocks. */
    int      bstart       = get_var_startblock(writer, evar);
    block_t *block        = writer->block;
    GMX_ASSERT(m < block[bstart].npoints, "Attempt to transfer AWH data to block for point index out of range");

    int      ncorrelation = awh_bias->forcecorr->corrmatrix[0].ncorr;

    /* Transfer the point data of this variable to the right block(s) */
    int b = bstart;
    switch (evar)
    {
        case evarMETA:
            switch (m)
            {
                case emetadataNBLOCK:
                    /* The number of subblocks per awh (needed by gmx_energy) */
                    block[b].data[m] = static_cast<double>(writer->nblock);
                    /* Note: a single subblock takes only a single type and we need doubles. */
                    break;
                case emetadataTARGETERROR:
                    /* The theoretical target error */
                    block[b].data[m] = awh_bias_params->error_initial*std::sqrt(awh_bias->histsize_initial/awh_bias->histsize);
                    break;
                case emetadataSCALEDSAMPLEWEIGHT:
                    /* The logarithm of the sample weight relative to a sample weight of 1 at the initial time.
                       In the normal case: this will increase in the initial stage and then stay at a constant value. */
                    block[b].data[m] = awh_bias->scaledSampleWeight;
                    break;
            }
            break;
        case evarCOORDVALUE:
            for (int d = 0; d < awh_bias->ndim; d++)
            {
                block[b].data[m] = awh_bias->grid->point[m].value[d];
                b++;
            }
            break;
        case evarPMF:
            block[b].data[m] = in_target_region(&awh_bias->coordpoint[m]) ? pmf[m] : 0;
            break;
        case evarBIAS:
            block[b].data[m] = in_target_region(&awh_bias->coordpoint[m]) ? calc_convolved_bias(awh_bias, awh_bias->grid->point[m].value) : 0;
            break;
        case evarVISITS:
            block[b].data[m] = awh_bias->coordpoint[m].visits_tot;
            break;
        case evarWEIGHTS:
            block[b].data[m] = awh_bias->coordpoint[m].weightsum_tot;
            break;
        case evarTARGET:
            block[b].data[m] =  awh_bias->coordpoint[m].target;
            break;
        case evarFORCECORRVOL:
            block[b].data[m] = get_correlation_volelem(awh_bias->forcecorr, m);
            break;
        case evarFRICTION:
            /* Store force correlation in units of friction, i.e. time/length^2 */
            for (int n = 0; n < ncorrelation; n++)
            {
                block[b].data[m] = get_correlation_timeintegral(awh_bias->forcecorr, m, n);
                b++;
            }
            break;
        case evarCORRTIME:
            /* Force correlation times */
            for (int n = 0; n < ncorrelation; n++)
            {
                block[b].data[m] = get_correlation_time(awh_bias->forcecorr, m, n);
                b++;
            }
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
static void prep_awh_bias_output(bias_writer_t *writer, const awh_bias_params_t *awh_bias_params, const awh_bias_t *awh_bias, const gmx_multisim_t *ms)
{
    float *pmf;

    /* Pack the AWH data into the writer data. */

    /* Evaluate the PMF for all points */
    pmf = writer->block[writer->var_to_block[evarPMF]].data;
    getPmf(awh_bias, ms, pmf);

    /* Pack the the data point by point. */

    /* The metadata has a fixed number of points so we do this variable separately.
       We could also switch the loop order below (to variables outer, points inner).
       I'm thinking this loop order performs better but it would need testing. */
    for (int i = 0; i < emetadataNR; i++)
    {
        transfer_variable_point_data_to_writer(evarMETA, writer, i, awh_bias_params, awh_bias, pmf);
    }
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        for (int evar = 0; evar < evarNR; evar++)
        {
            if (evar == evarMETA)
            {
                continue;
            }
            transfer_variable_point_data_to_writer(evar, writer, m, awh_bias_params, awh_bias, pmf);
        }
    }

    /* For looks of the output, normalize it */
    normalize_data(writer, awh_bias);
}

/* Prepare AWH output data for writing. */
void prep_awh_output(awh_energywriter_t *energywriter, const awh_params_t *awh_params,
                     const awh_t *awh, const gmx_multisim_t *ms)
{
    for (int k = 0; k < awh->nbias; k++)
    {
        prep_awh_bias_output(&energywriter->writer[k], &awh_params->awh_bias_params[k], &awh->awh_bias[k], ms);
    }

    energywriter->bPrint = TRUE;
}

/*! \brief
 * Write AWH bias data blocks to energy subblocks.
 *
 * \param[in,out] sub       Energy subblock.
 * \param[in] writer        Bias writer.
 */
static void write_awh_to_subblocks(t_enxsubblock *sub, const bias_writer_t *writer)
{
    for (int b = 0; b < writer->nblock; b++)
    {
        sub[b].type  = xdr_datatype_float;
        sub[b].nr    = writer->block[b].npoints;
        sub[b].fval  = writer->block[b].data;
    }
}

/* Fill the AWH data block of an energy frame with data (if there is any). */
void write_awh_to_frame(t_enxframe *frame, awh_energywriter_t *energywriter)
{
    if (energywriter->bPrint)
    {
        t_enxblock *awhenergyblock;

        /* Get the total number of energy subblocks that AWH needs */
        int nsub  = 0;
        for (int k = 0; k < energywriter->nwriter; k++)
        {
            nsub  += energywriter->writer[k].nblock;
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
        for (int k = 0; k < energywriter->nwriter; k++)
        {
            write_awh_to_subblocks(&(awhenergyblock->sub[energy_subblock_count]), &energywriter->writer[k]);
            energy_subblock_count += energywriter->writer[k].nblock;
        }

        /* Don't print again until a new dataset has been prepared */
        energywriter->bPrint = FALSE;
    }
}
