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

#include "correlation.h"

#include <assert.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "math.h"

/* Get the current number of blocks. */
int get_nblocks(const correlation_grid_t *corrgrid)
{
    int                      nblockdata     = corrgrid->corrmatrix[0].corr[0].nblockdata;
    correlation_blockdata_t *blockdata      = corrgrid->corrmatrix[0].corr[0].blockdata;
    double                   maxblocklength = blockdata[nblockdata - 1].blocklength;
    double                   minblocklength = blockdata[0].blocklength;

    /* If we have a finite block span we have a constant number of blocks, otherwise we are always adding more blocks (and we don't keep track of the number) */
    if (maxblocklength < GMX_DOUBLE_MAX)
    {
        return static_cast<int>(maxblocklength/minblocklength);
    }
    else
    {
        return -1;
    }
}

/* Get the current blocklength. */
double get_blocklength(const correlation_grid_t *corrgrid)
{
    /* Return the  minimum blocklength */
    return corrgrid->corrmatrix[0].corr[0].blockdata[0].blocklength;
}

/* Returns if correlation should be sampled at the given step. */
bool time_to_sample_correlation(const correlation_grid_t *corrgrid, gmx_int64_t step)
{
    return (corrgrid->nstsample > 0) && (step % corrgrid->nstsample == 0);
}

/*! \brief
 * Get the covariance.
 *
 * \param[in] corr      Correlation element.
 * \returns the covariance.
 */
static double get_covariance(const correlation_t *corr)
{
    double weight = corr->blockdata[corr->nblockdata - 1].sum_w;

    return weight > 0 ? corr->covariance/weight : 0;
}

/* Returns an element of the time integrated correlation matrix at a given point in the grid. */
double get_correlation_timeintegral(const correlation_grid_t *corrgrid, int pointindex, int corrindex)
{
    correlation_t *corr     = &corrgrid->corrmatrix[pointindex].corr[corrindex];
    double         weight   = corr->blockdata[corr->nblockdata - 1].sum_w;
    double         corrsum  = weight > 0  ? corr->blockdata[0].correlation_integral/weight : 0;
    double         covar    = get_covariance(corr);

    return 0.5*(corrsum + covar)*corrgrid->dt_sample;
}

/* Returns an element of the correlation time matrix at a given point in the correlation grid. */
double get_correlation_time(const correlation_grid_t *corrgrid, int pointindex, int corrindex)
{
    correlation_t *corr               = &corrgrid->corrmatrix[pointindex].corr[corrindex];
    double         covar              = get_covariance(corr);
    double         corr_timeintegral  = get_correlation_timeintegral(corrgrid, pointindex, corrindex);

    return covar > 0 ? corr_timeintegral/covar : 0;
}

/* Returns the volume element of the correlation metric. */
double get_correlation_volelem(const correlation_grid_t *corrgrid, int pointindex)
{
    double det, a, b, c;

    switch (corrgrid->corrmatrix[pointindex].ncorr)
    {
        case 1:
            /* 1-dimensional matrix: [a] */
            det = get_correlation_timeintegral(corrgrid, pointindex, 0);
            break;
        case 3:
            /* 2-dimensional matrix: [a b; b c] */
            a   = get_correlation_timeintegral(corrgrid, pointindex, 0);
            b   = get_correlation_timeintegral(corrgrid, pointindex, 1);
            c   = get_correlation_timeintegral(corrgrid, pointindex, 2);

            det = a*c - b*b;
            break;
        default:
            det = 0;
            /* meh */
    }

    /* Returns 0 if no data, not supported number of dims
       or not enough data to give a positive determinant (as it should be) */
    return det > 0 ? std::sqrt(det) : 0;
}


/*! \brief
 * Updates the covariance with new data.
 *
 * \param[in,out] covariance_ptr     Covariance to set.
 * \param[in] x                      Data of variable x
 * \param[in] y                      Data of variable y.
 * \param[in] w                      Weight of data.
 * \param[in] sum_wx_old             Sum w*x, before update.
 * \param[in] sum_wy_old             Sum w*y, before update.
 * \param[in] sum_w_old              Sum w, before update.
 */
static void update_covariance(double *covariance_ptr, double x, double y, double w,
                              double sum_wx_old, double sum_wy_old, double sum_w_old)
{
    double avg_x_old, avg_y_old, inv_sum_w_old;

    if (gmx_within_tol(sum_w_old, 0, GMX_DOUBLE_EPS))
    {
        /* If this is the first sample, the covariance is 0 */
        return;
    }

    inv_sum_w_old = 1./sum_w_old;

    avg_x_old = sum_wx_old*inv_sum_w_old;
    avg_y_old = sum_wy_old*inv_sum_w_old;

    (*covariance_ptr) += w*sum_w_old/(sum_w_old + w)*(x - avg_x_old)*(y - avg_y_old);
}

/*! \brief
 * Updates the block length by doubling.
 *
 * \param[in,out] blockdata      Block data array.
 * \param[in] nblockdata         Number of block data.
 */
static void double_blocklengths(correlation_blockdata_t *blockdata, int nblockdata)
{
    /* We need to shift the data so that a given blockdata gets the data for double the block length.
       The data for the shortest block length is not needed anymore. */

    for (int i = 0; i < nblockdata - 1; i++)
    {
        blockdata[i].blocklength          = blockdata[i + 1].blocklength;
        blockdata[i].correlation_integral = blockdata[i + 1].correlation_integral;
        blockdata[i].sum_wx               = blockdata[i + 1].sum_wx;
        blockdata[i].sum_wy               = blockdata[i + 1].sum_wy;
        blockdata[i].sum_w                = blockdata[i + 1].sum_w;
    }

    /* The blockdata which has 1 block is the same as the old one but with double the block length */
    blockdata[nblockdata - 1].blocklength *= 2.;
}

/*! \brief
 * Updates the block length such that data fits.
 *
 * \param[in,out] blockdata      Block data array.
 * \param[in] nblockdata         Number of block data.
 * \param[in] sampling_length    Sampling length of all data, in time or weight.
 */
static void update_blocklengths(correlation_blockdata_t *blockdata, int nblockdata, double sampling_length)
{
    /* How many times do we need to double the longest block length to fit the data? */
    int ndoublings = ceil_log2(sampling_length/blockdata[nblockdata - 1].blocklength);

    while (ndoublings > 0)
    {
        double_blocklengths(blockdata, nblockdata);
        ndoublings--;
    }
}

/*! \brief
 * Adds filled data block to correlation time integral.
 *
 * \param[in,out] blockdata      Block data array.
 * \param[in] tot_sum_wx         Sum w*x, all data.
 * \param[in] tot_sum_wy         Sum w*y, all data.
 * \param[in] tot_sum_w          Sum w, all data.
 */
static void add_block_to_correlation_integral(correlation_blockdata_t *blockdata,
                                              double tot_sum_wx, double tot_sum_wy, double tot_sum_w)
{
    double sum_dx_old, sum_dy_new, avg_x_old, avg_y_new;

    /* Need the old average, before the data of this block was added */
    if (!gmx_within_tol(tot_sum_w - blockdata->sum_w, 0, GMX_DOUBLE_EPS))
    {
        avg_x_old = (tot_sum_wx - blockdata->sum_wx)/(tot_sum_w - blockdata->sum_w);
    }
    else
    {
        /* This is the first block added */
        avg_x_old = 0;
    }

    avg_y_new = tot_sum_wy/tot_sum_w;

    sum_dx_old = blockdata->sum_wx - avg_x_old*blockdata->sum_w;
    sum_dy_new = blockdata->sum_wy - avg_y_new*blockdata->sum_w;

    blockdata->correlation_integral += sum_dx_old*sum_dy_new;

    /* Reset */
    blockdata->sum_wx = 0.;
    blockdata->sum_wy = 0.;
    blockdata->sum_w  = 0.;
}

/*! \brief
 * Adds filled data block to correlation time integral.
 *
 * \param[in] blocklength  Block length.
 * \param[in] length       Sampling length of all data, in time or weight.
 */
static int get_block_index(double blocklength, double length)
{
    return static_cast<int>(length/blocklength);
}

/*! \brief
 * Adds new data to block.
 *
 * \param[in,out] blockdata     Block data structure.
 * \param[in] x                 Data of variable x.
 * \param[in] y                 Data of variable y
 * \param[in] w                 Weight of data.
 */
static void add_data_to_block(correlation_blockdata_t *blockdata, double x, double y, double w)
{
    blockdata->sum_wx   += w*x;
    blockdata->sum_wy   += w*y;
    blockdata->sum_w    += w;
}

/*! \brief
 * Update correlation element with new data.
 *
 * \param[in,out] corr          Correlation element.
 * \param[in] x                 Data of variable x.
 * \param[in] y                 Data of variable y
 * \param[in] w                 Weight of data.
 * \param[in] sampling_length   Sampling length of all data, in time or weight.
 */
static void update_correlation(correlation_t *corr, double x, double y, double w, double sampling_length)
{
    correlation_blockdata_t *all_data;

    if (gmx_within_tol(w, 0, GMX_DOUBLE_EPS))
    {
        /* Nothing to add */
        return;
    }

    /* The last blockdata has only 1 block which contains all data so far */
    all_data = &corr->blockdata[corr->nblockdata - 1];

    /* The covariance is updated for every sample using the old, non-updated averages. */
    update_covariance(&corr->covariance, x, y, w, all_data->sum_wx, all_data->sum_wy, all_data->sum_w);

    /* Make sure the blocks are long enough to fit all data */
    update_blocklengths(corr->blockdata, corr->nblockdata, sampling_length);

    /* Store the data for each block length considered. First update the longest block which has all data since it's
       used for updating the correlation function for the other block lengths. */
    add_data_to_block(all_data, x, y, w);

    for (int i = 0; i < corr->nblockdata - 1; i++)
    {
        /* Find current block index given block length. */
        int block_index = get_block_index(corr->blockdata[i].blocklength, sampling_length);

        if (corr->blockdata[i].blockindex_prev >= 0 && corr->blockdata[i].blockindex_prev != block_index)
        {
            /* Changed block. Update correlation with old data before adding to new block. */
            add_block_to_correlation_integral(&corr->blockdata[i], all_data->sum_wx, all_data->sum_wy, all_data->sum_w);
        }

        /* Keep track of which block index data was last added to */
        corr->blockdata[i].blockindex_prev = block_index;

        /* Store the data */
        add_data_to_block(&corr->blockdata[i], x, y, w);
    }
}

/*! \brief
 * Get the total weight of the data in the correlation matrix.
 *
 * \param[in] corrmatrix    Correlation matrix.
 * \returns the weight of the added data.
 */
static double get_weight(const correlation_matrix_t *corrmatrix)
{
    /* All elements have the same weight */
    correlation_t *corr = &corrmatrix->corr[0];

    /* The last blockdata has only 1 block containing all data */
    return corr->blockdata[corr->nblockdata - 1].sum_w;
}

/* Adds a weighted data vector to one point in the correlation grid. */
void add_data_to_correlation_matrix(correlation_grid_t *corrgrid, int pointindex,
                                    double weight, const double *dimdata,
                                    double t)
{
    double                sampling_length;
    correlation_matrix_t *corrmatrix = &corrgrid->corrmatrix[pointindex];

    /*  The sampling length is measured either in the total (local) weight or the current time */
    sampling_length = corrgrid->bBlocklength_in_weight ? get_weight(corrmatrix) + weight : t;

    /* Add data for each correlation element in the matrix */
    for (int i = 0; i < corrmatrix->ncorr; i++)
    {
        correlation_t *corr = &corrmatrix->corr[i];
        update_correlation(corr, dimdata[corr->xy_index[0]], dimdata[corr->xy_index[1]],
                           weight, sampling_length);
    }
}

/*! \brief
 * Initialize block data.
 *
 * \param[in] b              Block data to initialize.
 * \param[in] blocklength    Block length.
 */
static void init_correlation_blockdata_t(correlation_blockdata_t *b, double blocklength)
{
    b->sum_wx  = 0;
    b->sum_wy  = 0;
    b->sum_w   = 0;

    b->blockindex_prev     = -1; /* no data added */
    b->blocklength         = blocklength;
}

/*! \brief
 * Initialize correlation element.
 *
 * \param[in,out] corr          Correlation matrix element.
 * \param[in] nblockdata        Number of block data structs.
 * \param[in,out] blockdata     Block data array.
 * \param[in] i1                First dimensional index of this matrix element.
 * \param[in] i2                Second dimensional index of this matrix element.
 * \param[in] blocklength       Block length.
 * \param[in] bConst_nblocks    True if a constant number of blocks should be used, otherwise increase when needed.
 */
static void init_correlation(correlation_t *corr, int nblockdata, correlation_blockdata_t *blockdata,
                             int i1, int i2,
                             double blocklength, bool bConst_nblocks)
{
    corr->xy_index[0]          = i1;
    corr->xy_index[1]          = i2;
    corr->covariance           = 0;
    corr->nblockdata           = nblockdata;
    corr->blockdata            = blockdata;

    if (bConst_nblocks)
    {
        int scaling = 1;

        for (int n = 0; n < corr->nblockdata; n++)
        {
            init_correlation_blockdata_t(&corr->blockdata[n], blocklength*scaling);
            scaling <<= 1; /* Double block length */
        }
    }
    else
    {
        /* Should have 2 different block lengths */
        init_correlation_blockdata_t(&corr->blockdata[0], blocklength);
        init_correlation_blockdata_t(&corr->blockdata[1], GMX_DOUBLE_MAX);
    }
}

/*! \brief
 * Return the number of block data structs needed for keeping a certain number of blocks.
 *
 * \param[in] nblocks  Number of blocks.
 * \returns the number of block data structs.
 */
static int get_nblockdata(int nblocks)
{
    int nblockdata = 0;

    while (nblocks > 0)
    {
        nblocks >>= 1; /* divide by 2 */
        nblockdata++;
    }
    return nblockdata;
}

/* Allocate, initialize and return a correlation grid struct. */
correlation_grid_t *init_correlation_grid(int npoints, int ndim,
                                          int nblocks, double blocklength,
                                          bool bBlocklength_in_weight,
                                          double dt_step,
                                          int nstsample)
{
    int                      nblockdata, nelem;
    bool                     bConst_nblocks;
    correlation_grid_t      *corrgrid       = NULL;
    correlation_t           *corr_temp      = NULL;
    correlation_blockdata_t *blockdata_temp = NULL;

    snew(corrgrid, 1);

    corrgrid->nstsample               = nstsample;
    corrgrid->dt_sample               = nstsample*dt_step;
    corrgrid->bBlocklength_in_weight  = bBlocklength_in_weight;

    /* Set the initial block length for the block averaging. The length doesn't really matter
       after the block length has been doubled a few times, as long as it's set small enough */
    if (corrgrid->bBlocklength_in_weight)
    {
        blocklength = blocklength > 0 ? blocklength : 1;
    }
    else
    {
        blocklength = blocklength > 0 ? blocklength : dt_step*nstsample;
    }

    /* Set the number of blocks. The number of blocks determines the current span of the data
       and how many different block lengths (nblockdata) we need to keep track of to be able to
       increase the block length later */
    if (nblocks == 0)
    {
        /* Default value */
        nblocks = 128;
    }
    if (nblocks > 0)
    {
        /* User value */
        nblockdata     = get_nblockdata(nblocks);
        bConst_nblocks = TRUE;
    }
    else
    {
        /* Do not keep a constant number of blocks and do not double the block length,
           instead just keep adding new blocks. We keep data of 2 different block lengths,
           the current one and the trivial 1-block partition with all the data. */
        nblockdata     = 2;
        bConst_nblocks = FALSE;
    }

    corrgrid->ncorrmatrix = npoints;
    snew(corrgrid->corrmatrix, corrgrid->ncorrmatrix);

    /* The correlation matrix is symmetric correlation. Only half of it needs to be stored. */
    nelem = ndim*(ndim + 1)/2;
    /* Allocate memory for the correlation matrix elements in one chunk and distribute */
    snew(corr_temp, corrgrid->ncorrmatrix*nelem);
    int corr_index = 0;
    for (int m = 0; m < corrgrid->ncorrmatrix; m++)
    {
        /* Initialize the correlation matrix */
        corrgrid->corrmatrix[m].ncorr = nelem;
        corrgrid->corrmatrix[m].corr  = &(corr_temp[corr_index]);
        corr_index                   += nelem;
    }
    corr_temp = NULL;

    /* Initialize the correlation elements. Allocate blockdata memory in one chunk and
       distribute it */
    snew(blockdata_temp, corrgrid->ncorrmatrix*nelem*nblockdata);
    int blockdata_index = 0;
    for (int m = 0; m < corrgrid->ncorrmatrix; m++)
    {
        int d = 0;
        for (int d1 = 0; d1 < ndim; d1++)
        {
            for (int d2 = 0; d2 <= d1; d2++)
            {
                init_correlation(&corrgrid->corrmatrix[m].corr[d],
                                 nblockdata, &(blockdata_temp[blockdata_index]),
                                 d1, d2,
                                 blocklength, bConst_nblocks);
                d++;
                blockdata_index += nblockdata;
            }
        }
    }
    blockdata_temp = NULL;

    return corrgrid;
}
