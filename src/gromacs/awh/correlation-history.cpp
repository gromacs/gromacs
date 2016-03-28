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

#include "correlation-history.h"

#include <assert.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/awh-correlation-history.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "correlation.h"

/*! \brief
 * Unpacks the correlation history elements into a linear array.
 *
 * \param[in,out] corr_flat      Flattened array of correlation elements.
 * \param[in] corrgrid_hist      Correlation grid with correlation elements to unpack.
 */
static void unpack_correlation_history(correlation_history_t *corr_flat, const correlation_grid_history_t *corrgrid_hist)
{
    int k = 0;
    for (int i = 0; i < corrgrid_hist->ncorrmatrix; i++)
    {
        for (int j = 0; j < corrgrid_hist->ncorr; j++)
        {
            corr_flat[k] = corrgrid_hist->corrmatrix[i].corr[j];
            k++;
        }
    }
}

/*! \brief
 * Unpacks the block data history elements into a linear array.
 *
 * \param[in,out] blockdata_flat      Flattened array of block data structs.
 * \param[in] corrgrid_hist           Correlation grid with block data structs to unpack.
 */
static void unpack_blockdata_history(correlation_blockdata_history_t *blockdata_flat, const correlation_grid_history_t *corrgrid_hist)
{
    int l = 0;
    for (int i = 0; i < corrgrid_hist->ncorrmatrix; i++)
    {
        for (int j = 0; j < corrgrid_hist->ncorr; j++)
        {
            for (int k = 0; k < corrgrid_hist->nblockdata; k++)
            {
                blockdata_flat[l] = corrgrid_hist->corrmatrix[i].corr[j].blockdata[k];
                l++;
            }
        }
    }
}

/*! \brief
 * Packs the correlation history elements from a linear array into the correlation grid struct.
 *
 * \param[in,out] corr_flat           Flattened array of correlation elements.
 * \param[in] corrgrid_hist           Correlation grid to pack correlation elements into.
 */
static void pack_correlation_history(correlation_history_t *corr_flat, correlation_grid_history_t *corrgrid_hist)
{
    int j = 0;
    for (int i  = 0; i < corrgrid_hist->ncorrmatrix; i++)
    {
        corrgrid_hist->corrmatrix[i].corr = &(corr_flat[j]);
        j += corrgrid_hist->ncorr;
    }
}

/*! \brief
 * Packs the block data history structs from a linear array into the correlation grid struct.
 *
 * \param[in,out] blockdata_flat      Flattened array of block data structs.
 * \param[in] corrgrid_hist           Correlation grid to pack block data structs into.
 */
static void pack_blockdata_history(correlation_blockdata_history_t *blockdata_flat, correlation_grid_history_t *corrgrid_hist)
{
    int k = 0;
    for (int i = 0; i < corrgrid_hist->ncorrmatrix; i++)
    {
        for (int j = 0; j < corrgrid_hist->ncorr; j++)
        {
            corrgrid_hist->corrmatrix[i].corr[j].blockdata = &(blockdata_flat[k]);
            k += corrgrid_hist->nblockdata;
        }
    }
}

/*! \brief
 * Initialize and allocate for the correlation grid history when only the master rank has been initialized.
 *
 * \param[in] corrgrid_hist      Correlation grid to pack block data structs into.
 * \param[in] cr                 Struct for communication.
 */
static void broadcast_initialize_correlation_grid_history(correlation_grid_history_t *corrgrid_hist, const t_commrec *cr)
{
    int ncorr_flat, nblockdata_flat;
    correlation_history_t           *corr_flat      = NULL;
    correlation_blockdata_history_t *blockdata_flat = NULL;

    /* Directly broadcast and allocate memory for the top levels of the data structure */
    gmx_bcast(sizeof(correlation_grid_history_t), corrgrid_hist, cr);
    if (!MASTER(cr))
    {
        snew(corrgrid_hist->corrmatrix, corrgrid_hist->ncorrmatrix);
    }
    gmx_bcast(sizeof(correlation_matrix_history_t)*corrgrid_hist->ncorrmatrix, corrgrid_hist->corrmatrix, cr);

    /* To avoid broadcasting each lower level of the data structure, first flatten the data, then broadcast, then pack the data. */
    ncorr_flat      = corrgrid_hist->ncorrmatrix*corrgrid_hist->ncorr;
    nblockdata_flat = corrgrid_hist->ncorrmatrix*corrgrid_hist->ncorr*corrgrid_hist->nblockdata;
    snew(corr_flat, ncorr_flat);
    snew(blockdata_flat, nblockdata_flat);

    if (MASTER(cr))
    {
        /* Unpack data of master to linear array */
        unpack_correlation_history(corr_flat, corrgrid_hist);
        unpack_blockdata_history(blockdata_flat, corrgrid_hist);
    }

    /* Broadcast */
    gmx_bcast(sizeof(correlation_history_t)*ncorr_flat, corr_flat, cr);
    gmx_bcast(sizeof(correlation_blockdata_history_t)*nblockdata_flat, blockdata_flat, cr);

    if (MASTER(cr))
    {
        sfree(blockdata_flat);
        sfree(corr_flat);
    }
    else
    {
        /* Pack data */
        pack_correlation_history(corr_flat, corrgrid_hist);
        pack_blockdata_history(blockdata_flat, corrgrid_hist);
        corr_flat      = NULL;
        blockdata_flat = NULL;
    }
}

/*  Allocate a correlation grid history when restarting from a checkpoint. */
correlation_grid_history_t *init_correlation_grid_history_from_checkpoint(correlation_grid_history_t *corrgrid_hist_in, const t_commrec *cr)
{

    correlation_grid_history_t *corrgrid_hist = NULL;

    /* Only master has read the checkpoint and only master has an initialized the history.
       Non-master ranks need to get the their history allocated and initialized by broadcasting from the master rank. */
    if (PAR(cr))
    {
        GMX_RELEASE_ASSERT((MASTER(cr) && corrgrid_hist_in != NULL) ||
                           (!MASTER(cr) && corrgrid_hist_in == NULL),
                           "(Non-)master rank(s) should (not) have been initialized before restoring AWH correlation history");

        if (!MASTER(cr))
        {
            snew(corrgrid_hist, 1);
        }
        else
        {
            corrgrid_hist = corrgrid_hist_in;
        }

        broadcast_initialize_correlation_grid_history(corrgrid_hist, cr);
    }

    return !MASTER(cr) ? corrgrid_hist : NULL;
}

/*! \brief
 * Allocate and initialize correlation grid history.
 *
 * \param[in,out] corrgrid_hist  Correlation grid history for master rank.
 * \param[in] ncorrmatrix        Number of correlation matrices in grid.
 * \param[in] ncorr              Number of correlation elements in each matrix.
 * \param[in] nblockdata         Number of block data structs needed for each correlation element.
 */
void init_correlation_grid_history(correlation_grid_history_t *corrgrid_hist, int ncorrmatrix, int ncorr, int nblockdata)
{
    correlation_history_t           *corr_flat      = NULL;
    correlation_blockdata_history_t *blockdata_flat = NULL;

    corrgrid_hist->ncorrmatrix = ncorrmatrix;
    corrgrid_hist->ncorr       = ncorr;
    corrgrid_hist->nblockdata  = nblockdata;

    snew(corrgrid_hist->corrmatrix, ncorrmatrix);

    /* Allocate memory for correlation elements and block data in one block and distribute */
    snew(corr_flat, ncorrmatrix*ncorr);
    snew(blockdata_flat, ncorrmatrix*ncorr*nblockdata);
    pack_correlation_history(corr_flat, corrgrid_hist);
    pack_blockdata_history(blockdata_flat, corrgrid_hist);

    corr_flat      = NULL;
    blockdata_flat = NULL;
}

/* Allocate a correlation grid history when restarting from a checkpoint. */
correlation_grid_history_t *init_correlation_grid_history_from_state(const correlation_grid_t *corrgrid)
{
    correlation_grid_history_t *corrgrid_hist;

    snew(corrgrid_hist, 1);

    int ncorrmatrix = corrgrid->ncorrmatrix;
    int ncorr       = corrgrid->corrmatrix[0].ncorr;
    int nblockdata  = corrgrid->corrmatrix[0].corr[0].nblockdata;

    init_correlation_grid_history(corrgrid_hist, ncorrmatrix, ncorr, nblockdata);

    return corrgrid_hist;
}

/* Update the correlation grid history for checkpointing. */
void update_correlation_grid_history(correlation_grid_history_t *corrgrid_hist, const correlation_grid_t *corrgrid)
{
    GMX_RELEASE_ASSERT(corrgrid_hist->corrmatrix != NULL, "AWH force correlation matrix not initialized when updating history");

    /* Matrix for each grid point */
    for (int m = 0; m < corrgrid->ncorrmatrix; m++)
    {
        correlation_matrix_history_t *corrmatrix_hist = &corrgrid_hist->corrmatrix[m];
        correlation_matrix_t         *corrmatrix      = &corrgrid->corrmatrix[m];

        /* Correlation elements for each matrix */
        for (int k = 0; k < corrmatrix->ncorr; k++)
        {
            correlation_history_t *corr_hist = &corrmatrix_hist->corr[k];
            correlation_t         *corr      = &corrmatrix->corr[k];

            corr_hist->covariance = corr->covariance;

            /* Blockdatas for each correlation element */
            for (int l = 0; l < corr->nblockdata; l++)
            {
                correlation_blockdata_history_t *blockdata_hist = &corr_hist->blockdata[l];
                correlation_blockdata_t         *blockdata      = &corr->blockdata[l];

                blockdata_hist->blocklength          = blockdata->blocklength;
                blockdata_hist->correlation_integral = blockdata->correlation_integral;
                blockdata_hist->sum_wx               = blockdata->sum_wx;
                blockdata_hist->sum_wy               = blockdata->sum_wy;
                blockdata_hist->sum_w                = blockdata->sum_w;
                blockdata_hist->blockindex_prev      = blockdata->blockindex_prev;
            }
        }
    }
}

/* Restores the correlation grid state from the correlation grid history. */
void restore_correlation_grid_state_from_history(const correlation_grid_history_t *corrgrid_hist, correlation_grid_t *corrgrid)
{
    /* Matrix for each grid point */
    for (int m = 0; m < corrgrid->ncorrmatrix; m++)
    {
        correlation_matrix_history_t *corrmatrix_hist = &corrgrid_hist->corrmatrix[m];
        correlation_matrix_t         *corrmatrix      = &corrgrid->corrmatrix[m];

        /* Correlation elements for each matrix */
        for (int k = 0; k < corrmatrix->ncorr; k++)
        {
            correlation_history_t *corr_hist = &corrmatrix_hist->corr[k];
            correlation_t         *corr      = &corrmatrix->corr[k];

            corr->covariance = corr_hist->covariance;

            /* Blockdatas for each correlation element */
            for (int l = 0; l < corr->nblockdata; l++)
            {
                correlation_blockdata_history_t *blockdata_hist = &corr_hist->blockdata[l];
                correlation_blockdata_t         *blockdata      = &corr->blockdata[l];

                blockdata->blocklength          = blockdata_hist->blocklength;
                blockdata->correlation_integral = blockdata_hist->correlation_integral;
                blockdata->sum_wx               = blockdata_hist->sum_wx;
                blockdata->sum_wy               = blockdata_hist->sum_wy;
                blockdata->sum_w                = blockdata_hist->sum_w;
                blockdata->blockindex_prev      = blockdata_hist->blockindex_prev;
            }
        }
    }
}
