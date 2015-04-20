/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "mdebin_bar.h"

#include <float.h>
#include <math.h>
#include <string.h>

#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/mdebin.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* reset the delta_h list to prepare it for new values */
static void mde_delta_h_reset(t_mde_delta_h *dh)
{
    dh->ndh     = 0;
    dh->written = FALSE;
}

/* initialize the delta_h list */
static void mde_delta_h_init(t_mde_delta_h *dh, int nbins,
                             double dx, unsigned int  ndhmax,
                             int type, int derivative, int nlambda,
                             double *lambda)
{
    int i;

    dh->type       = type;
    dh->derivative = derivative;
    dh->lambda     = lambda;
    dh->nlambda    = nlambda;

    snew(dh->lambda, nlambda);
    for (i = 0; i < nlambda; i++)
    {
        dh->lambda[i] = lambda[i];
    }


    snew(dh->subblock_meta_d, dh->nlambda+1);

    dh->ndhmax = ndhmax+2;
    for (i = 0; i < 2; i++)
    {
        dh->bin[i] = NULL;
    }

    snew(dh->dh, dh->ndhmax);
    snew(dh->dhf, dh->ndhmax);

    if (nbins <= 0 || dx < GMX_REAL_EPS*10)
    {
        dh->nhist = 0;
    }
    else
    {
        int i;
        /* pre-allocate the histogram */
        dh->nhist = 2; /* energies and derivatives histogram */
        dh->dx    = dx;
        dh->nbins = nbins;
        for (i = 0; i < dh->nhist; i++)
        {
            snew(dh->bin[i], dh->nbins);
        }
    }
    mde_delta_h_reset(dh);
}

/* Add a value to the delta_h list */
static void mde_delta_h_add_dh(t_mde_delta_h *dh, double delta_h)
{
    if (dh->ndh >= dh->ndhmax)
    {
        gmx_incons("delta_h array not big enough!");
    }
    dh->dh[dh->ndh] = delta_h;
    dh->ndh++;
}

/* construct histogram with index hi */
static void mde_delta_h_make_hist(t_mde_delta_h *dh, int hi, gmx_bool invert)
{
    double       min_dh = FLT_MAX;
    double       max_dh = -FLT_MAX;
    unsigned int i;
    double       max_dh_hist; /* maximum binnable dh value */
    double       min_dh_hist; /* minimum binnable dh value */
    double       dx = dh->dx;
    double       f;           /* energy mult. factor */

    /* by applying a -1 scaling factor on the energies we get the same as
       having a negative dx, but we don't need to fix the min/max values
       beyond inverting x0 */
    f = invert ? -1 : 1;

    /* first find min and max */
    for (i = 0; i < dh->ndh; i++)
    {
        if (f*dh->dh[i] < min_dh)
        {
            min_dh = f*dh->dh[i];
        }
        if (f*dh->dh[i] > max_dh)
        {
            max_dh = f*dh->dh[i];
        }
    }

    /* reset the histogram */
    for (i = 0; i < dh->nbins; i++)
    {
        dh->bin[hi][i] = 0;
    }
    dh->maxbin[hi] = 0;

    /* The starting point of the histogram is the lowest value found:
       that value has the highest contribution to the free energy.

       Get this start value in number of histogram dxs from zero,
       as an integer.*/

    dh->x0[hi] = (gmx_int64_t)floor(min_dh/dx);

    min_dh_hist = (dh->x0[hi])*dx;
    max_dh_hist = (dh->x0[hi] + dh->nbins + 1)*dx;

    /* and fill the histogram*/
    for (i = 0; i < dh->ndh; i++)
    {
        unsigned int bin;

        /* Determine the bin number. If it doesn't fit into the histogram,
           add it to the last bin.
           We check the max_dh_int range because converting to integers
           might lead to overflow with unpredictable results.*/
        if ( (f*dh->dh[i] >= min_dh_hist) && (f*dh->dh[i] <= max_dh_hist ) )
        {
            bin = (unsigned int)( (f*dh->dh[i] - min_dh_hist)/dx );
        }
        else
        {
            bin = dh->nbins-1;
        }
        /* double-check here because of possible round-off errors*/
        if (bin >= dh->nbins)
        {
            bin = dh->nbins-1;
        }
        if (bin > dh->maxbin[hi])
        {
            dh->maxbin[hi] = bin;
        }

        dh->bin[hi][bin]++;
    }

    /* make sure we include a bin with 0 if we didn't use the full
       histogram width. This can then be used as an indication that
       all the data was binned. */
    if (dh->maxbin[hi] < dh->nbins-1)
    {
        dh->maxbin[hi] += 1;
    }
}


void mde_delta_h_handle_block(t_mde_delta_h *dh, t_enxblock *blk)
{
    /* first check which type we should use: histogram or raw data */
    if (dh->nhist == 0)
    {
        int i;

        /* We write raw data.
           Raw data consists of 3 subblocks: an int metadata block
           with type and derivative index, a foreign lambda block
           and and the data itself */
        add_subblocks_enxblock(blk, 3);

        blk->id = enxDH;

        /* subblock 1 */
        dh->subblock_meta_i[0] = dh->type;       /* block data type */
        dh->subblock_meta_i[1] = dh->derivative; /* derivative direction if
                                                    applicable (in indices
                                                    starting from first coord in
                                                    the main delta_h_coll) */
        blk->sub[0].nr         = 2;
        blk->sub[0].type       = xdr_datatype_int;
        blk->sub[0].ival       = dh->subblock_meta_i;

        /* subblock 2 */
        for (i = 0; i < dh->nlambda; i++)
        {
            dh->subblock_meta_d[i] = dh->lambda[i];
        }
        blk->sub[1].nr   = dh->nlambda;
        blk->sub[1].type = xdr_datatype_double;
        blk->sub[1].dval = dh->subblock_meta_d;

        /* subblock 3 */
        /* check if there's actual data to be written. */
        /*if (dh->ndh > 1)*/
        if (dh->ndh > 0)
        {
            unsigned int i;

            blk->sub[2].nr = dh->ndh;
/* For F@H for now. */
#undef GMX_DOUBLE
#ifndef GMX_DOUBLE
            blk->sub[2].type = xdr_datatype_float;
            for (i = 0; i < dh->ndh; i++)
            {
                dh->dhf[i] = (float)dh->dh[i];
            }
            blk->sub[2].fval = dh->dhf;
#else
            blk->sub[2].type = xdr_datatype_double;
            blk->sub[2].dval = dh->dh;
#endif
            dh->written = TRUE;
        }
        else
        {
            blk->sub[2].nr = 0;
#ifndef GMX_DOUBLE
            blk->sub[2].type = xdr_datatype_float;
            blk->sub[2].fval = NULL;
#else
            blk->sub[2].type = xdr_datatype_double;
            blk->sub[2].dval = NULL;
#endif
        }
    }
    else
    {
        int nhist_written = 0;
        int i;
        int k;

        /* TODO histogram metadata */
        /* check if there's actual data to be written. */
        if (dh->ndh > 1)
        {
            gmx_bool prev_complete = FALSE;
            /* Make the histogram(s) */
            for (i = 0; i < dh->nhist; i++)
            {
                if (!prev_complete)
                {
                    /* the first histogram is always normal, and the
                       second one is always reverse */
                    mde_delta_h_make_hist(dh, i, i == 1);
                    nhist_written++;
                    /* check whether this histogram contains all data: if the
                       last bin is 0, it does */
                    if (dh->bin[i][dh->nbins-1] == 0)
                    {
                        prev_complete = TRUE;
                    }
                    if (!dh->derivative)
                    {
                        prev_complete = TRUE;
                    }
                }
            }
            dh->written = TRUE;
        }

        /* A histogram consists of 2, 3 or 4 subblocks:
           the foreign lambda value + histogram spacing, the starting point,
           and the histogram data (0, 1 or 2 blocks). */
        add_subblocks_enxblock(blk, nhist_written+2);
        blk->id = enxDHHIST;

        /* subblock 1: the lambda value + the histogram spacing */
        if (dh->nlambda == 1)
        {
            /* for backward compatibility */
            dh->subblock_meta_d[0] = dh->lambda[0];
        }
        else
        {
            dh->subblock_meta_d[0] = -1;
            for (i = 0; i < dh->nlambda; i++)
            {
                dh->subblock_meta_d[2+i] = dh->lambda[i];
            }
        }
        dh->subblock_meta_d[1] = dh->dx;
        blk->sub[0].nr         = 2+ ((dh->nlambda > 1) ? dh->nlambda : 0);
        blk->sub[0].type       = xdr_datatype_double;
        blk->sub[0].dval       = dh->subblock_meta_d;

        /* subblock 2: the starting point(s) as a long integer */
        dh->subblock_meta_l[0] = nhist_written;
        dh->subblock_meta_l[1] = dh->type; /*dh->derivative ? 1 : 0;*/
        k = 2;
        for (i = 0; i < nhist_written; i++)
        {
            dh->subblock_meta_l[k++] = dh->x0[i];
        }
        /* append the derivative data */
        dh->subblock_meta_l[k++] = dh->derivative;

        blk->sub[1].nr   = nhist_written+3;
        blk->sub[1].type = xdr_datatype_int64;
        blk->sub[1].lval = dh->subblock_meta_l;

        /* subblock 3 + 4 : the histogram data */
        for (i = 0; i < nhist_written; i++)
        {
            blk->sub[i+2].nr   = dh->maxbin[i]+1; /* it's +1 because size=index+1
                                                     in C */
            blk->sub[i+2].type = xdr_datatype_int;
            blk->sub[i+2].ival = dh->bin[i];
        }
    }
}

/* initialize the collection*/
void mde_delta_h_coll_init(t_mde_delta_h_coll *dhc, const t_inputrec *ir)
{
    int       i, j, n;
    double    lambda;
    double   *lambda_vec;
    int       ndhmax = ir->nstenergy/ir->nstcalcenergy;
    t_lambda *fep    = ir->fepvals;

    dhc->temperature    = ir->opts.ref_t[0]; /* only store system temperature */
    dhc->start_time     = 0.;
    dhc->delta_time     = ir->delta_t*ir->fepvals->nstdhdl;
    dhc->start_time_set = FALSE;

    /* this is the compatibility lambda value. If it is >=0, it is valid,
       and there is either an old-style lambda or a slow growth simulation. */
    dhc->start_lambda = ir->fepvals->init_lambda;
    /* for continuous change of lambda values */
    dhc->delta_lambda = ir->fepvals->delta_lambda*ir->fepvals->nstdhdl;

    if (dhc->start_lambda < 0)
    {
        /* create the native lambda vectors */
        dhc->lambda_index = fep->init_fep_state;
        dhc->n_lambda_vec = 0;
        for (i = 0; i < efptNR; i++)
        {
            if (fep->separate_dvdl[i])
            {
                dhc->n_lambda_vec++;
            }
        }
        snew(dhc->native_lambda_vec, dhc->n_lambda_vec);
        snew(dhc->native_lambda_components, dhc->n_lambda_vec);
        j = 0;
        for (i = 0; i < efptNR; i++)
        {
            if (fep->separate_dvdl[i])
            {
                dhc->native_lambda_components[j] = i;
                if (fep->init_fep_state >= 0 &&
                    fep->init_fep_state < fep->n_lambda)
                {
                    dhc->native_lambda_vec[j] =
                        fep->all_lambda[i][fep->init_fep_state];
                }
                else
                {
                    dhc->native_lambda_vec[j] = -1;
                }
                j++;
            }
        }
    }
    else
    {
        /* don't allocate the meta-data subblocks for lambda vectors */
        dhc->native_lambda_vec        = NULL;
        dhc->n_lambda_vec             = 0;
        dhc->native_lambda_components = 0;
        dhc->lambda_index             = -1;
    }
    /* allocate metadata subblocks */
    snew(dhc->subblock_d, 5 + dhc->n_lambda_vec);
    snew(dhc->subblock_i, 1 + dhc->n_lambda_vec);

    /* now decide which data to write out */
    dhc->nlambda     = 0;
    dhc->ndhdl       = 0;
    dhc->dh_expanded = NULL;
    dhc->dh_energy   = NULL;
    dhc->dh_pv       = NULL;

    /* total number of raw data point collections in the sample */
    dhc->ndh = 0;

    {
        gmx_bool bExpanded           = FALSE;
        gmx_bool bEnergy             = FALSE;
        gmx_bool bPV                 = FALSE;
        int      n_lambda_components = 0;

        /* first count the number of states */

        /* add the dhdl's */
        if (fep->dhdl_derivatives == edhdlderivativesYES)
        {
            for (i = 0; i < efptNR; i++)
            {
                if (ir->fepvals->separate_dvdl[i])
                {
                    dhc->ndh   += 1;
                    dhc->ndhdl += 1;
                }
            }
        }
        /* add the lambdas */
        dhc->nlambda = ir->fepvals->lambda_stop_n - ir->fepvals->lambda_start_n;
        dhc->ndh    += dhc->nlambda;
        /* another compatibility check */
        if (dhc->start_lambda < 0)
        {
            /* include one more for the specification of the state, by lambda or
               fep_state*/
            if (ir->expandedvals->elmcmove > elmcmoveNO)
            {
                dhc->ndh += 1;
                bExpanded = TRUE;
            }
            /* whether to print energies */
            if (ir->fepvals->edHdLPrintEnergy != edHdLPrintEnergyNO)
            {
                dhc->ndh += 1;
                bEnergy   = TRUE;
            }
            if (ir->epc > epcNO)
            {
                dhc->ndh += 1;  /* include pressure-volume work */
                bPV       = TRUE;
            }
        }
        /* allocate them */
        snew(dhc->dh, dhc->ndh);

        /* now initialize them */
        /* the order, for now, must match that of the dhdl.xvg file because of
           how g_energy -odh is implemented */
        n = 0;
        if (bExpanded)
        {
            dhc->dh_expanded = dhc->dh+n;
            mde_delta_h_init(dhc->dh+n, ir->fepvals->dh_hist_size,
                             ir->fepvals->dh_hist_spacing, ndhmax,
                             dhbtEXPANDED, 0, 0, NULL);
            n++;
        }
        if (bEnergy)
        {
            dhc->dh_energy = dhc->dh+n;
            mde_delta_h_init(dhc->dh+n, ir->fepvals->dh_hist_size,
                             ir->fepvals->dh_hist_spacing, ndhmax,
                             dhbtEN, 0, 0, NULL);
            n++;
        }
        /* add the dhdl's */
        n_lambda_components = 0;
        if (fep->dhdl_derivatives == edhdlderivativesYES)
        {
            dhc->dh_dhdl = dhc->dh + n;
            for (i = 0; i < efptNR; i++)
            {
                if (ir->fepvals->separate_dvdl[i])
                {
                    /* we give it init_lambda for compatibility */
                    mde_delta_h_init(dhc->dh+n, ir->fepvals->dh_hist_size,
                                     ir->fepvals->dh_hist_spacing, ndhmax,
                                     dhbtDHDL, n_lambda_components, 1,
                                     &(fep->init_lambda));
                    n++;
                    n_lambda_components++;
                }
            }
        }
        else
        {
            for (i = 0; i < efptNR; i++)
            {
                if (ir->fepvals->separate_dvdl[i])
                {
                    n_lambda_components++; /* count the components */
                }
            }

        }
        /* add the lambdas */
        dhc->dh_du = dhc->dh + n;
        snew(lambda_vec, n_lambda_components);
        for (i = ir->fepvals->lambda_start_n; i < ir->fepvals->lambda_stop_n; i++)
        {
            int k = 0;

            for (j = 0; j < efptNR; j++)
            {
                if (ir->fepvals->separate_dvdl[j])
                {
                    lambda_vec[k++] = fep->all_lambda[j][i];
                }
            }

            mde_delta_h_init(dhc->dh+n, ir->fepvals->dh_hist_size,
                             ir->fepvals->dh_hist_spacing, ndhmax,
                             dhbtDH, 0, n_lambda_components, lambda_vec);
            n++;
        }
        sfree(lambda_vec);
        if (bPV)
        {
            dhc->dh_pv = dhc->dh+n;
            mde_delta_h_init(dhc->dh+n, ir->fepvals->dh_hist_size,
                             ir->fepvals->dh_hist_spacing, ndhmax,
                             dhbtPV, 0, 0, NULL);
            n++;
        }
    }
}

/* add a bunch of samples - note fep_state is double to allow for better data storage */
void mde_delta_h_coll_add_dh(t_mde_delta_h_coll *dhc,
                             double              fep_state,
                             double              energy,
                             double              pV,
                             double             *dhdl,
                             double             *foreign_dU,
                             double              time)
{
    int i;

    if (!dhc->start_time_set)
    {
        dhc->start_time_set = TRUE;
        dhc->start_time     = time;
    }

    for (i = 0; i < dhc->ndhdl; i++)
    {
        mde_delta_h_add_dh(dhc->dh_dhdl+i, dhdl[i]);
    }
    for (i = 0; i < dhc->nlambda; i++)
    {
        mde_delta_h_add_dh(dhc->dh_du+i, foreign_dU[i]);
    }
    if (dhc->dh_pv != NULL)
    {
        mde_delta_h_add_dh(dhc->dh_pv, pV);
    }
    if (dhc->dh_energy != NULL)
    {
        mde_delta_h_add_dh(dhc->dh_energy, energy);
    }
    if (dhc->dh_expanded != NULL)
    {
        mde_delta_h_add_dh(dhc->dh_expanded, fep_state);
    }

}

/* write the metadata associated with all the du blocks, and call
   handle_block to write out all the du blocks */
void mde_delta_h_coll_handle_block(t_mde_delta_h_coll *dhc,
                                   t_enxframe *fr, int nblock)
{
    int         i;
    t_enxblock *blk;

    /* add one block with one subblock as the collection's own data */
    nblock++;
    add_blocks_enxframe(fr, nblock);
    blk = fr->block + (nblock-1);

    /* only allocate lambda vector component blocks if they must be written out
       for backward compatibility */
    if (dhc->native_lambda_components != NULL)
    {
        add_subblocks_enxblock(blk, 2);
    }
    else
    {
        add_subblocks_enxblock(blk, 1);
    }

    dhc->subblock_d[0] = dhc->temperature;  /* temperature */
    dhc->subblock_d[1] = dhc->start_time;   /* time of first sample */
    dhc->subblock_d[2] = dhc->delta_time;   /* time difference between samples */
    dhc->subblock_d[3] = dhc->start_lambda; /* old-style lambda at starttime */
    dhc->subblock_d[4] = dhc->delta_lambda; /* lambda diff. between samples */
    /* set the lambda vector components if they exist */
    if (dhc->native_lambda_components != NULL)
    {
        for (i = 0; i < dhc->n_lambda_vec; i++)
        {
            dhc->subblock_d[5+i] = dhc->native_lambda_vec[i];
        }
    }
    blk->id          = enxDHCOLL;
    blk->sub[0].nr   = 5 + dhc->n_lambda_vec;
    blk->sub[0].type = xdr_datatype_double;
    blk->sub[0].dval = dhc->subblock_d;

    if (dhc->native_lambda_components != NULL)
    {
        dhc->subblock_i[0] = dhc->lambda_index;
        /* set the lambda vector component IDs if they exist */
        dhc->subblock_i[1] = dhc->n_lambda_vec;
        for (i = 0; i < dhc->n_lambda_vec; i++)
        {
            dhc->subblock_i[i+2] = dhc->native_lambda_components[i];
        }
        blk->sub[1].nr   = 2 + dhc->n_lambda_vec;
        blk->sub[1].type = xdr_datatype_int;
        blk->sub[1].ival = dhc->subblock_i;
    }

    for (i = 0; i < dhc->ndh; i++)
    {
        nblock++;
        add_blocks_enxframe(fr, nblock);
        blk = fr->block + (nblock-1);

        mde_delta_h_handle_block(dhc->dh+i, blk);
    }
}

/* reset the data for a new round */
void mde_delta_h_coll_reset(t_mde_delta_h_coll *dhc)
{
    int i;
    for (i = 0; i < dhc->ndh; i++)
    {
        if (dhc->dh[i].written)
        {
            /* we can now throw away the data */
            mde_delta_h_reset(dhc->dh + i);
        }
    }
    dhc->start_time_set = FALSE;
}

/* set the energyhistory variables to save state */
void mde_delta_h_coll_update_energyhistory(t_mde_delta_h_coll *dhc,
                                           energyhistory_t    *enerhist)
{
    int i;
    if (!enerhist->dht)
    {
        snew(enerhist->dht, 1);
        snew(enerhist->dht->ndh, dhc->ndh);
        snew(enerhist->dht->dh, dhc->ndh);
        enerhist->dht->nndh = dhc->ndh;
    }
    else
    {
        if (enerhist->dht->nndh != dhc->ndh)
        {
            gmx_incons("energy history number of delta_h histograms != inputrec's number");
        }
    }
    for (i = 0; i < dhc->ndh; i++)
    {
        enerhist->dht->dh[i]  = dhc->dh[i].dh;
        enerhist->dht->ndh[i] = dhc->dh[i].ndh;
    }
    enerhist->dht->start_time   = dhc->start_time;
    enerhist->dht->start_lambda = dhc->start_lambda;
}



/* restore the variables from an energyhistory */
void mde_delta_h_coll_restore_energyhistory(t_mde_delta_h_coll *dhc,
                                            energyhistory_t    *enerhist)
{
    int          i;
    unsigned int j;

    if (dhc && !enerhist->dht)
    {
        gmx_incons("No delta_h histograms in energy history");
    }
    if (enerhist->dht->nndh != dhc->ndh)
    {
        gmx_incons("energy history number of delta_h histograms != inputrec's number");
    }

    for (i = 0; i < enerhist->dht->nndh; i++)
    {
        dhc->dh[i].ndh = enerhist->dht->ndh[i];
        for (j = 0; j < dhc->dh[i].ndh; j++)
        {
            dhc->dh[i].dh[j] = enerhist->dht->dh[i][j];
        }
    }
    dhc->start_time = enerhist->dht->start_time;
    if (enerhist->dht->start_lambda_set)
    {
        dhc->start_lambda = enerhist->dht->start_lambda;
    }
    if (dhc->dh[0].ndh > 0)
    {
        dhc->start_time_set = TRUE;
    }
    else
    {
        dhc->start_time_set = FALSE;
    }
}
