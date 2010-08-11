/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <float.h>
#include "typedefs.h"
#include "string2.h"
#include "gmx_fatal.h"
#include "mdebin.h"
#include "smalloc.h"
#include "enxio.h"
#include "gmxfio.h"
#include "mdebin_bar.h"

/* reset the delta_h list to prepare it for new values */
static void mde_delta_h_reset(t_mde_delta_h *dh)
{
    dh->ndh=0;
    dh->written=FALSE;
}

/* initialize the delta_h list */
static void mde_delta_h_init(t_mde_delta_h *dh, int nbins, 
                             double spacing, unsigned int  ndhmax, 
                             double foreign_lambda)
{
    dh->lambda=foreign_lambda;

    dh->nbins=nbins;
    dh->spacing=spacing;
    dh->ndhmax=ndhmax+2;

    snew(dh->dh, ndhmax);
    if ( nbins < 0 || spacing<GMX_REAL_EPS*10 )
    {
        dh->write_hist=FALSE;
    }
    else
    {
        dh->write_hist=TRUE;
        /* pre-allocate the histogram */
        snew(dh->hist, dh->nbins); 
    }
    mde_delta_h_reset(dh);
}

/* free the contents of the delta_h list */
static void mde_delta_h_free(t_mde_delta_h *dh)
{
    sfree(dh->dh);
    dh->dh=NULL;
}

/* Add a value to the delta_h list */
static void mde_delta_h_add_dh(t_mde_delta_h *dh, double delta_h, double time)
{
    if (dh->ndh >= dh->ndhmax)
    {
        gmx_incons("delta_h array not big enough!");
    }
    dh->dh[dh->ndh]=delta_h;
    dh->ndh++;
}

static void mde_delta_h_make_hist(t_mde_delta_h *dh)
{ 
    double min_dh = FLT_MAX;
    double max_dh = -FLT_MAX;
    unsigned int i;
    double max_dh_int;

    /* first find min and max */
    for(i=0;i<dh->ndh;i++)
    {
        if (dh->dh[i] < min_dh)
            min_dh=dh->dh[i];
        if (dh->dh[i] > max_dh)
            max_dh=dh->dh[i];
    }
    
    /* reset the histogram */
    for(i=0;i<dh->nbins;i++)
    {
        dh->hist[i]=0;
    }
    dh->maxbin=0;

    /* The starting point of the histogram is the lowest value found: 
       that value has the highest contribution to the free energy. 

       Get this start value in number of histogram spacings from zero, 
       as an integer.*/
    dh->start = (gmx_large_int_t)(min_dh/dh->spacing);

    max_dh_int=dh->start + (dh->nbins*dh->spacing);

    /* and fill the histogram*/
    for(i=0;i<dh->ndh;i++)
    {
        unsigned int bin;

        /* Determine the bin number. If it doesn't fit into the histogram, 
           add it to the last bin. 
           We check the max_dh_int range because converting to integers 
           might lead to overflow with unpredictable results.*/
        if (dh->dh[i] <= max_dh_int )
        {
            bin = (unsigned int)( (dh->dh[i] - min_dh)/dh->spacing );
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
        if (bin > dh->maxbin)
        {
            dh->maxbin = bin;
        }

        dh->hist[bin]++;
    }
    /* make sure we include a bin with 0 if we didn't use the full 
       histogram width */
    if (dh->maxbin < dh->nbins-1)
        dh->maxbin += 1;
}


void mde_delta_h_handle_block(t_mde_delta_h *dh, t_enxblock *blk)
{
    /* first check which type we should use: if a histogram
       actually turns out to be bigger than raw data, just 
       write raw data */
    if (dh->write_hist)
    {
        /* A histogram consists of 3 subblocks: the foreign labmda value +
           histogram spacing, the starting point, and the histogram data. */
        add_subblocks_enxblock(blk, 3);
        blk->id=enxDHHIST;

        /* check if there's actual data to be written. */
        if (dh->ndh > 1)
        {
            /* Make the histogram */
            mde_delta_h_make_hist(dh);
            dh->written=TRUE;
        }

        /* subblock 1: the foreign lambda value + the histogram spacing */
        dh->subblock_d[0]=dh->lambda;
        dh->subblock_d[1]=dh->spacing;
        blk->sub[0].nr=2;
        blk->sub[0].type=xdr_datatype_double;
        blk->sub[0].dval=dh->subblock_d;

        /* subblock 2: the starting point as a long integer */
        blk->sub[1].nr=1;
        blk->sub[1].type=xdr_datatype_large_int;
        blk->sub[1].lval=&(dh->start);

        /* subblock 3: the histogram data */
        if (dh->ndh > 1)
        {
            blk->sub[2].nr=dh->maxbin+1; /* it's +1 because size=index+1 in C */
            blk->sub[2].type=xdr_datatype_int;
            blk->sub[2].ival=dh->hist;
        }
        else
        {
            blk->sub[2].nr=0;
            blk->sub[2].type=xdr_datatype_int;
            blk->sub[2].ival=NULL;
        }
    }
    else
    {
        /* the histogram is bigger, we write raw data.
           Raw data consists of 2 subblocks: a block with the 
           the foreign lambda, and the data itself */
        add_subblocks_enxblock(blk, 2);

        blk->id=enxDH;


        /* subblock 1 */
        dh->subblock_d[0]=dh->lambda;
        blk->sub[0].nr=1;
        blk->sub[0].type=xdr_datatype_double;
        blk->sub[0].dval=dh->subblock_d;

        /* subblock 2 */
        /* check if there's actual data to be written. */
        if (dh->ndh > 1)
        {
            blk->sub[1].nr=dh->ndh;
#ifndef GMX_DOUBLE
            blk->sub[1].type=xdr_datatype_float;
            blk->sub[1].fval=dh->dh;
#else
            blk->sub[1].type=xdr_datatype_double;
            blk->sub[1].dval=dh->dh;
#endif
            dh->written=TRUE;
        }
        else
        {
            blk->sub[1].nr=0;
#ifndef GMX_DOUBLE
            blk->sub[1].type=xdr_datatype_float;
#else
            blk->sub[1].type=xdr_datatype_double;
#endif
            blk->sub[1].dval=NULL;
        }
    }
}

/* initialize the collection*/
void mde_delta_h_coll_init(t_mde_delta_h_coll *dhc, 
                           double temp,  
                           double native_lambda, 
                           int table_size, 
                           double table_spacing,
                           unsigned int ndhmax, 
                           int n_dh,
                           double *flambda)
{
    int i; 

    dhc->temp=temp;
    dhc->lambda=native_lambda;
    dhc->starttime=dhc->endtime=0.;
    dhc->starttime_set=FALSE;

    snew(dhc->dh, n_dh);
    dhc->ndh=n_dh;

    for(i=0;i<n_dh;i++)
    {
        mde_delta_h_init(dhc->dh + i, table_size, table_spacing, ndhmax, 
                         flambda[i] );
    }
}

/* add a bunch of samples */
void mde_delta_h_coll_add_dh(t_mde_delta_h_coll *dhc, double *U, double time)
{
    int i;

    if (!dhc->starttime_set)
    {
        dhc->starttime_set=TRUE;
        dhc->starttime=time;
    }
    dhc->endtime=time;
    for(i=0;i<dhc->ndh;i++)
    {
        mde_delta_h_add_dh(dhc->dh + i, U[i+1] - U[0], time);
    }
}

/* write the data associated with all the du blocks, but not the blocks 
   themselves. */
void mde_delta_h_coll_handle_block(t_mde_delta_h_coll *dhc,
                                   t_enxframe *fr, int nblock)
{
    int i;
    t_enxblock *blk;

    /* add one block with one subblock as the collection's own data */
    nblock++;
    add_blocks_enxframe(fr, nblock);
    blk=fr->block + (nblock-1);

    add_subblocks_enxblock(blk, 1);

    dhc->subblock_d[0] = dhc->temp;
    dhc->subblock_d[1] = dhc->lambda;
    dhc->subblock_d[2] = dhc->starttime;
    dhc->subblock_d[3] = dhc->endtime;

    blk->id=enxDHCOLL;
    blk->sub[0].nr=4;
    blk->sub[0].type=xdr_datatype_double;
    blk->sub[0].dval=dhc->subblock_d;

    for(i=0;i<dhc->ndh;i++)
    {
        nblock++;
        add_blocks_enxframe(fr, nblock);
        blk=fr->block + (nblock-1);

        mde_delta_h_handle_block(dhc->dh+i, blk);
    }
}

/* reset the data for a new round */
void mde_delta_h_coll_reset(t_mde_delta_h_coll *dhc)
{
    int i;
    for(i=0;i<dhc->ndh;i++)
    {
        if (dhc->dh[i].written)
        {
            /* we can now throw away the data */
            mde_delta_h_reset(dhc->dh + i);
        }
    }
    dhc->starttime_set=FALSE;
}

/* set the energyhistory variables to save state */
void mde_delta_h_coll_update_energyhistory(t_mde_delta_h_coll *dhc,
                                           energyhistory_t *enerhist)
{
    int i;
    if (!enerhist->dht)
    {
        snew(enerhist->dht, 1);
        snew(enerhist->dht->ndh, dhc->ndh);
        snew(enerhist->dht->dh, dhc->ndh);
        enerhist->dht->nndh=dhc->ndh;

        /* these don't change during the simulation */
        for(i=0;i<dhc->ndh;i++)
        {
            enerhist->dht->dh[i] = dhc->dh[i].dh;
        }
    }
    else
    {
        if (enerhist->dht->nndh != dhc->ndh)
            gmx_incons("energy history number of delta_h histograms != inputrec's number");
    }
    for(i=0;i<dhc->ndh;i++)
    {
        enerhist->dht->ndh[i] = dhc->dh[i].ndh;
    }
    enerhist->dht->starttime=dhc->starttime;
}



/* restore the variables from an energyhistory */
void mde_delta_h_coll_restore_energyhistory(t_mde_delta_h_coll *dhc,
                                            energyhistory_t *enerhist)
{
    int i;
    unsigned int j;

    if (dhc && !enerhist->dht)
        gmx_incons("No delta_h histograms in energy history");
    if (enerhist->dht->nndh != dhc->ndh)
        gmx_incons("energy history number of delta_h histograms != inputrec's number");

    for(i=0;i<enerhist->dht->nndh;i++)
    {
        dhc->dh[i].ndh=enerhist->dht->ndh[i];
        for(j=0;j<dhc->dh[i].ndh;j++)
        {
            dhc->dh[i].dh[j] = enerhist->dht->dh[i][j];
        }
    }
    dhc->starttime=enerhist->dht->starttime;
    if (dhc->dh[0].ndh > 0)
        dhc->starttime_set=TRUE;
    else
        dhc->starttime_set=FALSE;
}




