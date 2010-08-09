/*
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _mdebin_bar_h
#define _mdebin_bar_h

#ifdef __cplusplus
extern "C" {
#endif

/* The functions & data structures here describe writing 
   energy differences (or their histogram )for use with g_bar */

typedef struct 
{
    real *dh; /* the raw energy differences */
    unsigned int ndh; /* number of data points */
    unsigned int ndhmax; /* the maximum number of points */

    int *hist; /* the histogram values */
    unsigned int nbins; /* the number of bins for the histogram to create */
    gmx_large_int_t start; /* the starting point in units of spacing of the 
                              histogram */
    double spacing; /* the histogram spacing in kJ/mol */
    unsigned int maxbin; /* highest bin with data */

    double lambda; /* the 'foreign' lambda value associated with this delta H */
    bool write_hist; /* whether to write histograms or raw data */
    bool written; /* whether this data has already been written out */

    double subblock_d[4]; /* data for an mdebin subblock for I/O. */
    gmx_large_int_t subblock_l[2]; /* data for an mdebin subblock for I/O.  */
} t_mde_delta_h;

/* the type definition is in mdebin.h */
struct t_mde_delta_h_coll
{
    t_mde_delta_h *dh; /* the delta hs */
    int ndh; /* the number of delta_h structures */

    double lambda; /* the native lambda associated with the free energy 
                      calculations */
    double temp; /* the temperature */

    double starttime; /* start time of the current dh collection */
    double endtime; /* end time of the current dh collection */
    bool starttime_set; /* whether the start time has been set */

    double subblock_d[4]; /* data for writing an mdebin subblock for I/O */
};



/* initialize a collection of delta h histograms/sets 
    dhc = the collection
    temp = temperature
    native_lambda = the native lambda (we're assuminng this doesn't change)
    dh_table_size = the histogram size as given by the .mdp option
    table_spacing = the histogram spacing (dx) as given by the .mdp option
    ndhmax = the maximum size of the du buffer
    n_dh = the number of foreign lambdas (.mdp option)
    flambda = the array of foreign lambdas (.mdp option) */
void mde_delta_h_coll_init(t_mde_delta_h_coll *dhc,
                           double temp,
                           double native_lambda,
                           int dh_table_size,
                           double table_spacing,
                           unsigned int ndhmax,
                           int n_dh,
                           double *flambda);

/* add a bunch of samples to the delta_h collection
    dhc = the collection
    U = the array with energies: from enerd->enerpart_lambda.
    time = the current simulation time. */
void mde_delta_h_coll_add_dh(t_mde_delta_h_coll *dhc, double *U, double time);

/* write the data associated with the du blocks collection as a collection
    of mdebin blocks.
    dhc = the collection
    fr = the enxio frame
    nblock = the current number of blocks */
void mde_delta_h_coll_handle_block(t_mde_delta_h_coll *dhc,   
                                   t_enxframe *fr, int nblock);


/* reset the collection of delta_h buffers for a new round of 
   data gathering */
void mde_delta_h_coll_reset(t_mde_delta_h_coll *dhc);


/* set the energyhistory variables to save state */
void mde_delta_h_coll_update_energyhistory(t_mde_delta_h_coll *dhc, 
                                           energyhistory_t *enerhist);

/* restore the variables from an energyhistory */
void mde_delta_h_coll_restore_energyhistory(t_mde_delta_h_coll *dhc, 
                                            energyhistory_t *enerhist);


#ifdef __cplusplus
}
#endif

#endif	/* _mdebin_bar_h */

