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

/* Data for one foreign lambda, or derivative. */
typedef struct 
{
    real *dh; /* the raw energy differences */
    unsigned int ndh; /* number of data points */
    unsigned int ndhmax; /* the maximum number of points */

    int nhist; /* the number of histograms. There can either be
                  0 (for no histograms)
                  1 (for 'foreign lambda' histograms)
                  2 (for derivative histograms: there's
                     a 'forward' and 'backward' histogram
                     containing the minimum and maximum
                     values, respectively). */
    int *bin[2]; /* the histogram(s) */
    double dx; /* the histogram spacing in kJ/mol. This is the
                  same for the two histograms */
    unsigned int nbins; /* the number of bins */
    gmx_large_int_t x0[2]; /* the starting point in units of spacing 
                              of the histogram */
    unsigned int maxbin[2]; /* highest bin number with data */

    gmx_bool derivative; /* whether this delta_h contains derivatives */
    double lambda; /* the 'foreign' lambda value associated with this delta H */
    gmx_bool written; /* whether this data has already been written out */

    double subblock_d[4]; /* data for an mdebin subblock for I/O. */
    gmx_large_int_t subblock_l[4]; /* data for an mdebin subblock for I/O.  */
    int subblock_i[4]; /* data for an mdebin subblock for I/O.  */
} t_mde_delta_h;

/* the type definition is in mdebin.h */
struct t_mde_delta_h_coll
{
    t_mde_delta_h *dh; /* the delta hs */
    int ndh; /* the number of delta_h structures */
    int ndhdl; /* number of derivative delta_hs */

    double start_time; /* start time of the current dh collection */
    double delta_time; /* time difference between samples */
    gmx_bool start_time_set; /* whether the start time has been set */

    double start_lambda; /* the native lambda associated with the free energy 
                           calculations (at the time of the first sample) */
    double delta_lambda; /* lambda difference between samples */

    double temp; /* the temperature */
    double subblock_d[5]; /* data for writing an mdebin subblock for I/O */
};



/* initialize a collection of delta h histograms/sets 
    dhc = the collection
    ir = the input record */
void mde_delta_h_coll_init(t_mde_delta_h_coll *dhc,
                           const t_inputrec *ir);

/* add a bunch of samples to the delta_h collection
    dhc = the collection
    dhdl = the hamiltonian derivative
    U = the array with energies: from enerd->enerpart_lambda.
    time = the current simulation time. */
void mde_delta_h_coll_add_dh(t_mde_delta_h_coll *dhc, 
                             double dhdl,
                             double *U, double time,
                             double native_lambda);

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

