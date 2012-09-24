/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef GMX_OMP_NTHREADS
#define GMX_OMP_NTHREADS

/*! Enum values corresponding to multithreaded algorithmic modules. */
typedef enum module_nth
{
    /* Default is meant to be used in OMP regions outside the named
     * algorithmic modules listed below. */
    emntDefault, emntDomdec, emntPairsearch, emntNonbonded,
    emntBonded, emntPME,  emntUpdate, emntLINCS, emntSETTLE,
    emntNR
} module_nth_t;

/*! Initializes the per-module thread count. It is compatible with tMPI, 
 *  thread-safety is ensured (for the features available with tMPI). 
 *  This function should caled only once during the initialization of mdrun. */
void gmx_omp_nthreads_init(FILE *fplog, t_commrec *cr,
                           int omp_nthreads_req,
                           int omp_nthreads_pme_req,
                           gmx_bool bCurrNodePMEOnly,
                           gmx_bool bFullOmpSupport);

/*! Detect the maximum number of cores per node. This function should be called
 *  before thread-MPI is initialized so the Intel OpenMP doesn't outsmart us
 *  (by returning omp_get_max_threads=1 when nt=max virt. logical cores. */
void gmx_omp_nthreads_detecthw();

/*! Returns the number of threads to be used in the given module m. */
int gmx_omp_nthreads_get(int mod);

/*! Read the OMP_NUM_THREADS env. var. and check against the value set on the command line. */
void gmx_omp_nthreads_read_env(int *nthreads_omp);

#endif /* GMX_OMP_NTHREADS */
