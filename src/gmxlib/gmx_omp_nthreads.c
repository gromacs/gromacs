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
 * Copyright (c) 2001-2010, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "macros.h"
#include "gmx_omp_nthreads.h"

#ifdef GMX_OPENMP
#include <omp.h>
#endif

/*! Structure with the number of threads for each OpenMP multi-threaded
 *  algorithmic module in mdrun. */
typedef struct 
{
    int nth_pp;             /*! Global number of threads per MPI process/tMPI thread. */
    int nth[emntNR];        /*! Number of threads for each module, indexed with module_nth_t */
    gmx_bool initialized;   /*! TRUE if initialized. */
} omp_module_nthreads_t;

/*! Names of environment variables to set the per module number of threads.
 *
 *  Indexed with the values of module_nth_t.
 * */
static const char *modth_env_var[emntNR] =
{ 
    "GMX_DEFAULT_NUM_THREADS should never be set",
    "GMX_DOMDEC_NUM_THREADS", "GMX_PAIRSEARCH_NUM_THREADS",
    "GMX_NONBONDED_NUM_THREADS", "GMX_BONDED_NUM_THREADS",
    "GMX_PME_NUM_THREADS", "GMX_UPDATE_NUM_THREADS",
    "GMX_LINCS_NUM_THREADS", "GMX_SETTLE_NUM_THREADS"
};


/*! Number of threads for each algorithmic module.
 *
 *  File-scope global variable that gets set once in \init_module_nthreads
 *  and queried via gmx_omp_nthreads_get.
 *
 *  All fields are initialized to 0 which should result in erros if
 *  the init call is omitted
 * */
static omp_module_nthreads_t modth = { 0, {0, 0, 0, 0, 0, 0, 0, 0}, FALSE};


/*! Determine the number of threads for module \mod.
 *
 *  \mod takes values form the above enum and maps these to the corresponding
 *  value in modth_env_var. */
static int pick_module_nthreads(int mod)
{
    char *env;
    int  nth;
    char sbuf[STRLEN];
   
    /* the default should never be set through a GMX_OMP_ environment variable,
     * as it's always equal with nth_pp */
    if (mod == emntDefault)
    {
        return modth.nth[emntDefault];
    }

    if ((env = getenv(modth_env_var[mod])) != NULL)
    {
#ifndef GMX_OPENMP
        sprintf(sbuf, "%s is set, but mdrun is compiled without OpenMP!\n");
        gmx_warning(sbuf);
#endif
        sscanf(env, "%d", &nth);

        fprintf(stderr, "\n%s=%d set, using this value instead of the default %d\n",
                modth_env_var[mod], nth, modth.nth_pp);
    }
    else 
    {
        nth = modth.nth_pp;
    }

    return nth;
}

void gmx_omp_nthreads_init(t_commrec *cr)
{
    int  nth, omp_maxth, gmx_maxth, nppn;
    char *env;

    /* just return if it has already been initialized */
    if (modth.initialized)
    {
        return;
    }

    /* Set the number of threads to:
     * - 1 if not compiled with OpenMP on
     * - OMP_NUM_THREADS if the env. var is set (overrides everything!);
     * - take the max number of available or permitted (through the
     *   GMX_MAX_THREADS env var) threads and distribute them on the
     *   processes/tMPI threads.
     * */
#ifndef GMX_OPENMP
    nth = 1;
#else
    /* TODO OMP_NUM_THREADS and GMX_MAX_THREADS should be read only on the
     * master node and bcasted to the rest of the nodes */
    if (getenv("OMP_NUM_THREADS") != NULL)
    {
        nth = omp_get_max_threads();
    }
    else
    {
        /* max available threads per node */
        omp_maxth = omp_get_max_threads();
        /* max permitted threads per node set through env var */
        if ((env = getenv("GMX_MAX_THREADS")) != NULL)
        {
            sscanf(env, "%d",&gmx_maxth);
            nth = min(gmx_maxth, omp_maxth);

            if (nth != gmx_maxth && MASTER(cr))
            {
                fprintf(stderr, "\nOverriding GMX_MAX_THREADS=%d, "
                        "the maximum number of available threads is %d\n\n",
                         gmx_maxth, nth);
            }
        }
        else
        {
            nth = omp_maxth;
        }

        /* get the number of processes per node */
#ifdef GMX_THREADS
        nppn = cr->nnodes;
#else
#ifdef GMX_MPI
        /* FIXME: this doesn't work with separate PME nodes */
        MPI_Comm c_intra;
        MPI_Comm_split(MPI_COMM_WORLD, gmx_hostname_num() , gmx_node_rank(), &c_intra);
        MPI_Comm_size(c_intra, &nppn);
        MPI_Comm_free(&c_intra);
#else
        /* neither MPI nor tMPI */
        nppn = 1;
#endif /* GMX_MPI */
#endif /* GMX_THREADS */

        /* divide the threads within the MPI processes/tMPI threads */
        nth /= nppn;

        /* set the number of threads globally*/
        omp_set_num_threads(nth);
    }
#endif

    /* now we have the dafult, set it */
    modth.nth_pp = modth.nth[emntDefault] = nth;

    if (PAR(cr))
    {
        if (MASTER(cr))
        {
#ifdef GMX_THREADS
            fprintf(stderr, "Using %d OpenMP threads per tMPI thread\n", modth.nth_pp);
#else
            fprintf(stderr, "Using %d OpenMP threads per MPI process\n", modth.nth_pp);
#endif /* GMX_THREADS */
        }
    }
    else
    {
        fprintf(stderr, "Using %d OpenMP threads\n", modth.nth_pp);
    }

    /* now set the per-module values */
#ifdef GMX_THREADS
    /* This is required for tMPI thread-safety. It is not thread-safe with
     * multiple simulations, but that's anyway not supported by tMPI. */
    if (SIMMASTER(cr))
#endif        
    {
        modth.nth[emntDomdec]     = pick_module_nthreads(emntDomdec);
        modth.nth[emntPairsearch] = pick_module_nthreads(emntPairsearch);
        modth.nth[emntNonbonded]  = pick_module_nthreads(emntNonbonded);
        modth.nth[emntBonded]     = pick_module_nthreads(emntBonded);
        modth.nth[emntPME]        = pick_module_nthreads(emntPME);
        modth.nth[emntUpdate]     = pick_module_nthreads(emntUpdate);
        modth.nth[emntLincs]      = pick_module_nthreads(emntLincs);
        modth.nth[emntSettle]     = pick_module_nthreads(emntSettle);
    }
#ifdef GMX_THREADS
    if (PAR(cr))
    {
        MPI_Barrier(cr->mpi_comm_mysim);
    }
#endif

    modth.initialized = TRUE;
}

int gmx_omp_nthreads_get(int mod)
{
    if (mod < 0 || mod >= emntNR)
    {
        /* invalid module queried */
        return -1;
    }
    else
    {
        return modth.nth[mod];
    }
}
