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

#ifdef GMX_OPENMP
#include <omp.h>
#endif

/*! Structure which holds the number of threads for each OpenMP multi-threaded
 *  algorithmic part/module of mdrun. */
typedef struct 
{
    int domdec;
    int pairsearch;
    int nonbonded;
    int bonded;
    int pme;
    int update;
    int lincs;
    int settle;
    int gmxdefault; /* it is meant to be used in omp regions outside the
                       algorithmic parts listed below */
    gmx_bool initialized;
} omp_module_nthreads_t;

/*! Number of threads for each algorithmic module.
 *
 * All fields are initialized to 0 which should result in erros if 
 * the init call is omitted */
static omp_module_nthreads_t mod_nth = {0, 0, 0, 0, 0, 0, 0, 0, FALSE};

/*! Enum values corresponding to multithreaded algorithmic modules. */
enum 
{ 
    emodGMXDefault, emodDomdec, emodPairsearch, emodNonbonded,
    emodBonded, emodPME,  emodUpdate, emodLincs, emodSettle,
    emodNR 
};

/*! Names of environment variables through which the number of threads of 
 *  the corresponding multithreaded algorithmic modules can be set from the 
 *  command line. */
static const char *mod_nth_env_var[emodNR] = 
{ 
    "GMX_OMP_DEFAULT_NTHREADS should never be set",  /* as it has to be = OMP_NUM_THREADS */
    "GMX_OMP_DOMDEC_NTHREADS", "GMX_OMP_PAIRSEARCH_NTHREADS", "GMX_OMP_NONBONDED_NTHREADS",
    "GMX_OMP_BONDED_NTHREADS", "GMX_OMP_PME_NTHREADS",
    "GMX_OMP_UPDATE_NTHREADS", "GMX_OMP_LINCS_NTHREADS", "GMX_OMP_SETTLE_NTHREADS"
};

/*! Determine the number of threads for the given module. The \module variable 
 *  takes values form the above enum and maps these to the corrreponding value 
 *  in mod_nth_env_var. */
static int pick_module_nthreads(int module)
{
    char *env;
    int  nth, omp_maxth;
    char sbuf[STRLEN];
   
    /* set max nthreads to 1 if not using OpenMP */
#ifdef GMX_OPENMP
    omp_maxth = omp_get_max_threads();
#else
    omp_maxth = 1;
#endif

    /* the default should never be set through a GMX_OMP_ environment variable */
    if (module != emodGMXDefault && (env = getenv(mod_nth_env_var[module])) != NULL)
    {
#ifndef GMX_OPENMP
        sprintf(sbuf, "mdrun compiled without OpenMP, but %s is set!\n");
        gmx_warning(sbuf);
#endif
        sscanf(env, "%d", &nth);
    }
    else 
    {
        nth = omp_maxth;
    }

    return nth;
}

void init_module_nthreads(t_commrec *cr)
{
    /* return if struct has already been initialized */
    if (mod_nth.initialized)
    {
        return;
    }

#ifdef GMX_THREADS
    /* This is required for tMPI thread-safety. It is not thread-safe with
     * multiple simulations, but that's anyway not supported by tMPI. */
    if (SIMMASTER(cr))
#endif        
    {
        mod_nth.gmxdefault  = pick_module_nthreads(emodGMXDefault);
        mod_nth.domdec      = pick_module_nthreads(emodDomdec);
        mod_nth.pairsearch  = pick_module_nthreads(emodPairsearch);
        mod_nth.nonbonded   = pick_module_nthreads(emodNonbonded);
        mod_nth.bonded      = pick_module_nthreads(emodBonded);
        mod_nth.pme         = pick_module_nthreads(emodPME);
        mod_nth.update      = pick_module_nthreads(emodUpdate);
        mod_nth.lincs       = pick_module_nthreads(emodLincs);
        mod_nth.settle      = pick_module_nthreads(emodSettle);
    }
#ifdef GMX_THREADS
    if (PAR(cr))
    {
        MPI_Barrier(cr->mpi_comm_mysim);
    }
#endif

    mod_nth.initialized = TRUE;
}

int gmx_omp_get_default_nthreads()
{
    return mod_nth.gmxdefault;
}

int gmx_omp_get_domdec_nthreads()
{
    return mod_nth.domdec;
}

int gmx_omp_get_pairsearch_nthreads()
{
    return mod_nth.pairsearch;
}

int gmx_omp_get_nonbonded_nthreads()
{
    return mod_nth.nonbonded;
}

int gmx_omp_get_bonded_nthreads()
{
    return mod_nth.bonded;
}

int gmx_omp_get_pme_nthreads()
{
    return mod_nth.pme;
}

int gmx_omp_get_update_nthreads()
{
    return mod_nth.update;
}

int gmx_omp_get_lincs_nthreads()
{
    return mod_nth.lincs;
}

int gmx_omp_get_settle_nthreads()
{
    return mod_nth.settle;
}
