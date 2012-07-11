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

#include <stdlib.h>
#include <stdio.h>

#include "config.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "tables.h"
/* FIXME get rid of this dirty hack which is here just to avoid pulling in mpi.h */
#define _fcdata_h
#define _iteratedconstraints_h
#define _ifunc_h
#include "typedefs.h"
#undef  _fcdata_h
#undef  _iteratedconstraints_h
#undef  _ifunc_h
#include "types/nb_verlet.h"
#include "types/interaction_const.h"
#include "types/force_flags.h"

#include "nbnxn_cuda_types.h"
#include "../../gmxlib/cuda_tools/cudautils.cuh"
#include "nbnxn_cuda_data_mgmt.h"
#include "pmalloc_cuda.h"

#define USE_CUDA_EVENT_BLOCKING_SYNC FALSE  /* makes the CPU thread block */

/* coulomb force talble size chosen such that it fits along the non-bonded 
   parameters in the texture cache */
#define EWALD_COULOMB_FORCE_TABLE_SIZE (1536)

#define MY_PI               (3.1415926535897932384626433832795)
#define TWO_OVER_SQRT_PI    (2.0/sqrt(MY_PI))
    
#define TIME_GPU_TRANSFERS 1

#define NUM_NB_KERNELS 12

static void nbnxn_cuda_clear_e_fshift(nbnxn_cuda_ptr_t cu_nb);

/*! Names of default CUDA nonbonded kernel functions with mangling
    (used to be development k3). */
static const char * const nb_default_knames[NUM_NB_KERNELS] =
{
    "_Z10k_nbnxn_rf11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z13k_nbnxn_ewald11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z14k_nbnxn_cutoff11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z15k_nbnxn_rf_ener11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z16k_nbnxn_rf_prune11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z18k_nbnxn_ewald_ener11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z19k_nbnxn_ewald_prune11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z19k_nbnxn_cutoff_ener11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z20k_nbnxn_cutoff_prune11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z21k_nbnxn_rf_ener_prune11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z24k_nbnxn_ewald_ener_prune11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z25k_nbnxn_cutoff_ener_prune11cu_atomdata10cu_nbparam8cu_plisti"
};

/*! Names of CUDA legacy nonbonded kernel functions with mangling (to be used
    with CUDA 3.2/4.0 (used to be development k2). */
static const char * const nb_legacy_knames[NUM_NB_KERNELS] =
{
    "_Z17k_nbnxn_rf_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z20k_nbnxn_ewald_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z21k_nbnxn_cutoff_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z22k_nbnxn_rf_ener_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z23k_nbnxn_rf_prune_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z25k_nbnxn_ewald_ener_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z26k_nbnxn_ewald_prune_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z26k_nbnxn_cutoff_ener_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z27k_nbnxn_cutoff_prune_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z28k_nbnxn_rf_ener_prune_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z31k_nbnxn_ewald_ener_prune_legacy11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z32k_nbnxn_cutoff_ener_prune_legacy11cu_atomdata10cu_nbparam8cu_plisti"
};

/*! Names of old/unused nonbonded kernel fucntions with mangling (deprecated,
    used to be development k1). */
static const char * const nb_old_knames[NUM_NB_KERNELS] =
{
    "_Z14k_nbnxn_rf_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z17k_nbnxn_ewald_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z18k_nbnxn_cutoff_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z19k_nbnxn_rf_ener_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z20k_nbnxn_rf_prune_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z22k_nbnxn_ewald_ener_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z23k_nbnxn_ewald_prune_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z23k_nbnxn_cutoff_ener_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z24k_nbnxn_cutoff_prune_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z25k_nbnxn_rf_ener_prune_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z28k_nbnxn_ewald_ener_prune_old11cu_atomdata10cu_nbparam8cu_plisti",
    "_Z29k_nbnxn_cutoff_ener_prune_old11cu_atomdata10cu_nbparam8cu_plisti"
};


/*! Dummy kernel used for sanity check. */
__device__ __global__ void k_empty_test(){}

/*! Tabulates the Ewald Coulomb force and initializes the size/scale 
    and the table GPU array. If called with an already allocated table,
    it just re-uploads the table.
 */
static void init_ewald_coulomb_force_table(cu_nbparam_t *nbp)
{
    float       *ftmp, *coul_tab;
    int         tabsize;
    double      tabscale;
    cudaError_t stat;

    tabsize     = EWALD_COULOMB_FORCE_TABLE_SIZE;
    /* Subtract 2 iso 1 to avoid access out of range due to rounding */
    tabscale    = (tabsize - 2) / sqrt(nbp->rcoulomb_sq);

    pmalloc((void**)&ftmp, tabsize*sizeof(*ftmp));

    table_spline3_fill_ewald_lr(ftmp, NULL, tabsize, tableformatF,
                                1/tabscale, nbp->ewald_beta);

    /* If the table pointer == NULL the table is generated the first time =>
       the array pointer will be saved to nbparam and the texture is bound.
     */
    coul_tab = nbp->coulomb_tab;
    if (coul_tab == NULL)
    {
        stat = cudaMalloc((void **)&coul_tab, tabsize*sizeof(*coul_tab));
        CU_RET_ERR(stat, "cudaMalloc failed on coul_tab");

        nbp->coulomb_tab = coul_tab;
        cu_bind_texture("tex_coulomb_tab", coul_tab, tabsize*sizeof(*coul_tab));
    }

    cu_copy_H2D(coul_tab, ftmp, tabsize*sizeof(*coul_tab));

    nbp->coulomb_tab_size     = tabsize;
    nbp->coulomb_tab_scale    = tabscale;

    pfree(ftmp);
}


/*! Initilizes the atomdata structure first time, it only gets filled at 
    pair-search. */
static void init_atomdata_first(cu_atomdata_t *ad, int ntypes)
{
    cudaError_t stat;

    ad->ntypes  = ntypes;
    stat = cudaMalloc((void**)&ad->shift_vec, SHIFTS*sizeof(*ad->shift_vec));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->shift_vec"); 
    ad->shift_vec_uploaded = FALSE;

    stat = cudaMalloc((void**)&ad->fshift, SHIFTS*sizeof(*ad->fshift));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->fshift");

    stat = cudaMalloc((void**)&ad->e_lj, sizeof(*ad->e_lj));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->e_lj");
    stat = cudaMalloc((void**)&ad->e_el, sizeof(*ad->e_el));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->e_el");

    /* initilize to NULL poiters to data that is not allocated here and will
       need reallocation in nbnxn_cuda_init_atomdata */
    ad->xq = NULL;
    ad->f  = NULL;

    /* size -1 indicates that the repective array hasn't been initialized yet */
    ad->natoms = -1;
    ad->nalloc = -1;
}

/*! Initilizes the nonbonded parameter data structure. */
static void init_nbparam(cu_nbparam_t *nbp,
                           const interaction_const_t *ic,
                           const nonbonded_verlet_t *nbv)
{  
    cudaError_t stat;
    int         ntypes, nnbfp; 

    ntypes  = nbv->grp[0].nbat->ntype;
    
    nbp->ewald_beta = ic->ewaldcoeff;
    nbp->sh_ewald   = ic->sh_ewald;
    nbp->epsfac     = ic->epsfac;
    nbp->two_k_rf   = 2.0 * ic->k_rf;
    nbp->c_rf       = ic->c_rf;
    nbp->rvdw_sq    = ic->rvdw * ic->rvdw;
    nbp->rcoulomb_sq= ic->rcoulomb * ic->rcoulomb;
    nbp->rlist_sq   = ic->rlist * ic->rlist;
    nbp->sh_invrc6  = ic->sh_invrc6;

    if (ic->eeltype == eelCUT)
    {
        nbp->eeltype = cu_eelCUT;
    }
    else if (EEL_RF(ic->eeltype))
    {                
        nbp->eeltype = cu_eelRF;
    }
    else if ((EEL_PME(ic->eeltype) || ic->eeltype==eelEWALD))
    {
        nbp->eeltype = cu_eelEWALD;
    }
    else 
    {
        gmx_fatal(FARGS, "The requested electrostatics type is not implemented in the CUDA GPU accelerated kernels!");
    }

    /* generate table for PME */
    if (nbp->eeltype == cu_eelEWALD)
    {
        nbp->coulomb_tab = NULL;
        init_ewald_coulomb_force_table(nbp);
    }

    nnbfp = 2*ntypes*ntypes;
    stat = cudaMalloc((void **)&nbp->nbfp, nnbfp*sizeof(*nbp->nbfp));
    CU_RET_ERR(stat, "cudaMalloc failed on nbp->nbfp"); 
    cu_copy_H2D(nbp->nbfp, nbv->grp[0].nbat->nbfp, nnbfp*sizeof(*nbp->nbfp));
    cu_bind_texture("tex_nbfp", nbp->nbfp, nnbfp*sizeof(*nbp->nbfp));
}

void reset_gpu_rlist_ewaldtab(nbnxn_cuda_ptr_t cu_nb,
                              const interaction_const_t *ic)
{
    cu_nbparam_t *nbp = cu_nb->nbparam;

    nbp->rlist_sq       = ic->rlist * ic->rlist;
    nbp->rcoulomb_sq    = ic->rcoulomb * ic->rcoulomb;
    nbp->ewald_beta     = ic->ewaldcoeff;

    init_ewald_coulomb_force_table(cu_nb->nbparam);
}

/*! Initilizes the pair list data structure. */
static void init_plist(cu_plist_t *pl)
{
    /* initilize to NULL poiters to data that is not allocated here and will
       need reallocation in nbnxn_cuda_init_pairlist */
    pl->sci     = NULL;
    pl->cj4     = NULL;
    pl->excl    = NULL;
    
    /* size -1 indicates that the repective array hasn't been initialized yet */
    pl->na_c        = -1;
    pl->nsci        = -1;
    pl->sci_nalloc  = -1;
    pl->ncj4        = -1;
    pl->cj4_nalloc  = -1;
    pl->nexcl       = -1;
    pl->excl_nalloc = -1;
    pl->do_prune    = FALSE;
}

/*! Initilizes the timer data structure. */
static void init_timers(cu_timers_t *t, gmx_bool bDomDec)
{
    cudaError_t stat;
    int eventflags = ( USE_CUDA_EVENT_BLOCKING_SYNC ? cudaEventBlockingSync: cudaEventDefault );

    stat = cudaEventCreateWithFlags(&(t->start_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_atdat failed");
    stat = cudaEventCreateWithFlags(&(t->stop_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_atdat failed");

    /* The non-local counters/stream (second in the array) are needed only with DD. */
    for (int i = 0; i <= bDomDec ? 1 : 0; i++)
    {
        stat = cudaEventCreateWithFlags(&(t->start_nb_k[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_k failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_k[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_k failed");


        stat = cudaEventCreateWithFlags(&(t->start_pl_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_pl_h2d failed");
        stat = cudaEventCreateWithFlags(&(t->stop_pl_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_pl_h2d failed");

        stat = cudaEventCreateWithFlags(&(t->start_nb_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_h2d failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_h2d failed");

        stat = cudaEventCreateWithFlags(&(t->start_nb_d2h[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_d2h failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_d2h[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_d2h failed");
    }
}

/*! Initilizes the timings data structure. */
static void init_timings(wallclock_gpu_t *t)
{
    int i, j;

    t->nb_h2d_t = 0.0;
    t->nb_d2h_t = 0.0;
    t->nb_c    = 0;
    t->pl_h2d_t = 0.0;
    t->pl_h2d_c = 0;
    for (i = 0; i < 2; i++)
    {
        for(j = 0; j < 2; j++)
        {
            t->ktime[i][j].t = 0.0;
            t->ktime[i][j].c = 0;
        }
    }
}

void nbnxn_cuda_init(FILE *fplog,
                     nbnxn_cuda_ptr_t *p_cu_nb,
                     gmx_bool bDomDec)
{
    cudaError_t stat;
    nbnxn_cuda_ptr_t  nb;
    char sbuf[STRLEN];
    bool bOldKernel, bLegacyKernel, bDefaultKernel, bCUDA40, bCUDA32;

    if (p_cu_nb == NULL) return;

    snew(nb, 1);
    snew(nb->dev_info, 1);
    snew(nb->atdat, 1);
    snew(nb->nbparam, 1);
    snew(nb->plist[eintLocal], 1);
    if (bDomDec)
    {
        snew(nb->plist[eintNonlocal], 1);
    }

    nb->dd_run = bDomDec;

    snew(nb->timers, 1);
    snew(nb->timings, 1);

    /* init nbst */
    pmalloc((void**)&nb->nbst.e_lj, sizeof(*nb->nbst.e_lj));
    pmalloc((void**)&nb->nbst.e_el, sizeof(*nb->nbst.e_el));
    pmalloc((void**)&nb->nbst.fshift, SHIFTS * sizeof(*nb->nbst.fshift));

    init_plist(nb->plist[eintLocal]);

    /* local/non-local GPU streams */
    stat = cudaStreamCreate(&nb->stream[eintLocal]);
    CU_RET_ERR(stat, "cudaStreamCreate on stream[eintLocal] failed");
    if (bDomDec)
    {
        init_plist(nb->plist[eintNonlocal]);
        stat = cudaStreamCreate(&nb->stream[eintNonlocal]);
        CU_RET_ERR(stat, "cudaStreamCreate on stream[eintNonlocal] failed");
    }

    /* init events for sychronization (timing disabled for performance reasons!) */
    stat = cudaEventCreateWithFlags(&nb->nonlocal_done, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on nonlocal_done failed");
    stat = cudaEventCreateWithFlags(&nb->misc_ops_done, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on misc_ops_one failed");

    /* init device info */
    /* FIXME this should not be done here! */
    stat = cudaGetDevice(&nb->dev_info->dev_id);
    CU_RET_ERR(stat, "cudaGetDevice failed");
    stat = cudaGetDeviceProperties(&nb->dev_info->dev_prop, nb->dev_info->dev_id);
    CU_RET_ERR(stat, "cudaGetDeviceProperties failed");

    /* If ECC is enabled cudaStreamSynchronize introduces huge idling so we'll 
       switch to the (atmittedly fragile) memory polling waiting to preserve 
       performance. Event-timing also needs to be disabled. */

    gmx_bool bStreamSync    = getenv("GMX_CUDA_STREAMSYNC") != NULL;
    gmx_bool bNoStreamSync  = getenv("GMX_NO_CUDA_STREAMSYNC") != NULL;
    if (bStreamSync && bNoStreamSync)
    {
        gmx_fatal(FARGS, "Conflicting environment variables: both GMX_CUDA_STREAMSYNC and GMX_NO_CUDA_STREAMSYNC defined");
    }

    /* if we can't yield while poll-waiting, we need to use cudaStreamSynchronize */
    if (nb->dev_info->dev_prop.ECCEnabled == 1)
    {
        if (bStreamSync)
        {
            nb->use_stream_sync = TRUE;

            sprintf(sbuf, "NOTE: cudaStreamSynchronize-based waiting forced by GMX_CUDA_STREAMSYNC, but ECC is turned\n"
                    "      on which generally causes considerable performance loss");
            fprintf(stderr, "\n%s\n", sbuf);
            if (fplog)
            {
                fprintf(fplog, "\n%s\n", sbuf);
            }
        }
        else
        {
#ifdef CAN_CUTHREAD_YIELD
            nb->use_stream_sync = FALSE;

            sprintf(sbuf, "NOTE: running on a GPU with ECC on; will not use cudaStreamSynchronize-based\n"
                    "      waiting as it generally causes performance loss when used with ECC.");
            fprintf(stderr, "\n%s\n", sbuf);
            if (fplog)
            {
                fprintf(fplog, "\n%s\n", sbuf);
            }
#else
            nb->use_stream_sync = TRUE;

            sprintf(sbuf, "Can't do thread yield as sleep(0)/Sleep(0) (posix/win) is not available.\n"
                        "         Will use the standard cudaStreamSynchronize-based waiting.\n"
                        "         Note that with ECC enabled this causes considerable performance loss.");
            gmx_warning("%s\n", sbuf);
            if (fplog)
            {
                fprintf(fplog, "\n%s\n", sbuf);
            }
#endif /* CAN_CUTHREAD_YIELD */
        }
    }
    else
    {
        if (bNoStreamSync)
        {
            nb->use_stream_sync = FALSE;

            sprintf(sbuf, "NOTE: running on a GPU with no ECC, but cudaStreamSynchronize-based waiting\n"
                    "forced off by GMX_NO_CUDA_STREAMSYNC");
            fprintf(stderr, "\n%s\n", sbuf);
            if (fplog)
            {
                fprintf(fplog, "\n%s\n", sbuf);
            }
        }
        else
        {
            /* no/off ECC, cudaStreamSynchronize not turned off by env. var. */
            nb->use_stream_sync = TRUE;
        }
    }

    /* CUDA timing disabled as event timers don't work:
       - with multiple streams = domain-decomposition;
       - with the waiting hack, without cudaStreamSynchronize;
       - when turned off by GMX_DISABLE_CUDA_TIMING.
     */
    nb->do_time = (!bDomDec && nb->use_stream_sync && (getenv("GMX_DISABLE_CUDA_TIMING") == NULL));

    if (nb->do_time)
    {
        init_timers(nb->timers, bDomDec);
        init_timings(nb->timings);
    }

    /* FIXME move this out in a function
       Decide which kernel version to run based on:
       - GPU SM version
       - non-bonded kernel selector environment variables (both new and old form for backward-compatibility).
    */
    /* old kernel (former k1) */
    bOldKernel     = (getenv("GMX_NB_K1") != NULL) || (getenv("GMX_CU_NB_OLD") != NULL);
    /* legacy kernel (former k2), kept for now for backward compatiblity,
       faster than the default with  CUDA 3.2/4.0 (TODO: on Kepler?). */
    bLegacyKernel  = (getenv("GMX_NB_K2") != NULL) || (getenv("GMX_CU_NB_LEGACY") != NULL);
    /* default kernel (former k3). */
    bDefaultKernel = (getenv("GMX_NB_K3") != NULL) || (getenv("GMX_CU_NB_DEFAULT") != NULL);

    if ((unsigned)(bOldKernel + bLegacyKernel + bDefaultKernel) > 1)
    {
        /* TODO remove in release? */
        gmx_fatal(FARGS, "Multiple CUDA non-bonded kernels requested; to manually pick a kernel set only one \n"
                  "of the following environment variables: \n"
                  "GMX_CU_NB_DEFAULT, GMX_CU_NB_LEGACY [GMX_CU_NB_OLD deprecated!]");
    }

    /* set the kernel type for the current GPU */
    if (bOldKernel)
    {
        nb->kernel_ver = eNbnxnCuKOld;
        /* TODO remove in release? */
        gmx_warning("The old (and deprecated) CUDA non-bonded kernels (former k1) were selected.\n"
                    "      This kernel is always slower than the default and legacy kernels.");
    }
    else
    {
        char sbuf[STRLEN];

        bCUDA32 = bCUDA40 = false;
#if CUDA_VERSION == 3200
        bCUDA32 = true;
        sprintf(sbuf, "3.2");
#elif CUDA_VERSION == 4000
        bCUDA40 = false;
        sprintf(sbuf, "4.0");
#endif

        /* default is default ;) */
        nb->kernel_ver = eNbnxnCuKDefault;

        /* TODO: add warning to the sanity checks for the case when code compiled
                 with old nvcc is used on Kepler. */
        if (bCUDA32 || bCUDA40)
        {
            /* use legacy kernel unless something else is forced by an env. var */
            if (bDefaultKernel)
            {
                fprintf(stderr,
                        "\nNOTE: CUDA %s compilation detected; with this compiler version the legacy\n"
                        "      non-bonded kenels perform best. However, the default kernels were\n"
                        "      selected by the GMX_CU_NB_DEFAULT environment variable."
                        "      Consider using using the legacy kernels or upgrade your CUDA tookit.\n",
                        sbuf);
            }
            else
            {
                nb->kernel_ver = eNbnxnCuKLegacy;
            }
        }
        else
        {
            /* issue not if the non-default kernel is forced by an env. var */
            if (bLegacyKernel)
            {
                fprintf(stderr, "\nNOTE: Non-default non-bonded CUDA kernels were selected by the GMX_CU_NB_LEGACY\n"
                        "      env var. Consider using using the default kernels which should be faster!\n");

                nb->kernel_ver = eNbnxnCuKLegacy;
            }
        }
    }


    for (int i = 0; i < NUM_NB_KERNELS; i++)
    {
        /* Legacy kernel 16/48 kB Shared/L1 */
        stat = cudaFuncSetCacheConfig(nb_legacy_knames[i], cudaFuncCachePreferL1);
        CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");

        if (nb->dev_info->dev_prop.major >= 3)
        {
            /* Default kernel on sm 3.x 48/16 kB Shared/L1 */
            stat = cudaFuncSetCacheConfig(nb_default_knames[i], cudaFuncCachePreferShared);
        }
        else
        {
            /* On Fermi prefer L1 gives 2% higher performance */
            /* Default kernel on sm_2.x 16/48 kB Shared/L1 */
            stat = cudaFuncSetCacheConfig(nb_default_knames[i], cudaFuncCachePreferL1);
        }
        CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");

        /* Old kernel 48/16 kB Shared/L1 (TODO deprecated, remove) */
        stat = cudaFuncSetCacheConfig(nb_old_knames[i],  cudaFuncCachePreferShared);
        CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
    }

    *p_cu_nb = nb;

    if (debug)
    {
        fprintf(debug, "Initialized CUDA data structures.\n");
    }


    /* TODO: move this to gpu_utils module */
    k_empty_test<<<1, 512>>>();
    CU_LAUNCH_ERR_SYNC("dummy test kernel");
}

void nbnxn_cuda_init_const(nbnxn_cuda_ptr_t cu_nb,
                           const interaction_const_t *ic,
                           const nonbonded_verlet_t *nbv)
{
    init_atomdata_first(cu_nb->atdat, nbv->grp[0].nbat->ntype);
    init_nbparam(cu_nb->nbparam, ic, nbv);

    /* clear energy and shift force outputs */
    nbnxn_cuda_clear_e_fshift(cu_nb);
}

void nbnxn_cuda_init_pairlist(nbnxn_cuda_ptr_t cu_nb,
                              const nbnxn_pairlist_t *h_plist,
                              int iloc)
{
    char         sbuf[STRLEN];
    cudaError_t  stat;
    gmx_bool     do_time    = cu_nb->do_time;
    cudaStream_t stream     = cu_nb->stream[iloc];
    cu_plist_t  *d_plist    = cu_nb->plist[iloc];

    if (d_plist->na_c < 0)
    {
        d_plist->na_c = h_plist->na_ci;
    }
    else
    {
        if (d_plist->na_c != h_plist->na_ci)
        {
            sprintf(sbuf, "In cu_init_plist: the #atoms per cell has changed (from %d to %d)",
                    d_plist->na_c, h_plist->na_ci);
            gmx_incons(sbuf);
        }
    }

    if (do_time)
    {
        stat = cudaEventRecord(cu_nb->timers->start_pl_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    cu_realloc_buffered((void **)&d_plist->sci, h_plist->sci, sizeof(*d_plist->sci),
                         &d_plist->nsci, &d_plist->sci_nalloc,
                         h_plist->nsci,
                         stream, TRUE);

    cu_realloc_buffered((void **)&d_plist->cj4, h_plist->cj4, sizeof(*d_plist->cj4),
                         &d_plist->ncj4, &d_plist->cj4_nalloc,
                         h_plist->ncj4,
                         stream, TRUE);

    cu_realloc_buffered((void **)&d_plist->excl, h_plist->excl, sizeof(*d_plist->excl),
                         &d_plist->nexcl, &d_plist->excl_nalloc,
                         h_plist->nexcl,
                         stream, TRUE);

    if (do_time)
    {
        stat = cudaEventRecord(cu_nb->timers->stop_pl_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* need to prune the pair list during the next step */
    d_plist->do_prune = TRUE;
}

void nbnxn_cuda_upload_shiftvec(nbnxn_cuda_ptr_t cu_nb,
                                const nbnxn_atomdata_t *nbatom)
{
    cu_atomdata_t *adat = cu_nb->atdat;
    cudaStream_t  ls    = cu_nb->stream[eintLocal];

    /* only if we have a dynamic box */
    if (nbatom->dynamic_box || !adat->shift_vec_uploaded)
    {
        cu_copy_H2D_async(adat->shift_vec, nbatom->shift_vec, 
                          SHIFTS * sizeof(*adat->shift_vec), ls);
        adat->shift_vec_uploaded = TRUE;
    }
}

/*! Clears the first natoms_clear elements of the GPU nonbonded force output array. */
static void nbnxn_cuda_clear_f(nbnxn_cuda_ptr_t cu_nb, int natoms_clear)
{
    cudaError_t   stat;
    cu_atomdata_t *adat = cu_nb->atdat;
    cudaStream_t  ls    = cu_nb->stream[eintLocal];

    stat = cudaMemsetAsync(adat->f, 0, natoms_clear * sizeof(*adat->f), ls);
    CU_RET_ERR(stat, "cudaMemsetAsync on f falied");
}

/*! Clears nonbonded shift force output array and energy outputs on the GPU. */
static void nbnxn_cuda_clear_e_fshift(nbnxn_cuda_ptr_t cu_nb)
{
    cudaError_t   stat;
    cu_atomdata_t *adat = cu_nb->atdat;
    cudaStream_t  ls    = cu_nb->stream[eintLocal];

    stat = cudaMemsetAsync(adat->fshift, 0, SHIFTS * sizeof(*adat->fshift), ls);
    CU_RET_ERR(stat, "cudaMemsetAsync on fshift falied");
    stat = cudaMemsetAsync(adat->e_lj, 0, sizeof(*adat->e_lj), ls);
    CU_RET_ERR(stat, "cudaMemsetAsync on e_lj falied");
    stat = cudaMemsetAsync(adat->e_el, 0, sizeof(*adat->e_el), ls);
    CU_RET_ERR(stat, "cudaMemsetAsync on e_el falied");
}

void nbnxn_cuda_clear_outputs(nbnxn_cuda_ptr_t cu_nb, int flags)
{
    nbnxn_cuda_clear_f(cu_nb, cu_nb->atdat->natoms + 1);
    /* clear shift force array and energies if the outputs were 
       used in the current step */
    if (flags & GMX_FORCE_VIRIAL)
    {
        nbnxn_cuda_clear_e_fshift(cu_nb);
    }
}

void nbnxn_cuda_init_atomdata(nbnxn_cuda_ptr_t cu_nb,
                              const nbnxn_atomdata_t *nbat)
{
    cudaError_t   stat;
    int           nalloc, natoms;
    gmx_bool      realloced;
    gmx_bool      do_time   = cu_nb->do_time;
    cu_timers_t   *timers   = cu_nb->timers;
    cu_atomdata_t *d_atdat  = cu_nb->atdat;
    cudaStream_t  ls        = cu_nb->stream[eintLocal];

    natoms = nbat->natoms;
    realloced = FALSE;

    if (do_time)
    {
        /* time async copy */
        stat = cudaEventRecord(timers->start_atdat, ls);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* need to reallocate if we have to copy more atoms than the amount of space
       available and only allocate if we haven't initilzed yet, i.e d_atdat->natoms == -1 */
    if (natoms + 1 > d_atdat->nalloc)
    {
        nalloc = over_alloc_small(natoms);
    
        /* free up first if the arrays have already been initialized */
        if (d_atdat->nalloc != -1)
        {
            cu_free_buffered(d_atdat->f, &d_atdat->natoms, &d_atdat->nalloc);
            cu_free_buffered(d_atdat->xq);
            cu_free_buffered(d_atdat->atom_types);
        }
        
        /* Add one element to f for polling-based wait */
        stat = cudaMalloc((void **)&d_atdat->f, (nalloc + 1)*sizeof(*d_atdat->f));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->f");
        stat = cudaMalloc((void **)&d_atdat->xq, nalloc*sizeof(*d_atdat->xq));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->xq");

        stat = cudaMalloc((void **)&d_atdat->atom_types, nalloc*sizeof(*d_atdat->atom_types));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->atom_types");

        d_atdat->nalloc = nalloc;
        realloced = TRUE;
    }

    d_atdat->natoms = natoms;
    d_atdat->natoms_local = nbat->natoms_local;

    /* need to clear GPU f output if realloc happened */
    if (realloced)
    {
        nbnxn_cuda_clear_f(cu_nb, nalloc);
    }

    cu_copy_H2D_async(d_atdat->atom_types, nbat->type,
                      natoms*sizeof(*d_atdat->atom_types), ls);

    if (do_time)
    {
        stat = cudaEventRecord(timers->stop_atdat, ls);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

void nbnxn_cuda_free(FILE *fplog, nbnxn_cuda_ptr_t cu_nb, gmx_bool bDomDec)
{
    cudaError_t     stat;
    cu_atomdata_t   *atdat;
    cu_nbparam_t    *nbparam;
    cu_plist_t      *plist, *plist_nl;
    cu_timers_t     *timers;

    if (cu_nb == NULL) return;

    atdat       = cu_nb->atdat;
    nbparam     = cu_nb->nbparam;
    plist       = cu_nb->plist[eintLocal];
    plist_nl    = cu_nb->plist[eintNonlocal];
    timers      = cu_nb->timers;

    if (nbparam->eeltype == cu_eelEWALD)
    {
        cu_unbind_texture("tex_coulomb_tab");
        cu_free_buffered(nbparam->coulomb_tab, &nbparam->coulomb_tab_size);
    }

    stat = cudaEventDestroy(cu_nb->nonlocal_done);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->nonlocal_done");
    stat = cudaEventDestroy(cu_nb->misc_ops_done);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->misc_ops_done");

    if (cu_nb->do_time)
    {
        stat = cudaEventDestroy(timers->start_atdat);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_atdat");
        stat = cudaEventDestroy(timers->stop_atdat);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_atdat");

        /* The non-local counters/stream (second in the array) are needed only with DD. */
        for (int i = 0; i <= bDomDec ? 1 : 0; i++)
        {
            stat = cudaEventDestroy(timers->start_nb_k[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_k");
            stat = cudaEventDestroy(timers->stop_nb_k[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_k");

            stat = cudaEventDestroy(timers->start_pl_h2d[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_pl_h2d");
            stat = cudaEventDestroy(timers->stop_pl_h2d[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_pl_h2d");

            stat = cudaStreamDestroy(cu_nb->stream[i]);
            CU_RET_ERR(stat, "cudaStreamDestroy failed on stream");

            stat = cudaEventDestroy(timers->start_nb_h2d[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_h2d");
            stat = cudaEventDestroy(timers->stop_nb_h2d[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_h2d");

            stat = cudaEventDestroy(timers->start_nb_d2h[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_d2h");
            stat = cudaEventDestroy(timers->stop_nb_d2h[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_d2h");
        }
    }

    cu_unbind_texture("tex_nbfp");
    cu_free_buffered(nbparam->nbfp);

    stat = cudaFree(atdat->shift_vec);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atdat->shift_vec");
    stat = cudaFree(atdat->fshift);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atdat->fshift");

    stat = cudaFree(atdat->e_lj);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atdat->e_lj");
    stat = cudaFree(atdat->e_el);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atdat->e_el");

    cu_free_buffered(atdat->f, &atdat->natoms, &atdat->nalloc);
    cu_free_buffered(atdat->xq);
    cu_free_buffered(atdat->atom_types, &atdat->ntypes);

    cu_free_buffered(plist->sci, &plist->nsci, &plist->sci_nalloc);
    cu_free_buffered(plist->cj4, &plist->ncj4, &plist->cj4_nalloc);
    cu_free_buffered(plist->excl, &plist->nexcl, &plist->excl_nalloc);
    if (bDomDec)
    {
        cu_free_buffered(plist_nl->sci, &plist_nl->nsci, &plist_nl->sci_nalloc);
        cu_free_buffered(plist_nl->cj4, &plist_nl->ncj4, &plist_nl->cj4_nalloc);
        cu_free_buffered(plist_nl->excl, &plist_nl->nexcl, &plist->excl_nalloc);
    }

    stat = cudaThreadExit();
    CU_RET_ERR(stat, "cudaThreadExit failed");

    if (debug)
    {
        fprintf(debug, "Cleaned up CUDA data structures.\n");
    }
}

void cu_synchstream_atdat(nbnxn_cuda_ptr_t cu_nb, int iloc)
{
    cudaError_t stat;
    cudaStream_t stream = cu_nb->stream[iloc];

    stat = cudaStreamWaitEvent(stream, cu_nb->timers->stop_atdat, 0);
    CU_RET_ERR(stat, "cudaStreamWaitEvent failed");
}

wallclock_gpu_t * nbnxn_cuda_get_timings(nbnxn_cuda_ptr_t cu_nb)
{
    return (cu_nb != NULL && cu_nb->do_time) ? cu_nb->timings : NULL;
}

void nbnxn_cuda_reset_timings(nbnxn_cuda_ptr_t cu_nb)
{
    if (cu_nb->do_time)
    {
        init_timings(cu_nb->timings);
    }
}

int nbnxn_cuda_min_ci_balanced(nbnxn_cuda_ptr_t cu_nb)
{
    return cu_nb != NULL ? 
        GPU_MIN_CI_BALANCED_FACTOR*cu_nb->dev_info->dev_prop.multiProcessorCount : 0;
}


/****** FIXME Old stuff, mostly deprecated, remove before release  *****/

/* Upload asynchronously to the GPU the coordinate+charge array.
 * XXX not used  
 */
void cu_move_xq(nbnxn_cuda_ptr_t cu_nb, const nbnxn_atomdata_t *nbat, int aloc)
{
    int iloc = -1; 

    /* determine interaction locality from atom locality 
       (needed for indexing timers/streams) */
    if (LOCAL_A(aloc))
    {
        iloc = eintLocal;
    }
    else if (NONLOCAL_A(aloc))
    {
        iloc = eintNonlocal;
    }
    else 
    {
        gmx_incons("Invalid atom locality passed (valid here is only local or nonlocal)");
    }

    cu_atomdata_t   *d_nbat = cu_nb->atdat;
    cudaStream_t    stream  = cu_nb->stream[iloc];

    cu_copy_H2D_async(d_nbat->xq, nbat->x,
                        d_nbat->natoms * sizeof(*d_nbat->xq), stream);
}

/*! Waits until the atom data gets copied to the GPU and times the transfer.
 *  XXX not used  
 */
void cu_wait_atdat(nbnxn_cuda_ptr_t cu_nb)
{
    float t;
    cu_wait_event_time(cu_nb->timers->stop_atdat, cu_nb->timers->start_atdat, &t);
    cu_nb->timings->pl_h2d_t += t;
}
