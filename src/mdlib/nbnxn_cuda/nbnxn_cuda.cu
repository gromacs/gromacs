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

#if defined(_MSVC)
#include <limits>
#endif

#include "smalloc.h"
#include "types/simple.h" 
#include "types/nbnxn_pairlist.h"
#include "types/nb_verlet.h"
#include "types/ishift.h"
#include "types/force_flags.h"

#include "nbnxn_cuda_types.h"
#include "../../gmxlib/cuda_tools/cudautils.cuh"
#include "nbnxn_cuda.h"
#include "nbnxn_cuda_data_mgmt.h"
#include "pmalloc_cuda.h"

/* we can do yield which meaning that we have all the required posix/Win stuff */
#ifdef CAN_CUTHREAD_YIELD
#if (defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__
/* On Windows we assume that we always have windows.h */
#include <windows.h>
#else
/* Posix -- if we got here it means that unistd.h is available, but we'll
   check again to catch future stupid modifications. */
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#endif
#endif /* CAN_CUTHREAD_YIELD */


/***** The kernels come here *****/

#define CLUSTER_SIZE            (NBNXN_GPU_CLUSTER_SIZE)

#include "../../gmxlib/cuda_tools/vectype_ops.cuh"
#include "nbnxn_cuda_kernel_utils.cuh"

/* Generate all combinations of kernels through multiple inclusion:
   F, F + E, F + prune, F + E + prune. */
/** Force only **/
#include "nbnxn_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "nbnxn_cuda_kernels.cuh"
#undef CALC_ENERGIES

/*** Pair-list pruning kernels ***/
/** Force only **/
#define PRUNE_NBL
#include "nbnxn_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "nbnxn_cuda_kernels.cuh"
#undef CALC_ENERGIES
#undef PRUNE_NBL

/*! Nonbonded kernel function pointer type */
typedef void (*p_k_nbnxn)(const cu_atomdata_t,
                          const cu_nbparam_t,
                          const cu_plist_t,
                          gmx_bool);

/*********************************/

/* XXX always/never run the energy/pruning kernels -- only for benchmarking purposes */
static gmx_bool always_ener  = (getenv("GMX_GPU_ALWAYS_ENER") != NULL);
static gmx_bool never_ener   = (getenv("GMX_GPU_NEVER_ENER") != NULL);
static gmx_bool always_prune = (getenv("GMX_GPU_ALWAYS_PRUNE") != NULL);

/* Non-bonded kernel selector environment variables */
/* old kernel (former k1) */
static gmx_bool bRunOldKernel     = (getenv("GMX_NB_K1") != NULL) || (getenv("GMX_CU_NB_OLD") != NULL);
/* legacy kernel (former #2) -- kept for now for backward compatiblity. */
static gmx_bool bRunLegacyKernel  = (getenv("GMX_NB_K2") != NULL) || (getenv("GMX_CU_NB_LEGACY") != NULL);
/* default kernel (former #3). */
static gmx_bool bRunDefaultKernel = (getenv("GMX_NB_K3") != NULL) || (getenv("GMX_CU_NB_DEFAULT") != NULL);


/* Single byte of the bit-pattern used to fill the last element of the force
   array when using polling-based waiting. */
static unsigned char gpu_sync_signal_byte    = 0xAA;
/* Bit-pattern used to fill the last element of the force with when using
   polling-based waiting. */
static unsigned int  gpu_sync_signal_pattern = 0xAAAAAAAA;


/*! Returns the number of blocks to be used  for the nonbonded GPU kernel. */
static inline int calc_nb_blocknr(int nwork_units)
{
    int retval = (nwork_units <= GRID_MAX_DIM ? nwork_units : GRID_MAX_DIM);
    if (retval != nwork_units)
    {
        gmx_fatal(FARGS, "Watch out, the number of nonbonded work units exceeds ",
                  "the maximum grid size (%d > %d)!",
                  nwork_units, GRID_MAX_DIM);
    }
    return retval;
}

/*! Selects the kernel version to execute at the current step and 
 *  returns a function pointer to it. 
 */
static inline p_k_nbnxn select_nbnxn_kernel(int eeltype, gmx_bool doEne, 
                                            gmx_bool doPrune)
{
    p_k_nbnxn k = NULL;

    if (bRunOldKernel + bRunLegacyKernel + bRunDefaultKernel > 1)
    {
        /* TODO remove in release */
        gmx_fatal(FARGS, "Multiple GPU NB kernels requested at the same time; "
                  "set only one of the GMX_NB_K[123] env vars!");
    }

    /* select which kernel will be used */
    switch (eeltype)
    {
        case cu_eelCUT:
            if (bRunOldKernel)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_cutoff_old :
                        k_nbnxn_cutoff_prune_old;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_cutoff_ener_old :
                        k_nbnxn_cutoff_ener_prune_old;
                }
            }
            else if (bRunLegacyKernel)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_cutoff_legacy :
                                   k_nbnxn_cutoff_prune_legacy;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_cutoff_ener_legacy :
                                   k_nbnxn_cutoff_ener_prune_legacy;
                }
            }
            else
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_cutoff :
                                   k_nbnxn_cutoff_prune;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_cutoff_ener :
                                   k_nbnxn_cutoff_ener_prune;
                }
            }
            break;

        case cu_eelRF:
            if (bRunOldKernel)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_rf_old :
                        k_nbnxn_rf_prune_old;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_rf_ener_old :
                        k_nbnxn_rf_ener_prune_old;
                }
            }

            if (bRunLegacyKernel)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_rf_legacy :
                                   k_nbnxn_rf_prune_legacy;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_rf_ener_legacy :
                                   k_nbnxn_rf_ener_prune_legacy;
                }
            }
            else if (bRunDefaultKernel)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_rf :
                                   k_nbnxn_rf_prune;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_rf_ener :
                                   k_nbnxn_rf_ener_prune;
                }
            }
            break;

        case cu_eelEWALD:
            if (bRunOldKernel)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_ewald_old :
                        k_nbnxn_ewald_prune_old;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_ewald_ener_old :
                        k_nbnxn_ewald_ener_prune_old;
                }
            }
            else if (bRunLegacyKernel)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_ewald_legacy :
                                   k_nbnxn_ewald_prune_legacy;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_ewald_ener_legacy:
                                   k_nbnxn_ewald_ener_prune_legacy;
                }
            }
            else
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_ewald :
                                   k_nbnxn_ewald_prune;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_ewald_ener :
                                   k_nbnxn_ewald_ener_prune;
                }
            }
            break;

        default:
            gmx_incons("The provided electrostatics type does not exist in the "
                       "CUDA implementation!");
    }
    return k;
}

/*! As we execute nonbonded workload in separate streams, before launching 
   the kernel we need to make sure that he following operations have completed:
   - atomata allocation and related H2D transfers (every nstlist step);
   - pair list H2D transfer (every nstlist step);
   - shift vector H2D transfer (every nstlist step);
   - force (+shift force and energy) output clearing (every step);

   In CUDA all operations issued in stream 0 precede everything launched in 
   other streams. As all the above operations are launched in sream 0 (default), 
   while all nonbonded calculations are done in the respective local/non-local 
   stream, no explicit synchronization is required.
 */
void nbnxn_cuda_launch_kernel(nbnxn_cuda_ptr_t cu_nb,
                              const nbnxn_atomdata_t *nbatom,
                              int flags,
                              int iloc)
{
    cudaError_t stat;
    int adat_begin, adat_len;  /* local/nonlocal offset and length used for xq and f */

    cu_atomdata_t   *adat   = cu_nb->atdat;
    cu_nbparam_t    *nbp    = cu_nb->nbparam;
    cu_plist_t      *plist  = cu_nb->plist[iloc];
    cu_timers_t     *t      = cu_nb->timers;
    cudaStream_t    stream  = cu_nb->stream[iloc];

    /* kernel lunch-related stuff */
    p_k_nbnxn   nb_kernel = NULL; /* fn pointer to the nonbonded kernel */
    int         shmem;
    int         nb_blocks = calc_nb_blocknr(plist->nsci);
    dim3        dim_block(CLUSTER_SIZE, CLUSTER_SIZE, 1);
    dim3        dim_grid(nb_blocks, 1, 1);

    gmx_bool calc_ener   = flags & GMX_FORCE_VIRIAL;
    gmx_bool calc_fshift = flags & GMX_FORCE_VIRIAL;
    gmx_bool do_time     = cu_nb->do_time;

    /* turn energy calculation always on/off (for debugging/testing only) */
    calc_ener = (calc_ener || always_ener) && !never_ener; 

    /* don't launch the kernel if there is no work to do */
    if (plist->nsci == 0)
    {
        return;
    }

    /* calculate the atom data index range based on locality */
    if (LOCAL_I(iloc))
    {
        adat_begin  = 0;
        adat_len    = adat->natoms_local;
    }
    else
    {
        adat_begin  = adat->natoms_local;
        adat_len    = adat->natoms - adat->natoms_local;
    }

    if (debug)
    {
        fprintf(debug, "GPU launch configuration:\n\tThread block: %dx%dx%d\n\t"
                "Grid: %dx%d\n\t#Cells/Subcells: %d/%d (%d)\n",
                dim_block.x, dim_block.y, dim_block.z,
                dim_grid.x, dim_grid.y, plist->nsci*NSUBCELL,
                NSUBCELL, plist->na_c);
    }

    /* When we get here all misc operations issues in the local stream are done,
       so we record that in the local stream and wait for it in the nonlocal one. */
    if (cu_nb->dd_run)
    {
        if (iloc == eintLocal)
        {
            stat = cudaEventRecord(cu_nb->misc_ops_done, stream);
            CU_RET_ERR(stat, "cudaEventRecord on misc_ops_done failed");
        }
        else
        {
            stat = cudaStreamWaitEvent(stream, cu_nb->misc_ops_done, 0);
            CU_RET_ERR(stat, "cudaStreamWaitEvent on misc_ops_done failed");
        }
    }

    /* beginning of timed HtoD section */
    if (do_time)
    {
        stat = cudaEventRecord(t->start_nb_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* HtoD x, q */
    cu_copy_H2D_async(adat->xq + adat_begin, nbatom->x + adat_begin * 4,
                      adat_len * sizeof(*adat->xq), stream); 

    if (do_time)
    {
        stat = cudaEventRecord(t->stop_nb_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* beginning of timed nonbonded calculation section */
    if (do_time)
    {
        stat = cudaEventRecord(t->start_nb_k[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* size of shmem (force-buffers/xq/atom type preloading) */
    if (bRunOldKernel)
    {
        shmem = (1 + NSUBCELL) * CLUSTER_SIZE * CLUSTER_SIZE * 3 * sizeof(float);
    }
    else if (bRunLegacyKernel)
    {
        /* i-atom x+q in shared memory */
        shmem =  NSUBCELL * CLUSTER_SIZE * sizeof(float4);
        /* force reduction buffers in shared memory */
        shmem += CLUSTER_SIZE * CLUSTER_SIZE * 3 * sizeof(float);
    }
    else
    {
        /* NOTE: with the default kernel on sm3.0 we need shmem only for pre-loading */
        /* i-atom x+q in shared memory */
        shmem  = NSUBCELL * CLUSTER_SIZE * sizeof(float4);
#ifdef IATYPE_SHMEM
        /* i-atom types in shared memory */
        shmem += NSUBCELL * CLUSTER_SIZE * sizeof(int);
#endif
#if __CUDA_ARCH__ < 300
        /* force reduction buffers in shared memory */
        shmem += CLUSTER_SIZE * CLUSTER_SIZE * 3 * sizeof(float);
#endif
    }

    /* get the pointer to the kernel flavor we need to use */
    nb_kernel = select_nbnxn_kernel(nbp->eeltype, calc_ener,
                                    plist->do_prune || always_prune);

    nb_kernel<<<dim_grid, dim_block, shmem, stream>>>(*adat, *nbp, *plist,
                                                      calc_fshift);
    CU_LAUNCH_ERR("k_calc_nb");

    if (do_time)
    {
        stat = cudaEventRecord(t->stop_nb_k[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

void nbnxn_cuda_launch_cpyback(nbnxn_cuda_ptr_t cu_nb,
                               const nbnxn_atomdata_t *nbatom,
                               int flags,
                               int aloc)
{
    cudaError_t stat;
    int adat_begin, adat_len, adat_end;  /* local/nonlocal offset and length used for xq and f */
    int iloc = -1;
    float3 *signal_field;

    /* determine interaction locality from atom locality */
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
        char stmp[STRLEN];
        sprintf(stmp, "Invalid atom locality passed (%d); valid here is only "
                "local (%d) or nonlocal (%d)", aloc, eatLocal, eatNonlocal);
        gmx_incons(stmp);
    }

    cu_atomdata_t   *adat   = cu_nb->atdat;
    cu_timers_t     *t      = cu_nb->timers;
    gmx_bool        do_time = cu_nb->do_time;
    cudaStream_t    stream  = cu_nb->stream[iloc];

    gmx_bool calc_ener   = flags & GMX_FORCE_VIRIAL;
    gmx_bool calc_fshift = flags & GMX_FORCE_VIRIAL;

    /* don't launch copy-back if there was no work to do */
    if (cu_nb->plist[iloc]->nsci == 0)
    {
        return;
    }

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(aloc))
    {
        adat_begin  = 0;
        adat_len    = adat->natoms_local;
        adat_end    = cu_nb->atdat->natoms_local;
    }
    else
    {
        adat_begin  = adat->natoms_local;
        adat_len    = adat->natoms - adat->natoms_local;
        adat_end    = cu_nb->atdat->natoms;
    }

    /* beginning of timed D2H section */
    if (do_time)
    {
        stat = cudaEventRecord(t->start_nb_d2h[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    if (!cu_nb->use_stream_sync)
    {
        /* Set the z component of the extra element in the GPU force array
           (separately for l/nl parts) to a specific bit-pattern that we can check 
           for while waiting for the transfer to be done. */
        signal_field = adat->f + adat_end;
        stat = cudaMemsetAsync(&signal_field->z, gpu_sync_signal_byte,
                               sizeof(signal_field->z), stream);
        CU_RET_ERR(stat, "cudaMemsetAsync of the signal_field failed");

        /* For safety reasons set a few (5%) forces to NaN. This way even if the
           polling "hack" fails with some future NVIDIA driver we'll get a crash. */
        for (int i = adat_begin; i < 3*adat_end + 2; i += adat_len/20)
        {
#ifdef NAN
            nbatom->out[0].f[i] = NAN;
#else
#  ifdef _MSVC
            if (numeric_limits<float>::has_quiet_NaN)
            {
                nbatom->out[0].f[i] = numeric_limits<float>::quiet_NaN();
            }
            else
#  endif
            {
                nbatom->out[0].f[i] = GMX_REAL_MAX;
            }
#endif
        }

        /* Clear the bits in the force array (on the CPU) that we are going to poll. */
        nbatom->out[0].f[adat_end*3 + 2] = 0;
    }

    /* With DD the local D2H transfer can only start after the non-local 
       has been launched. */
    if (iloc == eintLocal && cu_nb->dd_run)
    {
        stat = cudaStreamWaitEvent(stream, cu_nb->nonlocal_done, 0);
        CU_RET_ERR(stat, "cudaStreamWaitEvent on nonlocal_done failed");
    }

    /* DtoH f */
    cu_copy_D2H_async(nbatom->out[0].f + adat_begin * 3, adat->f + adat_begin, 
                      (adat_len + 1) * sizeof(*adat->f), stream);

    /* After the non-local D2H is launched the nonlocal_done event can be
       recorded which signals that the local D2H can proceed. This event is not
       placed after the non-local kernel because we first need the non-local
       data back first. */
    if (iloc == eintNonlocal)
    {
        stat = cudaEventRecord(cu_nb->nonlocal_done, stream);
        CU_RET_ERR(stat, "cudaEventRecord on nonlocal_done failed");
    }

    /* only transfer energies in the local stream */
    if (LOCAL_I(iloc))
    {
        /* DtoH fshift */
        if (calc_fshift)
        {
            cu_copy_D2H_async(cu_nb->nbst.fshift, adat->fshift,
                              SHIFTS * sizeof(*cu_nb->nbst.fshift), stream);
        }

        /* DtoH energies */
        if (calc_ener)
        {
            cu_copy_D2H_async(cu_nb->nbst.e_lj, adat->e_lj,
                              sizeof(*cu_nb->nbst.e_lj), stream);
            cu_copy_D2H_async(cu_nb->nbst.e_el, adat->e_el,
                              sizeof(*cu_nb->nbst.e_el), stream);
        }
    }

    if (do_time)
    {
        stat = cudaEventRecord(t->stop_nb_d2h[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

static inline void cuthread_yield(void)
{
#ifndef CAN_CUTHREAD_YIELD
    gmx_incons("cuthread_yield called, but win/posix sleep is not available!")
#else
#if (defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__
        /* Windows */
        Sleep(0);
#else
        /* Posix */
        sleep(0);
#endif
#endif /* CAN_CUTHREAD_YIELD */
}

void nbnxn_cuda_wait_gpu(nbnxn_cuda_ptr_t cu_nb,
                         const nbnxn_atomdata_t *nbatom,
                         int flags, int aloc,
                         float *e_lj, float *e_el, rvec *fshift)
{
    cudaError_t stat;
    int i, adat_end, iloc = -1;
    volatile unsigned int *signal_bytes;

    /* determine interaction locality from atom locality */
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
        char stmp[STRLEN];
        sprintf(stmp, "Invalid atom locality passed (%d); valid here is only "
                "local (%d) or nonlocal (%d)", aloc, eatLocal, eatNonlocal);
        gmx_incons(stmp);
    }

    cu_plist_t      *plist   = cu_nb->plist[iloc];
    cu_timers_t     *timers  = cu_nb->timers;
    wallclock_gpu_t *timings = cu_nb->timings;
    nb_staging      nbst     = cu_nb->nbst;

    gmx_bool    calc_ener   = flags & GMX_FORCE_VIRIAL;
    gmx_bool    calc_fshift = flags & GMX_FORCE_VIRIAL;

    /* turn energy calculation always on/off (for debugging/testing only) */
    calc_ener = (calc_ener || always_ener) && !never_ener; 

    /* don't launch wait/update timers & counters if there was no work to do

       NOTE: if timing with multiple GPUs (streams) becomes possible, the
       counters could end up being inconsistent due to not being incremented
       on some of the nodes! */
    if (cu_nb->plist[iloc]->nsci == 0)
    {
        return;
    }

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(aloc))
    {
        adat_end = cu_nb->atdat->natoms_local;
    }
    else
    {
        adat_end = cu_nb->atdat->natoms;
    }

    if (cu_nb->use_stream_sync)
    {
        stat = cudaStreamSynchronize(cu_nb->stream[iloc]);
        CU_RET_ERR(stat, "cudaStreamSynchronize failed in cu_blockwait_nb");
    }
    else 
    {
        /* Busy-wait until we get the signalling pattern set in last 4-bytes 
           of the l/nl float vector. */
        signal_bytes    = (unsigned int*)&nbatom->out[0].f[adat_end*3 + 2];
        while (*signal_bytes != gpu_sync_signal_pattern)
        {
            /* back off, otherwise we get stuck */
            cuthread_yield();
        }
    }

    /* timing data accumulation */
    if (cu_nb->do_time)
    {
        /* only increase counter once (at local F wait) */
        if (LOCAL_I(iloc))
        {
            timings->nb_c++;
            timings->ktime[plist->do_prune ? 1 : 0][calc_ener ? 1 : 0].c += 1;
        }

        /* kernel timings */
        timings->ktime[plist->do_prune ? 1 : 0][calc_ener ? 1 : 0].t +=
            cu_event_elapsed(timers->start_nb_k[iloc], timers->stop_nb_k[iloc]);

        /* X/q H2D and F D2H timings */
        timings->nb_h2d_t += cu_event_elapsed(timers->start_nb_h2d[iloc],
                                                 timers->stop_nb_h2d[iloc]);
        timings->nb_d2h_t += cu_event_elapsed(timers->start_nb_d2h[iloc],
                                                 timers->stop_nb_d2h[iloc]);

        /* only count atdat and pair-list H2D at pair-search step */
        if (plist->do_prune)
        {
            /* atdat transfer timing (add only once, at local F wait) */
            if (LOCAL_A(aloc))
            {
                timings->pl_h2d_c++;
                timings->pl_h2d_t += cu_event_elapsed(timers->start_atdat,
                                                         timers->stop_atdat);
            }

            timings->pl_h2d_t += cu_event_elapsed(timers->start_pl_h2d[iloc],
                                                     timers->stop_pl_h2d[iloc]);
        }
    }

    /* add up enegies and shift forces (only once at local F wait) */
    if (LOCAL_I(iloc))
    {
        if (calc_ener)
        {
            *e_lj += *nbst.e_lj;
            *e_el += *nbst.e_el;
        }

        if (calc_fshift)
        {
            for (i = 0; i < SHIFTS; i++)
            {
                fshift[i][0] += nbst.fshift[i].x;
                fshift[i][1] += nbst.fshift[i].y;
                fshift[i][2] += nbst.fshift[i].z;
            }
        }
    }

    /* turn off pruning (doesn't matter if this is pair-search step or not) */
    plist->do_prune = FALSE;
}
