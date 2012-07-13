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
#include <assert.h>

#if defined(_MSVC)
#include <limits>
#endif

#include "smalloc.h"
#include "types/simple.h" 
#include "types/nbnxn_pairlist.h"
#include "types/nb_verlet.h"
#include "types/ishift.h"
#include "types/force_flags.h"

#ifdef TMPI_ATOMICS
#include "thread_mpi/atomic.h"
#endif

#include "nbnxn_cuda_types.h"
#include "../../gmxlib/cuda_tools/cudautils.cuh"
#include "nbnxn_cuda.h"
#include "nbnxn_cuda_data_mgmt.h"
#include "pmalloc_cuda.h"


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
typedef void (*nbnxn_cu_k_ptr_t)(const cu_atomdata_t,
                                 const cu_nbparam_t,
                                 const cu_plist_t,
                                 gmx_bool);

/*********************************/

/* XXX always/never run the energy/pruning kernels -- only for benchmarking purposes */
static gmx_bool always_ener  = (getenv("GMX_GPU_ALWAYS_ENER") != NULL);
static gmx_bool never_ener   = (getenv("GMX_GPU_NEVER_ENER") != NULL);
static gmx_bool always_prune = (getenv("GMX_GPU_ALWAYS_PRUNE") != NULL);


/* Bit-pattern used for polling-based GPU synchronization. It is used as a float
 * and corresponds to having the exponent set to the maximum (127 -- single
 * precision) and the matissa to 0.
 */
static unsigned int poll_wait_pattern = (0x7FU << 23);

/*! Returns the number of blocks to be used for the nonbonded GPU kernel. */
static inline int calc_nb_kernel_nblock(int nwork_units, cu_dev_info_t *dinfo)
{
    int max_grid_x_size;

    assert(dinfo);

    max_grid_x_size = dinfo->dev_prop.maxGridSize[0];

    /* do we exceed the grid x dimension limit? */
    if (nwork_units > max_grid_x_size)
    {
        gmx_fatal(FARGS, "Watch out system too large to simulate!\n"
                  "The number of nonbonded work units (=number of super-clusters) exceeds the"
                  "maximum grid size in x dimension (%d > %d)!", nwork_units, max_grid_x_size);
    }

    return nwork_units;
}


/* Constant arrays listing all kernel function pointers and enabling selection
   of a kernel in an elegant manner. */

static const int nEnergyKernelTypes = 2; /* 0 - no nergy, 1 - energy */
static const int nPruneKernelTypes  = 2; /* 0 - no prune, 1 - prune */

/* Default kernels */
static const nbnxn_cu_k_ptr_t
nb_cu_k_default[eelCuNR][nEnergyKernelTypes][nPruneKernelTypes] =
{
    { { k_nbnxn_ewald,              k_nbnxn_ewald_prune },
      { k_nbnxn_ewald_ener,         k_nbnxn_ewald_ener_prune } },
    { { k_nbnxn_rf,                 k_nbnxn_rf_prune },
      { k_nbnxn_rf_ener,            k_nbnxn_rf_ener_prune } },
    { { k_nbnxn_ewald,              k_nbnxn_ewald_prune },
      { k_nbnxn_cutoff_ener,        k_nbnxn_cutoff_ener_prune } },
};

/* Legacy kernels */
static const nbnxn_cu_k_ptr_t
nb_cu_k_legacy[eelCuNR][nEnergyKernelTypes][nPruneKernelTypes] =
{
    { { k_nbnxn_ewald_legacy,       k_nbnxn_ewald_prune_legacy },
      { k_nbnxn_ewald_ener_legacy,  k_nbnxn_ewald_ener_prune_legacy } },
    { { k_nbnxn_rf_legacy,          k_nbnxn_rf_prune_legacy },
      { k_nbnxn_rf_ener_legacy,     k_nbnxn_rf_ener_prune_legacy } },
    { { k_nbnxn_ewald_legacy,       k_nbnxn_ewald_prune_legacy },
      { k_nbnxn_cutoff_ener_legacy, k_nbnxn_cutoff_ener_prune_legacy } },
};

/* Old kernels */
static const nbnxn_cu_k_ptr_t
nb_cu_k_old[eelCuNR][nEnergyKernelTypes][nPruneKernelTypes] =
{
    { { k_nbnxn_ewald_old,          k_nbnxn_ewald_prune_old },
      { k_nbnxn_ewald_ener_old,     k_nbnxn_ewald_ener_prune_old } },
    { { k_nbnxn_rf_old,             k_nbnxn_rf_prune_old },
      { k_nbnxn_rf_ener_old,        k_nbnxn_rf_ener_prune_old } },
    { { k_nbnxn_ewald_old,          k_nbnxn_ewald_prune_old },
      { k_nbnxn_cutoff_ener_old,    k_nbnxn_cutoff_ener_prune_old } },
};

/*! Return a pointer to the kernel version to be executed at the current step. */
static inline nbnxn_cu_k_ptr_t select_nbnxn_kernel(int kver, int eeltype,
                                                    gmx_bool doEne, gmx_bool doPrune)
{
    assert(kver < eNbnxnCuKNR);
    assert(eeltype < eelCuNR);

    if (NBNXN_KVER_OLD(kver))
    {
        return nb_cu_k_old[eeltype][doEne][doPrune];
    }
    else if (NBNXN_KVER_LEGACY(kver))
    {
        return nb_cu_k_legacy[eeltype][doEne][doPrune];
    }
    else
    {
        return nb_cu_k_default[eeltype][doEne][doPrune];
    }
}

/*! Calculates the amount of shared memory required for kernel version in use. */
static inline int calc_shmem_required(int kver)
{
    int shmem;

    /* size of shmem (force-buffers/xq/atom type preloading) */
    if (NBNXN_KVER_OLD(kver))
    {
        shmem = (1 + NSUBCELL) * CLUSTER_SIZE * CLUSTER_SIZE * 3 * sizeof(float);
    }
    else if (NBNXN_KVER_LEGACY(kver))
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

    return shmem;
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
    /* CUDA kernel launch-related stuff */
    int  shmem, nblock;
    dim3 dim_block, dim_grid;
    nbnxn_cu_k_ptr_t nb_kernel = NULL; /* fn pointer to the nonbonded kernel */

    cu_atomdata_t   *adat   = cu_nb->atdat;
    cu_nbparam_t    *nbp    = cu_nb->nbparam;
    cu_plist_t      *plist  = cu_nb->plist[iloc];
    cu_timers_t     *t      = cu_nb->timers;
    cudaStream_t    stream  = cu_nb->stream[iloc];

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

    /* get the pointer to the kernel flavor we need to use */
    nb_kernel = select_nbnxn_kernel(cu_nb->kernel_ver, nbp->eeltype, calc_ener,
                                    plist->do_prune || always_prune);

    /* kernel launch config */
    nblock    = calc_nb_kernel_nblock(plist->nsci, cu_nb->dev_info);
    dim_block = dim3(CLUSTER_SIZE, CLUSTER_SIZE, 1);
    dim_grid  = dim3(nblock, 1, 1);
    shmem     = calc_shmem_required(cu_nb->kernel_ver);

    if (debug)
    {
        fprintf(debug, "GPU launch configuration:\n\tThread block: %dx%dx%d\n\t"
                "Grid: %dx%d\n\t#Cells/Subcells: %d/%d (%d)\n",
                dim_block.x, dim_block.y, dim_block.z,
                dim_grid.x, dim_grid.y, plist->nsci*NSUBCELL,
                NSUBCELL, plist->na_c);
    }

    nb_kernel<<<dim_grid, dim_block, shmem, stream>>>(*adat, *nbp, *plist, calc_fshift);
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

        /* Set the last four bytes of the force array to a bit pattern
           which can't be the result of the force calculation:
           max exponent (127) and zero mantissa. */
        *(unsigned int*)&nbatom->out[0].f[adat_end*3 - 1] = poll_wait_pattern;
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
                      (adat_len)*sizeof(*adat->f), stream);

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

/* Atomic compare-exchange operation on unsigned values. It is used in
 * polling wait for the GPU.
 */
static inline bool atomic_cas(volatile unsigned int *ptr,
                              unsigned int oldval,
                              unsigned int newval)
{
    assert(ptr);

#ifdef TMPI_ATOMICS
    return tMPI_Atomic_cas((tMPI_Atomic_t *)ptr, oldval, newval);
#else
    gmx_incons("Atomic operations not available, atomic_cas() should not have been called!");
    return TRUE;
#endif
}

void nbnxn_cuda_wait_gpu(nbnxn_cuda_ptr_t cu_nb,
                         const nbnxn_atomdata_t *nbatom,
                         int flags, int aloc,
                         float *e_lj, float *e_el, rvec *fshift)
{
    cudaError_t stat;
    int i, adat_end, iloc = -1;
    volatile unsigned int *poll_word;

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
        /* Busy-wait until we get the signalling pattern set in last byte
         * of the l/nl float vector. This pattern corresponds to a floating
         * point number which can't be the result of the force calculation
         * (maximum, 127 exponent and 0 mantissa).
         * The polling uses atomic compare-exchange.
         */
        poll_word = (volatile unsigned int*)&nbatom->out[0].f[adat_end*3 - 1];
        while (atomic_cas(poll_word, poll_wait_pattern, poll_wait_pattern)) {}
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
