#include "stdlib.h"

#include "smalloc.h"
#include "force.h"
#include "types/simple.h" 
#include "types/nbnxn_pairlist.h"
#include "types/nb_verlet.h"

#include "nbnxn_cuda_types.h"
#include "cudautils.cuh"
#include "nbnxn_cuda.h"
#include "nbnxn_cuda_data_mgmt.h"
#include "pmalloc_cuda.h"

#define CLUSTER_SIZE            (GPU_NS_CLUSTER_SIZE)

#include "vectype_ops.cuh"
#include "nbnxn_cuda_kernel_utils.cuh"

/* Generate all combinations of kernels through multiple inclusion:
   F, F + E, F + prune, F + E + prune. */
/** Force only **/
#include "nbnxn_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "nbnxn_cuda_kernels.cuh"
#undef CALC_ENERGIES

/*** Neighborlist pruning kernels ***/
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
                          const cu_nb_params_t, 
                          const cu_nblist_t,
                          gmx_bool /*calc virial*/);

/* XXX
    if GMX_GPU_ENE env var set it always runs the energy kernel unless the 
    GMX_GPU_NO_ENE env var is set, case in which it never runs the energy kernel.     
    --> only for benchmarking purposes */
static gmx_bool always_ener  = (getenv("GMX_GPU_ALWAYS_ENER") != NULL);
static gmx_bool never_ener   = (getenv("GMX_GPU_NEVER_ENER") != NULL);
static gmx_bool always_prune = (getenv("GMX_GPU_ALWAYS_PRUNE") != NULL);

/*! Returns the number of blocks to be used  for the nonbonded GPU kernel. */
static inline int calc_nb_blocknr(int nwork_units)
{
    int retval = (nwork_units <= GRID_MAX_DIM ? nwork_units : GRID_MAX_DIM);
    if (retval != nwork_units)
    {
        gmx_fatal(FARGS, "Watch out, the number of nonbonded work units exceeds the maximum grid size (%d > %d)!",
                nwork_units, GRID_MAX_DIM);
    }
    return retval;
}

/*! Selects the kernel version to execute at the current step and 
 *  returns a function pointer to it. 
 */
static inline p_k_nbnxn select_nbnxn_kernel(int eeltype, gmx_bool doEne, 
                                            gmx_bool doPrune, gmx_bool doKernel2)
{
    p_k_nbnxn k = NULL;

    /* select which kernel will be used */
    switch (eeltype)
    {
        case cu_eelCUT:
            if (!doKernel2)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_cutoff_1 :
                                   k_nbnxn_cutoff_prune_1;
                }
                else 
                {
                    k = !doPrune ? k_nbnxn_cutoff_ener_1 :
                                   k_nbnxn_cutoff_ener_prune_1;
                }
            }
            else 
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_cutoff_2 :
                                   k_nbnxn_cutoff_prune_2;
                }
                else 
                {
                    k = !doPrune ? k_nbnxn_cutoff_ener_2 :
                                   k_nbnxn_cutoff_ener_prune_2;
                }
            }
            break;

        case cu_eelRF:
            if (!doKernel2)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_rf_1 :
                                   k_nbnxn_rf_prune_1;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_rf_ener_1 :
                                   k_nbnxn_rf_ener_prune_1;
                }
            }
            else 
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_rf_2 :
                                   k_nbnxn_rf_prune_2;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_rf_ener_2 :
                                   k_nbnxn_rf_ener_prune_2;
                }
            }
            break;

        case cu_eelEWALD:
            if (!doKernel2)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_ewald_1 :
                                   k_nbnxn_ewald_prune_1;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_ewald_ener_1:
                                   k_nbnxn_ewald_ener_prune_1;
                }
            }
            else 
            {
                if (!doEne)
                {
                    k = !doPrune ? k_nbnxn_ewald_2 :
                                   k_nbnxn_ewald_prune_2;
                }
                else
                {
                    k = !doPrune ? k_nbnxn_ewald_ener_2 :
                                   k_nbnxn_ewald_ener_prune_2;
                }
            }
            break;

        default: 
            gmx_incons("The provided electrostatics type does not exist in the CUDA implementation!");
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
void nbnxn_cuda_launch_kernel(cu_nonbonded_t cu_nb,
                              const nbnxn_atomdata_t *nbatom,
                              int flags,
                              int iloc)
{
    cudaError_t stat;
    int adat_begin, adat_len;  /* local/nonlocal offset and length used for xq and f */

    cu_atomdata_t   *adat   = cu_nb->atomdata;
    cu_nb_params_t  *nbp    = cu_nb->nb_params;
    cu_nblist_t     *nblist = cu_nb->nblist[iloc];
    cu_timers_t     *timers = cu_nb->timers;
    cudaStream_t    stream  = cu_nb->stream[iloc];

    cudaEvent_t start_nb_k  = timers->start_nb_k[iloc];
    cudaEvent_t stop_nb_k   = timers->stop_nb_k[iloc];
    cudaEvent_t start_h2d   = timers->start_nb_h2d[iloc];
    cudaEvent_t stop_h2d    = timers->stop_nb_h2d[iloc];

    /* kernel lunch stuff */
    p_k_nbnxn   nb_kernel = NULL; /* fn pointer to the nonbonded kernel */
    int         shmem;
    int         nb_blocks = calc_nb_blocknr(nblist->nsci);
    dim3        dim_block(CLUSTER_SIZE, CLUSTER_SIZE, 1);
    dim3        dim_grid(nb_blocks, 1, 1);

    gmx_bool calc_ener   = flags & GMX_FORCE_VIRIAL;
    gmx_bool calc_fshift = flags & GMX_FORCE_VIRIAL;
    gmx_bool do_time     = cu_nb->do_time;

    /* turn energy calculation always on/off (for debugging/testing only) */
    calc_ener = (calc_ener || always_ener) && !never_ener; 

    static gmx_bool doKernel2 = (getenv("GMX_NB_K2") != NULL);        

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
                dim_grid.x, dim_grid.y, nblist->nsci*NSUBCELL,
                NSUBCELL, nblist->na_c);
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
        stat = cudaEventRecord(start_h2d, stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* HtoD x, q */
    cu_copy_H2D_async(adat->xq + adat_begin, nbatom->x + adat_begin * 4,
                      adat_len * sizeof(*adat->xq), stream); 

    if (do_time)
    {
        stat = cudaEventRecord(stop_h2d, stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* beginning of timed nonbonded calculation section */
    if (do_time)
    {
        stat = cudaEventRecord(start_nb_k, stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* size of force buffers in shmem */
    shmem = !doKernel2 ?
        (1 + NSUBCELL) * CLUSTER_SIZE * CLUSTER_SIZE * 3 * sizeof(float) :
        CLUSTER_SIZE * CLUSTER_SIZE * 3 * sizeof(float);

    /* get the pointer to the kernel we need & launch it */
    nb_kernel = select_nbnxn_kernel(nbp->eeltype, calc_ener,
                                    nblist->prune_nbl || always_prune,
                                    doKernel2);
    nb_kernel<<<dim_grid, dim_block, shmem, stream>>>(*adat, *nbp, *nblist,
                                                      calc_fshift);
    CU_LAUNCH_ERR("k_calc_nb");

    if (do_time)
    {
        stat = cudaEventRecord(stop_nb_k, stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

void nbnxn_cuda_launch_cpyback(cu_nonbonded_t cu_nb,
                               const nbnxn_atomdata_t *nbatom,
                               int flags,
                               int aloc)
{
    cudaError_t stat;
    int adat_begin, adat_len, adat_last;  /* local/nonlocal offset and length used for xq and f */
    int iloc = -1;
    float4 *signal_field;

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
        printf("aloc = %d\n", aloc);
        gmx_incons("Invalid atom locality passed (valid here is only local or nonlocal)");
    }

    cu_atomdata_t   *adat   = cu_nb->atomdata;
    cu_timers_t     *t      = cu_nb->timers;
    gmx_bool        do_time = cu_nb->do_time;
    cudaStream_t    stream  = cu_nb->stream[iloc];

    gmx_bool calc_ener   = flags & GMX_FORCE_VIRIAL;
    gmx_bool calc_fshift = flags & GMX_FORCE_VIRIAL;

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(aloc))
    {
        adat_begin  = 0;
        adat_len    = adat->natoms_local;
        adat_last   = cu_nb->atomdata->natoms_local - 1;
    }
    else
    {
        adat_begin  = adat->natoms_local;
        adat_len    = adat->natoms - adat->natoms_local;
        adat_last   = cu_nb->atomdata->natoms - 1;
    }

    /* beginning of timed D2H section */
    if (do_time)
    {
        stat = cudaEventRecord(t->start_nb_d2h[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    if (!cu_nb->use_stream_sync)
    {
        /* XXX: could move this to cu_clear_nb_f_out() 
             -> it's just that nbnxn_atomdata_t is not available there */

        /* Set the 4th unsused component of the last element in the GPU force array
           (separately for l/nl parts) to a specific bit-pattern that we can check 
           for while waiting for the transfer to be done. */
        signal_field = adat->f + adat_last;
        stat = cudaMemsetAsync(&signal_field->w, 0xAA, sizeof(signal_field->z), stream);
        CU_RET_ERR(stat, "cudaMemsetAsync of the signal_field failed");

        /* Clear the bits in the force array (on the CPU) that we are going to poll. */
        nbatom->out[0].f[adat_last*4 + 3] = 0;
    }

    /* With DD the local D2H transfer can only start after the non-local 
       has been launched. */
    if (iloc == eintLocal && cu_nb->dd_run)
    {
        stat = cudaStreamWaitEvent(stream, cu_nb->nonlocal_done, 0);
        CU_RET_ERR(stat, "cudaStreamWaitEvent on nonlocal_done failed");
    }

    /* DtoH f */
    cu_copy_D2H_async(nbatom->out[0].f + adat_begin * 4, adat->f + adat_begin, 
                      adat_len * sizeof(*adat->f), stream);

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
        /* DtoH f_shift */
        if (calc_fshift)
        {
            cu_copy_D2H_async(cu_nb->tmpdata.f_shift, adat->f_shift,
                              SHIFTS * sizeof(*cu_nb->tmpdata.f_shift), stream);
        }

        /* DtoH energies */
        if (calc_ener)
        {
            cu_copy_D2H_async(cu_nb->tmpdata.e_lj, adat->e_lj, 
                              sizeof(*cu_nb->tmpdata.e_lj), stream);
            cu_copy_D2H_async(cu_nb->tmpdata.e_el, adat->e_el, 
                              sizeof(*cu_nb->tmpdata.e_el), stream);
        }
    }

    if (do_time)
    {
        stat = cudaEventRecord(t->stop_nb_d2h[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

void nbnxn_cuda_wait_gpu(cu_nonbonded_t cu_nb,
                         const nbnxn_atomdata_t *nbatom,
                         int flags, int aloc,
                         float *e_lj, float *e_el, rvec *fshift)
{
    cudaError_t stat;
    int i, adat_last, iloc = -1;
    unsigned int *signal_bytes;
    unsigned int signal_pattern;

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
        gmx_incons("Invalid atom locality passed (valid here is only local or nonlocal)");
    }

    cu_nblist_t     *nblist  = cu_nb->nblist[iloc];
    cu_timers_t     *timers  = cu_nb->timers;
    cu_timings_t    *timings = cu_nb->timings;
    nb_tmp_data     td       = cu_nb->tmpdata;

    cudaEvent_t start_nb_k      = timers->start_nb_k[iloc];
    cudaEvent_t stop_nb_k       = timers->stop_nb_k[iloc];
    cudaEvent_t start_h2d       = timers->start_nb_h2d[iloc];
    cudaEvent_t stop_h2d        = timers->stop_nb_h2d[iloc];
    cudaEvent_t start_d2h       = timers->start_nb_d2h[iloc];
    cudaEvent_t stop_d2h        = timers->stop_nb_d2h[iloc];
    cudaEvent_t start_nbl_h2d   = timers->start_nbl_h2d[iloc];
    cudaEvent_t stop_nbl_h2d    = timers->stop_nbl_h2d[iloc];

    gmx_bool    calc_ener   = flags & GMX_FORCE_VIRIAL;
    gmx_bool    calc_fshift = flags & GMX_FORCE_VIRIAL;

    /* turn energy calculation always on/off (for debugging/testing only) */
    calc_ener = (calc_ener || always_ener) && !never_ener; 

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(aloc))
    {
        adat_last = cu_nb->atomdata->natoms_local - 1;
    }
    else
    {
        adat_last = cu_nb->atomdata->natoms - 1;
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
        signal_bytes    = (unsigned int*)&nbatom->out[0].f[adat_last*4 + 3];
        signal_pattern  = 0xAAAAAAAA;
        while (*signal_bytes != signal_pattern) 
        {
            /* back off, otherwise we get stuck. why??? */
            usleep(0); 
        }
    }

    /* timing data accumulation */
    if (cu_nb->do_time)
    {
        /* only increase counter once (at local F wait) */
        if (LOCAL_I(iloc))
        {
            timings->nb_count++;
            timings->k_time[nblist->prune_nbl ? 1 : 0][calc_ener ? 1 : 0].c += 1;
        }

        /* kernel timings */
        timings->k_time[nblist->prune_nbl ? 1 : 0][calc_ener ? 1 : 0].t +=
            cu_event_elapsed(start_nb_k, stop_nb_k);

        /* X/q H2D and F D2H timings */
        timings->nb_h2d_time += cu_event_elapsed(start_h2d, stop_h2d);
        timings->nb_d2h_time += cu_event_elapsed(start_d2h, stop_d2h);

        /* only count atomdata and nbl H2D at neighbor search step */
        if (nblist->prune_nbl)
        {
            /* atomdata transfer timing (add only once, at local F wait) */
            if (LOCAL_A(aloc))
            {
                timings->nbl_h2d_count++;
                timings->nbl_h2d_time +=
                    cu_event_elapsed(timers->start_atdat, timers->stop_atdat);
            }

            timings->nbl_h2d_time +=
                cu_event_elapsed(start_nbl_h2d, stop_nbl_h2d);
        }
    }
   
    /* turn off pruning (doesn't matter if this is neighbor search step or not) */
    nblist->prune_nbl = FALSE;

    /* add up enegies and shift forces (only once at local F wait) */
    if (LOCAL_I(iloc))
    {
        if (calc_ener)
        {
            *e_lj += *td.e_lj;
            *e_el += *td.e_el;
        }

        if (calc_fshift)
        {
            for (i = 0; i < SHIFTS; i++)
            {
                fshift[i][0] += td.f_shift[i].x;
                fshift[i][1] += td.f_shift[i].y;
                fshift[i][2] += td.f_shift[i].z;
            }
        }
    }
}
