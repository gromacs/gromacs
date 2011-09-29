#include "stdlib.h"

#include "smalloc.h"
#include "force.h"
#include "types/simple.h" 
#include "types/nblist_box.h"
#include "nb_verlet.h"

#include "cutypedefs.h"
#include "cudautils.h"
#include "cuda_nb.h"
#include "cuda_data_mgmt.h"
#include "cupmalloc.h"

#define CELL_SIZE               (GPU_NS_CELL_SIZE)
#define CELL_SIZE_POW2_EXPONENT (3) /* NOTE: change this together with GPU_NS_CELL_SIZE !*/
#define NB_DEFAULT_THREADS      (CELL_SIZE * CELL_SIZE)

#include "cutype_utils.cuh"
#include "nb_kernel_utils.cuh"

/* Generate all combinations of force and energy-calculation and/or pruning kernels. */
/** Force only kernels **/
#include "nb_kernels.cuh"
/** Force & energy kernels **/
#define CALC_ENERGIES
#include "nb_kernels.cuh"
#undef CALC_ENERGIES

/*** Neighborlist pruning kernels ***/
/** Force only kernels **/
#define PRUNE_NBL
#include "nb_kernels.cuh"
/** Force & energy kernels **/
#define CALC_ENERGIES
#include "nb_kernels.cuh"
#undef CALC_ENERGIES
#undef PRUNE_NBL

/*! nonbonded kernel function pointer type */
typedef void (*p_k_calc_nb)(const cu_atomdata_t,
                            const cu_nb_params_t, 
                            const cu_nblist_t,
                            gmx_bool /*calc virial*/);

/* XXX
    if GMX_GPU_ENE env var set it always runs the energy kernel unless the 
    GMX_GPU_NO_ENE env var is set, case in which it never runs the energy kernel.     
    --> only for benchmarking purposes */
static gmx_bool alwaysE = (getenv("GMX_GPU_ALWAYS_ENE") != NULL); 
static gmx_bool neverE  = (getenv("GMX_GPU_NEVER_ENE") != NULL);

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

/*! Selects the kernel version (force / energy / pruning) to execute and 
 * returns a function pointer to it. 
 */
static inline p_k_calc_nb select_nb_kernel(int eeltype, gmx_bool doEne, 
                                           gmx_bool doPrune, gmx_bool doKernel2)
{
    p_k_calc_nb k = NULL;

    /* select which kernel will be used */
    switch (eeltype)
    {
        case cu_eelCUT:
            if (!doKernel2)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_calc_nb_cutoff_forces_1 : 
                                   k_calc_nb_cutoff_forces_prunenbl_1;                                  
                }
                else 
                {
                    k = !doPrune ? k_calc_nb_cutoff_forces_energies_1 :
                                   k_calc_nb_cutoff_forces_energies_prunenbl_1;
                }
            }
            else 
            {
                if (!doEne)
                {
                    k = !doPrune ? k_calc_nb_cutoff_forces_2 :
                                   k_calc_nb_cutoff_forces_prunenbl_2;
                }
                else 
                {
                    k = !doPrune ? k_calc_nb_cutoff_forces_energies_2 :
                                   k_calc_nb_cutoff_forces_energies_prunenbl_2;
                }
            }
            break;

        case cu_eelRF:
            if (!doKernel2)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_calc_nb_RF_forces_1 :
                                   k_calc_nb_RF_forces_prunenbl_1;
                }
                else
                {
                    k = !doPrune ? k_calc_nb_RF_forces_energies_1 :
                                   k_calc_nb_RF_forces_energies_prunenbl_1;
                }
            }
            else 
            {
                if (!doEne)
                {
                    k = !doPrune ? k_calc_nb_RF_forces_2 :
                                   k_calc_nb_RF_forces_prunenbl_2;
                }
                else
                {
                    k = !doPrune ? k_calc_nb_RF_forces_energies_2 :
                                   k_calc_nb_RF_forces_energies_prunenbl_2;
                }
            }
            break;

        case cu_eelEWALD:
            if (!doKernel2)
            {
                if (!doEne)
                {
                    k = !doPrune ? k_calc_nb_ewald_forces_1 :
                                   k_calc_nb_ewald_forces_prunenbl_1;
                }
                else
                {
                    k = !doPrune ? k_calc_nb_ewald_forces_energies_1:
                                   k_calc_nb_ewald_forces_energies_prunenbl_1;
                }
            }
            else 
            {
                if (!doEne)
                {
                    k = !doPrune ? k_calc_nb_ewald_forces_2 :
                                   k_calc_nb_ewald_forces_prunenbl_2;
                }
                else
                {
                    k = !doPrune ? k_calc_nb_ewald_forces_energies_2 :
                                   k_calc_nb_ewald_forces_energies_prunenbl_2;
                }
            }
            break;

        default: 
            gmx_incons("The provided electrostatics type does not exist in the  CUDA implementation!");
    }
    return k;
}

/*!  Launch asynchronously the nonbonded force calculations. 

    This consists of the following (async) steps launched in the default stream 0: 
   - initilize to zero force output;
   - upload x and q;
   - upload shift vector;
   - launch kernel;
   - download forces/energies.
    
    Timing is done using:
    - start_nb/stop_nb events for total execution time;
    - start_nb_h2d/stop_nb_h2d and start_nb_h2d/stop_nb_h2d event for 
    the CPU->GPU and GPU->CPU transfers, respectively.
 */
void cu_stream_nb(cu_nonbonded_t cu_nb,
                  const gmx_nb_atomdata_t *nbatom,
                  int flags,
                  int iloc)
{
    cu_atomdata_t   *adat       = cu_nb->atomdata;
    cu_nb_params_t  *nb_params  = cu_nb->nb_params;
    cu_nblist_t     *nblist     = cu_nb->nblist[iloc];
    cu_timers_t     *timers     = cu_nb->timers;
    cudaStream_t    stream      = timers->nbstream[iloc];

    cudaEvent_t start_nb   = timers->start_nb[iloc];
    cudaEvent_t stop_nb    = timers->stop_nb[iloc];
    cudaEvent_t start_h2d  = timers->start_nb_h2d[iloc];
    cudaEvent_t stop_h2d   = timers->stop_nb_h2d[iloc];

    int adat_begin, adat_len;  /* local/nonlocal offset and length used for xq and f */

    /* kernel lunch stuff */
    p_k_calc_nb nb_kernel = NULL; /* fn pointer to the nonbonded kernel */
    int         shmem;
    int         nb_blocks = calc_nb_blocknr(nblist->nci);
    dim3        dim_block(CELL_SIZE, CELL_SIZE, 1);
    dim3        dim_grid(nb_blocks, 1, 1);

    gmx_bool calc_ene    = flags & GMX_FORCE_VIRIAL;
    gmx_bool calc_fshift = flags & GMX_FORCE_VIRIAL;
    gmx_bool time_trans  = timers->time_transfers;

    static gmx_bool doKernel2 = (getenv("GMX_NB_K2") != NULL);        
    static gmx_bool doAlwaysNsPrune = (getenv("GMX_GPU_ALWAYS_NS_PRUNE") != NULL);

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

    /* XXX debugging code, remove it */
    calc_ene = (calc_ene || alwaysE) && !neverE; 

    if (debug)
    {
        fprintf(debug, "GPU launch configuration:\n\tThread block: %dx%dx%d\n\tGrid: %dx%d\n\t#Cells/Subcells: %d/%d (%d)\n",         
        dim_block.x, dim_block.y, dim_block.z, dim_grid.x, dim_grid.y, nblist->nci*NSUBCELL, 
        NSUBCELL, nblist->naps);
    }
    
    /* wait for the atomdata trasfer to be finished */
    cu_synchstream_atomdata(cu_nb, iloc);
    // FIXME fix this mess - all this waiting is not necessary!
    cudaError_t stat = cudaStreamWaitEvent(stream, timers->stop_clear, 0);
    CU_RET_ERR(stat, "cudaEventElapsedTime on stop_clear failed in cu_blockwait_nb");

    /* beginning of timed HtoD section */
    if (time_trans)
    {
        cudaEventRecord(start_h2d, stream);
    }

    /* HtoD x, q */
    upload_cudata_async(adat->xq + adat_begin, nbatom->x + adat_begin * 4,
                        adat_len * sizeof(*adat->xq), stream); 

    if (time_trans)
    {
        cudaEventRecord(stop_h2d, stream);
    }

    /* beginning of timed nonbonded calculation section */
    cudaEventRecord(start_nb, stream);

    /* launch async nonbonded calculations */        
    /* size of force buffers in shmem */
     shmem = !doKernel2 ?
                (1 + NSUBCELL) * CELL_SIZE * CELL_SIZE * 3 * sizeof(float) :
                CELL_SIZE * CELL_SIZE * 3 * sizeof(float);
     
    nb_kernel = select_nb_kernel(nb_params->eeltype, calc_ene, 
                                 nblist->prune_nbl || doAlwaysNsPrune, doKernel2);
    nb_kernel<<<dim_grid, dim_block, shmem, stream>>>(*adat, *nb_params, *nblist, 
                                                 calc_fshift);

    CU_LAUNCH_ERR("k_calc_nb");

    cudaEventRecord(stop_nb, stream);
}

void cu_copyback_nb_data(cu_nonbonded_t cu_nb,
                  const gmx_nb_atomdata_t *nbatom,
                  int flags,
                  int aloc)
{
    int adat_begin, adat_len;  /* local/nonlocal offset and length used for xq and f */
    int iloc = -1;

    gmx_bool calc_ene    = flags & GMX_FORCE_VIRIAL;
    gmx_bool calc_fshift = flags & GMX_FORCE_VIRIAL;

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
        printf("aloc = %d\n", aloc);
        gmx_incons("Invalid atom locality passed (valid here is only local or nonlocal)");
    }

    cu_atomdata_t   *adat   = cu_nb->atomdata;
    cu_timers_t     *timers = cu_nb->timers;
    cudaStream_t    stream  = timers->nbstream[iloc];
    cudaEvent_t     start   = timers->start_nb_d2h[iloc];
    cudaEvent_t     stop    = timers->stop_nb_d2h[iloc];
    gmx_bool        time_trans  = timers->time_transfers;

    if (LOCAL_A(aloc))
    {
        adat_begin  = 0;
        adat_len    = adat->natoms_local;
    }
    else
    {
        adat_begin  = adat->natoms_local;
        adat_len    = adat->natoms - adat->natoms_local;
    }

    /* beginning of timed D2H section */
    if (time_trans)
    {
        cudaEventRecord(start, stream);
    }

    /* DtoH f */
    download_cudata_async(nbatom->out[0].f + adat_begin * 4, adat->f + adat_begin, 
                          adat_len * sizeof(*adat->f), stream);

    /* only transfer energies in the local stream */
    if (LOCAL_I(iloc))
    {
        /* DtoH f_shift */
        if (calc_fshift)
        {
            download_cudata_async(cu_nb->tmpdata.f_shift, adat->f_shift,
                    SHIFTS * sizeof(*cu_nb->tmpdata.f_shift), stream);
        }

        /* DtoH energies */
        if (calc_ene)
        {
            download_cudata_async(cu_nb->tmpdata.e_lj, adat->e_lj, sizeof(*cu_nb->tmpdata.e_lj), stream);
            download_cudata_async(cu_nb->tmpdata.e_el, adat->e_el, sizeof(*cu_nb->tmpdata.e_el), stream);
        }
    }

    if (time_trans)
    {
        cudaEventRecord(stop, stream);
    }
}

/*! Blocking wait for the asynchrounously launched nonbonded calculations to finish. */
void cu_blockwait_nb(cu_nonbonded_t cu_nb,
                     int flags,
                     int aloc,
                     float *e_lj, float *e_el, rvec *fshift)
{
    int i, iloc = -1;
    cudaError_t s;
    gmx_bool    calc_ene   = flags & GMX_FORCE_VIRIAL;
    gmx_bool    calc_fshift = flags & GMX_FORCE_VIRIAL;

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

    cu_nblist_t     *nblist  = cu_nb->nblist[iloc];
    cu_timers_t     *timers  = cu_nb->timers;
    cu_timings_t    *timings = cu_nb->timings;
    nb_tmp_data     td       = cu_nb->tmpdata;
    cudaStream_t    stream   = timers->nbstream[iloc];

    cudaEvent_t start_nb        = timers->start_nb[iloc];
    cudaEvent_t stop_nb         = timers->stop_nb[iloc];
    cudaEvent_t start_h2d       = timers->start_nb_h2d[iloc];
    cudaEvent_t stop_h2d        = timers->stop_nb_h2d[iloc];
    cudaEvent_t start_d2h       = timers->start_nb_d2h[iloc];
    cudaEvent_t stop_d2h        = timers->stop_nb_d2h[iloc];
    cudaEvent_t start_nbl_h2d   = timers->start_nbl_h2d[iloc];
    cudaEvent_t stop_nbl_h2d    = timers->stop_nbl_h2d[iloc];

    s = cudaStreamSynchronize(stream);
    CU_RET_ERR(s, "cudaStreamSynchronize failed in cu_blockwait_nb");

    /* only increase counter once (at local F wait) */
    if (LOCAL_I(iloc))
    {
        timings->nb_count++;
        timings->k_time[nblist->prune_nbl ? 1 : 0][calc_ene ? 1 : 0].c += 1;
    }
    
    /* accumulate kernel timings */
    timings->k_time[nblist->prune_nbl ? 1 : 0][calc_ene ? 1 : 0].t +=
        cu_event_elapsed(start_nb, stop_nb);

    // FIXME remove the time_transfers conditional 
    if (timers->time_transfers)
    {
        /* X/q H2D and F D2H timings */
        timings->nb_h2d_time += cu_event_elapsed(start_h2d, stop_h2d);
        timings->nb_d2h_time += cu_event_elapsed(start_d2h, stop_d2h);

        /* only count atomdata and nbl H2D if it's neighbor search step */
        if (nblist->prune_nbl)
        {
            /* only add atomdata transfer time once (at local F wait) */
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
   
    /* turn off neighborlist pruning */
    nblist->prune_nbl = FALSE;

    /* add up enegies and shift forces (only once at local F wait) */
    if (LOCAL_I(iloc))
    {
        /* XXX debugging code, remove this */
        calc_ene = (calc_ene || alwaysE) && !neverE; 

        if (calc_ene)
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
