#include "stdlib.h"

#include "smalloc.h"
#include "force.h"
#include "types/simple.h" 
#include "types/nblist_box.h"

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
                  // gmx_bool calc_ene,
                  int flags,
                  gmx_bool sync)
{
    cu_atomdata_t   *adat = cu_nb->atomdata;
    cu_nb_params_t  *nb_params = cu_nb->nb_params;
    cu_nblist_t     *nblist = cu_nb->nblist;
    cu_timers_t     *timers = cu_nb->timers;

    int     shmem; 
    int     nb_blocks = calc_nb_blocknr(nblist->nci);
    dim3    dim_block(CELL_SIZE, CELL_SIZE, 1); 
    dim3    dim_grid(nb_blocks, 1, 1); 

    gmx_bool calc_ene   = flags & GMX_FORCE_VIRIAL;
    gmx_bool calc_fshift = flags & GMX_FORCE_VIRIAL;
    gmx_bool time_trans = timers->time_transfers; 

    p_k_calc_nb nb_kernel = NULL; /* fn pointer to the nonbonded kernel */

    static gmx_bool doKernel2 = (getenv("GMX_NB_K2") != NULL);        
    static gmx_bool doAlwaysNsPrune = (getenv("GMX_GPU_ALWAYS_NS_PRUNE") != NULL);

    /* XXX debugging code, remove it */
    calc_ene = (calc_ene || alwaysE) && !neverE; 

    if (debug)
    {
        fprintf(debug, "GPU launch configuration:\n\tThread block: %dx%dx%d\n\tGrid: %dx%d\n\t#Cells/Subcells: %d/%d (%d)\n",         
        dim_block.x, dim_block.y, dim_block.z, dim_grid.x, dim_grid.y, nblist->nci*NSUBCELL, 
        NSUBCELL, nblist->naps);
    }
    
    /* beginning of timed nonbonded calculation section */
    cudaEventRecord(timers->start_nb, 0);

    /* beginning of timed HtoD section */
    if (time_trans)
    {
        cudaEventRecord(timers->start_nb_h2d, 0);
    }

    /* 0 the force output array */
    cudaMemsetAsync(adat->f, 0, adat->natoms * sizeof(*adat->f), 0);

    /* HtoD x, q */    
    upload_cudata_async(adat->xq, nbatom->x, adat->natoms * sizeof(*adat->xq), 0);

    /* HtoD shift vec if we have a dynamic box */
    if (nbatom->dynamic_box || !adat->shift_vec_copied)
    {
        upload_cudata_async(adat->shift_vec, nbatom->shift_vec, SHIFTS * sizeof(*adat->shift_vec), 0);
        adat->shift_vec_copied = TRUE;
    }

    /* set the shift force output to 0 */
    if (calc_fshift)
    {
        cudaMemsetAsync(adat->f_shift, 0, SHIFTS * sizeof(*adat->f_shift), 0);
    }

    /* set energy outputs to 0 */
    if (calc_ene)
    {
        cudaMemsetAsync(adat->e_lj, 0, sizeof(*adat->e_lj), 0);
        cudaMemsetAsync(adat->e_el, 0, sizeof(*adat->e_el), 0);
    }

    if (time_trans)
    {
        cudaEventRecord(timers->stop_nb_h2d, 0);
    }

    /* launch async nonbonded calculations */        
    /* size of force buffers in shmem */
     shmem = !doKernel2 ?
                (1 + NSUBCELL) * CELL_SIZE * CELL_SIZE * 3 * sizeof(float) :
                CELL_SIZE * CELL_SIZE * 3 * sizeof(float);
     
    nb_kernel = select_nb_kernel(nb_params->eeltype, calc_ene, 
                                 nblist->prune_nbl || doAlwaysNsPrune, doKernel2);
    nb_kernel<<<dim_grid, dim_block, shmem, 0>>>(*adat, *nb_params, *nblist, 
                                                 calc_fshift);

    if (sync)
    {
        CU_LAUNCH_ERR_SYNC("k_calc_nb");
    }
    else
    {
        CU_LAUNCH_ERR("k_calc_nb");
    }
   
    /* beginning of timed D2H section */
    if (time_trans)
    {
        cudaEventRecord(timers->start_nb_d2h, 0);
    }

    /* DtoH f */
    download_cudata_async(nbatom->f, adat->f, adat->natoms*sizeof(*adat->f), 0);

    /* DtoH f_shift */
    if (calc_fshift)
    {
        download_cudata_async(cu_nb->tmpdata.f_shift, adat->f_shift, 
                              SHIFTS * sizeof(*cu_nb->tmpdata.f_shift), 0);
    }

    /* DtoH energies */
    if (calc_ene)
    {
        download_cudata_async(cu_nb->tmpdata.e_lj, adat->e_lj, sizeof(*cu_nb->tmpdata.e_lj), 0);
        download_cudata_async(cu_nb->tmpdata.e_el, adat->e_el, sizeof(*cu_nb->tmpdata.e_el), 0);
    }

    if (time_trans)
    {        
        cudaEventRecord(timers->stop_nb_d2h, 0);
    }

    cudaEventRecord(timers->stop_nb, 0);
}

/*! Blocking wait for the asynchrounously launched nonbonded calculations to finish. */
void cu_blockwait_nb(cu_nonbonded_t cu_nb, int flags, 
                     float *e_lj, float *e_el, rvec *fshift)
{    
    cudaError_t     s;
    int             i;
    float           t_tot, t;
    gmx_bool        calc_ene   = flags & GMX_FORCE_VIRIAL;
    gmx_bool        calc_fshift = flags & GMX_FORCE_VIRIAL;

    cu_timers_t     *timers  = cu_nb->timers;
    cu_timings_t    *timings = cu_nb->timings;
    nb_tmp_data     td = cu_nb->tmpdata;    

    cu_blockwait_event(timers->stop_nb, timers->start_nb, &t_tot);
    timings->nb_count++;
    
    if (timers->time_transfers)
    {        
        s = cudaEventElapsedTime(&t, timers->start_nb_h2d, timers->stop_nb_h2d);
        CU_RET_ERR(s, "cudaEventElapsedTime failed in cu_blockwait_nb");
        timings->nb_h2d_time += t;
        t_tot -= t;
        
        s = cudaEventElapsedTime(&t, timers->start_nb_d2h, timers->stop_nb_d2h);
        CU_RET_ERR(s, "cudaEventElapsedTime failed in cu_blockwait_nb");    
        timings->nb_d2h_time += t;
        t_tot -= t;
    }

    timings->k_time[cu_nb->nblist->prune_nbl ? 1 : 0][calc_ene ? 1 : 0].t += t_tot;
    timings->k_time[cu_nb->nblist->prune_nbl ? 1 : 0][calc_ene ? 1 : 0].c += 1;
   
    /* turn off neighborlist pruning */
    cu_nb->nblist->prune_nbl = FALSE;

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

/*! Checks if the nonbonded calculation has finished. */
gmx_bool cu_checkstat_nb(cu_nonbonded_t cu_nb, float *time)
{
    cudaError_t stat; 
    cu_timers_t *timers = cu_nb->timers;

    time = NULL;
    stat = cudaEventQuery(timers->stop_nb);

    /* we're done, let's calculate times*/
    if (stat == cudaSuccess)
    {
        stat = cudaEventElapsedTime(time, timers->start_nb, timers->stop_nb);
        CU_RET_ERR(stat, "cudaEventElapsedTime on start_nb and stop_nb failed");
    }
    else 
    {
        /* do we have an error? */
        if (stat != cudaErrorNotReady) 
        {
            CU_RET_ERR(stat, "the execution of the nonbonded calculations has failed");
        }
    }
    
    return (stat == cudaSuccess ? TRUE : FALSE);
}
