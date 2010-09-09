#include "stdlib.h"

#include "smalloc.h"

#include "types/simple.h" 
#include "types/nblist_box.h"
#include "cutypedefs.h"
#include "cudautils.h"

#include "gpu_nb.h"
#include "gpu_data.h"
#include "cupmalloc.h"

#define CELL_SIZE               (GPU_NS_CELL_SIZE)
#define CELL_SIZE_POW2_EXPONENT (3) /* NOTE: change this togather with GPU_NS_CELL_SIZE !*/
#define NB_DEFAULT_THREADS      (CELL_SIZE * CELL_SIZE)

#include "gpu_nb_kernel_utils.h"

/*****  Generate kernels  *****/

/* Force-type defines (EL_*) and FUNCTION_NAME are undef-d at the end 
   of the kernels in gpu_nb_kernels.h
*/

/*** Force only kernels ***/
/* Cut-Off */
#define EL_CUTOFF
#define FUNCTION_NAME(x, y) x##_cutoff_##y
#include "gpu_nb_kernels.h"
/* Reaction-Field */
#define EL_RF
#define FUNCTION_NAME(x, y) x##_RF_##y
#include "gpu_nb_kernels.h"
/* Ewald */
#define EL_EWALD
#define FUNCTION_NAME(x, y) x##_ewald_##y
#include "gpu_nb_kernels.h"

/*** Force & energy kernels ***/
#define CALC_ENERGIES

/* Cut-Off */
#define EL_CUTOFF
#define FUNCTION_NAME(x, y) x##_cutoff_##y
#include "gpu_nb_kernels.h"
/* Reaction-Field */
#define EL_RF
#define FUNCTION_NAME(x, y) x##_RF_##y
#include "gpu_nb_kernels.h"
/* Ewald */
#define EL_EWALD
#define FUNCTION_NAME(x, y) x##_ewald_##y
#include "gpu_nb_kernels.h"

#undef CALC_ENERGIES

/* TODO clean this up! */
/****** Prune neighborlist ******/
#define PRUNE_NBL
/*** Force only kernels ***/
/* Cut-Off */
#define EL_CUTOFF
#define FUNCTION_NAME(x, y) x##_cutoff_##y
#include "gpu_nb_kernels.h"
/* Reaction-Field */
#define EL_RF
#define FUNCTION_NAME(x, y) x##_RF_##y
#include "gpu_nb_kernels.h"
/* Ewald */
#define EL_EWALD
#define FUNCTION_NAME(x, y) x##_ewald_##y
#include "gpu_nb_kernels.h"

/*** Force & energy kernels ***/
#define CALC_ENERGIES

/* Cut-Off */
#define EL_CUTOFF
#define FUNCTION_NAME(x, y) x##_cutoff_##y
#include "gpu_nb_kernels.h"
/* Reaction-Field */
#define EL_RF
#define FUNCTION_NAME(x, y) x##_RF_##y
#include "gpu_nb_kernels.h"
/* Ewald */
#define EL_EWALD
#define FUNCTION_NAME(x, y) x##_ewald_##y
#include "gpu_nb_kernels.h"

#undef CALC_ENERGIES
#undef PRUNE_NBL


/* XXX
    if GMX_GPU_ENE env var set it always runs the energy kernel unless the 
    GMX_GPU_NO_ENE env var is set, case in which it never runs the energy kernel.     
    --> only for benchmarking purposes */
static gmx_bool alwaysE = (getenv("GMX_GPU_ALWAYS_ENE") != NULL); 
static gmx_bool neverE  = (getenv("GMX_GPU_NEVER_ENE") != NULL);

/* based on the number of work units, return the number of blocks to be used 
   for the nonbonded GPU kernel */
inline int calc_nb_blocknr(int nwork_units)
{
    int retval = (nwork_units <= GRID_MAX_DIM ? nwork_units : GRID_MAX_DIM);
    if (retval != nwork_units)
    {
        gmx_fatal(FARGS, "Watch out, the number of nonbonded work units exceeds the maximum grid size (%d > %d)!",
                nwork_units, GRID_MAX_DIM);
    }
    return retval;
}

/*  Launch asynchronously the nonbonded force calculations. 

    This consists of the following (async) steps launched in the default stream 0: 
   - initilize to zero force output
   - upload x and q
   - upload shift vector
   - launch kernel
   - download forces
    
    Timing is done using the start_nb and stop_nb events.
 */
void cu_stream_nb(t_cudata d_data,
                  const gmx_nb_atomdata_t *nbatom,                                    
                  gmx_bool calc_ene,
                  gmx_bool sync)
{
    int     shmem; 
    int     nb_blocks = calc_nb_blocknr(d_data->nci);
    dim3    dim_block(CELL_SIZE, CELL_SIZE, 1); 
    dim3    dim_grid(nb_blocks, 1, 1); 
    gmx_bool time_trans = d_data->time_transfers; 

    /* fn pointer to the nonbonded kernel */
    /* force-only */
    void    (*p_k_calc_nb_f)(
                const gmx_nbl_ci_t * /*nbl_ci*/,
                const gmx_nbl_sj4_t * /*nbl_sj4*/,
                const gmx_nbl_excl_t * /*excl*/,
                const int * /*atom_types*/,
                int /*ntypes*/,
                const float4 * /*xq*/,
                const float * /*nbfp*/,
                const float3 * /*shift_vec*/,
                float /*two_k_rf*/,
                float /*cutoff_sq*/,
                float /*coulomb_tab_scale*/,
                float4 * /*f*/) = NULL;
    /* force & energy */
    void    (*p_k_calc_nb_fe)(
                const gmx_nbl_ci_t * /*nbl_ci*/,
                const gmx_nbl_sj4_t * /*nbl_sj4*/,
                const gmx_nbl_excl_t * /*excl*/,
                const int * /*atom_types*/,
                int /*ntypes*/,
                const float4 * /*xq*/,
                const float * /*nbfp*/,
                const float3 * /*shift_vec*/,
                float /*two_k_rf*/,
                float /*cutoff_sq*/,
                float /*coulomb_tab_scale*/,                
                float /*beta*/,
                float /*c_rf*/,
                float * /*e_lj*/,
                float * /*e_el*/,
                float4 * /*f*/) = NULL;

    /* TODO clean this up! */
    /****** Prune neighborlist ******/
    /* force-only */
    void    (*p_k_calc_nb_f_pnbl)(
                const gmx_nbl_ci_t * /*nbl_ci*/,
                gmx_nbl_sj4_t * /*nbl_sj4*/,
                const gmx_nbl_excl_t * /*excl*/,
                const int * /*atom_types*/,
                int /*ntypes*/,
                const float4 * /*xq*/,
                const float * /*nbfp*/,
                const float3 * /*shift_vec*/,
                float /*two_k_rf*/,
                float /*cutoff_sq*/,
                float /*coulomb_tab_scale*/,
                float /*rlist_sq*/,
                float4 * /*f*/) = NULL;
    /* force & energy */
    void    (*p_k_calc_nb_fe_pnbl)(
                const gmx_nbl_ci_t * /*nbl_ci*/,
                gmx_nbl_sj4_t * /*nbl_sj4*/,
                const gmx_nbl_excl_t * /*excl*/,
                const int * /*atom_types*/,
                int /*ntypes*/,
                const float4 * /*xq*/,
                const float * /*nbfp*/,
                const float3 * /*shift_vec*/,
                float /*two_k_rf*/,
                float /*cutoff_sq*/,
                float /*coulomb_tab_scale*/,                
                float /*rlist_sq*/,
                float /*beta*/,
                float /*c_rf*/,
                float * /*e_lj*/,
                float * /*e_el*/,
                float4 * /*f*/) = NULL;


    static gmx_bool doKernel2 = (getenv("GMX_NB_K2") != NULL);        
    static gmx_bool doAlwaysNsPrune = (getenv("GMX_GPU_ALWAYS_NS_PRUNE") != NULL);
    
    /* XXX debugging code, remove this */
    calc_ene = (calc_ene || alwaysE) && !neverE; 

    /* size of force buffers in shmem */
    if (!doKernel2)
    {
        shmem =  (1 + NSUBCELL) * CELL_SIZE * CELL_SIZE * 3 * sizeof(float);
    }
    else 
    {    
        shmem =  CELL_SIZE * CELL_SIZE * 3 * sizeof(float);
    }

    /* select which kernel will be used */
    switch (d_data->eeltype)
    {
        case cu_eelCUT:
            if (!doKernel2)
            {
                p_k_calc_nb_f   = k_calc_nb_cutoff_forces_1n;
                p_k_calc_nb_fe  = k_calc_nb_cutoff_forces_energies_1n;
                p_k_calc_nb_f_pnbl  = k_calc_nb_cutoff_forces_prunenbl_1n;
                p_k_calc_nb_fe_pnbl = k_calc_nb_cutoff_forces_energies_prunenbl_1n;
            }
            else 
            {
                p_k_calc_nb_f   = k_calc_nb_cutoff_forces_2n;
                p_k_calc_nb_fe  = k_calc_nb_cutoff_forces_energies_2n;
                p_k_calc_nb_f_pnbl  = k_calc_nb_cutoff_forces_prunenbl_2n;
                p_k_calc_nb_fe_pnbl = k_calc_nb_cutoff_forces_energies_prunenbl_2n;
            }
            break;
        case cu_eelRF:
            if (!doKernel2)
            {
                p_k_calc_nb_f   = k_calc_nb_RF_forces_1n;
                p_k_calc_nb_fe  = k_calc_nb_RF_forces_energies_1n;
                p_k_calc_nb_f_pnbl  = k_calc_nb_RF_forces_prunenbl_1n;
                p_k_calc_nb_fe_pnbl = k_calc_nb_RF_forces_energies_prunenbl_1n;
            }
            else 
            {
                p_k_calc_nb_f   = k_calc_nb_RF_forces_2n;
                p_k_calc_nb_fe  = k_calc_nb_RF_forces_energies_2n;
                p_k_calc_nb_f_pnbl  = k_calc_nb_RF_forces_prunenbl_2n;
                p_k_calc_nb_fe_pnbl = k_calc_nb_RF_forces_energies_prunenbl_2n;
            }
            break;
        case cu_eelEWALD:
            if (!doKernel2)
            {
                p_k_calc_nb_f   = k_calc_nb_ewald_forces_1n;
                p_k_calc_nb_fe  = k_calc_nb_ewald_forces_energies_1n;
                p_k_calc_nb_f_pnbl  = k_calc_nb_ewald_forces_prunenbl_1n;
                p_k_calc_nb_fe_pnbl = k_calc_nb_ewald_forces_energies_prunenbl_1n;
            }
            else 
            {
                p_k_calc_nb_f   = k_calc_nb_ewald_forces_2n;
                p_k_calc_nb_fe  = k_calc_nb_ewald_forces_energies_2n; 
                p_k_calc_nb_f_pnbl  = k_calc_nb_ewald_forces_prunenbl_2n;
                p_k_calc_nb_fe_pnbl = k_calc_nb_ewald_forces_energies_prunenbl_2n;
            }
            break;
        default: 
        {
            gmx_incons("The provided electrostatics type does not exist in the  CUDA implementation!");
        }
    }   

    if (debug)
    {
        fprintf(debug, "GPU launch configuration:\n\tThread block: %dx%dx%d\n\tGrid: %dx%d\n\t#Cells/Subcells: %d/%d (%d)\n",         
        dim_block.x, dim_block.y, dim_block.z, dim_grid.x, dim_grid.y, d_data->nci*NSUBCELL, 
        NSUBCELL, d_data->naps);
    }
    
    /* beginning of timed nonbonded calculation section */
    cudaEventRecord(d_data->start_nb, 0);

    /* beginning of timed HtoD section */
    if (time_trans)
    {
        cudaEventRecord(d_data->start_nb_h2d, 0);
    }

    /* 0 the force output array */
    cudaMemsetAsync(d_data->f, 0, d_data->natoms * sizeof(*d_data->f), 0);

    /* HtoD x, q */    
    upload_cudata_async(d_data->xq, nbatom->x, d_data->natoms * sizeof(*d_data->xq), 0);

    /* HtoD shift vec if we have a dynamic box */
    if (nbatom->dynamic_box || !d_data->shift_vec_copied)
    {
        upload_cudata_async(d_data->shift_vec, nbatom->shift_vec, SHIFTS * sizeof(*d_data->shift_vec), 0);   
        d_data->shift_vec_copied = TRUE;
    }
    
    if (time_trans)
    {
        cudaEventRecord(d_data->stop_nb_h2d, 0);
    }

    /* launch async nonbonded calculations */        
    if (!calc_ene)
    {
        if (!(d_data->prune_nbl || doAlwaysNsPrune))
        {
            p_k_calc_nb_f<<<dim_grid, dim_block, shmem, 0>>>(
                    d_data->ci,
                    d_data->sj4,
                    d_data->excl,
                    d_data->atom_types,
                    d_data->ntypes,
                    d_data->xq,
                    d_data->nbfp,
                    d_data->shift_vec,
                    d_data->two_k_rf,
                    d_data->cutoff_sq,
                    d_data->coulomb_tab_scale,
                    d_data->f);
        } 
        else 
        {
            p_k_calc_nb_f_pnbl<<<dim_grid, dim_block, shmem, 0>>>(
                    d_data->ci,
                    d_data->sj4,
                    d_data->excl,
                    d_data->atom_types,
                    d_data->ntypes,
                    d_data->xq,
                    d_data->nbfp,
                    d_data->shift_vec,
                    d_data->two_k_rf,
                    d_data->cutoff_sq,
                    d_data->coulomb_tab_scale,
                    d_data->rlist_sq,
                    d_data->f);

        }
    } 
    else 
    {
        /* set energy outputs to 0 */
        cudaMemsetAsync(d_data->e_lj, 0.0f, sizeof(*d_data->e_lj), 0);
        cudaMemsetAsync(d_data->e_el, 0.0f, sizeof(*d_data->e_el), 0);
        if (!(d_data->prune_nbl /*|| doAlwaysNsPrune */))
        {  
            p_k_calc_nb_fe<<<dim_grid, dim_block, shmem, 0>>>(
                    d_data->ci,
                    d_data->sj4,
                    d_data->excl,
                    d_data->atom_types,
                    d_data->ntypes,
                    d_data->xq,
                    d_data->nbfp,
                    d_data->shift_vec,
                    d_data->two_k_rf,
                    d_data->cutoff_sq,
                    d_data->coulomb_tab_scale,
                    d_data->ewald_beta,
                    d_data->c_rf,
                    d_data->e_lj, 
                    d_data->e_el,
                    d_data->f);       
        }
        else 
        {
            p_k_calc_nb_fe_pnbl<<<dim_grid, dim_block, shmem, 0>>>(
                    d_data->ci,
                    d_data->sj4,
                    d_data->excl,
                    d_data->atom_types,
                    d_data->ntypes,
                    d_data->xq,
                    d_data->nbfp,
                    d_data->shift_vec,
                    d_data->two_k_rf,
                    d_data->cutoff_sq,
                    d_data->coulomb_tab_scale,
                    d_data->rlist_sq,
                    d_data->ewald_beta,
                    d_data->c_rf,
                    d_data->e_lj, 
                    d_data->e_el,
                    d_data->f);

        }
    }

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
        cudaEventRecord(d_data->start_nb_d2h, 0);
    }

    /* DtoH f */
    download_cudata_async(nbatom->f, d_data->f, d_data->natoms*sizeof(*d_data->f), 0);
    /* DtoH energies */
    if (calc_ene)
    {
        download_cudata_async(d_data->tmpdata.e_lj, d_data->e_lj, sizeof(*d_data->e_lj), 0);
        download_cudata_async(d_data->tmpdata.e_el, d_data->e_el, sizeof(*d_data->e_el), 0);
    }

    if (time_trans)
    {        
        cudaEventRecord(d_data->stop_nb_d2h, 0);
    }

    cudaEventRecord(d_data->stop_nb, 0);

    /* turn off neighborlist pruning */
    // d_data->prune_nbl = FALSE;
}

/* Blocking wait for the asynchrounously launched nonbonded calculations to finish. */
void cu_blockwait_nb(t_cudata d_data, gmx_bool calc_ene, 
                     float *e_lj, float *e_el)
{    
    cudaError_t s;
    float t_tot, t;

    cu_blockwait_event(d_data->stop_nb, d_data->start_nb, &t_tot);
    d_data->timings.nb_count++;
    
    if (d_data->time_transfers)
    {        
        s = cudaEventElapsedTime(&t, d_data->start_nb_h2d, d_data->stop_nb_h2d);
        CU_RET_ERR(s, "cudaEventElapsedTime failed in cu_blockwait_nb");
        d_data->timings.nb_h2d_time += t;
        t_tot -= t;
        
        s = cudaEventElapsedTime(&t, d_data->start_nb_d2h, d_data->stop_nb_d2h);
        CU_RET_ERR(s, "cudaEventElapsedTime failed in cu_blockwait_nb");    
        d_data->timings.nb_d2h_time += t;
        t_tot -= t;
    }

    d_data->timings.k_time[d_data->prune_nbl ? 1 : 0][calc_ene ? 1 : 0].t += t_tot;
    d_data->timings.k_time[d_data->prune_nbl ? 1 : 0][calc_ene ? 1 : 0].c += 1;
   
    /* turn off neighborlist pruning */
    d_data->prune_nbl = FALSE;

    /* XXX debugging code, remove this */
    calc_ene = (calc_ene || alwaysE) && !neverE; 

    if (calc_ene)
    {
        *e_lj += *d_data->tmpdata.e_lj;
        *e_el += *d_data->tmpdata.e_el;
    }
}

/* Blocking wait for the asynchrounously launched nonbonded calculations to finish. */
void cu_blockwait_nb_OLD(t_cudata d_data, float *time)
{    
    cudaError_t stat;

    // stat = cudaStreamSynchronize(d_data->nb_stream);    
    stat = cudaEventSynchronize(d_data->stop_nb);
    CU_RET_ERR(stat, "the async execution of nonbonded calculations has failed"); 

    stat = cudaEventElapsedTime(time, d_data->start_nb, d_data->stop_nb);
    CU_RET_ERR(stat, "cudaEventElapsedTime on start_nb and stop_nb failed");
}

/* Check if the nonbonded calculation has finished. */
gmx_bool cu_checkstat_nb(t_cudata d_data, float *time)
{
    cudaError_t stat; 
    
    time = NULL;
    stat = cudaEventQuery(d_data->stop_nb);

    /* we're done, let's calculate times*/
    if (stat == cudaSuccess)
    {
        stat = cudaEventElapsedTime(time, d_data->start_nb, d_data->stop_nb);
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


/* XXX:  remove, not used anyomore! */
void cu_do_nb(t_cudata d_data, rvec shift_vec[]) 
{
#if 0
    int     nb_blocks = calc_nb_blocknr(d_data->nci);
    dim3    dim_block(CELL_SIZE, CELL_SIZE, 1); 
    dim3    dim_grid(nb_blocks, 1, 1); 
    int     shmem = (1 + NSUBCELL) * CELL_SIZE * CELL_SIZE * sizeof(float4); /* force buffer */

    if (debug)
    {
        printf("~> Thread block: %dx%dx%d\n~> Grid: %dx%d\n~> #SubCell pairs: %d (%d)\n", 
            dim_block.x, dim_block.y, dim_block.z, dim_grid.x, dim_grid.y, d_data->nsi, 
            d_data->naps);
    }

    /* set the forces to 0 */
    cudaMemset(d_data->f, 0, d_data->natoms*sizeof(*d_data->f));

    /* upload shift vec */
    upload_cudata(d_data->shift_vec, shift_vec, SHIFTS*sizeof(*d_data->shift_vec));   

    /* sync nonbonded calculations */      
    k_calc_nb_1<<<dim_grid, dim_block, shmem>>>(d_data->ci,
                                                  d_data->sj, 
                                                  d_data->si,
                                                  d_data->atom_types, 
                                                  d_data->ntypes, 
                                                  d_data->xq, 
                                                  d_data->nbfp,
                                                  d_data->shift_vec,
                                                  d_data->ewald_beta,
                                                  d_data->f);
    CU_LAUNCH_ERR_SYNC("k_calc_nb");
#endif 
}
